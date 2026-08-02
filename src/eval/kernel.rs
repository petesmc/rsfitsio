//! Arithmetic and comparison kernels over [`ColumnarValue`].
//!
//! Each kernel takes two values and returns one. Compare with the old engine:
//! `Do_BinOp_lng`, `Do_BinOp_dbl`, `Do_BinOp_log`, `Do_BinOp_str` and
//! `Do_BinOp_bit` were five near-identical 300-line functions, chosen by which
//! `DoOp` the parser installed, each opening with the same scalar-versus-vector
//! preamble and each looping with raw pointers into a shared arena.
//!
//! Here the sort dispatch happens once, in [`arith`] and [`compare`], and the
//! scalar case is a variant of the input type rather than a branch inside every
//! loop.

use super::value::{Array, ArrayData, ColumnarValue, Scalar, ValueError};
use crate::c_types::c_long;
use crate::eval_defs::ValueSort;

type Res = Result<ColumnarValue, ValueError>;

/// The arithmetic operators, in the sense `eval.y` gave them.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Arith {
    Add,
    Sub,
    Mul,
    Div,
    /// C's `%`: truncated, so the sign follows the dividend.
    Mod,
    Pow,
}

/// The comparison operators.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Compare {
    Eq,
    Ne,
    Gt,
    Lt,
    Gte,
    Lte,
    /// `~`, equal to within `APPROX`.
    Approx,
}

/// `eval.y`'s tolerance for `~`.
const APPROX: f64 = 1.0e-7;

/// The sort an arithmetic result takes: the wider of the two operands, with
/// booleans promoted to integers.
///
/// This is `PROMOTE` from `eval.y`, which worked because the sort tags were
/// ordered; [`ValueSort`] derives `Ord` over the same order.
fn promote(a: ValueSort, b: ValueSort) -> ValueSort {
    match a.max(b) {
        ValueSort::Boolean => ValueSort::Long,
        other => other,
    }
}

/// Apply `f` element by element, broadcasting scalars.
///
/// A result element is null when either input element is null; `f` is not
/// called for those, which is what stops a division by a null propagating a
/// spurious floating-point exception.
fn zip_numeric(
    lhs: &ColumnarValue,
    rhs: &ColumnarValue,
    out: ValueSort,
    op: &'static str,
    f: impl Fn(f64, f64) -> Option<f64>,
) -> Res {
    let len = lhs.broadcast_len(rhs)?;
    let l = lhs.numeric(op)?;
    let r = rhs.numeric(op)?;

    let Some(n) = len else {
        /* both sides are scalars, so the result is one too */
        let (lv, ln) = l.get(0);
        let (rv, rn) = r.get(0);
        if ln || rn {
            return Ok(ColumnarValue::Null(out));
        }
        return Ok(match f(lv, rv) {
            Some(v) => ColumnarValue::Scalar(scalar_of(out, v)),
            None => ColumnarValue::Null(out),
        });
    };

    let mut nulls = vec![false; n];
    let mut any_null = false;
    match out {
        ValueSort::Double => {
            let mut data = vec![0.0f64; n];
            for i in 0..n {
                let (lv, ln) = l.get(i);
                let (rv, rn) = r.get(i);
                match if ln || rn { None } else { f(lv, rv) } {
                    Some(v) => data[i] = v,
                    None => {
                        nulls[i] = true;
                        any_null = true;
                    }
                }
            }
            Ok(array(ArrayData::Double(data), nulls, any_null))
        }
        _ => {
            let mut data = vec![0 as c_long; n];
            for i in 0..n {
                let (lv, ln) = l.get(i);
                let (rv, rn) = r.get(i);
                match if ln || rn { None } else { f(lv, rv) } {
                    Some(v) => data[i] = v as c_long,
                    None => {
                        nulls[i] = true;
                        any_null = true;
                    }
                }
            }
            Ok(array(ArrayData::Long(data), nulls, any_null))
        }
    }
}

fn array(data: ArrayData, nulls: Vec<bool>, any_null: bool) -> ColumnarValue {
    ColumnarValue::Array(if any_null {
        Array::with_nulls(data, nulls)
    } else {
        Array::new(data)
    })
}

fn scalar_of(sort: ValueSort, v: f64) -> Scalar {
    match sort {
        ValueSort::Double => Scalar::Double(v),
        _ => Scalar::Long(v as c_long),
    }
}

/// `lhs OP rhs` for the arithmetic operators.
pub(crate) fn arith(op: Arith, lhs: &ColumnarValue, rhs: &ColumnarValue) -> Res {
    let (lt, rt) = (lhs.sort(), rhs.sort());
    if !numeric_sort(lt) || !numeric_sort(rt) {
        return Err(ValueError::Incompatible(name(op), lt, rt));
    }
    let out = promote(lt, rt);

    match op {
        Arith::Add => zip_numeric(lhs, rhs, out, "+", |a, b| Some(a + b)),
        Arith::Sub => zip_numeric(lhs, rhs, out, "-", |a, b| Some(a - b)),
        Arith::Mul => zip_numeric(lhs, rhs, out, "*", |a, b| Some(a * b)),
        /* integer division by zero yields a null rather than a trap, which is
        what the engine did via its undef flag */
        Arith::Div => zip_numeric(lhs, rhs, out, "/", move |a, b| {
            if b == 0.0 && out != ValueSort::Double {
                None
            } else {
                Some(a / b)
            }
        }),
        Arith::Mod => zip_numeric(lhs, rhs, out, "%", move |a, b| {
            if b == 0.0 {
                None
            } else if out == ValueSort::Double {
                Some(a % b)
            } else {
                Some(((a as c_long) % (b as c_long)) as f64)
            }
        }),
        Arith::Pow => zip_numeric(lhs, rhs, out, "**", |a, b| Some(a.powf(b))),
    }
}

fn numeric_sort(t: ValueSort) -> bool {
    t.is_expr() || t == ValueSort::Boolean
}

fn name(op: Arith) -> &'static str {
    match op {
        Arith::Add => "+",
        Arith::Sub => "-",
        Arith::Mul => "*",
        Arith::Div => "/",
        Arith::Mod => "%",
        Arith::Pow => "**",
    }
}

/// `lhs OP rhs` for the comparison operators, over numeric operands.
pub(crate) fn compare(op: Compare, lhs: &ColumnarValue, rhs: &ColumnarValue) -> Res {
    let (lt, rt) = (lhs.sort(), rhs.sort());
    if !numeric_sort(lt) || !numeric_sort(rt) {
        return Err(ValueError::Incompatible("comparison", lt, rt));
    }
    let len = lhs.broadcast_len(rhs)?;
    let l = lhs.numeric("comparison")?;
    let r = rhs.numeric("comparison")?;

    let test = |a: f64, b: f64| match op {
        Compare::Eq => a == b,
        Compare::Ne => a != b,
        Compare::Gt => a > b,
        Compare::Lt => a < b,
        Compare::Gte => a >= b,
        Compare::Lte => a <= b,
        /* the engine's approximate equality: within APPROX of the larger
        magnitude, so the tolerance scales with the values */
        Compare::Approx => (a - b).abs() < APPROX * a.abs().max(b.abs()).max(f64::MIN_POSITIVE),
    };

    let Some(n) = len else {
        let (lv, ln) = l.get(0);
        let (rv, rn) = r.get(0);
        return Ok(if ln || rn {
            ColumnarValue::Null(ValueSort::Boolean)
        } else {
            ColumnarValue::Scalar(Scalar::Boolean(test(lv, rv)))
        });
    };

    let mut data = vec![false; n];
    let mut nulls = vec![false; n];
    let mut any_null = false;
    for i in 0..n {
        let (lv, ln) = l.get(i);
        let (rv, rn) = r.get(i);
        if ln || rn {
            nulls[i] = true;
            any_null = true;
        } else {
            data[i] = test(lv, rv);
        }
    }
    Ok(array(ArrayData::Boolean(data), nulls, any_null))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn s_long(v: c_long) -> ColumnarValue {
        ColumnarValue::Scalar(Scalar::Long(v))
    }
    fn s_dbl(v: f64) -> ColumnarValue {
        ColumnarValue::Scalar(Scalar::Double(v))
    }
    fn a_long(v: &[c_long]) -> ColumnarValue {
        ColumnarValue::Array(Array::new(ArrayData::Long(v.to_vec())))
    }
    fn a_dbl(v: &[f64]) -> ColumnarValue {
        ColumnarValue::Array(Array::new(ArrayData::Double(v.to_vec())))
    }

    fn longs(v: &ColumnarValue) -> Vec<c_long> {
        match v {
            ColumnarValue::Array(a) => match a.data() {
                ArrayData::Long(d) => d.clone(),
                other => panic!("not a long array: {other:?}"),
            },
            other => panic!("not an array: {other:?}"),
        }
    }
    fn bools(v: &ColumnarValue) -> Vec<bool> {
        match v {
            ColumnarValue::Array(a) => match a.data() {
                ArrayData::Boolean(d) => d.clone(),
                other => panic!("not a boolean array: {other:?}"),
            },
            other => panic!("not an array: {other:?}"),
        }
    }

    #[test]
    fn scalar_plus_scalar_stays_a_scalar() {
        let r = arith(Arith::Add, &s_long(2), &s_long(3)).unwrap();
        assert_eq!(r, ColumnarValue::Scalar(Scalar::Long(5)));
        assert!(r.is_scalar());
    }

    #[test]
    fn a_scalar_broadcasts_over_an_array() {
        let r = arith(Arith::Mul, &a_long(&[1, 2, 3]), &s_long(10)).unwrap();
        assert_eq!(longs(&r), vec![10, 20, 30]);
        /* and the other way round */
        let r = arith(Arith::Mul, &s_long(10), &a_long(&[1, 2, 3])).unwrap();
        assert_eq!(longs(&r), vec![10, 20, 30]);
    }

    #[test]
    fn promotion_follows_the_sort_order() {
        /* long + double is double */
        let r = arith(Arith::Add, &s_long(1), &s_dbl(0.5)).unwrap();
        assert_eq!(r, ColumnarValue::Scalar(Scalar::Double(1.5)));
        /* boolean + long is long, not boolean */
        let r = arith(
            Arith::Add,
            &ColumnarValue::Scalar(Scalar::Boolean(true)),
            &s_long(1),
        )
        .unwrap();
        assert_eq!(r, ColumnarValue::Scalar(Scalar::Long(2)));
    }

    #[test]
    fn integer_division_and_modulo_by_zero_give_null() {
        let r = arith(Arith::Div, &a_long(&[6, 6]), &a_long(&[3, 0])).unwrap();
        assert_eq!(longs(&r)[0], 2);
        match &r {
            ColumnarValue::Array(a) => {
                assert!(!a.is_null(0));
                assert!(a.is_null(1), "6/0 must be null");
            }
            other => panic!("{other:?}"),
        }
        let r = arith(Arith::Div, &s_long(1), &s_long(0)).unwrap();
        assert_eq!(r, ColumnarValue::Null(ValueSort::Long));
    }

    #[test]
    fn float_division_by_zero_is_infinity_not_null() {
        let r = arith(Arith::Div, &s_dbl(1.0), &s_dbl(0.0)).unwrap();
        assert_eq!(r, ColumnarValue::Scalar(Scalar::Double(f64::INFINITY)));
    }

    #[test]
    fn modulo_follows_c_truncation() {
        /* C's % takes the sign of the dividend: -7 % 3 == -1 */
        let r = arith(Arith::Mod, &s_long(-7), &s_long(3)).unwrap();
        assert_eq!(r, ColumnarValue::Scalar(Scalar::Long(-1)));
        let r = arith(Arith::Mod, &s_long(7), &s_long(3)).unwrap();
        assert_eq!(r, ColumnarValue::Scalar(Scalar::Long(1)));
    }

    #[test]
    fn nulls_propagate_through_arithmetic() {
        let lhs = ColumnarValue::Array(Array::with_nulls(
            ArrayData::Long(vec![1, 2, 3]),
            vec![false, true, false],
        ));
        let r = arith(Arith::Add, &lhs, &s_long(10)).unwrap();
        assert_eq!(longs(&r)[0], 11);
        assert_eq!(longs(&r)[2], 13);
        match &r {
            ColumnarValue::Array(a) => {
                assert!(!a.is_null(0));
                assert!(a.is_null(1), "null in, null out");
                assert!(!a.is_null(2));
            }
            other => panic!("{other:?}"),
        }
    }

    #[test]
    fn a_null_scalar_poisons_the_whole_result() {
        let r = arith(
            Arith::Add,
            &a_long(&[1, 2, 3]),
            &ColumnarValue::Null(ValueSort::Long),
        )
        .unwrap();
        match &r {
            ColumnarValue::Array(a) => assert!((0..3).all(|i| a.is_null(i))),
            other => panic!("{other:?}"),
        }
    }

    #[test]
    fn comparisons_yield_booleans() {
        let r = compare(Compare::Gt, &a_long(&[1, 5, 9]), &s_long(4)).unwrap();
        assert_eq!(bools(&r), vec![false, true, true]);
        let r = compare(Compare::Eq, &s_long(3), &s_long(3)).unwrap();
        assert_eq!(r, ColumnarValue::Scalar(Scalar::Boolean(true)));
    }

    #[test]
    fn approximate_equality_scales_with_magnitude() {
        /* the engine's ~ is a relative tolerance, so a big pair compares equal
        at a difference that a small pair does not */
        let close = compare(Compare::Approx, &s_dbl(1.0e9), &s_dbl(1.0e9 + 1.0)).unwrap();
        assert_eq!(close, ColumnarValue::Scalar(Scalar::Boolean(true)));
        let far = compare(Compare::Approx, &s_dbl(1.0), &s_dbl(2.0)).unwrap();
        assert_eq!(far, ColumnarValue::Scalar(Scalar::Boolean(false)));
    }

    #[test]
    fn mismatched_array_lengths_are_an_error() {
        let e = arith(Arith::Add, &a_long(&[1, 2]), &a_long(&[1, 2, 3]));
        assert_eq!(e, Err(ValueError::LengthMismatch(2, 3)));
    }

    #[test]
    fn strings_are_not_arithmetic_operands() {
        let s = ColumnarValue::Scalar(Scalar::Str(b"a".to_vec()));
        assert_eq!(
            arith(Arith::Add, &s, &s_long(1)),
            Err(ValueError::Incompatible(
                "+",
                ValueSort::String,
                ValueSort::Long
            ))
        );
    }

    #[test]
    fn power_matches_the_engines_right_associative_use() {
        let r = arith(Arith::Pow, &s_long(2), &s_long(3)).unwrap();
        assert_eq!(r, ColumnarValue::Scalar(Scalar::Long(8)));
        let r = arith(Arith::Pow, &a_dbl(&[4.0, 9.0]), &s_dbl(0.5)).unwrap();
        match &r {
            ColumnarValue::Array(a) => match a.data() {
                ArrayData::Double(d) => assert_eq!(d, &vec![2.0, 3.0]),
                other => panic!("{other:?}"),
            },
            other => panic!("{other:?}"),
        }
    }
}
