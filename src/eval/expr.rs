//! The expression tree that replaces the `Node` arena.
//!
//! A node in the arena owned its output buffer and reached its children by
//! index, which is what forced raw pointers: `&mut Nodes[this]` together with
//! `&Nodes[child]` is not something the borrow checker will allow. Here a
//! parent owns its children and evaluation *returns* a value, so the two are
//! never borrowed at once and nothing needs a pointer.
//!
//! Evaluation is per batch of rows, matching how `ffiter` drives the engine.

use super::kernel::{self, Arith, Compare};
use super::value::{Array, ArrayData, ColumnarValue, Scalar, ValueError};
use crate::c_types::{c_char, c_long};
use crate::eval_defs::{ParseData, ValueSort};

/// A unary operator.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Unary {
    /// Arithmetic negation.
    Neg,
    /// Logical negation, defined for booleans only.
    Not,
    /// `(int)`
    ToLong,
    /// `(float)` / `(double)`
    ToDouble,
    /// An implicit widening inserted by the parser's `PROMOTE`.
    ToBoolean,
}

/// Boolean connectives, which are not comparisons and short-circuit nothing —
/// the engine evaluates both sides because it works a column at a time.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Logic {
    And,
    Or,
    Eq,
    Ne,
}

#[derive(Clone, Debug, PartialEq)]
pub(crate) enum Expr {
    Literal(Scalar),
    /// A constant that is undefined for every row: `#NULL`.
    Null(ValueSort),
    /// A column, by index into `ParseData::varData`.
    Column(usize),
    /// The one-based row number: `#ROW`.
    RowNumber,
    Unary {
        op: Unary,
        arg: Box<Expr>,
    },
    Arith {
        op: Arith,
        lhs: Box<Expr>,
        rhs: Box<Expr>,
    },
    Compare {
        op: Compare,
        lhs: Box<Expr>,
        rhs: Box<Expr>,
    },
    Logic {
        op: Logic,
        lhs: Box<Expr>,
        rhs: Box<Expr>,
    },
    /// `cond ? then : els`, over numeric or boolean branches.
    IfThenElse {
        cond: Box<Expr>,
        then: Box<Expr>,
        els: Box<Expr>,
    },
}

/// One column's data for the rows in the current batch.
#[derive(Clone, Debug)]
pub(crate) struct ColumnBatch {
    pub(crate) data: ArrayData,
    pub(crate) nulls: Vec<bool>,
    /// Elements per row, so a vector column can be reshaped.
    pub(crate) nelem: c_long,
}

/// The rows an [`Expr`] is evaluated over.
pub(crate) struct Batch {
    pub(crate) columns: Vec<ColumnBatch>,
    pub(crate) n_rows: c_long,
    /// Table row number of the first row in the batch, one-based.
    pub(crate) first_row: c_long,
}

impl Batch {
    /// Copy the current chunk out of `ParseData::varData`.
    ///
    /// This is the one place that touches the iterator's raw buffers. It copies
    /// rather than borrows: a borrowed view would have to carry a lifetime
    /// through every kernel, and the copy is a memcpy per column per chunk
    /// against per-element work downstream. Making it a borrow is a later
    /// optimisation, not a correctness question.
    pub(crate) fn gather(lParse: &ParseData, first_row: c_long, n_rows: c_long) -> Batch {
        let row_offset = first_row - lParse.firstDataRow;
        let mut columns = Vec::with_capacity(lParse.nCols as usize);

        for var in lParse.varData.iter().take(lParse.nCols as usize) {
            let count = (var.nelem * n_rows).max(0) as usize;
            let offset = (var.nelem * row_offset).max(0) as usize;

            let data = unsafe {
                match var.dtype {
                    ValueSort::Long => ArrayData::Long(
                        core::slice::from_raw_parts(var.data.cast::<c_long>().add(offset), count)
                            .to_vec(),
                    ),
                    ValueSort::Double => ArrayData::Double(
                        core::slice::from_raw_parts(var.data.cast::<f64>().add(offset), count)
                            .to_vec(),
                    ),
                    ValueSort::Boolean => ArrayData::Boolean(
                        core::slice::from_raw_parts(var.data.cast::<c_char>().add(offset), count)
                            .iter()
                            .map(|&b| b != 0)
                            .collect(),
                    ),
                    /* strings and bit strings are not ported yet; lowering
                    refuses any expression that reaches for one */
                    _ => ArrayData::Long(Vec::new()),
                }
            };

            let nulls = match &var.undef {
                Some(u) if u.len() >= offset + count => {
                    u[offset..offset + count].iter().map(|&v| v != 0).collect()
                }
                _ => Vec::new(),
            };

            columns.push(ColumnBatch {
                data,
                nulls,
                nelem: var.nelem,
            });
        }

        Batch {
            columns,
            n_rows,
            first_row,
        }
    }
}

impl Expr {
    /// The sort this expression produces, known without evaluating it.
    pub(crate) fn sort(&self, batch_sorts: &dyn Fn(usize) -> ValueSort) -> ValueSort {
        match self {
            Expr::Literal(s) => s.sort(),
            Expr::Null(t) => *t,
            Expr::Column(i) => batch_sorts(*i),
            Expr::RowNumber => ValueSort::Long,
            Expr::Unary { op, arg } => match op {
                Unary::Not => ValueSort::Boolean,
                Unary::ToLong => ValueSort::Long,
                Unary::ToDouble => ValueSort::Double,
                Unary::ToBoolean => ValueSort::Boolean,
                Unary::Neg => arg.sort(batch_sorts),
            },
            Expr::Arith { lhs, rhs, .. } => {
                let (a, b) = (lhs.sort(batch_sorts), rhs.sort(batch_sorts));
                match a.max(b) {
                    ValueSort::Boolean => ValueSort::Long,
                    other => other,
                }
            }
            Expr::Compare { .. } | Expr::Logic { .. } => ValueSort::Boolean,
            Expr::IfThenElse { then, els, .. } => then.sort(batch_sorts).max(els.sort(batch_sorts)),
        }
    }

    pub(crate) fn evaluate(&self, batch: &Batch) -> Result<ColumnarValue, ValueError> {
        match self {
            Expr::Literal(s) => Ok(ColumnarValue::Scalar(s.clone())),
            Expr::Null(t) => Ok(ColumnarValue::Null(*t)),

            Expr::Column(i) => {
                let col = &batch.columns[*i];
                Ok(ColumnarValue::Array(Array::with_nulls(
                    col.data.clone(),
                    col.nulls.clone(),
                )))
            }

            Expr::RowNumber => {
                let n = batch.n_rows.max(0) as usize;
                let rows: Vec<c_long> = (0..n).map(|i| batch.first_row + i as c_long).collect();
                Ok(ColumnarValue::Array(Array::new(ArrayData::Long(rows))))
            }

            Expr::Unary { op, arg } => {
                let v = arg.evaluate(batch)?;
                unary(*op, &v)
            }

            Expr::Arith { op, lhs, rhs } => {
                let l = lhs.evaluate(batch)?;
                let r = rhs.evaluate(batch)?;
                kernel::arith(*op, &l, &r)
            }

            Expr::Compare { op, lhs, rhs } => {
                let l = lhs.evaluate(batch)?;
                let r = rhs.evaluate(batch)?;
                kernel::compare(*op, &l, &r)
            }

            Expr::Logic { op, lhs, rhs } => {
                let l = lhs.evaluate(batch)?;
                let r = rhs.evaluate(batch)?;
                logic(*op, &l, &r)
            }

            Expr::IfThenElse { cond, then, els } => {
                let c = cond.evaluate(batch)?;
                let t = then.evaluate(batch)?;
                let e = els.evaluate(batch)?;
                if_then_else(&c, &t, &e)
            }
        }
    }
}

fn unary(op: Unary, v: &ColumnarValue) -> Result<ColumnarValue, ValueError> {
    match op {
        Unary::Neg => kernel::arith(Arith::Sub, &ColumnarValue::Scalar(Scalar::Long(0)), v),
        Unary::Not => {
            if v.sort() != ValueSort::Boolean {
                return Err(ValueError::BadSort("!", v.sort()));
            }
            /* !x is x == false */
            logic(Logic::Eq, v, &ColumnarValue::Scalar(Scalar::Boolean(false)))
        }
        Unary::ToLong => convert(v, ValueSort::Long),
        Unary::ToDouble => convert(v, ValueSort::Double),
        Unary::ToBoolean => convert(v, ValueSort::Boolean),
    }
}

/// Widen or narrow a numeric value to `to`, preserving nulls.
fn convert(v: &ColumnarValue, to: ValueSort) -> Result<ColumnarValue, ValueError> {
    let input = v.numeric("cast")?;
    let Some(n) = v.len() else {
        let (x, null) = input.get(0);
        return Ok(if null {
            ColumnarValue::Null(to)
        } else {
            ColumnarValue::Scalar(match to {
                ValueSort::Double => Scalar::Double(x),
                ValueSort::Boolean => Scalar::Boolean(x != 0.0),
                _ => Scalar::Long(x as c_long),
            })
        });
    };

    let mut nulls = vec![false; n];
    let mut any = false;
    let data = match to {
        ValueSort::Double => {
            let mut d = vec![0.0; n];
            for (i, slot) in d.iter_mut().enumerate() {
                let (x, null) = input.get(i);
                *slot = x;
                nulls[i] = null;
                any |= null;
            }
            ArrayData::Double(d)
        }
        ValueSort::Boolean => {
            let mut d = vec![false; n];
            for (i, slot) in d.iter_mut().enumerate() {
                let (x, null) = input.get(i);
                *slot = x != 0.0;
                nulls[i] = null;
                any |= null;
            }
            ArrayData::Boolean(d)
        }
        _ => {
            let mut d = vec![0 as c_long; n];
            for (i, slot) in d.iter_mut().enumerate() {
                let (x, null) = input.get(i);
                *slot = x as c_long;
                nulls[i] = null;
                any |= null;
            }
            ArrayData::Long(d)
        }
    };
    Ok(ColumnarValue::Array(if any {
        Array::with_nulls(data, nulls)
    } else {
        Array::new(data)
    }))
}

/// Read a boolean operand element by element, broadcasting a scalar.
fn bool_at(v: &ColumnarValue, i: usize) -> Result<(bool, bool), ValueError> {
    Ok(match v {
        ColumnarValue::Scalar(Scalar::Boolean(b)) => (*b, false),
        ColumnarValue::Null(ValueSort::Boolean) => (false, true),
        ColumnarValue::Array(a) => match a.data() {
            ArrayData::Boolean(d) => (d[i], a.is_null(i)),
            other => return Err(ValueError::BadSort("boolean operator", other.sort())),
        },
        other => return Err(ValueError::BadSort("boolean operator", other.sort())),
    })
}

fn logic(op: Logic, lhs: &ColumnarValue, rhs: &ColumnarValue) -> Result<ColumnarValue, ValueError> {
    if lhs.sort() != ValueSort::Boolean || rhs.sort() != ValueSort::Boolean {
        return Err(ValueError::Incompatible(
            "boolean operator",
            lhs.sort(),
            rhs.sort(),
        ));
    }
    let apply = |a: bool, b: bool| match op {
        Logic::And => a && b,
        Logic::Or => a || b,
        Logic::Eq => a == b,
        Logic::Ne => a != b,
    };

    let Some(n) = lhs.broadcast_len(rhs)? else {
        let (a, an) = bool_at(lhs, 0)?;
        let (b, bn) = bool_at(rhs, 0)?;
        return Ok(if an || bn {
            ColumnarValue::Null(ValueSort::Boolean)
        } else {
            ColumnarValue::Scalar(Scalar::Boolean(apply(a, b)))
        });
    };

    let mut data = vec![false; n];
    let mut nulls = vec![false; n];
    let mut any = false;
    for i in 0..n {
        let (a, an) = bool_at(lhs, if lhs.is_scalar() { 0 } else { i })?;
        let (b, bn) = bool_at(rhs, if rhs.is_scalar() { 0 } else { i })?;
        if an || bn {
            nulls[i] = true;
            any = true;
        } else {
            data[i] = apply(a, b);
        }
    }
    Ok(ColumnarValue::Array(if any {
        Array::with_nulls(ArrayData::Boolean(data), nulls)
    } else {
        Array::new(ArrayData::Boolean(data))
    }))
}

/// `cond ? then : els`.
///
/// Both branches are evaluated -- the engine works a column at a time, so
/// there is nothing to short-circuit -- and a null condition gives a null
/// result, matching `Do_Func`'s `ifthenelse_fct`.
fn if_then_else(
    cond: &ColumnarValue,
    then: &ColumnarValue,
    els: &ColumnarValue,
) -> Result<ColumnarValue, ValueError> {
    if cond.sort() != ValueSort::Boolean {
        return Err(ValueError::BadSort("?:", cond.sort()));
    }
    /* unlike arithmetic, `?:` over two booleans stays boolean */
    let out = then.sort().max(els.sort());

    let n = match (cond.len(), then.broadcast_len(els)?) {
        (None, None) => None,
        (Some(a), None) | (None, Some(a)) => Some(a),
        (Some(a), Some(b)) if a == b => Some(a),
        (Some(a), Some(b)) => return Err(ValueError::LengthMismatch(a, b)),
    };

    let pick = |i: usize| -> Result<(f64, bool), ValueError> {
        let (c, cn) = bool_at(cond, if cond.is_scalar() { 0 } else { i })?;
        if cn {
            return Ok((0.0, true));
        }
        let side = if c { then } else { els };
        let input = side.numeric("?:")?;
        Ok(input.get(if side.is_scalar() { 0 } else { i }))
    };

    let Some(n) = n else {
        let (v, null) = pick(0)?;
        return Ok(if null {
            ColumnarValue::Null(out)
        } else if out == ValueSort::Boolean {
            ColumnarValue::Scalar(Scalar::Boolean(v != 0.0))
        } else if out == ValueSort::Double {
            ColumnarValue::Scalar(Scalar::Double(v))
        } else {
            ColumnarValue::Scalar(Scalar::Long(v as c_long))
        });
    };

    let mut nulls = vec![false; n];
    let mut any = false;
    let mut vals = vec![0.0f64; n];
    for i in 0..n {
        let (v, null) = pick(i)?;
        vals[i] = v;
        nulls[i] = null;
        any |= null;
    }
    let data = match out {
        ValueSort::Double => ArrayData::Double(vals),
        ValueSort::Boolean => ArrayData::Boolean(vals.iter().map(|&v| v != 0.0).collect()),
        _ => ArrayData::Long(vals.iter().map(|&v| v as c_long).collect()),
    };
    Ok(ColumnarValue::Array(if any {
        Array::with_nulls(data, nulls)
    } else {
        Array::new(data)
    }))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn batch(cols: Vec<ColumnBatch>, n_rows: c_long) -> Batch {
        Batch {
            columns: cols,
            n_rows,
            first_row: 1,
        }
    }

    fn long_col(v: &[c_long]) -> ColumnBatch {
        ColumnBatch {
            data: ArrayData::Long(v.to_vec()),
            nulls: Vec::new(),
            nelem: 1,
        }
    }

    fn lit(v: c_long) -> Box<Expr> {
        Box::new(Expr::Literal(Scalar::Long(v)))
    }

    fn longs(v: &ColumnarValue) -> Vec<c_long> {
        match v {
            ColumnarValue::Array(a) => match a.data() {
                ArrayData::Long(d) => d.clone(),
                other => panic!("{other:?}"),
            },
            other => panic!("{other:?}"),
        }
    }
    fn bools(v: &ColumnarValue) -> Vec<bool> {
        match v {
            ColumnarValue::Array(a) => match a.data() {
                ArrayData::Boolean(d) => d.clone(),
                other => panic!("{other:?}"),
            },
            other => panic!("{other:?}"),
        }
    }

    #[test]
    fn a_column_evaluates_to_its_batch() {
        let b = batch(vec![long_col(&[7, -3, 10])], 3);
        let e = Expr::Column(0);
        assert_eq!(longs(&e.evaluate(&b).unwrap()), vec![7, -3, 10]);
    }

    #[test]
    fn nested_arithmetic_folds_left_to_right() {
        /* (col * 2) + 1 */
        let b = batch(vec![long_col(&[1, 2, 3])], 3);
        let e = Expr::Arith {
            op: Arith::Add,
            lhs: Box::new(Expr::Arith {
                op: Arith::Mul,
                lhs: Box::new(Expr::Column(0)),
                rhs: lit(2),
            }),
            rhs: lit(1),
        };
        assert_eq!(longs(&e.evaluate(&b).unwrap()), vec![3, 5, 7]);
    }

    #[test]
    fn a_wholly_constant_expression_stays_scalar() {
        /* 2 + 3 never materialises an array, however many rows there are */
        let b = batch(vec![], 1000);
        let e = Expr::Arith {
            op: Arith::Add,
            lhs: lit(2),
            rhs: lit(3),
        };
        let v = e.evaluate(&b).unwrap();
        assert!(v.is_scalar(), "got {v:?}");
        assert_eq!(v, ColumnarValue::Scalar(Scalar::Long(5)));
    }

    #[test]
    fn row_number_counts_from_the_batch_start() {
        let b = Batch {
            columns: vec![],
            n_rows: 3,
            first_row: 11,
        };
        assert_eq!(
            longs(&Expr::RowNumber.evaluate(&b).unwrap()),
            vec![11, 12, 13]
        );
    }

    #[test]
    fn comparison_and_logic_compose() {
        /* col > 2 && col < 10 */
        let b = batch(vec![long_col(&[1, 5, 20])], 3);
        let e = Expr::Logic {
            op: Logic::And,
            lhs: Box::new(Expr::Compare {
                op: Compare::Gt,
                lhs: Box::new(Expr::Column(0)),
                rhs: lit(2),
            }),
            rhs: Box::new(Expr::Compare {
                op: Compare::Lt,
                lhs: Box::new(Expr::Column(0)),
                rhs: lit(10),
            }),
        };
        assert_eq!(bools(&e.evaluate(&b).unwrap()), vec![false, true, false]);
    }

    #[test]
    fn negation_and_casts() {
        let b = batch(vec![long_col(&[1, -2, 3])], 3);
        let neg = Expr::Unary {
            op: Unary::Neg,
            arg: Box::new(Expr::Column(0)),
        };
        assert_eq!(longs(&neg.evaluate(&b).unwrap()), vec![-1, 2, -3]);

        /* (int) truncates toward zero, as C does */
        let trunc = Expr::Unary {
            op: Unary::ToLong,
            arg: Box::new(Expr::Literal(Scalar::Double(-2.7))),
        };
        assert_eq!(
            trunc.evaluate(&b).unwrap(),
            ColumnarValue::Scalar(Scalar::Long(-2))
        );
    }

    #[test]
    fn not_requires_a_boolean() {
        let b = batch(vec![long_col(&[1])], 1);
        let e = Expr::Unary {
            op: Unary::Not,
            arg: Box::new(Expr::Column(0)),
        };
        assert!(matches!(
            e.evaluate(&b),
            Err(ValueError::BadSort("!", ValueSort::Long))
        ));
    }

    #[test]
    fn conditional_picks_per_row_and_propagates_a_null_condition() {
        let b = batch(vec![long_col(&[1, 5, 9])], 3);
        let e = Expr::IfThenElse {
            cond: Box::new(Expr::Compare {
                op: Compare::Gt,
                lhs: Box::new(Expr::Column(0)),
                rhs: lit(4),
            }),
            then: lit(100),
            els: lit(0),
        };
        assert_eq!(longs(&e.evaluate(&b).unwrap()), vec![0, 100, 100]);

        let n = Expr::IfThenElse {
            cond: Box::new(Expr::Null(ValueSort::Boolean)),
            then: lit(1),
            els: lit(2),
        };
        assert_eq!(
            n.evaluate(&b).unwrap(),
            ColumnarValue::Null(ValueSort::Long)
        );
    }

    #[test]
    fn nulls_in_a_column_survive_a_whole_tree() {
        let mut col = long_col(&[1, 2, 3]);
        col.nulls = vec![false, true, false];
        let b = batch(vec![col], 3);
        let e = Expr::Arith {
            op: Arith::Add,
            lhs: Box::new(Expr::Arith {
                op: Arith::Mul,
                lhs: Box::new(Expr::Column(0)),
                rhs: lit(10),
            }),
            rhs: lit(1),
        };
        match e.evaluate(&b).unwrap() {
            ColumnarValue::Array(a) => {
                assert!(!a.is_null(0));
                assert!(a.is_null(1), "the null must survive two operations");
                assert!(!a.is_null(2));
            }
            other => panic!("{other:?}"),
        }
    }

    #[test]
    fn deep_nesting_does_not_overflow() {
        /* left-nested, so evaluation recurses to the full depth */
        let mut e = Expr::Literal(Scalar::Long(0));
        for _ in 0..500 {
            e = Expr::Arith {
                op: Arith::Add,
                lhs: Box::new(e),
                rhs: lit(1),
            };
        }
        let b = batch(vec![], 1);
        assert_eq!(
            e.evaluate(&b).unwrap(),
            ColumnarValue::Scalar(Scalar::Long(500))
        );
    }
}
