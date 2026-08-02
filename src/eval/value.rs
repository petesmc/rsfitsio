//! Values flowing through the expression evaluator.
//!
//! This is the first piece of the replacement for the `Node` arena described in
//! `PARSER_MIGRATION.md`. Two ideas from the old engine are collapsed here:
//!
//! * **Scalar versus vector.** Every `Do_*` routine hand-rolled "if my operand
//!   is a scalar constant do it once, otherwise loop" — 57 sites across the
//!   file. That distinction is [`ColumnarValue`], made once.
//!
//! * **The null mask.** `Node::value.undef` was a `*mut c_char` pointing into
//!   the *tail of the data allocation*, so the two could not be separated. Here
//!   an [`Array`] owns its own `nulls`.
//!
//! Nothing here indexes with raw pointers: an array is a `Vec`, and every
//! kernel works over slices.

use crate::c_types::{c_char, c_long};
use crate::eval_defs::ValueSort;

/// A single value, with no null of its own — a scalar that is null is
/// represented as [`ColumnarValue::Scalar`] with `null` set.
#[derive(Clone, Debug, PartialEq)]
pub(crate) enum Scalar {
    Boolean(bool),
    Long(c_long),
    Double(f64),
    /// A character string, without a terminator.
    Str(Vec<u8>),
    /// A bit string as ASCII `'0'` / `'1'`, as the engine stores them.
    Bits(Vec<u8>),
}

impl Scalar {
    pub(crate) fn sort(&self) -> ValueSort {
        match self {
            Scalar::Boolean(_) => ValueSort::Boolean,
            Scalar::Long(_) => ValueSort::Long,
            Scalar::Double(_) => ValueSort::Double,
            Scalar::Str(_) => ValueSort::String,
            Scalar::Bits(_) => ValueSort::Bits,
        }
    }

    /// The numeric value, widened. Booleans count as 0 or 1, matching the
    /// engine's promotion.
    pub(crate) fn as_f64(&self) -> Option<f64> {
        Some(match self {
            Scalar::Boolean(b) => f64::from(u8::from(*b)),
            Scalar::Long(v) => *v as f64,
            Scalar::Double(v) => *v,
            _ => return None,
        })
    }

    pub(crate) fn as_i64(&self) -> Option<c_long> {
        Some(match self {
            Scalar::Boolean(b) => c_long::from(*b),
            Scalar::Long(v) => *v,
            Scalar::Double(v) => *v as c_long,
            _ => return None,
        })
    }
}

/// The elements of an [`Array`], one `Vec` per sort.
#[derive(Clone, Debug, PartialEq)]
pub(crate) enum ArrayData {
    Boolean(Vec<bool>),
    Long(Vec<c_long>),
    Double(Vec<f64>),
    /// Fixed-width strings, one entry per element.
    Str(Vec<Vec<u8>>),
    /// Fixed-width bit strings, one entry per element.
    Bits(Vec<Vec<u8>>),
}

impl ArrayData {
    pub(crate) fn len(&self) -> usize {
        match self {
            ArrayData::Boolean(v) => v.len(),
            ArrayData::Long(v) => v.len(),
            ArrayData::Double(v) => v.len(),
            ArrayData::Str(v) => v.len(),
            ArrayData::Bits(v) => v.len(),
        }
    }

    pub(crate) fn sort(&self) -> ValueSort {
        match self {
            ArrayData::Boolean(_) => ValueSort::Boolean,
            ArrayData::Long(_) => ValueSort::Long,
            ArrayData::Double(_) => ValueSort::Double,
            ArrayData::Str(_) => ValueSort::String,
            ArrayData::Bits(_) => ValueSort::Bits,
        }
    }
}

/// A column of values with a null flag per element.
///
/// `nulls` is either empty — meaning nothing is null, the common case — or
/// exactly as long as the data. Keeping it empty avoids touching a second
/// allocation for the many expressions that cannot produce a null.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct Array {
    data: ArrayData,
    nulls: Vec<bool>,
}

impl Array {
    pub(crate) fn new(data: ArrayData) -> Self {
        Array {
            data,
            nulls: Vec::new(),
        }
    }

    pub(crate) fn with_nulls(data: ArrayData, nulls: Vec<bool>) -> Self {
        debug_assert!(
            nulls.is_empty() || nulls.len() == data.len(),
            "null mask must be empty or match the data length"
        );
        Array { data, nulls }
    }

    pub(crate) fn data(&self) -> &ArrayData {
        &self.data
    }

    pub(crate) fn len(&self) -> usize {
        self.data.len()
    }

    pub(crate) fn sort(&self) -> ValueSort {
        self.data.sort()
    }

    /// Whether element `i` is undefined.
    pub(crate) fn is_null(&self, i: usize) -> bool {
        self.nulls.get(i).copied().unwrap_or(false)
    }

    pub(crate) fn has_nulls(&self) -> bool {
        self.nulls.iter().any(|&n| n)
    }

    pub(crate) fn nulls(&self) -> &[bool] {
        &self.nulls
    }
}

/// What a sub-expression evaluates to for one batch of rows.
///
/// A `Scalar` is a value that is the same for every row — a literal, or a
/// sub-expression the parser folded. Keeping it distinct from a full `Array`
/// is what lets a kernel avoid materialising `n` copies of a constant, which
/// is the optimisation the old engine spelled out by hand in every operation.
#[derive(Clone, Debug, PartialEq)]
pub(crate) enum ColumnarValue {
    Scalar(Scalar),
    /// A scalar that is undefined for every row.
    Null(ValueSort),
    Array(Array),
}

impl ColumnarValue {
    pub(crate) fn sort(&self) -> ValueSort {
        match self {
            ColumnarValue::Scalar(s) => s.sort(),
            ColumnarValue::Null(t) => *t,
            ColumnarValue::Array(a) => a.sort(),
        }
    }

    /// Number of elements, or `None` for a scalar, which has no length of its
    /// own and adopts whatever the other operand has.
    pub(crate) fn len(&self) -> Option<usize> {
        match self {
            ColumnarValue::Array(a) => Some(a.len()),
            _ => None,
        }
    }

    pub(crate) fn is_scalar(&self) -> bool {
        !matches!(self, ColumnarValue::Array(_))
    }

    /// The length a binary operation over `self` and `other` produces, or
    /// `None` when both are scalars.
    ///
    /// Two arrays must agree; the engine's `Test_Dims` allowed a scalar
    /// against anything, which is the same rule.
    pub(crate) fn broadcast_len(&self, other: &ColumnarValue) -> Result<Option<usize>, ValueError> {
        match (self.len(), other.len()) {
            (None, None) => Ok(None),
            (Some(n), None) | (None, Some(n)) => Ok(Some(n)),
            (Some(a), Some(b)) if a == b => Ok(Some(a)),
            (Some(a), Some(b)) => Err(ValueError::LengthMismatch(a, b)),
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) enum ValueError {
    /// Two arrays of different lengths met in one operation.
    LengthMismatch(usize, usize),
    /// An operation was applied to a sort it is not defined for.
    BadSort(&'static str, ValueSort),
    /// Two operands whose sorts have no common operation.
    Incompatible(&'static str, ValueSort, ValueSort),
}

/// Helper for reading a numeric operand element by element, whether it is a
/// scalar broadcast over the batch or a real array.
///
/// This is the piece that removes the `vector1 > 1 ? buf[elem] : buf[row]`
/// branching the old kernels repeated in every loop.
#[derive(Debug)]
pub(crate) enum NumericInput<'a> {
    Scalar(f64, bool),
    Long(&'a [c_long], &'a [bool]),
    Double(&'a [f64], &'a [bool]),
    Boolean(&'a [bool], &'a [bool]),
}

impl NumericInput<'_> {
    pub(crate) fn get(&self, i: usize) -> (f64, bool) {
        match self {
            NumericInput::Scalar(v, null) => (*v, *null),
            NumericInput::Long(v, n) => (v[i] as f64, n.get(i).copied().unwrap_or(false)),
            NumericInput::Double(v, n) => (v[i], n.get(i).copied().unwrap_or(false)),
            NumericInput::Boolean(v, n) => (
                f64::from(u8::from(v[i])),
                n.get(i).copied().unwrap_or(false),
            ),
        }
    }
}

impl ColumnarValue {
    /// View this value as a numeric operand, or fail if it is not numeric.
    pub(crate) fn numeric(&self, op: &'static str) -> Result<NumericInput<'_>, ValueError> {
        match self {
            ColumnarValue::Scalar(s) => match s.as_f64() {
                Some(v) => Ok(NumericInput::Scalar(v, false)),
                None => Err(ValueError::BadSort(op, s.sort())),
            },
            ColumnarValue::Null(t) if t.is_expr() || *t == ValueSort::Boolean => {
                Ok(NumericInput::Scalar(0.0, true))
            }
            ColumnarValue::Null(t) => Err(ValueError::BadSort(op, *t)),
            ColumnarValue::Array(a) => Ok(match a.data() {
                ArrayData::Long(v) => NumericInput::Long(v, a.nulls()),
                ArrayData::Double(v) => NumericInput::Double(v, a.nulls()),
                ArrayData::Boolean(v) => NumericInput::Boolean(v, a.nulls()),
                other => return Err(ValueError::BadSort(op, other.sort())),
            }),
        }
    }
}

/// Convenience for a `c_char`-based null mask, which is how the old engine and
/// the FITS iterator both spell it.
pub(crate) fn nulls_from_chars(undef: &[c_char]) -> Vec<bool> {
    undef.iter().map(|&u| u != 0).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn long_array(v: &[c_long]) -> ColumnarValue {
        ColumnarValue::Array(Array::new(ArrayData::Long(v.to_vec())))
    }

    #[test]
    fn scalars_have_no_length_and_broadcast() {
        let s = ColumnarValue::Scalar(Scalar::Long(3));
        let a = long_array(&[1, 2, 3]);
        assert_eq!(s.len(), None);
        assert_eq!(a.len(), Some(3));
        assert_eq!(s.broadcast_len(&a), Ok(Some(3)));
        assert_eq!(a.broadcast_len(&s), Ok(Some(3)));
        assert_eq!(s.broadcast_len(&s), Ok(None));
    }

    #[test]
    fn arrays_of_different_lengths_do_not_broadcast() {
        let a = long_array(&[1, 2, 3]);
        let b = long_array(&[1, 2]);
        assert_eq!(a.broadcast_len(&b), Err(ValueError::LengthMismatch(3, 2)));
    }

    #[test]
    fn an_empty_null_mask_means_nothing_is_null() {
        let a = Array::new(ArrayData::Long(vec![1, 2, 3]));
        assert!(!a.has_nulls());
        for i in 0..3 {
            assert!(!a.is_null(i));
        }
    }

    #[test]
    fn null_masks_track_elements() {
        let a = Array::with_nulls(ArrayData::Long(vec![1, 2, 3]), vec![false, true, false]);
        assert!(a.has_nulls());
        assert!(!a.is_null(0));
        assert!(a.is_null(1));
        assert!(!a.is_null(2));
    }

    #[test]
    fn numeric_input_reads_scalars_and_arrays_alike() {
        let s = ColumnarValue::Scalar(Scalar::Long(7));
        let n = s.numeric("+").unwrap();
        /* a scalar reads the same at every index */
        assert_eq!(n.get(0), (7.0, false));
        assert_eq!(n.get(99), (7.0, false));

        let a = ColumnarValue::Array(Array::with_nulls(
            ArrayData::Double(vec![1.5, 2.5]),
            vec![false, true],
        ));
        let n = a.numeric("+").unwrap();
        assert_eq!(n.get(0), (1.5, false));
        assert_eq!(n.get(1), (2.5, true));
    }

    #[test]
    fn booleans_are_numeric_but_strings_are_not() {
        let b = ColumnarValue::Scalar(Scalar::Boolean(true));
        assert_eq!(b.numeric("+").unwrap().get(0), (1.0, false));

        let s = ColumnarValue::Scalar(Scalar::Str(b"abc".to_vec()));
        assert!(matches!(
            s.numeric("+"),
            Err(ValueError::BadSort("+", ValueSort::String))
        ));
    }

    #[test]
    fn a_null_scalar_reports_itself_as_null_at_every_index() {
        let n = ColumnarValue::Null(ValueSort::Long);
        let input = n.numeric("+").unwrap();
        assert_eq!(input.get(0), (0.0, true));
        assert_eq!(input.get(5), (0.0, true));
    }

    #[test]
    fn null_masks_convert_from_the_engines_char_form() {
        assert_eq!(nulls_from_chars(&[0, 1, 0]), vec![false, true, false]);
    }
}
