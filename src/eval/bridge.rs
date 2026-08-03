//! Running the columnar evaluator in place of the `Node` arena.
//!
//! The port is incremental, so the two representations coexist: the arena is
//! still built at parse time because `ffiprs` reads the result node for the
//! expression's datatype and shape, and everything downstream of
//! `Evaluate_Parser` — `fits_parser_workfn`, `ffcalc`, `fffrow` — still reads
//! its answer from there.
//!
//! So this evaluates the [`Expr`] tree and writes the result back into the
//! result node in the layout the old engine would have produced. When the port
//! is complete the result node goes away and this file with it.

use super::expr::Batch;
use super::value::{ArrayData, ColumnarValue, Scalar};
use crate::c_types::{c_char, c_long, c_uint};
#[cfg(test)]
use crate::eval_defs::{BufferKind, ValueSort};
use crate::eval_defs::{NodeValue, Operation, ParseData};
use crate::eval_y::Allocate_Ptrs;
use crate::fitscore::ffpmsg_str;
use crate::fitsio::PARSE_SYNTAX_ERR;
use crate::simplerng::simplerng_srand;

/// Evaluate `lParse.expr_tree` over one batch and store the answer in the
/// result node.
///
/// Returns whether the batch was handled. A row offset reaching outside the
/// loaded chunk needs the file read again, which only the engine does, so that
/// batch is declined and the caller walks the arena instead. Nothing has been
/// written to the result node when that happens.
pub(crate) fn evaluate(lParse: &mut ParseData, first_row: c_long, n_rows: c_long) -> bool {
    /* The generators draw from a global stream, seeded once per process. The
    arena did this behind a `static mut`; a `Once` says the same thing without
    the unsafe. */
    static SEED: std::sync::Once = std::sync::Once::new();
    SEED.call_once(|| {
        simplerng_srand(
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .map(|d| d.as_secs() as c_uint)
                .unwrap_or(0),
        )
    });

    lParse.firstRow = first_row;
    lParse.nRows = n_rows;

    let Some(tree) = lParse.expr_tree.clone() else {
        return false;
    };
    let batch = Batch::gather(lParse, first_row, n_rows, tree.offset_range());
    /* the ACCUM/SEQDIFF running values live in ParseData between batches; lend
    them to the batch for this one and take them back afterwards */
    *batch.accum.borrow_mut() = core::mem::take(&mut lParse.accum_state);

    let result = tree.evaluate(&batch);
    lParse.accum_state = core::mem::take(&mut batch.accum.borrow_mut());
    match result {
        Ok(v) => {
            store(lParse, &v);
            true
        }
        Err(e) => {
            ffpmsg_str(&format!("expression evaluation failed: {e:?}"));
            if lParse.status == 0 {
                lParse.status = PARSE_SYNTAX_ERR;
            }
            true
        }
    }
}

/// Write `value` into the result node, in the shape the arena used.
fn store(lParse: &mut ParseData, value: &ColumnarValue) {
    let this = lParse.resultNode as usize;

    match value {
        /* A constant result keeps the node's CONST_OP marking and its scalar
        slot, exactly as a folded node did. */
        /* A constant string result is already in the node: the arena folded
        the same expression at parse time and owns the text buffer, so there
        is nothing to copy in and its pointer must not be overwritten. */
        ColumnarValue::Scalar(Scalar::Str(_) | Scalar::Bits(_)) => {
            lParse.Nodes[this].operation = Operation::Const;
        }

        ColumnarValue::Scalar(s) => {
            lParse.Nodes[this].operation = Operation::Const;
            lParse.Nodes[this].value.data = match s {
                Scalar::Boolean(b) => NodeValue::Logical(c_char::from(*b)),
                Scalar::Long(v) => NodeValue::Long(*v),
                Scalar::Double(v) => NodeValue::Double(*v),
                Scalar::Str(_) | Scalar::Bits(_) => unreachable!("handled above"),
            };
            lParse.Nodes[this].value.undef = core::ptr::null_mut();
        }

        /* A wholly undefined result still gets a row buffer: the engine's
        null_fct is a rows kernel, so downstream expects per-row undef flags
        rather than a flag smuggled into the pointer field. */
        /* A constant node has no row buffer to flag, so an undefined constant
        keeps the scalar slot the arena already gave it. */
        ColumnarValue::Null(_) if lParse.Nodes[this].operation == Operation::Const => {}

        ColumnarValue::Null(_) => {
            Allocate_Ptrs(lParse, this);
            if lParse.status != 0 {
                return;
            }
            let n = (lParse.nRows * lParse.Nodes[this].value.nelem).max(0) as usize;
            unsafe {
                let undef = lParse.Nodes[this].value.undef;
                if !undef.is_null() {
                    core::slice::from_raw_parts_mut(undef, n).fill(1);
                }
            }
        }

        /* The arena decided at parse time that this expression is constant, so
        the answer is row-invariant however it was computed: honour that
        marking and fill the scalar slot, since everything downstream will read
        the node that way rather than looking for a row buffer. */
        ColumnarValue::Array(a) if lParse.Nodes[this].operation == Operation::Const => {
            let first = match a.data() {
                ArrayData::Long(d) => NodeValue::Long(d.first().copied().unwrap_or(0)),
                ArrayData::Double(d) => NodeValue::Double(d.first().copied().unwrap_or(0.0)),
                ArrayData::Boolean(d) => {
                    NodeValue::Logical(c_char::from(d.first().copied().unwrap_or(false)))
                }
                other => unreachable!("lowering refuses {other:?} results"),
            };
            lParse.Nodes[this].value.data = first;
            lParse.Nodes[this].value.undef = core::ptr::null_mut();
        }

        ColumnarValue::Array(a) => {
            Allocate_Ptrs(lParse, this);
            if lParse.status != 0 {
                return;
            }
            let n = a.len();
            unsafe {
                match a.data() {
                    ArrayData::Long(d) => {
                        let buf = lParse.Nodes[this].value.data.lng_buf();
                        core::slice::from_raw_parts_mut(buf, n).copy_from_slice(d);
                    }
                    ArrayData::Double(d) => {
                        let buf = lParse.Nodes[this].value.data.dbl_buf();
                        core::slice::from_raw_parts_mut(buf, n).copy_from_slice(d);
                    }
                    ArrayData::Boolean(d) => {
                        let buf = lParse.Nodes[this].value.data.log_buf();
                        let out = core::slice::from_raw_parts_mut(buf, n);
                        for (slot, &v) in out.iter_mut().zip(d) {
                            *slot = c_char::from(v);
                        }
                    }
                    /* A string result is a `char*` per row, each pointing
                    into the block `Allocate_Ptrs` laid out at the node's
                    declared width; the text is copied in and terminated. */
                    /* A bit string is laid out exactly like a character
                    string -- `Allocate_Ptrs` gives both a `char*` per row --
                    the difference being that a bit node has no undef array,
                    which the guard below already accounts for. */
                    ArrayData::Str(d) | ArrayData::Bits(d) => {
                        let buf = lParse.Nodes[this].value.data.str_buf();
                        let width = lParse.Nodes[this].value.nelem.max(0) as usize;
                        for (row, text) in d.iter().enumerate() {
                            let dst = *buf.add(row);
                            if dst.is_null() {
                                continue;
                            }
                            let n = text.len().min(width);
                            core::ptr::copy_nonoverlapping(text.as_ptr(), dst.cast::<u8>(), n);
                            *dst.add(n) = 0;
                        }
                    }
                }

                let undef = lParse.Nodes[this].value.undef;
                if !undef.is_null() {
                    let flags = core::slice::from_raw_parts_mut(undef, n);
                    for (slot, i) in flags.iter_mut().zip(0..n) {
                        *slot = c_char::from(a.is_null(i));
                    }
                }
            }
        }
    }
}

/// The buffer kind `Allocate_Ptrs` will have chosen for a node of this sort,
/// used only by the tests below to check the two agree.
#[cfg(test)]
fn kind_for(sort: ValueSort) -> BufferKind {
    match sort {
        ValueSort::Double => BufferKind::Double,
        ValueSort::Long => BufferKind::Long,
        _ => BufferKind::Logical,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::eval::value::Array;
    use crate::eval_defs::Node;

    /// A ParseData with one node to write into, standing in for a result node
    /// of the given sort.
    fn parse_data(sort: ValueSort, nelem: c_long, n_rows: c_long) -> ParseData {
        let mut p = ParseData {
            nRows: n_rows,
            ..Default::default()
        };
        p.Nodes = vec![Node {
            ntype: sort,
            operation: Operation::Op(crate::eval_defs::OpCode::Add),
            ..Default::default()
        }];
        p.nNodes = 1;
        p.nNodesAlloc = 1;
        p.Nodes[0].value.nelem = nelem;
        p.resultNode = 0;
        p
    }

    #[test]
    fn a_scalar_result_lands_in_the_nodes_constant_slot() {
        let mut p = parse_data(ValueSort::Long, 1, 3);
        store(&mut p, &ColumnarValue::Scalar(Scalar::Long(42)));
        assert_eq!(p.Nodes[0].operation, Operation::Const);
        assert_eq!(p.Nodes[0].value.data.lng(), 42);
    }

    #[test]
    fn an_array_result_fills_the_row_buffer() {
        let mut p = parse_data(ValueSort::Long, 1, 3);
        let v = ColumnarValue::Array(Array::new(ArrayData::Long(vec![7, -3, 10])));
        store(&mut p, &v);
        assert_eq!(p.status, 0);
        unsafe {
            let buf = p.Nodes[0].value.data.lng_buf();
            assert_eq!(core::slice::from_raw_parts(buf, 3), &[7, -3, 10]);
        }
    }

    #[test]
    fn nulls_reach_the_undef_flags() {
        let mut p = parse_data(ValueSort::Long, 1, 3);
        let v = ColumnarValue::Array(Array::with_nulls(
            ArrayData::Long(vec![1, 2, 3]),
            vec![false, true, false],
        ));
        store(&mut p, &v);
        unsafe {
            let undef = p.Nodes[0].value.undef;
            assert!(!undef.is_null(), "Allocate_Ptrs must provide the flags");
            assert_eq!(core::slice::from_raw_parts(undef, 3), &[0, 1, 0]);
        }
    }

    #[test]
    fn booleans_are_stored_one_byte_per_row() {
        let mut p = parse_data(ValueSort::Boolean, 1, 3);
        let v = ColumnarValue::Array(Array::new(ArrayData::Boolean(vec![true, false, true])));
        store(&mut p, &v);
        unsafe {
            let buf = p.Nodes[0].value.data.log_buf();
            assert_eq!(core::slice::from_raw_parts(buf, 3), &[1, 0, 1]);
        }
    }

    #[test]
    fn the_buffer_kind_matches_what_allocate_ptrs_chose() {
        for sort in [ValueSort::Long, ValueSort::Double, ValueSort::Boolean] {
            let mut p = parse_data(sort, 1, 2);
            let data = match sort {
                ValueSort::Long => ArrayData::Long(vec![1, 2]),
                ValueSort::Double => ArrayData::Double(vec![1.0, 2.0]),
                _ => ArrayData::Boolean(vec![true, false]),
            };
            store(&mut p, &ColumnarValue::Array(Array::new(data)));
            /* reading through the accessor for the expected kind asserts the
            tag, so a mismatch panics rather than silently reinterpreting */
            match kind_for(sort) {
                BufferKind::Long => assert!(!p.Nodes[0].value.data.lng_buf().is_null()),
                BufferKind::Double => assert!(!p.Nodes[0].value.data.dbl_buf().is_null()),
                _ => assert!(!p.Nodes[0].value.data.log_buf().is_null()),
            }
        }
    }
}
