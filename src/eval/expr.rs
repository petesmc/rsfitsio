//! The expression tree that replaces the `Node` arena.
//!
//! A node in the arena owned its output buffer and reached its children by
//! index, which is what forced raw pointers: `&mut Nodes[this]` together with
//! `&Nodes[child]` is not something the borrow checker will allow. Here a
//! parent owns its children and evaluation *returns* a value, so the two are
//! never borrowed at once and nothing needs a pointer.
//!
//! Evaluation is per batch of rows, matching how `ffiter` drives the engine.

use core::cell::RefCell;

use super::bits;
use super::kernel::{self, Arith, Compare};
use super::regions;
use super::strings;
use super::value::{Array, ArrayData, ColumnarValue, NumericInput, Scalar, ValueError};
use crate::c_types::{c_char, c_long};
use crate::eval_defs::OpCode;
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

/// A function applied element by element.
///
/// Only the elementwise ones so far. The reductions -- `SUM`, `MEDIAN`,
/// `NELEM` and the one-argument `MIN`/`MAX` -- fold across the elements *of a
/// row* rather than mapping over them, so they need a different shape and
/// come later.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Func {
    /* one argument, double result */
    Sin,
    Cos,
    Tan,
    Asin,
    Acos,
    Atan,
    Sinh,
    Cosh,
    Tanh,
    Exp,
    Ln,
    Log10,
    Sqrt,
    Ceil,
    Floor,
    /// `floor(x + 0.5)`, which is what the C computes -- not Rust's `round`,
    /// which rounds half away from zero and so differs at `-2.5`.
    Round,
    /// One argument, keeping its sort.
    Abs,
    /* two arguments */
    Atan2,
    Min,
    Max,
    /// `ISNULL(x)`: true where `x` is undefined. Never null itself.
    IsNull,
    /// `DEFNULL(x, y)`: `y` wherever `x` is undefined.
    DefNull,
    /// `ELEMENTNUM(x)`: each element's 1-based position within its row.
    ElementNum,
    /// `AXISELEM(x, n)`: each element's 1-based index along axis `n`.
    AxisElem,
    /// `ARRAY(x, dims)`: `x` broadcast into an array of the given shape.
    Array,
    /// `ANGSEP(ra1, dec1, ra2, dec2)`, in degrees.
    AngSep,
    /// `NEAR(x, y, tol)`.
    Near,
    /// `CIRCLE(xcen, ycen, rad, x, y)`.
    Circle,
    /// `BOX(xcen, ycen, xwid, ywid, rot, x, y)`.
    Box,
    /// `ELLIPSE(xcen, ycen, xrad, yrad, rot, x, y)`.
    Ellipse,
    /// `STRSTR(a, b)`: the 1-based position of `b` in `a`, undefined where it
    /// does not occur.
    StrStr,
    /// `STRMID(s, pos, len)`: a 1-based substring.
    StrMid,
    /// `SETNULL(sentinel, x)`: `x`, undefined wherever it equals `sentinel`.
    SetNull,

    /* Reductions. These fold across the elements *of a row* rather than
    mapping over them, so an argument with `nelem` elements per row yields one
    element per row. Every one of them skips undefined elements. */
    /// `SUM`: null only when the whole row is undefined.
    Sum,
    /// `AVERAGE`: the mean of the defined elements.
    Average,
    /// `STDDEV`: the sample deviation, and never null -- fewer than two
    /// defined elements give 0.
    Stddev,
    /// `MEDIAN`: the lower of the two middle defined elements.
    Median,
    /// One-argument `MIN`.
    Min1,
    /// One-argument `MAX`.
    Max1,
    /// `NVALID`: how many elements of the row are defined. Never null.
    NValid,
}

impl Func {
    /// Whether this folds a row rather than mapping over it.
    fn is_reduction(self) -> bool {
        matches!(
            self,
            Func::Sum
                | Func::Average
                | Func::Stddev
                | Func::Median
                | Func::Min1
                | Func::Max1
                | Func::NValid
        )
    }
}

/// The carried value of one `ACCUM` or `SEQDIFF` node.
///
/// Both sorts are kept because the node's own sort decides which is live:
/// `ACCUM(INTCOL)` sums in `i64` and `ACCUM(FLOATCOL)` in `f64`, and routing
/// the integer case through `f64` would lose exactness past 2^53.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub(crate) struct AccumState {
    pub(crate) prev: f64,
    pub(crate) prev_i: c_long,
    /// Whether the previous element was undefined, which only `SEQDIFF` reads.
    pub(crate) undef: bool,
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
    /// `==` or `!=`, whose meaning depends on the operand sorts: a numeric
    /// comparison, or boolean equality. `eval.y` had separate productions for
    /// each; here the sorts are only known once the operands are evaluated.
    Equality {
        negated: bool,
        lhs: Box<Expr>,
        rhs: Box<Expr>,
    },
    Call {
        func: Func,
        args: Vec<Expr>,
    },
    /// `{ a, b, c }`: one row of the result holds each element's value for
    /// that row, so the result has as many elements per row as the literal
    /// has entries.
    Vector(Vec<Expr>),
    /// `base[i, j, ...]` with one subscript per axis, so the result is a
    /// single element per row.
    Deref {
        base: Box<Expr>,
        idx: Vec<Expr>,
        /// The operand's axis lengths, innermost first, resolved at lowering.
        naxes: Vec<usize>,
    },
    /// `COL{k}`: the column's value `k` rows away from the row being
    /// evaluated. `k` is a constant, which the parser enforces.
    Offset {
        col: usize,
        offset: c_long,
    },
    /// `ACCUM(x)` and `SEQDIFF(x)`: a running total, and the difference from
    /// the previous element. Both run over the row-major element sequence and
    /// carry across batches, so each owns a slot in [`Batch::accum`].
    Accum {
        id: usize,
        /// `SEQDIFF` rather than `ACCUM`.
        diff: bool,
        arg: Box<Expr>,
    },
    /// `M[i]` on a column with more axes than subscripts: one index against
    /// the outermost axis selects a whole slice rather than an element. The
    /// engine only does this for a constant index.
    Slice {
        base: Box<Expr>,
        /// One-based index along the outermost axis.
        index: c_long,
        /// Elements in the slice, the product of the remaining axes.
        run: usize,
        /// Length of the axis being indexed, for the range check.
        axis_len: usize,
        /// The slice's own axis lengths.
        naxes: Vec<usize>,
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
    /// Axis lengths within a row, innermost first.
    pub(crate) naxes: Vec<usize>,
}

/// The rows an [`Expr`] is evaluated over.
pub(crate) struct Batch {
    pub(crate) columns: Vec<ColumnBatch>,
    pub(crate) n_rows: c_long,
    /// The first row of the chunk the columns were loaded from, and how many
    /// it holds. A row offset can reach anywhere inside this window; past it
    /// the engine has to reload, which the arena still does.
    pub(crate) first_data_row: c_long,
    pub(crate) n_data_rows: c_long,
    /// Rows in the whole table, which bounds what a row offset can name at
    /// all: past either end the answer is undefined rather than a reload.
    pub(crate) total_rows: c_long,
    /// Table row number of the first row in the batch, one-based.
    pub(crate) first_row: c_long,
    /// The running values of the `ACCUM` and `SEQDIFF` nodes, which carry from
    /// one batch to the next. The engine keeps the same state in the constant
    /// subnode each of those operators is paired with; here it is handed in by
    /// the bridge and handed back after the batch, so `evaluate` stays a
    /// `&self` walk.
    pub(crate) accum: RefCell<Vec<AccumState>>,
}

impl Batch {
    /// Column `col` as seen from each row of the batch, shifted by `offset`
    /// rows -- so `offset == 0` is the plain column reference.
    ///
    /// A row naming something outside the table is undefined, which is how
    /// `Do_Offset` reports the rows before the first and after the last. A row
    /// inside the table but outside the loaded chunk needs the file read again,
    /// which only the engine does; that is [`ValueError::NeedsReload`] and the
    /// caller hands the whole batch back to the arena.
    fn shifted(&self, col: usize, offset: c_long) -> Result<ColumnarValue, ValueError> {
        let c = &self.columns[col];
        let nelem = c.nelem.max(1) as usize;
        /* a text column keeps one entry per row whatever its declared width */
        let stride = if matches!(c.data, ArrayData::Str(_) | ArrayData::Bits(_)) {
            1
        } else {
            nelem
        };
        let rows = self.n_rows.max(0) as usize;

        let mut picked: Vec<Option<usize>> = Vec::with_capacity(rows);
        for r in 0..rows {
            let target = self.first_row + r as c_long + offset;
            /* `totalRows` is 0 for an image with no NAXIS2, where it says
            nothing about the top end -- so it only bounds when it is set. A
            plain column reference has offset 0 and is always in range. */
            if target < 1 || (self.total_rows > 0 && target > self.total_rows) {
                picked.push(None);
                continue;
            }
            let within = target - self.first_data_row;
            if within < 0 || within >= self.n_data_rows {
                return Err(ValueError::NeedsReload);
            }
            picked.push(Some(within as usize));
        }

        /* the fast path: every row is in the chunk at a fixed shift, so the
        window can be sliced rather than gathered element by element */
        if picked.iter().all(Option::is_some) {
            let start = picked[0].unwrap() * stride;
            let len = rows * stride;
            let data = c.data.slice(start, len);
            let nulls = if c.nulls.len() >= start + len {
                c.nulls[start..start + len].to_vec()
            } else {
                Vec::new()
            };
            return Ok(ColumnarValue::Array(
                Array::with_nulls(data, nulls)
                    .with_nelem(nelem)
                    .with_naxes(c.naxes.clone()),
            ));
        }

        /* otherwise the rows off the ends of the table are undefined */
        let mut idx = Vec::with_capacity(rows * stride);
        for p in &picked {
            for e in 0..stride {
                idx.push(p.map(|row| row * stride + e));
            }
        }
        let src = Array::with_nulls(c.data.clone(), c.nulls.clone());
        Ok(ColumnarValue::Array(
            gather_at(&src, &idx)
                .with_nelem(nelem)
                .with_naxes(c.naxes.clone()),
        ))
    }

    /// Copy the current chunk out of `ParseData::varData`.
    ///
    /// This is the one place that touches the iterator's raw buffers. It copies
    /// rather than borrows: a borrowed view would have to carry a lifetime
    /// through every kernel, and the copy is a memcpy per column per chunk
    /// against per-element work downstream. Making it a borrow is a later
    /// optimisation, not a correctness question.
    pub(crate) fn gather(lParse: &ParseData, first_row: c_long, n_rows: c_long) -> Batch {
        /* The whole loaded chunk is gathered, not just the batch's slice of
        it, so a row offset can reach the neighbouring rows without going back
        to the file. `Expr::Column` slices the batch window back out. */
        let chunk_rows = lParse.nDataRows.max(n_rows);
        let mut columns = Vec::with_capacity(lParse.nCols as usize);

        for var in lParse.varData.iter().take(lParse.nCols as usize) {
            let count = (var.nelem * chunk_rows).max(0) as usize;
            let offset = 0usize;

            /* A column the iterator has not filled in for this pass -- the
            expression may not reference it at all -- has no buffer yet. */
            if var.data.is_null() || count == 0 {
                columns.push(ColumnBatch {
                    data: ArrayData::Long(Vec::new()),
                    nulls: Vec::new(),
                    nelem: var.nelem,
                    naxes: Vec::new(),
                });
                continue;
            }

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
                    /* a string or bit-string column is a `char*` per row,
                    not a value per element, so it is indexed by row */
                    ValueSort::String | ValueSort::Bits => {
                        let rows = var.data.cast::<*mut c_char>();
                        let n = chunk_rows.max(0) as usize;
                        let text: Vec<Vec<u8>> = (0..n)
                            .map(|r| {
                                let p = *rows.add(r);
                                if p.is_null() {
                                    Vec::new()
                                } else {
                                    let mut len = 0;
                                    while *p.add(len) != 0 {
                                        len += 1;
                                    }
                                    core::slice::from_raw_parts(p.cast::<u8>(), len).to_vec()
                                }
                            })
                            .collect();
                        if var.dtype == ValueSort::Bits {
                            ArrayData::Bits(text)
                        } else {
                            ArrayData::Str(text)
                        }
                    }
                }
            };

            /* the text sorts hold one entry per row, so their flags are
            indexed by row rather than by element */
            let text = matches!(var.dtype, ValueSort::String | ValueSort::Bits);
            let (off, len) = if text {
                (0, chunk_rows.max(0) as usize)
            } else {
                (offset, count)
            };
            let nulls = match &var.undef {
                Some(u) if u.len() >= off + len => {
                    u[off..off + len].iter().map(|&v| v != 0).collect()
                }
                _ => Vec::new(),
            };

            columns.push(ColumnBatch {
                data,
                nulls,
                /* A text column holds one entry per row; `nelem` carries its
                declared width instead, which STRMID measures against. */
                nelem: var.nelem,
                naxes: if text {
                    Vec::new()
                } else {
                    var.naxes[..(var.naxis.max(0) as usize).min(var.naxes.len())]
                        .iter()
                        .map(|&n| n.max(0) as usize)
                        .collect()
                },
            });
        }

        Batch {
            columns,
            n_rows,
            first_data_row: lParse.firstDataRow,
            n_data_rows: chunk_rows,
            total_rows: lParse.totalRows,
            first_row,
            accum: RefCell::new(Vec::new()),
        }
    }
}

impl Expr {
    /// The sort this expression produces, known without evaluating it.
    pub(crate) fn sort(&self, batch_sorts: &dyn Fn(usize) -> ValueSort) -> ValueSort {
        match self {
            Expr::Literal(s) => s.sort(),
            Expr::Null(t) => *t,
            Expr::Column(i) | Expr::Offset { col: i, .. } => batch_sorts(*i),
            Expr::RowNumber => ValueSort::Long,
            Expr::Unary { op, arg } => match op {
                /* `!` over a bit string complements it rather than testing it */
                Unary::Not => match arg.sort(batch_sorts) {
                    ValueSort::Bits => ValueSort::Bits,
                    _ => ValueSort::Boolean,
                },
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
            Expr::Compare { .. } | Expr::Logic { .. } | Expr::Equality { .. } => ValueSort::Boolean,
            Expr::IfThenElse { then, els, .. } => then.sort(batch_sorts).max(els.sort(batch_sorts)),
            Expr::Vector(items) => items
                .iter()
                .map(|e| e.sort(batch_sorts))
                .max()
                .unwrap_or(ValueSort::Long),
            Expr::Deref { base, .. } => base.sort(batch_sorts),
            Expr::Slice { base, .. } => base.sort(batch_sorts),
            Expr::Accum { arg, .. } => match arg.sort(batch_sorts) {
                ValueSort::Boolean => ValueSort::Long,
                other => other,
            },
            Expr::Call { func, args } => match func {
                Func::IsNull => ValueSort::Boolean,
                Func::NValid | Func::StrStr => ValueSort::Long,
                Func::AngSep => ValueSort::Double,
                Func::ElementNum | Func::AxisElem => ValueSort::Long,
                Func::Array => match args[0].sort(batch_sorts) {
                    ValueSort::Boolean => ValueSort::Long,
                    other => other,
                },
                Func::Near | Func::Circle | Func::Box | Func::Ellipse => ValueSort::Boolean,
                Func::StrMid => ValueSort::String,
                Func::Average | Func::Stddev => ValueSort::Double,
                Func::Sum | Func::Median | Func::Min1 | Func::Max1 => {
                    match args[0].sort(batch_sorts) {
                        ValueSort::Boolean => ValueSort::Long,
                        other => other,
                    }
                }
                Func::Abs | Func::SetNull => args[0].sort(batch_sorts),
                Func::Min | Func::Max | Func::DefNull => args
                    .iter()
                    .map(|a| a.sort(batch_sorts))
                    .max()
                    .unwrap_or(ValueSort::Long),
                _ => ValueSort::Double,
            },
        }
    }

    pub(crate) fn evaluate(&self, batch: &Batch) -> Result<ColumnarValue, ValueError> {
        match self {
            Expr::Literal(s) => Ok(ColumnarValue::Scalar(s.clone())),
            /* `#NULL` is a rows kernel in the engine -- it fills a buffer and
            sets every undef flag -- so it must be an array here too, or the
            caller never sees anynul. */
            Expr::Null(t) => {
                let n = batch.n_rows.max(0) as usize;
                if n == 0 {
                    return Ok(ColumnarValue::Null(*t));
                }
                /* the sort has to survive: an undefined *string* still has
                to read as text downstream, which `#SNULL` in a conditional
                depends on */
                let data = match t {
                    ValueSort::Double => ArrayData::Double(vec![0.0; n]),
                    ValueSort::Boolean => ArrayData::Boolean(vec![false; n]),
                    ValueSort::String => ArrayData::Str(vec![Vec::new(); n]),
                    ValueSort::Bits => ArrayData::Bits(vec![Vec::new(); n]),
                    _ => ArrayData::Long(vec![0; n]),
                };
                Ok(ColumnarValue::Array(Array::with_nulls(data, vec![true; n])))
            }

            /* The gathered column covers the whole chunk, so a plain
            reference takes the batch's window out of it -- a shift of zero. */
            Expr::Column(i) => batch.shifted(*i, 0),

            Expr::Offset { col, offset } => batch.shifted(*col, *offset),

            Expr::RowNumber => {
                let n = batch.n_rows.max(0) as usize;
                let rows: Vec<c_long> = (0..n).map(|i| batch.first_row + i as c_long).collect();
                Ok(ColumnarValue::Array(Array::new(ArrayData::Long(rows))))
            }

            Expr::Unary { op, arg } => {
                let v = arg.evaluate(batch)?;
                if *op == Unary::Not && is_bits(&v) {
                    bit_not(&v, batch.n_rows.max(0) as usize)
                } else {
                    unary(*op, &v)
                }
            }

            Expr::Arith { op, lhs, rhs } => {
                let l = lhs.evaluate(batch)?;
                let r = rhs.evaluate(batch)?;
                /* the parser has already rejected every mixed combination, so
                two bit strings here mean the bit kernels rather than the
                numeric ones */
                if is_bits(&l) && is_bits(&r) {
                    bit_arith(*op, &l, &r, batch.n_rows.max(0) as usize)
                } else if is_text(&l) && is_text(&r) {
                    str_arith(*op, &l, &r, batch.n_rows.max(0) as usize)
                } else {
                    kernel::arith(*op, &l, &r)
                }
            }

            Expr::Compare { op, lhs, rhs } => {
                let l = lhs.evaluate(batch)?;
                let r = rhs.evaluate(batch)?;
                if is_bits(&l) && is_bits(&r) {
                    bit_compare(*op, &l, &r, batch.n_rows.max(0) as usize)
                } else if is_text(&l) && is_text(&r) {
                    str_compare(*op, &l, &r, batch.n_rows.max(0) as usize)
                } else {
                    kernel::compare(*op, &l, &r)
                }
            }

            Expr::Logic { op, lhs, rhs } => {
                let l = lhs.evaluate(batch)?;
                let r = rhs.evaluate(batch)?;
                logic(*op, &l, &r)
            }

            Expr::Equality { negated, lhs, rhs } => {
                let l = lhs.evaluate(batch)?;
                let r = rhs.evaluate(batch)?;
                if is_bits(&l) && is_bits(&r) {
                    bit_compare(
                        if *negated { Compare::Ne } else { Compare::Eq },
                        &l,
                        &r,
                        batch.n_rows.max(0) as usize,
                    )
                } else if is_text(&l) && is_text(&r) {
                    str_compare(
                        if *negated { Compare::Ne } else { Compare::Eq },
                        &l,
                        &r,
                        batch.n_rows.max(0) as usize,
                    )
                } else if l.sort() == ValueSort::Boolean && r.sort() == ValueSort::Boolean {
                    logic(if *negated { Logic::Ne } else { Logic::Eq }, &l, &r)
                } else {
                    kernel::compare(if *negated { Compare::Ne } else { Compare::Eq }, &l, &r)
                }
            }

            Expr::Call { func, args } => {
                let vals: Result<Vec<_>, _> = args.iter().map(|a| a.evaluate(batch)).collect();
                call(*func, &vals?, batch.n_rows.max(0) as usize)
            }

            Expr::Deref { base, idx, naxes } => {
                let base = base.evaluate(batch)?;
                let idx: Result<Vec<_>, _> = idx.iter().map(|e| e.evaluate(batch)).collect();
                deref(&base, &idx?, naxes, batch.n_rows.max(0) as usize)
            }
            Expr::Accum { id, diff, arg } => {
                let v = arg.evaluate(batch)?;
                accumulate(*id, *diff, &v, batch)
            }

            Expr::Slice {
                base,
                index,
                run,
                axis_len,
                naxes,
            } => {
                let v = base.evaluate(batch)?;
                slice_of(
                    &v,
                    *index,
                    *run,
                    *axis_len,
                    naxes,
                    batch.n_rows.max(0) as usize,
                )
            }

            Expr::Vector(items) => {
                let vals: Result<Vec<_>, _> = items.iter().map(|e| e.evaluate(batch)).collect();
                vector(&vals?, batch.n_rows.max(0) as usize)
            }

            Expr::IfThenElse { cond, then, els } => {
                let c = cond.evaluate(batch)?;
                let t = then.evaluate(batch)?;
                let e = els.evaluate(batch)?;
                if is_text(&t) || is_text(&e) {
                    text_if_then_else(&c, &t, &e, batch.n_rows.max(0) as usize)
                } else {
                    if_then_else(&c, &t, &e)
                }
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
        let (x, null) = input.get(0, input.nelem());
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
                let (x, null) = input.get(i, input.nelem());
                *slot = x;
                nulls[i] = null;
                any |= null;
            }
            ArrayData::Double(d)
        }
        ValueSort::Boolean => {
            let mut d = vec![false; n];
            for (i, slot) in d.iter_mut().enumerate() {
                let (x, null) = input.get(i, input.nelem());
                *slot = x != 0.0;
                nulls[i] = null;
                any |= null;
            }
            ArrayData::Boolean(d)
        }
        _ => {
            let mut d = vec![0 as c_long; n];
            for (i, slot) in d.iter_mut().enumerate() {
                let (x, null) = input.get(i, input.nelem());
                *slot = x as c_long;
                nulls[i] = null;
                any |= null;
            }
            ArrayData::Long(d)
        }
    };
    let nelem = match v {
        ColumnarValue::Array(a) => a.nelem(),
        _ => 1,
    };
    Ok(ColumnarValue::Array(
        if any {
            Array::with_nulls(data, nulls)
        } else {
            Array::new(data)
        }
        .with_nelem(nelem),
    ))
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

    let nelem_of = |v: &ColumnarValue| match v {
        ColumnarValue::Array(a) => a.nelem(),
        _ => 1,
    };
    let (ln, rn_) = (nelem_of(lhs), nelem_of(rhs));
    let out_nelem = ln.max(rn_);
    let idx = |v: &ColumnarValue, own: usize, i: usize| {
        if v.is_scalar() {
            0
        } else if own == out_nelem {
            i
        } else {
            i / out_nelem.max(1)
        }
    };

    let mut data = vec![false; n];
    let mut nulls = vec![false; n];
    let mut any = false;
    for i in 0..n {
        let (a, an) = bool_at(lhs, idx(lhs, ln, i))?;
        let (b, bn) = bool_at(rhs, idx(rhs, rn_, i))?;
        if an || bn {
            nulls[i] = true;
            any = true;
        } else {
            data[i] = apply(a, b);
        }
    }
    Ok(ColumnarValue::Array(
        if any {
            Array::with_nulls(ArrayData::Boolean(data), nulls)
        } else {
            Array::new(ArrayData::Boolean(data))
        }
        .with_nelem(out_nelem),
    ))
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
        /* a per-row condition against per-element branches */
        (Some(a), Some(b)) if a != 0 && b != 0 && (a % b == 0 || b % a == 0) => Some(a.max(b)),
        (Some(a), Some(b)) => return Err(ValueError::LengthMismatch(a, b)),
    };

    /* the branches decide the shape; the condition is per row and repeats
    across the elements of that row, as `Do_Func`'s ifthenelse_fct did */
    let out_nelem = [then, els]
        .iter()
        .filter_map(|v| match v {
            ColumnarValue::Array(a) => Some(a.nelem()),
            _ => None,
        })
        .max()
        .unwrap_or(1);
    let cond_nelem = match cond {
        ColumnarValue::Array(a) => a.nelem(),
        _ => 1,
    };

    let pick = |i: usize| -> Result<(f64, bool), ValueError> {
        let ci = if cond.is_scalar() {
            0
        } else if cond_nelem == out_nelem {
            i
        } else {
            i / out_nelem.max(1)
        };
        let (c, cn) = bool_at(cond, ci)?;
        if cn {
            return Ok((0.0, true));
        }
        let side = if c { then } else { els };
        let input = side.numeric("?:")?;
        Ok(input.get(if side.is_scalar() { 0 } else { i }, out_nelem))
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
    Ok(ColumnarValue::Array(
        if any {
            Array::with_nulls(data, nulls)
        } else {
            Array::new(data)
        }
        .with_nelem(out_nelem),
    ))
}

/// Apply an elementwise function.
fn call(func: Func, args: &[ColumnarValue], rows: usize) -> Result<ColumnarValue, ValueError> {
    if func.is_reduction() {
        return reduce(func, &args[0]);
    }
    match func {
        Func::IsNull => is_null(&args[0]),
        Func::DefNull => def_null(&args[0], &args[1]),
        /* the parser writes SETNULL(sentinel, value); `New_Func` reordered
        them, so the value is the second argument here */
        Func::SetNull => set_null(&args[1], &args[0]),
        Func::Abs => map1_keep_sort(&args[0], |v| v.abs(), |v| v.wrapping_abs()),
        Func::Atan2 => kernel::zip_double(&args[0], &args[1], "ARCTAN2", |a, b| a.atan2(b)),
        Func::Min => kernel::zip_keep_sort(&args[0], &args[1], "MIN", f64::min, c_long::min),
        Func::Max => kernel::zip_keep_sort(&args[0], &args[1], "MAX", f64::max, c_long::max),
        Func::StrStr => str_str(&args[0], &args[1], rows),
        Func::ElementNum => Ok(element_num(&args[0], rows)),
        Func::AxisElem => axis_elem(&args[0], &args[1], rows),
        Func::Array => array_of(&args[0], &args[1], rows),
        Func::AngSep => nary(args, rows, "ANGSEP", |a| {
            NaryOut::Double(regions::angsep(a[0], a[1], a[2], a[3]))
        }),
        Func::Near => nary(args, rows, "NEAR", |a| {
            NaryOut::Boolean(regions::near(a[0], a[1], a[2]))
        }),
        Func::Circle => nary(args, rows, "CIRCLE", |a| {
            NaryOut::Boolean(regions::circle(a[0], a[1], a[2], a[3], a[4]))
        }),
        Func::Box => nary(args, rows, "BOX", |a| {
            NaryOut::Boolean(regions::saobox(a[0], a[1], a[2], a[3], a[4], a[5], a[6]))
        }),
        Func::Ellipse => nary(args, rows, "ELLIPSE", |a| {
            NaryOut::Boolean(regions::ellipse(a[0], a[1], a[2], a[3], a[4], a[5], a[6]))
        }),
        Func::StrMid => str_mid(&args[0], &args[1], &args[2], rows),
        _ => {
            let f: fn(f64) -> f64 = match func {
                Func::Sin => f64::sin,
                Func::Cos => f64::cos,
                Func::Tan => f64::tan,
                Func::Asin => f64::asin,
                Func::Acos => f64::acos,
                Func::Atan => f64::atan,
                Func::Sinh => f64::sinh,
                Func::Cosh => f64::cosh,
                Func::Tanh => f64::tanh,
                Func::Exp => f64::exp,
                Func::Ln => f64::ln,
                Func::Log10 => f64::log10,
                Func::Sqrt => f64::sqrt,
                Func::Ceil => f64::ceil,
                Func::Floor => f64::floor,
                Func::Round => |v: f64| (v + 0.5).floor(),
                _ => unreachable!("handled above"),
            };
            /* the engine turns a domain error into a null rather than
            letting a NaN through: sqrt and log check their argument and set
            the undef flag. */
            let domain: fn(f64) -> bool = match func {
                Func::Sqrt => |v| v < 0.0,
                Func::Ln | Func::Log10 => |v| v <= 0.0,
                Func::Asin | Func::Acos => |v| !(-1.0..=1.0).contains(&v),
                _ => |_| false,
            };
            kernel::map_double_checked(&args[0], "function", f, domain)
        }
    }
}

/// `cond ? a : b` over string branches, which the numeric kernel cannot carry.
/// An undefined condition leaves the row undefined, as does taking a branch
/// that is itself undefined there.
fn text_if_then_else(
    cond: &ColumnarValue,
    then: &ColumnarValue,
    els: &ColumnarValue,
    rows: usize,
) -> Result<ColumnarValue, ValueError> {
    let c = cond.numeric("?:")?;
    let (t, e) = (then.text("?:")?, els.text("?:")?);
    let mut out = Vec::with_capacity(rows);
    let mut nulls = Vec::with_capacity(rows);
    for i in 0..rows {
        let (flag, cnull) = c.get_i64(i, 1);
        if cnull {
            out.push(Vec::new());
            nulls.push(true);
            continue;
        }
        let (text, null) = if flag != 0 { t.get(i) } else { e.get(i) };
        out.push(text.to_vec());
        nulls.push(null);
    }
    Ok(ColumnarValue::Array(Array::with_nulls(
        ArrayData::Str(out),
        nulls,
    )))
}

/// `M[i]` where `i` indexes the outermost axis and the rest of the shape comes
/// through: each row yields the `run` elements starting at `run * (i - 1)`.
fn slice_of(
    v: &ColumnarValue,
    index: c_long,
    run: usize,
    axis_len: usize,
    naxes: &[usize],
    rows: usize,
) -> Result<ColumnarValue, ValueError> {
    if index < 1 || index as usize > axis_len {
        return Err(ValueError::OutOfRange);
    }
    let ColumnarValue::Array(a) = v else {
        return Ok(ColumnarValue::Null(v.sort()));
    };
    let src_nelem = a.nelem().max(1);
    let start = run * (index as usize - 1);
    let picked: Vec<Option<usize>> = (0..rows)
        .flat_map(|r| (0..run).map(move |e| Some(r * src_nelem + start + e)))
        .collect();
    Ok(ColumnarValue::Array(
        gather_at(a, &picked)
            .with_nelem(run)
            .with_naxes(naxes.to_vec()),
    ))
}

/// `ACCUM(x)` and `SEQDIFF(x)` over one batch, carrying the running value in
/// and out of [`Batch::accum`].
///
/// `ACCUM` skips undefined elements -- they add nothing -- and its result is
/// never undefined. `SEQDIFF` is undefined wherever the current element or the
/// one before it is, and still advances its previous value across a gap.
fn accumulate(
    id: usize,
    diff: bool,
    v: &ColumnarValue,
    batch: &Batch,
) -> Result<ColumnarValue, ValueError> {
    let sort = match v.sort() {
        ValueSort::Boolean => ValueSort::Long,
        other => other,
    };
    if !matches!(sort, ValueSort::Long | ValueSort::Double) {
        return Err(ValueError::BadSort("ACCUM", v.sort()));
    }
    let input = v.numeric("ACCUM")?;
    let (nelem, _) = shape_at_runtime(v);
    let n = batch.n_rows.max(0) as usize * nelem;

    let mut states = batch.accum.borrow_mut();
    let st = states
        .get_mut(id)
        .ok_or(ValueError::BadSort("ACCUM", sort))?;

    let mut longs = Vec::new();
    let mut doubles = Vec::new();
    let mut nulls = Vec::with_capacity(n);
    for i in 0..n {
        let (x, undef) = input.get(i, nelem);
        let xi = input.get_i64(i, nelem).0;
        if diff {
            /* a gap makes this element undefined, but the previous value still
            advances past it */
            let bad = undef || st.undef;
            nulls.push(bad);
            if sort == ValueSort::Double {
                doubles.push(if bad { 0.0 } else { x - st.prev });
            } else {
                longs.push(if bad { 0 } else { xi.wrapping_sub(st.prev_i) });
            }
            st.prev = x;
            st.prev_i = xi;
            st.undef = undef;
        } else {
            if !undef {
                st.prev += x;
                st.prev_i = st.prev_i.wrapping_add(xi);
            }
            nulls.push(false);
            if sort == ValueSort::Double {
                doubles.push(st.prev);
            } else {
                longs.push(st.prev_i);
            }
        }
    }
    let data = if sort == ValueSort::Double {
        ArrayData::Double(doubles)
    } else {
        ArrayData::Long(longs)
    };
    Ok(ColumnarValue::Array(
        Array::with_nulls(data, nulls).with_nelem(nelem),
    ))
}

/// The shape of an evaluated argument: its elements per row and axis lengths.
fn shape_at_runtime(v: &ColumnarValue) -> (usize, Vec<usize>) {
    match v {
        ColumnarValue::Array(a) => (a.nelem().max(1), a.naxes()),
        _ => (1, vec![1]),
    }
}

/// `ELEMENTNUM(x)`: 1..nelem repeated for each row. Never undefined.
fn element_num(v: &ColumnarValue, rows: usize) -> ColumnarValue {
    let (nelem, _) = shape_at_runtime(v);
    let out: Vec<c_long> = (0..rows * nelem)
        .map(|i| (i % nelem) as c_long + 1)
        .collect();
    ColumnarValue::Array(Array::new(ArrayData::Long(out)).with_nelem(nelem))
}

/// `AXISELEM(x, n)`: each element's 1-based index along axis `n`.
///
/// The axes are innermost first, so axis 1 cycles fastest: over a 2x3 column
/// it runs 1,2,1,2,1,2 while axis 2 runs 1,1,2,2,3,3.
fn axis_elem(
    v: &ColumnarValue,
    axis: &ColumnarValue,
    rows: usize,
) -> Result<ColumnarValue, ValueError> {
    let (nelem, naxes) = shape_at_runtime(v);
    let n = axis.numeric("AXISELEM")?.get_i64(0, 1).0.max(1) as usize;
    /* the stride of an axis is the product of the ones inside it */
    let stride: usize = naxes.iter().take(n - 1).product::<usize>().max(1);
    let len = naxes.get(n - 1).copied().unwrap_or(1).max(1);
    let out: Vec<c_long> = (0..rows * nelem)
        .map(|i| ((i % nelem) / stride % len) as c_long + 1)
        .collect();
    Ok(ColumnarValue::Array(
        Array::new(ArrayData::Long(out)).with_nelem(nelem),
    ))
}

/// `ARRAY(x, dims)`: `x` repeated into an array of the given shape, one such
/// array per row. `dims` is a constant count or a vector of axis lengths.
fn array_of(
    v: &ColumnarValue,
    dims: &ColumnarValue,
    rows: usize,
) -> Result<ColumnarValue, ValueError> {
    let d = dims.numeric("ARRAY")?;
    let (dim_count, _) = shape_at_runtime(dims);
    let naxes: Vec<usize> = (0..dim_count)
        .map(|i| d.get_i64(i, dim_count).0.max(0) as usize)
        .collect();
    let nelem: usize = naxes.iter().product::<usize>().max(1);

    let src = v.numeric("ARRAY")?;
    let sort = match v.sort() {
        ValueSort::Boolean => ValueSort::Long,
        other => other,
    };
    let mut nulls = Vec::with_capacity(rows * nelem);
    let mut doubles = Vec::new();
    let mut longs = Vec::new();
    for i in 0..rows * nelem {
        /* A scalar source fills every element; a source that is already an
        array is re-dimensioned, so ARRAY(MATRIX,6) lays its six elements out
        flat rather than repeating the first. */
        let (x, null) = src.get(i, nelem);
        nulls.push(null);
        if sort == ValueSort::Double {
            doubles.push(x);
        } else {
            longs.push(src.get_i64(i, nelem).0);
        }
    }
    let data = if sort == ValueSort::Double {
        ArrayData::Double(doubles)
    } else {
        ArrayData::Long(longs)
    };
    Ok(ColumnarValue::Array(
        Array::with_nulls(data, nulls)
            .with_nelem(nelem)
            .with_naxes(naxes),
    ))
}

/// What an n-ary numeric kernel produces for one element.
enum NaryOut {
    Double(f64),
    Boolean(bool),
}

/// A function of several numeric arguments, applied elementwise.
///
/// These all share the rule that a row is undefined if *any* argument is
/// undefined there, and they all broadcast a scalar argument across the batch,
/// which is what the engine's `vector[i] > 1 ? buf[elem] : buf[row]` does.
fn nary(
    args: &[ColumnarValue],
    rows: usize,
    name: &'static str,
    f: impl Fn(&[f64]) -> NaryOut,
) -> Result<ColumnarValue, ValueError> {
    let inputs: Vec<NumericInput> = args
        .iter()
        .map(|v| v.numeric(name))
        .collect::<Result<_, _>>()?;

    /* the widest argument sets how many elements each row has */
    let nelem = args
        .iter()
        .filter_map(|v| match v {
            ColumnarValue::Array(a) => Some(a.nelem()),
            _ => None,
        })
        .max()
        .unwrap_or(1);
    let all_const = args.iter().all(|v| matches!(v, ColumnarValue::Scalar(_)));
    let n = if all_const { 1 } else { rows * nelem };

    let mut doubles = Vec::new();
    let mut bools = Vec::new();
    let mut nulls = Vec::with_capacity(n);
    let mut vals = vec![0.0; inputs.len()];
    for i in 0..n {
        let mut null = false;
        for (slot, input) in vals.iter_mut().zip(&inputs) {
            let (v, u) = input.get(i, nelem);
            *slot = v;
            null |= u;
        }
        nulls.push(null);
        /* a null row still needs a slot, and the engine leaves it untouched
        rather than computing over undefined inputs */
        match if null { NaryOut::Double(0.0) } else { f(&vals) } {
            NaryOut::Double(v) => doubles.push(v),
            NaryOut::Boolean(v) => bools.push(v),
        }
    }
    let data = if doubles.is_empty() {
        ArrayData::Boolean(bools)
    } else {
        ArrayData::Double(doubles)
    };
    if all_const {
        return Ok(if nulls[0] {
            ColumnarValue::Null(match data {
                ArrayData::Boolean(_) => ValueSort::Boolean,
                _ => ValueSort::Double,
            })
        } else {
            ColumnarValue::Scalar(match data {
                ArrayData::Boolean(v) => Scalar::Boolean(v[0]),
                ArrayData::Double(v) => Scalar::Double(v[0]),
                _ => unreachable!("nary produces only doubles and booleans"),
            })
        });
    }
    Ok(ColumnarValue::Array(
        Array::with_nulls(data, nulls).with_nelem(nelem),
    ))
}

/// `STRSTR(a, b)`. Absence is a *null*, not a zero, so a caller that wants a
/// number writes `DEFNULL(STRSTR(...), -1)`.
fn str_str(a: &ColumnarValue, b: &ColumnarValue, rows: usize) -> Result<ColumnarValue, ValueError> {
    let (l, r) = (a.text("STRSTR")?, b.text("STRSTR")?);
    let both_const = matches!((a, b), (ColumnarValue::Scalar(_), ColumnarValue::Scalar(_)));
    let n = if both_const { 1 } else { rows };
    let mut out = Vec::with_capacity(n);
    let mut nulls = Vec::with_capacity(n);
    for i in 0..n {
        let ((x, nx), (y, ny)) = (l.get(i), r.get(i));
        match if nx || ny { None } else { strings::find(x, y) } {
            Some(p) => {
                out.push(p);
                nulls.push(false);
            }
            None => {
                out.push(0);
                nulls.push(true);
            }
        }
    }
    if both_const {
        /* `str_pos_const` and `str_pos_rows` disagree: over rows a miss is a
        null, but folded over two constants it is a plain 0. */
        return Ok(ColumnarValue::Scalar(Scalar::Long(out[0])));
    }
    Ok(ColumnarValue::Array(Array::with_nulls(
        ArrayData::Long(out),
        nulls,
    )))
}

/// The width `STRMID` measures a source against: a column's declared width,
/// or a literal's own length when there is no column behind it.
fn src_width(v: &ColumnarValue, text: &[u8]) -> usize {
    match v {
        ColumnarValue::Array(a) if a.nelem() > 0 => a.nelem(),
        _ => text.len(),
    }
}

/// `STRMID(s, pos, len)`. A zero position is undefined, as is a row whose
/// source, position or length is.
fn str_mid(
    s: &ColumnarValue,
    pos: &ColumnarValue,
    len: &ColumnarValue,
    rows: usize,
) -> Result<ColumnarValue, ValueError> {
    let src = s.text("STRMID")?;
    let (p, n) = (pos.numeric("STRMID")?, len.numeric("STRMID")?);
    let mut out = Vec::with_capacity(rows);
    let mut nulls = Vec::with_capacity(rows);
    for i in 0..rows {
        let (text, ntext) = src.get(i);
        let (pos_v, npos) = p.get_i64(i, 1);
        let (len_v, nlen) = n.get_i64(i, 1);
        /* `pos == 0` asks for the null string, which the engine reports as
        undefined rather than as an empty one */
        if ntext || npos || nlen || pos_v <= 0 {
            out.push(Vec::new());
            nulls.push(true);
            continue;
        }
        out.push(strings::mid(
            text,
            src_width(s, text),
            pos_v as usize,
            len_v.max(0) as usize,
        ));
        nulls.push(false);
    }
    Ok(ColumnarValue::Array(Array::with_nulls(
        ArrayData::Str(out),
        nulls,
    )))
}

/// Whether a value is a bit string, which decides between the numeric and the
/// bit kernels for the operators the two share.
fn is_bits(v: &ColumnarValue) -> bool {
    match v {
        ColumnarValue::Scalar(Scalar::Bits(_)) => true,
        ColumnarValue::Null(s) => *s == ValueSort::Bits,
        ColumnarValue::Array(a) => matches!(a.data(), ArrayData::Bits(_)),
        _ => false,
    }
}

/// Whether a value is a string, which decides between the numeric and the
/// string kernels for the operators the two share.
fn is_text(v: &ColumnarValue) -> bool {
    match v {
        ColumnarValue::Scalar(Scalar::Str(_)) => true,
        ColumnarValue::Null(s) => *s == ValueSort::String,
        ColumnarValue::Array(a) => matches!(a.data(), ArrayData::Str(_)),
        _ => false,
    }
}

/// `+` over two strings concatenates them; nothing else is defined.
fn str_arith(
    op: Arith,
    lhs: &ColumnarValue,
    rhs: &ColumnarValue,
    rows: usize,
) -> Result<ColumnarValue, ValueError> {
    if op != Arith::Add {
        return Err(ValueError::BadSort("string operator", ValueSort::String));
    }
    let (l, r) = (lhs.text("string operator")?, rhs.text("string operator")?);
    if let (ColumnarValue::Scalar(_), ColumnarValue::Scalar(_)) = (lhs, rhs) {
        return Ok(ColumnarValue::Scalar(Scalar::Str(strings::concat(
            l.get(0).0,
            r.get(0).0,
        ))));
    }
    let mut out = Vec::with_capacity(rows);
    let mut nulls = Vec::with_capacity(rows);
    for i in 0..rows {
        let ((a, na), (b, nb)) = (l.get(i), r.get(i));
        out.push(strings::concat(a, b));
        nulls.push(na || nb);
    }
    Ok(ColumnarValue::Array(Array::with_nulls(
        ArrayData::Str(out),
        nulls,
    )))
}

/// The comparisons over two strings. Either operand being undefined makes the
/// answer undefined, as `Do_BinOp_str` does.
fn str_compare(
    op: Compare,
    lhs: &ColumnarValue,
    rhs: &ColumnarValue,
    rows: usize,
) -> Result<ColumnarValue, ValueError> {
    let (l, r) = (
        lhs.text("string comparison")?,
        rhs.text("string comparison")?,
    );
    let test = |a: &[u8], b: &[u8]| {
        let ord = strings::compare(a, b);
        match op {
            /* `~` over strings is a parse error, so Approx never runs */
            Compare::Eq | Compare::Approx => ord.is_eq(),
            Compare::Ne => ord.is_ne(),
            Compare::Lt => ord.is_lt(),
            Compare::Lte => ord.is_le(),
            Compare::Gt => ord.is_gt(),
            Compare::Gte => ord.is_ge(),
        }
    };
    if let (ColumnarValue::Scalar(_), ColumnarValue::Scalar(_)) = (lhs, rhs) {
        return Ok(ColumnarValue::Scalar(Scalar::Boolean(test(
            l.get(0).0,
            r.get(0).0,
        ))));
    }
    let mut out = Vec::with_capacity(rows);
    let mut nulls = Vec::with_capacity(rows);
    for i in 0..rows {
        let ((a, na), (b, nb)) = (l.get(i), r.get(i));
        out.push(test(a, b));
        nulls.push(na || nb);
    }
    Ok(ColumnarValue::Array(Array::with_nulls(
        ArrayData::Boolean(out),
        nulls,
    )))
}

/// `!` over a bit string complements each bit rather than testing the value.
fn bit_not(v: &ColumnarValue, rows: usize) -> Result<ColumnarValue, ValueError> {
    let t = v.text("bit not")?;
    if let ColumnarValue::Scalar(_) = v {
        return Ok(ColumnarValue::Scalar(Scalar::Bits(bits::not(t.get(0).0))));
    }
    let out: Vec<Vec<u8>> = (0..rows).map(|i| bits::not(t.get(i).0)).collect();
    Ok(ColumnarValue::Array(Array::new(ArrayData::Bits(out))))
}

/// `&`, `|` and `+` over two bit strings, per `Do_BinOp_bit`.
fn bit_arith(
    op: Arith,
    lhs: &ColumnarValue,
    rhs: &ColumnarValue,
    rows: usize,
) -> Result<ColumnarValue, ValueError> {
    let f: fn(&[u8], &[u8]) -> Vec<u8> = match op {
        Arith::BitAnd => bits::and,
        Arith::BitOr => bits::or,
        Arith::Add => bits::concat,
        _ => return Err(ValueError::BadSort("bit operator", ValueSort::Bits)),
    };
    let (l, r) = (lhs.text("bit operator")?, rhs.text("bit operator")?);

    /* two constants fold, as the engine folds them at parse time */
    if let (ColumnarValue::Scalar(_), ColumnarValue::Scalar(_)) = (lhs, rhs) {
        return Ok(ColumnarValue::Scalar(Scalar::Bits(f(
            l.get(0).0,
            r.get(0).0,
        ))));
    }
    let out: Vec<Vec<u8>> = (0..rows).map(|i| f(l.get(i).0, r.get(i).0)).collect();
    Ok(ColumnarValue::Array(Array::new(ArrayData::Bits(out))))
}

/// The comparisons over two bit strings. Bit strings carry no undef flags, so
/// the result is never null.
fn bit_compare(
    op: Compare,
    lhs: &ColumnarValue,
    rhs: &ColumnarValue,
    rows: usize,
) -> Result<ColumnarValue, ValueError> {
    let (l, r) = (lhs.text("bit comparison")?, rhs.text("bit comparison")?);
    let test = |a: &[u8], b: &[u8]| match op {
        Compare::Eq => bits::cmp_eq(a, b),
        Compare::Ne => !bits::cmp_eq(a, b),
        Compare::Lt => bits::cmp_ord(a, OpCode::Lt, b),
        Compare::Lte => bits::cmp_ord(a, OpCode::Lte, b),
        Compare::Gt => bits::cmp_ord(a, OpCode::Gt, b),
        Compare::Gte => bits::cmp_ord(a, OpCode::Gte, b),
        /* `BITS ~ BITS` is a parse error, so this never runs */
        Compare::Approx => bits::cmp_eq(a, b),
    };
    if let (ColumnarValue::Scalar(_), ColumnarValue::Scalar(_)) = (lhs, rhs) {
        return Ok(ColumnarValue::Scalar(Scalar::Boolean(test(
            l.get(0).0,
            r.get(0).0,
        ))));
    }
    let out: Vec<bool> = (0..rows).map(|i| test(l.get(i).0, r.get(i).0)).collect();
    Ok(ColumnarValue::Array(Array::new(ArrayData::Boolean(out))))
}

/// `base[i, j, ...]` with every axis subscripted, so each row yields one
/// element. `Do_Deref` folds the subscripts innermost-last:
/// `elem = naxes[i] * elem + idx[i] - 1`, and an out-of-range subscript is a
/// range error rather than a null.
fn deref(
    base: &ColumnarValue,
    idx: &[ColumnarValue],
    naxes: &[usize],
    rows: usize,
) -> Result<ColumnarValue, ValueError> {
    let ColumnarValue::Array(a) = base else {
        /* the parser rejects indexing a scalar, so this is a wholly null
        operand: it stays null, one element per row */
        return Ok(ColumnarValue::Null(base.sort()));
    };
    let subs: Vec<NumericInput> = idx
        .iter()
        .map(|v| v.numeric("subscript"))
        .collect::<Result<_, _>>()?;

    let mut picked = Vec::with_capacity(rows);
    for row in 0..rows {
        let mut elem: i64 = 0;
        for (axis, sub) in naxes.iter().zip(&subs).rev() {
            let (i, null) = sub.get_i64(row, 1);
            if null {
                /* an undefined subscript cannot pick an element */
                elem = -1;
                break;
            }
            if i < 1 || i > *axis as i64 {
                return Err(ValueError::OutOfRange);
            }
            elem = *axis as i64 * elem + i - 1;
        }
        picked.push(if elem < 0 {
            None
        } else {
            Some(row * a.nelem() + elem as usize)
        });
    }
    Ok(ColumnarValue::Array(gather_at(a, &picked)))
}

/// Take one element of `a` per row, `None` marking a row whose subscript was
/// undefined. The result holds one element per row.
fn gather_at(a: &Array, picked: &[Option<usize>]) -> Array {
    let nulls = picked
        .iter()
        .map(|p| p.is_none_or(|i| a.is_null(i)))
        .collect();
    let data = match a.data() {
        ArrayData::Long(d) => ArrayData::Long(pick(picked, d, 0)),
        ArrayData::Double(d) => ArrayData::Double(pick(picked, d, 0.0)),
        ArrayData::Boolean(d) => ArrayData::Boolean(pick(picked, d, false)),
        ArrayData::Str(d) => ArrayData::Str(pick(picked, d, Vec::new())),
        ArrayData::Bits(d) => ArrayData::Bits(pick(picked, d, Vec::new())),
    };
    Array::with_nulls(data, nulls)
}

fn pick<T: Clone>(picked: &[Option<usize>], src: &[T], fill: T) -> Vec<T> {
    picked
        .iter()
        .map(|p| p.and_then(|i| src.get(i)).cloned().unwrap_or(fill.clone()))
        .collect()
}

/// Build `{ a, b, c }`: each row gets one element per entry, in order.
fn vector(items: &[ColumnarValue], rows: usize) -> Result<ColumnarValue, ValueError> {
    if items.is_empty() {
        return Err(ValueError::BadSort("{}", ValueSort::Long));
    }
    /* the widest entry decides the sort, as `Close_Vec` did */
    let out = match items.iter().map(|v| v.sort()).max().unwrap() {
        ValueSort::Boolean if items.iter().all(|v| v.sort() == ValueSort::Boolean) => {
            ValueSort::Boolean
        }
        ValueSort::Boolean => ValueSort::Long,
        other => other,
    };

    /* an entry that is itself a vector contributes all of its elements */
    let widths: Vec<usize> = items
        .iter()
        .map(|v| match v {
            ColumnarValue::Array(a) => a.nelem(),
            _ => 1,
        })
        .collect();
    let out_nelem: usize = widths.iter().sum();

    let mut vals = vec![0.0f64; rows * out_nelem];
    let mut nulls = vec![false; rows * out_nelem];
    let mut any = false;
    for (item, (idx, width)) in items.iter().zip(
        widths
            .iter()
            .scan(0usize, |acc, w| {
                let start = *acc;
                *acc += w;
                Some((start, *w))
            })
            .collect::<Vec<_>>(),
    ) {
        let input = item.numeric("{}")?;
        for r in 0..rows {
            for k in 0..width {
                let (v, null) = if item.is_scalar() {
                    input.get(0, width)
                } else {
                    input.get(r * width + k, width)
                };
                let at = r * out_nelem + idx + k;
                vals[at] = v;
                nulls[at] = null;
                any |= null;
            }
        }
    }
    Ok(kernel::build(out, vals, nulls, any, out_nelem))
}

/// Fold each row of `v` down to a single element.
fn reduce(func: Func, v: &ColumnarValue) -> Result<ColumnarValue, ValueError> {
    /* Over a bit string a reduction counts bits rather than folding elements:
    SUM is the number of set bits and NVALID the width, every bit being
    defined. The rest are parse errors over a bit string. */
    if is_bits(v) {
        let t = v.text("reduction")?;
        return match func {
            Func::Sum => Ok(ColumnarValue::Scalar(Scalar::Long(bits::count_ones(
                t.get(0).0,
            )))),
            Func::NValid => Ok(ColumnarValue::Scalar(Scalar::Long(t.get(0).0.len() as i64))),
            _ => Err(ValueError::BadSort("reduction", ValueSort::Bits)),
        };
    }
    /* A string is one element per row, so the only reduction defined over it
    counts whether that element is there. */
    if is_text(v) {
        if func != Func::NValid {
            return Err(ValueError::BadSort("reduction", ValueSort::String));
        }
        let Some(n) = v.len() else {
            return Ok(ColumnarValue::Scalar(Scalar::Long(1)));
        };
        let out = (0..n)
            .map(|i| i64::from(!v.text("NVALID").map(|t| t.get(i).1).unwrap_or(true)))
            .collect();
        return Ok(ColumnarValue::Array(Array::new(ArrayData::Long(out))));
    }
    let out = match func {
        Func::Average | Func::Stddev => ValueSort::Double,
        Func::NValid => ValueSort::Long,
        _ => match v.sort() {
            ValueSort::Boolean => ValueSort::Long,
            other => other,
        },
    };
    let input = v.numeric("reduction")?;

    let Some(n) = v.len() else {
        /* A constant is a row of one defined element, which is what the
        engine's `*_const` kernels compute. */
        let (x, null) = input.get(0, 1);
        return Ok(match func {
            Func::NValid => ColumnarValue::Scalar(Scalar::Long(i64::from(!null))),
            Func::Stddev => ColumnarValue::Scalar(Scalar::Double(0.0)),
            _ if null => ColumnarValue::Null(out),
            Func::Average => ColumnarValue::Scalar(Scalar::Double(x)),
            _ => ColumnarValue::Scalar(kernel::scalar_of_sort(out, x)),
        });
    };

    let nelem = input.nelem().max(1);
    let rows = n / nelem;
    let mut vals = vec![0.0f64; rows];
    let mut nulls = vec![false; rows];
    let mut any = false;
    let mut defined: Vec<f64> = Vec::with_capacity(nelem);

    for r in 0..rows {
        defined.clear();
        for k in 0..nelem {
            let (x, null) = input.get(r * nelem + k, nelem);
            if !null {
                defined.push(x);
            }
        }
        let count = defined.len();
        let (val, null) = match func {
            Func::NValid => (count as f64, false),
            Func::Stddev => {
                if count > 1 {
                    let mean = defined.iter().sum::<f64>() / count as f64;
                    let ss: f64 = defined.iter().map(|x| (x - mean) * (x - mean)).sum();
                    ((ss / (count as f64 - 1.0)).sqrt(), false)
                } else {
                    (0.0, false)
                }
            }
            _ if count == 0 => (0.0, true),
            Func::Sum => (defined.iter().sum(), false),
            Func::Average => (defined.iter().sum::<f64>() / count as f64, false),
            Func::Min1 => (defined.iter().copied().fold(f64::INFINITY, f64::min), false),
            Func::Max1 => (
                defined.iter().copied().fold(f64::NEG_INFINITY, f64::max),
                false,
            ),
            Func::Median => {
                let mut scratch = defined.clone();
                (crate::eval_y::qselect_median(&mut scratch), false)
            }
            _ => unreachable!("not a reduction"),
        };
        vals[r] = val;
        nulls[r] = null;
        any |= null;
    }
    /* one element per row now */
    Ok(kernel::build(out, vals, nulls, any, 1))
}

/// `ABS`, which keeps its argument's sort rather than widening to double.
fn map1_keep_sort(
    v: &ColumnarValue,
    f: impl Fn(f64) -> f64,
    g: impl Fn(c_long) -> c_long,
) -> Result<ColumnarValue, ValueError> {
    if v.sort() == ValueSort::Double {
        kernel::map_double(v, "ABS", f)
    } else {
        kernel::map_long(v, "ABS", g)
    }
}

fn is_null(v: &ColumnarValue) -> Result<ColumnarValue, ValueError> {
    Ok(match v {
        ColumnarValue::Null(_) => ColumnarValue::Scalar(Scalar::Boolean(true)),
        ColumnarValue::Scalar(_) => ColumnarValue::Scalar(Scalar::Boolean(false)),
        ColumnarValue::Array(a) => {
            let flags: Vec<bool> = (0..a.len()).map(|i| a.is_null(i)).collect();
            ColumnarValue::Array(Array::new(ArrayData::Boolean(flags)).with_nelem(a.nelem()))
        }
    })
}

fn def_null(v: &ColumnarValue, fallback: &ColumnarValue) -> Result<ColumnarValue, ValueError> {
    let out = v.sort().max(fallback.sort());
    let a = v.numeric("DEFNULL")?;
    let b = fallback.numeric("DEFNULL")?;
    let out_nelem = a.nelem().max(b.nelem());

    let Some(n) = v.broadcast_len(fallback)? else {
        let (x, xn) = a.get(0, out_nelem);
        let (y, yn) = b.get(0, out_nelem);
        let (val, null) = if xn { (y, yn) } else { (x, false) };
        return Ok(if null {
            ColumnarValue::Null(out)
        } else {
            ColumnarValue::Scalar(kernel::scalar_of_sort(out, val))
        });
    };

    let mut vals = vec![0.0f64; n];
    let mut nulls = vec![false; n];
    let mut any = false;
    for i in 0..n {
        let (x, xn) = a.get(i, out_nelem);
        let (y, yn) = b.get(i, out_nelem);
        let (val, null) = if xn { (y, yn) } else { (x, false) };
        vals[i] = val;
        nulls[i] = null;
        any |= null;
    }
    Ok(kernel::build(out, vals, nulls, any, out_nelem))
}

fn set_null(v: &ColumnarValue, sentinel: &ColumnarValue) -> Result<ColumnarValue, ValueError> {
    let out = v.sort();
    let a = v.numeric("SETNULL")?;
    let s = sentinel.numeric("SETNULL")?;
    let out_nelem = a.nelem();

    let Some(n) = v.len() else {
        /* Bug-compatible with the engine: `set_null_const` copies the value
        and never consults the sentinel, so `SETNULL(1,1)` is 1 rather than
        null. CFITSIO does the same -- see the corpus. */
        let (x, xn) = a.get(0, out_nelem);
        return Ok(if xn {
            ColumnarValue::Null(out)
        } else {
            ColumnarValue::Scalar(kernel::scalar_of_sort(out, x))
        });
    };

    let mut vals = vec![0.0f64; n];
    let mut nulls = vec![false; n];
    let mut any = false;
    let (t, _) = s.get(0, 1);
    for i in 0..n {
        let (x, xn) = a.get(i, out_nelem);
        let null = xn || x == t;
        vals[i] = if null { 0.0 } else { x };
        nulls[i] = null;
        any |= null;
    }
    Ok(kernel::build(out, vals, nulls, any, out_nelem))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn batch(cols: Vec<ColumnBatch>, n_rows: c_long) -> Batch {
        Batch {
            columns: cols,
            n_rows,
            first_data_row: 1,
            n_data_rows: n_rows,
            total_rows: n_rows,
            first_row: 1,
            accum: RefCell::new(vec![AccumState::default(); 4]),
        }
    }

    fn str_col(v: &[&str], width: c_long) -> ColumnBatch {
        ColumnBatch {
            data: ArrayData::Str(v.iter().map(|s| s.as_bytes().to_vec()).collect()),
            nulls: Vec::new(),
            nelem: width,
            naxes: Vec::new(),
        }
    }

    fn slit(s: &str) -> Expr {
        Expr::Literal(Scalar::Str(s.as_bytes().to_vec()))
    }

    #[test]
    fn strings_compare_and_concatenate_over_a_column() {
        let b = batch(vec![str_col(&["alpha", "beta"], 10)], 2);
        let e = Expr::Equality {
            negated: false,
            lhs: Box::new(Expr::Column(0)),
            rhs: Box::new(slit("alpha")),
        };
        let ColumnarValue::Array(a) = e.evaluate(&b).unwrap() else {
            panic!("expected an array")
        };
        assert_eq!(a.data(), &ArrayData::Boolean(vec![true, false]));

        let e = Expr::Arith {
            op: Arith::Add,
            lhs: Box::new(Expr::Column(0)),
            rhs: Box::new(slit("!")),
        };
        let ColumnarValue::Array(a) = e.evaluate(&b).unwrap() else {
            panic!("expected an array")
        };
        assert_eq!(
            a.data(),
            &ArrayData::Str(vec![b"alpha!".to_vec(), b"beta!".to_vec()])
        );
    }

    #[test]
    fn strstr_reports_a_miss_as_a_null_over_rows() {
        let b = batch(vec![str_col(&["alpha", "beta"], 10)], 2);
        let e = Expr::Call {
            func: Func::StrStr,
            args: vec![Expr::Column(0), slit("ta")],
        };
        let ColumnarValue::Array(a) = e.evaluate(&b).unwrap() else {
            panic!("expected an array")
        };
        /* 'ta' is in 'beta' at 3 and absent from 'alpha' */
        assert!(a.is_null(0), "a miss is undefined, not zero");
        assert!(!a.is_null(1));
        assert_eq!(a.data(), &ArrayData::Long(vec![0, 3]));
    }

    #[test]
    fn strstr_folded_over_two_constants_reports_a_miss_as_zero() {
        /* str_pos_const differs from str_pos_rows here, and callers depend on
        the constant form being a plain 0 */
        let b = batch(vec![], 2);
        let e = Expr::Call {
            func: Func::StrStr,
            args: vec![slit("abc"), slit("z")],
        };
        assert_eq!(
            e.evaluate(&b).unwrap(),
            ColumnarValue::Scalar(Scalar::Long(0))
        );
    }

    #[test]
    fn strmid_measures_against_the_declared_width() {
        let b = batch(vec![str_col(&["alpha", "beta"], 10)], 2);
        let e = Expr::Call {
            func: Func::StrMid,
            args: vec![
                Expr::Column(0),
                Expr::Literal(Scalar::Long(1)),
                Expr::Literal(Scalar::Long(3)),
            ],
        };
        let ColumnarValue::Array(a) = e.evaluate(&b).unwrap() else {
            panic!("expected an array")
        };
        assert_eq!(
            a.data(),
            &ArrayData::Str(vec![b"alp".to_vec(), b"bet".to_vec()])
        );
    }

    #[test]
    fn a_zero_position_makes_strmid_undefined() {
        let b = batch(vec![str_col(&["alpha"], 10)], 1);
        let e = Expr::Call {
            func: Func::StrMid,
            args: vec![
                Expr::Column(0),
                Expr::Literal(Scalar::Long(0)),
                Expr::Literal(Scalar::Long(3)),
            ],
        };
        let ColumnarValue::Array(a) = e.evaluate(&b).unwrap() else {
            panic!("expected an array")
        };
        assert!(a.is_null(0));
    }

    #[test]
    fn a_conditional_can_pick_between_strings() {
        let b = batch(vec![str_col(&["alpha", "beta"], 10)], 2);
        let e = Expr::IfThenElse {
            cond: Box::new(Expr::Equality {
                negated: false,
                lhs: Box::new(Expr::Column(0)),
                rhs: Box::new(slit("alpha")),
            }),
            then: Box::new(slit("big")),
            els: Box::new(slit("small")),
        };
        let ColumnarValue::Array(a) = e.evaluate(&b).unwrap() else {
            panic!("expected an array")
        };
        assert_eq!(
            a.data(),
            &ArrayData::Str(vec![b"big".to_vec(), b"small".to_vec()])
        );
    }

    /// A batch whose loaded chunk is wider than the rows being evaluated, so
    /// an offset can reach outside the batch but stay inside the chunk.
    fn offset_batch(first_row: c_long, n_rows: c_long, total: c_long) -> Batch {
        Batch {
            columns: vec![long_col(&[10, 20, 30, 40, 50])],
            n_rows,
            first_data_row: 1,
            n_data_rows: 5,
            total_rows: total,
            first_row,
            accum: RefCell::new(Vec::new()),
        }
    }

    #[test]
    fn a_zero_offset_is_the_plain_column() {
        let b = offset_batch(2, 3, 5);
        let e = Expr::Offset { col: 0, offset: 0 };
        assert_eq!(longs(&e.evaluate(&b).unwrap()), vec![20, 30, 40]);
        assert_eq!(
            longs(&Expr::Column(0).evaluate(&b).unwrap()),
            vec![20, 30, 40]
        );
    }

    #[test]
    fn an_offset_reads_the_neighbouring_rows() {
        let b = offset_batch(2, 3, 5);
        let back = Expr::Offset { col: 0, offset: -1 };
        assert_eq!(longs(&back.evaluate(&b).unwrap()), vec![10, 20, 30]);
        let fwd = Expr::Offset { col: 0, offset: 1 };
        assert_eq!(longs(&fwd.evaluate(&b).unwrap()), vec![30, 40, 50]);
    }

    #[test]
    fn a_row_off_either_end_of_the_table_is_undefined() {
        let b = offset_batch(1, 3, 5);
        let back = Expr::Offset { col: 0, offset: -1 };
        let ColumnarValue::Array(a) = back.evaluate(&b).unwrap() else {
            panic!("expected an array")
        };
        assert!(a.is_null(0), "row 0 is before the table");
        assert!(!a.is_null(1));

        /* and past the last row, with the chunk covering the whole table */
        let b = offset_batch(3, 3, 5);
        let fwd = Expr::Offset { col: 0, offset: 2 };
        let ColumnarValue::Array(a) = fwd.evaluate(&b).unwrap() else {
            panic!("expected an array")
        };
        assert!(!a.is_null(0), "row 5 is the last one");
        assert!(a.is_null(1), "row 6 is past the table");
        assert!(a.is_null(2));
    }

    #[test]
    fn a_row_outside_the_loaded_chunk_asks_for_a_reload() {
        /* the chunk holds rows 3..5 of a 10-row table, so reaching back to
        row 2 is a row the batch cannot serve and the arena must */
        let b = Batch {
            columns: vec![long_col(&[30, 40, 50])],
            n_rows: 3,
            first_data_row: 3,
            n_data_rows: 3,
            total_rows: 10,
            first_row: 3,
            accum: RefCell::new(Vec::new()),
        };
        let e = Expr::Offset { col: 0, offset: -1 };
        assert!(matches!(e.evaluate(&b), Err(ValueError::NeedsReload)));
    }

    #[test]
    fn an_unset_total_row_count_does_not_null_the_column() {
        /* an image with no NAXIS2 leaves totalRows at 0, which says nothing
        about the top end */
        let b = offset_batch(1, 3, 0);
        let ColumnarValue::Array(a) = Expr::Column(0).evaluate(&b).unwrap() else {
            panic!("expected an array")
        };
        assert!((0..3).all(|i| !a.is_null(i)));
        assert_eq!(
            longs(&Expr::Column(0).evaluate(&b).unwrap()),
            vec![10, 20, 30]
        );
    }

    /// Evaluate `e` over two consecutive batches, sharing one accumulator
    /// vector the way the bridge lends it between calls.
    fn over_two_batches(e: &Expr, first: &[c_long], second: &[c_long]) -> Vec<c_long> {
        let mut state = vec![AccumState::default(); 1];
        let mut out = Vec::new();
        for (i, rows) in [first, second].iter().enumerate() {
            let mut b = batch(vec![long_col(rows)], rows.len() as c_long);
            b.first_row = 1 + i as c_long * first.len() as c_long;
            b.total_rows = (first.len() + second.len()) as c_long;
            b.n_data_rows = rows.len() as c_long;
            b.first_data_row = b.first_row;
            *b.accum.borrow_mut() = core::mem::take(&mut state);
            out.extend(longs(&e.evaluate(&b).unwrap()));
            state = core::mem::take(&mut b.accum.borrow_mut());
        }
        out
    }

    #[test]
    fn accum_runs_a_total_across_the_batch_boundary() {
        let e = Expr::Accum {
            id: 0,
            diff: false,
            arg: Box::new(Expr::Column(0)),
        };
        /* the second batch must continue from 6, not restart */
        assert_eq!(
            over_two_batches(&e, &[1, 2, 3], &[4, 5]),
            vec![1, 3, 6, 10, 15]
        );
    }

    #[test]
    fn seqdiff_remembers_the_last_element_of_the_previous_batch() {
        let e = Expr::Accum {
            id: 0,
            diff: true,
            arg: Box::new(Expr::Column(0)),
        };
        /* the first element differences against zero, and 10 against 3 */
        assert_eq!(
            over_two_batches(&e, &[1, 2, 3], &[10, 20]),
            vec![1, 1, 1, 7, 10]
        );
    }

    #[test]
    fn accum_skips_undefined_elements_and_is_never_null() {
        let mut b = batch(vec![long_col(&[7, 0, 10])], 3);
        b.columns[0].nulls = vec![false, true, false];
        let e = Expr::Accum {
            id: 0,
            diff: false,
            arg: Box::new(Expr::Column(0)),
        };
        let ColumnarValue::Array(a) = e.evaluate(&b).unwrap() else {
            panic!("expected an array")
        };
        assert_eq!(a.data(), &ArrayData::Long(vec![7, 7, 17]));
        assert!((0..3).all(|i| !a.is_null(i)), "ACCUM is never undefined");
    }

    #[test]
    fn seqdiff_is_undefined_on_both_sides_of_a_gap() {
        let mut b = batch(vec![long_col(&[1, 0, 10])], 3);
        b.columns[0].nulls = vec![false, true, false];
        let e = Expr::Accum {
            id: 0,
            diff: true,
            arg: Box::new(Expr::Column(0)),
        };
        let ColumnarValue::Array(a) = e.evaluate(&b).unwrap() else {
            panic!("expected an array")
        };
        assert!(!a.is_null(0));
        assert!(a.is_null(1), "the gap itself");
        assert!(a.is_null(2), "and the element after it");
    }

    fn matrix_batch() -> Batch {
        let mut b = batch(vec![matrix_col()], 2);
        b.total_rows = 2;
        b.n_data_rows = 2;
        b.first_data_row = 1;
        b
    }

    #[test]
    fn elementnum_counts_within_each_row() {
        let b = matrix_batch();
        let e = Expr::Call {
            func: Func::ElementNum,
            args: vec![Expr::Column(0)],
        };
        assert_eq!(
            longs(&e.evaluate(&b).unwrap()),
            vec![1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 6]
        );
    }

    #[test]
    fn axiselem_cycles_the_innermost_axis_fastest() {
        let b = matrix_batch();
        let axis = |n: c_long| Expr::Call {
            func: Func::AxisElem,
            args: vec![Expr::Column(0), Expr::Literal(Scalar::Long(n))],
        };
        /* a 2x3 column: axis 1 runs 1,2 repeatedly, axis 2 changes every two */
        assert_eq!(
            longs(&axis(1).evaluate(&b).unwrap())[..6],
            [1, 2, 1, 2, 1, 2]
        );
        assert_eq!(
            longs(&axis(2).evaluate(&b).unwrap())[..6],
            [1, 1, 2, 2, 3, 3]
        );
    }

    #[test]
    fn array_fills_a_scalar_and_redimensions_a_vector() {
        let b = matrix_batch();
        /* a scalar reaches every element */
        let e = Expr::Call {
            func: Func::Array,
            args: vec![
                Expr::Literal(Scalar::Long(3)),
                Expr::Literal(Scalar::Long(4)),
            ],
        };
        assert_eq!(
            longs(&e.evaluate(&b).unwrap()),
            vec![3, 3, 3, 3, 3, 3, 3, 3]
        );

        /* a source that is already an array is laid out flat, not repeated */
        let e = Expr::Call {
            func: Func::Array,
            args: vec![Expr::Column(0), Expr::Literal(Scalar::Long(6))],
        };
        assert_eq!(longs(&e.evaluate(&b).unwrap())[..6], [1, 2, 3, 4, 5, 6]);
    }

    #[test]
    fn array_takes_its_shape_from_a_vector_of_dimensions() {
        let b = matrix_batch();
        let e = Expr::Call {
            func: Func::Array,
            args: vec![
                Expr::Literal(Scalar::Long(0)),
                Expr::Vector(vec![
                    Expr::Literal(Scalar::Long(2)),
                    Expr::Literal(Scalar::Long(3)),
                    Expr::Literal(Scalar::Long(1)),
                ]),
            ],
        };
        let ColumnarValue::Array(a) = e.evaluate(&b).unwrap() else {
            panic!("expected an array")
        };
        assert_eq!(a.nelem(), 6);
        assert_eq!(a.naxes(), vec![2, 3, 1]);
    }

    /// A 2x3 column holding 1..6 in the first row and 11..16 in the second,
    /// laid out innermost axis first as the engine stores it.
    fn matrix_col() -> ColumnBatch {
        ColumnBatch {
            data: ArrayData::Long((1..=6).chain(11..=16).collect()),
            nulls: Vec::new(),
            nelem: 6,
            naxes: vec![2, 3],
        }
    }

    #[test]
    fn a_full_subscript_picks_one_element_per_row() {
        let b = batch(vec![matrix_col()], 2);
        /* M[1,1] is the first element of each row, M[2,3] the last */
        for (idx, want) in [([1, 1], [1, 11]), ([2, 3], [6, 16])] {
            let e = Expr::Deref {
                base: Box::new(Expr::Column(0)),
                idx: idx
                    .iter()
                    .map(|&i| Expr::Literal(Scalar::Long(i)))
                    .collect(),
                naxes: vec![2, 3],
            };
            let got = e.evaluate(&b).unwrap();
            let ColumnarValue::Array(a) = got else {
                panic!("expected an array")
            };
            assert_eq!(a.data(), &ArrayData::Long(want.to_vec()));
            assert_eq!(a.len(), 2, "one element per row");
        }
    }

    #[test]
    fn one_subscript_on_a_two_dimensional_column_selects_a_slice() {
        let b = matrix_batch();
        /* a 2x3 column: M[1] is its first pair, M[3] its last */
        let slice = |i: c_long| Expr::Slice {
            base: Box::new(Expr::Column(0)),
            index: i,
            run: 2,
            axis_len: 3,
            naxes: vec![2],
        };
        let ColumnarValue::Array(a) = slice(1).evaluate(&b).unwrap() else {
            panic!("expected an array")
        };
        assert_eq!(a.data(), &ArrayData::Long(vec![1, 2, 11, 12]));
        assert_eq!(a.nelem(), 2, "the slice keeps the remaining axis");
        assert_eq!(longs(&slice(3).evaluate(&b).unwrap()), vec![5, 6, 15, 16]);
    }

    #[test]
    fn a_slice_index_is_checked_against_the_outermost_axis() {
        let b = matrix_batch();
        let slice = |i: c_long| Expr::Slice {
            base: Box::new(Expr::Column(0)),
            index: i,
            run: 2,
            axis_len: 3,
            naxes: vec![2],
        };
        for bad in [0, 4, -1] {
            assert!(matches!(
                slice(bad).evaluate(&b),
                Err(ValueError::OutOfRange)
            ));
        }
    }

    #[test]
    fn a_subscript_outside_the_axis_is_a_range_error() {
        let b = batch(vec![matrix_col()], 2);
        for idx in [[3, 1], [1, 4], [0, 0]] {
            let e = Expr::Deref {
                base: Box::new(Expr::Column(0)),
                idx: idx
                    .iter()
                    .map(|&i| Expr::Literal(Scalar::Long(i)))
                    .collect(),
                naxes: vec![2, 3],
            };
            assert!(
                matches!(e.evaluate(&b), Err(ValueError::OutOfRange)),
                "{idx:?} is out of a 2x3 column"
            );
        }
    }

    #[test]
    fn an_undefined_subscript_yields_a_null() {
        let b = batch(vec![matrix_col()], 2);
        let e = Expr::Deref {
            base: Box::new(Expr::Column(0)),
            idx: vec![Expr::Literal(Scalar::Long(1)), Expr::Null(ValueSort::Long)],
            naxes: vec![2, 3],
        };
        let ColumnarValue::Array(a) = e.evaluate(&b).unwrap() else {
            panic!("expected an array")
        };
        assert!((0..2).all(|i| a.is_null(i)));
    }

    #[test]
    fn a_subscript_of_a_vector_literal_reads_that_entry() {
        let b = batch(vec![], 2);
        let e = Expr::Deref {
            base: Box::new(Expr::Vector(vec![
                Expr::Literal(Scalar::Long(7)),
                Expr::Literal(Scalar::Long(9)),
            ])),
            idx: vec![Expr::Literal(Scalar::Long(2))],
            naxes: vec![2],
        };
        let ColumnarValue::Array(a) = e.evaluate(&b).unwrap() else {
            panic!("expected an array")
        };
        assert_eq!(a.data(), &ArrayData::Long(vec![9, 9]));
    }

    fn long_col(v: &[c_long]) -> ColumnBatch {
        ColumnBatch {
            data: ArrayData::Long(v.to_vec()),
            nulls: Vec::new(),
            nelem: 1,
            naxes: Vec::new(),
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
            first_data_row: 11,
            n_data_rows: 3,
            total_rows: 13,
            first_row: 11,
            accum: RefCell::new(Vec::new()),
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

        /* a null condition gives a null for every row: `#NULL` is a rows
        kernel in the engine, not a folded constant */
        let n = Expr::IfThenElse {
            cond: Box::new(Expr::Null(ValueSort::Boolean)),
            then: lit(1),
            els: lit(2),
        };
        match n.evaluate(&b).unwrap() {
            ColumnarValue::Array(a) => assert!((0..3).all(|i| a.is_null(i))),
            other => panic!("{other:?}"),
        }
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
    fn nesting_up_to_the_parsers_limit_evaluates() {
        /* Evaluation recurses once per level, so the depth the parser admits
        is the depth this has to survive. `grammar::MAX_DEPTH` is 100; if that
        is ever raised, this is the test that says whether the evaluator can
        take it. Dropping the tree recurses too, so the margin covers both. */
        const DEPTH: usize = 100;
        let mut e = Expr::Literal(Scalar::Long(0));
        for _ in 0..DEPTH {
            e = Expr::Arith {
                op: Arith::Add,
                lhs: Box::new(e),
                rhs: lit(1),
            };
        }
        let b = batch(vec![], 1);
        assert_eq!(
            e.evaluate(&b).unwrap(),
            ColumnarValue::Scalar(Scalar::Long(DEPTH as c_long))
        );
    }
}
