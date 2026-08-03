//! `Ast` to [`Expr`], for the subset the new evaluator covers.
//!
//! The port is incremental: anything not yet handled returns
//! [`Unsupported`], and the caller keeps the `Node` arena for that expression.
//! That is what lets the new evaluator go in behind a flag and grow one kernel
//! family at a time, with `tests/test_eval_corpus.rs` proving each step.
//!
//! Name resolution has already happened — `parser::resolve` ran before the
//! parse — so a name here is a column index or a folded constant.

use super::expr::{Expr, Func, GtiIntervals, Logic, Unary};
use super::kernel::{Arith, Compare};
use super::value::Scalar;
use crate::c_types::c_long;
use crate::eval_defs::{ColumnSort, ParserValue, ValueSort};
use crate::parser::ast::{Ast, AstKind, BinOp, UnOp};
use crate::parser::resolve::Resolutions;
use crate::parser::token::CallKind;

/// Why an expression could not be lowered to an [`Expr`].
///
/// Carries the construct's name so that a corpus divergence points straight at
/// what is still missing rather than at "something".
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct Unsupported(pub(crate) &'static str);

type Res = Result<Expr, Unsupported>;

/// The constant an offset names, when it is written plainly enough to fold
/// here. Anything else -- `INTCOL{1+1}` -- is left to the arena, which already
/// folded it at parse time.
fn const_long(ast: &Ast) -> Option<c_long> {
    match &ast.kind {
        AstKind::Long(v) => Some(*v),
        AstKind::Unary { op: UnOp::Neg, arg } => const_long(arg).map(|v| -v),
        AstKind::Unary {
            op: UnOp::Plus,
            arg,
        } => const_long(arg),
        /* the arena folded these at parse time, so folding the same way here
        keeps `INTCOL{1+1}` on the new path rather than falling back */
        AstKind::Binary { op, lhs, rhs } => {
            let (a, b) = (const_long(lhs)?, const_long(rhs)?);
            match op {
                BinOp::Add => a.checked_add(b),
                BinOp::Sub => a.checked_sub(b),
                BinOp::Mul => a.checked_mul(b),
                _ => None,
            }
        }
        _ => None,
    }
}

/// What lowering needs to know about the table's columns: the shape of each,
/// for the shape functions, and its sort, for the places where an operator's
/// meaning depends on it (`+` concatenates text but adds numbers).
pub(crate) struct Columns {
    /// What each GTI call read while the arena was built, keyed by the call's
    /// byte offset. The file navigation happens once, there.
    pub(crate) gti: crate::parser::lower::GtiLoads,
    /// Hands out a slot to each `ACCUM`/`SEQDIFF` as it is lowered, so the
    /// running values can live in one flat vector that survives between
    /// batches. Interior mutability keeps `lower` taking `&Columns`.
    pub(crate) accums: core::cell::Cell<usize>,
    /// Element count and axis lengths per column. The count is not always the
    /// product of the axes: a string or bit-string column is one entry per row
    /// laid out on one axis, and its `nelem` is the declared *width*.
    pub(crate) shapes: Vec<(c_long, Vec<c_long>)>,
    pub(crate) sorts: Vec<ValueSort>,
}

impl Columns {
    fn sort(&self, i: usize) -> ValueSort {
        self.sorts.get(i).copied().unwrap_or(ValueSort::Long)
    }
}

/// An expression's shape: how many elements each row holds, and how those are
/// laid out across axes.
///
/// The arena computes this while building nodes and the shape functions read
/// it back off the result node, folding to a constant at parse time. Lowering
/// works from the `Ast`, so it recomputes the same thing from the column
/// cols -- the rules are `New_BinOp`'s (the non-scalar operand's shape wins)
/// plus the string and bit-string `+`, which concatenates and so *adds* the
/// two widths.
#[derive(Clone, Debug, PartialEq, Eq)]
struct Shape {
    nelem: c_long,
    naxes: Vec<c_long>,
}

impl Shape {
    fn scalar() -> Shape {
        Shape {
            nelem: 1,
            naxes: vec![1],
        }
    }

    fn flat(n: c_long) -> Shape {
        Shape {
            nelem: n,
            naxes: vec![n],
        }
    }

    fn naxis(&self) -> c_long {
        self.naxes.len().max(1) as c_long
    }

    /// The axis `n` names, one-based and clamped at both ends the way
    /// `NAXES(MATRIX,0)` and `NAXES(MATRIX,1)` both give the first axis.
    fn axis(&self, n: c_long) -> c_long {
        let i = (n.max(1) - 1) as usize;
        self.naxes.get(i).copied().unwrap_or(1)
    }

    /// The shape a binary operator gives: whichever operand is not a scalar.
    fn combine(a: &Shape, b: &Shape) -> Shape {
        if a.nelem == 1 { b.clone() } else { a.clone() }
    }
}

/// The shape of a lowered expression, or `None` when it cannot be determined
/// statically -- in which case the shape functions over it fall back.
fn shape_of(e: &Expr, cols: &Columns) -> Option<Shape> {
    let col_shape = |i: &usize| -> Option<Shape> {
        let (nelem, naxes) = cols.shapes.get(*i)?;
        Some(Shape {
            nelem: (*nelem).max(1),
            naxes: if naxes.is_empty() {
                vec![(*nelem).max(1)]
            } else {
                naxes.clone()
            },
        })
    };
    Some(match e {
        Expr::Literal(Scalar::Str(v)) | Expr::Literal(Scalar::Bits(v)) => {
            Shape::flat(v.len() as c_long)
        }
        Expr::Literal(_) | Expr::Null(_) | Expr::RowNumber => Shape::scalar(),
        Expr::Column(i) | Expr::Offset { col: i, .. } => col_shape(i)?,
        /* a vector literal expands any entry that is itself a vector, so the
        element counts add rather than the entries being counted */
        Expr::Vector(items) => Shape::flat(
            items
                .iter()
                .map(|i| shape_of(i, cols).map(|s| s.nelem))
                .sum::<Option<c_long>>()?,
        ),
        /* a fully-indexed subscript picks one element */
        Expr::Deref { .. } => Shape::scalar(),
        Expr::Slice { run, naxes, .. } => Shape {
            nelem: *run as c_long,
            naxes: naxes.iter().map(|&n| n as c_long).collect(),
        },
        Expr::Unary { arg, .. } => shape_of(arg, cols)?,
        /* one answer per element of the time expression */
        Expr::Gti { time, .. } => shape_of(time, cols)?,
        /* ACCUM and SEQDIFF are elementwise, so the shape carries through */
        Expr::Accum { arg, .. } => shape_of(arg, cols)?,
        Expr::Arith { op, lhs, rhs } => {
            let (a, b) = (shape_of(lhs, cols)?, shape_of(rhs, cols)?);
            /* `+` over text concatenates, so the widths add */
            if *op == Arith::Add && is_text_expr(lhs, cols) && is_text_expr(rhs, cols) {
                Shape::flat(a.nelem + b.nelem)
            } else {
                Shape::combine(&a, &b)
            }
        }
        Expr::Compare { lhs, rhs, .. } | Expr::Equality { lhs, rhs, .. } => {
            Shape::combine(&shape_of(lhs, cols)?, &shape_of(rhs, cols)?)
        }
        Expr::Logic { lhs, rhs, .. } => {
            Shape::combine(&shape_of(lhs, cols)?, &shape_of(rhs, cols)?)
        }
        Expr::IfThenElse { then, els, .. } => {
            Shape::combine(&shape_of(then, cols)?, &shape_of(els, cols)?)
        }
        /* the reductions collapse a row to one element; the elementwise
        functions keep their argument's shape */
        Expr::Call { func, args } => match func {
            Func::Sum
            | Func::Average
            | Func::Stddev
            | Func::Median
            | Func::NValid
            | Func::Min1
            | Func::Max1
            | Func::StrStr => Shape::scalar(),
            /* a generator with no argument yields one value per row */
            Func::Random | Func::RandomN if args.is_empty() => Shape::scalar(),
            /* ARRAY's shape is its second argument, which the parser requires
            to be constant, so it is read off the lowered dims directly */
            Func::Array => {
                let dims = const_dims(args.get(1)?)?;
                Shape {
                    nelem: dims.iter().product::<c_long>().max(1),
                    naxes: dims,
                }
            }
            _ => shape_of(args.first()?, cols)?,
        },
    })
}

/// `NELEM`, `NAXIS` and `NAXES` answered from the argument's shape.
///
/// `Ok(None)` means this is not one of them. Any other failure -- an argument
/// that will not lower, a shape that is not statically known, an axis index
/// that is not a constant -- returns `Err` so the expression falls back.
fn fold_shape_fn(
    name: &[u8],
    args: &[Ast],
    names: &Resolutions,
    cols: &Columns,
) -> Result<Option<Expr>, Unsupported> {
    let arity = match name {
        b"NELEM" | b"NAXIS" => 1,
        b"NAXES" => 2,
        _ => return Ok(None),
    };
    if args.len() != arity {
        /* a wrong arity is the arena's error to report */
        return Err(Unsupported("shape function with the wrong arity"));
    }
    let arg = lower(&args[0], names, cols)?;
    let shape = shape_of(&arg, cols).ok_or(Unsupported("shape function on an unknown shape"))?;

    let v = match name {
        b"NELEM" => shape.nelem,
        b"NAXIS" => shape.naxis(),
        _ => {
            let n = const_long(&args[1]).ok_or(Unsupported("NAXES axis is not a constant"))?;
            shape.axis(n)
        }
    };
    Ok(Some(Expr::Literal(Scalar::Long(v))))
}

/// The axis lengths `ARRAY`'s second argument names: a single count, or a
/// vector literal of counts. Anything else is not a shape known here.
fn const_dims(e: &Expr) -> Option<Vec<c_long>> {
    match e {
        Expr::Literal(Scalar::Long(n)) => Some(vec![*n]),
        Expr::Literal(Scalar::Double(n)) => Some(vec![*n as c_long]),
        Expr::Vector(items) => items
            .iter()
            .map(|i| match i {
                Expr::Literal(Scalar::Long(n)) => Some(*n),
                Expr::Literal(Scalar::Double(n)) => Some(*n as c_long),
                _ => None,
            })
            .collect(),
        _ => None,
    }
}

/// Whether an expression is textual, which decides whether `+` concatenates.
fn is_text_expr(e: &Expr, cols: &Columns) -> bool {
    match e {
        Expr::Literal(Scalar::Str(_) | Scalar::Bits(_)) => true,
        Expr::Column(i) | Expr::Offset { col: i, .. } => {
            matches!(cols.sort(*i), ValueSort::String | ValueSort::Bits)
        }
        Expr::Arith { lhs, .. } => is_text_expr(lhs, cols),
        Expr::IfThenElse { then, .. } => is_text_expr(then, cols),
        Expr::Call { func, args } => {
            *func == Func::StrMid || (*func == Func::DefNull && is_text_expr(&args[0], cols))
        }
        _ => false,
    }
}

/// Lower a whole expression, or report the first construct not yet ported.
pub(crate) fn lower(ast: &Ast, names: &Resolutions, cols: &Columns) -> Res {
    match &ast.kind {
        AstKind::Long(v) => Ok(Expr::Literal(Scalar::Long(*v))),
        AstKind::Double(v) => Ok(Expr::Literal(Scalar::Double(*v))),
        AstKind::Boolean(v) => Ok(Expr::Literal(Scalar::Boolean(*v))),
        AstKind::RowRef => Ok(Expr::RowNumber),
        AstKind::NullRef => Ok(Expr::Null(ValueSort::Long)),

        AstKind::Ident(_) | AstKind::Keyword(_) => match names.get(&ast.at) {
            Some(ParserValue::Column { index, sort }) => match sort {
                /* the string and bit-string columns need their own kernels */
                ColumnSort::Numeric | ColumnSort::Boolean => Ok(Expr::Column(*index as usize)),
                ColumnSort::String => Ok(Expr::Column(*index as usize)),
                ColumnSort::Bits => Ok(Expr::Column(*index as usize)),
            },
            Some(ParserValue::Long(v)) => Ok(Expr::Literal(Scalar::Long(*v))),
            Some(ParserValue::Double(v)) => Ok(Expr::Literal(Scalar::Double(*v))),
            Some(ParserValue::Boolean(v)) => Ok(Expr::Literal(Scalar::Boolean(*v))),
            /* a keyword's string value is NUL-terminated, and the text is
            what precedes the terminator */
            Some(ParserValue::Str(v)) => Ok(Expr::Literal(Scalar::Str(
                v.iter()
                    .take_while(|&&c| c != 0)
                    .map(|&c| c as u8)
                    .collect(),
            ))),
            None => Err(Unsupported("unresolved name")),
        },

        AstKind::Unary { op, arg } => {
            let inner = Box::new(lower(arg, names, cols)?);
            Ok(Expr::Unary {
                op: match op {
                    UnOp::Neg => Unary::Neg,
                    UnOp::Not => Unary::Not,
                    UnOp::IntCast => Unary::ToLong,
                    UnOp::FltCast => Unary::ToDouble,
                    /* `+x` is the identity once the sort has been checked, and
                    the parser has already done that */
                    UnOp::Plus => return Ok(*inner),
                },
                arg: inner,
            })
        }

        AstKind::Binary { op, lhs, rhs } => {
            let l = Box::new(lower(lhs, names, cols)?);
            let r = Box::new(lower(rhs, names, cols)?);
            if let Some(a) = arith_of(*op) {
                return Ok(Expr::Arith {
                    op: a,
                    lhs: l,
                    rhs: r,
                });
            }
            if let Some(c) = compare_of(*op) {
                return Ok(Expr::Compare {
                    op: c,
                    lhs: l,
                    rhs: r,
                });
            }
            match op {
                /* `==` and `!=` are numeric comparison or boolean equality
                depending on the operands, which lowering cannot see */
                BinOp::Eq | BinOp::Ne => Ok(Expr::Equality {
                    negated: *op == BinOp::Ne,
                    lhs: l,
                    rhs: r,
                }),
                BinOp::And => Ok(Expr::Logic {
                    op: Logic::And,
                    lhs: l,
                    rhs: r,
                }),
                BinOp::Or => Ok(Expr::Logic {
                    op: Logic::Or,
                    lhs: l,
                    rhs: r,
                }),
                /* On integers these are bitwise; on bit strings they are set
                operations, which need their own kernel. The arena lowering has
                already rejected every other combination, so reaching here with
                a bit string means the operands were bit strings -- and those
                are refused earlier, at the literal or the column. */
                BinOp::BitAnd => Ok(Expr::Arith {
                    op: Arith::BitAnd,
                    lhs: l,
                    rhs: r,
                }),
                BinOp::BitOr => Ok(Expr::Arith {
                    op: Arith::BitOr,
                    lhs: l,
                    rhs: r,
                }),
                BinOp::BitXor => Ok(Expr::Arith {
                    op: Arith::BitXor,
                    lhs: l,
                    rhs: r,
                }),
                _ => Err(Unsupported("operator")),
            }
        }

        AstKind::Ternary { cond, then, els } => Ok(Expr::IfThenElse {
            cond: Box::new(lower(cond, names, cols)?),
            then: Box::new(lower(then, names, cols)?),
            els: Box::new(lower(els, names, cols)?),
        }),

        /* `x = lo : hi` desugars to `lo <= x && x <= hi`, as `eval.y` did */
        AstKind::Range { val, lo, hi } => {
            let v = lower(val, names, cols)?;
            Ok(Expr::Logic {
                op: Logic::And,
                lhs: Box::new(Expr::Compare {
                    op: Compare::Lte,
                    lhs: Box::new(lower(lo, names, cols)?),
                    rhs: Box::new(v.clone()),
                }),
                rhs: Box::new(Expr::Compare {
                    op: Compare::Lte,
                    lhs: Box::new(v),
                    rhs: Box::new(lower(hi, names, cols)?),
                }),
            })
        }

        AstKind::Str(s) => Ok(Expr::Literal(Scalar::Str(s.clone()))),
        AstKind::BitStr(s) => Ok(Expr::Literal(Scalar::Bits(s.clone()))),
        /* `#SNULL` is the undefined string, which the arena builds as a null
        of string sort */
        AstKind::SNullRef => Ok(Expr::Null(ValueSort::String)),
        /* The parser has already checked the offset is a constant integer,
        so it is folded here rather than evaluated per row. */
        AstKind::Offset { name, off } => {
            let Some(ParserValue::Column { index, sort }) = names.get(&ast.at) else {
                return Err(Unsupported("row offset of a non-column"));
            };
            /* `shifted` reads a text column the same way as a numeric one,
            one entry per row, so a string offset needs nothing extra */
            let _ = (name, sort);
            match const_long(off) {
                Some(offset) => Ok(Expr::Offset {
                    col: *index as usize,
                    offset,
                }),
                None => Err(Unsupported("row offset that is not a plain constant")),
            }
        }
        /* Subscripting needs the operand's `naxes`, not just its element
        count: with a 2x3 column, `M[1,1]` picks an element but `M[1]` selects
        a whole slice, so even the single-subscript form cannot be lowered
        without the shape. Threading `naxes` through the value model is a
        design addition rather than another kernel. */
        /* A subscript needs the operand's shape, not just its element count:
        on a 2x3 column `M[1,1]` picks an element but `M[1]` selects a whole
        slice. The shape is known here — for a column from the table, for a
        vector literal from its length — so only the fully-indexed form is
        lowered and a slice still falls back. */
        AstKind::Deref { base, idx } => {
            let naxes = match &base.kind {
                AstKind::Ident(_) | AstKind::Keyword(_) => match names.get(&base.at) {
                    Some(ParserValue::Column { index, .. }) => cols
                        .shapes
                        .get(*index as usize)
                        .map(|(_, naxes)| naxes.iter().map(|&n| n as usize).collect::<Vec<_>>())
                        .ok_or(Unsupported("subscript of unknown column"))?,
                    _ => return Err(Unsupported("subscript of a scalar")),
                },
                AstKind::Vector(items) => vec![items.len()],
                _ => return Err(Unsupported("subscript of a computed value")),
            };
            /* One index against a shape with more axes selects a slice
            rather than an element -- the engine only does that for a constant
            index, which is what `Do_Deref`'s `allConst && nDims == 1` path
            handles, so anything else still falls back. */
            if naxes.len() != idx.len() {
                if idx.len() != 1 || naxes.is_empty() {
                    return Err(Unsupported("subscript slice"));
                }
                let Some(index) = const_long(&idx[0]) else {
                    return Err(Unsupported("subscript slice with a computed index"));
                };
                let inner = &naxes[..naxes.len() - 1];
                return Ok(Expr::Slice {
                    base: Box::new(lower(base, names, cols)?),
                    index,
                    run: inner.iter().product::<usize>().max(1),
                    axis_len: *naxes.last().unwrap(),
                    naxes: inner.to_vec(),
                });
            }
            let base = lower(base, names, cols)?;
            let idx: Result<Vec<Expr>, Unsupported> =
                idx.iter().map(|e| lower(e, names, cols)).collect();
            Ok(Expr::Deref {
                base: Box::new(base),
                idx: idx?,
                naxes,
            })
        }
        AstKind::Vector(items) => {
            let lowered: Result<Vec<Expr>, Unsupported> =
                items.iter().map(|e| lower(e, names, cols)).collect();
            Ok(Expr::Vector(lowered?))
        }
        AstKind::Call { kind, name, args } => {
            /* ISNULL lexes as a BFUNCTION; the other names in that class are
            the region tests, and GTI/STRSTR each need their own kernels */
            /* GTIFILTER and GTIFIND read their intervals while the arena was
            built; the lowering picks them up rather than reading the file
            again. GTIOVERLAP and the region filters are not ported yet. */
            if matches!(kind, CallKind::GtiFilter | CallKind::GtiFind) {
                let Some(data) = cols.gti.get(&ast.at) else {
                    return Err(Unsupported("GTI call with no recorded load"));
                };
                /* the time is the second argument when there is one, and the
                column `New_GTI` resolved otherwise */
                let time = match args.get(1) {
                    Some(a) => lower(a, names, cols)?,
                    None => match data.time_col {
                        Some(c) => Expr::Column(c as usize),
                        None => return Err(Unsupported("GTI call with no time column")),
                    },
                };
                return Ok(Expr::Gti {
                    find: *kind == CallKind::GtiFind,
                    time: Box::new(time),
                    intervals: alloc::rc::Rc::new(GtiIntervals {
                        start: data.start.clone(),
                        stop: data.stop.clone(),
                        ordered: data.ordered,
                    }),
                });
            }
            let admissible = *kind == CallKind::Function
                || (*kind == CallKind::BFunction
                    && matches!(
                        name.as_slice(),
                        b"ISNULL" | b"NEAR" | b"CIRCLE" | b"BOX" | b"ELLIPSE"
                    ))
                /* STRSTR is the only IFunction */
                || *kind == CallKind::IFunction;
            if !admissible {
                return Err(Unsupported("function call"));
            }
            /* The shape functions are answers about the expression, not about
            the data, so they fold to a constant here exactly as the arena
            folds them at parse time. */
            if let Some(folded) = fold_shape_fn(name, args, names, cols)? {
                return Ok(folded);
            }
            /* ACCUM and SEQDIFF carry a running value between rows and
            batches, so each gets a slot rather than being a plain kernel */
            if matches!(name.as_slice(), b"ACCUM" | b"SEQDIFF") {
                if args.len() != 1 {
                    return Err(Unsupported("ACCUM with the wrong arity"));
                }
                let arg = lower(&args[0], names, cols)?;
                /* the boolean and bit-string forms take a different path in
                the engine, and one of them is the known type confusion, so
                they are left with the arena */
                let sort = arg.sort(&|i| cols.sort(i));
                if !matches!(sort, ValueSort::Long | ValueSort::Double) {
                    return Err(Unsupported("ACCUM over a non-numeric argument"));
                }
                let id = cols.accums.get();
                cols.accums.set(id + 1);
                return Ok(Expr::Accum {
                    id,
                    diff: name.as_slice() == b"SEQDIFF",
                    arg: Box::new(arg),
                });
            }
            /* MIN and MAX map elementwise with two arguments and reduce with
            one, so the arity picks the kernel */
            let (func, arity) = match (name.as_slice(), args.len()) {
                (b"MIN", 1) => (Func::Min1, 1),
                (b"MAX", 1) => (Func::Max1, 1),
                /* the generators take an optional argument that only sets the
                shape, so the arity follows what was written */
                (b"RANDOM", n @ (0 | 1)) => (Func::Random, n),
                (b"RANDOMN", n @ (0 | 1)) => (Func::RandomN, n),
                _ => match func_of(name) {
                    Some(pair) => pair,
                    None => return Err(Unsupported("function call")),
                },
            };
            if args.len() != arity {
                return Err(Unsupported("function call"));
            }
            let lowered: Result<Vec<Expr>, Unsupported> =
                args.iter().map(|a| lower(a, names, cols)).collect();
            Ok(Expr::Call {
                func,
                args: lowered?,
            })
        }
    }
}

/// The elementwise functions and their arity.
///
/// Names arrive upper-cased from the lexer. `MIN` and `MAX` appear here in
/// their two-argument form only; with one argument they reduce across a row,
/// which needs a different kernel.
fn func_of(name: &[u8]) -> Option<(Func, usize)> {
    Some(match name {
        b"SIN" => (Func::Sin, 1),
        b"COS" => (Func::Cos, 1),
        b"TAN" => (Func::Tan, 1),
        b"ARCSIN" | b"ASIN" => (Func::Asin, 1),
        b"ARCCOS" | b"ACOS" => (Func::Acos, 1),
        b"ARCTAN" | b"ATAN" => (Func::Atan, 1),
        b"SINH" => (Func::Sinh, 1),
        b"COSH" => (Func::Cosh, 1),
        b"TANH" => (Func::Tanh, 1),
        b"EXP" => (Func::Exp, 1),
        b"LOG" => (Func::Ln, 1),
        b"LOG10" => (Func::Log10, 1),
        b"SQRT" => (Func::Sqrt, 1),
        b"CEIL" => (Func::Ceil, 1),
        b"FLOOR" => (Func::Floor, 1),
        b"ROUND" => (Func::Round, 1),
        b"ABS" => (Func::Abs, 1),
        b"ARCTAN2" => (Func::Atan2, 2),
        b"MIN" => (Func::Min, 2),
        b"MAX" => (Func::Max, 2),
        b"DEFNULL" => (Func::DefNull, 2),
        b"STRSTR" => (Func::StrStr, 2),
        b"ANGSEP" => (Func::AngSep, 4),
        b"RANDOMP" => (Func::RandomP, 1),
        b"ELEMENTNUM" => (Func::ElementNum, 1),
        b"AXISELEM" => (Func::AxisElem, 2),
        b"ARRAY" => (Func::Array, 2),
        b"NEAR" => (Func::Near, 3),
        b"CIRCLE" => (Func::Circle, 5),
        b"BOX" => (Func::Box, 7),
        b"ELLIPSE" => (Func::Ellipse, 7),
        b"STRMID" => (Func::StrMid, 3),
        b"SETNULL" => (Func::SetNull, 2),
        b"ISNULL" => (Func::IsNull, 1),
        b"SUM" => (Func::Sum, 1),
        b"AVERAGE" => (Func::Average, 1),
        b"STDDEV" => (Func::Stddev, 1),
        b"MEDIAN" => (Func::Median, 1),
        b"NVALID" => (Func::NValid, 1),
        _ => return None,
    })
}

fn arith_of(op: BinOp) -> Option<Arith> {
    Some(match op {
        BinOp::Add => Arith::Add,
        BinOp::Sub => Arith::Sub,
        BinOp::Mul => Arith::Mul,
        BinOp::Div => Arith::Div,
        BinOp::Mod => Arith::Mod,
        BinOp::Pow => Arith::Pow,
        _ => return None,
    })
}

fn compare_of(op: BinOp) -> Option<Compare> {
    Some(match op {
        BinOp::Gt => Compare::Gt,
        BinOp::Lt => Compare::Lt,
        BinOp::Gte => Compare::Gte,
        BinOp::Lte => Compare::Lte,
        BinOp::Approx => Compare::Approx,
        _ => return None,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::{grammar, lexer};

    fn ast(src: &str) -> Ast {
        let toks = lexer::tokenize(format!("{src}\n").as_bytes()).expect("lex");
        grammar::parse(&toks).expect("parse").expect("non-empty")
    }

    fn try_lower(src: &str) -> Res {
        lower(
            &ast(src),
            &Resolutions::new(),
            &Columns {
                gti: Default::default(),
                accums: core::cell::Cell::new(0),
                shapes: Vec::new(),
                sorts: Vec::new(),
            },
        )
    }

    #[test]
    fn arithmetic_and_comparison_lower() {
        assert!(try_lower("1 + 2 * 3").is_ok());
        assert!(try_lower("1 > 2").is_ok());
        assert!(try_lower("1 = 2 : 3").is_ok());
        assert!(try_lower("T && F").is_ok());
        assert!(try_lower("T ? 1 : 2").is_ok());
        assert!(try_lower("-(1)").is_ok());
        assert!(try_lower("(int)2.5").is_ok());
        assert!(try_lower("#ROW").is_ok());
    }

    #[test]
    fn unary_plus_disappears() {
        /* `+x` has no node of its own once the sort check has happened */
        assert_eq!(try_lower("+1").unwrap(), Expr::Literal(Scalar::Long(1)));
    }

    #[test]
    fn a_range_becomes_two_comparisons() {
        let e = try_lower("2 = 1 : 3").unwrap();
        match e {
            Expr::Logic { op: Logic::And, .. } => {}
            other => panic!("expected a conjunction, got {other:?}"),
        }
    }

    #[test]
    fn the_elementwise_functions_lower() {
        for src in [
            "SIN(1)",
            "sqrt(4)",
            "ABS(-1)",
            "ROUND(2.5)",
            "ARCTAN2(1,1)",
            "MIN(1,2)",
            "MAX(1,2)",
            "DEFNULL(1,2)",
            "SETNULL(1,2)",
            "ISNULL(1)",
        ] {
            assert!(try_lower(src).is_ok(), "src: {src}");
        }
    }

    #[test]
    fn the_reductions_lower() {
        /* MIN and MAX pick a kernel by arity: one argument reduces a row,
        two map elementwise */
        for src in [
            "SUM(1)",
            "AVERAGE(1)",
            "STDDEV(1)",
            "MEDIAN(1)",
            "NVALID(1)",
            "MIN(1)",
            "MAX(1)",
        ] {
            assert!(try_lower(src).is_ok(), "src: {src}");
        }
    }

    #[test]
    fn vector_literals_lower() {
        for src in ["{1,2}", "{1,2,3}", "{1.5,2}", "{T,F}"] {
            assert!(try_lower(src).is_ok(), "src: {src}");
        }
    }

    #[test]
    fn the_bitwise_operators_lower() {
        for src in ["1 & 2", "1 | 2", "1 ^^ 2"] {
            assert!(try_lower(src).is_ok(), "src: {src}");
        }
    }

    #[test]
    fn bit_strings_lower() {
        /* the operators over them are the bit kernels, chosen at evaluation
        from the operand sorts */
        for src in ["b101 == b110", "b101 & b110", "!b101", "b101 > b110"] {
            assert!(try_lower(src).is_ok(), "src: {src}");
        }
    }

    #[test]
    fn strings_lower() {
        for src in [
            "'a' == 'b'",
            "'a' + 'b'",
            "STRSTR('abc','b')",
            "STRMID('abc',1,2)",
        ] {
            assert!(try_lower(src).is_ok(), "src: {src}");
        }
    }

    #[test]
    fn the_generators_lower_with_or_without_an_argument() {
        for src in [
            "RANDOM()",
            "RANDOM(1)",
            "RANDOMN()",
            "RANDOMN(1)",
            "RANDOMP(1)",
        ] {
            assert!(try_lower(src).is_ok(), "src: {src}");
        }
    }

    #[test]
    fn the_undefined_string_lowers_as_a_string() {
        /* not as a long, or a conditional with a string branch cannot read it */
        assert_eq!(try_lower("#SNULL"), Ok(Expr::Null(ValueSort::String)));
    }

    #[test]
    fn the_shape_functions_fold_to_a_constant() {
        /* with no columns declared every argument is a scalar, which is the
        one shape these tests can assert without a table */
        assert_eq!(try_lower("NELEM(1)"), Ok(Expr::Literal(Scalar::Long(1))));
        assert_eq!(try_lower("NAXIS(1)"), Ok(Expr::Literal(Scalar::Long(1))));
        assert_eq!(try_lower("NAXES(1,1)"), Ok(Expr::Literal(Scalar::Long(1))));
        /* a string literal's element count is its length */
        assert_eq!(
            try_lower("NELEM('abc')"),
            Ok(Expr::Literal(Scalar::Long(3)))
        );
        assert_eq!(
            try_lower("NELEM({1,2,3})"),
            Ok(Expr::Literal(Scalar::Long(3)))
        );
        /* text `+` concatenates, so the widths add rather than combine */
        assert_eq!(
            try_lower("NELEM('ab' + 'cde')"),
            Ok(Expr::Literal(Scalar::Long(5)))
        );
        /* an axis index is clamped at the bottom, as NAXES(x,0) shows */
        assert_eq!(
            try_lower("NAXES({1,2},0)"),
            Ok(Expr::Literal(Scalar::Long(2)))
        );
    }

    #[test]
    fn unported_constructs_name_themselves() {
        for (src, want) in [
            /* GTIOVERLAP and the region filters are not ported yet */
            ("GTIOVERLAP('f',1,2)", "function call"),
            ("REGFILTER('x.reg')", "function call"),
        ] {
            assert_eq!(try_lower(src), Err(Unsupported(want)), "src: {src}");
        }
    }

    #[test]
    fn an_unsupported_leaf_stops_the_whole_expression() {
        /* the fallback is per expression, not per node */
        assert_eq!(
            try_lower("1 + REGFILTER('x.reg')"),
            Err(Unsupported("function call"))
        );
    }
}
