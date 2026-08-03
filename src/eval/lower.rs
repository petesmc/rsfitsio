//! `Ast` to [`Expr`], for the subset the new evaluator covers.
//!
//! The port is incremental: anything not yet handled returns
//! [`Unsupported`], and the caller keeps the `Node` arena for that expression.
//! That is what lets the new evaluator go in behind a flag and grow one kernel
//! family at a time, with `tests/test_eval_corpus.rs` proving each step.
//!
//! Name resolution has already happened — `parser::resolve` ran before the
//! parse — so a name here is a column index or a folded constant.

use super::expr::{Expr, Func, Logic, Unary};
use super::kernel::{Arith, Compare};
use super::value::Scalar;
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

/// Lower a whole expression, or report the first construct not yet ported.
pub(crate) fn lower(ast: &Ast, names: &Resolutions, shapes: &[Vec<usize>]) -> Res {
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
                ColumnSort::String => Err(Unsupported("string column")),
                ColumnSort::Bits => Err(Unsupported("bit-string column")),
            },
            Some(ParserValue::Long(v)) => Ok(Expr::Literal(Scalar::Long(*v))),
            Some(ParserValue::Double(v)) => Ok(Expr::Literal(Scalar::Double(*v))),
            Some(ParserValue::Boolean(v)) => Ok(Expr::Literal(Scalar::Boolean(*v))),
            Some(ParserValue::Str(_)) => Err(Unsupported("string keyword")),
            None => Err(Unsupported("unresolved name")),
        },

        AstKind::Unary { op, arg } => {
            let inner = Box::new(lower(arg, names, shapes)?);
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
            let l = Box::new(lower(lhs, names, shapes)?);
            let r = Box::new(lower(rhs, names, shapes)?);
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
            cond: Box::new(lower(cond, names, shapes)?),
            then: Box::new(lower(then, names, shapes)?),
            els: Box::new(lower(els, names, shapes)?),
        }),

        /* `x = lo : hi` desugars to `lo <= x && x <= hi`, as `eval.y` did */
        AstKind::Range { val, lo, hi } => {
            let v = lower(val, names, shapes)?;
            Ok(Expr::Logic {
                op: Logic::And,
                lhs: Box::new(Expr::Compare {
                    op: Compare::Lte,
                    lhs: Box::new(lower(lo, names, shapes)?),
                    rhs: Box::new(v.clone()),
                }),
                rhs: Box::new(Expr::Compare {
                    op: Compare::Lte,
                    lhs: Box::new(v),
                    rhs: Box::new(lower(hi, names, shapes)?),
                }),
            })
        }

        AstKind::Str(_) => Err(Unsupported("string literal")),
        AstKind::BitStr(_) => Err(Unsupported("bit-string literal")),
        AstKind::SNullRef => Err(Unsupported("#SNULL")),
        AstKind::Offset { .. } => Err(Unsupported("row offset")),
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
                    Some(ParserValue::Column { index, .. }) => shapes
                        .get(*index as usize)
                        .cloned()
                        .ok_or(Unsupported("subscript of unknown column"))?,
                    _ => return Err(Unsupported("subscript of a scalar")),
                },
                AstKind::Vector(items) => vec![items.len()],
                _ => return Err(Unsupported("subscript of a computed value")),
            };
            if naxes.len() != idx.len() {
                /* a partial index selects a slice, which needs a shape-aware
                gather rather than one element per row */
                return Err(Unsupported("subscript slice"));
            }
            let base = lower(base, names, shapes)?;
            let idx: Result<Vec<Expr>, Unsupported> =
                idx.iter().map(|e| lower(e, names, shapes)).collect();
            Ok(Expr::Deref {
                base: Box::new(base),
                idx: idx?,
                naxes,
            })
        }
        AstKind::Vector(items) => {
            let lowered: Result<Vec<Expr>, Unsupported> =
                items.iter().map(|e| lower(e, names, shapes)).collect();
            Ok(Expr::Vector(lowered?))
        }
        AstKind::Call { kind, name, args } => {
            /* ISNULL lexes as a BFUNCTION; the other names in that class are
            the region tests, and GTI/STRSTR each need their own kernels */
            let admissible = *kind == CallKind::Function
                || (*kind == CallKind::BFunction && name.as_slice() == b"ISNULL");
            if !admissible {
                return Err(Unsupported("function call"));
            }
            /* MIN and MAX map elementwise with two arguments and reduce with
            one, so the arity picks the kernel */
            let (func, arity) = match (name.as_slice(), args.len()) {
                (b"MIN", 1) => (Func::Min1, 1),
                (b"MAX", 1) => (Func::Max1, 1),
                _ => match func_of(name) {
                    Some(pair) => pair,
                    None => return Err(Unsupported("function call")),
                },
            };
            if args.len() != arity {
                return Err(Unsupported("function call"));
            }
            let lowered: Result<Vec<Expr>, Unsupported> =
                args.iter().map(|a| lower(a, names, shapes)).collect();
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
        lower(&ast(src), &Resolutions::new(), &[])
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
    fn unported_constructs_name_themselves() {
        for (src, want) in [
            ("'a'", "string literal"),
            ("b101", "bit-string literal"),
            ("NELEM(1)", "function call"),
            ("RANDOM()", "function call"),
            ("GTIFILTER()", "function call"),
            ("#SNULL", "#SNULL"),
        ] {
            assert_eq!(try_lower(src), Err(Unsupported(want)), "src: {src}");
        }
    }

    #[test]
    fn an_unsupported_leaf_stops_the_whole_expression() {
        /* the fallback is per expression, not per node */
        assert_eq!(try_lower("1 + NELEM(1)"), Err(Unsupported("function call")));
    }
}
