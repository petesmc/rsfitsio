//! `Ast` to [`Expr`], for the subset the new evaluator covers.
//!
//! The port is incremental: anything not yet handled returns
//! [`Unsupported`], and the caller keeps the `Node` arena for that expression.
//! That is what lets the new evaluator go in behind a flag and grow one kernel
//! family at a time, with `tests/test_eval_corpus.rs` proving each step.
//!
//! Name resolution has already happened — `parser::resolve` ran before the
//! parse — so a name here is a column index or a folded constant.

use super::expr::{Expr, Logic, Unary};
use super::kernel::{Arith, Compare};
use super::value::Scalar;
use crate::eval_defs::{ColumnSort, ParserValue, ValueSort};
use crate::parser::ast::{Ast, AstKind, BinOp, UnOp};
use crate::parser::resolve::Resolutions;

/// Why an expression could not be lowered to an [`Expr`].
///
/// Carries the construct's name so that a corpus divergence points straight at
/// what is still missing rather than at "something".
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct Unsupported(pub(crate) &'static str);

type Res = Result<Expr, Unsupported>;

/// Lower a whole expression, or report the first construct not yet ported.
pub(crate) fn lower(ast: &Ast, names: &Resolutions) -> Res {
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
            let inner = Box::new(lower(arg, names)?);
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
            let l = Box::new(lower(lhs, names)?);
            let r = Box::new(lower(rhs, names)?);
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
                /* `&`, `|` and `^^` are bitwise on integers and set operations
                on bit strings; neither kernel exists yet */
                BinOp::BitAnd | BinOp::BitOr | BinOp::BitXor => {
                    Err(Unsupported("bitwise operator"))
                }
                _ => Err(Unsupported("operator")),
            }
        }

        AstKind::Ternary { cond, then, els } => Ok(Expr::IfThenElse {
            cond: Box::new(lower(cond, names)?),
            then: Box::new(lower(then, names)?),
            els: Box::new(lower(els, names)?),
        }),

        /* `x = lo : hi` desugars to `lo <= x && x <= hi`, as `eval.y` did */
        AstKind::Range { val, lo, hi } => {
            let v = lower(val, names)?;
            Ok(Expr::Logic {
                op: Logic::And,
                lhs: Box::new(Expr::Compare {
                    op: Compare::Lte,
                    lhs: Box::new(lower(lo, names)?),
                    rhs: Box::new(v.clone()),
                }),
                rhs: Box::new(Expr::Compare {
                    op: Compare::Lte,
                    lhs: Box::new(v),
                    rhs: Box::new(lower(hi, names)?),
                }),
            })
        }

        AstKind::Str(_) => Err(Unsupported("string literal")),
        AstKind::BitStr(_) => Err(Unsupported("bit-string literal")),
        AstKind::SNullRef => Err(Unsupported("#SNULL")),
        AstKind::Offset { .. } => Err(Unsupported("row offset")),
        AstKind::Deref { .. } => Err(Unsupported("subscript")),
        AstKind::Vector(_) => Err(Unsupported("vector literal")),
        AstKind::Call { .. } => Err(Unsupported("function call")),
    }
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
        lower(&ast(src), &Resolutions::new())
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
    fn unported_constructs_name_themselves() {
        for (src, want) in [
            ("'a'", "string literal"),
            ("b101", "bit-string literal"),
            ("1 & 2", "bitwise operator"),
            ("{1,2}", "vector literal"),
            ("SIN(1)", "function call"),
            ("#SNULL", "#SNULL"),
        ] {
            assert_eq!(try_lower(src), Err(Unsupported(want)), "src: {src}");
        }
    }

    #[test]
    fn an_unsupported_leaf_stops_the_whole_expression() {
        /* the fallback is per expression, not per node */
        assert_eq!(try_lower("1 + SIN(1)"), Err(Unsupported("function call")));
    }
}
