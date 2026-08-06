//! The expression evaluator.
//!
//! `PARSER_MIGRATION.md` §9 and §10 set out the migration this completes. The
//! parser hands back an `Ast`; this module lowers it to an [`expr::Expr`] and
//! evaluates that -- a value model where a batch of rows is an owned array
//! with its own null mask, and each operation is a kernel over two of them.
//!
//! There is no longer a `Node` arena, and no fallback: every expression the
//! library accepts is lowered here, and [`lower::lower`] is what decides
//! whether it is accepted at all.

pub(crate) mod bits;
pub(crate) mod bridge;
pub(crate) mod expr;
pub(crate) mod gti;
pub(crate) mod kernel;
pub(crate) mod lower;
pub(crate) mod regions;
pub(crate) mod strings;
pub(crate) mod value;
