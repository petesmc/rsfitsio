//! The replacement expression evaluator.
//!
//! `PARSER_MIGRATION.md` §9 sets out the plan: the parser now hands back an
//! `Ast`, which gives a seam between the front end and the `Node` arena that
//! `eval_y` still evaluates. This module is the second half of that migration —
//! a value model where a batch of rows is an owned array with its own null
//! mask, and each operation is a kernel over two of them.
//!
//! **This is the evaluator.** Every expression the parser accepts is lowered
//! here and evaluated here, unless it hits one of the cases the lowering
//! declines -- a bit-valued result, a random generator, a row offset reaching
//! outside the loaded chunk -- which still go to the `Node` arena in `eval_y`.
//! `PARSER_MIGRATION.md` section 10.3 lists them and says why each is where it
//! is.

pub(crate) mod bits;
pub(crate) mod bridge;
pub(crate) mod expr;
pub(crate) mod gti;
pub(crate) mod kernel;
pub(crate) mod lower;
pub(crate) mod regions;
pub(crate) mod strings;
pub(crate) mod value;
