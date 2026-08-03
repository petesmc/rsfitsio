//! The replacement expression evaluator.
//!
//! `PARSER_MIGRATION.md` §9 sets out the plan: the parser now hands back an
//! `Ast`, which gives a seam between the front end and the `Node` arena that
//! `eval_y` still evaluates. This module is the second half of that migration —
//! a value model where a batch of rows is an owned array with its own null
//! mask, and each operation is a kernel over two of them.
//!
//! **Status: foundation only.** The value layer and the arithmetic and
//! comparison kernels are here and tested; the `Ast -> Expr` lowering and the
//! remaining kernels are not, so nothing calls into this yet. It is kept in the
//! tree rather than on a branch because the pieces are independently
//! reviewable and the tests are meaningful on their own.
#![allow(dead_code)]

#[cfg(feature = "new-eval")]
pub(crate) mod bridge;
pub(crate) mod expr;
pub(crate) mod kernel;
pub(crate) mod lower;
pub(crate) mod value;
