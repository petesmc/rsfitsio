//! The replacement expression evaluator.
//!
//! `PARSER_MIGRATION.md` §9 sets out the plan: the parser now hands back an
//! `Ast`, which gives a seam between the front end and the `Node` arena that
//! `eval_y` still evaluates. This module is the second half of that migration —
//! a value model where a batch of rows is an owned array with its own null
//! mask, and each operation is a kernel over two of them.
//!
//! **Status: in progress.** The value layer, the kernels, the `Ast -> Expr`
//! lowering and the bridge back into the arena's result node are here. Behind
//! `--features new-eval` most of the corpus evaluates through this; anything
//! the lowering cannot yet express returns `Unsupported` and falls back to the
//! arena, per expression. `PARSER_MIGRATION.md` section 10.3 tracks what is
//! left and why.
#![allow(dead_code)]

pub(crate) mod bits;
#[cfg(feature = "new-eval")]
pub(crate) mod bridge;
pub(crate) mod expr;
pub(crate) mod kernel;
pub(crate) mod lower;
pub(crate) mod strings;
pub(crate) mod value;
