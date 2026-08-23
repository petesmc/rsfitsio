//! Vendored pieces of a C standard library, and of Rust's own `std::io`.
//!
//! The port needs `printf`/`scanf` formatting and buffered stream I/O without
//! calling libc for them. Rather than write those from scratch, this module
//! carries the relevant parts of the Redox project's relibc
//! (<https://gitlab.redox-os.org/redox-os/relibc>) and of the Rust project's
//! `std::io`.
//!
//! **Vendored code.** It is kept close to upstream so it can be re-synced, so
//! it is documented in upstream's style rather than converted to this crate's;
//! see the plan's decision on vendored code. Copyright notices are preserved in
//! the files that carry them.
#![warn(missing_docs)]

pub mod header;
// Vendored from the Rust project's std::io (see the copyright headers in
// io/mod.rs and io/buffered.rs). Kept byte-identical to upstream so it can be
// re-synced, so its lints are allowed here rather than fixed in place.
#[allow(clippy::manual_clear, clippy::redundant_guards)]
pub mod io;
pub mod platform;
