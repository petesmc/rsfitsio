//! Small Rust-side helpers with no CFITSIO counterpart.
//!
//! These exist to express in safe Rust what the C did with raw allocation and
//! `FILE*` I/O: fallible boxing, an over-aligned byte buffer, the registry that
//! carries a `Vec`'s layout across a C-ABI boundary, a `Read`/`Write` adapter
//! over a C stream, and the test fixtures.
#![warn(missing_docs)]

pub mod aligned;
pub mod boxed;
pub mod cfile;
pub mod raw_owned;
#[cfg(test)]
pub mod testhelpers;
pub mod vec_raw_parts;
