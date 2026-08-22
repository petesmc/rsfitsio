//! `errno` values used by the vendored code.
// Vendored code, documented in upstream's style.
#![allow(missing_docs)]

use crate::c_types::c_int;

pub const EILSEQ: c_int = 84; /* Illegal byte sequence */
