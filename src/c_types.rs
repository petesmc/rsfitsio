//! The C scalar types the port is written against.
//!
//! Re-exported from `libc` in one place so the transpiled code can say
//! `c_int` / `c_long` / `LONGLONG` exactly as the original CFITSIO C does.
//! Their widths are target-dependent — `c_long` is 64-bit on Linux and macOS
//! but 32-bit on Windows — which is why the port uses them rather than `i32` /
//! `i64`.
#![warn(missing_docs)]

pub use libc::{
    FILE, c_char, c_double, c_float, c_int, c_long, c_longlong, c_schar, c_short, c_uchar, c_uint,
    c_ulong, c_ulonglong, c_ushort, c_void, intmax_t, off_t, ptrdiff_t, size_t, ssize_t, uintmax_t,
};
