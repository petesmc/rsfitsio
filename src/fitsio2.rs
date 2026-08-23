//! Internal constants and helpers shared across the library.
//!
//! Ported from CFITSIO's `fitsio2.h`, the private companion to `fitsio.h`. It
//! holds the buffer sizes, machine-type flags, IEEE special-value tests and
//! the floating-point bounds used by the datatype conversion routines.
#![warn(missing_docs)]

use core::ffi::c_char;
use core::mem;

use crate::c_types::{c_int, c_long, c_short, size_t};

use crate::wrappers::{strcmp_safe, strncmp_safe};

/// Flag used when writing images, meaning the value is too large for the
/// narrower parameter it would otherwise be passed in.
pub const USE_LARGE_VALUE: c_int = -99;

/// Size of the data buffer, in bytes.
pub const DBUFFSIZE: u64 = 28800;

/// Maximum number of FITS files that can be open at once.
///
/// CFITSIO allocates `NMAXFILES * 80` bytes up front, and each file that is
/// opened uses a further [`NIOBUF`](crate::fitsio::NIOBUF) × 2880 bytes.
pub const NMAXFILES: usize = 10000;

/// Minimum transfer size for direct reads and writes, bypassing the I/O
/// buffers. Must be at least 8640.
pub const MINDIRECT: i64 = 8640;

/*   it is useful to identify certain specific types of machines   */
/// Machine type: uses non-byteswapped IEEE formats.
pub const NATIVE: u64 = 0;
/// Machine type: any other type of machine.
pub const OTHERTYPE: u64 = 1;
/// Machine type: IBM PC. Used in `drvrfile.rs` to work around a bug on PCs.
pub const IBMPC: u64 = 5;
/// Machine type: Cray, which requires a special NaN test algorithm.
pub const CRAY: u64 = 6;

// VAXVMS (3) and ALPHAVMS (4) are deliberately absent, as are FLOATTYPE,
// GFLOAT and IEEEFLOAT. rustc has no VMS target and dropped VAX and Alpha
// support entirely, so CFITSIO_MACHINE could never take those values and
// every branch guarded on them was unreachable. The affected code in
// fitsio2.h / buffers.c / drvrfile.c / drvrmem.c / getcol*.c is:
//
//   - the G_FLOAT <-> IEEE conversions in ffgr4b/ffgr8b/ffpr4b/ffpr8b
//     (buffers.c), which were transpiled as todo!() and are now gone
//   - the `rfm=fix, mrs=2880, ctx=stm` fopen record-format arguments in
//     file_openfile/file_create (drvrfile.c) and mem_createfile (drvrmem.c)
//   - the MSB offset guards, which in C read
//       `if (BYTESWAPPED && CFITSIO_MACHINE != VAXVMS
//            && CFITSIO_MACHINE != ALPHAVMS)`
//     and are now just `if BYTESWAPPED`
//
// iraffits.rs likewise dropped its `#[cfg(target_os = "vms")]` arms, which
// searched for ']' and ':' path separators instead of '/'.
// Restore from the C if a VMS target ever appears.

/*  assume all other machine uses the same IEEE formats as used in FITS files */
/*  e.g., Macs fall into this category  */

/// The machine type this build targets. Every supported target uses the same
/// IEEE formats as FITS files do.
pub const CFITSIO_MACHINE: u64 = NATIVE;
/// Whether this machine's byte order differs from the big-endian order FITS
/// stores numbers in, so values must be swapped on the way in and out.
pub const BYTESWAPPED: bool = true;

/// Width of a C `long` in bits: 64 on LP64 targets such as Linux and macOS, 32
/// on Windows.
pub const LONGSIZE: u64 = c_long::BITS as u64;

/// Reaching the end of the file is not an error.
pub const IGNORE_EOF: c_int = 1;
/// Reaching the end of the file should be reported as an error.
pub const REPORT_EOF: c_int = 0;
/// Marks a data value as undefined.
pub const DATA_UNDEFINED: i64 = -1;
/// Marks the null value of a column as not having been defined.
pub const NULL_UNDEFINED: c_int = 1234554321;
/// Marks an ASCII table column as having no defined null value.
pub const ASCII_NULL_UNDEFINED: c_char = 1;

/// Compares two nul-terminated strings, testing the first character inline
/// before falling back to a full comparison.
///
/// # Parameters
///
/// * `a` — (I) first string
/// * `b` — (I) second string
///
/// # Returns
///
/// Negative, zero or positive as `a` sorts before, equal to, or after `b`.
pub fn FSTRCMP(a: &[c_char], b: &[c_char]) -> c_int {
    if a[0] < b[0] {
        -1
    } else if a[0] > b[0] {
        1
    } else {
        strcmp_safe(a, b)
    }
}

/// Compares the first `n` characters of two nul-terminated strings, testing the
/// first character inline before falling back to a full comparison.
///
/// # Parameters
///
/// * `a` — (I) first string
/// * `b` — (I) second string
/// * `n` — (I) maximum number of characters to compare
///
/// # Returns
///
/// Negative, zero or positive as `a` sorts before, equal to, or after `b`.
pub fn FSTRNCMP(a: &[c_char], b: &[c_char], n: size_t) -> c_int {
    if a[0] < b[0] {
        -1
    } else if a[0] > b[0] {
        1
    } else {
        strncmp_safe(a, b, n)
    }
}

/// Mask over bits 1–8 of a `float`: all set on a NaN, all clear on an
/// underflow or zero.
pub const FNANMASK: u64 = 0x7F80;

/// Mask over bits 1–11 of a `double`: all set on a NaN, all clear on an
/// underflow or zero.
pub const DNANMASK: u64 = 0x7FF0;

/* these functions work for both big and little endian machines */
/* that use the IEEE floating point format for internal numbers */

/* These functions tests whether the float value is a reserved IEEE     */
/* value such as a Not-a-Number (NaN), or underflow, overflow, or       */
/* infinity.   The functions returns 1 if the value is a NaN, overflow  */
/* or infinity; it returns 2 if the value is an denormalized underflow  */
/* value; otherwise it returns 0. fnan tests floats, dnan tests doubles */
/*
pub const fnan: u64 =(L) \
      ( (L & FNANMASK) == FNANMASK ?  1 : (L & FNANMASK) == 0 ? 2 : 0)

pub const dnan: u64 =(L) \
      ( (L & DNANMASK) == DNANMASK ?  1 : (L & DNANMASK) == 0 ? 2 : 0)
 */

/// Tests whether a `float` is a reserved IEEE value.
///
/// Works on both big- and little-endian machines that use the IEEE floating
/// point format internally.
///
/// # Parameters
///
/// * `L` — (I) the high-order 16 bits of the float
///
/// # Returns
///
/// `1` for a NaN, overflow or infinity, `2` for a denormalized underflow value,
/// `0` otherwise.
pub fn fnan(L: c_short) -> c_int {
    if (L & FNANMASK as c_short) == FNANMASK as c_short {
        1
    } else if (L & FNANMASK as c_short) == 0 {
        2
    } else {
        0
    }
}

/// Tests whether a `double` is a reserved IEEE value.
///
/// Works on both big- and little-endian machines that use the IEEE floating
/// point format internally.
///
/// # Parameters
///
/// * `L` — (I) the high-order 16 bits of the double
///
/// # Returns
///
/// `1` for a NaN, overflow or infinity, `2` for a denormalized underflow value,
/// `0` otherwise.
pub fn dnan(L: c_short) -> c_int {
    if (L & DNANMASK as c_short) == DNANMASK as c_short {
        1
    } else if (L & DNANMASK as c_short) == 0 {
        2
    } else {
        0
    }
}

/// Largest `double` that fits in a `signed char`.
pub const DSCHAR_MAX: f64 = 127.49;
/// Smallest `double` that fits in a `signed char`.
pub const DSCHAR_MIN: f64 = -128.49;
/// Largest `double` that fits in an `unsigned char`.
pub const DUCHAR_MAX: f64 = 255.49;
/// Smallest `double` that fits in an `unsigned char`.
pub const DUCHAR_MIN: f64 = -0.49;
/// Largest `double` that fits in an `unsigned short`.
pub const DUSHRT_MAX: f64 = 65535.49;
/// Smallest `double` that fits in an `unsigned short`.
pub const DUSHRT_MIN: f64 = -0.49;
/// Largest `double` that fits in a `short`.
pub const DSHRT_MAX: f64 = 32767.49;
/// Smallest `double` that fits in a `short`.
pub const DSHRT_MIN: f64 = -32768.49;

/* These bounds depend on the width of a C `long`, which is 32 bits on Windows
 * (LLP64) and 64 bits on most Unix targets (LP64).  CFITSIO selects them via
 * `#if LONGSIZE == 32` in fitsio2.h; we mirror that by branching on the size of
 * `c_long` at const-eval time so the overflow checks stay correct on Windows. */
/// Largest `double` that fits in a `long`. Depends on the width of a C `long`.
pub const DLONG_MAX: f64 = if mem::size_of::<c_long>() == 4 {
    2147483647.49 /* max double value that fits in a long */
} else {
    9.223_372_036_854_775E18 /* max double value  long */
};
/// Smallest `double` that fits in a `long`. Depends on the width of a C `long`.
pub const DLONG_MIN: f64 = if mem::size_of::<c_long>() == 4 {
    -2147483648.49 /* min double value that fits in a long */
} else {
    -9.223_372_036_854_775E18 /* min double value  long */
};
/// Largest `double` that fits in an `unsigned long`. Depends on the width of a
/// C `long`.
pub const DULONG_MAX: f64 = if mem::size_of::<c_long>() == 4 {
    4294967295.49 /* max double that fits in a unsigned long */
} else {
    1.844_674_407_370_955E19 /* max double value  ulong */
};

/// Smallest `double` that fits in an `unsigned long`.
pub const DULONG_MIN: f64 = -0.49;
/// Largest `double` that fits in an unsigned 64-bit integer.
pub const DULONGLONG_MAX: f64 = 18446744073709551615.;
/// Smallest `double` that fits in an unsigned 64-bit integer.
pub const DULONGLONG_MIN: f64 = -0.49;
/// Largest `double` that fits in a signed 64-bit integer.
pub const DLONGLONG_MAX: f64 = 9.223_372_036_854_776E18;
/// Smallest `double` that fits in a signed 64-bit integer.
pub const DLONGLONG_MIN: f64 = -9.223_372_036_854_776E18;
/// Largest `double` that fits in an unsigned 4-byte integer.
pub const DUINT_MAX: f64 = 4294967295.49;
/// Smallest `double` that fits in an unsigned 4-byte integer.
pub const DUINT_MIN: f64 = -0.49;
/// Largest `double` that fits in a signed 4-byte integer.
pub const DINT_MAX: f64 = 2147483647.49;
/// Smallest `double` that fits in a signed 4-byte integer.
pub const DINT_MIN: f64 = -2147483648.49;

/// Largest unsigned 64-bit integer.
pub const UINT64_MAX: u64 = 18446744073709551615;

/// Largest unsigned 32-bit integer.
pub const UINT32_MAX: u32 = 4294967295;
/// Largest signed 32-bit integer.
pub const INT32_MAX: i32 = 2147483647;
/// Smallest signed 32-bit integer.
pub const INT32_MIN: i32 = -INT32_MAX - 1;

/// Value written for a null pixel in a tile-compressed integer image.
pub const COMPRESS_NULL_VALUE: i32 = -2147483647;
/// Length of the random-number sequence used when quantizing real numbers. Do
/// not change: the dithering is reproducible only for this length.
pub const N_RANDOM: usize = 10000;
