use std::ffi::c_char;

use crate::c_types::{c_int, c_long, c_short, size_t};

use crate::wrappers::{strcmp_safe, strncmp_safe};

pub const USE_LARGE_VALUE: c_int = -99; /* flag used when writing images */

pub const DBUFFSIZE: u64 = 28800; /* size of data buffer in bytes */

pub const NMAXFILES: usize = 10000; /* maximum number of FITS files that can be opened */
/* CFITSIO will allocate (NMAXFILES * 80) bytes of memory */
/* plus each file that is opened will use NIOBUF * 2880 bytes of memeory */
/* where NIOBUF is defined in fitio.h and has a default value of 40 */

pub const MINDIRECT: i64 = 8640; /* minimum size for direct reads and writes */
/* MINDIRECT must have a value >= 8640 */

/*   it is useful to identify certain specific types of machines   */
pub const NATIVE: u64 = 0; /* machine that uses non-byteswapped IEEE formats */
pub const OTHERTYPE: u64 = 1; /* any other type of machine */
pub const VAXVMS: u64 = 3; /* uses an odd floating point format */
pub const ALPHAVMS: u64 = 4; /* uses an odd floating point format */
pub const IBMPC: u64 = 5; /* used in drvrfile.c to work around a bug on PCs */
pub const CRAY: u64 = 6; /* requires a special NaN test algorithm */

pub const FLOATTYPE: u64 = IEEEFLOAT;
pub const GFLOAT: u64 = 1; /* used for VMS */
pub const IEEEFLOAT: u64 = 2; /* used for VMS */

/*  assume all other machine uses the same IEEE formats as used in FITS files */
/*  e.g., Macs fall into this category  */

pub const CFITSIO_MACHINE: u64 = NATIVE;
pub const BYTESWAPPED: bool = true;

/*  assume longs are 8 bytes long, unless previously set otherwise */
pub const LONGSIZE: u64 = c_long::BITS as u64;

/*       end of block that determine long size and byte swapping        */
/* ==================================================================== */

pub const IGNORE_EOF: c_int = 1;
pub const REPORT_EOF: c_int = 0;
pub const DATA_UNDEFINED: i64 = -1;
pub const NULL_UNDEFINED: c_int = 1234554321;
pub const ASCII_NULL_UNDEFINED: c_char = 1; /* indicate no defined null value */

pub fn FSTRCMP(a: &[c_char], b: &[c_char]) -> c_int {
    if a[0] < b[0] {
        -1
    } else if a[0] > b[0] {
        1
    } else {
        strcmp_safe(a, b)
    }
}

pub fn FSTRNCMP(a: &[c_char], b: &[c_char], n: size_t) -> c_int {
    if a[0] < b[0] {
        -1
    } else if a[0] > b[0] {
        1
    } else {
        strncmp_safe(a, b, n)
    }
}

pub const FNANMASK: u64 = 0x7F80; /* mask bits 1 - 8; all set on NaNs */
/* all 0 on underflow  or 0. */

pub const DNANMASK: u64 = 0x7FF0; /* mask bits 1 - 11; all set on NaNs */
/* all 0 on underflow  or 0. */

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

pub fn fnan(L: c_short) -> c_int {
    if (L & FNANMASK as c_short) == FNANMASK as c_short {
        1
    } else if (L & FNANMASK as c_short) == 0 {
        2
    } else {
        0
    }
}

pub fn dnan(L: c_short) -> c_int {
    if (L & DNANMASK as c_short) == DNANMASK as c_short {
        1
    } else if (L & DNANMASK as c_short) == 0 {
        2
    } else {
        0
    }
}

pub const DSCHAR_MAX: f64 = 127.49; /* max double value that fits in an signed char */
pub const DSCHAR_MIN: f64 = -128.49; /* min double value that fits in an signed char */
pub const DUCHAR_MAX: f64 = 255.49; /* max double value that fits in an unsigned char */
pub const DUCHAR_MIN: f64 = -0.49; /* min double value that fits in an unsigned char */
pub const DUSHRT_MAX: f64 = 65535.49; /* max double value that fits in a unsigned short*/
pub const DUSHRT_MIN: f64 = -0.49; /* min double value that fits in an unsigned short */
pub const DSHRT_MAX: f64 = 32767.49; /* max double value that fits in a short */
pub const DSHRT_MIN: f64 = -32768.49; /* min double value that fits in a short */

pub const DLONG_MAX: f64 = 9.223_372_036_854_775E18; /* max double value  long */
pub const DLONG_MIN: f64 = -9.223_372_036_854_775E18; /* min double value  long */
pub const DULONG_MAX: f64 = 1.844_674_407_370_955E19; /* max double value  ulong */

pub const DULONG_MIN: f64 = -0.49; /* min double value that fits in an unsigned long */
pub const DULONGLONG_MAX: f64 = 18446744073709551615.; /* max unsigned  longlong */
pub const DULONGLONG_MIN: f64 = -0.49;
pub const DLONGLONG_MAX: f64 = 9.223_372_036_854_776E18; /* max double value  longlong */
pub const DLONGLONG_MIN: f64 = -9.223_372_036_854_776E18; /* min double value  longlong */
pub const DUINT_MAX: f64 = 4294967295.49; /* max dbl that fits in a unsigned 4-byte int */
pub const DUINT_MIN: f64 = -0.49; /* min dbl that fits in an unsigned 4-byte int */
pub const DINT_MAX: f64 = 2147483647.49; /* max double value that fits in a 4-byte int */
pub const DINT_MIN: f64 = -2147483648.49; /* min double value that fits in a 4-byte int */

pub const UINT64_MAX: u64 = 18446744073709551615; /* max unsigned 64-bit integer */

pub const UINT32_MAX: u32 = 4294967295; /* max unsigned 32-bit integer */
pub const INT32_MAX: i32 = 2147483647; /* max 32-bit integer */
pub const INT32_MIN: i32 = -INT32_MAX - 1; /* min 32-bit integer */

pub const COMPRESS_NULL_VALUE: i32 = -2147483647;
pub const N_RANDOM: usize = 10000; /* DO NOT CHANGE THIS;  used when quantizing real numbers */
