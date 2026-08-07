/* Transpiled from cfitsio/utilities/fverify.h, plus the small amount of
   infrastructure that C gets for free from <stdio.h>/<ctype.h>/<stdlib.h>.

   The C globals declared here (errmes, comm, prhead, testdata, ...) are file
   statics / externs in the C.  Buffers that are only ever "printf into it, then
   immediately consume it" (errmes, comm, temp) become locals at each use site,
   which is semantically identical and avoids global mutable state.  The truly
   persistent globals (the flags and counters) live in `flags` below.  */

use std::cell::RefCell;
use std::fmt::Display;
use std::io::Write;
use std::sync::atomic::{AtomicI32, AtomicI64, Ordering};

use rsfitsio::c_types::{c_char, c_int, c_long, c_ulong};
use rsfitsio::fitsio::{FLEN_KEYWORD, FLEN_VALUE, LONGLONG};

pub const MAXERRORS: c_int = 200;
pub const MAXWRNS: c_int = 200;

/* static char errmes[256]; */
pub const ERRMES_LEN: usize = 256;
/* static char comm[FLEN_FILENAME+6]; */
pub const COMM_LEN: usize = rsfitsio::fitsio::FLEN_FILENAME + 6;

/********************************
*				*
*       Keywords 		*
*				*
********************************/

#[derive(Clone, Copy, PartialEq, Eq, Debug, Default)]
pub enum kwdtyp {
    STR_KEY, /* string   key */
    LOG_KEY, /* Logical key */
    INT_KEY, /* Integer key */
    FLT_KEY, /* Float key   */
    CMI_KEY, /* Complex integer key */
    CMF_KEY, /* Complex float key */
    COM_KEY, /* history, comment, "  ", and end */
    #[default]
    UNKNOWN, /* Unknown types */
}

/* error number masks of  the keyword test */
pub const BAD_STR: c_ulong = 0x0001;
pub const NO_TRAIL_QUOTE: c_ulong = 0x0002;
pub const BAD_NUM: c_ulong = 0x0004;
pub const LOWCASE_EXPO: c_ulong = 0x0008;
pub const NO_TRAIL_PAREN: c_ulong = 0x0010;
pub const NO_COMMA: c_ulong = 0x0020;
pub const TOO_MANY_COMMA: c_ulong = 0x0040;
pub const BAD_REAL: c_ulong = 0x0080;
pub const BAD_IMG: c_ulong = 0x0100;
pub const BAD_LOGICAL: c_ulong = 0x0200;
pub const NO_START_SLASH: c_ulong = 0x0400;
pub const BAD_COMMENT: c_ulong = 0x0800;
pub const UNKNOWN_TYPE: c_ulong = 0x1000;

/* Number of possible WCS descriptions to check.*/
/* 1 for the primary + 26 for [A-Z] suffix. */
pub const NWCSDESCR: usize = 27;

/* keyword structure */
#[derive(Clone)]
pub struct FitsKey {
    pub kname: [c_char; FLEN_KEYWORD], /* fits keyword name */
    pub ktype: kwdtyp,                 /* fits keyword type */
    pub kvalue: [c_char; FLEN_VALUE],  /* fits keyword name */
    pub kindex: c_int,                 /* position at the header */
    pub goodkey: c_int,                /* good keyword flag (=1 good)*/
}

impl Default for FitsKey {
    fn default() -> Self {
        Self {
            kname: [0; FLEN_KEYWORD],
            ktype: kwdtyp::UNKNOWN,
            kvalue: [0; FLEN_VALUE],
            kindex: 0,
            goodkey: 0,
        }
    }
}

/********************************
*				*
*       Headers  		*
*				*
********************************/
pub struct FitsHdu {
    pub hdutype: c_int,                /* hdutype */
    pub hdunum: c_int,                 /* hdunum  */
    pub isgroup: c_int,                /* random group flag */
    pub istilecompressed: c_int,       /* tile compressed image */
    pub gcount: c_int,                 /* gcount  */
    pub pcount: LONGLONG,              /* pcount  */
    pub bitpix: c_int,                 /* pix number */
    pub naxis: c_int,                  /* number of the axis,used for image array*/
    pub naxes: Vec<LONGLONG>,          /* dimension of each axis,used for image array*/
    pub ncols: c_int,                  /* number of the columns, used for image only*/
    pub extname: [c_char; FLEN_VALUE], /* EXTENSION NAME */
    pub extver: c_int,                 /* extension version */
    pub datamax: Vec<Vec<c_char>>,     /* strings for the maximum of the data in a column */
    pub datamin: Vec<Vec<c_char>>,     /* strings for the minimum of the data in a column */
    pub tnull: Vec<Vec<c_char>>,       /* number of NULL values */
    pub nkeys: c_int,                  /* number of keys */
    pub tkeys: c_int,                  /* total of the keys tested*/
    pub heap: c_int,                   /* heap */
    pub kwds: Vec<FitsKey>,            /* keywords list starting from the
                                       last NAXISn keyword. The array
                                       is sorted in the ascending alphabetical
                                       order of keyword names. The last keyword END
                                       and commentary keywords are  excluded.
                                       The total number of element, tkey, is
                                       nkeys - 4 - naxis - ncomm. */
    pub use_longstr: c_int, /* flag indicates that the long string
                            convention is used */
}

impl Default for FitsHdu {
    fn default() -> Self {
        Self {
            hdutype: 0,
            hdunum: 0,
            isgroup: 0,
            istilecompressed: 0,
            gcount: 0,
            pcount: 0,
            bitpix: 0,
            naxis: 0,
            naxes: Vec::new(),
            ncols: 0,
            extname: [0; FLEN_VALUE],
            extver: 0,
            datamax: Vec::new(),
            datamin: Vec::new(),
            tnull: Vec::new(),
            nkeys: 0,
            tkeys: 0,
            heap: 0,
            kwds: Vec::new(),
            use_longstr: 0,
        }
    }
}

pub struct ColName {
    pub name: Vec<c_char>, /* column name */
    pub index: c_int,      /* column index */
}

/********************************
*				*
*       Files   		*
*				*
********************************/
#[derive(Clone)]
pub struct HduName {
    pub hdutype: c_int,                /* hdutype */
    pub hdunum: c_int,                 /* hdunum  */
    pub extname: [c_char; FLEN_VALUE], /* extension name, used for extension*/
    pub extver: c_int,                 /* extension version, used for extension */
    pub errnum: c_int,                 /* number of errors in this hdu */
    pub wrnno: c_int,                  /* number of warnning in this hdu */
}

impl Default for HduName {
    fn default() -> Self {
        Self {
            hdutype: 0,
            hdunum: 0,
            extname: [0; FLEN_VALUE],
            extver: 0,
            errnum: 0,
            wrnno: 0,
        }
    }
}

/*===========================================================================
 *  The `extern int' globals from fverify.h / ftverify.c.
 *==========================================================================*/

macro_rules! int_global {
    ($get:ident, $set:ident, $store:ident, $init:expr) => {
        static $store: AtomicI32 = AtomicI32::new($init);
        #[inline]
        pub fn $get() -> c_int {
            $store.load(Ordering::Relaxed)
        }
        #[inline]
        pub fn $set(v: c_int) {
            $store.store(v, Ordering::Relaxed);
        }
    };
}

int_global!(err_report, set_err_report, ERR_REPORT, 0);
int_global!(prhead, set_prhead, PRHEAD, 0);
int_global!(prstat, set_prstat, PRSTAT, 1);
int_global!(testdata, set_testdata, TESTDATA, 1);
int_global!(testcsum, set_testcsum, TESTCSUM, 1);
int_global!(testfill, set_testfill, TESTFILL, 1);
int_global!(heasarc_conv, set_heasarc_conv, HEASARC_CONV, 1);
int_global!(testhierarch, set_testhierarch, TESTHIERARCH, 0);
int_global!(totalhdu, set_totalhdu, TOTALHDU, 0);

/* long totalerr, totalwrn;  (ftverify.c) */
static TOTALERR: AtomicI64 = AtomicI64::new(0);
static TOTALWRN: AtomicI64 = AtomicI64::new(0);

pub fn totalerr() -> i64 {
    TOTALERR.load(Ordering::Relaxed)
}
pub fn totalwrn() -> i64 {
    TOTALWRN.load(Ordering::Relaxed)
}
pub fn add_totalerr(v: i64) {
    TOTALERR.fetch_add(v, Ordering::Relaxed);
}
pub fn add_totalwrn(v: i64) {
    TOTALWRN.fetch_add(v, Ordering::Relaxed);
}

/*===========================================================================
 *  `FILE *out'
 *
 *  In the standalone fitsverify the report stream is always stdout, stderr or
 *  NULL, but ftverify_work() also has a branch that fopen()s a log file, so
 *  the File variant is carried too.  `Out' stays Copy (like the C FILE*) by
 *  keeping the actual handle in a thread-local.
 *==========================================================================*/

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Out {
    Null,
    Stdout,
    Stderr,
    File,
}

thread_local! {
    static OUTFILE: RefCell<Option<std::io::BufWriter<std::fs::File>>> =
        const { RefCell::new(None) };
}

pub fn out_set_file(f: std::fs::File) {
    OUTFILE.with(|o| *o.borrow_mut() = Some(std::io::BufWriter::new(f)));
}

pub fn out_close_file() {
    OUTFILE.with(|o| {
        if let Some(f) = o.borrow_mut().as_mut() {
            let _ = f.flush();
        }
        *o.borrow_mut() = None;
    });
}

/* fwrite() of raw bytes to `out'; a no-op for the NULL stream. */
pub fn fwrite_out(out: Out, b: &[u8]) {
    match out {
        Out::Null => {}
        Out::Stdout => {
            let so = std::io::stdout();
            let mut l = so.lock();
            let _ = l.write_all(b);
        }
        Out::Stderr => {
            let se = std::io::stderr();
            let mut l = se.lock();
            let _ = l.write_all(b);
        }
        Out::File => OUTFILE.with(|o| {
            if let Some(f) = o.borrow_mut().as_mut() {
                let _ = f.write_all(b);
            }
        }),
    }
}

pub fn fflush_out(out: Out) {
    match out {
        Out::Stdout => {
            let _ = std::io::stdout().flush();
        }
        Out::Stderr => {
            let _ = std::io::stderr().flush();
        }
        Out::File => OUTFILE.with(|o| {
            if let Some(f) = o.borrow_mut().as_mut() {
                let _ = f.flush();
            }
        }),
        Out::Null => {}
    }
}

/*===========================================================================
 *  printf-alikes.
 *
 *  FITS cards may hold arbitrary bytes and fitsverify prints them back in its
 *  diagnostics, so every message is built as a byte string rather than going
 *  through Rust's UTF-8 `String'.  `spf!' is sprintf(), `pf!' is printf(), and
 *  each argument is one printf conversion:
 *
 *      C:     sprintf(errmes, "Keyword #%d, %s is duplicated.", kpos, kname);
 *      Rust:  spf!(errmes; "Keyword #", kpos, ", ", CS(&kname), " is duplicated.");
 *==========================================================================*/

pub trait Put {
    fn put(&self, v: &mut Vec<u8>);
}

impl Put for &str {
    fn put(&self, v: &mut Vec<u8>) {
        v.extend_from_slice(self.as_bytes());
    }
}

macro_rules! put_via_display {
    ($($t:ty),*) => { $(
        impl Put for $t {
            fn put(&self, v: &mut Vec<u8>) { v.extend_from_slice(self.to_string().as_bytes()); }
        }
    )* };
}
put_via_display!(i8, u8, i16, u16, i32, u32, i64, u64, isize, usize);

/* %s of a NUL-terminated c_char buffer */
pub struct CS<'a>(pub &'a [c_char]);

impl Put for CS<'_> {
    fn put(&self, v: &mut Vec<u8>) {
        for &c in self.0 {
            if c == 0 {
                break;
            }
            v.push(c as u8);
        }
    }
}

/* %c */
pub struct CHR(pub c_char);

impl Put for CHR {
    fn put(&self, v: &mut Vec<u8>) {
        v.push(self.0 as u8);
    }
}

/* %[-]<width>[.<prec>]s of a NUL-terminated c_char buffer.
   `width' < 0 means left justified, as in "%-20s". */
pub struct CSW<'a>(pub &'a [c_char], pub i32, pub Option<usize>);

impl Put for CSW<'_> {
    fn put(&self, v: &mut Vec<u8>) {
        let mut b: Vec<u8> = Vec::new();
        CS(self.0).put(&mut b);
        if let Some(p) = self.2 {
            b.truncate(p);
        }
        pad(v, &b, self.1);
    }
}

/* %[-]<width>[.<prec>]s of a Rust &str */
pub struct SW<'a>(pub &'a str, pub i32, pub Option<usize>);

impl Put for SW<'_> {
    fn put(&self, v: &mut Vec<u8>) {
        let mut b: Vec<u8> = self.0.as_bytes().to_vec();
        if let Some(p) = self.2 {
            b.truncate(p);
        }
        pad(v, &b, self.1);
    }
}

/* %[-]<width>d */
pub struct DW<T: Display>(pub T, pub i32);

impl<T: Display> Put for DW<T> {
    fn put(&self, v: &mut Vec<u8>) {
        pad(v, self.0.to_string().as_bytes(), self.1);
    }
}

fn pad(v: &mut Vec<u8>, b: &[u8], width: i32) {
    let w = width.unsigned_abs() as usize;
    if b.len() >= w {
        v.extend_from_slice(b);
    } else if width < 0 {
        v.extend_from_slice(b);
        v.extend(std::iter::repeat_n(b' ', w - b.len()));
    } else {
        v.extend(std::iter::repeat_n(b' ', w - b.len()));
        v.extend_from_slice(b);
    }
}

/* Formats the arguments and stores the result NUL-terminated in `dst'.
   C's sprintf() would run off the end of these fixed buffers; we truncate. */
#[macro_export]
macro_rules! spf {
    ($dst:expr; $($x:expr),* $(,)?) => {{
        let mut __v: Vec<u8> = Vec::new();
        $( $crate::common::Put::put(&$x, &mut __v); )*
        $crate::common::set_cstr(&mut $dst[..], &__v);
    }};
}

/* Formats the arguments and appends them, NUL-terminated, to `dst' (strcat). */
#[macro_export]
macro_rules! scat {
    ($dst:expr; $($x:expr),* $(,)?) => {{
        let mut __v: Vec<u8> = Vec::new();
        $( $crate::common::Put::put(&$x, &mut __v); )*
        $crate::common::cat_cstr(&mut $dst[..], &__v);
    }};
}

/* printf() */
#[macro_export]
macro_rules! pf {
    ($out:expr; $($x:expr),* $(,)?) => {{
        let mut __v: Vec<u8> = Vec::new();
        $( $crate::common::Put::put(&$x, &mut __v); )*
        $crate::common::fwrite_out($out, &__v);
    }};
}

pub fn set_cstr(dst: &mut [c_char], src: &[u8]) {
    let n = std::cmp::min(src.len(), dst.len() - 1);
    for i in 0..n {
        dst[i] = src[i] as c_char;
    }
    dst[n] = 0;
}

pub fn cat_cstr(dst: &mut [c_char], src: &[u8]) {
    let l = cstrlen(dst);
    let n = std::cmp::min(src.len(), dst.len() - 1 - l);
    for i in 0..n {
        dst[l + i] = src[i] as c_char;
    }
    dst[l + n] = 0;
}

/* strlen() that tolerates an unterminated buffer (returns its length). */
pub fn cstrlen(s: &[c_char]) -> usize {
    s.iter().position(|&c| c == 0).unwrap_or(s.len())
}

/* The bytes of a NUL-terminated c_char buffer. */
pub fn cbytes(s: &[c_char]) -> &[u8] {
    let n = cstrlen(s);
    // SAFETY: c_char and u8 have the same size and alignment.
    unsafe { std::slice::from_raw_parts(s.as_ptr().cast::<u8>(), n) }
}

/*===========================================================================
 *  <string.h>
 *
 *  rsfitsio's wrappers.rs has strcpy_safe/strncpy_safe/strcmp_safe/strchr_safe,
 *  but they are not drop-in here: strlen_safe/strcpy_safe panic on a buffer
 *  with no NUL, strncpy_safe panics when the source is shorter than n,
 *  strchr_safe scans the whole slice rather than stopping at the NUL, and the
 *  strcmp_safe family subtracts signed c_char rather than glibc's unsigned
 *  char.  fitsverify is fed deliberately malformed headers, so it needs the
 *  non-panicking, glibc-shaped variants below.
 *==========================================================================*/

/* strcpy(); truncates instead of overrunning `dst' the way the C would. */
pub fn strcpy_c(dst: &mut [c_char], src: &[c_char]) {
    let n = std::cmp::min(cstrlen(src), dst.len() - 1);
    dst[..n].copy_from_slice(&src[..n]);
    dst[n] = 0;
}

/* strncpy(): at most n bytes, NUL-padding a short source. */
pub fn strncpy_c(dst: &mut [c_char], src: &[c_char], n: usize) {
    let n = std::cmp::min(n, dst.len());
    let mut i = 0;
    while i < n && i < src.len() && src[i] != 0 {
        dst[i] = src[i];
        i += 1;
    }
    while i < n {
        dst[i] = 0;
        i += 1;
    }
}

/* strcmp() / strncmp() with glibc's unsigned-char comparison, tolerating a
   buffer that runs to its end without a NUL. */
pub fn strcmp_c(a: &[c_char], b: &[c_char]) -> c_int {
    strncmp_c(a, b, usize::MAX)
}

pub fn strncmp_c(a: &[c_char], b: &[c_char], n: usize) -> c_int {
    let mut i = 0usize;
    while i < n {
        let ca = if i < a.len() { a[i] as u8 } else { 0 };
        let cb = if i < b.len() { b[i] as u8 } else { 0 };
        if ca != cb {
            return ca as c_int - cb as c_int;
        }
        if ca == 0 {
            return 0;
        }
        i += 1;
    }
    0
}

/* strchr() as an offset; stops at the NUL, as C does. */
pub fn chr_pos(s: &[c_char], c: u8) -> Option<usize> {
    let n = cstrlen(s);
    s[..n].iter().position(|&x| x as u8 == c)
}

/* strchr() != NULL */
pub fn chr_in(s: &[c_char], c: u8) -> bool {
    chr_pos(s, c).is_some()
}

/* A NUL-terminated owned copy (what the C aliases with a `char *'). */
pub fn cvec(s: &[c_char]) -> Vec<c_char> {
    let n = cstrlen(s);
    let mut v = Vec::with_capacity(n + 1);
    v.extend_from_slice(&s[..n]);
    v.push(0);
    v
}

/*===========================================================================
 *  <ctype.h>.  These mirror glibc in the "C" locale, where every class is
 *  false for the negative values a signed `char' can hold.
 *==========================================================================*/

#[inline]
pub fn isprint_c(c: c_char) -> bool {
    (0x20..=0x7e).contains(&(c as u8))
}

#[inline]
pub fn isspace_c(c: c_char) -> bool {
    matches!(c as u8, b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r')
}

#[inline]
pub fn isdigit_c(c: c_char) -> bool {
    (b'0'..=b'9').contains(&(c as u8))
}

#[inline]
pub fn isupper_c(c: c_char) -> bool {
    (b'A'..=b'Z').contains(&(c as u8))
}

#[inline]
pub fn isascii_c(c: c_char) -> bool {
    (c as u8) < 128
}

#[inline]
pub fn toupper_c(c: c_char) -> c_char {
    if (b'a'..=b'z').contains(&(c as u8)) {
        (c as u8 - 32) as c_char
    } else {
        c
    }
}

/*===========================================================================
 *  <stdlib.h> string conversions.
 *
 *  fitsverify tests strtol()'s results against LONG_MAX / LONG_MIN to detect
 *  malformed values, so the saturating overflow behaviour has to be kept.
 *==========================================================================*/

/* strtol(s, &end, 10); returns (value, index one past the last char consumed) */
pub fn strtol_c(s: &[c_char]) -> (c_long, usize) {
    let mut i = 0usize;
    let n = cstrlen(s);

    while i < n && isspace_c(s[i]) {
        i += 1;
    }

    let mut neg = false;
    if i < n && (s[i] as u8 == b'+' || s[i] as u8 == b'-') {
        neg = s[i] as u8 == b'-';
        i += 1;
    }

    let start = i;
    let mut acc: i128 = 0;
    let mut overflow = false;
    while i < n && isdigit_c(s[i]) {
        if !overflow {
            acc = acc * 10 + (s[i] as u8 - b'0') as i128;
            if acc > (c_long::MAX as i128) + 1 {
                overflow = true;
            }
        }
        i += 1;
    }

    if i == start {
        /* no digits: strtol() returns 0 and endptr == the original string */
        return (0, 0);
    }

    let v = if neg { -acc } else { acc };
    /* C's strtol() saturates towards the sign of the value it was parsing. */
    let v = if overflow {
        if neg { c_long::MIN } else { c_long::MAX }
    } else if v > c_long::MAX as i128 {
        c_long::MAX
    } else if v < c_long::MIN as i128 {
        c_long::MIN
    } else {
        v as c_long
    };
    (v, i)
}

/* strtoll(s, NULL, 10) */
pub fn strtoll_c(s: &[c_char]) -> LONGLONG {
    strtol_c(s).0 as LONGLONG
}

/* strtod(s, NULL) / atof(s) */
pub fn strtod_c(s: &[c_char]) -> f64 {
    let b = cbytes(s);
    let mut i = 0usize;
    while i < b.len() && (b[i] as char).is_whitespace() {
        i += 1;
    }
    let start = i;
    if i < b.len() && (b[i] == b'+' || b[i] == b'-') {
        i += 1;
    }
    while i < b.len() && b[i].is_ascii_digit() {
        i += 1;
    }
    if i < b.len() && b[i] == b'.' {
        i += 1;
        while i < b.len() && b[i].is_ascii_digit() {
            i += 1;
        }
    }
    if i < b.len() && (b[i] == b'e' || b[i] == b'E' || b[i] == b'd' || b[i] == b'D') {
        let save = i;
        i += 1;
        if i < b.len() && (b[i] == b'+' || b[i] == b'-') {
            i += 1;
        }
        if i < b.len() && b[i].is_ascii_digit() {
            while i < b.len() && b[i].is_ascii_digit() {
                i += 1;
            }
        } else {
            i = save;
        }
    }
    if i == start {
        return 0.0;
    }
    /* FITS allows 'D'/'d' as the exponent marker; Rust's parser wants 'e'. */
    let mut t = String::from_utf8_lossy(&b[start..i]).into_owned();
    if t.contains(['d', 'D']) {
        t = t.replace(['d', 'D'], "e");
    }
    t.parse::<f64>().unwrap_or(0.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn to_buf<const N: usize>(s: &str) -> [c_char; N] {
        let mut b = [0 as c_char; N];
        for (i, &c) in s.as_bytes().iter().enumerate() {
            b[i] = c as c_char;
        }
        b
    }

    /* strtol() saturates at LONG_MAX/LONG_MIN on overflow, and fitsverify
       tests its results against those to spot malformed TDISPn widths, so the
       saturation has to survive the port. */
    #[test]
    fn test_strtol_c() {
        assert_eq!(strtol_c(&to_buf::<32>("123abc")), (123, 3));
        assert_eq!(strtol_c(&to_buf::<32>("  -45")), (-45, 5));
        assert_eq!(strtol_c(&to_buf::<32>("+7")), (7, 2));
        assert_eq!(strtol_c(&to_buf::<32>("0")), (0, 1));
        /* no conversion: C leaves endptr == nptr and returns 0 */
        assert_eq!(strtol_c(&to_buf::<32>("")), (0, 0));
        assert_eq!(strtol_c(&to_buf::<32>("abc")), (0, 0));
        assert_eq!(strtol_c(&to_buf::<32>("+")), (0, 0));
        /* overflow saturates rather than wrapping */
        assert_eq!(strtol_c(&to_buf::<32>("99999999999999999999")).0, c_long::MAX);
        assert_eq!(strtol_c(&to_buf::<32>("-99999999999999999999")).0, c_long::MIN);
        /* the endptr is what the PC/CD alt-suffix test keys off */
        assert_eq!(strtol_c(&to_buf::<32>("12_3")), (12, 2));
    }

    #[test]
    fn test_strtod_c() {
        assert_eq!(strtod_c(&to_buf::<32>("1.5")), 1.5);
        assert_eq!(strtod_c(&to_buf::<32>("-0.25")), -0.25);
        assert_eq!(strtod_c(&to_buf::<32>("1.5E2")), 150.0);
        /* FITS allows 'D' as the exponent marker */
        assert_eq!(strtod_c(&to_buf::<32>("1.5D2")), 150.0);
        assert_eq!(strtod_c(&to_buf::<32>("abc")), 0.0);
        assert_eq!(strtod_c(&to_buf::<32>("")), 0.0);
        /* a trailing exponent marker with no digits is not consumed */
        assert_eq!(strtod_c(&to_buf::<32>("2.0E")), 2.0);
    }

    /* glibc's C locale reports false for every class on the negative values a
       signed char can hold, and c_char is unsigned on some targets, so these
       must not depend on the platform's char signedness. */
    #[test]
    fn test_ctype_helpers() {
        assert!(!isprint_c(0x1f as c_char));
        assert!(isprint_c(0x20 as c_char));
        assert!(isprint_c(0x7e as c_char));
        assert!(!isprint_c(0x7f as c_char));
        assert!(!isprint_c(0x80u8 as c_char));
        assert!(!isprint_c(0xffu8 as c_char));

        assert!(isascii_c(0x7f as c_char));
        assert!(!isascii_c(0x80u8 as c_char));

        assert!(isspace_c(b' ' as c_char));
        assert!(isspace_c(b'\t' as c_char));
        assert!(!isspace_c(b'x' as c_char));
        assert!(!isspace_c(0xa0u8 as c_char));

        assert!(isdigit_c(b'0' as c_char));
        assert!(isdigit_c(b'9' as c_char));
        assert!(!isdigit_c(b'a' as c_char));

        assert!(isupper_c(b'A' as c_char));
        assert!(!isupper_c(b'a' as c_char));

        assert_eq!(toupper_c(b'a' as c_char), b'A' as c_char);
        assert_eq!(toupper_c(b'Z' as c_char), b'Z' as c_char);
        assert_eq!(toupper_c(b'_' as c_char), b'_' as c_char);
    }

    #[test]
    fn test_string_helpers() {
        let mut d = [0 as c_char; 8];
        strcpy_c(&mut d, &to_buf::<16>("abc"));
        assert_eq!(cstrlen(&d), 3);
        /* strcpy_c truncates where the C would overrun the destination */
        strcpy_c(&mut d, &to_buf::<32>("abcdefghijkl"));
        assert_eq!(cstrlen(&d), 7);

        /* strncpy NUL-pads a short source */
        let mut e = [b'x' as c_char; 8];
        strncpy_c(&mut e, &to_buf::<16>("ab"), 5);
        assert_eq!(&e[..5], &[b'a' as c_char, b'b' as c_char, 0, 0, 0]);
        assert_eq!(e[5], b'x' as c_char);

        assert_eq!(strcmp_c(&to_buf::<16>("abc"), &to_buf::<16>("abc")), 0);
        assert!(strcmp_c(&to_buf::<16>("abc"), &to_buf::<16>("abd")) < 0);
        assert!(strcmp_c(&to_buf::<16>("abd"), &to_buf::<16>("abc")) > 0);
        assert_eq!(strncmp_c(&to_buf::<16>("abcZ"), &to_buf::<16>("abcY"), 3), 0);

        /* strchr stops at the NUL, so bytes past it are not found */
        let mut s = to_buf::<16>("abc");
        s[5] = b'z' as c_char;
        assert_eq!(chr_pos(&s, b'b'), Some(1));
        assert_eq!(chr_pos(&s, b'z'), None);
        assert!(chr_in(&s, b'c'));
        assert!(!chr_in(&s, b'q'));

        assert_eq!(cvec(&to_buf::<16>("hi")), vec![104, 105, 0]);
    }

    /* The message builder has to pass bytes through untouched: fitsverify
       echoes malformed cards, which are routinely not valid UTF-8. */
    #[test]
    fn test_spf_is_byte_exact() {
        let mut b = [0 as c_char; 64];
        spf!(b; "n=", 42, " c=", CHR(0xfeu8 as c_char), "!");
        assert_eq!(cbytes(&b), b"n=42 c=\xfe!");

        /* %-4d / %-20s / %.8s style conversions */
        let mut w = [0 as c_char; 64];
        spf!(w; "[", DW(7, -4), "][", DW(7, 4), "]");
        assert_eq!(cbytes(&w), b"[7   ][   7]");

        let name = to_buf::<16>("ABCDEFGHIJ");
        let mut p = [0 as c_char; 64];
        spf!(p; CSW(&name, 0, Some(8)), "|", CSW(&name, -12, None), "|");
        assert_eq!(cbytes(&p), b"ABCDEFGH|ABCDEFGHIJ  |");

        /* sprintf into a fixed buffer truncates rather than overrunning it */
        let mut t = [0 as c_char; 5];
        spf!(t; "abcdefgh");
        assert_eq!(cbytes(&t), b"abcd");

        /* strcat */
        let mut c = [0 as c_char; 16];
        spf!(c; "ab");
        scat!(c; "cd", 9);
        assert_eq!(cbytes(&c), b"abcd9");
    }
}
