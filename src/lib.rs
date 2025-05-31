#![feature(vec_into_raw_parts)]
#![feature(allocator_api)]
// These are required for printf in Relibc
#![feature(c_variadic)]
#![allow(
    non_camel_case_types,
    non_snake_case,
    dead_code,
    unused_variables,
    unused_assignments,
    clippy::missing_safety_doc,
    unreachable_code,
    clippy::too_many_arguments,
    clippy::needless_range_loop,
    clippy::manual_range_contains
)]
#![deny(deprecated)]

pub mod c_types;

pub mod aliases;
pub mod buffers;
pub mod checksum;
pub mod edithdu;
pub mod fitsio;
pub mod fitsio2;
pub mod scalnull;

//#[cfg(feature = "bzip2")]
//mod bzip2;

pub mod cfileio;
pub mod drvrfile;
pub mod drvrgsiftp;
pub mod drvrmem;
pub mod drvrnet;

#[cfg(feature = "shared_mem")]
pub mod drvrsmem;

pub mod editcol;
pub mod eval_defs;
pub mod eval_f;
pub mod eval_l;
pub mod eval_y;
pub mod fits_hcompress;
pub mod fits_hdecompress;
pub mod fitscore;
pub mod getcol;
pub mod getcolb;
pub mod getcold;
pub mod getcole;
pub mod getcoli;
pub mod getcolj;
pub mod getcolk;
pub mod getcoll;
pub mod getcols;
pub mod getcolsb;
pub mod getcolui;
pub mod getcoluj;
pub mod getcoluk;
pub mod getkey;
pub mod group;
pub mod grparser;
pub mod histo;
pub mod imcompress;
pub mod iraffits;
pub mod modkey;
pub mod putcol;
pub mod putcolb;
pub mod putcold;
pub mod putcole;
pub mod putcoli;
pub mod putcolj;
pub mod putcolk;
pub mod putcoll;
pub mod putcols;
pub mod putcolsb;
pub mod putcolu;
pub mod putcolui;
pub mod putcoluj;
pub mod putcoluk;
pub mod putkey;
pub mod quantize;
pub mod region;
pub mod relibc;
pub mod simplerng;
pub mod swapproc;
pub mod testhelpers;
pub mod wcssub;
pub mod wcsutil;
pub mod wrappers;
pub mod zcompress;
pub mod zuncompress;

use std::{
    ffi::{CStr, CString},
    marker::PhantomData,
    str::FromStr,
    sync::{Mutex, MutexGuard},
};

use bytemuck::cast_slice;
use fitsio::{
    FLEN_VALUE, LONGLONG, TBIT, TBYTE, TCOMPLEX, TDBLCOMPLEX, TDOUBLE, TFLOAT, TINT, TLOGICAL,
    TLONG, TLONGLONG, TSBYTE, TSHORT, TSTRING, TUINT, TULONG, TULONGLONG, TUSHORT, ULONGLONG,
};

use crate::c_types::*;

pub(crate) static MUTEX_LOCK: Mutex<bool> = Mutex::new(false);

pub(crate) fn FFLOCK<'a>() -> MutexGuard<'a, bool> {
    MUTEX_LOCK.lock().unwrap()
}

pub(crate) fn FFUNLOCK(p: MutexGuard<'_, bool>) {
    drop(p);
}

pub trait ToRaw {
    fn as_raw_mut(&mut self) -> *mut Self;
}

pub trait AsMutPtr<T> {
    fn as_mut_ptr(&self) -> *mut T;
}

impl<T> AsMutPtr<T> for Option<&mut [T]> {
    fn as_mut_ptr(&self) -> *mut T {
        match self {
            Some(v) => v.as_ptr() as *mut T, // UNSAFE
            None => std::ptr::null_mut(),
        }
    }
}

#[inline(always)]
pub fn bb(n: u8) -> c_char {
    n as c_char
}

#[macro_export]
macro_rules! int_snprintf {
    ($dst:expr, $len:expr, $($arg:tt)*) => {
        {
            let s = format!($($arg)*);
            let s_bytes = s.as_bytes();
            let mut s_len = s_bytes.len();

            s_len = cmp::min($len-1, s_len);

            let w = cast_slice_mut::<c_char, u8>(&mut $dst[..s_len]);
            w.copy_from_slice(&s_bytes[..s_len]);
            $dst[s_len] = 0; // null-terminate

            s_len as isize
        }
    };
}

#[macro_export]
macro_rules! slice_to_str {
    ($e:expr) => {
        CStr::from_bytes_until_nul(cast_slice($e))
            .unwrap()
            .to_str()
            .unwrap()
    };
}

#[macro_export]
macro_rules! cs {
    ($e: expr) => {
        cast_slice($e.to_bytes_with_nul())
    };
}

#[macro_export]
macro_rules! nullable_slice_cstr {
    ($e: ident) => {
        let $e: Option<&[c_char]> = match $e.is_null() {
            true => None,
            false => Some(cast_slice(CStr::from_ptr($e).to_bytes_with_nul())),
        };
    };
}

#[macro_export]
macro_rules! nullable_slice_cstr_mut {
    ($e: ident) => {
        let mut $e: Option<&mut [c_char]> = match $e.is_null() {
            true => None,
            false => {
                let _c = CStr::from_ptr($e).to_bytes_with_nul();
                let _l = _c.len();

                Some(slice::from_raw_parts_mut($e, _l))
            }
        };
    };
}

#[macro_export]
macro_rules! raw_to_slice {
    ($e: ident) => {
        let $e: &[c_char] = cast_slice(CStr::from_ptr($e).to_bytes_with_nul());
    };
}

pub(crate) struct TKeywords<'a> {
    tfields: c_int,              /* I - number of columns in the table           */
    ttype: *const *const c_char, /* I - name of each column                      */
    tform: *const *const c_char, /* I - value of TFORMn keyword for each column  */
    tunit: *const *const c_char, /* I - value of TUNITn keyword for each column  */
    marker: PhantomData<&'a ()>,
}

impl<'a> TKeywords<'a> {
    pub fn new(
        tfields: c_int,
        ttype: *const *const c_char,
        tform: *const *const c_char,
        tunit: *const *const c_char,
    ) -> Self {
        TKeywords {
            tfields,
            ttype,
            tform,
            tunit,
            marker: PhantomData,
        }
    }

    pub unsafe fn tkeywords_to_vecs(
        &'a self,
    ) -> (
        Vec<Option<&'a [c_char]>>,
        Vec<&'a [c_char]>,
        Option<Vec<Option<&'a [c_char]>>>,
    ) {
        unsafe {
            // Convert ttype to rust types
            let ttype = core::slice::from_raw_parts(self.ttype, self.tfields as usize);
            let mut v_ttype = Vec::new();

            for item in ttype {
                let ttype_item = if item.is_null() {
                    None
                } else {
                    Some(core::slice::from_raw_parts(*item, FLEN_VALUE))
                };
                v_ttype.push(ttype_item);
            }

            // Convert tform to rust types
            let tform = core::slice::from_raw_parts(self.tform, self.tfields as usize);
            let mut v_tform = Vec::new();

            for item in tform {
                let tform_item = core::slice::from_raw_parts(*item, FLEN_VALUE);
                v_tform.push(tform_item);
            }

            // Convert tunit to rust types
            let mut v_tunit = Vec::new();
            let out_tunit = if self.tunit.is_null() {
                None
            } else {
                let tunit = core::slice::from_raw_parts(self.tunit, self.tfields as usize);

                for item in tunit {
                    let tunit_item = if item.is_null() {
                        None
                    } else {
                        Some(core::slice::from_raw_parts(*item, FLEN_VALUE))
                    };
                    v_tunit.push(tunit_item);
                }
                Some(v_tunit)
            };

            (v_ttype, v_tform, out_tunit)
        }
    }
}

pub(crate) fn calculate_subsection_length(blc: &[c_long], trc: &[c_long], inc: &[c_long]) -> usize {
    assert!(blc.len() == trc.len() && blc.len() == inc.len());

    let len = blc.len();
    let mut acc: usize = 1;
    for ii in 0..len {
        acc *= ((trc[ii] - blc[ii]) / inc[ii] + 1) as usize; // WARNING: This also could be wrong
    }
    acc
}

pub(crate) fn calculate_subsection_length_unit(blc: &[c_long], trc: &[c_long]) -> usize {
    assert!(blc.len() == trc.len());

    let len = blc.len();
    let mut acc: usize = 1;
    for ii in 0..len {
        acc *= ((trc[ii] - blc[ii]) + 1) as usize; // WARNING: This also could be wrong
    }
    acc
}

pub(crate) fn vecs_to_slices<T>(vecs: &[Vec<T>]) -> Vec<&[T]> {
    vecs.iter().map(Vec::as_slice).collect()
}

pub(crate) fn vecs_to_slices_mut<T>(vecs: &mut [Vec<T>]) -> Vec<&mut [T]> {
    vecs.iter_mut().map(Vec::as_mut_slice).collect()
}

#[derive(Debug, PartialEq, Clone, Copy)]
pub enum NullCheckType {
    None = 0,         // No null checking
    SetPixel = 1,     // Set undefined pixels as a null value
    SetNullArray = 2, // Set a null array = 1 for undefined pixels
}

#[derive(Debug, PartialEq, Clone)]
pub enum NullValue {
    Float(f32),
    Double(f64),
    Long(c_long),
    ULong(c_ulong),
    LONGLONG(LONGLONG),
    ULONGLONG(ULONGLONG),
    Int(c_int),
    UInt(c_uint),
    Short(c_short),
    UShort(c_ushort),
    Byte(i8),
    UByte(c_uchar),
    Logical(c_char),
    String(CString),
}

impl NullValue {
    pub fn get_value_as_f64(&self) -> f64 {
        match self {
            NullValue::Float(v) => *v as f64,
            NullValue::Double(v) => *v,
            NullValue::Long(v) => *v as f64,
            NullValue::ULong(v) => *v as f64,
            NullValue::LONGLONG(v) => *v as f64,
            NullValue::ULONGLONG(v) => *v as f64,
            NullValue::Int(v) => *v as f64,
            NullValue::UInt(v) => *v as f64,
            NullValue::Short(v) => *v as f64,
            NullValue::UShort(v) => *v as f64,
            NullValue::Byte(v) => *v as f64,
            NullValue::UByte(v) => *v as f64,
            NullValue::Logical(v) => *v as f64,
            _ => 0.0,
        }
    }

    pub fn from_raw_ptr(datatype: c_int, value: *const c_void) -> Option<Self> {
        if value.is_null() {
            return None;
        }

        match datatype {
            TFLOAT => Some(NullValue::Float(unsafe { *(value as *const f32) })),
            TDOUBLE => Some(NullValue::Double(unsafe { *(value as *const f64) })),
            TLONG => Some(NullValue::Long(unsafe { *(value as *const c_long) })),
            TULONG => Some(NullValue::ULong(unsafe { *(value as *const c_ulong) })),
            TLONGLONG => Some(NullValue::LONGLONG(unsafe { *(value as *const LONGLONG) })),
            TULONGLONG => Some(NullValue::ULONGLONG(unsafe {
                *(value as *const ULONGLONG)
            })),
            TINT => Some(NullValue::Int(unsafe { *(value as *const c_int) })),
            TUINT => Some(NullValue::UInt(unsafe { *(value as *const c_uint) })),
            TSHORT => Some(NullValue::Short(unsafe { *(value as *const c_short) })),
            TUSHORT => Some(NullValue::UShort(unsafe { *(value as *const c_ushort) })),
            TBYTE => Some(NullValue::UByte(unsafe { *(value as *const c_uchar) })),
            TSBYTE => Some(NullValue::Byte(unsafe { *(value as *const i8) })),
            TLOGICAL => Some(NullValue::Logical(unsafe { *(value as *const c_char) })),
            TSTRING => {
                let cstr = unsafe { CStr::from_ptr(value as *const c_char) };
                Some(NullValue::String(cstr.to_owned()))
            }
            _ => None, // Don't panic here, will be handled by the caller
        }
    }
}

#[derive(Debug, PartialEq)]
pub enum KeywordDatatype<'a> {
    TBYTE(&'a c_uchar),
    TSBYTE(&'a c_char),
    TSHORT(&'a c_short),
    TUSHORT(&'a c_ushort),
    TINT(&'a c_int),
    TUINT(&'a c_uint),
    TLONG(&'a c_long),
    TULONG(&'a c_ulong),
    TFLOAT(&'a f32),
    TDOUBLE(&'a f64),
    TSTRING(&'a [c_char]),
    TLOGICAL(&'a c_int),
    TCOMPLEX(&'a [f32; 2]),
    TDBLCOMPLEX(&'a [f64; 2]),
    TULONGLONG(&'a ULONGLONG),
    TLONGLONG(&'a LONGLONG),
    INVALID(c_int),
}

impl KeywordDatatype<'_> {
    pub fn from_datatype(datatype: c_int, value: *const c_void) -> Self {
        match datatype {
            TBYTE => KeywordDatatype::TBYTE(unsafe { &*(value as *const c_uchar) }),
            TSBYTE => KeywordDatatype::TSBYTE(unsafe { &*(value as *const c_char) }),
            TSHORT => KeywordDatatype::TSHORT(unsafe { &*(value as *const c_short) }),
            TUSHORT => KeywordDatatype::TUSHORT(unsafe { &*(value as *const c_ushort) }),
            TINT => KeywordDatatype::TINT(unsafe { &*(value as *const c_int) }),
            TUINT => KeywordDatatype::TUINT(unsafe { &*(value as *const c_uint) }),
            TLONG => KeywordDatatype::TLONG(unsafe { &*(value as *const c_long) }),
            TULONG => KeywordDatatype::TULONG(unsafe { &*(value as *const c_ulong) }),
            TFLOAT => KeywordDatatype::TFLOAT(unsafe { &*(value as *const f32) }),
            TDOUBLE => KeywordDatatype::TDOUBLE(unsafe { &*(value as *const f64) }),
            TSTRING => KeywordDatatype::TSTRING(unsafe {
                cast_slice(CStr::from_ptr(value as *const c_char).to_bytes_with_nul())
            }),
            TLOGICAL => KeywordDatatype::TLOGICAL(unsafe { &*(value as *const c_int) }),
            TCOMPLEX => KeywordDatatype::TCOMPLEX(unsafe { &*(value as *const [f32; 2]) }),
            TDBLCOMPLEX => KeywordDatatype::TDBLCOMPLEX(unsafe { &*(value as *const [f64; 2]) }),
            TULONGLONG => KeywordDatatype::TULONGLONG(unsafe { &*(value as *const ULONGLONG) }),
            TLONGLONG => KeywordDatatype::TLONGLONG(unsafe { &*(value as *const LONGLONG) }),
            _ => KeywordDatatype::INVALID(datatype),
        }
    }

    pub fn to_datatype_code(&self) -> c_int {
        match self {
            KeywordDatatype::TBYTE(_) => TBYTE,
            KeywordDatatype::TSBYTE(_) => TSBYTE,
            KeywordDatatype::TSHORT(_) => TSHORT,
            KeywordDatatype::TUSHORT(_) => TUSHORT,
            KeywordDatatype::TINT(_) => TINT,
            KeywordDatatype::TUINT(_) => TUINT,
            KeywordDatatype::TLONG(_) => TLONG,
            KeywordDatatype::TULONG(_) => TULONG,
            KeywordDatatype::TFLOAT(_) => TFLOAT,
            KeywordDatatype::TDOUBLE(_) => TDOUBLE,
            KeywordDatatype::TSTRING(_) => TSTRING,
            KeywordDatatype::TLOGICAL(_) => TLOGICAL,
            KeywordDatatype::TCOMPLEX(_) => TCOMPLEX,
            KeywordDatatype::TDBLCOMPLEX(_) => TDBLCOMPLEX,
            KeywordDatatype::TULONGLONG(_) => TULONGLONG,
            KeywordDatatype::TLONGLONG(_) => TLONGLONG,
            KeywordDatatype::INVALID(x) => *x,
        }
    }
}

#[derive(Debug, PartialEq)]
pub enum KeywordDatatypeMut<'a> {
    TBYTE(&'a mut c_uchar),
    TSBYTE(&'a mut c_char),
    TSHORT(&'a mut c_short),
    TUSHORT(&'a mut c_ushort),
    TINT(&'a mut c_int),
    TUINT(&'a mut c_uint),
    TLONG(&'a mut c_long),
    TULONG(&'a mut c_ulong),
    TFLOAT(&'a mut f32),
    TDOUBLE(&'a mut f64),
    TSTRING(&'a mut [c_char; FLEN_VALUE]),
    TLOGICAL(&'a mut c_int),
    TCOMPLEX(&'a mut [f32; 2]),
    TDBLCOMPLEX(&'a mut [f64; 2]),
    TULONGLONG(&'a mut ULONGLONG),
    TLONGLONG(&'a mut LONGLONG),
    INVALID(c_int),
}

impl KeywordDatatypeMut<'_> {
    pub fn from_datatype(datatype: c_int, value: *mut c_void) -> Self {
        match datatype {
            TBYTE => KeywordDatatypeMut::TBYTE(unsafe { &mut *(value as *mut c_uchar) }),
            TSBYTE => KeywordDatatypeMut::TSBYTE(unsafe { &mut *(value as *mut c_char) }),
            TSHORT => KeywordDatatypeMut::TSHORT(unsafe { &mut *(value as *mut c_short) }),
            TUSHORT => KeywordDatatypeMut::TUSHORT(unsafe { &mut *(value as *mut c_ushort) }),
            TINT => KeywordDatatypeMut::TINT(unsafe { &mut *(value as *mut c_int) }),
            TUINT => KeywordDatatypeMut::TUINT(unsafe { &mut *(value as *mut c_uint) }),
            TLONG => KeywordDatatypeMut::TLONG(unsafe { &mut *(value as *mut c_long) }),
            TULONG => KeywordDatatypeMut::TULONG(unsafe { &mut *(value as *mut c_ulong) }),
            TFLOAT => KeywordDatatypeMut::TFLOAT(unsafe { &mut *(value as *mut f32) }),
            TDOUBLE => KeywordDatatypeMut::TDOUBLE(unsafe { &mut *(value as *mut f64) }),
            TSTRING => {
                KeywordDatatypeMut::TSTRING(unsafe { &mut *(value as *mut [c_char; FLEN_VALUE]) })
            }
            TLOGICAL => KeywordDatatypeMut::TLOGICAL(unsafe { &mut *(value as *mut c_int) }),
            TCOMPLEX => KeywordDatatypeMut::TCOMPLEX(unsafe { &mut *(value as *mut [f32; 2]) }),
            TDBLCOMPLEX => {
                KeywordDatatypeMut::TDBLCOMPLEX(unsafe { &mut *(value as *mut [f64; 2]) })
            }
            TULONGLONG => {
                KeywordDatatypeMut::TULONGLONG(unsafe { &mut *(value as *mut ULONGLONG) })
            }
            TLONGLONG => KeywordDatatypeMut::TLONGLONG(unsafe { &mut *(value as *mut LONGLONG) }),
            _ => KeywordDatatypeMut::INVALID(datatype),
        }
    }

    pub fn to_datatype_code(&self) -> c_int {
        match self {
            KeywordDatatypeMut::TBYTE(_) => TBYTE,
            KeywordDatatypeMut::TSBYTE(_) => TSBYTE,
            KeywordDatatypeMut::TSHORT(_) => TSHORT,
            KeywordDatatypeMut::TUSHORT(_) => TUSHORT,
            KeywordDatatypeMut::TINT(_) => TINT,
            KeywordDatatypeMut::TUINT(_) => TUINT,
            KeywordDatatypeMut::TLONG(_) => TLONG,
            KeywordDatatypeMut::TULONG(_) => TULONG,
            KeywordDatatypeMut::TFLOAT(_) => TFLOAT,
            KeywordDatatypeMut::TDOUBLE(_) => TDOUBLE,
            KeywordDatatypeMut::TSTRING(_) => TSTRING,
            KeywordDatatypeMut::TLOGICAL(_) => TLOGICAL,
            KeywordDatatypeMut::TCOMPLEX(_) => TCOMPLEX,
            KeywordDatatypeMut::TDBLCOMPLEX(_) => TDBLCOMPLEX,
            KeywordDatatypeMut::TULONGLONG(_) => TULONGLONG,
            KeywordDatatypeMut::TLONGLONG(_) => TLONGLONG,
            KeywordDatatypeMut::INVALID(x) => *x,
        }
    }
}

pub(crate) fn bytes_per_datatype(datatype: c_int) -> Option<usize> {
    match datatype {
        TBIT => Some(1),
        TBYTE => Some(1),
        TLOGICAL => Some(1),
        TSBYTE => Some(1),
        TUSHORT => Some(2),
        TSHORT => Some(2),
        TUINT => Some(4),
        TINT => Some(4),
        TULONG => Some(std::mem::size_of::<c_ulong>()),
        TLONG => Some(std::mem::size_of::<c_long>()),
        TULONGLONG => Some(8),
        TLONGLONG => Some(8),
        TFLOAT => Some(4),
        TDOUBLE => Some(8),
        TCOMPLEX => Some(8),
        TDBLCOMPLEX => Some(16),
        TSTRING => Some(std::mem::size_of::<usize>()), // pointer size
        _ => None,
    }
}

//https://stackoverflow.com/questions/65601579/parse-an-integer-ignoring-any-non-numeric-suffix
fn atoi<F: FromStr>(input: &str) -> Result<F, <F as FromStr>::Err> {
    let input = input.trim();
    let i = input
        .find(|c: char| !c.is_numeric() && c != '-' && c != '+')
        .unwrap_or(input.len());
    input[..i].parse::<F>()
}

// https://stackoverflow.com/questions/65264069/alignment-of-floating-point-numbers-printed-in-scientific-notation
fn fmt_f64(num: f64, precision: usize, exp_pad: usize) -> String {
    let mut num = format!("{num:.precision$E}");
    // Safe to `unwrap` as `num` is guaranteed to contain `'e'`
    let exp = num.split_off(num.find('E').unwrap());

    let (sign, exp) = if exp.starts_with("E-") {
        ('-', &exp[2..])
    } else {
        ('+', &exp[1..])
    };
    num.push_str(&format!("E{sign}{exp:0>exp_pad$}"));

    num.to_string()
}

const WB_MODE: *const c_char = c"wb".as_ptr().cast::<c_char>();
const RB_MODE: *const c_char = c"rb".as_ptr().cast::<c_char>();

#[cfg(not(target_os = "windows"))]
unsafe extern "C" {

    pub unsafe static mut stdin: *mut FILE;

    pub unsafe static mut stdout: *mut FILE;

    pub unsafe static mut stderr: *mut FILE;
}

#[cfg(windows)]
unsafe extern "C" {
    pub unsafe fn __acrt_iob_func(idx: libc::c_uint) -> *mut FILE;
}

#[macro_export]
#[cfg(not(target_os = "windows"))]
macro_rules! STDIN {
    () => {
        $crate::stdin
    };
}

#[macro_export]
#[cfg(windows)]
macro_rules! STDIN {
    () => {
        $crate::__acrt_iob_func(0)
    };
}

#[macro_export]
#[cfg(not(target_os = "windows"))]
macro_rules! STDOUT {
    () => {
        $crate::stdout
    };
}

#[macro_export]
#[cfg(windows)]
macro_rules! STDOUT {
    () => {
        $crate::__acrt_iob_func(1)
    };
}

#[macro_export]
#[cfg(not(target_os = "windows"))]
macro_rules! STDERR {
    () => {
        $crate::stderr
    };
}

#[macro_export]
#[cfg(windows)]
macro_rules! STDERR {
    () => {
        $crate::__acrt_iob_func(2)
    };
}

#[cfg(test)]
mod tests {
    use std::{ffi::CString, slice};

    use crate::{aliases::ffclos_safer, cs};
    use bytemuck::cast_slice;
    use cfileio::ffinit_safer;

    use putkey::ffcrim_safer;
    use tempfile::Builder;

    use crate::{
        fitsio::{USHORT_IMG, fitsfile},
        testhelpers::with_temp_file,
    };

    use super::*;

    use crate::aliases::safer::{fits_update_key, fits_write_img};

    #[test]
    fn test_write_image() {
        unsafe {
            with_temp_file(|filename| {
                let bitpix = USHORT_IMG;
                let naxis = 2;
                const NAXES: [c_long; 2] = [300, 200];
                let mut storage: [[u16; NAXES[0] as usize]; NAXES[1] as usize] =
                    [[0; NAXES[0] as usize]; NAXES[1] as usize];
                let mut fptr: Option<Box<fitsfile>> = None;
                let mut status = 0;

                let tempfile = Builder::new()
                    .prefix("my-temporary-note")
                    .suffix(".fits")
                    .tempfile()
                    .unwrap();

                let tdir = Builder::new().prefix("rsfitsio-").tempdir().unwrap();
                let abc = Builder::new().prefix("prefix").tempfile();
                let tdir_path = tdir.path();
                let filename = tdir_path.join("test.fits");

                let filename_str = filename.to_str().expect("cannot create string filename");
                let filename_cstr = CString::new(filename_str).unwrap();

                status = ffinit_safer(
                    &mut fptr,
                    cast_slice(filename_cstr.as_bytes_with_nul()),
                    &mut status,
                );
                assert_eq!(status, 0);

                let mut fptr = fptr.unwrap();

                status = ffcrim_safer(&mut fptr, bitpix, naxis, &NAXES, &mut status);
                assert_eq!(status, 0);

                for jj in 0..(NAXES[1] as usize) {
                    for ii in 0..(NAXES[0] as usize) {
                        storage[jj][ii] = (ii + jj) as u16;
                    }
                }

                let fpixel = 1; /* first pixel to write      */
                let nelements = NAXES[0] * NAXES[1]; /* number of pixels to write */

                /* write the array of unsigned integers to the FITS file */
                let s = slice::from_raw_parts(
                    storage.as_ptr() as *mut u16,
                    (NAXES[0] * NAXES[1]) as usize,
                );
                fits_write_img(
                    &mut fptr,
                    TUSHORT,
                    fpixel,
                    nelements as LONGLONG,
                    cast_slice(s),
                    &mut status,
                );
                assert_eq!(status, 0);

                /* write another optional keyword to the header */
                /* Note that the ADDRESS of the value is passed in the routine */
                let exposure = 1500;
                fits_update_key(
                    &mut fptr,
                    KeywordDatatype::TLONG(&exposure),
                    cs!(c"EXPOSURE"),
                    Some(cs!(c"Total Exposure Time")),
                    &mut status,
                );
                assert_eq!(status, 0);

                ffclos_safer(fptr, &mut status); /* close the file */
                assert_eq!(status, 0);
            });
        }
    }
}
