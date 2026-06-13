use std::ffi::CStr;
use std::{f32, f64};

use bytemuck::cast_slice;
use tempfile::Builder;

use crate::aliases::rust_api::{fits_read_key_str, fits_write_imghdr};
use crate::c_types::{c_char, c_int};
use crate::cfileio::ffinit_safe;
use crate::cs;
use crate::fitsio::{BYTE_IMG, FLEN_FILENAME, FLEN_VALUE, fitsfile};

/// Function to allow access to a temporary file
pub fn with_temp_file<F>(callback: F)
where
    F: for<'a> Fn(&'a str),
{
    let tdir = Builder::new().prefix("rsfitsio-").tempdir().unwrap();
    let tdir_path = tdir.path();
    let filename = tdir_path.join("test.fits");

    let filename_str = filename.to_str().expect("cannot create string filename");
    callback(filename_str);
}

/// Copy a `&str` into a NUL-terminated `[c_char; FLEN_FILENAME]` buffer.
pub fn to_buf(s: &str) -> [c_char; FLEN_FILENAME] {
    let mut buf = [0 as c_char; FLEN_FILENAME];
    for (i, &b) in s.as_bytes().iter().enumerate() {
        buf[i] = b as c_char;
    }
    buf
}

/// Read a NUL-terminated `[c_char]` buffer back into a `&str`.
pub fn from_buf(buf: &[c_char]) -> &str {
    CStr::from_bytes_until_nul(cast_slice(buf))
        .unwrap()
        .to_str()
        .unwrap()
}

/// Create an in-memory FITS file with an empty `BYTE_IMG` primary HDU.
pub fn new_mem_file(status: &mut c_int) -> Option<Box<fitsfile>> {
    let mut fptr: Option<Box<fitsfile>> = None;
    ffinit_safe(&mut fptr, cs!(c"mem://"), status);
    assert_eq!(*status, 0, "ffinit failed");
    let f = fptr.as_deref_mut().unwrap();
    fits_write_imghdr(f, BYTE_IMG, 0, &[], status);
    assert_eq!(*status, 0, "ffphps failed");
    fptr
}

/// Read a string-valued keyword and return it as an owned `String`.
pub fn read_str_key(f: &mut fitsfile, key: &[c_char], status: &mut c_int) -> String {
    let mut buf = [0 as c_char; FLEN_VALUE];
    fits_read_key_str(f, key, &mut buf, None, status);
    from_buf(&buf).to_string()
}

/// Build a `path` + `ext` (e.g. "[GROUPING]") NUL-terminated filename buffer.
pub fn path_with_ext(path: &str, ext: &str) -> [c_char; FLEN_FILENAME] {
    let mut s = String::from(path);
    s.push_str(ext);
    to_buf(&s)
}

/// Helper function for float comparisons
pub fn floats_close_f32(a: f32, b: f32) -> bool {
    (a - b).abs() < f32::EPSILON
}

pub fn floats_close_f64(a: f64, b: f64) -> bool {
    (a - b).abs() < f64::EPSILON
}
