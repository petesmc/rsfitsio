//! The bzip2 entry points used by the `BZIP2_1` tile compressor.
//!
//! Re-exported from `libbz2-rs-sys` behind the default `bzip2` feature. BZIP2
//! is not publicly supported by CFITSIO and exists for test purposes only.
#![warn(missing_docs)]

use crate::c_types::{c_int, c_void, FILE};

pub use bzip2_sys::{
    BZ_CONFIG_ERROR, BZ_DATA_ERROR, BZ_IO_ERROR, BZ_MEM_ERROR, BZ_OK, BZ_STREAM_END,
    BZ_UNEXPECTED_EOF, 
};

pub type BZFILE = c_void;


unsafe extern {

    pub fn BZ2_bzReadOpen(
        bzerror: *mut c_int,
        f: *mut FILE,
        verbosity: c_int,
        small: c_int,
        unused: *mut c_void,
        nUnused: c_int,
    ) -> *mut BZFILE;

    pub fn BZ2_bzReadClose(bzerror: *mut c_int, b: *mut BZFILE);

    pub fn BZ2_bzRead(bzerror: *mut c_int, b: *mut BZFILE, buf: *mut c_void, len: c_int) -> c_int;
}
