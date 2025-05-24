#![allow(deprecated)]

use std::{ffi::CStr, process::ExitCode};

use bytemuck::cast_slice;
use rsfitsio::c_types::{c_char, c_ulong};
use rsfitsio::{
    aliases::{ffcmsg, fits_open_file},
    cfileio::ffclos,
    fitsio::{READWRITE, fitsfile},
    wrappers::strcpy_safe,
};

fn strcpy_cstr(dest: &mut [c_char], src: &CStr) {
    strcpy_safe(dest, cast_slice(src.to_bytes_with_nul()));
}

pub fn main() -> ExitCode {
    unsafe {
        let mut fptr: Option<Box<fitsfile>> = None;
        let mut tblname: [c_char; 40] = [0; 40];

        let mut status = 0;
        strcpy_cstr(&mut tblname, c"Test-ASCII");

        println!("CFITSIO TESTPROG\n");
        println!("Try opening then closing a nonexistent file:");
        fits_open_file(&mut fptr, c"tq123x.kjl".as_ptr(), READWRITE, &mut status);
        println!(
            "  ffopen fptr, status  = {} {} (expect an error)",
            fptr.as_deref_mut().unwrap() as *mut _ as c_ulong,
            status
        );
        ffclos(fptr.unwrap().as_mut(), &mut status);
        println!("  ffclos status = {status}\n");
        ffcmsg();
        status = 0;
    }

    ExitCode::from(0)
}
