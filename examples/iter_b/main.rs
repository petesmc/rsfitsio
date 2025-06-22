#![allow(deprecated)]

use core::slice;
use std::ffi::CStr;
use std::ffi::CString;
use std::path::PathBuf;
use std::process::ExitCode;
use std::process::exit;
use std::ptr;

use libc::{c_char, c_int, c_long, strcpy};
use rsfitsio::STDERR;
use rsfitsio::aliases::c_api::*;
use rsfitsio::fitsio::{
    BINARY_TBL, FALSE, INPUT_OUTPUT_COL, READWRITE, TLOGICAL, TRUE, TSTRING, fitsfile, iteratorCol,
};
use rsfitsio::putcol::{fits_iter_get_array, fits_iter_get_datatype, fits_iter_set_by_name};

/*
  This program illustrates how to use the CFITSIO iterator function.
  It simply prints out the values in a character string and a logical
  type column in a table, and toggles the value in the logical column
  so that T -> F and F -> T.
*/
pub fn main() -> ExitCode {
    let mut fptr: Option<Box<fitsfile>> = None;
    let mut cols: [iteratorCol; 2] = unsafe { std::mem::zeroed() };
    let n_cols: c_int = 2; /* number of columns */
    let rows_per_loop: c_long = 0; /* use default optimum number of rows */
    let offset: c_long = 0; /* process all the rows */
    let mut status: c_int = 0;

    let mut d = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    d.push("examples/iter_b/iter_b.fit");
    let filename = d.to_str().unwrap(); /* name of rate FITS file */

    /* open the file and move to the correct extension */
    if unsafe {
        fits_open_file(
            &mut fptr,
            CString::new(filename).unwrap().as_ptr(),
            READWRITE,
            &mut status,
        )
    } != 0
    {
        printerror(status);
    }

    if let Some(ref mut fptr_box) = fptr {
        let extname = c"iter_test";
        if unsafe {
            fits_movnam_hdu(
                fptr_box.as_mut(),
                BINARY_TBL,
                extname.as_ptr(),
                0,
                &mut status,
            )
        } != 0
        {
            printerror(status);
        }
    }

    /* define input column structure members for the iterator function */
    if let Some(ref mut fptr_box) = fptr {
        let avalue_name = c"Avalue";
        let lvalue_name = c"Lvalue";

        unsafe {
            fits_iter_set_by_name(
                &mut cols[0],
                fptr_box.as_mut(),
                avalue_name.as_ptr(),
                TSTRING,
                INPUT_OUTPUT_COL as c_int,
            );
            fits_iter_set_by_name(
                &mut cols[1],
                fptr_box.as_mut(),
                lvalue_name.as_ptr(),
                TLOGICAL,
                INPUT_OUTPUT_COL as c_int,
            );
        }
    }

    /* apply the function to each row of the table */
    println!("Calling iterator function...{status}");

    unsafe {
        fits_iterate_data(
            n_cols,
            cols.as_mut_ptr(),
            offset,
            rows_per_loop,
            str_iter,
            ptr::null_mut(),
            &mut status,
        );
    }

    /* all done */
    if let Some(fptr_box) = fptr {
        if unsafe { fits_close_file(Some(fptr_box), &mut status) } != 0 {
            printerror(status);
        }
    }

    if status != 0 {
        printerror(status);
    }

    ExitCode::from(status as u8)
}
/*--------------------------------------------------------------------------*/
/// Sample iterator function that demonstrates string and logical column processing.
extern "C" fn str_iter(
    totalrows: c_long,
    _offset: c_long,
    firstrow: c_long,
    nrows: c_long,
    ncols: c_int,
    cols: *mut iteratorCol,
    _user_strct: *mut std::os::raw::c_void,
) -> c_int {
    // Static variables to preserve values between calls - using unsafe static
    static mut STRINGVALS: *mut *mut c_char = ptr::null_mut();
    static mut LOGICALVALS: *mut c_char = ptr::null_mut();

    unsafe {
        let cols = slice::from_raw_parts_mut(cols, ncols as usize);

        /*--------------------------------------------------------*/
        /*  Initialization procedures: execute on the first call  */
        /*--------------------------------------------------------*/
        if firstrow == 1 {
            if ncols != 2 {
                return -1; /* number of columns incorrect */
            }

            if fits_iter_get_datatype(&mut cols[0]) != TSTRING
                || fits_iter_get_datatype(&mut cols[1]) != TLOGICAL
            {
                return -2; /* bad data type */
            }

            /* assign the input pointers to the appropriate arrays */
            STRINGVALS = fits_iter_get_array(&mut cols[0]) as *mut *mut c_char;
            LOGICALVALS = fits_iter_get_array(&mut cols[1]) as *mut c_char;

            println!("Total rows, No. rows = {totalrows} {nrows}");
        }

        /*------------------------------------------*/
        /*  Main loop: process all the rows of data */
        /*------------------------------------------*/

        /*  NOTE: 1st element of array is the null pixel value!  */
        /*  Loop from 1 to nrows, not 0 to nrows - 1.  */

        let stringvals_slice = slice::from_raw_parts_mut(STRINGVALS, nrows as usize + 1);
        let logicalvals_slice = slice::from_raw_parts_mut(LOGICALVALS, nrows as usize + 1);

        for ii in 1..=nrows {
            let string_ptr = stringvals_slice[ii as usize];
            let logical_val = logicalvals_slice[ii as usize] as c_int;

            // Convert the C string to a Rust string for printing
            let string_val = if !string_ptr.is_null() {
                CStr::from_ptr(string_ptr).to_string_lossy()
            } else {
                "null".into()
            };

            println!("{string_val} {logical_val}");

            if logical_val != 0 {
                logicalvals_slice[ii as usize] = FALSE as c_char;
                let new_string = c"changed to false";
                strcpy(string_ptr, new_string.as_ptr());
            } else {
                logicalvals_slice[ii as usize] = TRUE as c_char;
                let new_string = c"changed to true";
                strcpy(string_ptr, new_string.as_ptr());
            }
        }

        /*-------------------------------------------------------*/
        /*  Clean up procedures:  after processing all the rows  */
        /*-------------------------------------------------------*/

        if firstrow + nrows - 1 == totalrows { /* no action required in this case */ }
    }

    0 /* return successful status */
}

/*--------------------------------------------------------------------------*/
/// Print out cfitsio error messages and exit program
fn printerror(status: c_int) {
    if status != 0 {
        unsafe {
            fits_report_error(STDERR!(), status); /* print error report */
        }
        exit(status); /* terminate the program, returning error status */
    }
}
