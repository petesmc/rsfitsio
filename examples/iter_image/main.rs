#![allow(deprecated)]

use std::process::ExitCode;
use std::process::exit;
use std::ptr;

use libc::{c_int, c_long};
use rsfitsio::STDERR;
use rsfitsio::aliases::c_api::*;
use rsfitsio::fitsio::{INPUT_OUTPUT_COL, READWRITE, fitsfile, iteratorCol};
use rsfitsio::putcol::{
    fits_iter_get_array, fits_iter_set_datatype, fits_iter_set_file, fits_iter_set_iotype,
};

/*
  This program illustrates how to use the CFITSIO iterator function.
  It reads and modifies the input 'iter_image.fit' image file by setting
  all the pixel values to zero (DESTROYING THE ORIGINAL IMAGE!!!)
*/
pub fn main() -> ExitCode {
    let mut fptr: Option<Box<fitsfile>> = None;
    let mut cols: [iteratorCol; 3] = unsafe { std::mem::zeroed() }; /* structure used by the iterator function */
    let n_cols: c_int = 1;
    let rows_per_loop: c_long = 0; /* use default optimum number of rows */
    let offset: c_long = 0; /* process all the rows */

    let mut status: c_int = 0;
    let filename = c"iter_image.fit"; /* name of rate FITS file */

    /* open file */
    if unsafe { fits_open_file(&mut fptr, filename.as_ptr(), READWRITE, &mut status) } != 0 {
        printerror(status);
    }

    /* define input column structure members for the iterator function */
    if let Some(ref mut fptr_box) = fptr {
        unsafe {
            fits_iter_set_file(&mut cols[0], fptr_box.as_mut());
            fits_iter_set_iotype(&mut cols[0], INPUT_OUTPUT_COL as c_int);
            fits_iter_set_datatype(&mut cols[0], 0);
        }
    }

    /* apply the rate function to each row of the table */
    println!("Calling iterator function...{status}");

    unsafe {
        fits_iterate_data(
            n_cols,
            cols.as_mut_ptr(),
            offset,
            rows_per_loop,
            zero_image,
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
/// Sample iterator function that zeros out image pixels.
/// This demonstrates how to process FITS image data using the iterator.
extern "C" fn zero_image(
    _totalrows: c_long,
    _offset: c_long,
    firstrow: c_long,
    nrows: c_long,
    ncols: c_int,
    cols: *mut iteratorCol,
    _user_strct: *mut std::os::raw::c_void,
) -> c_int {
    // Static variables to preserve values between calls - using unsafe static
    static mut COUNTS: *mut c_int = ptr::null_mut();

    unsafe {
        /*--------------------------------------------------------*/
        /*  Initialization procedures: execute on the first call  */
        /*--------------------------------------------------------*/
        if firstrow == 1 {
            if ncols != 1 {
                return -1; /* number of columns incorrect */
            }

            /* assign the input pointers to the appropriate arrays */
            COUNTS = fits_iter_get_array(cols.offset(0)) as *mut c_int;
        }

        /*--------------------------------------------*/
        /*  Main loop: process all the rows of data */
        /*--------------------------------------------*/

        /*  NOTE: 1st element of array is the null pixel value!  */
        /*  Loop from 1 to nrows, not 0 to nrows - 1.  */

        for ii in 1..=nrows {
            *COUNTS.offset(ii as isize) = 1;
        }

        println!("firstrows, nrows = {firstrow} {nrows}");
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
