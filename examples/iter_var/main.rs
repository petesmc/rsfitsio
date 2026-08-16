#![allow(deprecated)]

use std::process::ExitCode;
use std::process::exit;
use std::ptr;

use libc::{c_int, c_long};
use rsfitsio::STDERR;
use rsfitsio::aliases::c_api::*;
use rsfitsio::fitsio::{BINARY_TBL, INPUT_COL, READWRITE, fitsfile, iteratorCol};
use rsfitsio::putcol::{
    fits_iter_get_array, fits_iter_get_datatype, fits_iter_get_repeat, fits_iter_set_by_name,
};

/*
  This program illustrates how to use the CFITSIO iterator function.
  It reads the input 'vari.fits' file, moves to the binary
  table in the "COMPRESSED_IMAGE" extension, and prints the
  float-values in the column COMPRESSED_DATA.
*/
pub fn main() -> ExitCode {
    let mut fptr: Option<Box<fitsfile>> = None;
    let mut cols: [iteratorCol; 3] = Default::default(); /* structure used by the iterator function */
    let n_cols: c_int = 1; /* number of columns */
    let rows_per_loop: c_long = 0; /* use default optimum number of rows */
    let offset: c_long = 0; /* process all the rows */

    let mut status: c_int = 0;
    let filename = c"vari.fits"; /* name of rate FITS file */

    /* open file */
    if unsafe { fits_open_file(&raw mut fptr, filename.as_ptr(), READWRITE, &raw mut status) } != 0
    {
        printerror(status);
    }

    /* move to the desired binary table extension */
    if let Some(ref mut fptr_box) = fptr {
        let extname = c"COMPRESSED_IMAGE";
        if unsafe {
            fits_movnam_hdu(
                fptr_box.as_mut(),
                BINARY_TBL,
                extname.as_ptr(),
                0,
                &raw mut status,
            )
        } != 0
        {
            printerror(status);
        }
    }

    /* define input column structure members for the iterator function */
    if let Some(ref mut fptr_box) = fptr {
        let compressed_data_name = c"COMPRESSED_DATA";

        unsafe {
            fits_iter_set_by_name(
                &raw mut cols[0],
                fptr_box.as_mut(),
                compressed_data_name.as_ptr(),
                0,
                INPUT_COL as c_int,
            );
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
            flux_rate,
            ptr::null_mut(),
            &raw mut status,
        );
    }

    /* all done */
    if let Some(fptr_box) = fptr
        && unsafe { fits_close_file(Some(fptr_box), &raw mut status) } != 0
    {
        printerror(status);
    }

    if status != 0 {
        printerror(status);
    }

    ExitCode::from(status as u8)
}

/*--------------------------------------------------------------------------*/
/// Sample iterator function that prints the values (assumed to beof type float).
extern "C" fn flux_rate(
    _totalrows: c_long,
    _offset: c_long,
    firstrow: c_long,
    nrows: c_long,
    _ncols: c_int,
    cols: *mut iteratorCol,
    _user_strct: *mut std::os::raw::c_void,
) -> c_int {
    // Static variables to preserve values between calls - using unsafe static
    static mut COUNTS: *mut f32 = ptr::null_mut();

    unsafe {
        /*--------------------------------------------------------*/
        /*  Initialization procedures: execute on the first call  */
        /*--------------------------------------------------------*/
        if firstrow == 1 {
            println!(
                "Datatype of column = {}",
                fits_iter_get_datatype(cols.offset(0))
            );

            /* assign the input pointers to the appropriate arrays */
            COUNTS = fits_iter_get_array(cols.offset(0)) as *mut f32;
        }

        /*--------------------------------------------*/
        /*  Main loop: process all the rows of data */
        /*--------------------------------------------*/

        /*  NOTE: 1st element of array is the null pixel value!  */
        /*  Loop from 1 to nrows, not 0 to nrows - 1.  */

        for ii in 1..=nrows {
            let repeat = fits_iter_get_repeat(cols.offset(0));
            let first_element = if !COUNTS.is_null() {
                *COUNTS.add(ii as usize)
            } else {
                0.0
            };
            println!("repeat = {repeat}, {first_element}");
        }
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
