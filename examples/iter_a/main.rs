#![allow(deprecated)]

use core::slice;
use std::ffi::CStr;
use std::ffi::CString;
use std::path::PathBuf;
use std::process::ExitCode;
use std::process::exit;
use std::ptr;

use libc::{c_float, c_int, c_long};
use rsfitsio::STDERR;
use rsfitsio::aliases::c_api::*;
use rsfitsio::fitsio::{
    BINARY_TBL, DOUBLENULLVALUE, INPUT_COL, OUTPUT_COL, READWRITE, TFLOAT, TLONG, fitsfile,
    iteratorCol,
};
use rsfitsio::putcol::{fits_iter_get_array, fits_iter_get_datatype, fits_iter_set_by_name};

/*
  This program illustrates how to use the CFITSIO iterator function.
  It reads and modifies the input 'iter_a.fit' file by computing a
  value for the 'rate' column as a function of the values in the other
  'counts' and 'time' columns.
*/
pub fn main() -> ExitCode {
    let mut fptr: Option<Box<fitsfile>> = None;
    let mut cols: [iteratorCol; 3] = unsafe { std::mem::zeroed() }; /* structure used by the iterator function */
    let n_cols: c_int = 3; /* number of columns */
    let rows_per_loop: c_long = 0; /* use default optimum number of rows */
    let offset: c_long = 0; /* process all the rows */

    let mut status: c_int = 0;

    let mut d = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    d.push("examples/iter_a/iter_a.fit");
    let filename = d.to_str().unwrap(); /* name of rate FITS file */

    /* open file */
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

    /* move to the desired binary table extension */
    if let Some(ref mut fptr_box) = fptr {
        let extname = c"RATE";
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
        let counts_name = c"COUNTS";
        let time_name = c"TIME";
        let rate_name = c"RATE";

        unsafe {
            fits_iter_set_by_name(
                &mut cols[0],
                fptr_box.as_mut(),
                counts_name.as_ptr(),
                TLONG,
                INPUT_COL as c_int,
            );
            fits_iter_set_by_name(
                &mut cols[1],
                fptr_box.as_mut(),
                time_name.as_ptr(),
                TFLOAT,
                INPUT_COL as c_int,
            );
            fits_iter_set_by_name(
                &mut cols[2],
                fptr_box.as_mut(),
                rate_name.as_ptr(),
                TFLOAT,
                OUTPUT_COL as c_int,
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
/// Sample iterator function that calculates the output flux 'rate' column
/// by dividing the input 'counts' by the 'time' column.
/// It also applies a constant deadtime correction factor if the 'deadtime'
/// keyword exists.  Finally, this creates or updates the 'LIVETIME'
/// keyword with the sum of all the individual integration times.
extern "C" fn flux_rate(
    totalrows: c_long,
    _offset: c_long,
    firstrow: c_long,
    nrows: c_long,
    ncols: c_int,
    cols: *mut iteratorCol,
    _user_strct: *mut std::os::raw::c_void,
) -> c_int {
    let mut status: c_int = 0;

    // Static variables to preserve values between calls - using unsafe static
    static mut COUNTS: *mut c_long = ptr::null_mut();
    static mut INTERVAL: *mut c_float = ptr::null_mut();
    static mut RATE: *mut c_float = ptr::null_mut();
    static mut DEADTIME: c_float = 1.0;
    static mut LIVETIME: c_float = 0.0;

    unsafe {
        let cols = slice::from_raw_parts_mut(cols, ncols as usize);

        /*--------------------------------------------------------*/
        /*  Initialization procedures: execute on the first call  */
        /*--------------------------------------------------------*/
        if firstrow == 1 {
            if ncols != 3 {
                return -1; /* number of columns incorrect */
            }

            if fits_iter_get_datatype(&mut cols[0]) != TLONG
                || fits_iter_get_datatype(&mut cols[1]) != TFLOAT
                || fits_iter_get_datatype(&mut cols[2]) != TFLOAT
            {
                return -2; /* bad data type */
            }

            /* assign the input pointers to the appropriate arrays */
            COUNTS = fits_iter_get_array(&mut cols[0]) as *mut c_long;
            INTERVAL = fits_iter_get_array(&mut cols[1]) as *mut c_float;
            RATE = fits_iter_get_array(&mut cols[2]) as *mut c_float;

            LIVETIME = 0.0; /* initialize the total integration time */

            /* try to get the deadtime keyword value */
            let deadtime_key = c"DEADTIME";
            let mut null_comment = [0u8; 1];

            let read_status = fits_read_key(
                (cols[0]).fptr,
                TFLOAT,
                deadtime_key.as_ptr(),
                &raw mut DEADTIME as *mut std::os::raw::c_void,
                null_comment.as_mut_ptr() as *mut libc::c_char,
                &mut status,
            );

            if read_status != 0 {
                DEADTIME = 1.0; /* default deadtime if keyword doesn't exist */
                status = 0; /* reset status */
            } else if DEADTIME < 0.0 || DEADTIME > 1.0 {
                return -1; /* bad deadtime value */
            }

            println!("deadtime = {:.6}", ptr::read(&raw const DEADTIME));
        }

        /*--------------------------------------------*/
        /*  Main loop: process all the rows of data */
        /*--------------------------------------------*/

        /*  NOTE: 1st element of array is the null pixel value!  */
        /*  Loop from 1 to nrows, not 0 to nrows - 1.  */

        /* this version tests for null values */
        *RATE.offset(0) = DOUBLENULLVALUE as c_float; /* define the value that represents null */

        for ii in 1..=nrows {
            let counts_val = *COUNTS.offset(ii as isize);
            let counts_null = *COUNTS.offset(0);
            let interval_val = *INTERVAL.offset(ii as isize);

            if counts_val == counts_null {
                /*  undefined counts value? */
                *RATE.offset(ii as isize) = DOUBLENULLVALUE as c_float;
            } else if interval_val > 0.0 {
                *RATE.offset(ii as isize) = (counts_val as c_float) / interval_val / DEADTIME;
                LIVETIME += interval_val; /* accumulate total integration time */
            } else {
                return -2; /* bad integration time */
            }
        }

        /*-------------------------------------------------------*/
        /*  Clean up procedures:  after processing all the rows  */
        /*-------------------------------------------------------*/

        if firstrow + nrows - 1 == totalrows {
            /*  update the LIVETIME keyword value */
            let livetime_key = c"LIVETIME";
            let livetime_comment = c"total integration time";

            fits_update_key(
                (cols[0]).fptr,
                TFLOAT,
                livetime_key.as_ptr(),
                &raw const LIVETIME as *const std::os::raw::c_void,
                livetime_comment.as_ptr(),
                &mut status,
            );
            println!("livetime = {:.6}", ptr::read(&raw const LIVETIME));
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
