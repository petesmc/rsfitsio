#![allow(deprecated)]

use core::slice;
use std::ffi::CString;
use std::fs::remove_file;
use std::path::PathBuf;
use std::process::ExitCode;
use std::process::exit;
use std::ptr;

use libc::{c_int, c_long};
use rsfitsio::STDERR;
use rsfitsio::aliases::c_api::*;
use rsfitsio::fitsio::{
    BINARY_TBL, INPUT_COL, LONG_IMG, OUTPUT_COL, READONLY, TLONG, fitsfile, iteratorCol,
};
use rsfitsio::putcol::{fits_iter_get_array, fits_iter_get_datatype, fits_iter_set_by_name};

/*
    This example program illustrates how to use the CFITSIO iterator function.

    This program creates a 2D histogram of the X and Y columns of an event
    list.  The 'main' routine just creates the empty new image, then executes
    the 'writehisto' work function by calling the CFITSIO iterator function.

    'writehisto' opens the FITS event list that contains the X and Y columns.
    It then calls a second work function, calchisto, (by recursively calling
    the CFITSIO iterator function) which actually computes the 2D histogram.
*/

/*   Globally defined parameters */
static XSIZE: c_long = 480; /* size of the histogram image */
static YSIZE: c_long = 480;
static XBINSIZE: c_long = 32;
static YBINSIZE: c_long = 32;

pub fn main() -> ExitCode {
    let mut fptr: Option<Box<fitsfile>> = None;
    let mut cols: [iteratorCol; 1] = Default::default();
    let n_cols: c_int = 1; /* number of columns */
    let mut status: c_int = 0;
    let n_per_loop: c_long = -1; /* force whole array to be passed at one time */
    let offset: c_long = 0; /* don't skip over any pixels */
    let naxes: [c_long; 2] = [XSIZE, YSIZE];
    let filename = c"histoimg.fit"; /* name of FITS image */

    let _ = remove_file(filename.to_str().unwrap()); /* delete previous version of the file if it exists */

    /* create new output image */
    if unsafe { fits_create_file(&mut fptr, filename.as_ptr(), &mut status) } != 0 {
        printerror(status);
    }

    /* create primary HDU */
    if let Some(ref mut fptr_box) = fptr
        && unsafe {
            fits_create_img(
                fptr_box.as_mut(),
                LONG_IMG,
                2,
                naxes.as_ptr().cast_mut(),
                &mut status,
            )
        } != 0
    {
        printerror(status);
    }

    /* define input column structure members for the iterator function */
    if let Some(ref mut fptr_box) = fptr {
        let empty_name = c" ";
        unsafe {
            fits_iter_set_by_name(
                &mut cols[0],
                fptr_box.as_mut(),
                empty_name.as_ptr(),
                TLONG,
                OUTPUT_COL as c_int,
            );
        }
    }

    /* execute the function to create and write the 2D histogram */
    println!("Calling writehisto iterator work function... {status}");

    unsafe {
        fits_iterate_data(
            n_cols,
            cols.as_mut_ptr(),
            offset,
            n_per_loop,
            writehisto,
            ptr::null_mut(),
            &mut status,
        );
    }

    /* all done; close the file */
    if let Some(fptr_box) = fptr
        && unsafe { fits_close_file(Some(fptr_box), &mut status) } != 0
    {
        printerror(status);
    }

    if status != 0 {
        printerror(status);
    } else {
        println!("Program completed successfully.");
    }

    ExitCode::from(status as u8)
}

/*--------------------------------------------------------------------------*/
/// Iterator work function that writes out the 2D histogram.
/// The histogram values are calculated by another work function, calchisto.
/// This routine is executed only once since nvalues was forced to = totaln.
extern "C" fn writehisto(
    totaln: c_long,
    _offset: c_long,
    _firstn: c_long,
    nvalues: c_long,
    narrays: c_int,
    histo: *mut iteratorCol,
    _user_pointer: *mut std::os::raw::c_void,
) -> c_int {
    let mut tblptr: Option<Box<fitsfile>> = None;
    let mut cols: [iteratorCol; 2] = Default::default();
    let n_cols: c_int = 2; /* number of columns */
    let mut status: c_int = 0;
    let rows_per_loop: c_long = 0; /* take default number of rows per iteration */
    let rowoffset: c_long = 0;

    let mut d = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    d.push("examples/iter_c/iter_c.fit");
    let filename = d.to_str().unwrap(); /* name of FITS table */

    unsafe {
        /* do sanity checking of input values */
        if totaln != nvalues {
            return -1; /* whole image must be passed at one time */
        }

        if narrays != 1 {
            return -2; /* number of images is incorrect */
        }

        if fits_iter_get_datatype(histo.offset(0)) != TLONG {
            return -3; /* input array has wrong data type */
        }

        /* assign the FITS array pointer to the global histogram pointer */
        let histogram = fits_iter_get_array(histo.offset(0)) as *mut c_long;

        /* open the file and move to the table containing the X and Y columns */

        if fits_open_file(
            &mut tblptr,
            CString::new(filename).unwrap().as_ptr(),
            READONLY,
            &mut status,
        ) != 0
        {
            return status;
        }

        if let Some(ref mut tblptr_box) = tblptr {
            let extname = c"EVENTS";
            if fits_movnam_hdu(
                tblptr_box.as_mut(),
                BINARY_TBL,
                extname.as_ptr(),
                0,
                &mut status,
            ) != 0
            {
                return status;
            }
        }

        /* define input column structure members for the iterator function */
        if let Some(ref mut tblptr_box) = tblptr {
            let x_name = c"X";
            let y_name = c"Y";

            fits_iter_set_by_name(
                &mut cols[0],
                tblptr_box.as_mut(),
                x_name.as_ptr(),
                TLONG,
                INPUT_COL as c_int,
            );
            fits_iter_set_by_name(
                &mut cols[1],
                tblptr_box.as_mut(),
                y_name.as_ptr(),
                TLONG,
                INPUT_COL as c_int,
            );
        }

        /* calculate the histogram */
        println!("Calling calchisto iterator work function... {status}");

        fits_iterate_data(
            n_cols,
            cols.as_mut_ptr(),
            rowoffset,
            rows_per_loop,
            calchisto,
            histogram as *mut std::os::raw::c_void,
            &mut status,
        );

        /* all done */
        if let Some(tblptr_box) = tblptr {
            let _ = fits_close_file(Some(tblptr_box), &mut status);
        }
    }

    status
}

/*--------------------------------------------------------------------------*/
/// Iterator work function that calculates values for the 2D histogram.
extern "C" fn calchisto(
    _totalrows: c_long,
    _offset: c_long,
    firstrow: c_long,
    nrows: c_long,
    ncols: c_int,
    cols: *mut iteratorCol,
    user_pointer: *mut std::os::raw::c_void,
) -> c_int {
    // Static variables to preserve values between calls - using unsafe static
    static mut XCOL: *mut c_long = ptr::null_mut();
    static mut YCOL: *mut c_long = ptr::null_mut();
    static mut HISTOGRAM: *mut c_long = ptr::null_mut();

    unsafe {
        let cols = slice::from_raw_parts_mut(cols, ncols as usize);

        /*--------------------------------------------------------*/
        /*  Initialization procedures: execute on the first call  */
        /*--------------------------------------------------------*/
        if firstrow == 1 {
            /* do sanity checking of input values */
            if ncols != 2 {
                return -3; /* number of arrays is incorrect */
            }

            if fits_iter_get_datatype(&mut cols[0]) != TLONG
                || fits_iter_get_datatype(&mut cols[1]) != TLONG
            {
                return -4; /* wrong datatypes */
            }

            /* assign the input array pointers to the X and Y arrays */
            XCOL = fits_iter_get_array(&mut cols[0]) as *mut c_long;
            YCOL = fits_iter_get_array(&mut cols[1]) as *mut c_long;
            HISTOGRAM = user_pointer as *mut c_long;

            /* initialize the histogram image pixels = 0 */
            for ii in 0..=(XSIZE * YSIZE) {
                *HISTOGRAM.offset(ii as isize) = 0;
            }
        }

        /*------------------------------------------------------------------*/
        /*  Main loop: increment the 2D histogram at position of each event */
        /*------------------------------------------------------------------*/

        for ii in 1..=nrows {
            let xbin = *XCOL.offset(ii as isize) / XBINSIZE;
            let ybin = *YCOL.offset(ii as isize) / YBINSIZE;

            let ihisto = (ybin * XSIZE) + xbin + 1;
            if ihisto >= 0 && ihisto < (XSIZE * YSIZE) {
                *HISTOGRAM.offset(ihisto as isize) += 1;
            }
        }
    }

    0
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
