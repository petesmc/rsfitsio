#![allow(deprecated)]

use std::ffi::{CStr, CString};
use std::fs::remove_file;
use std::process::ExitCode;
use std::process::exit;
use std::ptr;

use libc::{c_char, c_float, c_int, c_long, c_ushort, free, malloc};
use rsfitsio::STDERR;
use rsfitsio::aliases::c_api::*;
use rsfitsio::fitsio::{
    ASCII_TBL, BINARY_TBL, CASEINSEN, END_OF_FILE, FLEN_CARD, FLEN_VALUE, LONGLONG, READONLY,
    READWRITE, TFLOAT, TLONG, TSTRING, TUSHORT, USHORT_IMG, fitsfile,
};

pub fn main() -> ExitCode {
    /*************************************************************************
       This is a simple main program that calls the following routines:

        writeimage    - write a FITS primary array image
        writeascii    - write a FITS ASCII table extension
        writebintable - write a FITS binary table extension
        copyhdu       - copy a header/data unit from one FITS file to another
        selectrows    - copy selected row from one HDU to another
        readheader    - read and print the header keywords in every extension
        readimage     - read a FITS image and compute the min and max value
        readtable     - read columns of data from ASCII and binary tables

    **************************************************************************/

    writeimage();
    writeascii();
    writebintable();
    copyhdu();
    selectrows();
    readheader();
    readimage();
    readtable();

    println!("\nAll the cfitsio cookbook routines ran successfully.");

    ExitCode::from(0)
}

/*--------------------------------------------------------------------------*/
/// Create a FITS primary array containing a 2-D image
fn writeimage() {
    let mut fptr: Option<Box<fitsfile>> = None;
    let mut status: c_int = 0;
    let fpixel: c_long = 1;

    /* initialize FITS image parameters */
    let filename = "atestfil.fit"; /* name for new FITS file */
    let bitpix: c_int = USHORT_IMG; /* 16-bit unsigned short pixel values       */
    let naxis: c_long = 2; /* 2-dimensional image                            */
    let naxes: [c_long; 2] = [300, 200]; /* image is 300 pixels wide by 200 rows */
    let nelements: c_long = naxes[0] * naxes[1];

    /* allocate memory for the whole image */
    let array_size = (naxes[0] * naxes[1]) as usize;
    let mut array: Vec<c_ushort> = vec![0; array_size];

    let _ = remove_file(filename); /* Delete old file if it already exists */

    let filename_cstr = CString::new(filename).unwrap();

    if unsafe { fits_create_file(&raw mut fptr, filename_cstr.as_ptr(), &raw mut status) } != 0 {
        printerror(status); /* call printerror if error occurs */
    }

    /* write the required keywords for the primary array image.     */
    /* Since bitpix = USHORT_IMG, this will cause cfitsio to create */
    /* a FITS image with BITPIX = 16 (signed short integers) with   */
    /* BSCALE = 1.0 and BZERO = 32768.  This is the convention that */
    /* FITS uses to store unsigned integers.  Note that the BSCALE  */
    /* and BZERO keywords will be automatically written by cfitsio  */
    /* in this case.                                                */

    if let Some(ref mut fptr_box) = fptr
        && unsafe {
            fits_create_img(
                fptr_box.as_mut(),
                bitpix,
                naxis as c_int,
                naxes.as_ptr().cast_mut(),
                &raw mut status,
            )
        } != 0
    {
        printerror(status);
    }

    /* initialize the values in the image with a linear ramp function */
    for jj in 0..naxes[1] {
        for ii in 0..naxes[0] {
            let index = (jj * naxes[0] + ii) as usize;
            array[index] = (ii + jj) as c_ushort;
        }
    }

    /* write the array of unsigned integers to the FITS file */
    if let Some(ref mut fptr_box) = fptr
        && unsafe {
            fits_write_img(
                fptr_box.as_mut(),
                TUSHORT,
                fpixel as LONGLONG,
                nelements as LONGLONG,
                array.as_ptr().cast(),
                &raw mut status,
            )
        } != 0
    {
        printerror(status);
    }

    /* write another optional keyword to the header */
    /* Note that the ADDRESS of the value is passed in the routine */
    let exposure: c_long = 1500;
    if let Some(ref mut fptr_box) = fptr {
        let exposure_key = CString::new("EXPOSURE").unwrap();
        let exposure_comment = CString::new("Total Exposure Time").unwrap();
        if unsafe {
            fits_update_key(
                fptr_box.as_mut(),
                TLONG,
                exposure_key.as_ptr(),
                core::ptr::from_ref::<c_long>(&exposure).cast(),
                exposure_comment.as_ptr(),
                &raw mut status,
            )
        } != 0
        {
            printerror(status);
        }
    }

    if let Some(fptr_box) = fptr
        && unsafe { fits_close_file(Some(fptr_box), &raw mut status) } != 0
    {
        /* close the file */
        printerror(status);
    }
}

/*--------------------------------------------------------------------------*/
/// Create an ASCII table extension containing 3 columns and 6 rows
fn writeascii() {
    let mut fptr: Option<Box<fitsfile>> = None;
    let mut status: c_int = 0;
    let firstrow: LONGLONG = 1;
    let firstelem: LONGLONG = 1;

    let tfields: c_int = 3; /* table will have 3 columns */
    let nrows: LONGLONG = 6; /* table will have 6 rows    */

    let filename = "atestfil.fit"; /* name for new FITS file */
    let extname = "PLANETS_ASCII"; /* extension name */

    /* define the name, datatype, and physical units for the 3 columns */
    let ttype_strs = [
        CString::new("Planet").unwrap(),
        CString::new("Diameter").unwrap(),
        CString::new("Density").unwrap(),
    ];
    let tform_strs = [
        CString::new("a8").unwrap(),
        CString::new("I6").unwrap(),
        CString::new("F4.2").unwrap(),
    ];
    let tunit_strs = [
        CString::new("").unwrap(),
        CString::new("km").unwrap(),
        CString::new("g/cm").unwrap(),
    ];

    let ttype: Vec<*const c_char> = ttype_strs.iter().map(|s| s.as_ptr()).collect();
    let tform: Vec<*const c_char> = tform_strs.iter().map(|s| s.as_ptr()).collect();
    let tunit: Vec<*const c_char> = tunit_strs.iter().map(|s| s.as_ptr()).collect();

    /* define the name diameter, and density of each planet */
    let planet_strs = [
        CString::new("Mercury").unwrap(),
        CString::new("Venus").unwrap(),
        CString::new("Earth").unwrap(),
        CString::new("Mars").unwrap(),
        CString::new("Jupiter").unwrap(),
        CString::new("Saturn").unwrap(),
    ];
    let planet: Vec<*const c_char> = planet_strs.iter().map(|s| s.as_ptr()).collect();
    let diameter: [c_long; 6] = [4880, 12112, 12742, 6800, 143000, 121000];
    let density: [c_float; 6] = [5.1, 5.3, 5.52, 3.94, 1.33, 0.69];

    let filename_cstr = CString::new(filename).unwrap();
    let extname_cstr = CString::new(extname).unwrap();

    /* open with write access the FITS file containing a primary array */
    if unsafe {
        fits_open_file(
            &raw mut fptr,
            filename_cstr.as_ptr(),
            READWRITE,
            &raw mut status,
        )
    } != 0
    {
        printerror(status);
    }

    /* append a new empty ASCII table onto the FITS file */
    if let Some(ref mut fptr_box) = fptr
        && unsafe {
            fits_create_tbl(
                fptr_box.as_mut(),
                ASCII_TBL,
                nrows,
                tfields,
                ttype.as_ptr(),
                tform.as_ptr(),
                tunit.as_ptr(),
                extname_cstr.as_ptr(),
                &raw mut status,
            )
        } != 0
    {
        printerror(status);
    }

    /* write names to the first column (character strings) */
    /* write diameters to the second column (longs) */
    /* write density to the third column (floats) */

    if let Some(ref mut fptr_box) = fptr {
        unsafe {
            fits_write_col(
                fptr_box.as_mut(),
                TSTRING,
                1,
                firstrow,
                firstelem,
                nrows,
                planet.as_ptr().cast(),
                &raw mut status,
            );
            fits_write_col(
                fptr_box.as_mut(),
                TLONG,
                2,
                firstrow,
                firstelem,
                nrows,
                diameter.as_ptr().cast(),
                &raw mut status,
            );
            fits_write_col(
                fptr_box.as_mut(),
                TFLOAT,
                3,
                firstrow,
                firstelem,
                nrows,
                density.as_ptr().cast(),
                &raw mut status,
            );
        }
    }

    if let Some(fptr_box) = fptr
        && unsafe { fits_close_file(Some(fptr_box), &raw mut status) } != 0
    {
        /* close the FITS file */
        printerror(status);
    }
}

/*--------------------------------------------------------------------------*/
/// Create a binary table extension containing 3 columns and 6 rows
fn writebintable() {
    let mut fptr: Option<Box<fitsfile>> = None;
    let mut status: c_int = 0;
    let mut hdutype: c_int = 0;
    let firstrow: LONGLONG = 1;
    let firstelem: LONGLONG = 1;

    let tfields: c_int = 3; /* table will have 3 columns */
    let nrows: LONGLONG = 6; /* table will have 6 rows    */

    let filename = "atestfil.fit"; /* name for new FITS file */
    let extname = "PLANETS_Binary"; /* extension name */

    /* define the name, datatype, and physical units for the 3 columns */
    let ttype_strs = [
        CString::new("Planet").unwrap(),
        CString::new("Diameter").unwrap(),
        CString::new("Density").unwrap(),
    ];
    let tform_strs = [
        CString::new("8a").unwrap(),
        CString::new("1J").unwrap(),
        CString::new("1E").unwrap(),
    ];
    let tunit_strs = [
        CString::new("").unwrap(),
        CString::new("km").unwrap(),
        CString::new("g/cm").unwrap(),
    ];

    let ttype: Vec<*const c_char> = ttype_strs.iter().map(|s| s.as_ptr()).collect();
    let tform: Vec<*const c_char> = tform_strs.iter().map(|s| s.as_ptr()).collect();
    let tunit: Vec<*const c_char> = tunit_strs.iter().map(|s| s.as_ptr()).collect();

    /* define the name diameter, and density of each planet */
    let planet_strs = [
        CString::new("Mercury").unwrap(),
        CString::new("Venus").unwrap(),
        CString::new("Earth").unwrap(),
        CString::new("Mars").unwrap(),
        CString::new("Jupiter").unwrap(),
        CString::new("Saturn").unwrap(),
    ];
    let planet: Vec<*const c_char> = planet_strs.iter().map(|s| s.as_ptr()).collect();
    let diameter: [c_long; 6] = [4880, 12112, 12742, 6800, 143000, 121000];
    let density: [c_float; 6] = [5.1, 5.3, 5.52, 3.94, 1.33, 0.69];

    let filename_cstr = CString::new(filename).unwrap();
    let extname_cstr = CString::new(extname).unwrap();

    /* open the FITS file containing a primary array and an ASCII table */
    if unsafe {
        fits_open_file(
            &raw mut fptr,
            filename_cstr.as_ptr(),
            READWRITE,
            &raw mut status,
        )
    } != 0
    {
        printerror(status);
    }

    if let Some(ref mut fptr_box) = fptr
        && unsafe { fits_movabs_hdu(fptr_box.as_mut(), 2, &raw mut hdutype, &raw mut status) } != 0
    {
        /* move to 2nd HDU */
        printerror(status);
    }

    /* append a new empty binary table onto the FITS file */
    if let Some(ref mut fptr_box) = fptr
        && unsafe {
            fits_create_tbl(
                fptr_box.as_mut(),
                BINARY_TBL,
                nrows,
                tfields,
                ttype.as_ptr(),
                tform.as_ptr(),
                tunit.as_ptr(),
                extname_cstr.as_ptr(),
                &raw mut status,
            )
        } != 0
    {
        printerror(status);
    }

    /* write names to the first column (character strings) */
    /* write diameters to the second column (longs) */
    /* write density to the third column (floats) */

    if let Some(ref mut fptr_box) = fptr {
        unsafe {
            fits_write_col(
                fptr_box.as_mut(),
                TSTRING,
                1,
                firstrow,
                firstelem,
                nrows,
                planet.as_ptr().cast(),
                &raw mut status,
            );
            fits_write_col(
                fptr_box.as_mut(),
                TLONG,
                2,
                firstrow,
                firstelem,
                nrows,
                diameter.as_ptr().cast(),
                &raw mut status,
            );
            fits_write_col(
                fptr_box.as_mut(),
                TFLOAT,
                3,
                firstrow,
                firstelem,
                nrows,
                density.as_ptr().cast(),
                &raw mut status,
            );
        }
    }

    if let Some(fptr_box) = fptr
        && unsafe { fits_close_file(Some(fptr_box), &raw mut status) } != 0
    {
        /* close the FITS file */
        printerror(status);
    }
}

/*--------------------------------------------------------------------------*/
/// Copy the 1st and 3rd HDUs from the input file to a new FITS file
fn copyhdu() {
    let mut infptr: Option<Box<fitsfile>> = None;
    let mut outfptr: Option<Box<fitsfile>> = None;

    let infilename = "atestfil.fit"; /* name for existing FITS file   */
    let outfilename = "btestfil.fit"; /* name for new FITS file        */

    let mut status: c_int = 0;
    let morekeys: c_int = 0;
    let mut hdutype: c_int = 0;

    let _ = remove_file(outfilename); /* Delete old file if it already exists */

    let infilename_cstr = CString::new(infilename).unwrap();
    let outfilename_cstr = CString::new(outfilename).unwrap();

    /* open the existing FITS file */
    if unsafe {
        fits_open_file(
            &raw mut infptr,
            infilename_cstr.as_ptr(),
            READONLY,
            &raw mut status,
        )
    } != 0
    {
        printerror(status);
    }

    if unsafe { fits_create_file(&raw mut outfptr, outfilename_cstr.as_ptr(), &raw mut status) }
        != 0
    {
        /*create FITS file*/
        printerror(status); /* call printerror if error occurs */
    }

    /* copy the primary array from the input file to the output file */
    if let (Some(infptr_box), Some(outfptr_box)) = (&mut infptr, &mut outfptr)
        && unsafe {
            fits_copy_hdu(
                infptr_box.as_mut(),
                outfptr_box.as_mut(),
                morekeys,
                &raw mut status,
            )
        } != 0
    {
        printerror(status);
    }

    /* move to the 3rd HDU in the input file */
    if let Some(ref mut infptr_box) = infptr
        && unsafe { fits_movabs_hdu(infptr_box.as_mut(), 3, &raw mut hdutype, &raw mut status) }
            != 0
    {
        printerror(status);
    }

    /* copy 3rd HDU from the input file to the output file (to 2nd HDU) */
    if let (Some(infptr_box), Some(outfptr_box)) = (&mut infptr, &mut outfptr)
        && unsafe {
            fits_copy_hdu(
                infptr_box.as_mut(),
                outfptr_box.as_mut(),
                morekeys,
                &raw mut status,
            )
        } != 0
    {
        printerror(status);
    }

    let mut close_status = 0;
    if let Some(outfptr_box) = outfptr
        && unsafe { fits_close_file(Some(outfptr_box), &raw mut close_status) } != 0
    {
        printerror(close_status);
    }
    if let Some(infptr_box) = infptr
        && unsafe { fits_close_file(Some(infptr_box), &raw mut status) } != 0
    {
        /* close files */
        printerror(status);
    }
}

/*--------------------------------------------------------------------------*/
/// select rows from an input table and copy them to the output table
fn selectrows() {
    let mut infptr: Option<Box<fitsfile>> = None;
    let mut outfptr: Option<Box<fitsfile>> = None;

    let mut card = [0u8; FLEN_CARD];
    let mut status: c_int = 0;
    let mut hdutype: c_int = 0;
    let mut nkeys: c_int = 0;
    let mut keypos: c_int = 0;
    let mut nfound: c_int = 0;
    let mut colnum: c_int = 0;
    let mut anynulls: c_int = 0;
    let mut naxes: [c_long; 2] = [0; 2];
    let frow: LONGLONG = 1;
    let felem: LONGLONG = 1;
    let mut noutrows: LONGLONG;
    let nullval: c_float = -99.0;
    let mut density: [c_float; 6] = [0.0; 6];

    let infilename = "atestfil.fit"; /* name for existing FITS file   */
    let outfilename = "btestfil.fit"; /* name for new FITS file        */

    let infilename_cstr = CString::new(infilename).unwrap();
    let outfilename_cstr = CString::new(outfilename).unwrap();

    /* open the existing FITS files */
    if unsafe {
        fits_open_file(
            &raw mut infptr,
            infilename_cstr.as_ptr(),
            READONLY,
            &raw mut status,
        )
    } != 0
        || unsafe {
            fits_open_file(
                &raw mut outfptr,
                outfilename_cstr.as_ptr(),
                READWRITE,
                &raw mut status,
            )
        } != 0
    {
        printerror(status);
    }

    /* move to the 3rd HDU in the input file (a binary table in this case) */
    if let Some(ref mut infptr_box) = infptr
        && unsafe { fits_movabs_hdu(infptr_box.as_mut(), 3, &raw mut hdutype, &raw mut status) }
            != 0
    {
        printerror(status);
    }

    if hdutype != BINARY_TBL {
        println!("Error: expected to find a binary table in this HDU");
        return;
    }

    /* move to the last (2rd) HDU in the output file */
    if let Some(ref mut outfptr_box) = outfptr
        && unsafe { fits_movabs_hdu(outfptr_box.as_mut(), 2, &raw mut hdutype, &raw mut status) }
            != 0
    {
        printerror(status);
    }

    /* create new extension in the output file */
    if let Some(ref mut outfptr_box) = outfptr
        && unsafe { fits_create_hdu(outfptr_box.as_mut(), &raw mut status) } != 0
    {
        printerror(status);
    }

    /* get number of keywords */
    if let Some(ref mut infptr_box) = infptr
        && unsafe {
            fits_get_hdrpos(
                infptr_box.as_mut(),
                &raw mut nkeys,
                &raw mut keypos,
                &raw mut status,
            )
        } != 0
    {
        printerror(status);
    }

    /* copy all the keywords from the input to the output extension */
    for ii in 1..=nkeys {
        if let (Some(infptr_box), Some(outfptr_box)) = (&mut infptr, &mut outfptr) {
            unsafe {
                fits_read_record(
                    infptr_box.as_mut(),
                    ii,
                    card.as_mut_ptr().cast::<c_char>(),
                    &raw mut status,
                );
                fits_write_record(
                    outfptr_box.as_mut(),
                    card.as_ptr().cast::<c_char>(),
                    &raw mut status,
                );
            }
        }
    }

    /* read the NAXIS1 and NAXIS2 keyword to get table size */
    if let Some(ref mut infptr_box) = infptr {
        let naxis_key = CString::new("NAXIS").unwrap();
        if unsafe {
            fits_read_keys_lng(
                infptr_box.as_mut(),
                naxis_key.as_ptr(),
                1,
                2,
                naxes.as_mut_ptr(),
                &raw mut nfound,
                &raw mut status,
            )
        } != 0
        {
            printerror(status);
        }
    }

    /* find which column contains the DENSITY values */
    if let Some(ref mut infptr_box) = infptr {
        let density_key = CString::new("density").unwrap();
        if unsafe {
            fits_get_colnum(
                infptr_box.as_mut(),
                CASEINSEN as c_int,
                density_key.as_ptr(),
                &raw mut colnum,
                &raw mut status,
            )
        } != 0
        {
            printerror(status);
        }
    }

    /* read the DENSITY column values */
    if let Some(ref mut infptr_box) = infptr
        && unsafe {
            fits_read_col(
                infptr_box.as_mut(),
                TFLOAT,
                colnum,
                frow,
                felem,
                naxes[1] as LONGLONG,
                core::ptr::from_ref::<c_float>(&nullval).cast(),
                density.as_mut_ptr().cast(),
                &raw mut anynulls,
                &raw mut status,
            )
        } != 0
    {
        printerror(status);
    }

    /* allocate buffer large enough for 1 row of the table */
    let buffer: *mut u8 = unsafe { malloc(naxes[0] as usize).cast::<u8>() };

    /* If the density is less than 3.0, copy the row to the output table */
    noutrows = 0;
    for irow in 1..=naxes[1] {
        if density[(irow - 1) as usize] < 3.0 {
            noutrows += 1;
            if let (Some(infptr_box), Some(outfptr_box)) = (&mut infptr, &mut outfptr) {
                unsafe {
                    fits_read_tblbytes(
                        infptr_box.as_mut(),
                        irow as LONGLONG,
                        1,
                        naxes[0] as LONGLONG,
                        buffer,
                        &raw mut status,
                    );
                    fits_write_tblbytes(
                        outfptr_box.as_mut(),
                        noutrows,
                        1,
                        naxes[0] as LONGLONG,
                        buffer,
                        &raw mut status,
                    );
                }
            }
        }
    }

    /* update the NAXIS2 keyword with the correct number of rows */
    if let Some(ref mut outfptr_box) = outfptr {
        let naxis2_key = CString::new("NAXIS2").unwrap();
        if unsafe {
            fits_update_key(
                outfptr_box.as_mut(),
                TLONG,
                naxis2_key.as_ptr(),
                core::ptr::from_ref::<c_long>(&(noutrows as c_long)).cast(),
                ptr::null(),
                &raw mut status,
            )
        } != 0
        {
            printerror(status);
        }
    }

    unsafe { free(buffer.cast()) };

    let mut close_status = 0;
    if let Some(outfptr_box) = outfptr
        && unsafe { fits_close_file(Some(outfptr_box), &raw mut close_status) } != 0
    {
        printerror(close_status);
    }
    if let Some(infptr_box) = infptr
        && unsafe { fits_close_file(Some(infptr_box), &raw mut status) } != 0
    {
        printerror(status);
    }
}

/*--------------------------------------------------------------------------*/
/// Print out all the header keywords in all extensions of a FITS file
fn readheader() {
    let mut fptr: Option<Box<fitsfile>> = None;

    let mut status: c_int = 0;
    let mut nkeys: c_int = 0;
    let mut keypos: c_int = 0;
    let mut hdutype: c_int = 0;
    let mut ii: c_int;
    let filename = "atestfil.fit"; /* name of existing FITS file   */
    let mut card = [0u8; FLEN_CARD]; /* standard string lengths defined in fitsioc.h */

    let filename_cstr = CString::new(filename).unwrap();

    if unsafe {
        fits_open_file(
            &raw mut fptr,
            filename_cstr.as_ptr(),
            READONLY,
            &raw mut status,
        )
    } != 0
    {
        printerror(status);
    }

    /* attempt to move to next HDU, until we get an EOF error */
    ii = 1;
    while let Some(ref mut fptr_box) = fptr {
        if unsafe { fits_movabs_hdu(fptr_box.as_mut(), ii, &raw mut hdutype, &raw mut status) } != 0
        {
            break;
        }

        /* get no. of keywords */
        if unsafe {
            fits_get_hdrpos(
                fptr_box.as_mut(),
                &raw mut nkeys,
                &raw mut keypos,
                &raw mut status,
            )
        } != 0
        {
            printerror(status);
        }

        println!("Header listing for HDU #{ii}:");
        for jj in 1..=nkeys {
            if unsafe {
                fits_read_record(
                    fptr_box.as_mut(),
                    jj,
                    card.as_mut_ptr().cast::<c_char>(),
                    &raw mut status,
                )
            } != 0
            {
                printerror(status);
            }

            let card_str = unsafe { CStr::from_ptr(card.as_ptr().cast::<c_char>()) };
            println!("{}", card_str.to_string_lossy()); /* print the keyword card */
        }
        println!("END\n"); /* terminate listing with END */
        ii += 1;
    }

    if status == END_OF_FILE {
        /* status values are defined in fitsioc.h */
        status = 0; /* got the expected EOF error; reset = 0  */
    } else {
        printerror(status); /* got an unexpected error                */
    }

    if let Some(fptr_box) = fptr
        && unsafe { fits_close_file(Some(fptr_box), &raw mut status) } != 0
    {
        printerror(status);
    }
}

/*--------------------------------------------------------------------------*/
/// Read a FITS image and determine the minimum and maximum pixel values
fn readimage() {
    let mut fptr: Option<Box<fitsfile>> = None;
    let mut status: c_int = 0;
    let mut nfound: c_int = 0;
    let mut anynull: c_int = 0;
    let mut naxes: [c_long; 2] = [0; 2];
    let mut fpixel: LONGLONG;
    let mut nbuffer: LONGLONG;
    let mut npixels: LONGLONG;

    const BUFFSIZE: usize = 1000;
    let mut datamin: c_float = 1.0e30;
    let mut datamax: c_float = -1.0e30;
    let nullval: c_float = 0.0;
    let mut buffer: [c_float; BUFFSIZE] = [0.0; BUFFSIZE];
    let filename = "atestfil.fit"; /* name of existing FITS file   */

    let filename_cstr = CString::new(filename).unwrap();

    if unsafe {
        fits_open_file(
            &raw mut fptr,
            filename_cstr.as_ptr(),
            READONLY,
            &raw mut status,
        )
    } != 0
    {
        printerror(status);
    }

    /* read the NAXIS1 and NAXIS2 keyword to get image size */
    if let Some(ref mut fptr_box) = fptr {
        let naxis_key = CString::new("NAXIS").unwrap();
        if unsafe {
            fits_read_keys_lng(
                fptr_box.as_mut(),
                naxis_key.as_ptr(),
                1,
                2,
                naxes.as_mut_ptr(),
                &raw mut nfound,
                &raw mut status,
            )
        } != 0
        {
            printerror(status);
        }
    }

    npixels = (naxes[0] * naxes[1]) as LONGLONG; /* number of pixels in the image */
    fpixel = 1;

    while npixels > 0 {
        nbuffer = npixels;
        if npixels > BUFFSIZE as LONGLONG {
            nbuffer = BUFFSIZE as LONGLONG; /* read as many pixels as will fit in buffer */
        }

        /* Note that even though the FITS images contains unsigned integer */
        /* pixel values (or more accurately, signed integer pixels with    */
        /* a bias of 32768),  this routine is reading the values into a    */
        /* float array.   Cfitsio automatically performs the datatype      */
        /* conversion in cases like this.                                  */

        if let Some(ref mut fptr_box) = fptr
            && unsafe {
                fits_read_img(
                    fptr_box.as_mut(),
                    TFLOAT,
                    fpixel,
                    nbuffer,
                    core::ptr::from_ref::<c_float>(&nullval).cast(),
                    buffer.as_mut_ptr().cast(),
                    &raw mut anynull,
                    &raw mut status,
                )
            } != 0
        {
            printerror(status);
        }

        for ii in 0..nbuffer {
            if buffer[ii as usize] < datamin {
                datamin = buffer[ii as usize];
            }

            if buffer[ii as usize] > datamax {
                datamax = buffer[ii as usize];
            }
        }
        npixels -= nbuffer; /* increment remaining number of pixels */
        fpixel += nbuffer; /* next pixel to be read in image */
    }

    println!("\nMin and max image pixels =  {datamin:.0}, {datamax:.0}");

    if let Some(fptr_box) = fptr
        && unsafe { fits_close_file(Some(fptr_box), &raw mut status) } != 0
    {
        printerror(status);
    }
}

/*--------------------------------------------------------------------------*/
/// read and print data values from an ASCII or binary table
fn readtable() {
    let mut fptr: Option<Box<fitsfile>> = None;
    let mut status: c_int = 0;
    let mut hdutype: c_int = 0;
    let mut nfound: c_int = 0;
    let mut anynull: c_int = 0;
    let frow: LONGLONG = 1;
    let felem: LONGLONG = 1;
    let nelem: LONGLONG = 6;
    let longnull: c_long = 0;
    let mut dia: [c_long; 6] = [0; 6];
    let floatnull: c_float = 0.0;
    let mut den: [c_float; 6] = [0.0; 6];
    let strnull = CString::new(" ").unwrap();

    let filename = "atestfil.fit"; /* name of existing FITS file   */

    let filename_cstr = CString::new(filename).unwrap();

    if unsafe {
        fits_open_file(
            &raw mut fptr,
            filename_cstr.as_ptr(),
            READONLY,
            &raw mut status,
        )
    } != 0
    {
        printerror(status);
    }

    /* allocate space for the column labels */
    let mut ttype_ptrs: Vec<*mut c_char> = Vec::new();
    for _i in 0..3 {
        let ptr = unsafe { malloc(FLEN_VALUE).cast::<c_char>() };
        ttype_ptrs.push(ptr);
    }

    /* allocate space for string column values */
    let mut name_ptrs: Vec<*mut c_char> = Vec::new();
    for _i in 0..6 {
        let ptr = unsafe { malloc(10).cast::<c_char>() };
        name_ptrs.push(ptr);
    }

    for hdunum in 2..=3 {
        /*read ASCII, then binary table */
        /* move to the HDU */
        if let Some(ref mut fptr_box) = fptr
            && unsafe {
                fits_movabs_hdu(fptr_box.as_mut(), hdunum, &raw mut hdutype, &raw mut status)
            } != 0
        {
            printerror(status);
        }

        if hdutype == ASCII_TBL {
            println!("\nReading ASCII table in HDU {hdunum}:");
        } else if hdutype == BINARY_TBL {
            println!("\nReading binary table in HDU {hdunum}:");
        } else {
            println!("Error: this HDU is not an ASCII or binary table");
            printerror(status);
        }

        /* read the column names from the TTYPEn keywords */
        if let Some(ref mut fptr_box) = fptr {
            let ttype_key = CString::new("TTYPE").unwrap();
            unsafe {
                fits_read_keys_str(
                    fptr_box.as_mut(),
                    ttype_key.as_ptr(),
                    1,
                    3,
                    ttype_ptrs.as_mut_ptr(),
                    &raw mut nfound,
                    &raw mut status,
                );
            }
        }

        let ttype0 = unsafe { CStr::from_ptr(ttype_ptrs[0]) }.to_string_lossy();
        let ttype1 = unsafe { CStr::from_ptr(ttype_ptrs[1]) }.to_string_lossy();
        let ttype2 = unsafe { CStr::from_ptr(ttype_ptrs[2]) }.to_string_lossy();
        println!(" Row  {ttype0:>10} {ttype1:>10} {ttype2:>10}");

        /*  read the columns */
        if let Some(ref mut fptr_box) = fptr {
            unsafe {
                fits_read_col(
                    fptr_box.as_mut(),
                    TSTRING,
                    1,
                    frow,
                    felem,
                    nelem,
                    strnull.as_ptr().cast(),
                    name_ptrs.as_mut_ptr().cast(),
                    &raw mut anynull,
                    &raw mut status,
                );
                fits_read_col(
                    fptr_box.as_mut(),
                    TLONG,
                    2,
                    frow,
                    felem,
                    nelem,
                    core::ptr::from_ref::<c_long>(&longnull).cast(),
                    dia.as_mut_ptr().cast(),
                    &raw mut anynull,
                    &raw mut status,
                );
                fits_read_col(
                    fptr_box.as_mut(),
                    TFLOAT,
                    3,
                    frow,
                    felem,
                    nelem,
                    core::ptr::from_ref::<c_float>(&floatnull).cast(),
                    den.as_mut_ptr().cast(),
                    &raw mut anynull,
                    &raw mut status,
                );
            }
        }

        for ii in 0..6 {
            let name_str = unsafe { CStr::from_ptr(name_ptrs[ii as usize]) }.to_string_lossy();
            println!(
                "{:5} {:>10} {:>10} {:>10.2}",
                ii + 1,
                name_str,
                dia[ii as usize],
                den[ii as usize]
            );
        }
    }

    /* free the memory for the column labels */
    for ptr in ttype_ptrs {
        unsafe { free(ptr.cast()) };
    }

    /* free the memory for the string column */
    for ptr in name_ptrs {
        unsafe { free(ptr.cast()) };
    }

    if let Some(fptr_box) = fptr
        && unsafe { fits_close_file(Some(fptr_box), &raw mut status) } != 0
    {
        printerror(status);
    }
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
