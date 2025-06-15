#![allow(deprecated)]

use std::ffi::CString;
use std::fs::{File, remove_file};
use std::io::{BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::process::ExitCode;
use std::process::exit;
use std::time::{Instant, SystemTime};

use libc::{c_char, c_int, c_long, c_short};
use rsfitsio::aliases::c_api::{
    fits_close_file, fits_create_file, fits_create_img, fits_create_tbl, fits_flush_file,
    fits_get_errstatus, fits_get_rowsize, fits_movabs_hdu, fits_movrel_hdu, fits_read_col_lng,
    fits_read_errmsg, fits_read_img_lng, fits_write_col_lng, fits_write_img_lng,
    fits_write_img_sht,
};
use rsfitsio::fitsio::{ASCII_TBL, BINARY_TBL, FLEN_ERRMSG, FLEN_STATUS, fitsfile};

/* size of the image */
const XSIZE: usize = 3000;
const YSIZE: usize = 3000;

/* size of data buffer */
const SHTSIZE: usize = 20000;

/* no. of rows in binary table */
const BROWS: usize = 2500000;

/* no. of rows in ASCII table */
const AROWS: usize = 400000;

pub fn main() -> ExitCode {
    /*************************************************************************
    This program tests the speed of writing/reading FITS files with cfitsio
    **************************************************************************/

    let mut fptr: Option<Box<fitsfile>> = None;
    let mut status: c_int = 0;
    let filename = "speedcc.fit";
    let mut buffer = [0u8; 2880];
    let rawloop: c_long = (XSIZE * YSIZE / 720) as c_long;
    let mut rate: f32;
    let mut size: f32;
    let mut elapcpu: f32 = 0.0;
    let mut cpufrac: f32;
    let mut elapse: f64 = 0.0;

    let mut sarray: [c_long; SHTSIZE] = [0; SHTSIZE];
    let ssarray: [c_short; SHTSIZE] = [0; SHTSIZE];

    let tbegin = SystemTime::now();

    let _ = remove_file(filename); /* Delete old file if it already exists */

    let mut diskfile = File::options()
        .create(true)
        .truncate(true)
        .read(true)
        .write(true)
        .open(filename)
        .expect("Failed to create file");

    print!("                                                ");
    println!(" SIZE / ELAPSE(%CPU) = RATE");
    print!("RAW fwrite (2880 bytes/loop)...                 ");
    marktime(&mut status);

    let mut buff_writer = BufWriter::new(&diskfile);
    for _ in 0..rawloop {
        if buff_writer.write_all(&buffer).is_err() {
            println!("write error");
        }
    }

    drop(buff_writer); // Ensure all data is flushed to disk

    gettime(&mut elapse, &mut elapcpu, &mut status);

    cpufrac = ((elapcpu as f64) / elapse * 100.0) as f32;
    size = (2880.0 * rawloop as f32) / 1000000.0;
    rate = (size as f64 / elapse) as f32;
    println!(" {size:4.1}MB/{elapse:6.3}s({cpufrac:3.0}) = {rate:5.2}MB/s");

    /* read back the binary records */
    diskfile.seek(SeekFrom::Start(0)).expect("Failed to seek");

    print!("RAW fread  (2880 bytes/loop)...                 ");
    marktime(&mut status);

    let mut buff_reader = BufReader::new(&diskfile);
    for _ in 0..rawloop {
        match buff_reader.read_exact(&mut buffer) {
            Ok(_) => {}
            Err(e) => {
                dbg!(e);
                println!("read error");
            }
        }
    }

    drop(buff_reader); // Ensure all data is flushed from the buffer

    gettime(&mut elapse, &mut elapcpu, &mut status);

    cpufrac = ((elapcpu as f64) / elapse * 100.0) as f32;
    size = (2880.0 * rawloop as f32) / 1000000.0;
    rate = (size as f64 / elapse) as f32;
    println!(" {size:4.1}MB/{elapse:6.3}s({cpufrac:3.0}) = {rate:5.2}MB/s");

    drop(diskfile);
    let _ = remove_file(filename);

    status = 0;

    let filename_cstr = CString::new(filename).unwrap();
    if unsafe { fits_create_file(&mut fptr, filename_cstr.as_ptr(), &mut status) } != 0 {
        /* create new FITS file */
        printerror(status);
    }

    if let Some(ref mut fptr_box) = fptr {
        if writesimage(fptr_box, &ssarray, &mut status) != 0 {
            printerror(status);
        }
    }

    if let Some(fptr_box) = fptr {
        if unsafe { fits_close_file(Some(fptr_box), &mut status) } != 0 {
            printerror(status);
        }
    }

    let _ = remove_file(filename); /* Delete old file if it already exists */
    fptr = None;

    if unsafe { fits_create_file(&mut fptr, filename_cstr.as_ptr(), &mut status) } != 0 {
        /* create new FITS file */
        printerror(status);
    }

    if let Some(ref mut fptr_box) = fptr {
        if writeimage(fptr_box, &sarray, &mut status) != 0 {
            printerror(status);
        }
    }

    if let Some(ref mut fptr_box) = fptr {
        if writebintable(fptr_box, &sarray, &mut status) != 0 {
            printerror(status);
        }
    }

    if let Some(ref mut fptr_box) = fptr {
        if writeasctable(fptr_box, &sarray, &mut status) != 0 {
            printerror(status);
        }
    }

    if let Some(ref mut fptr_box) = fptr {
        if readimage(fptr_box, &mut sarray, &mut status) != 0 {
            printerror(status);
        }
    }

    if let Some(ref mut fptr_box) = fptr {
        if readbtable(fptr_box, &mut sarray, &mut status) != 0 {
            printerror(status);
        }
    }

    if let Some(ref mut fptr_box) = fptr {
        if readatable(fptr_box, &mut sarray, &mut status) != 0 {
            printerror(status);
        }
    }

    if let Some(fptr_box) = fptr {
        if unsafe { fits_close_file(Some(fptr_box), &mut status) } != 0 {
            printerror(status);
        }
    }

    let elapsed_time = tbegin.elapsed().unwrap();
    let total_elapse = elapsed_time.as_secs_f64();
    println!("Total elapsed time = {total_elapse:.3}s, status = {status}");

    ExitCode::from(0)
}

/*--------------------------------------------------------------------------*/
/// Write the primary array containing a 2-D image
fn writeimage(fptr: &mut fitsfile, sarray: &[c_long], status: &mut c_int) -> c_int {
    let mut ii: c_long;

    let mut elapcpu: f32 = 0.0;
    
    let mut elapse: f64 = 0.0;

    /* initialize FITS image parameters */
    let bitpix: c_int = 32; /* 32-bit  signed integer pixel values       */
    let naxis: c_long = 2; /* 2-dimensional image                            */
    let naxes: [c_long; 2] = [XSIZE as c_long, YSIZE as c_long]; /* image size */

    /* write the required keywords for the primary array image */
    if unsafe {
        fits_create_img(
            fptr,
            bitpix,
            naxis as c_int,
            naxes.as_ptr() as *mut c_long,
            status,
        )
    } != 0
    {
        printerror(*status);
    }

    print!("\nWrite {XSIZE}x{YSIZE} I*4 image, {SHTSIZE} pixels/loop:   ");
    marktime(status);

    let nremain: c_long = (XSIZE * YSIZE) as c_long;
    ii = 1;
    while ii <= nremain {
        unsafe {
            fits_write_img_lng(fptr, 0, ii, SHTSIZE as c_long, sarray.as_ptr(), status);
        }
        ii += SHTSIZE as c_long;
    }

    unsafe { fits_flush_file(fptr, status) }; /* flush all buffers to disk */

    gettime(&mut elapse, &mut elapcpu, status);

    let cpufrac: f32 = ((elapcpu as f64) / elapse * 100.0) as f32;
    let size: f32 = (XSIZE as f32 * 4.0 * YSIZE as f32) / 1000000.0;
    let rate: f32 = (size as f64 / elapse) as f32;
    println!(" {size:4.1}MB/{elapse:6.3}s({cpufrac:3.0}) = {rate:5.2}MB/s");

    *status
}

/*--------------------------------------------------------------------------*/
/// Write the primary array containing a 2-D image
fn writesimage(fptr: &mut fitsfile, ssarray: &[c_short], status: &mut c_int) -> c_int {
    let mut ii: c_long;

    let mut elapcpu: f32 = 0.0;
    
    let mut elapse: f64 = 0.0;

    /* initialize FITS image parameters */
    let bitpix: c_int = 16; /* 16-bit  signed integer pixel values       */
    let naxis: c_long = 2; /* 2-dimensional image                            */
    let naxes: [c_long; 2] = [XSIZE as c_long, YSIZE as c_long]; /* image size */

    /* write the required keywords for the primary array image */
    if unsafe {
        fits_create_img(
            fptr,
            bitpix,
            naxis as c_int,
            naxes.as_ptr() as *mut c_long,
            status,
        )
    } != 0
    {
        printerror(*status);
    }

    print!("\nWrite {XSIZE}x{YSIZE} I*2 image, {SHTSIZE} pixels/loop:   ");
    marktime(status);

    let nremain: c_long = (XSIZE * YSIZE) as c_long;
    ii = 1;
    while ii <= nremain {
        unsafe {
            fits_write_img_sht(fptr, 0, ii, SHTSIZE as c_long, ssarray.as_ptr(), status);
        }
        ii += SHTSIZE as c_long;
    }

    unsafe { fits_flush_file(fptr, status) }; /* flush all buffers to disk */

    gettime(&mut elapse, &mut elapcpu, status);

    let cpufrac: f32 = ((elapcpu as f64) / elapse * 100.0) as f32;
    let size: f32 = (XSIZE as f32 * 2.0 * YSIZE as f32) / 1000000.0;
    let rate: f32 = (size as f64 / elapse) as f32;
    println!(" {size:4.1}MB/{elapse:6.3}s({cpufrac:3.0}) = {rate:5.2}MB/s");

    *status
}

/*--------------------------------------------------------------------------*/
/// Create a binary table extension containing 3 columns
fn writebintable(fptr: &mut fitsfile, sarray: &[c_long], status: &mut c_int) -> c_int {
    let tfields: c_int = 2;
    let mut nremain: c_long;
    let mut ntodo: c_long;
    let mut firstrow: c_long = 1;
    let firstelem: c_long = 1;
    let mut nrows: c_long = 0;

    let mut elapcpu: f32 = 0.0;
    
    let mut elapse: f64 = 0.0;

    let extname = CString::new("Speed_Test").unwrap(); /* extension name */

    /* define the name, datatype, and physical units for the columns */
    let ttype_strs = [
        CString::new("first").unwrap(),
        CString::new("second").unwrap(),
    ];
    let tform_strs = [CString::new("1J").unwrap(), CString::new("1J").unwrap()];
    let tunit_strs = [CString::new(" ").unwrap(), CString::new(" ").unwrap()];

    let ttype: Vec<*const c_char> = ttype_strs.iter().map(|s| s.as_ptr()).collect();
    let tform: Vec<*const c_char> = tform_strs.iter().map(|s| s.as_ptr()).collect();
    let tunit: Vec<*const c_char> = tunit_strs.iter().map(|s| s.as_ptr()).collect();

    /* append a new empty binary table onto the FITS file */
    if unsafe {
        fits_create_tbl(
            fptr,
            BINARY_TBL,
            BROWS as c_long,
            tfields,
            ttype.as_ptr(),
            tform.as_ptr(),
            tunit.as_ptr(),
            extname.as_ptr() as *const c_char,
            status,
        )
    } != 0
    {
        printerror(*status);
    }

    /* get table row size and optimum number of rows to write per loop */
    unsafe { fits_get_rowsize(fptr, &mut nrows, status) };
    nrows = minvalue(nrows, SHTSIZE as c_long);
    nremain = BROWS as c_long;

    print!("Write {BROWS:7}row x {tfields}col bintable {nrows:4} rows/loop:");
    marktime(status);

    while nremain > 0 {
        ntodo = minvalue(nrows, nremain);
        unsafe {
            fits_write_col_lng(fptr, 1, firstrow, firstelem, ntodo, sarray.as_ptr(), status);
            fits_write_col_lng(fptr, 2, firstrow, firstelem, ntodo, sarray.as_ptr(), status);
        }
        firstrow += ntodo;
        nremain -= ntodo;
    }

    unsafe { fits_flush_file(fptr, status) }; /* flush all buffers to disk */

    gettime(&mut elapse, &mut elapcpu, status);

    let cpufrac: f32 = ((elapcpu as f64) / elapse * 100.0) as f32;
    let size: f32 = (BROWS as f32 * 8.0) / 1000000.0;
    let rate: f32 = (size as f64 / elapse) as f32;
    println!(" {size:4.1}MB/{elapse:6.3}s({cpufrac:3.0}) = {rate:5.2}MB/s");

    *status
}

/*--------------------------------------------------------------------------*/
/// Create an ASCII table extension containing 2 columns
fn writeasctable(fptr: &mut fitsfile, sarray: &[c_long], status: &mut c_int) -> c_int {
    let tfields: c_int = 2;
    let mut nremain: c_long;
    let mut ntodo: c_long;
    let mut firstrow: c_long = 1;
    let firstelem: c_long = 1;
    let mut nrows: c_long = 0;

    let mut elapcpu: f32 = 0.0;
    
    let mut elapse: f64 = 0.0;

    let extname = CString::new("Speed_Test").unwrap(); /* extension name */

    /* define the name, datatype, and physical units for the columns */
    let ttype_strs = [
        CString::new("first").unwrap(),
        CString::new("second").unwrap(),
    ];
    let tform_strs = [CString::new("I6").unwrap(), CString::new("I6").unwrap()];
    let tunit_strs = [CString::new(" ").unwrap(), CString::new(" ").unwrap()];

    let ttype: Vec<*const c_char> = ttype_strs.iter().map(|s| s.as_ptr()).collect();
    let tform: Vec<*const c_char> = tform_strs.iter().map(|s| s.as_ptr()).collect();
    let tunit: Vec<*const c_char> = tunit_strs.iter().map(|s| s.as_ptr()).collect();

    /* append a new empty ASCII table onto the FITS file */
    if unsafe {
        fits_create_tbl(
            fptr,
            ASCII_TBL,
            AROWS as c_long,
            tfields,
            ttype.as_ptr(),
            tform.as_ptr(),
            tunit.as_ptr(),
            extname.as_ptr() as *const c_char,
            status,
        )
    } != 0
    {
        printerror(*status);
    }

    /* get table row size and optimum number of rows to write per loop */
    unsafe { fits_get_rowsize(fptr, &mut nrows, status) };
    nrows = minvalue(nrows, SHTSIZE as c_long);
    nremain = AROWS as c_long;

    print!("Write {AROWS:7}row x {tfields}col asctable {nrows:4} rows/loop:");
    marktime(status);

    while nremain > 0 {
        ntodo = minvalue(nrows, nremain);
        unsafe {
            fits_write_col_lng(fptr, 1, firstrow, firstelem, ntodo, sarray.as_ptr(), status);
            fits_write_col_lng(fptr, 2, firstrow, firstelem, ntodo, sarray.as_ptr(), status);
        }
        firstrow += ntodo;
        nremain -= ntodo;
    }

    unsafe { fits_flush_file(fptr, status) }; /* flush all buffers to disk */

    gettime(&mut elapse, &mut elapcpu, status);

    let cpufrac: f32 = ((elapcpu as f64) / elapse * 100.0) as f32;
    let size: f32 = (AROWS as f32 * 13.0) / 1000000.0;
    let rate: f32 = (size as f64 / elapse) as f32;
    println!(" {size:4.1}MB/{elapse:6.3}s({cpufrac:3.0}) = {rate:5.2}MB/s");

    *status
}

/*--------------------------------------------------------------------------*/
/// Read a FITS image
fn readimage(fptr: &mut fitsfile, sarray: &mut [c_long], status: &mut c_int) -> c_int {
    let mut anynull: c_int = 0;
    let mut hdutype: c_int = 0;

    let mut ii: c_long;
    let longnull: c_long = 0;

    let mut elapcpu: f32 = 0.0;
    
    let mut elapse: f64 = 0.0;

    /* move to the primary array */
    if unsafe { fits_movabs_hdu(fptr, 1, &mut hdutype, status) } != 0 {
        printerror(*status);
    }

    print!("\nRead back image                                 ");
    marktime(status);

    let nremain: c_long = (XSIZE * YSIZE) as c_long;
    ii = 1;
    while ii <= nremain {
        unsafe {
            fits_read_img_lng(
                fptr,
                0,
                ii,
                SHTSIZE as c_long,
                longnull,
                sarray.as_mut_ptr(),
                &mut anynull,
                status,
            );
        }
        ii += SHTSIZE as c_long;
    }

    gettime(&mut elapse, &mut elapcpu, status);

    let cpufrac: f32 = ((elapcpu as f64) / elapse * 100.0) as f32;
    let size: f32 = (XSIZE as f32 * 4.0 * YSIZE as f32) / 1000000.0;
    let rate: f32 = (size as f64 / elapse) as f32;
    println!(" {size:4.1}MB/{elapse:6.3}s({cpufrac:3.0}) = {rate:5.2}MB/s");

    *status
}

/*--------------------------------------------------------------------------*/
/// read and print data values from the binary table
fn readbtable(fptr: &mut fitsfile, sarray: &mut [c_long], status: &mut c_int) -> c_int {
    let mut hdutype: c_int = 0;
    let mut anynull: c_int = 0;
    let mut nremain: c_long;
    let mut ntodo: c_long;
    let mut firstrow: c_long = 1;
    let firstelem: c_long = 1;
    let mut nrows: c_long = 0;
    let lnull: c_long = 0;

    let mut elapcpu: f32 = 0.0;
    
    let mut elapse: f64 = 0.0;

    /* move to the table */
    if unsafe { fits_movrel_hdu(fptr, 1, &mut hdutype, status) } != 0 {
        printerror(*status);
    }

    /* get table row size and optimum number of rows to read per loop */
    unsafe { fits_get_rowsize(fptr, &mut nrows, status) };
    nrows = minvalue(nrows, SHTSIZE as c_long);

    /*  read the columns */
    nremain = BROWS as c_long;

    print!("Read back BINTABLE                              ");
    marktime(status);

    while nremain > 0 {
        ntodo = minvalue(nrows, nremain);
        unsafe {
            fits_read_col_lng(
                fptr,
                1,
                firstrow,
                firstelem,
                ntodo,
                lnull,
                sarray.as_mut_ptr(),
                &mut anynull,
                status,
            );
            fits_read_col_lng(
                fptr,
                2,
                firstrow,
                firstelem,
                ntodo,
                lnull,
                sarray.as_mut_ptr(),
                &mut anynull,
                status,
            );
        }
        firstrow += ntodo;
        nremain -= ntodo;
    }

    gettime(&mut elapse, &mut elapcpu, status);

    let cpufrac: f32 = ((elapcpu as f64) / elapse * 100.0) as f32;
    let size: f32 = (BROWS as f32 * 8.0) / 1000000.0;
    let rate: f32 = (size as f64 / elapse) as f32;
    println!(" {size:4.1}MB/{elapse:6.3}s({cpufrac:3.0}) = {rate:5.2}MB/s");

    *status
}

/*--------------------------------------------------------------------------*/
/// read and print data values from an ASCII or binary table
fn readatable(fptr: &mut fitsfile, sarray: &mut [c_long], status: &mut c_int) -> c_int {
    let mut hdutype: c_int = 0;
    let mut anynull: c_int = 0;
    let mut nremain: c_long;
    let mut ntodo: c_long;
    let mut firstrow: c_long = 1;
    let firstelem: c_long = 1;
    let mut nrows: c_long = 0;
    let lnull: c_long = 0;

    let mut elapcpu: f32 = 0.0;
    
    let mut elapse: f64 = 0.0;

    /* move to the table */
    if unsafe { fits_movrel_hdu(fptr, 1, &mut hdutype, status) } != 0 {
        printerror(*status);
    }

    /* get table row size and optimum number of rows to read per loop */
    unsafe { fits_get_rowsize(fptr, &mut nrows, status) };
    nrows = minvalue(nrows, SHTSIZE as c_long);

    /*  read the columns */
    nremain = AROWS as c_long;

    print!("Read back ASCII Table                           ");
    marktime(status);

    while nremain > 0 {
        ntodo = minvalue(nrows, nremain);
        unsafe {
            fits_read_col_lng(
                fptr,
                1,
                firstrow,
                firstelem,
                ntodo,
                lnull,
                sarray.as_mut_ptr(),
                &mut anynull,
                status,
            );
            fits_read_col_lng(
                fptr,
                2,
                firstrow,
                firstelem,
                ntodo,
                lnull,
                sarray.as_mut_ptr(),
                &mut anynull,
                status,
            );
        }
        firstrow += ntodo;
        nremain -= ntodo;
    }

    gettime(&mut elapse, &mut elapcpu, status);

    let cpufrac: f32 = ((elapcpu as f64) / elapse * 100.0) as f32;
    let size: f32 = (AROWS as f32 * 13.0) / 1000000.0;
    let rate: f32 = (size as f64 / elapse) as f32;
    println!(" {size:4.1}MB/{elapse:6.3}s({cpufrac:3.0}) = {rate:5.2}MB/s");

    *status
}

/*--------------------------------------------------------------------------*/
/// Print out cfitsio error messages and exit program
fn printerror(status: c_int) {
    let mut status_str = [0u8; FLEN_STATUS];
    let mut errmsg = [0u8; FLEN_ERRMSG];

    if status != 0 {
        eprintln!("\n*** Error occurred during program execution ***");
    }

    unsafe {
        fits_get_errstatus(status, status_str.as_mut_ptr() as *mut c_char);
    } /* get the error description */
    let status_str = String::from_utf8_lossy(&status_str);
    eprintln!("\nstatus = {status}: {status_str}");

    /* get first message; null if stack is empty */
    if unsafe { fits_read_errmsg(errmsg.as_mut_ptr() as *mut c_char) } != 0 {
        eprintln!("\nError message stack:");
        let msg = String::from_utf8_lossy(&errmsg);
        eprintln!(" {msg}");

        while unsafe { fits_read_errmsg(errmsg.as_mut_ptr() as *mut c_char) } != 0 {
            /* get remaining messages */
            let msg = String::from_utf8_lossy(&errmsg);
            eprintln!(" {msg}");
        }
    }

    exit(status); /* terminate the program, returning error status */
}

static mut START_TIME: Option<Instant> = None;

fn marktime(_status: &mut c_int) -> c_int {
    unsafe {
        START_TIME = Some(Instant::now());
    }
    0
}

fn gettime(elapse: &mut f64, elapcpu: &mut f32, _status: &mut c_int) -> c_int {
    unsafe {
        if let Some(start) = START_TIME {
            *elapse = start.elapsed().as_secs_f64();
            *elapcpu = *elapse as f32; // CPU time approximation
        }
    }
    0
}

fn minvalue(a: c_long, b: c_long) -> c_long {
    if a < b { a } else { b }
}
