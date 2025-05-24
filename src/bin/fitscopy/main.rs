#![allow(deprecated)]

use std::ffi::CString;
use std::process::ExitCode;

use rsfitsio::fitsio::READONLY;

use rsfitsio::STDERR;

#[cfg(windows)]
use rsfitsio::__acrt_iob_func;

use rsfitsio::{aliases::fits_open_file, fitsio::fitsfile};

use rsfitsio::cfileio::ffclos as fits_close_file;
use rsfitsio::cfileio::ffinit as fits_create_file;
use rsfitsio::cfileio::ffrprt as fits_report_error;
use rsfitsio::edithdu::ffcpfl as fits_copy_file;

pub fn main() -> ExitCode {
    /* FITS file pointers defined in fitsio.h */
    let mut infptr: Option<Box<fitsfile>> = None;
    let mut outfptr: Option<Box<fitsfile>> = None;
    let mut status = 0; /* status must always be initialized = 0  */

    let args = std::env::args();

    if args.len() != 3 {
        println!("Usage:  fitscopy inputfile outputfile");
        println!();
        println!("Copy an input file to an output file, optionally filtering");
        println!("the file in the process.  This seemingly simple program can");
        println!("apply powerful filters which transform the input file as");
        println!("it is being copied.  Filters may be used to extract a");
        println!("subimage from a larger image, select rows from a table,");
        println!("filter a table with a GTI time extension or a SAO region file,");
        println!("create or delete columns in a table, create an image by");
        println!("binning (histogramming) 2 table columns, and convert IRAF");
        println!("format *.imh or raw binary data files into FITS images.");
        println!("See the CFITSIO User's Guide for a complete description of");
        println!("the Extended File Name filtering syntax.");
        println!();
        println!("Examples:");
        println!();
        println!("fitscopy in.fit out.fit                   (simple file copy)");
        println!("fitscopy - -                              (stdin to stdout)");
        println!("fitscopy in.fit[11:50,21:60] out.fit      (copy a subimage)");
        println!("fitscopy iniraf.imh out.fit               (IRAF image to FITS)");
        println!("fitscopy in.dat[i512,512] out.fit         (raw array to FITS)");
        println!("fitscopy in.fit[events][pi>35] out.fit    (copy rows with pi>35)");
        println!("fitscopy in.fit[events][bin X,Y] out.fit  (bin an image) ");
        println!("fitscopy in.fit[events][col x=.9*y] out.fit        (new x column)");
        println!("fitscopy in.fit[events][gtifilter()] out.fit       (time filter)");
        println!("fitscopy in.fit[2][regfilter(\"pow.reg\")] out.fit (spatial filter)");
        println!();
        println!("Note that it may be necessary to enclose the input file name");
        println!("in single quote characters on the Unix command line.");

        return ExitCode::from(0);
    }

    let args: Vec<String> = args.collect();
    let infile = CString::new(args[1].clone()).unwrap();
    let outfile = CString::new(args[2].clone()).unwrap();

    unsafe {
        /* Open the input file */
        if fits_open_file(&mut infptr, infile.as_ptr(), READONLY, &mut status) == 0 {
            /* Create the output file */
            if fits_create_file(&mut outfptr, outfile.as_ptr(), &mut status) == 0 {
                /* copy the previous, current, and following HDUs */
                fits_copy_file(
                    infptr.as_deref_mut().unwrap(),
                    outfptr.as_deref_mut().unwrap(),
                    1,
                    1,
                    1,
                    &mut status,
                );

                let mut outfptr = outfptr.unwrap();
                fits_close_file(outfptr.as_mut(), &mut status);
            }

            let mut infptr = infptr.unwrap();
            fits_close_file(infptr.as_mut(), &mut status);
        }

        /* if error occured, print out error message */
        if status != 0 {
            fits_report_error(STDERR!(), status);
        }
    }

    ExitCode::from(status as u8)
}
