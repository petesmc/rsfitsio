/* Transpiled from cfitsio/utilities/ftverify.c and fitsverify.c.

ftverify.c #includes fitsverify.c for the STANDALONE build, so both live
here: `main_ftverify' holds ftverify_work()/update_parfile() and the PIL
stubs, and main() below is fitsverify.c's main().  */

// kwdtyp, FitsHdu and friends keep their fitsverify C spellings.
#![allow(non_camel_case_types, non_snake_case, non_upper_case_globals)]
#![allow(clippy::upper_case_acronyms)]
// Several C-derived declarations are only referenced from one arm of a
// platform/feature branch, or exist purely to mirror the C's shape.
#![allow(dead_code)]
// The C has dead stores that a faithful transpile keeps: C-style declarations
// initialised and then immediately reassigned (`int status = 0; ... status =
// f();`), loop counters it bumps but never reads again, and fvrf_key.c's
// `temp1` buffer, which is filled and then unused. Keeping them preserves the
// line-by-line correspondence with the original.
#![allow(unused_assignments, unused_variables)]
// Same C-shape allowances the library takes (see src/lib.rs): the C declares
// its locals at the top of a function and indexes its arrays with an explicit
// counter, and its integer widths are target-dependent, so a cast that looks
// like a no-op on this target is load-bearing on another.
#![allow(
    clippy::needless_late_init,
    clippy::needless_range_loop,
    clippy::manual_range_contains,
    clippy::unnecessary_cast
)]

use std::process::ExitCode;

mod common;
mod fvrf_data;
mod fvrf_file;
mod fvrf_head;
mod fvrf_key;
mod fvrf_misc;

use bytemuck::cast_slice;
use rsfitsio::c_types::c_char;
use rsfitsio::cs;
use rsfitsio::fitsio::FLEN_FILENAME;

use crate::common::*;

/*---------------------------------------------------------------------------*/
/*  ftverify.c                                                               */
/*---------------------------------------------------------------------------*/
pub(crate) mod main_ftverify {
    use std::io::BufRead;

    use rsfitsio::aliases::rust_api::fits_get_version;
    use rsfitsio::c_types::{c_char, c_int};
    use rsfitsio::fitsio::{FILE_NOT_CREATED, FILE_NOT_OPENED, FLEN_FILENAME};

    use crate::common::*;
    use crate::fvrf_file::{get_total_err, get_total_warn};
    use crate::fvrf_head::{leave_early, verify_fits};
    use crate::fvrf_misc::*;
    use crate::{pf, spf};

    const MAXMSG: usize = 256;

    const PIL_LINESIZE: usize = 1024;

    /*---------------------------------------------------------------------------*/
    /* call work function to verify that infile conforms to the FITS
    standard and write report to the output file */
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn ftverify_work(
        infile: &[c_char],  /* I - Input file name (Fits) */
        outfile: &[c_char], /* I - Output file name (ASCII) */
        _prehead: c_int,
        prstat: c_int,
        errreport: &[c_char],
        _testdata: c_int,
        _testcsum: c_int,
        _testfill: c_int,
        heasarc_conv: c_int,
        testhierarch: c_int,
    ) -> c_int {
        let mut outfptr: Out = Out::Null;
        let mut list: Option<std::io::BufReader<std::fs::File>> = None;
        let mut status: c_int = 0;
        let mut vfstatus: c_int;
        let mut filestatus: c_int;
        let mut p: usize;
        let mut task: [c_char; 80] = [0; 80];
        let mut tversion: [c_char; 80] = [0; 80];
        let mut fversion: f32 = 0.0;
        let mut i: c_int;
        let mut nerrs: c_int;
        let mut nwarns: c_int;
        let mut msg: [c_char; MAXMSG] = [0; MAXMSG];
        let mut comm: [c_char; COMM_LEN] = [0; COMM_LEN];
        let mut infile: Vec<c_char> = infile.to_vec();
        infile.resize(FLEN_FILENAME.max(infile.len()), 0);

        /* determine 'Severe error", "Error", or "Warning" report level */
        if errreport[0] as u8 == b's' || errreport[0] as u8 == b'S' {
            set_err_report(2);
        }
        if errreport[0] as u8 == b'e' || errreport[0] as u8 == b'E' {
            set_err_report(1);
        }

        p = 0;
        if infile[p] as u8 == b'@' {
            p += 1;
            match std::fs::File::open(cstr_to_os(&infile[p..])) {
                Ok(f) => list = Some(std::io::BufReader::new(f)),
                Err(_) => {
                    pf!(Out::Stderr; "Cannot open the list file: ", CS(&infile[p..]), ".");
                    leave_early(Out::Null);
                    return FILE_NOT_OPENED;
                }
            }
        }

        /* headas_clobberfile(outfile) -- delete existing file if clobber=YES */
        headas_clobberfile(outfile);
        let op = outfile;

        /* test if writing output log to a disk file */
        if prstat != 0
            && cstrlen(op) != 0
            && !ceq(op, b"STDOUT")
            && !ceq(op, b"STDERR")
            && std::fs::File::open(cstr_to_os(op)).is_ok()
        {
            spf!(msg; "Clobber is not set. Cannot overwrite the file", CS(op));
            status = FILE_NOT_CREATED;
            HD_ERROR_THROW(&msg, status);
            leave_early(Out::Null);
            return status;
        }

        if prstat != 0 && (cstrlen(op) == 0 || ceq(op, b"STDOUT")) {
            outfptr = Out::Stdout;
        } else if prstat != 0 && (cstrlen(op) == 0 || ceq(op, b"STDERR")) {
            outfptr = Out::Stderr;
        } else if prstat == 0 {
            outfptr = Out::Null;
        } else {
            match std::fs::File::create(cstr_to_os(op)) {
                Ok(f) => {
                    out_set_file(f);
                    outfptr = Out::File;
                }
                Err(_) => {
                    pf!(Out::Stderr;
                        "Error open output file ", CS(outfile), ". Using stdout instead.");
                    outfptr = Out::Stdout;
                }
            }
        }

        wrtout_str(outfptr, " ");
        fits_get_version(&mut fversion);
        get_toolname(&mut task);
        get_toolversion(&mut tversion);
        let fversion_int = ((fversion as f64 * 10000.0) + 0.5) as c_int;
        let fversion_major = fversion_int / 10000;
        let fversion_minor = (fversion_int % 10000) / 100;
        let fversion_patch = fversion_int % 100;
        spf!(comm;
            CS(&task), " ", CS(&tversion), " (CFITSIO V", fversion_major, ".", fversion_minor,
            ".", fversion_patch, ")");
        wrtsep(outfptr, b' ' as c_char, &comm, 60);
        i = 0;
        while comm[i as usize] != 0 {
            comm[i as usize] = b'-' as c_char;
            i += 1;
        }
        wrtsep(outfptr, b' ' as c_char, &comm, 60);
        wrtout_str(outfptr, " ");
        match err_report() {
            2 => {
                spf!(comm;
                    "Caution: Only checking for the most severe FITS format errors.");
                wrtout(outfptr, &comm);
            }
            1 => {}
            _ => {}
        }

        if heasarc_conv != 0 {
            spf!(comm; "HEASARC conventions are being checked.");
            wrtout(outfptr, &comm);
        }

        if testhierarch != 0 {
            spf!(comm; "ESO HIERARCH keywords are being checked.");
            wrtout(outfptr, &comm);
        }

        /* process each file */
        match list.as_mut() {
            None => {
                vfstatus = verify_fits(&mut infile, outfptr);

                if outfptr == Out::Null {
                    /* print one-line file summary */

                    /* verify_fits returns a non-zero status for catastrophic
                     * file I/O problems (an abort), and in this case total_err
                     * is not updated via close_report(), so we need to set
                     * nerrs accordingly for the one-line file summary. */
                    nerrs = if vfstatus != 0 { 1 } else { get_total_err() };

                    nwarns = get_total_warn();

                    filestatus = if (nerrs + nwarns) > 0 { 1 } else { 0 };
                    if filestatus != 0 {
                        if err_report() != 0 {
                            pf!(Out::Stdout;
                                "verification FAILED: ", CSW(&infile, -20, None), ", ", nerrs,
                                " errors\n");
                        } else {
                            pf!(Out::Stdout;
                                "verification FAILED: ", CSW(&infile, -20, None), ", ", nwarns,
                                " warnings and ", nerrs, " errors\n");
                        }
                    } else {
                        pf!(Out::Stdout;
                            "verification OK: ", CSW(&infile, -20, None), "\n");
                    }
                }
            }
            Some(list) => {
                let mut line = String::new();
                loop {
                    line.clear();
                    /* fgets(infile, FLEN_FILENAME, list) */
                    let n = match list.read_line(&mut line) {
                        Ok(0) | Err(_) => break,
                        Ok(n) => n,
                    };
                    let _ = n;
                    infile.iter_mut().for_each(|c| *c = 0);
                    set_cstr(&mut infile, line.as_bytes());

                    vfstatus = verify_fits(&mut infile, outfptr);

                    if outfptr == Out::Null {
                        /* print one-line file summary */
                        nerrs = if vfstatus != 0 { 1 } else { get_total_err() };

                        nwarns = get_total_warn();

                        filestatus = if (nerrs + nwarns) > 0 { 1 } else { 0 };
                        if filestatus != 0 {
                            if err_report() != 0 {
                                pf!(Out::Stdout;
                                    "verification FAILED: ", CSW(&infile, -20, None), ", ",
                                    nerrs, " errors\n");
                            } else {
                                pf!(Out::Stdout;
                                    "verification FAILED: ", CSW(&infile, -20, None), ", ",
                                    nwarns, " warnings and ", nerrs, " errors\n");
                            }
                        } else {
                            pf!(Out::Stdout;
                                "verification OK: ", CSW(&infile, -20, None), "\n");
                        }
                    }

                    for _i in 1..3 {
                        wrtout_str(outfptr, " ");
                    }
                }
            }
        }

        /* close the output file  */
        if outfptr != Out::Stdout && outfptr != Out::Null {
            out_close_file();
        }

        status
    }

    /******************************************************************************
     * Function
     *      update_parfile
     *
     * DESCRIPTION:
     *      Update the numerrs and numwrns parameters in the parfile.
     *
     *******************************************************************************/
    pub(crate) fn update_parfile(nerr: c_int, nwrn: c_int) {
        let mut status: c_int;
        let mut parname: [c_char; 32] = [0; 32];

        add_totalerr(nerr as i64);
        add_totalwrn(nwrn as i64);
        /* write the total accumulated total warnings and errors to the
        parfile */
        spf!(parname; "numwrns");
        status = PILPutInt(&parname, totalwrn() as c_int);
        if status != 0 {
            pf!(Out::Stderr; "Error to update the numwrns keyword.\n");
            status = 0;
        }
        spf!(parname; "numerrs");
        status = PILPutInt(&parname, totalerr() as c_int);
        if status != 0 {
            pf!(Out::Stderr; "Error to update the numerrs keyword.\n");
            status = 0;
        }
        let _ = status;
    }

    /*------------------------------------------------------------------
      The following are all dummy stub routines for functions that are
      only needed when ftverify is built in the HEADAS environment.
    --------------------------------------------------------------------*/

    fn PILGetFname(_parname: &[c_char], _filename: &mut [c_char]) -> c_int {
        0
    }

    fn PILGetString(_parname: &[c_char], _stringname: &mut [c_char]) -> c_int {
        0
    }

    fn PILGetBool(_parname: &[c_char], _intvalue: &mut c_int) -> c_int {
        0
    }

    fn PILPutInt(_parname: &[c_char], _intvalue: c_int) -> c_int {
        0
    }

    fn set_toolname(_taskname: &[c_char]) {}

    fn set_toolversion(_taskname: &[c_char]) {}

    fn get_toolname(taskname: &mut [c_char]) {
        spf!(taskname; "fitsverify");
    }

    fn get_toolversion(taskvers: &mut [c_char]) {
        spf!(taskvers; "4.22");
    }

    fn headas_clobberfile(_filename: &[c_char]) {}

    fn HD_ERROR_THROW(_msg: &[c_char], _status: c_int) {}

    /* An OS path from a NUL-terminated c_char buffer. */
    fn cstr_to_os(s: &[c_char]) -> std::ffi::OsString {
        #[cfg(unix)]
        {
            use std::os::unix::ffi::OsStringExt;
            std::ffi::OsString::from_vec(cbytes(s).to_vec())
        }
        #[cfg(not(unix))]
        {
            std::ffi::OsString::from(String::from_utf8_lossy(cbytes(s)).into_owned())
        }
    }
}

use main_ftverify::ftverify_work;

/*---------------------------------------------------------------------------*/
/*  fitsverify.c                                                             */
/*---------------------------------------------------------------------------*/
pub(crate) fn main() -> ExitCode {
    let mut status: c_int_alias = 0;
    let mut invalid = 0;
    let mut file1 = 0usize;
    let mut errormode: [c_char; 2] = [b'w' as c_char, 0];

    let argv: Vec<Vec<c_char>> = std::env::args_os()
        .map(|a| {
            let mut v: Vec<c_char> = os_bytes(&a).iter().map(|&b| b as c_char).collect();
            v.push(0);
            v
        })
        .collect();
    let argc = argv.len();

    if argc == 2 && ceq(&argv[1], b"-h") {
        print!(
            "{}",
            "fitsverify -- Verify that the input files conform to the FITS Standard.\n\
             \n\
             USAGE:   fitsverify filename ...  - verify one or more FITS files\n\
             \x20                                   (may use wildcard characters)\n\
             \x20  or    fitsverify @filelist.txt - verify a list of FITS files\n"
        );
        println!("      ");
        println!("   Optional flags:");
        println!("          -H  test ESO HIERARCH keywords");
        println!("          -l  list all header keywords");
        println!("          -q  quiet; print one-line pass/fail summary per file");
        println!("          -e  only test for error conditions (ignore warnings)");
        println!(" ");
        println!("   fitsverify exits with a status equal to the number of errors + warnings.");
        println!("        ");
        println!("EXAMPLES:");
        println!("     fitsverify -l m101.fits    - produce a detailed verificaton report of");
        println!("                                  a single file, including a keyword listing");
        println!("     fitsverify -q *.fits *.fit - verify all files with .fits or .fit");
        println!("                                  extensions, writing a 1-line pass/fail");
        println!("                                  message for each file");
        println!(" ");
        println!("DESCRIPTION:");
        println!("    ");
        println!("    This task reads one or more input FITS files and verifies that the");
        println!("    files conform to the specifications of the FITS Standard, Definition");
        println!("    of the Flexible Image Transport System (FITS), Version 3.0, available");
        println!("    online  at http://fits.gsfc.nasa.gov/.  The input filename template may");
        println!("    contain wildcard characters, in which case all matching files will be ");
        println!("    tested.  Alternatively, the name of an ASCII text file containing a list");
        println!("    of file names, one per line, may be entered preceded by an '@' character.");
        println!("    The following error or warning conditions will be reported:");
        println!("    ");
        println!("    ERROR CONDITIONS");
        println!("    ");
        println!("     - Mandatory keyword not present or out of order");
        println!("     - Mandatory keyword has wrong datatype or illegal value");
        println!("     - END header keyword is not present");
        println!("     - Sum of table column widths is inconsistent with NAXIS1 value");
        println!("     - BLANK keyword present in image with floating-point datatype");
        println!("     - TNULLn keyword present for floating-point binary table column");
        println!("     - Bit column has non-zero fill bits or is not left adjusted ");
        println!("     - ASCII TABLE column contains illegal value inconsistent with TFORMn");
        println!("     - Address to a variable length array not within the data heap ");
        println!("     - Extraneous bytes in the FITS file following the last HDU    ");
        println!("     - Mandatory keyword values not expressed in fixed format");
        println!("     - Mandatory keyword duplicated elsewhere in the header");
        println!("     - Header contains illegal ASCII character (not ASCII 32 - 126)");
        println!("     - Keyword name contains illegal character");
        println!("     - Keyword value field has illegal format");
        println!("     - Value and comment fields not separated by a slash character");
        println!("     - END keyword not filled with blanks in columns 9 - 80");
        println!("     - Reserved keyword with wrong datatype or illegal value");
        println!("     - XTENSION keyword in the primary array");
        println!("     - Column related keyword (TFIELDS, TTYPEn,TFORMn, etc.) in an image");
        println!("     - SIMPLE, EXTEND, or BLOCKED keyword in any extension");
        println!("     - BSCALE, BZERO, BUNIT, BLANK, DATAMAX, DATAMIN keywords in a table");
        println!("     - Table WCS keywords (TCTYPn, TCRPXn, TCRVLn, etc.) in an image");
        println!("     - TDIMn or THEAP keyword in an ASCII table ");
        println!("     - TBCOLn keyword in a Binary table");
        println!("     - THEAP keyword in a binary table that has PCOUNT = 0 ");
        println!("     - XTENSION, TFORMn, TDISPn or TDIMn value contains leading space(s)");
        println!("     - WCSAXES keyword appears after other WCS keywords");
        println!("     - Index of any WCS keyword (CRPIXn, CRVALn, etc.) greater than ");
        println!("       value of WCSAXES");
        println!("     - Index of any table column descriptor keyword (TTYPEn, TFORMn,");
        println!("       etc.) greater than value of TFIELDS");
        println!("     - TSCALn or TZEROn present for an ASCII, logical, or Bit column");
        println!("     - TDISPn value is inconsistent with the column datatype ");
        println!("     - Length of a variable length array greater than the maximum ");
        println!("       length as given by the TFORMn keyword");
        println!("     - ASCII table floating-point column value does not have decimal point(*)");
        println!("     - ASCII table numeric column value has embedded space character");
        println!("     - Logical column contains illegal value not equal to 'T', 'F', or 0");
        println!("     - Character string column contains non-ASCII text character");
        println!("     - Header fill bytes not all blanks");
        println!("     - Data fill bytes not all blanks in ASCII tables or all zeros ");
        println!("       in any other type of HDU ");
        println!("     - Gaps between defined ASCII table columns contain characters with");
        println!("       ASCII value > 127");
        println!("    ");
        println!("    WARNING CONDITIONS");
        println!("    ");
        println!("     - SIMPLE = F");
        println!("     - Presence of deprecated keywords BLOCKED or EPOCH");
        println!("     - 2 HDUs have identical EXTNAME, EXTVER, and EXTLEVEL values");
        println!("     - BSCALE or TSCALn value = 0.");
        println!("     - BLANK OR TNULLn value exceeds the legal range");
        println!("     - TFORMn has 'rAw' format and r is not a multiple of w");
        println!("     - DATE = 'dd/mm/yy' and yy is less than 10 (Y2K problem?)");
        println!("     - Index of any WCS keyword (CRPIXn, CRVALn, etc.) greater than");
        println!("       value of NAXIS, if the WCSAXES keyword is not present");
        println!("     - Duplicated keyword (except COMMENT, HISTORY, blank, etc.)");
        println!("     - Column name (TTYPEn) does not exist or contains characters ");
        println!("       other than letter, digit and underscore");
        println!("     - Calculated checksum inconsistent with CHECKSUM or DATASUM keyword");
        println!("        ");
        println!("    This is the stand alone version of the FTOOLS 'fverify' program.  It is");
        println!("    maintained by the HEASARC at NASA/GSFC.  Any comments about this program");
        println!("    should be submitted to http://heasarc.gsfc.nasa.gov/cgi-bin/ftoolshelp");

        return ExitCode::from(0);
    }

    set_prhead(0); /* don't print header by default */
    set_prstat(1); /* print HDU summary by default */
    set_testhierarch(0); /* do not test ESO HIERARCH keywords by default */

    /* check for flags on the command line */
    for ii in 1..argc {
        if argv[ii][0] as u8 != b'-' || ceq(&argv[ii], b"-") {
            file1 = ii;
            break;
        }

        if ceq(&argv[ii], b"-l") {
            set_prhead(1);
        } else if ceq(&argv[ii], b"-H") {
            set_testhierarch(1);
        } else if ceq(&argv[ii], b"-e") {
            errormode[0] = b'e' as c_char;
            errormode[1] = 0;
        } else if ceq(&argv[ii], b"-q") {
            set_prstat(0);
        } else {
            invalid = 1;
        }
    }

    if invalid != 0 || argc == 1 || file1 == 0 {
        /*  invalid input, so print brief help message */

        println!();
        println!("fitsverify - test if the input file(s) conform to the FITS format.");
        println!();
        println!("Usage:  fitsverify filename ...   or   fitsverify @filelist.txt");
        println!();
        println!("  where 'filename' is a filename template (with optional wildcards), and");
        println!("        'filelist.txt' is an ASCII text file with a list of");
        println!("         FITS file names, one per line.");
        println!();
        println!("   Optional flags:");
        println!("          -H  test ESO HIERARCH keywords");
        println!("          -l  list all header keywords");
        println!("          -q  quiet; print one-line pass/fail summary per file");
        println!("          -e  only test for error conditions; don't issue warnings");
        println!();
        println!("Help:   fitsverify -h");
        return ExitCode::from(0);
    }

    /*
         call work function to verify that infile conforms to the FITS
         standard and write report to the output file.
    */
    for ii in file1..argc {
        status = ftverify_work(
            &argv[ii],      /* name of file to verify */
            cs!(c"STDOUT"), /* write report to this stream */
            prhead(),       /* print listing of header keywords? */
            prstat(),       /* print detailed summary report */
            &errormode,     /* report errors only, or errors and warnings */
            1,              /* test the data  */
            1,              /* test checksum, if checksum keywords are present */
            1,              /* test data fill areas (should contain all zeros */
            0,              /* do not test for conformance with HEASARC convensions */
            /*    that are not required by the FITS Standard */
            testhierarch(), /* test format of ESO HIERARCH keywords? */
        );

        if status != 0 {
            return ExitCode::from(status as u8);
        }
    }

    if (totalerr() + totalwrn()) > 255 {
        ExitCode::from(255)
    } else {
        ExitCode::from((totalerr() + totalwrn()) as u8)
    }
}

type c_int_alias = rsfitsio::c_types::c_int;

fn os_bytes(s: &std::ffi::OsStr) -> Vec<u8> {
    #[cfg(unix)]
    {
        use std::os::unix::ffi::OsStrExt;
        s.as_bytes().to_vec()
    }
    #[cfg(not(unix))]
    {
        s.to_string_lossy().into_owned().into_bytes()
    }
}

#[allow(dead_code)]
fn _unused() -> usize {
    FLEN_FILENAME
}
