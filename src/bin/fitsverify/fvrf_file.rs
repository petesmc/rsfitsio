/* Transpiled from cfitsio/utilities/fvrf_file.c

   `static HduName **hduname' becomes a thread-local Vec<HduName>; the C's
   malloc/calloc/free bookkeeping is dropped in favour of RAII.  */

use std::cell::{Cell, RefCell};

use bytemuck::cast_slice;
use rsfitsio::aliases::rust_api::{
    fits_clear_errmsg, fits_get_hduaddrll, fits_movrel_hdu,
};
use rsfitsio::buffers::ffmbyt_safe;
use rsfitsio::c_types::{c_char, c_int};
use rsfitsio::cs;
use rsfitsio::fitsio::{ASCII_TBL, BINARY_TBL, END_OF_FILE, FLEN_VALUE, IMAGE_HDU, LONGLONG, fitsfile};

use crate::common::*;
use crate::fvrf_misc::*;
use crate::main_ftverify::update_parfile;
use crate::{scat, spf};

thread_local! {
    static HDUNAME: RefCell<Vec<HduName>> = const { RefCell::new(Vec::new()) };
}
thread_local! { static TOTAL_ERR: Cell<c_int> = const { Cell::new(1) }; } /* initialzed to 1 in case fail to open file */
thread_local! { static TOTAL_WARN: Cell<c_int> = const { Cell::new(0) }; }

pub(crate) fn get_total_warn() -> c_int {
    TOTAL_WARN.with(Cell::get)
}

pub(crate) fn get_total_err() -> c_int {
    TOTAL_ERR.with(Cell::get)
}

/* Get the total hdu number and allocate the memory for hdu array */
fn init_hduname() {
    /* allocate memories for the hdu structure array  */
    HDUNAME.with(|h| {
        let mut h = h.borrow_mut();
        h.clear();
        for _i in 0..totalhdu() {
            let mut e = HduName::default();
            e.hdutype = -1;
            e.errnum = 0;
            e.wrnno = 0;
            e.extname[0] = 0;
            e.extver = 0;
            h.push(e);
        }
    });
}

/* set the hduname memeber hdutype, extname, extver */
pub(crate) fn set_hduname(
    hdunum: c_int,             /* hdu number */
    hdutype: c_int,            /* hdutype */
    extname: Option<&[c_char]>, /* extension name */
    extver: c_int,             /* extension version */
) {
    let i = (hdunum - 1) as usize;
    HDUNAME.with(|h| {
        let mut h = h.borrow_mut();
        h[i].hdutype = hdutype;
        match extname {
            Some(e) => {
                let n = std::cmp::min(cstrlen(e), FLEN_VALUE - 1);
                h[i].extname[..n].copy_from_slice(&e[..n]);
                h[i].extname[n] = 0;
            }
            None => h[i].extname[0] = 0,
        }
        h[i].extver = extver;
    });
}

/* get the total errors and total warnings in this hdu */
pub(crate) fn set_hduerr(hdunum: c_int /* hdu number */) {
    let i = (hdunum - 1) as usize;
    let mut errnum = 0;
    let mut wrnno = 0;
    num_err_wrn(&mut errnum, &mut wrnno);
    HDUNAME.with(|h| {
        let mut h = h.borrow_mut();
        h[i].errnum = errnum;
        h[i].wrnno = wrnno;
    });
    reset_err_wrn(); /* reset the error and warning counter */
}

/* set the basic information for hduname structure */
pub(crate) fn set_hdubasic(hdunum: c_int, hdutype: c_int) {
    set_hduname(hdunum, hdutype, None, 0);
    set_hduerr(hdunum);
}

/* test to see whether the two extension having the same name */
/* return 1: identical 0: different */
pub(crate) fn test_hduname(
    hdunum1: c_int, /* index of first hdu */
    hdunum2: c_int, /* index of second hdu */
) -> c_int {
    HDUNAME.with(|h| {
        let h = h.borrow();
        let p1 = &h[(hdunum1 - 1) as usize];
        let p2 = &h[(hdunum2 - 1) as usize];
        if cstrlen(&p1.extname) == 0 || cstrlen(&p2.extname) == 0 {
            return 0;
        }
        if cbytes(&p1.extname) == cbytes(&p2.extname)
            && p1.hdutype == p2.hdutype
            && p2.extver == p1.extver
            && hdunum1 != hdunum2
        {
            return 1;
        }
        0
    })
}

/* Added the error numbers */
fn total_errors(toterr: &mut c_int, totwrn: &mut c_int) {
    let mut ierr = 0;
    let mut iwrn = 0;
    *toterr = 0;
    *totwrn = 0;

    if totalhdu() == 0 {
        /* this means the file couldn't be opened */
        *toterr = 1;
        return;
    }

    HDUNAME.with(|h| {
        let h = h.borrow();
        for i in 0..totalhdu() as usize {
            *toterr += h[i].errnum;
            *totwrn += h[i].wrnno;
        }
    });
    /*check the end of file errors */
    num_err_wrn(&mut ierr, &mut iwrn);
    *toterr += ierr;
    *totwrn += iwrn;
}

/* print the extname, exttype, extver, errnum and wrnno in a  table */
fn hdus_summary(out: Out) {
    let mut ierr = 0;
    let mut iwrn = 0;
    let mut temp: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut temp1: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut comm: [c_char; COMM_LEN] = [0; COMM_LEN];

    wrtsep_str(out, b'+' as c_char, " Error Summary  ", 60);
    wrtout_str(out, " ");
    spf!(comm; " HDU#  Name (version)       Type             Warnings  Errors");
    wrtout(out, &comm);

    let (w0, e0) = HDUNAME.with(|h| {
        let h = h.borrow();
        (h[0].wrnno, h[0].errnum)
    });
    spf!(comm;
        " 1                          Primary Array    ", format!("{:<4}", w0), "      ", format!("{:<4}", e0), "  ");
    wrtout(out, &comm);
    for i in 2..=totalhdu() {
        let (extname, extver, wrnno, errnum, hdutype) = HDUNAME.with(|h| {
            let h = h.borrow();
            let p = &h[(i - 1) as usize];
            (p.extname, p.extver, p.wrnno, p.errnum, p.hdutype)
        });
        spf!(temp; CS(&extname));
        if extver != 0 && extver != -999 {
            spf!(temp1; " (", extver, ")");
            scat!(temp; CS(&temp1));
        }
        match hdutype {
            IMAGE_HDU => {
                spf!(comm; " ", format!("{:<5}", i), " ", CSW(&temp, -20, None), " Image Array      ",
                    format!("{:<4}", wrnno), "      ", format!("{:<4}", errnum), "  ");
                wrtout(out, &comm);
            }
            ASCII_TBL => {
                spf!(comm; " ", format!("{:<5}", i), " ", CSW(&temp, -20, None), " ASCII Table      ",
                    format!("{:<4}", wrnno), "      ", format!("{:<4}", errnum), "  ");
                wrtout(out, &comm);
            }
            BINARY_TBL => {
                spf!(comm; " ", format!("{:<5}", i), " ", CSW(&temp, -20, None), " Binary Table     ",
                    format!("{:<4}", wrnno), "      ", format!("{:<4}", errnum), "  ");
                wrtout(out, &comm);
            }
            _ => {
                spf!(comm; " ", format!("{:<5}", i), " ", CSW(&temp, -20, None), " Unknown HDU      ",
                    format!("{:<4}", wrnno), "      ", format!("{:<4}", errnum), "  ");
                wrtout(out, &comm);
            }
        }
    }
    /* check the end of file */
    num_err_wrn(&mut ierr, &mut iwrn);
    if iwrn != 0 || ierr != 0 {
        spf!(comm; " End-of-file ", format!("{:<30}", ""), "  ", format!("{:<4}", iwrn), "      ",
            format!("{:<4}", ierr), "  ");
        wrtout(out, &comm);
    }
    wrtout_str(out, " ");
}

fn destroy_hduname() {
    HDUNAME.with(|h| h.borrow_mut().clear());
}

/* Routine to test the extra bytes at the end of file */
pub(crate) fn test_end(infits: &mut fitsfile, out: Out) {
    let mut status: c_int = 0;
    let mut headstart: LONGLONG = 0;
    let mut datastart: LONGLONG = 0;
    let mut dataend: LONGLONG = 0;
    let mut hdutype: c_int = 0;
    let mut errmes: [c_char; ERRMES_LEN] = [0; ERRMES_LEN];

    /* check whether there are any HDU left */
    fits_movrel_hdu(infits, 1, Some(&mut hdutype), &mut status);
    if status == 0 {
        wrtout_str(out, "< End-of-File >");
        spf!(errmes; "There are extraneous HDU(s) beyond the end of last HDU.");
        wrterr(out, &errmes, 2);
        wrtout_str(out, " ");
        return;
    }

    if status != END_OF_FILE {
        wrtserr_str(out, "Bad HDU? ", &mut status, 2);
        return;
    }

    status = 0;
    fits_clear_errmsg();
    if fits_get_hduaddrll(
        infits,
        Some(&mut headstart),
        Some(&mut datastart),
        Some(&mut dataend),
        &mut status,
    ) != 0
    {
        wrtferr_str(out, "", &mut status, 1);
    }

    /* try to move to the last byte of this extension.  */
    if ffmbyt_safe(infits, dataend - 1, 0, &mut status) != 0 {
        spf!(errmes;
            "Error trying to read last byte of the file at byte ", dataend as i64, ".");
        wrterr(out, &errmes, 2);
        wrtout_str(out, "< End-of-File >");
        wrtout_str(out, " ");
        return;
    }

    /* try to move to what would be the first byte of the next extension.
      If successfull, we have a problem... */

    ffmbyt_safe(infits, dataend, 0, &mut status);
    if status == 0 {
        wrtout_str(out, "< End-of-File >");
        spf!(errmes; "File has extra byte(s) after last HDU at byte ", dataend as i64, ".");
        wrterr(out, &errmes, 2);
        wrtout_str(out, " ");
    }
}

/******************************************************************************
* Function
*      init_report
*
*
* DESCRIPTION:
*      Initialize the fverify report
*
*******************************************************************************/
pub(crate) fn init_report(
    out: Out,               /* output file */
    _rootnam: &[c_char],    /* input file name */
) {
    let mut comm: [c_char; COMM_LEN] = [0; COMM_LEN];
    spf!(comm; "\n", totalhdu(), " Header-Data Units in this file.");
    wrtout(out, &comm);
    wrtout_str(out, " ");

    reset_err_wrn();
    init_hduname();
}

/******************************************************************************
* Function
*      close_report
*
*
* DESCRIPTION:
*      Close the fverify report
*
*******************************************************************************/
pub(crate) fn close_report(out: Out /* output file */) {
    let mut numerrs: c_int = 0; /* number of the errors         */
    let mut numwrns: c_int = 0; /* number of the warnings       */
    let mut comm: [c_char; COMM_LEN] = [0; COMM_LEN];

    /* print out a summary of all the hdus */
    if prstat() != 0 {
        hdus_summary(out);
    }
    total_errors(&mut numerrs, &mut numwrns);

    TOTAL_WARN.with(|c| c.set(numwrns));
    TOTAL_ERR.with(|c| c.set(numerrs));

    /* get the total number of errors and warnnings */
    spf!(comm;
        "**** Verification found ", numwrns, " warning(s) and ", numerrs, " error(s). ****");
    wrtout(out, &comm);

    update_parfile(numerrs, numwrns);
    /* destroy the hdu name */
    destroy_hduname();
}

#[allow(dead_code)]
fn _unused() -> &'static [c_char] {
    cs!(c"")
}
