//! ****************************************************************************
//! Function
//! wrtout: print messages in the streams of stdout and out.
//! wrterr: print erro messages  in the streams of stderr and out.
//! wrtferr: print cfitsio erro messages in the streams of stderr and out.
//! wrtwrn: print warning messages in the streams of stdout and out.
//! wrtsep: print seperators.
//! num_err_wrn: Return the number of errors and  warnings.
//!
//! *****************************************************************************
#![warn(missing_docs)]

/* Transpiled from cfitsio/utilities/fvrf_misc.c */

use std::cell::Cell;
use std::cmp::Ordering;

use rsfitsio::aliases::rust_api::{fits_clear_errmsg, fits_get_errstatus, fits_read_errmsg};
use rsfitsio::c_types::{c_char, c_int};
use rsfitsio::fitsio::FLEN_ERRMSG;
use rsfitsio::wrappers::{strcat_safe, strlen_safe};

use crate::common::*;
use crate::fvrf_file::close_report;
use crate::{pf, scat, spf};

thread_local! { static NWRNS: Cell<c_int> = const { Cell::new(0) }; }
thread_local! { static NERRS: Cell<c_int> = const { Cell::new(0) }; }
/* static char temp[512]; -- a scratch buffer, made local at each use site */
const TEMP_LEN: usize = 512;

pub(crate) fn num_err_wrn(num_err: &mut c_int, num_wrn: &mut c_int) {
    *num_wrn = NWRNS.with(Cell::get);
    *num_err = NERRS.with(Cell::get);
}

pub(crate) fn reset_err_wrn() {
    NWRNS.with(|c| c.set(0));
    NERRS.with(|c| c.set(0));
}

pub(crate) fn wrtout(out: Out, mess: &[c_char]) -> c_int {
    if out != Out::Null {
        pf!(out; CS(mess), "\n");
    }
    if out == Out::Stdout {
        fflush_out(Out::Stdout);
    }
    0
}

/* wrtout() for a plain Rust literal */
pub(crate) fn wrtout_str(out: Out, mess: &str) -> c_int {
    if out != Out::Null {
        pf!(out; mess, "\n");
    }
    if out == Out::Stdout {
        fflush_out(Out::Stdout);
    }
    0
}

pub(crate) fn wrtwrn(out: Out, mess: &[c_char], isheasarc: c_int) -> c_int {
    if err_report() != 0 {
        return 0;
    } /* Don't print the warnings */
    if heasarc_conv() == 0 && isheasarc != 0 {
        return 0;
    } /* heasarc warnings  but with
    heasarc convention turns off */
    let nwrns = NWRNS.with(|c| {
        let v = c.get();
        c.set(v + 1);
        v
    }) + 1;

    let mut temp: [c_char; TEMP_LEN] = [0; TEMP_LEN];
    spf!(temp; "*** Warning: ", CS(mess));
    if isheasarc != 0 {
        scat!(temp; " (HEASARC Convention)");
    }
    print_fmt(out, &temp, 13);
    /*    if(nwrns > MAXWRNS ) {
         fprintf(stderr,"??? Too many Warnings! I give up...\n");

    }  */
    nwrns
}

pub(crate) fn wrtwrn_str(out: Out, mess: &str, isheasarc: c_int) -> c_int {
    let mut b: [c_char; TEMP_LEN] = [0; TEMP_LEN];
    spf!(b; mess);
    wrtwrn(out, &b, isheasarc)
}

pub(crate) fn wrterr(out: Out, mess: &[c_char], severity: c_int) -> c_int {
    if severity < err_report() {
        fits_clear_errmsg();
        return 0;
    }
    let nerrs = NERRS.with(|c| {
        let v = c.get();
        c.set(v + 1);
        v
    }) + 1;

    let mut temp: [c_char; TEMP_LEN] = [0; TEMP_LEN];
    spf!(temp; "*** Error:   ", CS(mess));
    if out != Out::Null {
        if out != Out::Stdout && out != Out::Stderr {
            print_fmt(out, &temp, 13);
        }
        /*
           if ERR2OUT is defined, then error messages will be sent to the
           stdout stream rather than to stderr
        */
        /* #ifdef ERR2OUT ... #else */
        print_fmt(Out::Stderr, &temp, 13);
    }

    if nerrs > MAXERRORS {
        pf!(Out::Stderr; "??? Too many Errors! I give up...\n");
        close_report(out);
        fflush_out(Out::Stdout);
        fflush_out(Out::Stderr);
        std::process::exit(1);
    }
    fits_clear_errmsg();
    nerrs
}

pub(crate) fn wrterr_str(out: Out, mess: &str, severity: c_int) -> c_int {
    let mut b: [c_char; TEMP_LEN] = [0; TEMP_LEN];
    spf!(b; mess);
    wrterr(out, &b, severity)
}

/* construct an error message: mess + cfitsio error */
pub(crate) fn wrtferr(out: Out, mess: &[c_char], status: &mut c_int, severity: c_int) -> c_int {
    let mut ttemp: [c_char; 255] = [0; 255];

    if severity < err_report() {
        fits_clear_errmsg();
        return 0;
    }
    let nerrs = NERRS.with(|c| {
        let v = c.get();
        c.set(v + 1);
        v
    }) + 1;

    let mut temp: [c_char; TEMP_LEN] = [0; TEMP_LEN];
    spf!(temp; "*** Error:   ", CS(mess));
    fits_get_errstatus(*status, &mut ttemp);
    strcat_safe(&mut temp, &ttemp);
    if out != Out::Null {
        if out != Out::Stdout && out != Out::Stderr {
            print_fmt(out, &temp, 13);
        }
        /* #ifdef ERR2OUT ... #else */
        print_fmt(Out::Stderr, &temp, 13);
    }

    *status = 0;
    fits_clear_errmsg();
    if nerrs > MAXERRORS {
        pf!(Out::Stderr; "??? Too many Errors! I give up...\n");
        close_report(out);
        fflush_out(Out::Stdout);
        fflush_out(Out::Stderr);
        std::process::exit(1);
    }
    nerrs
}

pub(crate) fn wrtferr_str(out: Out, mess: &str, status: &mut c_int, severity: c_int) -> c_int {
    let mut b: [c_char; TEMP_LEN] = [0; TEMP_LEN];
    spf!(b; mess);
    wrtferr(out, &b, status, severity)
}

/* dump the cfitsio stack */
pub(crate) fn wrtserr(out: Out, mess: &[c_char], status: &mut c_int, severity: c_int) -> c_int {
    /* char* errfmt = "             %.67s\n"; */
    let mut i;
    /* C declares tmp[20][80] but then prints tmp[nstack] with nstack possibly
    == 20; one extra row keeps that in bounds. */
    let mut tmp: [[c_char; FLEN_ERRMSG]; 21] = [[0; FLEN_ERRMSG]; 21];
    let mut nstack = 0usize;

    if severity < err_report() {
        fits_clear_errmsg();
        return 0;
    }
    let nerrs = NERRS.with(|c| {
        let v = c.get();
        c.set(v + 1);
        v
    }) + 1;

    let mut temp: [c_char; TEMP_LEN] = [0; TEMP_LEN];
    spf!(temp; "*** Error:   ", CS(mess), "(from CFITSIO error stack:)");
    while nstack < 20 {
        tmp[nstack][0] = 0;
        let r = fits_read_errmsg(&mut tmp[nstack]);
        if r == 0 && tmp[nstack][0] == 0 {
            break;
        }
        nstack += 1;
    }

    if out != Out::Null {
        if out != Out::Stdout && out != Out::Stderr {
            print_fmt(out, &temp, 13);
            i = 0;
            while i <= nstack {
                pf!(out; "             ", CSW(&tmp[i], 0, Some(67)), "\n");
                i += 1;
            }
        }

        /* #ifdef ERR2OUT ... #else */
        print_fmt(Out::Stderr, &temp, 13);
        i = 0;
        while i <= nstack {
            pf!(Out::Stderr; "             ", CSW(&tmp[i], 0, Some(67)), "\n");
            i += 1;
        }
    }

    *status = 0;
    fits_clear_errmsg();
    if nerrs > MAXERRORS {
        pf!(Out::Stderr; "??? Too many Errors! I give up...\n");
        close_report(out);
        fflush_out(Out::Stdout);
        fflush_out(Out::Stderr);
        std::process::exit(1);
    }
    nerrs
}

pub(crate) fn wrtserr_str(out: Out, mess: &str, status: &mut c_int, severity: c_int) -> c_int {
    let mut b: [c_char; TEMP_LEN] = [0; TEMP_LEN];
    spf!(b; mess);
    wrtserr(out, &b, status, severity)
}

/* Print output of messages in a 80 character record.
Continue lines are aligned. */
pub(crate) fn print_fmt(out: Out, temp: &[c_char], nprompt: c_int) {
    let mut j: usize;
    let clen: usize;
    let mut tmp: [c_char; 81] = [0; 81];
    /* The C builds a `static char cont_fmt[80]' of `nprompt' blanks followed by
    "%.67s\n" the first time it sees a new nprompt.  Every call site in
    fitsverify passes nprompt = 13, so the prompt is simply re-emitted here. */

    if out == Out::Null {
        return;
    }

    let slen = cstrlen(temp);
    let mut i: isize = slen as isize - 80;
    if i <= 0 {
        pf!(out; CSW(temp, 0, Some(80)), "\n");
    } else {
        let mut p = 0usize;
        clen = (80 - nprompt) as usize;
        strncpy_c(&mut tmp, &temp[p..], 80);
        tmp[80] = 0;
        if isprint_c(temp[p + 79]) && isprint_c(temp[p + 80]) && temp[p + 80] != 0 {
            j = 79;
            while temp[p + j] as u8 != b' ' && j > 0 {
                j -= 1;
            }
            p += j;
            while temp[p] as u8 == b' ' {
                p += 1;
            }
            tmp[j] = 0;
        } else if temp[p + 80] as u8 == b' ' {
            j = 80;
            while temp[p + j] as u8 == b' ' {
                j += 1;
            }
            p += j;
        } else {
            p += 80;
        }
        pf!(out; CSW(&tmp, 0, Some(80)), "\n");
        while temp[p] != 0 && i > 0 {
            strncpy_c(&mut tmp, &temp[p..], clen);
            tmp[clen] = 0;
            i = cstrlen(&temp[p..]) as isize - clen as isize;
            if i > 0
                && isprint_c(temp[p + clen - 1])
                && isprint_c(temp[p + clen])
                && temp[p + clen] != 0
            {
                j = clen;
                while temp[p + j] as u8 != b' ' && j > 0 {
                    j -= 1;
                }
                p += j;
                while temp[p] as u8 == b' ' {
                    p += 1;
                }
                tmp[j] = 0;
            } else if i > 0 && temp[p + clen] as u8 == b' ' {
                j = clen;
                while temp[p + j] as u8 == b' ' {
                    j += 1;
                }
                p += j;
            } else if i > 0 {
                p += clen;
            }
            /* fprintf(out, cont_fmt, tmp) */
            let mut v: Vec<u8> = Vec::new();
            v.extend(std::iter::repeat_n(b' ', nprompt as usize));
            CSW(&tmp, 0, Some(67)).put(&mut v);
            v.push(b'\n');
            fwrite_out(out, &v);
        }
    }
    if out == Out::Stdout {
        fflush_out(Out::Stdout);
    }
}

pub(crate) fn print_fmt_str(out: Out, s: &str, nprompt: c_int) {
    let mut b: [c_char; TEMP_LEN] = [0; TEMP_LEN];
    spf!(b; s);
    print_fmt(out, &b, nprompt);
}

/* print a line of char fill with string title in the middle */
pub(crate) fn wrtsep(out: Out, fill: c_char, title: &[c_char], nchar: c_int) {
    let ntitle: usize;
    let line: &mut Vec<c_char>;
    let mut p: usize;
    let first_end: usize;
    let mut i;
    let mut nchar = nchar;

    ntitle = strlen_safe(title);
    if ntitle as c_int > nchar {
        nchar = ntitle as c_int;
    }
    if nchar <= 0 {
        return;
    }
    let nchar = nchar as usize;
    let mut linebuf: Vec<c_char> = vec![0; nchar + 1];
    line = &mut linebuf;
    p = 0;
    if ntitle < 1 {
        i = 0;
        while i < nchar {
            line[p] = fill;
            p += 1;
            i += 1;
        }
        line[p] = 0;
    } else {
        first_end = (nchar - ntitle) / 2;
        i = 0;
        while i < first_end {
            line[p] = fill;
            p += 1;
            i += 1;
        }
        line[p] = 0;
        strcat_safe(line, title);
        p += ntitle;
        i = first_end + ntitle;
        while i < nchar {
            line[p] = fill;
            p += 1;
            i += 1;
        }
        line[p] = 0;
    }
    if out != Out::Null {
        pf!(out; CS(line), "\n");
    }
    if out == Out::Stdout {
        fflush_out(out);
    }
}

pub(crate) fn wrtsep_str(out: Out, fill: c_char, title: &str, nchar: c_int) {
    let mut b: [c_char; TEMP_LEN] = [0; TEMP_LEN];
    spf!(b; title);
    wrtsep(out, fill, &b, nchar);
}

/* comparison function for the FitsKey structure array */
pub(crate) fn compkey(key1: &FitsKey, key2: &FitsKey) -> Ordering {
    ccmp(&key1.kname, &key2.kname)
}

/* comparison function for the colname structure array */
pub(crate) fn compcol(col1: &ColName, col2: &ColName) -> Ordering {
    ccmp(&col1.name, &col2.name)
}

/* comparison function for the string pattern maching.
Equal when the pattern is a prefix of the candidate -- this is what lets
key_match() find every CRPIXn from the pattern "CRPIX". */
pub(crate) fn compstrp(pattern: &[c_char], candidate: &[c_char]) -> Ordering {
    let p = cbytes(pattern);
    let q = cbytes(candidate);
    if q.starts_with(p) {
        Ordering::Equal
    } else {
        p.cmp(q)
    }
}

/* comparison function for the string exact maching*/
pub(crate) fn compstre(str1: &[c_char], str2: &[c_char]) -> Ordering {
    ccmp(str1, str2)
}

#[cfg(test)]
mod tests {
    use super::*;
    use rsfitsio::fitsio::FLEN_KEYWORD;

    fn kw(s: &str) -> [c_char; FLEN_KEYWORD] {
        let mut b = [0 as c_char; FLEN_KEYWORD];
        for (i, &c) in s.as_bytes().iter().enumerate() {
            b[i] = c as c_char;
        }
        b
    }

    /* compstrp is the prefix comparator bsearch() uses for indexed keywords:
    it reports equal when the pattern is a prefix of the array element. */
    #[test]
    fn test_compstrp() {
        assert_eq!(compstrp(&kw("CRPIX"), &kw("CRPIX1")), Ordering::Equal);
        assert_eq!(compstrp(&kw("CRPIX"), &kw("CRPIX")), Ordering::Equal);
        assert_eq!(compstrp(&kw("CRPIX"), &kw("CRPIX12A")), Ordering::Equal);
        /* candidate shorter than the pattern is not a match, and sorts before it */
        assert_eq!(compstrp(&kw("CRPIX"), &kw("CRP")), Ordering::Greater);
        /* ordinary ordering either side */
        assert_eq!(compstrp(&kw("CRPIX"), &kw("CS")), Ordering::Less);
        assert_eq!(compstrp(&kw("CRPIX"), &kw("BITPIX")), Ordering::Greater);
    }

    #[test]
    fn test_compstre() {
        assert_eq!(compstre(&kw("NAXIS"), &kw("NAXIS")), Ordering::Equal);
        /* exact matching, so a longer element is not equal */
        assert_eq!(compstre(&kw("NAXIS"), &kw("NAXIS1")), Ordering::Less);
        assert_eq!(compstre(&kw("NAXIS1"), &kw("NAXIS")), Ordering::Greater);
    }

    #[test]
    fn test_compkey_and_compcol() {
        let a = FitsKey {
            kname: kw("AAA"),
            ..Default::default()
        };
        let b = FitsKey {
            kname: kw("AAB"),
            ..Default::default()
        };
        assert_eq!(compkey(&a, &b), Ordering::Less);
        assert_eq!(compkey(&a, &a.clone()), Ordering::Equal);

        let c1 = ColName {
            name: cvec(&kw("ALPHA")),
            index: 1,
        };
        let c2 = ColName {
            name: cvec(&kw("BETA")),
            index: 2,
        };
        assert_eq!(compcol(&c1, &c2), Ordering::Less);
    }
}
