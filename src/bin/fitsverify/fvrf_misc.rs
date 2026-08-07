/******************************************************************************
* Function
*      wrtout: print messages in the streams of stdout and out.
*      wrterr: print erro messages  in the streams of stderr and out.
*      wrtferr: print cfitsio erro messages in the streams of stderr and out.
*      wrtwrn: print warning messages in the streams of stdout and out.
*      wrtsep: print seperators.
*      num_err_wrn: Return the number of errors and  warnings.
*
*******************************************************************************/
/* Transpiled from cfitsio/utilities/fvrf_misc.c */

use std::sync::atomic::{AtomicI32, Ordering};

use rsfitsio::aliases::rust_api::{fits_clear_errmsg, fits_get_errstatus, fits_read_errmsg};
use rsfitsio::c_types::{c_char, c_int};
use rsfitsio::fitsio::FLEN_ERRMSG;
use rsfitsio::wrappers::{strcat_safe, strlen_safe};

use crate::common::*;
use crate::fvrf_file::close_report;
use crate::{pf, scat, spf};

static NWRNS: AtomicI32 = AtomicI32::new(0);
static NERRS: AtomicI32 = AtomicI32::new(0);
/* static char temp[512]; -- a scratch buffer, made local at each use site */
const TEMP_LEN: usize = 512;

pub fn num_err_wrn(num_err: &mut c_int, num_wrn: &mut c_int) {
    *num_wrn = NWRNS.load(Ordering::Relaxed);
    *num_err = NERRS.load(Ordering::Relaxed);
}

pub fn reset_err_wrn() {
    NWRNS.store(0, Ordering::Relaxed);
    NERRS.store(0, Ordering::Relaxed);
}

pub fn wrtout(out: Out, mess: &[c_char]) -> c_int {
    if out != Out::Null {
        pf!(out; CS(mess), "\n");
    }
    if out == Out::Stdout {
        fflush_out(Out::Stdout);
    }
    0
}

/* wrtout() for a plain Rust literal */
pub fn wrtout_str(out: Out, mess: &str) -> c_int {
    if out != Out::Null {
        pf!(out; mess, "\n");
    }
    if out == Out::Stdout {
        fflush_out(Out::Stdout);
    }
    0
}

pub fn wrtwrn(out: Out, mess: &[c_char], isheasarc: c_int) -> c_int {
    if err_report() != 0 {
        return 0;
    } /* Don't print the warnings */
    if heasarc_conv() == 0 && isheasarc != 0 {
        return 0;
    } /* heasarc warnings  but with
      heasarc convention turns off */
    let nwrns = NWRNS.fetch_add(1, Ordering::Relaxed) + 1;

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

pub fn wrtwrn_str(out: Out, mess: &str, isheasarc: c_int) -> c_int {
    let mut b: [c_char; TEMP_LEN] = [0; TEMP_LEN];
    spf!(b; mess);
    wrtwrn(out, &b, isheasarc)
}

pub fn wrterr(out: Out, mess: &[c_char], severity: c_int) -> c_int {
    if severity < err_report() {
        fits_clear_errmsg();
        return 0;
    }
    let nerrs = NERRS.fetch_add(1, Ordering::Relaxed) + 1;

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

pub fn wrterr_str(out: Out, mess: &str, severity: c_int) -> c_int {
    let mut b: [c_char; TEMP_LEN] = [0; TEMP_LEN];
    spf!(b; mess);
    wrterr(out, &b, severity)
}

/* construct an error message: mess + cfitsio error */
pub fn wrtferr(out: Out, mess: &[c_char], status: &mut c_int, severity: c_int) -> c_int {
    let mut ttemp: [c_char; 255] = [0; 255];

    if severity < err_report() {
        fits_clear_errmsg();
        return 0;
    }
    let nerrs = NERRS.fetch_add(1, Ordering::Relaxed) + 1;

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

pub fn wrtferr_str(out: Out, mess: &str, status: &mut c_int, severity: c_int) -> c_int {
    let mut b: [c_char; TEMP_LEN] = [0; TEMP_LEN];
    spf!(b; mess);
    wrtferr(out, &b, status, severity)
}

/* dump the cfitsio stack */
pub fn wrtserr(out: Out, mess: &[c_char], status: &mut c_int, severity: c_int) -> c_int {
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
    let nerrs = NERRS.fetch_add(1, Ordering::Relaxed) + 1;

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

pub fn wrtserr_str(out: Out, mess: &str, status: &mut c_int, severity: c_int) -> c_int {
    let mut b: [c_char; TEMP_LEN] = [0; TEMP_LEN];
    spf!(b; mess);
    wrtserr(out, &b, status, severity)
}

/* Print output of messages in a 80 character record.
    Continue lines are aligned. */
pub fn print_fmt(out: Out, temp: &[c_char], nprompt: c_int) {
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

pub fn print_fmt_str(out: Out, s: &str, nprompt: c_int) {
    let mut b: [c_char; TEMP_LEN] = [0; TEMP_LEN];
    spf!(b; s);
    print_fmt(out, &b, nprompt);
}

/* print a line of char fill with string title in the middle */
pub fn wrtsep(out: Out, fill: c_char, title: &[c_char], nchar: c_int) {
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

pub fn wrtsep_str(out: Out, fill: c_char, title: &str, nchar: c_int) {
    let mut b: [c_char; TEMP_LEN] = [0; TEMP_LEN];
    spf!(b; title);
    wrtsep(out, fill, &b, nchar);
}

/* comparison function for the FitsKey structure array */
pub fn compkey(key1: &FitsKey, key2: &FitsKey) -> c_int {
    strncmp_c(&key1.kname, &key2.kname, FLEN_KEYWORD_C)
}

const FLEN_KEYWORD_C: usize = rsfitsio::fitsio::FLEN_KEYWORD;

/* comparison function for the colname structure array */
pub fn compcol(col1: &ColName, col2: &ColName) -> c_int {
    strcmp_c(&col1.name, &col2.name)
}

/* comparison function for the string pattern maching*/
pub fn compstrp(str1: &[c_char], str2: &[c_char]) -> c_int {
    let mut p = 0usize;
    let mut q = 0usize;
    while str2[q] == str1[p] && str2[q] != 0 {
        p += 1;
        q += 1;
        if str1[p] == 0 {
            return 0; /* str2 is longer than str1, but
                      matched */
        }
    }
    (str1[p] as u8) as c_int - (str2[q] as u8) as c_int
}

/* comparison function for the string exact maching*/
pub fn compstre(str1: &[c_char], str2: &[c_char]) -> c_int {
    strcmp_c(str1, str2)
}


#[cfg(test)]
mod tests {
    use super::*;

    fn kw(s: &str) -> [c_char; FLEN_KEYWORD_C] {
        let mut b = [0 as c_char; FLEN_KEYWORD_C];
        for (i, &c) in s.as_bytes().iter().enumerate() {
            b[i] = c as c_char;
        }
        b
    }

    /* compstrp is the prefix comparator bsearch() uses for indexed keywords:
       it reports equal when the pattern is a prefix of the array element. */
    #[test]
    fn test_compstrp() {
        assert_eq!(compstrp(&kw("CRPIX"), &kw("CRPIX1")), 0);
        assert_eq!(compstrp(&kw("CRPIX"), &kw("CRPIX")), 0);
        assert_eq!(compstrp(&kw("CRPIX"), &kw("CRPIX12A")), 0);
        /* element shorter than the pattern sorts before it: 'I' - '\0' */
        assert_eq!(compstrp(&kw("CRPIX"), &kw("CRP")), b'I' as c_int);
        /* ordinary ordering either side */
        assert!(compstrp(&kw("CRPIX"), &kw("CS")) < 0);
        assert!(compstrp(&kw("CRPIX"), &kw("BITPIX")) > 0);
    }

    #[test]
    fn test_compstre() {
        assert_eq!(compstre(&kw("NAXIS"), &kw("NAXIS")), 0);
        /* exact matching, so a longer element is not equal */
        assert!(compstre(&kw("NAXIS"), &kw("NAXIS1")) < 0);
        assert!(compstre(&kw("NAXIS1"), &kw("NAXIS")) > 0);
    }

    #[test]
    fn test_compkey_and_compcol() {
        let a = FitsKey { kname: kw("AAA"), ..Default::default() };
        let b = FitsKey { kname: kw("AAB"), ..Default::default() };
        assert!(compkey(&a, &b) < 0);
        assert_eq!(compkey(&a, &a.clone()), 0);

        let c1 = ColName { name: cvec(&kw("ALPHA")), index: 1 };
        let c2 = ColName { name: cvec(&kw("BETA")), index: 2 };
        assert!(compcol(&c1, &c2) < 0);
    }
}
