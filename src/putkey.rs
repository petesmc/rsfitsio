/*  This file, putkey.rs, contains routines that write keywords to          */
/*  a FITS header.                                                         */

/*  The FITSIO software was written by William Pence at the High Energy    */
/*  Astrophysic Science Archive Research Center (HEASARC) at the NASA      */
/*  Goddard Space Flight Center.                                           */

use alloc::ffi::CString;
use core::ffi::CStr;
use core::num::ParseIntError;
use core::{cmp, mem};
use core::{slice, str};
use std::fs::File;
use std::io::BufRead;

use chrono::{DateTime, Datelike, Utc};

use crate::c_types::*;

use bytemuck::{cast_slice, cast_slice_mut};

use crate::fitscore::{
    ffbnfm_safe, ffbnfmll_safe, ffcrhd_safe, ffgabc_safe, ffgthd_safe, ffiblk, ffkeyn_safe,
    ffmahd_safe, ffmkky_safe, ffpmsg_slice, ffpmsg_str, ffrdef_safe, fftkey_safe, ffupch_safe,
    fits_strncasecmp,
};
use crate::getkey::ffgkys_safe;
use crate::imcompress::imcomp_init_table;
use crate::modkey::{ffdkey_safe, ffirec_safe, ffmnam_safe, ffucrd_safe};
use crate::relibc::header::stdio::snprintf_f64_decim;
use crate::{KeywordDatatype, fitsio2::*};
use crate::{TKeywords, wrappers::*};
use crate::{atoi, int_snprintf};
use crate::{bb, cs, parse_c_int};
use crate::{buffers::*, nullable_slice_cstr, raw_to_slice};
use crate::{fitsio::*, fmt_f64};

/*--------------------------------------------------------------------------*/
/// create an IMAGE extension following the current HDU. If the
/// current HDU is empty (contains no header keywords), then simply
/// write the required image (or primary array) keywords to the current
/// HDU.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffcrim(
    fptr: *mut fitsfile,  /* I - FITS file pointer           */
    bitpix: c_int,        /* I - bits per pixel              */
    naxis: c_int,         /* I - number of axes in the array */
    naxes: *const c_long, /* I - size of each axis           */
    status: *mut c_int,   /* IO - error status               */
) -> c_int {
    unsafe {
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);

        let naxes = match naxes.is_null() {
            true => &[],
            false => slice::from_raw_parts(naxes, naxis as usize),
        };

        ffcrim_safe(fptr, bitpix, naxis, naxes, status)
    }
}

/*--------------------------------------------------------------------------*/
/// create an IMAGE extension following the current HDU. If the
/// current HDU is empty (contains no header keywords), then simply
/// write the required image (or primary array) keywords to the current
/// HDU.
pub fn ffcrim_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer           */
    bitpix: c_int,       /* I - bits per pixel              */
    naxis: c_int,        /* I - number of axes in the array */
    naxes: &[c_long],    /* I - size of each axis           */
    status: &mut c_int,  /* IO - error status               */
) -> c_int {
    if *status > 0 {
        return *status;
    }

    if fptr.HDUposition != fptr.Fptr.curhdu {
        ffmahd_safe(fptr, (fptr.HDUposition) + 1, None, status);
    }

    /* create new extension if current header is not empty */
    let h = fptr.Fptr.get_headstart_as_slice();
    if fptr.Fptr.headend != h[fptr.Fptr.curhdu as usize] {
        ffcrhd_safe(fptr, status);
    }

    /* write the required header keywords */
    ffphpr_safe(
        fptr,
        TRUE as c_int,
        bitpix,
        naxis,
        naxes,
        0,
        1,
        TRUE as c_int,
        status,
    );

    *status
}

/*--------------------------------------------------------------------------*/
/// create an IMAGE extension following the current HDU. If the
/// current HDU is empty (contains no header keywords), then simply
/// write the required image (or primary array) keywords to the current
/// HDU.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffcrimll(
    fptr: *mut fitsfile,    /* I - FITS file pointer           */
    bitpix: c_int,          /* I - bits per pixel              */
    naxis: c_int,           /* I - number of axes in the array */
    naxes: *const LONGLONG, /* I - size of each axis           */
    status: *mut c_int,     /* IO - error status               */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        let naxes = slice::from_raw_parts(naxes, naxis as usize);

        ffcrimll_safe(fptr, bitpix, naxis, naxes, status)
    }
}

/*--------------------------------------------------------------------------*/
/// create an IMAGE extension following the current HDU. If the
/// current HDU is empty (contains no header keywords), then simply
/// write the required image (or primary array) keywords to the current
/// HDU.
pub fn ffcrimll_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer           */
    bitpix: c_int,       /* I - bits per pixel              */
    naxis: c_int,        /* I - number of axes in the array */
    naxes: &[LONGLONG],  /* I - size of each axis           */
    status: &mut c_int,  /* IO - error status               */
) -> c_int {
    if *status > 0 {
        return *status;
    }

    if fptr.HDUposition != fptr.Fptr.curhdu {
        ffmahd_safe(fptr, (fptr.HDUposition) + 1, None, status);
    }

    /* create new extension if current header is not empty */
    let headstart = fptr.Fptr.get_headstart_as_slice();
    if fptr.Fptr.headend != headstart[fptr.Fptr.curhdu as usize] {
        ffcrhd_safe(fptr, status);
    }

    /* write the required header keywords */
    ffphprll_safe(
        fptr,
        TRUE as c_int,
        bitpix,
        naxis,
        naxes,
        0,
        1,
        TRUE as c_int,
        status,
    );

    *status
}

/*--------------------------------------------------------------------------*/
/// Create a table extension in a FITS file.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffcrtb(
    fptr: *mut fitsfile,         /* I - FITS file pointer                        */
    tbltype: c_int,              /* I - type of table to create                  */
    naxis2: LONGLONG,            /* I - number of rows in the table              */
    tfields: c_int,              /* I - number of columns in the table           */
    ttype: *const *const c_char, /* I - name of each column                      */
    tform: *const *const c_char, /* I - value of TFORMn keyword for each column  */
    tunit: *const *const c_char, /* I - value of TUNITn keyword for each column  */
    extnm: *const c_char,        /* I - value of EXTNAME keyword, if any         */
    status: *mut c_int,          /* IO - error status                            */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        let tkeywords = TKeywords::new(tfields, ttype, tform, tunit);
        let (v_ttype, v_tform, v_tunit) = tkeywords.tkeywords_to_vecs();

        nullable_slice_cstr!(extnm);

        ffcrtb_safe(
            fptr,
            tbltype,
            naxis2,
            tfields,
            &v_ttype,
            &v_tform,
            v_tunit.as_deref(),
            extnm,
            status,
        )
    }
}

/*--------------------------------------------------------------------------*/
/// Create a table extension in a FITS file.
pub fn ffcrtb_safe(
    fptr: &mut fitsfile,         /* I - FITS file pointer                        */
    tbltype: c_int,              /* I - type of table to create                  */
    naxis2: LONGLONG,            /* I - number of rows in the table              */
    tfields: c_int,              /* I - number of columns in the table           */
    ttype: &[Option<&[c_char]>], /* I - name of each column                      */
    tform: &[&[c_char]],         /* I - value of TFORMn keyword for each column  */
    tunit: Option<&[Option<&[c_char]>]>, /* I - value of TUNITn keyword for each column  */
    extnm: Option<&[c_char]>,    /* I - value of EXTNAME keyword, if any         */
    status: &mut c_int,          /* IO - error status                            */
) -> c_int {
    let naxis1: LONGLONG = 0;

    if *status > 0 {
        return *status;
    }

    if fptr.HDUposition != fptr.Fptr.curhdu {
        ffmahd_safe(fptr, (fptr.HDUposition) + 1, None, status);
    }

    /* create new extension if current header is not empty */
    let headstart = fptr.Fptr.get_headstart_as_slice();
    if fptr.Fptr.headend != headstart[fptr.Fptr.curhdu as usize] {
        ffcrhd_safe(fptr, status);
    }

    if fptr.Fptr.curhdu == 0 {
        /* have to create dummy primary array */
        ffcrim_safe(fptr, 16, 0, &[], status);
        ffcrhd_safe(fptr, status);
    }

    if tbltype == BINARY_TBL {
        /* write the required header keywords. This will write PCOUNT = 0 */
        ffphbn_safe(fptr, naxis2, tfields, ttype, tform, tunit, extnm, 0, status);
    } else if tbltype == ASCII_TBL {
        /* write the required header keywords */
        /* default values for naxis1 and tbcol will be calculated */
        ffphtb_safe(
            fptr, naxis1, naxis2, tfields, ttype, None, tform, tunit, extnm, status,
        );
    } else {
        *status = NOT_TABLE;
    }
    *status
}

/*--------------------------------------------------------------------------*/
/// write STANDARD set of required primary header keywords
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffphps(
    fptr: *mut fitsfile,  /* I - FITS file pointer                        */
    bitpix: c_int,        /* I - number of bits per data value pixel      */
    naxis: c_int,         /* I - number of axes in the data array         */
    naxes: *const c_long, /* I - length of each data axis                 */
    status: *mut c_int,   /* IO - error status                            */
) -> c_int {
    unsafe {
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);

        let naxes = if naxes.is_null() {
            &[]
        } else {
            slice::from_raw_parts(naxes, naxis as usize)
        };

        ffphps_safe(fptr, bitpix, naxis, naxes, status)
    }
}

/*--------------------------------------------------------------------------*/
/// write STANDARD set of required primary header keywords
pub fn ffphps_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer                        */
    bitpix: c_int,       /* I - number of bits per data value pixel      */
    naxis: c_int,        /* I - number of axes in the data array         */
    naxes: &[c_long],    /* I - length of each data axis                 */
    status: &mut c_int,  /* IO - error status                            */
) -> c_int {
    let simple: c_int = 1; /* does file conform to FITS standard? 1/0  */
    let pcount: LONGLONG = 0; /* number of group parameters (usually 0)   */
    let gcount: LONGLONG = 1; /* number of random groups (usually 1 or 0) */
    let extend: c_int = 1; /* may FITS file have extensions?           */

    ffphpr_safe(
        fptr, simple, bitpix, naxis, naxes, pcount, gcount, extend, status,
    );

    *status
}

/*--------------------------------------------------------------------------*/
/// Write STANDARD set of required primary header keywords
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffphpsll(
    fptr: *mut fitsfile,    /* I - FITS file pointer                        */
    bitpix: c_int,          /* I - number of bits per data value pixel      */
    naxis: c_int,           /* I - number of axes in the data array         */
    naxes: *const LONGLONG, /* I - length of each data axis                 */
    status: *mut c_int,     /* IO - error status                            */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        let naxes = slice::from_raw_parts(naxes, naxis as usize);

        ffphpsll_safe(fptr, bitpix, naxis, naxes, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Write STANDARD set of required primary header keywords
pub fn ffphpsll_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer                        */
    bitpix: c_int,       /* I - number of bits per data value pixel      */
    naxis: c_int,        /* I - number of axes in the data array         */
    naxes: &[LONGLONG],  /* I - length of each data axis                 */
    status: &mut c_int,  /* IO - error status                            */
) -> c_int {
    let simple: c_int = 1; /* does file conform to FITS standard? 1/0  */
    let pcount: LONGLONG = 0; /* number of group parameters (usually 0)   */
    let gcount: LONGLONG = 1; /* number of random groups (usually 1 or 0) */
    let extend: c_int = 1; /* may FITS file have extensions?           */

    ffphprll_safe(
        fptr, simple, bitpix, naxis, naxes, pcount, gcount, extend, status,
    );
    *status
}

/*--------------------------------------------------------------------------*/
/// write required primary header keywords
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffphpr(
    fptr: *mut fitsfile,  /* I - FITS file pointer                        */
    simple: c_int,        /* I - does file conform to FITS standard? 1/0  */
    bitpix: c_int,        /* I - number of bits per data value pixel      */
    naxis: c_int,         /* I - number of axes in the data array         */
    naxes: *const c_long, /* I - length of each data axis                 */
    pcount: LONGLONG,     /* I - number of group parameters (usually 0)   */
    gcount: LONGLONG,     /* I - number of random groups (usually 1 or 0) */
    extend: c_int,        /* I - may FITS file have extensions?           */
    status: *mut c_int,   /* IO - error status                            */
) -> c_int {
    unsafe {
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);

        let naxes = if naxes.is_null() {
            &[]
        } else {
            slice::from_raw_parts(naxes, naxis as usize)
        };

        ffphpr_safe(
            fptr, simple, bitpix, naxis, naxes, pcount, gcount, extend, status,
        )
    }
}

/*--------------------------------------------------------------------------*/
/// write required primary header keywords
pub fn ffphpr_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer                        */
    simple: c_int,       /* I - does file conform to FITS standard? 1/0  */
    bitpix: c_int,       /* I - number of bits per data value pixel      */
    naxis: c_int,        /* I - number of axes in the data array         */
    naxes: &[c_long],    /* I - length of each data axis                 */
    pcount: LONGLONG,    /* I - number of group parameters (usually 0)   */
    gcount: LONGLONG,    /* I - number of random groups (usually 1 or 0) */
    extend: c_int,       /* I - may FITS file have extensions?           */
    status: &mut c_int,  /* IO - error status                            */
) -> c_int {
    let mut ii: usize = 0;
    let mut naxesll = [0 as LONGLONG; 20];

    if naxis > 0 {
        while (ii < naxis as usize) && (ii < 20) {
            naxesll[ii] = naxes[ii] as LONGLONG;
            ii += 1;
        }
    }

    ffphprll_safe(
        fptr, simple, bitpix, naxis, &naxesll, pcount, gcount, extend, status,
    );

    *status
}

/*--------------------------------------------------------------------------*/
/// write required primary header keywords
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffphprll(
    fptr: *mut fitsfile,    /* I - FITS file pointer                        */
    simple: c_int,          /* I - does file conform to FITS standard? 1/0  */
    bitpix: c_int,          /* I - number of bits per data value pixel      */
    naxis: c_int,           /* I - number of axes in the data array         */
    naxes: *const LONGLONG, /* I - length of each data axis                 */
    pcount: LONGLONG,       /* I - number of group parameters (usually 0)   */
    gcount: LONGLONG,       /* I - number of random groups (usually 1 or 0) */
    extend: c_int,          /* I - may FITS file have extensions?           */
    status: *mut c_int,     /* IO - error status                            */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let naxes = if naxes.is_null() {
            &[]
        } else {
            slice::from_raw_parts(naxes, naxis as usize)
        };

        ffphprll_safe(
            fptr, simple, bitpix, naxis, naxes, pcount, gcount, extend, status,
        )
    }
}

/*--------------------------------------------------------------------------*/
/// write required primary header keywords
pub fn ffphprll_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer                        */
    simple: c_int,       /* I - does file conform to FITS standard? 1/0  */
    bitpix: c_int,       /* I - number of bits per data value pixel      */
    naxis: c_int,        /* I - number of axes in the data array         */
    naxes: &[LONGLONG],  /* I - length of each data axis                 */
    pcount: LONGLONG,    /* I - number of group parameters (usually 0)   */
    gcount: LONGLONG,    /* I - number of random groups (usually 1 or 0) */
    extend: c_int,       /* I - may FITS file have extensions?           */
    status: &mut c_int,  /* IO - error status                            */
) -> c_int {
    let mut ii: usize = 0;
    let mut tnaxes: [c_long; 20] = [0; 20];
    let mut name: [c_char; FLEN_KEYWORD] = [0; FLEN_KEYWORD];
    let mut comm: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];
    let mut message: [c_char; FLEN_ERRMSG] = [0; FLEN_ERRMSG];
    let mut card: [c_char; FLEN_CARD] = [0; FLEN_CARD];

    if *status > 0 {
        return *status;
    }

    if fptr.HDUposition != fptr.Fptr.curhdu {
        ffmahd_safe(fptr, (fptr.HDUposition) + 1, None, status);
    }

    let h = fptr.Fptr.get_headstart_as_slice();
    let headend = h[fptr.Fptr.curhdu as usize];

    if fptr.Fptr.headend != headend {
        *status = HEADER_NOT_EMPTY;
        return *status;
    }

    if naxis != 0 {
        /* never try to compress a null image */
        if fptr.Fptr.request_compress_type > 0 {
            ii = 0;
            while ii < naxis as usize {
                tnaxes[ii] = naxes[ii] as c_long;
                ii += 1;
            }
            /* write header for a compressed image */
            imcomp_init_table(fptr, bitpix, naxis, &tnaxes, true, status);

            return *status;
        }
    }

    if fptr.Fptr.curhdu == 0 {
        /* write primary array header */
        if simple > 0 {
            strcpy_safe(&mut comm, cs!(c"file does conform to FITS standard"));
        } else {
            strcpy_safe(&mut comm, cs!(c"file does not conform to FITS standard"));
        }
        ffpkyl_safe(fptr, cs!(c"SIMPLE"), simple, Some(&comm), status);
    } else {
        /* write IMAGE extension header */
        strcpy_safe(&mut comm, cs!(c"IMAGE extension"));
        ffpkys_safe(fptr, cs!(c"XTENSION"), cs!(c"IMAGE"), Some(&comm), status);
    }

    let mut longbitpix = bitpix;

    /* test for the 3 special cases that represent unsigned integers */
    if longbitpix == USHORT_IMG {
        longbitpix = SHORT_IMG;
    } else if longbitpix == ULONG_IMG {
        longbitpix = LONG_IMG;
    } else if longbitpix == ULONGLONG_IMG {
        longbitpix = LONGLONG_IMG;
    } else if longbitpix == SBYTE_IMG {
        longbitpix = BYTE_IMG;
    }
    if longbitpix != BYTE_IMG
        && longbitpix != SHORT_IMG
        && longbitpix != LONG_IMG
        && longbitpix != LONGLONG_IMG
        && longbitpix != FLOAT_IMG
        && longbitpix != DOUBLE_IMG
    {
        int_snprintf!(
            &mut message,
            FLEN_ERRMSG,
            "Illegal value for BITPIX keyword: {}",
            bitpix,
        );
        ffpmsg_slice(&message);
        *status = BAD_BITPIX;
        return *status;
    }

    strcpy_safe(&mut comm, cs!(c"number of bits per data pixel"));
    if ffpkyj_safe(
        fptr,
        cs!(c"BITPIX"),
        LONGLONG::from(longbitpix),
        Some(&comm),
        status,
    ) > 0
    {
        return *status;
    }
    if naxis < 0 || naxis > 999 {
        int_snprintf!(
            &mut message,
            FLEN_ERRMSG,
            "Illegal value for NAXIS keyword: {}",
            naxis,
        );
        ffpmsg_slice(&message);
        *status = BAD_NAXIS;
        return *status;
    }

    strcpy_safe(&mut comm, cs!(c"number of data axes"));
    ffpkyj_safe(
        fptr,
        cs!(c"NAXIS"),
        LONGLONG::from(naxis),
        Some(&comm),
        status,
    );

    strcpy_safe(&mut comm, cs!(c"length of data axis "));
    ii = 0;
    while ii < naxis as usize {
        if naxes[ii] < 0 {
            int_snprintf!(
                &mut message,
                FLEN_ERRMSG,
                "Illegal negative value for NAXIS{} keyword: {:.0}",
                ii + 1,
                naxes[ii] as f64,
            );
            ffpmsg_slice(&message);
            *status = BAD_NAXES;
            return *status;
        }

        int_snprintf!(&mut comm[20..], FLEN_COMMENT - 20, "{}", ii + 1,);
        ffkeyn_safe(cs!(c"NAXIS"), ii as c_int + 1, &mut name, status);
        ffpkyj_safe(fptr, &name, naxes[ii], Some(&comm), status);
        ii += 1;
    }

    if fptr.Fptr.curhdu == 0 {
        /* the primary array */

        if extend > 0 {
            /* only write EXTEND keyword if value = true */
            strcpy_safe(&mut comm, cs!(c"FITS dataset may contain extensions"));
            ffpkyl_safe(fptr, cs!(c"EXTEND"), extend, Some(&comm), status);
        }

        if pcount < 0 {
            ffpmsg_str("pcount value is less than 0");
            *status = BAD_PCOUNT;
            return *status;
        } else if gcount < 1 {
            ffpmsg_str("gcount value is less than 1");
            *status = BAD_GCOUNT;
            return *status;
        } else if pcount > 0 || gcount > 1 {
            /* only write these keyword if non-standard values */
            strcpy_safe(&mut comm, cs!(c"random group records are present"));
            ffpkyl_safe(fptr, cs!(c"GROUPS"), 1, Some(&comm), status);

            strcpy_safe(&mut comm, cs!(c"number of random group parameters"));
            ffpkyj_safe(fptr, cs!(c"PCOUNT"), pcount, Some(&comm), status);

            strcpy_safe(&mut comm, cs!(c"number of random groups"));
            ffpkyj_safe(fptr, cs!(c"GCOUNT"), gcount, Some(&comm), status);
        }

        /* write standard block of self-documentating comments */
        ffprec_safe(
            fptr,
            cs!(
                c"COMMENT   FITS (Flexible Image Transport System) format is defined in 'Astronomy"
            ),
            status,
        );
        ffprec_safe(
            fptr,
            cs!(c"COMMENT   and Astrophysics', volume 376, page 359; bibcode: 2001A&A...376..359H"),
            status,
        );
    } else
    /* an IMAGE extension */
    {
        /* image extension; cannot have random groups */
        if pcount != 0 {
            ffpmsg_str("image extensions must have pcount = 0");
            *status = BAD_PCOUNT;
        } else if gcount != 1 {
            ffpmsg_str("image extensions must have gcount = 1");
            *status = BAD_GCOUNT;
        } else {
            strcpy_safe(&mut comm, cs!(c"required keyword; must = 0"));
            ffpkyj_safe(fptr, cs!(c"PCOUNT"), 0, Some(&comm), status);

            strcpy_safe(&mut comm, cs!(c"required keyword; must = 1"));
            ffpkyj_safe(fptr, cs!(c"GCOUNT"), 1, Some(&comm), status);
        }
    }

    /* Write the BSCALE and BZERO keywords, if an unsigned integer image */
    if bitpix == USHORT_IMG {
        strcpy_safe(
            &mut comm,
            cs!(c"offset data range to that of unsigned short"),
        );
        ffpkyg_safe(fptr, cs!(c"BZERO"), 32768., 0, Some(&comm), status);
        strcpy_safe(&mut comm, cs!(c"default scaling factor"));
        ffpkyg_safe(fptr, cs!(c"BSCALE"), 1.0, 0, Some(&comm), status);
    } else if bitpix == ULONG_IMG {
        strcpy_safe(
            &mut comm,
            cs!(c"offset data range to that of unsigned long"),
        );
        ffpkyg_safe(fptr, cs!(c"BZERO"), 2147483648., 0, Some(&comm), status);
        strcpy_safe(&mut comm, cs!(c"default scaling factor"));
        ffpkyg_safe(fptr, cs!(c"BSCALE"), 1.0, 0, Some(&comm), status);
    } else if bitpix == ULONGLONG_IMG {
        strcpy_safe(
            &mut card,
            cs!(
                c"BZERO   =  9223372036854775808 / offset data range to that of unsigned long long"
            ),
        );
        ffprec_safe(fptr, &card, status);
        strcpy_safe(&mut comm, cs!(c"default scaling factor"));
        ffpkyg_safe(fptr, cs!(c"BSCALE"), 1.0, 0, Some(&comm), status);
    } else if bitpix == SBYTE_IMG {
        strcpy_safe(&mut comm, cs!(c"offset data range to that of signed byte"));
        ffpkyg_safe(fptr, cs!(c"BZERO"), -128., 0, Some(&comm), status);
        strcpy_safe(&mut comm, cs!(c"default scaling factor"));
        ffpkyg_safe(fptr, cs!(c"BSCALE"), 1.0, 0, Some(&comm), status);
    }
    *status
}

/*--------------------------------------------------------------------------*/
/// Put required Header keywords into the ASCII TaBle:
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffphtb(
    fptr: *mut fitsfile,         /* I - FITS file pointer                        */
    naxis1: LONGLONG,            /* I - width of row in the table                */
    naxis2: LONGLONG,            /* I - number of rows in the table              */
    tfields: c_int,              /* I - number of columns in the table           */
    ttype: *const *const c_char, /* I - name of each column                      */
    tbcol: *const c_long,        /* I - byte offset in row to each column        */
    tform: *const *const c_char, /* I - value of TFORMn keyword for each column  */
    tunit: *const *const c_char, /* I - value of TUNITn keyword for each column  */
    extnmx: *const c_char,       /* I - value of EXTNAME keyword, if any         */
    status: *mut c_int,          /* IO - error status                            */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        let tbcol = if tbcol.is_null() {
            None
        } else {
            Some(slice::from_raw_parts(tbcol, tfields as usize))
        };

        let tkeywords = TKeywords::new(tfields, ttype, tform, tunit);
        let (v_ttype, v_tform, v_tunit) = tkeywords.tkeywords_to_vecs();

        nullable_slice_cstr!(extnmx);

        ffphtb_safe(
            fptr,
            naxis1,
            naxis2,
            tfields,
            &v_ttype,
            tbcol,
            &v_tform,
            v_tunit.as_deref(),
            extnmx,
            status,
        )
    }
}

/*--------------------------------------------------------------------------*/
/// Put required Header keywords into the ASCII TaBle:
pub fn ffphtb_safe(
    fptr: &mut fitsfile,         /* I - FITS file pointer                        */
    naxis1: LONGLONG,            /* I - width of row in the table                */
    naxis2: LONGLONG,            /* I - number of rows in the table              */
    tfields: c_int,              /* I - number of columns in the table           */
    ttype: &[Option<&[c_char]>], /* I - name of each column                      */
    tbcol: Option<&[c_long]>,    /* I - byte offset in row to each column        */
    tform: &[&[c_char]],         /* I - value of TFORMn keyword for each column  */
    tunit: Option<&[Option<&[c_char]>]>, /* I - value of TUNITn keyword for each column  */
    extnmx: Option<&[c_char]>,   /* I - value of EXTNAME keyword, if any         */
    status: &mut c_int,          /* IO - error status                            */
) -> c_int {
    let mut ncols: c_int = 0;

    let mut rowlen: c_long = 0; /* must be 'long' because it is passed to ffgabc */
    let mut tfmt: [c_char; 30] = [0; 30];
    let mut name: [c_char; FLEN_KEYWORD] = [0; FLEN_KEYWORD];
    let mut comm: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];
    let mut extnm: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];

    let mut v = Vec::new();
    let tbcol_slice: &[c_long];

    if fptr.HDUposition != fptr.Fptr.curhdu {
        ffmahd_safe(fptr, (fptr.HDUposition) + 1, None, status);
    }

    let headstart = fptr.Fptr.get_headstart_as_slice();

    if *status > 0 {
        return *status;
    } else if fptr.Fptr.headend != headstart[fptr.Fptr.curhdu as usize] {
        *status = HEADER_NOT_EMPTY;
        return *status;
    } else if naxis1 < 0 {
        *status = NEG_WIDTH;
        return *status;
    } else if naxis2 < 0 {
        *status = NEG_ROWS;
        return *status;
    } else if tfields < 0 || tfields > 999 {
        *status = BAD_TFIELDS;
        return *status;
    }

    extnm[0] = 0;
    if let Some(extnmx) = extnmx {
        strncat_safe(&mut extnm, extnmx, FLEN_VALUE - 1);
    }
    rowlen = naxis1 as c_long;

    if tbcol.is_none() || (naxis1 == 0 && tfields != 0) {
        /* spacing not defined? */

        /* allocate mem for tbcol; malloc can have problems allocating small */
        /* arrays, so allocate at least 20 bytes */

        ncols = cmp::max(5, tfields);

        let tmp_size = mem::size_of::<c_long>() * ncols as usize;
        if v.try_reserve_exact(tmp_size).is_ok() {
            v.resize(tmp_size, 0);

            /* calculate width of a row and starting position of each column. */
            /* Each column will be separated by 1 blank space */
            ffgabc_safe(tfields, tform, 1, &mut rowlen, &mut v, status);
        }
        tbcol_slice = &v;
    } else if let Some(tbcol) = tbcol
        && !tbcol.is_empty()
        && tbcol[0] == 0
    {
        // REPEAT all above since can't chain

        /* spacing not defined? */

        /* allocate mem for tbcol; malloc can have problems allocating small */
        /* arrays, so allocate at least 20 bytes */

        ncols = cmp::max(5, tfields);

        let tmp_size = mem::size_of::<c_long>() * ncols as usize;
        if v.try_reserve_exact(tmp_size).is_ok() {
            v.resize(tmp_size, 0);

            /* calculate width of a row and starting position of each column. */
            /* Each column will be separated by 1 blank space */
            ffgabc_safe(tfields, tform, 1, &mut rowlen, &mut v, status);
        }
        tbcol_slice = &v;
    } else if let Some(tbcol) = tbcol {
        tbcol_slice = tbcol;
    } else {
        /* unreachable: the first arm above already handles tbcol == None */
        tbcol_slice = &v;
    }

    ffpkys_safe(
        fptr,
        cs!(c"XTENSION"),
        cs!(c"TABLE"),
        Some(cs!(c"ASCII table extension")),
        status,
    );
    ffpkyj_safe(
        fptr,
        cs!(c"BITPIX"),
        8,
        Some(cs!(c"8-bit ASCII characters")),
        status,
    );
    ffpkyj_safe(
        fptr,
        cs!(c"NAXIS"),
        2,
        Some(cs!(c"2-dimensional ASCII table")),
        status,
    );
    ffpkyj_safe(
        fptr,
        cs!(c"NAXIS1"),
        rowlen as LONGLONG,
        Some(cs!(c"width of table in characters")),
        status,
    );
    ffpkyj_safe(
        fptr,
        cs!(c"NAXIS2"),
        naxis2,
        Some(cs!(c"number of rows in table")),
        status,
    );
    ffpkyj_safe(
        fptr,
        cs!(c"PCOUNT"),
        0,
        Some(cs!(c"no group parameters (required keyword)")),
        status,
    );
    ffpkyj_safe(
        fptr,
        cs!(c"GCOUNT"),
        1,
        Some(cs!(c"one data group (required keyword)")),
        status,
    );
    ffpkyj_safe(
        fptr,
        cs!(c"TFIELDS"),
        tfields.into(),
        Some(cs!(c"number of fields in each row")),
        status,
    );

    for ii in 0..(tfields as usize) {
        /* loop over every column */
        //for (ii = 0; ii < tfields; ii++) {
        if let Some(ttype_item) = ttype[ii] {
            /* optional TTYPEn keyword */
            int_snprintf!(&mut comm, FLEN_COMMENT, "label for field {:3}", ii + 1,);

            ffkeyn_safe(cs!(c"TTYPE"), (ii + 1) as c_int, &mut name, status);
            ffpkys_safe(fptr, &name, ttype_item, Some(&comm), status);
        }

        if tbcol_slice[ii] < 1 || tbcol_slice[ii] > rowlen {
            *status = BAD_TBCOL;
        }

        int_snprintf!(
            &mut comm,
            FLEN_COMMENT,
            "beginning column of field {:3}",
            ii + 1,
        );
        ffkeyn_safe(cs!(c"TBCOL"), (ii + 1) as c_int, &mut name, status);
        ffpkyj_safe(
            fptr,
            &name,
            tbcol_slice[ii] as LONGLONG,
            Some(&comm),
            status,
        );

        if strlen_safe(tform[ii]) > 29 {
            ffpmsg_str("Error: ASCII table TFORM code is too long (ffphtb)");
            *status = BAD_TFORM;
            break;
        }
        strcpy_safe(&mut tfmt, tform[ii]); /* required TFORMn keyword */
        ffupch_safe(&mut tfmt);
        ffkeyn_safe(cs!(c"TFORM"), (ii + 1) as c_int, &mut name, status);
        ffpkys_safe(
            fptr,
            &name,
            &tfmt,
            Some(cs!(c"Fortran-77 format of field")),
            status,
        );

        if let Some(tunit) = tunit
            && let Some(tunit_item) = tunit[ii]
            && tunit_item[0] != 0
        {
            /* optional TUNITn keyword */

            ffkeyn_safe(cs!(c"TUNIT"), (ii + 1) as c_int, &mut name, status);
            ffpkys_safe(
                fptr,
                &name,
                tunit_item,
                Some(cs!(c"physical unit of field")),
                status,
            );
        }

        if *status > 0 {
            break;
        } /* abort loop on error */
    }

    if extnm[0] > 0 {
        /* optional EXTNAME keyword */
        ffpkys_safe(
            fptr,
            cs!(c"EXTNAME"),
            &extnm,
            Some(cs!(c"name of this ASCII table extension")),
            status,
        );
    }
    if *status > 0 {
        ffpmsg_str("Failed to write ASCII table header keywords (ffphtb)");
    }

    *status
}

/*--------------------------------------------------------------------------*/
/// Put required Header keywords into the Binary Table:
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffphbn(
    fptr: *mut fitsfile,         /* I - FITS file pointer                        */
    naxis2: LONGLONG,            /* I - number of rows in the table              */
    tfields: c_int,              /* I - number of columns in the table           */
    ttype: *const *const c_char, /* I - name of each column                      */
    tform: *const *const c_char, /* I - value of TFORMn keyword for each column  */
    tunit: *const *const c_char, /* I - value of TUNITn keyword for each column  */
    extnmx: *const c_char,       /* I - value of EXTNAME keyword, if any         */
    pcount: LONGLONG,            /* I - size of the variable length heap area    */
    status: *mut c_int,          /* IO - error status                            */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        let tkeywords = TKeywords::new(tfields, ttype, tform, tunit);
        let (v_ttype, v_tform, v_tunit) = tkeywords.tkeywords_to_vecs();

        nullable_slice_cstr!(extnmx);

        ffphbn_safe(
            fptr,
            naxis2,
            tfields,
            &v_ttype,
            &v_tform,
            v_tunit.as_deref(),
            extnmx,
            pcount,
            status,
        )
    }
}

/*--------------------------------------------------------------------------*/
/// Put required Header keywords into the Binary Table:
#[allow(clippy::if_same_then_else)]
// C dispatch chain: distinct conditions deliberately share an action.
pub fn ffphbn_safe(
    fptr: &mut fitsfile,         /* I - FITS file pointer                        */
    naxis2: LONGLONG,            /* I - number of rows in the table              */
    tfields: c_int,              /* I - number of columns in the table           */
    ttype: &[Option<&[c_char]>], /* I - name of each column                      */
    tform: &[&[c_char]],         /* I - value of TFORMn keyword for each column  */
    tunit: Option<&[Option<&[c_char]>]>, /* I - value of TUNITn keyword for each column  */
    extnmx: Option<&[c_char]>,   /* I - value of EXTNAME keyword, if any         */
    pcount: LONGLONG,            /* I - size of the variable length heap area    */
    status: &mut c_int,          /* IO - error status                            */
) -> c_int {
    let mut datatype: c_int = 0;
    let mut iread: c_int = 0;
    let mut repeat: c_long = 0;
    let mut width: c_long = 0;
    let mut naxis1: LONGLONG = 0;
    let mut tfmt: [c_char; 30] = [0; 30];
    let mut name: [c_char; FLEN_KEYWORD] = [0; FLEN_KEYWORD];
    let mut comm: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];
    let mut extnm: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut card: [c_char; FLEN_CARD] = [0; FLEN_CARD];

    if *status > 0 {
        return *status;
    }

    if fptr.HDUposition != fptr.Fptr.curhdu {
        ffmahd_safe(fptr, (fptr.HDUposition) + 1, None, status);
    }

    let h = fptr.Fptr.get_headstart_as_slice();

    if fptr.Fptr.headend != h[fptr.Fptr.curhdu as usize] {
        *status = HEADER_NOT_EMPTY;
        return *status;
    } else if naxis2 < 0 {
        *status = NEG_ROWS;
        return *status;
    } else if pcount < 0 {
        *status = BAD_PCOUNT;
        return *status;
    } else if tfields < 0 || tfields > 999 {
        *status = BAD_TFIELDS;
        return *status;
    }

    extnm[0] = 0;

    if let Some(extnmx) = extnmx {
        strncat_safe(&mut extnm, extnmx, FLEN_VALUE - 1);
    }
    ffpkys_safe(
        fptr,
        cs!(c"XTENSION"),
        cs!(c"BINTABLE"),
        Some(cs!(c"binary table extension")),
        status,
    );
    ffpkyj_safe(fptr, cs!(c"BITPIX"), 8, Some(cs!(c"8-bit bytes")), status);
    ffpkyj_safe(
        fptr,
        cs!(c"NAXIS"),
        2,
        Some(cs!(c"2-dimensional binary table")),
        status,
    );

    naxis1 = 0;
    for ii in 0..(tfields as usize) {
        /* sum the width of each field */
        let tform_item: &[c_char] = tform[ii];

        ffbnfm_safe(
            tform_item,
            Some(&mut datatype),
            Some(&mut repeat),
            Some(&mut width),
            status,
        );

        if datatype == TSTRING {
            naxis1 += repeat as LONGLONG; /* one byte per char */
        } else if datatype == TBIT {
            naxis1 += ((repeat + 7) / 8) as LONGLONG;
        } else if datatype > 0 {
            naxis1 += repeat as LONGLONG * (LONGLONG::from(datatype) / 10);
        } else if tform_item[0] == bb(b'P')
            || tform_item[1] == bb(b'P')
            || tform_item[0] == bb(b'p')
            || tform_item[1] == bb(b'p')
        {
            /* this is a 'P' variable length descriptor (neg. datatype) */
            naxis1 += 8;
        } else {
            /* this is a 'Q' variable length descriptor (neg. datatype) */
            naxis1 += 16;
        }
        if *status > 0 {
            break; /* abort loop on error */
        }
    }

    ffpkyj_safe(
        fptr,
        cs!(c"NAXIS1"),
        naxis1,
        Some(cs!(c"width of table in bytes")),
        status,
    );
    ffpkyj_safe(
        fptr,
        cs!(c"NAXIS2"),
        naxis2,
        Some(cs!(c"number of rows in table")),
        status,
    );

    /*
    the initial value of PCOUNT (= size of the variable length array heap)
    should always be zero.  If any variable length data is written, then
    the value of PCOUNT will be updated when the HDU is closed
    */
    ffpkyj_safe(
        fptr,
        cs!(c"PCOUNT"),
        0,
        Some(cs!(c"size of special data area")),
        status,
    );
    ffpkyj_safe(
        fptr,
        cs!(c"GCOUNT"),
        1,
        Some(cs!(c"one data group (required keyword)")),
        status,
    );
    ffpkyj_safe(
        fptr,
        cs!(c"TFIELDS"),
        LONGLONG::from(tfields),
        Some(cs!(c"number of fields in each row")),
        status,
    );

    /* sum the width of each field */
    for ii in 0..(tfields as usize) {
        let tform_item: &[c_char] = tform[ii];

        /* loop over every column */
        if let Some(ttype_item) = ttype[ii] {
            /* optional TTYPEn keyword */

            int_snprintf!(&mut comm, FLEN_COMMENT, "label for field {:3}", ii + 1,);
            ffkeyn_safe(cs!(c"TTYPE"), (ii + 1) as c_int, &mut name, status);
            ffpkys_safe(fptr, &name, ttype_item, Some(&comm), status);
        }

        if strlen_safe(tform_item) > 29 {
            ffpmsg_str("Error: BIN table TFORM code is too long (ffphbn)");
            *status = BAD_TFORM;
            break;
        }
        strcpy_safe(&mut tfmt, tform_item); /* required TFORMn keyword */
        ffupch_safe(&mut tfmt);

        ffkeyn_safe(cs!(c"TFORM"), (ii + 1) as c_int, &mut name, status);
        strcpy_safe(&mut comm, cs!(c"data format of field"));

        ffbnfm_safe(
            &tfmt,
            Some(&mut datatype),
            Some(&mut repeat),
            Some(&mut width),
            status,
        );

        if datatype == TSTRING {
            strcat_safe(&mut comm, cs!(c": ASCII Character"));

            /* Do sanity check to see if an ASCII table format was used,  */
            /* e.g., 'A8' instead of '8A', or a bad unit width eg '8A9'.  */
            /* Don't want to return an error status, so write error into  */
            /* the keyword comment.  */

            let cptr = strchr_safe(&tfmt, bb(b'A'));

            if let Some(cptr) = cptr {
                let c = cptr + 1;

                // iread = sscanf_ld(&tfmt[c..], cs!(c"%ld"), &mut width);
                let tmp: Result<c_long, ParseIntError> =
                    atoi(str::from_utf8(cast_slice(&tfmt[c..])).unwrap());

                match tmp {
                    Ok(x) => {
                        width = x;
                        iread = 1;
                    }
                    Err(_) => {
                        iread = 0;
                    }
                }
            }

            if iread == 1 && (width > repeat) {
                if repeat == 1 {
                    strcpy_safe(
                        &mut comm,
                        cs!(c"ERROR??  USING ASCII TABLE SYNTAX BY MISTAKE??"),
                    );
                } else {
                    strcpy_safe(
                        &mut comm,
                        cs!(c"rAw FORMAT ERROR! UNIT WIDTH w > COLUMN WIDTH r"),
                    );
                }
            }
        } else if datatype == TBIT {
            strcat_safe(&mut comm, cs!(c": BIT"));
        } else if datatype == TBYTE {
            strcat_safe(&mut comm, cs!(c": BYTE"));
        } else if datatype == TLOGICAL {
            strcat_safe(&mut comm, cs!(c": 1-byte LOGICAL"));
        } else if datatype == TSHORT {
            strcat_safe(&mut comm, cs!(c": 2-byte INTEGER"));
        } else if datatype == TUSHORT {
            strcat_safe(&mut comm, cs!(c": 2-byte INTEGER"));
        } else if datatype == TLONG {
            strcat_safe(&mut comm, cs!(c": 4-byte INTEGER"));
        } else if datatype == TLONGLONG {
            strcat_safe(&mut comm, cs!(c": 8-byte INTEGER"));
        } else if datatype == TULONG {
            strcat_safe(&mut comm, cs!(c": 4-byte INTEGER"));
        } else if datatype == TULONGLONG {
            strcat_safe(&mut comm, cs!(c": 8-byte INTEGER"));
        } else if datatype == TFLOAT {
            strcat_safe(&mut comm, cs!(c": 4-byte REAL"));
        } else if datatype == TDOUBLE {
            strcat_safe(&mut comm, cs!(c": 8-byte DOUBLE"));
        } else if datatype == TCOMPLEX {
            strcat_safe(&mut comm, cs!(c": COMPLEX"));
        } else if datatype == TDBLCOMPLEX {
            strcat_safe(&mut comm, cs!(c": DOUBLE COMPLEX"));
        } else if datatype < 0 {
            strcat_safe(&mut comm, cs!(c": variable length array"));
        }
        if datatype.abs() == TSBYTE
        /* signed bytes */
        {
            /* Replace the 'S' with an 'B' in the TFORMn code */
            let mut ci = 0;
            while tfmt[ci] != bb(b'S') {
                ci += 1;
            }
            tfmt[ci] = bb(b'B');
            ffpkys_safe(fptr, &name, &tfmt, Some(&comm), status);

            /* write the TZEROn and TSCALn keywords */
            ffkeyn_safe(cs!(c"TZERO"), (ii + 1) as c_int, &mut name, status);
            strcpy_safe(&mut comm, cs!(c"offset for signed bytes"));

            ffpkyg_safe(fptr, &name, -128., 0, Some(&comm), status);

            ffkeyn_safe(cs!(c"TSCAL"), (ii + 1) as c_int, &mut name, status);
            strcpy_safe(&mut comm, cs!(c"data are not scaled"));
            ffpkyg_safe(fptr, &name, 1.0, 0, Some(&comm), status);
        } else if datatype.abs() == TUSHORT {
            /* Replace the 'U' with an 'I' in the TFORMn code */
            let mut ci = 0;
            while tfmt[ci] != bb(b'U') {
                ci += 1;
            }
            tfmt[ci] = bb(b'I');
            ffpkys_safe(fptr, &name, &tfmt, Some(&comm), status);

            /* write the TZEROn and TSCALn keywords */
            ffkeyn_safe(cs!(c"TZERO"), (ii + 1) as c_int, &mut name, status);
            strcpy_safe(&mut comm, cs!(c"offset for unsigned integers"));

            ffpkyg_safe(fptr, &name, 32768., 0, Some(&comm), status);

            ffkeyn_safe(cs!(c"TSCAL"), (ii + 1) as c_int, &mut name, status);
            strcpy_safe(&mut comm, cs!(c"data are not scaled"));
            ffpkyg_safe(fptr, &name, 1.0, 0, Some(&comm), status);
        } else if datatype.abs() == TULONG {
            /* Replace the 'V' with an 'J' in the TFORMn code */
            let mut ci = 0;
            while tfmt[ci] != bb(b'V') {
                ci += 1;
            }
            tfmt[ci] = bb(b'J');
            ffpkys_safe(fptr, &name, &tfmt, Some(&comm), status);

            /* write the TZEROn and TSCALn keywords */
            ffkeyn_safe(cs!(c"TZERO"), (ii + 1) as c_int, &mut name, status);
            strcpy_safe(&mut comm, cs!(c"offset for unsigned integers"));

            ffpkyg_safe(fptr, &name, 2147483648., 0, Some(&comm), status);

            ffkeyn_safe(cs!(c"TSCAL"), (ii + 1) as c_int, &mut name, status);
            strcpy_safe(&mut comm, cs!(c"data are not scaled"));
            ffpkyg_safe(fptr, &name, 1.0, 0, Some(&comm), status);
        } else if datatype.abs() == TULONGLONG {
            /* Replace the 'W' with an 'K' in the TFORMn code */
            let mut ci = 0;
            while tfmt[ci] != bb(b'W') {
                ci += 1;
            }
            tfmt[ci] = bb(b'K');
            ffpkys_safe(fptr, &name, &tfmt, Some(&comm), status);

            /* write the TZEROn and TSCALn keywords */
            ffkeyn_safe(cs!(c"TZERO"), ii as c_int + 1, &mut card, status);
            strcat_safe(&mut card, cs!(c"     ")); /* make sure name is >= 8 chars long */

            card[8] = 0;
            strcat_safe(
                &mut card,
                cs!(c"=  9223372036854775808 / offset for unsigned integers"),
            );
            ffprec_safe(fptr, &card, status);

            ffkeyn_safe(cs!(c"TSCAL"), (ii + 1) as c_int, &mut name, status);
            strcpy_safe(&mut comm, cs!(c"data are not scaled"));
            ffpkyg_safe(fptr, &name, 1., 0, Some(&comm), status);
        } else {
            ffpkys_safe(fptr, &name, &tfmt, Some(&comm), status);
        }

        if let Some(tunit) = tunit
            && let Some(tunit_item) = tunit[ii]
            && tunit_item[0] != 0
        {
            /* optional TUNITn keyword */

            ffkeyn_safe(cs!(c"TUNIT"), (ii + 1) as c_int, &mut name, status);
            ffpkys_safe(
                fptr,
                &name,
                tunit_item,
                Some(cs!(c"physical unit of field")),
                status,
            );
        }

        if *status > 0 {
            break; /* abort loop on error */
        }
    }

    if extnm[0] > 0 {
        /* optional EXTNAME keyword */
        ffpkys_safe(
            fptr,
            cs!(c"EXTNAME"),
            &extnm,
            Some(cs!(c"name of this binary table extension")),
            status,
        );
    }
    if *status > 0 {
        ffpmsg_str("Failed to write binary table header keywords (ffphbn)");
    }
    *status
}

/*--------------------------------------------------------------------------*/
/// Put required Header keywords into a conforming extension:
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffphext(
    fptr: *mut fitsfile,      /* I - FITS file pointer                       */
    xtensionx: *const c_char, /* I - value for the XTENSION keyword          */
    bitpix: c_int,            /* I - value for the BIXPIX keyword            */
    naxis: c_int,             /* I - value for the NAXIS keyword             */
    naxes: *const c_long,     /* I - value for the NAXISn keywords           */
    pcount: LONGLONG,         /* I - value for the PCOUNT keyword            */
    gcount: LONGLONG,         /* I - value for the GCOUNT keyword            */
    status: *mut c_int,       /* IO - error status                           */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        raw_to_slice!(xtensionx);
        let naxes = slice::from_raw_parts(naxes, naxis as usize);

        ffphext_safe(
            fptr, xtensionx, bitpix, naxis, naxes, pcount, gcount, status,
        )
    }
}

/*--------------------------------------------------------------------------*/
/// Put required Header keywords into a conforming extension:
pub fn ffphext_safe(
    fptr: &mut fitsfile,  /* I - FITS file pointer                       */
    xtensionx: &[c_char], /* I - value for the XTENSION keyword          */
    bitpix: c_int,        /* I - value for the BIXPIX keyword            */
    naxis: c_int,         /* I - value for the NAXIS keyword             */
    naxes: &[c_long],     /* I - value for the NAXISn keywords           */
    pcount: LONGLONG,     /* I - value for the PCOUNT keyword            */
    gcount: LONGLONG,     /* I - value for the GCOUNT keyword            */
    status: &mut c_int,   /* IO - error status                           */
) -> c_int {
    let mut message: [c_char; FLEN_ERRMSG] = [0; FLEN_ERRMSG];
    let mut comm: [c_char; 81] = [0; 81];
    let mut name: [c_char; 20] = [0; 20];
    let mut xtension: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];

    if fptr.HDUposition != fptr.Fptr.curhdu {
        ffmahd_safe(fptr, (fptr.HDUposition) + 1, None, status);
    }

    if *status > 0 {
        return *status;
    } else {
        let h = fptr.Fptr.get_headstart_as_slice();
        if fptr.Fptr.headend != h[fptr.Fptr.curhdu as usize] {
            *status = HEADER_NOT_EMPTY;
            return *status;
        }
    }

    if naxis < 0 || naxis > 999 {
        int_snprintf!(
            &mut message,
            FLEN_ERRMSG,
            "Illegal value for NAXIS keyword: {}",
            naxis
        );
        ffpmsg_slice(&message);
        *status = BAD_NAXIS;
        return *status;
    }

    xtension[0] = 0;
    strncat_safe(&mut xtension, xtensionx, FLEN_VALUE - 1);

    ffpkys_safe(
        fptr,
        cs!(c"XTENSION"),
        &xtension,
        Some(cs!(c"extension type")),
        status,
    );
    ffpkyj_safe(
        fptr,
        cs!(c"BITPIX"),
        bitpix as LONGLONG,
        Some(cs!(c"number of bits per data pixel")),
        status,
    );
    ffpkyj_safe(
        fptr,
        cs!(c"NAXIS"),
        naxis as LONGLONG,
        Some(cs!(c"number of data axes")),
        status,
    );

    strcpy_safe(&mut comm, cs!(c"length of data axis "));
    for ii in 0..(naxis as usize) {
        if naxes[ii] < 0 {
            int_snprintf!(
                &mut message,
                FLEN_ERRMSG,
                "Illegal negative value for NAXIS{} keyword: {}",
                ii + 1,
                naxes[ii]
            );
            ffpmsg_slice(&message);
            *status = BAD_NAXES;
            return *status;
        }

        int_snprintf!(&mut comm[20..], 61, "{}", ii + 1);
        ffkeyn_safe(cs!(c"NAXIS"), (ii + 1) as c_int, &mut name, status);
        ffpkyj_safe(fptr, &name, naxes[ii] as LONGLONG, Some(&comm), status);
    }

    ffpkyj_safe(fptr, cs!(c"PCOUNT"), pcount, Some(cs!(c" ")), status);
    ffpkyj_safe(fptr, cs!(c"GCOUNT"), gcount, Some(cs!(c" ")), status);

    if *status > 0 {
        ffpmsg_str("Failed to write extension header keywords (ffphext)");
    }

    *status
}

/*-------------------------------------------------------------------------*/
/// write a keyword record (80 bytes long) to the end of the header
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffprec(
    fptr: *mut fitsfile, /* I - FITS file pointer        */
    card: *const c_char, /* I - string to be written     */
    status: *mut c_int,  /* IO - error status            */
) -> c_int {
    unsafe {
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);

        raw_to_slice!(card);

        ffprec_safe(fptr, card, status)
    }
}

/*-------------------------------------------------------------------------*/
/// write a keyword record (80 bytes long) to the end of the header
pub fn ffprec_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer        */
    card: &[c_char],     /* I - string to be written     */
    status: &mut c_int,  /* IO - error status            */
) -> c_int {
    let mut tcard: [c_char; FLEN_CARD] = [0; FLEN_CARD];

    let mut ii: usize;
    let nblocks: c_long;
    let mut keylength: c_int;

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    if fptr.HDUposition != fptr.Fptr.curhdu {
        ffmahd_safe(fptr, (fptr.HDUposition) + 1, None, status);
    }

    if (fptr.Fptr.datastart - fptr.Fptr.headend) == 80 {
        /* no room */
        nblocks = 1;
        if ffiblk(fptr, nblocks, 0, status) > 0 {
            /* insert 2880-byte block */
            return *status;
        }
    }

    strncpy_safe(&mut tcard, card, 80);
    tcard[80] = 0;

    let len: usize = strlen_safe(&tcard);

    /* silently replace any illegal characters with a space */
    ii = 0;
    while ii < len {
        if tcard[ii] < b' ' as c_char || tcard[ii] > 126 {
            tcard[ii] = b' ' as c_char;
        }
        ii += 1;
    }
    ii = len;
    while ii < 80 {
        /* fill card with spaces if necessary */
        tcard[ii] = b' ' as c_char;
        ii += 1;
    }
    keylength = strcspn_safe(&tcard, cs!(c"=")) as c_int; /* support for free-format keywords */
    if keylength == 80 {
        keylength = 8;
    }
    /* test for the common commentary keywords which by definition have 8-char names */
    if fits_strncasecmp(cs!(c"COMMENT "), &tcard, 8) == 0
        || fits_strncasecmp(cs!(c"HISTORY "), &tcard, 8) == 0
        || fits_strncasecmp(cs!(c"        "), &tcard, 8) == 0
        || fits_strncasecmp(cs!(c"CONTINUE"), &tcard, 8) == 0
    {
        keylength = 8;
    }
    ii = 0;
    while ii < keylength as usize {
        /* make sure keyword name is uppercase */
        tcard[ii] = toupper(tcard[ii]);
        ii += 1;
    }

    fftkey_safe(&tcard, status); /* test keyword name contains legal chars */

    /*  no need to do this any more, since any illegal characters have been removed
    fftrec(tcard, status);  */
    /* test rest of keyword for legal chars */

    ffmbyt_safe(fptr, fptr.Fptr.headend, IGNORE_EOF, status); /* move to end */

    ffpbyt(fptr, 80, cast_slice(&tcard), status); /* write the 80 byte card */

    if *status <= 0 {
        fptr.Fptr.headend += 80; /* update end-of-header position */
    }
    *status
}

/*--------------------------------------------------------------------------*/
/// Write (put) a null-valued keyword and comment into the FITS header.  
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffpkyu(
    fptr: *mut fitsfile,    /* I - FITS file pointer        */
    keyname: *const c_char, /* I - name of keyword to write */
    comm: *const c_char,    /* I - keyword comment          */
    status: *mut c_int,     /* IO - error status            */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        raw_to_slice!(keyname);
        nullable_slice_cstr!(comm);

        ffpkyu_safe(fptr, keyname, comm, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Write (put) a null-valued keyword and comment into the FITS header.  
pub fn ffpkyu_safe(
    fptr: &mut fitsfile,     /* I - FITS file pointer        */
    keyname: &[c_char],      /* I - name of keyword to write */
    comm: Option<&[c_char]>, /* I - keyword comment          */
    status: &mut c_int,      /* IO - error status            */
) -> c_int {
    let mut valstring: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut card: [c_char; FLEN_CARD] = [0; FLEN_CARD];

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    strcpy_safe(&mut valstring, cs!(c" ")); /* create a dummy value string */
    ffmkky_safe(keyname, &valstring, comm, &mut card, status); /* construct the keyword */
    ffprec_safe(fptr, &card, status);

    *status
}

/*--------------------------------------------------------------------------*/
/// Write (put) the keyword, value and comment into the FITS header.
/// The value string will be truncated at 68 characters which is the
/// maximum length that will fit on a single FITS keyword.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffpkys(
    fptr: *mut fitsfile,    /* I - FITS file pointer        */
    keyname: *const c_char, /* I - name of keyword to write */
    value: *const c_char,   /* I - keyword value            */
    comm: *const c_char,    /* I - keyword comment          */
    status: *mut c_int,     /* IO - error status            */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);
        raw_to_slice!(keyname);
        raw_to_slice!(value);

        let comm: Option<&[c_char]> = match comm.is_null() {
            true => None,
            false => Some(cast_slice(CStr::from_ptr(comm).to_bytes_with_nul())),
        };

        ffpkys_safe(fptr, keyname, value, comm, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Write (put) the keyword, value and comment into the FITS header.
/// The value string will be truncated at 68 characters which is the
/// maximum length that will fit on a single FITS keyword.
pub fn ffpkys_safe(
    fptr: &mut fitsfile,     /* I - FITS file pointer        */
    keyname: &[c_char],      /* I - name of keyword to write */
    value: &[c_char],        /* I - keyword value            */
    comm: Option<&[c_char]>, /* I - keyword comment          */
    status: &mut c_int,      /* IO - error status            */
) -> c_int {
    let mut valstring: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut card: [c_char; FLEN_CARD] = [0; FLEN_CARD];

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    ffs2c(value, &mut valstring, status); /* put quotes around the string */

    ffmkky_safe(keyname, &valstring, comm, &mut card, status); /* construct the keyword */
    ffprec_safe(fptr, &card, status);

    *status
}

/*--------------------------------------------------------------------------*/
/// Write (put) the keyword, value and comment into the FITS header.
/// Writes an integer keyword value.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffpkyuj(
    fptr: *mut fitsfile,    /* I - FITS file pointer        */
    keyname: *const c_char, /* I - name of keyword to write */
    value: ULONGLONG,       /* I - keyword value            */
    comm: *const c_char,    /* I - keyword comment          */
    status: *mut c_int,     /* IO - error status            */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        raw_to_slice!(keyname);
        nullable_slice_cstr!(comm);

        ffpkyuj_safe(fptr, keyname, value, comm, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Write (put) the keyword, value and comment into the FITS header.
/// Writes an integer keyword value.
pub fn ffpkyuj_safe(
    fptr: &mut fitsfile,     /* I - FITS file pointer        */
    keyname: &[c_char],      /* I - name of keyword to write */
    value: ULONGLONG,        /* I - keyword value            */
    comm: Option<&[c_char]>, /* I - keyword comment          */
    status: &mut c_int,      /* IO - error status            */
) -> c_int {
    let mut valstring: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut card: [c_char; FLEN_CARD] = [0; FLEN_CARD];

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    ffu2c(value, &mut valstring, status); /* convert to formatted string */
    ffmkky_safe(keyname, &valstring, comm, &mut card, status); /* construct the keyword*/
    ffprec_safe(fptr, &card, status); /* write the keyword*/

    *status
}

/*--------------------------------------------------------------------------*/
/// Write (put) the keyword, value and comment into the FITS header.
/// Writes a fixed float keyword value.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffpkyf(
    fptr: *mut fitsfile,    /* I - FITS file pointer                   */
    keyname: *const c_char, /* I - name of keyword to write            */
    value: f32,             /* I - keyword value                       */
    decim: c_int,           /* I - number of decimal places to display */
    comm: *const c_char,    /* I - keyword comment                     */
    status: *mut c_int,     /* IO - error status                       */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        raw_to_slice!(keyname);
        nullable_slice_cstr!(comm);

        ffpkyf_safe(fptr, keyname, value, decim, comm, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Write (put) the keyword, value and comment into the FITS header.
/// Writes a fixed float keyword value.
pub fn ffpkyf_safe(
    fptr: &mut fitsfile,     /* I - FITS file pointer                   */
    keyname: &[c_char],      /* I - name of keyword to write            */
    value: f32,              /* I - keyword value                       */
    decim: c_int,            /* I - number of decimal places to display */
    comm: Option<&[c_char]>, /* I - keyword comment                     */
    status: &mut c_int,      /* IO - error status                       */
) -> c_int {
    let mut valstring: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut card: [c_char; FLEN_CARD] = [0; FLEN_CARD];

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    ffr2f(value, decim, &mut valstring, status); /* convert to formatted string */
    ffmkky_safe(keyname, &valstring, comm, &mut card, status); /* construct the keyword*/
    ffprec_safe(fptr, &card, status); /* write the keyword*/

    *status
}

/*--------------------------------------------------------------------------*/
/// Write (put) the keyword, value and comment into the FITS header.
/// Writes an exponential float keyword value.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffpkye(
    fptr: *mut fitsfile,    /* I - FITS file pointer                   */
    keyname: *const c_char, /* I - name of keyword to write            */
    value: f32,             /* I - keyword value                       */
    decim: c_int,           /* I - number of decimal places to display */
    comm: *const c_char,    /* I - keyword comment                     */
    status: *mut c_int,     /* IO - error status                       */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        raw_to_slice!(keyname);
        nullable_slice_cstr!(comm);

        ffpkye_safe(fptr, keyname, value, decim, comm, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Write (put) the keyword, value and comment into the FITS header.
/// Writes an exponential float keyword value.
pub fn ffpkye_safe(
    fptr: &mut fitsfile,     /* I - FITS file pointer                   */
    keyname: &[c_char],      /* I - name of keyword to write            */
    value: f32,              /* I - keyword value                       */
    decim: c_int,            /* I - number of decimal places to display */
    comm: Option<&[c_char]>, /* I - keyword comment                     */
    status: &mut c_int,      /* IO - error status                       */
) -> c_int {
    let mut valstring: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut card: [c_char; FLEN_CARD] = [0; FLEN_CARD];

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    ffr2e(value, decim, &mut valstring, status); /* convert to formatted string */
    ffmkky_safe(keyname, &valstring, comm, &mut card, status); /* construct the keyword*/
    ffprec_safe(fptr, &card, status); /* write the keyword*/

    *status
}

/*--------------------------------------------------------------------------*/
/// Write (put) the keyword, value and comment into the FITS header.
/// Writes a fixed double keyword value.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffpkyg(
    fptr: *mut fitsfile,    /* I - FITS file pointer                   */
    keyname: *const c_char, /* I - name of keyword to write            */
    value: f64,             /* I - keyword value                       */
    decim: c_int,           /* I - number of decimal places to display */
    comm: *const c_char,    /* I - keyword comment                     */
    status: *mut c_int,     /* IO - error status                       */
) -> c_int {
    unsafe {
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);

        raw_to_slice!(keyname);
        nullable_slice_cstr!(comm);

        ffpkyg_safe(fptr, keyname, value, decim, comm, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Write (put) the keyword, value and comment into the FITS header.
/// Writes a fixed double keyword value.
pub fn ffpkyg_safe(
    fptr: &mut fitsfile,     /* I - FITS file pointer                   */
    keyname: &[c_char],      /* I - name of keyword to write            */
    value: f64,              /* I - keyword value                       */
    decim: c_int,            /* I - number of decimal places to display */
    comm: Option<&[c_char]>, /* I - keyword comment                     */
    status: &mut c_int,      /* IO - error status                       */
) -> c_int {
    let mut valstring: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut card: [c_char; FLEN_CARD] = [0; FLEN_CARD];

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    ffd2f(value, decim, &mut valstring, status); /* convert to formatted string */
    ffmkky_safe(keyname, &valstring, comm, &mut card, status); /* construct the keyword*/
    ffprec_safe(fptr, &card, status); /* write the keyword*/

    *status
}

/*--------------------------------------------------------------------------*/
/// Write (put) the keyword, value and comment into the FITS header.
/// Writes an exponential double keyword value.*/
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffpkyd(
    fptr: *mut fitsfile,    /* I - FITS file pointer                   */
    keyname: *const c_char, /* I - name of keyword to write            */
    value: f64,             /* I - keyword value                       */
    decim: c_int,           /* I - number of decimal places to display */
    comm: *const c_char,    /* I - keyword comment                     */
    status: *mut c_int,     /* IO - error status                       */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        raw_to_slice!(keyname);
        nullable_slice_cstr!(comm);

        ffpkyd_safe(fptr, keyname, value, decim, comm, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Write (put) the keyword, value and comment into the FITS header.
/// Writes an exponential double keyword value.*/
pub fn ffpkyd_safe(
    fptr: &mut fitsfile,     /* I - FITS file pointer                   */
    keyname: &[c_char],      /* I - name of keyword to write            */
    value: f64,              /* I - keyword value                       */
    decim: c_int,            /* I - number of decimal places to display */
    comm: Option<&[c_char]>, /* I - keyword comment                     */
    status: &mut c_int,      /* IO - error status                       */
) -> c_int {
    let mut valstring: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut card: [c_char; FLEN_CARD] = [0; FLEN_CARD];

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    ffd2e(value, decim, &mut valstring, status); /* convert to formatted string */
    ffmkky_safe(keyname, &valstring, comm, &mut card, status); /* construct the keyword*/
    ffprec_safe(fptr, &card, status); /* write the keyword*/

    *status
}

/*--------------------------------------------------------------------------*/
/// Write (put) the keyword, value and comment into the FITS header.
/// Writes an complex float keyword value. Format = (realvalue, imagvalue)
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffpkyc(
    fptr: *mut fitsfile,    /* I - FITS file pointer                   */
    keyname: *const c_char, /* I - name of keyword to write            */
    value: *const [f32; 2], /* I - keyword value (real, imaginary)     */
    decim: c_int,           /* I - number of decimal places to display */
    comm: *const c_char,    /* I - keyword comment                     */
    status: *mut c_int,     /* IO - error status                       */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let value = value.as_ref().expect(NULL_MSG);

        raw_to_slice!(keyname);
        nullable_slice_cstr!(comm);

        ffpkyc_safe(fptr, keyname, value, decim, comm, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Write (put) the keyword, value and comment into the FITS header.
/// Writes an complex float keyword value. Format = (realvalue, imagvalue)
pub fn ffpkyc_safe(
    fptr: &mut fitsfile,     /* I - FITS file pointer                   */
    keyname: &[c_char],      /* I - name of keyword to write            */
    value: &[f32; 2],        /* I - keyword value (real, imaginary)     */
    decim: c_int,            /* I - number of decimal places to display */
    comm: Option<&[c_char]>, /* I - keyword comment                     */
    status: &mut c_int,      /* IO - error status                       */
) -> c_int {
    let mut valstring: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut tmpstring: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut card: [c_char; FLEN_CARD] = [0; FLEN_CARD];

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    strcpy_safe(&mut valstring, cs!(c"("));
    ffr2e(value[0], decim, &mut tmpstring, status); /* convert to string */
    if strlen_safe(&valstring) + strlen_safe(&tmpstring) + 2 > FLEN_VALUE - 1 {
        ffpmsg_str("Error converting complex to string (ffpkyc)");
        *status = BAD_F2C;
        return *status;
    }

    strcat_safe(&mut valstring, &tmpstring);
    strcat_safe(&mut valstring, cs!(c", "));
    ffr2e(value[1], decim, &mut tmpstring, status); /* convert to string */

    if strlen_safe(&valstring) + strlen_safe(&tmpstring) + 1 > FLEN_VALUE - 1 {
        ffpmsg_str("Error converting complex to string (ffpkyc)");
        *status = BAD_F2C;
        return *status;
    }

    strcat_safe(&mut valstring, &tmpstring);
    strcat_safe(&mut valstring, cs!(c")"));

    ffmkky_safe(keyname, &valstring, comm, &mut card, status); /* construct the keyword*/
    ffprec_safe(fptr, &card, status); /* write the keyword*/

    *status
}

/*--------------------------------------------------------------------------*/
/// Write (put) the keyword, value and comment into the FITS header.
/// Writes an complex double keyword value. Format = (realvalue, imagvalue)
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffpkym(
    fptr: *mut fitsfile,    /* I - FITS file pointer                   */
    keyname: *const c_char, /* I - name of keyword to write            */
    value: *const [f64; 2], /* I - keyword value (real, imaginary)     */
    decim: c_int,           /* I - number of decimal places to display */
    comm: *const c_char,    /* I - keyword comment                     */
    status: *mut c_int,     /* IO - error status                       */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let value = value.as_ref().expect(NULL_MSG);

        raw_to_slice!(keyname);
        nullable_slice_cstr!(comm);

        ffpkym_safe(fptr, keyname, value, decim, comm, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Write (put) the keyword, value and comment into the FITS header.
/// Writes an complex double keyword value. Format = (realvalue, imagvalue)
pub fn ffpkym_safe(
    fptr: &mut fitsfile,     /* I - FITS file pointer                   */
    keyname: &[c_char],      /* I - name of keyword to write            */
    value: &[f64; 2],        /* I - keyword value (real, imaginary)     */
    decim: c_int,            /* I - number of decimal places to display */
    comm: Option<&[c_char]>, /* I - keyword comment                     */
    status: &mut c_int,      /* IO - error status                       */
) -> c_int {
    let mut valstring: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut tmpstring: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut card: [c_char; FLEN_CARD] = [0; FLEN_CARD];

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    strcpy_safe(&mut valstring, cs!(c"("));
    ffd2e(value[0], decim, &mut tmpstring, status); /* convert to string */
    if strlen_safe(&valstring) + strlen_safe(&tmpstring) + 2 > FLEN_VALUE - 1 {
        ffpmsg_str("Error converting complex to string (ffpkym)");
        *status = BAD_F2C;
        return *status;
    }

    strcat_safe(&mut valstring, &tmpstring);
    strcat_safe(&mut valstring, cs!(c", "));
    ffd2e(value[1], decim, &mut tmpstring, status); /* convert to string */

    if strlen_safe(&valstring) + strlen_safe(&tmpstring) + 1 > FLEN_VALUE - 1 {
        ffpmsg_str("Error converting complex to string (ffpkym)");
        *status = BAD_F2C;
        return *status;
    }

    strcat_safe(&mut valstring, &tmpstring);
    strcat_safe(&mut valstring, cs!(c")"));

    ffmkky_safe(keyname, &valstring, comm, &mut card, status); /* construct the keyword*/
    ffprec_safe(fptr, &card, status); /* write the keyword*/

    *status
}

/*--------------------------------------------------------------------------*/
/// Write (put) the keyword, value and comment into the FITS header.
/// Writes an complex float keyword value. Format = (realvalue, imagvalue)
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffpkfc(
    fptr: *mut fitsfile,    /* I - FITS file pointer                   */
    keyname: *const c_char, /* I - name of keyword to write            */
    value: *const [f32; 2], /* I - keyword value (real, imaginary)     */
    decim: c_int,           /* I - number of decimal places to display */
    comm: *const c_char,    /* I - keyword comment                     */
    status: *mut c_int,     /* IO - error status                       */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let value = value.as_ref().expect(NULL_MSG);

        raw_to_slice!(keyname);
        nullable_slice_cstr!(comm);

        ffpkfc_safe(fptr, keyname, value, decim, comm, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Write (put) the keyword, value and comment into the FITS header.
/// Writes an complex float keyword value. Format = (realvalue, imagvalue)
pub fn ffpkfc_safe(
    fptr: &mut fitsfile,     /* I - FITS file pointer                   */
    keyname: &[c_char],      /* I - name of keyword to write            */
    value: &[f32; 2],        /* I - keyword value (real, imaginary)     */
    decim: c_int,            /* I - number of decimal places to display */
    comm: Option<&[c_char]>, /* I - keyword comment                     */
    status: &mut c_int,      /* IO - error status                       */
) -> c_int {
    let mut valstring: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut tmpstring: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut card: [c_char; FLEN_CARD] = [0; FLEN_CARD];

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    strcpy_safe(&mut valstring, cs!(c"("));
    ffr2f(value[0], decim, &mut tmpstring, status); /* convert to string */
    if strlen_safe(&valstring) + strlen_safe(&tmpstring) + 2 > FLEN_VALUE - 1 {
        ffpmsg_str("Error converting complex to string (ffpkfc)");
        *status = BAD_F2C;
        return *status;
    }

    strcat_safe(&mut valstring, &tmpstring);
    strcat_safe(&mut valstring, cs!(c", "));
    ffr2f(value[1], decim, &mut tmpstring, status); /* convert to string */
    if strlen_safe(&valstring) + strlen_safe(&tmpstring) + 1 > FLEN_VALUE - 1 {
        ffpmsg_str("Error converting complex to string (ffpkfc)");
        *status = BAD_F2C;
        return *status;
    }

    strcat_safe(&mut valstring, &tmpstring);
    strcat_safe(&mut valstring, cs!(c")"));

    ffmkky_safe(keyname, &valstring, comm, &mut card, status); /* construct the keyword*/
    ffprec_safe(fptr, &card, status); /* write the keyword*/

    *status
}

/*--------------------------------------------------------------------------*/
/// Write (put) the keyword, value and comment into the FITS header.
/// Writes an complex double keyword value. Format = (realvalue, imagvalue)
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffpkfm(
    fptr: *mut fitsfile,    /* I - FITS file pointer                   */
    keyname: *const c_char, /* I - name of keyword to write            */
    value: *const [f64; 2], /* I - keyword value (real, imaginary)     */
    decim: c_int,           /* I - number of decimal places to display */
    comm: *const c_char,    /* I - keyword comment                     */
    status: *mut c_int,     /* IO - error status                       */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let value = value.as_ref().expect(NULL_MSG);

        raw_to_slice!(keyname);
        nullable_slice_cstr!(comm);

        ffpkfm_safe(fptr, keyname, value, decim, comm, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Write (put) the keyword, value and comment into the FITS header.
/// Writes an complex double keyword value. Format = (realvalue, imagvalue)
pub fn ffpkfm_safe(
    fptr: &mut fitsfile,     /* I - FITS file pointer                   */
    keyname: &[c_char],      /* I - name of keyword to write            */
    value: &[f64; 2],        /* I - keyword value (real, imaginary)     */
    decim: c_int,            /* I - number of decimal places to display */
    comm: Option<&[c_char]>, /* I - keyword comment                     */
    status: &mut c_int,      /* IO - error status                       */
) -> c_int {
    let mut valstring: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut tmpstring: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut card: [c_char; FLEN_CARD] = [0; FLEN_CARD];

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    strcpy_safe(&mut valstring, cs!(c"("));
    ffd2f(value[0], decim, &mut tmpstring, status); /* convert to string */
    if strlen_safe(&valstring) + strlen_safe(&tmpstring) + 1 > FLEN_VALUE - 1 {
        ffpmsg_str("Error converting complex to string (ffpkfm)");
        *status = BAD_F2C;
        return *status;
    }

    strcat_safe(&mut valstring, &tmpstring);
    strcat_safe(&mut valstring, cs!(c", "));
    ffd2f(value[1], decim, &mut tmpstring, status); /* convert to string */
    if strlen_safe(&valstring) + strlen_safe(&tmpstring) + 1 > FLEN_VALUE - 1 {
        ffpmsg_str("Error converting complex to string (ffpkfm)");
        *status = BAD_F2C;
        return *status;
    }

    strcat_safe(&mut valstring, &tmpstring);
    strcat_safe(&mut valstring, cs!(c")"));

    ffmkky_safe(keyname, &valstring, comm, &mut card, status); /* construct the keyword*/
    ffprec_safe(fptr, &card, status); /* write the keyword*/

    *status
}

/*--------------------------------------------------------------------------*/
/// Write (put) the keyword, value and comment into the FITS header.
///
/// This routine is a modified version of ffpkys which supports the
/// HEASARC long string convention and can write arbitrarily long string
/// keyword values.  The value is continued over multiple keywords that
/// have the name CONTINUE without an equal sign in column 9 of the card.
/// This routine also supports simple string keywords which are less than
/// 75 characters in length.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffpkls(
    fptr: *mut fitsfile,    /* I - FITS file pointer        */
    keyname: *const c_char, /* I - name of keyword to write */
    value: *const c_char,   /* I - keyword value            */
    comm: *const c_char,    /* I - keyword comment          */
    status: *mut c_int,     /* IO - error status            */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        raw_to_slice!(keyname);
        raw_to_slice!(value);
        nullable_slice_cstr!(comm);

        if *status > 0 {
            return *status;
        }

        ffpkls_safe(fptr, keyname, value, comm, status) // call the safe version
    }
}

/*--------------------------------------------------------------------------*/
/// Write (put) the keyword, value and comment into the FITS header.
///
/// This routine is a modified version of ffpkys which supports the
/// HEASARC long string convention and can write arbitrarily long string
/// keyword values.  The value is continued over multiple keywords that
/// have the name CONTINUE without an equal sign in column 9 of the card.
/// This routine also supports simple string keywords which are less than
/// 75 characters in length.
pub fn ffpkls_safe(
    fptr: &mut fitsfile,     /* I - FITS file pointer        */
    keyname: &[c_char],      /* I - name of keyword to write */
    value: &[c_char],        /* I - keyword value            */
    comm: Option<&[c_char]>, /* I - keyword comment          */
    status: &mut c_int,      /* IO - error status            */
) -> c_int {
    if *status > 0 {
        return *status;
    }

    fits_make_longstr_key_util(fptr, keyname, value, comm, -1, status);

    *status
}

/*--------------------------------------------------------------------------*/
/// Write (put) the keyword, value and comment into the FITS header.
///
/// This routine is a modified version of ffpkys which supports the
/// HEASARC long string convention and can write arbitrarily long string
/// keyword values.  The value is continued over multiple keywords that
/// have the name CONTINUE without an equal sign in column 9 of the card.
/// This routine also supports simple string keywords which are less than
/// 75 characters in length.
pub fn fits_make_longstr_key_util(
    fptr: &mut fitsfile,     /* I - FITS file pointer        */
    keyname: &[c_char],      /* I - name of keyword to write */
    value: &[c_char],        /* I - keyword value            */
    comm: Option<&[c_char]>, /* I - keyword comment          */
    mut position: c_int,     /* I - position to insert (-1 for end) */
    status: &mut c_int,      /* IO - error status            */
) -> c_int {
    let mut valstring: [c_char; FLEN_CARD] = [0; FLEN_CARD];
    let mut commstring: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];

    let mut card: [c_char; FLEN_CARD] = [0; FLEN_CARD];
    let mut tmpkeyname: [c_char; FLEN_CARD] = [0; FLEN_CARD]; /* give tmpkeyname same size restriction as in ffmkky */
    let mut tstring: [c_char; FLEN_CARD] = [0; FLEN_CARD];

    let mut nocomment = false;
    let mut tstatus = -1;

    let mut addline = true;

    let mut spaceForComments = 0;
    let mut processingComment = false;
    let mut nblanks = 0;
    let mut allInOne = false;

    let mut vlen = 0;
    let mut nchar = 0;
    let mut nquote = 0;

    /* This setting is arbitrary */
    let fixedSpaceForComments = 50;

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    let mut remainval = strlen_safe(value);
    let mut remaincom = if let Some(comm) = comm {
        strlen_safe(comm)
    } else {
        0
    };

    tmpkeyname[0] = 0;
    let mut ci = 0; //keyname
    while keyname[ci] == bb(b' ') {
        /* skip over leading spaces in name */
        ci += 1;
    }

    let mut cptr = &keyname[ci..];

    strncpy_safe(&mut tmpkeyname, cptr, FLEN_KEYWORD - 1);
    tmpkeyname[FLEN_KEYWORD - 1] = 0;

    let mut namelen = strlen_safe(&tmpkeyname);
    if namelen != 0 {
        /* skip trailing spaces in name */
        ci = namelen - 1; //tmpkeyname
        while tmpkeyname[ci] == bb(b' ') {
            tmpkeyname[ci] = 0;
            ci -= 1;
        }

        cptr = &tmpkeyname[ci..];

        namelen = strlen_safe(&tmpkeyname);
    }

    /* First determine final length of keyword.  ffmkky may prepend
    "HIERARCH " to it, and we need to determine that now using the
    same criteria as ffmkky. */

    let mut finalnamelen = 0;
    let mut maxvalchars = 0;

    if namelen <= 8 && (fftkey_safe(cptr, &mut tstatus) <= 0) {
        /* This a normal 8-character FITS keyword. ffmkky
        will pad it to 8 if necessary, and add "= ". */
        finalnamelen = 10;
        /* 2 additional chars are needed for opening/closing quotes. */
        maxvalchars = (FLEN_CARD - 1) - finalnamelen - 2;
    } else {
        if namelen != 0
            && ((FSTRNCMP(&tmpkeyname, cs!(c"HIERARCH "), 9) == 0)
                || (FSTRNCMP(&tmpkeyname, cs!(c"hierarch "), 9) == 0))
        {
            /* We have an explicitly marked long keyword, so HIERARCH
            will not be prepended.  However it can then have
            " = " or "= ", depending on size of value string.
            For now, assume "= ".

            If we're here, must have 75 > namelen > 9. */
            finalnamelen = namelen + 2;
        } else {
            /* ffmkky is going to prepend "HIERARCH " to the keyword, and " = " or "= ". */
            finalnamelen = namelen + 11;
            if finalnamelen > FLEN_CARD - 1 {
                ffpmsg_str("The following keyword is too long to fit on a card in ffpkls:");
                ffpmsg_slice(keyname);
                *status = BAD_KEYCHAR;
                return *status;
            }
        }
        maxvalchars = (FLEN_CARD - 1) - finalnamelen - 2;
    }

    let mut contin = false;
    let mut next = 0; /* pointer to next character to write */

    while addline {
        if processingComment {
            let comm = comm.unwrap(); // safe to unwrap because we checked it earlier

            if remaincom > (fixedSpaceForComments - 3) {
                strcpy_safe(&mut valstring, cs!(c"'&'"));
                nblanks = (FLEN_CARD - 1) - fixedSpaceForComments - 13;
                valstring[3..(3 + nblanks)].fill(32);
                valstring[nblanks + 3] = 0;
            } else {
                strcpy_safe(&mut valstring, cs!(c"''"));
                nblanks = (FLEN_CARD - 1) - fixedSpaceForComments - 12;
                valstring[2..(2 + nblanks)].fill(32);
                valstring[nblanks + 2] = 0;
            }

            nchar = cmp::min(remaincom, fixedSpaceForComments - 3);
            strncpy_safe(&mut commstring, &comm[next..], nchar);
            commstring[nchar] = 0;
            next += nchar;
            remaincom -= nchar;
        } else {
            vlen = strlen_safe(&value[next..]);
            nquote = 0;
            let mut ichar = 0;
            while ichar < vlen && (ichar + nquote) < maxvalchars {
                if value[next + ichar] == bb(b'\'') {
                    nquote += 1;
                }
                ichar += 1;
            }
            /* Note that (ichar+nquote) can be 1 greater than maxvalchars
            if last processed char is a quote.  Therefore do this check: */
            nchar = cmp::min(ichar, maxvalchars - nquote);

            tstring[0] = 0;
            strncat_safe(&mut tstring, &value[next..], nchar); /* copy string to temp buff */
            /* expand quotes, and put quotes around the string */
            if contin {
                ffs2c_nopad(&tstring, &mut valstring, status);
                vlen = strlen_safe(&valstring);
                spaceForComments = (FLEN_CARD - 1) - (10 + vlen);
            } else {
                ffs2c(&tstring, &mut valstring, status);
                vlen = strlen_safe(&valstring);
                spaceForComments = (FLEN_CARD - 1) - (finalnamelen + vlen);
            }

            /* Check for simplest case where everything fits on first line.*/
            if !contin
                && (remainval == nchar)
                && (finalnamelen + vlen + remaincom + 3 < FLEN_CARD)
                && remaincom < fixedSpaceForComments - 3
            {
                allInOne = true;
            }

            if !allInOne {
                /* There are 2 situations which require overwriting the last char of
                valstring with a continue symbol '&' */
                if spaceForComments == 0 && (remaincom != 0 || (remainval > nchar)) {
                    nchar -= 1; /* outputting one less character now */

                    if valstring[vlen - 2] != bb(b'\'') {
                        valstring[vlen - 2] = bb(b'&'); /*  overwrite last char with &  */
                    } else {
                        /* last char was a pair of single quotes, so over write both */
                        valstring[vlen - 3] = bb(b'&');
                        valstring[vlen - 1] = 0;
                    }
                } else if (spaceForComments != 0 && nchar < remainval)
                    || (remaincom != 0
                        && (spaceForComments < fixedSpaceForComments
                            || spaceForComments - 3 < remaincom
                            || remaincom > fixedSpaceForComments - 3))
                {
                    /* Cases where '&' should be appended to valstring rather than
                    overwritten.  This would mostly be due to the inclusion
                    of a comment string requiring additional lines.  But there's
                    also the obscure case where the last character that can
                    fit happened to be a single quote.  Since this was removed
                    with the earlier 'nchar = minvlaue()' test, the valstring
                    must be continued even though it's one space short of filling
                    this line.  We then append it with a '&'. */

                    valstring[vlen - 1] = bb(b'&');
                    valstring[vlen] = bb(b'\'');
                    valstring[vlen + 1] = 0;
                    vlen += 1;
                }
            }

            if allInOne {
                nocomment = false;
                /* The allInOne test ensures that comm length will
                fit within FLEN_CARD buffer size */
                if let Some(comm) = comm {
                    strcpy_safe(&mut commstring, comm);
                } else {
                    commstring[0] = 0;
                }
                /* Ensure that loop exits after this iteration */
                remainval = 0;
                remaincom = 0;
            } else if remainval > nchar {
                nocomment = true;
                remainval -= nchar;
                next += nchar;
                maxvalchars = (FLEN_CARD - 1) - 12;
            } else {
                /* We've reached the end of val input.  Now switch to writing
                comment (if any).  This block can only be reached once. */

                /* Do not write comments on this line if fewer than
                fixedSpaceForComments are available for the comment string
                and " / ". */
                nocomment = true;
                remainval = 0;
                next = 0;
                processingComment = true;
                if remaincom != 0 && spaceForComments >= fixedSpaceForComments {
                    let comm = comm.unwrap();
                    nocomment = false;
                    nchar = cmp::min(remaincom, fixedSpaceForComments - 3);
                    strncpy_safe(&mut commstring, comm, nchar);
                    commstring[nchar] = 0;
                    next = nchar;
                    remaincom -= nchar;
                }
            }
        } /* end if processing valstring and not comment */

        if contin {
            /* This is a CONTINUEd keyword */

            if nocomment {
                ffmkky_safe(cs!(c"CONTINUE"), &valstring, None, &mut card, status);
            /* make keyword w/o comment */
            } else {
                ffmkky_safe(
                    cs!(c"CONTINUE"),
                    &valstring,
                    Some(&commstring),
                    &mut card,
                    status,
                ); /* make keyword */
            }
            strncpy_safe(&mut card[8..], cs!(c"   "), 2); /* overwrite the '=' */
        } else if nocomment {
            ffmkky_safe(keyname, &valstring, None, &mut card, status); /* make keyword */
        } else {
            ffmkky_safe(keyname, &valstring, Some(&commstring), &mut card, status);
            /* make keyword */
        }

        if position < 0 {
            ffprec_safe(fptr, &card, status); /* write the keyword */
        } else {
            ffirec_safe(fptr, position, &card, status); /* insert the keyword */
            position += 1;
        }

        contin = true;
        nocomment = false;
        addline = remainval > 0 || remaincom > 0;
    }
    *status
}

/*--------------------------------------------------------------------------*/
/// Write the LONGSTRN keyword and a series of related COMMENT keywords
/// which document that this FITS header may contain long string keyword
/// values which are continued over multiple keywords using the HEASARC
/// long string keyword convention.  If the LONGSTRN keyword already exists
/// then this routine simple returns without doing anything.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffplsw(
    fptr: *mut fitsfile, /* I - FITS file pointer  */
    status: *mut c_int,  /* IO - error status       */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        ffplsw_safe(fptr, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Write the LONGSTRN keyword and a series of related COMMENT keywords
/// which document that this FITS header may contain long string keyword
/// values which are continued over multiple keywords using the HEASARC
/// long string keyword convention.  If the LONGSTRN keyword already exists
/// then this routine simple returns without doing anything.
pub fn ffplsw_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer  */
    status: &mut c_int,  /* IO - error status       */
) -> c_int {
    let mut valstring: [c_char; FLEN_CARD] = [0; FLEN_CARD];
    let mut comm: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];
    let mut tstatus = 0;

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    tstatus = 0;
    if ffgkys_safe(
        fptr,
        cs!(c"LONGSTRN"),
        &mut valstring,
        Some(&mut comm),
        &mut tstatus,
    ) == 0
    {
        return *status; /* keyword already exists, so just return */
    }

    ffpkys_safe(
        fptr,
        cs!(c"LONGSTRN"),
        cs!(c"OGIP 1.0"),
        Some(cs!(c"The HEASARC Long String Convention may be used.")),
        status,
    );

    ffpcom_safe(
        fptr,
        cs!(c"  This FITS file may contain long string keyword values that are"),
        status,
    );

    ffpcom_safe(
        fptr,
        cs!(c"  continued over multiple keywords.  The HEASARC convention uses the &"),
        status,
    );

    ffpcom_safe(
        fptr,
        cs!(c"  character at the end of each substring which is then continued"),
        status,
    );

    ffpcom_safe(
        fptr,
        cs!(c"  on the next keyword which has the name CONTINUE."),
        status,
    );

    *status
}

/*--------------------------------------------------------------------------*/
/// Write (put) the keyword, value and comment into the FITS header.
/// Values equal to 0 will result in a False FITS keyword; any other
/// non-zero value will result in a True FITS keyword.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffpkyl(
    fptr: *mut fitsfile,    /* I - FITS file pointer        */
    keyname: *const c_char, /* I - name of keyword to write */
    value: c_int,           /* I - keyword value            */
    comm: *const c_char,    /* I - keyword comment          */
    status: *mut c_int,     /* IO - error status            */
) -> c_int {
    unsafe {
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);

        raw_to_slice!(keyname);
        nullable_slice_cstr!(comm);

        ffpkyl_safe(fptr, keyname, value, comm, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Write (put) the keyword, value and comment into the FITS header.
/// Values equal to 0 will result in a False FITS keyword; any other
/// non-zero value will result in a True FITS keyword.
pub fn ffpkyl_safe(
    fptr: &mut fitsfile,     /* I - FITS file pointer        */
    keyname: &[c_char],      /* I - name of keyword to write */
    value: c_int,            /* I - keyword value            */
    comm: Option<&[c_char]>, /* I - keyword comment          */
    status: &mut c_int,      /* IO - error status            */
) -> c_int {
    let mut valstring: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut card: [c_char; FLEN_CARD] = [0; FLEN_CARD];

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }
    ffl2c(value, &mut valstring, status); /* convert to formatted string */
    ffmkky_safe(keyname, &valstring, comm, &mut card, status); /* construct the keyword*/
    ffprec_safe(fptr, &card, status); /* write the keyword*/

    *status
}

/*--------------------------------------------------------------------------*/
/// Write (put) the keyword, value and comment into the FITS header.
/// Writes an integer keyword value.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffpkyj(
    fptr: *mut fitsfile,    /* I - FITS file pointer        */
    keyname: *const c_char, /* I - name of keyword to write */
    value: LONGLONG,        /* I - keyword value            */
    comm: *const c_char,    /* I - keyword comment          */
    status: *mut c_int,     /* IO - error status            */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        raw_to_slice!(keyname);
        nullable_slice_cstr!(comm);

        ffpkyj_safe(fptr, keyname, value, comm, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Write (put) the keyword, value and comment into the FITS header.
/// Writes an integer keyword value.
pub fn ffpkyj_safe(
    fptr: &mut fitsfile,     /* I - FITS file pointer        */
    keyname: &[c_char],      /* I - name of keyword to write */
    value: LONGLONG,         /* I - keyword value            */
    comm: Option<&[c_char]>, /* I - keyword comment          */
    status: &mut c_int,      /* IO - error status            */
) -> c_int {
    let mut valstring: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut card: [c_char; FLEN_CARD] = [0; FLEN_CARD];

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    ffi2c(value, &mut valstring, status); /* convert to formatted string */
    ffmkky_safe(keyname, &valstring, comm, &mut card, status); /* construct the keyword*/
    ffprec_safe(fptr, &card, status); /* write the keyword*/

    *status
}

/*--------------------------------------------------------------------------*/
/// Write (put) a 'triple' precision keyword where the integer and
/// fractional parts of the value are passed in separate parameters to
/// increase the total amount of numerical precision.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffpkyt(
    fptr: *mut fitsfile,    /* I - FITS file pointer        */
    keyname: *const c_char, /* I - name of keyword to write */
    intval: c_long,         /* I - integer part of value    */
    fraction: f64,          /* I - fractional part of value */
    comm: *const c_char,    /* I - keyword comment          */
    status: *mut c_int,     /* IO - error status            */
) -> c_int {
    unsafe {
        raw_to_slice!(keyname);
        nullable_slice_cstr!(comm);

        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        ffpkyt_safe(fptr, keyname, intval, fraction, comm, status)
    }
}
/*--------------------------------------------------------------------------*/
/// Write (put) a 'triple' precision keyword where the integer and
/// fractional parts of the value are passed in separate parameters to
/// increase the total amount of numerical precision.
pub fn ffpkyt_safe(
    fptr: &mut fitsfile,     /* I - FITS file pointer        */
    keyname: &[c_char],      /* I - name of keyword to write */
    intval: c_long,          /* I - integer part of value    */
    fraction: f64,           /* I - fractional part of value */
    comm: Option<&[c_char]>, /* I - keyword comment          */
    status: &mut c_int,      /* IO - error status            */
) -> c_int {
    let mut valstring: [c_char; FLEN_CARD] = [0; FLEN_CARD];
    let mut card: [c_char; FLEN_CARD] = [0; FLEN_CARD];
    let mut fstring: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    if fraction > 1.0 || fraction < 0.0 {
        ffpmsg_str("fraction must be between 0. and 1. (ffpkyt)");
        *status = BAD_F2C;
        return *status;
    }

    ffi2c(intval as LONGLONG, &mut valstring, status); /* convert integer to string */
    ffd2f(fraction, 16, &mut fstring, status); /* convert to 16 decimal string */

    let cptr = strchr_safe(&fstring, bb(b'.')).unwrap(); /* find the decimal point */

    if strlen_safe(&valstring) + strlen_safe(&fstring[cptr..]) > FLEN_VALUE - 1 {
        ffpmsg_str("converted numerical string too long");
        *status = BAD_F2C;
        return *status;
    }
    strcat_safe(&mut valstring, &fstring[cptr..]); /* append the fraction to the integer */

    ffmkky_safe(keyname, &valstring, comm, &mut card, status); /* construct the keyword*/
    ffprec_safe(fptr, &card, status); /* write the keyword*/

    *status
}

/*-----------------------------------------------------------------*/
/// Write 1 or more COMMENT keywords.  If the comment string is too
/// long to fit on a single keyword (72 chars) then it will automatically
/// be continued on multiple CONTINUE keywords.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffpcom(
    fptr: *mut fitsfile, /* I - FITS file pointer   */
    comm: *const c_char, /* I - comment string      */
    status: *mut c_int,  /* IO - error status       */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);
        raw_to_slice!(comm);

        ffpcom_safe(fptr, comm, status)
    }
}

/*-----------------------------------------------------------------*/
/// Write 1 or more COMMENT keywords.  If the comment string is too
/// long to fit on a single keyword (72 chars) then it will automatically
/// be continued on multiple CONTINUE keywords.
pub fn ffpcom_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer   */
    comm: &[c_char],     /* I - comment string      */
    status: &mut c_int,  /* IO - error status       */
) -> c_int {
    let mut card: [c_char; FLEN_CARD] = [0; FLEN_CARD];

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    let mut len = strlen_safe(comm) as isize;
    let mut ii = 0;

    while len > 0 {
        strcpy_safe(&mut card, cs!(c"COMMENT "));
        strncat_safe(&mut card, &comm[ii..], 72);
        ffprec_safe(fptr, &card, status);
        ii += 72;
        len -= 72;
    }

    *status
}

/*-----------------------------------------------------------------*/
/// Write 1 or more HISTORY keywords.  If the history string is too
/// long to fit on a single keyword (72 chars) then it will automatically
/// be continued on multiple HISTORY keywords.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffphis(
    fptr: *mut fitsfile,    /* I - FITS file pointer  */
    history: *const c_char, /* I - history string     */
    status: *mut c_int,     /* IO - error status      */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);
        raw_to_slice!(history);

        ffphis_safe(fptr, history, status)
    }
}

/*-----------------------------------------------------------------*/
/// Write 1 or more HISTORY keywords.  If the history string is too
/// long to fit on a single keyword (72 chars) then it will automatically
/// be continued on multiple HISTORY keywords.
pub fn ffphis_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer  */
    history: &[c_char],  /* I - history string     */
    status: &mut c_int,  /* IO - error status      */
) -> c_int {
    let mut card: [c_char; FLEN_CARD] = [0; FLEN_CARD];

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    let mut len = strlen_safe(history) as isize;
    let mut ii = 0;

    while len > 0 {
        strcpy_safe(&mut card, cs!(c"HISTORY "));
        strncat_safe(&mut card, &history[ii..], 72);
        ffprec_safe(fptr, &card, status);
        ii += 72;
        len -= 72;
    }

    *status
}

/*-----------------------------------------------------------------*/
/// Write the DATE keyword into the FITS header.  If the keyword already
/// exists then the date will simply be updated in the existing keyword.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffpdat(
    fptr: *mut fitsfile, /* I - FITS file pointer  */
    status: *mut c_int,  /* IO - error status      */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        ffpdat_safe(fptr, status)
    }
}

/*-----------------------------------------------------------------*/
/// Write the DATE keyword into the FITS header.  If the keyword already
/// exists then the date will simply be updated in the existing keyword.
pub fn ffpdat_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer  */
    status: &mut c_int,  /* IO - error status      */
) -> c_int {
    let mut timeref = 0;
    let mut date: [c_char; 20] = [0; 20];
    let mut tmzone: [c_char; 10] = [0; 10];
    let mut card: [c_char; FLEN_CARD] = [0; FLEN_CARD];

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    ffgstm_safe(&mut date, Some(&mut timeref), status);

    if timeref > 0 {
        /* GMT not available on this machine */
        strcpy_safe(&mut tmzone, cs!(c" Local"));
    } else {
        strcpy_safe(&mut tmzone, cs!(c" UT"));
    }

    strcpy_safe(&mut card, cs!(c"DATE    = '"));
    strcat_safe(&mut card, &date);
    strcat_safe(
        &mut card,
        cs!(c"' / file creation date (YYYY-MM-DDThh:mm:ss"),
    );
    strcat_safe(&mut card, &tmzone);
    strcat_safe(&mut card, cs!(c")"));

    ffucrd_safe(fptr, cs!(c"DATE"), &card, status);

    *status
}

/*--------------------------------------------------------------------------*/
/// Write (put) an indexed array of keywords with index numbers between
/// NSTART and (NSTART + NKEY -1) inclusive.  Writes string keywords.
/// The value strings will be truncated at 68 characters, and the HEASARC
/// long string keyword convention is not supported by this routine.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffpkns(
    fptr: *mut fitsfile,         /* I - FITS file pointer                    */
    keyroot: *const c_char,      /* I - root name of keywords to write       */
    nstart: c_int,               /* I - starting index number                */
    nkey: c_int,                 /* I - number of keywords to write          */
    value: *const *const c_char, /* I - array of pointers to keyword values  */
    comm: *const *const c_char,  /* I - array of pointers to keyword comment */
    status: *mut c_int,          /* IO - error status                        */
) -> c_int {
    unsafe {
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);

        raw_to_slice!(keyroot);

        let comm = if comm.is_null() {
            &[]
        } else {
            slice::from_raw_parts(comm, nkey as usize)
        };

        let value = slice::from_raw_parts(value, nkey as usize);

        ffpkns_safe(fptr, keyroot, nstart, nkey, value, comm, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Write (put) an indexed array of keywords with index numbers between
/// NSTART and (NSTART + NKEY -1) inclusive.  Writes string keywords.
/// The value strings will be truncated at 68 characters, and the HEASARC
/// long string keyword convention is not supported by this routine.
pub fn ffpkns_safe(
    fptr: &mut fitsfile,     /* I - FITS file pointer                    */
    keyroot: &[c_char],      /* I - root name of keywords to write       */
    nstart: c_int,           /* I - starting index number                */
    nkey: c_int,             /* I - number of keywords to write          */
    value: &[*const c_char], /* I - array of pointers to keyword values  */
    comm: &[*const c_char],  /* I - array of pointers to keyword comment */
    status: &mut c_int,      /* IO - error status                        */
) -> c_int {
    let mut keyname: [c_char; FLEN_KEYWORD] = [0; FLEN_KEYWORD];
    let mut tcomment: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    /* check if first comment string is to be repeated for all the keywords */
    /* by looking to see if the last non-blank character is a '&' char      */

    let mut repeat = false;

    if !comm.is_empty() {
        let comm_item = unsafe { cast_slice(CStr::from_ptr(comm[0]).to_bytes_with_nul()) };

        let mut len = strlen_safe(comm_item);

        while len > 0 && comm_item[len - 1] == bb(b' ') {
            len -= 1; /* ignore trailing blanks */
        }

        if len > 0 && comm_item[len - 1] == bb(b'&') {
            len = cmp::min(len, FLEN_COMMENT);
            tcomment[0] = 0;
            strncat_safe(&mut tcomment, comm_item, len - 1); /* don't copy the final '&' char */
            repeat = true;
        }
    } else {
        repeat = true;
        tcomment[0] = 0;
    }

    let mut ii: usize = 0;
    let mut jj = nstart;
    while ii < nkey as usize {
        ffkeyn_safe(keyroot, jj, &mut keyname, status);

        let value_item = unsafe { cast_slice(CStr::from_ptr(value[ii]).to_bytes_with_nul()) };

        if repeat {
            ffpkys_safe(fptr, &keyname, value_item, Some(&tcomment), status);
        } else {
            let c: Option<&[c_char]> = match comm[ii].is_null() {
                true => None,
                false => Some(unsafe { cast_slice(CStr::from_ptr(comm[ii]).to_bytes_with_nul()) }),
            };

            ffpkys_safe(fptr, &keyname, value_item, c, status);
        }

        if *status > 0 {
            return *status;
        }

        ii += 1;
        jj += 1;
    }
    *status
}

/*--------------------------------------------------------------------------*/
/// Write (put) an indexed array of keywords with index numbers between
/// NSTART and (NSTART + NKEY -1) inclusive.  Writes logical keywords
/// Values equal to zero will be written as a False FITS keyword value; any
/// other non-zero value will result in a True FITS keyword.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffpknl(
    fptr: *mut fitsfile,        /* I - FITS file pointer                    */
    keyroot: *const c_char,     /* I - root name of keywords to write       */
    nstart: c_int,              /* I - starting index number                */
    nkey: c_int,                /* I - number of keywords to write          */
    value: *const c_int,        /* I - array of keyword values              */
    comm: *const *const c_char, /* I - array of pointers to keyword comment */
    status: *mut c_int,         /* IO - error status                        */
) -> c_int {
    unsafe {
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);

        raw_to_slice!(keyroot);

        let comm = if comm.is_null() {
            &[]
        } else {
            slice::from_raw_parts(comm, nkey as usize)
        };

        let value = slice::from_raw_parts(value, nkey as usize);

        ffpknl_safe(fptr, keyroot, nstart, nkey, value, comm, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Write (put) an indexed array of keywords with index numbers between
/// NSTART and (NSTART + NKEY -1) inclusive.  Writes logical keywords
/// Values equal to zero will be written as a False FITS keyword value; any
/// other non-zero value will result in a True FITS keyword.
pub fn ffpknl_safe(
    fptr: &mut fitsfile,    /* I - FITS file pointer                    */
    keyroot: &[c_char],     /* I - root name of keywords to write       */
    nstart: c_int,          /* I - starting index number                */
    nkey: c_int,            /* I - number of keywords to write          */
    value: &[c_int],        /* I - array of keyword values              */
    comm: &[*const c_char], /* I - array of pointers to keyword comment */
    status: &mut c_int,     /* IO - error status                        */
) -> c_int {
    let mut keyname: [c_char; FLEN_KEYWORD] = [0; FLEN_KEYWORD];
    let mut tcomment: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    /* check if first comment string is to be repeated for all the keywords */
    /* by looking to see if the last non-blank character is a '&' char      */

    let mut repeat = false;

    if !comm.is_empty() {
        let comm_item = unsafe { cast_slice(CStr::from_ptr(comm[0]).to_bytes_with_nul()) };

        let mut len = strlen_safe(comm_item);

        while len > 0 && comm_item[len - 1] == bb(b' ') {
            len -= 1; /* ignore trailing blanks */
        }

        if len > 0 && comm_item[len - 1] == bb(b'&') {
            len = cmp::min(len, FLEN_COMMENT);
            tcomment[0] = 0;
            strncat_safe(&mut tcomment, comm_item, len - 1); /* don't copy the final '&' char */
            repeat = true;
        }
    } else {
        repeat = true;
        tcomment[0] = 0;
    }

    let mut ii: usize = 0;
    let mut jj = nstart;
    while ii < nkey as usize {
        ffkeyn_safe(keyroot, jj, &mut keyname, status);

        if repeat {
            ffpkyl_safe(fptr, &keyname, value[ii], Some(&tcomment), status);
        } else {
            let comm_item = if comm[ii].is_null() {
                None
            } else {
                Some(unsafe { cast_slice(CStr::from_ptr(comm[ii]).to_bytes_with_nul()) })
            };
            ffpkyl_safe(fptr, &keyname, value[ii], comm_item, status);
        }
        if *status > 0 {
            return *status;
        }
        ii += 1;
        jj += 1;
    }

    *status
}

/*--------------------------------------------------------------------------*/
/// Write (put) an indexed array of keywords with index numbers between
/// NSTART and (NSTART + NKEY -1) inclusive.  Write integer keywords
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffpknj(
    fptr: *mut fitsfile,        /* I - FITS file pointer                    */
    keyroot: *const c_char,     /* I - root name of keywords to write       */
    nstart: c_int,              /* I - starting index number                */
    nkey: c_int,                /* I - number of keywords to write          */
    value: *const c_long,       /* I - array of keyword values              */
    comm: *const *const c_char, /* I - array of pointers to keyword comment */
    status: *mut c_int,         /* IO - error status                        */
) -> c_int {
    unsafe {
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);

        raw_to_slice!(keyroot);

        let comm = if comm.is_null() {
            &[]
        } else {
            slice::from_raw_parts(comm, nkey as usize)
        };

        let value = slice::from_raw_parts(value, nkey as usize);

        ffpknj_safe(fptr, keyroot, nstart, nkey, value, comm, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Write (put) an indexed array of keywords with index numbers between
/// NSTART and (NSTART + NKEY -1) inclusive.  Write integer keywords
pub fn ffpknj_safe(
    fptr: &mut fitsfile,    /* I - FITS file pointer                    */
    keyroot: &[c_char],     /* I - root name of keywords to write       */
    nstart: c_int,          /* I - starting index number                */
    nkey: c_int,            /* I - number of keywords to write          */
    value: &[c_long],       /* I - array of keyword values              */
    comm: &[*const c_char], /* I - array of pointers to keyword comment */
    status: &mut c_int,     /* IO - error status                        */
) -> c_int {
    let mut keyname: [c_char; FLEN_KEYWORD] = [0; FLEN_KEYWORD];
    let mut tcomment: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    /* check if first comment string is to be repeated for all the keywords */
    /* by looking to see if the last non-blank character is a '&' char      */

    let mut repeat = false;

    if !comm.is_empty() {
        let comm_item = unsafe { cast_slice(CStr::from_ptr(comm[0]).to_bytes_with_nul()) };

        let mut len = strlen_safe(comm_item);

        while len > 0 && comm_item[len - 1] == bb(b' ') {
            len -= 1; /* ignore trailing blanks */
        }

        if len > 0 && comm_item[len - 1] == bb(b'&') {
            len = cmp::min(len, FLEN_COMMENT);
            tcomment[0] = 0;
            strncat_safe(&mut tcomment, comm_item, len - 1); /* don't copy the final '&' char */
            repeat = true;
        }
    } else {
        repeat = true;
        tcomment[0] = 0;
    }

    let mut ii: usize = 0;
    let mut jj = nstart;
    while ii < nkey as usize {
        ffkeyn_safe(keyroot, jj, &mut keyname, status);

        if repeat {
            ffpkyj_safe(
                fptr,
                &keyname,
                value[ii] as LONGLONG,
                Some(&tcomment),
                status,
            );
        } else {
            let comm_item = if comm[ii].is_null() {
                None
            } else {
                Some(unsafe { cast_slice(CStr::from_ptr(comm[ii]).to_bytes_with_nul()) })
            };
            ffpkyj_safe(fptr, &keyname, value[ii] as LONGLONG, comm_item, status);
        }
        if *status > 0 {
            return *status;
        }
        ii += 1;
        jj += 1;
    }

    *status
}

/*--------------------------------------------------------------------------*/
/// Write (put) an indexed array of keywords with index numbers between
/// NSTART and (NSTART + NKEY -1) inclusive.  Writes fixed float values.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffpknf(
    fptr: *mut fitsfile,        /* I - FITS file pointer                    */
    keyroot: *const c_char,     /* I - root name of keywords to write       */
    nstart: c_int,              /* I - starting index number                */
    nkey: c_int,                /* I - number of keywords to write          */
    value: *const f32,          /* I - array of keyword values              */
    decim: c_int,               /* I - number of decimals to display        */
    comm: *const *const c_char, /* I - array of pointers to keyword comment */
    status: *mut c_int,         /* IO - error status                        */
) -> c_int {
    unsafe {
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);

        raw_to_slice!(keyroot);

        let comm = if comm.is_null() {
            &[]
        } else {
            slice::from_raw_parts(comm, nkey as usize)
        };

        let value = slice::from_raw_parts(value, nkey as usize);

        ffpknf_safe(fptr, keyroot, nstart, nkey, value, decim, comm, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Write (put) an indexed array of keywords with index numbers between
/// NSTART and (NSTART + NKEY -1) inclusive.  Writes fixed float values.
pub fn ffpknf_safe(
    fptr: &mut fitsfile,    /* I - FITS file pointer                    */
    keyroot: &[c_char],     /* I - root name of keywords to write       */
    nstart: c_int,          /* I - starting index number                */
    nkey: c_int,            /* I - number of keywords to write          */
    value: &[f32],          /* I - array of keyword values              */
    decim: c_int,           /* I - number of decimals to display        */
    comm: &[*const c_char], /* I - array of pointers to keyword comment */
    status: &mut c_int,     /* IO - error status                        */
) -> c_int {
    let mut keyname: [c_char; FLEN_KEYWORD] = [0; FLEN_KEYWORD];
    let mut tcomment: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    /* check if first comment string is to be repeated for all the keywords */
    /* by looking to see if the last non-blank character is a '&' char      */

    let mut repeat = false;

    if !comm.is_empty() {
        let comm_item = unsafe { cast_slice(CStr::from_ptr(comm[0]).to_bytes_with_nul()) };

        let mut len = strlen_safe(comm_item);

        while len > 0 && comm_item[len - 1] == bb(b' ') {
            len -= 1; /* ignore trailing blanks */
        }

        if len > 0 && comm_item[len - 1] == bb(b'&') {
            len = cmp::min(len, FLEN_COMMENT);
            tcomment[0] = 0;
            strncat_safe(&mut tcomment, comm_item, len - 1); /* don't copy the final '&' char */
            repeat = true;
        }
    } else {
        repeat = true;
        tcomment[0] = 0;
    }

    let mut ii: usize = 0;
    let mut jj = nstart;
    while ii < nkey as usize {
        ffkeyn_safe(keyroot, jj, &mut keyname, status);

        if repeat {
            ffpkyf_safe(fptr, &keyname, value[ii], decim, Some(&tcomment), status);
        } else {
            let comm_item = if comm[ii].is_null() {
                None
            } else {
                Some(unsafe { cast_slice(CStr::from_ptr(comm[ii]).to_bytes_with_nul()) })
            };
            ffpkyf_safe(fptr, &keyname, value[ii], decim, comm_item, status);
        }
        if *status > 0 {
            return *status;
        }
        ii += 1;
        jj += 1;
    }

    *status
}

/*--------------------------------------------------------------------------*/
/// Write (put) an indexed array of keywords with index numbers between
/// NSTART and (NSTART + NKEY -1) inclusive.  Writes exponential float values.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffpkne(
    fptr: *mut fitsfile,        /* I - FITS file pointer                    */
    keyroot: *const c_char,     /* I - root name of keywords to write       */
    nstart: c_int,              /* I - starting index number                */
    nkey: c_int,                /* I - number of keywords to write          */
    value: *const f32,          /* I - array of keyword values              */
    decim: c_int,               /* I - number of decimals to display        */
    comm: *const *const c_char, /* I - array of pointers to keyword comment */
    status: *mut c_int,         /* IO - error status                        */
) -> c_int {
    unsafe {
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);

        raw_to_slice!(keyroot);

        let comm = if comm.is_null() {
            &[]
        } else {
            slice::from_raw_parts(comm, nkey as usize)
        };

        let value = slice::from_raw_parts(value, nkey as usize);

        ffpkne_safe(fptr, keyroot, nstart, nkey, value, decim, comm, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Write (put) an indexed array of keywords with index numbers between
/// NSTART and (NSTART + NKEY -1) inclusive.  Writes exponential float values.
pub fn ffpkne_safe(
    fptr: &mut fitsfile,    /* I - FITS file pointer                    */
    keyroot: &[c_char],     /* I - root name of keywords to write       */
    nstart: c_int,          /* I - starting index number                */
    nkey: c_int,            /* I - number of keywords to write          */
    value: &[f32],          /* I - array of keyword values              */
    decim: c_int,           /* I - number of decimals to display        */
    comm: &[*const c_char], /* I - array of pointers to keyword comment */
    status: &mut c_int,     /* IO - error status                        */
) -> c_int {
    let mut keyname: [c_char; FLEN_KEYWORD] = [0; FLEN_KEYWORD];
    let mut tcomment: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    /* check if first comment string is to be repeated for all the keywords */
    /* by looking to see if the last non-blank character is a '&' char      */

    let mut repeat = false;

    if !comm.is_empty() {
        let comm_item = unsafe { cast_slice(CStr::from_ptr(comm[0]).to_bytes_with_nul()) };

        let mut len = strlen_safe(comm_item);

        while len > 0 && comm_item[len - 1] == bb(b' ') {
            len -= 1; /* ignore trailing blanks */
        }

        if len > 0 && comm_item[len - 1] == bb(b'&') {
            len = cmp::min(len, FLEN_COMMENT);
            tcomment[0] = 0;
            strncat_safe(&mut tcomment, comm_item, len - 1); /* don't copy the final '&' char */
            repeat = true;
        }
    } else {
        repeat = true;
        tcomment[0] = 0;
    }

    let mut ii: usize = 0;
    let mut jj = nstart;
    while ii < nkey as usize {
        ffkeyn_safe(keyroot, jj, &mut keyname, status);

        if repeat {
            ffpkye_safe(fptr, &keyname, value[ii], decim, Some(&tcomment), status);
        } else {
            let comm_item = if comm[ii].is_null() {
                None
            } else {
                Some(unsafe { cast_slice(CStr::from_ptr(comm[ii]).to_bytes_with_nul()) })
            };
            ffpkye_safe(fptr, &keyname, value[ii], decim, comm_item, status);
        }
        if *status > 0 {
            return *status;
        }
        ii += 1;
        jj += 1;
    }

    *status
}

/*--------------------------------------------------------------------------*/
/// Write (put) an indexed array of keywords with index numbers between
/// NSTART and (NSTART + NKEY -1) inclusive.  Writes fixed double values.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffpkng(
    fptr: *mut fitsfile,        /* I - FITS file pointer                    */
    keyroot: *const c_char,     /* I - root name of keywords to write       */
    nstart: c_int,              /* I - starting index number                */
    nkey: c_int,                /* I - number of keywords to write          */
    value: *const f64,          /* I - array of keyword values              */
    decim: c_int,               /* I - number of decimals to display        */
    comm: *const *const c_char, /* I - array of pointers to keyword comment */
    status: *mut c_int,         /* IO - error status                        */
) -> c_int {
    unsafe {
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);

        raw_to_slice!(keyroot);

        let comm = if comm.is_null() {
            &[]
        } else {
            slice::from_raw_parts(comm, nkey as usize)
        };

        let value = slice::from_raw_parts(value, nkey as usize);

        ffpkng_safe(fptr, keyroot, nstart, nkey, value, decim, comm, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Write (put) an indexed array of keywords with index numbers between
/// NSTART and (NSTART + NKEY -1) inclusive.  Writes fixed double values.
pub fn ffpkng_safe(
    fptr: &mut fitsfile,    /* I - FITS file pointer                    */
    keyroot: &[c_char],     /* I - root name of keywords to write       */
    nstart: c_int,          /* I - starting index number                */
    nkey: c_int,            /* I - number of keywords to write          */
    value: &[f64],          /* I - array of keyword values              */
    decim: c_int,           /* I - number of decimals to display        */
    comm: &[*const c_char], /* I - array of pointers to keyword comment */
    status: &mut c_int,     /* IO - error status                        */
) -> c_int {
    let mut keyname: [c_char; FLEN_KEYWORD] = [0; FLEN_KEYWORD];
    let mut tcomment: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    /* check if first comment string is to be repeated for all the keywords */
    /* by looking to see if the last non-blank character is a '&' char      */

    let mut repeat = false;

    if !comm.is_empty() {
        let comm_item = unsafe { cast_slice(CStr::from_ptr(comm[0]).to_bytes_with_nul()) };

        let mut len = strlen_safe(comm_item);

        while len > 0 && comm_item[len - 1] == bb(b' ') {
            len -= 1; /* ignore trailing blanks */
        }

        if len > 0 && comm_item[len - 1] == bb(b'&') {
            len = cmp::min(len, FLEN_COMMENT);
            tcomment[0] = 0;
            strncat_safe(&mut tcomment, comm_item, len - 1); /* don't copy the final '&' char */
            repeat = true;
        }
    } else {
        repeat = true;
        tcomment[0] = 0;
    }

    let mut ii: usize = 0;
    let mut jj = nstart;
    while ii < nkey as usize {
        ffkeyn_safe(keyroot, jj, &mut keyname, status);

        if repeat {
            ffpkyg_safe(fptr, &keyname, value[ii], decim, Some(&tcomment), status);
        } else {
            let comm_item = if comm[ii].is_null() {
                None
            } else {
                Some(unsafe { cast_slice(CStr::from_ptr(comm[ii]).to_bytes_with_nul()) })
            };
            ffpkyg_safe(fptr, &keyname, value[ii], decim, comm_item, status);
        }
        if *status > 0 {
            return *status;
        }
        ii += 1;
        jj += 1;
    }

    *status
}

/*--------------------------------------------------------------------------*/
/// Write (put) an indexed array of keywords with index numbers between
/// NSTART and (NSTART + NKEY -1) inclusive.  Writes exponential double values.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffpknd(
    fptr: *mut fitsfile,        /* I - FITS file pointer                    */
    keyroot: *const c_char,     /* I - root name of keywords to write       */
    nstart: c_int,              /* I - starting index number                */
    nkey: c_int,                /* I - number of keywords to write          */
    value: *const f64,          /* I - array of keyword values              */
    decim: c_int,               /* I - number of decimals to display        */
    comm: *const *const c_char, /* I - array of pointers to keyword comment */
    status: *mut c_int,         /* IO - error status                        */
) -> c_int {
    unsafe {
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);

        raw_to_slice!(keyroot);

        let comm = if comm.is_null() {
            &[]
        } else {
            slice::from_raw_parts(comm, nkey as usize)
        };

        let value = slice::from_raw_parts(value, nkey as usize);

        ffpknd_safe(fptr, keyroot, nstart, nkey, value, decim, comm, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Write (put) an indexed array of keywords with index numbers between
/// NSTART and (NSTART + NKEY -1) inclusive.  Writes exponential double values.
pub fn ffpknd_safe(
    fptr: &mut fitsfile,    /* I - FITS file pointer                    */
    keyroot: &[c_char],     /* I - root name of keywords to write       */
    nstart: c_int,          /* I - starting index number                */
    nkey: c_int,            /* I - number of keywords to write          */
    value: &[f64],          /* I - array of keyword values              */
    decim: c_int,           /* I - number of decimals to display        */
    comm: &[*const c_char], /* I - array of pointers to keyword comment */
    status: &mut c_int,     /* IO - error status                        */
) -> c_int {
    let mut keyname: [c_char; FLEN_KEYWORD] = [0; FLEN_KEYWORD];
    let mut tcomment: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    /* check if first comment string is to be repeated for all the keywords */
    /* by looking to see if the last non-blank character is a '&' char      */

    let mut repeat = false;

    if !comm.is_empty() {
        let comm_item = unsafe { cast_slice(CStr::from_ptr(comm[0]).to_bytes_with_nul()) };

        let mut len = strlen_safe(comm_item);

        while len > 0 && comm_item[len - 1] == bb(b' ') {
            len -= 1; /* ignore trailing blanks */
        }

        if len > 0 && comm_item[len - 1] == bb(b'&') {
            len = cmp::min(len, FLEN_COMMENT);
            tcomment[0] = 0;
            strncat_safe(&mut tcomment, comm_item, len - 1); /* don't copy the final '&' char */
            repeat = true;
        }
    } else {
        repeat = true;
        tcomment[0] = 0;
    }

    let mut ii: usize = 0;
    let mut jj = nstart;
    while ii < nkey as usize {
        ffkeyn_safe(keyroot, jj, &mut keyname, status);

        if repeat {
            ffpkyd_safe(fptr, &keyname, value[ii], decim, Some(&tcomment), status);
        } else {
            let comm_item = if comm[ii].is_null() {
                None
            } else {
                Some(unsafe { cast_slice(CStr::from_ptr(comm[ii]).to_bytes_with_nul()) })
            };

            ffpkyd_safe(fptr, &keyname, value[ii], decim, comm_item, status);
        }
        if *status > 0 {
            return *status;
        }
        ii += 1;
        jj += 1;
    }

    *status
}

/*--------------------------------------------------------------------------*/
// Write (put) the keyword, value and comment into the FITS header.
// Writes a keyword value with the datatype specified by the 2nd argument.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffpky(
    fptr: *mut fitsfile,    /* I - FITS file pointer        */
    datatype: c_int,        /* I - datatype of the value    */
    keyname: *const c_char, /* I - name of keyword to write */
    value: *const c_void,   /* I - keyword value            */
    comm: *const c_char,    /* I - keyword comment          */
    status: *mut c_int,     /* IO - error status            */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        raw_to_slice!(keyname);
        nullable_slice_cstr!(comm);

        let datatype_with_data = KeywordDatatype::from_datatype(datatype, value.cast_mut());

        ffpky_safe(fptr, datatype_with_data, keyname, comm, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Write (put) the keyword, value and comment into the FITS header.
/// Writes a keyword value with the datatype specified by the 2nd argument.
///
/// Heavily modified to use safe Rust types and idioms.
pub fn ffpky_safe(
    fptr: &mut fitsfile,       /* I - FITS file pointer        */
    datatype: KeywordDatatype, /* I - datatype of the value    */
    keyname: &[c_char],        /* I - name of keyword to write */
    comm: Option<&[c_char]>,   /* I - keyword comment          */
    status: &mut c_int,        /* IO - error status            */
) -> c_int {
    let mut errmsg: [c_char; FLEN_ERRMSG] = [0; FLEN_ERRMSG];

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    match datatype {
        KeywordDatatype::TSTRING(value) => {
            ffpkys_safe(fptr, keyname, value, comm, status);
        }
        KeywordDatatype::TBYTE(value) => {
            ffpkyj_safe(fptr, keyname, LONGLONG::from(*(value)), comm, status);
        }
        KeywordDatatype::TSBYTE(value) => {
            ffpkyj_safe(fptr, keyname, LONGLONG::from(*(value)), comm, status);
        }
        KeywordDatatype::TUSHORT(value) => {
            ffpkyj_safe(fptr, keyname, LONGLONG::from(*(value)), comm, status);
        }
        KeywordDatatype::TSHORT(value) => {
            ffpkyj_safe(fptr, keyname, LONGLONG::from(*(value)), comm, status);
        }
        KeywordDatatype::TUINT(value) => {
            ffpkyg_safe(fptr, keyname, f64::from(*(value)), 0, comm, status);
        }
        KeywordDatatype::TINT(value) => {
            ffpkyj_safe(fptr, keyname, LONGLONG::from(*(value)), comm, status);
        }
        KeywordDatatype::TLOGICAL(value) => {
            ffpkyl_safe(fptr, keyname, *(value), comm, status);
        }
        KeywordDatatype::TULONG(value) => {
            ffpkyuj_safe(fptr, keyname, (*(value)) as ULONGLONG, comm, status);
        }
        KeywordDatatype::TULONGLONG(value) => {
            ffpkyuj_safe(fptr, keyname, (*(value)) as ULONGLONG, comm, status);
        }
        KeywordDatatype::TLONG(value) => {
            ffpkyj_safe(fptr, keyname, (*(value)) as LONGLONG, comm, status);
        }
        KeywordDatatype::TLONGLONG(value) => {
            ffpkyj_safe(fptr, keyname, *(value), comm, status);
        }
        KeywordDatatype::TFLOAT(value) => {
            ffpkye_safe(fptr, keyname, *(value), -7, comm, status);
        }
        KeywordDatatype::TDOUBLE(value) => {
            ffpkyd_safe(fptr, keyname, *(value), -15, comm, status);
        }
        KeywordDatatype::TCOMPLEX(value) => {
            ffpkyc_safe(fptr, keyname, value, -7, comm, status);
        }
        KeywordDatatype::TDBLCOMPLEX(value) => {
            ffpkym_safe(fptr, keyname, value, -15, comm, status);
        }
        _ => {
            int_snprintf!(
                &mut errmsg,
                FLEN_ERRMSG,
                "Bad keyword datatype code: {} (ffpky)",
                datatype.to_datatype_code(),
            );
            ffpmsg_slice(&errmsg);
            *status = BAD_DATATYPE;
        }
    }

    *status
}

/*-------------------------------------------------------------------------*/
/// read keywords from template file and append to the FITS file
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffpktp(
    fptr: *mut fitsfile,     /* I - FITS file pointer       */
    filename: *const c_char, /* I - name of template file   */
    status: *mut c_int,      /* IO - error status           */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);
        raw_to_slice!(filename);

        ffpktp_safe(fptr, filename, status)
    }
}

/*-------------------------------------------------------------------------*/
/// read keywords from template file and append to the FITS file
pub fn ffpktp_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer       */
    filename: &[c_char], /* I - name of template file   */
    status: &mut c_int,  /* IO - error status           */
) -> c_int {
    let mut card: [c_char; FLEN_CARD] = [0; FLEN_CARD];
    let mut template: [c_char; 161] = [0; 161];
    let mut keyname: [c_char; FLEN_KEYWORD] = [0; FLEN_KEYWORD];
    let mut newname: [c_char; FLEN_KEYWORD] = [0; FLEN_KEYWORD];
    let mut keytype: c_int = 0;

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    let filename_str = CStr::from_bytes_until_nul(cast_slice(filename))
        .unwrap()
        .to_str()
        .unwrap();
    let diskfile = File::options().read(true).open(filename_str);

    if diskfile.is_err() {
        ffpmsg_str("ffpktp could not open the following template file:");
        ffpmsg_slice(filename);
        *status = FILE_NOT_OPENED;
        return *status;
    }

    let diskfile = diskfile.unwrap();

    let bufreader = std::io::BufReader::with_capacity(161, diskfile).lines();

    /* get next template line */
    for template_line in bufreader.map_while(Result::ok) {
        let nread = template_line.len();

        if nread == 0 {
            /* end of file */
            break;
        }

        let to_copy = cmp::min(nread, 161);
        template[..to_copy].copy_from_slice(cast_slice(&template_line.as_bytes()[..to_copy])); /* copy the line */

        template[to_copy] = 0; /* make sure string is terminated */
        template[160] = 0; /* make sure string is terminated */

        if ffgthd_safe(&template, &mut card, &mut keytype, status) > 0 {
            /* parse template */
            break;
        }

        strncpy_safe(&mut keyname, &card, 8);
        keyname[8] = 0;

        if keytype == -2 {
            /* rename the card */

            strncpy_safe(&mut newname, &card[40..], 8);
            newname[8] = 0;

            ffmnam_safe(fptr, &keyname, &newname, status);
        } else if keytype == -1 {
            /* delete the card */
            ffdkey_safe(fptr, &keyname, status);
        } else if keytype == 0 {
            /* update the card */
            ffucrd_safe(fptr, &keyname, &card, status);
        } else if keytype == 1 {
            /* append the card */
            ffprec_safe(fptr, &card, status);
        } else {
            /* END card; stop here */
            break;
        }
    }

    *status
}

/*--------------------------------------------------------------------------*/
/// write the TDIMnnn keyword describing the dimensionality of a column
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffptdm(
    fptr: *mut fitsfile,  /* I - FITS file pointer                        */
    colnum: c_int,        /* I - column number                            */
    naxis: c_int,         /* I - number of axes in the data array         */
    naxes: *const c_long, /* I - length of each data axis                 */
    status: *mut c_int,   /* IO - error status                            */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let naxes = slice::from_raw_parts(naxes, naxis as usize);

        ffptdm_safe(fptr, colnum, naxis, naxes, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Write the TDIMnnn keyword describing the dimensionality of a column
pub fn ffptdm_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer                        */
    colnum: c_int,       /* I - column number                            */
    naxis: c_int,        /* I - number of axes in the data array         */
    naxes: &[c_long],    /* I - length of each data axis                 */
    status: &mut c_int,  /* IO - error status                            */
) -> c_int {
    let mut keyname: [c_char; FLEN_KEYWORD] = [0; FLEN_KEYWORD];
    let mut tdimstr: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut comm: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];
    let mut value: [c_char; 80] = [0; 80];
    let mut message: [c_char; FLEN_ERRMSG] = [0; FLEN_ERRMSG];
    let mut totalpix: c_long = 1;
    let mut repeat: c_long = 0;

    if *status > 0 {
        return *status;
    }

    if colnum < 1 || colnum > 999 {
        ffpmsg_str("column number is out of range 1 - 999 (ffptdm)");
        *status = BAD_COL_NUM;
        return *status;
    }

    if naxis < 1 {
        ffpmsg_str("naxis is less than 1 (ffptdm)");
        *status = BAD_DIMEN;
        return *status;
    }

    /* reset position to the correct HDU if necessary */
    if fptr.HDUposition != fptr.Fptr.curhdu {
        ffmahd_safe(fptr, (fptr.HDUposition) + 1, None, status);
    } else if fptr.Fptr.datastart == DATA_UNDEFINED as LONGLONG && ffrdef_safe(fptr, status) > 0 {
        /* rescan header */
        return *status;
    }

    if fptr.Fptr.hdutype != BINARY_TBL {
        ffpmsg_str("Error: The TDIMn keyword is only allowed in BINTABLE extensions (ffptdm)");
        *status = NOT_BTABLE;
        return *status;
    }

    strcpy_safe(&mut tdimstr, cs!(c"(")); /* start constructing the TDIM value */

    for ii in 0..(naxis as usize) {
        if ii > 0 {
            strcat_safe(&mut tdimstr, cs!(c",")); /* append the comma separator */
        }

        if naxes[ii] < 0 {
            ffpmsg_str("one or more TDIM values are less than 0 (ffptdm)");
            *status = BAD_TDIM;
            return *status;
        }

        int_snprintf!(&mut value, 80, "{}", naxes[ii]);
        /* This will either be followed by a ',' or ')'. */
        if strlen_safe(&tdimstr) + strlen_safe(&value) + 1 > FLEN_VALUE - 1 {
            ffpmsg_str("TDIM string too long (ffptdm)");
            *status = BAD_TDIM;
            return *status;
        }
        strcat_safe(&mut tdimstr, &value); /* append the axis size */

        totalpix *= naxes[ii];
    }

    let colptr = fptr.Fptr.tableptr; /* point to first column structure */
    let c = unsafe { slice::from_raw_parts_mut(colptr, fptr.Fptr.tfield as usize) };
    let ci = (colnum - 1) as usize; /* point to the specified column number */

    if (c[ci].trepeat as c_long) != totalpix {
        /* There is an apparent inconsistency between TDIMn and TFORMn. */
        /* The colptr->trepeat value may be out of date, so re-read     */
        /* the TFORMn keyword to be sure.                               */

        ffkeyn_safe(cs!(c"TFORM"), colnum, &mut keyname, status); /* construct TFORMn name  */
        ffgkys_safe(fptr, &keyname, &mut value, None, status); /* read TFORMn keyword    */
        ffbnfm_safe(&value, None, Some(&mut repeat), None, status); /* parse the repeat count */

        if *status > 0 || repeat != totalpix {
            int_snprintf!(
                &mut message,
                FLEN_ERRMSG,
                "column vector length, {}, does not equal TDIMn array size, {}",
                c[ci].trepeat as c_long,
                totalpix,
            );
            ffpmsg_slice(&message);
            *status = BAD_TDIM;
            return *status;
        }
    }

    strcat_safe(&mut tdimstr, cs!(c")")); /* append the closing parenthesis */

    strcpy_safe(&mut comm, cs!(c"size of the multidimensional array"));
    ffkeyn_safe(cs!(c"TDIM"), colnum, &mut keyname, status); /* construct TDIMn name */
    ffpkys_safe(fptr, &keyname, &tdimstr, Some(&comm), status); /* write the keyword */
    *status
}

/*--------------------------------------------------------------------------*/
/// Write the TDIMnnn keyword describing the dimensionality of a column
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffptdmll(
    fptr: *mut fitsfile,    /* I - FITS file pointer                      */
    colnum: c_int,          /* I - column number                            */
    naxis: c_int,           /* I - number of axes in the data array         */
    naxes: *const LONGLONG, /* I - length of each data axis               */
    status: *mut c_int,     /* IO - error status                            */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        let naxes = slice::from_raw_parts(naxes, naxis as usize);

        ffptdmll_safe(fptr, colnum, naxis, naxes, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Write the TDIMnnn keyword describing the dimensionality of a column
pub fn ffptdmll_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer                      */
    colnum: c_int,       /* I - column number                            */
    naxis: c_int,        /* I - number of axes in the data array         */
    naxes: &[LONGLONG],  /* I - length of each data axis                 */
    status: &mut c_int,  /* IO - error status                            */
) -> c_int {
    let mut keyname: [c_char; FLEN_KEYWORD] = [0; FLEN_KEYWORD];
    let mut tdimstr: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut comm: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];
    let mut value: [c_char; 80] = [0; 80];
    let mut message: [c_char; FLEN_ERRMSG] = [0; FLEN_ERRMSG];
    let mut totalpix: LONGLONG = 1;
    let mut repeat: LONGLONG = 0;

    if *status > 0 {
        return *status;
    }

    if colnum < 1 || colnum > 999 {
        ffpmsg_str("column number is out of range 1 - 999 (ffptdm)");
        *status = BAD_COL_NUM;
        return *status;
    }

    if naxis < 1 {
        ffpmsg_str("naxis is less than 1 (ffptdm)");
        *status = BAD_DIMEN;
        return *status;
    }

    /* reset position to the correct HDU if necessary */
    if fptr.HDUposition != fptr.Fptr.curhdu {
        ffmahd_safe(fptr, (fptr.HDUposition) + 1, None, status);
    } else if fptr.Fptr.datastart == DATA_UNDEFINED as LONGLONG && ffrdef_safe(fptr, status) > 0 {
        /* rescan header */
        return *status;
    }

    if fptr.Fptr.hdutype != BINARY_TBL {
        ffpmsg_str("Error: The TDIMn keyword is only allowed in BINTABLE extensions (ffptdm)");
        *status = NOT_BTABLE;
        return *status;
    }

    strcpy_safe(&mut tdimstr, cs!(c"(")); /* start constructing the TDIM value */

    for ii in 0..(naxis as usize) {
        if ii > 0 {
            strcat_safe(&mut tdimstr, cs!(c",")); /* append the comma separator */
        }

        if naxes[ii] < 0 {
            ffpmsg_str("one or more TDIM values are less than 0 (ffptdm)");
            *status = BAD_TDIM;
            return *status;
        }

        /* cast to double because the 64-bit int conversion character in */
        /* sprintf is platform dependent ( %lld, %ld, %I64d )            */

        int_snprintf!(&mut value, 80, "{}", naxes[ii]);

        if strlen_safe(&tdimstr) + strlen_safe(&value) + 1 > FLEN_VALUE - 1 {
            ffpmsg_str("TDIM string too long (ffptdmll)");
            *status = BAD_TDIM;
            return *status;
        }
        strcat_safe(&mut tdimstr, &value); /* append the axis size */

        totalpix *= naxes[ii];
    }

    let colptr = fptr.Fptr.tableptr; /* point to first column structure */
    let c = unsafe { slice::from_raw_parts_mut(colptr, fptr.Fptr.tfield as usize) };
    let ci = (colnum - 1) as usize; /* point to the specified column number */

    if c[ci].trepeat != totalpix {
        /* There is an apparent inconsistency between TDIMn and TFORMn. */
        /* The colptr->trepeat value may be out of date, so re-read     */
        /* the TFORMn keyword to be sure.                               */

        ffkeyn_safe(cs!(c"TFORM"), colnum, &mut keyname, status); /* construct TFORMn name  */
        ffgkys_safe(fptr, &keyname, &mut value, None, status); /* read TFORMn keyword    */
        ffbnfmll_safe(&value, None, Some(&mut repeat), None, status); /* parse the repeat count */

        if *status > 0 || repeat != totalpix {
            int_snprintf!(
                &mut message,
                FLEN_ERRMSG,
                "column vector length, {}, does not equal TDIMn array size, {}",
                c[ci].trepeat,
                totalpix,
            );
            ffpmsg_slice(&message);
            *status = BAD_TDIM;
            return *status;
        }
    }

    strcat_safe(&mut tdimstr, cs!(c")")); /* append the closing parenthesis */

    strcpy_safe(&mut comm, cs!(c"size of the multidimensional array"));
    ffkeyn_safe(cs!(c"TDIM"), colnum, &mut keyname, status); /* construct TDIMn name */
    ffpkys_safe(fptr, &keyname, &tdimstr, Some(&comm), status); /* write the keyword */
    *status
}

/*-----------------------------------------------------------------*/
/// Returns the current date and time in format 'yyyy-mm-ddThh:mm:ss'.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgstm(
    timestr: *mut c_char, /* O  - returned system date and time string  */
    timeref: *mut c_int,  /* O - GMT = 0, Local time = 1  */
    status: *mut c_int,   /* IO - error status      */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let timestr = timestr.as_mut().expect(NULL_MSG);
        let timeref = timeref.as_mut();

        let timestr: &mut [c_char; 20] = slice::from_raw_parts_mut(timestr, 20).try_into().unwrap();

        ffgstm_safe(timestr, timeref, status)
    }
}

/*-----------------------------------------------------------------*/
/// Returns the current date and time in format 'yyyy-mm-ddThh:mm:ss'.
pub fn ffgstm_safe(
    timestr: &mut [c_char; 20], /* O  - returned system date and time string  */
    timeref: Option<&mut c_int>, /* O - GMT = 0, Local time = 1  */
    status: &mut c_int,         /* IO - error status      */
) -> c_int {
    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    if let Some(timeref) = timeref {
        *timeref = 0; /* returning GMT */
        //*timeref = 1; /* returning local time */
    }

    let utc: DateTime<Utc> = Utc::now();

    let dt_str = CString::new(utc.format("%Y-%m-%dT%H:%M:%S").to_string()).unwrap();

    strncpy_safe(timestr, cast_slice(dt_str.as_bytes_with_nul()), 20);

    *status
}

/*-----------------------------------------------------------------*/
/// Construct a date character string
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffdt2s(
    year: c_int,          /* I - year (0 - 9999)           */
    month: c_int,         /* I - month (1 - 12)            */
    day: c_int,           /* I - day (1 - 31)              */
    datestr: *mut c_char, /* O - date string: "YYYY-MM-DD" */
    status: *mut c_int,   /* IO - error status             */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);

        let datestr: &mut [c_char; 11] = slice::from_raw_parts_mut(datestr, 11).try_into().unwrap();

        ffdt2s_safe(year, month, day, datestr, status)
    }
}

/*-----------------------------------------------------------------*/
/// Construct a date character string
pub fn ffdt2s_safe(
    year: c_int,                /* I - year (0 - 9999)           */
    month: c_int,               /* I - month (1 - 12)            */
    day: c_int,                 /* I - day (1 - 31)              */
    datestr: &mut [c_char; 11], /* O - date string: "YYYY-MM-DD" */
    status: &mut c_int,         /* IO - error status             */
) -> c_int {
    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    datestr[0] = 0;

    if ffverifydate_safe(year, month, day, status) > 0 {
        ffpmsg_str("invalid date (ffdt2s)");
        return *status;
    }

    if year >= 1900 && year <= 1998 {
        /* use old 'dd/mm/yy' format */
        int_snprintf!(datestr, 9, "{:02}/{:02}/{:02}", day, month, year - 1900);
    } else {
        /* use the new 'YYYY-MM-DD' format */
        int_snprintf!(datestr, 11, "{:04}-{:02}-{:02}", year, month, day);
    }

    *status
}

/*-----------------------------------------------------------------*/
/// Parse a date character string into year, month, and day values
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffs2dt(
    datestr: *const c_char, /* I - date string: "YYYY-MM-DD" or "dd/mm/yy" */
    year: *mut c_int,       /* O - year (0 - 9999)                         */
    month: *mut c_int,      /* O - month (1 - 12)                          */
    day: *mut c_int,        /* O - day (1 - 31)                            */
    status: *mut c_int,     /* IO - error status                           */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let year = year.as_mut();
        let month = month.as_mut();
        let day = day.as_mut();
        nullable_slice_cstr!(datestr);

        ffs2dt_safe(datestr.map(cast_slice), year, month, day, status)
    }
}

/*-----------------------------------------------------------------*/
/// Parse a date character string into year, month, and day values
pub fn ffs2dt_safe(
    datestr: Option<&[c_char]>, /* I - date string: "YYYY-MM-DD" or "dd/mm/yy" */
    mut year: Option<&mut c_int>, /* O - year (0 - 9999)                         */
    mut month: Option<&mut c_int>, /* O - month (1 - 12)                          */
    mut day: Option<&mut c_int>, /* O - day (1 - 31)                            */
    status: &mut c_int,         /* IO - error status                           */
) -> c_int {
    let lyear: c_int;
    let lmonth: c_int;
    let lday: c_int;

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    if let Some(y) = year.as_deref_mut() {
        *y = 0;
    }
    if let Some(m) = month.as_deref_mut() {
        *m = 0;
    }
    if let Some(d) = day.as_deref_mut() {
        *d = 0;
    }

    let datestr = match datestr {
        Some(d) if !d.is_empty() => d,
        _ => {
            ffpmsg_str("error: null input date string (ffs2dt)");
            *status = BAD_DATE;
            return *status;
        }
    };

    let slen: usize = strlen_safe(datestr);

    if slen == 8 && datestr[2] == bb(b'/') && datestr[5] == bb(b'/') {
        if isdigit_safe(datestr[0])
            && isdigit_safe(datestr[1])
            && isdigit_safe(datestr[3])
            && isdigit_safe(datestr[4])
            && isdigit_safe(datestr[6])
            && isdigit_safe(datestr[7])
        {
            /* this is an old format string: "dd/mm/yy" */
            lyear = parse_c_int(&datestr[6..]) + 1900;
            lmonth = parse_c_int(&datestr[3..]);
            lday = parse_c_int(datestr);

            if let Some(y) = year.as_deref_mut() {
                *y = lyear;
            }
            if let Some(m) = month.as_deref_mut() {
                *m = lmonth;
            }
            if let Some(d) = day.as_deref_mut() {
                *d = lday;
            }
        } else {
            ffpmsg_str("input date string has illegal format (ffs2dt):");
            ffpmsg_slice(datestr);
            *status = BAD_DATE;
            return *status;
        }
    } else if slen >= 10 && datestr[4] == bb(b'-') && datestr[7] == bb(b'-') {
        if isdigit_safe(datestr[0])
            && isdigit_safe(datestr[1])
            && isdigit_safe(datestr[2])
            && isdigit_safe(datestr[3])
            && isdigit_safe(datestr[5])
            && isdigit_safe(datestr[6])
            && isdigit_safe(datestr[8])
            && isdigit_safe(datestr[9])
        {
            if slen > 10 && datestr[10] != bb(b'T') {
                ffpmsg_str("input date string has illegal format (ffs2dt):");
                ffpmsg_slice(datestr);
                *status = BAD_DATE;
                return *status;
            }

            /* this is a new format string: "yyyy-mm-dd" */
            lyear = parse_c_int(datestr);
            lmonth = parse_c_int(&datestr[5..]);
            lday = parse_c_int(&datestr[8..]);

            if let Some(y) = year {
                *y = lyear;
            }
            if let Some(m) = month {
                *m = lmonth;
            }
            if let Some(d) = day {
                *d = lday;
            }
        } else {
            ffpmsg_str("input date string has illegal format (ffs2dt):");
            ffpmsg_slice(datestr);
            *status = BAD_DATE;
            return *status;
        }
    } else {
        ffpmsg_str("input date string has illegal format (ffs2dt):");
        ffpmsg_slice(datestr);
        *status = BAD_DATE;
        return *status;
    }

    if ffverifydate_safe(lyear, lmonth, lday, status) > 0 {
        ffpmsg_str("invalid date (ffs2dt)");
    }

    *status
}

/*-----------------------------------------------------------------*/
/// Construct a date and time character string
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn fftm2s(
    year: c_int,          /* I - year (0 - 9999)           */
    month: c_int,         /* I - month (1 - 12)            */
    day: c_int,           /* I - day (1 - 31)              */
    hour: c_int,          /* I - hour (0 - 23)             */
    minute: c_int,        /* I - minute (0 - 59)           */
    second: c_double,     /* I - second (0. - 60.9999999)  */
    decimals: c_int,      /* I - number of decimal points to write      */
    datestr: *mut c_char, /* O - date string: "YYYY-MM-DDThh:mm:ss.ddd" or "hh:mm:ss.ddd" if year, month day = 0 */
    status: *mut c_int,   /* IO - error status             */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let datestr = datestr.as_mut().expect(NULL_MSG);

        let datestr: &mut [c_char; 24] = slice::from_raw_parts_mut(datestr, 24).try_into().unwrap();

        fftm2s_safe(
            year, month, day, hour, minute, second, decimals, datestr, status,
        )
    }
}

/*-----------------------------------------------------------------*/
/// Construct a date and time character string
pub fn fftm2s_safe(
    year: c_int,            /* I - year (0 - 9999)           */
    month: c_int,           /* I - month (1 - 12)            */
    day: c_int,             /* I - day (1 - 31)              */
    hour: c_int,            /* I - hour (0 - 23)             */
    minute: c_int,          /* I - minute (0 - 59)           */
    second: c_double,       /* I - second (0. - 60.9999999)  */
    decimals: c_int,        /* I - number of decimal points to write      */
    datestr: &mut [c_char], /* O - date string: "YYYY-MM-DDThh:mm:ss.ddd" or "hh:mm:ss.ddd" if year, month day = 0 */
    status: &mut c_int,     /* IO - error status             */
) -> c_int {
    let width: c_int;
    let mut errmsg: [c_char; FLEN_ERRMSG] = [0; FLEN_ERRMSG];

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    datestr[0] = 0;

    if (year != 0 || month != 0 || day != 0) && ffverifydate_safe(year, month, day, status) > 0 {
        ffpmsg_str("invalid date (fftm2s)");
        return *status;
    }

    if hour < 0 || hour > 23 {
        int_snprintf!(
            &mut errmsg,
            FLEN_ERRMSG,
            "input hour value is out of range 0 - 23: {} (fftm2s)",
            hour
        );
        ffpmsg_slice(&errmsg);
        *status = BAD_DATE;
        return *status;
    } else if minute < 0 || minute > 59 {
        int_snprintf!(
            &mut errmsg,
            FLEN_ERRMSG,
            "input minute value is out of range 0 - 59: {} (fftm2s)",
            minute
        );
        ffpmsg_slice(&errmsg);
        *status = BAD_DATE;
        return *status;
    } else if second < 0. || second >= 61. {
        int_snprintf!(
            &mut errmsg,
            FLEN_ERRMSG,
            "input second value is out of range 0 - 60.999: {} (fftm2s)",
            second
        );
        ffpmsg_slice(&errmsg);
        *status = BAD_DATE;
        return *status;
    } else if decimals > 25 {
        int_snprintf!(
            &mut errmsg,
            FLEN_ERRMSG,
            "input decimals value is out of range 0 - 25: {} (fftm2s)",
            decimals
        );
        ffpmsg_slice(&errmsg);
        *status = BAD_DATE;
        return *status;
    }

    if decimals == 0 {
        width = 2;
    } else {
        width = decimals + 3;
    }

    if decimals < 0 {
        /* a negative decimals value means return only the date, not time */
        int_snprintf!(
            datestr,
            datestr.len(),
            "{:04}-{:02}-{:02}",
            year,
            month,
            day
        );
    } else if year == 0 && month == 0 && day == 0 {
        /* return only the time, not the date */
        int_snprintf!(
            datestr,
            datestr.len(),
            "{:02}:{:02}:{:0w$.p$}",
            hour,
            minute,
            second,
            w = width as usize,
            p = decimals as usize
        );
    } else {
        /* return both the time and date */
        int_snprintf!(
            datestr,
            datestr.len(),
            "{:04}-{:02}-{:02}T{:02}:{:02}:{:0w$.p$}",
            year,
            month,
            day,
            hour,
            minute,
            second,
            w = width as usize,
            p = decimals as usize
        );
    }

    *status
}

/*-----------------------------------------------------------------*/
/// Parse a date character string into date and time values
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffs2tm(
    datestr: *const c_char, /* I - date string: "YYYY-MM-DD"    */
    /*     or "YYYY-MM-DDThh:mm:ss.ddd" */
    /*     or "dd/mm/yy"                */
    year: *mut c_int,      /* O - year (0 - 9999)              */
    month: *mut c_int,     /* O - month (1 - 12)               */
    day: *mut c_int,       /* O - day (1 - 31)                 */
    hour: *mut c_int,      /* O - hour (0 - 23)                */
    minute: *mut c_int,    /* O - minute (0 - 59)              */
    second: *mut c_double, /* O - second (0. - 60.9999999)     */
    status: *mut c_int,    /* IO - error status                */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let year = year.as_mut();
        let month = month.as_mut();
        let day = day.as_mut();
        let hour = hour.as_mut();
        let minute = minute.as_mut();
        let second = second.as_mut();

        nullable_slice_cstr!(datestr);

        ffs2tm_safe(
            datestr.map(cast_slice),
            year,
            month,
            day,
            hour,
            minute,
            second,
            status,
        )
    }
}

/*-----------------------------------------------------------------*/
/// Parse a date character string into date and time values
pub fn ffs2tm_safe(
    datestr: Option<&[c_char]>, /* I - date string: "YYYY-MM-DD"    */
    /*     or "YYYY-MM-DDThh:mm:ss.ddd" */
    /*     or "dd/mm/yy"                */
    mut year: Option<&mut c_int>,  /* O - year (0 - 9999)              */
    mut month: Option<&mut c_int>, /* O - month (1 - 12)               */
    mut day: Option<&mut c_int>,   /* O - day (1 - 31)                 */
    mut hour: Option<&mut c_int>,  /* O - hour (0 - 23)                */
    mut minute: Option<&mut c_int>, /* O - minute (0 - 59)              */
    mut second: Option<&mut c_double>, /* O - second (0. - 60.9999999)     */
    status: &mut c_int,            /* IO - error status                */
) -> c_int {
    let slen: usize;
    let mut errmsg: [c_char; FLEN_ERRMSG] = [0; FLEN_ERRMSG];

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    if let Some(y) = year.as_deref_mut() {
        *y = 0;
    }
    if let Some(m) = month.as_deref_mut() {
        *m = 0;
    }
    if let Some(d) = day.as_deref_mut() {
        *d = 0;
    }
    if let Some(h) = hour.as_deref_mut() {
        *h = 0;
    }
    if let Some(min) = minute.as_deref_mut() {
        *min = 0;
    }
    if let Some(s) = second.as_deref_mut() {
        *s = 0.;
    }

    let datestr = match datestr {
        Some(d) if !d.is_empty() => d,
        _ => {
            ffpmsg_str("error: null input date string (ffs2tm)");
            *status = BAD_DATE;
            return *status;
        }
    };

    if datestr.len() > 2
        && (datestr[2] == bb(b'/') || (datestr.len() > 4 && datestr[4] == bb(b'-')))
    {
        /*  Parse the year, month, and date */
        if ffs2dt_safe(Some(datestr), year, month, day, status) > 0 {
            return *status;
        }

        slen = strlen_safe(datestr);
        if slen == 8 || slen == 10 {
            return *status; /* OK, no time fields */
        } else if slen < 19 {
            ffpmsg_str("input date string has illegal format:");
            ffpmsg_slice(datestr);
            *status = BAD_DATE;
            return *status;
        } else if datestr[10] == bb(b'T') {
            if datestr[13] == bb(b':') && datestr[16] == bb(b':') {
                if isdigit_safe(datestr[11])
                    && isdigit_safe(datestr[12])
                    && isdigit_safe(datestr[14])
                    && isdigit_safe(datestr[15])
                    && isdigit_safe(datestr[17])
                    && isdigit_safe(datestr[18])
                {
                    if slen > 19 && datestr[19] != bb(b'.') {
                        ffpmsg_str("input date string has illegal format:");
                        ffpmsg_slice(datestr);
                        *status = BAD_DATE;
                        return *status;
                    }

                    /* this is a new format string: "yyyy-mm-ddThh:mm:ss.dddd" */
                    if let Some(h) = hour.as_deref_mut() {
                        *h = parse_c_int(&datestr[11..]);
                    }
                    if let Some(min) = minute.as_deref_mut() {
                        *min = parse_c_int(&datestr[14..]);
                    }
                    if let Some(s) = second.as_deref_mut() {
                        *s = atof_safe(&datestr[17..]);
                    }
                } else {
                    ffpmsg_str("input date string has illegal format:");
                    ffpmsg_slice(datestr);
                    *status = BAD_DATE;
                    return *status;
                }
            } else {
                ffpmsg_str("input date string has illegal format:");
                ffpmsg_slice(datestr);
                *status = BAD_DATE;
                return *status;
            }
        }
    } else
    /* no date fields */
    if datestr.len() > 5 && datestr[2] == bb(b':') && datestr[5] == bb(b':')
    /* time string */
    {
        if datestr.len() > 7
            && isdigit_safe(datestr[0])
            && isdigit_safe(datestr[1])
            && isdigit_safe(datestr[3])
            && isdigit_safe(datestr[4])
            && isdigit_safe(datestr[6])
            && isdigit_safe(datestr[7])
        {
            /* this is a time string: "hh:mm:ss.dddd" */
            if let Some(h) = hour.as_deref_mut() {
                *h = parse_c_int(datestr);
            }
            if let Some(min) = minute.as_deref_mut() {
                *min = parse_c_int(&datestr[3..]);
            }
            if let Some(s) = second.as_deref_mut() {
                *s = atof_safe(&datestr[6..]);
            }
        } else {
            ffpmsg_str("input date string has illegal format:");
            ffpmsg_slice(datestr);
            *status = BAD_DATE;
            return *status;
        }
    } else {
        ffpmsg_str("input date string has illegal format:");
        ffpmsg_slice(datestr);
        *status = BAD_DATE;
        return *status;
    }

    if let Some(h) = hour.as_deref()
        && (*h < 0 || *h > 23)
    {
        int_snprintf!(
            errmsg,
            FLEN_ERRMSG,
            "hour value is out of range 0 - 23: {} (ffs2tm)",
            *h
        );
        ffpmsg_slice(&errmsg);
        *status = BAD_DATE;
        return *status;
    }

    if let Some(min) = minute.as_deref()
        && (*min < 0 || *min > 59)
    {
        int_snprintf!(
            errmsg,
            FLEN_ERRMSG,
            "minute value is out of range 0 - 59: {} (ffs2tm)",
            *min
        );
        ffpmsg_slice(&errmsg);
        *status = BAD_DATE;
        return *status;
    }

    if let Some(s) = second.as_deref()
        && (*s < 0. || *s >= 61.)
    {
        int_snprintf!(
            errmsg,
            FLEN_ERRMSG,
            "second value is out of range 0 - 60.9999: {} (ffs2tm)",
            *s
        );
        ffpmsg_slice(&errmsg);
        *status = BAD_DATE;
        return *status;
    }

    *status
}

/*--------------------------------------------------------------------------*/
/// his routine is included for backward compatibility with the Fortran FITSIO library.
/// Get current System DaTe (GMT if available)
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgsdt(
    day: *mut c_int,
    month: *mut c_int,
    year: *mut c_int,
    status: *mut c_int,
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let day = day.as_mut().expect(NULL_MSG);
        let month = month.as_mut().expect(NULL_MSG);
        let year = year.as_mut().expect(NULL_MSG);

        ffgsdt_safe(day, month, year, status)
    }
}

/*--------------------------------------------------------------------------*/
/// This routine is included for backward compatibility with the Fortran FITSIO library.
/// Get current System DaTe (GMT if available)
pub fn ffgsdt_safe(
    day: &mut c_int,
    month: &mut c_int,
    year: &mut c_int,
    status: &mut c_int,
) -> c_int {
    let utc: DateTime<Utc> = Utc::now();
    let utc = utc.date_naive();

    *day = utc.day() as c_int;
    *month = utc.month() as c_int;
    *year = utc.year() as c_int;

    *status
}

/*--------------------------------------------------------------------------*/
/// convert value to a null-terminated formatted string.
pub(crate) fn ffi2c(
    ival: LONGLONG,      /* I - value to be converted to a string */
    cval: &mut [c_char], /* O - character string representation of the value */
    status: &mut c_int,  /* IO - error status */
) -> c_int {
    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }
    cval[0] = 0;

    if int_snprintf!(cval, cval.len() - 1, "{}", ival) < 0 {
        ffpmsg_str("Error in ffi2c converting integer to string");
        *status = BAD_I2C;
    }
    *status
}

/*--------------------------------------------------------------------------*/
/// convert  value to a null-terminated formatted string.
pub(crate) fn ffu2c(
    ival: ULONGLONG,     /* I - value to be converted to a string */
    cval: &mut [c_char], /* O - character string representation of the value */
    status: &mut c_int,  /* IO - error status */
) -> c_int {
    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    cval[0] = 0;

    if int_snprintf!(cval, cval.len(), "{}", ival) < 0 {
        ffpmsg_str("Error in ffu2c converting integer to string");
        *status = BAD_I2C;
    }

    *status
}

/*--------------------------------------------------------------------------*/
/// Convert logical value to a null-terminated formatted string.  If the
/// input value == 0, then the output character is the letter F, else
/// the output character is the letter T.  The output string is null terminated.
pub(crate) fn ffl2c(
    lval: c_int,         /* I - value to be converted to a string */
    cval: &mut [c_char], /* O - character string representation of the value */
    status: &mut c_int,  /* IO - error status ) */
) -> c_int {
    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    if lval > 0 {
        strcpy_safe(cval, cs!(c"T"));
    } else {
        strcpy_safe(cval, cs!(c"F"));
    }
    *status
}

/*--------------------------------------------------------------------------*/
/// convert an input string to a quoted string. Leading spaces
/// are significant.  FITS string keyword values must be at least
/// 8 chars long so pad out string with spaces if necessary.
/// (*** This 8 char requirement is now obsolete.  See ffs2c_nopad
/// for an alternative ***)
/// Example:   km/s ==> 'km/s    '
/// Single quote characters in the input string will be replace by
/// two single quote characters. e.g., o'brian ==> 'o''brian'
pub(crate) fn ffs2c(
    instr: &[c_char],      /* I - null terminated input string  */
    outstr: &mut [c_char], /* O - null terminated quoted output string */
    status: &mut c_int,    /* IO - error status */
) -> c_int {
    let mut len = 0;

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    if instr.is_empty() {
        /* a null input pointer?? */

        strcpy_safe(outstr, cs!(c"''")); /* a null FITS string */
        return *status;
    }

    outstr[0] = bb(b'\''); /* start output string with a quote */

    len = strlen_safe(instr);
    if len > 68 {
        len = 68; /* limit input string to 68 chars */
    }
    let mut ii = 0;
    let mut jj = 1;
    while ii < len && jj < 69 {
        outstr[jj] = instr[ii]; /* copy each char from input to output */
        if instr[ii] == bb(b'\'') {
            jj += 1;
            outstr[jj] = bb(b'\''); /* duplicate any apostrophies in the input */
        }
        ii += 1;
        jj += 1;
    }
    while jj < 9 {
        /* pad string so it is at least 8 chars long */
        outstr[jj] = bb(b' ');
        jj += 1;
    }
    if jj == 70 {
        /* only occurs if the last char of string was a quote */
        outstr[69] = 0;
    } else {
        outstr[jj] = bb(b'\''); /* append closing quote character */
        outstr[jj + 1] = 0; /* terminate the string */
    }

    *status
}

/*--------------------------------------------------------------------------*/
/// This performs identically to ffs2c except that it won't pad output
/// strings to make them a minimum of 8 chars long.  The requirement
/// that FITS keyword string values be 8 characters is now obsolete
/// (except for "XTENSION" keyword), but for backwards compatibility we'll
/// keep ffs2c the way it is.  A better solution would be to add another
/// argument to ffs2c for 'pad' or 'nopad', but it is called from many other
/// places in Heasoft outside of CFITSIO.  
pub(crate) fn ffs2c_nopad(
    instr: &[c_char],      /* I - null terminated input string  */
    outstr: &mut [c_char], /* O - null terminated quoted output string */
    status: &mut c_int,    /* IO - error status */
) -> c_int {
    //size_t len, ii, jj;

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    if strlen_safe(instr) == 0 {
        /* a null input pointer?? */
        strcpy_safe(outstr, cs!(c"''")); /* a null FITS string */
        return *status;
    }

    outstr[0] = bb(b'\''); /* start output string with a quote */

    let mut len = strlen_safe(instr);
    if len > 68 {
        len = 68; /* limit input string to 68 chars */
    }

    let mut jj = 1;
    let mut ii = 0;
    while ii < len && jj < 69 {
        outstr[jj] = instr[ii]; /* copy each char from input to output */
        if instr[ii] == bb(b'\'') {
            jj += 1;
            outstr[jj] = bb(b'\''); /* duplicate any apostrophies in the input */
        }
        ii += 1;
        jj += 1;
    }

    if jj == 70 {
        /* only occurs if the last char of string was a quote */
        outstr[69] = 0;
    } else {
        outstr[jj] = bb(b'\''); /* append closing quote character */
        outstr[jj + 1] = 0; /* terminate the string */
    }

    *status
}

/*--------------------------------------------------------------------------*/
/// convert float value to a null-terminated F format string
pub(crate) fn ffr2f(
    fval: f32,           /* I - value to be converted to a string */
    decim: c_int,        /* I - number of decimal places to display */
    cval: &mut [c_char], /* O - character string representation of the value */
    status: &mut c_int,  /* IO - error status */
) -> c_int {
    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    cval[0] = 0;

    if decim < 0 {
        ffpmsg_str("Error in ffr2f:  no. of decimal places < 0");
        *status = BAD_DECIM;
        return *status;
    }

    if int_snprintf!(cval, FLEN_VALUE, "{:.*}", decim as usize, f64::from(fval),) < 0 {
        ffpmsg_str("Error in ffr2f converting float to string");
        *status = BAD_F2C;
    }

    /* replace comma with a period (e.g. in French locale) */
    let cptr = strchr_safe(cval, bb(b','));
    if let Some(cptr) = cptr {
        cval[cptr] = bb(b'.');
    }

    /* test if output string is 'NaN', 'INDEF', or 'INF' */
    if strchr_safe(cval, bb(b'N')).is_some() {
        ffpmsg_str("Error in ffr2f: float value is a NaN or INDEF");
        *status = BAD_F2C;
    }

    *status
}

/*--------------------------------------------------------------------------*/
/// convert float value to a null-terminated exponential format string
pub(crate) fn ffr2e(
    fval: f32,           /* I - value to be converted to a string */
    decim: c_int,        /* I - number of decimal places to display */
    cval: &mut [c_char], /* O - character string representation of the value */
    status: &mut c_int,  /* IO - error status */
) -> c_int {
    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    cval[0] = 0;

    if decim < 0 {
        /* use G format if decim is negative */
        // if int_snprintf!(
        //     cval,
        //     FLEN_VALUE,
        //     "{:.*}",
        //     (-decim) as usize,
        //     fval as f64, // GPoint(fval as f64),
        // ) < 0
        if snprintf_f64_decim(
            cval,
            FLEN_VALUE,
            cast_slice(c"%.*G".to_bytes_with_nul()),
            -decim,
            f64::from(fval),
        ) < 0
        {
            ffpmsg_str("Error in ffr2e converting float to string");
            *status = BAD_F2C;
        } else {
            /* test if E format was used, and there is no displayed decimal */
            if strchr_safe(cval, bb(b'.')).is_none()
                && strchr_safe(cval, bb(b',')).is_none()
                && strchr_safe(cval, bb(b'E')).is_some()
            {
                /* reformat value with a decimal point and single zero */
                if int_snprintf!(cval, FLEN_VALUE, "{:.1E}", f64::from(fval),) < 0 {
                    ffpmsg_str("Error in ffr2e converting float to string");
                    *status = BAD_F2C;
                }

                /* convert French locale comma to a decimal point.*/
                let comma = strchr_safe(cval, bb(b','));
                if let Some(comma) = comma {
                    cval[comma] = bb(b'.');
                }

                return *status;
            }
        }
    } else if int_snprintf!(
        cval,
        FLEN_VALUE,
        "{}",
        fmt_f64(f64::from(fval), decim as usize, 2)
    ) < 0
    {
        ffpmsg_str("Error in ffr2e converting float to string");
        *status = BAD_F2C;
    }

    if *status <= 0 {
        /* replace comma with a period (e.g. in French locale) */
        let cptr = strchr_safe(cval, bb(b','));
        if let Some(cptr) = cptr {
            cval[cptr] = bb(b'.');
        }

        /* test if output string is 'NaN', 'INDEF', or 'INF' */
        if strchr_safe(cval, bb(b'N')).is_some() {
            ffpmsg_str("Error in ffr2e: float value is a NaN or INDEF");
            *status = BAD_F2C;
        } else if strchr_safe(cval, bb(b'.')).is_none()
            && strchr_safe(cval, bb(b'E')).is_none()
            && strlen_safe(cval) < FLEN_VALUE - 1
        {
            /* add decimal point if necessary to distinquish from integer */
            strcat_safe(cval, cs!(c"."));
        }
    }

    *status
}

/*--------------------------------------------------------------------------*/
/// convert double value to a null-terminated F format string
pub(crate) fn ffd2f(
    dval: f64,           /* I - value to be converted to a string */
    decim: c_int,        /* I - number of decimal places to display */
    cval: &mut [c_char], /* O - character string representation of the value */
    status: &mut c_int,  /* IO - error status */
) -> c_int {
    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }
    cval[0] = 0;

    if decim < 0 {
        ffpmsg_str("Error in ffd2f:  no. of decimal places < 0");
        *status = BAD_DECIM;
        return *status;
    }

    if int_snprintf!(cval, FLEN_VALUE, "{:.*}", decim as usize, dval,) < 0 {
        ffpmsg_str("Error in ffd2f converting double to string");
        *status = BAD_F2C;
    }

    /* replace comma with a period (e.g. in French locale) */
    let cptr = strchr_safe(cval, bb(b','));
    if let Some(cptr) = cptr {
        cval[cptr] = bb(b'.');
    }

    /* test if output string is 'NaN', 'INDEF', or 'INF' */
    if strchr_safe(cval, bb(b'N')).is_some() {
        ffpmsg_str("Error in ffd2f: double value is a NaN or INDEF");
        *status = BAD_F2C;
    }

    *status
}

/*--------------------------------------------------------------------------*/
/// convert double value to a null-terminated exponential format string.
pub(crate) fn ffd2e(
    dval: f64,           /* I - value to be converted to a string */
    decim: c_int,        /* I - number of decimal places to display */
    cval: &mut [c_char], /* O - character string representation of the value */
    status: &mut c_int,  /* IO - error status */
) -> c_int {
    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    cval[0] = 0;

    if decim < 0 {
        /* use G format if decim is negative */
        // if int_snprintf!(cval, FLEN_VALUE, "{:.*}", (-decim) as usize, dval /* GPoint(dval) */) < 0 {
        if snprintf_f64_decim(
            cval,
            FLEN_VALUE,
            cast_slice(c"%.*G".to_bytes_with_nul()),
            -decim,
            dval,
        ) < 0
        {
            ffpmsg_str("Error in ffd2e converting float to string");
            *status = BAD_F2C;
        } else {
            /* test if E format was used, and there is no displayed decimal */
            if strchr_safe(cval, bb(b'.')).is_none()
                && strchr_safe(cval, bb(b',')).is_none()
                && strchr_safe(cval, bb(b'E')).is_some()
            {
                /* reformat value with a decimal point and single zero */
                if int_snprintf!(cval, FLEN_VALUE, "{:.1E}", dval) < 0 {
                    ffpmsg_str("Error in ffd2e converting float to string");
                    *status = BAD_F2C;
                }

                /* convert French locale comma to a decimal point.*/
                let comma = strchr_safe(cval, bb(b','));
                if let Some(comma) = comma {
                    cval[comma] = bb(b'.');
                }

                return *status;
            }
        }
    } else if int_snprintf!(cval, FLEN_VALUE, "{}", fmt_f64(dval, decim as usize, 2)) < 0 {
        ffpmsg_str("Error in ffd2e converting float to string");
        *status = BAD_F2C;
    }

    if *status <= 0 {
        /* replace comma with a period (e.g. in French locale) */
        let cptr = strchr_safe(cval, bb(b','));
        if let Some(cptr) = cptr {
            cval[cptr] = bb(b'.');
        }

        /* test if output string is 'NaN', 'INDEF', or 'INF' */
        if strchr_safe(cval, bb(b'N')).is_some() {
            ffpmsg_str("Error in ffd2e: double value is a NaN or INDEF");
            *status = BAD_F2C;
        } else if strchr_safe(cval, bb(b'.')).is_none()
            && strchr_safe(cval, bb(b'E')).is_none()
            && strlen_safe(cval) < FLEN_VALUE - 1
        {
            /* add decimal point if necessary to distinquish from integer */
            strcat_safe(cval, cs!(c"."));
        }
    }

    *status
}

/*--------------------------------------------------------------------------*/
/// Verify that the specified date is valid
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffverifydate(
    year: c_int,        /* I - year */
    month: c_int,       /* I - month (1-12) */
    day: c_int,         /* I - day (1-31) */
    status: *mut c_int, /* IO - error status */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        ffverifydate_safe(year, month, day, status)
    }
}

/// Verify that the specified date is valid (safe version)
pub fn ffverifydate_safe(
    year: c_int,        /* I - year */
    month: c_int,       /* I - month (1-12) */
    day: c_int,         /* I - day (1-31) */
    status: &mut c_int, /* IO - error status */
) -> c_int {
    const NDAYS: [c_int; 13] = [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
    let mut errmsg: [c_char; FLEN_ERRMSG] = [0; FLEN_ERRMSG];

    if year < 0 || year > 9999 {
        int_snprintf!(
            errmsg,
            FLEN_ERRMSG,
            "input year value = {} is out of range 0 - 9999",
            year
        );
        ffpmsg_slice(&errmsg);
        *status = BAD_DATE;
        return *status;
    } else if month < 1 || month > 12 {
        int_snprintf!(
            errmsg,
            FLEN_ERRMSG,
            "input month value = {} is out of range 1 - 12",
            month
        );
        ffpmsg_slice(&errmsg);
        *status = BAD_DATE;
        return *status;
    }

    if NDAYS[month as usize] == 31 {
        if day < 1 || day > 31 {
            int_snprintf!(
                errmsg,
                FLEN_ERRMSG,
                "input day value = {} is out of range 1 - 31 for month {}",
                day,
                month
            );
            ffpmsg_slice(&errmsg);
            *status = BAD_DATE;
            return *status;
        }
    } else if NDAYS[month as usize] == 30 {
        if day < 1 || day > 30 {
            int_snprintf!(
                errmsg,
                FLEN_ERRMSG,
                "input day value = {} is out of range 1 - 30 for month {}",
                day,
                month
            );
            ffpmsg_slice(&errmsg);
            *status = BAD_DATE;
            return *status;
        }
    } else if day < 1 || day > 28 {
        if day == 29 {
            /* year is a leap year if it is divisible by 4 but not by 100,
               except years divisible by 400 are leap years
            */
            if (year % 4 == 0 && year % 100 != 0) || year % 400 == 0 {
                return *status;
            }

            int_snprintf!(
                errmsg,
                FLEN_ERRMSG,
                "input day value = {} is out of range 1 - 28 for February {} (not leap year)",
                day,
                year
            );
            ffpmsg_slice(&errmsg);
        } else {
            int_snprintf!(
                errmsg,
                FLEN_ERRMSG,
                "input day value = {} is out of range 1 - 28 (or 29) for February",
                day
            );
            ffpmsg_slice(&errmsg);
        }

        *status = BAD_DATE;
        return *status;
    }
    *status
}

/*--------------------------------------------------------------------------*/
/// Write (put) an indexed array of keywords with index numbers between
/// NSTART and (NSTART + NKEY -1) inclusive.  Write integer keywords
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffpknjj(
    fptr: *mut fitsfile,        /* I - FITS file pointer */
    keyroot: *const c_char,     /* I - root name of keywords */
    nstart: c_int,              /* I - starting index number */
    nkey: c_int,                /* I - number of keywords to write */
    value: *const LONGLONG,     /* I - array of keyword values */
    comm: *const *const c_char, /* I - array of keyword comments */
    status: *mut c_int,         /* IO - error status */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let value = slice::from_raw_parts(value, nkey as usize);

        let comm = match comm.is_null() {
            true => None,
            false => Some(slice::from_raw_parts(comm, nkey as usize)),
        };

        let mut v_comm = Vec::new();

        if let Some(comm) = comm {
            for item in comm {
                let comm_item = slice::from_raw_parts(*item, FLEN_COMMENT);
                v_comm.push(comm_item);
            }
        }

        raw_to_slice!(keyroot);

        ffpknjj_safe(
            fptr,
            keyroot,
            nstart,
            nkey,
            value,
            if comm.is_some() { Some(&v_comm) } else { None },
            status,
        )
    }
}

/*--------------------------------------------------------------------------*/
/// Write (put) an indexed array of keywords with index numbers between
/// NSTART and (NSTART + NKEY -1) inclusive.  Write integer keywords
pub fn ffpknjj_safe(
    fptr: &mut fitsfile,        /* I - FITS file pointer */
    keyroot: &[c_char],         /* I - root name of keywords */
    nstart: c_int,              /* I - starting index number */
    nkey: c_int,                /* I - number of keywords to write */
    value: &[LONGLONG],         /* I - array of keyword values */
    comm: Option<&[&[c_char]]>, /* I - array of keyword comments */
    status: &mut c_int,         /* IO - error status */
) -> c_int {
    let mut keyname: [c_char; FLEN_KEYWORD] = [0; FLEN_KEYWORD];
    let mut tcomment: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];
    let mut ii: c_int;
    let mut jj: c_int;
    let mut repeat: c_int;
    let mut len: usize;

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }
    /* check if first comment string is to be repeated for all the keywords */
    /* by looking to see if the last non-blank character is a '&' char      */

    repeat = 0;

    if let Some(comm) = comm {
        len = strlen_safe(comm[0]);

        while len > 0 && comm[0][len - 1] == bb(b' ') {
            len -= 1; /* ignore trailing blanks */
        }

        if len > 0 && comm[0][len - 1] == bb(b'&') {
            len = cmp::min(len, FLEN_COMMENT);
            tcomment[0] = 0;
            strncat_safe(&mut tcomment, comm[0], len - 1); /* don't copy the final '&' char */
            repeat = 1;
        }
    } else {
        repeat = 1;
        tcomment[0] = 0;
    }

    let mut jj = nstart;
    for ii in 0..nkey as usize {
        ffkeyn_safe(keyroot, jj, &mut keyname, status);
        if repeat != 0 {
            ffpkyj_safe(fptr, &keyname, value[ii], Some(&tcomment), status);
        } else {
            ffpkyj_safe(fptr, &keyname, value[ii], Some(comm.unwrap()[ii]), status);
        }

        if *status > 0 {
            return *status;
        }
        jj += 1;
    }
    *status
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::KeywordDatatype;
    use crate::aliases::rust_api::*;
    use crate::fitsio::{
        BAD_COL_NUM, BAD_DATATYPE, BAD_DATE, BAD_DIMEN, BAD_F2C, BAD_TDIM, BINARY_TBL, BYTE_IMG,
        FILE_NOT_OPENED, FLEN_VALUE, LONGLONG, LONGLONG_IMG, NOT_BTABLE, READONLY, SBYTE_IMG,
        SHORT_IMG, ULONG_IMG, ULONGLONG_IMG, USHORT_IMG, fitsfile,
    };
    use crate::helpers::testhelpers::{from_buf, to_buf, with_temp_file};
    use libc::{c_char, c_double, c_int, c_long};

    /// Make a NUL-terminated `Vec<c_char>` from a `&str`.
    fn cc(s: &str) -> Vec<c_char> {
        let mut v: Vec<c_char> = s.bytes().map(|b| b as c_char).collect();
        v.push(0);
        v
    }

    #[test]
    fn test_ffverifydate_safe() {
        let mut status = 0;

        // Test valid dates for each month
        assert_eq!(ffverifydate_safe(2024, 1, 15, &mut status), 0);
        assert_eq!(status, 0);

        assert_eq!(ffverifydate_safe(2024, 3, 31, &mut status), 0);
        assert_eq!(status, 0);

        assert_eq!(ffverifydate_safe(2024, 4, 30, &mut status), 0);
        assert_eq!(status, 0);

        // Test leap year - 2024 is a leap year
        status = 0;
        assert_eq!(ffverifydate_safe(2024, 2, 29, &mut status), 0);
        assert_eq!(status, 0);

        // Test non-leap year - 2023 is not a leap year
        status = 0;
        assert_eq!(ffverifydate_safe(2023, 2, 29, &mut status), BAD_DATE);
        assert_eq!(status, BAD_DATE);

        // Test leap year rules - 2000 is a leap year (divisible by 400)
        status = 0;
        assert_eq!(ffverifydate_safe(2000, 2, 29, &mut status), 0);
        assert_eq!(status, 0);

        // Test leap year rules - 1900 is not a leap year (divisible by 100 but not 400)
        status = 0;
        assert_eq!(ffverifydate_safe(1900, 2, 29, &mut status), BAD_DATE);
        assert_eq!(status, BAD_DATE);

        // Test invalid year
        status = 0;
        assert_eq!(ffverifydate_safe(-1, 1, 1, &mut status), BAD_DATE);
        assert_eq!(status, BAD_DATE);

        status = 0;
        assert_eq!(ffverifydate_safe(10000, 1, 1, &mut status), BAD_DATE);
        assert_eq!(status, BAD_DATE);

        // Test invalid month
        status = 0;
        assert_eq!(ffverifydate_safe(2024, 0, 1, &mut status), BAD_DATE);
        assert_eq!(status, BAD_DATE);

        status = 0;
        assert_eq!(ffverifydate_safe(2024, 13, 1, &mut status), BAD_DATE);
        assert_eq!(status, BAD_DATE);

        // Test invalid day for 31-day month
        status = 0;
        assert_eq!(ffverifydate_safe(2024, 1, 0, &mut status), BAD_DATE);
        assert_eq!(status, BAD_DATE);

        status = 0;
        assert_eq!(ffverifydate_safe(2024, 1, 32, &mut status), BAD_DATE);
        assert_eq!(status, BAD_DATE);

        // Test invalid day for 30-day month
        status = 0;
        assert_eq!(ffverifydate_safe(2024, 4, 31, &mut status), BAD_DATE);
        assert_eq!(status, BAD_DATE);

        // Test invalid day for February (non-leap year)
        status = 0;
        assert_eq!(ffverifydate_safe(2023, 2, 30, &mut status), BAD_DATE);
        assert_eq!(status, BAD_DATE);
    }

    #[test]
    fn test_ffs2dt_safe() {
        let mut status = 0;
        let mut year = 0;
        let mut month = 0;
        let mut day = 0;

        // Test old format: "dd/mm/yy"
        let date_old = cs!(c"15/03/99");
        status = 0;
        assert_eq!(
            ffs2dt_safe(
                Some(date_old),
                Some(&mut year),
                Some(&mut month),
                Some(&mut day),
                &mut status
            ),
            0
        );
        assert_eq!(status, 0);
        assert_eq!(year, 1999);
        assert_eq!(month, 3);
        assert_eq!(day, 15);

        // Test new format: "yyyy-mm-dd"
        let date_new = cs!(c"2024-03-15");
        status = 0;
        assert_eq!(
            ffs2dt_safe(
                Some(date_new),
                Some(&mut year),
                Some(&mut month),
                Some(&mut day),
                &mut status
            ),
            0
        );
        assert_eq!(status, 0);
        assert_eq!(year, 2024);
        assert_eq!(month, 3);
        assert_eq!(day, 15);

        // Test new format with time extension: "yyyy-mm-ddThh:mm:ss"
        let date_time = cs!(c"2024-03-15T12:30:45");
        status = 0;
        assert_eq!(
            ffs2dt_safe(
                Some(date_time),
                Some(&mut year),
                Some(&mut month),
                Some(&mut day),
                &mut status
            ),
            0
        );
        assert_eq!(status, 0);
        assert_eq!(year, 2024);
        assert_eq!(month, 3);
        assert_eq!(day, 15);

        // Test invalid old format
        let invalid_old = cs!(c"1a/03/99");
        status = 0;
        assert_eq!(
            ffs2dt_safe(
                Some(invalid_old),
                Some(&mut year),
                Some(&mut month),
                Some(&mut day),
                &mut status
            ),
            BAD_DATE
        );
        assert_eq!(status, BAD_DATE);

        // Test invalid new format
        let invalid_new = cs!(c"202a-03-15");
        status = 0;
        assert_eq!(
            ffs2dt_safe(
                Some(invalid_new),
                Some(&mut year),
                Some(&mut month),
                Some(&mut day),
                &mut status
            ),
            BAD_DATE
        );
        assert_eq!(status, BAD_DATE);

        // Test invalid format with wrong separator after date
        let invalid_sep = cs!(c"2024-03-15X12:30:45");
        status = 0;
        assert_eq!(
            ffs2dt_safe(
                Some(invalid_sep),
                Some(&mut year),
                Some(&mut month),
                Some(&mut day),
                &mut status
            ),
            BAD_DATE
        );
        assert_eq!(status, BAD_DATE);

        // Test empty string
        let empty: &[c_char] = &[];
        status = 0;
        assert_eq!(
            ffs2dt_safe(
                Some(empty),
                Some(&mut year),
                Some(&mut month),
                Some(&mut day),
                &mut status
            ),
            BAD_DATE
        );

        // Test None (null pointer)
        status = 0;
        assert_eq!(
            ffs2dt_safe(
                None,
                Some(&mut year),
                Some(&mut month),
                Some(&mut day),
                &mut status
            ),
            BAD_DATE
        );
        assert_eq!(status, BAD_DATE);

        // Test date validation through ffverifydate_safe
        let invalid_date = cs!(c"2023-02-29"); // Not a leap year
        status = 0;
        assert_eq!(
            ffs2dt_safe(
                Some(invalid_date),
                Some(&mut year),
                Some(&mut month),
                Some(&mut day),
                &mut status
            ),
            BAD_DATE
        );
        assert_eq!(status, BAD_DATE);
    }

    #[test]
    fn test_ffs2tm_safe() {
        let mut status = 0;
        let mut year = 0;
        let mut month = 0;
        let mut day = 0;
        let mut hour = 0;
        let mut minute = 0;
        let mut second = 0.0;

        // Test full datetime: "yyyy-mm-ddThh:mm:ss.ddd"
        let datetime = cs!(c"2024-03-15T14:30:45.123");
        status = 0;
        assert_eq!(
            ffs2tm_safe(
                Some(datetime),
                Some(&mut year),
                Some(&mut month),
                Some(&mut day),
                Some(&mut hour),
                Some(&mut minute),
                Some(&mut second),
                &mut status
            ),
            0
        );
        assert_eq!(status, 0);
        assert_eq!(year, 2024);
        assert_eq!(month, 3);
        assert_eq!(day, 15);
        assert_eq!(hour, 14);
        assert_eq!(minute, 30);
        assert!(second >= 45.123 - 0.001 && second <= 45.123 + 0.001);

        // Test date only: "yyyy-mm-dd"
        let date_only = cs!(c"2024-03-15");
        status = 0;
        year = 0;
        month = 0;
        day = 0;
        hour = 0;
        minute = 0;
        second = 0.0;
        assert_eq!(
            ffs2tm_safe(
                Some(date_only),
                Some(&mut year),
                Some(&mut month),
                Some(&mut day),
                Some(&mut hour),
                Some(&mut minute),
                Some(&mut second),
                &mut status
            ),
            0
        );
        assert_eq!(status, 0);
        assert_eq!(year, 2024);
        assert_eq!(month, 3);
        assert_eq!(day, 15);
        assert_eq!(hour, 0);
        assert_eq!(minute, 0);
        assert_eq!(second, 0.0);

        // Test old date format: "dd/mm/yy"
        let old_date = cs!(c"15/03/99");
        status = 0;
        year = 0;
        month = 0;
        day = 0;
        hour = 0;
        minute = 0;
        second = 0.0;
        assert_eq!(
            ffs2tm_safe(
                Some(old_date),
                Some(&mut year),
                Some(&mut month),
                Some(&mut day),
                Some(&mut hour),
                Some(&mut minute),
                Some(&mut second),
                &mut status
            ),
            0
        );
        assert_eq!(status, 0);
        assert_eq!(year, 1999);
        assert_eq!(month, 3);
        assert_eq!(day, 15);

        // Test time only: "hh:mm:ss.ddd"
        let time_only = cs!(c"14:30:45.678");
        status = 0;
        year = 0;
        month = 0;
        day = 0;
        hour = 0;
        minute = 0;
        second = 0.0;
        assert_eq!(
            ffs2tm_safe(
                Some(time_only),
                Some(&mut year),
                Some(&mut month),
                Some(&mut day),
                Some(&mut hour),
                Some(&mut minute),
                Some(&mut second),
                &mut status
            ),
            0
        );
        assert_eq!(status, 0);
        assert_eq!(year, 0);
        assert_eq!(month, 0);
        assert_eq!(day, 0);
        assert_eq!(hour, 14);
        assert_eq!(minute, 30);
        assert!(second >= 45.678 - 0.001 && second <= 45.678 + 0.001);

        // Test invalid hour
        let invalid_hour = cs!(c"24:30:45");
        status = 0;
        assert_eq!(
            ffs2tm_safe(
                Some(invalid_hour),
                Some(&mut year),
                Some(&mut month),
                Some(&mut day),
                Some(&mut hour),
                Some(&mut minute),
                Some(&mut second),
                &mut status
            ),
            BAD_DATE
        );
        assert_eq!(status, BAD_DATE);

        // Test invalid minute
        let invalid_minute = cs!(c"14:60:45");
        status = 0;
        assert_eq!(
            ffs2tm_safe(
                Some(invalid_minute),
                Some(&mut year),
                Some(&mut month),
                Some(&mut day),
                Some(&mut hour),
                Some(&mut minute),
                Some(&mut second),
                &mut status
            ),
            BAD_DATE
        );
        assert_eq!(status, BAD_DATE);

        // Test invalid second (>= 61)
        let invalid_second = cs!(c"14:30:61.0");
        status = 0;
        assert_eq!(
            ffs2tm_safe(
                Some(invalid_second),
                Some(&mut year),
                Some(&mut month),
                Some(&mut day),
                Some(&mut hour),
                Some(&mut minute),
                Some(&mut second),
                &mut status
            ),
            BAD_DATE
        );
        assert_eq!(status, BAD_DATE);

        // Test valid leap second (60.999...)
        let leap_second = cs!(c"14:30:60.5");
        status = 0;
        assert_eq!(
            ffs2tm_safe(
                Some(leap_second),
                Some(&mut year),
                Some(&mut month),
                Some(&mut day),
                Some(&mut hour),
                Some(&mut minute),
                Some(&mut second),
                &mut status
            ),
            0
        );
        assert_eq!(status, 0);
        assert!(second >= 60.5 - 0.001 && second <= 60.5 + 0.001);

        // Test invalid format
        let invalid_format = cs!(c"not-a-date");
        status = 0;
        assert_eq!(
            ffs2tm_safe(
                Some(invalid_format),
                Some(&mut year),
                Some(&mut month),
                Some(&mut day),
                Some(&mut hour),
                Some(&mut minute),
                Some(&mut second),
                &mut status
            ),
            BAD_DATE
        );
        assert_eq!(status, BAD_DATE);

        // Test empty string
        let empty: &[c_char] = &[];
        status = 0;
        assert_eq!(
            ffs2tm_safe(
                Some(empty),
                Some(&mut year),
                Some(&mut month),
                Some(&mut day),
                Some(&mut hour),
                Some(&mut minute),
                Some(&mut second),
                &mut status
            ),
            BAD_DATE
        );

        // Test None (null pointer)
        status = 0;
        assert_eq!(
            ffs2tm_safe(
                None,
                Some(&mut year),
                Some(&mut month),
                Some(&mut day),
                Some(&mut hour),
                Some(&mut minute),
                Some(&mut second),
                &mut status
            ),
            BAD_DATE
        );
        assert_eq!(status, BAD_DATE);
    }

    #[test]
    fn test_ffdt2s_safe() {
        let mut status = 0;
        let mut datestr = [0; 11];

        // Test new format (YYYY-MM-DD) for year 2024
        status = 0;
        datestr = [0; 11];
        assert_eq!(ffdt2s_safe(2024, 3, 15, &mut datestr, &mut status), 0);
        assert_eq!(status, 0);
        let result_str = CStr::from_bytes_until_nul(cast_slice(&datestr))
            .unwrap()
            .to_str()
            .unwrap();
        assert_eq!(result_str, "2024-03-15");

        // Test old format (dd/mm/yy) for year 1950
        status = 0;
        datestr = [0; 11];
        assert_eq!(ffdt2s_safe(1950, 12, 25, &mut datestr, &mut status), 0);
        assert_eq!(status, 0);
        let result_str = CStr::from_bytes_until_nul(cast_slice(&datestr))
            .unwrap()
            .to_str()
            .unwrap();
        assert_eq!(result_str, "25/12/50");

        // Test edge case: year 1900 (should use old format)
        status = 0;
        datestr = [0; 11];
        assert_eq!(ffdt2s_safe(1900, 1, 1, &mut datestr, &mut status), 0);
        assert_eq!(status, 0);
        let result_str = CStr::from_bytes_until_nul(cast_slice(&datestr))
            .unwrap()
            .to_str()
            .unwrap();
        assert_eq!(result_str, "01/01/00");

        // Test edge case: year 1998 (should use old format)
        status = 0;
        datestr = [0; 11];
        assert_eq!(ffdt2s_safe(1998, 12, 31, &mut datestr, &mut status), 0);
        assert_eq!(status, 0);
        let result_str = CStr::from_bytes_until_nul(cast_slice(&datestr))
            .unwrap()
            .to_str()
            .unwrap();
        assert_eq!(result_str, "31/12/98");

        // Test edge case: year 1999 (should use new format)
        status = 0;
        datestr = [0; 11];
        assert_eq!(ffdt2s_safe(1999, 6, 15, &mut datestr, &mut status), 0);
        assert_eq!(status, 0);
        let result_str = CStr::from_bytes_until_nul(cast_slice(&datestr))
            .unwrap()
            .to_str()
            .unwrap();
        assert_eq!(result_str, "1999-06-15");

        // Test year 0 (should use new format)
        status = 0;
        datestr = [0; 11];
        assert_eq!(ffdt2s_safe(0, 1, 1, &mut datestr, &mut status), 0);
        assert_eq!(status, 0);
        let result_str = CStr::from_bytes_until_nul(cast_slice(&datestr))
            .unwrap()
            .to_str()
            .unwrap();
        assert_eq!(result_str, "0000-01-01");

        // Test year 9999 (should use new format)
        status = 0;
        datestr = [0; 11];
        assert_eq!(ffdt2s_safe(9999, 12, 31, &mut datestr, &mut status), 0);
        assert_eq!(status, 0);
        let result_str = CStr::from_bytes_until_nul(cast_slice(&datestr))
            .unwrap()
            .to_str()
            .unwrap();
        assert_eq!(result_str, "9999-12-31");

        // Test invalid date: month > 12
        status = 0;
        datestr = [0; 11];
        assert_eq!(
            ffdt2s_safe(2024, 13, 15, &mut datestr, &mut status),
            BAD_DATE
        );
        assert_eq!(status, BAD_DATE);

        // Test invalid date: day > 31
        status = 0;
        datestr = [0; 11];
        assert_eq!(
            ffdt2s_safe(2024, 3, 32, &mut datestr, &mut status),
            BAD_DATE
        );
        assert_eq!(status, BAD_DATE);

        // Test invalid date: month = 0
        status = 0;
        datestr = [0; 11];
        assert_eq!(
            ffdt2s_safe(2024, 0, 15, &mut datestr, &mut status),
            BAD_DATE
        );
        assert_eq!(status, BAD_DATE);

        // Test invalid date: day = 0
        status = 0;
        datestr = [0; 11];
        assert_eq!(ffdt2s_safe(2024, 3, 0, &mut datestr, &mut status), BAD_DATE);
        assert_eq!(status, BAD_DATE);

        // Test invalid date: year > 9999
        status = 0;
        datestr = [0; 11];
        assert_eq!(
            ffdt2s_safe(10000, 3, 15, &mut datestr, &mut status),
            BAD_DATE
        );
        assert_eq!(status, BAD_DATE);

        // Test leap year February 29
        status = 0;
        datestr = [0; 11];
        assert_eq!(ffdt2s_safe(2024, 2, 29, &mut datestr, &mut status), 0);
        assert_eq!(status, 0);
        let result_str = CStr::from_bytes_until_nul(cast_slice(&datestr))
            .unwrap()
            .to_str()
            .unwrap();
        assert_eq!(result_str, "2024-02-29");

        // Test single digit values with padding
        status = 0;
        datestr = [0; 11];
        assert_eq!(ffdt2s_safe(2024, 1, 5, &mut datestr, &mut status), 0);
        assert_eq!(status, 0);
        let result_str = CStr::from_bytes_until_nul(cast_slice(&datestr))
            .unwrap()
            .to_str()
            .unwrap();
        assert_eq!(result_str, "2024-01-05");

        // Test status inheritance
        status = 123; // Pre-existing error status
        datestr = [0; 11];
        assert_eq!(ffdt2s_safe(2024, 3, 15, &mut datestr, &mut status), 123);
        assert_eq!(status, 123);
    }

    #[test]
    fn test_putkey_string() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            assert_eq!(status, 0, "ffinit failed");
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            assert_eq!(status, 0, "ffphps failed");
            fits_write_key_str(
                f.as_deref_mut().unwrap(),
                &cc("STRKEY"),
                &cc("test value"),
                Some(&cc("string keyword")),
                &mut status,
            );
            assert_eq!(status, 0, "ffpkys failed");
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            assert_eq!(status, 0, "ffopen failed");
            let mut sval = [0 as c_char; FLEN_VALUE];
            fits_read_key_str(
                f.as_deref_mut().unwrap(),
                &cc("STRKEY"),
                &mut sval,
                None,
                &mut status,
            );
            assert_eq!(status, 0, "ffgkys failed");
            assert_eq!(from_buf(&sval), "test value");
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_putkey_long() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            assert_eq!(status, 0, "ffinit failed");
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            assert_eq!(status, 0, "ffphps failed");
            fits_write_key_lng(
                f.as_deref_mut().unwrap(),
                &cc("LONGKEY"),
                123456789,
                Some(&cc("long keyword")),
                &mut status,
            );
            assert_eq!(status, 0, "ffpkyj failed");
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            assert_eq!(status, 0, "ffopen failed");
            let mut lval: c_long = 0;
            fits_read_key(
                f.as_deref_mut().unwrap(),
                crate::KeywordDatatypeMut::TLONG(&mut lval),
                &cc("LONGKEY"),
                None,
                &mut status,
            );
            assert_eq!(status, 0, "ffgky failed");
            assert_eq!(lval, 123456789);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_putkey_double() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            assert_eq!(status, 0, "setup failed");
            fits_write_key_dbl(
                f.as_deref_mut().unwrap(),
                &cc("DBLKEY"),
                3.14159265358979,
                15,
                Some(&cc("double keyword")),
                &mut status,
            );
            assert_eq!(status, 0, "ffpkyd failed");
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut dval: c_double = 0.0;
            fits_read_key(
                f.as_deref_mut().unwrap(),
                crate::KeywordDatatypeMut::TDOUBLE(&mut dval),
                &cc("DBLKEY"),
                None,
                &mut status,
            );
            assert_eq!(status, 0, "ffgky failed");
            assert!(dval >= 3.14 && dval <= 3.15);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_putkey_logical() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            assert_eq!(status, 0, "setup failed");
            fits_write_key_log(
                f.as_deref_mut().unwrap(),
                &cc("BOOLKEY"),
                1,
                Some(&cc("logical keyword")),
                &mut status,
            );
            assert_eq!(status, 0, "ffpkyl failed");
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut lval: c_int = 0;
            fits_read_key_log(
                f.as_deref_mut().unwrap(),
                &cc("BOOLKEY"),
                &mut lval,
                None,
                &mut status,
            );
            assert_eq!(status, 0, "ffgkyl failed");
            assert_eq!(lval, 1);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_putkey_float() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            assert_eq!(status, 0, "setup failed");
            fits_write_key_flt(
                f.as_deref_mut().unwrap(),
                &cc("FLTKEY"),
                2.71828,
                6,
                Some(&cc("float keyword")),
                &mut status,
            );
            assert_eq!(status, 0, "ffpkye failed");
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut fval: f32 = 0.0;
            fits_read_key(
                f.as_deref_mut().unwrap(),
                crate::KeywordDatatypeMut::TFLOAT(&mut fval),
                &cc("FLTKEY"),
                None,
                &mut status,
            );
            assert_eq!(status, 0, "ffgky failed");
            assert!(fval >= 2.71 && fval <= 2.72);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_putkey_generic_types() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);

            let sbval: c_schar = -42;
            let usval: u16 = 12345;
            let uival: u32 = 3000000000;
            let llval: LONGLONG = 9876543210;
            let fval: f32 = 1.5;
            let dval: f64 = 2.5;
            let cval: [f32; 2] = [1.0, 2.0];
            let mval: [f64; 2] = [3.0, 4.0];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            assert_eq!(status, 0, "setup failed");
            let g = f.as_deref_mut().unwrap();

            fits_write_key(
                g,
                KeywordDatatype::TSBYTE(&sbval),
                &cc("SBKEY"),
                Some(&cc("signed byte")),
                &mut status,
            );
            fits_write_key(
                g,
                KeywordDatatype::TUSHORT(&usval),
                &cc("USKEY"),
                Some(&cc("unsigned short")),
                &mut status,
            );
            fits_write_key(
                g,
                KeywordDatatype::TUINT(&uival),
                &cc("UIKEY"),
                Some(&cc("unsigned int")),
                &mut status,
            );
            fits_write_key(
                g,
                KeywordDatatype::TLONGLONG(&llval),
                &cc("LLKEY"),
                Some(&cc("longlong")),
                &mut status,
            );
            fits_write_key(
                g,
                KeywordDatatype::TFLOAT(&fval),
                &cc("FKEY"),
                Some(&cc("float")),
                &mut status,
            );
            fits_write_key(
                g,
                KeywordDatatype::TDOUBLE(&dval),
                &cc("DKEY"),
                Some(&cc("double")),
                &mut status,
            );
            fits_write_key(
                g,
                KeywordDatatype::TCOMPLEX(&cval),
                &cc("CKEY"),
                Some(&cc("complex")),
                &mut status,
            );
            fits_write_key(
                g,
                KeywordDatatype::TDBLCOMPLEX(&mval),
                &cc("MKEY"),
                Some(&cc("dblcomplex")),
                &mut status,
            );
            assert_eq!(status, 0, "ffpky failed");

            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_putkey_unsigned_long() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);

            let ulval: c_ulong = 4000000000;
            let ullval: u64 = 18000000000000000000;

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            assert_eq!(status, 0, "setup failed");
            let g = f.as_deref_mut().unwrap();

            fits_write_key(
                g,
                KeywordDatatype::TULONG(&ulval),
                &cc("ULKEY"),
                Some(&cc("unsigned long")),
                &mut status,
            );
            fits_write_key(
                g,
                KeywordDatatype::TULONGLONG(&ullval),
                &cc("ULLKEY"),
                Some(&cc("ulonglong")),
                &mut status,
            );
            assert_eq!(status, 0, "ffpky failed");

            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_putkey_null() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            assert_eq!(status, 0, "setup failed");
            fits_write_key_null(
                f.as_deref_mut().unwrap(),
                &cc("NULLKEY"),
                Some(&cc("undefined value")),
                &mut status,
            );
            assert_eq!(status, 0, "ffpkyu failed");
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_putkey_units() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            assert_eq!(status, 0, "setup failed");
            fits_write_key_lng(
                f.as_deref_mut().unwrap(),
                &cc("EXPTIME"),
                300,
                Some(&cc("Exposure time")),
                &mut status,
            );
            fits_write_key_unit(
                f.as_deref_mut().unwrap(),
                &cc("EXPTIME"),
                &cc("seconds"),
                &mut status,
            );
            assert_eq!(status, 0, "ffpunt failed");
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut units = [0 as c_char; FLEN_VALUE];
            fits_read_key_unit(
                f.as_deref_mut().unwrap(),
                &cc("EXPTIME"),
                &mut units,
                &mut status,
            );
            assert_eq!(status, 0, "ffgunt failed");
            assert_eq!(from_buf(&units), "seconds");
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_putkey_fixfmt() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);

            let dval: f64 = 123.456789;
            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            assert_eq!(status, 0, "setup failed");
            fits_write_key_fixdbl(
                f.as_deref_mut().unwrap(),
                &cc("GKEY"),
                dval,
                10,
                Some(&cc("G format")),
                &mut status,
            );
            assert_eq!(status, 0, "ffpkyg failed");
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_putkey_fixfmt_float() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);

            let dval: f64 = 987.654321;
            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            assert_eq!(status, 0, "setup failed");
            // ffpkyf takes f32 value
            fits_write_key_fixflt(
                f.as_deref_mut().unwrap(),
                &cc("FKEY"),
                dval as f32,
                6,
                Some(&cc("F format")),
                &mut status,
            );
            assert_eq!(status, 0, "ffpkyf failed");
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_putkey_complex() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);

            let cval: [f32; 2] = [1.5, -2.5];
            let mval: [f64; 2] = [3.14159, -2.71828];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            assert_eq!(status, 0, "setup failed");
            fits_write_key_cmp(
                f.as_deref_mut().unwrap(),
                &cc("CFLT"),
                &cval,
                6,
                Some(&cc("complex float")),
                &mut status,
            );
            fits_write_key_dblcmp(
                f.as_deref_mut().unwrap(),
                &cc("CDBL"),
                &mval,
                14,
                Some(&cc("complex double")),
                &mut status,
            );
            assert_eq!(status, 0, "ffpkyc/ffpkym failed");
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_putkey_complex_fixed() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);

            let cval: [f32; 2] = [100.0, 200.0];
            let mval: [f64; 2] = [300.0, 400.0];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            assert_eq!(status, 0, "setup failed");
            fits_write_key_fixcmp(
                f.as_deref_mut().unwrap(),
                &cc("FCFLT"),
                &cval,
                4,
                Some(&cc("fixed complex float")),
                &mut status,
            );
            fits_write_key_fixdblcmp(
                f.as_deref_mut().unwrap(),
                &cc("FCDBL"),
                &mval,
                8,
                Some(&cc("fixed complex double")),
                &mut status,
            );
            assert_eq!(status, 0, "ffpkfc/ffpkfm failed");
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_putkey_comment_history() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            assert_eq!(status, 0, "setup failed");
            fits_write_comment(
                f.as_deref_mut().unwrap(),
                &cc("This is a comment keyword"),
                &mut status,
            );
            fits_write_history(
                f.as_deref_mut().unwrap(),
                &cc("This is a history keyword"),
                &mut status,
            );
            assert_eq!(status, 0, "ffpcom/ffphis failed");
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_putkey_date() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            assert_eq!(status, 0, "setup failed");
            fits_write_date(f.as_deref_mut().unwrap(), &mut status);
            assert_eq!(status, 0, "ffpdat failed");
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_update_key() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            fits_write_key_lng(
                f.as_deref_mut().unwrap(),
                &cc("UPDKEY"),
                100,
                Some(&cc("original")),
                &mut status,
            );
            fits_update_key_lng(
                f.as_deref_mut().unwrap(),
                &cc("UPDKEY"),
                200,
                Some(&cc("updated")),
                &mut status,
            );
            assert_eq!(status, 0, "ffukyj failed");
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut lval: c_long = 0;
            fits_read_key(
                f.as_deref_mut().unwrap(),
                crate::KeywordDatatypeMut::TLONG(&mut lval),
                &cc("UPDKEY"),
                None,
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(lval, 200);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_putkey_indexed_arrays() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);

            let svals = [cc("val1"), cc("val2"), cc("val3")];
            let svals_ptrs: Vec<*const c_char> = svals.iter().map(|v| v.as_ptr()).collect();
            let lvals: [c_int; 3] = [1, 0, 1];
            let jvals: [c_long; 3] = [100, 200, 300];
            let jjvals: [LONGLONG; 3] = [1000000000, 2000000000, 3000000000];
            let fvals: [f32; 3] = [1.1, 2.2, 3.3];
            let evals: [f32; 3] = [1.0e10, 2.0e20, 3.0e30];
            let gvals: [f64; 3] = [1.111, 2.222, 3.333];
            let dvals: [f64; 3] = [1.0e100, 2.0e200, 3.0e-100];
            let comm = [cc("comment &"), cc("comment 2"), cc("comment 3")];
            let comm_ptrs: Vec<*const c_char> = comm.iter().map(|v| v.as_ptr()).collect();
            let comm_norep = [cc("comm1"), cc("comm2"), cc("comm3")];
            let comm_norep_ptrs: Vec<*const c_char> =
                comm_norep.iter().map(|v| v.as_ptr()).collect();

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            assert_eq!(status, 0, "setup failed");
            let g = f.as_deref_mut().unwrap();

            // Test with repeated comment (ending with &)
            fits_write_keys_str(g, &cc("SKEY"), 1, 3, &svals_ptrs, &comm_ptrs, &mut status);
            fits_write_keys_log(g, &cc("LKEY"), 1, 3, &lvals, &comm_ptrs, &mut status);
            fits_write_keys_lng(g, &cc("JKEY"), 1, 3, &jvals, &comm_ptrs, &mut status);
            // Test with NULL comment via ffpknjj (no public alias -> call safe)
            ffpknjj_safe(g, &cc("JJKEY"), 1, 3, &jjvals, None, &mut status);
            fits_write_keys_fixflt(g, &cc("FKEY"), 1, 3, &fvals, 2, &comm_ptrs, &mut status);
            fits_write_keys_flt(g, &cc("EKEY"), 1, 3, &evals, 3, &comm_ptrs, &mut status);
            fits_write_keys_fixdbl(g, &cc("GKEY"), 1, 3, &gvals, 4, &comm_ptrs, &mut status);
            fits_write_keys_dbl(g, &cc("DKEY"), 1, 3, &dvals, 5, &comm_ptrs, &mut status);

            // Test with non-repeated comments
            fits_write_keys_log(g, &cc("LKEY2"), 1, 3, &lvals, &comm_norep_ptrs, &mut status);
            fits_write_keys_lng(g, &cc("JKEY2"), 1, 3, &jvals, &comm_norep_ptrs, &mut status);

            assert_eq!(status, 0, "indexed array keywords failed");
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_putkey_long_string() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);

            let longstr = "This is a very long string value that exceeds the normal \
                FITS keyword value length limit of 68 characters and requires continuation";

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            fits_write_key_longstr(
                f.as_deref_mut().unwrap(),
                &cc("LONGSTR"),
                &cc(longstr),
                Some(&cc("long string test")),
                &mut status,
            );
            assert_eq!(status, 0, "ffpkls failed");
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut result: *mut c_char = core::ptr::null_mut();
            fits_read_key_longstr(
                f.as_deref_mut().unwrap(),
                &cc("LONGSTR"),
                &mut result,
                None,
                &mut status,
            );
            assert_eq!(status, 0, "ffgkls failed");
            assert!(!result.is_null());
            let got = unsafe { core::ffi::CStr::from_ptr(result) }
                .to_str()
                .unwrap()
                .to_string();
            assert_eq!(got, longstr);
            unsafe { fits_free_memory(result.cast::<libc::c_void>(), &mut status) };
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_putkey_template() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);

            // Create a template file alongside the fits file.
            let tpl_path = std::path::Path::new(filename).with_file_name("test_putkey.tpl");
            std::fs::write(
                &tpl_path,
                "TPLKEY1 = 'template value' / added by template\n\
                 TPLKEY2 = 42 / integer from template\n\
                 END\n",
            )
            .unwrap();
            let tpl_buf = to_buf(tpl_path.to_str().unwrap());

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            fits_write_key_template(f.as_deref_mut().unwrap(), &tpl_buf, &mut status);
            assert_eq!(status, 0, "ffpktp failed");
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut sval = [0 as c_char; FLEN_VALUE];
            fits_read_key_str(
                f.as_deref_mut().unwrap(),
                &cc("TPLKEY1"),
                &mut sval,
                None,
                &mut status,
            );
            assert_eq!(status, 0, "ffgkys failed");
            assert_eq!(from_buf(&sval), "template value");
            fits_close_file(f.take().unwrap(), &mut status);

            let _ = std::fs::remove_file(&tpl_path);
        });
    }

    /// Mirrors test_putkey_datetime in ~/code/cfitsio/tests/test_putkey.c
    #[test]
    fn test_putkey_datetime() {
        let mut status: c_int = 0;
        let mut yr: c_int = 0;
        let mut mon: c_int = 0;
        let mut day: c_int = 0;
        let mut hr: c_int = 0;
        let mut min: c_int = 0;
        let mut sec: c_double = 0.0;
        let mut datestr = [0 as c_char; 11];
        let mut timestr = [0 as c_char; 30];

        /* ffdt2s - convert date to string (new format) */
        fits_date2str(2024, 1, 15, &mut datestr, &mut status);
        assert_eq!(status, 0);
        assert_eq!(from_buf(&datestr), "2024-01-15");

        /* ffdt2s - old format date (year 1900-1998) */
        fits_date2str(1995, 6, 15, &mut datestr, &mut status);
        assert_eq!(status, 0);
        assert_eq!(from_buf(&datestr), "15/06/95");

        /* ffs2dt - convert string to date */
        fits_str2date(
            Some(&cc("2024-06-20")),
            Some(&mut yr),
            Some(&mut mon),
            Some(&mut day),
            &mut status,
        );
        assert_eq!(status, 0);
        assert_eq!((yr, mon, day), (2024, 6, 20));

        /* old-style date format */
        fits_str2date(
            Some(&cc("15/03/99")),
            Some(&mut yr),
            Some(&mut mon),
            Some(&mut day),
            &mut status,
        );
        assert_eq!(status, 0);
        assert_eq!((yr, mon, day), (1999, 3, 15));

        /* fftm2s - convert date-time to string */
        fits_time2str(2024, 7, 4, 12, 30, 45.5, 3, &mut timestr, &mut status);
        assert_eq!(status, 0);
        assert_eq!(&from_buf(&timestr)[..19], "2024-07-04T12:30:45");

        /* fftm2s with 0 decimals */
        fits_time2str(2024, 8, 15, 10, 20, 30.0, 0, &mut timestr, &mut status);
        assert_eq!(status, 0);
        assert_eq!(&from_buf(&timestr)[..19], "2024-08-15T10:20:30");

        /* fftm2s with negative decimals (date only) */
        fits_time2str(2024, 9, 1, 0, 0, 0.0, -1, &mut timestr, &mut status);
        assert_eq!(status, 0);
        assert_eq!(from_buf(&timestr), "2024-09-01");

        /* fftm2s with time only (year=month=day=0) */
        fits_time2str(0, 0, 0, 14, 30, 45.5, 3, &mut timestr, &mut status);
        assert_eq!(status, 0);
        assert_eq!(&from_buf(&timestr)[..8], "14:30:45");

        /* ffs2tm - convert string to date-time */
        fits_str2time(
            Some(&cc("2024-12-25T08:15:30.123")),
            Some(&mut yr),
            Some(&mut mon),
            Some(&mut day),
            Some(&mut hr),
            Some(&mut min),
            Some(&mut sec),
            &mut status,
        );
        assert_eq!(status, 0);
        assert_eq!((yr, mon, day), (2024, 12, 25));
        assert_eq!((hr, min), (8, 15));
        assert!((sec - 30.123).abs() < 1e-6, "second = {sec}");
    }

    #[test]
    fn test_putkey_tdim() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);

            let naxes: [c_long; 3] = [10, 20, 30];
            let ttype = [Some(cc("DATA"))];
            let ttype_ref: Vec<Option<&[c_char]>> = ttype.iter().map(|o| o.as_deref()).collect();
            let tform_v = [cc("6000J")];
            let tform_ref: Vec<&[c_char]> = tform_v.iter().map(|v| v.as_slice()).collect();

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);

            fits_create_tbl(
                f.as_deref_mut().unwrap(),
                BINARY_TBL,
                0,
                1,
                &ttype_ref,
                &tform_ref,
                None,
                None,
                &mut status,
            );
            assert_eq!(status, 0, "ffcrtb failed");

            fits_write_tdim(f.as_deref_mut().unwrap(), 1, 3, &naxes, &mut status);
            assert_eq!(status, 0, "ffptdm failed");

            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_putkey_time() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            fits_write_key_triple(
                f.as_deref_mut().unwrap(),
                &cc("TIMEKEY"),
                12345,
                0.6789,
                Some(&cc("time value")),
                &mut status,
            );
            assert_eq!(status, 0, "ffpkyt failed");
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_putkey_rawrec() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);

            let card =
                "RAWKEY  = 'raw value'          / written as raw record                         ";
            fits_write_record(f.as_deref_mut().unwrap(), &cc(card), &mut status);
            assert_eq!(status, 0, "ffprec failed");

            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_putkey_verifydate() {
        let mut status: c_int = 0;

        // Valid dates
        ffverifydate_safe(2024, 1, 1, &mut status);
        assert_eq!(status, 0);
        ffverifydate_safe(2024, 12, 31, &mut status);
        assert_eq!(status, 0);
        ffverifydate_safe(2000, 2, 29, &mut status); // leap year
        assert_eq!(status, 0);

        // Boundary conditions
        ffverifydate_safe(1999, 1, 1, &mut status);
        assert_eq!(status, 0);
        ffverifydate_safe(9999, 12, 31, &mut status);
        assert_eq!(status, 0);
    }

    #[test]
    fn test_putkey_sysdate() {
        let mut status: c_int = 0;
        let mut day: c_int = 0;
        let mut month: c_int = 0;
        let mut year: c_int = 0;

        fits_get_system_date(&mut day, &mut month, &mut year, &mut status);
        assert_eq!(status, 0, "ffgsdt failed");
        assert!(year >= 2020 && year <= 2100);
        assert!(month >= 1 && month <= 12);
        assert!(day >= 1 && day <= 31);
    }

    #[test]
    fn test_putkey_insert() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);

            fits_write_key_lng(
                f.as_deref_mut().unwrap(),
                &cc("KEY1"),
                1,
                Some(&cc("first key")),
                &mut status,
            );
            fits_write_key_lng(
                f.as_deref_mut().unwrap(),
                &cc("KEY3"),
                3,
                Some(&cc("third key")),
                &mut status,
            );

            let card =
                "KEY2    =                    2 / second key                                     ";
            fits_insert_record(f.as_deref_mut().unwrap(), 4, &cc(card), &mut status);
            assert_eq!(status, 0, "ffirec failed");

            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_putkey_header() {
        // ffphpr (own temp file: ffinit cannot overwrite without '!')
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 2] = [100, 100];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_grphdr(
                f.as_deref_mut().unwrap(),
                1,
                16,
                2,
                &naxes,
                0,
                1,
                1,
                &mut status,
            );
            assert_eq!(status, 0, "ffphpr failed");
            fits_close_file(f.take().unwrap(), &mut status);
        });

        // ffphprll
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);
            let naxesll: [LONGLONG; 2] = [200, 200];

            let mut f2: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f2, &name, &mut status);
            fits_write_grphdrll(
                f2.as_deref_mut().unwrap(),
                1,
                16,
                2,
                &naxesll,
                0,
                1,
                1,
                &mut status,
            );
            assert_eq!(status, 0, "ffphprll failed");
            fits_close_file(f2.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_putkey_table_headers() {
        // ffphtb (ASCII table header) -- own temp file
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);

            let ttype = [Some(cc("COL1")), Some(cc("COL2"))];
            let ttype_ref: Vec<Option<&[c_char]>> = ttype.iter().map(|o| o.as_deref()).collect();
            let tform_ascii_v = [cc("A10"), cc("I5")];
            let tform_ascii: Vec<&[c_char]> = tform_ascii_v.iter().map(|v| v.as_slice()).collect();
            let tunit = [Some(cc("unit1")), Some(cc("unit2"))];
            let tunit_ref: Vec<Option<&[c_char]>> = tunit.iter().map(|o| o.as_deref()).collect();
            let tbcol: [c_long; 2] = [1, 11];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            fits_create_hdu(f.as_deref_mut().unwrap(), &mut status);
            fits_write_atblhdr(
                f.as_deref_mut().unwrap(),
                20,
                10,
                2,
                &ttype_ref,
                Some(&tbcol),
                &tform_ascii,
                Some(&tunit_ref),
                Some(&cc("TESTTBL")),
                &mut status,
            );
            assert_eq!(status, 0, "ffphtb failed");
            fits_close_file(f.take().unwrap(), &mut status);
        });

        // ffphbn (binary table header) -- own temp file
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);

            let ttype = [Some(cc("COL1")), Some(cc("COL2"))];
            let ttype_ref: Vec<Option<&[c_char]>> = ttype.iter().map(|o| o.as_deref()).collect();
            let tform_bin_v = [cc("10A"), cc("1J")];
            let tform_bin: Vec<&[c_char]> = tform_bin_v.iter().map(|v| v.as_slice()).collect();
            let tunit = [Some(cc("unit1")), Some(cc("unit2"))];
            let tunit_ref: Vec<Option<&[c_char]>> = tunit.iter().map(|o| o.as_deref()).collect();

            let mut f2: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f2, &name, &mut status);
            fits_write_imghdr(f2.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            fits_create_hdu(f2.as_deref_mut().unwrap(), &mut status);
            fits_write_btblhdr(
                f2.as_deref_mut().unwrap(),
                10,
                2,
                &ttype_ref,
                &tform_bin,
                Some(&tunit_ref),
                Some(&cc("BINTBL")),
                0,
                &mut status,
            );
            assert_eq!(status, 0, "ffphbn failed");
            fits_close_file(f2.take().unwrap(), &mut status);
        });
    }

    /// Mirrors test_putkey_exthead in ~/code/cfitsio/tests/test_putkey.c
    #[test]
    fn test_putkey_exthead() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 2] = [50, 50];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            fits_create_hdu(f.as_deref_mut().unwrap(), &mut status);
            fits_write_exthdr(
                f.as_deref_mut().unwrap(),
                &cc("IMAGE"),
                16,
                2,
                &naxes,
                0,
                1,
                &mut status,
            );
            assert_eq!(status, 0, "ffphext failed");
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0);

            /* the C stops there; check the keywords actually landed */
            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut hdutype: c_int = 0;
            fits_movabs_hdu(
                f.as_deref_mut().unwrap(),
                2,
                Some(&mut hdutype),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(hdutype, IMAGE_HDU);

            let mut bitpix: c_int = 0;
            let mut naxis: c_int = 0;
            let mut naxes_out = [0 as c_long; 2];
            fits_get_img_type(f.as_deref_mut().unwrap(), &mut bitpix, &mut status);
            fits_get_img_dim(f.as_deref_mut().unwrap(), &mut naxis, &mut status);
            fits_get_img_size(f.as_deref_mut().unwrap(), 2, &mut naxes_out, &mut status);
            assert_eq!(status, 0);
            assert_eq!(bitpix, SHORT_IMG);
            assert_eq!(naxis, 2);
            assert_eq!(naxes_out, [50, 50]);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_putkey_updates() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            let g = f.as_deref_mut().unwrap();

            // string
            fits_write_key_str(
                g,
                &cc("STRKEY"),
                &cc("original"),
                Some(&cc("string comment")),
                &mut status,
            );
            fits_update_key_str(
                g,
                &cc("STRKEY"),
                &cc("updated"),
                Some(&cc("new comment")),
                &mut status,
            );

            // float (E format)
            fits_write_key_flt(
                g,
                &cc("FLTKEY"),
                1.5,
                3,
                Some(&cc("float comment")),
                &mut status,
            );
            fits_update_key_flt(
                g,
                &cc("FLTKEY"),
                2.5,
                3,
                Some(&cc("updated float")),
                &mut status,
            );

            // double (D format)
            fits_write_key_dbl(
                g,
                &cc("DBLKEY"),
                1.5,
                6,
                Some(&cc("double comment")),
                &mut status,
            );
            fits_update_key_dbl(
                g,
                &cc("DBLKEY"),
                2.5,
                6,
                Some(&cc("updated double")),
                &mut status,
            );

            // fixed format G
            fits_write_key_fixdbl(g, &cc("GKEY"), 100.5, 4, Some(&cc("G format")), &mut status);
            fits_update_key_fixdbl(
                g,
                &cc("GKEY"),
                200.5,
                4,
                Some(&cc("updated G")),
                &mut status,
            );

            // fixed format F
            fits_write_key_fixflt(g, &cc("FKEY"), 100.5, 3, Some(&cc("F format")), &mut status);
            fits_update_key_fixflt(
                g,
                &cc("FKEY"),
                200.5,
                3,
                Some(&cc("updated F")),
                &mut status,
            );

            assert_eq!(status, 0, "updates failed");
            fits_close_file(f.take().unwrap(), &mut status);

            // Verify updates
            fits_open_file(&mut f, &name, READONLY, &mut status);
            let g = f.as_deref_mut().unwrap();
            let mut sval = [0 as c_char; FLEN_VALUE];
            fits_read_key_str(g, &cc("STRKEY"), &mut sval, None, &mut status);
            assert_eq!(from_buf(&sval), "updated");
            let mut dval: c_double = 0.0;
            fits_read_key(
                g,
                crate::KeywordDatatypeMut::TDOUBLE(&mut dval),
                &cc("DBLKEY"),
                None,
                &mut status,
            );
            assert!(dval >= 2.4 && dval <= 2.6);
            let mut fval: f32 = 0.0;
            fits_read_key(
                g,
                crate::KeywordDatatypeMut::TFLOAT(&mut fval),
                &cc("FLTKEY"),
                None,
                &mut status,
            );
            assert!(fval >= 2.4 && fval <= 2.6);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_putkey_image_ext() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);

            let naxesll: [LONGLONG; 2] = [10, 10];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            fits_create_imgll(
                f.as_deref_mut().unwrap(),
                SHORT_IMG,
                2,
                &naxesll,
                &mut status,
            );
            assert_eq!(status, 0, "ffcrimll failed");
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_putkey_special_images() {
        let img_test = |bitpix: c_int| {
            with_temp_file(|filename| {
                let mut status: c_int = 0;
                let name = to_buf(filename);
                let naxes: [c_long; 2] = [10, 10];
                let mut f: Option<Box<fitsfile>> = None;
                fits_create_file(&mut f, &name, &mut status);
                fits_write_grphdr(
                    f.as_deref_mut().unwrap(),
                    1,
                    bitpix,
                    2,
                    &naxes,
                    0,
                    1,
                    1,
                    &mut status,
                );
                assert_eq!(status, 0, "ffphpr failed for bitpix {bitpix}");
                fits_close_file(f.take().unwrap(), &mut status);
            });
        };
        img_test(USHORT_IMG);
        img_test(ULONG_IMG);
        img_test(SBYTE_IMG);
        img_test(ULONGLONG_IMG);
    }

    #[test]
    fn test_putkey_longstrn() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);

            fits_write_key_longwarn(f.as_deref_mut().unwrap(), &mut status);
            // Calling again should return without writing (keyword exists)
            fits_write_key_longwarn(f.as_deref_mut().unwrap(), &mut status);
            assert_eq!(status, 0, "ffplsw failed");

            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_putkey_random_groups() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);

            let naxes: [c_long; 2] = [10, 10];
            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            // pcount=2, gcount=3
            fits_write_grphdr(
                f.as_deref_mut().unwrap(),
                1,
                16,
                2,
                &naxes,
                2,
                3,
                1,
                &mut status,
            );
            assert_eq!(status, 0, "ffphpr failed");
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_putkey_bad_datatype() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);

            let value: c_int = 42;
            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);

            // Use an invalid datatype code
            status = 0;
            fits_write_key(
                f.as_deref_mut().unwrap(),
                KeywordDatatype::from_datatype(9999, core::ptr::from_ref::<c_int>(&value).cast::<libc::c_void>()),
                &cc("BADTYPE"),
                Some(&cc("bad datatype")),
                &mut status,
            );
            assert_eq!(status, BAD_DATATYPE);

            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_putkey_template_notfound() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);

            status = 0;
            fits_write_key_template(
                f.as_deref_mut().unwrap(),
                &cc("nonexistent_template.tpl"),
                &mut status,
            );
            assert_eq!(status, FILE_NOT_OPENED);

            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_putkey_date_errors() {
        let mut status: c_int;
        let mut yr: c_int = 0;
        let mut mon: c_int = 0;
        let mut day: c_int = 0;
        let mut hr: c_int = 0;
        let mut min: c_int = 0;
        let mut sec: c_double = 0.0;

        // invalid date format for ffs2dt
        status = 0;
        fits_str2date(
            Some(&cc("invalid")),
            Some(&mut yr),
            Some(&mut mon),
            Some(&mut day),
            &mut status,
        );
        assert_eq!(status, BAD_DATE);

        // null input
        status = 0;
        fits_str2date(
            None,
            Some(&mut yr),
            Some(&mut mon),
            Some(&mut day),
            &mut status,
        );
        assert_eq!(status, BAD_DATE);

        // too short string
        status = 0;
        fits_str2date(
            Some(&cc("12")),
            Some(&mut yr),
            Some(&mut mon),
            Some(&mut day),
            &mut status,
        );
        assert_eq!(status, BAD_DATE);

        // invalid old format (wrong separator)
        status = 0;
        fits_str2date(
            Some(&cc("12-34-56")),
            Some(&mut yr),
            Some(&mut mon),
            Some(&mut day),
            &mut status,
        );
        assert_eq!(status, BAD_DATE);

        // ffs2tm with null input
        status = 0;
        fits_str2time(
            None,
            Some(&mut yr),
            Some(&mut mon),
            Some(&mut day),
            Some(&mut hr),
            Some(&mut min),
            Some(&mut sec),
            &mut status,
        );
        assert_eq!(status, BAD_DATE);

        // ffs2tm with invalid format
        status = 0;
        fits_str2time(
            Some(&cc("not-a-date")),
            Some(&mut yr),
            Some(&mut mon),
            Some(&mut day),
            Some(&mut hr),
            Some(&mut min),
            Some(&mut sec),
            &mut status,
        );
        assert_eq!(status, BAD_DATE);

        // ffs2tm with invalid time format
        status = 0;
        fits_str2time(
            Some(&cc("12:34")),
            Some(&mut yr),
            Some(&mut mon),
            Some(&mut day),
            Some(&mut hr),
            Some(&mut min),
            Some(&mut sec),
            &mut status,
        );
        assert_eq!(status, BAD_DATE);

        // ffs2tm with invalid date-time separator
        status = 0;
        fits_str2time(
            Some(&cc("2024-01-01X12:00:00")),
            Some(&mut yr),
            Some(&mut mon),
            Some(&mut day),
            Some(&mut hr),
            Some(&mut min),
            Some(&mut sec),
            &mut status,
        );
        assert_eq!(status, BAD_DATE);
    }

    /// Mirrors test_putkey_time_errors in ~/code/cfitsio/tests/test_putkey.c
    #[test]
    fn test_putkey_time_errors() {
        let mut timestr = [0 as c_char; 40];

        /* each out-of-range field must be rejected with BAD_DATE */
        let cases: [(c_int, c_int, c_int, c_int, c_int, c_double, &str); 6] = [
            (10000, 1, 1, 0, 0, 0.0, "year > 9999"),
            (2024, 13, 1, 0, 0, 0.0, "month > 12"),
            (2024, 1, 32, 0, 0, 0.0, "day > 31"),
            (2024, 1, 1, 24, 0, 0.0, "hour > 23"),
            (2024, 1, 1, 0, 60, 0.0, "minute > 59"),
            (2024, 1, 1, 0, 0, 61.0, "second >= 61"),
        ];

        for (year, month, day, hour, minute, second, what) in cases {
            let mut status: c_int = 0;
            fits_time2str(
                year,
                month,
                day,
                hour,
                minute,
                second,
                0,
                &mut timestr,
                &mut status,
            );
            assert_eq!(status, BAD_DATE, "fftm2s accepted {what}");
        }
    }

    #[test]
    fn test_putkey_verifydate_errors() {
        let mut status: c_int;

        // Invalid year
        status = 0;
        ffverifydate_safe(10000, 1, 1, &mut status);
        assert_eq!(status, BAD_DATE);

        // Invalid month high
        status = 0;
        ffverifydate_safe(2024, 13, 1, &mut status);
        assert_eq!(status, BAD_DATE);

        // Invalid month low
        status = 0;
        ffverifydate_safe(2024, 0, 1, &mut status);
        assert_eq!(status, BAD_DATE);

        // Invalid day high
        status = 0;
        ffverifydate_safe(2024, 1, 32, &mut status);
        assert_eq!(status, BAD_DATE);

        // Invalid day low
        status = 0;
        ffverifydate_safe(2024, 1, 0, &mut status);
        assert_eq!(status, BAD_DATE);

        // Invalid Feb 29 in non-leap year
        status = 0;
        ffverifydate_safe(2023, 2, 29, &mut status);
        assert_eq!(status, BAD_DATE);

        // Valid Feb 29 in leap year
        status = 0;
        ffverifydate_safe(2024, 2, 29, &mut status);
        assert_eq!(status, 0);
    }

    #[test]
    fn test_putkey_dt2s_errors() {
        let mut status: c_int = 0;
        let mut datestr = [0 as c_char; 11];

        // Invalid date
        status = 0;
        fits_date2str(2024, 13, 1, &mut datestr, &mut status);
        assert_eq!(status, BAD_DATE);
    }

    #[test]
    fn test_putkey_time_keyword_errors() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);

            // fraction out of range (negative)
            status = 0;
            fits_write_key_triple(
                f.as_deref_mut().unwrap(),
                &cc("BADTIME1"),
                12345,
                -0.5,
                Some(&cc("bad fraction")),
                &mut status,
            );
            assert_eq!(status, BAD_F2C);

            // fraction out of range (> 1)
            status = 0;
            fits_write_key_triple(
                f.as_deref_mut().unwrap(),
                &cc("BADTIME2"),
                12345,
                1.5,
                Some(&cc("bad fraction")),
                &mut status,
            );
            assert_eq!(status, BAD_F2C);

            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_putkey_tdim_errors() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);

            let naxes: [c_long; 3] = [10, 20, 30];
            let ttype = [Some(cc("DATA"))];
            let ttype_ref: Vec<Option<&[c_char]>> = ttype.iter().map(|o| o.as_deref()).collect();
            let tform_v = [cc("6000J")];
            let tform_ref: Vec<&[c_char]> = tform_v.iter().map(|v| v.as_slice()).collect();

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_create_tbl(
                f.as_deref_mut().unwrap(),
                BINARY_TBL,
                0,
                1,
                &ttype_ref,
                &tform_ref,
                None,
                None,
                &mut status,
            );
            assert_eq!(status, 0, "ffcrtb failed");

            // column number out of range (0)
            status = 0;
            fits_write_tdim(f.as_deref_mut().unwrap(), 0, 3, &naxes, &mut status);
            assert_eq!(status, BAD_COL_NUM);

            // column number out of range (1000)
            status = 0;
            fits_write_tdim(f.as_deref_mut().unwrap(), 1000, 3, &naxes, &mut status);
            assert_eq!(status, BAD_COL_NUM);

            // naxis < 1
            status = 0;
            fits_write_tdim(f.as_deref_mut().unwrap(), 1, 0, &naxes, &mut status);
            assert_eq!(status, BAD_DIMEN);

            // negative axis value
            let bad_naxes: [c_long; 3] = [10, -20, 30];
            status = 0;
            fits_write_tdim(f.as_deref_mut().unwrap(), 1, 3, &bad_naxes, &mut status);
            assert_eq!(status, BAD_TDIM);

            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_putkey_tdim_nottable() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);

            let naxes: [c_long; 2] = [10, 10];
            let tdim: [c_long; 2] = [5, 20];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_grphdr(
                f.as_deref_mut().unwrap(),
                1,
                16,
                2,
                &naxes,
                0,
                1,
                1,
                &mut status,
            );

            // Try to write TDIM on image HDU - should fail
            status = 0;
            fits_write_tdim(f.as_deref_mut().unwrap(), 1, 2, &tdim, &mut status);
            assert_eq!(status, NOT_BTABLE);

            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_putkey_very_long_string() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);

            let longstr: String = "x".repeat(400);

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            fits_write_key_longstr(
                f.as_deref_mut().unwrap(),
                &cc("VERYLONG"),
                &cc(&longstr),
                Some(&cc("very long string")),
                &mut status,
            );
            assert_eq!(status, 0, "ffpkls failed");
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut result: *mut c_char = core::ptr::null_mut();
            fits_read_key_longstr(
                f.as_deref_mut().unwrap(),
                &cc("VERYLONG"),
                &mut result,
                None,
                &mut status,
            );
            assert_eq!(status, 0, "ffgkls failed");
            assert!(!result.is_null());
            let got = unsafe { core::ffi::CStr::from_ptr(result) }
                .to_str()
                .unwrap()
                .to_string();
            assert_eq!(got, longstr);
            unsafe { fits_free_memory(result.cast::<libc::c_void>(), &mut status) };
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_putkey_long_string_quotes() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);

            let longstr = "This is a string with 'quotes' and more 'quotes' \
                that needs to be properly escaped when written";

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            fits_write_key_longstr(
                f.as_deref_mut().unwrap(),
                &cc("QUOTSTR"),
                &cc(longstr),
                Some(&cc("string with quotes")),
                &mut status,
            );
            assert_eq!(status, 0, "ffpkls failed");
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut result: *mut c_char = core::ptr::null_mut();
            fits_read_key_longstr(
                f.as_deref_mut().unwrap(),
                &cc("QUOTSTR"),
                &mut result,
                None,
                &mut status,
            );
            assert_eq!(status, 0, "ffgkls failed");
            assert!(!result.is_null());
            let got = unsafe { core::ffi::CStr::from_ptr(result) }
                .to_str()
                .unwrap()
                .to_string();
            assert_eq!(got, longstr);
            unsafe { fits_free_memory(result.cast::<libc::c_void>(), &mut status) };
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_putkey_header_longlong() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);

            let naxes: [c_long; 2] = [20, 20];
            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_grphdr(
                f.as_deref_mut().unwrap(),
                1,
                LONGLONG_IMG,
                2,
                &naxes,
                0,
                1,
                1,
                &mut status,
            );
            assert_eq!(status, 0, "ffphpr failed");
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_putkey_bintable_vla() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);

            let ttype = [Some(cc("COL1")), Some(cc("COL2"))];
            let ttype_ref: Vec<Option<&[c_char]>> = ttype.iter().map(|o| o.as_deref()).collect();
            let tform_v = [cc("1J"), cc("PJ")]; // variable length array
            let tform_ref: Vec<&[c_char]> = tform_v.iter().map(|v| v.as_slice()).collect();

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            fits_create_hdu(f.as_deref_mut().unwrap(), &mut status);
            fits_write_btblhdr(
                f.as_deref_mut().unwrap(),
                10,
                2,
                &ttype_ref,
                &tform_ref,
                None,
                Some(&cc("TESTTBL")),
                100,
                &mut status,
            );
            assert_eq!(status, 0, "ffphbn failed");
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_putkey_asciitable_multi() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);

            let ttype = [Some(cc("NAME")), Some(cc("VALUE")), Some(cc("COUNT"))];
            let ttype_ref: Vec<Option<&[c_char]>> = ttype.iter().map(|o| o.as_deref()).collect();
            let tform_v = [cc("A20"), cc("E15.7"), cc("I10")];
            let tform_ref: Vec<&[c_char]> = tform_v.iter().map(|v| v.as_slice()).collect();
            let tunit = [Some(cc("text")), Some(cc("meters")), Some(cc("items"))];
            let tunit_ref: Vec<Option<&[c_char]>> = tunit.iter().map(|o| o.as_deref()).collect();
            let tbcol: [c_long; 3] = [1, 21, 37];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            fits_create_hdu(f.as_deref_mut().unwrap(), &mut status);
            fits_write_atblhdr(
                f.as_deref_mut().unwrap(),
                50,
                100,
                3,
                &ttype_ref,
                Some(&tbcol),
                &tform_ref,
                Some(&tunit_ref),
                Some(&cc("TESTTBL")),
                &mut status,
            );
            assert_eq!(status, 0, "ffphtb failed");
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    /// Mirrors test_putkey_extheader_3d in ~/code/cfitsio/tests/test_putkey.c
    #[test]
    fn test_putkey_extheader_3d() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 3] = [10, 20, 5];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            fits_create_hdu(f.as_deref_mut().unwrap(), &mut status);
            fits_write_exthdr(
                f.as_deref_mut().unwrap(),
                &cc("IMAGE"),
                32,
                3,
                &naxes,
                0,
                1,
                &mut status,
            );
            assert_eq!(status, 0, "ffphext failed");
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut bitpix: c_int = 0;
            let mut naxis: c_int = 0;
            let mut naxes_out = [0 as c_long; 3];
            fits_get_img_type(f.as_deref_mut().unwrap(), &mut bitpix, &mut status);
            fits_get_img_dim(f.as_deref_mut().unwrap(), &mut naxis, &mut status);
            fits_get_img_size(f.as_deref_mut().unwrap(), 3, &mut naxes_out, &mut status);
            assert_eq!(status, 0);
            assert_eq!(bitpix, LONG_IMG);
            assert_eq!(naxis, 3);
            assert_eq!(naxes_out, [10, 20, 5]);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_putkey_empty_string() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);

            // Test writing empty string
            fits_write_key_str(
                f.as_deref_mut().unwrap(),
                &cc("EMPTY"),
                &cc(""),
                Some(&cc("empty string test")),
                &mut status,
            );
            assert_eq!(status, 0, "ffpkys failed");
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut sval = [0 as c_char; FLEN_VALUE];
            fits_read_key_str(
                f.as_deref_mut().unwrap(),
                &cc("EMPTY"),
                &mut sval,
                None,
                &mut status,
            );
            assert_eq!(status, 0, "ffgkys failed");
            assert_eq!(from_buf(&sval), "");
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    /// Mirrors test_putkey_timeonly in ~/code/cfitsio/tests/test_putkey.c
    #[test]
    fn test_putkey_timeonly() {
        let mut status: c_int = 0;
        let mut timestr = [0 as c_char; 40];

        /* year = month = day = 0 selects the time-only format hh:mm:ss.ddd */
        fits_time2str(0, 0, 0, 14, 30, 45.5, 3, &mut timestr, &mut status);
        assert_eq!(status, 0);

        let got = from_buf(&timestr);
        assert!(got.len() >= 8, "too short: {got:?}");
        assert_eq!(got.as_bytes()[2], b':');
        assert_eq!(got.as_bytes()[5], b':');
        assert_eq!(got, "14:30:45.500");
    }

    #[test]
    fn test_putkey_date_edge_cases() {
        let mut status: c_int = 0;
        let mut yr: c_int = 0;
        let mut mon: c_int = 0;
        let mut day: c_int = 0;
        let mut datestr = [0 as c_char; 11];

        // century boundary (year 2000)
        fits_date2str(2000, 1, 1, &mut datestr, &mut status);
        assert_eq!(status, 0);
        assert_eq!(from_buf(&datestr), "2000-01-01");

        // old format with single digit day
        fits_str2date(
            Some(&cc("01/02/98")),
            Some(&mut yr),
            Some(&mut mon),
            Some(&mut day),
            &mut status,
        );
        assert_eq!(status, 0);
        assert_eq!(yr, 1998);
        assert_eq!(mon, 2);
        assert_eq!(day, 1);
    }

    #[test]
    fn test_putkey_complex_special() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);

            let cval: [f32; 2] = [1.5, 2.5];
            let mval: [f64; 2] = [3.5, 4.5];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            let g = f.as_deref_mut().unwrap();

            // very small precision
            fits_write_key_cmp(
                g,
                &cc("CKEY1"),
                &cval,
                0,
                Some(&cc("complex low precision")),
                &mut status,
            );
            fits_write_key_dblcmp(
                g,
                &cc("MKEY1"),
                &mval,
                0,
                Some(&cc("dblcomplex low precision")),
                &mut status,
            );

            // high precision
            fits_write_key_cmp(
                g,
                &cc("CKEY2"),
                &cval,
                14,
                Some(&cc("complex high precision")),
                &mut status,
            );
            fits_write_key_dblcmp(
                g,
                &cc("MKEY2"),
                &mval,
                14,
                Some(&cc("dblcomplex high precision")),
                &mut status,
            );

            // fixed format complex
            fits_write_key_fixcmp(
                g,
                &cc("FCKEY"),
                &cval,
                0,
                Some(&cc("fixed complex low prec")),
                &mut status,
            );
            fits_write_key_fixdblcmp(
                g,
                &cc("FMKEY"),
                &mval,
                0,
                Some(&cc("fixed dblcomplex low prec")),
                &mut status,
            );

            assert_eq!(status, 0, "complex special failed");
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_putkey_old_date() {
        let mut status: c_int = 0;
        let mut datestr = [0 as c_char; 11];

        // pre-1999 date (should use DD/MM/YY format)
        fits_date2str(1985, 3, 15, &mut datestr, &mut status);
        assert_eq!(status, 0);
        assert_eq!(from_buf(&datestr), "15/03/85");

        // 1998 (last year for old format)
        fits_date2str(1998, 12, 31, &mut datestr, &mut status);
        assert_eq!(status, 0);
        assert_eq!(from_buf(&datestr), "31/12/98");

        // 1999 (first year for new format)
        fits_date2str(1999, 1, 1, &mut datestr, &mut status);
        assert_eq!(status, 0);
        assert_eq!(from_buf(&datestr), "1999-01-01");
    }

    #[test]
    fn test_putkey_update_unsigned() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);

            let mut ulval: c_ulong = 3000000000;

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            let g = f.as_deref_mut().unwrap();

            // Write and update unsigned long
            fits_write_key(
                g,
                KeywordDatatype::TULONG(&ulval),
                &cc("ULKEY"),
                Some(&cc("unsigned long")),
                &mut status,
            );
            ulval = 4000000000;
            fits_update_key(
                g,
                KeywordDatatype::TULONG(&ulval),
                &cc("ULKEY"),
                Some(&cc("updated unsigned long")),
                &mut status,
            );
            assert_eq!(status, 0, "ffuky failed");

            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_putkey_array_null_comment() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);

            let svals = [cc("alpha"), cc("beta"), cc("gamma")];
            let svals_ptrs: Vec<*const c_char> = svals.iter().map(|v| v.as_ptr()).collect();
            let lvals: [c_int; 3] = [1, 0, 1];
            let jvals: [c_long; 3] = [100, 200, 300];
            let fvals: [f32; 3] = [1.1, 2.2, 3.3];
            let dvals: [f64; 3] = [4.4, 5.5, 6.6];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            let g = f.as_deref_mut().unwrap();

            // Test array keywords with NULL comments (empty pointer slice)
            let no_comm: [*const c_char; 0] = [];
            fits_write_keys_str(g, &cc("SKEY"), 1, 3, &svals_ptrs, &no_comm, &mut status);
            fits_write_keys_log(g, &cc("LKEY"), 1, 3, &lvals, &no_comm, &mut status);
            fits_write_keys_lng(g, &cc("JKEY"), 1, 3, &jvals, &no_comm, &mut status);
            fits_write_keys_fixflt(g, &cc("FKEY"), 1, 3, &fvals, 2, &no_comm, &mut status);
            fits_write_keys_dbl(g, &cc("DKEY"), 1, 3, &dvals, 4, &no_comm, &mut status);

            assert_eq!(status, 0, "array null comment failed");
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_putkey_array_nonrepeat_comment() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);

            let svals = [cc("one"), cc("two"), cc("three")];
            let svals_ptrs: Vec<*const c_char> = svals.iter().map(|v| v.as_ptr()).collect();
            let lvals: [c_int; 3] = [0, 1, 0];
            let jvals: [c_long; 3] = [1000, 2000, 3000];
            let evals: [f32; 3] = [1.0e10, 2.0e10, 3.0e10];
            let gvals: [f64; 3] = [1.111, 2.222, 3.333];
            // Non-repeating comments (no trailing &)
            let comm = [cc("comment one"), cc("comment two"), cc("comment three")];
            let comm_ptrs: Vec<*const c_char> = comm.iter().map(|v| v.as_ptr()).collect();

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            let g = f.as_deref_mut().unwrap();

            fits_write_keys_str(g, &cc("NRSKEY"), 1, 3, &svals_ptrs, &comm_ptrs, &mut status);
            fits_write_keys_log(g, &cc("NRLKEY"), 1, 3, &lvals, &comm_ptrs, &mut status);
            fits_write_keys_lng(g, &cc("NRJKEY"), 1, 3, &jvals, &comm_ptrs, &mut status);
            fits_write_keys_flt(g, &cc("NREKEY"), 1, 3, &evals, 3, &comm_ptrs, &mut status);
            fits_write_keys_fixdbl(g, &cc("NRGKEY"), 1, 3, &gvals, 4, &comm_ptrs, &mut status);

            assert_eq!(status, 0, "array nonrepeat comment failed");
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_putkey_longstr_null_comment() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);

            let longstr: String = "a".repeat(100);

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            fits_write_key_longstr(
                f.as_deref_mut().unwrap(),
                &cc("LSNULL"),
                &cc(&longstr),
                None,
                &mut status,
            );
            assert_eq!(status, 0, "ffpkls failed");
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut result: *mut c_char = core::ptr::null_mut();
            fits_read_key_longstr(
                f.as_deref_mut().unwrap(),
                &cc("LSNULL"),
                &mut result,
                None,
                &mut status,
            );
            assert_eq!(status, 0, "ffgkls failed");
            assert!(!result.is_null());
            let got = unsafe { core::ffi::CStr::from_ptr(result) }
                .to_str()
                .unwrap()
                .to_string();
            assert_eq!(got, longstr);
            unsafe { fits_free_memory(result.cast::<libc::c_void>(), &mut status) };
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_putkey_hierarch() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);

            // Explicit HIERARCH keyword
            fits_write_key_str(
                f.as_deref_mut().unwrap(),
                &cc("HIERARCH MY.LONG.KEYWORD"),
                &cc("hier value"),
                Some(&cc("hierarch test")),
                &mut status,
            );

            // Long keyword that will have HIERARCH prepended
            fits_write_key_str(
                f.as_deref_mut().unwrap(),
                &cc("MY.EXTRA.LONG.KEYWORD.NAME"),
                &cc("auto value"),
                Some(&cc("auto hierarch")),
                &mut status,
            );

            assert_eq!(status, 0, "hierarch failed");
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_putkey_leading_spaces() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);

            // Keyword with trailing space - should be trimmed
            fits_write_key_str(
                f.as_deref_mut().unwrap(),
                &cc("SPACEKEY "),
                &cc("test value"),
                Some(&cc("space test")),
                &mut status,
            );
            assert_eq!(status, 0, "ffpkys failed");
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut sval = [0 as c_char; FLEN_VALUE];
            fits_read_key_str(
                f.as_deref_mut().unwrap(),
                &cc("SPACEKEY"),
                &mut sval,
                None,
                &mut status,
            );
            assert_eq!(status, 0, "ffgkys failed");
            assert_eq!(from_buf(&sval), "test value");
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_putkey_float_formats() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            let g = f.as_deref_mut().unwrap();

            // zero values
            let mut fval: f32 = 0.0;
            let mut dval: f64 = 0.0;
            fits_write_key_flt(
                g,
                &cc("FZERO"),
                fval,
                6,
                Some(&cc("zero float E")),
                &mut status,
            );
            fits_write_key_fixflt(
                g,
                &cc("FZEROF"),
                fval,
                3,
                Some(&cc("zero float F")),
                &mut status,
            );
            fits_write_key_dbl(
                g,
                &cc("DZERO"),
                dval,
                10,
                Some(&cc("zero double E")),
                &mut status,
            );
            fits_write_key_fixdbl(
                g,
                &cc("DZEROG"),
                dval,
                6,
                Some(&cc("zero double G")),
                &mut status,
            );

            // very small values
            fval = 1.23e-30;
            dval = 4.56e-200;
            fits_write_key_flt(
                g,
                &cc("FSMALL"),
                fval,
                6,
                Some(&cc("small float E")),
                &mut status,
            );
            fits_write_key_dbl(
                g,
                &cc("DSMALL"),
                dval,
                10,
                Some(&cc("small double E")),
                &mut status,
            );

            // very large values
            fval = 9.87e30;
            dval = 1.23e200;
            fits_write_key_flt(
                g,
                &cc("FLARGE"),
                fval,
                6,
                Some(&cc("large float E")),
                &mut status,
            );
            fits_write_key_dbl(
                g,
                &cc("DLARGE"),
                dval,
                10,
                Some(&cc("large double E")),
                &mut status,
            );

            assert_eq!(status, 0, "float formats failed");
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_putkey_array_trailing_blanks() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);

            let jvals: [c_long; 3] = [10, 20, 30];
            // Comments with trailing blanks followed by &
            let comm = [cc("comment with spaces   &"), cc("second"), cc("third")];
            let comm_ptrs: Vec<*const c_char> = comm.iter().map(|v| v.as_ptr()).collect();

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);

            fits_write_keys_lng(
                f.as_deref_mut().unwrap(),
                &cc("TBKEY"),
                1,
                3,
                &jvals,
                &comm_ptrs,
                &mut status,
            );
            assert_eq!(status, 0, "ffpknj failed");
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_putkey_complex_too_long() {
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let name = to_buf(filename);

            let cval: [f32; 2] = [1.5, 2.5];
            let mval: [f64; 2] = [3.5, 4.5];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);

            // Use very high precision to trigger BAD_F2C error
            status = 0;
            fits_write_key_cmp(
                f.as_deref_mut().unwrap(),
                &cc("CKEY"),
                &cval,
                40,
                Some(&cc("too long")),
                &mut status,
            );
            assert_eq!(status, BAD_F2C);

            status = 0;
            fits_write_key_dblcmp(
                f.as_deref_mut().unwrap(),
                &cc("MKEY"),
                &mval,
                40,
                Some(&cc("too long")),
                &mut status,
            );
            assert_eq!(status, BAD_F2C);

            status = 0;
            fits_write_key_fixcmp(
                f.as_deref_mut().unwrap(),
                &cc("FCKEY"),
                &cval,
                40,
                Some(&cc("too long")),
                &mut status,
            );
            assert_eq!(status, BAD_F2C);

            status = 0;
            fits_write_key_fixdblcmp(
                f.as_deref_mut().unwrap(),
                &cc("FMKEY"),
                &mval,
                40,
                Some(&cc("too long")),
                &mut status,
            );
            assert_eq!(status, BAD_F2C);

            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }
}
