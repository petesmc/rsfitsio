/*  This file, putkey.c, contains routines that write keywords to          */
/*  a FITS header.                                                         */

/*  The FITSIO software was written by William Pence at the High Energy    */
/*  Astrophysic Science Archive Research Center (HEASARC) at the NASA      */
/*  Goddard Space Flight Center.                                           */

use core::{slice, str};
use std::ffi::{CStr, CString};
use std::fs::File;
use std::io::BufRead;
use std::num::ParseIntError;
use std::{cmp, mem, ptr};

use chrono::{DateTime, Utc};

use crate::c_types::*;

use bytemuck::{cast_slice, cast_slice_mut};

use crate::fitscore::{
    ffbnfm_safe, ffcrhd_safer, ffgabc_safe, ffgthd_safe, ffiblk, ffkeyn_safe, ffmahd_safe,
    ffmkky_safe, ffpmsg_slice, ffpmsg_str, ffrdef_safe, fftkey_safe, ffupch_safe, fits_strncasecmp,
};
use crate::getkey::ffgkys_safe;
use crate::imcompress::imcomp_init_table;
use crate::modkey::{ffdkey_safe, ffirec_safe, ffmnam_safe, ffucrd_safe};
use crate::relibc::header::stdio::snprintf_f64_decim;
use crate::{KeywordDatatype, fitsio2::*};
use crate::{TKeywords, wrappers::*};
use crate::{atoi, int_snprintf};
use crate::{bb, cs};
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
        let naxes = slice::from_raw_parts(naxes, naxis as usize);
        ffcrim_safer(fptr, bitpix, naxis, naxes, status)
    }
}

/*--------------------------------------------------------------------------*/
/// create an IMAGE extension following the current HDU. If the
/// current HDU is empty (contains no header keywords), then simply
/// write the required image (or primary array) keywords to the current
/// HDU.
pub unsafe fn ffcrim_safer(
    fptr: &mut fitsfile, /* I - FITS file pointer           */
    bitpix: c_int,       /* I - bits per pixel              */
    naxis: c_int,        /* I - number of axes in the array */
    naxes: &[c_long],    /* I - size of each axis           */
    status: &mut c_int,  /* IO - error status               */
) -> c_int {
    unsafe {
        if *status > 0 {
            return *status;
        }

        if fptr.HDUposition != fptr.Fptr.curhdu {
            ffmahd_safe(fptr, (fptr.HDUposition) + 1, None, status);
        }

        /* create new extension if current header is not empty */
        let h = fptr.Fptr.get_headstart_as_slice();
        if fptr.Fptr.headend != h[fptr.Fptr.curhdu as usize] {
            ffcrhd_safer(fptr, status);
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

        ffcrimll_safer(fptr, bitpix, naxis, naxes, status)
    }
}

/*--------------------------------------------------------------------------*/
/// create an IMAGE extension following the current HDU. If the
/// current HDU is empty (contains no header keywords), then simply
/// write the required image (or primary array) keywords to the current
/// HDU.
pub unsafe fn ffcrimll_safer(
    fptr: &mut fitsfile, /* I - FITS file pointer           */
    bitpix: c_int,       /* I - bits per pixel              */
    naxis: c_int,        /* I - number of axes in the array */
    naxes: &[LONGLONG],  /* I - size of each axis           */
    status: &mut c_int,  /* IO - error status               */
) -> c_int {
    unsafe {
        if *status > 0 {
            return *status;
        }

        if fptr.HDUposition != fptr.Fptr.curhdu {
            ffmahd_safe(fptr, (fptr.HDUposition) + 1, None, status);
        }

        /* create new extension if current header is not empty */
        let headstart = fptr.Fptr.get_headstart_as_slice();
        if fptr.Fptr.headend != headstart[fptr.Fptr.curhdu as usize] {
            ffcrhd_safer(fptr, status);
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

        ffcrtb_safer(
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
pub unsafe fn ffcrtb_safer(
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
    unsafe {
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
            ffcrhd_safer(fptr, status);
        }

        if fptr.Fptr.curhdu == 0 {
            /* have to create dummy primary array */
            ffcrim_safer(fptr, 16, 0, &[], status);
            ffcrhd_safer(fptr, status);
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
    todo!()
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
        longbitpix as LONGLONG,
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
    ffpkyj_safe(fptr, cs!(c"NAXIS"), naxis as LONGLONG, Some(&comm), status);

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
    let ii: c_int = 0;
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
    } else {
        tbcol_slice = tbcol.unwrap();
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

        if let Some(tunit) = tunit {
            if let Some(tunit_item) = tunit[ii]
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
    let ii: c_int = 0;
    let mut datatype: c_int = 0;
    let mut iread: c_int = 0;
    let mut repeat: c_long = 0;
    let mut width: c_long = 0;
    let mut naxis1: LONGLONG = 0;
    let mut tfmt: [c_char; 30] = [0; 30];
    let mut name: [c_char; FLEN_KEYWORD] = [0; FLEN_KEYWORD];
    let mut comm: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];
    let mut extnm: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let cptr: *mut c_char = ptr::null_mut();
    let mut card: [c_char; FLEN_CARD] = [0; FLEN_CARD];
    let colptr: *mut tcolumn = ptr::null_mut();

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
            naxis1 += repeat as LONGLONG * (datatype as LONGLONG / 10);
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
        tfields as LONGLONG,
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

            if cptr.is_some() {
                let c = cptr.unwrap() + 1;

                // iread = sscanf(tfmt[c..].as_ptr(), c"%ld".as_ptr(), &width);
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

        if let Some(tunit) = tunit {
            if let Some(tunit_item) = tunit[ii]
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
    todo!();
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
        let tcard: [c_char; FLEN_CARD] = [0; FLEN_CARD];
        let len: usize = 0;
        let ii: usize = 0;
        let nblocks: c_long = 0;
        let keylength: c_int = 0;

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
    let mut len: usize = 0;
    let mut ii: usize = 0;
    let mut nblocks: c_long = 0;
    let mut keylength: c_int = 0;

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

    len = strlen_safe(&tcard);

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

        fits_make_longstr_key_util(fptr, keyname, value, comm, -1, status);

        *status
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

    let commlen = 0;
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
        let mut keyname: [c_char; FLEN_KEYWORD] = [0; FLEN_KEYWORD];
        let mut tcomment: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];

        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        raw_to_slice!(keyroot);

        let comm = if comm.is_null() {
            &[]
        } else {
            slice::from_raw_parts(comm, nkey as usize)
        };

        let value = slice::from_raw_parts(value, nkey as usize);

        if *status > 0 {
            /* inherit input status value if > 0 */
            return *status;
        }

        /* check if first comment string is to be repeated for all the keywords */
        /* by looking to see if the last non-blank character is a '&' char      */

        let mut repeat = false;

        if !comm.is_empty() {
            let comm_item = cast_slice(CStr::from_ptr(comm[0]).to_bytes_with_nul());

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

            let value_item = cast_slice(CStr::from_ptr(value[ii]).to_bytes_with_nul());

            if repeat {
                ffpkys_safe(fptr, &keyname, value_item, Some(&tcomment), status);
            } else {
                let c: Option<&[c_char]> = match comm[ii].is_null() {
                    true => None,
                    false => Some(cast_slice(CStr::from_ptr(comm[ii]).to_bytes_with_nul())),
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
        let mut keyname: [c_char; FLEN_KEYWORD] = [0; FLEN_KEYWORD];
        let mut tcomment: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];

        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        raw_to_slice!(keyroot);

        let comm = if comm.is_null() {
            &[]
        } else {
            slice::from_raw_parts(comm, nkey as usize)
        };

        let value = slice::from_raw_parts(value, nkey as usize);

        if *status > 0 {
            /* inherit input status value if > 0 */
            return *status;
        }

        /* check if first comment string is to be repeated for all the keywords */
        /* by looking to see if the last non-blank character is a '&' char      */

        let mut repeat = false;

        if !comm.is_empty() {
            let comm_item = cast_slice(CStr::from_ptr(comm[0]).to_bytes_with_nul());

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
                    Some(cast_slice(CStr::from_ptr(comm[ii]).to_bytes_with_nul()))
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
        let mut keyname: [c_char; FLEN_KEYWORD] = [0; FLEN_KEYWORD];
        let mut tcomment: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];

        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        raw_to_slice!(keyroot);

        let comm = if comm.is_null() {
            &[]
        } else {
            slice::from_raw_parts(comm, nkey as usize)
        };

        let value = slice::from_raw_parts(value, nkey as usize);

        if *status > 0 {
            /* inherit input status value if > 0 */
            return *status;
        }

        /* check if first comment string is to be repeated for all the keywords */
        /* by looking to see if the last non-blank character is a '&' char      */

        let mut repeat = false;

        if !comm.is_empty() {
            let comm_item = cast_slice(CStr::from_ptr(comm[0]).to_bytes_with_nul());

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
                    Some(cast_slice(CStr::from_ptr(comm[ii]).to_bytes_with_nul()))
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
        let mut keyname: [c_char; FLEN_KEYWORD] = [0; FLEN_KEYWORD];
        let mut tcomment: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];

        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        raw_to_slice!(keyroot);

        let comm = if comm.is_null() {
            &[]
        } else {
            slice::from_raw_parts(comm, nkey as usize)
        };

        let value = slice::from_raw_parts(value, nkey as usize);

        if *status > 0 {
            /* inherit input status value if > 0 */
            return *status;
        }

        /* check if first comment string is to be repeated for all the keywords */
        /* by looking to see if the last non-blank character is a '&' char      */

        let mut repeat = false;

        if !comm.is_empty() {
            let comm_item = cast_slice(CStr::from_ptr(comm[0]).to_bytes_with_nul());

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
                    Some(cast_slice(CStr::from_ptr(comm[ii]).to_bytes_with_nul()))
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
        let mut keyname: [c_char; FLEN_KEYWORD] = [0; FLEN_KEYWORD];
        let mut tcomment: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];

        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        raw_to_slice!(keyroot);

        let comm = if comm.is_null() {
            &[]
        } else {
            slice::from_raw_parts(comm, nkey as usize)
        };

        let value = slice::from_raw_parts(value, nkey as usize);

        if *status > 0 {
            /* inherit input status value if > 0 */
            return *status;
        }

        /* check if first comment string is to be repeated for all the keywords */
        /* by looking to see if the last non-blank character is a '&' char      */

        let mut repeat = false;

        if !comm.is_empty() {
            let comm_item = cast_slice(CStr::from_ptr(comm[0]).to_bytes_with_nul());

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
                    Some(cast_slice(CStr::from_ptr(comm[ii]).to_bytes_with_nul()))
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
        let mut keyname: [c_char; FLEN_KEYWORD] = [0; FLEN_KEYWORD];
        let mut tcomment: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];

        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        raw_to_slice!(keyroot);

        let comm = if comm.is_null() {
            &[]
        } else {
            slice::from_raw_parts(comm, nkey as usize)
        };

        let value = slice::from_raw_parts(value, nkey as usize);

        if *status > 0 {
            /* inherit input status value if > 0 */
            return *status;
        }

        /* check if first comment string is to be repeated for all the keywords */
        /* by looking to see if the last non-blank character is a '&' char      */

        let mut repeat = false;

        if !comm.is_empty() {
            let comm_item = cast_slice(CStr::from_ptr(comm[0]).to_bytes_with_nul());

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
                    Some(cast_slice(CStr::from_ptr(comm[ii]).to_bytes_with_nul()))
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
        let mut keyname: [c_char; FLEN_KEYWORD] = [0; FLEN_KEYWORD];
        let mut tcomment: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];

        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        raw_to_slice!(keyroot);

        let comm = if comm.is_null() {
            &[]
        } else {
            slice::from_raw_parts(comm, nkey as usize)
        };

        let value = slice::from_raw_parts(value, nkey as usize);

        if *status > 0 {
            /* inherit input status value if > 0 */
            return *status;
        }

        /* check if first comment string is to be repeated for all the keywords */
        /* by looking to see if the last non-blank character is a '&' char      */

        let mut repeat = false;

        if !comm.is_empty() {
            let comm_item = cast_slice(CStr::from_ptr(comm[0]).to_bytes_with_nul());

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
                    Some(cast_slice(CStr::from_ptr(comm[ii]).to_bytes_with_nul()))
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

        let datatype_with_data = KeywordDatatype::from_datatype(datatype, value as *mut _);

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
            ffpkyj_safe(fptr, keyname, (*(value)) as LONGLONG, comm, status);
        }
        KeywordDatatype::TSBYTE(value) => {
            ffpkyj_safe(fptr, keyname, (*(value)) as LONGLONG, comm, status);
        }
        KeywordDatatype::TUSHORT(value) => {
            ffpkyj_safe(fptr, keyname, (*(value)) as LONGLONG, comm, status);
        }
        KeywordDatatype::TSHORT(value) => {
            ffpkyj_safe(fptr, keyname, (*(value)) as LONGLONG, comm, status);
        }
        KeywordDatatype::TUINT(value) => {
            ffpkyg_safe(fptr, keyname, (*(value)) as f64, 0, comm, status);
        }
        KeywordDatatype::TINT(value) => {
            ffpkyj_safe(fptr, keyname, (*(value)) as LONGLONG, comm, status);
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
    let mut slen: usize = 0;

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
        slen = strlen_safe(&template); /* get string length */

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
        let mut keyname: [c_char; FLEN_KEYWORD] = [0; FLEN_KEYWORD];
        let mut tdimstr: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
        let mut comm: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];
        let mut value: [c_char; 80] = [0; 80];
        let mut message: [c_char; FLEN_ERRMSG] = [0; FLEN_ERRMSG];
        let ii: c_int = 0;
        let mut totalpix: c_long = 1;
        let mut repeat: c_long = 0;

        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        let naxes = slice::from_raw_parts(naxes, naxis as usize);

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
        } else if fptr.Fptr.datastart == DATA_UNDEFINED as LONGLONG && ffrdef_safe(fptr, status) > 0
        {
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
        let c = slice::from_raw_parts_mut(colptr, fptr.Fptr.tfield as usize);
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
    todo!();
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
        let datestr = datestr.as_mut().expect(NULL_MSG);

        let datestr: &mut [c_char; 11] = slice::from_raw_parts_mut(datestr, 11).try_into().unwrap();

        ffdt2s_safe(year, month, day, datestr, status)
    }
}

/*-----------------------------------------------------------------*/
/// Construct a date character string
pub fn ffdt2s_safe(
    year: c_int,            /* I - year (0 - 9999)           */
    month: c_int,           /* I - month (1 - 12)            */
    day: c_int,             /* I - day (1 - 31)              */
    datestr: &[c_char; 11], /* O - date string: "YYYY-MM-DD" */
    status: *mut c_int,     /* IO - error status             */
) -> c_int {
    todo!()
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
        let year = year.as_mut().expect(NULL_MSG);
        let month = month.as_mut().expect(NULL_MSG);
        let day = day.as_mut().expect(NULL_MSG);
        raw_to_slice!(datestr);

        ffs2dt_safe(cast_slice(datestr), year, month, day, status)
    }
}

/*-----------------------------------------------------------------*/
/// Parse a date character string into year, month, and day values
pub fn ffs2dt_safe(
    datestr: &[c_char], /* I - date string: "YYYY-MM-DD" or "dd/mm/yy" */
    year: &mut c_int,   /* O - year (0 - 9999)                         */
    month: &mut c_int,  /* O - month (1 - 12)                          */
    day: &mut c_int,    /* O - day (1 - 31)                            */
    status: &mut c_int, /* IO - error status                           */
) -> c_int {
    todo!();
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
    todo!();
}

/*-----------------------------------------------------------------*/
/// Parse a date character string into date and time values
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffs2tm(
    datestr: *const c_char, /* I - date string: "YYYY-MM-DD"    */
    /*     or "YYYY-MM-DDThh:mm:ss.ddd" */
    /*     or "dd/mm/yy"                */
    year: *mut c_int,        /* O - year (0 - 9999)              */
    month: *mut c_int,       /* O - month (1 - 12)               */
    day: *mut c_int,         /* O - day (1 - 31)                 */
    hour: *const c_int,      /* I - hour (0 - 23)                */
    minute: *const c_int,    /* I - minute (0 - 59)              */
    second: *const c_double, /* I - second (0. - 60.9999999)     */
    status: *mut c_int,      /* IO - error status                */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let year = year.as_mut().expect(NULL_MSG);
        let month = month.as_mut().expect(NULL_MSG);
        let day = day.as_mut().expect(NULL_MSG);
        let hour = hour.as_ref().expect(NULL_MSG);
        let minute = minute.as_ref().expect(NULL_MSG);
        let second = second.as_ref().expect(NULL_MSG);

        raw_to_slice!(datestr);

        ffs2tm_safe(
            cast_slice(datestr),
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
    datestr: &[c_char], /* I - date string: "YYYY-MM-DD"    */
    /*     or "YYYY-MM-DDThh:mm:ss.ddd" */
    /*     or "dd/mm/yy"                */
    year: &mut c_int,   /* O - year (0 - 9999)              */
    month: &mut c_int,  /* O - month (1 - 12)               */
    day: &mut c_int,    /* O - day (1 - 31)                 */
    hour: &c_int,       /* I - hour (0 - 23)                */
    minute: &c_int,     /* I - minute (0 - 59)              */
    second: &c_double,  /* I - second (0. - 60.9999999)     */
    status: &mut c_int, /* IO - error status                */
) -> c_int {
    todo!();
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
/// his routine is included for backward compatibility with the Fortran FITSIO library.
/// Get current System DaTe (GMT if available)
pub fn ffgsdt_safe(
    day: &mut c_int,
    month: &mut c_int,
    year: &mut c_int,
    status: &mut c_int,
) -> c_int {
    todo!();
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

    if int_snprintf!(cval, FLEN_VALUE, "{:.*}", decim as usize, fval as f64,) < 0 {
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
            fval as f64,
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
                if int_snprintf!(cval, FLEN_VALUE, "{:.1E}", fval as f64,) < 0 {
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
        fmt_f64(fval as f64, decim as usize, 2)
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
    let cptr: *mut c_char = ptr::null_mut();

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
