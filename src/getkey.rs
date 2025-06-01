/*  This file, getkey.c, contains routines that read keywords from         */
/*  a FITS header.                                                         */

/*  The FITSIO software was written by William Pence at the High Energy    */
/*  Astrophysic Science Archive Research Center (HEASARC) at the NASA      */
/*  Goddard Space Flight Center.                                           */

use core::slice;
use std::ffi::CStr;
use std::{cmp, ptr};

use crate::aliases::ALLOCATIONS;
use crate::c_types::{c_char, c_int, c_long, c_short, c_uint, c_ulong, c_ushort, c_void};
use crate::helpers::vec_raw_parts::vec_into_raw_parts;

use bytemuck::{cast_slice, cast_slice_mut};

use crate::buffers::*;
use crate::fitscore::{
    ffc2d, ffc2dd, ffc2i, ffc2ii, ffc2j, ffc2jj, ffc2l, ffc2ll, ffc2r, ffc2s, ffc2uj, ffcmps_safe,
    ffghadll_safe, ffgidm_safe, ffgidt_safe, ffgiszll_safe, ffkeyn_safe, ffmahd_safe, ffpmsg_slice,
    ffpmsg_str, ffpsvc_safe, fftrec_safe, ffxmsg_safer, fits_is_compressed_image_safe,
};
use crate::raw_to_slice;
use crate::wrappers::*;
use crate::{KeywordDatatypeMut, bb, cs};
use crate::{fitsio::*, int_snprintf};
use crate::{fitsio2::*, slice_to_str};

/*--------------------------------------------------------------------------*/
/// returns the number of existing keywords (not counting the END keyword)
/// and the number of more keyword that will fit in the current header
/// without having to insert more FITS blocks.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffghsp(
    fptr: *mut fitsfile, /* I - FITS file pointer                     */
    nexist: *mut c_int,  /* O - number of existing keywords in header */
    nmore: *mut c_int,   /* O - how many more keywords will fit       */
    status: *mut c_int,  /* IO - error status                         */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        let nexist = nexist.as_mut();
        let nmore = nmore.as_mut();

        ffghsp_safe(fptr, nexist, nmore, status)
    }
}

/*--------------------------------------------------------------------------*/
/// returns the number of existing keywords (not counting the END keyword)
/// and the number of more keyword that will fit in the current header
/// without having to insert more FITS blocks.
pub fn ffghsp_safe(
    fptr: &mut fitsfile,        /* I - FITS file pointer                     */
    nexist: Option<&mut c_int>, /* O - number of existing keywords in header */
    nmore: Option<&mut c_int>,  /* O - how many more keywords will fit       */
    status: &mut c_int,         /* IO - error status                         */
) -> c_int {
    if *status > 0 {
        return *status;
    }

    if fptr.HDUposition != fptr.Fptr.curhdu {
        ffmahd_safe(fptr, (fptr.HDUposition) + 1, None, status);
    }

    if let Some(nexist) = nexist {
        *nexist = (((fptr.Fptr.headend)
            - (fptr.Fptr.get_headstart_as_slice()[fptr.Fptr.curhdu as usize]))
            / 80) as c_int;
    }

    if fptr.Fptr.datastart == DATA_UNDEFINED as LONGLONG {
        if let Some(nmore) = nmore {
            *nmore = -1; /* data not written yet, so room for any keywords */
        }
    } else {
        /* calculate space available between the data and the END card */
        if let Some(nmore) = nmore {
            *nmore = ((fptr.Fptr.datastart - fptr.Fptr.headend) / 80 - 1) as c_int;
        }
    }

    *status
}

/*--------------------------------------------------------------------------*/
/// Return the number of existing keywords and the position of the next
/// keyword that will be read.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffghps(
    fptr: *mut fitsfile,  /* I - FITS file pointer                     */
    nexist: *mut c_int,   /* O - number of existing keywords in header */
    position: *mut c_int, /* O - position of next keyword to be read   */
    status: *mut c_int,   /* IO - error status                         */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let nexist = nexist.as_mut();
        let position = position.as_mut();

        ffghps_safe(fptr, nexist, position, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Return the number of existing keywords and the position of the next
/// keyword that will be read.
pub fn ffghps_safe(
    fptr: &mut fitsfile,          /* I - FITS file pointer                     */
    nexist: Option<&mut c_int>,   /* O - number of existing keywords in header */
    position: Option<&mut c_int>, /* O - position of next keyword to be read   */
    status: &mut c_int,           /* IO - error status                         */
) -> c_int {
    if *status > 0 {
        return *status;
    }

    if fptr.HDUposition != fptr.Fptr.curhdu {
        ffmahd_safe(fptr, (fptr.HDUposition) + 1, None, status);
    }

    let headstart = fptr.Fptr.get_headstart_as_slice();

    if let Some(nexist) = nexist {
        *nexist = (((fptr.Fptr.headend) - (headstart[fptr.Fptr.curhdu as usize])) / 80) as c_int;
    }

    if let Some(position) = position {
        *position =
            (((fptr.Fptr.nextkey) - (headstart[fptr.Fptr.curhdu as usize])) / 80 + 1) as c_int;
    }
    *status
}

/*--------------------------------------------------------------------------*/
/// Function returns the position of the first null character (ASCII 0), if
/// any, in the current header.  Null characters are illegal, but the other
/// CFITSIO routines that read the header will not detect this error, because
/// the null gets interpreted as a normal end of string character.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffnchk(
    fptr: *mut fitsfile, /* I - FITS file pointer                     */
    status: *mut c_int,  /* IO - error status                         */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        ffnchk_safe(fptr, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Function returns the position of the first null character (ASCII 0), if
/// any, in the current header.  Null characters are illegal, but the other
/// CFITSIO routines that read the header will not detect this error, because
/// the null gets interpreted as a normal end of string character.
pub(crate) fn ffnchk_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer                     */
    status: &mut c_int,  /* IO - error status                         */
) -> c_int {
    let mut bytepos: LONGLONG = 0;
    let mut length = 0;
    let mut nullpos: c_int = 0;
    let mut block: [c_char; (IOBUFLEN + 1) as usize] = [0; (IOBUFLEN + 1) as usize];

    if *status > 0 {
        return *status;
    }

    if fptr.HDUposition != fptr.Fptr.curhdu {
        ffmahd_safe(fptr, (fptr.HDUposition) + 1, None, status);
    }

    let headstart = fptr.Fptr.get_headstart_as_slice();

    if fptr.Fptr.datastart == DATA_UNDEFINED as LONGLONG {
        return 0; /* Don't check a file that is just being created.  */
        /* It cannot contain nulls since CFITSIO wrote it. */
    }

    /* calculate number of blocks in the header */
    let nblock = (fptr.Fptr.datastart - headstart[fptr.Fptr.curhdu as usize]) / IOBUFLEN;

    bytepos = headstart[fptr.Fptr.curhdu as usize];
    ffmbyt_safe(fptr, bytepos, REPORT_EOF, status); /* move to read pos. */

    block[IOBUFLEN as usize] = 0;

    for ii in 0..(nblock as usize) {
        if ffgbyt(fptr, IOBUFLEN, cast_slice_mut(&mut block), status) > 0 {
            return 0; /* read error of some sort */
        }

        length = strlen_safe(&block);
        if length != IOBUFLEN as usize {
            nullpos = ((ii * IOBUFLEN as usize) + length + 1) as c_int;
            return nullpos;
        }
    }

    0
}

/*--------------------------------------------------------------------------*/
/// Move pointer to the specified absolute keyword position.  E.g. this keyword
/// will then be read by the next call to ffgnky.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffmaky(
    fptr: *mut fitsfile, /* I - FITS file pointer                    */
    nrec: c_int,         /* I - one-based keyword number to move to  */
    status: *mut c_int,  /* IO - error status                        */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        ffmaky_safe(fptr, nrec, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Move pointer to the specified absolute keyword position.  E.g. this keyword
/// will then be read by the next call to ffgnky.
pub fn ffmaky_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer                    */
    nrec: c_int,         /* I - one-based keyword number to move to  */
    status: &mut c_int,  /* IO - error status                        */
) -> c_int {
    if fptr.HDUposition != fptr.Fptr.curhdu {
        ffmahd_safe(fptr, (fptr.HDUposition) + 1, None, status);
    }

    let headstart = fptr.Fptr.get_headstart_as_slice();
    fptr.Fptr.nextkey = headstart[fptr.Fptr.curhdu as usize] + ((nrec as LONGLONG - 1) * 80);
    *status
}

/*--------------------------------------------------------------------------*/
/// move pointer to the specified keyword position relative to the current
/// position.  E.g. this keyword  will then be read by the next call to ffgnky.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffmrky(
    fptr: *mut fitsfile, /* I - FITS file pointer                   */
    nmove: c_int,        /* I - relative number of keywords to move */
    status: *mut c_int,  /* IO - error status                       */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        ffmrky_safe(fptr, nmove, status)
    }
}

/*--------------------------------------------------------------------------*/
/// move pointer to the specified keyword position relative to the current
/// position.  E.g. this keyword  will then be read by the next call to ffgnky.
pub fn ffmrky_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer                   */
    nmove: c_int,        /* I - relative number of keywords to move */
    status: &mut c_int,  /* IO - error status                       */
) -> c_int {
    if fptr.HDUposition != fptr.Fptr.curhdu {
        ffmahd_safe(fptr, (fptr.HDUposition) + 1, None, status);
    }

    fptr.Fptr.nextkey += nmove as LONGLONG * 80;

    *status
}

/*--------------------------------------------------------------------------*/
/// Read the next keyword from the header - used internally by cfitsio
pub(crate) fn ffgnky(
    fptr: &mut fitsfile, /* I - FITS file pointer     */
    card: &mut [c_char], /* O - card string           */
    status: &mut c_int,  /* IO - error status         */
) -> c_int {
    let mut message: [c_char; FLEN_ERRMSG] = [0; FLEN_ERRMSG];

    if *status > 0 {
        return *status;
    }

    if fptr.HDUposition != fptr.Fptr.curhdu {
        ffmahd_safe(fptr, (fptr.HDUposition) + 1, None, status);
    }

    card[0] = 0; /* make sure card is terminated, even affer read error */

    /*
    Check that nextkey points to a legal keyword position.  Note that headend
    is the current end of the header, i.e., the position where a new keyword
    would be appended, however, if there are more than 1 FITS block worth of
    blank keywords at the end of the header (36 keywords per 2880 byte block)
    then the actual physical END card must be located at a starting position
    which is just 2880 bytes prior to the start of the data unit.
    */
    let bytepos = fptr.Fptr.nextkey;
    let endhead = cmp::max(fptr.Fptr.headend, fptr.Fptr.datastart - 2880);

    /* nextkey must be < endhead and > than  headstart */
    let headstart = fptr.Fptr.get_headstart_as_slice();
    if bytepos > endhead || bytepos < headstart[fptr.Fptr.curhdu as usize] {
        let nrec = ((bytepos - headstart[fptr.Fptr.curhdu as usize]) / 80 + 1) as c_int;
        int_snprintf!(
            &mut message,
            FLEN_ERRMSG,
            "Cannot get keyword number {}.  It does not exist.",
            nrec,
        );
        ffpmsg_slice(&message);

        *status = KEY_OUT_BOUNDS;
        return *status;
    }

    ffmbyt_safe(fptr, bytepos, REPORT_EOF, status); /* move to read pos. */

    card[80] = 0; /* make sure card is terminate, even if ffgbyt fails */

    if ffgbyt(fptr, 80, cast_slice_mut(card), status) <= 0 {
        fptr.Fptr.nextkey += 80; /* increment pointer to next keyword */

        /* strip off trailing blanks with terminated string */
        let mut jj: isize = 79; //isize because of >= comparison
        while jj >= 0 && card[jj as usize] == bb(b' ') {
            jj -= 1;
        }
        card[(jj + 1) as usize] = 0;
    }
    *status
}

/*--------------------------------------------------------------------------*/
/// Return the next keyword that matches one of the names in inclist
/// but does not match any of the names in exclist.  The search
/// goes from the current position to the end of the header, only.
/// Wild card characters may be used in the name lists ('*', '?' and '#').
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgnxk(
    fptr: *mut fitsfile,           /* I - FITS file pointer              */
    inclist: *const *const c_char, /* I - list of included keyword names */
    ninc: c_int,                   /* I - number of names in inclist     */
    exclist: *const *const c_char, /* I - list of excluded keyword names */
    nexc: c_int,                   /* I - number of names in exclist     */
    card: *mut c_char,             /* O - first matching keyword         */
    status: *mut c_int,            /* IO - error status                  */
) -> c_int {
    unsafe {
        let mut casesn: c_int = 0;
        let mut is_match: c_int = 0;
        let mut exact: c_int = 0;
        let mut namelen: c_int = 0;

        let mut keybuf: [c_char; FLEN_CARD] = [0; FLEN_CARD];
        let mut keyname: [c_char; FLEN_KEYWORD] = [0; FLEN_KEYWORD];

        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        let nexc = nexc as usize;
        let ninc = ninc as usize;

        let inclist = slice::from_raw_parts(inclist, ninc);
        let exclist = slice::from_raw_parts(exclist, nexc);
        let card = slice::from_raw_parts_mut(card, FLEN_CARD);

        card[0] = 0;

        if *status > 0 {
            return *status;
        }

        casesn = FALSE as c_int;

        /* get next card, and return with an error if hit end of header */
        while ffgcrd_safe(fptr, cs!(c"*"), &mut keybuf, status) <= 0 {
            ffgknm_safe(&keybuf, &mut keyname, &mut namelen, status); /* get the keyword name */

            /* does keyword match any names in the include list? */
            for ii in 0..ninc {
                let inclist_ii = cast_slice(CStr::from_ptr(inclist[ii]).to_bytes_with_nul());

                ffcmps_safe(inclist_ii, &keyname, casesn, &mut is_match, &mut exact);
                if is_match != 0 {
                    /* does keyword match any names in the exclusion list? */
                    let mut jj = 0;
                    while jj < nexc {
                        let exclist_ii =
                            cast_slice(CStr::from_ptr(exclist[jj]).to_bytes_with_nul());
                        ffcmps_safe(exclist_ii, &keyname, casesn, &mut is_match, &mut exact);
                        if is_match != 0 {
                            break;
                        }
                        jj += 1;
                    }

                    if jj >= nexc {
                        /* not in exclusion list, so return this keyword */
                        strcat_safe(card, &keybuf);
                        return *status;
                    }
                }
            }
        }
        *status
    }
}

/*--------------------------------------------------------------------------*/
/// Read (get) the keyword value and comment from the FITS header.
/// Reads a keyword value with the datatype specified by the 2nd argument.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgky(
    fptr: *mut fitsfile,    /* I - FITS file pointer        */
    datatype: c_int,        /* I - datatype of the value    */
    keyname: *const c_char, /* I - name of keyword to read  */
    value: *mut c_void,     /* O - keyword value            */
    comm: *mut c_char,      /* O - keyword comment          */
    status: *mut c_int,     /* IO - error status            */
) -> c_int {
    unsafe {
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);
        raw_to_slice!(keyname);

        let comm: Option<&mut [c_char; FLEN_COMMENT]> = match comm.is_null() {
            true => None,
            false => Some(
                slice::from_raw_parts_mut(comm, FLEN_COMMENT)
                    .try_into()
                    .unwrap(),
            ),
        };

        let datatype_with_value = KeywordDatatypeMut::from_datatype(datatype, value);

        ffgky_safe(fptr, datatype_with_value, keyname, comm, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Read (get) the keyword value and comment from the FITS header.
/// Reads a keyword value with the datatype specified by the 2nd argument.
///
/// This has been modified heavily from the original for safety.
pub fn ffgky_safe(
    fptr: &mut fitsfile,                           /* I - FITS file pointer        */
    datatype: KeywordDatatypeMut,                  /* I - datatype of the value    */
    keyname: &[c_char],                            /* I - name of keyword to read  */
    mut comm: Option<&mut [c_char; FLEN_COMMENT]>, /* O - keyword comment          */
    status: &mut c_int,                            /* IO - error status            */
) -> c_int {
    let mut longval: LONGLONG = 0;
    let mut ulongval: ULONGLONG = 0;
    let doubleval: f64 = 0.0;

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    match datatype {
        KeywordDatatypeMut::TSTRING(value) => {
            ffgkys_safe(fptr, keyname, value, comm, status);
        }
        KeywordDatatypeMut::TBYTE(value) => {
            if ffgkyjj_safe(fptr, keyname, &mut longval, comm, status) <= 0 {
                if longval > u8::MAX as LONGLONG || longval < 0 {
                    *status = NUM_OVERFLOW;
                } else {
                    *(value) = longval as u8;
                };
            };
        }
        KeywordDatatypeMut::TSBYTE(value) => {
            if ffgkyjj_safe(fptr, keyname, &mut longval, comm, status) <= 0 {
                if longval > 127i64 || longval < -128 {
                    *status = NUM_OVERFLOW;
                } else {
                    *(value) = longval as c_char;
                };
            };
        }
        KeywordDatatypeMut::TUSHORT(value) => {
            if ffgkyjj_safe(fptr, keyname, &mut longval, comm, status) <= 0 {
                if longval > (c_ushort::MAX as LONGLONG) || longval < 0 {
                    *status = NUM_OVERFLOW;
                } else {
                    *(value) = longval as c_ushort;
                };
            };
        }
        KeywordDatatypeMut::TSHORT(value) => {
            if ffgkyjj_safe(fptr, keyname, &mut longval, comm, status) <= 0 {
                if longval > c_short::MAX as LONGLONG || longval < c_short::MIN as LONGLONG {
                    *status = NUM_OVERFLOW;
                } else {
                    *(value) = longval as c_short;
                };
            };
        }
        KeywordDatatypeMut::TUINT(value) => {
            if ffgkyjj_safe(fptr, keyname, &mut longval, comm, status) <= 0 {
                if longval > (c_uint::MAX as LONGLONG) || longval < 0 {
                    *status = NUM_OVERFLOW;
                } else {
                    *(value) = longval as c_uint;
                };
            };
        }
        KeywordDatatypeMut::TINT(value) => {
            if ffgkyjj_safe(fptr, keyname, &mut longval, comm, status) <= 0 {
                if longval > c_int::MAX as LONGLONG || longval < c_int::MIN as LONGLONG {
                    *status = NUM_OVERFLOW;
                } else {
                    *(value) = longval as c_int;
                };
            };
        }
        KeywordDatatypeMut::TLOGICAL(value) => {
            ffgkyl_safe(fptr, keyname, value, comm, status);
        }
        KeywordDatatypeMut::TULONG(value) => {
            if ffgkyujj_safe(fptr, keyname, &mut ulongval, comm, status) <= 0 {
                if ulongval > c_ulong::MAX as ULONGLONG {
                    *status = NUM_OVERFLOW;
                } else {
                    *(value) = ulongval as c_ulong;
                };
            };
        }
        KeywordDatatypeMut::TLONG(value) => {
            if ffgkyjj_safe(fptr, keyname, &mut longval, comm.as_deref_mut(), status) <= 0 {
                if longval > LONG_MAX as LONGLONG || longval < LONG_MIN as LONGLONG {
                    *status = NUM_OVERFLOW;
                } else {
                    *(value) = longval as c_long;
                };
            }
        }
        KeywordDatatypeMut::TULONGLONG(value) => {
            ffgkyujj_safe(fptr, keyname, value, comm, status);
        }
        KeywordDatatypeMut::TLONGLONG(value) => {
            ffgkyjj_safe(fptr, keyname, value, comm, status);
        }
        KeywordDatatypeMut::TFLOAT(value) => {
            ffgkye_safe(fptr, keyname, value, comm, status);
        }
        KeywordDatatypeMut::TDOUBLE(value) => {
            ffgkyd_safe(fptr, keyname, value, comm, status);
        }
        KeywordDatatypeMut::TCOMPLEX(value) => {
            ffgkyc_safe(fptr, keyname, value, comm, status);
        }
        KeywordDatatypeMut::TDBLCOMPLEX(value) => {
            ffgkym_safe(fptr, keyname, value, comm, status);
        }
        _ => {
            *status = BAD_DATATYPE;
        }
    }

    *status
}

/*--------------------------------------------------------------------------*/
/// Read (get) the named keyword, returning the keyword value and comment.
///
/// The value is just the literal string of characters in the value field
/// of the keyword.  In the case of a string valued keyword, the returned
/// value includes the leading and closing quote characters.  The value may be
/// up to 70 characters long, and the comment may be up to 72 characters long.
/// If the keyword has no value (no equal sign in column 9) then a null value
/// is returned.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgkey(
    fptr: *mut fitsfile,    /* I - FITS file pointer        */
    keyname: *const c_char, /* I - name of keyword to read  */
    keyval: *mut c_char,    /* O - keyword value            */
    comm: *mut c_char,      /* O - keyword comment          */
    status: *mut c_int,     /* IO - error status            */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let keyval = slice::from_raw_parts_mut(keyval, FLEN_VALUE);

        raw_to_slice!(keyname);

        let comm: Option<&mut [c_char; FLEN_COMMENT]> = match comm.is_null() {
            true => None,
            false => Some(
                slice::from_raw_parts_mut(comm, FLEN_COMMENT)
                    .try_into()
                    .unwrap(),
            ),
        };

        ffgkey_safe(fptr, keyname, keyval, comm, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Read (get) the named keyword, returning the keyword value and comment.
///
/// The value is just the literal string of characters in the value field
/// of the keyword.  In the case of a string valued keyword, the returned
/// value includes the leading and closing quote characters.  The value may be
/// up to 70 characters long, and the comment may be up to 72 characters long.
/// If the keyword has no value (no equal sign in column 9) then a null value
/// is returned.
pub fn ffgkey_safe(
    fptr: &mut fitsfile,                       /* I - FITS file pointer        */
    keyname: &[c_char],                        /* I - name of keyword to read  */
    keyval: &mut [c_char],                     /* O - keyword value            */
    comm: Option<&mut [c_char; FLEN_COMMENT]>, /* O - keyword comment          */
    status: &mut c_int,                        /* IO - error status            */
) -> c_int {
    let mut card: [c_char; FLEN_CARD] = [0; FLEN_CARD];

    keyval[0] = 0;

    let mut comm = comm;
    if let Some(c) = comm.as_deref_mut() {
        c[0] = 0;
    }

    if *status > 0 {
        return *status;
    }

    if ffgcrd_safe(fptr, keyname, &mut card, status) > 0 {
        /* get the 80-byte card */
        return *status;
    }

    ffpsvc_safe(&card, keyval, comm, status); /* parse the value and comment */

    *status
}

/*--------------------------------------------------------------------------*/
/// Read (get) the nrec-th keyword, returning the entire keyword card up to
/// 80 characters long.  The first keyword in the header has nrec = 1, not 0.
/// The returned card value is null terminated with any trailing blank
/// characters removed.  If nrec = 0, then this routine simply moves the
/// current header pointer to the top of the header.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgrec(
    fptr: *mut fitsfile, /* I - FITS file pointer          */
    nrec: c_int,         /* I - number of keyword to read  */
    card: *mut c_char,   /* O - keyword card               */
    status: *mut c_int,  /* IO - error status              */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        if *status > 0 {
            return *status;
        }

        let card = if card.is_null() {
            None
        } else {
            Some(slice::from_raw_parts_mut(card, FLEN_CARD))
        };

        ffgrec_safe(fptr, nrec, card, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Read (get) the nrec-th keyword, returning the entire keyword card up to
/// 80 characters long.  The first keyword in the header has nrec = 1, not 0.
/// The returned card value is null terminated with any trailing blank
/// characters removed.  If nrec = 0, then this routine simply moves the
/// current header pointer to the top of the header.
pub fn ffgrec_safe(
    fptr: &mut fitsfile,         /* I - FITS file pointer          */
    nrec: c_int,                 /* I - number of keyword to read  */
    card: Option<&mut [c_char]>, /* O - keyword card               */
    status: &mut c_int,          /* IO - error status              */
) -> c_int {
    if *status > 0 {
        return *status;
    }
    if nrec == 0 {
        ffmaky_safe(fptr, 1, status); /* simply move to beginning of header */

        if let Some(card) = card {
            card[0] = 0; /* and return null card */
        };
    } else if nrec > 0 {
        let card = card.unwrap(); // WARNING: Original code doesn't check for null ptr here
        ffmaky_safe(fptr, nrec, status);
        ffgnky(fptr, card, status);
    }
    *status
}

/*--------------------------------------------------------------------------*/
/// Read (get) the named keyword, returning the entire keyword card up to
/// 80 characters long.  
/// The returned card value is null terminated with any trailing blank
/// characters removed.
///
/// If the input name contains wild cards ('?' matches any single char
/// and '*' matches any sequence of chars, # matches any string of decimal
/// digits) then the search ends once the end of header is reached and does
/// not automatically resume from the top of the header.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgcrd(
    fptr: *mut fitsfile, /* I - FITS file pointer        */
    name: *const c_char, /* I - name of keyword to read  */
    card: *mut c_char,   /* O - keyword card             */
    status: *mut c_int,  /* IO - error status            */
) -> c_int {
    unsafe {
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);
        let name: &[c_char] = cast_slice(CStr::from_ptr(name).to_bytes_with_nul());
        let card: &mut [c_char; FLEN_CARD] = slice::from_raw_parts_mut(card, FLEN_CARD)
            .try_into()
            .unwrap();

        ffgcrd_safe(fptr, name, card, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Read (get) the named keyword, returning the entire keyword card up to
/// 80 characters long.  
/// The returned card value is null terminated with any trailing blank
/// characters removed.
///
/// If the input name contains wild cards ('?' matches any single char
/// and '*' matches any sequence of chars, # matches any string of decimal
/// digits) then the search ends once the end of header is reached and does
/// not automatically resume from the top of the header.
pub fn ffgcrd_safe(
    fptr: &mut fitsfile,            /* I - FITS file pointer        */
    name: &[c_char],                /* I - name of keyword to read  */
    card: &mut [c_char; FLEN_CARD], /* O - keyword card             */
    status: &mut c_int,             /* IO - error status            */
) -> c_int {
    let mut nkeys: isize = 0;
    let mut nextkey: isize = 0;
    let mut ntodo: isize = 0;
    let mut namelen: usize = 0;
    let mut namelen_limit: c_int = 0;
    let mut namelenminus1: c_int = 0;
    let mut cardlen: usize = 0;
    let mut ii: usize = 0;
    let mut wild: c_int = 0;
    let mut match_item: c_int = 0;
    let mut exact: c_int = 0;
    let mut hier: c_int = 0;
    let mut keyname: [c_char; FLEN_KEYWORD] = [0; FLEN_KEYWORD];
    let mut cardname: [c_char; FLEN_KEYWORD] = [0; FLEN_KEYWORD];
    let mut ptr1: usize;
    let mut ptr2: usize;
    let mut gotstar: Option<usize>;

    if *status > 0 {
        return *status;
    }

    keyname[0] = 0;

    while name[ii] == bb(b' ') {
        /* skip leading blanks in name */
        ii += 1;
    }
    strncat_safe(&mut keyname, &name[ii..], FLEN_KEYWORD - 1);

    namelen = strlen_safe(&keyname);

    while namelen > 0 && keyname[namelen - 1] == bb(b' ') {
        namelen -= 1; /* ignore trailing blanks in name */
    }
    keyname[namelen] = 0; /* terminate the name */

    let mut ii = 0;
    while ii < namelen {
        keyname[ii] = toupper(keyname[ii]); /*  make upper case  */
        ii += 1;
    }
    if FSTRNCMP(cs!(c"HIERARCH"), &keyname, 8) == 0 {
        if namelen == 8 {
            /* special case: just looking for any HIERARCH keyword */
            hier = 1;
        } else {
            /* ignore the leading HIERARCH and look for the 'real' name */
            /* starting with first non-blank character following HIERARCH */
            ptr1 = 0;
            ptr2 = 8;
            {
                while keyname[ptr2] == bb(b' ') {
                    ptr2 += 1;
                }
            }
            namelen = 0;
            while keyname[ptr2] > 0 {
                keyname[ptr1] = keyname[ptr2];
                ptr1 += 1;
                ptr2 += 1;
                namelen += 1;
            }
            keyname[ptr1] = 0;
        }
    }

    /* does input name contain wild card chars?  ('?',  '*', or '#') */
    /* wild cards are currently not supported with HIERARCH keywords */

    namelen_limit = namelen as c_int;
    let gotstar = strchr_safe(&keyname, bb(b'*'));
    if namelen < 9 && strchr_safe(&keyname, bb(b'?')).is_some()
        || (gotstar.is_some() || strchr_safe(&keyname, bb(b'#')).is_some())
    {
        wild = 1;

        /* if we found a '*' wild card in the name, there might be */
        /* more than one.  Support up to 2 '*' in the template. */
        /* Thus we need to compare keywords whose names have at least */
        /* namelen - 2 characters.                                   */
        if gotstar.is_some() {
            namelen_limit -= 2;
        }
    } else {
        wild = 0;
    }

    let mut nk = 0;
    let mut n_k = 0;
    ffghps_safe(fptr, Some(&mut n_k), Some(&mut nk), status); /* get no. keywords and position */

    nkeys = n_k as isize;
    nextkey = nk as isize;

    namelenminus1 = cmp::max(namelen - 1, 1) as c_int;
    ntodo = nkeys + 1 - nextkey; /* first, read from next keyword to end */

    let mut jj: usize = 0;
    while jj < 2 {
        let mut kk = 0;
        while kk < ntodo {
            ffgnky(fptr, card, status); /* get next keyword */

            if hier > 0 {
                if FSTRNCMP(cs!(c"HIERARCH"), card, 8) == 0 {
                    return *status; /* found a HIERARCH keyword */
                }
            } else {
                let mut clen = 0;
                ffgknm_safe(card, &mut cardname, &mut clen, status); /* get the keyword name */
                cardlen = clen as usize;

                if (cardlen as c_int) >= namelen_limit {
                    /* can't match if card < name */

                    /* if there are no wild cards, lengths must be the same */
                    if !(wild == 0 && cardlen != namelen) {
                        let mut ii: usize = 0;
                        while ii < cardlen {
                            /* make sure keyword is in uppercase */
                            if cardname[ii] > 96 {
                                /* This assumes the ASCII character set in which */
                                /* upper case characters start at ASCII(97)  */
                                /* Timing tests showed that this is 20% faster */
                                /* than calling the isupper function.          */

                                cardname[ii] = toupper(cardname[ii]); /* make upper case */
                            }
                            ii += 1;
                        }

                        if wild > 0 {
                            ffcmps_safe(&keyname, &cardname, 1, &mut match_item, &mut exact);
                            if match_item > 0 {
                                return *status; /* found a matching keyword */
                            }
                        } else if keyname[namelenminus1 as usize]
                            == cardname[namelenminus1 as usize]
                        {
                            /* test the last character of the keyword name first, on */
                            /* the theory that it is less likely to match then the first */
                            /* character since many keywords begin with 'T', for example */

                            if FSTRNCMP(&keyname, &cardname, namelenminus1 as usize) == 0 {
                                return *status; /* found the matching keyword */
                            }
                        } else if namelen == 0 && cardlen == 0 {
                            /* matched a blank keyword */
                            return *status;
                        }
                    }
                }
            }
            kk += 1;
        }

        if wild > 0 || jj == 1 {
            break; /* stop at end of header if template contains wildcards */
        }
        ffmaky_safe(fptr, 1, status); /* reset pointer to beginning of header */
        ntodo = nextkey - 1; /* number of keyword to read */

        jj += 1;
    }

    *status = KEY_NO_EXIST; /* couldn't find the keyword */
    *status
}

/*--------------------------------------------------------------------------*/
/// Read (get) the next keyword record that contains the input character string,
/// returning the entire keyword card up to 80 characters long.
/// The returned card value is null terminated with any trailing blank
/// characters removed.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgstr(
    fptr: *mut fitsfile,   /* I - FITS file pointer        */
    string: *const c_char, /* I - string to match  */
    card: *mut c_char,     /* O - keyword card             */
    status: *mut c_int,    /* IO - error status            */
) -> c_int {
    unsafe {
        let mut nkeys: c_int = 0;
        let mut nextkey: c_int = 0;
        let mut ntodo: c_int = 0;

        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        let card = slice::from_raw_parts_mut(card, FLEN_CARD);
        raw_to_slice!(string);

        if *status > 0 {
            return *status;
        }

        let stringlen = strlen_safe(string);
        if stringlen > 80 {
            *status = KEY_NO_EXIST; /* matching string is too long to exist */
            return *status;
        }

        ffghps_safe(fptr, Some(&mut nkeys), Some(&mut nextkey), status); /* get no. keywords and position */
        ntodo = nkeys - nextkey + 1; /* first, read from next keyword to end */

        for jj in 0..2 {
            for kk in 0..(ntodo as usize) {
                ffgnky(fptr, card, status); /* get next keyword */
                if strstr_safe(card, string).is_some() {
                    return *status; /* found the matching string */
                }
            }

            ffmaky_safe(fptr, 1, status); /* reset pointer to beginning of header */
            ntodo = nextkey - 1; /* number of keyword to read */
        }

        *status = KEY_NO_EXIST; /* couldn't find the keyword */
        *status
    }
}

/*--------------------------------------------------------------------------*/
/// Return the name of the keyword, and the name length.  This supports the
/// ESO HIERARCH convention where keyword names may be > 8 characters long.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgknm(
    card: *const c_char, /* I - keyword card                   */
    name: *mut c_char,   /* O - name of the keyword            */
    length: *mut c_int,  /* O - length of the keyword name     */
    status: *mut c_int,  /* IO - error status                  */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let length = length.as_mut().expect(NULL_MSG);

        let card = slice::from_raw_parts(card, FLEN_CARD).try_into().unwrap();
        let name = slice::from_raw_parts_mut(name, FLEN_KEYWORD)
            .try_into()
            .unwrap();

        ffgknm_safe(card, name, length, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Return the name of the keyword, and the name length.  This supports the
/// ESO HIERARCH convention where keyword names may be > 8 characters long.
pub fn ffgknm_safe(
    card: &[c_char; FLEN_CARD],        /* I - keyword card                   */
    name: &mut [c_char; FLEN_KEYWORD], /* O - name of the keyword            */
    length: &mut c_int,                /* O - length of the keyword name     */
    status: &mut c_int,                /* IO - error status                  */
) -> c_int {
    let mut ptr1: usize = 0;

    let mut ii: usize;
    let namelength = FLEN_KEYWORD - 1;

    name[0] = 0;
    *length = 0;

    /* support for ESO HIERARCH keywords; find the '=' */
    if FSTRNCMP(card, cs!(c"HIERARCH "), 9) == 0 {
        let ptr2 = strchr_safe(card, bb(b'=')); /* no value indicator ??? */

        if ptr2.is_none() {
            /* this probably indicates an error, so just return FITS name */
            strcat_safe(name, cs!(c"HIERARCH"));
            *length = 8;
            return *status;
        }

        let ptr2 = ptr2.unwrap();

        /* find the start and end of the HIERARCH name */
        ptr1 = 9;

        while card[ptr1] == bb(b' ') {
            /* skip spaces */
            ptr1 += 1;
        }

        strncat_safe(name, &card[ptr1..], ptr2 - ptr1);

        let mut ii = ptr2 - ptr1;

        while ii > 0 && name[ii - 1] == bb(b' ') {
            /* remove trailing spaces */
            ii -= 1;
        }

        name[ii] = 0;
        *length = ii as c_int;
    } else {
        let ii = 0;
        for ii in 0..namelength {
            /* look for string terminator, or a blank */
            if card[ii] != bb(b' ') && card[ii] != bb(b'=') && card[ii] != 0 {
                name[ii] = card[ii];
            } else {
                name[ii] = 0;
                *length = ii as c_int;
                return *status;
            };
        }
        /* if we got here, keyword is namelength characters long */
        name[namelength] = 0;
        *length = namelength as c_int;
    }
    *status
}

/*--------------------------------------------------------------------------*/
/// Read (get) the units string from the comment field of the existing
/// keyword. This routine uses a local FITS convention (not defined in the
/// official FITS standard) in which the units are enclosed in
/// square brackets following the '/' comment field delimiter, e.g.:
///
/// KEYWORD =                   12 / [kpc] comment string goes here
/// WARNING: No definition in spec of length of units
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgunt(
    fptr: *mut fitsfile,    /* I - FITS file pointer         */
    keyname: *const c_char, /* I - name of keyword to read   */
    unit: *mut c_char,      /* O - keyword units             */
    status: *mut c_int,     /* IO - error status             */
) -> c_int {
    unsafe {
        let mut valstring: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
        let mut comm: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];

        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        let unit = slice::from_raw_parts_mut(unit, FLEN_COMMENT); //WARNING: See above

        raw_to_slice!(keyname);

        if *status > 0 {
            return *status;
        }

        ffgkey_safe(fptr, keyname, &mut valstring, Some(&mut comm), status); /* read the keyword */

        if comm[0] == bb(b'[') {
            let loc = strchr_safe(&comm, bb(b']')); /*  find the closing bracket */
            if let Some(loc) = loc {
                comm[loc] = 0; /*  terminate the string */
            }

            strcpy_safe(unit, &comm[1..]); /*  copy the string */
        } else {
            unit[0] = 0;
        }

        *status
    }
}

/*--------------------------------------------------------------------------*/
/// Get KeYword with a String value:
/// Read (get) a simple string valued keyword.  The returned value may be up to
/// 68 chars long ( + 1 null terminator char).  The routine does not support the
/// HEASARC convention for continuing long string values over multiple keywords.
/// The ffgkls routine may be used to read long continued strings. The returned
/// comment string may be up to 69 characters long (including null terminator).
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgkys(
    fptr: *mut fitsfile,    /* I - FITS file pointer         */
    keyname: *const c_char, /* I - name of keyword to read   */
    value: *mut c_char,     /* O - keyword value             */
    comm: *mut c_char,      /* O - keyword comment           */
    status: *mut c_int,     /* IO - error status             */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        raw_to_slice!(keyname);

        let c: Option<&mut [c_char; FLEN_COMMENT]> = if comm.is_null() {
            None
        } else {
            Some(
                slice::from_raw_parts_mut(comm, FLEN_COMMENT)
                    .try_into()
                    .unwrap(),
            )
        };

        let value = slice::from_raw_parts_mut(value, 69);

        ffgkys_safe(fptr, keyname, value, c, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Get KeYword with a String value:
/// Read (get) a simple string valued keyword.  The returned value may be up to
/// 68 chars long ( + 1 null terminator char).  The routine does not support the
/// HEASARC convention for continuing long string values over multiple keywords.
/// The ffgkls routine may be used to read long continued strings. The returned
/// comment string may be up to 69 characters long (including null terminator).
pub(crate) fn ffgkys_safe(
    fptr: &mut fitsfile,                       /* I - FITS file pointer         */
    keyname: &[c_char],                        /* I - name of keyword to read   */
    value: &mut [c_char],                      /* O - keyword value             */
    comm: Option<&mut [c_char; FLEN_COMMENT]>, /* O - keyword comment           */
    status: &mut c_int,                        /* IO - error status             */
) -> c_int {
    let mut valstring: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];

    if *status > 0 {
        return *status;
    }

    ffgkey_safe(fptr, keyname, &mut valstring, comm, status); /* read the keyword */
    value[0] = 0;
    ffc2s(&valstring, value, status); /* remove quotes from string */

    *status
}

/*--------------------------------------------------------------------------*/
/// Get the length of the keyword value string.
/// This routine explicitly supports the CONTINUE convention for long string values.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgksl(
    fptr: *mut fitsfile,    /* I - FITS file pointer             */
    keyname: *const c_char, /* I - name of keyword to read       */
    length: *mut c_int,     /* O - length of the string value    */
    status: *mut c_int,     /* IO - error status                 */
) -> c_int {
    unsafe {
        let mut dummy = 0;

        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let length = length.as_mut().expect(NULL_MSG);

        raw_to_slice!(keyname);

        if *status > 0 {
            return *status;
        }

        ffgkcsl_safe(fptr, keyname, length, &mut dummy, status);

        *status
    }
}

/*--------------------------------------------------------------------------*/
/// Get the length of the keyword value string and comment string.
/// This routine explicitly supports the CONTINUE convention for long string values.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgkcsl(
    fptr: *mut fitsfile,    /* I - FITS file pointer             */
    keyname: *const c_char, /* I - name of keyword to read       */
    length: *mut c_int,     /* O - length of the string value    */
    comlength: *mut c_int,  /* O - length of comment string      */
    status: *mut c_int,     /* IO - error status                 */
) -> c_int {
    unsafe {
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);
        let length = length.as_mut().expect(NULL_MSG);
        let comlength = comlength.as_mut().expect(NULL_MSG);

        raw_to_slice!(keyname);

        if *status > 0 {
            return *status;
        }

        let dummy_vallen: c_int = 0;
        let dummy_comlen: c_int = 0;

        ffglkut(
            fptr, keyname, 0, 0, 0, None, length, None, comlength, status,
        );

        *status
    }
}

/*--------------------------------------------------------------------------*/
/// Get the length of the keyword value string and comment string.
/// This routine explicitly supports the CONTINUE convention for long string values.
pub(crate) fn ffgkcsl_safe(
    fptr: &mut fitsfile,   /* I - FITS file pointer             */
    keyname: &[c_char],    /* I - name of keyword to read       */
    length: &mut c_int,    /* O - length of the string value    */
    comlength: &mut c_int, /* O - length of comment string      */
    status: &mut c_int,    /* IO - error status                 */
) -> c_int {
    if *status > 0 {
        return *status;
    }

    let dummy_vallen: c_int = 0;
    let dummy_comlen: c_int = 0;

    ffglkut(
        fptr, keyname, 0, 0, 0, None, length, None, comlength, status,
    );

    *status
}

/*--------------------------------------------------------------------------*/
/// This is the original routine for reading long string keywords that use
/// the CONTINUE keyword convention.  In 2016 a new routine called
/// ffgsky / fits_read_string_key was added, which may provide a more
/// convenient user interface  for most applications.
///
/// Get Keyword with possible Long String value:
/// Read (get) the named keyword, returning the value and comment.
/// The returned value string may be arbitrarily long (by using the HEASARC
/// convention for continuing long string values over multiple keywords) so
/// this routine allocates the required memory for the returned string value.
/// It is up to the calling routine to free the memory once it is finished
/// with the value string.  The returned comment string may be up to 69
/// characters long.
// #[deprecated]
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgkls(
    fptr: *mut fitsfile,     /* I - FITS file pointer         */
    keyname: *const c_char,  /* I - name of keyword to read       */
    value: *mut *mut c_char, /* O - pointer to keyword value      */
    comm: *mut c_char,       /* O - keyword comment (may be NULL) */
    status: *mut c_int,      /* IO - error status             */
) -> c_int {
    unsafe {
        let mut valstring: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
        let mut nextcomm: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];
        let mut card: [c_char; FLEN_CARD] = [0; FLEN_CARD];

        let mut addCommDelim = false;
        let mut keynum = 0;

        let mut commspace = 0;
        let mut len;

        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        raw_to_slice!(keyname);

        let mut c: Option<&mut [c_char; FLEN_COMMENT]> = if comm.is_null() {
            None
        } else {
            Some(
                slice::from_raw_parts_mut(comm, FLEN_COMMENT)
                    .try_into()
                    .unwrap(),
            )
        };

        if *status > 0 {
            return *status;
        }

        *value = ptr::null_mut(); /* initialize a null pointer in case of error */

        card[0] = 0;
        if let Some(comm) = c.as_deref_mut() {
            comm[0] = 0;
        }

        ffgcrd_safe(fptr, keyname, &mut card, status);

        if *status > 0 {
            return *status;
        }

        if strlen_safe(&card) < FLEN_CARD - 1 {
            addCommDelim = true;
        }

        ffpsvc_safe(&card, &mut valstring, c, status);

        if *status > 0 {
            return *status;
        }

        let mut c = if comm.is_null() {
            None
        } else {
            Some(slice::from_raw_parts_mut(comm, 69))
        };

        if let Some(c) = c.as_deref_mut() {
            /* remaining space in comment string */
            commspace = FLEN_COMMENT - 1 - strlen_safe(c);
        }

        if valstring[0] == 0 {
            /* null value string? */
            // HEAP ALLOCATION

            /* allocate and return a null string */
            let mut v: Vec<c_char> = vec![0; 1];
            v[0] = 0;
        } else {
            /* allocate space,  plus 1 for null */
            // HEAP ALLOCATION
            let mut v: Vec<c_char> = vec![0; strlen_safe(&valstring) + 1];

            ffc2s(&valstring, &mut v, status); /* convert string to value */
            len = strlen_safe(&v);

            /* If last character is a & then value may be continued on next keyword */
            let mut contin = true;
            while contin {
                if len > 0 && v[len - 1] == bb(b'&') {
                    /*  is last char an ampersand?  */
                    valstring[0] = 0;
                    nextcomm[0] = 0;
                    ffgcnt(fptr, &mut valstring, Some(&mut nextcomm), status);

                    if valstring[0] != 0 || nextcomm[0] != 0 {
                        /* If either valstring or nextcom is filled, this must be a CONTINUE line */
                        v[len - 1] = 0; /* erase the trailing & char */

                        if valstring[0] != 0 {
                            len += strlen_safe(&valstring) - 1;
                            // HEAP ALLOCATION
                            v.resize(len + 1, 0); /* increase size */
                            strcat_safe(&mut v, &valstring); /* append the continued chars */
                        }

                        if (commspace > 0) && (nextcomm[0] != 0) {
                            /* If in here, input 'comm' cannot be 0 */
                            /* concantenate comment strings (if any) */

                            let c = c.as_deref_mut().unwrap();

                            if strlen_safe(c) != 0 && addCommDelim {
                                strcat_safe(c, cs!(c" "));
                                commspace -= 1;
                            }

                            strncat_safe(c, &nextcomm, commspace);
                            commspace = FLEN_COMMENT - 1 - strlen_safe(c);
                        }

                        /* Determine if a space delimiter is needed for next
                        comment concatenation (if any).  Assume it is if card length
                        of the most recently read keyword is less than max.
                        keynum is 1-based. */
                        ffghps_safe(fptr, Some(&mut 0), Some(&mut keynum), status);
                        ffgrec_safe(fptr, keynum - 1, Some(&mut card), status);
                        addCommDelim = strlen_safe(&card) < FLEN_CARD - 1;
                    } else {
                        contin = false;
                    }
                } else {
                    contin = false;
                }
            }

            let (p, l, c) = vec_into_raw_parts(v);
            ALLOCATIONS.lock().unwrap().insert(p as usize, (l, c));
            *value = p;
        }
        *status
    }
}

/*--------------------------------------------------------------------------*/
/// Read and return the value of the specified string-valued keyword.
///
/// This new routine was added in 2016 to provide a more convenient user
/// interface than the older ffgkls routine.
///
/// Read a string keyword, returning up to 'naxchars' characters of the value
/// starting with the 'firstchar' character.
/// The input 'value' string must be allocated at least 1 char bigger to
/// allow for the terminating null character.
///
/// This routine may be used to read continued string keywords that use
/// the CONTINUE keyword convention, as well as normal string keywords
/// that are contained within a single header record.
///
/// This routine differs from the ffkls routine in that it does not
/// internally allocate memory for the returned value string, and consequently
/// the calling routine does not need to call fffree to free the memory.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgsky(
    fptr: *mut fitsfile,    /* I - FITS file pointer             */
    keyname: *const c_char, /* I - name of keyword to read       */
    firstchar: c_int,       /* I - first character of string to return */
    maxchar: c_int,         /* I - maximum length of string to return */
    /*    (string will be null terminated)  */
    value: *mut c_char,   /* O - pointer to keyword value      */
    valuelen: *mut c_int, /* O - total length of the keyword value string */
    /*     The returned 'value' string may only */
    /*     contain a piece of the total string, depending */
    /*     on the value of firstchar and maxchar */
    comm: *mut c_char,  /* O - keyword comment (may be NULL) */
    status: *mut c_int, /* IO - error status                 */
) -> c_int {
    unsafe {
        let mut valstring: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
        let mut nextcomm: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];
        let mut card: [c_char; FLEN_CARD] = [0; FLEN_CARD];

        let mut commspace = 0;
        let mut addCommDelim = false;
        let mut keynum = 0;
        let mut len = 0;

        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        raw_to_slice!(keyname);

        // comm is an output, only know size is limited to 69 chars including null
        let mut c: Option<&mut [c_char; FLEN_COMMENT]> = if comm.is_null() {
            None
        } else {
            Some(
                slice::from_raw_parts_mut(comm, FLEN_COMMENT)
                    .try_into()
                    .unwrap(),
            )
        };

        if *status > 0 {
            return *status;
        }

        *value = 0;
        if !valuelen.is_null() {
            *valuelen = 0;
        }

        card[0] = 0;
        if let Some(comm) = c.as_deref_mut() {
            comm[0] = 0;
        }

        ffgcrd_safe(fptr, keyname, &mut card, status);

        if *status > 0 {
            return *status;
        }

        if strlen_safe(&card) < FLEN_CARD - 1 {
            addCommDelim = true;
        }

        ffpsvc_safe(&card, &mut valstring, c.as_deref_mut(), status);

        if *status > 0 {
            return *status;
        }

        if !comm.is_null() {
            /* remaining space in comment string */
            let c = slice::from_raw_parts(comm, FLEN_COMMENT);
            commspace = FLEN_COMMENT - 1 - strlen_safe(c);
        }

        let mut tempstring: Vec<c_char>;

        if valstring[0] == 0 {
            /* null value string? */
            tempstring = vec![0; 1]; /* allocate and return a null string */
        } else {
            /* allocate space,  plus 1 for null */
            tempstring = vec![0; strlen_safe(&valstring) + 1];

            ffc2s(&valstring, tempstring.as_mut_slice(), status); /* convert string to value */
            len = strlen_safe(&tempstring);

            /* If last character is a & then value may be continued on next keyword */
            let mut contin = true;
            while contin && *status <= 0 {
                if len != 0 && tempstring[len - 1] == bb(b'&') {
                    /*  is last char an ampersand?  */
                    valstring[0] = 0;
                    nextcomm[0] = 0;
                    ffgcnt(fptr, &mut valstring, Some(&mut nextcomm), status);

                    if valstring[0] != 0 || nextcomm[0] != 0 {
                        /* If either valstring or nextcom is filled, this must be a CONTINUE line */

                        tempstring[len - 1] = 0; /* erase the trailing & char */

                        if valstring[0] != 0 {
                            len += strlen_safe(&valstring) - 1;
                            tempstring.resize(len + 1, 0); /* increase size */
                            strcat_safe(tempstring.as_mut_slice(), &valstring); /* append the continued chars */
                        }

                        if (commspace > 0) && (nextcomm[0] != 0) {
                            /* If in here, input 'comm' cannot be 0 */
                            /* concantenate comment strings (if any) */

                            let c = c.as_deref_mut().unwrap();

                            if strlen_safe(c) != 0 && addCommDelim {
                                strcat_safe(c, cs!(c" "));
                                commspace -= 1;
                            }
                            strncat_safe(c, &nextcomm, commspace);
                            commspace = FLEN_COMMENT - 1 - strlen_safe(c);
                        }
                        /* Determine if a space delimiter is needed for next
                        comment concatenation (if any).  Assume it is if card length
                        of the most recently read keyword is less than max.
                        keynum is 1-based. */
                        ffghps_safe(fptr, Some(&mut 0), Some(&mut keynum), status);
                        ffgrec_safe(fptr, keynum - 1, Some(&mut card), status);
                        addCommDelim = strlen_safe(&card) < FLEN_CARD - 1;
                    } else {
                        contin = false;
                    }
                } else {
                    contin = false;
                }
            }
        }

        len = strlen_safe(&tempstring);
        if firstchar <= len as c_int {
            let value = slice::from_raw_parts_mut(value, maxchar as usize);

            strncat_safe(
                value,
                &tempstring[(firstchar as usize - 1)..],
                maxchar as usize,
            );
        }

        if !valuelen.is_null() {
            *valuelen = len as c_int; /* total length of the keyword value */
        }

        *status
    }
}

/*--------------------------------------------------------------------------*/
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgskyc(
    fptr: *mut fitsfile,    /* I - FITS file pointer             */
    keyname: *const c_char, /* I - name of keyword to read       */
    firstchar: c_int,       /* I - first character of string to return */
    maxchar: c_int,         /* I - maximum length of string to return */
    /*    (string will be null terminated)  */
    maxcomchar: c_int,    /* I - maximum length of comment to return */
    value: *mut c_char,   /* O - pointer to keyword value      */
    valuelen: *mut c_int, /* O - total length of the keyword value string */
    /*     The returned 'value' string may only */
    /*     contain a piece of the total string, depending */
    /*     on the value of firstchar and maxchar */
    comm: *mut c_char,  /* O - keyword comment (may be NULL) */
    comlen: *mut c_int, /* O - total length of the comment string */
    status: *mut c_int, /* IO - error status                 */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let valuelen = valuelen.as_mut().expect(NULL_MSG);
        let comlen = comlen.as_mut().expect(NULL_MSG);

        raw_to_slice!(keyname);

        let comm = match comm.is_null() {
            true => None,
            false => Some(slice::from_raw_parts_mut(comm, 1 + maxcomchar as usize)),
        };

        let value = match value.is_null() {
            true => None,
            false => Some(slice::from_raw_parts_mut(value, 1 + maxchar as usize)),
        };

        if *status > 0 {
            return *status;
        }

        ffglkut(
            fptr, keyname, firstchar, maxchar, maxcomchar, value, valuelen, comm, comlen, status,
        );

        *status
    }
}

/*--------------------------------------------------------------------------*/
pub(crate) fn ffglkut(
    fptr: &mut fitsfile, /* I - FITS file pointer             */
    keyname: &[c_char],  /* I - name of keyword to read       */
    firstchar: c_int,    /* I - first character of string to return */
    maxvalchar: c_int,   /* I - maximum length of string to return */
    /*    (string will be null terminated)  */
    maxcomchar: c_int,                /* I - maximum length of comment to return */
    mut value: Option<&mut [c_char]>, /* O - pointer to keyword value (may be NULL) */
    valuelen: &mut c_int,             /* O - total length of the keyword value string */
    /*     The returned 'value' string may only */
    /*     contain a piece of the total string, depending */
    /*     on the value of firstchar and maxvalchar */
    mut comm: Option<&mut [c_char]>, /* O - keyword comment (may be NULL) */
    comlen: &mut c_int,              /* O - total length of comment string */
    status: &mut c_int,              /* IO - error status                 */
) -> c_int {
    let mut valstring: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut comstring: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];
    let mut card: [c_char; FLEN_CARD] = [0; FLEN_CARD];
    let mut dynValStr: Vec<c_char> = Vec::new();
    let mut dynComStr: Vec<c_char> = Vec::new();
    let mut contin = false;
    let mut addCommDelim = false;
    let mut keynum = 0;
    let mut lenOnly = false;
    let mut savedKeyPos = 1;
    let mut len = 0;
    let mut lenc = 0;

    if *status > 0 {
        return *status;
    }

    if maxvalchar == 0 && maxcomchar == 0 {
        lenOnly = true;
    }

    if let Some(value) = value.as_deref_mut() {
        value[0] = 0;
    }
    if let Some(comm) = comm.as_deref_mut() {
        comm[0] = 0;
    }
    /* If lenOnly, 'value' and 'comm' should not be accessed after this point.*/

    *valuelen = 0;
    *comlen = 0;
    card[0] = 0;
    valstring[0] = 0;
    comstring[0] = 0;

    ffgcrd_safe(fptr, keyname, &mut card, status);
    if *status > 0 {
        return *status;
    }
    ffpsvc_safe(&card, &mut valstring, Some(&mut comstring), status);
    if *status > 0 {
        return *status;
    }
    if strlen_safe(&card) < FLEN_CARD - 1 && comstring[0] != 0 {
        addCommDelim = true;
    }

    /* If called in lenOnly mode, there's a good chance the user will soon call
    this again to read the value string.  Therefore we'll save and later restore
    the original keyword position. */
    if lenOnly {
        ffghps_safe(fptr, Some(&mut 0), Some(&mut savedKeyPos), status);
        if savedKeyPos > 1 {
            savedKeyPos -= 1;
        }
    }

    if valstring[0] == 0
    /* null value string? */
    {
        dynValStr.push(0); /* allocate and return a null string */
        dynComStr.push(0);
    } else {
        /* allocate space,  plus 1 for null */
        dynValStr.resize(strlen_safe(&valstring) + 1, 0);

        ffc2s(&valstring, &mut dynValStr, status); /* convert string to value */
        len = strlen_safe(&dynValStr);

        dynComStr.resize(strlen_safe(&comstring) + 1, 0);
        dynComStr[0] = 0;
        strcpy_safe(&mut dynComStr, &comstring);
        lenc = strlen_safe(&dynComStr);

        /* If last character is a & then value may be continued on next keyword */
        contin = true;
        while contin && *status <= 0 {
            if len != 0 && dynValStr[len - 1] == bb(b'&')
            /*  is last char an ampersand?  */
            {
                valstring[0] = 0;
                comstring[0] = 0;
                ffgcnt(fptr, &mut valstring, Some(&mut comstring), status);
                if valstring[0] != 0 || comstring[0] != 0
                /* If either valstring or comstring
                is filled, this must be a CONTINUE line */
                {
                    dynValStr[len - 1] = 0; /* erase the trailing & char */
                    len -= 1;
                    if valstring[0] != 0 {
                        len += strlen_safe(&valstring);
                        dynValStr.resize(len + 1, 0); /* increase size */
                        strcat_safe(&mut dynValStr, &valstring); /* append the continued chars */
                    }
                    if comstring[0] != 0 {
                        /* concantenate comment strings */
                        if addCommDelim {
                            lenc += strlen_safe(&comstring) + 1;
                            dynComStr.resize(lenc + 1, 0);
                            strcat_safe(&mut dynComStr, cs!(c" "));
                            strcat_safe(&mut dynComStr, &comstring);
                        } else {
                            lenc += strlen_safe(&comstring);
                            dynComStr.resize(lenc + 1, 0);
                            strcat_safe(&mut dynComStr, &comstring);
                        }
                    }
                    /* Determine if a space delimiter is needed for next
                    comment concatenation (if any).  Assume it is if card length
                    of the most recently read keyword is less than max.
                    keynum is 1-based. */
                    ffghps_safe(fptr, Some(&mut 0), Some(&mut keynum), status);
                    ffgrec_safe(fptr, keynum - 1, Some(&mut card), status);
                    addCommDelim = (strlen_safe(&card) < FLEN_CARD - 1) && lenc != 0;
                } else {
                    contin = false;
                }
            } else {
                contin = false;
            }
        }
    }

    /* Resetting len and lenc really shouldn't be necessary here, but
    just to make sure ... */
    len = strlen_safe(&dynValStr);
    lenc = strlen_safe(&dynComStr);
    *valuelen = len as c_int;
    *comlen = lenc as c_int;
    if lenOnly {
        ffmaky_safe(fptr, savedKeyPos, status);
    } else {
        if let Some(value) = value {
            if firstchar <= len as c_int {
                strncat_safe(
                    value,
                    &dynValStr[(firstchar as usize - 1)..],
                    maxvalchar as usize,
                );
            }
        }
        if let Some(comm) = comm {
            strncat_safe(comm, &dynComStr, maxcomchar as usize);
        }
    }

    *status
}

/*--------------------------------------------------------------------------*/
/// Free the memory that was previously allocated by CFITSIO,
/// such as by ffgkls or fits_hdr2str
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn fffree(
    value: *mut c_void, /* I - pointer to keyword value  */
    status: *mut c_int, /* IO - error status             */
) -> c_int {
    unsafe {
        if *status > 0 {
            return *status;
        }

        if !value.is_null() {
            // HEAP DEALLOCATION

            let mut alloc_lock = ALLOCATIONS.lock().unwrap();
            let alloc = alloc_lock.remove(&(value as usize));
            if let Some((l, c)) = alloc {
                // HEAP DEALLOCATION
                let _ = Vec::from_raw_parts(value, l, c);
            } else {
                let _ = Vec::from_raw_parts(value, 1, 1);
            }
        }
        *status
    }
}

/*--------------------------------------------------------------------------*/
/// Attempt to read the next keyword, returning the string value
/// if it is a continuation of the previous string keyword value.
/// This uses the HEASARC convention for continuing long string values
/// over multiple keywords.  Each continued string is terminated with a
/// backslash character, and the continuation follows on the next keyword
/// which must have the name CONTINUE without an equal sign in column 9
/// of the card.  If the next card is not a continuation, then the returned
/// value string will be null.
pub(crate) fn ffgcnt(
    fptr: &mut fitsfile,                       /* I - FITS file pointer         */
    value: &mut [c_char],                      /* O - continued string value    */
    comm: Option<&mut [c_char; FLEN_COMMENT]>, /* O - continued comment string  */
    status: &mut c_int,                        /* IO - error status             */
) -> c_int {
    let mut tstatus = 0;
    let mut card: [c_char; FLEN_CARD] = [0; FLEN_CARD];
    let mut strval: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];

    if *status > 0 {
        return *status;
    }

    tstatus = 0;
    value[0] = 0;

    if ffgnky(fptr, &mut card, &mut tstatus) > 0 {
        /*  read next keyword  */
        return *status; /*  hit end of header  */
    }

    if strncmp_safe(&card, cs!(c"CONTINUE  "), 10) == 0 {
        /* a continuation card? */

        strncpy_safe(&mut card, cs!(c"D2345678=  "), 10); /* overwrite a dummy keyword name */
        ffpsvc_safe(&card, &mut strval, comm, &mut tstatus); /*  get the string value & comment */

        ffc2s(&strval, value, &mut tstatus); /* remove the surrounding quotes */

        if tstatus > 0 {
            /*  return null if error status was returned  */
            value[0] = 0;
        }
    } else {
        ffmrky_safe(fptr, -1, status); /* reset the keyword pointer */
    }

    *status
}

/*--------------------------------------------------------------------------*/
/// Read (get) the named keyword, returning the value and comment.
/// The returned value = 1 if the keyword is true, else = 0 if false.
/// The comment may be up to 69 characters long.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgkyl(
    fptr: *mut fitsfile,    /* I - FITS file pointer         */
    keyname: *const c_char, /* I - name of keyword to read   */
    value: *mut c_int,      /* O - keyword value             */
    comm: *mut c_char,      /* O - keyword comment           */
    status: *mut c_int,     /* IO - error status             */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let value = value.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        raw_to_slice!(keyname);

        // comm is an output, only know size is limited to 69 chars including null
        let comm: Option<&mut [c_char; FLEN_COMMENT]> = if comm.is_null() {
            None
        } else {
            Some(
                slice::from_raw_parts_mut(comm, FLEN_COMMENT)
                    .try_into()
                    .unwrap(),
            )
        };

        ffgkyl_safe(fptr, keyname, value, comm, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Read (get) the named keyword, returning the value and comment.
/// The returned value = 1 if the keyword is true, else = 0 if false.
/// The comment may be up to 69 characters long.
pub fn ffgkyl_safe(
    fptr: &mut fitsfile,                       /* I - FITS file pointer         */
    keyname: &[c_char],                        /* I - name of keyword to read   */
    value: &mut c_int,                         /* O - keyword value             */
    comm: Option<&mut [c_char; FLEN_COMMENT]>, /* O - keyword comment           */
    status: &mut c_int,                        /* IO - error status             */
) -> c_int {
    let mut valstring: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];

    if *status > 0 {
        return *status;
    }

    ffgkey_safe(fptr, keyname, &mut valstring, comm, status); /* read the keyword */
    ffc2l(&valstring, value, status); /* convert string to value */

    *status
}

/*--------------------------------------------------------------------------*/
///Read (get) the named keyword, returning the value and comment.
///
///The value will be implicitly converted to a (long) integer if it not
///already of this datatype.  The comment may be up to 69 characters long.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgkyj(
    fptr: *mut fitsfile,    /* I - FITS file pointer         */
    keyname: *const c_char, /* I - name of keyword to read   */
    value: *mut c_long,     /* O - keyword value             */
    comm: *mut c_char,      /* O - keyword comment           */
    status: *mut c_int,     /* IO - error status             */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let value = value.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        raw_to_slice!(keyname);

        // comm is an output, only know size is limited to 69 chars including null
        let comm: Option<&mut [c_char; FLEN_COMMENT]> = if comm.is_null() {
            None
        } else {
            Some(
                slice::from_raw_parts_mut(comm, FLEN_COMMENT)
                    .try_into()
                    .unwrap(),
            )
        };

        ffgkyj_safe(fptr, keyname, value, comm, status)
    }
}

/*--------------------------------------------------------------------------*/
///Read (get) the named keyword, returning the value and comment.
///
///The value will be implicitly converted to a (long) integer if it not
///already of this datatype.  The comment may be up to 69 characters long.
pub fn ffgkyj_safe(
    fptr: &mut fitsfile,                       /* I - FITS file pointer         */
    keyname: &[c_char],                        /* I - name of keyword to read   */
    value: &mut c_long,                        /* O - keyword value             */
    comm: Option<&mut [c_char; FLEN_COMMENT]>, /* O - keyword comment           */
    status: &mut c_int,                        /* IO - error status             */
) -> c_int {
    let mut valstring: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];

    if *status > 0 {
        return *status;
    }

    ffgkey_safe(fptr, keyname, &mut valstring, comm, status); /* read the keyword */
    ffc2i(&valstring, value, status); /* convert string to value */

    *status
}

/*--------------------------------------------------------------------------*/
/// Read (get) the named keyword, returning the value and comment.
///
/// The value will be implicitly converted to a (LONGLONG) integer if it not
/// already of this datatype.  The comment may be up to 69 characters long.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgkyjj(
    fptr: *mut fitsfile,    /* I - FITS file pointer         */
    keyname: *const c_char, /* I - name of keyword to read   */
    value: *mut LONGLONG,   /* O - keyword value             */
    comm: *mut c_char,      /* O - keyword comment           */
    status: *mut c_int,     /* IO - error status             */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let value = value.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        raw_to_slice!(keyname);

        // comm is an output, only know size is limited to 69 chars including null
        let comm: Option<&mut [c_char; FLEN_COMMENT]> = if comm.is_null() {
            None
        } else {
            Some(
                slice::from_raw_parts_mut(comm, FLEN_COMMENT)
                    .try_into()
                    .unwrap(),
            )
        };

        ffgkyjj_safe(fptr, keyname, value, comm, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Read (get) the named keyword, returning the value and comment.
///
/// The value will be implicitly converted to a (LONGLONG) integer if it not
/// already of this datatype.  The comment may be up to 69 characters long.
pub fn ffgkyjj_safe(
    fptr: &mut fitsfile,                       /* I - FITS file pointer         */
    keyname: &[c_char],                        /* I - name of keyword to read   */
    value: &mut LONGLONG,                      /* O - keyword value             */
    comm: Option<&mut [c_char; FLEN_COMMENT]>, /* O - keyword comment           */
    status: &mut c_int,                        /* IO - error status             */
) -> c_int {
    let mut valstring: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];

    if *status > 0 {
        return *status;
    }

    ffgkey_safe(fptr, keyname, &mut valstring, comm, status); /* read the keyword */
    ffc2j(&valstring, value, status); /* convert string to value */

    *status
}

/*--------------------------------------------------------------------------*/
/// Read (get) the named keyword, returning the value and comment.
///
/// The value will be implicitly converted to a (ULONGLONG) integer if it not
/// already of this datatype.  The comment may be up to 69 characters long.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgkyujj(
    fptr: *mut fitsfile,
    keyname: *const c_char,
    value: *mut u64,
    comm: *mut c_char,
    status: *mut c_int,
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let value = value.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        raw_to_slice!(keyname);

        // comm is an output, only know size is limited to 69 chars including null
        let comm: Option<&mut [c_char; FLEN_COMMENT]> = if comm.is_null() {
            None
        } else {
            Some(
                slice::from_raw_parts_mut(comm, FLEN_COMMENT)
                    .try_into()
                    .unwrap(),
            )
        };

        ffgkyujj_safe(fptr, keyname, value, comm, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Read (get) the named keyword, returning the value and comment.
///
/// The value will be implicitly converted to a (ULONGLONG) integer if it not
/// already of this datatype.  The comment may be up to 69 characters long.
pub fn ffgkyujj_safe(
    fptr: &mut fitsfile,
    keyname: &[c_char],
    value: &mut u64,
    comm: Option<&mut [c_char; FLEN_COMMENT]>,
    status: &mut c_int,
) -> c_int {
    let mut valstring: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];

    if *status > 0 {
        return *status;
    }

    ffgkey_safe(fptr, keyname, &mut valstring, comm, status); /* read the keyword */
    ffc2uj(&valstring, value, status); /* convert string to value */

    *status
}

/*--------------------------------------------------------------------------*/
/// Read (get) the named keyword, returning the value and comment.
/// The value will be implicitly converted to a float if it not
/// already of this datatype.  The comment may be up to 69 characters long.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgkye(
    fptr: *mut fitsfile,    /* I - FITS file pointer         */
    keyname: *const c_char, /* I - name of keyword to read   */
    value: *mut f32,        /* O - keyword value             */
    comm: *mut c_char,      /* O - keyword comment           */
    status: *mut c_int,     /* IO - error status             */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let value = value.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        raw_to_slice!(keyname);

        // comm is an output, only know size is limited to 69 chars including null
        let comm: Option<&mut [c_char; FLEN_COMMENT]> = if comm.is_null() {
            None
        } else {
            Some(
                slice::from_raw_parts_mut(comm, FLEN_COMMENT)
                    .try_into()
                    .unwrap(),
            )
        };

        ffgkye_safe(fptr, keyname, value, comm, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Read (get) the named keyword, returning the value and comment.
/// The value will be implicitly converted to a float if it not
/// already of this datatype.  The comment may be up to 69 characters long.
pub fn ffgkye_safe(
    fptr: &mut fitsfile,                       /* I - FITS file pointer         */
    keyname: &[c_char],                        /* I - name of keyword to read   */
    value: &mut f32,                           /* O - keyword value             */
    comm: Option<&mut [c_char; FLEN_COMMENT]>, /* O - keyword comment           */
    status: &mut c_int,                        /* IO - error status             */
) -> c_int {
    let mut valstring: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];

    if *status > 0 {
        return *status;
    }

    ffgkey_safe(fptr, keyname, &mut valstring, comm, status); /* read the keyword */
    ffc2r(&valstring, value, status); /* convert string to value */

    *status
}

/*--------------------------------------------------------------------------*/
/// Read (get) the named keyword, returning the value and comment.
/// The value will be implicitly converted to a double if it not
/// already of this datatype.  The comment may be up to 69 characters long.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgkyd(
    fptr: *mut fitsfile,    /* I - FITS file pointer         */
    keyname: *const c_char, /* I - name of keyword to read   */
    value: *mut f64,        /* O - keyword value             */
    comm: *mut c_char,      /* O - keyword comment           */
    status: *mut c_int,     /* IO - error status             */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let value = value.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        raw_to_slice!(keyname);

        // comm is an output, only know size is limited to 69 chars including null
        let comm: Option<&mut [c_char; FLEN_COMMENT]> = if comm.is_null() {
            None
        } else {
            Some(
                slice::from_raw_parts_mut(comm, FLEN_COMMENT)
                    .try_into()
                    .unwrap(),
            )
        };

        ffgkyd_safe(fptr, keyname, value, comm, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Read (get) the named keyword, returning the value and comment.
/// The value will be implicitly converted to a double if it not
/// already of this datatype.  The comment may be up to 69 characters long.
pub fn ffgkyd_safe(
    fptr: &mut fitsfile,                       /* I - FITS file pointer         */
    keyname: &[c_char],                        /* I - name of keyword to read   */
    value: &mut f64,                           /* O - keyword value             */
    comm: Option<&mut [c_char; FLEN_COMMENT]>, /* O - keyword comment           */
    status: &mut c_int,                        /* IO - error status             */
) -> c_int {
    let mut valstring: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];

    if *status > 0 {
        return *status;
    }

    ffgkey_safe(fptr, keyname, &mut valstring, comm, status); /* read the keyword */
    ffc2d(&valstring, value, status); /* convert string to value */

    *status
}

/*--------------------------------------------------------------------------*/
/// Read (get) the named keyword, returning the value and comment.
/// The keyword must have a complex value. No implicit data conversion
/// will be performed.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgkyc(
    fptr: *mut fitsfile,    /* I - FITS file pointer         */
    keyname: *const c_char, /* I - name of keyword to read   */
    value: *mut [f32; 2],   /* O - keyword value (real,imag) */
    comm: *mut c_char,      /* O - keyword comment           */
    status: *mut c_int,     /* IO - error status             */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let value = value.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        raw_to_slice!(keyname);

        // comm is an output, only know size is limited to 69 chars including null
        let comm: Option<&mut [c_char; FLEN_COMMENT]> = if comm.is_null() {
            None
        } else {
            Some(
                slice::from_raw_parts_mut(comm, FLEN_COMMENT)
                    .try_into()
                    .unwrap(),
            )
        };

        ffgkyc_safe(fptr, keyname, value, comm, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Read (get) the named keyword, returning the value and comment.
/// The keyword must have a complex value. No implicit data conversion
/// will be performed.
pub fn ffgkyc_safe(
    fptr: &mut fitsfile,                       /* I - FITS file pointer         */
    keyname: &[c_char],                        /* I - name of keyword to read   */
    value: &mut [f32; 2],                      /* O - keyword value (real,imag) */
    comm: Option<&mut [c_char; FLEN_COMMENT]>, /* O - keyword comment           */
    status: &mut c_int,                        /* IO - error status             */
) -> c_int {
    let mut valstring: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut message: [c_char; FLEN_ERRMSG] = [0; FLEN_ERRMSG];
    let mut len = 0;

    if *status > 0 {
        return *status;
    }

    ffgkey_safe(fptr, keyname, &mut valstring, comm, status); /* read the keyword */

    if valstring[0] != bb(b'(') {
        /* test that this is a complex keyword */

        int_snprintf!(
            &mut message,
            FLEN_ERRMSG,
            "keyword {} does not have a complex value (ffgkyc):",
            slice_to_str!(&keyname),
        );
        ffpmsg_slice(&message);
        ffpmsg_slice(&valstring);
        *status = BAD_C2F;
        return *status;
    }

    valstring[0] = bb(b' '); /* delete the opening parenthesis */
    len = strcspn_safe(&valstring, cs!(c")"));
    valstring[len] = 0; /* delete the closing parenthesis */

    len = strcspn_safe(&valstring, cs!(c","));
    valstring[len] = 0;

    ffc2r(&valstring, &mut value[0], status); /* convert the real part */
    ffc2r(&valstring[(len + 1)..], &mut value[1], status); /* convert imag. part */
    *status
}

/*--------------------------------------------------------------------------*/
/// Read (get) the named keyword, returning the value and comment.
/// The keyword must have a complex value. No implicit data conversion
/// will be performed.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgkym(
    fptr: *mut fitsfile,    /* I - FITS file pointer         */
    keyname: *const c_char, /* I - name of keyword to read   */
    value: *mut [f64; 2],   /* O - keyword value (real,imag) */
    comm: *mut c_char,      /* O - keyword comment           */
    status: *mut c_int,     /* IO - error status             */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let value = value.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        raw_to_slice!(keyname);

        // comm is an output, only know size is limited to 69 chars including null
        let comm: Option<&mut [c_char; FLEN_COMMENT]> = if comm.is_null() {
            None
        } else {
            Some(
                slice::from_raw_parts_mut(comm, FLEN_COMMENT)
                    .try_into()
                    .unwrap(),
            )
        };

        ffgkym_safe(fptr, keyname, value, comm, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Read (get) the named keyword, returning the value and comment.
/// The keyword must have a complex value. No implicit data conversion
/// will be performed.
pub fn ffgkym_safe(
    fptr: &mut fitsfile,                       /* I - FITS file pointer         */
    keyname: &[c_char],                        /* I - name of keyword to read   */
    value: &mut [f64; 2],                      /* O - keyword value (real,imag) */
    comm: Option<&mut [c_char; FLEN_COMMENT]>, /* O - keyword comment           */
    status: &mut c_int,                        /* IO - error status             */
) -> c_int {
    let mut valstring: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut message: [c_char; FLEN_ERRMSG] = [0; FLEN_ERRMSG];
    let mut len = 0;

    if *status > 0 {
        return *status;
    }

    ffgkey_safe(fptr, keyname, &mut valstring, comm, status); /* read the keyword */

    if valstring[0] != bb(b'(') {
        /* test that this is a complex keyword */

        int_snprintf!(
            &mut message,
            FLEN_ERRMSG,
            "keyword {} does not have a complex value (ffgkym):",
            slice_to_str!(&keyname),
        );
        ffpmsg_slice(&message);
        ffpmsg_slice(&valstring);
        *status = BAD_C2D;
        return *status;
    }

    valstring[0] = bb(b' '); /* delete the opening parenthesis */
    len = strcspn_safe(&valstring, cs!(c")"));
    valstring[len] = 0; /* delete the closing parenthesis */

    len = strcspn_safe(&valstring, cs!(c","));
    valstring[len] = 0;

    ffc2d(&valstring, &mut value[0], status); /* convert the real part */
    ffc2d(&valstring[(len + 1)..], &mut value[1], status); /* convert the imag. part */

    *status
}

/*--------------------------------------------------------------------------*/
/// Read (get) the named keyword, returning the value and comment.
///
/// The integer and fractional parts of the value are returned in separate
/// variables, to allow more numerical precision to be passed.  This
/// effectively passes a 'triple' precision value, with a 4-byte integer
/// and an 8-byte fraction.  The comment may be up to 69 characters long.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgkyt(
    fptr: *mut fitsfile,    /* I - FITS file pointer         */
    keyname: *const c_char, /* I - name of keyword to read   */
    ivalue: *mut c_long,    /* O - integer part of keyword value     */
    fraction: *mut f64,     /* O - fractional part of keyword value  */
    comm: *mut c_char,      /* O - keyword comment           */
    status: *mut c_int,     /* IO - error status             */
) -> c_int {
    unsafe {
        let mut valstring: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];

        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let ivalue = ivalue.as_mut().expect(NULL_MSG);
        let fraction = fraction.as_mut().expect(NULL_MSG);

        raw_to_slice!(keyname);

        // comm is an output, only know size is limited to 69 chars including null
        let c: Option<&mut [c_char; FLEN_COMMENT]> = if comm.is_null() {
            None
        } else {
            Some(
                slice::from_raw_parts_mut(comm, FLEN_COMMENT)
                    .try_into()
                    .unwrap(),
            )
        };

        if *status > 0 {
            return *status;
        }

        ffgkey_safe(fptr, keyname, &mut valstring, c, status); /* read the keyword */

        /*  read the entire value string as a double, to get the integer part */
        ffc2d(&valstring, fraction, status);

        *ivalue = *fraction as c_long;

        *fraction -= *ivalue as f64;

        /* see if we need to read the fractional part again with more precision */
        /* look for decimal point, without an exponential E or D character */

        let loc = strchr_safe(&valstring, bb(b'.'));
        if let Some(loc) = loc {
            if strchr_safe(&valstring, bb(b'E')).is_none()
                && strchr_safe(&valstring, bb(b'D')).is_none()
            {
                ffc2d(&valstring[loc..], fraction, status);
            }
        }

        *status
    }
}

/*--------------------------------------------------------------------------*/
/// Read (get) the nkey-th keyword returning the keyword name, value and comment.
///
/// The value is just the literal string of characters in the value field
/// of the keyword.  In the case of a string valued keyword, the returned
/// value includes the leading and closing quote characters.  The value may be
/// up to 70 characters long, and the comment may be up to 72 characters long.
/// If the keyword has no value (no equal sign in column 9) then a null value
/// is returned.  If comm = NULL, then do not return the comment string.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgkyn(
    fptr: *mut fitsfile,  /* I - FITS file pointer             */
    nkey: c_int,          /* I - number of the keyword to read */
    keyname: *mut c_char, /* O - name of the keyword           */
    value: *mut c_char,   /* O - keyword value                 */
    comm: *mut c_char,    /* O - keyword comment               */
    status: *mut c_int,   /* IO - error status                 */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        let keyname = slice::from_raw_parts_mut(keyname, FLEN_KEYWORD)
            .try_into()
            .unwrap();
        let value = slice::from_raw_parts_mut(value, FLEN_VALUE)
            .try_into()
            .unwrap();

        let comm: Option<&mut [c_char; FLEN_COMMENT]> = match comm.is_null() {
            true => None,
            false => Some(
                slice::from_raw_parts_mut(comm, FLEN_COMMENT)
                    .try_into()
                    .unwrap(),
            ),
        };

        ffgkyn_safe(fptr, nkey, keyname, value, comm, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Read (get) the nkey-th keyword returning the keyword name, value and comment.
///
/// The value is just the literal string of characters in the value field
/// of the keyword.  In the case of a string valued keyword, the returned
/// value includes the leading and closing quote characters.  The value may be
/// up to 70 characters long, and the comment may be up to 72 characters long.
/// If the keyword has no value (no equal sign in column 9) then a null value
/// is returned.  If comm = NULL, then do not return the comment string.
pub fn ffgkyn_safe(
    fptr: &mut fitsfile,                       /* I - FITS file pointer             */
    nkey: c_int,                               /* I - number of the keyword to read */
    keyname: &mut [c_char; FLEN_KEYWORD],      /* O - name of the keyword           */
    value: &mut [c_char; FLEN_VALUE],          /* O - keyword value                 */
    comm: Option<&mut [c_char; FLEN_COMMENT]>, /* O - keyword comment               */
    status: &mut c_int,                        /* IO - error status                 */
) -> c_int {
    let mut card: [c_char; FLEN_CARD] = [0; FLEN_CARD];
    let mut sbuff: [c_char; FLEN_CARD] = [0; FLEN_CARD];
    let mut namelen = 0;

    if *status > 0 {
        return *status;
    }

    if ffgrec_safe(fptr, nkey, Some(&mut card), status) > 0 {
        /* get the 80-byte card */
        return *status;
    }

    ffgknm_safe(&card, keyname, &mut namelen, status); /* get the keyword name */

    //println!("keyword name len={}", namelen);
    //std::io::stdout().write_all(cast_slice(&keyname));

    if ffpsvc_safe(&card, value, comm, status) > 0 {
        /* parse value and comment */

        return *status;
    }

    if fftrec_safe(keyname, status) > 0 {
        /* test keyword name; catches no END */

        int_snprintf!(
            &mut sbuff,
            FLEN_CARD,
            "Name of keyword no. {} contains illegal character(s): {}",
            nkey,
            slice_to_str!(keyname),
        );
        ffpmsg_slice(&sbuff);

        if nkey % 36 == 0 {
            /* test if at beginning of 36-card FITS record */
            ffpmsg_str("  (This may indicate a missing END keyword).");
        };
    }
    *status
}

/*--------------------------------------------------------------------------*/
/// Read (get) an indexed array of keywords with index numbers between
/// NSTART and (NSTART + NMAX -1) inclusive.  
/// This routine does NOT support the HEASARC long string convention.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgkns(
    fptr: *mut fitsfile,     /* I - FITS file pointer                    */
    keyname: *const c_char,  /* I - root name of keywords to read        */
    nstart: c_int,           /* I - starting index number                */
    nmax: c_int,             /* I - maximum number of keywords to return */
    value: *mut *mut c_char, /* O - array of pointers to keyword values  */
    nfound: *mut c_int,      /* O - number of values that were returned  */
    status: *mut c_int,      /* IO - error status                        */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let nfound = nfound.as_mut().expect(NULL_MSG);
        raw_to_slice!(keyname);

        let value = slice::from_raw_parts_mut(value, nmax as usize);
        let mut v_value = Vec::new();

        for item in value {
            let value_item = slice::from_raw_parts_mut(*item, FLEN_VALUE);
            v_value.push(value_item);
        }

        ffgkns_safe(fptr, keyname, nstart, nmax, &mut v_value, nfound, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Read (get) an indexed array of keywords with index numbers between
/// NSTART and (NSTART + NMAX -1) inclusive.  
/// This routine does NOT support the HEASARC long string convention.
pub fn ffgkns_safe(
    fptr: &mut fitsfile,         /* I - FITS file pointer                    */
    keyname: &[c_char],          /* I - root name of keywords to read        */
    nstart: c_int,               /* I - starting index number                */
    nmax: c_int,                 /* I - maximum number of keywords to return */
    value: &mut [&mut [c_char]], /* O - array of pointers to keyword values  */
    nfound: &mut c_int,          /* O - number of values that were returned  */
    status: &mut c_int,          /* IO - error status                        */
) -> c_int {
    let mut nend: c_int = 0;
    let mut lenroot: usize = 0;
    let ii: c_int = 0;
    let mut nkeys: c_int = 0;
    let mut mkeys: c_int = 0;
    let mut tstatus: c_int = 0;
    let mut undefinedval: bool;
    let mut ival: c_long = 0;
    let mut keyroot: [c_char; FLEN_KEYWORD] = [0; FLEN_KEYWORD];
    let mut keyindex: [c_char; 8] = [0; 8];
    let mut card: [c_char; FLEN_CARD] = [0; FLEN_CARD];
    let mut svalue: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut comm: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];

    if *status > 0 {
        return *status;
    }

    *nfound = 0;
    nend = nstart + nmax - 1;

    keyroot[0] = 0;

    strncat_safe(&mut keyroot, keyname, FLEN_KEYWORD - 1);

    lenroot = strlen_safe(&keyroot);

    if lenroot == 0 {
        /*  root must be at least 1 char long  */
        return *status;
    }

    for ii in 0..(lenroot) {
        /*  make sure upper case  */
        keyroot[ii] = toupper(keyroot[ii]);
    }

    ffghps_safe(fptr, Some(&mut nkeys), Some(&mut mkeys), status); /*  get the number of keywords  */

    undefinedval = false;
    for ii in 3..(nkeys as usize) {
        if ffgrec_safe(fptr, ii as c_int, Some(&mut card), status) > 0 {
            /*  get next keyword  */
            return *status;
        }
        if strncmp_safe(&keyroot, &card, lenroot) == 0 {
            /* see if keyword matches */
            keyindex[0] = 0;
            let equalssign = strchr_safe(&card, bb(b'='));
            if equalssign.is_none() {
                continue; /* keyword has no value */
            }
            let equalssign = equalssign.unwrap();

            if equalssign - lenroot > 7 {
                *status = BAD_KEYCHAR;
                return *status;
            }

            /*  copy suffix  */
            strncat_safe(&mut keyindex, &card[lenroot..], equalssign - lenroot);

            tstatus = 0;
            if ffc2ii(&keyindex, &mut ival, &mut tstatus) <= 0 {
                /*  test suffix  */
                if ival <= nend as c_long && ival >= nstart as c_long {
                    let v = &mut value[(ival - nstart as c_long) as usize];

                    ffpsvc_safe(&card, &mut svalue, Some(&mut comm), status); /*  parse the value */

                    ffc2s(&svalue, v, status); /* convert */
                    if ival - nstart as c_long + 1 > (*nfound) as c_long {
                        *nfound = (ival - nstart as c_long + 1) as c_int; /*  max found */
                    }

                    if *status == VALUE_UNDEFINED {
                        undefinedval = true;
                        *status = 0; /* reset status to read remaining values */
                    }
                }
            }
        }
    }
    if undefinedval && (*status <= 0) {
        *status = VALUE_UNDEFINED; /* report at least 1 value undefined */
    }
    *status
}

/*--------------------------------------------------------------------------*/
/// Read (get) an indexed array of keywords with index numbers between
/// NSTART and (NSTART + NMAX -1) inclusive.  
/// The returned value = 1 if the keyword is true, else = 0 if false.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgknl(
    fptr: *mut fitsfile,    /* I - FITS file pointer                    */
    keyname: *const c_char, /* I - root name of keywords to read        */
    nstart: c_int,          /* I - starting index number                */
    nmax: c_int,            /* I - maximum number of keywords to return */
    value: *mut c_int,      /* O - array of keyword values              */
    nfound: *mut c_int,     /* O - number of values that were returned  */
    status: *mut c_int,     /* IO - error status                        */
) -> c_int {
    unsafe {
        //int nend, lenroot, ii, nkeys, mkeys, tstatus, undefinedval;
        //long ival;
        let mut nkeys: c_int = 0;
        let mut mkeys: c_int = 0;
        let mut ival: c_long = 0;

        let mut keyroot: [c_char; FLEN_KEYWORD] = [0; FLEN_KEYWORD];
        let mut keyindex: [c_char; 8] = [0; 8];
        let mut card: [c_char; FLEN_CARD] = [0; FLEN_CARD];
        let mut svalue: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
        let mut comm: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];

        let nmax = nmax as usize;
        let nstart = nstart as usize;

        // *equalssign;

        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);
        raw_to_slice!(keyname);

        let value = slice::from_raw_parts_mut(value, nmax);

        if *status > 0 {
            return *status;
        }

        *nfound = 0;
        let nend = nstart + nmax - 1;

        keyroot[0] = 0;
        strncat_safe(&mut keyroot, keyname, FLEN_KEYWORD - 1);

        let lenroot = strlen_safe(&keyroot);

        if lenroot == 0 {
            /*  root must be at least 1 char long  */
            return *status;
        }

        for ii in 0..(lenroot) {
            /*  make sure upper case  */
            keyroot[ii] = toupper(keyroot[ii]);
        }

        ffghps_safe(fptr, Some(&mut nkeys), Some(&mut mkeys), status); /*  get the number of keywords  */

        ffmaky_safe(fptr, 3, status); /* move to 3rd keyword (skip 1st 2 keywords) */

        let mut undefinedval = false;
        for ii in 3..=(nkeys as usize) {
            if ffgnky(fptr, &mut card, status) > 0 {
                /*  get next keyword  */
                return *status;
            }

            if strncmp_safe(&keyroot, &card, lenroot) == 0
            /* see if keyword matches */
            {
                keyindex[0] = 0;
                let equalssign = strchr_safe(&card, bb(b'='));

                if equalssign.is_none() {
                    continue; /* keyword has no value */
                }
                let equalssign = equalssign.unwrap();
                if equalssign - lenroot > 7 {
                    *status = BAD_KEYCHAR;
                    return *status;
                }
                strncat_safe(&mut keyindex, &card[lenroot..], equalssign - lenroot); /*  copy suffix  */

                let mut tstatus = 0;
                if ffc2ii(&keyindex, &mut ival, &mut tstatus) <= 0 {
                    /*  test suffix  */
                    let ival = ival as usize;
                    if ival <= nend && ival >= nstart {
                        ffpsvc_safe(&card, &mut svalue, Some(&mut comm), status); /*  parse the value */
                        ffc2l(&svalue, &mut value[ival - nstart], status); /* convert*/
                        if (ival - nstart + 1) as c_int > *nfound {
                            *nfound = (ival - nstart + 1) as c_int; /*  max found */
                        }

                        if *status == VALUE_UNDEFINED {
                            undefinedval = true;
                            *status = 0; /* reset status to read remaining values */
                        }
                    }
                }
            }
        }
        if undefinedval && (*status <= 0) {
            *status = VALUE_UNDEFINED; /* report at least 1 value undefined */
        }

        *status
    }
}

/*--------------------------------------------------------------------------*/
/// Read (get) an indexed array of keywords with index numbers between
/// NSTART and (NSTART + NMAX -1) inclusive.  
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgknj(
    fptr: *mut fitsfile,    /* I - FITS file pointer                    */
    keyname: *const c_char, /* I - root name of keywords to read        */
    nstart: c_int,          /* I - starting index number                */
    nmax: c_int,            /* I - maximum number of keywords to return */
    value: *mut c_long,     /* O - array of keyword values              */
    nfound: *mut c_int,     /* O - number of values that were returned  */
    status: *mut c_int,     /* IO - error status                        */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let nfound = nfound.as_mut().expect(NULL_MSG);

        raw_to_slice!(keyname);

        let value = slice::from_raw_parts_mut(value, nmax as usize);

        ffgknj_safe(fptr, keyname, nstart, nmax, value, nfound, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Read (get) an indexed array of keywords with index numbers between
/// NSTART and (NSTART + NMAX -1) inclusive.  
pub fn ffgknj_safe(
    fptr: &mut fitsfile,  /* I - FITS file pointer                    */
    keyname: &[c_char],   /* I - root name of keywords to read        */
    nstart: c_int,        /* I - starting index number                */
    nmax: c_int,          /* I - maximum number of keywords to return */
    value: &mut [c_long], /* O - array of keyword values              */
    nfound: &mut c_int,   /* O - number of values that were returned  */
    status: &mut c_int,   /* IO - error status                        */
) -> c_int {
    let mut nkeys = 0;
    let mut mkeys = 0;
    let mut ival: c_long = 0;

    let mut keyroot: [c_char; FLEN_KEYWORD] = [0; FLEN_KEYWORD];
    let mut keyindex: [c_char; 8] = [0; 8];
    let mut card: [c_char; FLEN_CARD] = [0; FLEN_CARD];
    let mut svalue: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut comm: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];

    let nmax = nmax as usize;
    let nstart = nstart as usize;

    // *equalssign;

    if *status > 0 {
        return *status;
    }

    *nfound = 0;
    let nend = nstart + nmax - 1;

    keyroot[0] = 0;
    strncat_safe(&mut keyroot, keyname, FLEN_KEYWORD - 1);

    let lenroot = strlen_safe(&keyroot);

    if lenroot == 0 {
        /*  root must be at least 1 char long  */
        return *status;
    }
    for ii in 0..lenroot {
        /*  make sure upper case  */
        keyroot[ii] = toupper(keyroot[ii]);
    }

    ffghps_safe(fptr, Some(&mut nkeys), Some(&mut mkeys), status); /*  get the number of keywords  */
    ffmaky_safe(fptr, 3, status); /* move to 3rd keyword (skip 1st 2 keywords) */

    let mut undefinedval = false;

    for ii in 3..=(nkeys) {
        if ffgnky(fptr, &mut card, status) > 0 {
            /*  get next keyword  */
            return *status;
        }

        if strncmp_safe(&keyroot, &card, lenroot) == 0 {
            /* see if keyword matches */
            keyindex[0] = 0;
            let equalssign = strchr_safe(&card, bb(b'='));
            if equalssign.is_none() {
                continue; /* keyword has no value */
            }
            let equalssign = equalssign.unwrap();
            if equalssign - lenroot > 7 {
                *status = BAD_KEYCHAR;
                return *status;
            }

            strncat_safe(&mut keyindex, &card[lenroot..], equalssign - lenroot); /*  copy suffix  */

            let mut tstatus = 0;
            if ffc2ii(&keyindex, &mut ival, &mut tstatus) <= 0 {
                /*  test suffix  */
                let ival: usize = ival.try_into().unwrap();
                if ival <= nend && ival >= nstart {
                    ffpsvc_safe(&card, &mut svalue, Some(&mut comm), status); /*  parse the value */

                    ffc2i(&svalue, &mut value[ival - nstart], status); /* convert */

                    if (ival - nstart + 1) as c_int > *nfound {
                        /*  max found */
                        *nfound = (ival - nstart + 1) as c_int;
                    }

                    if *status == VALUE_UNDEFINED {
                        undefinedval = true;
                        *status = 0; /* reset status to read remaining values */
                    };
                };
            };
        };
    }
    if undefinedval && (*status <= 0) {
        *status = VALUE_UNDEFINED;
    }
    *status
}

/*--------------------------------------------------------------------------*/
/// Read (get) an indexed array of keywords with index numbers between
/// NSTART and (NSTART + NMAX -1) inclusive.  
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgknjj(
    fptr: *mut fitsfile,    /* I - FITS file pointer                    */
    keyname: *const c_char, /* I - root name of keywords to read        */
    nstart: c_int,          /* I - starting index number                */
    nmax: c_int,            /* I - maximum number of keywords to return */
    value: *mut LONGLONG,   /* O - array of keyword values              */
    nfound: *mut c_int,     /* O - number of values that were returned  */
    status: *mut c_int,     /* IO - error status                        */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let nfound = nfound.as_mut().expect(NULL_MSG);

        raw_to_slice!(keyname);

        let value = slice::from_raw_parts_mut(value, nmax as usize);

        ffgknjj_safe(fptr, keyname, nstart, nmax, value, nfound, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Read (get) an indexed array of keywords with index numbers between
/// NSTART and (NSTART + NMAX -1) inclusive.  
pub fn ffgknjj_safe(
    fptr: &mut fitsfile,    /* I - FITS file pointer                    */
    keyname: &[c_char],     /* I - root name of keywords to read        */
    nstart: c_int,          /* I - starting index number                */
    nmax: c_int,            /* I - maximum number of keywords to return */
    value: &mut [LONGLONG], /* O - array of keyword values              */
    nfound: &mut c_int,     /* O - number of values that were returned  */
    status: &mut c_int,     /* IO - error status                        */
) -> c_int {
    let mut nkeys = 0;
    let mut mkeys = 0;
    let mut ival: c_long = 0;

    let mut keyroot: [c_char; FLEN_KEYWORD] = [0; FLEN_KEYWORD];
    let mut keyindex: [c_char; 8] = [0; 8];
    let mut card: [c_char; FLEN_CARD] = [0; FLEN_CARD];
    let mut svalue: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut comm: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];

    let nmax = nmax as usize;
    let nstart = nstart as usize;

    // *equalssign;
    if *status > 0 {
        return *status;
    }

    *nfound = 0;
    let nend = nstart + nmax - 1;

    keyroot[0] = 0;
    strncat_safe(&mut keyroot, keyname, FLEN_KEYWORD - 1);

    let lenroot = strlen_safe(&keyroot);

    if lenroot == 0 {
        /*  root must be at least 1 char long  */
        return *status;
    }
    for ii in 0..lenroot {
        /*  make sure upper case  */
        keyroot[ii] = toupper(keyroot[ii]);
    }

    ffghps_safe(fptr, Some(&mut nkeys), Some(&mut mkeys), status); /*  get the number of keywords  */
    ffmaky_safe(fptr, 3, status); /* move to 3rd keyword (skip 1st 2 keywords) */

    let mut undefinedval = false;

    for ii in 3..=(nkeys) {
        if ffgnky(fptr, &mut card, status) > 0 {
            /*  get next keyword  */
            return *status;
        }

        if strncmp_safe(&keyroot, &card, lenroot) == 0 {
            /* see if keyword matches */
            keyindex[0] = 0;
            let equalssign = strchr_safe(&card, bb(b'='));
            if equalssign.is_none() {
                continue; /* keyword has no value */
            }
            let equalssign = equalssign.unwrap();
            if equalssign - lenroot > 7 {
                *status = BAD_KEYCHAR;
                return *status;
            }

            strncat_safe(&mut keyindex, &card[lenroot..], equalssign - lenroot); /*  copy suffix  */

            let mut tstatus = 0;
            if ffc2ii(&keyindex, &mut ival, &mut tstatus) <= 0 {
                /*  test suffix  */
                let ival: usize = ival.try_into().unwrap();
                if ival <= nend && ival >= nstart {
                    ffpsvc_safe(&card, &mut svalue, Some(&mut comm), status); /*  parse the value */

                    ffc2j(&svalue, &mut value[ival - nstart], status); /* convert */

                    if (ival - nstart + 1) as c_int > *nfound {
                        /*  max found */
                        *nfound = (ival - nstart + 1) as c_int;
                    }

                    if *status == VALUE_UNDEFINED {
                        undefinedval = true;
                        *status = 0; /* reset status to read remaining values */
                    };
                };
            };
        };
    }
    if undefinedval && (*status <= 0) {
        *status = VALUE_UNDEFINED;
    }
    *status
}

/*--------------------------------------------------------------------------*/
/// Read (get) an indexed array of keywords with index numbers between
/// NSTART and (NSTART + NMAX -1) inclusive.  
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgkne(
    fptr: *mut fitsfile,    /* I - FITS file pointer                    */
    keyname: *const c_char, /* I - root name of keywords to read        */
    nstart: c_int,          /* I - starting index number                */
    nmax: c_int,            /* I - maximum number of keywords to return */
    value: *mut f32,        /* O - array of keyword values              */
    nfound: *mut c_int,     /* O - number of values that were returned  */
    status: *mut c_int,     /* IO - error status                        */
) -> c_int {
    unsafe {
        let mut nkeys = 0;
        let mut mkeys = 0;
        let mut ival: c_long = 0;

        let mut keyroot: [c_char; FLEN_KEYWORD] = [0; FLEN_KEYWORD];
        let mut keyindex: [c_char; 8] = [0; 8];
        let mut card: [c_char; FLEN_CARD] = [0; FLEN_CARD];
        let mut svalue: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
        let mut comm: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];

        let nmax = nmax as usize;
        let nstart = nstart as usize;

        // *equalssign;

        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);
        raw_to_slice!(keyname);

        let value = slice::from_raw_parts_mut(value, nmax);

        if *status > 0 {
            return *status;
        }

        *nfound = 0;
        let nend = nstart + nmax - 1;

        keyroot[0] = 0;
        strncat_safe(&mut keyroot, keyname, FLEN_KEYWORD - 1);

        let lenroot = strlen_safe(&keyroot);

        if lenroot == 0 {
            /*  root must be at least 1 char long  */
            return *status;
        }
        for ii in 0..lenroot {
            /*  make sure upper case  */
            keyroot[ii] = toupper(keyroot[ii]);
        }

        ffghps_safe(fptr, Some(&mut nkeys), Some(&mut mkeys), status); /*  get the number of keywords  */
        ffmaky_safe(fptr, 3, status); /* move to 3rd keyword (skip 1st 2 keywords) */

        let mut undefinedval = false;

        for ii in 3..=(nkeys) {
            if ffgnky(fptr, &mut card, status) > 0 {
                /*  get next keyword  */
                return *status;
            }

            if strncmp_safe(&keyroot, &card, lenroot) == 0 {
                /* see if keyword matches */
                keyindex[0] = 0;
                let equalssign = strchr_safe(&card, bb(b'='));
                if equalssign.is_none() {
                    continue; /* keyword has no value */
                }
                let equalssign = equalssign.unwrap();
                if equalssign - lenroot > 7 {
                    *status = BAD_KEYCHAR;
                    return *status;
                }

                strncat_safe(&mut keyindex, &card[lenroot..], equalssign - lenroot); /*  copy suffix  */

                let mut tstatus = 0;
                if ffc2ii(&keyindex, &mut ival, &mut tstatus) <= 0 {
                    /*  test suffix  */
                    let ival: usize = ival.try_into().unwrap();
                    if ival <= nend && ival >= nstart {
                        ffpsvc_safe(&card, &mut svalue, Some(&mut comm), status); /*  parse the value */

                        ffc2r(&svalue, &mut value[ival - nstart], status); /* convert */

                        if (ival - nstart + 1) as c_int > *nfound {
                            /*  max found */
                            *nfound = (ival - nstart + 1) as c_int;
                        }

                        if *status == VALUE_UNDEFINED {
                            undefinedval = true;
                            *status = 0; /* reset status to read remaining values */
                        };
                    };
                };
            };
        }
        if undefinedval && (*status <= 0) {
            *status = VALUE_UNDEFINED;
        }
        *status
    }
}

/*--------------------------------------------------------------------------*/
/// Read (get) an indexed array of keywords with index numbers between
/// NSTART and (NSTART + NMAX -1) inclusive.  
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgknd(
    fptr: *mut fitsfile,    /* I - FITS file pointer                    */
    keyname: *const c_char, /* I - root name of keywords to read        */
    nstart: c_int,          /* I - starting index number                */
    nmax: c_int,            /* I - maximum number of keywords to return */
    value: *mut f64,        /* O - array of keyword values              */
    nfound: *mut c_int,     /* O - number of values that were returned  */
    status: *mut c_int,     /* IO - error status                        */
) -> c_int {
    unsafe {
        let mut nkeys = 0;
        let mut mkeys = 0;
        let mut ival: c_long = 0;

        let mut keyroot: [c_char; FLEN_KEYWORD] = [0; FLEN_KEYWORD];
        let mut keyindex: [c_char; 8] = [0; 8];
        let mut card: [c_char; FLEN_CARD] = [0; FLEN_CARD];
        let mut svalue: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
        let mut comm: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];

        let nmax = nmax as usize;
        let nstart = nstart as usize;

        // *equalssign;

        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);
        raw_to_slice!(keyname);

        let value = slice::from_raw_parts_mut(value, nmax);

        if *status > 0 {
            return *status;
        }

        *nfound = 0;
        let nend = nstart + nmax - 1;

        keyroot[0] = 0;
        strncat_safe(&mut keyroot, keyname, FLEN_KEYWORD - 1);

        let lenroot = strlen_safe(&keyroot);

        if lenroot == 0 {
            /*  root must be at least 1 char long  */
            return *status;
        }
        for ii in 0..lenroot {
            /*  make sure upper case  */
            keyroot[ii] = toupper(keyroot[ii]);
        }

        ffghps_safe(fptr, Some(&mut nkeys), Some(&mut mkeys), status); /*  get the number of keywords  */
        ffmaky_safe(fptr, 3, status); /* move to 3rd keyword (skip 1st 2 keywords) */

        let mut undefinedval = false;

        for ii in 3..=(nkeys) {
            if ffgnky(fptr, &mut card, status) > 0 {
                /*  get next keyword  */
                return *status;
            }

            if strncmp_safe(&keyroot, &card, lenroot) == 0 {
                /* see if keyword matches */
                keyindex[0] = 0;
                let equalssign = strchr_safe(&card, bb(b'='));
                if equalssign.is_none() {
                    continue; /* keyword has no value */
                }
                let equalssign = equalssign.unwrap();
                if equalssign - lenroot > 7 {
                    *status = BAD_KEYCHAR;
                    return *status;
                }

                strncat_safe(&mut keyindex, &card[lenroot..], equalssign - lenroot); /*  copy suffix  */

                let mut tstatus = 0;
                if ffc2ii(&keyindex, &mut ival, &mut tstatus) <= 0 {
                    /*  test suffix  */
                    let ival: usize = ival.try_into().unwrap();
                    if ival <= nend && ival >= nstart {
                        /* is index within range? */

                        ffpsvc_safe(&card, &mut svalue, Some(&mut comm), status); /*  parse the value */

                        ffc2d(&svalue, &mut value[ival - nstart], status); /* convert */

                        if (ival - nstart + 1) as c_int > *nfound {
                            /*  max found */
                            *nfound = (ival - nstart + 1) as c_int;
                        }

                        if *status == VALUE_UNDEFINED {
                            undefinedval = true;
                            *status = 0; /* reset status to read remaining values */
                        };
                    };
                };
            };
        }
        if undefinedval && (*status <= 0) {
            *status = VALUE_UNDEFINED;
        }
        *status
    }
}

/*--------------------------------------------------------------------------*/
/// read and parse the TDIMnnn keyword to get the dimensionality of a column
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgtdm(
    fptr: *mut fitsfile, /* I - FITS file pointer                        */
    colnum: c_int,       /* I - number of the column to read             */
    maxdim: c_int,       /* I - maximum no. of dimensions to read;       */
    naxis: *mut c_int,   /* O - number of axes in the data array         */
    naxes: *mut c_long,  /* O - length of each data axis                 */
    status: *mut c_int,  /* IO - error status                            */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let naxis = naxis.as_mut().expect(NULL_MSG);

        let naxes = slice::from_raw_parts_mut(naxes, maxdim as usize);

        ffgtdm_safe(fptr, colnum, maxdim, naxis, naxes, status)
    }
}

/*--------------------------------------------------------------------------*/
/// read and parse the TDIMnnn keyword to get the dimensionality of a column
pub fn ffgtdm_safe(
    fptr: &mut fitsfile,  /* I - FITS file pointer                        */
    colnum: c_int,        /* I - number of the column to read             */
    maxdim: c_int,        /* I - maximum no. of dimensions to read;       */
    naxis: &mut c_int,    /* O - number of axes in the data array         */
    naxes: &mut [c_long], /* O - length of each data axis                 */
    status: &mut c_int,   /* IO - error status                            */
) -> c_int {
    let mut tstatus = 0;
    let mut keyname: [c_char; FLEN_KEYWORD] = [0; FLEN_KEYWORD];
    let mut tdimstr: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];

    if *status > 0 {
        return *status;
    }

    ffkeyn_safe(cs!(c"TDIM"), colnum, &mut keyname, status); /* construct keyword name */

    ffgkys_safe(fptr, &keyname, &mut tdimstr, None, &mut tstatus); /* try reading keyword */

    ffdtdm_safe(fptr, &tdimstr, colnum, maxdim, naxis, naxes, status); /* decode it */

    *status
}

/*--------------------------------------------------------------------------*/
/// read and parse the TDIMnnn keyword to get the dimensionality of a column
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgtdmll(
    fptr: *mut fitsfile,  /* I - FITS file pointer                        */
    colnum: c_int,        /* I - number of the column to read             */
    maxdim: c_int,        /* I - maximum no. of dimensions to read;       */
    naxis: *mut c_int,    /* O - number of axes in the data array         */
    naxes: *mut LONGLONG, /* O - length of each data axis                 */
    status: *mut c_int,   /* IO - error status                            */
) -> c_int {
    unsafe {
        let mut tstatus = 0;
        let mut keyname: [c_char; FLEN_KEYWORD] = [0; FLEN_KEYWORD];
        let mut tdimstr: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];

        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let naxis = naxis.as_mut().expect(NULL_MSG);

        let naxes = slice::from_raw_parts_mut(naxes, maxdim as usize);

        if *status > 0 {
            return *status;
        }

        ffkeyn_safe(cs!(c"TDIM"), colnum, &mut keyname, status); /* construct keyword name */

        ffgkys_safe(fptr, &keyname, &mut tdimstr, None, &mut tstatus); /* try reading keyword */

        ffdtdmll_safe(fptr, &tdimstr, colnum, maxdim, naxis, naxes, status); /* decode it */

        *status
    }
}

/*--------------------------------------------------------------------------*/
/// Decode the TDIMnnn keyword to get the dimensionality of a column.
/// Check that the value is legal and consistent with the TFORM value.
/// If colnum = 0, then the validity checking is disabled.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffdtdm(
    fptr: *mut fitsfile,    /* I - FITS file pointer                        */
    tdimstr: *const c_char, /* I - TDIMn keyword value string. e.g. (10,10) */
    colnum: c_int,          /* I - number of the column             */
    maxdim: c_int,          /* I - maximum no. of dimensions to read;       */
    naxis: *mut c_int,      /* O - number of axes in the data array         */
    naxes: *mut c_long,     /* O - length of each data axis                 */
    status: *mut c_int,     /* IO - error status                            */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let naxis = naxis.as_mut().expect(NULL_MSG);

        raw_to_slice!(tdimstr);

        let naxes = slice::from_raw_parts_mut(naxes, maxdim as usize);

        ffdtdm_safe(fptr, tdimstr, colnum, maxdim, naxis, naxes, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Decode the TDIMnnn keyword to get the dimensionality of a column.
/// Check that the value is legal and consistent with the TFORM value.
/// If colnum = 0, then the validity checking is disabled.
pub fn ffdtdm_safe(
    fptr: &mut fitsfile,  /* I - FITS file pointer                        */
    tdimstr: &[c_char],   /* I - TDIMn keyword value string. e.g. (10,10) */
    colnum: c_int,        /* I - number of the column             */
    maxdim: c_int,        /* I - maximum no. of dimensions to read;       */
    naxis: &mut c_int,    /* O - number of axes in the data array         */
    naxes: &mut [c_long], /* O - length of each data axis                 */
    status: &mut c_int,   /* IO - error status                            */
) -> c_int {
    let mut lastloc = 0;
    let mut dimsize = 0;
    let mut totalpix = 1;

    let mut message: [c_char; FLEN_ERRMSG] = [0; FLEN_ERRMSG];

    if *status > 0 {
        return *status;
    }

    if colnum != 0 {
        if fptr.HDUposition != fptr.Fptr.curhdu {
            ffmahd_safe(fptr, (fptr.HDUposition) + 1, None, status);
        }

        if colnum < 1 || colnum > fptr.Fptr.tfield {
            *status = BAD_COL_NUM;
            return *status;
        }

        let c = fptr.Fptr.get_tableptr_as_slice(); /* set pointer to first column */
        let ci = colnum as usize - 1; /* increment to the correct column */

        if tdimstr[0] == 0 {
            /* TDIMn keyword doesn't exist? */
            *naxis = 1; /* default = 1 dimensional */
            if maxdim > 0 {
                naxes[0] = c[ci].trepeat as c_long; /* default length = repeat */
            }

            return *status;
        }
    }

    *naxis = 0;

    let mut loc = strchr_safe(tdimstr, bb(b'(')); /* find the opening quote */
    if loc.is_none() {
        int_snprintf!(
            &mut message,
            FLEN_ERRMSG,
            "Illegal dimensions format: {}",
            slice_to_str!(&tdimstr),
        );
        *status = BAD_TDIM;
        return *status;
    }

    let mut carry_loc = 0;
    while let Some(loc_inner) = loc {
        let mut loc_inner = carry_loc + loc_inner;
        loc_inner += 1;

        let mut loc_offset = 0;

        /* read size of next dimension */
        // dimsize = strtol_safe(&tdimstr[loc_inner..], &mut loc_offset, 10);
        let (r, p) = strtol_safe(&tdimstr[loc_inner..]).unwrap();
        dimsize = r;
        loc_offset = p;

        loc_inner += loc_offset;

        if *naxis < maxdim {
            naxes[(*naxis) as usize] = dimsize;
        }

        if dimsize < 0 {
            ffpmsg_str("one or more dimension are less than 0 (ffdtdm)");
            ffpmsg_slice(tdimstr);
            *status = BAD_TDIM;
            return *status;
        }

        totalpix *= dimsize;
        (*naxis) += 1;
        lastloc = loc_inner;
        carry_loc = loc_inner;
        let comma_loc = strchr_safe(&tdimstr[loc_inner..], bb(b',')); /* look for comma before next dimension */
        loc = comma_loc;
    }

    let loc = strchr_safe(&tdimstr[lastloc..], bb(b')')); /* check for the closing quote */
    if loc.is_none() {
        int_snprintf!(
            &mut message,
            FLEN_ERRMSG,
            "Illegal dimensions format: {}",
            slice_to_str!(&tdimstr),
        );
        *status = BAD_TDIM;
        return *status;
    }

    if colnum != 0 {
        let c = fptr.Fptr.get_tableptr_as_slice(); /* set pointer to first column */
        let ci = colnum as usize - 1; /* increment to the correct column */

        if (c[ci].tdatatype > 0) && (c[ci].trepeat != totalpix as LONGLONG) {
            int_snprintf!(
                &mut message,
                FLEN_ERRMSG,
                "column vector length, {}, does not equal TDIMn array size, {}",
                c[ci].trepeat,
                totalpix,
            );
            ffpmsg_slice(&message);
            ffpmsg_slice(tdimstr);
            *status = BAD_TDIM;
            return *status;
        }
    }
    *status
}

/*--------------------------------------------------------------------------*/
/// Decode the TDIMnnn keyword to get the dimensionality of a column.
/// Check that the value is legal and consistent with the TFORM value.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffdtdmll(
    fptr: *mut fitsfile,    /* I - FITS file pointer                        */
    tdimstr: *const c_char, /* I - TDIMn keyword value string. e.g. (10,10) */
    colnum: c_int,          /* I - number of the column             */
    maxdim: c_int,          /* I - maximum no. of dimensions to read;       */
    naxis: *mut c_int,      /* O - number of axes in the data array         */
    naxes: *mut LONGLONG,   /* O - length of each data axis                 */
    status: *mut c_int,     /* IO - error status                            */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let naxis = naxis.as_mut().expect(NULL_MSG);

        raw_to_slice!(tdimstr);

        let naxes = slice::from_raw_parts_mut(naxes, maxdim as usize);

        ffdtdmll_safe(fptr, tdimstr, colnum, maxdim, naxis, naxes, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Decode the TDIMnnn keyword to get the dimensionality of a column.
/// Check that the value is legal and consistent with the TFORM value.
pub fn ffdtdmll_safe(
    fptr: &mut fitsfile,    /* I - FITS file pointer                        */
    tdimstr: &[c_char],     /* I - TDIMn keyword value string. e.g. (10,10) */
    colnum: c_int,          /* I - number of the column             */
    maxdim: c_int,          /* I - maximum no. of dimensions to read;       */
    naxis: &mut c_int,      /* O - number of axes in the data array         */
    naxes: &mut [LONGLONG], /* O - length of each data axis                 */
    status: &mut c_int,     /* IO - error status                            */
) -> c_int {
    let mut lastloc = 0;
    let mut dimsize = 0;
    let mut totalpix = 1;
    let mut doublesize: f64 = 0.0;

    let mut message: [c_char; FLEN_ERRMSG] = [0; FLEN_ERRMSG];

    if *status > 0 {
        return *status;
    }

    if fptr.HDUposition != fptr.Fptr.curhdu {
        ffmahd_safe(fptr, (fptr.HDUposition) + 1, None, status);
    }

    if colnum < 1 || colnum > fptr.Fptr.tfield {
        *status = BAD_COL_NUM;
        return *status;
    }

    let c = fptr.Fptr.get_tableptr_as_slice(); /* set pointer to first column */
    let ci = colnum as usize - 1; /* increment to the correct column */

    if tdimstr[0] == 0 {
        /* TDIMn keyword doesn't exist? */
        *naxis = 1; /* default = 1 dimensional */
        if maxdim > 0 {
            naxes[0] = c[ci].trepeat as LONGLONG; /* default length = repeat */
        }
    } else {
        *naxis = 0;

        let mut loc = strchr_safe(tdimstr, bb(b'(')); /* find the opening quote */
        if loc.is_none() {
            int_snprintf!(
                &mut message,
                FLEN_ERRMSG,
                "Illegal TDIM keyword value: {}",
                slice_to_str!(&tdimstr),
            );
            *status = BAD_TDIM;
            return *status;
        }

        while let Some(mut loc_inner) = loc {
            loc_inner += 1;

            /* Read value as a double because the string to 64-bit int function is  */
            /* platform dependent (strtoll, strtol, _atoI64).  This still gives     */
            /* about 48 bits of precision, which is plenty for this purpose.        */

            doublesize = strtod_safe(&tdimstr[loc_inner..], &mut loc_inner);
            dimsize = (doublesize + 0.1) as LONGLONG;

            if *naxis < maxdim {
                naxes[(*naxis) as usize] = dimsize;
            }

            if dimsize < 0 {
                ffpmsg_str("one or more TDIM values are less than 0 (ffdtdmll)");
                ffpmsg_slice(tdimstr);
                *status = BAD_TDIM;
                return *status;
            }

            totalpix *= dimsize;
            (*naxis) += 1;
            lastloc = loc_inner;
            loc = strchr_safe(&tdimstr[loc_inner..], bb(b',')); /* look for comma before next dimension */
        }

        let loc = strchr_safe(&tdimstr[lastloc..], bb(b')')); /* check for the closing quote */
        if loc.is_none() {
            int_snprintf!(
                &mut message,
                FLEN_ERRMSG,
                "Illegal TDIM keyword value: {}",
                slice_to_str!(&tdimstr),
            );
            *status = BAD_TDIM;
            return *status;
        }

        if (c[ci].tdatatype > 0) && (c[ci].trepeat != totalpix) {
            int_snprintf!(
                &mut message,
                FLEN_ERRMSG,
                "column vector length, {:.0}, does not equal TDIMn array size, {:.0}",
                c[ci].trepeat as f64,
                totalpix as f64,
            );
            ffpmsg_slice(&message);
            ffpmsg_slice(tdimstr);
            *status = BAD_TDIM;
            return *status;
        }
    }
    *status
}

/*--------------------------------------------------------------------------*/
/// Get keywords from the Header of the Primary array:
/// Check that the keywords conform to the FITS standard and return the
/// parameters which determine the size and structure of the primary array
/// or IMAGE extension.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffghpr(
    fptr: *mut fitsfile, /* I - FITS file pointer                        */
    maxdim: c_int,       /* I - maximum no. of dimensions to read;       */
    simple: *mut c_int,  /* O - does file conform to FITS standard? 1/0  */
    bitpix: *mut c_int,  /* O - number of bits per data value pixel      */
    naxis: *mut c_int,   /* O - number of axes in the data array         */
    naxes: *mut c_long,  /* O - length of each data axis                 */
    pcount: *mut c_long, /* O - number of group parameters (usually 0)   */
    gcount: *mut c_long, /* O - number of random groups (usually 1 or 0) */
    extend: *mut c_int,  /* O - may FITS file haave extensions?          */
    status: *mut c_int,  /* IO - error status                            */
) -> c_int {
    unsafe {
        let mut idummy: c_int = 0;
        let mut lldummy: LONGLONG = 0;
        let mut ddummy1: f64 = 0.0;
        let mut ddummy2: f64 = 0.0;
        let mut tnaxes: [LONGLONG; 99] = [0; 99];

        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let simple = simple.as_mut().expect(NULL_MSG);
        let bitpix = bitpix.as_mut().expect(NULL_MSG);
        let mut naxis = naxis.as_mut();
        let pcount = pcount.as_mut().expect(NULL_MSG);
        let gcount = gcount.as_mut().expect(NULL_MSG);
        let extend = extend.as_mut().expect(NULL_MSG);

        ffgphd(
            fptr,
            maxdim,
            simple,
            bitpix,
            naxis.as_deref_mut(),
            &mut tnaxes,
            pcount,
            gcount,
            extend,
            &mut ddummy1,
            &mut ddummy2,
            &mut lldummy,
            &mut idummy,
            status,
        );

        if let Some(naxis) = naxis
            && !naxes.is_null()
        {
            let naxes = slice::from_raw_parts_mut(naxes, *naxis as usize);

            let mut ii = 0;
            while (ii < *naxis) && (ii < maxdim) {
                naxes[ii as usize] = tnaxes[ii as usize] as c_long;
                ii += 1;
            }
        } else if !naxes.is_null() {
            let naxes = slice::from_raw_parts_mut(naxes, maxdim as usize);
            for ii in 0..(maxdim as usize) {
                naxes[ii] = tnaxes[ii] as c_long;
            }
        }

        *status
    }
}

/*--------------------------------------------------------------------------*/
/// Get keywords from the Header of the PRimary array:
/// Check that the keywords conform to the FITS standard and return the
/// parameters which determine the size and structure of the primary array
/// or IMAGE extension.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffghprll(
    fptr: *mut fitsfile,  /* I - FITS file pointer                        */
    maxdim: c_int,        /* I - maximum no. of dimensions to read;       */
    simple: *mut c_int,   /* O - does file conform to FITS standard? 1/0  */
    bitpix: *mut c_int,   /* O - number of bits per data value pixel      */
    naxis: *mut c_int,    /* O - number of axes in the data array         */
    naxes: *mut LONGLONG, /* O - length of each data axis                 */
    pcount: *mut c_long,  /* O - number of group parameters (usually 0)   */
    gcount: *mut c_long,  /* O - number of random groups (usually 1 or 0) */
    extend: *mut c_int,   /* O - may FITS file haave extensions?          */
    status: *mut c_int,   /* IO - error status                            */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let simple = simple.as_mut().expect(NULL_MSG);
        let bitpix = bitpix.as_mut().expect(NULL_MSG);
        let naxis = naxis.as_mut();
        let pcount = pcount.as_mut().expect(NULL_MSG);
        let gcount = gcount.as_mut().expect(NULL_MSG);
        let extend = extend.as_mut().expect(NULL_MSG);

        let naxes = slice::from_raw_parts_mut(naxes, maxdim as usize);

        ffghprll_safe(
            fptr, maxdim, simple, bitpix, naxis, naxes, pcount, gcount, extend, status,
        )
    }
}

/*--------------------------------------------------------------------------*/
/// Get keywords from the Header of the PRimary array:
/// Check that the keywords conform to the FITS standard and return the
/// parameters which determine the size and structure of the primary array
/// or IMAGE extension.
pub fn ffghprll_safe(
    fptr: &mut fitsfile,       /* I - FITS file pointer                        */
    maxdim: c_int,             /* I - maximum no. of dimensions to read;       */
    simple: &mut c_int,        /* O - does file conform to FITS standard? 1/0  */
    bitpix: &mut c_int,        /* O - number of bits per data value pixel      */
    naxis: Option<&mut c_int>, /* O - number of axes in the data array         */
    naxes: &mut [LONGLONG],    /* O - length of each data axis                 */
    pcount: &mut c_long,       /* O - number of group parameters (usually 0)   */
    gcount: &mut c_long,       /* O - number of random groups (usually 1 or 0) */
    extend: &mut c_int,        /* O - may FITS file haave extensions?          */
    status: &mut c_int,        /* IO - error status                            */
) -> c_int {
    let mut idummy: c_int = 0;
    let mut lldummy: LONGLONG = 0;
    let mut ddummy1: f64 = 0.0;
    let mut ddummy2: f64 = 0.0;

    ffgphd(
        fptr,
        maxdim,
        simple,
        bitpix,
        naxis,
        naxes,
        pcount,
        gcount,
        extend,
        &mut ddummy1,
        &mut ddummy2,
        &mut lldummy,
        &mut idummy,
        status,
    );

    *status
}

/*--------------------------------------------------------------------------*/
/// Get keywords from the Header of the ASCII TaBle:
/// Check that the keywords conform to the FITS standard and return the
/// parameters which describe the table.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffghtb(
    fptr: *mut fitsfile,     /* I - FITS file pointer                        */
    maxfield: c_int,         /* I - maximum no. of columns to read;          */
    naxis1: *mut c_long,     /* O - length of table row in bytes             */
    naxis2: *mut c_long,     /* O - number of rows in the table              */
    tfields: *mut c_int,     /* O - number of columns in the table           */
    ttype: *mut *mut c_char, /* O - name of each column                      */
    tbcol: *mut c_long,      /* O - byte offset in row to each column        */
    tform: *mut *mut c_char, /* O - value of TFORMn keyword for each column  */
    tunit: *mut *mut c_char, /* O - value of TUNITn keyword for each column  */
    extnm: *mut c_char,      /* O - value of EXTNAME keyword, if any         */
    status: *mut c_int,      /* IO - error status                            */
) -> c_int {
    unsafe {
        let mut maxf: c_int = 0;
        let mut nfound: c_int = 0;
        let mut tstatus: c_int = 0;
        let mut fields: c_long = 0;
        let mut name: [c_char; FLEN_KEYWORD] = [0; FLEN_KEYWORD];
        let mut value: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
        let mut comm: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];
        let mut xtension: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
        let mut message: [c_char; FLEN_ERRMSG] = [0; FLEN_ERRMSG];
        let mut llnaxis1: LONGLONG = 0;
        let mut llnaxis2: LONGLONG = 0;
        let mut pcount: LONGLONG = 0;

        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        if *status > 0 {
            return *status;
        }

        /* read the first keyword of the extension */
        ffgkyn_safe(fptr, 1, &mut name, &mut value, Some(&mut comm), status);

        if strcmp_safe(&name, cs!(c"XTENSION")) == 0 {
            if ffc2s(&value, &mut xtension, status) > 0 {
                /* get the value string */
                ffpmsg_str("Bad value string for XTENSION keyword:");
                ffpmsg_slice(&value);
                return *status;
            }

            /* allow the quoted string value to begin in any column and */
            /* allow any number of trailing blanks before the closing quote */
            /* first char must be a quote */
            if (value[0] != bb(b'\'')) || (strcmp_safe(&xtension, cs!(c"TABLE")) != 0) {
                int_snprintf!(
                    &mut message,
                    FLEN_ERRMSG,
                    "This is not a TABLE extension: {}",
                    slice_to_str!(&value),
                );
                ffpmsg_slice(&message);
                *status = NOT_ATABLE;
                return *status;
            }
        } else {
            /* error: 1st keyword of extension != XTENSION */

            int_snprintf!(
                &mut message,
                FLEN_ERRMSG,
                "First keyword of the extension is not XTENSION: {}",
                slice_to_str!(&name),
            );
            ffpmsg_slice(&message);
            *status = NO_XTENSION;
            return *status;
        }

        if ffgttb(
            fptr,
            &mut llnaxis1,
            &mut llnaxis2,
            &mut pcount,
            &mut fields,
            status,
        ) > 0
        {
            return *status;
        }

        if let Some(naxis1) = naxis1.as_mut() {
            *naxis1 = llnaxis1 as c_long;
        }

        if let Some(naxis2) = naxis2.as_mut() {
            *naxis2 = llnaxis2 as c_long;
        }

        if pcount != 0 {
            int_snprintf!(
                &mut message,
                FLEN_ERRMSG,
                "PCOUNT = {:.0} is illegal in ASCII table; must = 0",
                pcount as f64,
            );
            ffpmsg_slice(&message);
            *status = BAD_PCOUNT;
            return *status;
        }

        if let Some(tfields) = tfields.as_mut() {
            *tfields = fields as c_int;
        }

        if maxfield < 0 {
            maxf = fields as c_int;
        } else {
            maxf = cmp::min(maxfield, fields as c_int);
        }

        if maxf > 0 {
            for ii in 0..(maxf as usize) {
                /* initialize optional keyword values */
                if !ttype.is_null() {
                    let ttype = slice::from_raw_parts_mut(ttype, maxf as usize);
                    *(ttype[ii]) = 0;
                }

                if !tunit.is_null() {
                    let tunit_arr = slice::from_raw_parts_mut(tunit, maxf as usize);
                    (*(tunit_arr[ii])) = 0;
                }
            }

            if !ttype.is_null() {
                let ttype = slice::from_raw_parts_mut(ttype, maxf as usize);
                let mut v_ttype = Vec::new();

                for item in ttype {
                    let ttype_item = slice::from_raw_parts_mut(*item, FLEN_VALUE);
                    v_ttype.push(ttype_item);
                }
                ffgkns_safe(
                    fptr,
                    cs!(c"TTYPE"),
                    1,
                    maxf,
                    &mut v_ttype,
                    &mut nfound,
                    status,
                );
            }

            if !tunit.is_null() {
                let tunit_arr = slice::from_raw_parts_mut(tunit, maxf as usize);
                let mut v_tunit = Vec::new();

                for item in tunit_arr {
                    let tunit_item = slice::from_raw_parts_mut(*item, FLEN_VALUE);
                    v_tunit.push(tunit_item);
                }
                ffgkns_safe(
                    fptr,
                    cs!(c"TUNIT"),
                    1,
                    maxf,
                    &mut v_tunit,
                    &mut nfound,
                    status,
                );
            }

            if *status > 0 {
                return *status;
            }

            if !tbcol.is_null() {
                let tbcol = slice::from_raw_parts_mut(tbcol, maxf as usize);

                ffgknj_safe(fptr, cs!(c"TBCOL"), 1, maxf, tbcol, &mut nfound, status);

                if *status > 0 || nfound != maxf {
                    ffpmsg_str(
                        "Required TBCOL keyword(s) not found in ASCII table header (ffghtb).",
                    );
                    *status = NO_TBCOL;
                    return *status;
                }
            }

            if !tform.is_null() {
                let tform = slice::from_raw_parts_mut(tform, maxf as usize);
                let mut v_tform = Vec::new();

                for item in tform {
                    let tform_item = slice::from_raw_parts_mut(*item, FLEN_VALUE);
                    v_tform.push(tform_item);
                }
                ffgkns_safe(
                    fptr,
                    cs!(c"TFORM"),
                    1,
                    maxf,
                    &mut v_tform,
                    &mut nfound,
                    status,
                );

                if *status > 0 || nfound != maxf {
                    ffpmsg_str(
                        "Required TFORM keyword(s) not found in ASCII table header (ffghtb).",
                    );
                    *status = NO_TFORM;
                    return *status;
                }
            }
        }

        if !extnm.is_null() {
            let extnm = slice::from_raw_parts_mut(extnm, 69);
            extnm[0] = 0;

            tstatus = *status;
            ffgkys_safe(fptr, cs!(c"EXTNAME"), extnm, Some(&mut comm), status);

            if *status == KEY_NO_EXIST {
                *status = tstatus; /* keyword not required, so ignore error */
            }
        }

        *status
    }
}

/*--------------------------------------------------------------------------*/
/// Get keywords from the Header of the ASCII TaBle:
/// Check that the keywords conform to the FITS standard and return the
/// parameters which describe the table.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffghtbll(
    fptr: *mut fitsfile,     /* I - FITS file pointer                        */
    maxfield: c_int,         /* I - maximum no. of columns to read;          */
    naxis1: *mut LONGLONG,   /* O - length of table row in bytes             */
    naxis2: *mut LONGLONG,   /* O - number of rows in the table              */
    tfields: *mut c_int,     /* O - number of columns in the table           */
    ttype: *mut *mut c_char, /* O - name of each column                      */
    tbcol: *mut LONGLONG,    /* O - byte offset in row to each column        */
    tform: *mut *mut c_char, /* O - value of TFORMn keyword for each column  */
    tunit: *mut *mut c_char, /* O - value of TUNITn keyword for each column  */
    extnm: *mut c_char,      /* O - value of EXTNAME keyword, if any         */
    status: *mut c_int,      /* IO - error status                            */
) -> c_int {
    unsafe {
        let mut maxf: c_int = 0;
        let mut nfound: c_int = 0;
        let mut tstatus: c_int = 0;
        let mut fields: c_long = 0;
        let mut name: [c_char; FLEN_KEYWORD] = [0; FLEN_KEYWORD];
        let mut value: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
        let mut comm: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];
        let mut xtension: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
        let mut message: [c_char; FLEN_ERRMSG] = [0; FLEN_ERRMSG];
        let mut llnaxis1: LONGLONG = 0;
        let mut llnaxis2: LONGLONG = 0;
        let mut pcount: LONGLONG = 0;

        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        if *status > 0 {
            return *status;
        }

        /* read the first keyword of the extension */
        ffgkyn_safe(fptr, 1, &mut name, &mut value, Some(&mut comm), status);

        if strcmp_safe(&name, cs!(c"XTENSION")) == 0 {
            if ffc2s(&value, &mut xtension, status) > 0 {
                /* get the value string */
                ffpmsg_str("Bad value string for XTENSION keyword:");
                ffpmsg_slice(&value);
                return *status;
            }

            /* allow the quoted string value to begin in any column and */
            /* allow any number of trailing blanks before the closing quote */
            /* first char must be a quote */
            if (value[0] != bb(b'\'')) || (strcmp_safe(&xtension, cs!(c"TABLE")) != 0) {
                int_snprintf!(
                    &mut message,
                    FLEN_ERRMSG,
                    "This is not a TABLE extension: {}",
                    slice_to_str!(&value),
                );
                ffpmsg_slice(&message);
                *status = NOT_ATABLE;
                return *status;
            }
        } else
        /* error: 1st keyword of extension != XTENSION */
        {
            int_snprintf!(
                &mut message,
                FLEN_ERRMSG,
                "First keyword of the extension is not XTENSION: {}",
                slice_to_str!(&name),
            );
            ffpmsg_slice(&message);
            *status = NO_XTENSION;
            return *status;
        }

        if ffgttb(
            fptr,
            &mut llnaxis1,
            &mut llnaxis2,
            &mut pcount,
            &mut fields,
            status,
        ) > 0
        {
            return *status;
        }

        if let Some(naxis1) = naxis1.as_mut() {
            *naxis1 = llnaxis1;
        }

        if let Some(naxis2) = naxis2.as_mut() {
            *naxis2 = llnaxis2;
        }

        if pcount != 0 {
            int_snprintf!(
                &mut message,
                FLEN_ERRMSG,
                "PCOUNT = {:.0} is illegal in ASCII table; must = 0",
                pcount as f64,
            );
            ffpmsg_slice(&message);
            *status = BAD_PCOUNT;
            return *status;
        }

        if let Some(tfields) = tfields.as_mut() {
            *tfields = fields as c_int;
        }

        if maxfield < 0 {
            maxf = fields as c_int;
        } else {
            maxf = cmp::min(maxfield, fields as c_int);
        }

        if maxf > 0 {
            for ii in 0..(maxf as usize) {
                /* initialize optional keyword values */
                if !ttype.is_null() {
                    let ttype = slice::from_raw_parts_mut(
                        ttype as *mut [c_char; FLEN_VALUE],
                        maxf as usize,
                    );
                    ttype[ii][0] = 0;
                }

                if !tunit.is_null() {
                    let tunit = slice::from_raw_parts_mut(
                        tunit as *mut [c_char; FLEN_VALUE],
                        maxf as usize,
                    );
                    tunit[ii][0] = 0;
                }
            }

            if !ttype.is_null() {
                let ttype = slice::from_raw_parts_mut(ttype, maxf as usize);
                let mut v_ttype = Vec::new();

                for item in ttype {
                    let ttype_item = slice::from_raw_parts_mut(*item, FLEN_VALUE);
                    v_ttype.push(ttype_item);
                }

                ffgkns_safe(
                    fptr,
                    cs!(c"TTYPE"),
                    1,
                    maxf,
                    &mut v_ttype,
                    &mut nfound,
                    status,
                );
            }

            if !tunit.is_null() {
                let tunit_arr = slice::from_raw_parts_mut(tunit, maxf as usize);
                let mut v_tunit = Vec::new();

                for item in tunit_arr {
                    let tunit_item = slice::from_raw_parts_mut(*item, FLEN_VALUE);
                    v_tunit.push(tunit_item);
                }
                ffgkns_safe(
                    fptr,
                    cs!(c"TUNIT"),
                    1,
                    maxf,
                    &mut v_tunit,
                    &mut nfound,
                    status,
                );
            }

            if *status > 0 {
                return *status;
            }

            if !tbcol.is_null() {
                let tbcol = slice::from_raw_parts_mut(tbcol, maxf as usize);

                ffgknjj_safe(fptr, cs!(c"TBCOL"), 1, maxf, tbcol, &mut nfound, status);

                if *status > 0 || nfound != maxf {
                    ffpmsg_str(
                        "Required TBCOL keyword(s) not found in ASCII table header (ffghtbll).",
                    );
                    *status = NO_TBCOL;
                    return *status;
                }
            }

            if !tform.is_null() {
                let tform = slice::from_raw_parts_mut(tform, maxf as usize);
                let mut v_tform = Vec::new();

                for item in tform {
                    let tform_item = slice::from_raw_parts_mut(*item, FLEN_VALUE);
                    v_tform.push(tform_item);
                }
                ffgkns_safe(
                    fptr,
                    cs!(c"TFORM"),
                    1,
                    maxf,
                    &mut v_tform,
                    &mut nfound,
                    status,
                );

                if *status > 0 || nfound != maxf {
                    ffpmsg_str(
                        "Required TFORM keyword(s) not found in ASCII table header (ffghtbll).",
                    );
                    *status = NO_TFORM;
                    return *status;
                }
            }
        }

        if !extnm.is_null() {
            let extnm = slice::from_raw_parts_mut(extnm, 69);
            extnm[0] = 0;

            tstatus = *status;
            ffgkys_safe(fptr, cs!(c"EXTNAME"), extnm, Some(&mut comm), status);

            if *status == KEY_NO_EXIST {
                *status = tstatus; /* keyword not required, so ignore error */
            }
        }

        *status
    }
}

/*--------------------------------------------------------------------------*/
/// Get keywords from the Header of the BiNary table:
/// Check that the keywords conform to the FITS standard and return the
/// parameters which describe the table.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffghbn(
    fptr: *mut fitsfile,     /* I - FITS file pointer                        */
    maxfield: c_int,         /* I - maximum no. of columns to read;          */
    naxis2: *mut c_long,     /* O - number of rows in the table              */
    tfields: *mut c_int,     /* O - number of columns in the table           */
    ttype: *mut *mut c_char, /* O - name of each column                      */
    tform: *mut *mut c_char, /* O - TFORMn value for each column             */
    tunit: *mut *mut c_char, /* O - TUNITn value for each column             */
    extnm: *mut c_char,      /* O - value of EXTNAME keyword, if any         */
    pcount: *mut c_long,     /* O - value of PCOUNT keyword                  */
    status: *mut c_int,      /* IO - error status                            */
) -> c_int {
    unsafe {
        let mut maxf: c_int = 0;
        let mut nfound: c_int = 0;
        let mut tstatus: c_int = 0;
        let mut fields: c_long = 0;
        let mut name: [c_char; FLEN_KEYWORD] = [0; FLEN_KEYWORD];
        let mut value: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
        let mut comm: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];
        let mut xtension: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
        let mut message: [c_char; FLEN_ERRMSG] = [0; FLEN_ERRMSG];
        let mut naxis1ll: LONGLONG = 0;
        let mut naxis2ll: LONGLONG = 0;
        let mut pcountll: LONGLONG = 0;

        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        if *status > 0 {
            return *status;
        }

        /* read the first keyword of the extension */
        ffgkyn_safe(fptr, 1, &mut name, &mut value, Some(&mut comm), status);

        if strcmp_safe(&name, cs!(c"XTENSION")) == 0 {
            if ffc2s(&value, &mut xtension, status) > 0 {
                /* get the value string */
                ffpmsg_str("Bad value string for XTENSION keyword:");
                ffpmsg_slice(&value);
                return *status;
            }

            /* allow the quoted string value to begin in any column and */
            /* allow any number of trailing blanks before the closing quote */
            /* first char must be a quote */
            if (value[0] != bb(b'\''))
                || (strcmp_safe(&xtension, cs!(c"BINTABLE")) != 0)
                    && (strcmp_safe(&xtension, cs!(c"A3DTABLE")) != 0)
                    && (strcmp_safe(&xtension, cs!(c"3DTABLE")) != 0)
            {
                int_snprintf!(
                    &mut message,
                    FLEN_ERRMSG,
                    "This is not a BINTABLE extension: {}",
                    slice_to_str!(&value),
                );
                ffpmsg_slice(&message);
                *status = NOT_BTABLE;
                return *status;
            }
        } else
        /* error: 1st keyword of extension != XTENSION */
        {
            int_snprintf!(
                &mut message,
                FLEN_ERRMSG,
                "First keyword of the extension is not XTENSION: {}",
                slice_to_str!(&name),
            );
            ffpmsg_slice(&message);
            *status = NO_XTENSION;
            return *status;
        }

        if ffgttb(
            fptr,
            &mut naxis1ll,
            &mut naxis2ll,
            &mut pcountll,
            &mut fields,
            status,
        ) > 0
        {
            return *status;
        }

        if let Some(naxis2) = naxis2.as_mut() {
            *naxis2 = naxis2ll as c_long;
        }

        if let Some(pcount) = pcount.as_mut() {
            *pcount = pcountll as c_long;
        }

        if let Some(tfields) = tfields.as_mut() {
            *tfields = fields as c_int;
        }

        if maxfield < 0 {
            maxf = fields as c_int;
        } else {
            maxf = cmp::min(maxfield, fields as c_int);
        }

        if maxf > 0 {
            for ii in 0..(maxf as usize) {
                /* initialize optional keyword values */
                if !ttype.is_null() {
                    let ttype = slice::from_raw_parts_mut(ttype, maxf as usize);
                    *(ttype[ii]) = 0;
                }

                if !tunit.is_null() {
                    let tunit = slice::from_raw_parts_mut(tunit, maxf as usize);
                    *(tunit[ii]) = 0;
                }
            }

            if !ttype.is_null() {
                let ttype = slice::from_raw_parts_mut(ttype, maxf as usize);
                let mut v_ttype = Vec::new();

                for item in ttype {
                    let ttype_item = slice::from_raw_parts_mut(*item, FLEN_VALUE);
                    v_ttype.push(ttype_item);
                }

                ffgkns_safe(
                    fptr,
                    cs!(c"TTYPE"),
                    1,
                    maxf,
                    &mut v_ttype,
                    &mut nfound,
                    status,
                );
            }

            if !tunit.is_null() {
                let tunit_arr = slice::from_raw_parts_mut(tunit, maxf as usize);
                let mut v_tunit = Vec::new();

                for item in tunit_arr {
                    let tunit_item = slice::from_raw_parts_mut(*item, FLEN_VALUE);
                    v_tunit.push(tunit_item);
                }
                ffgkns_safe(
                    fptr,
                    cs!(c"TUNIT"),
                    1,
                    maxf,
                    &mut v_tunit,
                    &mut nfound,
                    status,
                );
            }

            if *status > 0 {
                return *status;
            }

            if !tform.is_null() {
                let tform = slice::from_raw_parts_mut(tform, maxf as usize);
                let mut v_tform = Vec::new();

                for item in tform {
                    let tform_item = slice::from_raw_parts_mut(*item, FLEN_VALUE);
                    v_tform.push(tform_item);
                }
                ffgkns_safe(
                    fptr,
                    cs!(c"TFORM"),
                    1,
                    maxf,
                    &mut v_tform,
                    &mut nfound,
                    status,
                );

                if *status > 0 || nfound != maxf {
                    ffpmsg_str(
                        "Required TFORM keyword(s) not found in binary table header (ffghbn).",
                    );
                    *status = NO_TFORM;
                    return *status;
                }
            }
        }

        if !extnm.is_null() {
            let extnm = slice::from_raw_parts_mut(extnm, 69);
            extnm[0] = 0;

            tstatus = *status;
            ffgkys_safe(fptr, cs!(c"EXTNAME"), extnm, Some(&mut comm), status);

            if *status == KEY_NO_EXIST {
                *status = tstatus; /* keyword not required, so ignore error */
            }
        }

        *status
    }
}

/*--------------------------------------------------------------------------*/
/// Get keywords from the Header of the BiNary table:
/// Check that the keywords conform to the FITS standard and return the
/// parameters which describe the table.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffghbnll(
    fptr: *mut fitsfile,     /* I - FITS file pointer                        */
    maxfield: c_int,         /* I - maximum no. of columns to read;          */
    naxis2: *mut LONGLONG,   /* O - number of rows in the table              */
    tfields: *mut c_int,     /* O - number of columns in the table           */
    ttype: *mut *mut c_char, /* O - name of each column                      */
    tform: *mut *mut c_char, /* O - TFORMn value for each column             */
    tunit: *mut *mut c_char, /* O - TUNITn value for each column             */
    extnm: *mut c_char,      /* O - value of EXTNAME keyword, if any         */
    pcount: *mut LONGLONG,   /* O - value of PCOUNT keyword                  */
    status: *mut c_int,      /* IO - error status                            */
) -> c_int {
    unsafe {
        let mut maxf: c_int = 0;
        let mut nfound: c_int = 0;
        let mut tstatus: c_int = 0;
        let mut fields: c_long = 0;
        let mut name: [c_char; FLEN_KEYWORD] = [0; FLEN_KEYWORD];
        let mut value: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
        let mut comm: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];
        let mut xtension: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
        let mut message: [c_char; FLEN_ERRMSG] = [0; FLEN_ERRMSG];
        let mut naxis1ll: LONGLONG = 0;
        let mut naxis2ll: LONGLONG = 0;
        let mut pcountll: LONGLONG = 0;

        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        if *status > 0 {
            return *status;
        }

        /* read the first keyword of the extension */
        ffgkyn_safe(fptr, 1, &mut name, &mut value, Some(&mut comm), status);

        if strcmp_safe(&name, cs!(c"XTENSION")) == 0 {
            if ffc2s(&value, &mut xtension, status) > 0 {
                /* get the value string */
                ffpmsg_str("Bad value string for XTENSION keyword:");
                ffpmsg_slice(&value);
                return *status;
            }

            /* allow the quoted string value to begin in any column and */
            /* allow any number of trailing blanks before the closing quote */
            /* first char must be a quote */
            if (value[0] != bb(b'\''))
                || (strcmp_safe(&xtension, cs!(c"BINTABLE")) != 0)
                    && (strcmp_safe(&xtension, cs!(c"A3DTABLE")) != 0)
                    && (strcmp_safe(&xtension, cs!(c"3DTABLE")) != 0)
            {
                int_snprintf!(
                    &mut message,
                    FLEN_ERRMSG,
                    "This is not a BINTABLE extension: {}",
                    slice_to_str!(&value),
                );
                ffpmsg_slice(&message);
                *status = NOT_BTABLE;
                return *status;
            }
        } else
        /* error: 1st keyword of extension != XTENSION */
        {
            int_snprintf!(
                &mut message,
                FLEN_ERRMSG,
                "First keyword of the extension is not XTENSION: {}",
                slice_to_str!(&name),
            );
            ffpmsg_slice(&message);
            *status = NO_XTENSION;
            return *status;
        }

        if ffgttb(
            fptr,
            &mut naxis1ll,
            &mut naxis2ll,
            &mut pcountll,
            &mut fields,
            status,
        ) > 0
        {
            return *status;
        }

        if let Some(naxis2) = naxis2.as_mut() {
            *naxis2 = naxis2ll;
        }

        if let Some(pcount) = pcount.as_mut() {
            *pcount = pcountll;
        }

        if let Some(tfields) = tfields.as_mut() {
            *tfields = fields as c_int;
        }

        if maxfield < 0 {
            maxf = fields as c_int;
        } else {
            maxf = cmp::min(maxfield, fields as c_int);
        }

        if maxf > 0 {
            for ii in 0..(maxf as usize) {
                /* initialize optional keyword values */
                if !ttype.is_null() {
                    let ttype = slice::from_raw_parts_mut(
                        ttype as *mut [c_char; FLEN_VALUE],
                        maxf as usize,
                    );
                    ttype[ii][0] = 0;
                }

                if !tunit.is_null() {
                    let tunit = slice::from_raw_parts_mut(
                        tunit as *mut [c_char; FLEN_VALUE],
                        maxf as usize,
                    );
                    tunit[ii][0] = 0;
                }
            }

            if !ttype.is_null() {
                let ttype = slice::from_raw_parts_mut(ttype, maxf as usize);
                let mut v_ttype = Vec::new();

                for item in ttype {
                    let ttype_item = slice::from_raw_parts_mut(*item, FLEN_VALUE);
                    v_ttype.push(ttype_item);
                }
                ffgkns_safe(
                    fptr,
                    cs!(c"TTYPE"),
                    1,
                    maxf,
                    &mut v_ttype,
                    &mut nfound,
                    status,
                );
            }

            if !tunit.is_null() {
                let tunit_arr = slice::from_raw_parts_mut(tunit, maxf as usize);
                let mut v_tunit = Vec::new();

                for item in tunit_arr {
                    let tunit_item = slice::from_raw_parts_mut(*item, FLEN_VALUE);
                    v_tunit.push(tunit_item);
                }
                ffgkns_safe(
                    fptr,
                    cs!(c"TUNIT"),
                    1,
                    maxf,
                    &mut v_tunit,
                    &mut nfound,
                    status,
                );
            }

            if *status > 0 {
                return *status;
            }

            if !tform.is_null() {
                let tform = slice::from_raw_parts_mut(tform, maxf as usize);
                let mut v_tform = Vec::new();

                for item in tform {
                    let tform_item = slice::from_raw_parts_mut(*item, FLEN_VALUE);
                    v_tform.push(tform_item);
                }
                ffgkns_safe(
                    fptr,
                    cs!(c"TFORM"),
                    1,
                    maxf,
                    &mut v_tform,
                    &mut nfound,
                    status,
                );

                if *status > 0 || nfound != maxf {
                    ffpmsg_str(
                        "Required TFORM keyword(s) not found in binary table header (ffghbnll).",
                    );
                    *status = NO_TFORM;
                    return *status;
                }
            }
        }

        if !extnm.is_null() {
            let extnm = slice::from_raw_parts_mut(extnm, 69);
            extnm[0] = 0;

            tstatus = *status;
            ffgkys_safe(fptr, cs!(c"EXTNAME"), extnm, Some(&mut comm), status);

            if *status == KEY_NO_EXIST {
                *status = tstatus; /* keyword not required, so ignore error */
            }
        }

        *status
    }
}

/*--------------------------------------------------------------------------*/
/// Get the Primary HeaDer parameters.  Check that the keywords conform to
/// the FITS standard and return the parameters which determine the size and
/// structure of the primary array or IMAGE extension.
pub(crate) fn ffgphd(
    fptr: &mut fitsfile,       /* I - FITS file pointer                        */
    maxdim: c_int,             /* I - maximum no. of dimensions to read;       */
    simple: &mut c_int,        /* O - does file conform to FITS standard? 1/0  */
    bitpix: &mut c_int,        /* O - number of bits per data value pixel      */
    naxis: Option<&mut c_int>, /* O - number of axes in the data array         */
    naxes: &mut [LONGLONG],    /* O - length of each data axis                 */
    pcount: &mut c_long,       /* O - number of group parameters (usually 0)   */
    gcount: &mut c_long,       /* O - number of random groups (usually 1 or 0) */
    extend: &mut c_int,        /* O - may FITS file have extensions?          */
    bscale: &mut f64,          /* O - array pixel linear scaling factor        */
    bzero: &mut f64,           /* O - array pixel linear scaling zero point    */
    blank: &mut LONGLONG,      /* O - value used to represent undefined pixels */
    nspace: &mut c_int,        /* O - number of blank keywords prior to END    */
    status: &mut c_int,        /* IO - error status                            */
) -> c_int {
    let mut ii = 0;
    let mut nextkey = 0;
    let mut longbitpix: c_long = 0;
    let mut longnaxis: c_long = 0;
    let mut axislen: LONGLONG = 0;
    let mut message: [c_char; FLEN_ERRMSG] = [0; FLEN_ERRMSG];
    let mut keyword: [c_char; FLEN_KEYWORD] = [0; FLEN_KEYWORD];
    let mut card: [c_char; FLEN_CARD] = [0; FLEN_CARD];
    let mut name: [c_char; FLEN_KEYWORD] = [0; FLEN_KEYWORD];
    let mut value: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut comm: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];
    let mut xtension: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];

    if *status > 0 {
        return *status;
    }

    if fptr.HDUposition != fptr.Fptr.curhdu {
        ffmahd_safe(fptr, (fptr.HDUposition) + 1, None, status);
    }

    //if !simple.is_null() {
    *simple = 1;
    //}

    let mut unknown = 0;

    /*--------------------------------------------------------------------*/
    /*  Get 1st keyword of HDU and test whether it is SIMPLE or XTENSION  */
    /*--------------------------------------------------------------------*/

    ffgkyn_safe(fptr, 1, &mut name, &mut value, Some(&mut comm), status);

    if fptr.Fptr.curhdu == 0 {
        /* Is this the beginning of the FITS file? */
        if strcmp_safe(&name, cs!(c"SIMPLE")) == 0 {
            if value[0] == bb(b'F') {
                //if !simple.is_null() {
                *simple = 0; /* not a simple FITS file */
            // };
            } else if value[0] != bb(b'T') {
                *status = BAD_SIMPLE;
                return *status;
            };
        } else {
            int_snprintf!(
                &mut message,
                FLEN_ERRMSG,
                "First keyword of the file is not SIMPLE: {}",
                slice_to_str!(&name),
            );
            ffpmsg_slice(&message);

            *status = NO_SIMPLE;
            return *status;
        };
    } else {
        /* not beginning of the file, so presumably an IMAGE extension */
        /* or it could be a compressed image in a binary table */

        if strcmp_safe(&name, cs!(c"XTENSION")) == 0 {
            if ffc2s(&value, &mut xtension, status) > 0 {
                /* get the value string */
                ffpmsg_str("Bad value string for XTENSION keyword:");
                ffpmsg_slice(&value);
                return *status;
            }

            /* allow the quoted string value to begin in any column and */
            /* allow any number of trailing blanks before the closing quote */
            /* first char must be a quote */
            if (value[0] != bb(b'\''))
                || (strcmp_safe(&xtension, cs!(c"IMAGE")) != 0
                    && strcmp_safe(&xtension, cs!(c"IUEIMAGE")) != 0)
            {
                /* unknown type of extension; press on anyway */
                unknown = 1;
                int_snprintf!(
                    &mut message,
                    FLEN_ERRMSG,
                    "This is not an IMAGE extension: {}",
                    slice_to_str!(&value),
                );
                ffpmsg_slice(&message);
            };
        } else {
            /* error: 1st keyword of extension != XTENSION */
            int_snprintf!(
                &mut message,
                FLEN_ERRMSG,
                "First keyword of the extension is not XTENSION: {}",
                slice_to_str!(&name),
            );
            ffpmsg_slice(&message);

            *status = NO_XTENSION;
            return *status;
        };
    }

    if unknown != 0 && fptr.Fptr.compressimg != 0 {
        /* this is a compressed image, so read ZBITPIX, ZNAXIS keywords */

        unknown = 0; /* reset flag */
        ffxmsg_safer(3, Some(&mut message)); /* clear previous spurious error message */
        // if !bitpix.is_null() {
        ffgidt_safe(fptr, bitpix, status); /* get bitpix value */
        if *status > 0 {
            ffpmsg_str("Error reading BITPIX value of compressed image");
            return *status;
        };
        //}

        if let Some(naxis) = naxis {
            ffgidm_safe(fptr, naxis, status); /* get NAXIS value */
            if *status > 0 {
                ffpmsg_str("Error reading NAXIS value of compressed image");
                return *status;
            };
        }

        //if !naxes.is_null() {
        ffgiszll_safe(fptr, maxdim, naxes, status); /* get NAXISn value */
        if *status > 0 {
            ffpmsg_str("Error reading NAXISn values of compressed image");
            return *status;
        };
        //}
        nextkey = 9; /* skip required table keywords in the following search */
    } else {
        /*----------------------------------------------------------------*/
        /*  Get 2nd keyword;  test whether it is BITPIX with legal value  */
        /*----------------------------------------------------------------*/

        /* BITPIX = 2nd keyword */
        ffgkyn_safe(fptr, 2, &mut name, &mut value, Some(&mut comm), status);
        if strcmp_safe(&name, cs!(c"BITPIX")) != 0 {
            int_snprintf!(
                &mut message,
                FLEN_ERRMSG,
                "Second keyword of the extension is not BITPIX: {}",
                slice_to_str!(&name),
            );
            ffpmsg_slice(&message);

            *status = NO_BITPIX;
            return *status;
        }
        if ffc2ii(&value, &mut longbitpix, status) > 0 {
            int_snprintf!(
                &mut message,
                FLEN_ERRMSG,
                "Value of BITPIX keyword is not an integer: {}",
                slice_to_str!(&value),
            );
            ffpmsg_slice(&message);

            *status = BAD_BITPIX;
            return *status;
        } else if longbitpix != BYTE_IMG as c_long
            && longbitpix != SHORT_IMG as c_long
            && longbitpix != LONG_IMG as c_long
            && longbitpix != LONGLONG_IMG as c_long
            && longbitpix != FLOAT_IMG as c_long
            && longbitpix != DOUBLE_IMG as c_long
        {
            int_snprintf!(
                &mut message,
                FLEN_ERRMSG,
                "Illegal value for BITPIX keyword: {}",
                slice_to_str!(&value),
            );
            ffpmsg_slice(&message);

            *status = BAD_BITPIX;
            return *status;
        }
        // if !bitpix.is_null() {
        *bitpix = longbitpix as c_int; /* do explicit type conversion */
        //}

        /*---------------------------------------------------------------*/
        /*  Get 3rd keyword;  test whether it is NAXIS with legal value  */
        /*---------------------------------------------------------------*/
        ffgtkn(fptr, 3, cs!(c"NAXIS"), &mut longnaxis, status);
        if *status == BAD_ORDER {
            *status = NO_NAXIS;
            return *status;
        } else if *status == NOT_POS_INT || longnaxis > 999 {
            int_snprintf!(
                &mut message,
                FLEN_ERRMSG,
                "NAXIS = {} is illegal",
                longnaxis,
            );
            ffpmsg_slice(&message);

            *status = BAD_NAXIS;
            return *status;
        } else if let Some(naxis) = naxis {
            *naxis = longnaxis as c_int; /* do explicit type conversion */
        }

        /*---------------------------------------------------------*/
        /*  Get the next NAXISn keywords and test for legal values */
        /*---------------------------------------------------------*/

        ii = 0;
        nextkey = 4;
        while ii < longnaxis {
            ffkeyn_safe(cs!(c"NAXIS"), ii as c_int + 1, &mut keyword, status);
            ffgtknjj(fptr, 4 + ii as c_int, &keyword, &mut axislen, status);
            if *status == BAD_ORDER {
                *status = NO_NAXES;
                return *status;
            } else if *status == NOT_POS_INT {
                *status = BAD_NAXES;
                return *status;
            } else if (ii as c_int) < maxdim {
                //if !naxes.is_null() {
                naxes[ii as usize] = axislen;
                //}
            };
            {
                ii += 1;
                nextkey += 1
            }
        }
    }

    /*---------------------------------------------------------*/
    /*  now look for other keywords of interest:               */
    /*  BSCALE, BZERO, BLANK, PCOUNT, GCOUNT, EXTEND, and END  */
    /*---------------------------------------------------------*/

    /*  initialize default values in case keyword is not present */

    //if !bscale.is_null() {
    *bscale = 1.0;
    //}
    // if !bzero.is_null() {
    *bzero = 0.0;
    //}
    // if !pcount.is_null() {
    *pcount = 0;
    //}
    //if !gcount.is_null() {
    *gcount = 1;
    //}
    //if !extend.is_null() {
    *extend = 0;
    //}
    //if !blank.is_null() {
    *blank = NULL_UNDEFINED as LONGLONG; /* no default null value for BITPIX=8,16,32 */
    //}

    *nspace = 0;
    let mut found_end = 0;
    let tstatus = *status;
    let mut namelen = 0;

    while found_end == 0 {
        /* get next keyword */
        /* don't use ffgkyn here because it trys to parse the card to read */
        /* the value string, thus failing to read the file just because of */
        /* minor syntax errors in optional keywords.                       */

        if ffgrec_safe(fptr, nextkey, Some(&mut card), status) > 0 {
            /* get the 80-byte card */

            if *status == KEY_OUT_BOUNDS {
                found_end = 1; /* simply hit the end of the header */
                *status = tstatus; /* reset error status */
            } else {
                ffpmsg_str("Failed to find the END keyword in header (ffgphd).");
            };
        } else {
            /* got the next keyword without error */

            /* get the keyword name */
            ffgknm_safe(&card, &mut name, &mut namelen, status);

            if fftrec_safe(&name, status) > 0 {
                /* test keyword name; catches no END */
                int_snprintf!(
                    &mut message,
                    FLEN_ERRMSG,
                    "Name of keyword no. {} contains illegal character(s): {}",
                    nextkey,
                    slice_to_str!(&name),
                );
                ffpmsg_slice(&message);

                if nextkey % 36 == 0 {
                    /* test if at beginning of 36-card record */
                    ffpmsg_str("  (This may indicate a missing END keyword).");
                };
            }

            if strcmp_safe(&name, cs!(c"BSCALE")) == 0 {
                /* todo && !bscale.is_null() */
                *nspace = 0; /* reset count of blank keywords */

                /* parse value and comment */
                ffpsvc_safe(&card, &mut value, Some(&mut comm), status);

                if ffc2dd(&value, bscale, status) > 0 {
                    /* convert to double */

                    /* reset error status and continue, but still issue warning */
                    *status = tstatus;
                    *bscale = 1.0;
                    int_snprintf!(
                        &mut message,
                        FLEN_ERRMSG,
                        "Error reading BSCALE keyword value as a double: {}",
                        slice_to_str!(&value),
                    );
                    ffpmsg_slice(&message);
                };
            } else if strcmp_safe(&name, cs!(c"BZERO")) == 0 {
                /* todo && !bzero.is_null() */
                *nspace = 0; /* reset count of blank keywords */

                /* parse value and comment */
                ffpsvc_safe(&card, &mut value, Some(&mut comm), status);

                if ffc2dd(&value, bzero, status) > 0 {
                    /* convert to double */

                    /* reset error status and continue, but still issue warning */
                    *status = tstatus;
                    *bzero = 0.0;
                    int_snprintf!(
                        &mut message,
                        FLEN_ERRMSG,
                        "Error reading BZERO keyword value as a double: {}",
                        slice_to_str!(&value),
                    );
                    ffpmsg_slice(&message);
                };
            } else if strcmp_safe(&name, cs!(c"BLANK")) == 0 {
                /* todo && !blank.is_null() */
                *nspace = 0; /* reset count of blank keywords */

                /* parse value and comment */
                ffpsvc_safe(&card, &mut value, Some(&mut comm), status);

                if ffc2jj(&value, blank, status) > 0 {
                    /* convert to LONGLONG */

                    /* reset error status and continue, but still issue warning */
                    *status = tstatus;
                    *blank = NULL_UNDEFINED as LONGLONG;
                    int_snprintf!(
                        &mut message,
                        FLEN_ERRMSG,
                        "Error reading BLANK keyword value as an integer: {}",
                        slice_to_str!(&value),
                    );
                    ffpmsg_slice(&message);
                };
            } else if strcmp_safe(&name, cs!(c"PCOUNT")) == 0 {
                /* todo && !pcount.is_null() */
                *nspace = 0; /* reset count of blank keywords */

                /* parse value and comment */
                ffpsvc_safe(&card, &mut value, Some(&mut comm), status);

                if ffc2ii(&value, pcount, status) > 0 {
                    /* convert to long */
                    int_snprintf!(
                        &mut message,
                        FLEN_ERRMSG,
                        "Error reading PCOUNT keyword value as an integer: {}",
                        slice_to_str!(&value),
                    );
                    ffpmsg_slice(&message);
                };
            } else if strcmp_safe(&name, cs!(c"GCOUNT")) == 0 {
                /* todo && !gcount.is_null() */
                *nspace = 0; /* reset count of blank keywords */

                /* parse value and comment */
                ffpsvc_safe(&card, &mut value, Some(&mut comm), status);

                if ffc2ii(&value, gcount, status) > 0 {
                    /* convert to long */
                    int_snprintf!(
                        &mut message,
                        FLEN_ERRMSG,
                        "Error reading GCOUNT keyword value as an integer: {}",
                        slice_to_str!(&value),
                    );
                    ffpmsg_slice(&message);
                };
            } else if strcmp_safe(&name, cs!(c"EXTEND")) == 0
            /* todo && !extend.is_null() */
            {
                *nspace = 0; /* reset count of blank keywords */

                /* parse value and comment */
                ffpsvc_safe(&card, &mut value, Some(&mut comm), status);

                if ffc2ll(&value, extend, status) > 0 {
                    /* convert to logical */

                    /* reset error status and continue, but still issue warning */
                    *status = tstatus;
                    *extend = 0;
                    int_snprintf!(
                        &mut message,
                        FLEN_ERRMSG,
                        "Error reading EXTEND keyword value as a logical: {}",
                        slice_to_str!(&value),
                    );
                    ffpmsg_slice(&message);
                };
            } else if strcmp_safe(&name, cs!(c"END")) == 0 {
                found_end = 1;
            } else if card[0] == 0 {
                *nspace += 1; /* this is a blank card in the header */
            } else {
                /* reset count of blank keywords immediately before the END keyword to zero   */
                *nspace = 0;
            };
        }
        if *status > 0 {
            /* exit on error after writing error message */
            if fptr.Fptr.curhdu == 0 {
                ffpmsg_str("Failed to read the required primary array header keywords.");
            } else {
                ffpmsg_str("Failed to read the required image extension header keywords.");
            }
            return *status;
        };
        nextkey += 1
    }
    if unknown != 0 {
        *status = NOT_IMAGE;
    }
    *status
}

/*--------------------------------------------------------------------------*/
/// Get and Test TaBle;
/// Test that this is a legal ASCII or binary table and get some keyword values.
/// We assume that the calling routine has already tested the 1st keyword
/// of the extension to ensure that this is really a table extension.
pub(crate) fn ffgttb(
    fptr: &mut fitsfile,   /* I - FITS file pointer*/
    rowlen: &mut LONGLONG, /* O - length of a table row, in bytes */
    nrows: &mut LONGLONG,  /* O - number of rows in the table */
    pcount: &mut LONGLONG, /* O - value of PCOUNT keyword */
    tfields: &mut c_long,  /* O - number of fields in the table */
    status: &mut c_int,    /* IO - error status    */
) -> c_int {
    if *status > 0 {
        return *status;
    }

    if fftkyn(fptr, 2, cs!(c"BITPIX"), cs!(c"8"), status) == BAD_ORDER {
        /* 2nd keyword */
        *status = NO_BITPIX; /* keyword not BITPIX */
        return *status;
    } else if *status == NOT_POS_INT {
        *status = BAD_BITPIX; /* value != 8 */
        return *status;
    }

    if fftkyn(fptr, 3, cs!(c"NAXIS"), cs!(c"2"), status) == BAD_ORDER {
        /* 3rd keyword */
        *status = NO_NAXIS; /* keyword not NAXIS */
        return *status;
    } else if *status == NOT_POS_INT {
        *status = BAD_NAXIS; /* value != 2 */
        return *status;
    }

    if ffgtknjj(fptr, 4, cs!(c"NAXIS1"), rowlen, status) == BAD_ORDER {
        /* 4th key */
        *status = NO_NAXES; /* keyword not NAXIS1 */
        return *status;
    } else if *status == NOT_POS_INT {
        *status = BAD_NAXES; /* bad NAXIS1 value */
        return *status;
    }

    if ffgtknjj(fptr, 5, cs!(c"NAXIS2"), nrows, status) == BAD_ORDER {
        /* 5th key */
        *status = NO_NAXES; /* keyword not NAXIS2 */
        return *status;
    } else if *status == NOT_POS_INT {
        *status = BAD_NAXES; /* bad NAXIS2 value */
        return *status;
    }

    if ffgtknjj(fptr, 6, cs!(c"PCOUNT"), pcount, status) == BAD_ORDER {
        /* 6th key */
        *status = NO_PCOUNT; /* keyword not PCOUNT */
        return *status;
    } else if *status == NOT_POS_INT {
        *status = BAD_PCOUNT; /* bad PCOUNT value */
        return *status;
    }

    if fftkyn(fptr, 7, cs!(c"GCOUNT"), cs!(c"1"), status) == BAD_ORDER {
        /* 7th keyword */
        *status = NO_GCOUNT; /* keyword not GCOUNT */
        return *status;
    } else if *status == NOT_POS_INT {
        *status = BAD_GCOUNT; /* value != 1 */
        return *status;
    }

    if ffgtkn(fptr, 8, cs!(c"TFIELDS"), tfields, status) == BAD_ORDER {
        /* 8th key*/
        *status = NO_TFIELDS; /* keyword not TFIELDS */
        return *status;
    } else if *status == NOT_POS_INT || *tfields > 999 {
        *status = BAD_TFIELDS; /* bad TFIELDS value */
        return *status;
    }

    if *status > 0 {
        ffpmsg_str("Error reading required keywords in the table header (FTGTTB).");
    }
    *status
}

/*--------------------------------------------------------------------------*/
/// Test that keyword number NUMKEY has the expected name and get the
/// integer value of the keyword.  Return an error if the keyword
/// name does not match the input name, or if the value of the
/// keyword is not a positive integer.
pub(crate) fn ffgtkn(
    fptr: &mut fitsfile, /* I - FITS file pointer              */
    numkey: c_int,       /* I - number of the keyword to read  */
    name: &[c_char],     /* I - expected name of the keyword   */
    value: &mut c_long,  /* O - integer value of the keyword   */
    status: &mut c_int,  /* IO - error status                  */
) -> c_int {
    let mut keyname: [c_char; FLEN_KEYWORD] = [0; FLEN_KEYWORD]; /* incorrect keyword name */
    let mut valuestring: [c_char; FLEN_VALUE] = [0; FLEN_VALUE]; /* convert to integer */
    let mut comm: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];
    let mut message: [c_char; FLEN_ERRMSG] = [0; FLEN_ERRMSG];

    if *status > 0 {
        return *status;
    }

    keyname[0] = 0;
    valuestring[0] = 0;
    if ffgkyn_safe(
        fptr,
        numkey,
        &mut keyname,
        &mut valuestring,
        Some(&mut comm),
        status,
    ) <= 0
    {
        if strncmp_safe(&keyname, name, FLEN_KEYWORD) != 0 {
            *status = BAD_ORDER; /* incorrect keyword name */
        } else {
            /* convert to integer */
            ffc2ii(&valuestring, value, status);
            if *status > 0 || *value < 0 {
                *status = NOT_POS_INT;
            };
        }
        if *status > 0 {
            int_snprintf!(
                &mut message,
                FLEN_ERRMSG,
                "ffgtkn found unexpected keyword or value for keyword no. {}.",
                numkey,
            );
            ffpmsg_slice(&message);

            int_snprintf!(
                &mut message,
                FLEN_ERRMSG,
                " Expected positive integer keyword {}, but instead",
                slice_to_str!(&name),
            );
            ffpmsg_slice(&message);

            int_snprintf!(
                &mut message,
                FLEN_ERRMSG,
                " found keyword {} with value {}",
                slice_to_str!(&keyname),
                slice_to_str!(&valuestring),
            );
            ffpmsg_slice(&message);
        };
    }
    *status
}

/*--------------------------------------------------------------------------*/
/// Test that keyword number NUMKEY has the expected name and get the
/// integer value of the keyword.  Return an error if the keyword
/// name does not match the input name, or if the value of the
/// keyword is not a positive integer.
pub(crate) fn ffgtknjj(
    fptr: &mut fitsfile,  /* I - FITS file pointer              */
    numkey: c_int,        /* I - number of the keyword to read  */
    name: &[c_char],      /* I - expected name of the keyword   */
    value: &mut LONGLONG, /* O - integer value of the keyword   */
    status: &mut c_int,   /* IO - error status                  */
) -> c_int {
    let mut keyname: [c_char; FLEN_KEYWORD] = [0; FLEN_KEYWORD];
    let mut valuestring: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut comm: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];
    let mut message: [c_char; FLEN_ERRMSG] = [0; FLEN_ERRMSG];

    if *status > 0 {
        return *status;
    }

    keyname[0] = 0;
    valuestring[0] = 0;

    if ffgkyn_safe(
        fptr,
        numkey,
        &mut keyname,
        &mut valuestring,
        Some(&mut comm),
        status,
    ) <= 0
    {
        if strcmp_safe(&keyname, name) != 0 {
            *status = BAD_ORDER; /* incorrect keyword name */
        } else {
            /* convert to integer */
            ffc2jj(&valuestring, value, status);
            if *status > 0 || *value < 0 {
                *status = NOT_POS_INT;
            };
        }

        if *status > 0 {
            int_snprintf!(
                &mut message,
                FLEN_ERRMSG,
                "ffgtknjj found unexpected keyword or value for keyword no. {}.",
                numkey,
            );
            ffpmsg_slice(&message);

            int_snprintf!(
                &mut message,
                FLEN_ERRMSG,
                " Expected positive integer keyword {}, but instead",
                slice_to_str!(&name),
            );
            ffpmsg_slice(&message);

            int_snprintf!(
                &mut message,
                FLEN_ERRMSG,
                " found keyword {} with value {}",
                slice_to_str!(&keyname),
                slice_to_str!(&valuestring),
            );
            ffpmsg_slice(&message);
        };
    }
    *status
}

/*--------------------------------------------------------------------------*/
/// Test that keyword number NUMKEY has the expected name and the
/// expected value string.
pub(crate) fn fftkyn(
    fptr: &mut fitsfile, /* I - FITS file pointer              */
    numkey: c_int,       /* I - number of the keyword to read  */
    name: &[c_char],     /* I - expected name of the keyword   */
    value: &[c_char],    /* I - expected value of the keyword  */
    status: &mut c_int,  /* IO - error status                  */
) -> c_int {
    let mut keyname: [c_char; FLEN_KEYWORD] = [0; FLEN_KEYWORD];
    let mut valuestring: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut comm: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];
    let mut message: [c_char; FLEN_ERRMSG] = [0; FLEN_ERRMSG];

    if *status > 0 {
        return *status;
    }
    keyname[0] = 0;
    valuestring[0] = 0;

    if ffgkyn_safe(
        fptr,
        numkey,
        &mut keyname,
        &mut valuestring,
        Some(&mut comm),
        status,
    ) <= 0
    {
        if strcmp_safe(&keyname, name) > 0 {
            *status = BAD_ORDER; /* incorrect keyword name */
        }
        if strcmp_safe(value, &valuestring) > 0 {
            *status = NOT_POS_INT; /* incorrect keyword value */
        }
    }

    if *status > 0 {
        int_snprintf!(
            &mut message,
            FLEN_ERRMSG,
            "fftkyn found unexpected keyword or value for keyword no. {}.",
            numkey,
        );
        ffpmsg_slice(&message);

        int_snprintf!(
            &mut message,
            FLEN_ERRMSG,
            " Expected keyword {} with value {}, but",
            slice_to_str!(&name),
            slice_to_str!(&value),
        );
        ffpmsg_slice(&message);

        int_snprintf!(
            &mut message,
            FLEN_ERRMSG,
            " found keyword {} with value {}",
            slice_to_str!(&keyname),
            slice_to_str!(&valuestring),
        );
        ffpmsg_slice(&message);
    }

    *status
}

/*--------------------------------------------------------------------------*/
/// Read header keywords into a long string of chars.  This routine allocates
/// memory for the string, so the calling routine must eventually free the
/// memory when it is not needed any more.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffh2st(
    fptr: *mut fitsfile,      /* I - FITS file pointer           */
    header: *mut *mut c_char, /* O - returned header string      */
    status: *mut c_int,       /* IO - error status               */
) -> c_int {
    // FFI WRAPPER
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let header = header.as_mut().expect(NULL_MSG);

        ffh2st_safe(fptr, header, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Read header keywords into a long string of chars.  This routine allocates
/// memory for the string, so the calling routine must eventually free the
/// memory when it is not needed any more.
pub(crate) fn ffh2st_safe(
    fptr: &mut fitsfile,      /* I - FITS file pointer           */
    header: &mut *mut c_char, /* O - returned header string      */
    status: &mut c_int,       /* IO - error status               */
) -> c_int {
    let mut nkeys: c_int = 0;
    let mut nrec: c_long = 0;
    let mut headstart: LONGLONG = 0;

    if *status > 0 {
        return *status;
    }

    /* get number of keywords in the header (doesn't include END) */
    if ffghsp_safe(fptr, Some(&mut nkeys), None, status) > 0 {
        return *status;
    }

    nrec = nkeys as c_long / 36 + 1;

    /* allocate memory for all the keywords (multiple of 2880 bytes) */
    // HEAP ALLOCATION
    let mut hdr: Vec<c_char> = Vec::new();
    if hdr
        .try_reserve_exact((nrec as LONGLONG * IOBUFLEN + 1) as usize)
        .is_err()
    {
        *status = MEMORY_ALLOCATION;
        ffpmsg_str("failed to allocate memory to hold all the header keywords");
        return *status;
    } else {
        hdr.resize((nrec as LONGLONG * IOBUFLEN + 1) as usize, 0);
    }

    ffghadll_safe(fptr, Some(&mut headstart), None, None, status); /* get header address */
    ffmbyt_safe(fptr, headstart, REPORT_EOF, status); /* move to header */
    ffgbyt(
        fptr,
        nrec as LONGLONG * IOBUFLEN,
        cast_slice_mut(&mut hdr),
        status,
    ); /* copy header */
    hdr[(nrec as LONGLONG * IOBUFLEN) as usize] = 0;

    let (header_ptr, l, c) = vec_into_raw_parts(hdr);
    ALLOCATIONS
        .lock()
        .unwrap()
        .insert(header_ptr as usize, (l, c));

    *header = header_ptr;

    *status
}

/*--------------------------------------------------------------------------*/
/// Read header keywords into a long string of chars.  This routine allocates
/// memory for the string, so the calling routine must eventually free the
/// memory when it is not needed any more.  If exclude_comm is TRUE, then all
/// the COMMENT, HISTORY, and <blank> keywords will be excluded from the output
/// string of keywords.  Any other list of keywords to be excluded may be
/// specified with the exclist parameter.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffhdr2str(
    fptr: *mut fitsfile,           /* I - FITS file pointer                    */
    exclude_comm: c_int,           /* I - if TRUE, exclude commentary keywords */
    exclist: *const *const c_char, /* I - list of excluded keyword names       */
    nexc: c_int,                   /* I - number of names in exclist           */
    header: *mut *mut c_char,      /* O - returned header string               */
    nkeys: *mut c_int,             /* O - returned number of 80-char keywords  */
    status: *mut c_int,            /* IO - error status                        */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let nkeys = nkeys.as_mut().expect(NULL_MSG);
        let header = header.as_mut().expect(NULL_MSG);
        let exclist = slice::from_raw_parts(exclist, nexc as usize);

        let mut v_exclist = Vec::new();

        for item in exclist {
            let exclist_item: &[c_char] = cast_slice(CStr::from_ptr(*item).to_bytes_with_nul());
            v_exclist.push(exclist_item);
        }

        ffhdr2str_safe(fptr, exclude_comm, &v_exclist, nexc, header, nkeys, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Read header keywords into a long string of chars.  This routine allocates
/// memory for the string, so the calling routine must eventually free the
/// memory when it is not needed any more.  If exclude_comm is TRUE, then all
/// the COMMENT, HISTORY, and <blank> keywords will be excluded from the output
/// string of keywords.  Any other list of keywords to be excluded may be
/// specified with the exclist parameter.
pub(crate) fn ffhdr2str_safe(
    fptr: &mut fitsfile,      /* I - FITS file pointer                    */
    exclude_comm: c_int,      /* I - if TRUE, exclude commentary keywords */
    exclist: &[&[c_char]],    /* I - list of excluded keyword names       */
    nexc: c_int,              /* I - number of names in exclist           */
    header: &mut *mut c_char, /* O - returned header string               */
    nkeys: &mut c_int,        /* O - returned number of 80-char keywords  */
    status: &mut c_int,       /* IO - error status                        */
) -> c_int {
    let mut casesn: c_int = 0;
    let mut matched: c_int = 0;
    let mut exact: c_int = 0;
    let mut totkeys: c_int = 0;

    let mut keybuf: [c_char; 162] = [0; 162];
    let mut keyname: [c_char; FLEN_KEYWORD] = [0; FLEN_KEYWORD];

    let nexc = nexc as usize;

    *nkeys = 0;

    if *status > 0 {
        return *status;
    }

    /* get number of keywords in the header (doesn't include END) */
    if ffghsp_safe(fptr, Some(&mut totkeys), None, status) > 0 {
        return *status;
    }

    /* allocate memory for all the keywords */
    /* (will reallocate it later to minimize the memory size) */
    // HEAP ALLOCATION
    let mut hdr: Vec<c_char> = Vec::new();
    if hdr
        .try_reserve_exact(((totkeys + 1) * 80 + 1) as usize)
        .is_err()
    {
        *status = MEMORY_ALLOCATION;
        ffpmsg_str("failed to allocate memory to hold all the header keywords");
        return *status;
    } else {
        hdr.resize(((totkeys + 1) * 80 + 1) as usize, 0);
    }

    let mut hi = 0; //Index into *header / h
    casesn = FALSE as c_int;

    /* read every keyword */
    for ii in 1..=(totkeys as usize) {
        ffgrec_safe(fptr, ii as c_int, Some(&mut keybuf), status);
        /* pad record with blanks so that it is at least 80 chars long */
        strcat_safe(
            &mut keybuf,
            cs!(
                c"                                                                                "
            ),
        );

        keyname[0] = 0;
        strncat_safe(&mut keyname, &keybuf, 8); /* copy the keyword name */

        if exclude_comm != 0
            && (FSTRCMP(cs!(c"COMMENT "), &keyname) == 0
                || FSTRCMP(cs!(c"HISTORY "), &keyname) == 0
                || FSTRCMP(cs!(c"        "), &keyname) == 0)
        {
            continue; /* skip this commentary keyword */
        }

        /* does keyword match any names in the exclusion list? */
        let mut jj: usize = 0;

        while jj < nexc {
            ffcmps_safe(exclist[jj], &keyname, casesn, &mut matched, &mut exact);
            if matched != 0 {
                break;
            }
            jj += 1;
        }

        if jj == nexc {
            /* not in exclusion list, add this keyword to the string */
            strcpy_safe(&mut hdr[hi..], &keybuf);
            hi += 80;
            (*nkeys) += 1;
        }
    }

    /* add the END keyword */
    strcpy_safe(
        &mut hdr[hi..],
        cs!(c"END                                                                             "),
    );
    hi += 80;
    (*nkeys) += 1;

    hdr[hi] = 0; /* terminate the header string */

    /* minimize the allocated memory */
    hdr.resize(((*nkeys * 80) + 1) as usize, 0);
    hdr.shrink_to_fit();

    let (header_ptr, l, c) = vec_into_raw_parts(hdr);
    ALLOCATIONS
        .lock()
        .unwrap()
        .insert(header_ptr as usize, (l, c));

    *header = header_ptr;

    *status
}

/*--------------------------------------------------------------------------*/
/// Same as ffhdr2str, except that if the input HDU is a tile compressed image
/// (stored in a binary table) then it will first convert that header back
/// to that of a normal uncompressed FITS image before concatenating the header
/// keyword records.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffcnvthdr2str(
    fptr: *mut fitsfile,           /* I - FITS file pointer                    */
    exclude_comm: c_int,           /* I - if TRUE, exclude commentary keywords */
    exclist: *const *const c_char, /* I - list of excluded keyword names       */
    nexc: c_int,                   /* I - number of names in exclist           */
    header: *mut *mut c_char,      /* O - returned header string               */
    nkeys: *mut c_int,             /* O - returned number of 80-char keywords  */
    status: *mut c_int,            /* IO - error status                        */
) -> c_int {
    unsafe {
        let mut tempfptr: fitsfile;

        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let nkeys = nkeys.as_mut().expect(NULL_MSG);
        let header = header.as_mut().expect(NULL_MSG);
        let exclist = slice::from_raw_parts(exclist, nexc as usize);

        let mut v_exclist = Vec::new();

        for item in exclist {
            let exclist_item: &[c_char] = cast_slice(CStr::from_ptr(*item).to_bytes_with_nul());
            v_exclist.push(exclist_item);
        }

        if *status > 0 {
            return *status;
        }

        if fits_is_compressed_image_safe(fptr, status) != 0 {
            /* this is a tile compressed image, so need to make an uncompressed */
            /* copy of the image header in memory before concatenating the keywords */
            todo!();
        /*
        if (fits_create_file(&tempfptr, "mem://", status) > 0) {
            return(*status);
        }


        if (fits_img_decompress_header(fptr, tempfptr, status) > 0) {
         fits_delete_file(tempfptr, status);
         return(*status);
        }

        ffhdr2str(tempfptr, exclude_comm, exclist, nexc, header, nkeys, status);
        fits_close_file(tempfptr, status);
        */
        } else {
            ffhdr2str_safe(fptr, exclude_comm, &v_exclist, nexc, header, nkeys, status);
        }

        *status
    }
}
