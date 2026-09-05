//! Routines that create, copy, insert and delete whole HDUs.
//!
//! Inserting or deleting an HDU shifts every following HDU in the file, which
//! is why these routines rewrite the file's HDU offsets as they go. Creating an
//! image or table extension here writes only the required keywords; the data
//! unit is filled in by the column and image I/O routines.
//!
//! Ported from CFITSIO's `edithdu.c`, written by William Pence at the High
//! Energy Astrophysics Science Archive Research Center (HEASARC), NASA Goddard
//! Space Flight Center.
#![warn(missing_docs)]

use crate::c_types::{FILE, c_char, c_int, c_long};
use crate::helpers::cfile::CFile;
use crate::{
    BL, TKeywords,
    buffers::*,
    cs,
    editcol::ffcprw_safe,
    fitscore::{
        ffbnfm_safe, ffcmsg_safe, ffcrhd_safe, ffdblk, ffgabc_safe, ffgext, ffghadll_safe,
        ffghdn_safe, ffgidm_safe, ffhdef_safe, ffiblk, ffkeyn_safe, ffmahd_safe, ffpdfl,
        ffpmsg_slice, ffpmsg_str, ffrdef_safe, ffrhdu_safe, ffwend,
    },
    fitsio::*,
    fitsio2::*,
    getkey::{ffgcrd_safe, ffghsp_safe, ffgkyj_safe, ffgrec_safe},
    int_snprintf,
    modkey::{ffdkey_safe, ffikyj_safe, ffukyj_safe},
    nullable_slice_cstr,
    putkey::{
        ffcrim_safe, ffcrimll_safe, ffcrtb_safe, ffphbn_safe, ffphpr_safe, ffphprll_safe,
        ffphtb_safe, ffpkyj_safe, ffpkyl_safe, ffpkys_safe, ffprec_safe,
    },
    wrappers::{strcpy_safe, strncat_safe},
};
use bytemuck::{cast_slice, cast_slice_mut};
use core::slice;
use std::io::Write;

use core::{cmp, ffi::CStr};

/// Copy the current HDU from `infptr` and append it to the end of the file
/// associated with `outfptr`.
///
/// Space may be reserved in the output header for `morekeys` additional
/// keywords.
///
/// # Parameters
///
/// * `infptr`   — (I) FITS file pointer to input file
/// * `outfptr`  — (I) FITS file pointer to output file
/// * `morekeys` — (I) reserve space in output header
/// * `status`   — (IO) error status
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffcopy(
    infptr: *mut fitsfile,
    outfptr: *mut fitsfile,
    morekeys: c_int,
    status: *mut c_int,
) -> c_int {
    // FFI WRAPPER
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);

        if *status > 0 {
            return *status;
        }

        if core::ptr::eq(infptr, outfptr) {
            *status = SAME_FILE;
            return *status;
        }

        let infptr = infptr.as_mut().expect(NULL_MSG);
        let outfptr = outfptr.as_mut().expect(NULL_MSG);

        ffcopy_safe(infptr, outfptr, morekeys, status)
    }
}

/// Copy the current HDU from `infptr` and append it to the end of the file
/// associated with `outfptr`.
///
/// Space may be reserved in the output header for `morekeys` additional
/// keywords.
///
/// # Parameters
///
/// * `infptr`   — (I) FITS file pointer to input file
/// * `outfptr`  — (I) FITS file pointer to output file
/// * `morekeys` — (I) reserve space in output header
/// * `status`   — (IO) error status
pub fn ffcopy_safe(
    infptr: &mut fitsfile,
    outfptr: &mut fitsfile,
    morekeys: c_int,
    status: &mut c_int,
) -> c_int {
    let mut nspace: c_int = 0;

    if *status > 0 {
        return *status;
    }

    if ffcphd_safe(infptr, outfptr, status) > 0 {
        /* copy the header keywords */
        return *status;
    }

    if morekeys > 0 {
        ffhdef_safe(outfptr, morekeys, status); /* reserve space for more keywords */
    } else {
        if ffghsp_safe(infptr, None, Some(&mut nspace), status) > 0 {
            /* get existing space */
            return *status;
        }
        if nspace > 0 {
            ffhdef_safe(outfptr, nspace, status); /* preserve same amount of space */
            if nspace >= 35 {
                /* There is at least 1 full empty FITS block in the header. */
                /* Physically write the END keyword at the beginning of the */
                /* last block to preserve this extra space now rather than */
                /* later.  This is needed by the stream: driver which cannot */
                /* seek back to the header to write the END keyword later. */

                ffwend(outfptr, status);
            }
        }
    }

    ffcpdt_safe(infptr, outfptr, status); /* now copy the data unit */

    *status
}

/// Copy all or part of the input file to the output file.
///
/// # Parameters
///
/// * `infptr`    — (I) FITS file pointer to input file
/// * `outfptr`   — (I) FITS file pointer to output file
/// * `previous`  — (I) copy any previous HDUs?
/// * `current`   — (I) copy the current HDU?
/// * `following` — (I) copy any following HDUs?
/// * `status`    — (IO) error status
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffcpfl(
    infptr: *mut fitsfile,
    outfptr: *mut fitsfile,
    previous: c_int,
    current: c_int,
    following: c_int,
    status: *mut c_int,
) -> c_int {
    // FFI WRAPPER
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);

        if *status > 0 {
            return *status;
        }

        if core::ptr::eq(infptr, outfptr) {
            *status = SAME_FILE;
            return *status;
        }

        let infptr = infptr.as_mut().expect(NULL_MSG);
        let outfptr = outfptr.as_mut().expect(NULL_MSG);

        ffcpfl_safe(infptr, outfptr, previous, current, following, status)
    }
}

/// Copy all or part of the input file to the output file.
///
/// # Parameters
///
/// * `infptr`    — (I) FITS file pointer to input file
/// * `outfptr`   — (I) FITS file pointer to output file
/// * `previous`  — (I) copy any previous HDUs?
/// * `current`   — (I) copy the current HDU?
/// * `following` — (I) copy any following HDUs?
/// * `status`    — (IO) error status
pub fn ffcpfl_safe(
    infptr: &mut fitsfile,
    outfptr: &mut fitsfile,
    previous: c_int,
    current: c_int,
    following: c_int,
    status: &mut c_int,
) -> c_int {
    let mut hdunum = 0;

    if *status > 0 {
        return *status;
    }

    ffghdn_safe(infptr, &mut hdunum);

    if previous != 0 {
        /* copy any previous HDUs */
        for ii in 1..(hdunum) {
            ffmahd_safe(infptr, ii, None, status);
            ffcopy_safe(infptr, outfptr, 0, status);
        }
    }

    if current != 0 && (*status <= 0) {
        /* copy current HDU */
        ffmahd_safe(infptr, hdunum, None, status);
        ffcopy_safe(infptr, outfptr, 0, status);
    }

    if following != 0 && (*status <= 0) {
        /* copy any remaining HDUs */
        let mut ii = hdunum + 1;
        loop {
            if ffmahd_safe(infptr, ii, None, status) != 0 {
                /* reset expected end of file status */
                if *status == END_OF_FILE {
                    *status = 0;
                }
                break;
            }

            if ffcopy_safe(infptr, outfptr, 0, status) != 0 {
                break; /* quit on unexpected error */
            }

            ii += 1;
        }
    }

    ffmahd_safe(infptr, hdunum, None, status); /* restore initial position */
    *status
}

/// Copy the header keywords from `infptr` to `outfptr`.
///
/// # Parameters
///
/// * `infptr`  — (I) FITS file pointer to input file
/// * `outfptr` — (I) FITS file pointer to output file
/// * `status`  — (IO) error status
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffcphd(
    infptr: *mut fitsfile,
    outfptr: *mut fitsfile,
    status: *mut c_int,
) -> c_int {
    // FFI WRAPPER
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);

        if *status > 0 {
            return *status;
        }

        if core::ptr::eq(infptr, outfptr) {
            *status = SAME_FILE;
            return *status;
        }

        let infptr = infptr.as_mut().expect(NULL_MSG);
        let outfptr = outfptr.as_mut().expect(NULL_MSG);

        ffcphd_safe(infptr, outfptr, status)
    }
}

/// Copy the header keywords from `infptr` to `outfptr`.
///
/// # Parameters
///
/// * `infptr`  — (I) FITS file pointer to input file
/// * `outfptr` — (I) FITS file pointer to output file
/// * `status`  — (IO) error status
pub fn ffcphd_safe(infptr: &mut fitsfile, outfptr: &mut fitsfile, status: &mut c_int) -> c_int {
    let mut nkeys: c_int = 0;
    let mut inPrim: c_int = 0;
    let mut outPrim: c_int = 0;
    let mut naxis: c_long = 0;
    let naxes: [c_long; 1] = [0; 1];
    let mut comm: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];

    if *status > 0 {
        return *status;
    }

    /* set the input pointer to the correct HDU */
    if infptr.HDUposition != infptr.Fptr.curhdu {
        ffmahd_safe(infptr, (infptr.HDUposition) + 1, None, status);
    }

    if ffghsp_safe(infptr, Some(&mut nkeys), None, status) > 0 {
        /* get no. of keywords */
        return *status;
    }

    /* create a memory buffer to hold the header records */
    let mut tmpbuff2: Vec<c_char> = Vec::new();
    vec![0 as c_char; nkeys as usize * FLEN_CARD];

    if tmpbuff2
        .try_reserve_exact(nkeys as usize * FLEN_CARD)
        .is_err()
    {
        *status = MEMORY_ALLOCATION;
        return *status;
    } else {
        tmpbuff2.resize(nkeys as usize * FLEN_CARD, 0);
    }

    /* read all of the header records in the input HDU */
    for ii in 0..(nkeys as usize) {
        ffgrec_safe(
            infptr,
            ii as c_int + 1,
            Some(&mut tmpbuff2[(ii * FLEN_CARD)..]),
            status,
        );
    }

    if infptr.HDUposition == 0 {
        /* set flag if this is the Primary HDU */
        inPrim = 1;
    }

    /* if input is an image hdu, get the number of axes */
    naxis = -1; /* negative if HDU is a table */

    if infptr.Fptr.hdutype == IMAGE_HDU {
        ffgkyj_safe(infptr, cs!(c"NAXIS"), &mut naxis, None, status);
    }

    /* set the output pointer to the correct HDU */
    if outfptr.HDUposition != outfptr.Fptr.curhdu {
        ffmahd_safe(outfptr, (outfptr.HDUposition) + 1, None, status);
    }

    /* check if output header is empty; if not create new empty HDU */
    let headstart = outfptr.Fptr.get_headstart_as_slice()[outfptr.Fptr.curhdu as usize];
    if outfptr.Fptr.headend != headstart {
        ffcrhd_safe(outfptr, status);
    }
    if outfptr.HDUposition == 0 {
        if naxis < 0 {
            /* the input HDU is a table, so we have to create */
            /* a dummy Primary array before copying it to the output */
            ffcrim_safe(outfptr, 8, 0, &naxes, status);
            ffcrhd_safe(outfptr, status); /* create new empty HDU */
        } else {
            /* set flag that this is the Primary HDU */
            outPrim = 1;
        }
    }

    if *status > 0 {
        /* check for errors before proceeding */
        return *status;
    }

    if inPrim == 1 && outPrim == 0 {
        /* copying from primary array to image extension */
        strcpy_safe(&mut comm, cs!(c"IMAGE extension"));
        ffpkys_safe(
            outfptr,
            cs!(c"XTENSION"),
            cs!(c"IMAGE"),
            Some(&comm),
            status,
        );

        /* copy BITPIX through NAXISn keywords */
        for ii in 1..(3 + naxis as usize) {
            let card = &tmpbuff2[(ii * FLEN_CARD)..];
            ffprec_safe(outfptr, card, status);
        }

        strcpy_safe(&mut comm, cs!(c"number of random group parameters"));
        ffpkyj_safe(outfptr, cs!(c"PCOUNT"), 0, Some(&comm), status);

        strcpy_safe(&mut comm, cs!(c"number of random groups"));
        ffpkyj_safe(outfptr, cs!(c"GCOUNT"), 1, Some(&comm), status);

        /* copy remaining keywords, excluding EXTEND, and reference COMMENT keywords */
        for ii in (3 + naxis as usize)..(nkeys as usize) {
            let card = &tmpbuff2[(ii * FLEN_CARD)..];
            if FSTRNCMP(card, cs!(c"EXTEND  "), 8) > 0
                && FSTRNCMP(
                    card,
                    cs!(c"COMMENT   FITS (Flexible Image Transport System) format is"),
                    58,
                ) > 0
                && FSTRNCMP(
                    card,
                    cs!(c"COMMENT   and Astrophysics', volume 376, page 3"),
                    47,
                ) > 0
            {
                ffprec_safe(outfptr, card, status);
            }
        }
    } else if inPrim == 0 && outPrim == 1 {
        /* copying between image extension and primary array */
        strcpy_safe(&mut comm, cs!(c"file does conform to FITS standard"));
        ffpkyl_safe(outfptr, cs!(c"SIMPLE"), TRUE as c_int, Some(&comm), status);

        /* copy BITPIX through NAXISn keywords */
        for ii in 1..(3 + naxis as usize) {
            let card = &tmpbuff2[(ii * FLEN_CARD)..];
            ffprec_safe(outfptr, card, status);
        }

        /* add the EXTEND keyword */
        strcpy_safe(&mut comm, cs!(c"FITS dataset may contain extensions"));
        ffpkyl_safe(outfptr, cs!(c"EXTEND"), TRUE as c_int, Some(&comm), status);

        /* write standard block of self-documentating comments */
        ffprec_safe(
            outfptr,
            cs!(
                c"COMMENT   FITS (Flexible Image Transport System) format is defined in 'Astronomy"
            ),
            status,
        );
        ffprec_safe(
            outfptr,
            cs!(c"COMMENT   and Astrophysics', volume 376, page 359; bibcode: 2001A&A...376..359H"),
            status,
        );

        /* copy remaining keywords, excluding pcount, gcount */
        for ii in (3 + naxis as usize)..(nkeys as usize) {
            let card = &tmpbuff2[(ii * FLEN_CARD)..];
            if FSTRNCMP(card, cs!(c"PCOUNT  "), 8) > 0 && FSTRNCMP(card, cs!(c"GCOUNT  "), 8) > 0 {
                ffprec_safe(outfptr, card, status);
            }
        }
    } else {
        /* input and output HDUs are same type; simply copy all keywords */
        for ii in 0..(nkeys as usize) {
            let card = &tmpbuff2[(ii * FLEN_CARD)..];
            ffprec_safe(outfptr, card, status);
        }
    }

    *status
}

/// Copy the table structure from an existing table HDU, but only
/// Copy the structure of an open table to a new table, optionally copying a
/// limited range of rows.
///
/// All header keywords from the input table are copied directly, but `NAXIS2`
/// and `PCOUNT` are set to their correct values. Useful when a task will filter
/// rows before transferring them, so that a pristine output table with zero
/// rows is wanted to start with. `nrows` may be 0. The first row of a table is
/// row 1.
///
/// # Parameters
///
/// * `infptr`   — (I) FITS file pointer to input file
/// * `outfptr`  — (I) FITS file pointer to output file
/// * `firstrow` — (I) number of first row to copy (1 based)
/// * `nrows`    — (I) number of rows to copy
/// * `status`   — (IO) error status
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffcpht(
    infptr: *mut fitsfile,
    outfptr: *mut fitsfile,
    firstrow: LONGLONG,
    nrows: LONGLONG,
    status: *mut c_int,
) -> c_int {
    // FFI WRAPPER
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let infptr = infptr.as_mut().expect(NULL_MSG);
        let outfptr = outfptr.as_mut().expect(NULL_MSG);

        ffcpht_safe(infptr, outfptr, firstrow, nrows, status)
    }
}

/// Copy the table structure from an existing table HDU, but only
/// Copy the structure of an open table to a new table, optionally copying a
/// limited range of rows.
///
/// All header keywords from the input table are copied directly, but `NAXIS2`
/// and `PCOUNT` are set to their correct values. Useful when a task will filter
/// rows before transferring them, so that a pristine output table with zero
/// rows is wanted to start with. `nrows` may be 0. The first row of a table is
/// row 1.
///
/// # Parameters
///
/// * `infptr`   — (I) FITS file pointer to input file
/// * `outfptr`  — (I) FITS file pointer to output file
/// * `firstrow` — (I) number of first row to copy (1 based)
/// * `nrows`    — (I) number of rows to copy
/// * `status`   — (IO) error status
pub fn ffcpht_safe(
    infptr: &mut fitsfile,
    outfptr: &mut fitsfile,
    firstrow: LONGLONG,
    nrows: LONGLONG,
    status: &mut c_int,
) -> c_int {
    if *status > 0 {
        return *status;
    }

    /* Copy the header only */
    ffcphd_safe(infptr, outfptr, status);
    /* Note that we now have a copied header that describes the table,
    and that is the current header, but the original number of table
    rows and heap area sizes are still there. */

    /* Zero out the size-related keywords */
    if *status == 0 {
        ffukyj_safe(outfptr, cs!(c"NAXIS2"), 0, None, status); /* NAXIS2 = 0 */
        ffukyj_safe(outfptr, cs!(c"PCOUNT"), 0, None, status); /* PCOUNT = 0 */
        /* Update the internal structure variables within CFITSIO now
        that we have a valid table header */
        ffrdef_safe(outfptr, status);
    }

    /* OK now that we have a pristine HDU, copy the requested rows */
    if *status == 0 && nrows > 0 {
        ffcprw_safe(infptr, outfptr, firstrow, nrows, status);
    }

    *status
}

/// Copy the data unit, and not the header, from the current HDU of `infptr` to
/// the current HDU of `outfptr`.
///
/// This overwrites any data previously in the output HDU. A low level routine
/// used by [`ffcopy_safe`], but also useful to an application that wants to copy
/// the data from one file to another while modifying the header itself.
/// This will overwrite any data already in the outfptr CHDU.
///
/// # Parameters
///
/// * `infptr`  — (I) FITS file pointer to input file
/// * `outfptr` — (I) FITS file pointer to output file
/// * `status`  — (IO) error status
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffcpdt(
    infptr: *mut fitsfile,
    outfptr: *mut fitsfile,
    status: *mut c_int,
) -> c_int {
    // FFI WRAPPER
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);

        if *status > 0 {
            return *status;
        }

        if core::ptr::eq(infptr, outfptr) {
            *status = SAME_FILE;
            return *status;
        }

        let infptr = infptr.as_mut().expect(NULL_MSG);
        let outfptr = outfptr.as_mut().expect(NULL_MSG);

        ffcpdt_safe(infptr, outfptr, status)
    }
}

/// Copy the data unit, and not the header, from the current HDU of `infptr` to
/// the current HDU of `outfptr`.
///
/// This overwrites any data previously in the output HDU. A low level routine
/// used by [`ffcopy_safe`], but also useful to an application that wants to copy
/// the data from one file to another while modifying the header itself.
/// This will overwrite any data already in the outfptr CHDU.
///
/// # Parameters
///
/// * `infptr`  — (I) FITS file pointer to input file
/// * `outfptr` — (I) FITS file pointer to output file
/// * `status`  — (IO) error status
pub fn ffcpdt_safe(infptr: &mut fitsfile, outfptr: &mut fitsfile, status: &mut c_int) -> c_int {
    let mut nb: c_long = 0;
    let mut indatastart: LONGLONG = 0;
    let mut indataend: LONGLONG = 0;
    let mut outdatastart: LONGLONG = 0;
    let mut buffer: [c_char; IOBUFLEN as usize] = [0; IOBUFLEN as usize];

    if *status > 0 {
        return *status;
    }

    ffghadll_safe(
        infptr,
        None,
        Some(&mut indatastart),
        Some(&mut indataend),
        status,
    );

    ffghadll_safe(outfptr, None, Some(&mut outdatastart), None, status);

    /* Calculate the number of blocks to be copied  */
    nb = ((indataend - indatastart) / IOBUFLEN) as c_long;

    if nb > 0 {
        if core::ptr::eq(infptr.Fptr.as_ptr(), outfptr.Fptr.as_ptr()) {
            /* copying between 2 HDUs in the SAME file */
            for _ii in 0..(nb as usize) {
                ffmbyt_safe(infptr, indatastart, REPORT_EOF, status);
                ffgbyt(infptr, IOBUFLEN, cast_slice_mut(&mut buffer), status); /* read input block */

                ffmbyt_safe(outfptr, outdatastart, IGNORE_EOF, status);
                ffpbyt(outfptr, IOBUFLEN, cast_slice_mut(&mut buffer), status); /* write output block */

                indatastart += IOBUFLEN; /* move address */
                outdatastart += IOBUFLEN; /* move address */
            }
        } else {
            /* copying between HDUs in separate files */
            /* move to the initial copy position in each of the files */
            ffmbyt_safe(infptr, indatastart, REPORT_EOF, status);
            ffmbyt_safe(outfptr, outdatastart, IGNORE_EOF, status);

            for _ii in 0..(nb as usize) {
                ffgbyt(infptr, IOBUFLEN, cast_slice_mut(&mut buffer), status); /* read input block */
                ffpbyt(outfptr, IOBUFLEN, cast_slice_mut(&mut buffer), status); /* write output block */
            }
        }
    }
    *status
}

/// Write the data unit from the current HDU of `infptr` to the output file
/// stream.
///
/// # Parameters
///
/// * `infptr`    — (I) FITS file pointer to input file
/// * `outstream` — (I) stream to write HDU to
/// * `status`    — (IO) error status
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffwrhdu(
    infptr: *mut fitsfile,
    outstream: *mut FILE,
    status: *mut c_int,
) -> c_int {
    // FFI WRAPPER
    unsafe {
        let infptr = infptr.as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);

        ffwrhdu_safe(infptr, outstream, status)
    }
}

/// Write the current HDU to the output stream.
///
/// # Parameters
///
/// * `infptr`    — (I) FITS file pointer to input file
/// * `outstream` — (I) stream to write HDU to
/// * `status`    — (IO) error status
pub fn ffwrhdu_safe(infptr: &mut fitsfile, outstream: *mut FILE, status: &mut c_int) -> c_int {
    let mut hdustart: LONGLONG = 0;
    let mut hduend: LONGLONG = 0;
    let mut buffer: [c_char; IOBUFLEN as usize] = [0; IOBUFLEN as usize];

    if *status > 0 {
        return *status;
    }

    ffghadll_safe(infptr, Some(&mut hdustart), None, Some(&mut hduend), status);

    let nb = (hduend - hdustart) / BL!() as LONGLONG; /* number of blocks to copy */

    if nb > 0 {
        /* move to the start of the HDU */
        ffmbyt_safe(infptr, hdustart, REPORT_EOF, status);

        let mut outstream_cfile = CFile::from(outstream);

        for _ii in 0..(nb as usize) {
            ffgbyt(infptr, BL!(), cast_slice_mut(&mut buffer), status); /* read input block */
            let _ = outstream_cfile.write(cast_slice(&buffer[..BL!()])); /* write to output stream */
        }
    }
    *status
}

/// Insert an IMAGE extension following the current HDU.
///
/// # Parameters
///
/// * `fptr`   — (I) FITS file pointer
/// * `bitpix` — (I) bits per pixel
/// * `naxis`  — (I) number of axes in the array
/// * `naxes`  — (I) size of each axis
/// * `status` — (IO) error status
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffiimg(
    fptr: *mut fitsfile,
    bitpix: c_int,
    naxis: c_int,
    naxes: *const c_long,
    status: *mut c_int,
) -> c_int {
    // FFI WRAPPER
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        let naxes = slice::from_raw_parts(naxes, naxis as usize);

        ffiimg_safe(fptr, bitpix, naxis, naxes, status)
    }
}

/// Insert an IMAGE extension following the current HDU.
///
/// # Parameters
///
/// * `fptr`   — (I) FITS file pointer
/// * `bitpix` — (I) bits per pixel
/// * `naxis`  — (I) number of axes in the array
/// * `naxes`  — (I) size of each axis
/// * `status` — (IO) error status
pub fn ffiimg_safe(
    fptr: &mut fitsfile,
    bitpix: c_int,
    naxis: c_int,
    naxes: &[c_long],
    status: &mut c_int,
) -> c_int {
    let mut tnaxes: [LONGLONG; 99] = [0; 99];

    if *status > 0 {
        return *status;
    }

    if naxis > 99 {
        ffpmsg_str("NAXIS value is too large (>99)  (ffiimg)");
        *status = BAD_NAXIS;
        return *status;
    }

    let mut ii = 0;
    while (ii as c_int) < naxis {
        tnaxes[ii] = naxes[ii] as LONGLONG;
        ii += 1;
    }

    ffiimgll_safe(fptr, bitpix, naxis, &tnaxes, status);

    *status
}

/// Insert an IMAGE extension following the current HDU.
///
/// # Parameters
///
/// * `fptr`   — (I) FITS file pointer
/// * `bitpix` — (I) bits per pixel
/// * `naxis`  — (I) number of axes in the array
/// * `naxes`  — (I) size of each axis
/// * `status` — (IO) error status
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffiimgll(
    fptr: *mut fitsfile,
    bitpix: c_int,
    naxis: c_int,
    naxes: *const LONGLONG,
    status: *mut c_int,
) -> c_int {
    // FFI WRAPPER
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        let naxes = slice::from_raw_parts(naxes, naxis as usize);

        ffiimgll_safe(fptr, bitpix, naxis, naxes, status)
    }
}

/// Insert an IMAGE extension following the current HDU.
///
/// # Parameters
///
/// * `fptr`   — (I) FITS file pointer
/// * `bitpix` — (I) bits per pixel
/// * `naxis`  — (I) number of axes in the array
/// * `naxes`  — (I) size of each axis
/// * `status` — (IO) error status
pub fn ffiimgll_safe(
    fptr: &mut fitsfile,
    bitpix: c_int,
    naxis: c_int,
    naxes: &[LONGLONG],
    status: &mut c_int,
) -> c_int {
    let mut bytlen: c_int = 0;
    let mut nexthdu: c_int = 0;
    let mut maxhdu: c_int = 0;
    let mut onaxis: c_int = 0;

    let mut npixels: LONGLONG;
    let newstart: LONGLONG;

    let mut errmsg: [c_char; FLEN_ERRMSG] = [0; FLEN_ERRMSG];
    let mut card: [c_char; FLEN_CARD] = [0; FLEN_CARD];
    let mut naxiskey: [c_char; FLEN_KEYWORD] = [0; FLEN_KEYWORD];

    if *status > 0 {
        return *status;
    }

    if fptr.HDUposition != fptr.Fptr.curhdu {
        ffmahd_safe(fptr, (fptr.HDUposition) + 1, None, status);
    }

    maxhdu = fptr.Fptr.maxhdu;

    if *status != PREPEND_PRIMARY {
        /* if the current header is completely empty ...  */
        /* or, if we are at the end of the file, ... */
        let headstart = fptr.Fptr.get_headstart_as_slice();
        if (fptr.Fptr.headend == headstart[fptr.Fptr.curhdu as usize])
            || (((fptr.Fptr.curhdu) == maxhdu)
                && (headstart[(maxhdu + 1) as usize] >= fptr.Fptr.logfilesize))
        {
            /* then simply append new image extension */
            ffcrimll_safe(fptr, bitpix, naxis, naxes, status);
            return *status;
        }
    }

    if bitpix == 8 {
        bytlen = 1;
    } else if bitpix == 16 {
        bytlen = 2;
    } else if bitpix == 32 || bitpix == -32 {
        bytlen = 4;
    } else if bitpix == 64 || bitpix == -64 {
        bytlen = 8;
    } else {
        int_snprintf!(
            &mut errmsg,
            FLEN_ERRMSG,
            "Illegal value for BITPIX keyword: {}",
            bitpix,
        );
        ffpmsg_slice(&errmsg);
        *status = BAD_BITPIX; /* illegal bitpix value */
        return *status;
    }

    if naxis < 0 || naxis > 999 {
        int_snprintf!(
            &mut errmsg,
            FLEN_ERRMSG,
            "Illegal value for NAXIS keyword: {}",
            naxis,
        );
        ffpmsg_slice(&errmsg);
        *status = BAD_NAXIS;
        return *status;
    }

    for ii in 0..(naxis as usize) {
        if naxes[ii] < 0 {
            int_snprintf!(
                &mut errmsg,
                FLEN_ERRMSG,
                "Illegal value for NAXIS{} keyword: {}",
                ii + 1,
                naxes[ii] as c_long,
            );

            ffpmsg_slice(&errmsg);
            *status = BAD_NAXES;
            return *status;
        }
    }

    /* calculate number of pixels in the image */
    if naxis == 0 {
        npixels = 0;
    } else {
        npixels = naxes[0];
    }

    for ii in 1..(naxis as usize) {
        npixels *= naxes[ii];
    }

    let datasize: LONGLONG = npixels * LONGLONG::from(bytlen); /* size of image in bytes */
    let nblocks: c_long = (((datasize + (BL!() - 1)) / BL!()) + 1) as c_long; /* +1 for the header */

    if fptr.Fptr.writemode == READWRITE
    /* must have write access */
    {
        /* close the CHDU */
        ffrdef_safe(fptr, status); /* scan header to redefine structure */
        ffpdfl(fptr, status); /* insure correct data file values */
    } else {
        *status = READONLY_FILE;
        return *status;
    }

    if *status == PREPEND_PRIMARY {
        /* inserting a new primary array; the current primary */
        /* array must be transformed into an image extension. */

        *status = 0;
        ffmahd_safe(fptr, 1, None, status); /* move to the primary array */

        ffgidm_safe(fptr, &mut onaxis, status);
        if onaxis > 0 {
            ffkeyn_safe(cs!(c"NAXIS"), onaxis, &mut naxiskey, status);
        } else {
            strcpy_safe(&mut naxiskey, cs!(c"NAXIS"));
        }

        ffgcrd_safe(fptr, &naxiskey, &mut card, status); /* read last NAXIS keyword */

        ffikyj_safe(
            fptr,
            cs!(c"PCOUNT"),
            0,
            Some(cs!(c"required keyword")),
            status,
        ); /* add PCOUNT and */
        ffikyj_safe(
            fptr,
            cs!(c"GCOUNT"),
            1,
            Some(cs!(c"required keyword")),
            status,
        ); /* GCOUNT keywords */

        if *status > 0 {
            return *status;
        }

        if ffdkey_safe(fptr, cs!(c"EXTEND"), status) != 0 {
            /* delete the EXTEND keyword */
            *status = 0;
        }

        /* redefine internal structure for this HDU */
        ffrdef_safe(fptr, status);

        /* insert space for the primary array */
        if ffiblk(fptr, nblocks, -1, status) > 0 {
            /* insert the blocks */
            return *status;
        }

        nexthdu = 0; /* number of the new hdu */
        newstart = 0; /* starting addr of HDU */
    } else {
        let headstart = fptr.Fptr.get_headstart_as_slice();
        nexthdu = (fptr.Fptr.curhdu) + 1; /* number of the next (new) hdu */
        newstart = headstart[nexthdu as usize]; /* save starting addr of HDU */

        fptr.Fptr.hdutype = IMAGE_HDU; /* so that correct fill value is used */
        /* ffiblk also increments headstart for all following HDUs */
        if ffiblk(fptr, nblocks, 1, status) > 0 {
            /* insert the blocks */
            return *status;
        }
    }

    (fptr.Fptr.maxhdu) += 1; /* increment known number of HDUs in the file */

    let maxhdu = fptr.Fptr.maxhdu as usize;
    let curhdu = fptr.Fptr.curhdu as usize;
    let mut ii = maxhdu;
    let headstart = fptr.Fptr.get_headstart_as_mut_slice();
    while ii > curhdu {
        headstart[ii + 1] = headstart[ii]; /* incre start addr */
        ii -= 1;
    }

    if nexthdu == 0 {
        headstart[1] = (nblocks * BL!()) as LONGLONG; /* start of the old Primary array */
    }

    headstart[nexthdu as usize] = newstart; /* set starting addr of HDU */

    /* set default parameters for this new empty HDU */
    let headstart = fptr.Fptr.get_headstart_as_slice();
    let hs_item = headstart[nexthdu as usize];

    fptr.Fptr.curhdu = nexthdu; /* we are now located at the next HDU */
    fptr.HDUposition = nexthdu; /* we are now located at the next HDU */
    fptr.Fptr.nextkey = hs_item;
    fptr.Fptr.headend = hs_item;
    fptr.Fptr.datastart = (hs_item) + BL!();
    fptr.Fptr.hdutype = IMAGE_HDU; /* might need to be reset... */

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

    /* redefine internal structure for this HDU */
    ffrdef_safe(fptr, status);
    *status
}

/// Insert an ASCII table extension following the current HDU.
///
/// # Parameters
///
/// * `fptr`    — (I) FITS file pointer
/// * `naxis1`  — (I) width of row in the table
/// * `naxis2`  — (I) number of rows in the table
/// * `tfields` — (I) number of columns in the table
/// * `ttype`   — (I) name of each column
/// * `tbcol`   — (I) byte offset in row to each column
/// * `tform`   — (I) value of TFORMn keyword for each column
/// * `tunit`   — (I) value of TUNITn keyword for each column
/// * `extnmx`  — (I) value of EXTNAME keyword, if any
/// * `status`  — (IO) error status
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffitab(
    fptr: *mut fitsfile,
    naxis1: LONGLONG,
    naxis2: LONGLONG,
    tfields: c_int,
    ttype: *const *const c_char,
    tbcol: *const c_long,
    tform: *const *const c_char,
    tunit: *const *const c_char,
    extnmx: *const c_char,
    status: *mut c_int,
) -> c_int {
    // FFI WRAPPER
    unsafe {
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);

        nullable_slice_cstr!(extnmx);

        let tkeywords = TKeywords::new(tfields, ttype, tform, tunit);
        let (v_ttype, v_tform, v_tunit) = tkeywords.tkeywords_to_vecs();

        let tbcol_slice = if tbcol.is_null() {
            None
        } else {
            Some(slice::from_raw_parts(tbcol, cmp::max(5, tfields as usize)))
        };

        ffitab_safe(
            fptr,
            naxis1,
            naxis2,
            tfields,
            &v_ttype,
            tbcol_slice,
            &v_tform,
            v_tunit.as_deref(),
            extnmx,
            status,
        )
    }
}

/// Insert an ASCII table extension following the current HDU.
///
/// # Parameters
///
/// * `fptr`    — (I) FITS file pointer
/// * `naxis1`  — (I) width of row in the table
/// * `naxis2`  — (I) number of rows in the table
/// * `tfields` — (I) number of columns in the table
/// * `v_ttype` — (I) name of each column
/// * `tbcol`   — (I) byte offset in row to each column
/// * `v_tform` — (I) value of TFORMn keyword for each column
/// * `v_tunit` — (I) value of TUNITn keyword for each column
/// * `extnmx`  — (I) value of EXTNAME keyword, if any
/// * `status`  — (IO) error status
pub fn ffitab_safe(
    fptr: &mut fitsfile,
    naxis1: LONGLONG,
    naxis2: LONGLONG,
    tfields: c_int,
    v_ttype: &[Option<&[c_char]>],
    tbcol: Option<&[c_long]>,
    v_tform: &[&[c_char]],
    v_tunit: Option<&[Option<&[c_char]>]>,
    extnmx: Option<&[c_char]>,
    status: &mut c_int,
) -> c_int {
    let mut nunit: c_int = 0;

    let mut rowlen: c_long;

    let mut errmsg: [c_char; FLEN_ERRMSG] = [0; FLEN_ERRMSG];
    let mut extnm: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];

    if *status > 0 {
        return *status;
    }

    extnm[0] = 0;
    if let Some(extnmx) = extnmx {
        strncat_safe(&mut extnm, extnmx, FLEN_VALUE - 1);
    }

    if fptr.HDUposition != fptr.Fptr.curhdu {
        ffmahd_safe(fptr, (fptr.HDUposition) + 1, None, status);
    }

    let maxhdu: c_int = fptr.Fptr.maxhdu;

    let headstart = fptr.Fptr.get_headstart_as_slice();
    /* if the current header is completely empty or, if we are at the end of the file, ...  */
    if (fptr.Fptr.headend == headstart[fptr.Fptr.curhdu as usize])
        || (((fptr.Fptr.curhdu) == maxhdu)
            && (headstart[(maxhdu + 1) as usize] >= fptr.Fptr.logfilesize))
    {
        /* then simply append new image extension */
        ffcrtb_safe(
            fptr,
            ASCII_TBL,
            naxis2,
            tfields,
            v_ttype,
            v_tform,
            v_tunit,
            Some(&extnm),
            status,
        );
        return *status;
    }

    if naxis1 < 0 {
        *status = NEG_WIDTH;
        return *status;
    } else if naxis2 < 0 {
        *status = NEG_ROWS;
        return *status;
    } else if tfields < 0 || tfields > 999 {
        int_snprintf!(
            &mut errmsg,
            FLEN_ERRMSG,
            "Illegal value for TFIELDS keyword: {}",
            tfields,
        );
        ffpmsg_slice(&errmsg);
        *status = BAD_TFIELDS;
        return *status;
    }

    /* count number of optional TUNIT keywords to be written */
    for ii in 0..(tfields as usize) {
        if let Some(v_tunit) = v_tunit
            && v_tunit[ii].is_some()
        {
            nunit += 1;
        }
    }

    if extnm[0] != 0 {
        nunit += 1; /* add one for the EXTNAME keyword */
    }

    rowlen = naxis1 as c_long;

    /* spacing not defined? */
    let mut gotmem = false;
    let ncols = cmp::max(5, tfields as usize); /* arrays, so allocate at least 20 bytes */

    let mut tbcol = match tbcol.is_none() || (naxis1 == 0 && tfields != 0) {
        true => {
            gotmem = true;
            vec![0 as c_long; ncols]
        }
        false => {
            let x = tbcol.unwrap().to_vec();
            if x[0] == 0 {
                gotmem = true;
                vec![0 as c_long; ncols]
            } else {
                x
            }
        }
    };

    if gotmem {
        /* calculate width of a row and starting position of each column. */
        /* Each column will be separated by 1 blank space */
        ffgabc_safe(tfields, v_tform, 1, &mut rowlen, &mut tbcol, status);
    }

    let nhead: c_int = (9 + (3 * tfields) + nunit + 35) / 36; /* no. of header blocks */
    let datasize: LONGLONG = (rowlen as LONGLONG) * naxis2; /* size of table in bytes */
    let nblocks: c_long = (((datasize + (BL!() - 1)) / BL!()) + LONGLONG::from(nhead)) as c_long; /* size of HDU */

    if fptr.Fptr.writemode == READWRITE {
        /* must have write access */

        /* close the CHDU */
        ffrdef_safe(fptr, status); /* scan header to redefine structure */
        ffpdfl(fptr, status); /* insure correct data file values */
    } else {
        *status = READONLY_FILE;
        return *status;
    }
    let headstart = fptr.Fptr.get_headstart_as_slice();

    let nexthdu: c_int = (fptr.Fptr.curhdu) + 1; /* number of the next (new) hdu */
    let newstart: LONGLONG = headstart[nexthdu as usize]; /* save starting addr of HDU */

    fptr.Fptr.hdutype = ASCII_TBL; /* so that correct fill value is used */

    /* ffiblk also increments headstart for all following HDUs */
    if ffiblk(fptr, nblocks, 1, status) > 0 {
        /* insert the blocks */
        return *status;
    }

    (fptr.Fptr.maxhdu) += 1; /* increment known number of HDUs in the file */

    let maxhdu = fptr.Fptr.maxhdu as usize;
    let curhdu = fptr.Fptr.curhdu as usize;
    let headstart = fptr.Fptr.get_headstart_as_mut_slice();
    let mut ii = maxhdu;
    while ii > curhdu {
        headstart[ii + 1] = headstart[ii]; /* incre start addr */
        ii -= 1;
    }

    headstart[nexthdu as usize] = newstart; /* set starting addr of HDU */

    /* set default parameters for this new empty HDU */
    let headstart = fptr.Fptr.get_headstart_as_slice();
    let hs_item = headstart[nexthdu as usize];

    fptr.Fptr.curhdu = nexthdu; /* we are now located at the next HDU */
    fptr.HDUposition = nexthdu; /* we are now located at the next HDU */
    fptr.Fptr.nextkey = hs_item;
    fptr.Fptr.headend = hs_item;
    fptr.Fptr.datastart = (hs_item) + (LONGLONG::from(nhead) * BL!());
    fptr.Fptr.hdutype = ASCII_TBL; /* might need to be reset... */

    /* write the required header keywords */

    ffphtb_safe(
        fptr,
        rowlen as LONGLONG,
        naxis2,
        tfields,
        v_ttype,
        Some(&tbcol),
        v_tform,
        v_tunit,
        Some(&extnm),
        status,
    );

    /* redefine internal structure for this HDU */
    ffrdef_safe(fptr, status);
    *status
}

/// Insert a binary table extension following the current HDU.
///
/// # Parameters
///
/// * `fptr`    — (I) FITS file pointer
/// * `naxis2`  — (I) number of rows in the table
/// * `tfields` — (I) number of columns in the table
/// * `ttype`   — (I) name of each column
/// * `tform`   — (I) value of TFORMn keyword for each column
/// * `tunit`   — (I) value of TUNITn keyword for each column
/// * `extnmx`  — (I) value of EXTNAME keyword, if any
/// * `pcount`  — (I) size of special data area (heap)
/// * `status`  — (IO) error status
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffibin(
    fptr: *mut fitsfile,
    naxis2: LONGLONG,
    tfields: c_int,
    ttype: *const *const c_char,
    tform: *const *const c_char,
    tunit: *const *const c_char,
    extnmx: *const c_char,
    pcount: LONGLONG,
    status: *mut c_int,
) -> c_int {
    // FFI WRAPPER
    unsafe {
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);

        nullable_slice_cstr!(extnmx);

        let tkeywords = TKeywords::new(tfields, ttype, tform, tunit);
        let (v_ttype, v_tform, v_tunit) = tkeywords.tkeywords_to_vecs();

        ffibin_safe(
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

/// Insert a Binary table extension following the current HDU.
///
/// # Parameters
///
/// * `fptr`    — (I) FITS file pointer
/// * `naxis2`  — (I) number of rows in the table
/// * `tfields` — (I) number of columns in the table
/// * `ttype`   — (I) name of each column
/// * `tform`   — (I) value of TFORMn keyword for each column
/// * `tunit`   — (I) value of TUNITn keyword for each column
/// * `extnmx`  — (I) value of EXTNAME keyword, if any
/// * `pcount`  — (I) size of special data area (heap)
/// * `status`  — (IO) error status
pub fn ffibin_safe(
    fptr: &mut fitsfile,
    naxis2: LONGLONG,
    tfields: c_int,
    ttype: &[Option<&[c_char]>],
    tform: &[&[c_char]],
    tunit: Option<&[Option<&[c_char]>]>,
    extnmx: Option<&[c_char]>,
    pcount: LONGLONG,
    status: &mut c_int,
) -> c_int {
    let mut nunit: c_int = 0;

    let mut datacode: c_int = 0;
    let mut naxis1: LONGLONG = 0;

    let mut repeat: c_long = 0;
    let mut width: c_long = 0;

    let mut errmsg: [c_char; FLEN_ERRMSG] = [0; FLEN_ERRMSG];
    let mut extnm: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];

    if *status > 0 {
        return *status;
    }

    extnm[0] = 0;
    if let Some(extnmx) = extnmx {
        strncat_safe(&mut extnm, extnmx, FLEN_VALUE - 1);
    }

    if fptr.HDUposition != fptr.Fptr.curhdu {
        ffmahd_safe(fptr, (fptr.HDUposition) + 1, None, status);
    }

    let maxhdu: c_int = fptr.Fptr.maxhdu;

    let headstart = fptr.Fptr.get_headstart_as_slice();
    /* if the current header is completely empty ...  */
    /* or, if we are at the end of the file, ... */
    if (fptr.Fptr.headend == headstart[fptr.Fptr.curhdu as usize])
        || (((fptr.Fptr.curhdu) == maxhdu)
            && (headstart[(maxhdu + 1) as usize] >= fptr.Fptr.logfilesize))
    {
        /* then simply append new image extension */
        ffcrtb_safe(
            fptr,
            BINARY_TBL,
            naxis2,
            tfields,
            ttype,
            tform,
            tunit,
            Some(&extnm),
            status,
        );
        return *status;
    }

    if naxis2 < 0 {
        *status = NEG_ROWS;
        return *status;
    } else if tfields < 0 || tfields > 999 {
        int_snprintf!(
            &mut errmsg,
            FLEN_ERRMSG,
            "Illegal value for TFIELDS keyword: {}",
            tfields,
        );
        ffpmsg_slice(&errmsg);
        *status = BAD_TFIELDS;
        return *status;
    }

    /* count number of optional TUNIT keywords to be written */
    for ii in 0..(tfields as usize) {
        if let Some(v_tunit) = tunit
            && v_tunit[ii].is_some()
        {
            nunit += 1;
        }
    }

    if extnm[0] != 0 {
        nunit += 1; /* add one for the EXTNAME keyword */
    }

    let nhead: c_int = (9 + (2 * tfields) + nunit + 35) / 36; /* no. of header blocks */

    /* calculate total width of the table */
    for ii in 0..(tfields as usize) {
        let tform_item = tform[ii];
        ffbnfm_safe(
            tform_item,
            Some(&mut datacode),
            Some(&mut repeat),
            Some(&mut width),
            status,
        );

        if datacode == TBIT {
            naxis1 += (repeat as LONGLONG + 7) / 8;
        } else if datacode == TSTRING {
            naxis1 += repeat as LONGLONG;
        } else {
            naxis1 += (repeat * width) as LONGLONG;
        }
    }

    let datasize: LONGLONG = ((naxis1 as LONGLONG) * naxis2) + pcount; /* size of table in bytes */
    let nblocks: c_long = (((datasize + (BL!() - 1)) / BL!()) + LONGLONG::from(nhead)) as c_long; /* size of HDU */

    if fptr.Fptr.writemode == READWRITE {
        /* must have write access */

        /* close the CHDU */
        ffrdef_safe(fptr, status); /* scan header to redefine structure */
        ffpdfl(fptr, status); /* insure correct data file values */
    } else {
        *status = READONLY_FILE;
        return *status;
    }

    let headstart = fptr.Fptr.get_headstart_as_slice();

    let nexthdu: c_int = (fptr.Fptr.curhdu) + 1; /* number of the next (new) hdu */
    let newstart: LONGLONG = headstart[nexthdu as usize]; /* save starting addr of HDU */

    fptr.Fptr.hdutype = BINARY_TBL; /* so that correct fill value is used */

    /* ffiblk also increments headstart for all following HDUs */
    if ffiblk(fptr, nblocks, 1, status) > 0 {
        /* insert the blocks */
        return *status;
    }

    (fptr.Fptr.maxhdu) += 1; /* increment known number of HDUs in the file */

    let maxhdu = fptr.Fptr.maxhdu as usize;
    let curhdu = fptr.Fptr.curhdu as usize;
    let headstart = fptr.Fptr.get_headstart_as_mut_slice();
    let mut ii = maxhdu;
    while ii > curhdu {
        headstart[ii + 1] = headstart[ii]; /* incre start addr */
        ii -= 1;
    }

    headstart[nexthdu as usize] = newstart; /* set starting addr of HDU */

    /* set default parameters for this new empty HDU */
    let headstart = fptr.Fptr.get_headstart_as_slice();
    let hs_item = headstart[nexthdu as usize];

    fptr.Fptr.curhdu = nexthdu; /* we are now located at the next HDU */
    fptr.HDUposition = nexthdu; /* we are now located at the next HDU */
    fptr.Fptr.nextkey = hs_item;
    fptr.Fptr.headend = hs_item;
    fptr.Fptr.datastart = (hs_item) + (LONGLONG::from(nhead) * BL!());
    fptr.Fptr.hdutype = BINARY_TBL; /* might need to be reset... */

    /* write the required header keywords. This will write PCOUNT = 0 */
    /* so that the variable length data will be written at the right place */
    ffphbn_safe(
        fptr,
        naxis2,
        tfields,
        ttype,
        tform,
        tunit,
        Some(&extnm),
        pcount,
        status,
    );

    /* redefine internal structure for this HDU (with PCOUNT = 0) */
    ffrdef_safe(fptr, status);

    *status
}

/// Delete the CHDU.  If the CHDU is the primary array, then replace the HDU
/// with an empty primary array with no data.   Return the
/// type of the new CHDU after the old CHDU is deleted.
///
/// # Parameters
///
/// * `fptr`    — (I) FITS file pointer
/// * `hdutype` — (O) type of the new CHDU after deletion
/// * `status`  — (IO) error status
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffdhdu(
    fptr: *mut fitsfile,
    hdutype: *mut c_int,
    status: *mut c_int,
) -> c_int {
    // FFI WRAPPER
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let hdutype = hdutype.as_mut();

        ffdhdu_safe(fptr, hdutype, status)
    }
}

/// Delete the CHDU.  If the CHDU is the primary array, then replace the HDU
/// with an empty primary array with no data.   Return the
/// type of the new CHDU after the old CHDU is deleted.
///
/// # Parameters
///
/// * `fptr`    — (I) FITS file pointer
/// * `hdutype` — (O) type of the new CHDU after deletion
/// * `status`  — (IO) error status
pub fn ffdhdu_safe(fptr: &mut fitsfile, hdutype: Option<&mut c_int>, status: &mut c_int) -> c_int {
    let mut tmptype = 0;
    let mut nblocks: c_long = 0;
    let naxes: [c_long; 1] = [0; 1];

    if *status > 0 {
        return *status;
    }

    if fptr.HDUposition != fptr.Fptr.curhdu {
        ffmahd_safe(fptr, (fptr.HDUposition) + 1, None, status);
    }

    if fptr.Fptr.curhdu == 0 {
        /* replace primary array with null image */

        /* ignore any existing keywords */
        fptr.Fptr.headend = 0;
        fptr.Fptr.nextkey = 0;

        /* write default primary array header */
        ffphpr_safe(fptr, 1, 8, 0, &naxes, 0, 1, 1, status);

        /* calc number of blocks to delete (leave just 1 block) */
        let headstart = fptr.Fptr.get_headstart_as_slice();
        nblocks = ((headstart[(fptr.Fptr.curhdu + 1) as usize] - BL!()) / BL!()) as c_long;

        /* ffdblk also updates the starting address of all following HDUs */
        if nblocks > 0 && ffdblk(fptr, nblocks, status) > 0 {
            /* delete the HDU */
            return *status;
        }

        /* this might not be necessary, but is doesn't hurt */
        fptr.Fptr.datastart = DATA_UNDEFINED as LONGLONG;

        ffrdef_safe(fptr, status); /* reinitialize the primary array */
    } else {
        let headstart = fptr.Fptr.get_headstart_as_slice();

        /* calc number of blocks to delete */
        nblocks = ((headstart[(fptr.Fptr.curhdu + 1) as usize]
            - headstart[fptr.Fptr.curhdu as usize])
            / BL!()) as c_long;

        /* ffdblk also updates the starting address of all following HDUs */
        if ffdblk(fptr, nblocks, status) > 0 {
            /* delete the HDU */
            return *status;
        }

        /* delete the CHDU from the list of HDUs */
        let curhdu = 1 + fptr.Fptr.curhdu as usize;
        let maxhdu = fptr.Fptr.maxhdu as usize;
        let headstart = fptr.Fptr.get_headstart_as_mut_slice();
        for ii in curhdu..=maxhdu {
            headstart[ii] = headstart[ii + 1];
        }

        headstart[maxhdu + 1] = 0;
        (fptr.Fptr.maxhdu) -= 1; /* decrement the known number of HDUs */

        if ffrhdu_safe(fptr, Some(&mut tmptype), status) > 0 {
            /* initialize next HDU */

            /* failed (end of file?), so move back one HDU */
            *status = 0;
            ffcmsg_safe(); /* clear extraneous error messages */
            ffgext(fptr, (fptr.Fptr.curhdu) - 1, Some(&mut tmptype), status);
        }
    }

    if let Some(hdutype) = hdutype {
        *hdutype = tmptype;
    }

    *status
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::KeywordDatatypeMut;
    use crate::aliases::rust_api::*;
    use crate::fitsio::{
        BAD_BITPIX, BAD_NAXES, BAD_NAXIS, BAD_TFIELDS, BINARY_TBL, BYTE_IMG, DOUBLE_IMG,
        FLEN_VALUE, LONGLONG, LONGLONG_IMG, NEG_ROWS, NEG_WIDTH, READONLY, READONLY_FILE,
        READWRITE, SAME_FILE, SHORT_IMG, fitsfile,
    };
    use crate::helpers::testhelpers::{from_buf, to_buf, with_temp_file};
    use libc::{c_char, c_int, c_long};

    /// Make a NUL-terminated `Vec<c_char>` from a `&str`.
    fn cc(s: &str) -> Vec<c_char> {
        let mut v: Vec<c_char> = s.bytes().map(|b| b as c_char).collect();
        v.push(0);
        v
    }

    /// Derive a sibling path in the same temp dir as `base`.
    fn sibling(base: &str, name: &str) -> String {
        let p = std::path::Path::new(base);
        p.with_file_name(name).to_str().unwrap().to_string()
    }

    #[test]
    fn test_cpfl_copy_all() {
        // Test ffcpfl - copy all HDUs from input to output.
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let path = to_buf(filename);
            let path2_str = sibling(filename, "test_edithdu2.fits");
            let path2 = to_buf(&path2_str);
            let naxes: [c_long; 1] = [10];
            let ttype = [Some(cc("COL1"))];
            let ttype_ref: Vec<Option<&[c_char]>> = ttype.iter().map(|o| o.as_deref()).collect();
            let tform_v = [cc("1J")];
            let tform_ref: Vec<&[c_char]> = tform_v.iter().map(|v| v.as_slice()).collect();
            let data: [c_long; 3] = [100, 200, 300];

            // Create input file with primary + 2 extensions.
            let mut inf: Option<Box<fitsfile>> = None;
            fits_create_file(&mut inf, &path, &mut status);
            fits_write_imghdr(
                inf.as_deref_mut().unwrap(),
                BYTE_IMG,
                1,
                &naxes,
                &mut status,
            );
            fits_create_tbl(
                inf.as_deref_mut().unwrap(),
                BINARY_TBL,
                3,
                1,
                &ttype_ref,
                &tform_ref,
                None,
                Some(&cc("EXT1")),
                &mut status,
            );
            fits_write_col_lng(inf.as_deref_mut().unwrap(), 1, 1, 1, 3, &data, &mut status);
            fits_create_tbl(
                inf.as_deref_mut().unwrap(),
                BINARY_TBL,
                2,
                1,
                &ttype_ref,
                &tform_ref,
                None,
                Some(&cc("EXT2")),
                &mut status,
            );
            fits_close_file(inf.take().unwrap(), &mut status);
            assert_eq!(status, 0, "input setup failed");

            // Create empty output file.
            let mut out: Option<Box<fitsfile>> = None;
            fits_create_file(&mut out, &path2, &mut status);
            fits_write_imghdr(out.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            fits_close_file(out.take().unwrap(), &mut status);
            assert_eq!(status, 0, "output setup failed");

            // Open both and copy all HDUs using ffcpfl.
            fits_open_file(&mut inf, &path, READONLY, &mut status);
            fits_open_file(&mut out, &path2, READWRITE, &mut status);
            fits_movabs_hdu(inf.as_deref_mut().unwrap(), 2, None, &mut status);
            fits_copy_file(
                inf.as_deref_mut().unwrap(),
                out.as_deref_mut().unwrap(),
                1,
                1,
                1,
                &mut status,
            );
            assert_eq!(status, 0, "ffcpfl failed");
            fits_close_file(inf.take().unwrap(), &mut status);
            fits_close_file(out.take().unwrap(), &mut status);

            // Verify output has 4 HDUs (original empty + 3 copied).
            fits_open_file(&mut out, &path2, READONLY, &mut status);
            let mut numhdus: c_int = 0;
            fits_get_num_hdus(out.as_deref_mut().unwrap(), &mut numhdus, &mut status);
            assert_eq!(status, 0);
            assert_eq!(numhdus, 4);
            fits_close_file(out.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_cpfl_copy_current_only() {
        // Test ffcpfl - copy only current HDU.
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let path = to_buf(filename);
            let path2_str = sibling(filename, "test_edithdu2.fits");
            let path2 = to_buf(&path2_str);
            let naxes: [c_long; 1] = [10];
            let ttype = [Some(cc("COL1"))];
            let ttype_ref: Vec<Option<&[c_char]>> = ttype.iter().map(|o| o.as_deref()).collect();
            let tform_v = [cc("1J")];
            let tform_ref: Vec<&[c_char]> = tform_v.iter().map(|v| v.as_slice()).collect();
            let data: [c_long; 3] = [100, 200, 300];

            let mut inf: Option<Box<fitsfile>> = None;
            fits_create_file(&mut inf, &path, &mut status);
            fits_write_imghdr(
                inf.as_deref_mut().unwrap(),
                BYTE_IMG,
                1,
                &naxes,
                &mut status,
            );
            fits_create_tbl(
                inf.as_deref_mut().unwrap(),
                BINARY_TBL,
                3,
                1,
                &ttype_ref,
                &tform_ref,
                None,
                Some(&cc("FIRST")),
                &mut status,
            );
            fits_write_col_lng(inf.as_deref_mut().unwrap(), 1, 1, 1, 3, &data, &mut status);
            fits_create_tbl(
                inf.as_deref_mut().unwrap(),
                BINARY_TBL,
                2,
                1,
                &ttype_ref,
                &tform_ref,
                None,
                Some(&cc("SECOND")),
                &mut status,
            );
            fits_close_file(inf.take().unwrap(), &mut status);
            assert_eq!(status, 0);

            let mut out: Option<Box<fitsfile>> = None;
            fits_create_file(&mut out, &path2, &mut status);
            fits_write_imghdr(out.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            fits_close_file(out.take().unwrap(), &mut status);
            assert_eq!(status, 0);

            fits_open_file(&mut inf, &path, READONLY, &mut status);
            fits_open_file(&mut out, &path2, READWRITE, &mut status);
            fits_movabs_hdu(inf.as_deref_mut().unwrap(), 2, None, &mut status);
            fits_copy_file(
                inf.as_deref_mut().unwrap(),
                out.as_deref_mut().unwrap(),
                0,
                1,
                0,
                &mut status,
            );
            assert_eq!(status, 0);
            fits_close_file(inf.take().unwrap(), &mut status);
            fits_close_file(out.take().unwrap(), &mut status);

            fits_open_file(&mut out, &path2, READONLY, &mut status);
            let mut numhdus: c_int = 0;
            fits_get_num_hdus(out.as_deref_mut().unwrap(), &mut numhdus, &mut status);
            assert_eq!(numhdus, 2);
            fits_movabs_hdu(out.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut extname = [0 as c_char; FLEN_VALUE];
            fits_read_key(
                out.as_deref_mut().unwrap(),
                KeywordDatatypeMut::TSTRING(&mut extname),
                &cc("EXTNAME"),
                None,
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(from_buf(&extname), "FIRST");
            fits_close_file(out.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_cpfl_copy_following() {
        // Test ffcpfl - copy following HDUs only.
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let path = to_buf(filename);
            let path2_str = sibling(filename, "test_edithdu2.fits");
            let path2 = to_buf(&path2_str);
            let naxes: [c_long; 1] = [10];
            let ttype = [Some(cc("COL1"))];
            let ttype_ref: Vec<Option<&[c_char]>> = ttype.iter().map(|o| o.as_deref()).collect();
            let tform_v = [cc("1J")];
            let tform_ref: Vec<&[c_char]> = tform_v.iter().map(|v| v.as_slice()).collect();

            let mut inf: Option<Box<fitsfile>> = None;
            fits_create_file(&mut inf, &path, &mut status);
            fits_write_imghdr(
                inf.as_deref_mut().unwrap(),
                BYTE_IMG,
                1,
                &naxes,
                &mut status,
            );
            for nm in ["A", "B", "C"] {
                fits_create_tbl(
                    inf.as_deref_mut().unwrap(),
                    BINARY_TBL,
                    2,
                    1,
                    &ttype_ref,
                    &tform_ref,
                    None,
                    Some(&cc(nm)),
                    &mut status,
                );
            }
            fits_close_file(inf.take().unwrap(), &mut status);
            assert_eq!(status, 0);

            let mut out: Option<Box<fitsfile>> = None;
            fits_create_file(&mut out, &path2, &mut status);
            fits_write_imghdr(out.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            fits_close_file(out.take().unwrap(), &mut status);
            assert_eq!(status, 0);

            fits_open_file(&mut inf, &path, READONLY, &mut status);
            fits_open_file(&mut out, &path2, READWRITE, &mut status);
            fits_movabs_hdu(inf.as_deref_mut().unwrap(), 2, None, &mut status); // Move to A
            fits_copy_file(
                inf.as_deref_mut().unwrap(),
                out.as_deref_mut().unwrap(),
                0,
                0,
                1,
                &mut status,
            );
            assert_eq!(status, 0);
            fits_close_file(inf.take().unwrap(), &mut status);
            fits_close_file(out.take().unwrap(), &mut status);

            fits_open_file(&mut out, &path2, READONLY, &mut status);
            let mut numhdus: c_int = 0;
            fits_get_num_hdus(out.as_deref_mut().unwrap(), &mut numhdus, &mut status);
            assert_eq!(numhdus, 3);
            let mut extname = [0 as c_char; FLEN_VALUE];
            fits_movabs_hdu(out.as_deref_mut().unwrap(), 2, None, &mut status);
            fits_read_key(
                out.as_deref_mut().unwrap(),
                KeywordDatatypeMut::TSTRING(&mut extname),
                &cc("EXTNAME"),
                None,
                &mut status,
            );
            assert_eq!(from_buf(&extname), "B");
            fits_movabs_hdu(out.as_deref_mut().unwrap(), 3, None, &mut status);
            fits_read_key(
                out.as_deref_mut().unwrap(),
                KeywordDatatypeMut::TSTRING(&mut extname),
                &cc("EXTNAME"),
                None,
                &mut status,
            );
            assert_eq!(from_buf(&extname), "C");
            fits_close_file(out.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_cpht_copy_table_subset() {
        // Test ffcpht - copy table structure with limited row range.
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let path = to_buf(filename);
            let path2_str = sibling(filename, "test_edithdu2.fits");
            let path2 = to_buf(&path2_str);
            let ttype = [Some(cc("COL1")), Some(cc("COL2"))];
            let ttype_ref: Vec<Option<&[c_char]>> = ttype.iter().map(|o| o.as_deref()).collect();
            let tform_v = [cc("1J"), cc("1E")];
            let tform_ref: Vec<&[c_char]> = tform_v.iter().map(|v| v.as_slice()).collect();
            let jdata: [c_long; 5] = [10, 20, 30, 40, 50];
            let edata: [f32; 5] = [1.1, 2.2, 3.3, 4.4, 5.5];

            let mut inf: Option<Box<fitsfile>> = None;
            fits_create_file(&mut inf, &path, &mut status);
            fits_write_imghdr(inf.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            fits_create_tbl(
                inf.as_deref_mut().unwrap(),
                BINARY_TBL,
                5,
                2,
                &ttype_ref,
                &tform_ref,
                None,
                Some(&cc("DATA")),
                &mut status,
            );
            fits_write_col_lng(inf.as_deref_mut().unwrap(), 1, 1, 1, 5, &jdata, &mut status);
            fits_write_col_flt(inf.as_deref_mut().unwrap(), 2, 1, 1, 5, &edata, &mut status);
            fits_close_file(inf.take().unwrap(), &mut status);
            assert_eq!(status, 0);

            let mut out: Option<Box<fitsfile>> = None;
            fits_create_file(&mut out, &path2, &mut status);
            fits_write_imghdr(out.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            fits_close_file(out.take().unwrap(), &mut status);
            assert_eq!(status, 0);

            fits_open_file(&mut inf, &path, READONLY, &mut status);
            fits_open_file(&mut out, &path2, READWRITE, &mut status);
            fits_movabs_hdu(inf.as_deref_mut().unwrap(), 2, None, &mut status);
            fits_copy_hdutab(
                inf.as_deref_mut().unwrap(),
                out.as_deref_mut().unwrap(),
                2,
                3,
                &mut status,
            );
            assert_eq!(status, 0, "ffcpht failed");
            fits_close_file(inf.take().unwrap(), &mut status);
            fits_close_file(out.take().unwrap(), &mut status);

            fits_open_file(&mut out, &path2, READONLY, &mut status);
            fits_movabs_hdu(out.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut nrows: LONGLONG = 0;
            fits_get_num_rowsll(out.as_deref_mut().unwrap(), &mut nrows, &mut status);
            assert_eq!(nrows, 3);
            let mut jresult = [0 as c_long; 3];
            let mut eresult = [0.0f32; 3];
            fits_read_col_lng(
                out.as_deref_mut().unwrap(),
                1,
                1,
                1,
                3,
                0,
                &mut jresult,
                None,
                &mut status,
            );
            fits_read_col_flt(
                out.as_deref_mut().unwrap(),
                2,
                1,
                1,
                3,
                0.0,
                &mut eresult,
                None,
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(jresult[0], 20);
            assert_eq!(jresult[1], 30);
            assert_eq!(jresult[2], 40);
            assert!(eresult[0] >= 2.1 && eresult[0] <= 2.3);
            assert!(eresult[1] >= 3.2 && eresult[1] <= 3.4);
            assert!(eresult[2] >= 4.3 && eresult[2] <= 4.5);
            fits_close_file(out.take().unwrap(), &mut status);
        });
    }

    #[test]
    // ffwrhdu takes a C `FILE *`, as in the C, so exercising it needs a real
    // stdio stream from libc's fopen.  miri cannot call foreign functions, so
    // this test can never run under it -- unlike the rest of the suite, which
    // it is worth keeping noise-free so that real findings stand out.
    #[cfg_attr(miri, ignore)]
    fn test_wrhdu_write_to_stream() {
        // Test ffwrhdu - write HDU to output stream.
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let path = to_buf(filename);
            let path2_str = sibling(filename, "test_edithdu2.fits");
            let path2 = to_buf(&path2_str);
            let naxes: [c_long; 2] = [10, 10];
            let mut data = [0i16; 100];
            for (i, d) in data.iter_mut().enumerate() {
                *d = (i * 2) as i16;
            }

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &path, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), SHORT_IMG, 2, &naxes, &mut status);
            fits_write_img_sht(f.as_deref_mut().unwrap(), 1, 1, 100, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0);

            // Open and write HDU to a stream (file).
            fits_open_file(&mut f, &path, READONLY, &mut status);
            let cpath = alloc::ffi::CString::new(path2_str.as_str()).unwrap();
            let mode = alloc::ffi::CString::new("wb").unwrap();
            let stream = unsafe { libc::fopen(cpath.as_ptr(), mode.as_ptr()) };
            assert!(!stream.is_null());
            ffwrhdu_safe(
                f.as_deref_mut().unwrap(),
                stream.cast::<crate::c_types::FILE>(),
                &mut status,
            );
            assert_eq!(status, 0, "ffwrhdu failed");
            unsafe {
                libc::fclose(stream);
            }
            fits_close_file(f.take().unwrap(), &mut status);

            // Verify the stream output is a valid FITS file.
            let mut f2: Option<Box<fitsfile>> = None;
            fits_open_file(&mut f2, &path2, READONLY, &mut status);
            assert_eq!(status, 0, "stream output not a valid FITS file");
            let mut result = [0i16; 100];
            fits_read_img_sht(
                f2.as_deref_mut().unwrap(),
                1,
                1,
                100,
                0,
                &mut result,
                None,
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[0], 0);
            assert_eq!(result[50], 100);
            assert_eq!(result[99], 198);
            fits_close_file(f2.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_error_status_returns() {
        // Test that functions return early when status > 0.
        // Same-pointer ffcopy/ffcphd/ffcpdt use the raw C-ABI wrappers so the
        // same handle can be passed twice.
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let path = to_buf(filename);
            let naxes: [c_long; 1] = [10];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &path, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 1, &naxes, &mut status);
            assert_eq!(status, 0);

            let fp: *mut fitsfile = f.as_deref_mut().unwrap();
            status = 1;
            unsafe {
                #[allow(deprecated)]
                ffcopy(fp, fp, 0, &raw mut status);
            }
            assert_eq!(status, 1);

            unsafe {
                #[allow(deprecated)]
                ffcphd(fp, fp, &raw mut status);
            }
            assert_eq!(status, 1);

            unsafe {
                #[allow(deprecated)]
                ffcpdt(fp, fp, &raw mut status);
            }
            assert_eq!(status, 1);

            fits_write_hdu(
                f.as_deref_mut().unwrap(),
                core::ptr::null_mut(),
                &mut status,
            );
            assert_eq!(status, 1);

            fits_insert_img(f.as_deref_mut().unwrap(), 16, 1, &naxes, &mut status);
            assert_eq!(status, 1);

            let empty_ttype: Vec<Option<&[c_char]>> = vec![];
            let empty_tform: Vec<&[c_char]> = vec![];
            fits_insert_atbl(
                f.as_deref_mut().unwrap(),
                10,
                1,
                1,
                &empty_ttype,
                None,
                &empty_tform,
                None,
                None,
                &mut status,
            );
            assert_eq!(status, 1);

            fits_insert_btbl(
                f.as_deref_mut().unwrap(),
                1,
                1,
                &empty_ttype,
                &empty_tform,
                None,
                None,
                0,
                &mut status,
            );
            assert_eq!(status, 1);

            fits_delete_hdu(f.as_deref_mut().unwrap(), None, &mut status);
            assert_eq!(status, 1);

            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_same_file_error() {
        // Test SAME_FILE error when copying to same file pointer.
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let path = to_buf(filename);
            let naxes: [c_long; 1] = [10];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &path, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 1, &naxes, &mut status);
            assert_eq!(status, 0);

            let fp: *mut fitsfile = f.as_deref_mut().unwrap();
            unsafe {
                #[allow(deprecated)]
                ffcopy(fp, fp, 0, &raw mut status);
            }
            assert_eq!(status, SAME_FILE);
            status = 0;

            unsafe {
                #[allow(deprecated)]
                ffcphd(fp, fp, &raw mut status);
            }
            assert_eq!(status, SAME_FILE);
            status = 0;

            unsafe {
                #[allow(deprecated)]
                ffcpdt(fp, fp, &raw mut status);
            }
            assert_eq!(status, SAME_FILE);
            status = 0;

            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_insert_image() {
        // Test ffiimg - insert image extension.
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let path = to_buf(filename);
            let naxes: [c_long; 2] = [10, 10];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &path, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            fits_insert_img(f.as_deref_mut().unwrap(), SHORT_IMG, 2, &naxes, &mut status);
            assert_eq!(status, 0);
            let mut numhdus: c_int = 0;
            fits_get_num_hdus(f.as_deref_mut().unwrap(), &mut numhdus, &mut status);
            assert_eq!(numhdus, 2);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_insert_ascii_table() {
        // Test ffitab - insert ASCII table.
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let path = to_buf(filename);
            let ttype = [Some(cc("COL1"))];
            let ttype_ref: Vec<Option<&[c_char]>> = ttype.iter().map(|o| o.as_deref()).collect();
            let tform_v = [cc("F10.2")];
            let tform_ref: Vec<&[c_char]> = tform_v.iter().map(|v| v.as_slice()).collect();
            let tunit = [Some(cc("m"))];
            let tunit_ref: Vec<Option<&[c_char]>> = tunit.iter().map(|o| o.as_deref()).collect();

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &path, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            fits_insert_atbl(
                f.as_deref_mut().unwrap(),
                80,
                1,
                1,
                &ttype_ref,
                None,
                &tform_ref,
                Some(&tunit_ref),
                Some(&cc("TEST")),
                &mut status,
            );
            assert_eq!(status, 0);
            let mut numhdus: c_int = 0;
            fits_get_num_hdus(f.as_deref_mut().unwrap(), &mut numhdus, &mut status);
            assert_eq!(numhdus, 2);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_insert_binary_table() {
        // Test ffibin - insert binary table.
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let path = to_buf(filename);
            let ttype = [Some(cc("COL1")), Some(cc("COL2"))];
            let ttype_ref: Vec<Option<&[c_char]>> = ttype.iter().map(|o| o.as_deref()).collect();
            let tform_v = [cc("1J"), cc("1E")];
            let tform_ref: Vec<&[c_char]> = tform_v.iter().map(|v| v.as_slice()).collect();

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &path, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            fits_insert_btbl(
                f.as_deref_mut().unwrap(),
                5,
                2,
                &ttype_ref,
                &tform_ref,
                None,
                Some(&cc("TABLE")),
                0,
                &mut status,
            );
            assert_eq!(status, 0);
            let mut numhdus: c_int = 0;
            fits_get_num_hdus(f.as_deref_mut().unwrap(), &mut numhdus, &mut status);
            assert_eq!(numhdus, 2);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_delete_hdu() {
        // Test ffdhdu - delete current HDU.
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let path3_str = sibling(filename, "test_edithdu3.fits");
            let path3 = to_buf(&path3_str);
            let ttype = [Some(cc("COL1"))];
            let ttype_ref: Vec<Option<&[c_char]>> = ttype.iter().map(|o| o.as_deref()).collect();
            let tform_v = [cc("1J")];
            let tform_ref: Vec<&[c_char]> = tform_v.iter().map(|v| v.as_slice()).collect();

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &path3, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            fits_create_tbl(
                f.as_deref_mut().unwrap(),
                BINARY_TBL,
                1,
                1,
                &ttype_ref,
                &tform_ref,
                None,
                Some(&cc("TABLE1")),
                &mut status,
            );
            fits_create_tbl(
                f.as_deref_mut().unwrap(),
                BINARY_TBL,
                1,
                1,
                &ttype_ref,
                &tform_ref,
                None,
                Some(&cc("TABLE2")),
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0);

            fits_open_file(&mut f, &path3, READWRITE, &mut status);
            let mut numhdus: c_int = 0;
            fits_get_num_hdus(f.as_deref_mut().unwrap(), &mut numhdus, &mut status);
            assert_eq!(numhdus, 3);

            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut hdutype: c_int = 0;
            fits_delete_hdu(f.as_deref_mut().unwrap(), Some(&mut hdutype), &mut status);
            assert_eq!(status, 0);

            fits_get_num_hdus(f.as_deref_mut().unwrap(), &mut numhdus, &mut status);
            assert_eq!(numhdus, 2);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_bad_naxis() {
        // Test BAD_NAXIS error with too many dimensions.
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let path = to_buf(filename);
            let naxes = [2 as c_long; 200];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &path, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            assert_eq!(status, 0);

            fits_insert_img(
                f.as_deref_mut().unwrap(),
                SHORT_IMG,
                200,
                &naxes,
                &mut status,
            );
            assert_eq!(status, BAD_NAXIS);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_bad_bitpix() {
        // Test BAD_BITPIX error with invalid bitpix value in ffiimgll.
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let path = to_buf(filename);
            let naxes: [c_long; 1] = [10];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &path, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            assert_eq!(status, 0);

            fits_insert_img(f.as_deref_mut().unwrap(), 99, 1, &naxes, &mut status);
            assert_eq!(status, BAD_BITPIX);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_bad_naxes() {
        // Test BAD_NAXES error with negative axis size in ffiimgll.
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let path = to_buf(filename);
            let naxes: [c_long; 1] = [-10];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &path, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            assert_eq!(status, 0);

            fits_insert_img(f.as_deref_mut().unwrap(), SHORT_IMG, 1, &naxes, &mut status);
            assert_eq!(status, BAD_NAXES);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_insert_64bit_image() {
        // Test ffiimg with 64-bit bitpix values in ffiimgll.
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let path = to_buf(filename);
            let naxes: [c_long; 1] = [5];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &path, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);

            fits_insert_img(
                f.as_deref_mut().unwrap(),
                LONGLONG_IMG,
                1,
                &naxes,
                &mut status,
            );
            assert_eq!(status, 0);
            let mut hdunum: c_int = 0;
            fits_get_hdu_num(f.as_deref_mut().unwrap(), &mut hdunum);
            assert_eq!(hdunum, 2);

            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            fits_insert_img(
                f.as_deref_mut().unwrap(),
                DOUBLE_IMG,
                1,
                &naxes,
                &mut status,
            );
            assert_eq!(status, 0);
            fits_get_hdu_num(f.as_deref_mut().unwrap(), &mut hdunum);
            assert_eq!(hdunum, 3);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_copy_with_morekeys() {
        // Test ffcopy with morekeys parameter to reserve header space.
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let path = to_buf(filename);
            let path2_str = sibling(filename, "test_edithdu2.fits");
            let path2 = to_buf(&path2_str);
            let naxes: [c_long; 1] = [10];

            let mut inf: Option<Box<fitsfile>> = None;
            fits_create_file(&mut inf, &path, &mut status);
            fits_write_imghdr(
                inf.as_deref_mut().unwrap(),
                BYTE_IMG,
                1,
                &naxes,
                &mut status,
            );
            fits_close_file(inf.take().unwrap(), &mut status);

            let mut out: Option<Box<fitsfile>> = None;
            fits_create_file(&mut out, &path2, &mut status);
            fits_write_imghdr(out.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            fits_close_file(out.take().unwrap(), &mut status);
            assert_eq!(status, 0);

            fits_open_file(&mut inf, &path, READONLY, &mut status);
            fits_open_file(&mut out, &path2, READWRITE, &mut status);
            fits_copy_hdu(
                inf.as_deref_mut().unwrap(),
                out.as_deref_mut().unwrap(),
                10,
                &mut status,
            );
            assert_eq!(status, 0);
            fits_close_file(inf.take().unwrap(), &mut status);
            fits_close_file(out.take().unwrap(), &mut status);
        });
    }

    /// Copying between two fitsfile handles that share one FITSfile -- the same
    /// file opened twice -- used to abort with `double free detected` and was
    /// ignored for it.
    ///
    /// The second free came from `FITSfile::drop`'s guessed-length fallbacks:
    /// once the registry entry had been consumed, the `else` branch rebuilt a
    /// `Vec` from an assumed length and freed the block again. `free_registered`
    /// has no such fallback -- an unregistered pointer is left alone -- and the
    /// six tile* pointers, which `TILE_STRUCTS` owns, are no longer freed there
    /// at all.
    #[test]
    fn test_copy_within_same_file() {
        // Test copying HDU data within the same file.
        // Open the file twice with different fitsfile pointers.
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let path = to_buf(filename);
            let naxes: [c_long; 1] = [100];
            let mut data = [0u8; 100];
            for (i, d) in data.iter_mut().enumerate() {
                *d = i as u8;
            }

            let mut inf: Option<Box<fitsfile>> = None;
            fits_create_file(&mut inf, &path, &mut status);
            fits_write_imghdr(
                inf.as_deref_mut().unwrap(),
                BYTE_IMG,
                1,
                &naxes,
                &mut status,
            );
            fits_write_img_byt(inf.as_deref_mut().unwrap(), 1, 1, 100, &data, &mut status);
            fits_close_file(inf.take().unwrap(), &mut status);
            assert_eq!(status, 0);

            // Open the same file twice - two fitsfile pointers sharing the Fptr.
            let mut out: Option<Box<fitsfile>> = None;
            fits_open_file(&mut inf, &path, READWRITE, &mut status);
            fits_open_file(&mut out, &path, READWRITE, &mut status);
            fits_movabs_hdu(inf.as_deref_mut().unwrap(), 1, None, &mut status);
            fits_copy_hdu(
                inf.as_deref_mut().unwrap(),
                out.as_deref_mut().unwrap(),
                0,
                &mut status,
            );
            assert_eq!(status, 0);
            let mut numhdus: c_int = 0;
            fits_get_num_hdus(out.as_deref_mut().unwrap(), &mut numhdus, &mut status);
            assert_eq!(numhdus, 2);
            fits_close_file(inf.take().unwrap(), &mut status);
            fits_close_file(out.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_itab_negative_width() {
        // Test ffitab with negative width.
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let path = to_buf(filename);
            let ttype = [Some(cc("COL1"))];
            let ttype_ref: Vec<Option<&[c_char]>> = ttype.iter().map(|o| o.as_deref()).collect();
            let tform_v = [cc("F10.2")];
            let tform_ref: Vec<&[c_char]> = tform_v.iter().map(|v| v.as_slice()).collect();

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &path, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);

            fits_insert_atbl(
                f.as_deref_mut().unwrap(),
                -10,
                1,
                1,
                &ttype_ref,
                None,
                &tform_ref,
                None,
                None,
                &mut status,
            );
            assert_eq!(status, NEG_WIDTH);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_itab_negative_rows() {
        // Test ffitab with negative rows.
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let path = to_buf(filename);
            let ttype = [Some(cc("COL1"))];
            let ttype_ref: Vec<Option<&[c_char]>> = ttype.iter().map(|o| o.as_deref()).collect();
            let tform_v = [cc("F10.2")];
            let tform_ref: Vec<&[c_char]> = tform_v.iter().map(|v| v.as_slice()).collect();

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &path, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);

            fits_insert_atbl(
                f.as_deref_mut().unwrap(),
                80,
                -5,
                1,
                &ttype_ref,
                None,
                &tform_ref,
                None,
                None,
                &mut status,
            );
            assert_eq!(status, NEG_ROWS);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_itab_bad_tfields() {
        // Test ffitab with bad tfields.
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let path = to_buf(filename);
            let empty_ttype: Vec<Option<&[c_char]>> = vec![];
            let empty_tform: Vec<&[c_char]> = vec![];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &path, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);

            fits_insert_atbl(
                f.as_deref_mut().unwrap(),
                80,
                1,
                -5,
                &empty_ttype,
                None,
                &empty_tform,
                None,
                None,
                &mut status,
            );
            assert_eq!(status, BAD_TFIELDS);
            status = 0;

            fits_insert_atbl(
                f.as_deref_mut().unwrap(),
                80,
                1,
                1000,
                &empty_ttype,
                None,
                &empty_tform,
                None,
                None,
                &mut status,
            );
            assert_eq!(status, BAD_TFIELDS);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_ibin_negative_rows() {
        // Test ffibin with negative rows.
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let path = to_buf(filename);
            let ttype = [Some(cc("COL1"))];
            let ttype_ref: Vec<Option<&[c_char]>> = ttype.iter().map(|o| o.as_deref()).collect();
            let tform_v = [cc("1J")];
            let tform_ref: Vec<&[c_char]> = tform_v.iter().map(|v| v.as_slice()).collect();

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &path, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);

            fits_insert_btbl(
                f.as_deref_mut().unwrap(),
                -5,
                1,
                &ttype_ref,
                &tform_ref,
                None,
                None,
                0,
                &mut status,
            );
            assert_eq!(status, NEG_ROWS);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_ibin_bad_tfields() {
        // Test ffibin with bad tfields.
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let path = to_buf(filename);
            let empty_ttype: Vec<Option<&[c_char]>> = vec![];
            let empty_tform: Vec<&[c_char]> = vec![];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &path, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);

            fits_insert_btbl(
                f.as_deref_mut().unwrap(),
                1,
                -5,
                &empty_ttype,
                &empty_tform,
                None,
                None,
                0,
                &mut status,
            );
            assert_eq!(status, BAD_TFIELDS);
            status = 0;

            fits_insert_btbl(
                f.as_deref_mut().unwrap(),
                1,
                1000,
                &empty_ttype,
                &empty_tform,
                None,
                None,
                0,
                &mut status,
            );
            assert_eq!(status, BAD_TFIELDS);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_insert_byte_image() {
        // Test ffiimg with BYTE_IMG bitpix.
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let path = to_buf(filename);
            let naxes: [c_long; 1] = [5];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &path, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), SHORT_IMG, 0, &[], &mut status);

            fits_insert_img(f.as_deref_mut().unwrap(), BYTE_IMG, 1, &naxes, &mut status);
            assert_eq!(status, 0);
            let mut hdunum: c_int = 0;
            fits_get_hdu_num(f.as_deref_mut().unwrap(), &mut hdunum);
            assert_eq!(hdunum, 2);
            let mut bitpix: c_int = 0;
            fits_read_key(
                f.as_deref_mut().unwrap(),
                KeywordDatatypeMut::TINT(&mut bitpix),
                &cc("BITPIX"),
                None,
                &mut status,
            );
            assert_eq!(bitpix, 8);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_insert_naxis_zero() {
        // Test ffiimg with naxis=0 (null image).
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let path = to_buf(filename);

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &path, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);

            fits_insert_img(f.as_deref_mut().unwrap(), SHORT_IMG, 0, &[], &mut status);
            assert_eq!(status, 0);
            let mut hdunum: c_int = 0;
            fits_get_hdu_num(f.as_deref_mut().unwrap(), &mut hdunum);
            assert_eq!(hdunum, 2);
            let mut naxis: c_int = 0;
            fits_read_key(
                f.as_deref_mut().unwrap(),
                KeywordDatatypeMut::TINT(&mut naxis),
                &cc("NAXIS"),
                None,
                &mut status,
            );
            assert_eq!(naxis, 0);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_copy_with_large_header_space() {
        // Test ffcopy with source having large empty header space (ffwend path).
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let path = to_buf(filename);
            let path2_str = sibling(filename, "test_edithdu2.fits");
            let path2 = to_buf(&path2_str);
            let naxes: [c_long; 1] = [10];

            let mut inf: Option<Box<fitsfile>> = None;
            fits_create_file(&mut inf, &path, &mut status);
            fits_write_imghdr(
                inf.as_deref_mut().unwrap(),
                BYTE_IMG,
                1,
                &naxes,
                &mut status,
            );
            // Reserve 100 additional keywords - creates multiple empty blocks.
            fits_set_hdrsize(inf.as_deref_mut().unwrap(), 100, &mut status);
            fits_close_file(inf.take().unwrap(), &mut status);
            assert_eq!(status, 0);

            let mut out: Option<Box<fitsfile>> = None;
            fits_create_file(&mut out, &path2, &mut status);
            fits_write_imghdr(out.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            fits_close_file(out.take().unwrap(), &mut status);
            assert_eq!(status, 0);

            fits_open_file(&mut inf, &path, READONLY, &mut status);
            fits_open_file(&mut out, &path2, READWRITE, &mut status);
            fits_copy_hdu(
                inf.as_deref_mut().unwrap(),
                out.as_deref_mut().unwrap(),
                0,
                &mut status,
            );
            assert_eq!(status, 0);
            fits_close_file(inf.take().unwrap(), &mut status);
            fits_close_file(out.take().unwrap(), &mut status);

            fits_open_file(&mut out, &path2, READONLY, &mut status);
            let mut numhdus: c_int = 0;
            fits_get_num_hdus(out.as_deref_mut().unwrap(), &mut numhdus, &mut status);
            assert!(numhdus >= 2);
            fits_close_file(out.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_negative_naxis() {
        // Test ffiimg with negative naxis.
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let path = to_buf(filename);
            let naxes: [c_long; 1] = [10];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &path, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);

            fits_insert_img(
                f.as_deref_mut().unwrap(),
                SHORT_IMG,
                -1,
                &naxes,
                &mut status,
            );
            assert_eq!(status, BAD_NAXIS);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_naxis_over_limit() {
        // Test ffiimg with naxis > 999.
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let path = to_buf(filename);
            let naxes = [1 as c_long; 1000];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &path, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);

            fits_insert_img(
                f.as_deref_mut().unwrap(),
                SHORT_IMG,
                1000,
                &naxes,
                &mut status,
            );
            assert_eq!(status, BAD_NAXIS);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_cpfl_same_file_error() {
        // Test ffcpfl with same file (should fail).
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let path = to_buf(filename);
            let naxes: [c_long; 1] = [10];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &path, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 1, &naxes, &mut status);
            assert_eq!(status, 0);

            let fp: *mut fitsfile = f.as_deref_mut().unwrap();
            unsafe {
                #[allow(deprecated)]
                ffcpfl(fp, fp, 1, 1, 1, &raw mut status);
            }
            assert_eq!(status, SAME_FILE);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_cpfl_status_error() {
        // Test ffcpfl with initial error status.
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let path = to_buf(filename);
            let naxes: [c_long; 1] = [10];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &path, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 1, &naxes, &mut status);
            assert_eq!(status, 0);

            let fp: *mut fitsfile = f.as_deref_mut().unwrap();
            status = 1;
            unsafe {
                #[allow(deprecated)]
                ffcpfl(fp, fp, 1, 1, 1, &raw mut status);
            }
            assert_eq!(status, 1);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_insert_image_readonly() {
        // Test ffiimg on readonly file.
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let path = to_buf(filename);
            let naxes: [c_long; 1] = [10];
            let empty_ttype: Vec<Option<&[c_char]>> = vec![];
            let empty_tform: Vec<&[c_char]> = vec![];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &path, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 1, &naxes, &mut status);
            fits_create_tbl(
                f.as_deref_mut().unwrap(),
                BINARY_TBL,
                1,
                0,
                &empty_ttype,
                &empty_tform,
                None,
                None,
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0);

            fits_open_file(&mut f, &path, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 1, None, &mut status);
            fits_insert_img(f.as_deref_mut().unwrap(), SHORT_IMG, 1, &naxes, &mut status);
            assert_eq!(status, READONLY_FILE);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_insert_ascii_table_readonly() {
        // Test ffitab on readonly file.
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let path = to_buf(filename);
            let naxes: [c_long; 1] = [10];
            let ttype = [Some(cc("COL1"))];
            let ttype_ref: Vec<Option<&[c_char]>> = ttype.iter().map(|o| o.as_deref()).collect();
            let tform_v = [cc("F10.2")];
            let tform_ref: Vec<&[c_char]> = tform_v.iter().map(|v| v.as_slice()).collect();
            let empty_ttype: Vec<Option<&[c_char]>> = vec![];
            let empty_tform: Vec<&[c_char]> = vec![];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &path, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 1, &naxes, &mut status);
            fits_create_tbl(
                f.as_deref_mut().unwrap(),
                BINARY_TBL,
                1,
                0,
                &empty_ttype,
                &empty_tform,
                None,
                None,
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0);

            fits_open_file(&mut f, &path, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 1, None, &mut status);
            fits_insert_atbl(
                f.as_deref_mut().unwrap(),
                80,
                1,
                1,
                &ttype_ref,
                None,
                &tform_ref,
                None,
                None,
                &mut status,
            );
            assert_eq!(status, READONLY_FILE);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_insert_binary_table_readonly() {
        // Test ffibin on readonly file.
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let path = to_buf(filename);
            let naxes: [c_long; 1] = [10];
            let ttype = [Some(cc("COL1"))];
            let ttype_ref: Vec<Option<&[c_char]>> = ttype.iter().map(|o| o.as_deref()).collect();
            let tform_v = [cc("1J")];
            let tform_ref: Vec<&[c_char]> = tform_v.iter().map(|v| v.as_slice()).collect();
            let empty_ttype: Vec<Option<&[c_char]>> = vec![];
            let empty_tform: Vec<&[c_char]> = vec![];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &path, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 1, &naxes, &mut status);
            fits_create_tbl(
                f.as_deref_mut().unwrap(),
                BINARY_TBL,
                1,
                0,
                &empty_ttype,
                &empty_tform,
                None,
                None,
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0);

            fits_open_file(&mut f, &path, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 1, None, &mut status);
            fits_insert_btbl(
                f.as_deref_mut().unwrap(),
                1,
                1,
                &ttype_ref,
                &tform_ref,
                None,
                None,
                0,
                &mut status,
            );
            assert_eq!(status, READONLY_FILE);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_cpht_status_error() {
        // Test ffcpht with status already set.
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let path = to_buf(filename);
            let path2_str = sibling(filename, "test_edithdu2.fits");
            let path2 = to_buf(&path2_str);
            let naxes: [c_long; 1] = [10];

            let mut inf: Option<Box<fitsfile>> = None;
            fits_create_file(&mut inf, &path, &mut status);
            fits_write_imghdr(
                inf.as_deref_mut().unwrap(),
                BYTE_IMG,
                1,
                &naxes,
                &mut status,
            );
            fits_close_file(inf.take().unwrap(), &mut status);

            let mut out: Option<Box<fitsfile>> = None;
            fits_create_file(&mut out, &path2, &mut status);
            fits_write_imghdr(out.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            fits_close_file(out.take().unwrap(), &mut status);
            assert_eq!(status, 0);

            fits_open_file(&mut inf, &path, READONLY, &mut status);
            fits_open_file(&mut out, &path2, READWRITE, &mut status);

            status = 1;
            fits_copy_hdutab(
                inf.as_deref_mut().unwrap(),
                out.as_deref_mut().unwrap(),
                1,
                1,
                &mut status,
            );
            assert_eq!(status, 1);
            status = 0;
            fits_close_file(inf.take().unwrap(), &mut status);
            fits_close_file(out.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_copy_table_to_empty_primary() {
        // Test copying a table to an empty output file (create dummy primary).
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let path = to_buf(filename);
            let ttype = [Some(cc("COL1"))];
            let ttype_ref: Vec<Option<&[c_char]>> = ttype.iter().map(|o| o.as_deref()).collect();
            let tform_v = [cc("1J")];
            let tform_ref: Vec<&[c_char]> = tform_v.iter().map(|v| v.as_slice()).collect();

            let mut inf: Option<Box<fitsfile>> = None;
            fits_create_file(&mut inf, &path, &mut status);
            fits_write_imghdr(inf.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            fits_create_tbl(
                inf.as_deref_mut().unwrap(),
                BINARY_TBL,
                3,
                1,
                &ttype_ref,
                &tform_ref,
                None,
                Some(&cc("DATA")),
                &mut status,
            );
            fits_close_file(inf.take().unwrap(), &mut status);
            assert_eq!(status, 0);

            // Memory file to avoid automatic primary array creation.
            let mut out: Option<Box<fitsfile>> = None;
            let memname = to_buf("mem://test_empty");
            fits_create_file(&mut out, &memname, &mut status);
            assert_eq!(status, 0);

            fits_open_file(&mut inf, &path, READONLY, &mut status);
            fits_movabs_hdu(inf.as_deref_mut().unwrap(), 2, None, &mut status);
            fits_copy_header(
                inf.as_deref_mut().unwrap(),
                out.as_deref_mut().unwrap(),
                &mut status,
            );
            assert_eq!(status, 0);
            fits_close_file(inf.take().unwrap(), &mut status);
            fits_close_file(out.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_iimgll_status_error() {
        // Test ffiimgll with status already set.
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let path = to_buf(filename);
            let llnaxes: [LONGLONG; 1] = [10];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &path, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);

            status = 1;
            fits_insert_imgll(
                f.as_deref_mut().unwrap(),
                SHORT_IMG,
                1,
                &llnaxes,
                &mut status,
            );
            assert_eq!(status, 1);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_delete_last_hdu() {
        // Test ffdhdu on last HDU making primary empty.
        with_temp_file(|filename| {
            let mut status: c_int = 0;
            let path3_str = sibling(filename, "test_edithdu3.fits");
            let path3 = to_buf(&path3_str);
            let naxes: [c_long; 1] = [10];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &path3, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 1, &naxes, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0);

            fits_open_file(&mut f, &path3, READWRITE, &mut status);
            let mut numhdus: c_int = 0;
            fits_get_num_hdus(f.as_deref_mut().unwrap(), &mut numhdus, &mut status);
            assert_eq!(numhdus, 1);

            fits_movabs_hdu(f.as_deref_mut().unwrap(), 1, None, &mut status);
            fits_delete_hdu(f.as_deref_mut().unwrap(), None, &mut status);
            // After deleting the only HDU, should get a minimal primary.
            assert_eq!(status, 0);

            fits_get_num_hdus(f.as_deref_mut().unwrap(), &mut numhdus, &mut status);
            assert_eq!(numhdus, 1);
            let mut naxis: c_int = 0;
            fits_read_key(
                f.as_deref_mut().unwrap(),
                KeywordDatatypeMut::TINT(&mut naxis),
                &cc("NAXIS"),
                None,
                &mut status,
            );
            assert_eq!(naxis, 0); // Primary should now be empty.
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }
}
