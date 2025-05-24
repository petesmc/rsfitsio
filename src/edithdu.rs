use crate::c_types::{FILE, c_char, c_int, c_long};
use crate::{
    BL, TKeywords,
    buffers::*,
    cs,
    editcol::ffcprw_safe,
    fitscore::{
        ffbnfm_safe, ffcmsg_safe, ffcrhd_safer, ffdblk, ffgabc_safe, ffgext, ffghadll_safe,
        ffghdn_safe, ffgidm_safe, ffhdef_safe, ffiblk, ffkeyn_safe, ffmahd_safe, ffpdfl,
        ffpmsg_slice, ffpmsg_str, ffrdef_safe, ffrhdu_safer, ffwend,
    },
    fitsio::*,
    fitsio2::*,
    getkey::{ffgcrd_safe, ffghsp_safe, ffgkyj_safe, ffgrec_safe},
    int_snprintf,
    modkey::{ffdkey_safe, ffikyj_safe, ffukyj_safe},
    nullable_slice_cstr,
    putkey::{
        ffcrim_safer, ffcrimll_safer, ffcrtb_safer, ffphbn_safe, ffphpr_safe, ffphprll_safe,
        ffphtb_safe, ffpkyj_safe, ffpkyl_safe, ffpkys_safe, ffprec_safe,
    },
    wrappers::{strcpy_safe, strncat_safe},
};
use bytemuck::{cast_slice, cast_slice_mut};
use core::slice;
use cstr::cstr;
use libc::fwrite;
use std::{
    cmp,
    ffi::CStr,
    ptr::{self},
};

/*--------------------------------------------------------------------------*/
/// copy the CHDU from infptr to the CHDU of outfptr.
/// This will also allocate space in the output header for MOREKY keywords
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffcopy(
    infptr: *mut fitsfile,  /* I - FITS file pointer to input file  */
    outfptr: *mut fitsfile, /* I - FITS file pointer to output file */
    morekeys: c_int,        /* I - reserve space in output header   */
    status: *mut c_int,     /* IO - error status     */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);

        if std::ptr::eq(infptr, outfptr) {
            *status = SAME_FILE;
            return *status;
        }

        let infptr = infptr.as_mut().expect(NULL_MSG);
        let outfptr = outfptr.as_mut().expect(NULL_MSG);

        ffcopy_safer(infptr, outfptr, morekeys, status)
    }
}

/*--------------------------------------------------------------------------*/
/// copy the CHDU from infptr to the CHDU of outfptr.
/// This will also allocate space in the output header for MOREKY keywords
pub unsafe fn ffcopy_safer(
    infptr: &mut fitsfile,  /* I - FITS file pointer to input file  */
    outfptr: &mut fitsfile, /* I - FITS file pointer to output file */
    morekeys: c_int,        /* I - reserve space in output header   */
    status: &mut c_int,     /* IO - error status     */
) -> c_int {
    unsafe {
        let mut nspace: c_int = 0;

        if *status > 0 {
            return *status;
        }

        if ffcphd_safer(infptr, outfptr, status) > 0 {
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
}

/*--------------------------------------------------------------------------*/
/// copy all or part of the input file to the output file.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffcpfl(
    infptr: *mut fitsfile,  /* I - FITS file pointer to input file  */
    outfptr: *mut fitsfile, /* I - FITS file pointer to output file */
    previous: c_int,        /* I - copy any previous HDUs?   */
    current: c_int,         /* I - copy the current HDU?     */
    following: c_int,       /* I - copy any following HDUs?   */
    status: *mut c_int,     /* IO - error status     */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);

        if *status > 0 {
            return *status;
        }

        if std::ptr::eq(infptr, outfptr) {
            *status = SAME_FILE;
            return *status;
        }

        let infptr = infptr.as_mut().expect(NULL_MSG);
        let outfptr = outfptr.as_mut().expect(NULL_MSG);

        ffcpfl_safer(infptr, outfptr, previous, current, following, status)
    }
}

/*--------------------------------------------------------------------------*/
/// copy all or part of the input file to the output file.
pub unsafe fn ffcpfl_safer(
    infptr: &mut fitsfile,  /* I - FITS file pointer to input file  */
    outfptr: &mut fitsfile, /* I - FITS file pointer to output file */
    previous: c_int,        /* I - copy any previous HDUs?   */
    current: c_int,         /* I - copy the current HDU?     */
    following: c_int,       /* I - copy any following HDUs?   */
    status: &mut c_int,     /* IO - error status     */
) -> c_int {
    unsafe {
        let mut hdunum = 0;

        if *status > 0 {
            return *status;
        }

        ffghdn_safe(infptr, &mut hdunum);

        if previous != 0 {
            /* copy any previous HDUs */
            for ii in 1..(hdunum) {
                ffmahd_safe(infptr, ii, None, status);
                ffcopy_safer(infptr, outfptr, 0, status);
            }
        }

        if current != 0 && (*status <= 0) {
            /* copy current HDU */
            ffmahd_safe(infptr, hdunum, None, status);
            ffcopy_safer(infptr, outfptr, 0, status);
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

                if ffcopy_safer(infptr, outfptr, 0, status) != 0 {
                    break; /* quit on unexpected error */
                }

                ii += 1;
            }
        }

        ffmahd_safe(infptr, hdunum, None, status); /* restore initial position */
        *status
    }
}

/*--------------------------------------------------------------------------*/
/// copy the header keywords from infptr to outfptr.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffcphd(
    infptr: *mut fitsfile,  /* I - FITS file pointer to input file  */
    outfptr: *mut fitsfile, /* I - FITS file pointer to output file */
    status: *mut c_int,     /* IO - error status     */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);

        if *status > 0 {
            return *status;
        }

        if std::ptr::eq(infptr, outfptr) {
            *status = SAME_FILE;
            return *status;
        }

        let infptr = infptr.as_mut().expect(NULL_MSG);
        let outfptr = outfptr.as_mut().expect(NULL_MSG);

        ffcphd_safer(infptr, outfptr, status)
    }
}

/*--------------------------------------------------------------------------*/
/// copy the header keywords from infptr to outfptr.
pub(crate) unsafe fn ffcphd_safer(
    infptr: &mut fitsfile,  /* I - FITS file pointer to input file  */
    outfptr: &mut fitsfile, /* I - FITS file pointer to output file */
    status: &mut c_int,     /* IO - error status     */
) -> c_int {
    unsafe {
        let mut nkeys: c_int = 0;
        let ii: c_int = 0;
        let mut inPrim: c_int = 0;
        let mut outPrim: c_int = 0;
        let mut naxis: c_long = 0;
        let naxes: [c_long; 1] = [0; 1];
        let card: *mut c_char = ptr::null_mut();
        let mut comm: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];
        let tmpbuff2: *mut c_char = ptr::null_mut();

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
            ffgkyj_safe(infptr, cs!("NAXIS"), &mut naxis, None, status);
        }
        /* set the output pointer to the correct HDU */
        if outfptr.HDUposition != outfptr.Fptr.curhdu {
            ffmahd_safe(outfptr, (outfptr.HDUposition) + 1, None, status);
        }
        /* check if output header is empty; if not create new empty HDU */
        let headstart = outfptr.Fptr.get_headstart_as_slice()[outfptr.Fptr.curhdu as usize];
        if outfptr.Fptr.headend != headstart {
            ffcrhd_safer(outfptr, status);
        }
        if outfptr.HDUposition == 0 {
            if naxis < 0 {
                /* the input HDU is a table, so we have to create */
                /* a dummy Primary array before copying it to the output */
                ffcrim_safer(outfptr, 8, 0, &naxes, status);
                ffcrhd_safer(outfptr, status); /* create new empty HDU */
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
            strcpy_safe(&mut comm, cs!("IMAGE extension"));
            ffpkys_safe(outfptr, cs!("XTENSION"), cs!("IMAGE"), Some(&comm), status);

            /* copy BITPIX through NAXISn keywords */
            for ii in 1..(3 + naxis as usize) {
                let card = &tmpbuff2[(ii * FLEN_CARD)..];
                ffprec_safe(outfptr, card, status);
            }

            strcpy_safe(&mut comm, cs!("number of random group parameters"));
            ffpkyj_safe(outfptr, cs!("PCOUNT"), 0, Some(&comm), status);

            strcpy_safe(&mut comm, cs!("number of random groups"));
            ffpkyj_safe(outfptr, cs!("GCOUNT"), 1, Some(&comm), status);

            /* copy remaining keywords, excluding EXTEND, and reference COMMENT keywords */
            for ii in (3 + naxis as usize)..(nkeys as usize) {
                let card = &tmpbuff2[(ii * FLEN_CARD)..];
                if FSTRNCMP(card, cs!("EXTEND  "), 8) > 0
                    && FSTRNCMP(
                        card,
                        cs!("COMMENT   FITS (Flexible Image Transport System) format is"),
                        58,
                    ) > 0
                    && FSTRNCMP(
                        card,
                        cs!("COMMENT   and Astrophysics', volume 376, page 3"),
                        47,
                    ) > 0
                {
                    ffprec_safe(outfptr, card, status);
                }
            }
        } else if inPrim == 0 && outPrim == 1 {
            /* copying between image extension and primary array */
            strcpy_safe(&mut comm, cs!("file does conform to FITS standard"));
            ffpkyl_safe(outfptr, cs!("SIMPLE"), TRUE as c_int, Some(&comm), status);

            /* copy BITPIX through NAXISn keywords */
            for ii in 1..(3 + naxis as usize) {
                let card = &tmpbuff2[(ii * FLEN_CARD)..];
                ffprec_safe(outfptr, card, status);
            }

            /* add the EXTEND keyword */
            strcpy_safe(&mut comm, cs!("FITS dataset may contain extensions"));
            ffpkyl_safe(outfptr, cs!("EXTEND"), TRUE as c_int, Some(&comm), status);

            /* write standard block of self-documentating comments */
            ffprec_safe(
                outfptr,
                cs!(
                    "COMMENT   FITS (Flexible Image Transport System) format is defined in 'Astronomy"
                ),
                status,
            );
            ffprec_safe(
                outfptr,
                cs!(
                    "COMMENT   and Astrophysics', volume 376, page 359; bibcode: 2001A&A...376..359H"
                ),
                status,
            );

            /* copy remaining keywords, excluding pcount, gcount */
            for ii in (3 + naxis as usize)..(nkeys as usize) {
                let card = &tmpbuff2[(ii * FLEN_CARD)..];
                if FSTRNCMP(card, cs!("PCOUNT  "), 8) > 0 && FSTRNCMP(card, cs!("GCOUNT  "), 8) > 0
                {
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
}

/*--------------------------------------------------------------------------*/
/// Copy the table structure from an existing table HDU, but only
/// copy a limited row range.  All header keywords from the input
/// table are copied directly, but NAXSI2 and PCOUNT are set to their
/// correct values.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffcpht(
    infptr: *mut fitsfile,  /* I - FITS file pointer to input file  */
    outfptr: *mut fitsfile, /* I - FITS file pointer to output file */
    firstrow: LONGLONG,     /* I - number of first row to copy (1 based)  */
    nrows: LONGLONG,        /* I - number of rows to copy  */
    status: *mut c_int,     /* IO - error status     */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let infptr = infptr.as_mut().expect(NULL_MSG);
        let outfptr = outfptr.as_mut().expect(NULL_MSG);

        ffcpht_safer(infptr, outfptr, firstrow, nrows, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Copy the table structure from an existing table HDU, but only
/// copy a limited row range.  All header keywords from the input
/// table are copied directly, but NAXSI2 and PCOUNT are set to their
/// correct values.
pub unsafe fn ffcpht_safer(
    infptr: &mut fitsfile,  /* I - FITS file pointer to input file  */
    outfptr: &mut fitsfile, /* I - FITS file pointer to output file */
    firstrow: LONGLONG,     /* I - number of first row to copy (1 based)  */
    nrows: LONGLONG,        /* I - number of rows to copy  */
    status: &mut c_int,     /* IO - error status     */
) -> c_int {
    unsafe {
        if *status > 0 {
            return *status;
        }

        /* Copy the header only */
        ffcphd_safer(infptr, outfptr, status);
        /* Note that we now have a copied header that describes the table,
        and that is the current header, but the original number of table
        rows and heap area sizes are still there. */

        /* Zero out the size-related keywords */
        if *status == 0 {
            ffukyj_safe(outfptr, cs!("NAXIS2"), 0, None, status); /* NAXIS2 = 0 */
            ffukyj_safe(outfptr, cs!("PCOUNT"), 0, None, status); /* PCOUNT = 0 */
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
}

/*--------------------------------------------------------------------------*/
/// copy the data unit from the CHDU of infptr to the CHDU of outfptr.
/// This will overwrite any data already in the outfptr CHDU.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffcpdt(
    infptr: *mut fitsfile,  /* I - FITS file pointer to input file  */
    outfptr: *mut fitsfile, /* I - FITS file pointer to output file */
    status: *mut c_int,     /* IO - error status     */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);

        if *status > 0 {
            return *status;
        }

        if std::ptr::eq(infptr, outfptr) {
            *status = SAME_FILE;
            return *status;
        }

        let infptr = infptr.as_mut().expect(NULL_MSG);
        let outfptr = outfptr.as_mut().expect(NULL_MSG);

        ffcpdt_safe(infptr, outfptr, status)
    }
}

/*--------------------------------------------------------------------------*/
/// copy the data unit from the CHDU of infptr to the CHDU of outfptr.
/// This will overwrite any data already in the outfptr CHDU.
pub fn ffcpdt_safe(
    infptr: &mut fitsfile,  /* I - FITS file pointer to input file  */
    outfptr: &mut fitsfile, /* I - FITS file pointer to output file */
    status: &mut c_int,     /* IO - error status     */
) -> c_int {
    let mut nb: c_long = 0;
    let ii: c_long = 0;
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
        if std::ptr::eq(infptr.Fptr.as_mut(), outfptr.Fptr.as_mut()) {
            /* copying between 2 HDUs in the SAME file */
            unreachable!(
                "Above ptr comparison prevents us from landing here. Matching original code"
            );
            for ii in 0..(nb as usize) {
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

            for ii in 0..(nb as usize) {
                ffgbyt(infptr, IOBUFLEN, cast_slice_mut(&mut buffer), status); /* read input block */
                ffpbyt(outfptr, IOBUFLEN, cast_slice_mut(&mut buffer), status); /* write output block */
            }
        }
    }
    *status
}

/*--------------------------------------------------------------------------*/
/// write the data unit from the CHDU of infptr to the output file stream
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffwrhdu(
    infptr: *mut fitsfile, /* I - FITS file pointer to input file  */
    outstream: *mut FILE,  /* I - stream to write HDU to */
    status: *mut c_int,    /* IO - error status     */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let infptr = infptr.as_mut().expect(NULL_MSG);

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

            for ii in 0..(nb as usize) {
                ffgbyt(infptr, BL!(), cast_slice_mut(&mut buffer), status); /* read input block */
                fwrite(buffer.as_ptr() as *mut _, 1, BL!(), outstream); /* write to output stream */
            }
        }
        *status
    }
}

/*--------------------------------------------------------------------------*/
/// insert an IMAGE extension following the current HDU
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffiimg(
    fptr: *mut fitsfile,  /* I - FITS file pointer           */
    bitpix: c_int,        /* I - bits per pixel              */
    naxis: c_int,         /* I - number of axes in the array */
    naxes: *const c_long, /* I - size of each axis           */
    status: *mut c_int,   /* IO - error status               */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        let naxes = slice::from_raw_parts(naxes, naxis as usize);

        ffiimg_safer(fptr, bitpix, naxis, naxes, status)
    }
}

/*--------------------------------------------------------------------------*/
/// insert an IMAGE extension following the current HDU
pub unsafe fn ffiimg_safer(
    fptr: &mut fitsfile, /* I - FITS file pointer           */
    bitpix: c_int,       /* I - bits per pixel              */
    naxis: c_int,        /* I - number of axes in the array */
    naxes: &[c_long],    /* I - size of each axis           */
    status: &mut c_int,  /* IO - error status               */
) -> c_int {
    unsafe {
        let mut tnaxes: [LONGLONG; 99] = [0; 99];

        if *status > 0 {
            return *status;
        }

        if naxis > 99 {
            ffpmsg_str("NAXIS value is too large (>99)  (ffiimg)");
            *status = BAD_NAXIS;
            return *status;
        }

        for i in 0..(naxis as usize) {
            tnaxes[i] = naxes[i] as LONGLONG;
        }

        ffiimgll_safer(fptr, bitpix, naxis, &tnaxes, status);

        *status
    }
}

/*--------------------------------------------------------------------------*/
/// insert an IMAGE extension following the current HDU
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffiimgll(
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

        ffiimgll_safer(fptr, bitpix, naxis, naxes, status)
    }
}

/*--------------------------------------------------------------------------*/
/// insert an IMAGE extension following the current HDU
pub unsafe fn ffiimgll_safer(
    fptr: &mut fitsfile, /* I - FITS file pointer           */
    bitpix: c_int,       /* I - bits per pixel              */
    naxis: c_int,        /* I - number of axes in the array */
    naxes: &[LONGLONG],  /* I - size of each axis           */
    status: &mut c_int,  /* IO - error status               */
) -> c_int {
    unsafe {
        let mut bytlen: c_int = 0;
        let mut nexthdu: c_int = 0;
        let mut maxhdu: c_int = 0;
        let ii: c_int = 0;
        let mut onaxis: c_int = 0;
        let mut nblocks: c_long = 0;
        let mut npixels: LONGLONG = 0;
        let mut newstart: LONGLONG = 0;
        let mut datasize: LONGLONG = 0;
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
                ffcrimll_safer(fptr, bitpix, naxis, naxes, status);
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

        datasize = npixels * bytlen as LONGLONG; /* size of image in bytes */
        nblocks = (((datasize + (BL!() - 1)) / BL!()) + 1) as c_long; /* +1 for the header */

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
                ffkeyn_safe(cs!("NAXIS"), onaxis, &mut naxiskey, status);
            } else {
                strcpy_safe(&mut naxiskey, cs!("NAXIS"));
            }

            ffgcrd_safe(fptr, &naxiskey, &mut card, status); /* read last NAXIS keyword */

            ffikyj_safe(
                fptr,
                cs!("PCOUNT"),
                0,
                Some(cs!("required keyword")),
                status,
            ); /* add PCOUNT and */
            ffikyj_safe(
                fptr,
                cs!("GCOUNT"),
                1,
                Some(cs!("required keyword")),
                status,
            ); /* GCOUNT keywords */

            if *status > 0 {
                return *status;
            }

            if ffdkey_safe(fptr, cs!("EXTEND"), status) != 0 {
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
}

/*--------------------------------------------------------------------------*/
/// insert an ASCII table extension following the current HDU
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffitab(
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
        let mut nexthdu: c_int = 0;
        let mut maxhdu: c_int = 0;
        let ii: c_int = 0;
        let mut nunit: c_int = 0;
        let mut nhead: c_int = 0;
        let ncols: c_int = 0;
        let gotmem: c_int = 0;
        let mut nblocks: c_long = 0;
        let mut rowlen: c_long = 0;
        let mut datasize: LONGLONG = 0;
        let mut newstart: LONGLONG = 0;
        let mut errmsg: [c_char; FLEN_ERRMSG] = [0; FLEN_ERRMSG];
        let mut extnm: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];

        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        if *status > 0 {
            return *status;
        }

        nullable_slice_cstr!(extnmx);

        let tkeywords = TKeywords::new(tfields, ttype, tform, tunit);
        let (v_ttype, v_tform, v_tunit) = tkeywords.tkeywords_to_vecs();

        extnm[0] = 0;
        if let Some(extnmx) = extnmx {
            strncat_safe(&mut extnm, extnmx, FLEN_VALUE - 1);
        }

        if fptr.HDUposition != fptr.Fptr.curhdu {
            ffmahd_safe(fptr, (fptr.HDUposition) + 1, None, status);
        }

        maxhdu = fptr.Fptr.maxhdu;

        let headstart = fptr.Fptr.get_headstart_as_slice();
        /* if the current header is completely empty or, if we are at the end of the file, ...  */
        if (fptr.Fptr.headend == headstart[fptr.Fptr.curhdu as usize])
            || (((fptr.Fptr.curhdu) == maxhdu)
                && (headstart[(maxhdu + 1) as usize] >= fptr.Fptr.logfilesize))
        {
            /* then simply append new image extension */
            ffcrtb_safer(
                fptr,
                ASCII_TBL,
                naxis2,
                tfields,
                &v_ttype,
                &v_tform,
                v_tunit.as_deref(),
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
        nunit = 0;

        for ii in 0..(tfields as usize) {
            if let Some(v_tunit) = v_tunit.as_deref()
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

        let mut tbcol = match tbcol.is_null() || (naxis1 == 0 && tfields != 0) {
            true => {
                gotmem = true;
                vec![0 as c_long; ncols]
            }
            false => {
                let x = slice::from_raw_parts(tbcol, ncols).to_vec();
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
            ffgabc_safe(tfields, &v_tform, 1, &mut rowlen, &mut tbcol, status);
        }

        nhead = (9 + (3 * tfields) + nunit + 35) / 36; /* no. of header blocks */
        datasize = (rowlen as LONGLONG) * naxis2; /* size of table in bytes */
        nblocks = (((datasize + (BL!() - 1)) / BL!()) + nhead as LONGLONG) as c_long; /* size of HDU */

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

        nexthdu = (fptr.Fptr.curhdu) + 1; /* number of the next (new) hdu */
        newstart = headstart[nexthdu as usize]; /* save starting addr of HDU */

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
        fptr.Fptr.datastart = (hs_item) + (nhead as LONGLONG * BL!());
        fptr.Fptr.hdutype = ASCII_TBL; /* might need to be reset... */

        /* write the required header keywords */

        ffphtb_safe(
            fptr,
            rowlen as LONGLONG,
            naxis2,
            tfields,
            &v_ttype,
            Some(&tbcol),
            &v_tform,
            v_tunit.as_deref(),
            Some(&extnm),
            status,
        );

        /* redefine internal structure for this HDU */
        ffrdef_safe(fptr, status);
        *status
    }
}

/*--------------------------------------------------------------------------*/
/// insert a Binary table extension following the current HDU
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffibin(
    fptr: *mut fitsfile,         /* I - FITS file pointer                        */
    naxis2: LONGLONG,            /* I - number of rows in the table              */
    tfields: c_int,              /* I - number of columns in the table           */
    ttype: *const *const c_char, /* I - name of each column                      */
    tform: *const *const c_char, /* I - value of TFORMn keyword for each column  */
    tunit: *const *const c_char, /* I - value of TUNITn keyword for each column  */
    extnmx: *const c_char,       /* I - value of EXTNAME keyword, if any         */
    pcount: LONGLONG,            /* I - size of special data area (heap)         */
    status: *mut c_int,          /* IO - error status                            */
) -> c_int {
    unsafe {
        let mut nexthdu: c_int = 0;
        let mut maxhdu: c_int = 0;
        let ii: c_int = 0;
        let mut nunit: c_int = 0;
        let mut nhead: c_int = 0;
        let mut datacode: c_int = 0;
        let mut naxis1: LONGLONG = 0;
        let mut nblocks: c_long = 0;
        let mut repeat: c_long = 0;
        let mut width: c_long = 0;
        let mut datasize: LONGLONG = 0;
        let mut newstart: LONGLONG = 0;
        let mut errmsg: [c_char; FLEN_ERRMSG] = [0; FLEN_ERRMSG];
        let mut extnm: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];

        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        nullable_slice_cstr!(extnmx);

        let tkeywords = TKeywords::new(tfields, ttype, tform, tunit);
        let (v_ttype, v_tform, v_tunit) = tkeywords.tkeywords_to_vecs();

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

        maxhdu = fptr.Fptr.maxhdu;

        let headstart = fptr.Fptr.get_headstart_as_slice();
        /* if the current header is completely empty ...  */
        if ( fptr.Fptr.headend == headstart[fptr.Fptr.curhdu as usize] )
        /* or, if we are at the end of the file, ... */
    ||  ( ((fptr.Fptr.curhdu) == maxhdu ) &&
       (headstart[(maxhdu + 1) as usize] >= fptr.Fptr.logfilesize ) )
        {
            /* then simply append new image extension */
            ffcrtb_safer(
                fptr,
                BINARY_TBL,
                naxis2,
                tfields,
                &v_ttype,
                &v_tform,
                v_tunit.as_deref(),
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
        nunit = 0;
        for ii in 0..(tfields as usize) {
            if let Some(v_tunit) = v_tunit.as_deref()
                && v_tunit[ii].is_some()
            {
                nunit += 1;
            }
        }

        if extnm[0] != 0 {
            nunit += 1; /* add one for the EXTNAME keyword */
        }

        nhead = (9 + (2 * tfields) + nunit + 35) / 36; /* no. of header blocks */

        /* calculate total width of the table */
        naxis1 = 0;
        for ii in 0..(tfields as usize) {
            let tform_item = v_tform[ii];
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

        datasize = ((naxis1 as LONGLONG) * naxis2) + pcount; /* size of table in bytes */
        nblocks = (((datasize + (BL!() - 1)) / BL!()) + nhead as LONGLONG) as c_long; /* size of HDU */

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

        nexthdu = (fptr.Fptr.curhdu) + 1; /* number of the next (new) hdu */
        newstart = headstart[nexthdu as usize]; /* save starting addr of HDU */

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
        fptr.Fptr.datastart = (hs_item) + (nhead as LONGLONG * BL!());
        fptr.Fptr.hdutype = BINARY_TBL; /* might need to be reset... */

        /* write the required header keywords. This will write PCOUNT = 0 */
        /* so that the variable length data will be written at the right place */
        ffphbn_safe(
            fptr,
            naxis2,
            tfields,
            &v_ttype,
            &v_tform,
            v_tunit.as_deref(),
            Some(&extnm),
            pcount,
            status,
        );

        /* redefine internal structure for this HDU (with PCOUNT = 0) */
        ffrdef_safe(fptr, status);

        *status
    }
}

/*--------------------------------------------------------------------------*/
/// Delete the CHDU.  If the CHDU is the primary array, then replace the HDU
/// with an empty primary array with no data.   Return the
/// type of the new CHDU after the old CHDU is deleted.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffdhdu(
    fptr: *mut fitsfile, /* I - FITS file pointer                   */
    hdutype: *mut c_int, /* O - type of the new CHDU after deletion */
    status: *mut c_int,  /* IO - error status                       */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let hdutype = hdutype.as_mut();

        ffdhdu_safer(fptr, hdutype, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Delete the CHDU.  If the CHDU is the primary array, then replace the HDU
/// with an empty primary array with no data.   Return the
/// type of the new CHDU after the old CHDU is deleted.
pub(crate) unsafe fn ffdhdu_safer(
    fptr: &mut fitsfile,         /* I - FITS file pointer                   */
    hdutype: Option<&mut c_int>, /* O - type of the new CHDU after deletion */
    status: &mut c_int,          /* IO - error status                       */
) -> c_int {
    unsafe {
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

            if ffrhdu_safer(fptr, Some(&mut tmptype), status) > 0 {
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
}
