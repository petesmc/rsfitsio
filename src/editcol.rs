/*  This file, editcol.c, contains the set of FITSIO routines that    */
/*  insert or delete rows or columns in a table or resize an image    */

/*  The FITSIO software was written by William Pence at the High Energy    */
/*  Astrophysic Science Archive Research Center (HEASARC) at the NASA      */
/*  Goddard Space Flight Center.                                           */

use core::slice;
use std::ffi::CStr;
use std::{cmp, mem};

use crate::c_types::{c_char, c_int, c_long, c_longlong, c_ulonglong};

use bytemuck::{cast_slice, cast_slice_mut};

use crate::fitscore::{
    ffasfm_safe, ffbnfm_safe, ffc2ii, ffcmph_safer, ffdblk, ffeqty_safe, ffgdesll_safe,
    ffgtcl_safe, ffiblk, ffkeyn_safe, ffmahd_safe, ffmkky_safe, ffnkey_safe, ffpdes_safe,
    ffrdef_safe, ffupch_safe,
};
use crate::fitscore::{ffpmsg_slice, ffpmsg_str};
use crate::getcold::{ffgcvd_safe, ffgcvm_safe};
use crate::getcole::ffgcvc_safe;
use crate::getcolj::ffgcvjj_safe;
use crate::getcoll::ffgcvl_safe;
use crate::getcols::ffgcvs_safe;
use crate::getcoluj::ffgcvujj_safe;
use crate::getkey::{
    ffghprll_safe, ffghsp_safe, ffgkey_safe, ffgkyj_safe, ffgkyjj_safe, ffgkys_safe, ffgrec_safe,
};
use crate::modkey::{
    ffdkey_safe, ffdrec_safe, ffikyj_safe, ffmcom_safe, ffmkyj_safe, ffmkys_safe, ffmrec_safe,
    ffuky_safe, ffukyg_safe,
};
use crate::{KeywordDatatype, fitsio::*, vecs_to_slices, vecs_to_slices_mut};
use crate::{fitsio2::*, int_snprintf, slice_to_str};

use crate::BL;
use crate::putcold::{ffpcld_safe, ffpclm_safe, ffpcnd_safe};
use crate::putcole::ffpclc_safe;
use crate::putcolj::ffpcljj_safe;
use crate::putcoll::ffpcnl_safe;
use crate::putcols::{ffpcls_safe, ffpcns_safe};
use crate::putcoluj::ffpclujj_safe;
use crate::putkey::{ffpkyg_safe, ffpkyj_safe, ffpkys_safe, ffprec_safe};
use crate::wrappers::*;
use crate::{bb, cs};
use crate::{buffers::*, raw_to_slice};

/*--------------------------------------------------------------------------*/
/// resize an existing primary array or IMAGE extension.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffrsim(
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

        ffrsim_safe(fptr, bitpix, naxis, naxes, status)
    }
}

/*--------------------------------------------------------------------------*/
/// resize an existing primary array or IMAGE extension.
pub(crate) fn ffrsim_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer           */
    bitpix: c_int,       /* I - bits per pixel              */
    naxis: c_int,        /* I - number of axes in the array */
    naxes: &[c_long],    /* I - size of each axis           */
    status: &mut c_int,  /* IO - error status               */
) -> c_int {
    let mut tnaxes: [LONGLONG; 99] = [0; 99];

    if *status > 0 {
        return *status;
    }

    let _to = cmp::min(naxis, 99) as usize;

    for ii in 0.._to {
        tnaxes[ii] = naxes[ii] as LONGLONG;
    }

    ffrsimll_safe(fptr, bitpix, naxis, &tnaxes, status);

    *status
}

/*--------------------------------------------------------------------------*/
/// resize an existing primary array or IMAGE extension.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffrsimll(
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

        ffrsimll_safe(fptr, bitpix, naxis, naxes, status)
    }
}

/*--------------------------------------------------------------------------*/
/// resize an existing primary array or IMAGE extension.
pub(crate) fn ffrsimll_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer           */
    bitpix: c_int,       /* I - bits per pixel              */
    naxis: c_int,        /* I - number of axes in the array */
    naxes: &[LONGLONG],  /* I - size of each axis           */
    status: &mut c_int,  /* IO - error status               */
) -> c_int {
    let mut simple: c_int = 0;
    let mut obitpix: c_int = 0;
    let mut onaxis: c_int = 0;
    let mut extend: c_int = 0;
    let mut nmodify: c_int = 0;
    let mut nblocks: c_long = 0;
    let mut longval: c_long = 0;
    let mut pcount: c_long = 0;
    let mut gcount: c_long = 0;
    let mut longbitpix: c_int = 0; // WARNING: The name of the variable indicates it should be a long but seems superfluous
    let mut onaxes: [LONGLONG; 99] = [0; 99];
    let mut newsize: LONGLONG = 0;
    let mut oldsize: LONGLONG = 0;
    let mut comment: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];
    let mut keyname: [c_char; FLEN_KEYWORD] = [0; FLEN_KEYWORD];
    let mut message: [c_char; FLEN_ERRMSG] = [0; FLEN_ERRMSG];

    if *status > 0 {
        return *status;
    }

    /* reset position to the correct HDU if necessary */
    if fptr.HDUposition != fptr.Fptr.curhdu {
        ffmahd_safe(fptr, (fptr.HDUposition) + 1, None, status);
    }
    /* rescan header if data structure is undefined */
    else if fptr.Fptr.datastart == DATA_UNDEFINED as LONGLONG && ffrdef_safe(fptr, status) > 0 {
        return *status;
    }

    /* get current image size parameters */
    if ffghprll_safe(
        fptr,
        99,
        &mut simple,
        &mut obitpix,
        Some(&mut onaxis),
        &mut onaxes,
        &mut pcount,
        &mut gcount,
        &mut extend,
        status,
    ) > 0
    {
        return *status;
    }

    longbitpix = bitpix;

    /* test for the 4 special cases that represent unsigned integers
    or signed bytes */
    if longbitpix == USHORT_IMG {
        longbitpix = SHORT_IMG;
    } else if longbitpix == ULONG_IMG {
        longbitpix = LONG_IMG;
    } else if longbitpix == SBYTE_IMG {
        longbitpix = BYTE_IMG;
    } else if longbitpix == ULONGLONG_IMG {
        longbitpix = LONGLONG_IMG;
    }

    /* test that the new values are legal */

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

    if naxis == 0 {
        newsize = 0;
    } else {
        newsize = 1;
    }

    for ii in 0..(naxis as usize) {
        if naxes[ii] < 0 {
            int_snprintf!(
                &mut message,
                FLEN_ERRMSG,
                "Illegal value for NAXIS{} keyword: {:.0}",
                ii + 1,
                (naxes[ii]) as f64,
            );
            ffpmsg_slice(&message);
            *status = BAD_NAXES;
            return *status;
        }

        newsize *= naxes[ii]; /* compute new image size, in pixels */
    }

    /* compute size of old image, in bytes */

    if onaxis == 0 {
        oldsize = 0;
    } else {
        oldsize = 1;
        for ii in 0..(onaxis as usize) {
            oldsize *= onaxes[ii];
        }
        oldsize =
            (oldsize + pcount as LONGLONG) * (gcount as LONGLONG) * (obitpix.abs() / 8) as LONGLONG;
    }

    oldsize = (oldsize + BL!() - 1) / BL!(); /* old size, in blocks */

    newsize =
        (newsize + pcount as LONGLONG) * (gcount as LONGLONG) * (longbitpix.abs() / 8) as LONGLONG;
    newsize = (newsize + BL!() - 1) / BL!(); /* new size, in blocks */

    if newsize > oldsize {
        /* have to insert new blocks for image */

        nblocks = (newsize - oldsize) as c_long;
        if ffiblk(fptr, nblocks, 1, status) > 0 {
            return *status;
        }
    } else if oldsize > newsize {
        /* have to delete blocks from image */

        nblocks = (oldsize - newsize) as c_long;
        if ffdblk(fptr, nblocks, status) > 0 {
            return *status;
        }
    }

    /* now update the header keywords */

    strcpy_safe(&mut comment, cs!(c"&")); /* special value to leave comments unchanged */

    if longbitpix != obitpix {
        /* update BITPIX value */
        ffmkyj_safe(
            fptr,
            cs!(c"BITPIX"),
            longbitpix as LONGLONG,
            Some(&comment),
            status,
        );
    }

    if naxis != onaxis {
        /* update NAXIS value */
        longval = naxis as c_long;
        ffmkyj_safe(
            fptr,
            cs!(c"NAXIS"),
            longval as LONGLONG,
            Some(&comment),
            status,
        );
    }

    /* modify the existing NAXISn keywords */
    nmodify = cmp::min(naxis, onaxis);
    for ii in 0..(nmodify as usize) {
        ffkeyn_safe(cs!(c"NAXIS"), (ii + 1) as c_int, &mut keyname, status);
        ffmkyj_safe(fptr, &keyname, naxes[ii], Some(&comment), status);
    }

    if naxis > onaxis
    /* insert additional NAXISn keywords */
    {
        strcpy_safe(&mut comment, cs!(c"length of data axis"));
        for ii in (onaxis as usize)..(naxis as usize) {
            ffkeyn_safe(cs!(c"NAXIS"), (ii + 1) as c_int, &mut keyname, status);
            ffikyj_safe(fptr, &keyname, naxes[ii], Some(&comment), status);
        }
    } else if onaxis > naxis
    /* delete old NAXISn keywords */
    {
        for ii in (naxis as usize)..(onaxis as usize) {
            ffkeyn_safe(cs!(c"NAXIS"), (ii + 1) as c_int, &mut keyname, status);
            ffdkey_safe(fptr, &keyname, status);
        }
    }

    /* Update the BSCALE and BZERO keywords, if an unsigned integer image
    or a signed byte image.  */
    if bitpix == USHORT_IMG {
        strcpy_safe(
            &mut comment,
            cs!(c"offset data range to that of unsigned short"),
        );
        ffukyg_safe(fptr, cs!(c"BZERO"), 32768.0, 0, Some(&comment), status);
        strcpy_safe(&mut comment, cs!(c"default scaling factor"));
        ffukyg_safe(fptr, cs!(c"BSCALE"), 1.0, 0, Some(&comment), status);
    } else if bitpix == ULONG_IMG {
        strcpy_safe(
            &mut comment,
            cs!(c"offset data range to that of unsigned long"),
        );
        ffukyg_safe(fptr, cs!(c"BZERO"), 2147483648.0, 0, Some(&comment), status);
        strcpy_safe(&mut comment, cs!(c"default scaling factor"));
        ffukyg_safe(fptr, cs!(c"BSCALE"), 1.0, 0, Some(&comment), status);
    } else if bitpix == ULONGLONG_IMG {
        strcpy_safe(
            &mut comment,
            cs!(c"offset data range to that of unsigned long long"),
        );
        ffukyg_safe(
            fptr,
            cs!(c"BZERO"),
            9223372036854775808.,
            0,
            Some(&comment),
            status,
        );
        strcpy_safe(&mut comment, cs!(c"default scaling factor"));
        ffukyg_safe(fptr, cs!(c"BSCALE"), 1.0, 0, Some(&comment), status);
    } else if bitpix == SBYTE_IMG {
        strcpy_safe(
            &mut comment,
            cs!(c"offset data range to that of signed byte"),
        );
        ffukyg_safe(fptr, cs!(c"BZERO"), -128., 0, Some(&comment), status);
        strcpy_safe(&mut comment, cs!(c"default scaling factor"));
        ffukyg_safe(fptr, cs!(c"BSCALE"), 1.0, 0, Some(&comment), status);
    }

    /* re-read the header, to make sure structures are updated */
    ffrdef_safe(fptr, status);
    *status
}

/*--------------------------------------------------------------------------*/
/// insert NROWS blank rows immediated after row firstrow (1 = first row).
/// Set firstrow = 0 to insert space at the beginning of the table.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffirow(
    fptr: *mut fitsfile, /* I - FITS file pointer                        */
    firstrow: LONGLONG, /* I - insert space AFTER this row, 0 = insert space at beginning of table              */
    nrows: LONGLONG,    /* I - number of rows to insert                 */
    status: *mut c_int, /* IO - error status                            */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        ffirow_safe(fptr, firstrow, nrows, status)
    }
}

/*--------------------------------------------------------------------------*/
/// insert NROWS blank rows immediated after row firstrow (1 = first row).
/// Set firstrow = 0 to insert space at the beginning of the table.
pub fn ffirow_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer                        */
    firstrow: LONGLONG, /* I - insert space AFTER this row, 0 = insert space at beginning of table              */
    nrows: LONGLONG,    /* I - number of rows to insert                 */
    status: &mut c_int, /* IO - error status                            */
) -> c_int {
    let mut naxis1: LONGLONG = 0;
    let mut naxis2: LONGLONG = 0;
    let mut datasize: LONGLONG = 0;
    let mut firstbyte: LONGLONG = 0;
    let mut nshift: LONGLONG = 0;
    let mut nbytes: LONGLONG = 0;
    let mut freespace: LONGLONG = 0;
    let mut nblock: c_long = 0;

    if *status > 0 {
        return *status;
    }

    /* reset position to the correct HDU if necessary */
    if fptr.HDUposition != fptr.Fptr.curhdu {
        ffmahd_safe(fptr, (fptr.HDUposition) + 1, None, status);
    }
    /* rescan header if data structure is undefined */
    else if fptr.Fptr.datastart == DATA_UNDEFINED as LONGLONG && ffrdef_safe(fptr, status) > 0 {
        return *status;
    }

    if fptr.Fptr.hdutype == IMAGE_HDU {
        ffpmsg_str("Can only add rows to TABLE or BINTABLE extension (ffirow)");
        *status = NOT_TABLE;
        return *status;
    }

    if nrows < 0 {
        *status = NEG_BYTES;
        return *status;
    } else if nrows == 0 {
        return *status; /* no op, so just return */
    }

    /* get the current size of the table */
    /* use internal structure since NAXIS2 keyword may not be up to date */
    naxis1 = fptr.Fptr.rowlength;
    naxis2 = fptr.Fptr.numrows;

    if firstrow > naxis2 {
        ffpmsg_str("Insert position greater than the number of rows in the table (ffirow)");
        *status = BAD_ROW_NUM;
        return *status;
    } else if firstrow < 0 {
        ffpmsg_str("Insert position is less than 0 (ffirow)");
        *status = BAD_ROW_NUM;
        return *status;
    }

    /* current data size */
    datasize = fptr.Fptr.heapstart + fptr.Fptr.heapsize;
    freespace = (((datasize + 2879) / 2880) * 2880) - datasize;
    nshift = naxis1 * nrows; /* no. of bytes to add to table */

    if (freespace - nshift) < 0
    /* not enough existing space? */
    {
        nblock = ((nshift - freespace + 2879) / 2880) as c_long; /* number of blocks */
        ffiblk(fptr, nblock, 1, status); /* insert the blocks */
    }

    firstbyte = naxis1 * firstrow; /* relative insert position */
    nbytes = datasize - firstbyte; /* no. of bytes to shift down */
    firstbyte += fptr.Fptr.datastart; /* absolute insert position */

    if nshift > 0 {
        /* nshift may be zero if naxis1 == naxis2 == 0 */
        ffshft_safe(fptr, firstbyte, nbytes, nshift, status); /* shift rows and heap */
    }

    /* update the heap starting address */
    fptr.Fptr.heapstart += nshift;

    /* update the THEAP keyword if it exists */
    let mut tstatus = 0;
    ffmkyj_safe(
        fptr,
        cs!(c"THEAP"),
        fptr.Fptr.heapstart,
        Some(cs!(c"&")),
        &mut tstatus,
    );

    /* update the NAXIS2 keyword */
    ffmkyj_safe(
        fptr,
        cs!(c"NAXIS2"),
        naxis2 + nrows,
        Some(cs!(c"&")),
        status,
    );
    (fptr.Fptr.numrows) += nrows;
    (fptr.Fptr.origrows) += nrows;

    *status
}

/*--------------------------------------------------------------------------*/
/// delete NROWS rows from table starting with firstrow (1 = first row of table).
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffdrow(
    fptr: *mut fitsfile, /* I - FITS file pointer                        */
    firstrow: LONGLONG,  /* I - first row to delete (1 = first)          */
    nrows: LONGLONG,     /* I - number of rows to delete                 */
    status: *mut c_int,  /* IO - error status                            */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        ffdrow_safe(fptr, firstrow, nrows, status)
    }
}

/*--------------------------------------------------------------------------*/
/// delete NROWS rows from table starting with firstrow (1 = first row of table).
pub(crate) fn ffdrow_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer                        */
    firstrow: LONGLONG,  /* I - first row to delete (1 = first)          */
    nrows: LONGLONG,     /* I - number of rows to delete                 */
    status: &mut c_int,  /* IO - error status                            */
) -> c_int {
    let mut naxis1: LONGLONG = 0;
    let mut naxis2: LONGLONG = 0;
    let mut datasize: LONGLONG = 0;
    let mut firstbyte: LONGLONG = 0;
    let mut nshift: LONGLONG = 0;
    let mut nbytes: LONGLONG = 0;
    let mut freespace: LONGLONG = 0;
    let mut nblock: c_long = 0;
    let mut comm: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];

    if *status > 0 {
        return *status;
    }

    if fptr.HDUposition != fptr.Fptr.curhdu {
        ffmahd_safe(fptr, (fptr.HDUposition) + 1, None, status);
    }
    /* rescan header if data structure is undefined */
    else if fptr.Fptr.datastart == DATA_UNDEFINED as LONGLONG && ffrdef_safe(fptr, status) > 0 {
        return *status;
    }

    if fptr.Fptr.hdutype == IMAGE_HDU {
        ffpmsg_str("Can only delete rows in TABLE or BINTABLE extension (ffdrow)");
        *status = NOT_TABLE;
        return *status;
    }

    if nrows < 0 {
        *status = NEG_BYTES;
        return *status;
    } else if nrows == 0 {
        return *status; /* no op, so just return */
    }

    ffgkyjj_safe(fptr, cs!(c"NAXIS1"), &mut naxis1, Some(&mut comm), status); /* get the current   */

    /* ffgkyj(fptr, "NAXIS2", &naxis2, comm, status);*/
    /* size of the table */

    /* the NAXIS2 keyword may not be up to date, so use the structure value */
    naxis2 = fptr.Fptr.numrows;

    if firstrow > naxis2 {
        ffpmsg_str("Delete position greater than the number of rows in the table (ffdrow)");
        *status = BAD_ROW_NUM;
        return *status;
    } else if firstrow < 1 {
        ffpmsg_str("Delete position is less than 1 (ffdrow)");
        *status = BAD_ROW_NUM;
        return *status;
    } else if firstrow + nrows - 1 > naxis2 {
        ffpmsg_str("No. of rows to delete exceeds size of table (ffdrow)");
        *status = BAD_ROW_NUM;
        return *status;
    }

    nshift = naxis1 * nrows; /* no. of bytes to delete from table */
    /* cur size of data */
    datasize = fptr.Fptr.heapstart + fptr.Fptr.heapsize;

    firstbyte = naxis1 * (firstrow + nrows - 1); /* relative del pos */
    nbytes = datasize - firstbyte; /* no. of bytes to shift up */
    firstbyte += fptr.Fptr.datastart; /* absolute delete position */

    ffshft_safe(fptr, firstbyte, nbytes, -nshift, status); /* shift data */

    freespace = (((datasize + 2879) / 2880) * 2880) - datasize;
    nblock = ((nshift + freespace) / 2880) as c_long; /* number of blocks */

    /* delete integral number blocks */
    if nblock > 0 {
        ffdblk(fptr, nblock, status);
    }

    /* update the heap starting address */
    fptr.Fptr.heapstart -= nshift;

    /* update the THEAP keyword if it exists */
    let mut tstatus = 0;
    ffmkyj_safe(
        fptr,
        cs!(c"THEAP"),
        fptr.Fptr.heapstart as LONGLONG,
        Some(cs!(c"&")),
        &mut tstatus,
    );

    /* update the NAXIS2 keyword */
    ffmkyj_safe(
        fptr,
        cs!(c"NAXIS2"),
        naxis2 - nrows,
        Some(cs!(c"&")),
        status,
    );
    (fptr.Fptr.numrows) -= nrows;
    (fptr.Fptr.origrows) -= nrows;

    /* Update the heap data, if any.  This will remove any orphaned data */
    /* that was only pointed to by the rows that have been deleted */
    unsafe {
        ffcmph_safer(fptr, status);
    }
    *status
}

/*--------------------------------------------------------------------------*/
/// Delete the ranges of rows from the table (1 = first row of table).
///
/// The 'ranges' parameter typically looks like:
/// '10-20, 30 - 40, 55' or '50-'
/// and gives a list of rows or row ranges separated by commas.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffdrrg(
    fptr: *mut fitsfile,   /* I - FITS file pointer to table               */
    ranges: *const c_char, /* I - ranges of rows to delete (1 = first)     */
    status: *mut c_int,    /* IO - error status                            */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        raw_to_slice!(ranges);

        ffdrrg_safe(fptr, ranges, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Delete the ranges of rows from the table (1 = first row of table).
///
/// The 'ranges' parameter typically looks like:
/// '10-20, 30 - 40, 55' or '50-'
/// and gives a list of rows or row ranges separated by commas.
pub(crate) fn ffdrrg_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer to table               */
    ranges: &[c_char],   /* I - ranges of rows to delete (1 = first)     */
    status: &mut c_int,  /* IO - error status                            */
) -> c_int {
    let nranges: c_int = 0;
    let mut nranges2: c_int = 0;
    let ii: c_int = 0;
    let mut minrow: *mut c_long;
    let mut maxrow: *mut c_long;
    let mut nrows: c_long = 0;
    let mut rowarray: *mut c_long;
    let jj: c_long = 0;
    let kk: c_long = 0;
    let mut naxis2: LONGLONG = 0;

    if *status > 0 {
        return *status;
    }

    if fptr.HDUposition != fptr.Fptr.curhdu {
        ffmahd_safe(fptr, (fptr.HDUposition) + 1, None, status);
    }
    /* rescan header if data structure is undefined */
    else if fptr.Fptr.datastart == DATA_UNDEFINED as LONGLONG && ffrdef_safe(fptr, status) > 0 {
        return *status;
    }

    if fptr.Fptr.hdutype == IMAGE_HDU {
        ffpmsg_str("Can only delete rows in TABLE or BINTABLE extension (ffdrrg)");
        *status = NOT_TABLE;
        return *status;
    }

    /* the NAXIS2 keyword may not be up to date, so use the structure value */
    naxis2 = fptr.Fptr.numrows;

    /* find how many ranges were specified ( = no. of commas in string + 1) */
    let mut cptr = 0;

    let mut nranges = 1;

    let mut c_found = strchr_safe(&ranges[cptr..], bb(b','));
    while c_found.is_some() {
        //for (nranges = 1; (cptr = strchr_safe(cptr, ',')); nranges++){
        cptr += 1;
        nranges += 1;
        c_found = strchr_safe(&ranges[cptr..], bb(b','));
    }

    let mut minrow: Vec<c_long> = Vec::new();
    let mut maxrow: Vec<c_long> = Vec::new();

    if minrow.try_reserve_exact(nranges).is_err() || minrow.try_reserve_exact(nranges).is_err() {
        *status = MEMORY_ALLOCATION;
        ffpmsg_str("failed to allocate memory for row ranges (ffdrrg)");
        return *status;
    } else {
        minrow.resize(nranges, 0);
        maxrow.resize(nranges, 0);
    }

    /* parse range list into array of range min and max values */
    ffrwrg_safe(
        ranges,
        naxis2,
        nranges as c_int,
        &mut nranges2,
        &mut minrow,
        &mut maxrow,
        status,
    );
    if *status > 0 || nranges2 == 0 {
        return *status;
    }

    /* determine total number or rows to delete */
    nrows = 0;
    for ii in 0..(nranges2 as usize) {
        nrows = nrows + maxrow[ii] - minrow[ii] + 1;
    }

    let mut rowarray: Vec<c_long> = Vec::new();
    if rowarray.try_reserve_exact(nrows as usize).is_err() {
        *status = MEMORY_ALLOCATION;
        ffpmsg_str("failed to allocate memory for row array (ffdrrg)");
        return *status;
    } else {
        rowarray.resize(nrows as usize, 0);
    }

    let mut kk = 0;
    for ii in 0..(nranges2 as usize) {
        for jj in (minrow[ii])..=(maxrow[ii]) {
            rowarray[kk] = jj;
            kk += 1;
        }
    }

    /* delete the rows */
    ffdrws_safe(fptr, &rowarray, nrows, status);

    *status
}

/*--------------------------------------------------------------------------*/
/// delete the list of rows from the table (1 = first row of table).
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffdrws(
    fptr: *mut fitsfile,   /* I - FITS file pointer                        */
    rownum: *const c_long, /* I - list of rows to delete (1 = first)       */
    nrows: c_long,         /* I - number of rows to delete                 */
    status: *mut c_int,    /* IO - error status                            */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        let rownum = slice::from_raw_parts(rownum, nrows as usize);

        ffdrws_safe(fptr, rownum, nrows, status)
    }
}

/*--------------------------------------------------------------------------*/
/// delete the list of rows from the table (1 = first row of table).
pub(crate) fn ffdrws_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer                        */
    rownum: &[c_long],   /* I - list of rows to delete (1 = first)       */
    nrows: c_long,       /* I - number of rows to delete                 */
    status: &mut c_int,  /* IO - error status                            */
) -> c_int {
    let mut naxis1: LONGLONG = 0;
    let mut naxis2: LONGLONG = 0;
    let mut insertpos: LONGLONG = 0;
    let mut nextrowpos: LONGLONG = 0;
    let ii: c_long = 0;
    let mut nextrow: LONGLONG = 0;
    let mut comm: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];
    let mut buffer: *mut u8;

    if *status > 0 {
        return *status;
    }

    if fptr.HDUposition != fptr.Fptr.curhdu {
        ffmahd_safe(fptr, (fptr.HDUposition) + 1, None, status);
    }

    /* rescan header if data structure is undefined */
    if fptr.Fptr.datastart == DATA_UNDEFINED as LONGLONG && ffrdef_safe(fptr, status) > 0 {
        return *status;
    }

    if fptr.Fptr.hdutype == IMAGE_HDU {
        ffpmsg_str("Can only delete rows in TABLE or BINTABLE extension (ffdrws)");
        *status = NOT_TABLE;
        return *status;
    }

    if nrows < 0 {
        *status = NEG_BYTES;
        return *status;
    } else if nrows == 0 {
        return *status; /* no op, so just return */
    }

    ffgkyjj_safe(fptr, cs!(c"NAXIS1"), &mut naxis1, Some(&mut comm), status); /* row width   */
    ffgkyjj_safe(fptr, cs!(c"NAXIS2"), &mut naxis2, Some(&mut comm), status); /* number of rows */

    /* check that input row list is in ascending order */
    for ii in 1..(nrows as usize) {
        if rownum[ii - 1] >= rownum[ii] {
            ffpmsg_str("row numbers are not in increasing order (ffdrws)");
            *status = BAD_ROW_NUM;
            return *status;
        }
    }

    if rownum[0] < 1 {
        ffpmsg_str("first row to delete is less than 1 (ffdrws)");
        *status = BAD_ROW_NUM;
        return *status;
    } else if rownum[(nrows - 1) as usize] as LONGLONG > naxis2 {
        ffpmsg_str("last row to delete exceeds size of table (ffdrws)");
        *status = BAD_ROW_NUM;
        return *status;
    }

    let mut buffer: Vec<u8> = Vec::new();

    /* buffer for one row */
    if buffer.try_reserve_exact(naxis1 as usize).is_err() {
        ffpmsg_str("malloc failed (ffdrws)");
        *status = MEMORY_ALLOCATION;
        return *status;
    } else {
        buffer.resize(naxis1 as usize, 0);
    }

    /* byte location to start of first row to delete, and the next row */
    insertpos = fptr.Fptr.datastart + ((rownum[0] as LONGLONG - 1) * naxis1);
    nextrowpos = insertpos + naxis1;
    nextrow = rownum[0] as LONGLONG + 1;

    /* work through the list of rows to delete */
    let mut ii: usize = 1;
    while ii < nrows as usize {
        if nextrow < rownum[ii] as LONGLONG {
            /* keep this row, so copy it to the new position */

            ffmbyt_safe(fptr, nextrowpos, REPORT_EOF, status);
            ffgbyt(fptr, naxis1, cast_slice_mut(&mut buffer), status); /* read the bytes */

            ffmbyt_safe(fptr, insertpos, IGNORE_EOF, status);
            ffpbyt(fptr, naxis1, cast_slice_mut(&mut buffer), status); /* write the bytes */

            if *status > 0 {
                ffpmsg_str("error while copying good rows in table (ffdrws)");
                return *status;
            }
            insertpos += naxis1;
        } else {
            /* skip over this row since it is in the list */
            ii += 1;
        }

        nextrow += 1;
        nextrowpos += naxis1;
    }

    /* finished with all the rows to delete; copy remaining rows */
    while nextrow <= naxis2 {
        ffmbyt_safe(fptr, nextrowpos, REPORT_EOF, status);
        ffgbyt(fptr, naxis1, cast_slice_mut(&mut buffer), status); /* read the bytes */

        ffmbyt_safe(fptr, insertpos, IGNORE_EOF, status);
        ffpbyt(fptr, naxis1, cast_slice_mut(&mut buffer), status); /* write the bytes */

        if *status > 0 {
            ffpmsg_str("failed to copy remaining rows in table (ffdrws)");
            return *status;
        }
        insertpos += naxis1;
        nextrowpos += naxis1;
        nextrow += 1;
    }

    /* now delete the empty rows at the end of the table */
    ffdrow_safe(
        fptr,
        naxis2 - (nrows as LONGLONG) + 1,
        nrows as LONGLONG,
        status,
    );

    /* Update the heap data, if any.  This will remove any orphaned data */
    /* that was only pointed to by the rows that have been deleted */
    unsafe {
        ffcmph_safer(fptr, status);
    }

    *status
}

/*--------------------------------------------------------------------------*/
/// delete the list of rows from the table (1 = first row of table).
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffdrwsll(
    fptr: *mut fitsfile,     /* I - FITS file pointer                        */
    rownum: *const LONGLONG, /* I - list of rows to delete (1 = first)       */
    nrows: LONGLONG,         /* I - number of rows to delete                 */
    status: *mut c_int,      /* IO - error status                            */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        let rownum = slice::from_raw_parts(rownum, nrows as usize);

        ffdrwsll_safe(fptr, rownum, nrows, status)
    }
}

/*--------------------------------------------------------------------------*/
/// delete the list of rows from the table (1 = first row of table).
pub(crate) fn ffdrwsll_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer                        */
    rownum: &[LONGLONG], /* I - list of rows to delete (1 = first)       */
    nrows: LONGLONG,     /* I - number of rows to delete                 */
    status: &mut c_int,  /* IO - error status                            */
) -> c_int {
    let mut insertpos: LONGLONG = 0;
    let mut nextrowpos: LONGLONG = 0;
    let mut naxis1: LONGLONG = 0;
    let mut naxis2: LONGLONG = 0;
    let ii: LONGLONG = 0;
    let mut nextrow: LONGLONG = 0;
    let mut comm: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];
    let mut buffer: *mut u8;

    if *status > 0 {
        return *status;
    }

    if fptr.HDUposition != fptr.Fptr.curhdu {
        ffmahd_safe(fptr, (fptr.HDUposition) + 1, None, status);
    }

    /* rescan header if data structure is undefined */
    if fptr.Fptr.datastart == DATA_UNDEFINED as LONGLONG && ffrdef_safe(fptr, status) > 0 {
        return *status;
    }

    if fptr.Fptr.hdutype == IMAGE_HDU {
        ffpmsg_str("Can only delete rows in TABLE or BINTABLE extension (ffdrws)");
        *status = NOT_TABLE;
        return *status;
    }

    if nrows < 0 {
        *status = NEG_BYTES;
        return *status;
    } else if nrows == 0 {
        return *status; /* no op, so just return */
    }

    ffgkyjj_safe(fptr, cs!(c"NAXIS1"), &mut naxis1, Some(&mut comm), status); /* row width   */
    ffgkyjj_safe(fptr, cs!(c"NAXIS2"), &mut naxis2, Some(&mut comm), status); /* number of rows */

    /* check that input row list is in ascending order */
    for ii in 1..(nrows as usize) {
        if rownum[ii - 1] >= rownum[ii] {
            ffpmsg_str("row numbers are not in increasing order (ffdrws)");
            *status = BAD_ROW_NUM;
            return *status;
        }
    }

    if rownum[0] < 1 {
        ffpmsg_str("first row to delete is less than 1 (ffdrws)");
        *status = BAD_ROW_NUM;
        return *status;
    } else if rownum[(nrows - 1) as usize] > naxis2 {
        ffpmsg_str("last row to delete exceeds size of table (ffdrws)");
        *status = BAD_ROW_NUM;
        return *status;
    }

    let mut buffer: Vec<u8> = Vec::new();

    /* buffer for one row */
    if buffer.try_reserve_exact(naxis1 as usize).is_err() {
        ffpmsg_str("malloc failed (ffdrwsll)");
        *status = MEMORY_ALLOCATION;
        return *status;
    } else {
        buffer.resize(naxis1 as usize, 0);
    }

    /* byte location to start of first row to delete, and the next row */
    insertpos = fptr.Fptr.datastart + ((rownum[0] - 1) * naxis1);
    nextrowpos = insertpos + naxis1;
    nextrow = rownum[0] + 1;

    /* work through the list of rows to delete */
    let mut ii: usize = 1;
    while ii < nrows as usize {
        if nextrow < rownum[ii] {
            /* keep this row, so copy it to the new position */

            ffmbyt_safe(fptr, nextrowpos, REPORT_EOF, status);
            ffgbyt(fptr, naxis1, cast_slice_mut(&mut buffer), status); /* read the bytes */

            ffmbyt_safe(fptr, insertpos, IGNORE_EOF, status);
            ffpbyt(fptr, naxis1, cast_slice_mut(&mut buffer), status); /* write the bytes */

            if *status > 0 {
                ffpmsg_str("error while copying good rows in table (ffdrws)");
                return *status;
            }
            insertpos += naxis1;
        } else {
            /* skip over this row since it is in the list */
            ii += 1;
        }
        nextrow += 1;
        nextrowpos += naxis1;
    }

    /* finished with all the rows to delete; copy remaining rows */
    while nextrow <= naxis2 {
        ffmbyt_safe(fptr, nextrowpos, REPORT_EOF, status);
        ffgbyt(fptr, naxis1, cast_slice_mut(&mut buffer), status); /* read the bytes */

        ffmbyt_safe(fptr, insertpos, IGNORE_EOF, status);
        ffpbyt(fptr, naxis1, cast_slice_mut(&mut buffer), status); /* write the bytes */

        if *status > 0 {
            ffpmsg_str("failed to copy remaining rows in table (ffdrws)");
            return *status;
        }
        insertpos += naxis1;
        nextrowpos += naxis1;
        nextrow += 1;
    }

    /* now delete the empty rows at the end of the table */
    ffdrow_safe(fptr, naxis2 - nrows + 1, nrows, status);

    /* Update the heap data, if any.  This will remove any orphaned data */
    /* that was only pointed to by the rows that have been deleted */
    unsafe {
        ffcmph_safer(fptr, status);
    }

    *status
}

/*--------------------------------------------------------------------------*/
/// parse the input list of row ranges, returning the number of ranges,
/// and the min and max row value in each range.
///
/// The only characters allowed in the input rowlist are
/// decimal digits, minus sign, and comma (and non-significant spaces)
///
/// Example:  
///
/// list = "10-20, 30-35,50"
///
/// would return numranges = 3, minrow[] = {10, 30, 50}, maxrow[] = {20, 35, 50}
///
/// error is returned if min value of range is > max value of range or if the
/// ranges are not monotonically increasing.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffrwrg(
    rowlist: *const c_char, /* I - list of rows and row ranges */
    maxrows: LONGLONG,      /* I - number of rows in the table */
    maxranges: c_int,       /* I - max number of ranges to be returned */
    numranges: *mut c_int,  /* O - number ranges returned */
    minrow: *mut c_long,    /* O - first row in each range */
    maxrow: *mut c_long,    /* O - last row in each range */
    status: *mut c_int,     /* IO - status value */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let numranges = numranges.as_mut().expect(NULL_MSG);
        let minrow = slice::from_raw_parts_mut(minrow, maxranges as usize);
        let maxrow = slice::from_raw_parts_mut(maxrow, maxranges as usize);

        raw_to_slice!(rowlist);

        ffrwrg_safe(
            rowlist, maxrows, maxranges, numranges, minrow, maxrow, status,
        )
    }
}

/*--------------------------------------------------------------------------*/
/// parse the input list of row ranges, returning the number of ranges,
/// and the min and max row value in each range.
///
/// The only characters allowed in the input rowlist are
/// decimal digits, minus sign, and comma (and non-significant spaces)
///
/// Example:  
///
/// list = "10-20, 30-35,50"
///
/// would return numranges = 3, minrow[] = {10, 30, 50}, maxrow[] = {20, 35, 50}
///
/// error is returned if min value of range is > max value of range or if the
/// ranges are not monotonically increasing.
pub(crate) fn ffrwrg_safe(
    rowlist: &[c_char],    /* I - list of rows and row ranges */
    maxrows: LONGLONG,     /* I - number of rows in the table */
    maxranges: c_int,      /* I - max number of ranges to be returned */
    numranges: &mut c_int, /* O - number ranges returned */
    minrow: &mut [c_long], /* O - first row in each range */
    maxrow: &mut [c_long], /* O - last row in each range */
    status: &mut c_int,    /* IO - status value */
) -> c_int {
    let mut next: *mut c_char;
    let mut minval: c_long = 0;
    let mut maxval: c_long = 0;

    if *status > 0 {
        return *status;
    }

    if maxrows <= 0 {
        *status = RANGE_PARSE_ERROR;
        ffpmsg_str("Input maximum range value is <= 0 (fits_parse_ranges)");
        return *status;
    }

    //next = rowlist;
    let mut ni = 0;

    *numranges = 0;

    while rowlist[ni] == bb(b' ') {
        ni += 1; /* skip spaces */
    }

    while rowlist[ni] != 0 {
        /* find min value of next range; rowlist[ni] must be '-' or a digit */
        if rowlist[ni] == bb(b'-') {
            minval = 1; /* implied minrow value = 1 */
        } else if isdigit_safe(rowlist[ni]) {
            minval = strtol_safe(&rowlist[ni..], &mut ni, 10);
        } else {
            *status = RANGE_PARSE_ERROR;
            ffpmsg_str("Syntax error in this row range list:");
            ffpmsg_slice(rowlist);
            return *status;
        }

        while rowlist[ni] == bb(b' ') {
            ni += 1; /* skip spaces */
        }

        /* find max value of next range; rowlist[ni] must be '-', or ',' */
        if rowlist[ni] == bb(b'-') {
            ni += 1;
            while rowlist[ni] == bb(b' ') {
                ni += 1; /* skip spaces */
            }

            if isdigit_safe(rowlist[ni]) {
                maxval = strtol_safe(&rowlist[ni..], &mut ni, 10);
            } else if rowlist[ni] == bb(b',') || rowlist[ni] == 0 {
                maxval = maxrows as c_long; /* implied max value */
            } else {
                *status = RANGE_PARSE_ERROR;
                ffpmsg_str("Syntax error in this row range list:");
                ffpmsg_slice(rowlist);
                return *status;
            }
        } else if rowlist[ni] == bb(b',') || rowlist[ni] == 0 {
            maxval = minval; /* only a single integer in this range */
        } else {
            *status = RANGE_PARSE_ERROR;
            ffpmsg_str("Syntax error in this row range list:");
            ffpmsg_slice(rowlist);
            return *status;
        }

        if *numranges + 1 > maxranges {
            *status = RANGE_PARSE_ERROR;
            ffpmsg_str("Overflowed maximum number of ranges (fits_parse_ranges)");
            return *status;
        }

        if minval < 1 {
            *status = RANGE_PARSE_ERROR;
            ffpmsg_str("Syntax error in this row range list: row number < 1");
            ffpmsg_slice(rowlist);
            return *status;
        }

        if maxval < minval {
            *status = RANGE_PARSE_ERROR;
            ffpmsg_str("Syntax error in this row range list: min > max");
            ffpmsg_slice(rowlist);
            return *status;
        }

        if *numranges > 0 && minval <= maxrow[(*numranges as usize) - 1] {
            *status = RANGE_PARSE_ERROR;
            ffpmsg_str("Syntax error in this row range list.  Range minimum is");
            ffpmsg_str("  less than or equal to previous range maximum");
            ffpmsg_slice(rowlist);
            return *status;
        }

        if minval as LONGLONG <= maxrows {
            /* ignore range if greater than maxrows */
            if maxval as LONGLONG > maxrows {
                maxval = maxrows as c_long;
            }

            minrow[(*numranges) as usize] = minval;
            maxrow[(*numranges) as usize] = maxval;

            (*numranges) += 1;
        }

        while rowlist[ni] == bb(b' ') {
            ni += 1; /* skip spaces */
        }
        if rowlist[ni] == bb(b',') {
            ni += 1;
            while rowlist[ni] == bb(b' ') {
                ni += 1; /* skip more spaces */
            }
        }
    }

    if *numranges == 0 {
        /* a null string was entered */
        minrow[0] = 1;
        maxrow[0] = maxrows as c_long;
        *numranges = 1;
    }

    *status
}

/*--------------------------------------------------------------------------*/
// parse the input list of row ranges, returning the number of ranges,
// and the min and max row value in each range.
//
// The only characters allowed in the input rowlist are
// decimal digits, minus sign, and comma (and non-significant spaces)
//
// Example:
//
// list = "10-20, 30-35,50"
//
// would return numranges = 3, minrow[] = {10, 30, 50}, maxrow[] = {20, 35, 50}
//
// error is returned if min value of range is > max value of range or if the
// ranges are not monotonically increasing.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffrwrgll(
    rowlist: *const c_char, /* I - list of rows and row ranges */
    maxrows: LONGLONG,      /* I - number of rows in the list */
    maxranges: c_int,       /* I - max number of ranges to be returned */
    numranges: *mut c_int,  /* O - number ranges returned */
    minrow: *mut LONGLONG,  /* O - first row in each range */
    maxrow: *mut LONGLONG,  /* O - last row in each range */
    status: *mut c_int,     /* IO - status value */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let numranges = numranges.as_mut().expect(NULL_MSG);
        let minrow = slice::from_raw_parts_mut(minrow, maxranges as usize);
        let maxrow = slice::from_raw_parts_mut(maxrow, maxranges as usize);

        raw_to_slice!(rowlist);

        ffrwrgll_safe(
            rowlist, maxrows, maxranges, numranges, minrow, maxrow, status,
        )
    }
}

/*--------------------------------------------------------------------------*/
// parse the input list of row ranges, returning the number of ranges,
// and the min and max row value in each range.
//
// The only characters allowed in the input rowlist are
// decimal digits, minus sign, and comma (and non-significant spaces)
//
// Example:
//
// list = "10-20, 30-35,50"
//
// would return numranges = 3, minrow[] = {10, 30, 50}, maxrow[] = {20, 35, 50}
//
// error is returned if min value of range is > max value of range or if the
// ranges are not monotonically increasing.
pub(crate) fn ffrwrgll_safe(
    rowlist: &[c_char],      /* I - list of rows and row ranges */
    maxrows: LONGLONG,       /* I - number of rows in the list */
    maxranges: c_int,        /* I - max number of ranges to be returned */
    numranges: &mut c_int,   /* O - number ranges returned */
    minrow: &mut [LONGLONG], /* O - first row in each range */
    maxrow: &mut [LONGLONG], /* O - last row in each range */
    status: &mut c_int,      /* IO - status value */
) -> c_int {
    let mut minval: LONGLONG = 0;
    let mut maxval: LONGLONG = 0;
    let mut dvalue: f64 = 0.0;

    if *status > 0 {
        return *status;
    }

    if maxrows <= 0 {
        *status = RANGE_PARSE_ERROR;
        ffpmsg_str("Input maximum range value is <= 0 (fits_parse_ranges)");
        return *status;
    }

    let mut ni = 0;
    //next = rowlist;
    *numranges = 0;

    while rowlist[ni] == bb(b' ') {
        ni += 1; /* skip spaces */
    }

    while rowlist[ni] != 0 {
        /* find min value of next range; rowlist[ni] must be '-' or a digit */
        if rowlist[ni] == bb(b'-') {
            minval = 1; /* implied minrow value = 1 */
        } else if isdigit_safe(rowlist[ni]) {
            /* read as a double, because the string to LONGLONG function */
            /* is platform dependent (strtoll, strtol, _atoI64)          */

            dvalue = strtod_safe(&rowlist[ni..], &mut ni);
            minval = (dvalue + 0.1) as LONGLONG;
        } else {
            *status = RANGE_PARSE_ERROR;
            ffpmsg_str("Syntax error in this row range list:");
            ffpmsg_slice(rowlist);
            return *status;
        }

        while rowlist[ni] == bb(b' ') {
            ni += 1; /* skip spaces */
        }

        /* find max value of next range; rowlist[ni] must be '-', or ',' */
        if rowlist[ni] == bb(b'-') {
            ni += 1;
            while rowlist[ni] == bb(b' ') {
                ni += 1; /* skip spaces */
            }

            if isdigit_safe(rowlist[ni]) {
                /* read as a double, because the string to LONGLONG function */
                /* is platform dependent (strtoll, strtol, _atoI64)          */

                dvalue = strtod_safe(&rowlist[ni..], &mut ni);
                maxval = (dvalue + 0.1) as LONGLONG;
            } else if rowlist[ni] == bb(b',') || rowlist[ni] == 0 {
                maxval = maxrows; /* implied max value */
            } else {
                *status = RANGE_PARSE_ERROR;
                ffpmsg_str("Syntax error in this row range list:");
                ffpmsg_slice(rowlist);
                return *status;
            }
        } else if rowlist[ni] == bb(b',') || rowlist[ni] == 0 {
            maxval = minval; /* only a single integer in this range */
        } else {
            *status = RANGE_PARSE_ERROR;
            ffpmsg_str("Syntax error in this row range list:");
            ffpmsg_slice(rowlist);
            return *status;
        }

        if *numranges + 1 > maxranges {
            *status = RANGE_PARSE_ERROR;
            ffpmsg_str("Overflowed maximum number of ranges (fits_parse_ranges)");
            return *status;
        }

        if minval < 1 {
            *status = RANGE_PARSE_ERROR;
            ffpmsg_str("Syntax error in this row range list: row number < 1");
            ffpmsg_slice(rowlist);
            return *status;
        }

        if maxval < minval {
            *status = RANGE_PARSE_ERROR;
            ffpmsg_str("Syntax error in this row range list: min > max");
            ffpmsg_slice(rowlist);
            return *status;
        }

        if *numranges > 0 && minval <= maxrow[(*numranges as usize) - 1] {
            *status = RANGE_PARSE_ERROR;
            ffpmsg_str("Syntax error in this row range list.  Range minimum is");
            ffpmsg_str("  less than or equal to previous range maximum");
            ffpmsg_slice(rowlist);
            return *status;
        }

        if minval <= maxrows {
            /* ignore range if greater than maxrows */
            if maxval > maxrows {
                maxval = maxrows;
            }

            minrow[(*numranges) as usize] = minval;
            maxrow[(*numranges) as usize] = maxval;

            (*numranges) += 1;
        }

        while rowlist[ni] == bb(b' ') {
            ni += 1; /* skip spaces */
        }
        if rowlist[ni] == bb(b',') {
            ni += 1;
            while rowlist[ni] == bb(b' ') {
                ni += 1; /* skip more spaces */
            }
        }
    }

    if *numranges == 0 {
        /* a null string was entered */
        minrow[0] = 1;
        maxrow[0] = maxrows;
        *numranges = 1;
    }

    *status
}

/*--------------------------------------------------------------------------*/
/// Insert a new column into an existing table at position numcol.  If
/// numcol is greater than the number of existing columns in the table
/// then the new column will be appended as the last column in the table.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn fficol(
    fptr: *mut fitsfile,  /* I - FITS file pointer                        */
    numcol: c_int,        /* I - position for new col. (1 = 1st)          */
    ttype: *const c_char, /* I - name of column (TTYPE keyword)           */
    tform: *const c_char, /* I - format of column (TFORM keyword)         */
    status: *mut c_int,   /* IO - error status                            */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        raw_to_slice!(ttype);
        raw_to_slice!(tform);

        fficol_safe(fptr, numcol, ttype, tform, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Insert a new column into an existing table at position numcol.  If
/// numcol is greater than the number of existing columns in the table
/// then the new column will be appended as the last column in the table.
pub(crate) fn fficol_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer                        */
    numcol: c_int,       /* I - position for new col. (1 = 1st)          */
    ttype: &[c_char],    /* I - name of column (TTYPE keyword)           */
    tform: &[c_char],    /* I - format of column (TFORM keyword)         */
    status: &mut c_int,  /* IO - error status                            */
) -> c_int {
    let name = ttype;
    let format = tform;

    fficls_safe(fptr, numcol, 1, &[name], &[format], status);
    *status
}

/*--------------------------------------------------------------------------*/
/// Insert 1 or more new columns into an existing table at position numcol.  If
/// fstcol is greater than the number of existing columns in the table
/// then the new column will be appended as the last column in the table.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn fficls(
    fptr: *mut fitsfile,         /* I - FITS file pointer                        */
    fstcol: c_int,               /* I - position for first new col. (1 = 1st)    */
    ncols: c_int,                /* I - number of columns to insert              */
    ttype: *const *const c_char, /* I - array of column names(TTYPE keywords)    */
    tform: *const *const c_char, /* I - array of formats of column (TFORM)       */
    status: *mut c_int,          /* IO - error status                            */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        let ttype = slice::from_raw_parts(ttype, ncols as usize);
        let tform = slice::from_raw_parts(tform, ncols as usize);

        let mut v_tform = Vec::new();

        for item in tform {
            let tform_item = slice::from_raw_parts(*item, FLEN_VALUE);
            v_tform.push(tform_item);
        }

        let mut v_ttype = Vec::new();

        for item in ttype {
            let tform_item = slice::from_raw_parts(*item, FLEN_VALUE);
            v_ttype.push(tform_item);
        }

        fficls_safe(fptr, fstcol, ncols, &v_ttype, &v_tform, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Insert 1 or more new columns into an existing table at position numcol.  If
/// fstcol is greater than the number of existing columns in the table
/// then the new column will be appended as the last column in the table.
pub(crate) fn fficls_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer                        */
    fstcol: c_int,       /* I - position for first new col. (1 = 1st)    */
    ncols: c_int,        /* I - number of columns to insert              */
    ttype: &[&[c_char]], /* I - array of column names(TTYPE keywords)    */
    tform: &[&[c_char]], /* I - array of formats of column (TFORM)       */
    status: &mut c_int,  /* IO - error status                            */
) -> c_int {
    let mut colnum: c_int = 0;
    let mut datacode: c_int = 0;
    let mut decims: c_int = 0;
    let mut tfields: c_int = 0;
    let mut tstatus: c_int = 0;
    let ii: c_int = 0;
    let mut datasize: LONGLONG = 0;
    let mut firstbyte: LONGLONG = 0;
    let mut nbytes: LONGLONG = 0;
    let mut nadd: LONGLONG = 0;
    let mut naxis1: LONGLONG = 0;
    let mut naxis2: LONGLONG = 0;
    let mut freespace: LONGLONG = 0;
    let mut tbcol: LONGLONG = 0;
    let mut firstcol: LONGLONG = 0;
    let mut delbyte: LONGLONG = 0;
    let mut nblock: c_long = 0;
    let mut width: c_long = 0;
    let mut repeat: c_long = 0;
    let mut tfm: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut keyname: [c_char; FLEN_KEYWORD] = [0; FLEN_KEYWORD];
    let mut comm: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];
    let mut cptr: *mut c_char;
    let mut card: [c_char; FLEN_CARD] = [0; FLEN_CARD];
    let mut colptr: *mut tcolumn;

    if *status > 0 {
        return *status;
    }

    if fptr.HDUposition != fptr.Fptr.curhdu {
        ffmahd_safe(fptr, (fptr.HDUposition) + 1, None, status);
    }
    /* rescan header if data structure is undefined */
    else if fptr.Fptr.datastart == DATA_UNDEFINED as LONGLONG && ffrdef_safe(fptr, status) > 0 {
        return *status;
    }

    if fptr.Fptr.hdutype == IMAGE_HDU {
        ffpmsg_str("Can only add columns to TABLE or BINTABLE extension (fficls)");
        *status = NOT_TABLE;
        return *status;
    }

    /*  is the column number valid?  */
    tfields = fptr.Fptr.tfield;
    if fstcol < 1 {
        *status = BAD_COL_NUM;
        return *status;
    } else if fstcol > tfields {
        colnum = tfields + 1; /* append as last column */
    } else {
        colnum = fstcol;
    }

    /* parse the tform value and calc number of bytes to add to each row */
    delbyte = 0;
    for ii in 0..(ncols as usize) {
        let t = tform[ii];

        if strlen_safe(t) > FLEN_VALUE - 1 {
            ffpmsg_str("Column format string too long (fficls)");
            *status = BAD_TFORM;
            return *status;
        }
        strcpy_safe(&mut tfm, t);
        ffupch_safe(&mut tfm); /* make sure format is in upper case */

        if fptr.Fptr.hdutype == ASCII_TBL {
            ffasfm_safe(
                &tfm,
                Some(&mut datacode),
                Some(&mut width),
                Some(&mut decims),
                status,
            );
            delbyte += width as LONGLONG + 1; /*  add one space between the columns */
        } else {
            ffbnfm_safe(
                &tfm,
                Some(&mut datacode),
                Some(&mut repeat),
                Some(&mut width),
                status,
            );

            if datacode < 0 {
                /* variable length array column */
                match strchr_safe(&tfm, bb(b'Q')) {
                    Some(x) => {
                        delbyte += 16;
                    }
                    None => {
                        delbyte += 8;
                    }
                }
            } else if datacode == 1 {
                /* bit column; round up  */
                delbyte += (repeat as LONGLONG + 7) / 8; /* to multiple of 8 bits */
            } else if datacode == 16 {
                /* ASCII string column */
                delbyte += repeat as LONGLONG;
            } else {
                /* numerical data type */
                delbyte += (datacode as LONGLONG / 10) * repeat as LONGLONG;
            }
        }
    }

    if *status > 0 {
        return *status;
    }

    /* get the current size of the table */
    /* use internal structure since NAXIS2 keyword may not be up to date */
    naxis1 = fptr.Fptr.rowlength;
    naxis2 = fptr.Fptr.numrows;

    /* current size of data */
    datasize = fptr.Fptr.heapstart + fptr.Fptr.heapsize;
    freespace = (((datasize + 2879) / 2880) * 2880) - datasize;
    nadd = delbyte * naxis2; /* no. of bytes to add to table */

    if (freespace - nadd) < 0
    /* not enough existing space? */
    {
        nblock = ((nadd - freespace + 2879) / 2880) as c_long; /* number of blocks  */
        if ffiblk(fptr, nblock, 1, status) > 0 {
            /* insert the blocks */
            return *status;
        }
    }

    /* shift heap down (if it exists) */
    if fptr.Fptr.heapsize > 0 {
        nbytes = fptr.Fptr.heapsize; /* no. of bytes to shift down */

        /* absolute heap pos */
        firstbyte = fptr.Fptr.datastart + fptr.Fptr.heapstart;

        if ffshft_safe(fptr, firstbyte, nbytes, nadd, status) > 0 {
            /* move heap */
            return *status;
        }
    }

    /* update the heap starting address */
    fptr.Fptr.heapstart += nadd;

    /* update the THEAP keyword if it exists */
    tstatus = 0;
    ffmkyj_safe(
        fptr,
        cs!(c"THEAP"),
        fptr.Fptr.heapstart,
        Some(cs!(c"&")),
        &mut tstatus,
    );

    /* calculate byte position in the row where to insert the new column */
    if colnum > tfields {
        firstcol = naxis1;
    } else {
        let colptr = fptr.Fptr.tableptr; /* point to first column structure */
        let c = fptr.Fptr.get_tableptr_as_slice();
        let ci = (colnum - 1) as usize; /* offset to the correct column */
        firstcol = c[ci].tbcol;
    }

    /* insert delbyte bytes in every row, at byte position firstcol */
    ffcins_safe(fptr, naxis1, naxis2, delbyte, firstcol, status);

    if fptr.Fptr.hdutype == ASCII_TBL {
        /* adjust the TBCOL values of the existing columns */
        for ii in 0..(tfields as usize) {
            ffkeyn_safe(cs!(c"TBCOL"), (ii + 1) as c_int, &mut keyname, status);
            ffgkyjj_safe(fptr, &keyname, &mut tbcol, Some(&mut comm), status);
            if tbcol > firstcol {
                tbcol += delbyte;
                ffmkyj_safe(fptr, &keyname, tbcol, Some(cs!(c"&")), status);
            }
        }
    }

    /* update the mandatory keywords */
    ffmkyj_safe(
        fptr,
        cs!(c"TFIELDS"),
        (tfields + ncols) as _,
        Some(cs!(c"&")),
        status,
    );
    ffmkyj_safe(
        fptr,
        cs!(c"NAXIS1"),
        naxis1 + delbyte,
        Some(cs!(c"&")),
        status,
    );

    /* increment the index value on any existing column keywords */
    if colnum <= tfields {
        ffkshf_safe(fptr, colnum, tfields, ncols, status);
    }

    /* add the required keywords for the new columns */
    for ii in 0..(ncols as usize) {
        let tf = tform[ii];
        let ttype_item = ttype[ii];

        strcpy_safe(&mut comm, cs!(c"label for field"));
        ffkeyn_safe(cs!(c"TTYPE"), colnum, &mut keyname, status);
        ffpkys_safe(fptr, &keyname, ttype_item, Some(&comm), status);

        strcpy_safe(&mut comm, cs!(c"format of field"));
        strcpy_safe(&mut tfm, tf);
        ffupch_safe(&mut tfm); /* make sure format is in upper case */
        ffkeyn_safe(cs!(c"TFORM"), colnum, &mut keyname, status);

        if (datacode.abs()) == TSBYTE {
            /* Replace the 'S' with an 'B' in the TFORMn code */
            let mut cptr = 0;
            while tfm[cptr] != bb(b'S') {
                cptr += 1;
            }

            tfm[cptr] = bb(b'B');
            ffpkys_safe(fptr, &keyname, &tfm, Some(&comm), status);

            /* write the TZEROn and TSCALn keywords */
            ffkeyn_safe(cs!(c"TZERO"), colnum, &mut keyname, status);
            strcpy_safe(&mut comm, cs!(c"offset for signed bytes"));

            ffpkyg_safe(fptr, &keyname, -128.0, 0, Some(&comm), status);

            ffkeyn_safe(cs!(c"TSCAL"), colnum, &mut keyname, status);
            strcpy_safe(&mut comm, cs!(c"data are not scaled"));
            ffpkyg_safe(fptr, &keyname, 1., 0, Some(&comm), status);
        } else if (datacode.abs()) == TUSHORT {
            /* Replace the 'U' with an 'I' in the TFORMn code */
            let mut cptr = 0;
            while tfm[cptr] != bb(b'U') {
                cptr += 1;
            }

            tfm[cptr] = bb(b'I');
            ffpkys_safe(fptr, &keyname, &tfm, Some(&comm), status);

            /* write the TZEROn and TSCALn keywords */
            ffkeyn_safe(cs!(c"TZERO"), colnum, &mut keyname, status);
            strcpy_safe(&mut comm, cs!(c"offset for unsigned integers"));

            ffpkyg_safe(fptr, &keyname, 32768., 0, Some(&comm), status);

            ffkeyn_safe(cs!(c"TSCAL"), colnum, &mut keyname, status);
            strcpy_safe(&mut comm, cs!(c"data are not scaled"));
            ffpkyg_safe(fptr, &keyname, 1., 0, Some(&comm), status);
        } else if (datacode.abs()) == TULONG {
            /* Replace the 'V' with an 'J' in the TFORMn code */
            let mut cptr = 0;
            while tfm[cptr] != bb(b'V') {
                cptr += 1;
            }

            tfm[cptr] = bb(b'J');
            ffpkys_safe(fptr, &keyname, &tfm, Some(&comm), status);

            /* write the TZEROn and TSCALn keywords */
            ffkeyn_safe(cs!(c"TZERO"), colnum, &mut keyname, status);
            strcpy_safe(&mut comm, cs!(c"offset for unsigned integers"));

            ffpkyg_safe(fptr, &keyname, 2147483648., 0, Some(&comm), status);

            ffkeyn_safe(cs!(c"TSCAL"), colnum, &mut keyname, status);
            strcpy_safe(&mut comm, cs!(c"data are not scaled"));
            ffpkyg_safe(fptr, &keyname, 1., 0, Some(&comm), status);
        } else if (datacode.abs()) == TULONGLONG {
            /* Replace the 'W' with an 'K' in the TFORMn code */
            let mut cptr = 0;
            while tfm[cptr] != bb(b'W') {
                cptr += 1;
            }

            tfm[cptr] = bb(b'K');
            ffpkys_safe(fptr, &keyname, &tfm, Some(&comm), status);

            /* write the TZEROn and TSCALn keywords */
            ffkeyn_safe(cs!(c"TZERO"), colnum, &mut card, status);
            strcat_safe(&mut card, cs!(c"     ")); /* make sure name is >= 8 chars long */
            card[8] = 0;
            strcat_safe(
                &mut card,
                cs!(c"=  9223372036854775808 / offset for unsigned integers"),
            );
            ffprec_safe(fptr, &card, status);

            ffkeyn_safe(cs!(c"TSCAL"), colnum, &mut keyname, status);
            strcpy_safe(&mut comm, cs!(c"data are not scaled"));
            ffpkyg_safe(fptr, &keyname, 1., 0, Some(&comm), status);
        } else {
            ffpkys_safe(fptr, &keyname, &tfm, Some(&comm), status);
        }

        if fptr.Fptr.hdutype == ASCII_TBL
        /* write the TBCOL keyword */
        {
            if colnum == tfields + 1 {
                tbcol = firstcol + 2; /* allow space between preceding col */
            } else {
                tbcol = firstcol + 1;
            }

            strcpy_safe(&mut comm, cs!(c"beginning column of field"));
            ffkeyn_safe(cs!(c"TBCOL"), colnum, &mut keyname, status);
            ffpkyj_safe(fptr, &keyname, tbcol, Some(&comm), status);

            /* increment the column starting position for the next column */
            ffasfm_safe(
                &tfm,
                Some(&mut datacode),
                Some(&mut width),
                Some(&mut decims),
                status,
            );
            firstcol += width as LONGLONG + 1; /*  add one space between the columns */
        }
        colnum += 1;
    }
    ffrdef_safe(fptr, status); /* initialize the new table structure */
    *status
}

/*--------------------------------------------------------------------------*/
/// Modify the vector length of a column in a binary table, larger or smaller.
/// E.g., change a column from TFORMn = '1E' to '20E'.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffmvec(
    fptr: *mut fitsfile, /* I - FITS file pointer                        */
    colnum: c_int,       /* I - position of col to be modified           */
    newveclen: LONGLONG, /* I - new vector length of column (TFORM)       */
    status: *mut c_int,  /* IO - error status                            */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        ffmvec_safe(fptr, colnum, newveclen, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Modify the vector length of a column in a binary table, larger or smaller.
/// E.g., change a column from TFORMn = '1E' to '20E'.
pub(crate) fn ffmvec_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer                        */
    colnum: c_int,       /* I - position of col to be modified           */
    newveclen: LONGLONG, /* I - new vector length of column (TFORM)       */
    status: &mut c_int,  /* IO - error status                            */
) -> c_int {
    let mut datacode: c_int = 0;
    let mut tfields: c_int = 0;
    let mut tstatus: c_int = 0;
    let mut datasize: LONGLONG = 0;
    let mut size: LONGLONG = 0;
    let mut firstbyte: LONGLONG = 0;
    let mut nbytes: LONGLONG = 0;
    let mut nadd: LONGLONG = 0;
    let mut ndelete: LONGLONG = 0;
    let mut naxis1: LONGLONG = 0;
    let mut naxis2: LONGLONG = 0;
    let mut firstcol: LONGLONG = 0;
    let mut freespace: LONGLONG = 0;
    let mut width: LONGLONG = 0;
    let mut delbyte: LONGLONG = 0;
    let mut repeat: LONGLONG = 0;
    let mut nblock: c_long = 0;
    let mut tfm: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut keyname: [c_char; FLEN_KEYWORD] = [0; FLEN_KEYWORD];
    let mut tcode: [c_char; 2] = [0; 2];
    let mut colptr: *mut tcolumn;

    if *status > 0 {
        return *status;
    }

    if fptr.HDUposition != fptr.Fptr.curhdu {
        ffmahd_safe(fptr, (fptr.HDUposition) + 1, None, status);
    }
    /* rescan header if data structure is undefined */
    else if fptr.Fptr.datastart == DATA_UNDEFINED as LONGLONG && ffrdef_safe(fptr, status) > 0 {
        return *status;
    }

    if fptr.Fptr.hdutype != BINARY_TBL {
        ffpmsg_str("Can only change vector length of a column in BINTABLE extension (ffmvec)");
        *status = NOT_TABLE;
        return *status;
    }

    /*  is the column number valid?  */
    tfields = fptr.Fptr.tfield;
    if colnum < 1 || colnum > tfields {
        *status = BAD_COL_NUM;
        return *status;
    }

    /* look up the current vector length and element width */

    let colptr = fptr.Fptr.tableptr; /* point to first column structure */
    let c = fptr.Fptr.get_tableptr_as_slice();
    let ci = (colnum - 1) as usize; /* offset to the correct column */

    datacode = c[ci].tdatatype; /* datatype of the column */
    repeat = c[ci].trepeat; /* field repeat count  */
    width = c[ci].twidth as LONGLONG; /*  width of a single element in chars */

    if datacode < 0 {
        ffpmsg_str("Can't modify vector length of variable length column (ffmvec)");
        *status = BAD_TFORM;
        return *status;
    }

    if repeat == newveclen {
        return *status; /* column already has the desired vector length */
    }

    if datacode == TSTRING {
        width = 1; /* width was equal to width of unit string */
    }

    naxis1 = fptr.Fptr.rowlength; /* current width of the table */
    naxis2 = fptr.Fptr.numrows;

    delbyte = (newveclen - repeat) * width; /* no. of bytes to insert */
    if datacode == TBIT {
        /* BIT column is a special case */
        delbyte = ((newveclen + 7) / 8) - ((repeat + 7) / 8);
    }

    if delbyte > 0
    /* insert space for more elements */
    {
        /* current size of data */
        datasize = fptr.Fptr.heapstart + fptr.Fptr.heapsize;
        freespace = (((datasize + 2879) / 2880) * 2880) - datasize;

        nadd = (delbyte as LONGLONG) * naxis2; /* no. of bytes to add to table */

        if (freespace - nadd) < 0 {
            /* not enough existing space? */

            nblock = ((nadd - freespace + 2879) / 2880) as c_long; /* number of blocks  */
            if ffiblk(fptr, nblock, 1, status) > 0 {
                /* insert the blocks */
                return *status;
            }
        }

        /* shift heap down (if it exists) */
        if fptr.Fptr.heapsize > 0 {
            nbytes = fptr.Fptr.heapsize; /* no. of bytes to shift down */

            /* absolute heap pos */
            firstbyte = fptr.Fptr.datastart + fptr.Fptr.heapstart;

            if ffshft_safe(fptr, firstbyte, nbytes, nadd, status) > 0 {
                /* move heap */
                return *status;
            }
        }

        /* update the heap starting address */
        fptr.Fptr.heapstart += nadd;

        /* update the THEAP keyword if it exists */
        tstatus = 0;
        ffmkyj_safe(
            fptr,
            cs!(c"THEAP"),
            fptr.Fptr.heapstart,
            Some(cs!(c"&")),
            &mut tstatus,
        );

        /* Must reset colptr before using it again.  fptr.Fptr.tableptr
        may have been reallocated down in ffbinit via the call to ffiblk above.*/
        let colptr = fptr.Fptr.tableptr; /* point to first column structure */
        let c = fptr.Fptr.get_tableptr_as_slice();
        let ci = (colnum - 1) as usize; /* offset to the correct column */

        firstcol = c[ci].tbcol + (repeat * width); /* insert position */

        /* insert delbyte bytes in every row, at byte position firstcol */
        ffcins_safe(fptr, naxis1, naxis2, delbyte, firstcol, status);
    } else if delbyte < 0 {
        /* current size of table */
        size = fptr.Fptr.heapstart + fptr.Fptr.heapsize;
        freespace = ((size + 2879) / 2880) * 2880 - size - ((delbyte as LONGLONG) * naxis2);
        nblock = (freespace / 2880) as c_long; /* number of empty blocks to delete */
        firstcol = c[ci].tbcol + (newveclen * width); /* delete position */

        /* delete elements from the vector */
        ffcdel_safe(fptr, naxis1, naxis2, -delbyte, firstcol, status);

        /* abs heap pos */
        firstbyte = fptr.Fptr.datastart + fptr.Fptr.heapstart;
        ndelete = (delbyte as LONGLONG) * naxis2; /* size of shift (negative) */

        /* shift heap up (if it exists) */
        if fptr.Fptr.heapsize > 0 {
            nbytes = fptr.Fptr.heapsize; /* no. of bytes to shift up */
            if ffshft_safe(fptr, firstbyte, nbytes, ndelete, status) > 0 {
                return *status;
            }
        }

        /* delete the empty  blocks at the end of the HDU */
        if nblock > 0 {
            ffdblk(fptr, nblock, status);
        }

        /* update the heap starting address */
        fptr.Fptr.heapstart += ndelete; /* ndelete is negative */

        /* update the THEAP keyword if it exists */
        tstatus = 0;
        ffmkyj_safe(
            fptr,
            cs!(c"THEAP"),
            fptr.Fptr.heapstart,
            Some(cs!(c"&")),
            &mut tstatus,
        );
    }

    /* construct the new TFORM keyword for the column */
    if datacode == TBIT {
        strcpy_safe(&mut tcode, cs!(c"X"));
    } else if datacode == TBYTE {
        strcpy_safe(&mut tcode, cs!(c"B"));
    } else if datacode == TLOGICAL {
        strcpy_safe(&mut tcode, cs!(c"L"));
    } else if datacode == TSTRING {
        strcpy_safe(&mut tcode, cs!(c"A"));
    } else if datacode == TSHORT {
        strcpy_safe(&mut tcode, cs!(c"I"));
    } else if datacode == TLONG {
        strcpy_safe(&mut tcode, cs!(c"J"));
    } else if datacode == TLONGLONG {
        strcpy_safe(&mut tcode, cs!(c"K"));
    } else if datacode == TFLOAT {
        strcpy_safe(&mut tcode, cs!(c"E"));
    } else if datacode == TDOUBLE {
        strcpy_safe(&mut tcode, cs!(c"D"));
    } else if datacode == TCOMPLEX {
        strcpy_safe(&mut tcode, cs!(c"C"));
    } else if datacode == TDBLCOMPLEX {
        strcpy_safe(&mut tcode, cs!(c"M"));
    }

    /* write as a double value because the LONGLONG conversion */

    int_snprintf!(
        &mut tfm,
        FLEN_VALUE,
        "{:.0}{}",
        newveclen as f64,
        slice_to_str!(&tcode),
    );

    ffkeyn_safe(cs!(c"TFORM"), colnum, &mut keyname, status); /* Keyword name */
    ffmkys_safe(fptr, &keyname, &tfm, Some(cs!(c"&")), status); /* modify TFORM keyword */

    ffmkyj_safe(
        fptr,
        cs!(c"NAXIS1"),
        naxis1 + delbyte,
        Some(cs!(c"&")),
        status,
    ); /* modify NAXIS1 */

    ffrdef_safe(fptr, status); /* reinitialize the new table structure */
    *status
}

/*--------------------------------------------------------------------------*/
/// copy a column from infptr and insert it in the outfptr table.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffcpcl(
    infptr: *mut fitsfile,  /* I - FITS file pointer to input file  */
    outfptr: *mut fitsfile, /* I - FITS file pointer to output file */
    incol: c_int,           /* I - number of input column   */
    outcol: c_int,          /* I - number for output column  */
    create_col: c_int,      /* I - create new col if TRUE, else overwrite */
    status: *mut c_int,     /* IO - error status     */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let infptr = infptr.as_mut().expect(NULL_MSG);
        let outfptr = outfptr.as_mut().expect(NULL_MSG);

        ffcpcl_safe(infptr, outfptr, incol, outcol, create_col, status)
    }
}

/*--------------------------------------------------------------------------*/
/// copy a column from infptr and insert it in the outfptr table.
pub(crate) fn ffcpcl_safe(
    infptr: &mut fitsfile,  /* I - FITS file pointer to input file  */
    outfptr: &mut fitsfile, /* I - FITS file pointer to output file */
    incol: c_int,           /* I - number of input column   */
    outcol: c_int,          /* I - number for output column  */
    create_col: c_int,      /* I - create new col if TRUE, else overwrite */
    status: &mut c_int,     /* IO - error status     */
) -> c_int {
    let mut tstatus: c_int = 0;
    let mut colnum: c_int = 0;
    let mut typecode: c_int = 0;
    let mut otypecode: c_int = 0;
    let mut etypecode: c_int = 0;
    let mut anynull: c_int = 0;
    let mut inHduType: c_int = 0;
    let mut outHduType: c_int = 0;
    let mut tfields: c_long = 0;
    let mut repeat: c_long = 0;
    let mut orepeat: c_long = 0;
    let mut width: c_long = 0;
    let mut owidth: c_long = 0;
    let mut nrows: c_long = 0;
    let mut outrows: c_long = 0;
    let mut inloop: c_long = 0;
    let mut outloop: c_long = 0;
    let mut maxloop: c_long = 0;
    let mut ndone: c_long = 0;
    let mut ntodo: c_long = 0;
    let mut npixels: c_long = 0;
    let mut firstrow: c_long = 0;
    let mut firstelem: c_long = 0;
    let ii: c_long = 0;
    let mut keyname: [c_char; FLEN_KEYWORD] = [0; FLEN_KEYWORD];
    let mut ttype: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut tform: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut ttype_comm: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];
    let mut tform_comm: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];
    let mut lvalues: Vec<c_char> = Vec::new();
    let mut nullflag: c_char = 0;
    let mut strarray: Vec<Vec<c_char>> = Vec::new();
    let nulstr: [c_char; 2] = [5, 0];
    let mut dnull: f64 = 0.0;
    let mut dvalues: Vec<f64> = Vec::new();
    let mut fnull: f32 = 0.0;
    let mut fvalues: Vec<f32> = Vec::new();
    let typecodes: [c_int; 1000] = [0; 1000];

    let icol: c_int = 0;
    let incol1: c_int = 0;
    let outcol1: c_int = 0;

    let mut jjvalues: Vec<c_longlong> = Vec::new();
    let mut ujjvalues: Vec<c_ulonglong> = Vec::new();

    let mut incol = incol; // shadow as mut

    if *status > 0 {
        return *status;
    }

    if infptr.HDUposition != infptr.Fptr.curhdu {
        ffmahd_safe(infptr, (infptr.HDUposition) + 1, None, status);
    } else if infptr.Fptr.datastart == DATA_UNDEFINED as LONGLONG {
        ffrdef_safe(infptr, status); /* rescan header */
    }
    inHduType = infptr.Fptr.hdutype;

    if outfptr.HDUposition != outfptr.Fptr.curhdu {
        ffmahd_safe(outfptr, (outfptr.HDUposition) + 1, None, status);
    } else if outfptr.Fptr.datastart == DATA_UNDEFINED as LONGLONG {
        ffrdef_safe(outfptr, status); /* rescan header */
    }
    outHduType = outfptr.Fptr.hdutype;

    if *status > 0 {
        return *status;
    }

    if inHduType == IMAGE_HDU || outHduType == IMAGE_HDU {
        ffpmsg_str("Can not copy columns to or from IMAGE HDUs (ffcpcl)");
        *status = NOT_TABLE;
        return *status;
    }

    if inHduType == BINARY_TBL && outHduType == ASCII_TBL {
        ffpmsg_str("Copying from Binary table to ASCII table is not supported (ffcpcl)");
        *status = NOT_BTABLE;
        return *status;
    }

    /* get the datatype and vector repeat length of the column */
    ffgtcl_safe(
        infptr,
        incol,
        Some(&mut typecode),
        Some(&mut repeat),
        Some(&mut width),
        status,
    );

    /* ... and equivalent type code */
    ffeqty_safe(infptr, incol, Some(&mut etypecode), None, None, status);

    if typecode < 0 {
        ffpmsg_str("Variable-length columns are not supported (ffcpcl)");
        *status = BAD_TFORM;
        return *status;
    }

    if create_col != 0 {
        /* insert new column in output table? */
        tstatus = 0;
        ffkeyn_safe(cs!(c"TTYPE"), incol, &mut keyname, &mut tstatus);
        ffgkys_safe(
            infptr,
            &keyname,
            &mut ttype,
            Some(&mut ttype_comm),
            &mut tstatus,
        );
        ffkeyn_safe(cs!(c"TFORM"), incol, &mut keyname, &mut tstatus);

        if ffgkys_safe(
            infptr,
            &keyname,
            &mut tform,
            Some(&mut tform_comm),
            &mut tstatus,
        ) != 0
        {
            ffpmsg_str("Could not find TTYPE and TFORM keywords in input table (ffcpcl)");
            *status = NO_TFORM;
            return *status;
        }

        if inHduType == ASCII_TBL && outHduType == BINARY_TBL {
            /* convert from ASCII table to BINARY table format string */
            if typecode == TSTRING {
                ffnkey_safe(width as _, cs!(c"A"), &mut tform, status);
            } else if typecode == TLONG {
                strcpy_safe(&mut tform, cs!(c"1J"));
            } else if typecode == TSHORT {
                strcpy_safe(&mut tform, cs!(c"1I"));
            } else if typecode == TFLOAT {
                strcpy_safe(&mut tform, cs!(c"1E"));
            } else if typecode == TDOUBLE {
                strcpy_safe(&mut tform, cs!(c"1D"));
            }
        }

        if ffgkyj_safe(outfptr, cs!(c"TFIELDS"), &mut tfields, None, &mut tstatus) != 0 {
            ffpmsg_str("Could not read TFIELDS keyword in output table (ffcpcl)");
            *status = NO_TFIELDS;
            return *status;
        }

        colnum = cmp::min((tfields + 1) as c_int, outcol); /* output col. number */

        /* create the empty column */
        if fficol_safe(outfptr, colnum, &ttype, &tform, status) > 0 {
            ffpmsg_str("Could not append new column to output file (ffcpcl)");
            return *status;
        }

        if std::ptr::eq(infptr.Fptr.as_mut(), outfptr.Fptr.as_mut())
            && (infptr.HDUposition == outfptr.HDUposition)
            && (colnum <= incol)
        {
            incol += 1; /* the input column has been shifted over */
        }

        /* copy the comment strings from the input file for TTYPE and TFORM */
        tstatus = 0;
        ffkeyn_safe(cs!(c"TTYPE"), colnum, &mut keyname, &mut tstatus);
        ffmcom_safe(outfptr, &keyname, Some(&ttype_comm), &mut tstatus);
        ffkeyn_safe(cs!(c"TFORM"), colnum, &mut keyname, &mut tstatus);
        ffmcom_safe(outfptr, &keyname, Some(&tform_comm), &mut tstatus);

        /* copy other column-related keywords if they exist */

        ffcpky_safe(infptr, outfptr, incol, colnum, cs!(c"TUNIT"), status);
        ffcpky_safe(infptr, outfptr, incol, colnum, cs!(c"TSCAL"), status);
        ffcpky_safe(infptr, outfptr, incol, colnum, cs!(c"TZERO"), status);
        ffcpky_safe(infptr, outfptr, incol, colnum, cs!(c"TDISP"), status);
        ffcpky_safe(infptr, outfptr, incol, colnum, cs!(c"TLMIN"), status);
        ffcpky_safe(infptr, outfptr, incol, colnum, cs!(c"TLMAX"), status);
        ffcpky_safe(infptr, outfptr, incol, colnum, cs!(c"TDIM"), status);

        /*  WCS keywords */
        ffcpky_safe(infptr, outfptr, incol, colnum, cs!(c"TCTYP"), status);
        ffcpky_safe(infptr, outfptr, incol, colnum, cs!(c"TCUNI"), status);
        ffcpky_safe(infptr, outfptr, incol, colnum, cs!(c"TCRVL"), status);
        ffcpky_safe(infptr, outfptr, incol, colnum, cs!(c"TCRPX"), status);
        ffcpky_safe(infptr, outfptr, incol, colnum, cs!(c"TCDLT"), status);
        ffcpky_safe(infptr, outfptr, incol, colnum, cs!(c"TCROT"), status);

        if inHduType == ASCII_TBL && outHduType == BINARY_TBL {
            /* binary tables only have TNULLn keyword for integer columns */
            if typecode == TLONG || typecode == TSHORT {
                /* check if null string is defined; replace with integer */
                ffkeyn_safe(cs!(c"TNULL"), incol, &mut keyname, &mut tstatus);
                if ffgkys_safe(infptr, &keyname, &mut ttype, None, &mut tstatus) <= 0 {
                    ffkeyn_safe(cs!(c"TNULL"), colnum, &mut keyname, &mut tstatus);
                    if typecode == TLONG {
                        ffpkyj_safe(
                            outfptr,
                            &keyname,
                            -9999999,
                            Some(cs!(c"Null value")),
                            status,
                        );
                    } else {
                        ffpkyj_safe(outfptr, &keyname, -32768, Some(cs!(c"Null value")), status);
                    }
                }
            }
        } else {
            ffcpky_safe(infptr, outfptr, incol, colnum, cs!(c"TNULL"), status);
        }

        /* rescan header to recognize the new keywords */
        if ffrdef_safe(outfptr, status) != 0 {
            return *status;
        }
    } else {
        colnum = outcol;
        /* get the datatype and vector repeat length of the output column */
        ffgtcl_safe(
            outfptr,
            outcol,
            Some(&mut otypecode),
            Some(&mut orepeat),
            Some(&mut owidth),
            status,
        );

        if orepeat != repeat {
            ffpmsg_str("Input and output vector columns must have same length (ffcpcl)");
            *status = BAD_TFORM;
            return *status;
        }
    }

    ffgkyj_safe(infptr, cs!(c"NAXIS2"), &mut nrows, None, status); /* no. of input rows */
    ffgkyj_safe(outfptr, cs!(c"NAXIS2"), &mut outrows, None, status); /* no. of output rows */
    nrows = cmp::min(nrows, outrows);

    if typecode == TBIT {
        repeat = (repeat + 7) / 8; /* convert from bits to bytes */
    } else if typecode == TSTRING && inHduType == BINARY_TBL {
        repeat /= width; /* convert from chars to unit strings */
    }

    /* get optimum number of rows to copy at one time */
    ffgrsz_safe(infptr, &mut inloop, status);
    ffgrsz_safe(outfptr, &mut outloop, status);

    /* adjust optimum number, since 2 tables are open at once */
    maxloop = cmp::min(inloop, outloop); /* smallest of the 2 tables */
    maxloop = cmp::max(1, maxloop / 2); /* at least 1 row */
    maxloop = cmp::min(maxloop, nrows); /* max = nrows to be copied */
    maxloop *= repeat; /* mult by no of elements in a row */

    /* allocate memory for arrays */
    if typecode == TLOGICAL {
        if lvalues.try_reserve_exact(maxloop as usize).is_err() {
            ffpmsg_str("malloc failed to get memory for logicals (ffcpcl)");
            *status = ARRAY_TOO_BIG;
            return *status;
        } else {
            lvalues.resize(maxloop as usize, 0);
        }
    } else if typecode == TSTRING {
        /* allocate array of pointers */
        strarray.reserve_exact(maxloop as usize);

        /* allocate space for each string */
        for ii in 0..(maxloop as usize) {
            let str_storage = vec![0; (width + 1) as usize];

            strarray.push(str_storage);
        }
    } else if typecode == TCOMPLEX {
        if fvalues.try_reserve_exact(2 * maxloop as usize).is_err() {
            ffpmsg_str("malloc failed to get memory for complex (ffcpcl)");
            *status = ARRAY_TOO_BIG;
            return *status;
        } else {
            fvalues.resize(2 * maxloop as usize, 0.0);
        }
        fnull = 0.0;
    } else if typecode == TDBLCOMPLEX {
        if dvalues.try_reserve_exact(2 * maxloop as usize).is_err() {
            ffpmsg_str("malloc failed to get memory for dbl complex (ffcpcl)");
            *status = ARRAY_TOO_BIG;
            return *status;
        } else {
            dvalues.resize(2 * maxloop as usize, 0.0);
        }
        dnull = 0.0;
    } else if typecode == TLONGLONG && etypecode == TULONGLONG {
        /* These are unsigned long-long ints that are not rescaled to floating point numbers */

        if ujjvalues.try_reserve_exact(maxloop as usize).is_err() {
            ffpmsg_str("malloc failed to get memory for unsigned long long int (ffcpcl)");
            *status = ARRAY_TOO_BIG;
            return *status;
        } else {
            ujjvalues.resize(maxloop as usize, 0);
        }
    } else if typecode == TLONGLONG && etypecode != TDOUBLE {
        /* These are long-long ints that are not rescaled to floating point numbers */

        if jjvalues.try_reserve_exact(maxloop as usize).is_err() {
            ffpmsg_str("malloc failed to get memory for long long int (ffcpcl)");
            *status = ARRAY_TOO_BIG;
            return *status;
        } else {
            jjvalues.resize(maxloop as usize, 0);
        }
    } else {
        /* other numerical datatype; read them all as doubles */

        if dvalues.try_reserve_exact(maxloop as usize).is_err() {
            ffpmsg_str("malloc failed to get memory for doubles (ffcpcl)");
            *status = ARRAY_TOO_BIG;
            return *status;
        } else {
            dvalues.resize(maxloop as usize, 0.0);
        }
        dnull = -9.99991999E31; /* use an unlikely value for nulls */
    }

    npixels = nrows * repeat; /* total no. of pixels to copy */
    ntodo = cmp::min(npixels, maxloop); /* no. to copy per iteration */
    ndone = 0; /* total no. of pixels that have been copied */

    while ntodo != 0 {
        /* iterate through the table */
        firstrow = ndone / repeat + 1;
        firstelem = ndone - ((firstrow - 1) * repeat) + 1;

        /* read from input table */
        if typecode == TLOGICAL {
            ffgcvl_safe(
                infptr,
                incol,
                firstrow as LONGLONG,
                firstelem as LONGLONG,
                ntodo as LONGLONG,
                0,
                &mut lvalues,
                Some(&mut anynull),
                status,
            );
        } else if typecode == TSTRING {
            ffgcvs_safe(
                infptr,
                incol,
                firstrow as LONGLONG,
                firstelem as LONGLONG,
                ntodo as LONGLONG,
                Some(&nulstr),
                &mut vecs_to_slices_mut(&mut strarray),
                Some(&mut anynull),
                status,
            );
        } else if typecode == TCOMPLEX {
            ffgcvc_safe(
                infptr,
                incol,
                firstrow as LONGLONG,
                firstelem as LONGLONG,
                ntodo as LONGLONG,
                fnull,
                &mut fvalues,
                Some(&mut anynull),
                status,
            );
        } else if typecode == TDBLCOMPLEX {
            ffgcvm_safe(
                infptr,
                incol,
                firstrow as LONGLONG,
                firstelem as LONGLONG,
                ntodo as LONGLONG,
                dnull,
                &mut dvalues,
                Some(&mut anynull),
                status,
            );

        /* Neither TULONGLONG nor TLONGLONG does null checking.  Whatever
        null value is in input table is transferred to output table
        without checking.  Since the TNULL value was copied, this
        should preserve null values */
        } else if typecode == TLONGLONG && etypecode == TULONGLONG {
            ffgcvujj_safe(
                infptr,
                incol,
                firstrow as LONGLONG,
                firstelem as LONGLONG,
                ntodo as LONGLONG,
                /*nulval*/ 0,
                &mut ujjvalues,
                Some(&mut anynull),
                status,
            );
        } else if typecode == TLONGLONG && etypecode != TDOUBLE {
            ffgcvjj_safe(
                infptr,
                incol,
                firstrow as LONGLONG,
                firstelem as LONGLONG,
                ntodo as LONGLONG,
                /*nulval*/ 0,
                &mut jjvalues,
                Some(&mut anynull),
                status,
            );
        } else {
            /* all numerical types */
            ffgcvd_safe(
                infptr,
                incol,
                firstrow as LONGLONG,
                firstelem as LONGLONG,
                ntodo as LONGLONG,
                dnull,
                &mut dvalues,
                Some(&mut anynull),
                status,
            );
        }

        if *status > 0 {
            ffpmsg_str("Error reading input copy of column (ffcpcl)");
            break;
        }

        /* write to output table */
        if typecode == TLOGICAL {
            nullflag = 2;

            ffpcnl_safe(
                outfptr,
                colnum,
                firstrow as LONGLONG,
                firstelem as LONGLONG,
                ntodo as LONGLONG,
                &lvalues,
                nullflag,
                status,
            );
        } else if typecode == TSTRING {
            if anynull != 0 {
                ffpcns_safe(
                    outfptr,
                    colnum,
                    firstrow as LONGLONG,
                    firstelem as LONGLONG,
                    ntodo as LONGLONG,
                    &vecs_to_slices(&strarray),
                    &nulstr,
                    status,
                );
            } else {
                ffpcls_safe(
                    outfptr,
                    colnum,
                    firstrow as LONGLONG,
                    firstelem as LONGLONG,
                    ntodo as LONGLONG,
                    &vecs_to_slices(&strarray),
                    status,
                );
            }
        } else if typecode == TCOMPLEX {
            /* doesn't support writing nulls */
            ffpclc_safe(
                outfptr,
                colnum,
                firstrow as LONGLONG,
                firstelem as LONGLONG,
                ntodo as LONGLONG,
                &fvalues,
                status,
            );
        } else if typecode == TDBLCOMPLEX {
            /* doesn't support writing nulls */
            ffpclm_safe(
                outfptr,
                colnum,
                firstrow as LONGLONG,
                firstelem as LONGLONG,
                ntodo as LONGLONG,
                &dvalues,
                status,
            );
        } else if typecode == TLONGLONG && etypecode == TULONGLONG {
            /* No null checking because we did none to read */
            ffpclujj_safe(
                outfptr,
                colnum,
                firstrow as LONGLONG,
                firstelem as LONGLONG,
                ntodo as LONGLONG,
                &ujjvalues,
                status,
            );
        } else if typecode == TLONGLONG && etypecode != TDOUBLE {
            /* No null checking because we did none to read */
            ffpcljj_safe(
                outfptr,
                colnum,
                firstrow as LONGLONG,
                firstelem as LONGLONG,
                ntodo as LONGLONG,
                &jjvalues,
                status,
            );
        } else
        /* all other numerical types */
        if anynull != 0 {
            ffpcnd_safe(
                outfptr,
                colnum,
                firstrow as LONGLONG,
                firstelem as LONGLONG,
                ntodo as LONGLONG,
                &dvalues,
                dnull,
                status,
            );
        } else {
            ffpcld_safe(
                outfptr,
                colnum,
                firstrow as LONGLONG,
                firstelem as LONGLONG,
                ntodo as LONGLONG,
                &dvalues,
                status,
            );
        }

        if *status > 0 {
            ffpmsg_str("Error writing output copy of column (ffcpcl)");
            break;
        }

        npixels -= ntodo;
        ndone += ntodo;
        ntodo = cmp::min(npixels, maxloop);
    }

    *status
}

/*--------------------------------------------------------------------------*/
/// Copy multiple columns from infptr and insert them in the outfptr
/// table.  Optimized for multiple-column case since it only expands the
/// output file once using fits_insert_cols() instead of calling
/// fits_insert_col() multiple times.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffccls(
    infptr: *mut fitsfile,  /* I - FITS file pointer to input file  */
    outfptr: *mut fitsfile, /* I - FITS file pointer to output file */
    incol: c_int,           /* I - number of first input column   */
    outcol: c_int,          /* I - number for first output column  */
    ncols: c_int,           /* I - number of columns to copy from input to output */
    create_col: c_int,      /* I - create new col if TRUE, else overwrite */
    status: *mut c_int,     /* IO - error status     */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let infptr = infptr.as_mut().expect(NULL_MSG);
        let outfptr = outfptr.as_mut().expect(NULL_MSG);

        ffccls_safe(infptr, outfptr, incol, outcol, ncols, create_col, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Copy multiple columns from infptr and insert them in the outfptr
/// table.  Optimized for multiple-column case since it only expands the
/// output file once using fits_insert_cols() instead of calling
/// fits_insert_col() multiple times.
pub(crate) fn ffccls_safe(
    infptr: &mut fitsfile,  /* I - FITS file pointer to input file  */
    outfptr: &mut fitsfile, /* I - FITS file pointer to output file */
    incol: c_int,           /* I - number of first input column   */
    outcol: c_int,          /* I - number for first output column  */
    ncols: c_int,           /* I - number of columns to copy from input to output */
    create_col: c_int,      /* I - create new col if TRUE, else overwrite */
    status: &mut c_int,     /* IO - error status     */
) -> c_int {
    let mut tstatus: c_int = 0;
    let mut colnum: c_int = 0;
    let mut typecode: c_int = 0;
    let mut otypecode: c_int = 0;
    let anynull: c_int = 0;
    let mut inHduType: c_int = 0;
    let mut outHduType: c_int = 0;
    let mut tfields: c_long = 0;
    let mut repeat: c_long = 0;
    let mut orepeat: c_long = 0;
    let mut width: c_long = 0;
    let mut owidth: c_long = 0;
    let nrows: c_long = 0;
    let outrows: c_long = 0;
    let inloop: c_long = 0;
    let outloop: c_long = 0;
    let maxloop: c_long = 0;
    let ndone: c_long = 0;
    let ntodo: c_long = 0;
    let npixels: c_long = 0;
    let firstrow: c_long = 0;
    let firstelem: c_long = 0;
    let ii: c_long = 0;
    let mut keyname: [c_char; FLEN_KEYWORD] = [0; FLEN_KEYWORD];
    let mut ttype: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut tform: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut ttype_comm: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];
    let mut tform_comm: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];
    let lvalues: Vec<c_char> = Vec::new();
    let nullflag: c_char = 0;
    let strarray: Vec<*mut c_char> = Vec::new();
    let nulstr: [c_char; 2] = [5, 0];
    let dnull: f64 = 0.0;
    let dvalues: Vec<f64> = Vec::new();
    let fnull: f32 = 0.0;
    let fvalues: Vec<f32> = Vec::new();
    let mut typecodes: [c_int; 1000] = [0; 1000];
    let mut ttypes: [[c_char; FLEN_CARD]; 1000] = [[0; FLEN_CARD]; 1000];
    let mut tforms: [[c_char; FLEN_CARD]; 1000] = [[0; FLEN_CARD]; 1000];

    let ikey: c_int = 0;
    let jkey = 0;
    let icol: c_int = 0;
    let incol1: c_int = 0;
    let outcol1: c_int = 0;

    if *status > 0 {
        return *status;
    }

    /* Do not allow more than internal array limit to be copied */
    if ncols > 1000 {
        *status = ARRAY_TOO_BIG;
        return *status;
    }

    if infptr.HDUposition != infptr.Fptr.curhdu {
        ffmahd_safe(infptr, (infptr.HDUposition) + 1, None, status);
    } else if infptr.Fptr.datastart == DATA_UNDEFINED as LONGLONG {
        ffrdef_safe(infptr, status); /* rescan header */
    }

    inHduType = infptr.Fptr.hdutype;

    if outfptr.HDUposition != outfptr.Fptr.curhdu {
        ffmahd_safe(outfptr, (outfptr.HDUposition) + 1, None, status);
    } else if outfptr.Fptr.datastart == DATA_UNDEFINED as LONGLONG {
        ffrdef_safe(outfptr, status); /* rescan header */
    }

    outHduType = outfptr.Fptr.hdutype;

    if *status > 0 {
        return *status;
    }

    if inHduType == IMAGE_HDU || outHduType == IMAGE_HDU {
        ffpmsg_str("Can not copy columns to or from IMAGE HDUs (ffccls)");
        *status = NOT_TABLE;
        return *status;
    }

    if (inHduType == BINARY_TBL && outHduType == ASCII_TBL)
        || (inHduType == ASCII_TBL && outHduType == BINARY_TBL)
    {
        ffpmsg_str("Copying between Binary and ASCII tables is not supported (ffccls)");
        *status = NOT_BTABLE;
        return *status;
    }

    /* Do not allow copying multiple columns in the same HDU because the
    permutations of possible overlapping copies is mind-bending */
    if std::ptr::eq(infptr.Fptr.as_mut(), outfptr.Fptr.as_mut())
        && (infptr.HDUposition == outfptr.HDUposition)
    {
        ffpmsg_str("Copying multiple columns in same HDU is not supported (ffccls)");
        *status = NOT_BTABLE;
        return *status;
    }

    /* Retrieve the number of columns in output file */
    tstatus = 0;
    if ffgkyj_safe(outfptr, cs!(c"TFIELDS"), &mut tfields, None, &mut tstatus) != 0 {
        ffpmsg_str("Could not read TFIELDS keyword in output table (ffccls)");
        *status = NO_TFIELDS;
        return *status;
    }

    colnum = cmp::min((tfields + 1) as c_int, outcol); /* output col. number */

    /* Collect data about input column (type, repeat, etc) */
    let mut incol1 = incol;
    let mut outcol1 = colnum;
    for icol in 0..(ncols as usize) {
        ffgtcl_safe(
            infptr,
            incol1,
            Some(&mut typecode),
            Some(&mut repeat),
            Some(&mut width),
            status,
        );

        if typecode < 0 {
            ffpmsg_str("Variable-length columns are not supported (ffccls)");
            *status = BAD_TFORM;
            return *status;
        }

        typecodes[icol] = typecode;

        tstatus = 0;
        ffkeyn_safe(cs!(c"TTYPE"), incol1, &mut keyname, &mut tstatus);
        ffgkys_safe(
            infptr,
            &keyname,
            &mut ttype,
            Some(&mut ttype_comm),
            &mut tstatus,
        );

        ffkeyn_safe(cs!(c"TFORM"), incol1, &mut keyname, &mut tstatus);

        if ffgkys_safe(
            infptr,
            &keyname,
            &mut tform,
            Some(&mut tform_comm),
            &mut tstatus,
        ) != 0
        {
            ffpmsg_str("Could not find TTYPE and TFORM keywords in input table (ffccls)");
            *status = NO_TFORM;
            return *status;
        }

        /* If creating columns, we need to save these values */
        if create_col != 0 {
            strcpy_safe(&mut tforms[icol], &tform);
            strcpy_safe(&mut ttypes[icol], &ttype);
        } else {
            /* If not creating columns, then check the datatype and vector
            repeat length of the output column */
            ffgtcl_safe(
                outfptr,
                outcol1,
                Some(&mut otypecode),
                Some(&mut orepeat),
                Some(&mut owidth),
                status,
            );

            if orepeat != repeat {
                ffpmsg_str("Input and output vector columns must have same length (ffccls)");
                *status = BAD_TFORM;
                return *status;
            }
        }
        incol1 += 1;
        outcol1 += 1;
    }

    /* Insert columns into output file and copy all meta-data
    keywords, if requested */
    if create_col != 0 {
        /* create the empty columns */
        let tforms_slice = tforms
            .iter()
            .map(|x| x.as_slice())
            .collect::<Vec<&[c_char]>>();
        let ttypes_slice = ttypes
            .iter()
            .map(|x| x.as_slice())
            .collect::<Vec<&[c_char]>>();
        if fficls_safe(outfptr, colnum, ncols, &ttypes_slice, &tforms_slice, status) > 0 {
            ffpmsg_str("Could not append new columns to output file (ffccls)");
            return *status;
        }

        /* Copy meta-data strings from input column to output */
        let mut incol1 = incol;
        let mut outcol1 = colnum;
        for icol in 0..(ncols as usize) {
            /* copy the comment strings from the input file for TTYPE and TFORM */
            ffkeyn_safe(cs!(c"TTYPE"), incol1, &mut keyname, status);
            ffgkys_safe(infptr, &keyname, &mut ttype, Some(&mut ttype_comm), status);
            ffkeyn_safe(cs!(c"TTYPE"), outcol1, &mut keyname, status);
            ffmcom_safe(outfptr, &keyname, Some(&ttype_comm), status);

            ffkeyn_safe(cs!(c"TFORM"), incol1, &mut keyname, status);
            ffgkys_safe(infptr, &keyname, &mut tform, Some(&mut tform_comm), status);
            ffkeyn_safe(cs!(c"TFORM"), outcol1, &mut keyname, status);
            ffmcom_safe(outfptr, &keyname, Some(&tform_comm), status);

            /* copy other column-related keywords if they exist */

            ffcpky_safe(infptr, outfptr, incol1, outcol1, cs!(c"TUNIT"), status);
            ffcpky_safe(infptr, outfptr, incol1, outcol1, cs!(c"TSCAL"), status);
            ffcpky_safe(infptr, outfptr, incol1, outcol1, cs!(c"TZERO"), status);
            ffcpky_safe(infptr, outfptr, incol1, outcol1, cs!(c"TDISP"), status);
            ffcpky_safe(infptr, outfptr, incol1, outcol1, cs!(c"TLMIN"), status);
            ffcpky_safe(infptr, outfptr, incol1, outcol1, cs!(c"TLMAX"), status);
            ffcpky_safe(infptr, outfptr, incol1, outcol1, cs!(c"TDIM"), status);

            /*  WCS keywords */
            ffcpky_safe(infptr, outfptr, incol1, outcol1, cs!(c"TCTYP"), status);
            ffcpky_safe(infptr, outfptr, incol1, outcol1, cs!(c"TCUNI"), status);
            ffcpky_safe(infptr, outfptr, incol1, outcol1, cs!(c"TCRVL"), status);
            ffcpky_safe(infptr, outfptr, incol1, outcol1, cs!(c"TCRPX"), status);
            ffcpky_safe(infptr, outfptr, incol1, outcol1, cs!(c"TCDLT"), status);
            ffcpky_safe(infptr, outfptr, incol1, outcol1, cs!(c"TCROT"), status);

            ffcpky_safe(infptr, outfptr, incol1, outcol1, cs!(c"TNULL"), status);

            incol1 += 1;
            outcol1 += 1;
        }

        /* rescan header to recognize the new keywords */
        if ffrdef_safe(outfptr, status) != 0 {
            return *status;
        }
    }

    /* Copy columns using standard ffcpcl(); do this in a loop because
    the I/O-intensive column expanding is done */
    let mut incol1 = incol;
    let mut outcol1 = colnum;
    for icol in 0..(ncols as usize) {
        ffcpcl_safe(infptr, outfptr, incol1, outcol1, 0, status);
        incol1 += 1;
        outcol1 += 1;
    }

    *status
}

/*--------------------------------------------------------------------------*/
/// copy consecutive set of rows from infptr and append it in the outfptr table.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffcprw(
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

        ffcprw_safe(infptr, outfptr, firstrow, nrows, status)
    }
}

/*--------------------------------------------------------------------------*/
/// copy consecutive set of rows from infptr and append it in the outfptr table.
pub(crate) fn ffcprw_safe(
    infptr: &mut fitsfile,  /* I - FITS file pointer to input file  */
    outfptr: &mut fitsfile, /* I - FITS file pointer to output file */
    firstrow: LONGLONG,     /* I - number of first row to copy (1 based)  */
    nrows: LONGLONG,        /* I - number of rows to copy  */
    status: &mut c_int,     /* IO - error status     */
) -> c_int {
    let mut current_block: LONGLONG;
    let mut innaxis1: LONGLONG = 0;
    let mut innaxis2: LONGLONG = 0;
    let mut outnaxis1: LONGLONG = 0;
    let mut outnaxis2: LONGLONG = 0;
    let ii: LONGLONG = 0;
    let mut jj: LONGLONG = 0;
    let icol: LONGLONG = 0;
    let mut iVarCol: LONGLONG = 0;
    let mut inPos: LONGLONG = 0;
    let mut outPos: LONGLONG = 0;
    let mut nVarBytes: LONGLONG = 0;
    let mut nVarAllocBytes: LONGLONG = 0;
    let mut buffer: *mut u8;
    let mut varColBuff: Vec<u8> = Vec::new();
    let mut nInVarCols: c_int = 0;
    let mut nOutVarCols: c_int = 0;
    let mut varColDiff: c_int = 0;
    let mut inVarCols: *mut c_int;
    let mut outVarCols: *mut c_int;
    let mut nNewBlocks: c_long = 0;
    let mut hrepeat: LONGLONG = 0;
    let mut hoffset: LONGLONG = 0;
    let mut colptr: *mut tcolumn;

    if *status > 0 {
        return *status;
    }

    if infptr.HDUposition != infptr.Fptr.curhdu {
        ffmahd_safe(infptr, (infptr.HDUposition) + 1, None, status);
    } else if infptr.Fptr.datastart == DATA_UNDEFINED as LONGLONG {
        ffrdef_safe(infptr, status); /* rescan header */
    }

    if outfptr.HDUposition != outfptr.Fptr.curhdu {
        ffmahd_safe(outfptr, (outfptr.HDUposition) + 1, None, status);
    } else if outfptr.Fptr.datastart == DATA_UNDEFINED as LONGLONG {
        ffrdef_safe(outfptr, status); /* rescan header */
    }

    if *status > 0 {
        return *status;
    }

    if infptr.Fptr.hdutype == IMAGE_HDU || outfptr.Fptr.hdutype == IMAGE_HDU {
        ffpmsg_str("Can not copy rows to or from IMAGE HDUs (ffcprw)");
        *status = NOT_TABLE;
        return *status;
    }

    if (infptr.Fptr.hdutype == BINARY_TBL && outfptr.Fptr.hdutype == ASCII_TBL)
        || (infptr.Fptr.hdutype == ASCII_TBL && outfptr.Fptr.hdutype == BINARY_TBL)
    {
        ffpmsg_str("Copying rows between Binary and ASCII tables is not supported (ffcprw)");
        *status = NOT_BTABLE;
        return *status;
    }

    ffgkyjj_safe(infptr, cs!(c"NAXIS1"), &mut innaxis1, None, status); /* width of input rows */
    ffgkyjj_safe(infptr, cs!(c"NAXIS2"), &mut innaxis2, None, status); /* no. of input rows */
    ffgkyjj_safe(outfptr, cs!(c"NAXIS1"), &mut outnaxis1, None, status); /* width of output rows */
    ffgkyjj_safe(outfptr, cs!(c"NAXIS2"), &mut outnaxis2, None, status); /* no. of output rows */

    if *status > 0 {
        return *status;
    }

    if outnaxis1 != innaxis1 {
        ffpmsg_str("Input and output tables do not have same width (ffcprw)");
        *status = BAD_ROW_WIDTH;
        return *status;
    }

    if firstrow + nrows - 1 > innaxis2 {
        ffpmsg_str("Not enough rows in input table to copy (ffcprw)");
        *status = BAD_ROW_NUM;
        return *status;
    }

    if infptr.Fptr.tfield != outfptr.Fptr.tfield {
        ffpmsg_str("Input and output tables do not have same number of columns (ffcprw)");
        *status = BAD_COL_NUM;
        return *status;
    }

    let mut buffer: Vec<u8> = Vec::new();

    /* allocate buffer to hold 1 row of data */
    if buffer.try_reserve_exact(innaxis1 as usize).is_err() {
        ffpmsg_str("Unable to allocate memory (ffcprw)");
        *status = MEMORY_ALLOCATION;
        return *status;
    } else {
        buffer.resize(innaxis1 as usize, 0);
    }

    let mut inVarCols: Vec<c_int> = vec![0; infptr.Fptr.tfield as usize];
    let mut outVarCols: Vec<c_int> = vec![0; outfptr.Fptr.tfield as usize];

    fffvcl_safe(
        infptr,
        &mut nInVarCols,
        Some(inVarCols.as_mut_slice()),
        status,
    );
    fffvcl_safe(
        outfptr,
        &mut nOutVarCols,
        Some(outVarCols.as_mut_slice()),
        status,
    );
    if nInVarCols != nOutVarCols {
        varColDiff = 1;
    } else {
        for ii in 0..(nInVarCols as usize) {
            if inVarCols[ii] != outVarCols[ii] {
                varColDiff = 1;
                break;
            }
        }
    }

    if varColDiff != 0 {
        ffpmsg_str("Input and output tables have different variable columns (ffcprw)");
        *status = BAD_COL_NUM;
        return *status;
    }

    jj = outnaxis2 + 1;
    if nInVarCols != 0 {
        ffirow_safe(outfptr, outnaxis2, nrows, status);
        for ii in (firstrow)..(firstrow + nrows) {
            ffgtbb_safe(infptr, ii, 1, innaxis1, &mut buffer, status);
            ffptbb_safe(outfptr, jj, 1, innaxis1, &buffer, status);
            /* Now make corrections for variable length columns */
            iVarCol = 0;

            let colptr = infptr.Fptr.tableptr; /* point to first column structure */
            let c = infptr.Fptr.get_tableptr_as_slice();

            for icol in 0..(infptr.Fptr.tfield as usize) {
                if iVarCol < nInVarCols as LONGLONG
                    && inVarCols[iVarCol as usize] == (icol + 1) as c_int
                {
                    /* Copy from a variable length column */

                    ffgdesll_safe(
                        infptr,
                        (icol + 1) as c_int,
                        ii,
                        Some(&mut hrepeat),
                        Some(&mut hoffset),
                        status,
                    );
                    /* If this is a bit column, hrepeat will be number of
                            bits, not bytes. If it is a string column, hrepeat
                    is the number of bytes, twidth is the max col width
                    and can be ignored.*/
                    let c = infptr.Fptr.get_tableptr_as_slice();
                    if c[icol].tdatatype == -TBIT {
                        nVarBytes = (hrepeat + 7) / 8;
                    } else if c[icol].tdatatype == -TSTRING {
                        nVarBytes = hrepeat;
                    } else {
                        nVarBytes = hrepeat
                            * c[icol].twidth as LONGLONG
                            * mem::size_of::<c_char>() as LONGLONG;
                    }
                    inPos = infptr.Fptr.datastart + infptr.Fptr.heapstart + hoffset;
                    outPos =
                        outfptr.Fptr.datastart + outfptr.Fptr.heapstart + outfptr.Fptr.heapsize;
                    ffmbyt_safe(infptr, inPos, REPORT_EOF, status);
                    /* If this is not the last HDU in the file, then check if */
                    /* extending the heap would overwrite the following header. */
                    /* If so, then have to insert more blocks. */
                    if (outfptr.Fptr.lasthdu) == 0 {
                        let out_headstart = outfptr.Fptr.get_headstart_as_slice();
                        if outPos + nVarBytes > out_headstart[(outfptr.Fptr.curhdu + 1) as usize] {
                            nNewBlocks = (((outPos + nVarBytes
                                - 1
                                - out_headstart[(outfptr.Fptr.curhdu + 1) as usize])
                                / 2880)
                                + 1) as c_long;
                            if ffiblk(outfptr, nNewBlocks, 1, status) > 0 {
                                ffpmsg_str(
                                    "Failed to extend the size of the variable length heap (ffcprw)",
                                );
                                return *status;
                            }
                        }
                    }

                    if nVarBytes != 0 {
                        if nVarBytes > nVarAllocBytes {
                            /* Grow the copy buffer to accomodate the new maximum size.
                            Note it is safe to call realloc() with null input pointer,
                            which is equivalent to malloc(). */

                            let additional = nVarBytes as usize - varColBuff.capacity();
                            if varColBuff.try_reserve_exact(additional).is_err() {
                                *status = MEMORY_ALLOCATION;
                                ffpmsg_str(
                                    "failed to allocate memory for variable column copy (ffcprw)",
                                );
                                return *status;
                            } else {
                                /* Record the new state */
                                nVarAllocBytes = nVarBytes;
                                varColBuff.resize(nVarBytes as usize, 0);
                            };
                        }
                        /* Copy date from input to output */
                        ffgbyt(infptr, nVarBytes, varColBuff.as_mut_slice(), status);
                        ffmbyt_safe(outfptr, outPos, IGNORE_EOF, status);
                        ffpbyt(outfptr, nVarBytes, varColBuff.as_mut_slice(), status);
                    }
                    ffpdes_safe(
                        outfptr,
                        (icol + 1) as c_int,
                        jj,
                        hrepeat,
                        outfptr.Fptr.heapsize,
                        status,
                    );
                    outfptr.Fptr.heapsize += nVarBytes;
                    iVarCol += 1;
                }
            }
            jj += 1;
        }
    } else {
        /* copy the rows, 1 at a time */
        for ii in firstrow..(firstrow + nrows) {
            ffgtbb_safe(infptr, ii, 1, innaxis1, &mut buffer, status);
            ffptbb_safe(outfptr, jj, 1, innaxis1, &buffer, status);
            jj += 1;
        }
    }
    outnaxis2 += nrows;
    ffuky_safe(
        outfptr,
        KeywordDatatype::TLONGLONG(&outnaxis2),
        cs!(c"NAXIS2"),
        None,
        status,
    );

    *status
}

/*--------------------------------------------------------------------------*/
/// copy consecutive set of rows from infptr and append it in the outfptr table.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffcpsr(
    infptr: *mut fitsfile,     /* I - FITS file pointer to input file  */
    outfptr: *mut fitsfile,    /* I - FITS file pointer to output file */
    firstrow: LONGLONG,        /* I - number of first row to copy (1 based)  */
    nrows: LONGLONG,           /* I - number of rows to copy  */
    row_status: *const c_char, /* I - quality list of rows to keep (1) or not keep (0) */
    status: *mut c_int,        /* IO - error status     */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let infptr = infptr.as_mut().expect(NULL_MSG);
        let outfptr = outfptr.as_mut().expect(NULL_MSG);

        let row_status = match row_status.is_null() {
            false => Some(slice::from_raw_parts(row_status, nrows as usize)),
            true => None,
        };

        ffcpsr_safe(infptr, outfptr, firstrow, nrows, row_status, status)
    }
}

/*--------------------------------------------------------------------------*/
/// copy consecutive set of rows from infptr and append it in the outfptr table.
pub(crate) fn ffcpsr_safe(
    infptr: &mut fitsfile,         /* I - FITS file pointer to input file  */
    outfptr: &mut fitsfile,        /* I - FITS file pointer to output file */
    firstrow: LONGLONG,            /* I - number of first row to copy (1 based)  */
    nrows: LONGLONG,               /* I - number of rows to copy  */
    row_status: Option<&[c_char]>, /* I - quality list of rows to keep (1) or not keep (0) */
    status: &mut c_int,            /* IO - error status     */
) -> c_int {
    let mut current_block: LONGLONG;
    let mut innaxis1: LONGLONG = 0;
    let mut innaxis2: LONGLONG = 0;
    let mut outnaxis1: LONGLONG = 0;
    let mut outnaxis2: LONGLONG = 0;
    let ii: LONGLONG = 0;
    let mut jj: LONGLONG = 0;
    let i0: LONGLONG = 0;
    let icol: LONGLONG = 0;
    let mut iVarCol: LONGLONG = 0;
    let mut inPos: LONGLONG = 0;
    let mut outPos: LONGLONG = 0;
    let mut nVarBytes: LONGLONG = 0;
    let mut nVarAllocBytes: LONGLONG = 0;

    let mut varColBuff: Vec<u8> = Vec::new();
    let mut nInVarCols: c_int = 0;
    let mut nOutVarCols: c_int = 0;
    let mut varColDiff: c_int = 0;
    let mut inVarCols: *mut c_int;
    let mut outVarCols: *mut c_int;
    let mut nNewBlocks: c_long = 0;
    let mut hrepeat: LONGLONG = 0;
    let mut hoffset: LONGLONG = 0;
    let mut colptr: *mut tcolumn;
    let mut n_good_rows: LONGLONG = nrows;

    if *status > 0 {
        return *status;
    }

    if infptr.HDUposition != infptr.Fptr.curhdu {
        ffmahd_safe(infptr, (infptr.HDUposition) + 1, None, status);
    } else if infptr.Fptr.datastart == DATA_UNDEFINED as LONGLONG {
        ffrdef_safe(infptr, status); /* rescan header */
    }

    if outfptr.HDUposition != outfptr.Fptr.curhdu {
        ffmahd_safe(outfptr, (outfptr.HDUposition) + 1, None, status);
    } else if outfptr.Fptr.datastart == DATA_UNDEFINED as LONGLONG {
        ffrdef_safe(outfptr, status); /* rescan header */
    }

    if *status > 0 {
        return *status;
    }

    if infptr.Fptr.hdutype == IMAGE_HDU || outfptr.Fptr.hdutype == IMAGE_HDU {
        ffpmsg_str("Can not copy rows to or from IMAGE HDUs (ffcprw)");
        *status = NOT_TABLE;
        return *status;
    }

    if (infptr.Fptr.hdutype == BINARY_TBL && outfptr.Fptr.hdutype == ASCII_TBL)
        || (infptr.Fptr.hdutype == ASCII_TBL && outfptr.Fptr.hdutype == BINARY_TBL)
    {
        ffpmsg_str("Copying rows between Binary and ASCII tables is not supported (ffcprw)");
        *status = NOT_BTABLE;
        return *status;
    }

    ffgkyjj_safe(infptr, cs!(c"NAXIS1"), &mut innaxis1, None, status); /* width of input rows */
    ffgkyjj_safe(infptr, cs!(c"NAXIS2"), &mut innaxis2, None, status); /* no. of input rows */
    ffgkyjj_safe(outfptr, cs!(c"NAXIS1"), &mut outnaxis1, None, status); /* width of output rows */
    ffgkyjj_safe(outfptr, cs!(c"NAXIS2"), &mut outnaxis2, None, status); /* no. of output rows */

    if *status > 0 {
        return *status;
    }

    if outnaxis1 != innaxis1 {
        ffpmsg_str("Input and output tables do not have same width (ffcprw)");
        *status = BAD_ROW_WIDTH;
        return *status;
    }

    if firstrow + nrows - 1 > innaxis2 {
        ffpmsg_str("Not enough rows in input table to copy (ffcprw)");
        *status = BAD_ROW_NUM;
        return *status;
    }

    if infptr.Fptr.tfield != outfptr.Fptr.tfield {
        ffpmsg_str("Input and output tables do not have same number of columns (ffcprw)");
        *status = BAD_COL_NUM;
        return *status;
    }

    /* allocate buffer to hold 1 row of data */
    let mut buffer: Vec<u8> = Vec::new();

    if buffer.try_reserve_exact(innaxis1 as usize).is_err() {
        ffpmsg_str("Unable to allocate memory (ffcprw)");
        *status = MEMORY_ALLOCATION;
        return *status;
    } else {
        buffer.resize(innaxis1 as usize, 0);
    }

    let mut inVarCols: Vec<c_int> = vec![0; infptr.Fptr.tfield as usize];
    let mut outVarCols: Vec<c_int> = vec![0; outfptr.Fptr.tfield as usize];
    fffvcl_safe(
        infptr,
        &mut nInVarCols,
        Some(inVarCols.as_mut_slice()),
        status,
    );
    fffvcl_safe(
        outfptr,
        &mut nOutVarCols,
        Some(outVarCols.as_mut_slice()),
        status,
    );
    if nInVarCols != nOutVarCols {
        varColDiff = 1;
    } else {
        for ii in 0..(nInVarCols as usize) {
            if inVarCols[ii] != outVarCols[ii] {
                varColDiff = 1;
                break;
            }
        }
    }

    if varColDiff != 0 {
        ffpmsg_str("Input and output tables have different variable columns (ffcprw)");
        *status = BAD_COL_NUM;
        return *status;
    }

    jj = outnaxis2 + 1;
    if nInVarCols != 0 {
        if let Some(row_status) = row_status {
            let mut n_good_rows = 0;
            for ii in 0..(nrows as usize) {
                if row_status[ii] != 0 {
                    n_good_rows += 1;
                }
            }
        }

        ffirow_safe(outfptr, outnaxis2, n_good_rows, status);
        let mut ii = firstrow;
        for i0 in 0..(nrows as usize) {
            /* Ignore rows with row_status[] == 0 */
            if let Some(row_status) = row_status {
                if row_status[i0] == 0 {
                    continue;
                }
            }

            ffgtbb_safe(infptr, ii, 1, innaxis1, &mut buffer, status);
            ffptbb_safe(outfptr, jj, 1, innaxis1, &buffer, status);
            /* Now make corrections for variable length columns */
            iVarCol = 0;
            let colptr = infptr.Fptr.tableptr; /* point to first column structure */
            let c = infptr.Fptr.get_tableptr_as_slice();

            for icol in 0..(infptr.Fptr.tfield as usize) {
                if iVarCol < (nInVarCols as LONGLONG)
                    && inVarCols[iVarCol as usize] == (icol + 1) as c_int
                {
                    /* Copy from a variable length column */

                    ffgdesll_safe(
                        infptr,
                        (icol + 1) as c_int,
                        ii,
                        Some(&mut hrepeat),
                        Some(&mut hoffset),
                        status,
                    );
                    /*
                    If this is a bit column, hrepeat will be number of bits, not bytes.
                    If it is a string column, hrepeat is the number of bytes, twidth is the max col width
                    and can be ignored.
                    */

                    let c = infptr.Fptr.get_tableptr_as_slice();
                    if c[icol].tdatatype == -TBIT {
                        nVarBytes = (hrepeat + 7) / 8;
                    } else if c[icol].tdatatype == -TSTRING {
                        nVarBytes = hrepeat;
                    } else {
                        nVarBytes = hrepeat
                            * c[icol].twidth as LONGLONG
                            * mem::size_of::<c_char>() as LONGLONG;
                    }
                    inPos = infptr.Fptr.datastart + infptr.Fptr.heapstart + hoffset;
                    outPos =
                        outfptr.Fptr.datastart + outfptr.Fptr.heapstart + outfptr.Fptr.heapsize;
                    ffmbyt_safe(infptr, inPos, REPORT_EOF, status);
                    /* If this is not the last HDU in the file, then check if */
                    /* extending the heap would overwrite the following header. */
                    /* If so, then have to insert more blocks. */
                    if (outfptr.Fptr.lasthdu) == 0 {
                        let out_headstart = outfptr.Fptr.get_headstart_as_slice();
                        if outPos + nVarBytes > out_headstart[(outfptr.Fptr.curhdu + 1) as usize] {
                            nNewBlocks = (((outPos + nVarBytes
                                - 1
                                - out_headstart[(outfptr.Fptr.curhdu + 1) as usize])
                                / BL!())
                                + 1) as c_long;
                            if ffiblk(outfptr, nNewBlocks, 1, status) > 0 {
                                ffpmsg_str(
                                    "Failed to extend the size of the variable length heap (ffcprw)",
                                );
                                return *status;
                            }
                        }
                    }
                    if nVarBytes != 0 {
                        if nVarBytes > nVarAllocBytes {
                            /* Grow the copy buffer to accomodate the new maximum size.
                            Note it is safe to call realloc() with null input pointer,
                            which is equivalent to malloc(). */

                            let additional = nVarBytes as usize - varColBuff.capacity();
                            if varColBuff.try_reserve_exact(additional).is_err() {
                                *status = MEMORY_ALLOCATION;
                                ffpmsg_str(
                                    "failed to allocate memory for variable column copy (ffcprw)",
                                );
                                return *status;
                            } else {
                                /* Record the new state */
                                nVarAllocBytes = nVarBytes;
                                varColBuff.resize(nVarBytes as usize, 0);
                            };
                        }
                        /* Copy date from input to output */
                        ffgbyt(infptr, nVarBytes, varColBuff.as_mut_slice(), status);
                        ffmbyt_safe(outfptr, outPos, IGNORE_EOF, status);
                        ffpbyt(outfptr, nVarBytes, varColBuff.as_mut_slice(), status);
                    }
                    ffpdes_safe(
                        outfptr,
                        (icol + 1) as c_int,
                        jj,
                        hrepeat,
                        outfptr.Fptr.heapsize,
                        status,
                    );
                    outfptr.Fptr.heapsize += nVarBytes;
                    iVarCol += 1;
                }
            }
            jj += 1;
            ii += 1;
        }
    } else {
        /* copy the rows, 1 at a time */
        n_good_rows = 0;
        let mut ii = firstrow;
        for i0 in 0..(nrows as usize) {
            /* Ignore rows with row_status[] == 0 */
            if let Some(row_status) = row_status {
                if row_status[i0] == 0 {
                    continue;
                }
            }

            ffgtbb_safe(infptr, ii, 1, innaxis1, &mut buffer, status);
            ffptbb_safe(outfptr, jj, 1, innaxis1, &buffer, status);
            n_good_rows += 1;
            jj += 1;
            ii += 1;
        }
    }
    outnaxis2 += n_good_rows;
    ffuky_safe(
        outfptr,
        KeywordDatatype::TLONGLONG(&outnaxis2),
        cs!(c"NAXIS2"),
        None,
        status,
    );

    *status
}

/*--------------------------------------------------------------------------*/
/// copy an indexed keyword from infptr to outfptr.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffcpky(
    infptr: *mut fitsfile,   /* I - FITS file pointer to input file  */
    outfptr: *mut fitsfile,  /* I - FITS file pointer to output file */
    incol: c_int,            /* I - input index number   */
    outcol: c_int,           /* I - output index number  */
    rootname: *const c_char, /* I - root name of the keyword to be copied */
    status: *mut c_int,      /* IO - error status     */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let infptr = infptr.as_mut().expect(NULL_MSG);
        let outfptr = outfptr.as_mut().expect(NULL_MSG);

        raw_to_slice!(rootname);

        ffcpky_safe(infptr, outfptr, incol, outcol, rootname, status)
    }
}

/*--------------------------------------------------------------------------*/
/// copy an indexed keyword from infptr to outfptr.
pub(crate) fn ffcpky_safe(
    infptr: &mut fitsfile,  /* I - FITS file pointer to input file  */
    outfptr: &mut fitsfile, /* I - FITS file pointer to output file */
    incol: c_int,           /* I - input index number   */
    outcol: c_int,          /* I - output index number  */
    rootname: &[c_char],    /* I - root name of the keyword to be copied */
    status: &mut c_int,     /* IO - error status     */
) -> c_int {
    let mut tstatus: c_int = 0;

    let mut keyname: [c_char; FLEN_KEYWORD] = [0; FLEN_KEYWORD];
    let mut value: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut comment: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];
    let mut card: [c_char; FLEN_CARD] = [0; FLEN_CARD];

    ffkeyn_safe(rootname, incol, &mut keyname, &mut tstatus);
    if ffgkey_safe(
        infptr,
        &keyname,
        &mut value,
        Some(&mut comment),
        &mut tstatus,
    ) <= 0
    {
        ffkeyn_safe(rootname, outcol, &mut keyname, &mut tstatus);
        ffmkky_safe(&keyname, &value, Some(&comment), &mut card, status);
        ffprec_safe(outfptr, &card, status);
    }
    *status
}

/*--------------------------------------------------------------------------*/
/// Delete a column from a table.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffdcol(
    fptr: *mut fitsfile, /* I - FITS file pointer                        */
    colnum: c_int,       /* I - column to delete (1 = 1st)               */
    status: *mut c_int,  /* IO - error status                            */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        ffdcol_safe(fptr, colnum, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Delete a column from a table.
pub(crate) fn ffdcol_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer                        */
    colnum: c_int,       /* I - column to delete (1 = 1st)               */
    status: &mut c_int,  /* IO - error status                            */
) -> c_int {
    let ii: c_int = 0;
    let mut tstatus: c_int = 0;
    let mut firstbyte: LONGLONG = 0;
    let mut size: LONGLONG = 0;
    let mut ndelete: LONGLONG = 0;
    let mut nbytes: LONGLONG = 0;
    let mut naxis1: LONGLONG = 0;
    let mut naxis2: LONGLONG = 0;
    let mut firstcol: LONGLONG = 0;
    let mut delbyte: LONGLONG = 0;
    let mut freespace: LONGLONG = 0;
    let mut tbcol: LONGLONG = 0;
    let mut nblock: c_long = 0;
    let mut nspace: c_long = 0;
    let mut keyname: [c_char; FLEN_KEYWORD] = [0; FLEN_KEYWORD];
    let mut comm: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];
    let mut colptr: *mut tcolumn;
    let mut nextcol: usize = 0;

    if *status > 0 {
        return *status;
    }

    if fptr.HDUposition != fptr.Fptr.curhdu {
        ffmahd_safe(fptr, (fptr.HDUposition) + 1, None, status);
    }
    /* rescan header if data structure is undefined */
    else if fptr.Fptr.datastart == DATA_UNDEFINED as LONGLONG && ffrdef_safe(fptr, status) > 0 {
        return *status;
    }

    if fptr.Fptr.hdutype == IMAGE_HDU {
        ffpmsg_str("Can only delete column from TABLE or BINTABLE extension (ffdcol)");
        *status = NOT_TABLE;
        return *status;
    }

    if colnum < 1 || colnum > fptr.Fptr.tfield {
        *status = BAD_COL_NUM;
        return *status;
    }

    let colptr = fptr.Fptr.tableptr; /* point to first column structure */
    let c = fptr.Fptr.get_tableptr_as_slice();
    let ci = (colnum - 1) as usize; /* offset to the correct column */

    firstcol = c[ci].tbcol; /* starting byte position of the column */

    /* use column width to determine how many bytes to delete in each row */
    if fptr.Fptr.hdutype == ASCII_TBL {
        delbyte = c[ci].twidth as LONGLONG; /* width of ASCII column */

        if colnum < fptr.Fptr.tfield
        /* check for space between next column */
        {
            nextcol = ci + 1;
            nspace = ((c[nextcol].tbcol) - (c[ci].tbcol) - delbyte) as c_long;
            if nspace > 0 {
                delbyte += 1;
            }
        } else if colnum > 1
        /* check for space between last 2 columns */
        {
            nextcol = ci - 1;
            nspace =
                ((c[ci].tbcol) - (c[nextcol].tbcol) - (c[nextcol].twidth as LONGLONG)) as c_long;
            if nspace > 0 {
                delbyte += 1;
                firstcol -= 1; /* delete the leading space */
            }
        }
    } else
    /* a binary table */
    if colnum < fptr.Fptr.tfield {
        nextcol = ci + 1;
        delbyte = (c[nextcol].tbcol) - (c[ci].tbcol);
    } else {
        delbyte = (fptr.Fptr.rowlength) - (c[ci].tbcol);
    }

    naxis1 = fptr.Fptr.rowlength; /* current width of the table */
    naxis2 = fptr.Fptr.numrows;

    /* current size of table */
    size = fptr.Fptr.heapstart + fptr.Fptr.heapsize;
    freespace = ((delbyte as LONGLONG) * naxis2) + ((size + 2879) / 2880) * 2880 - size;
    nblock = (freespace / 2880) as c_long; /* number of empty blocks to delete */

    ffcdel_safe(fptr, naxis1, naxis2, delbyte, firstcol, status); /* delete col */

    /* absolute heap position */
    firstbyte = fptr.Fptr.datastart + fptr.Fptr.heapstart;
    ndelete = (delbyte as LONGLONG) * naxis2; /* size of shift */

    /* shift heap up (if it exists) */
    if fptr.Fptr.heapsize > 0 {
        nbytes = fptr.Fptr.heapsize; /* no. of bytes to shift up */

        if ffshft_safe(fptr, firstbyte, nbytes, -ndelete, status) > 0 {
            /* mv heap */
            return *status;
        }
    }

    /* delete the empty  blocks at the end of the HDU */
    if nblock > 0 {
        ffdblk(fptr, nblock, status);
    }

    /* update the heap starting address */
    fptr.Fptr.heapstart -= ndelete;

    /* update the THEAP keyword if it exists */
    tstatus = 0;
    ffmkyj_safe(
        fptr,
        cs!(c"THEAP"),
        fptr.Fptr.heapstart as LONGLONG,
        Some(cs!(c"&")),
        &mut tstatus,
    );

    if fptr.Fptr.hdutype == ASCII_TBL {
        /* adjust the TBCOL values of the remaining columns */
        for ii in 1..=(fptr.Fptr.tfield) {
            ffkeyn_safe(cs!(c"TBCOL"), ii, &mut keyname, status);
            ffgkyjj_safe(fptr, &keyname, &mut tbcol, Some(&mut comm), status);
            if tbcol > firstcol {
                tbcol -= delbyte;
                ffmkyj_safe(fptr, &keyname, tbcol, Some(cs!(c"&")), status);
            }
        }
    }

    /* update the mandatory keywords */
    ffmkyj_safe(
        fptr,
        cs!(c"TFIELDS"),
        ((fptr.Fptr.tfield) - 1) as LONGLONG,
        Some(cs!(c"&")),
        status,
    );
    ffmkyj_safe(
        fptr,
        cs!(c"NAXIS1"),
        naxis1 - delbyte,
        Some(cs!(c"&")),
        status,
    );
    /*
    delete the index keywords starting with 'T' associated with the
    deleted column and subtract 1 from index of all higher keywords
    */
    ffkshf_safe(fptr, colnum, fptr.Fptr.tfield, -1, status);

    ffrdef_safe(fptr, status); /* initialize the new table structure */
    *status
}

/*--------------------------------------------------------------------------*/
/// Insert 'ninsert' bytes into each row of the table at position 'bytepos'.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffcins(
    fptr: *mut fitsfile, /* I - FITS file pointer                        */
    naxis1: LONGLONG,    /* I - width of the table, in bytes             */
    naxis2: LONGLONG,    /* I - number of rows in the table              */
    ninsert: LONGLONG,   /* I - number of bytes to insert in each row    */
    bytepos: LONGLONG,   /* I - rel. position in row to insert bytes     */
    status: *mut c_int,  /* IO - error status                            */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        ffcins_safe(fptr, naxis1, naxis2, ninsert, bytepos, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Insert 'ninsert' bytes into each row of the table at position 'bytepos'.
pub(crate) fn ffcins_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer                        */
    naxis1: LONGLONG,    /* I - width of the table, in bytes             */
    naxis2: LONGLONG,    /* I - number of rows in the table              */
    ninsert: LONGLONG,   /* I - number of bytes to insert in each row    */
    bytepos: LONGLONG,   /* I - rel. position in row to insert bytes     */
    status: &mut c_int,  /* IO - error status                            */
) -> c_int {
    let mut buffer: [u8; 10000] = [0; 10000];
    let mut cfill: u8 = 0;
    let mut newlen: LONGLONG = 0;
    let mut fbyte: LONGLONG = 0;
    let mut nbytes: LONGLONG = 0;
    let irow: LONGLONG = 0;
    let mut nseg: LONGLONG = 0;
    let ii: LONGLONG = 0;

    if *status > 0 {
        return *status;
    }

    if naxis2 == 0 {
        return *status; /* just return if there are 0 rows in the table */
    }

    /* select appropriate fill value */
    if fptr.Fptr.hdutype == ASCII_TBL {
        cfill = 32; /* ASCII tables use blank fill */
    } else {
        cfill = 0; /* primary array and binary tables use zero fill */
    }

    newlen = naxis1 + ninsert;

    if newlen <= 10000 {
        /*******************************************************************
        CASE #1: optimal case where whole new row fits in the work buffer
        *******************************************************************/
        for ii in 0..(ninsert as usize) {
            buffer[ii] = cfill; /* initialize buffer with fill value */
        }

        /* first move the trailing bytes (if any) in the last row */
        fbyte = bytepos + 1;
        nbytes = naxis1 - bytepos;
        /* If the last row hasn't yet been accessed in full, it's possible
           that logfilesize hasn't been updated to account for it (by way
           of an ffldrc call).  This could cause ffgtbb to return with an
           EOF error.  To prevent this, we must increase logfilesize here.
        */
        if fptr.Fptr.logfilesize < fptr.Fptr.datastart + fptr.Fptr.heapstart {
            fptr.Fptr.logfilesize =
                ((fptr.Fptr.datastart + fptr.Fptr.heapstart + (BL!() - 1)) / BL!()) * BL!();
        }

        ffgtbb_safe(
            fptr,
            naxis2,
            fbyte,
            nbytes,
            &mut buffer[(ninsert as usize)..],
            status,
        );
        fptr.Fptr.rowlength = newlen; /*  new row length */

        /* write the row (with leading fill bytes) in the new place */
        nbytes += ninsert;
        ffptbb_safe(fptr, naxis2, fbyte, nbytes, &buffer, status);
        fptr.Fptr.rowlength = naxis1; /* reset to orig. value */

        /*  now move the rest of the rows */
        for irow in (1..=(naxis2 as usize - 1)).rev() {
            /* read the row to be shifted (work backwards thru the table) */
            ffgtbb_safe(
                fptr,
                irow as LONGLONG,
                fbyte,
                naxis1,
                &mut buffer[(ninsert as usize)..],
                status,
            );
            fptr.Fptr.rowlength = newlen; /* new row length */

            /* write the row (with the leading fill bytes) in the new place */
            ffptbb_safe(fptr, irow as LONGLONG, fbyte, newlen, &buffer, status);
            fptr.Fptr.rowlength = naxis1; /* reset to orig value */
        }
    } else {
        /*****************************************************************
        CASE #2:  whole row doesn't fit in work buffer; move row in pieces
        ******************************************************************
        first copy the data, then go back and write fill into the new column
        start by copying the trailing bytes (if any) in the last row.     */

        nbytes = naxis1 - bytepos;
        nseg = (nbytes + 9999) / 10000;
        fbyte = (nseg - 1) * 10000 + bytepos + 1;
        nbytes = naxis1 - fbyte + 1;

        for ii in 0..(nseg as usize) {
            ffgtbb_safe(fptr, naxis2, fbyte, nbytes, &mut buffer, status);
            fptr.Fptr.rowlength = newlen; /* new row length */

            ffptbb_safe(fptr, naxis2, fbyte + ninsert, nbytes, &buffer, status);
            fptr.Fptr.rowlength = naxis1; /* reset to orig value */

            fbyte -= 10000;
            nbytes = 10000;
        }

        /* now move the rest of the rows */
        nseg = (naxis1 + 9999) / 10000;
        for irow in (1..=(naxis2 as usize - 1)).rev() {
            //for (irow = naxis2 - 1; irow > 0; irow--)

            fbyte = (nseg - 1) * 10000 + bytepos + 1;
            nbytes = naxis1 - (nseg - 1) * 10000;
            for ii in 0..(nseg as usize) {
                /* read the row to be shifted (work backwards thru the table) */
                ffgtbb_safe(fptr, irow as LONGLONG, fbyte, nbytes, &mut buffer, status);
                fptr.Fptr.rowlength = newlen; /* new row length */

                /* write the row in the new place */
                ffptbb_safe(
                    fptr,
                    irow as LONGLONG,
                    fbyte + ninsert,
                    nbytes,
                    &buffer,
                    status,
                );
                fptr.Fptr.rowlength = naxis1; /* reset to orig value */

                fbyte -= 10000;
                nbytes = 10000;
            }
        }

        /* now write the fill values into the new column */
        nbytes = cmp::min(ninsert, 10000);
        buffer.fill(cfill); /* initialize with fill value */

        nseg = (ninsert + 9999) / 10000;
        fptr.Fptr.rowlength = newlen; /* new row length */

        for irow in 1..=(naxis2 as usize) {
            fbyte = bytepos + 1;
            nbytes = ninsert - ((nseg - 1) * 10000);
            for ii in 0..(nseg as usize) {
                ffptbb_safe(fptr, irow as LONGLONG, fbyte, nbytes, &buffer, status);
                fbyte += nbytes;
                nbytes = 10000;
            }
        }
        fptr.Fptr.rowlength = naxis1; /* reset to orig value */
    }
    *status
}

/*--------------------------------------------------------------------------*/
/// delete 'ndelete' bytes from each row of the table at position 'bytepos'.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffcdel(
    fptr: *mut fitsfile, /* I - FITS file pointer                        */
    naxis1: LONGLONG,    /* I - width of the table, in bytes             */
    naxis2: LONGLONG,    /* I - number of rows in the table              */
    ndelete: LONGLONG,   /* I - number of bytes to delete in each row    */
    bytepos: LONGLONG,   /* I - rel. position in row to delete bytes     */
    status: *mut c_int,  /* IO - error status                            */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        ffcdel_safe(fptr, naxis1, naxis2, ndelete, bytepos, status)
    }
}

/*--------------------------------------------------------------------------*/
/// delete 'ndelete' bytes from each row of the table at position 'bytepos'.
pub(crate) fn ffcdel_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer                        */
    naxis1: LONGLONG,    /* I - width of the table, in bytes             */
    naxis2: LONGLONG,    /* I - number of rows in the table              */
    ndelete: LONGLONG,   /* I - number of bytes to delete in each row    */
    bytepos: LONGLONG,   /* I - rel. position in row to delete bytes     */
    status: &mut c_int,  /* IO - error status                            */
) -> c_int {
    let mut buffer: [u8; 10000] = [0; 10000];
    let mut i1: LONGLONG = 0;
    let mut i2: LONGLONG = 0;
    let ii: LONGLONG = 0;
    let irow: LONGLONG = 0;
    let mut nseg: LONGLONG = 0;
    let mut newlen: LONGLONG = 0;
    let mut remain: LONGLONG = 0;
    let mut nbytes: LONGLONG = 0;

    if *status > 0 {
        return *status;
    }

    if naxis2 == 0 {
        return *status; /* just return if there are 0 rows in the table */
    }

    newlen = naxis1 - ndelete;

    if newlen <= 10000 {
        /*******************************************************************
        CASE #1: optimal case where whole new row fits in the work buffer
        *******************************************************************/
        i1 = bytepos + 1;
        i2 = i1 + ndelete;
        for irow in 1..(naxis2 as usize) {
            ffgtbb_safe(fptr, irow as LONGLONG, i2, newlen, &mut buffer, status); /* read row */
            fptr.Fptr.rowlength = newlen; /* new row length */

            ffptbb_safe(fptr, irow as LONGLONG, i1, newlen, &buffer, status); /* write row */
            fptr.Fptr.rowlength = naxis1; /* reset to orig value */
        }

        /* now do the last row */
        remain = naxis1 - (bytepos + ndelete);

        if remain > 0 {
            ffgtbb_safe(fptr, naxis2, i2, remain, &mut buffer, status); /* read row */
            fptr.Fptr.rowlength = newlen; /* new row length */

            ffptbb_safe(fptr, naxis2, i1, remain, &buffer, status); /* write row */
            fptr.Fptr.rowlength = naxis1; /* reset to orig value */
        }
    } else {
        /*****************************************************************
        CASE #2:  whole row doesn't fit in work buffer; move row in pieces
        ******************************************************************/

        nseg = (newlen + 9999) / 10000;
        for irow in 1..(naxis2 as usize) {
            i1 = bytepos + 1;
            i2 = i1 + ndelete;

            nbytes = newlen - (nseg - 1) * 10000;
            for ii in 0..(nseg as usize) {
                ffgtbb_safe(fptr, irow as LONGLONG, i2, nbytes, &mut buffer, status); /* read bytes */
                fptr.Fptr.rowlength = newlen; /* new row length */

                ffptbb_safe(fptr, irow as LONGLONG, i1, nbytes, &buffer, status); /* rewrite bytes */
                fptr.Fptr.rowlength = naxis1; /* reset to orig value */

                i1 += nbytes;
                i2 += nbytes;
                nbytes = 10000;
            }
        }

        /* now do the last row */
        remain = naxis1 - (bytepos + ndelete);

        if remain > 0 {
            nseg = (remain + 9999) / 10000;
            i1 = bytepos + 1;
            i2 = i1 + ndelete;
            nbytes = remain - (nseg - 1) * 10000;
            for ii in 0..(nseg as usize) {
                ffgtbb_safe(fptr, naxis2, i2, nbytes, &mut buffer, status);
                fptr.Fptr.rowlength = newlen; /* new row length */

                ffptbb_safe(fptr, naxis2, i1, nbytes, &buffer, status); /* write row */
                fptr.Fptr.rowlength = naxis1; /* reset to orig value */

                i1 += nbytes;
                i2 += nbytes;
                nbytes = 10000;
            }
        }
    }
    *status
}

/*--------------------------------------------------------------------------*/
/// shift the index value on any existing column keywords
/// This routine will modify the name of any keyword that begins with 'T'
/// and has an index number in the range COLMIN - COLMAX, inclusive.
///
/// if incre is positive, then the index values will be incremented.
/// if incre is negative, then the kewords with index = COLMIN
/// will be deleted and the index of higher numbered keywords will
/// be decremented.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffkshf(
    fptr: *mut fitsfile, /* I - FITS file pointer                        */
    colmin: c_int,       /* I - starting col. to be incremented; 1 = 1st */
    colmax: c_int,       /* I - last column to be incremented            */
    incre: c_int,        /* I - shift index number by this amount        */
    status: *mut c_int,  /* IO - error status                            */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        ffkshf_safe(fptr, colmin, colmax, incre, status)
    }
}

/*--------------------------------------------------------------------------*/
/// shift the index value on any existing column keywords
/// This routine will modify the name of any keyword that begins with 'T'
/// and has an index number in the range COLMIN - COLMAX, inclusive.
///
/// if incre is positive, then the index values will be incremented.
/// if incre is negative, then the kewords with index = COLMIN
/// will be deleted and the index of higher numbered keywords will
/// be decremented.
pub(crate) fn ffkshf_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer                        */
    colmin: c_int,       /* I - starting col. to be incremented; 1 = 1st */
    colmax: c_int,       /* I - last column to be incremented            */
    incre: c_int,        /* I - shift index number by this amount        */
    status: &mut c_int,  /* IO - error status                            */
) -> c_int {
    let mut nkeys: c_int = 0;
    let mut nmore: c_int = 0;
    let nrec: c_int = 0;
    let mut tstatus: c_int = 0;
    let mut i1 = 0;
    let mut ivalue: c_long = 0;
    let mut rec: [c_char; FLEN_CARD] = [0; FLEN_CARD];
    let mut q: [c_char; FLEN_KEYWORD] = [0; FLEN_KEYWORD];
    let mut newkey: [c_char; FLEN_KEYWORD] = [0; FLEN_KEYWORD];

    ffghsp_safe(fptr, Some(&mut nkeys), Some(&mut nmore), status); /* get number of keywords */

    /* go thru header starting with the 9th keyword looking for 'TxxxxNNN' */

    let mut nrec = 9;
    while nrec <= nkeys {
        ffgrec_safe(fptr, nrec as c_int, Some(&mut rec), status);

        if rec[0] == bb(b'T') {
            i1 = 0;
            strncpy_safe(&mut q, &rec[1..], 4);
            if strncmp_safe(&q, cs!(c"BCOL"), 4) == 0
                || strncmp_safe(&q, cs!(c"FORM"), 4) == 0
                || strncmp_safe(&q, cs!(c"TYPE"), 4) == 0
                || strncmp_safe(&q, cs!(c"SCAL"), 4) == 0
                || strncmp_safe(&q, cs!(c"UNIT"), 4) == 0
                || strncmp_safe(&q, cs!(c"NULL"), 4) == 0
                || strncmp_safe(&q, cs!(c"ZERO"), 4) == 0
                || strncmp_safe(&q, cs!(c"DISP"), 4) == 0
                || strncmp_safe(&q, cs!(c"LMIN"), 4) == 0
                || strncmp_safe(&q, cs!(c"LMAX"), 4) == 0
                || strncmp_safe(&q, cs!(c"DMIN"), 4) == 0
                || strncmp_safe(&q, cs!(c"DMAX"), 4) == 0
                || strncmp_safe(&q, cs!(c"CTYP"), 4) == 0
                || strncmp_safe(&q, cs!(c"CRPX"), 4) == 0
                || strncmp_safe(&q, cs!(c"CRVL"), 4) == 0
                || strncmp_safe(&q, cs!(c"CDLT"), 4) == 0
                || strncmp_safe(&q, cs!(c"CROT"), 4) == 0
                || strncmp_safe(&q, cs!(c"CUNI"), 4) == 0
            {
                i1 = 5;
            } else if strncmp_safe(&rec, cs!(c"TDIM"), 4) == 0 {
                i1 = 4;
            }

            if i1 != 0 {
                /* try reading the index number suffix */
                q[0] = 0;
                strncat_safe(&mut q, &rec[i1..], 8 - i1);

                tstatus = 0;
                ffc2ii(&q, &mut ivalue, status);

                if tstatus == 0 && ivalue >= (colmin as c_long) && ivalue <= (colmax as c_long) {
                    if incre <= 0 && ivalue == colmin as c_long {
                        ffdrec_safe(fptr, nrec as c_int, status); /* delete keyword */
                        nkeys -= 1;
                        nrec -= 1;
                    } else {
                        ivalue += incre as c_long;
                        q[0] = 0;
                        strncat_safe(&mut q, &rec, i1);

                        ffkeyn_safe(&q, ivalue as c_int, &mut newkey, status);
                        /* NOTE: because of null termination, it is not
                        equivalent to use strcpy_safe() for the same calls */
                        strncpy_safe(&mut rec, cs!(c"        "), 8); /* erase old keyword name */
                        i1 = strlen_safe(&newkey);
                        strncpy_safe(&mut rec, &newkey, i1); /* overwrite new keyword name */
                        ffmrec_safe(fptr, nrec as c_int, &rec, status); /* modify the record */
                    }
                }
            }
        }
        nrec += 1;
    }
    *status
}

/*--------------------------------------------------------------------------*/
/// Internal function to identify which columns in a binary table are variable length.
///
/// The colnums array will be filled with nvarcols elements - the 1-based numbers
/// of all variable length columns in the table.  This ASSUMES calling function
/// has passed in a colnums array large enough to hold these (colnums==NULL also
/// allowed).
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn fffvcl(
    fptr: *mut fitsfile,  /* I - FITS file pointer                       */
    nvarcols: *mut c_int, /* O - Number of variable length columns found */
    colnums: *mut c_int,  /* O - 1-based variable column positions       */
    status: *mut c_int,   /* IO - error status                           */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let nvarcols = nvarcols.as_mut().expect(NULL_MSG);

        let i = if !fptr.Fptr.tableptr.is_null() {
            fptr.Fptr.tfield
        } else {
            0
        } as usize;

        if colnums.is_null() {
            fffvcl_safe(fptr, nvarcols, None, status)
        } else {
            let mut colnums = slice::from_raw_parts_mut(colnums, i);
            fffvcl_safe(fptr, nvarcols, Some(&mut colnums), status)
        }
    }
}

/*--------------------------------------------------------------------------*/
/// Internal function to identify which columns in a binary table are variable length.
/// The colnums array will be filled with nvarcols elements - the 1-based numbers
/// of all variable length columns in the table.  This ASSUMES calling function
/// has passed in a colnums array large enough to hold these (colnums==NULL also
/// allowed).
pub(crate) fn fffvcl_safe(
    fptr: &mut fitsfile,  /* I - FITS file pointer                       */
    nvarcols: &mut c_int, /* O - Number of variable length columns found */
    mut colnums: Option<&mut [c_int]>, /* O - 1-based variable column positions       */
    status: &mut c_int,   /* IO - error status                           */
) -> c_int {
    let mut tfields: c_int = 0;
    let icol: c_int = 0;
    let mut colptr: *mut tcolumn;

    *nvarcols = 0;

    if *status > 0 {
        return *status;
    }

    if fptr.Fptr.hdutype != BINARY_TBL {
        ffpmsg_str("Var-length column search can only be performed on Binary tables (fffvcl)");
        *status = NOT_BTABLE;
        return *status;
    }

    if !fptr.Fptr.tableptr.is_null() {
        tfields = fptr.Fptr.tfield;
        let colptr = fptr.Fptr.tableptr; /* point to first column structure */
        let c = fptr.Fptr.get_tableptr_as_slice();

        for (ci, icol) in (0..(tfields as usize)).enumerate() {
            /* Condition for variable length column: negative tdatatype */
            if c[ci].tdatatype < 0 {
                if let Some(colnums) = colnums.as_mut() {
                    colnums[*nvarcols as usize] = (icol + 1) as c_int;
                }
                *nvarcols += 1;
            }
        }
    }
    *status
}

/*--------------------------------------------------------------------------*/
/// Shift block of bytes by nshift bytes (positive or negative).
///
/// A positive nshift value moves the block down further in the file, while a
/// negative value shifts the block towards the beginning of the file.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffshft(
    fptr: *mut fitsfile, /* I - FITS file pointer                        */
    firstbyte: LONGLONG, /* I - position of first byte in block to shift */
    nbytes: LONGLONG,    /* I - size of block of bytes to shift          */
    nshift: LONGLONG,    /* I - size of shift in bytes (+ or -)          */
    status: *mut c_int,  /* IO - error status                            */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        ffshft_safe(fptr, firstbyte, nbytes, nshift, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Shift block of bytes by nshift bytes (positive or negative).
///
/// A positive nshift value moves the block down further in the file, while a
/// negative value shifts the block towards the beginning of the file.
pub fn ffshft_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer                        */
    firstbyte: LONGLONG, /* I - position of first byte in block to shift */
    nbytes: LONGLONG,    /* I - size of block of bytes to shift          */
    nshift: LONGLONG,    /* I - size of shift in bytes (+ or -)          */
    status: &mut c_int,  /* IO - error status                            */
) -> c_int {
    const SHFTBUFFSIZE: usize = 100000;
    let mut ntomov: LONGLONG = 0;
    let mut ptr: LONGLONG = 0;
    let mut ntodo: LONGLONG = 0;
    let mut buffer: [c_char; SHFTBUFFSIZE] = [0; SHFTBUFFSIZE];

    if *status > 0 {
        return *status;
    }

    ntodo = nbytes; /* total number of bytes to shift */

    if nshift > 0 {
        /* start at the end of the block and work backwards */
        ptr = firstbyte + nbytes;
    } else {
        /* start at the beginning of the block working forwards */
        ptr = firstbyte;
    }

    while ntodo != 0 {
        /* number of bytes to move at one time */
        ntomov = cmp::min(ntodo, SHFTBUFFSIZE as LONGLONG);

        if nshift > 0 {
            /* if moving block down ... */
            ptr -= ntomov;
        }

        /* move to position and read the bytes to be moved */

        ffmbyt_safe(fptr, ptr, REPORT_EOF, status);
        ffgbyt(fptr, ntomov, cast_slice_mut(&mut buffer), status);

        /* move by shift amount and write the bytes */
        ffmbyt_safe(fptr, ptr + nshift, IGNORE_EOF, status);
        if ffpbyt(fptr, ntomov, cast_slice(&buffer), status) > 0 {
            ffpmsg_str("Error while shifting block (ffshft)");
            return *status;
        }

        ntodo -= ntomov;
        if nshift < 0 {
            /* if moving block up ... */
            ptr += ntomov;
        }
    }

    /* now overwrite the old data with fill */
    if fptr.Fptr.hdutype == ASCII_TBL {
        buffer.fill(32); /* fill ASCII tables with spaces */
    } else {
        buffer.fill(0); /* fill other HDUs with zeros */
    }

    if nshift < 0 {
        ntodo = -nshift;
        /* point to the end of the shifted block */
        ptr = firstbyte + nbytes + nshift;
    } else {
        ntodo = nshift;
        /* point to original beginning of the block */
        ptr = firstbyte;
    }

    ffmbyt_safe(fptr, ptr, REPORT_EOF, status);

    while ntodo != 0 {
        ntomov = cmp::min(ntodo, SHFTBUFFSIZE as LONGLONG);
        ffpbyt(fptr, ntomov, cast_slice(&buffer), status);
        ntodo -= ntomov;
    }
    *status
}
