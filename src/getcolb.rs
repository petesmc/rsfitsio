/*  This file, getcolb.c, contains routines that read data elements from   */
/*  a FITS image or table, with unsigned char (unsigned byte) data type.   */

/*  The FITSIO software was written by William Pence at the High Energy    */
/*  Astrophysic Science Archive Research Center (HEASARC) at the NASA      */
/*  Goddard Space Flight Center.                                           */

use core::slice;
use std::ffi::CStr;
use std::{cmp, mem};

use crate::c_types::{c_char, c_int, c_long, c_short, c_void};

use bytemuck::{cast_slice, cast_slice_mut};

use crate::bb;
use crate::fitscore::{
    ffasfm_safe, ffgcprll, ffghdt_safe, ffmahd_safe, ffpmsg_slice, ffpmsg_str, ffrdef_safe,
    fits_is_compressed_image_safe,
};
use crate::fitsio2::*;
use crate::getcoll::ffgcll;
use crate::wrappers::*;
use crate::{NullCheckType, fitsio::*};
use crate::{buffers::*, calculate_subsection_length};
use crate::{int_snprintf, slice_to_str};

/*--------------------------------------------------------------------------*/
/// Read an array of values from the primary array. Data conversion
/// and scaling will be performed if necessary (e.g, if the datatype of
/// the FITS array is not the same as the array being read).
/// Undefined elements will be set equal to NULVAL, unless NULVAL=0
/// in which case no checking for undefined values will be performed.
/// ANYNUL is returned with a value of .true. if any pixels are undefined.
///
/// The primary array is represented as a binary table:
/// each group of the primary array is a row in the table,
/// where the first column contains the group parameters
/// and the second column contains the image itself.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgpvb(
    fptr: *mut fitsfile, /* I - FITS file pointer                       */
    group: c_long,       /* I - group to read (1 = 1st group)           */
    firstelem: LONGLONG, /* I - first vector element to read (1 = 1st)  */
    nelem: LONGLONG,     /* I - number of values to read                */
    nulval: u8,          /* I - value for undefined pixels          */
    array: *mut u8,      /* O - array of values that are returned   */
    anynul: *mut c_int,  /* O - set to 1 if any values are null; else 0 */
    status: *mut c_int,  /* IO - error status                           */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let anynul = anynul.as_mut();

        let array = slice::from_raw_parts_mut(array, nelem as usize);

        ffgpvb_safe(fptr, group, firstelem, nelem, nulval, array, anynul, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Read an array of values from the primary array. Data conversion
/// and scaling will be performed if necessary (e.g, if the datatype of
/// the FITS array is not the same as the array being read).
/// Undefined elements will be set equal to NULVAL, unless NULVAL=0
/// in which case no checking for undefined values will be performed.
/// ANYNUL is returned with a value of .true. if any pixels are undefined.
///
/// The primary array is represented as a binary table:
/// each group of the primary array is a row in the table,
/// where the first column contains the group parameters
/// and the second column contains the image itself.
pub fn ffgpvb_safe(
    fptr: &mut fitsfile,        /* I - FITS file pointer                       */
    group: c_long,              /* I - group to read (1 = 1st group)           */
    firstelem: LONGLONG,        /* I - first vector element to read (1 = 1st)  */
    nelem: LONGLONG,            /* I - number of values to read                */
    nulval: u8,                 /* I - value for undefined pixels          */
    array: &mut [u8],           /* O - array of values that are returned   */
    anynul: Option<&mut c_int>, /* O - set to 1 if any values are null; else 0 */
    status: &mut c_int,         /* IO - error status                           */
) -> c_int {
    let cdummy = 0;
    let nullcheck = NullCheckType::SetPixel;
    let mut nullvalue: u8 = 0;

    let mut dummy_nularray = vec![0; nelem as usize];

    if fits_is_compressed_image_safe(fptr, status) > 0 {
        /* this is a compressed image in a binary table */
        nullvalue = nulval; /* set local variable */

        todo!();
        //fits_read_compressed_pixels(fptr, TBYTE, firstelem, nelem, nullcheck, &nullvalue, cast_slice_mut(array), None, anynul, status);
        return *status;
    }

    let row = cmp::max(1, group);

    ffgclb(
        fptr,
        2,
        row as LONGLONG,
        firstelem,
        nelem,
        1,
        NullCheckType::SetPixel,
        nulval,
        array,
        &mut dummy_nularray,
        anynul,
        status,
    );
    *status
}

/*--------------------------------------------------------------------------*/
/// Read an array of values from the primary array. Data conversion
/// and scaling will be performed if necessary (e.g, if the datatype of
/// the FITS array is not the same as the array being read).
/// Any undefined pixels in the returned array will be set = 0 and the
/// corresponding nularray value will be set = 1.
/// ANYNUL is returned with a value of .true. if any pixels are undefined.
///
/// The primary array is represented as a binary table:
/// each group of the primary array is a row in the table,
/// where the first column contains the group parameters
/// and the second column contains the image itself.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgpfb(
    fptr: *mut fitsfile,   /* I - FITS file pointer                       */
    group: c_long,         /* I - group to read (1 = 1st group)           */
    firstelem: LONGLONG,   /* I - first vector element to read (1 = 1st)  */
    nelem: LONGLONG,       /* I - number of values to read                */
    array: *mut u8,        /* O - array of values that are returned   */
    nularray: *mut c_char, /* O - array of null pixel flags               */
    anynul: *mut c_int,    /* O - set to 1 if any values are null; else 0 */
    status: *mut c_int,    /* IO - error status                           */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let anynul = anynul.as_mut();

        let array = slice::from_raw_parts_mut(array, nelem as usize);
        let nularray = slice::from_raw_parts_mut(nularray, nelem as usize);

        ffgpfb_safe(
            fptr, group, firstelem, nelem, array, nularray, anynul, status,
        )
    }
}

/*--------------------------------------------------------------------------*/
/// Read an array of values from the primary array. Data conversion
/// and scaling will be performed if necessary (e.g, if the datatype of
/// the FITS array is not the same as the array being read).
/// Any undefined pixels in the returned array will be set = 0 and the
/// corresponding nularray value will be set = 1.
/// ANYNUL is returned with a value of .true. if any pixels are undefined.
///
/// The primary array is represented as a binary table:
/// each group of the primary array is a row in the table,
/// where the first column contains the group parameters
/// and the second column contains the image itself.
pub fn ffgpfb_safe(
    fptr: &mut fitsfile,        /* I - FITS file pointer                       */
    group: c_long,              /* I - group to read (1 = 1st group)           */
    firstelem: LONGLONG,        /* I - first vector element to read (1 = 1st)  */
    nelem: LONGLONG,            /* I - number of values to read                */
    array: &mut [u8],           /* O - array of values that are returned   */
    nularray: &mut [c_char],    /* O - array of null pixel flags               */
    anynul: Option<&mut c_int>, /* O - set to 1 if any values are null; else 0 */
    status: &mut c_int,         /* IO - error status                           */
) -> c_int {
    let nullcheck = NullCheckType::SetNullArray;

    if fits_is_compressed_image_safe(fptr, status) > 0 {
        /* this is a compressed image in a binary table */

        todo!();
        //fits_read_compressed_pixels(fptr, TBYTE, firstelem, nelem,          nullcheck, None, cast_slice_mut(array), nularray, anynul, status);
        return *status;
    }

    let row = cmp::max(1, group);

    ffgclb(
        fptr,
        2,
        row as LONGLONG,
        firstelem,
        nelem,
        1,
        NullCheckType::SetNullArray,
        0,
        array,
        cast_slice_mut(nularray),
        anynul,
        status,
    );
    *status
}

/*--------------------------------------------------------------------------*/
/// Read an entire 2-D array of values to the primary array. Data conversion
/// and scaling will be performed if necessary (e.g, if the datatype of the
/// FITS array is not the same as the array being read).  Any null
/// values in the array will be set equal to the value of nulval, unless
/// nulval = 0 in which case no null checking will be performed.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffg2db(
    fptr: *mut fitsfile, /* I - FITS file pointer                       */
    group: c_long,       /* I - group to read (1 = 1st group)           */
    nulval: u8,          /* set undefined pixels equal to this     */
    ncols: LONGLONG,     /* I - number of pixels in each row of array   */
    naxis1: LONGLONG,    /* I - FITS image NAXIS1 value                 */
    naxis2: LONGLONG,    /* I - FITS image NAXIS2 value                 */
    array: *mut u8,      /* O - array to be filled and returned    */
    anynul: *mut c_int,  /* O - set to 1 if any values are null; else 0 */
    status: *mut c_int,  /* IO - error status                           */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let anynul = anynul.as_mut();
        let array = slice::from_raw_parts_mut(array, (ncols * naxis2 * naxis2) as usize);

        ffg2db_safe(
            fptr, group, nulval, ncols, naxis1, naxis2, array, anynul, status,
        )
    }
}

/*--------------------------------------------------------------------------*/
/// Read an entire 2-D array of values to the primary array. Data conversion
/// and scaling will be performed if necessary (e.g, if the datatype of the
/// FITS array is not the same as the array being read).  Any null
/// values in the array will be set equal to the value of nulval, unless
/// nulval = 0 in which case no null checking will be performed.
pub fn ffg2db_safe(
    fptr: &mut fitsfile,        /* I - FITS file pointer                       */
    group: c_long,              /* I - group to read (1 = 1st group)           */
    nulval: u8,                 /* set undefined pixels equal to this     */
    ncols: LONGLONG,            /* I - number of pixels in each row of array   */
    naxis1: LONGLONG,           /* I - FITS image NAXIS1 value                 */
    naxis2: LONGLONG,           /* I - FITS image NAXIS2 value                 */
    array: &mut [u8],           /* O - array to be filled and returned    */
    anynul: Option<&mut c_int>, /* O - set to 1 if any values are null; else 0 */
    status: &mut c_int,         /* IO - error status                           */
) -> c_int {
    /* call the 3D reading routine, with the 3rd dimension = 1 */
    ffg3db_safe(
        fptr, group, nulval, ncols, naxis2, naxis1, naxis2, 1, array, anynul, status,
    );

    *status
}

/*--------------------------------------------------------------------------*/
/// Read an entire 3-D array of values to the primary array. Data conversion
/// and scaling will be performed if necessary (e.g, if the datatype of the
/// FITS array is not the same as the array being read).  Any null
/// values in the array will be set equal to the value of nulval, unless
/// nulval = 0 in which case no null checking will be performed.
///
/// The primary array is represented as a binary table:
/// each group of the primary array is a row in the table,
/// where the first column contains the group parameters
/// and the second column contains the image itself.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffg3db(
    fptr: *mut fitsfile, /* I - FITS file pointer                       */
    group: c_long,       /* I - group to read (1 = 1st group)           */
    nulval: u8,          /* set undefined pixels equal to this     */
    ncols: LONGLONG,     /* I - number of pixels in each row of array   */
    nrows: LONGLONG,     /* I - number of rows in each plane of array   */
    naxis1: LONGLONG,    /* I - FITS image NAXIS1 value                 */
    naxis2: LONGLONG,    /* I - FITS image NAXIS2 value                 */
    naxis3: LONGLONG,    /* I - FITS image NAXIS3 value                 */
    array: *mut u8,      /* O - array to be filled and returned    */
    anynul: *mut c_int,  /* O - set to 1 if any values are null; else 0 */
    status: *mut c_int,  /* IO - error status                           */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let anynul = anynul.as_mut();
        let array = slice::from_raw_parts_mut(array, (ncols * naxis2 * naxis3) as usize);

        ffg3db_safe(
            fptr, group, nulval, ncols, nrows, naxis1, naxis2, naxis3, array, anynul, status,
        )
    }
}

/*--------------------------------------------------------------------------*/
/// Read an entire 3-D array of values to the primary array. Data conversion
/// and scaling will be performed if necessary (e.g, if the datatype of the
/// FITS array is not the same as the array being read).  Any null
/// values in the array will be set equal to the value of nulval, unless
/// nulval = 0 in which case no null checking will be performed.
///
/// The primary array is represented as a binary table:
/// each group of the primary array is a row in the table,
/// where the first column contains the group parameters
/// and the second column contains the image itself.
pub fn ffg3db_safe(
    fptr: &mut fitsfile,            /* I - FITS file pointer                       */
    group: c_long,                  /* I - group to read (1 = 1st group)           */
    nulval: u8,                     /* set undefined pixels equal to this     */
    ncols: LONGLONG,                /* I - number of pixels in each row of array   */
    nrows: LONGLONG,                /* I - number of rows in each plane of array   */
    naxis1: LONGLONG,               /* I - FITS image NAXIS1 value                 */
    naxis2: LONGLONG,               /* I - FITS image NAXIS2 value                 */
    naxis3: LONGLONG,               /* I - FITS image NAXIS3 value                 */
    array: &mut [u8],               /* O - array to be filled and returned    */
    mut anynul: Option<&mut c_int>, /* O - set to 1 if any values are null; else 0 */
    status: &mut c_int,             /* IO - error status                           */
) -> c_int {
    let mut narray = 0;
    let mut nfits: LONGLONG = 0;
    let cdummy: c_char = 0;
    let nullcheck = NullCheckType::SetPixel;
    let inc: [c_long; 3] = [1; 3];
    let fpixel: [c_long; 3] = [1; 3];
    let mut lpixel: [c_long; 3] = [0; 3];
    let mut nullvalue: u8 = 0;

    let mut dummy_nularray = vec![0; (ncols * naxis2 * naxis3) as usize];

    if fits_is_compressed_image_safe(fptr, status) > 0 {
        /* this is a compressed image in a binary table */

        lpixel[0] = ncols as c_long;
        lpixel[1] = nrows as c_long;
        lpixel[2] = naxis3 as c_long;
        nullvalue = nulval; /* set local variable */
        todo!();
        //fits_read_compressed_img(fptr, TBYTE, fpixel, lpixel, inc,            nullcheck, &nullvalue, cast_slice_mut(array), None, anynul, status);
        return *status;
    }

    let tablerow = cmp::max(1, group);

    if ncols == naxis1 && nrows == naxis2 {
        /* arrays have same size? */

        /* all the image pixels are contiguous, so read all at once */
        ffgclb(
            fptr,
            2,
            tablerow as LONGLONG,
            1,
            naxis1 * naxis2 * naxis3,
            1,
            NullCheckType::SetPixel,
            nulval,
            array,
            &mut dummy_nularray,
            anynul.as_deref_mut(),
            status,
        );
        return *status;
    }

    if ncols < naxis1 || nrows < naxis2 {
        *status = BAD_DIMEN;
        return *status;
    }

    nfits = 1; /* next pixel in FITS image to read */
    narray = 0; /* next pixel in output array to be filled */

    /* loop over naxis3 planes in the data cube */
    for jj in 0..(naxis3 as usize) {
        /* loop over the naxis2 rows in the FITS image, */
        /* reading naxis1 pixels to each row            */

        for ii in 0..(naxis2 as usize) {
            if ffgclb(
                fptr,
                2,
                tablerow as LONGLONG,
                nfits,
                naxis1,
                1,
                NullCheckType::SetPixel,
                nulval,
                &mut array[narray..],
                &mut dummy_nularray,
                anynul.as_deref_mut(),
                status,
            ) > 0
            {
                return *status;
            }

            nfits += naxis1;
            narray += ncols as usize;
        }
        narray += ((nrows - naxis2) * ncols) as usize;
    }

    *status
}

/*--------------------------------------------------------------------------*/
/// Read a subsection of data values from an image or a table column.
/// This routine is set up to handle a maximum of nine dimensions.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgsvb(
    fptr: *mut fitsfile,  /* I - FITS file pointer                         */
    colnum: c_int,        /* I - number of the column to read (1 = 1st)    */
    naxis: c_int,         /* I - number of dimensions in the FITS array    */
    naxes: *const c_long, /* I - size of each dimension                    */
    blc: *const c_long,   /* I - 'bottom left corner' of the subsection    */
    trc: *const c_long,   /* I - 'top right corner' of the subsection      */
    inc: *const c_long,   /* I - increment to be applied in each dimension */
    nulval: u8,           /* I - value to set undefined pixels       */
    array: *mut u8,       /* O - array to be filled and returned     */
    anynul: *mut c_int,   /* O - set to 1 if any values are null; else 0   */
    status: *mut c_int,   /* IO - error status                             */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let anynul = anynul.as_mut();

        let naxes = slice::from_raw_parts(naxes, naxis as usize);
        let blc = slice::from_raw_parts(blc, naxis as usize);
        let trc = slice::from_raw_parts(trc, naxis as usize);
        let inc = slice::from_raw_parts(inc, naxis as usize);

        let total_nelem = calculate_subsection_length(blc, trc, inc);

        let array = slice::from_raw_parts_mut(array, total_nelem);

        ffgsvb_safe(
            fptr, colnum, naxis, naxes, blc, trc, inc, nulval, array, anynul, status,
        )
    }
}

/*--------------------------------------------------------------------------*/
/// Read a subsection of data values from an image or a table column.
/// This routine is set up to handle a maximum of nine dimensions.
pub fn ffgsvb_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer                         */
    colnum: c_int,       /* I - number of the column to read (1 = 1st)    */
    naxis: c_int,        /* I - number of dimensions in the FITS array    */
    naxes: &[c_long],    /* I - size of each dimension                    */
    blc: &[c_long],      /* I - 'bottom left corner' of the subsection    */
    trc: &[c_long],      /* I - 'top right corner' of the subsection      */
    inc: &[c_long],      /* I - increment to be applied in each dimension */
    nulval: u8,          /* I - value to set undefined pixels       */
    array: &mut [u8],    /* O - array to be filled and returned     */
    mut anynul: Option<&mut c_int>, /* O - set to 1 if any values are null; else 0   */
    status: &mut c_int,  /* IO - error status                             */
) -> c_int {
    let mut rstr: c_long = 0;
    let mut rstp: c_long = 0;
    let mut rinc: c_long = 0;
    let mut str: [c_long; 9] = [0; 9];
    let mut stp: [c_long; 9] = [0; 9];
    let mut incr: [c_long; 9] = [0; 9];
    let mut dir: [c_long; 9] = [0; 9];
    let mut nelem: c_long = 0;
    let mut nultyp = NullCheckType::None;
    let mut ninc: c_long = 0;
    let mut numcol: c_long = 0;
    let mut felem: LONGLONG = 0;
    let mut dsize: [LONGLONG; 10] = [0; 10];
    let mut blcll: [LONGLONG; 9] = [0; 9];
    let mut trcll: [LONGLONG; 9] = [0; 9];
    let mut hdutype: c_int = 0;
    let mut anyf: c_int = 0;
    let ldummy = 0;
    let mut msg: [c_char; FLEN_ERRMSG] = [0; FLEN_ERRMSG];
    let nullcheck = NullCheckType::SetPixel;
    let mut nullvalue: u8 = 0;

    let naxis = naxis as usize;

    if naxis < 1 || naxis > 9 {
        int_snprintf!(
            &mut msg,
            FLEN_ERRMSG,
            "NAXIS = {} in call to ffgsvb is out of range",
            naxis,
        );
        ffpmsg_slice(&msg);
        *status = BAD_DIMEN;
        return *status;
    }

    if fits_is_compressed_image_safe(fptr, status) > 0 {
        /* this is a compressed image in a binary table */

        for ii in 0..naxis {
            blcll[ii] = blc[ii] as LONGLONG;
            trcll[ii] = trc[ii] as LONGLONG;
        }

        nullvalue = nulval; /* set local variable */
        todo!();
        // fits_read_compressed_img(fptr, TBYTE, blcll, trcll, inc,            nullcheck, &nullvalue, cast_slice_mut(array), None, anynul, status);
        return *status;
    }

    /*
    if this is a primary array, then the input COLNUM parameter should
    be interpreted as the row number, and we will alway read the image
    data from column 2 (any group parameters are in column 1).
    */
    if ffghdt_safe(fptr, &mut hdutype, status) > 0 {
        return *status;
    }

    if hdutype == IMAGE_HDU {
        /* this is a primary array, or image extension */
        if colnum == 0 {
            rstr = 1;
            rstp = 1;
        } else {
            rstr = colnum as c_long;
            rstp = colnum as c_long;
        }
        rinc = 1;
        numcol = 2;
    } else {
        /* this is a table, so the row info is in the (naxis+1) elements */
        rstr = blc[naxis];
        rstp = trc[naxis];
        rinc = inc[naxis];
        numcol = colnum as c_long;
    }

    nultyp = NullCheckType::SetPixel;

    if let Some(anynul) = anynul.as_deref_mut() {
        *anynul = FALSE as c_int;
    }

    let mut i0 = 0;
    for ii in 0..9 {
        str[ii] = 1;
        stp[ii] = 1;
        incr[ii] = 1;
        dsize[ii] = 1;
        dir[ii] = 1;
    }

    for ii in 0..(naxis) {
        if trc[ii] < blc[ii] {
            if hdutype == IMAGE_HDU {
                dir[ii] = -1;
            } else {
                int_snprintf!(
                    &mut msg,
                    FLEN_ERRMSG,
                    "ffgsvb: illegal range specified for axis {}",
                    ii + 1,
                );
                ffpmsg_slice(&msg);
                *status = BAD_PIX_NUM;
                return *status;
            }
        }

        str[ii] = blc[ii];
        stp[ii] = trc[ii];
        incr[ii] = inc[ii];
        dsize[ii + 1] = dsize[ii] * naxes[ii] as LONGLONG;
        dsize[ii] *= dir[ii] as LONGLONG;
    }
    dsize[naxis] *= dir[naxis] as LONGLONG;

    if naxis == 1 && naxes[0] == 1 {
        /* This is not a vector column, so read all the rows at once */
        nelem = (rstp - rstr) / rinc + 1;
        ninc = rinc;
        rstp = rstr;
    } else {
        /* have to read each row individually, in all dimensions */
        nelem = (stp[0] * dir[0] - str[0] * dir[0]) / inc[0] + 1;
        ninc = incr[0] * dir[0];
    }

    for row in (rstr..=rstp).step_by(rinc as usize) {
        for i8 in ((str[8] * dir[8])..=(stp[8] * dir[8])).step_by(incr[8] as usize) {
            for i7 in ((str[7] * dir[7])..=(stp[7] * dir[7])).step_by(incr[7] as usize) {
                for i6 in ((str[6] * dir[6])..=(stp[6] * dir[6])).step_by(incr[6] as usize) {
                    for i5 in ((str[5] * dir[5])..=(stp[5] * dir[5])).step_by(incr[5] as usize) {
                        for i4 in ((str[4] * dir[4])..=(stp[4] * dir[4])).step_by(incr[4] as usize)
                        {
                            for i3 in
                                ((str[3] * dir[3])..=(stp[3] * dir[3])).step_by(incr[3] as usize)
                            {
                                for i2 in ((str[2] * dir[2])..=(stp[2] * dir[2]))
                                    .step_by(incr[2] as usize)
                                {
                                    for i1 in ((str[1] * dir[1])..=(stp[1] * dir[1]))
                                        .step_by(incr[1] as usize)
                                    {
                                        felem = (str[0] as LONGLONG)
                                            + (i1 - dir[1]) as LONGLONG * dsize[1]
                                            + (i2 - dir[2]) as LONGLONG * dsize[2]
                                            + (i3 - dir[3]) as LONGLONG * dsize[3]
                                            + (i4 - dir[4]) as LONGLONG * dsize[4]
                                            + (i5 - dir[5]) as LONGLONG * dsize[5]
                                            + (i6 - dir[6]) as LONGLONG * dsize[6]
                                            + (i7 - dir[7]) as LONGLONG * dsize[7]
                                            + (i8 - dir[8]) as LONGLONG * dsize[8];

                                        if ffgclb(
                                            fptr,
                                            numcol as c_int,
                                            row as LONGLONG,
                                            felem,
                                            nelem as LONGLONG,
                                            ninc,
                                            nultyp,
                                            nulval,
                                            &mut array[i0..],
                                            &mut [ldummy],
                                            Some(&mut anyf),
                                            status,
                                        ) > 0
                                        {
                                            return *status;
                                        }

                                        if anyf > 0 {
                                            if let Some(anynul) = anynul.as_deref_mut() {
                                                *anynul = TRUE as c_int;
                                            }
                                        }
                                        i0 += nelem as usize;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    *status
}

pub fn ffgsfb_safe(
    fptr: &mut fitsfile,    /* I - FITS file pointer                         */
    colnum: c_int,          /* I - number of the column to read (1 = 1st)    */
    naxis: c_int,           /* I - number of dimensions in the FITS array    */
    naxes: &[c_long],       /* I - size of each dimension                    */
    blc: &[c_long],         /* I - 'bottom left corner' of the subsection    */
    trc: &[c_long],         /* I - 'top right corner' of the subsection      */
    inc: &[c_long],         /* I - increment to be applied in each dimension */
    array: &mut [u8],       /* O - array to be filled and returned     */
    flagval: &mut [c_char], /* O - set to 1 if corresponding value is null   */
    mut anynul: Option<&mut c_int>, /* O - set to 1 if any values are null; else 0   */
    status: &mut c_int,     /* IO - error status                             */
) -> c_int {
    let mut rstr: c_long = 0;
    let mut rstp: c_long = 0;
    let mut rinc: c_long = 0;
    let mut str: [c_long; 9] = [0; 9];
    let mut stp: [c_long; 9] = [0; 9];
    let mut incr: [c_long; 9] = [0; 9];
    let dir: [c_long; 9] = [0; 9];
    let mut nelem: c_long = 0;
    let mut nultyp = NullCheckType::None;
    let mut ninc: c_long = 0;
    let mut numcol: c_long = 0;
    let mut felem: LONGLONG = 0;
    let mut dsize: [LONGLONG; 10] = [0; 10];
    let mut blcll: [LONGLONG; 9] = [0; 9];
    let mut trcll: [LONGLONG; 9] = [0; 9];
    let mut hdutype: c_int = 0;
    let anyf: c_int = 0;
    let ldummy: c_char = 0;
    let mut msg: [c_char; FLEN_ERRMSG] = [0; FLEN_ERRMSG];
    let nullcheck = NullCheckType::SetNullArray;
    let nullval: u8 = 0;

    let naxis = naxis as usize;

    if naxis < 1 || naxis > 9 {
        int_snprintf!(
            &mut msg,
            FLEN_ERRMSG,
            "NAXIS = {} in call to ffgsfb is out of range",
            naxis,
        );
        ffpmsg_slice(&msg);
        *status = BAD_DIMEN;
        return *status;
    }

    if fits_is_compressed_image_safe(fptr, status) > 0 {
        /* this is a compressed image in a binary table */

        for ii in 0..naxis {
            blcll[ii] = blc[ii] as LONGLONG;
            trcll[ii] = trc[ii] as LONGLONG;
        }

        todo!();
        // fits_read_compressed_img(fptr, TBYTE, blcll, trcll, inc,     nullcheck, None, cast_slice_mut(array), flagval, anynul, status);
        return *status;
    }

    /*
    if this is a primary array, then the input COLNUM parameter should
    be interpreted as the row number, and we will alway read the image
    data from column 2 (any group parameters are in column 1).
    */
    if ffghdt_safe(fptr, &mut hdutype, status) > 0 {
        return *status;
    }

    if hdutype == IMAGE_HDU {
        /* this is a primary array, or image extension */
        if colnum == 0 {
            rstr = 1;
            rstp = 1;
        } else {
            rstr = colnum as c_long;
            rstp = colnum as c_long;
        }
        rinc = 1;
        numcol = 2;
    } else {
        /* this is a table, so the row info is in the (naxis+1) elements */
        rstr = blc[naxis];
        rstp = trc[naxis];
        rinc = inc[naxis];
        numcol = colnum as c_long;
    }

    nultyp = NullCheckType::SetNullArray;

    if let Some(anynul) = anynul.as_deref_mut() {
        *anynul = FALSE as c_int;
    }

    let mut i0 = 0;
    for ii in 0..9 {
        str[ii] = 1;
        stp[ii] = 1;
        incr[ii] = 1;
        dsize[ii] = 1;
    }

    for ii in 0..(naxis) {
        if trc[ii] < blc[ii] {
            if hdutype == IMAGE_HDU {
                /*
                support negative strides for FITS primary arrays
                (not tables, because the FITSIO column handlers don't
                support it).
                */
                str[ii] = blc[ii];
                stp[ii] = trc[ii];
                incr[ii] = inc[ii];
                dsize[ii + 1] = dsize[ii] * naxes[ii] as LONGLONG;
            } else {
                int_snprintf!(
                    &mut msg,
                    FLEN_ERRMSG,
                    "ffgsfb: illegal range specified for axis {}",
                    ii + 1,
                );
                ffpmsg_slice(&msg);
                *status = BAD_PIX_NUM;
                return *status;
            }
        } else {
            str[ii] = blc[ii];
            stp[ii] = trc[ii];
            incr[ii] = inc[ii];
            dsize[ii + 1] = dsize[ii] * naxes[ii] as LONGLONG;
        }
    }
    nelem = 1;
    for ii in 0..naxis {
        nelem *= ((stp[ii] - str[ii]) / inc[ii]) + 1;
    }

    i0 = 0;

    if hdutype == IMAGE_HDU {
        felem = str[0];
    } else {
        felem = rstr + (str[0] - 1) * (rstp - rstr + 1) / rinc;
    }
    let mut hh = str[0];

    for jj in 1..naxis {
        felem += (str[jj] - 1) * dsize[jj];
        hh += (str[jj] - 1) * dsize[jj];
    }

    /* determine the number of pixels to process in each loop */
    if naxis == 1 {
        ninc = incr[0];
    } else {
        ninc = incr[0];
        for jj in 1..naxis {
            if (stp[jj] - str[jj]) / inc[jj] > 0 {
                ninc = incr[0];
            } else {
                ninc = cmp::min(ninc, nelem);
            }
        }
    }
    ninc = cmp::min(ninc, nelem);

    if hdutype == IMAGE_HDU {
        ffmbyt_safe(
            fptr,
            (felem - 1) * mem::size_of::<c_char>() as c_long,
            REPORT_EOF,
            status,
        );
        if *status > 0 {
            return *status;
        }
    }

    /* read the null value parameters */
    // TODO: Implement table null value reading
    // let mut tnull: c_long = 0;
    // if (hdutype != IMAGE_HDU)
    //     && ffgtnlll(fptr, numcol, &mut tnull, status) > 0
    // {
    //     *status = COL_NOT_FOUND;
    //     return *status;
    // }

    let anynul_int = 0;
    let nularray = vec![0u8; ninc as usize];

    i0 = 0;
    while nelem > 0 && *status <= 0 {
        /* read the next subset of pixels */
        if hdutype == IMAGE_HDU {
            if incr[0] != 1 {
                ffgbytoff(
                    fptr,
                    1,           // gsize
                    ninc,        // ngroups
                    incr[0] - 1, // offset
                    &mut array[i0 as usize..i0 as usize + ninc as usize],
                    status,
                );
            } else {
                ffgbyt(
                    fptr,
                    ninc as LONGLONG,
                    &mut array[i0 as usize..i0 as usize + ninc as usize],
                    status,
                );
            }
        } else {
            /* read from a table - this should use the nested loop structure from the original */
            /* For now, return an error indicating this is not implemented */
            ffpmsg_str("ffgsfb_safe: table reading not yet implemented");
            *status = COL_NOT_FOUND;
            return *status;
        }

        if *status != 0 {
            /* test for EOF */
            if *status == END_OF_FILE {
                if hdutype > 0 || (i0 + 1) >= nelem {
                    *status = 0; /* reading from image extension or */
                }
                /* single pixel table so ignore EOF */
                else {
                    return *status;
                }
            } else if *status == ARRAY_TOO_BIG {
                /* ignore this error message */
                *status = 0;
            } else {
                return *status;
            }
        }

        /* increment the counters for the next loop */
        i0 += ninc;
        let remain = nelem - i0;

        let nextelem = ninc;
        ninc = cmp::min(ninc, remain); /* don't exceed the maximum */

        if incr[0] == 1
            && naxis > 1
            && (felem + nextelem - 1 - dsize[1]) / dsize[1] != (felem - 1 - dsize[1]) / dsize[1]
        {
            /* we have reached the boundary between successive planes */
            felem += dsize[1] - (felem - 1 - dsize[1]) % dsize[1];

            /* recalculate the indices of the next element to read */
            for kk in 1..naxis {
                hh /= dsize[kk];
                if hh == ((str[kk] + ((stp[kk] - str[kk]) / inc[kk]) * inc[kk]) / naxes[kk - 1]) {
                    str[kk] += 1;
                    hh = 1;
                    for ll in kk + 1..naxis {
                        hh += (str[ll] - 1) / naxes[ll - 1];
                    }
                    if kk == naxis - 1 {
                        ninc = incr[0]; /* completed a row */
                    }
                    break;
                }
            }
        } else {
            felem += ninc;
        }

        nelem = remain;
    }

    // TODO: Implement null checking for images
    // if anynul.is_some() {
    //     if nultyp == NullCheckType::SetNullArray {
    //         /* this is an image, so check entire array for nulls */
    //         // Implement null checking logic here
    //     }
    // }

    if let Some(anynul) = anynul {
        *anynul = anynul_int;
    }

    *status
}

/*--------------------------------------------------------------------------*/
/// Read a subsection of data values from an image or a table column.
/// This routine is set up to handle a maximum of nine dimensions.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgsfb(
    fptr: *mut fitsfile,  /* I - FITS file pointer                         */
    colnum: c_int,        /* I - number of the column to read (1 = 1st)    */
    naxis: c_int,         /* I - number of dimensions in the FITS array    */
    naxes: *const c_long, /* I - size of each dimension                    */
    blc: *const c_long,   /* I - 'bottom left corner' of the subsection    */
    trc: *const c_long,   /* I - 'top right corner' of the subsection      */
    inc: *const c_long,   /* I - increment to be applied in each dimension */
    array: *mut u8,       /* O - array to be filled and returned     */
    flagval: *mut c_char, /* O - set to 1 if corresponding value is null   */
    anynul: *mut c_int,   /* O - set to 1 if any values are null; else 0   */
    status: *mut c_int,   /* IO - error status                             */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let anynul = anynul.as_mut();

        let naxes = slice::from_raw_parts(naxes, naxis as usize);
        let blc = slice::from_raw_parts(blc, naxis as usize);
        let trc = slice::from_raw_parts(trc, naxis as usize);
        let inc = slice::from_raw_parts(inc, naxis as usize);

        let total_nelem = calculate_subsection_length(blc, trc, inc);

        let array = slice::from_raw_parts_mut(array, total_nelem);
        let flagval = slice::from_raw_parts_mut(flagval, total_nelem);

        ffgsfb_safe(
            fptr, colnum, naxis, naxes, blc, trc, inc, array, flagval, anynul, status,
        )
    }
}

/*--------------------------------------------------------------------------*/
///Read an array of group parameters from the primary array. Data conversion
///and scaling will be performed if necessary (e.g, if the datatype of
///the FITS array is not the same as the array being read).
///
///The primary array is represented as a binary table:
///each group of the primary array is a row in the table,
///where the first column contains the group parameters
///and the second column contains the image itself.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffggpb(
    fptr: *mut fitsfile, /* I - FITS file pointer                       */
    group: c_long,       /* I - group to read (1 = 1st group)           */
    firstelem: c_long,   /* I - first vector element to read (1 = 1st)  */
    nelem: c_long,       /* I - number of values to read                */
    array: *mut u8,      /* O - array of values that are returned   */
    status: *mut c_int,  /* IO - error status                           */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let array = slice::from_raw_parts_mut(array, nelem as usize);

        ffggpb_safe(fptr, group, firstelem, nelem, array, status)
    }
}

/// Safe wrapper for ffggpb - Read group parameters as unsigned bytes
///
/// Read an array of group parameters from the primary array.
/// The primary array is represented as a binary table where each group
/// is a row, with the first column containing group parameters.
pub fn ffggpb_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer                       */
    group: c_long,       /* I - group to read (1 = 1st group)           */
    firstelem: c_long,   /* I - first vector element to read (1 = 1st)  */
    nelem: c_long,       /* I - number of values to read                */
    array: &mut [u8],    /* O - array of values that are returned   */
    status: &mut c_int,  /* IO - error status                           */
) -> c_int {
    let cdummy = 0;
    let mut anynul = 0;
    let mut dummy_nularray = vec![0; nelem as usize];

    let row = cmp::max(1, group);

    ffgclb(
        fptr,
        1,
        row as LONGLONG,
        firstelem as LONGLONG,
        nelem as LONGLONG,
        1,
        NullCheckType::SetPixel,
        cdummy,
        array,
        &mut dummy_nularray,
        Some(&mut anynul),
        status,
    );
    *status
}

/*--------------------------------------------------------------------------*/
/// Read an array of values from a column in the current FITS HDU. Automatic
/// datatype conversion will be performed if the datatype of the column does not
/// match the datatype of the array parameter. The output values will be scaled
/// by the FITS TSCALn and TZEROn values if these values have been defined.
/// Any undefined pixels will be set equal to the value of 'nulval' unless
/// nulval = 0 in which case no checks for undefined pixels will be made.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgcvb(
    fptr: *mut fitsfile, /* I - FITS file pointer                       */
    colnum: c_int,       /* I - number of column to read (1 = 1st col)  */
    firstrow: LONGLONG,  /* I - first row to read (1 = 1st row)         */
    firstelem: LONGLONG, /* I - first vector element to read (1 = 1st)  */
    nelem: LONGLONG,     /* I - number of values to read                */
    nulval: u8,          /* I - value for null pixels               */
    array: *mut u8,      /* O - array of values that are read       */
    anynul: *mut c_int,  /* O - set to 1 if any values are null; else 0 */
    status: *mut c_int,  /* IO - error status                           */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let anynul = anynul.as_mut();
        let array = slice::from_raw_parts_mut(array, nelem as usize);

        ffgcvb_safe(
            fptr,
            colnum,
            firstrow as LONGLONG,
            firstelem as LONGLONG,
            nelem as LONGLONG,
            nulval,
            array,
            anynul,
            status,
        )
    }
}

/*--------------------------------------------------------------------------*/
/// Read an array of values from a column in the current FITS HDU. Automatic
/// datatype conversion will be performed if the datatype of the column does not
/// match the datatype of the array parameter. The output values will be scaled
/// by the FITS TSCALn and TZEROn values if these values have been defined.
/// Any undefined pixels will be set equal to the value of 'nulval' unless
/// nulval = 0 in which case no checks for undefined pixels will be made.
pub fn ffgcvb_safe(
    fptr: &mut fitsfile,        /* I - FITS file pointer                       */
    colnum: c_int,              /* I - number of column to read (1 = 1st col)  */
    firstrow: LONGLONG,         /* I - first row to read (1 = 1st row)         */
    firstelem: LONGLONG,        /* I - first vector element to read (1 = 1st)  */
    nelem: LONGLONG,            /* I - number of values to read                */
    nulval: u8,                 /* I - value for null pixels               */
    array: &mut [u8],           /* O - array of values that are read       */
    anynul: Option<&mut c_int>, /* O - set to 1 if any values are null; else 0 */
    status: &mut c_int,         /* IO - error status                           */
) -> c_int {
    let cdummy = 0;

    let mut dummy_nularray = vec![0; (nelem) as usize];

    ffgclb(
        fptr,
        colnum,
        firstrow,
        firstelem,
        nelem,
        1,
        NullCheckType::SetPixel,
        nulval,
        array,
        &mut dummy_nularray,
        anynul,
        status,
    );
    *status
}

/*--------------------------------------------------------------------------*/
/// Read an array of values from a column in the current FITS HDU. Automatic
/// datatype conversion will be performed if the datatype of the column does not
/// match the datatype of the array parameter. The output values will be scaled
/// by the FITS TSCALn and TZEROn values if these values have been defined.
/// Nularray will be set = 1 if the corresponding array pixel is undefined,
/// otherwise nularray will = 0.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgcfb(
    fptr: *mut fitsfile,   /* I - FITS file pointer                       */
    colnum: c_int,         /* I - number of column to read (1 = 1st col)  */
    firstrow: LONGLONG,    /* I - first row to read (1 = 1st row)         */
    firstelem: LONGLONG,   /* I - first vector element to read (1 = 1st)  */
    nelem: LONGLONG,       /* I - number of values to read                */
    array: *mut u8,        /* O - array of values that are read       */
    nularray: *mut c_char, /* O - array of flags: 1 if null pixel; else 0 */
    anynul: *mut c_int,    /* O - set to 1 if any values are null; else 0 */
    status: *mut c_int,    /* IO - error status                           */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let anynul = anynul.as_mut();

        let array = slice::from_raw_parts_mut(array, nelem as usize);
        let nularray = slice::from_raw_parts_mut(nularray, nelem as usize);

        ffgcfb_safe(
            fptr, colnum, firstrow, firstelem, nelem, array, nularray, anynul, status,
        )
    }
}

/*--------------------------------------------------------------------------*/
/// Read an array of values from a column in the current FITS HDU. Automatic
/// datatype conversion will be performed if the datatype of the column does not
/// match the datatype of the array parameter. The output values will be scaled
/// by the FITS TSCALn and TZEROn values if these values have been defined.
/// Nularray will be set = 1 if the corresponding array pixel is undefined,
/// otherwise nularray will = 0.
pub fn ffgcfb_safe(
    fptr: &mut fitsfile,        /* I - FITS file pointer                       */
    colnum: c_int,              /* I - number of column to read (1 = 1st col)  */
    firstrow: LONGLONG,         /* I - first row to read (1 = 1st row)         */
    firstelem: LONGLONG,        /* I - first vector element to read (1 = 1st)  */
    nelem: LONGLONG,            /* I - number of values to read                */
    array: &mut [u8],           /* O - array of values that are read       */
    nularray: &mut [c_char],    /* O - array of flags: 1 if null pixel; else 0 */
    anynul: Option<&mut c_int>, /* O - set to 1 if any values are null; else 0 */
    status: &mut c_int,         /* IO - error status                           */
) -> c_int {
    let dummy: u8 = 0;

    ffgclb(
        fptr,
        colnum,
        firstrow,
        firstelem,
        nelem,
        1,
        NullCheckType::SetNullArray,
        dummy,
        array,
        cast_slice_mut(nularray),
        anynul,
        status,
    );
    *status
}

/*--------------------------------------------------------------------------*/
/// Read an array of values from a column in the current FITS HDU.
/// The column number may refer to a real column in an ASCII or binary table,
/// or it may refer be a virtual column in a 1 or more grouped FITS primary
/// array or image extension.  FITSIO treats a primary array as a binary table
/// with 2 vector columns: the first column contains the group parameters (often
/// with length = 0) and the second column contains the array of image pixels.
/// Each row of the table represents a group in the case of multigroup FITS
/// images.
///
/// The output array of values will be converted from the datatype of the column
/// and will be scaled by the FITS TSCALn and TZEROn values if necessary.
pub(crate) fn ffgclb(
    fptr: &mut fitsfile,   /* I - FITS file pointer                       */
    colnum: c_int,         /* I - number of column to read (1 = 1st col)  */
    firstrow: LONGLONG,    /* I - first row to read (1 = 1st row)         */
    firstelem: LONGLONG,   /* I - first vector element to read (1 = 1st)  */
    nelem: LONGLONG,       /* I - number of values to read                */
    elemincre: c_long,     /* I - pixel increment; e.g., 2 = every other  */
    nultyp: NullCheckType, /* I - null value handling code:               */
    /*     1: set undefined pixels = nulval        */
    /*     2: set nularray=1 for undefined pixels  */
    nulval: u8,                     /* I - value for null pixels if nultyp = 1 */
    array: &mut [u8],               /* O - array of values that are read       */
    nularray: &mut [c_char],        /* O - array of flags = 1 if nultyp = 2        */
    mut anynul: Option<&mut c_int>, /* O - set to 1 if any values are null; else 0 */
    status: &mut c_int,             /* IO - error status                           */
) -> c_int {
    let mut scale: f64 = 0.0;
    let mut zero: f64 = 0.0;
    let mut power: f64 = 1.0;
    let mut dtemp: f64 = 0.0;
    let mut tcode: c_int = 0;
    let mut maxelem2: c_int = 0;
    let mut hdutype: c_int = 0;
    let mut xcode: c_int = 0;
    let mut decimals: c_int = 0;
    let mut twidth: c_long = 0;
    let mut incre: c_long = 0;
    let mut ntodo: c_long = 0;
    let mut xwidth: c_long = 0;
    let mut convert: bool = false;
    let mut nulcheck = NullCheckType::None;
    let mut readcheck: c_int = 16; /* see note below on readcheck */
    let mut repeat: LONGLONG = 0;
    let mut startpos: LONGLONG = 0;
    let mut elemnum: LONGLONG = 0;
    let mut readptr: LONGLONG = 0;
    let mut tnull: LONGLONG = 0;
    let mut rowlen: LONGLONG = 0;
    let mut rownum: LONGLONG = 0;
    let mut remain: LONGLONG = 0;
    let mut next: usize = 0;
    let mut rowincre: LONGLONG = 0;
    let mut maxelem: LONGLONG = 0;
    let mut tform: [c_char; 20] = [0; 20];
    let mut message: [c_char; FLEN_ERRMSG] = [0; FLEN_ERRMSG];
    let mut snull: [c_char; 20] = [0; 20]; /*  the FITS null value if reading from ASCII table  */
    let mut buffer: [f64; DBUFFSIZE as usize / mem::size_of::<f64>()] =
        [0.0; DBUFFSIZE as usize / mem::size_of::<f64>()]; /* align cbuff on word boundary */

    let mut u: c_char = 0;

    if *status > 0 || nelem == 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    if let Some(anynul) = anynul.as_deref_mut() {
        *anynul = 0;
    }

    if nultyp == NullCheckType::SetNullArray {
        nularray[..nelem as usize].fill(0); /* initialize nullarray */
    }

    /*---------------------------------------------------*/
    /*  Check input and get parameters about the column: */
    /*---------------------------------------------------*/
    if elemincre < 0 {
        readcheck -= 1; /* don't do range checking in this case */
    }

    /* IMPORTANT NOTE: that the special case of using this subroutine
    to read bytes from a character column are handled internally
    by the call to ffgcprll() below.  It will adjust the effective
    *tcode, repeats, etc, to appear as a TBYTE column. */

    /* Note that readcheck = 16 is equivalent to readcheck = 0
    and readcheck = 15 is equivalent to readcheck = -1,
    but either of those settings allow TSTRINGS to be
    treated as TBYTE vectors, but with full error checking */

    ffgcprll(
        fptr,
        colnum,
        firstrow,
        firstelem,
        nelem,
        readcheck,
        &mut scale,
        &mut zero,
        &mut tform,
        &mut twidth,
        &mut tcode,
        &mut maxelem2,
        &mut startpos,
        &mut elemnum,
        &mut incre,
        &mut repeat,
        &mut rowlen,
        &mut hdutype,
        &mut tnull,
        &mut snull,
        status,
    );
    maxelem = maxelem2 as LONGLONG;

    /* special case */
    if tcode == TLOGICAL && elemincre == 1 {
        u = nulval as c_char;
        ffgcll(
            fptr,
            colnum,
            firstrow as LONGLONG,
            firstelem as LONGLONG,
            nelem as LONGLONG,
            nultyp,
            u,
            cast_slice_mut(array),
            nularray,
            anynul.as_deref_mut(),
            status,
        );

        return *status;
    }

    if *status > 0 {
        return *status;
    }

    incre *= elemincre; /* multiply incre to just get every nth pixel */

    if tcode == TSTRING && hdutype == ASCII_TBL
    /* setup for ASCII tables */
    {
        /* get the number of implied decimal places if no explicit decmal point */
        ffasfm_safe(
            &tform,
            Some(&mut xcode),
            Some(&mut xwidth),
            Some(&mut decimals),
            status,
        );
        for ii in 0..(decimals as usize) {
            power *= 10.0;
        }
    }
    /*------------------------------------------------------------------*/
    /*  Decide whether to check for null values in the input FITS file: */
    /*------------------------------------------------------------------*/
    nulcheck = nultyp; /* by default, check for null values in the FITS file */

    if nultyp == NullCheckType::SetPixel && nulval == 0 {
        nulcheck = NullCheckType::None; /* calling routine does not want to check for nulls */
    } else if tcode % 10 == 1 && tnull == NULL_UNDEFINED as LONGLONG {
        /* if reading an integer column, and  */
        /* if a null value is not defined,    */

        nulcheck = NullCheckType::None; /* then do not check for null values. */
    } else if tcode == TSHORT
        && (tnull > c_short::MAX as LONGLONG || tnull < c_short::MIN as LONGLONG)
    {
        nulcheck = NullCheckType::None; /* Impossible null value */
    } else if tcode == TBYTE && (tnull > 255 || tnull < 0) {
        nulcheck = NullCheckType::None; /* Impossible null value */
    } else if tcode == TSTRING && snull[0] == ASCII_NULL_UNDEFINED {
        nulcheck = NullCheckType::None;
    }

    /*----------------------------------------------------------------------*/
    /*  If FITS column and output data array have same datatype, then we do */
    /*  not need to use a temporary buffer to store intermediate datatype.  */
    /*----------------------------------------------------------------------*/
    convert = true;
    if tcode == TBYTE {
        /* Special Case:                        */
        /* no type convertion required, so read */
        /* data directly into output buffer.    */

        if nelem < INT32_MAX as LONGLONG {
            maxelem = nelem;
        } else {
            maxelem = INT32_MAX as LONGLONG;
        }

        if nulcheck == NullCheckType::None && scale == 1.0 && zero == 0.0 {
            convert = false; /* no need to scale data or find nulls */
        }
    }

    /*---------------------------------------------------------------------*/
    /*  Now read the pixels from the FITS column. If the column does not   */
    /*  have the same datatype as the output array, then we have to read   */
    /*  the raw values into a temporary buffer (of limited size).  In      */
    /*  the case of a vector colum read only 1 vector of values at a time  */
    /*  then skip to the next row if more values need to be read.          */
    /*  After reading the raw values, then call the fffXXYY routine to (1) */
    /*  test for undefined values, (2) convert the datatype if necessary,  */
    /*  and (3) scale the values by the FITS TSCALn and TZEROn linear      */
    /*  scaling parameters.                                                */
    /*---------------------------------------------------------------------*/
    remain = nelem; /* remaining number of values to read */
    next = 0; /* next element in array to be read   */
    rownum = 0; /* row number, relative to firstrow   */

    while remain != 0 {
        /* limit the number of pixels to read at one time to the number that
           will fit in the buffer or to the number of pixels that remain in
           the current vector, which ever is smaller.
        */
        ntodo = cmp::min(remain, maxelem) as c_long;
        if elemincre >= 0 {
            ntodo = cmp::min(
                ntodo as LONGLONG,
                (repeat - elemnum - 1) / (elemincre as LONGLONG) + 1,
            ) as c_long;
        } else {
            ntodo = cmp::min(ntodo as LONGLONG, elemnum / (-elemincre as LONGLONG) + 1) as c_long;
        }

        readptr = startpos
            + ((rownum as LONGLONG) * rowlen)
            + (elemnum * (incre / elemincre) as LONGLONG);

        match tcode {
            TBYTE => {
                /*
                // WARNING: In 'C' version, ffgi1b reads directly into array[next], so null values will be in the output as per the input
                // E.g. input [0, 99], nulltyp = 2
                // Expected array = output = [0, 0] because buffer is zero'd and the null of 99 is left out, so 0 in place
                // In C version, input is copied to output here, and then processed, so array = output= [0, 99]
                // Documentation not clear what should happen.

                // Below is the appropriate safe rust version with what i believe is the intended logic
                if convert {
                    ffgi1b(
                        fptr,
                        readptr,
                        ntodo,
                        incre,
                        cast_slice_mut(&mut buffer),
                        status,
                    );

                    fffi1i1(
                        cast_slice(&buffer),
                        ntodo,
                        scale,
                        zero,
                        nulcheck,
                        tnull as u8,
                        nulval,
                        &mut nularray[next..],
                        anynul.as_deref_mut(),
                        &mut array[next..],
                        status,
                    );
                } else {
                    ffgi1b(fptr, readptr, ntodo, incre, &mut array[next..], status);
                }
                */

                // The original C logic
                ffgi1b(fptr, readptr, ntodo, incre, &mut array[next..], status);

                if convert {
                    fffi1i1_inplace(
                        cast_slice_mut(&mut array[next..]),
                        ntodo,
                        scale,
                        zero,
                        nulcheck,
                        tnull as u8,
                        nulval,
                        &mut nularray[next..],
                        anynul.as_deref_mut(),
                        status,
                    );
                }
            }
            TSHORT => {
                ffgi2b(
                    fptr,
                    readptr,
                    ntodo,
                    incre,
                    cast_slice_mut(&mut buffer),
                    status,
                );
                fffi2i1(
                    cast_slice(&buffer),
                    ntodo,
                    scale,
                    zero,
                    nulcheck,
                    tnull as c_short,
                    nulval,
                    &mut nularray[next..],
                    anynul.as_deref_mut(),
                    &mut array[next..],
                    status,
                );
            }
            TLONG => {
                ffgi4b(
                    fptr,
                    readptr,
                    ntodo,
                    incre,
                    cast_slice_mut(&mut buffer),
                    status,
                );
                fffi4i1(
                    cast_slice(&buffer),
                    ntodo,
                    scale,
                    zero,
                    nulcheck,
                    tnull as INT32BIT,
                    nulval,
                    &mut nularray[next..],
                    anynul.as_deref_mut(),
                    &mut array[next..],
                    status,
                );
            }
            TLONGLONG => {
                ffgi8b(
                    fptr,
                    readptr,
                    ntodo,
                    incre,
                    cast_slice_mut(&mut buffer),
                    status,
                );
                fffi8i1(
                    cast_slice(&buffer),
                    ntodo,
                    scale,
                    zero,
                    nulcheck,
                    tnull,
                    nulval,
                    &mut nularray[next..],
                    anynul.as_deref_mut(),
                    &mut array[next..],
                    status,
                );
            }
            TFLOAT => {
                ffgr4b(
                    fptr,
                    readptr,
                    ntodo,
                    incre,
                    cast_slice_mut(&mut buffer),
                    status,
                );
                fffr4i1(
                    cast_slice(&buffer),
                    ntodo,
                    scale,
                    zero,
                    nulcheck,
                    nulval,
                    &mut nularray[next..],
                    anynul.as_deref_mut(),
                    &mut array[next..],
                    status,
                );
            }
            TDOUBLE => {
                ffgr8b(
                    fptr,
                    readptr,
                    ntodo,
                    incre,
                    cast_slice_mut(&mut buffer),
                    status,
                );
                fffr8i1(
                    cast_slice(&buffer),
                    ntodo,
                    scale,
                    zero,
                    nulcheck,
                    nulval,
                    &mut nularray[next..],
                    anynul.as_deref_mut(),
                    &mut array[next..],
                    status,
                );
            }
            TSTRING => {
                ffmbyt_safe(fptr, readptr, REPORT_EOF, status);

                if incre == twidth {
                    /* contiguous bytes */
                    ffgbyt(
                        fptr,
                        (ntodo * twidth) as LONGLONG,
                        cast_slice_mut(&mut buffer),
                        status,
                    );
                } else {
                    ffgbytoff(
                        fptr,
                        twidth,
                        ntodo,
                        incre - twidth,
                        cast_slice_mut(&mut buffer),
                        status,
                    );
                }

                /* interpret the string as an ASCII formated number */
                fffstri1(
                    cast_slice_mut(&mut buffer),
                    ntodo,
                    scale,
                    zero,
                    twidth,
                    power,
                    nulcheck,
                    &snull,
                    nulval,
                    &mut nularray[next..],
                    anynul.as_deref_mut(),
                    &mut array[next..],
                    status,
                );
            }
            _ => {
                /*  error trap for invalid column format */
                int_snprintf!(
                    &mut message,
                    FLEN_ERRMSG,
                    "Cannot read bytes from column {} which has format {}",
                    colnum,
                    slice_to_str!(&tform),
                );
                ffpmsg_slice(&message);
                if hdutype == ASCII_TBL {
                    *status = BAD_ATABLE_FORMAT;
                    return *status;
                } else {
                    *status = BAD_BTABLE_FORMAT;
                    return *status;
                }
            }
        } /* End of switch block */

        /*-------------------------*/
        /*  Check for fatal error  */
        /*-------------------------*/
        if *status > 0 {
            /* test for error during previous read operation */

            dtemp = next as f64;
            if hdutype > 0 {
                int_snprintf!(
                    &mut message,
                    FLEN_ERRMSG,
                    "Error reading elements {:.0} thru {:.0} from column {} (ffgclb).",
                    dtemp + 1.0,
                    dtemp + ntodo as f64,
                    colnum,
                );
            } else {
                int_snprintf!(
                    &mut message,
                    FLEN_ERRMSG,
                    "Error reading elements {:.0} thru {:.0} from image (ffgclb).",
                    dtemp + 1.0,
                    dtemp + ntodo as f64,
                );
            }

            ffpmsg_slice(&message);
            return *status;
        }

        /*--------------------------------------------*/
        /*  increment the counters for the next loop  */
        /*--------------------------------------------*/
        remain -= ntodo as LONGLONG;
        if remain != 0 {
            next += ntodo as usize;
            elemnum += (ntodo * elemincre) as LONGLONG;

            if elemnum >= repeat
            /* completed a row; start on later row */
            {
                rowincre = elemnum / repeat;
                rownum += rowincre;
                elemnum -= rowincre * repeat;
            } else if elemnum < 0
            /* completed a row; start on a previous row */
            {
                rowincre = (-elemnum - 1) / repeat + 1;
                rownum -= rowincre;
                elemnum += rowincre * repeat;
            }
        }
    } /*  End of main while Loop  */

    /*--------------------------------*/
    /*  check for numerical overflow  */
    /*--------------------------------*/
    if *status == OVERFLOW_ERR {
        ffpmsg_str("Numerical overflow during type conversion while reading FITS data.");
        *status = NUM_OVERFLOW;
    }

    *status
}

/*--------------------------------------------------------------------------*/
/// Read a stream of bytes from the current FITS HDU.  This primative routine is mainly
/// for reading non-standard "conforming" extensions and should not be used
/// for standard IMAGE, TABLE or BINTABLE extensions.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgextn(
    fptr: *mut fitsfile, /* I - FITS file pointer                        */
    offset: LONGLONG,    /* I - byte offset from start of extension data */
    nelem: LONGLONG,     /* I - number of elements to read               */
    buffer: *mut c_void, /* I - stream of bytes to read                  */
    status: *mut c_int,  /* IO - error status                            */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let buffer = slice::from_raw_parts_mut(buffer as *mut u8, nelem as usize);

        ffgextn_safe(fptr, offset, nelem, buffer, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Read a stream of bytes from the current FITS HDU.  This primative routine is mainly
/// for reading non-standard "conforming" extensions and should not be used
/// for standard IMAGE, TABLE or BINTABLE extensions.
pub fn ffgextn_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer                        */
    offset: LONGLONG,    /* I - byte offset from start of extension data */
    nelem: LONGLONG,     /* I - number of elements to read               */
    buffer: &mut [u8],   /* I - stream of bytes to read                  */
    status: &mut c_int,  /* IO - error status                            */
) -> c_int {
    if *status > 0 {
        /* inherit input status value if > 0 */
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
    /* move to write position */
    ffmbyt_safe(fptr, fptr.Fptr.datastart + offset, IGNORE_EOF, status);

    /* read the buffer */
    ffgbyt(fptr, nelem, buffer, status);

    *status
}

/*--------------------------------------------------------------------------*/
/// Copy input to output following reading of the input from a FITS file.
/// Check for null values and do datatype conversion and scaling if required.
/// The nullcheck code value determines how any null values in the input array
/// are treated.  A null value is an input pixel that is equal to tnull.  If
/// nullcheck = 0, then no checking for nulls is performed and any null values
/// will be transformed just like any other pixel.  If nullcheck = 1, then the
/// output pixel will be set = nullval if the corresponding input pixel is null.
/// If nullcheck = 2, then if the pixel is null then the corresponding value of
/// nullarray will be set to 1; the value of nullarray for non-null pixels
/// will = 0.  The anynul parameter will be set = 1 if any of the returned
/// pixels are null, otherwise anynul will be returned with a value = 0;
pub(crate) fn fffi1i1(
    input: &[u8],             /* I - array of values to be converted     */
    ntodo: c_long,            /* I - number of elements in the array     */
    scale: f64,               /* I - FITS TSCALn or BSCALE value         */
    zero: f64,                /* I - FITS TZEROn or BZERO  value         */
    nullcheck: NullCheckType, /* I - null checking code; 0 = don't check */
    /*     1:set null pixels = nullval         */
    /*     2: if null pixel, set nullarray = 1 */
    tnull: u8,                  /* I - value of FITS TNULLn keyword if any */
    nullval: u8,                /* I - set null pixels, if nullcheck = 1   */
    nullarray: &mut [c_char],   /* O - bad pixel array, if nullcheck = 2   */
    anynul: Option<&mut c_int>, /* O - set to 1 if any pixels are null     */
    output: &mut [u8],          /* O - array of converted pixels           */
    status: &mut c_int,         /* IO - error status                       */
) -> c_int {
    let mut dvalue: f64 = 0.0;

    if nullcheck == NullCheckType::None {
        /* no null checking required */

        if scale == 1.0 && zero == 0.0 {
            /* no scaling */

            /* this routine is normally not called in this case */
            output.copy_from_slice(input);
        } else {
            /* must scale the data */
            for ii in 0..(ntodo as usize) {
                dvalue = (input[ii] as f64) * scale + zero;

                if dvalue < DUCHAR_MIN {
                    *status = OVERFLOW_ERR;
                    output[ii] = 0;
                } else if dvalue > DUCHAR_MAX {
                    *status = OVERFLOW_ERR;
                    output[ii] = u8::MAX;
                } else {
                    output[ii] = dvalue as u8;
                }
            }
        }
    } else {
        /* must check for null values */
        let anynul = anynul.unwrap();

        if scale == 1.0 && zero == 0.0 {
            /* no scaling */

            for ii in 0..(ntodo as usize) {
                if input[ii] == tnull {
                    *anynul = 1;
                    if nullcheck == NullCheckType::SetPixel {
                        output[ii] = nullval;
                    } else {
                        nullarray[ii] = 1;
                    }
                } else {
                    output[ii] = input[ii];
                }
            }
        } else {
            /* must scale the data */

            for ii in 0..(ntodo as usize) {
                if input[ii] == tnull {
                    *anynul = 1;
                    if nullcheck == NullCheckType::SetPixel {
                        output[ii] = nullval;
                    } else {
                        nullarray[ii] = 1;
                    }
                } else {
                    dvalue = (input[ii] as f64) * scale + zero;

                    if dvalue < DUCHAR_MIN {
                        *status = OVERFLOW_ERR;
                        output[ii] = 0;
                    } else if dvalue > DUCHAR_MAX {
                        *status = OVERFLOW_ERR;
                        output[ii] = u8::MAX;
                    } else {
                        output[ii] = dvalue as u8;
                    }
                }
            }
        }
    }
    *status
}

/*--------------------------------------------------------------------------*/
/// Copy input to output following reading of the input from a FITS file.
/// Check for null values and do datatype conversion and scaling if required.
/// The nullcheck code value determines how any null values in the input array
/// are treated.  A null value is an input pixel that is equal to tnull.  If
/// nullcheck = 0, then no checking for nulls is performed and any null values
/// will be transformed just like any other pixel.  If nullcheck = 1, then the
/// output pixel will be set = nullval if the corresponding input pixel is null.
/// If nullcheck = 2, then if the pixel is null then the corresponding value of
/// nullarray will be set to 1; the value of nullarray for non-null pixels
/// will = 0.  The anynul parameter will be set = 1 if any of the returned
/// pixels are null, otherwise anynul will be returned with a value = 0;
pub(crate) fn fffi1i1_inplace(
    inout: &mut [u8],         /* I - array of values to be converted     */
    ntodo: c_long,            /* I - number of elements in the array     */
    scale: f64,               /* I - FITS TSCALn or BSCALE value         */
    zero: f64,                /* I - FITS TZEROn or BZERO  value         */
    nullcheck: NullCheckType, /* I - null checking code; 0 = don't check */
    /*     1:set null pixels = nullval         */
    /*     2: if null pixel, set nullarray = 1 */
    tnull: u8,                  /* I - value of FITS TNULLn keyword if any */
    nullval: u8,                /* I - set null pixels, if nullcheck = 1   */
    nullarray: &mut [c_char],   /* O - bad pixel array, if nullcheck = 2   */
    anynul: Option<&mut c_int>, /* O - set to 1 if any pixels are null     */
    status: &mut c_int,         /* IO - error status                       */
) -> c_int {
    let mut dvalue: f64 = 0.0;

    if nullcheck == NullCheckType::None {
        /* no null checking required */

        if scale == 1.0 && zero == 0.0 {
            /* no scaling */
            /* this routine is normally not called in this case */
            // no op
        } else {
            /* must scale the data */
            for ii in 0..(ntodo as usize) {
                dvalue = (inout[ii] as f64) * scale + zero;

                if dvalue < DUCHAR_MIN {
                    *status = OVERFLOW_ERR;
                    inout[ii] = 0;
                } else if dvalue > DUCHAR_MAX {
                    *status = OVERFLOW_ERR;
                    inout[ii] = u8::MAX;
                } else {
                    inout[ii] = dvalue as u8;
                }
            }
        }
    } else {
        /* must check for null values */
        let anynul = anynul.unwrap();

        if scale == 1.0 && zero == 0.0 {
            /* no scaling */

            for ii in 0..(ntodo as usize) {
                if inout[ii] == tnull {
                    *anynul = 1;
                    if nullcheck == NullCheckType::SetPixel {
                        inout[ii] = nullval;
                    } else {
                        nullarray[ii] = 1;
                    }
                }
            }
        } else {
            /* must scale the data */

            for ii in 0..(ntodo as usize) {
                if inout[ii] == tnull {
                    *anynul = 1;
                    if nullcheck == NullCheckType::SetPixel {
                        inout[ii] = nullval;
                    } else {
                        nullarray[ii] = 1;
                    }
                } else {
                    dvalue = (inout[ii] as f64) * scale + zero;

                    if dvalue < DUCHAR_MIN {
                        *status = OVERFLOW_ERR;
                        inout[ii] = 0;
                    } else if dvalue > DUCHAR_MAX {
                        *status = OVERFLOW_ERR;
                        inout[ii] = u8::MAX;
                    } else {
                        inout[ii] = dvalue as u8;
                    }
                }
            }
        }
    }
    *status
}

/*--------------------------------------------------------------------------*/
/// Copy input to output following reading of the input from a FITS file.
/// Check for null values and do datatype conversion and scaling if required.
/// The nullcheck code value determines how any null values in the input array
/// are treated.  A null value is an input pixel that is equal to tnull.  If
/// nullcheck = 0, then no checking for nulls is performed and any null values
/// will be transformed just like any other pixel.  If nullcheck = 1, then the
/// output pixel will be set = nullval if the corresponding input pixel is null.
/// If nullcheck = 2, then if the pixel is null then the corresponding value of
/// nullarray will be set to 1; the value of nullarray for non-null pixels
/// will = 0.  The anynul parameter will be set = 1 if any of the returned
/// pixels are null, otherwise anynul will be returned with a value = 0;
pub(crate) fn fffi2i1(
    input: &[c_short],        /* I - array of values to be converted     */
    ntodo: c_long,            /* I - number of elements in the array     */
    scale: f64,               /* I - FITS TSCALn or BSCALE value         */
    zero: f64,                /* I - FITS TZEROn or BZERO  value         */
    nullcheck: NullCheckType, /* I - null checking code; 0 = don't check */
    /*     1:set null pixels = nullval         */
    /*     2: if null pixel, set nullarray = 1 */
    tnull: c_short,             /* I - value of FITS TNULLn keyword if any */
    nullval: u8,                /* I - set null pixels, if nullcheck = 1   */
    nullarray: &mut [c_char],   /* O - bad pixel array, if nullcheck = 2   */
    anynul: Option<&mut c_int>, /* O - set to 1 if any pixels are null     */
    output: &mut [u8],          /* O - array of converted pixels           */
    status: &mut c_int,         /* IO - error status                       */
) -> c_int {
    let mut dvalue: f64 = 0.0;

    if nullcheck == NullCheckType::None {
        /* no null checking required */
        if scale == 1.0 && zero == 0.0 {
            /* no scaling */

            for ii in 0..(ntodo as usize) {
                if input[ii] < 0 {
                    *status = OVERFLOW_ERR;
                    output[ii] = 0;
                } else if input[ii] > u8::MAX as c_short {
                    *status = OVERFLOW_ERR;
                    output[ii] = u8::MAX;
                } else {
                    output[ii] = input[ii] as u8;
                }
            }
        } else {
            /* must scale the data */

            for ii in 0..(ntodo as usize) {
                dvalue = (input[ii] as f64) * scale + zero;

                if dvalue < DUCHAR_MIN {
                    *status = OVERFLOW_ERR;
                    output[ii] = 0;
                } else if dvalue > DUCHAR_MAX {
                    *status = OVERFLOW_ERR;
                    output[ii] = u8::MAX;
                } else {
                    output[ii] = dvalue as u8;
                }
            }
        }
    } else {
        /* must check for null values */
        let anynul = anynul.unwrap();

        if scale == 1.0 && zero == 0.0 {
            /* no scaling */

            for ii in 0..(ntodo as usize) {
                if input[ii] == tnull {
                    *anynul = 1;
                    if nullcheck == NullCheckType::SetPixel {
                        output[ii] = nullval;
                    } else {
                        nullarray[ii] = 1;
                    }
                } else if input[ii] < 0 {
                    *status = OVERFLOW_ERR;
                    output[ii] = 0;
                } else if input[ii] > u8::MAX as c_short {
                    *status = OVERFLOW_ERR;
                    output[ii] = u8::MAX;
                } else {
                    output[ii] = input[ii] as u8;
                }
            }
        } else {
            /* must scale the data */

            for ii in 0..(ntodo as usize) {
                if input[ii] == tnull {
                    *anynul = 1;
                    if nullcheck == NullCheckType::SetPixel {
                        output[ii] = nullval;
                    } else {
                        nullarray[ii] = 1;
                    }
                } else {
                    dvalue = (input[ii] as f64) * scale + zero;

                    if dvalue < DUCHAR_MIN {
                        *status = OVERFLOW_ERR;
                        output[ii] = 0;
                    } else if dvalue > DUCHAR_MAX {
                        *status = OVERFLOW_ERR;
                        output[ii] = u8::MAX;
                    } else {
                        output[ii] = dvalue as u8;
                    }
                }
            }
        }
    }
    *status
}

/*--------------------------------------------------------------------------*/
/// Copy input to output following reading of the input from a FITS file.
/// Check for null values and do datatype conversion and scaling if required.
/// The nullcheck code value determines how any null values in the input array
/// are treated.  A null value is an input pixel that is equal to tnull.  If
/// nullcheck = 0, then no checking for nulls is performed and any null values
/// will be transformed just like any other pixel.  If nullcheck = 1, then the
/// output pixel will be set = nullval if the corresponding input pixel is null.
/// If nullcheck = 2, then if the pixel is null then the corresponding value of
/// nullarray will be set to 1; the value of nullarray for non-null pixels
/// will = 0.  The anynul parameter will be set = 1 if any of the returned
/// pixels are null, otherwise anynul will be returned with a value = 0;
pub(crate) fn fffi4i1(
    input: &[INT32BIT],       /* I - array of values to be converted     */
    ntodo: c_long,            /* I - number of elements in the array     */
    scale: f64,               /* I - FITS TSCALn or BSCALE value         */
    zero: f64,                /* I - FITS TZEROn or BZERO  value         */
    nullcheck: NullCheckType, /* I - null checking code; 0 = don't check */
    /*     1:set null pixels = nullval         */
    /*     2: if null pixel, set nullarray = 1 */
    tnull: INT32BIT,            /* I - value of FITS TNULLn keyword if any */
    nullval: u8,                /* I - set null pixels, if nullcheck = 1   */
    nullarray: &mut [c_char],   /* O - bad pixel array, if nullcheck = 2   */
    anynul: Option<&mut c_int>, /* O - set to 1 if any pixels are null     */
    output: &mut [u8],          /* O - array of converted pixels           */
    status: &mut c_int,         /* IO - error status                       */
) -> c_int {
    let mut dvalue: f64 = 0.0;

    if nullcheck == NullCheckType::None {
        /* no null checking required */

        if scale == 1.0 && zero == 0.0 {
            /* no scaling */

            for ii in 0..(ntodo as usize) {
                if input[ii] < 0 {
                    *status = OVERFLOW_ERR;
                    output[ii] = 0;
                } else if input[ii] > u8::MAX as INT32BIT {
                    *status = OVERFLOW_ERR;
                    output[ii] = u8::MAX;
                } else {
                    output[ii] = input[ii] as u8;
                }
            }
        } else {
            /* must scale the data */

            for ii in 0..(ntodo as usize) {
                dvalue = (input[ii] as f64) * scale + zero;

                if dvalue < DUCHAR_MIN {
                    *status = OVERFLOW_ERR;
                    output[ii] = 0;
                } else if dvalue > DUCHAR_MAX {
                    *status = OVERFLOW_ERR;
                    output[ii] = u8::MAX;
                } else {
                    output[ii] = dvalue as u8;
                }
            }
        }
    } else {
        /* must check for null values */
        let anynul = anynul.unwrap();

        if scale == 1.0 && zero == 0.0 {
            /* no scaling */

            for ii in 0..(ntodo as usize) {
                if input[ii] == tnull {
                    *anynul = 1;
                    if nullcheck == NullCheckType::SetPixel {
                        output[ii] = nullval;
                    } else {
                        nullarray[ii] = 1;
                    }
                } else if input[ii] < 0 {
                    *status = OVERFLOW_ERR;
                    output[ii] = 0;
                } else if input[ii] > u8::MAX as INT32BIT {
                    *status = OVERFLOW_ERR;
                    output[ii] = u8::MAX;
                } else {
                    output[ii] = input[ii] as u8;
                }
            }
        } else {
            /* must scale the data */

            for ii in 0..(ntodo as usize) {
                if input[ii] == tnull {
                    *anynul = 1;
                    if nullcheck == NullCheckType::SetPixel {
                        output[ii] = nullval;
                    } else {
                        nullarray[ii] = 1;
                    }
                } else {
                    dvalue = (input[ii] as f64) * scale + zero;

                    if dvalue < DUCHAR_MIN {
                        *status = OVERFLOW_ERR;
                        output[ii] = 0;
                    } else if dvalue > DUCHAR_MAX {
                        *status = OVERFLOW_ERR;
                        output[ii] = u8::MAX;
                    } else {
                        output[ii] = dvalue as u8;
                    }
                }
            }
        }
    }
    *status
}

/*--------------------------------------------------------------------------*/
/// Copy input to output following reading of the input from a FITS file.
/// Check for null values and do datatype conversion and scaling if required.
/// The nullcheck code value determines how any null values in the input array
/// are treated.  A null value is an input pixel that is equal to tnull.  If
/// nullcheck = 0, then no checking for nulls is performed and any null values
/// will be transformed just like any other pixel.  If nullcheck = 1, then the
/// output pixel will be set = nullval if the corresponding input pixel is null.
/// If nullcheck = 2, then if the pixel is null then the corresponding value of
/// nullarray will be set to 1; the value of nullarray for non-null pixels
/// will = 0.  The anynul parameter will be set = 1 if any of the returned
/// pixels are null, otherwise anynul will be returned with a value = 0;
pub(crate) fn fffi8i1(
    input: &[LONGLONG],       /* I - array of values to be converted     */
    ntodo: c_long,            /* I - number of elements in the array     */
    scale: f64,               /* I - FITS TSCALn or BSCALE value         */
    zero: f64,                /* I - FITS TZEROn or BZERO  value         */
    nullcheck: NullCheckType, /* I - null checking code; 0 = don't check */
    /*     1:set null pixels = nullval         */
    /*     2: if null pixel, set nullarray = 1 */
    tnull: LONGLONG,            /* I - value of FITS TNULLn keyword if any */
    nullval: u8,                /* I - set null pixels, if nullcheck = 1   */
    nullarray: &mut [c_char],   /* O - bad pixel array, if nullcheck = 2   */
    anynul: Option<&mut c_int>, /* O - set to 1 if any pixels are null     */
    output: &mut [u8],          /* O - array of converted pixels           */
    status: &mut c_int,         /* IO - error status                       */
) -> c_int {
    let mut dvalue: f64 = 0.0;
    let mut ulltemp: ULONGLONG = 0;

    if nullcheck == NullCheckType::None {
        /* no null checking required */

        if scale == 1.0 && zero == 9223372036854775808. {
            /* The column we read contains unsigned long long values. */
            /* Instead of adding 9223372036854775808, it is more efficient */
            /* and more precise to just flip the sign bit with the XOR operator */

            for ii in 0..(ntodo as usize) {
                ulltemp = (input[ii] as ULONGLONG) ^ 0x8000000000000000;

                if ulltemp > u8::MAX as ULONGLONG {
                    *status = OVERFLOW_ERR;
                    output[ii] = u8::MAX;
                } else {
                    output[ii] = ulltemp as u8;
                }
            }
        } else if scale == 1.0 && zero == 0.0 {
            /* no scaling */

            for ii in 0..(ntodo as usize) {
                if input[ii] < 0 {
                    *status = OVERFLOW_ERR;
                    output[ii] = 0;
                } else if input[ii] > u8::MAX as LONGLONG {
                    *status = OVERFLOW_ERR;
                    output[ii] = u8::MAX;
                } else {
                    output[ii] = input[ii] as u8;
                }
            }
        } else {
            /* must scale the data */

            for ii in 0..(ntodo as usize) {
                dvalue = (input[ii] as f64) * scale + zero;

                if dvalue < DUCHAR_MIN {
                    *status = OVERFLOW_ERR;
                    output[ii] = 0;
                } else if dvalue > DUCHAR_MAX {
                    *status = OVERFLOW_ERR;
                    output[ii] = u8::MAX;
                } else {
                    output[ii] = dvalue as u8;
                }
            }
        }
    } else {
        /* must check for null values */
        let anynul = anynul.unwrap();

        if scale == 1.0 && zero == 9223372036854775808. {
            /* The column we read contains unsigned long long values. */
            /* Instead of adding 9223372036854775808, it is more efficient */
            /* and more precise to just flip the sign bit with the XOR operator */

            for ii in 0..(ntodo as usize) {
                if input[ii] == tnull {
                    *anynul = 1;
                    if nullcheck == NullCheckType::SetPixel {
                        output[ii] = nullval;
                    } else {
                        nullarray[ii] = 1;
                    }
                } else {
                    ulltemp = (input[ii] as ULONGLONG) ^ 0x8000000000000000;

                    if ulltemp > u8::MAX as ULONGLONG {
                        *status = OVERFLOW_ERR;
                        output[ii] = u8::MAX;
                    } else {
                        output[ii] = ulltemp as u8;
                    }
                }
            }
        } else if scale == 1.0 && zero == 0.0 {
            /* no scaling */

            for ii in 0..(ntodo as usize) {
                if input[ii] == tnull {
                    *anynul = 1;
                    if nullcheck == NullCheckType::SetPixel {
                        output[ii] = nullval;
                    } else {
                        nullarray[ii] = 1;
                    }
                } else if input[ii] < 0 {
                    *status = OVERFLOW_ERR;
                    output[ii] = 0;
                } else if input[ii] > u8::MAX as LONGLONG {
                    *status = OVERFLOW_ERR;
                    output[ii] = u8::MAX;
                } else {
                    output[ii] = input[ii] as u8;
                }
            }
        } else {
            /* must scale the data */

            for ii in 0..(ntodo as usize) {
                if input[ii] == tnull {
                    *anynul = 1;
                    if nullcheck == NullCheckType::SetPixel {
                        output[ii] = nullval;
                    } else {
                        nullarray[ii] = 1;
                    }
                } else {
                    dvalue = (input[ii] as f64) * scale + zero;

                    if dvalue < DUCHAR_MIN {
                        *status = OVERFLOW_ERR;
                        output[ii] = 0;
                    } else if dvalue > DUCHAR_MAX {
                        *status = OVERFLOW_ERR;
                        output[ii] = u8::MAX;
                    } else {
                        output[ii] = dvalue as u8;
                    }
                }
            }
        }
    }
    *status
}

/*--------------------------------------------------------------------------*/
/// Copy input to output following reading of the input from a FITS file.
/// Check for null values and do datatype conversion and scaling if required.
/// The nullcheck code value determines how any null values in the input array
/// are treated.  A null value is an input pixel that is equal to NaN.  If
/// nullcheck = 0, then no checking for nulls is performed and any null values
/// will be transformed just like any other pixel.  If nullcheck = 1, then the
/// output pixel will be set = nullval if the corresponding input pixel is null.
/// If nullcheck = 2, then if the pixel is null then the corresponding value of
/// nullarray will be set to 1; the value of nullarray for non-null pixels
/// will = 0.  The anynul parameter will be set = 1 if any of the returned
/// pixels are null, otherwise anynul will be returned with a value = 0;
pub(crate) fn fffr4i1(
    input: &[f32],            /* I - array of values to be converted     */
    ntodo: c_long,            /* I - number of elements in the array     */
    scale: f64,               /* I - FITS TSCALn or BSCALE value         */
    zero: f64,                /* I - FITS TZEROn or BZERO  value         */
    nullcheck: NullCheckType, /* I - null checking code; 0 = don't check */
    /*     1:set null pixels = nullval         */
    /*     2: if null pixel, set nullarray = 1 */
    nullval: u8,                /* I - set null pixels, if nullcheck = 1   */
    nullarray: &mut [c_char],   /* O - bad pixel array, if nullcheck = 2   */
    anynul: Option<&mut c_int>, /* O - set to 1 if any pixels are null     */
    output: &mut [u8],          /* O - array of converted pixels           */
    status: &mut c_int,         /* IO - error status                       */
) -> c_int {
    let mut dvalue: f64 = 0.0;
    let mut sptr = 0;
    let iret = 0;

    if nullcheck == NullCheckType::None {
        /* no null checking required */

        if scale == 1.0 && zero == 0.0 {
            /* no scaling */

            for ii in 0..(ntodo as usize) {
                if (input[ii] as f64) < DUCHAR_MIN {
                    *status = OVERFLOW_ERR;
                    output[ii] = 0;
                } else if (input[ii] as f64) > DUCHAR_MAX {
                    *status = OVERFLOW_ERR;
                    output[ii] = u8::MAX;
                } else {
                    output[ii] = input[ii] as u8;
                }
            }
        } else {
            /* must scale the data */

            for ii in 0..(ntodo as usize) {
                dvalue = (input[ii] as f64) * scale + zero;

                if dvalue < DUCHAR_MIN {
                    *status = OVERFLOW_ERR;
                    output[ii] = 0;
                } else if dvalue > DUCHAR_MAX {
                    *status = OVERFLOW_ERR;
                    output[ii] = u8::MAX;
                } else {
                    output[ii] = dvalue as u8;
                }
            }
        }
    } else {
        /* must check for null values */
        let anynul = anynul.unwrap();

        if BYTESWAPPED && CFITSIO_MACHINE != VAXVMS && CFITSIO_MACHINE != ALPHAVMS {
            sptr += 1; /* point to MSBs */
        }

        let shortBuffer: &[c_short] = cast_slice(input);

        if scale == 1.0 && zero == 0.0 {
            /* no scaling */

            let ii = 0;
            for ii in 0..(ntodo as usize) {
                let iret = fnan(shortBuffer[sptr]);
                if 0 != iret {
                    /* test for NaN or underflow */

                    if iret == 1 {
                        /* is it a NaN? */
                        *anynul = 1;
                        if nullcheck == NullCheckType::SetPixel {
                            output[ii] = nullval;
                        } else {
                            nullarray[ii] = 1;
                        }
                    } else {
                        /* it's an underflow */
                        output[ii] = 0;
                    }
                } else if (input[ii] as f64) < DUCHAR_MIN {
                    *status = OVERFLOW_ERR;
                    output[ii] = 0;
                } else if (input[ii] as f64) > DUCHAR_MAX {
                    *status = OVERFLOW_ERR;
                    output[ii] = u8::MAX;
                } else {
                    output[ii] = input[ii] as u8;
                }

                sptr += 2;
            }
        } else {
            /* must scale the data */

            let ii = 0;
            for ii in 0..(ntodo as usize) {
                let iret = fnan(shortBuffer[sptr]);
                if 0 != iret {
                    /* test for NaN or underflow */
                    if iret == 1 {
                        /* is it a NaN? */

                        *anynul = 1;
                        if nullcheck == NullCheckType::SetPixel {
                            output[ii] = nullval;
                        } else {
                            nullarray[ii] = 1;
                        }
                    } else {
                        /* it's an underflow */

                        if zero < DUCHAR_MIN {
                            *status = OVERFLOW_ERR;
                            output[ii] = 0;
                        } else if zero > DUCHAR_MAX {
                            *status = OVERFLOW_ERR;
                            output[ii] = u8::MAX;
                        } else {
                            output[ii] = zero as u8;
                        }
                    }
                } else {
                    dvalue = (input[ii] as f64) * scale + zero;

                    if dvalue < DUCHAR_MIN {
                        *status = OVERFLOW_ERR;
                        output[ii] = 0;
                    } else if dvalue > DUCHAR_MAX {
                        *status = OVERFLOW_ERR;
                        output[ii] = u8::MAX;
                    } else {
                        output[ii] = dvalue as u8;
                    }
                }

                sptr += 2;
            }
        }
    }
    *status
}

/*--------------------------------------------------------------------------*/
/// Copy input to output following reading of the input from a FITS file.
/// Check for null values and do datatype conversion and scaling if required.
/// The nullcheck code value determines how any null values in the input array
/// are treated.  A null value is an input pixel that is equal to NaN.  If
/// nullcheck = 0, then no checking for nulls is performed and any null values
/// will be transformed just like any other pixel.  If nullcheck = 1, then the
/// output pixel will be set = nullval if the corresponding input pixel is null.
/// If nullcheck = 2, then if the pixel is null then the corresponding value of
/// nullarray will be set to 1; the value of nullarray for non-null pixels
/// will = 0.  The anynul parameter will be set = 1 if any of the returned
/// pixels are null, otherwise anynul will be returned with a value = 0;
pub(crate) fn fffr8i1(
    input: &[f64],            /* I - array of values to be converted     */
    ntodo: c_long,            /* I - number of elements in the array     */
    scale: f64,               /* I - FITS TSCALn or BSCALE value         */
    zero: f64,                /* I - FITS TZEROn or BZERO  value         */
    nullcheck: NullCheckType, /* I - null checking code; 0 = don't check */
    /*     1:set null pixels = nullval         */
    /*     2: if null pixel, set nullarray = 1 */
    nullval: u8,                /* I - set null pixels, if nullcheck = 1   */
    nullarray: &mut [c_char],   /* O - bad pixel array, if nullcheck = 2   */
    anynul: Option<&mut c_int>, /* O - set to 1 if any pixels are null     */
    output: &mut [u8],          /* O - array of converted pixels           */
    status: &mut c_int,         /* IO - error status                       */
) -> c_int {
    let mut dvalue: f64 = 0.0;
    let mut sptr = 0;
    let mut iret = 0;

    if nullcheck == NullCheckType::None {
        /* no null checking required */

        if scale == 1.0 && zero == 0.0 {
            /* no scaling */

            for ii in 0..(ntodo as usize) {
                if input[ii] < DUCHAR_MIN {
                    *status = OVERFLOW_ERR;
                    output[ii] = 0;
                } else if input[ii] > DUCHAR_MAX {
                    *status = OVERFLOW_ERR;
                    output[ii] = u8::MAX;
                } else {
                    output[ii] = input[ii] as u8;
                }
            }
        } else {
            /* must scale the data */

            for ii in 0..(ntodo as usize) {
                dvalue = input[ii] * scale + zero;

                if dvalue < DUCHAR_MIN {
                    *status = OVERFLOW_ERR;
                    output[ii] = 0;
                } else if dvalue > DUCHAR_MAX {
                    *status = OVERFLOW_ERR;
                    output[ii] = u8::MAX;
                } else {
                    output[ii] = dvalue as u8;
                }
            }
        }
    } else {
        /* must check for null values */
        let anynul = anynul.unwrap();

        let shortBuffer: &[c_short] = cast_slice(input);

        if BYTESWAPPED && CFITSIO_MACHINE != VAXVMS && CFITSIO_MACHINE != ALPHAVMS {
            sptr += 3; /* point to MSBs */
        }

        if scale == 1.0 && zero == 0.0 {
            /* no scaling */

            for ii in 0..(ntodo as usize) {
                iret = dnan(shortBuffer[sptr]);
                if 0 != iret {
                    /* test for NaN or underflow */

                    if iret == 1 {
                        /* is it a NaN? */

                        *anynul = 1;
                        if nullcheck == NullCheckType::SetPixel {
                            output[ii] = nullval;
                        } else {
                            nullarray[ii] = 1;
                        }
                    } else {
                        /* it's an underflow */
                        output[ii] = 0;
                    }
                } else if input[ii] < DUCHAR_MIN {
                    *status = OVERFLOW_ERR;
                    output[ii] = 0;
                } else if input[ii] > DUCHAR_MAX {
                    *status = OVERFLOW_ERR;
                    output[ii] = u8::MAX;
                } else {
                    output[ii] = input[ii] as u8;
                }
                sptr += 4;
            }
        } else {
            /* must scale the data */

            for ii in 0..(ntodo as usize) {
                iret = dnan(shortBuffer[sptr]);
                if 0 != iret {
                    /* test for NaN or underflow */

                    if iret == 1 {
                        /* is it a NaN? */

                        *anynul = 1;
                        if nullcheck == NullCheckType::SetPixel {
                            output[ii] = nullval;
                        } else {
                            nullarray[ii] = 1;
                        }
                    } else {
                        /* it's an underflow */

                        if zero < DUCHAR_MIN {
                            *status = OVERFLOW_ERR;
                            output[ii] = 0;
                        } else if zero > DUCHAR_MAX {
                            *status = OVERFLOW_ERR;
                            output[ii] = u8::MAX;
                        } else {
                            output[ii] = zero as u8;
                        }
                    }
                } else {
                    dvalue = input[ii] * scale + zero;

                    if dvalue < DUCHAR_MIN {
                        *status = OVERFLOW_ERR;
                        output[ii] = 0;
                    } else if dvalue > DUCHAR_MAX {
                        *status = OVERFLOW_ERR;
                        output[ii] = u8::MAX;
                    } else {
                        output[ii] = dvalue as u8;
                    }
                }
                sptr += 4;
            }
        }
    }
    *status
}

/*--------------------------------------------------------------------------*/
/// Copy input to output following reading of the input from a FITS file. Check
/// for null values and do scaling if required. The nullcheck code value
/// determines how any null values in the input array are treated. A null
/// value is an input pixel that is equal to snull.  If nullcheck= 0, then
/// no special checking for nulls is performed.  If nullcheck = 1, then the
/// output pixel will be set = nullval if the corresponding input pixel is null.
/// If nullcheck = 2, then if the pixel is null then the corresponding value of
/// nullarray will be set to 1; the value of nullarray for non-null pixels
/// will = 0.  The anynul parameter will be set = 1 if any of the returned
/// pixels are null, otherwise anynul will be returned with a value = 0;
pub(crate) fn fffstri1(
    input: &mut [c_char],     /* I - array of values to be converted     */
    ntodo: c_long,            /* I - number of elements in the array     */
    scale: f64,               /* I - FITS TSCALn or BSCALE value         */
    zero: f64,                /* I - FITS TZEROn or BZERO  value         */
    twidth: c_long,           /* I - width of each substring of chars    */
    implipower: f64,          /* I - power of 10 of implied decimal      */
    nullcheck: NullCheckType, /* I - null checking code; 0 = don't check */
    /*     1:set null pixels = nullval         */
    /*     2: if null pixel, set nullarray = 1 */
    snull: &[c_char],               /* I - value of FITS null string, if any   */
    nullval: u8,                    /* I - set null pixels, if nullcheck = 1  */
    nullarray: &mut [c_char],       /* O - bad pixel array, if nullcheck = 2   */
    mut anynul: Option<&mut c_int>, /* O - set to 1 if any pixels are null     */
    output: &mut [u8],              /* O - array of converted pixels          */
    status: &mut c_int,             /* IO - error status                       */
) -> c_int {
    let mut dvalue: f64 = 0.0;
    let mut message: [c_char; FLEN_ERRMSG] = [0; FLEN_ERRMSG];
    let tempstore: c_char = 0;
    let chrzero: c_char = bb(b'0'); // 49
    let mut val: f64 = 0.0;
    let mut power: f64 = 0.0;
    let mut exponent: c_int = 0;
    let mut sign: c_int = 0;
    let mut esign: c_int = 0;
    let mut decpt: c_int = 0;

    let nullen = strlen_safe(snull);
    let mut cptr: usize = 0; /* pointer to start of input string */

    let twidth = twidth as usize;

    for ii in 0..(ntodo as usize) {
        let originalCptr = cptr;

        /* temporarily insert a null terminator at end of the string */
        let tpos = cptr + twidth;
        let tempstore = input[tpos];
        input[tpos] = 0;

        /* check if null value is defined, and if the    */
        /* column string is identical to the null string */
        if snull[0] != ASCII_NULL_UNDEFINED && strncmp_safe(snull, &input[cptr..], nullen) == 0 {
            if nullcheck != NullCheckType::None {
                let anynul = anynul.as_deref_mut().unwrap(); //Original code doesn't check if null pointer
                *anynul = 1;
                if nullcheck == NullCheckType::SetPixel {
                    output[ii] = nullval;
                } else {
                    nullarray[ii] = 1;
                }
            }
            cptr += twidth;
        } else {
            /* value is not the null value, so decode it */
            /* remove any embedded blank characters from the string */

            decpt = 0;
            sign = 1;
            val = 0.0;
            power = 1.;
            exponent = 0;
            esign = 1;

            while input[cptr] == bb(b' ') {
                /* skip leading blanks */
                cptr += 1;
            }

            if input[cptr] == bb(b'-') || input[cptr] == bb(b'+') {
                /* check for leading sign */

                if input[cptr] == bb(b'-') {
                    sign = -1;
                }

                cptr += 1;

                while input[cptr] == bb(b' ') {
                    /* skip blanks between sign and value */
                    cptr += 1;
                }
            }

            while input[cptr] >= bb(b'0') && input[cptr] <= bb(b'9') {
                val = val * 10.0 + (input[cptr] - chrzero) as f64; /* accumulate the value */
                cptr += 1;

                while input[cptr] == bb(b' ') {
                    /* skip embedded blanks in the value */
                    cptr += 1;
                }
            }

            if input[cptr] == bb(b'.') || input[cptr] == bb(b',') {
                /* check for decimal point */

                decpt = 1;
                cptr += 1;
                while input[cptr] == bb(b' ') {
                    /* skip any blanks */
                    cptr += 1;
                }

                while input[cptr] >= bb(b'0') && input[cptr] <= bb(b'9') {
                    val = val * 10.0 + (input[cptr] - chrzero) as f64; /* accumulate the value */
                    power *= 10.;
                    cptr += 1;

                    while input[cptr] == bb(b' ') {
                        /* skip embedded blanks in the value */
                        cptr += 1;
                    }
                }
            }

            if input[cptr] == bb(b'E') || input[cptr] == bb(b'D') {
                /* check for exponent */

                cptr += 1;
                while input[cptr] == bb(b' ') {
                    /* skip blanks */
                    cptr += 1;
                }

                if input[cptr] == bb(b'-') || input[cptr] == bb(b'+') {
                    /* check for exponent sign */

                    if input[cptr] == bb(b'-') {
                        esign = -1;
                    }

                    cptr += 1;

                    while input[cptr] == bb(b' ') {
                        /* skip blanks between sign and exp */
                        cptr += 1;
                    }
                }

                while input[cptr] >= bb(b'0') && input[cptr] <= bb(b'9') {
                    exponent = exponent * 10 + (input[cptr] - chrzero) as c_int; /* accumulate exp */
                    cptr += 1;

                    while input[cptr] == bb(b' ') {
                        /* skip embedded blanks */
                        cptr += 1;
                    }
                }
            }

            if input[cptr] != 0 {
                /* should end up at the null terminator */

                int_snprintf!(
                    &mut message,
                    FLEN_ERRMSG,
                    "Cannot read number from ASCII table",
                );
                ffpmsg_slice(&message);
                int_snprintf!(
                    &mut message,
                    FLEN_ERRMSG,
                    "Column field = {}.",
                    slice_to_str!(&input[originalCptr..]),
                );
                ffpmsg_slice(&message);
                /* restore the char that was overwritten by the null */
                input[tpos] = tempstore;
                *status = BAD_C2D;
                return *status;
            }

            if decpt == 0 {
                /* if no explicit decimal, use implied */
                power = implipower;
            }
            dvalue = (sign as f64 * val / power) * 10.0_f64.powi(esign * exponent);

            dvalue = dvalue * scale + zero; /* apply the scaling */

            if dvalue < DUCHAR_MIN {
                *status = OVERFLOW_ERR;
                output[ii] = 0;
            } else if dvalue > DUCHAR_MAX {
                *status = OVERFLOW_ERR;
                output[ii] = u8::MAX;
            } else {
                output[ii] = dvalue as u8;
            }
        }
        /* restore the char that was overwritten by the null */
        input[tpos] = tempstore;
    }
    *status
}
