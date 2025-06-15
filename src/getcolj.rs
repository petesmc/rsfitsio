/*  This file, getcolj.c, contains routines that read data elements from   */
/*  a FITS image or table, with long data type.                            */

/*  The FITSIO software was written by William Pence at the High Energy    */
/*  Astrophysic Science Archive Research Center (HEASARC) at the NASA      */
/*  Goddard Space Flight Center.                                           */

use core::slice;
use std::ffi::CStr;
use std::{cmp, mem};

use crate::c_types::{c_char, c_int, c_long, c_short};

use bytemuck::{cast_slice, cast_slice_mut};

use crate::bb;
use crate::fitscore::{
    ffasfm_safe, ffgcprll, ffghdt_safe, ffpmsg_slice, ffpmsg_str, fits_is_compressed_image_safe,
};
use crate::imcompress::fits_read_compressed_img;
use crate::wrappers::*;
use crate::{NullCheckType, fitsio::*};
use crate::{NullValue, fitsio2::*};
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
pub unsafe extern "C" fn ffgpvj(
    fptr: *mut fitsfile, /* I - FITS file pointer                       */
    group: c_long,       /* I - group to read (1 = 1st group)           */
    firstelem: LONGLONG, /* I - first vector element to read (1 = 1st)  */
    nelem: LONGLONG,     /* I - number of values to read                */
    nulval: c_long,      /* I - value for undefined pixels              */
    array: *mut c_long,  /* O - array of values that are returned       */
    anynul: *mut c_int,  /* O - set to 1 if any values are null; else 0 */
    status: *mut c_int,  /* IO - error status                           */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let anynul = anynul.as_mut();

        let array = slice::from_raw_parts_mut(array, nelem as usize);

        ffgpvj_safe(fptr, group, firstelem, nelem, nulval, array, anynul, status)
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
pub fn ffgpvj_safe(
    fptr: &mut fitsfile,        /* I - FITS file pointer                       */
    group: c_long,              /* I - group to read (1 = 1st group)           */
    firstelem: LONGLONG,        /* I - first vector element to read (1 = 1st)  */
    nelem: LONGLONG,            /* I - number of values to read                */
    nulval: c_long,             /* I - value for undefined pixels              */
    array: &mut [c_long],       /* O - array of values that are returned       */
    anynul: Option<&mut c_int>, /* O - set to 1 if any values are null; else 0 */
    status: &mut c_int,         /* IO - error status                           */
) -> c_int {
    let mut row: c_long = 0;
    let nullcheck = NullCheckType::SetPixel;
    let mut nullvalue: c_long = 0;

    let mut dummy_nularray = vec![0; nelem as usize];

    if fits_is_compressed_image_safe(fptr, status) > 0 {
        /* this is a compressed image in a binary table */
        nullvalue = nulval; /* set local variable */

        todo!();
        /*
        fits_read_compressed_pixels(
            fptr,
            TLONG,
            firstelem,
            nelem,
            nullcheck,
            &nullvalue,
            array,
            ptr::null_mut(),
            anynul,
            status,
        );
        */
        return *status;
    }

    row = cmp::max(1, group);

    ffgclj(
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
pub unsafe extern "C" fn ffgpfj(
    fptr: *mut fitsfile,   /* I - FITS file pointer                       */
    group: c_long,         /* I - group to read (1 = 1st group)           */
    firstelem: LONGLONG,   /* I - first vector element to read (1 = 1st)  */
    nelem: LONGLONG,       /* I - number of values to read                */
    array: *mut c_long,    /* O - array of values that are returned       */
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

        ffgpfj_safe(
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
pub fn ffgpfj_safe(
    fptr: &mut fitsfile,        /* I - FITS file pointer                       */
    group: c_long,              /* I - group to read (1 = 1st group)           */
    firstelem: LONGLONG,        /* I - first vector element to read (1 = 1st)  */
    nelem: LONGLONG,            /* I - number of values to read                */
    array: &mut [c_long],       /* O - array of values that are returned       */
    nularray: &mut [c_char],    /* O - array of null pixel flags               */
    anynul: Option<&mut c_int>, /* O - set to 1 if any values are null; else 0 */
    status: &mut c_int,         /* IO - error status                           */
) -> c_int {
    let mut row = 0;
    let nullcheck = NullCheckType::SetNullArray;

    if fits_is_compressed_image_safe(fptr, status) > 0 {
        /* this is a compressed image in a binary table */

        todo!();
        /*
        fits_read_compressed_pixels(
            fptr,
            TLONG,
            firstelem,
            nelem,
            nullcheck,
            ptr::null_mut(),
            array,
            nularray,
            anynul,
            status,
        );
        */
        return *status;
    }

    row = cmp::max(1, group);

    ffgclj(
        fptr,
        2,
        row as LONGLONG,
        firstelem,
        nelem,
        1,
        NullCheckType::SetNullArray,
        0,
        array,
        nularray,
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
pub unsafe extern "C" fn ffg2dj(
    fptr: *mut fitsfile, /* I - FITS file pointer                       */
    group: c_long,       /* I - group to read (1 = 1st group)           */
    nulval: c_long,      /* set undefined pixels equal to this          */
    ncols: LONGLONG,     /* I - number of pixels in each row of array   */
    naxis1: LONGLONG,    /* I - FITS image NAXIS1 value                 */
    naxis2: LONGLONG,    /* I - FITS image NAXIS2 value                 */
    array: *mut c_long,  /* O - array to be filled and returned         */
    anynul: *mut c_int,  /* O - set to 1 if any values are null; else 0 */
    status: *mut c_int,  /* IO - error status                           */
) -> c_int {
    unsafe {
        /* call the 3D reading routine, with the 3rd dimension = 1 */

        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let anynul = anynul.as_mut();

        let array = slice::from_raw_parts_mut(array, (ncols * naxis2) as usize);

        ffg2dj_safe(
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
pub fn ffg2dj_safe(
    fptr: &mut fitsfile,        /* I - FITS file pointer                       */
    group: c_long,              /* I - group to read (1 = 1st group)           */
    nulval: c_long,             /* set undefined pixels equal to this     */
    ncols: LONGLONG,            /* I - number of pixels in each row of array   */
    naxis1: LONGLONG,           /* I - FITS image NAXIS1 value                 */
    naxis2: LONGLONG,           /* I - FITS image NAXIS2 value                 */
    array: &mut [c_long],       /* O - array to be filled and returned    */
    anynul: Option<&mut c_int>, /* O - set to 1 if any values are null; else 0 */
    status: &mut c_int,         /* IO - error status                           */
) -> c_int {
    /* call the 3D reading routine, with the 3rd dimension = 1 */
    ffg3dj_safe(
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
pub unsafe extern "C" fn ffg3dj(
    fptr: *mut fitsfile, /* I - FITS file pointer                       */
    group: c_long,       /* I - group to read (1 = 1st group)           */
    nulval: c_long,      /* set undefined pixels equal to this          */
    ncols: LONGLONG,     /* I - number of pixels in each row of array   */
    nrows: LONGLONG,     /* I - number of rows in each plane of array   */
    naxis1: LONGLONG,    /* I - FITS image NAXIS1 value                 */
    naxis2: LONGLONG,    /* I - FITS image NAXIS2 value                 */
    naxis3: LONGLONG,    /* I - FITS image NAXIS3 value                 */
    array: *mut c_long,  /* O - array to be filled and returned         */
    anynul: *mut c_int,  /* O - set to 1 if any values are null; else 0 */
    status: *mut c_int,  /* IO - error status                           */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let anynul = anynul.as_mut();

        let array = slice::from_raw_parts_mut(array, (ncols * nrows) as usize);

        ffg3dj_safe(
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
pub fn ffg3dj_safe(
    fptr: &mut fitsfile,            /* I - FITS file pointer                       */
    group: c_long,                  /* I - group to read (1 = 1st group)           */
    nulval: c_long,                 /* set undefined pixels equal to this          */
    ncols: LONGLONG,                /* I - number of pixels in each row of array   */
    nrows: LONGLONG,                /* I - number of rows in each plane of array   */
    naxis1: LONGLONG,               /* I - FITS image NAXIS1 value                 */
    naxis2: LONGLONG,               /* I - FITS image NAXIS2 value                 */
    naxis3: LONGLONG,               /* I - FITS image NAXIS3 value                 */
    array: &mut [c_long],           /* O - array to be filled and returned         */
    mut anynul: Option<&mut c_int>, /* O - set to 1 if any values are null; else 0 */
    status: &mut c_int,             /* IO - error status                           */
) -> c_int {
    let mut tablerow = 0;

    let nullcheck = NullCheckType::SetPixel;
    let inc: [c_long; 3] = [1, 1, 1];
    let fpixel: [LONGLONG; 3] = [1, 1, 1];
    let mut nfits = 0;
    let mut narray = 0;
    let mut lpixel: [LONGLONG; 3] = [0; 3];
    let mut nullvalue: LONGLONG = 0;

    let mut dummy_nularray = vec![0; (ncols * nrows) as usize];

    if fits_is_compressed_image_safe(fptr, status) > 0 {
        /* this is a compressed image in a binary table */

        lpixel[0] = ncols as LONGLONG;
        lpixel[1] = nrows as LONGLONG;
        lpixel[2] = naxis3 as LONGLONG;
        nullvalue = nulval as LONGLONG; /* set local variable */

        fits_read_compressed_img(
            fptr,
            TLONG,
            &fpixel,
            &lpixel,
            &inc,
            nullcheck,
            &Some(NullValue::Long(nullvalue as c_long)),
            cast_slice_mut(array),
            None,
            anynul.as_deref_mut(),
            status,
        );

        return *status;
    }

    tablerow = cmp::max(1, group);

    if ncols == naxis1 && nrows == naxis2 {
        /* arrays have same size? */

        /* all the image pixels are contiguous, so read all at once */
        ffgclj(
            fptr,
            2,
            tablerow as LONGLONG,
            1,
            naxis1 * naxis2 * naxis3,
            1,
            NullCheckType::SetPixel,
            nulval,
            cast_slice_mut(array),
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
            if ffgclj(
                fptr,
                2,
                tablerow as LONGLONG,
                nfits,
                naxis1,
                1,
                NullCheckType::SetPixel,
                nulval,
                cast_slice_mut(&mut array[(narray as usize)..]),
                &mut dummy_nularray,
                anynul.as_deref_mut(),
                status,
            ) > 0
            {
                return *status;
            }

            nfits += naxis1;
            narray += ncols;
        }
        narray += (nrows - naxis2) * ncols;
    }

    *status
}

/*--------------------------------------------------------------------------*/
/// Read a subsection of data values from an image or a table column.
/// This routine is set up to handle a maximum of nine dimensions.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgsvj(
    fptr: *mut fitsfile,  /* I - FITS file pointer                         */
    colnum: c_int,        /* I - number of the column to read (1 = 1st)    */
    naxis: c_int,         /* I - number of dimensions in the FITS array    */
    naxes: *const c_long, /* I - size of each dimension                    */
    blc: *const c_long,   /* I - 'bottom left corner' of the subsection    */
    trc: *const c_long,   /* I - 'top right corner' of the subsection      */
    inc: *const c_long,   /* I - increment to be applied in each dimension */
    nulval: c_long,       /* I - value to set undefined pixels             */
    array: *mut c_long,   /* O - array to be filled and returned           */
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

        ffgsvj_safe(
            fptr, colnum, naxis, naxes, blc, trc, inc, nulval, array, anynul, status,
        )
    }
}

/*--------------------------------------------------------------------------*/
/// Read a subsection of data values from an image or a table column.
/// This routine is set up to handle a maximum of nine dimensions.
pub fn ffgsvj_safe(
    fptr: &mut fitsfile,  /* I - FITS file pointer                         */
    colnum: c_int,        /* I - number of the column to read (1 = 1st)    */
    naxis: c_int,         /* I - number of dimensions in the FITS array    */
    naxes: &[c_long],     /* I - size of each dimension                    */
    blc: &[c_long],       /* I - 'bottom left corner' of the subsection    */
    trc: &[c_long],       /* I - 'top right corner' of the subsection      */
    inc: &[c_long],       /* I - increment to be applied in each dimension */
    nulval: c_long,       /* I - value to set undefined pixels             */
    array: &mut [c_long], /* O - array to be filled and returned           */
    mut anynul: Option<&mut c_int>, /* O - set to 1 if any values are null; else 0   */
    status: &mut c_int,   /* IO - error status                             */
) -> c_int {
    let row: c_long = 0;
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
    let ldummy: c_char = 0;
    let mut msg: [c_char; FLEN_ERRMSG] = [0; FLEN_ERRMSG];
    let nullcheck = NullCheckType::SetPixel;
    let mut nullvalue: c_long = 0;

    let naxis: usize = naxis as usize;

    if naxis < 1 || naxis > 9 {
        int_snprintf!(
            &mut msg,
            FLEN_ERRMSG,
            "NAXIS = {} in call to ffgsvj is out of range",
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

        fits_read_compressed_img(
            fptr,
            TLONG,
            &blcll,
            &trcll,
            inc,
            nullcheck,
            &Some(NullValue::Long(nullvalue)),
            cast_slice_mut(array),
            None,
            anynul,
            status,
        );

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
                    "ffgsvj: illegal range specified for axis {}",
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

                                        if ffgclj(
                                            fptr,
                                            numcol as c_int,
                                            row as LONGLONG,
                                            felem as LONGLONG,
                                            nelem as LONGLONG,
                                            ninc,
                                            nultyp,
                                            nulval,
                                            &mut array[(i0 as usize)..],
                                            &mut [0],
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
                                        i0 += nelem;
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

/*--------------------------------------------------------------------------*/
/// Read a subsection of data values from an image or a table column.
/// This routine is set up to handle a maximum of nine dimensions.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgsfj(
    fptr: *mut fitsfile,  /* I - FITS file pointer                         */
    colnum: c_int,        /* I - number of the column to read (1 = 1st)    */
    naxis: c_int,         /* I - number of dimensions in the FITS array    */
    naxes: *const c_long, /* I - size of each dimension                    */
    blc: *const c_long,   /* I - 'bottom left corner' of the subsection    */
    trc: *const c_long,   /* I - 'top right corner' of the subsection      */
    inc: *const c_long,   /* I - increment to be applied in each dimension */
    array: *mut c_long,   /* O - array to be filled and returned           */
    flagval: *mut c_char, /* O - set to 1 if corresponding value is null   */
    anynul: *mut c_int,   /* O - set to 1 if any values are null; else 0   */
    status: *mut c_int,   /* IO - error status                             */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        let naxes = slice::from_raw_parts(naxes, naxis as usize);
        let blc = slice::from_raw_parts(blc, naxis as usize);
        let trc = slice::from_raw_parts(trc, naxis as usize);
        let inc = slice::from_raw_parts(inc, naxis as usize);

        let total_nelem = calculate_subsection_length(blc, trc, inc);
        let array = slice::from_raw_parts_mut(array, total_nelem);
        let flagval = slice::from_raw_parts_mut(flagval, total_nelem);

        let anynul = anynul.as_mut();

        ffgsfj_safe(
            fptr, colnum, naxis, naxes, blc, trc, inc, array, flagval, anynul, status,
        )
    }
}

pub fn ffgsfj_safe(
    fptr: &mut fitsfile,    /* I - FITS file pointer                         */
    colnum: c_int,          /* I - number of the column to read (1 = 1st)    */
    naxis: c_int,           /* I - number of dimensions in the FITS array    */
    naxes: &[c_long],       /* I - size of each dimension                    */
    blc: &[c_long],         /* I - 'bottom left corner' of the subsection    */
    trc: &[c_long],         /* I - 'top right corner' of the subsection      */
    inc: &[c_long],         /* I - increment to be applied in each dimension */
    array: &mut [c_long],   /* O - array to be filled and returned           */
    flagval: &mut [c_char], /* O - set to 1 if corresponding value is null   */
    mut anynul: Option<&mut c_int>, /* O - set to 1 if any values are null; else 0   */
    status: &mut c_int,     /* IO - error status                             */
) -> c_int {
    let row: c_long = 0;
    let mut rstr: c_long = 0;
    let mut rstp: c_long = 0;
    let mut rinc: c_long = 0;
    let mut str: [c_long; 9] = [0; 9];
    let mut stp: [c_long; 9] = [0; 9];
    let mut incr: [c_long; 9] = [0; 9];
    let mut dsize: [c_long; 10] = [0; 10];
    let mut blcll: [LONGLONG; 9] = [0; 9];
    let mut trcll: [LONGLONG; 9] = [0; 9];
    let mut felem: c_long = 0;
    let mut nelem: c_long = 0;
    let mut nultyp: NullCheckType = NullCheckType::None;
    let mut ninc: c_long = 0;
    let mut numcol: c_long = 0;
    let nulval: c_long = 0;
    let mut hdutype: c_int = 0;
    let mut anyf: c_int = 0;
    let mut msg: [c_char; FLEN_ERRMSG] = [0; FLEN_ERRMSG];
    let nullcheck = NullCheckType::SetNullArray;
    let naxis: usize = naxis as usize;

    if naxis < 1 || naxis > 9 {
        int_snprintf!(
            &mut msg,
            FLEN_ERRMSG,
            "NAXIS = {} in call to ffgsvj is out of range",
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

        fits_read_compressed_img(
            fptr,
            TLONG,
            &blcll,
            &trcll,
            inc,
            nullcheck,
            &None,
            cast_slice_mut(array),
            Some(flagval),
            anynul,
            status,
        );

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
            int_snprintf!(
                &mut msg,
                FLEN_ERRMSG,
                "ffgsvj: illegal range specified for axis {}",
                ii + 1,
            );
            ffpmsg_slice(&msg);
            *status = BAD_PIX_NUM;
            return *status;
        }

        str[ii] = blc[ii];
        stp[ii] = trc[ii];
        incr[ii] = inc[ii];
        dsize[ii + 1] = dsize[ii] * naxes[ii] as c_long;
    }

    if naxis == 1 && naxes[0] == 1 {
        /* This is not a vector column, so read all the rows at once */
        nelem = (rstp - rstr) / rinc + 1;
        ninc = rinc;
        rstp = rstr;
    } else {
        /* have to read each row individually, in all dimensions */
        nelem = (stp[0] - str[0]) / inc[0] + 1;
        ninc = incr[0];
    }

    for row in (rstr..=rstp).step_by(rinc as usize) {
        for i8 in ((str[8])..=(stp[8])).step_by(incr[8] as usize) {
            for i7 in ((str[7])..=(stp[7])).step_by(incr[7] as usize) {
                for i6 in ((str[6])..=(stp[6])).step_by(incr[6] as usize) {
                    for i5 in ((str[5])..=(stp[5])).step_by(incr[5] as usize) {
                        for i4 in ((str[4])..=(stp[4])).step_by(incr[4] as usize) {
                            for i3 in ((str[3])..=(stp[3])).step_by(incr[3] as usize) {
                                for i2 in ((str[2])..=(stp[2])).step_by(incr[2] as usize) {
                                    for i1 in ((str[1])..=(stp[1])).step_by(incr[1] as usize) {
                                        felem = (str[0] as c_long)
                                            + (i1 as c_long - 1) * dsize[1]
                                            + (i2 as c_long - 1) * dsize[2]
                                            + (i3 as c_long - 1) * dsize[3]
                                            + (i4 as c_long - 1) * dsize[4]
                                            + (i5 as c_long - 1) * dsize[5]
                                            + (i6 as c_long - 1) * dsize[6]
                                            + (i7 as c_long - 1) * dsize[7]
                                            + (i8 as c_long - 1) * dsize[8];

                                        if ffgclj(
                                            fptr,
                                            numcol as c_int,
                                            row as LONGLONG,
                                            felem as LONGLONG,
                                            nelem as LONGLONG,
                                            ninc,
                                            nultyp,
                                            nulval,
                                            &mut array[(i0 as usize)..],
                                            &mut flagval[(i0 as usize)..],
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
                                        i0 += nelem;
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

/*--------------------------------------------------------------------------*/
/// Read an array of group parameters from the primary array. Data conversion
/// and scaling will be performed if necessary (e.g, if the datatype of
/// the FITS array is not the same as the array being read).
///
/// The primary array is represented as a binary table:
/// each group of the primary array is a row in the table,
/// where the first column contains the group parameters
/// and the second column contains the image itself.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffggpj(
    fptr: *mut fitsfile, /* I - FITS file pointer                       */
    group: c_long,       /* I - group to read (1 = 1st group)           */
    firstelem: c_long,   /* I - first vector element to read (1 = 1st)  */
    nelem: c_long,       /* I - number of values to read                */
    array: *mut c_long,  /* O - array of values that are returned       */
    status: *mut c_int,  /* IO - error status                           */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let array = slice::from_raw_parts_mut(array, nelem as usize);

        ffggpj_safe(fptr, group, firstelem, nelem, array, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Read an array of group parameters from the primary array. Data conversion
/// and scaling will be performed if necessary (e.g, if the datatype of
/// the FITS array is not the same as the array being read).
///
/// The primary array is represented as a binary table:
/// each group of the primary array is a row in the table,
/// where the first column contains the group parameters
/// and the second column contains the image itself.
pub fn ffggpj_safe(
    fptr: &mut fitsfile,  /* I - FITS file pointer                       */
    group: c_long,        /* I - group to read (1 = 1st group)           */
    firstelem: c_long,    /* I - first vector element to read (1 = 1st)  */
    nelem: c_long,        /* I - number of values to read                */
    array: &mut [c_long], /* O - array of values that are returned       */
    status: &mut c_int,   /* IO - error status                           */
) -> c_int {
    let mut dummy_nularray = vec![0; nelem as usize];
    let mut anynul = 0;

    let row = cmp::max(1, group);

    ffgclj(
        fptr,
        2,
        row as LONGLONG,
        firstelem as LONGLONG,
        nelem as LONGLONG,
        1,
        NullCheckType::SetPixel,
        0,
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
pub unsafe extern "C" fn ffgcvj(
    fptr: *mut fitsfile, /* I - FITS file pointer                       */
    colnum: c_int,       /* I - number of column to read (1 = 1st col)  */
    firstrow: LONGLONG,  /* I - first row to read (1 = 1st row)         */
    firstelem: LONGLONG, /* I - first vector element to read (1 = 1st)  */
    nelem: LONGLONG,     /* I - number of values to read                */
    nulval: c_long,      /* I - value for null pixels                   */
    array: *mut c_long,  /* O - array of values that are read           */
    anynul: *mut c_int,  /* O - set to 1 if any values are null; else 0 */
    status: *mut c_int,  /* IO - error status                           */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let anynul = anynul.as_mut();

        let array = slice::from_raw_parts_mut(array, nelem as usize);

        ffgcvj_safe(
            fptr, colnum, firstrow, firstelem, nelem, nulval, array, anynul, status,
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
pub fn ffgcvj_safe(
    fptr: &mut fitsfile,        /* I - FITS file pointer                       */
    colnum: c_int,              /* I - number of column to read (1 = 1st col)  */
    firstrow: LONGLONG,         /* I - first row to read (1 = 1st row)         */
    firstelem: LONGLONG,        /* I - first vector element to read (1 = 1st)  */
    nelem: LONGLONG,            /* I - number of values to read                */
    nulval: c_long,             /* I - value for null pixels                   */
    array: &mut [c_long],       /* O - array of values that are read           */
    anynul: Option<&mut c_int>, /* O - set to 1 if any values are null; else 0 */
    status: &mut c_int,         /* IO - error status                           */
) -> c_int {
    let mut dummy_nularray = vec![0; nelem as usize];

    ffgclj(
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
pub unsafe extern "C" fn ffgcfj(
    fptr: *mut fitsfile,   /* I - FITS file pointer                       */
    colnum: c_int,         /* I - number of column to read (1 = 1st col)  */
    firstrow: LONGLONG,    /* I - first row to read (1 = 1st row)         */
    firstelem: LONGLONG,   /* I - first vector element to read (1 = 1st)  */
    nelem: LONGLONG,       /* I - number of values to read                */
    array: *mut c_long,    /* O - array of values that are read           */
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

        ffgcfj_safe(
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
pub fn ffgcfj_safe(
    fptr: &mut fitsfile,        /* I - FITS file pointer                       */
    colnum: c_int,              /* I - number of column to read (1 = 1st col)  */
    firstrow: LONGLONG,         /* I - first row to read (1 = 1st row)         */
    firstelem: LONGLONG,        /* I - first vector element to read (1 = 1st)  */
    nelem: LONGLONG,            /* I - number of values to read                */
    array: &mut [c_long],       /* O - array of values that are read           */
    nularray: &mut [c_char],    /* O - array of flags: 1 if null pixel; else 0 */
    anynul: Option<&mut c_int>, /* O - set to 1 if any values are null; else 0 */
    status: &mut c_int,         /* IO - error status                           */
) -> c_int {
    let dummy = 0;

    ffgclj(
        fptr,
        colnum,
        firstrow,
        firstelem,
        nelem,
        1,
        NullCheckType::SetNullArray,
        dummy,
        array,
        nularray,
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
pub(crate) fn ffgclj(
    fptr: &mut fitsfile,   /* I - FITS file pointer                       */
    colnum: c_int,         /* I - number of column to read (1 = 1st col)  */
    firstrow: LONGLONG,    /* I - first row to read (1 = 1st row)         */
    firstelem: LONGLONG,   /* I - first vector element to read (1 = 1st)  */
    nelem: LONGLONG,       /* I - number of values to read                */
    elemincre: c_long,     /* I - pixel increment; e.g., 2 = every other  */
    nultyp: NullCheckType, /* I - null value handling code:               */
    /*     1: set undefined pixels = nulval        */
    /*     2: set nularray=1 for undefined pixels  */
    nulval: c_long,                 /* I - value for null pixels if nultyp = 1     */
    array: &mut [c_long],           /* O - array of values that are read           */
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
    let ii: c_long = 0;
    let mut xwidth: c_long = 0;
    let mut ntodo: c_long = 0;
    let mut convert = false;
    let mut nulcheck = NullCheckType::None;
    let mut readcheck: c_int = 0;
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
    let mut cbuff: [f64; DBUFFSIZE as usize / mem::size_of::<f64>()] =
        [0.0; DBUFFSIZE as usize / mem::size_of::<f64>()]; /* align cbuff on word boundary */

    if *status > 0 || nelem == 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    if let Some(anynul) = anynul.as_deref_mut() {
        *anynul = 0;
    }

    if nultyp == NullCheckType::SetNullArray {
        nularray.fill(0); /* initialize nullarray */
    }
    /*---------------------------------------------------*/
    /*  Check input and get parameters about the column: */
    /*---------------------------------------------------*/
    if elemincre < 0 {
        readcheck = -1; /* don't do range checking in this case */
    }
    if ffgcprll(
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
    ) > 0
    {
        return *status;
    }
    maxelem = maxelem2 as LONGLONG;

    incre *= elemincre; /* multiply incre to just get every nth pixel */

    if tcode == TSTRING {
        /* setup for ASCII tables */
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
    nulcheck = nultyp; /* by default check for null values in the FITS file */

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

    if (tcode == TLONG) && (LONGSIZE == 32) {
        /* Special Case:                        */
        /* no type convertion required, so read */
        /* data directly into output buffer.    */

        if nelem < INT32_MAX as LONGLONG / 4 {
            maxelem = nelem;
        } else {
            maxelem = (INT32_MAX / 4) as LONGLONG;
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

    let buffer = &mut cbuff;

    while remain > 0 {
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

        readptr =
            startpos + (rownum as LONGLONG * rowlen) + (elemnum * (incre / elemincre) as LONGLONG);

        match tcode {
            TLONG => {
                if LONGSIZE == 32 {
                    if convert {
                        ffgi4b(fptr, readptr, ntodo, incre, cast_slice_mut(buffer), status);

                        fffi4i4(
                            cast_slice_mut(buffer),
                            ntodo,
                            scale,
                            zero,
                            nulcheck,
                            tnull as c_int,
                            nulval,
                            cast_slice_mut(&mut nularray[next..]),
                            anynul.as_deref_mut(),
                            &mut array[next..],
                            status,
                        );
                    } else {
                        ffgi4b(
                            fptr,
                            readptr,
                            ntodo,
                            incre,
                            cast_slice_mut(&mut array[next..]),
                            status,
                        );
                    }
                } else {
                    /* case where sizeof(long) = 8 */

                    if convert {
                        ffgi4b(fptr, readptr, ntodo, incre, cast_slice_mut(buffer), status);

                        fffi4i4(
                            cast_slice_mut(buffer),
                            ntodo,
                            scale,
                            zero,
                            nulcheck,
                            tnull as c_int,
                            nulval,
                            cast_slice_mut(&mut nularray[next..]),
                            anynul.as_deref_mut(),
                            cast_slice_mut(&mut array[next..]),
                            status,
                        );
                    } else {
                        ffgi4b(
                            fptr,
                            readptr,
                            ntodo,
                            incre,
                            cast_slice_mut(&mut array[next..]),
                            status,
                        );
                    }
                }
            }
            TLONGLONG => {
                ffgi8b(fptr, readptr, ntodo, incre, cast_slice_mut(buffer), status);

                fffi8i4(
                    cast_slice_mut(buffer),
                    ntodo,
                    scale,
                    zero,
                    nulcheck,
                    tnull,
                    nulval,
                    cast_slice_mut(&mut nularray[next..]),
                    anynul.as_deref_mut(),
                    cast_slice_mut(&mut array[next..]),
                    status,
                );
            }
            TBYTE => {
                ffgi1b(fptr, readptr, ntodo, incre, cast_slice_mut(buffer), status);

                fffi1i4(
                    cast_slice_mut(buffer),
                    ntodo,
                    scale,
                    zero,
                    nulcheck,
                    tnull as u8,
                    nulval,
                    cast_slice_mut(&mut nularray[next..]),
                    anynul.as_deref_mut(),
                    cast_slice_mut(&mut array[next..]),
                    status,
                );
            }
            TSHORT => {
                ffgi2b(fptr, readptr, ntodo, incre, cast_slice_mut(buffer), status);

                fffi2i4(
                    cast_slice_mut(buffer),
                    ntodo,
                    scale,
                    zero,
                    nulcheck,
                    tnull as c_short,
                    nulval,
                    cast_slice_mut(&mut nularray[next..]),
                    anynul.as_deref_mut(),
                    cast_slice_mut(&mut array[next..]),
                    status,
                );
            }
            TFLOAT => {
                ffgr4b(fptr, readptr, ntodo, incre, cast_slice_mut(buffer), status);

                fffr4i4(
                    cast_slice_mut(buffer),
                    ntodo,
                    scale,
                    zero,
                    nulcheck,
                    nulval,
                    cast_slice_mut(&mut nularray[next..]),
                    anynul.as_deref_mut(),
                    cast_slice_mut(&mut array[next..]),
                    status,
                );
            }
            TDOUBLE => {
                ffgr8b(fptr, readptr, ntodo, incre, buffer, status);

                fffr8i4(
                    cast_slice_mut(buffer),
                    ntodo,
                    scale,
                    zero,
                    nulcheck,
                    nulval,
                    cast_slice_mut(&mut nularray[next..]),
                    anynul.as_deref_mut(),
                    cast_slice_mut(&mut array[next..]),
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
                        cast_slice_mut(buffer),
                        status,
                    );
                } else {
                    ffgbytoff(
                        fptr,
                        twidth,
                        ntodo,
                        incre - twidth,
                        cast_slice_mut(buffer),
                        status,
                    );
                }

                fffstri4(
                    cast_slice_mut(buffer),
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
                    cast_slice_mut(&mut array[next..]),
                    status,
                );
            }
            _ => {
                /*  error trap for invalid column format */
                int_snprintf!(
                    &mut message,
                    FLEN_ERRMSG,
                    "Cannot read numbers from column {} which has format {}",
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
                    "Error reading elements {:.0} thru {:.0} from column {} (ffgclj).",
                    dtemp + 1.0,
                    dtemp + ntodo as f64,
                    colnum,
                );
            } else {
                int_snprintf!(
                    &mut message,
                    FLEN_ERRMSG,
                    "Error reading elements {:.0} thru {:.0} from image (ffgclj).",
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
        if remain > 0 {
            next += ntodo as usize;
            elemnum += (ntodo * elemincre) as LONGLONG;

            if elemnum >= repeat {
                /* completed a row; start on later row */
                rowincre = elemnum / repeat;
                rownum += rowincre;
                elemnum -= rowincre * repeat;
            } else if elemnum < 0 {
                /* completed a row; start on a previous row */
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
/// Copy input to output following reading of the input from a FITS file.
/// Check for null values and do datatype conversion and scaling if required.
/// The nullcheck code value determines how any null values in the input array
/// are treated.  A null value is an input pixel that is equal to tnull.  If
/// nullcheck = 0, then no checking for nulls is performed and any null values
/// will be transformed just like any other pixel.  If nullcheck = 1, then the
/// output pixel will be set = nullval if the corresponding input pixel is null.
/// If nullcheck = 2, then if the pixel is null then the corresponding value of
/// nullarray will be set to 1; the value of nullarray for non-null pixels
/// will = 0.  The anynull parameter will be set = 1 if any of the returned
/// pixels are null, otherwise anynull will be returned with a value = 0;
pub(crate) fn fffi1i4(
    input: &[u8],             /* I - array of values to be converted     */
    ntodo: c_long,            /* I - number of elements in the array     */
    scale: f64,               /* I - FITS TSCALn or BSCALE value         */
    zero: f64,                /* I - FITS TZEROn or BZERO  value         */
    nullcheck: NullCheckType, /* I - null checking code; 0 = don't check */
    /*     1:set null pixels = nullval         */
    /*     2: if null pixel, set nullarray = 1 */
    tnull: u8,                   /* I - value of FITS TNULLn keyword if any */
    nullval: c_long,             /* I - set null pixels, if nullcheck = 1   */
    nullarray: &mut [c_char],    /* I - bad pixel array, if nullcheck = 2   */
    anynull: Option<&mut c_int>, /* O - set to 1 if any pixels are null     */
    output: &mut [c_long],       /* O - array of converted pixels           */
    status: &mut c_int,          /* IO - error status                       */
) -> c_int {
    let mut dvalue: f64 = 0.0;

    if nullcheck == NullCheckType::None {
        /* no null checking required */
        if scale == 1.0 && zero == 0.0 {
            /* no scaling */
            for ii in 0..(ntodo as usize) {
                output[ii] = input[ii] as c_long; /* copy input to output */
            }
        } else {
            /* must scale the data */
            for ii in 0..(ntodo as usize) {
                dvalue = input[ii] as f64 * scale + zero;

                if dvalue < DLONG_MIN {
                    *status = OVERFLOW_ERR;
                    output[ii] = LONG_MIN;
                } else if dvalue > DLONG_MAX {
                    *status = OVERFLOW_ERR;
                    output[ii] = LONG_MAX;
                } else {
                    output[ii] = dvalue as c_long;
                }
            }
        }
    } else {
        /* must check for null values */
        let anynull = anynull.unwrap();

        if scale == 1.0 && zero == 0.0 {
            /* no scaling */
            for ii in 0..(ntodo as usize) {
                if input[ii] == tnull {
                    *anynull = 1;
                    if nullcheck == NullCheckType::SetPixel {
                        output[ii] = nullval;
                    } else {
                        nullarray[ii] = 1;
                    }
                } else {
                    output[ii] = input[ii] as c_long;
                }
            }
        } else {
            /* must scale the data */
            for ii in 0..(ntodo as usize) {
                if input[ii] == tnull {
                    *anynull = 1;
                    if nullcheck == NullCheckType::SetPixel {
                        output[ii] = nullval;
                    } else {
                        nullarray[ii] = 1;
                    }
                } else {
                    dvalue = input[ii] as f64 * scale + zero;

                    if dvalue < DLONG_MIN {
                        *status = OVERFLOW_ERR;
                        output[ii] = LONG_MIN;
                    } else if dvalue > DLONG_MAX {
                        *status = OVERFLOW_ERR;
                        output[ii] = LONG_MAX;
                    } else {
                        output[ii] = dvalue as c_long;
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
/// will = 0.  The anynull parameter will be set = 1 if any of the returned
/// pixels are null, otherwise anynull will be returned with a value = 0;
pub(crate) fn fffi2i4(
    input: &[c_short],        /* I - array of values to be converted     */
    ntodo: c_long,            /* I - number of elements in the array     */
    scale: f64,               /* I - FITS TSCALn or BSCALE value         */
    zero: f64,                /* I - FITS TZEROn or BZERO  value         */
    nullcheck: NullCheckType, /* I - null checking code; 0 = don't check */
    /*     1:set null pixels = nullval         */
    /*     2: if null pixel, set nullarray = 1 */
    tnull: c_short,              /* I - value of FITS TNULLn keyword if any */
    nullval: c_long,             /* I - set null pixels, if nullcheck = 1   */
    nullarray: &mut [c_char],    /* I - bad pixel array, if nullcheck = 2   */
    anynull: Option<&mut c_int>, /* O - set to 1 if any pixels are null     */
    output: &mut [c_long],       /* O - array of converted pixels           */
    status: &mut c_int,          /* IO - error status                       */
) -> c_int {
    let mut dvalue: f64 = 0.0;

    if nullcheck == NullCheckType::None {
        /* no null checking required */
        if scale == 1.0 && zero == 0.0 {
            /* no scaling */
            for ii in 0..(ntodo as usize) {
                output[ii] = input[ii] as c_long; /* copy input to output */
            }
        } else {
            /* must scale the data */
            for ii in 0..(ntodo as usize) {
                dvalue = input[ii] as f64 * scale + zero;

                if dvalue < DLONG_MIN {
                    *status = OVERFLOW_ERR;
                    output[ii] = LONG_MIN;
                } else if dvalue > DLONG_MAX {
                    *status = OVERFLOW_ERR;
                    output[ii] = LONG_MAX;
                } else {
                    output[ii] = dvalue as c_long;
                }
            }
        }
    } else {
        /* must check for null values */
        let anynull = anynull.unwrap();

        if scale == 1.0 && zero == 0.0 {
            /* no scaling */
            for ii in 0..(ntodo as usize) {
                if input[ii] == tnull {
                    *anynull = 1;
                    if nullcheck == NullCheckType::SetPixel {
                        output[ii] = nullval;
                    } else {
                        nullarray[ii] = 1;
                    }
                } else {
                    output[ii] = input[ii] as c_long;
                }
            }
        } else {
            /* must scale the data */
            for ii in 0..(ntodo as usize) {
                if input[ii] == tnull {
                    *anynull = 1;
                    if nullcheck == NullCheckType::SetPixel {
                        output[ii] = nullval;
                    } else {
                        nullarray[ii] = 1;
                    }
                } else {
                    dvalue = input[ii] as f64 * scale + zero;

                    if dvalue < DLONG_MIN {
                        *status = OVERFLOW_ERR;
                        output[ii] = LONG_MIN;
                    } else if dvalue > DLONG_MAX {
                        *status = OVERFLOW_ERR;
                        output[ii] = LONG_MAX;
                    } else {
                        output[ii] = dvalue as c_long;
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
/// will = 0.  The anynull parameter will be set = 1 if any of the returned
/// pixels are null, otherwise anynull will be returned with a value = 0;
pub(crate) fn fffi4i4(
    input: &[INT32BIT],       /* I - array of values to be converted     */
    ntodo: c_long,            /* I - number of elements in the array     */
    scale: f64,               /* I - FITS TSCALn or BSCALE value         */
    zero: f64,                /* I - FITS TZEROn or BZERO  value         */
    nullcheck: NullCheckType, /* I - null checking code; 0 = don't check */
    /*     1:set null pixels = nullval         */
    /*     2: if null pixel, set nullarray = 1 */
    tnull: INT32BIT,             /* I - value of FITS TNULLn keyword if any */
    nullval: c_long,             /* I - set null pixels, if nullcheck = 1   */
    nullarray: &mut [c_char],    /* I - bad pixel array, if nullcheck = 2   */
    anynull: Option<&mut c_int>, /* O - set to 1 if any pixels are null     */
    output: &mut [c_long],       /* O - array of converted pixels           */
    status: &mut c_int,          /* IO - error status                       */
) -> c_int {
    let mut dvalue: f64 = 0.0;

    if nullcheck == NullCheckType::None {
        /* no null checking required */
        if scale == 1.0 && zero == 0.0 {
            /* no scaling */
            for ii in 0..(ntodo as usize) {
                output[ii] = input[ii] as c_long; /* copy input to output */
            }
        } else {
            /* must scale the data */
            for ii in 0..(ntodo as usize) {
                dvalue = input[ii] as f64 * scale + zero;

                if dvalue < DLONG_MIN {
                    *status = OVERFLOW_ERR;
                    output[ii] = LONG_MIN;
                } else if dvalue > DLONG_MAX {
                    *status = OVERFLOW_ERR;
                    output[ii] = LONG_MAX;
                } else {
                    output[ii] = dvalue as c_long;
                }
            }
        }
    } else {
        /* must check for null values */
        let anynull = anynull.unwrap();

        if scale == 1.0 && zero == 0.0 {
            /* no scaling */
            for ii in 0..(ntodo as usize) {
                if input[ii] == tnull {
                    *anynull = 1;
                    if nullcheck == NullCheckType::SetPixel {
                        output[ii] = nullval;
                    } else {
                        nullarray[ii] = 1;
                    }
                } else {
                    output[ii] = input[ii] as c_long;
                }
            }
        } else {
            /* must scale the data */
            for ii in 0..(ntodo as usize) {
                if input[ii] == tnull {
                    *anynull = 1;
                    if nullcheck == NullCheckType::SetPixel {
                        output[ii] = nullval;
                    } else {
                        nullarray[ii] = 1;
                    }
                } else {
                    dvalue = input[ii] as f64 * scale + zero;

                    if dvalue < DLONG_MIN {
                        *status = OVERFLOW_ERR;
                        output[ii] = LONG_MIN;
                    } else if dvalue > DLONG_MAX {
                        *status = OVERFLOW_ERR;
                        output[ii] = LONG_MAX;
                    } else {
                        output[ii] = dvalue as c_long;
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
/// will = 0.  The anynull parameter will be set = 1 if any of the returned
/// pixels are null, otherwise anynull will be returned with a value = 0;
pub(crate) fn fffi8i4(
    input: &[LONGLONG],       /* I - array of values to be converted     */
    ntodo: c_long,            /* I - number of elements in the array     */
    scale: f64,               /* I - FITS TSCALn or BSCALE value         */
    zero: f64,                /* I - FITS TZEROn or BZERO  value         */
    nullcheck: NullCheckType, /* I - null checking code; 0 = don't check */
    /*     1:set null pixels = nullval         */
    /*     2: if null pixel, set nullarray = 1 */
    tnull: LONGLONG,             /* I - value of FITS TNULLn keyword if any */
    nullval: c_long,             /* I - set null pixels, if nullcheck = 1   */
    nullarray: &mut [c_char],    /* I - bad pixel array, if nullcheck = 2   */
    anynull: Option<&mut c_int>, /* O - set to 1 if any pixels are null     */
    output: &mut [c_long],       /* O - array of converted pixels           */
    status: &mut c_int,          /* IO - error status                       */
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
                ulltemp = ((input[ii] as ULONGLONG) ^ 0x8000000000000000) as ULONGLONG;

                if ulltemp > LONG_MAX as ULONGLONG {
                    *status = OVERFLOW_ERR;
                    output[ii] = LONG_MAX;
                } else {
                    output[ii] = ulltemp as c_long;
                }
            }
        } else if scale == 1.0 && zero == 0.0 {
            /* no scaling */
            for ii in 0..(ntodo as usize) {
                if input[ii] < LONG_MIN as LONGLONG {
                    *status = OVERFLOW_ERR;
                    output[ii] = LONG_MIN;
                } else if input[ii] > LONG_MAX as LONGLONG {
                    *status = OVERFLOW_ERR;
                    output[ii] = LONG_MAX;
                } else {
                    output[ii] = input[ii] as c_long;
                }
            }
        } else {
            /* must scale the data */
            for ii in 0..(ntodo as usize) {
                dvalue = input[ii] as f64 * scale + zero;

                if dvalue < DLONG_MIN {
                    *status = OVERFLOW_ERR;
                    output[ii] = LONG_MIN;
                } else if dvalue > DLONG_MAX {
                    *status = OVERFLOW_ERR;
                    output[ii] = LONG_MAX;
                } else {
                    output[ii] = dvalue as c_long;
                }
            }
        }
    } else {
        /* must check for null values */
        let anynull = anynull.unwrap();

        if scale == 1.0 && zero == 9223372036854775808. {
            /* The column we read contains unsigned long long values. */
            /* Instead of subtracting 9223372036854775808, it is more efficient */
            /* and more precise to just flip the sign bit with the XOR operator */

            for ii in 0..(ntodo as usize) {
                if input[ii] == tnull {
                    *anynull = 1;
                    if nullcheck == NullCheckType::SetPixel {
                        output[ii] = nullval;
                    } else {
                        nullarray[ii] = 1;
                    }
                } else {
                    ulltemp = ((input[ii] as ULONGLONG) ^ 0x8000000000000000) as ULONGLONG;

                    if ulltemp > LONG_MAX as ULONGLONG {
                        *status = OVERFLOW_ERR;
                        output[ii] = LONG_MAX;
                    } else {
                        output[ii] = ulltemp as c_long;
                    }
                }
            }
        } else if scale == 1.0 && zero == 0.0 {
            /* no scaling */
            for ii in 0..(ntodo as usize) {
                if input[ii] == tnull {
                    *anynull = 1;
                    if nullcheck == NullCheckType::SetPixel {
                        output[ii] = nullval;
                    } else {
                        nullarray[ii] = 1;
                    }
                } else if input[ii] < LONG_MIN as LONGLONG {
                    *status = OVERFLOW_ERR;
                    output[ii] = LONG_MIN;
                } else if input[ii] > LONG_MAX as LONGLONG {
                    *status = OVERFLOW_ERR;
                    output[ii] = LONG_MAX;
                } else {
                    output[ii] = input[ii] as c_long;
                }
            }
        } else {
            /* must scale the data */
            for ii in 0..(ntodo as usize) {
                if input[ii] == tnull {
                    *anynull = 1;
                    if nullcheck == NullCheckType::SetPixel {
                        output[ii] = nullval;
                    } else {
                        nullarray[ii] = 1;
                    }
                } else {
                    dvalue = input[ii] as f64 * scale + zero;

                    if dvalue < DLONG_MIN {
                        *status = OVERFLOW_ERR;
                        output[ii] = LONG_MIN;
                    } else if dvalue > DLONG_MAX {
                        *status = OVERFLOW_ERR;
                        output[ii] = LONG_MAX;
                    } else {
                        output[ii] = dvalue as c_long;
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
/// will = 0.  The anynull parameter will be set = 1 if any of the returned
/// pixels are null, otherwise anynull will be returned with a value = 0;
pub(crate) fn fffr4i4(
    input: &[f32],            /* I - array of values to be converted     */
    ntodo: c_long,            /* I - number of elements in the array     */
    scale: f64,               /* I - FITS TSCALn or BSCALE value         */
    zero: f64,                /* I - FITS TZEROn or BZERO  value         */
    nullcheck: NullCheckType, /* I - null checking code; 0 = don't check */
    /*     1:set null pixels = nullval         */
    /*     2: if null pixel, set nullarray = 1 */
    nullval: c_long,             /* I - set null pixels, if nullcheck = 1   */
    nullarray: &mut [c_char],    /* I - bad pixel array, if nullcheck = 2   */
    anynull: Option<&mut c_int>, /* O - set to 1 if any pixels are null     */
    output: &mut [c_long],       /* O - array of converted pixels           */
    status: &mut c_int,          /* IO - error status                       */
) -> c_int {
    let mut dvalue: f64 = 0.0;
    let mut sptr = 0;
    let iret = 0;

    if nullcheck == NullCheckType::None {
        /* no null checking required */
        if scale == 1.0 && zero == 0.0 {
            /* no scaling */
            for ii in 0..(ntodo as usize) {
                if (input[ii] as f64) < DLONG_MIN {
                    *status = OVERFLOW_ERR;
                    output[ii] = LONG_MIN;
                } else if input[ii] as f64 > DLONG_MAX {
                    *status = OVERFLOW_ERR;
                    output[ii] = LONG_MAX;
                } else {
                    output[ii] = input[ii] as c_long;
                }
            }
        } else {
            /* must scale the data */
            for ii in 0..(ntodo as usize) {
                dvalue = input[ii] as f64 * scale + zero;

                if dvalue < DLONG_MIN {
                    *status = OVERFLOW_ERR;
                    output[ii] = LONG_MIN;
                } else if dvalue > DLONG_MAX {
                    *status = OVERFLOW_ERR;
                    output[ii] = LONG_MAX;
                } else {
                    output[ii] = dvalue as c_long;
                }
            }
        }
    } else {
        /* must check for null values */
        let anynull = anynull.unwrap();

        //sptr = (short *) input;

        if BYTESWAPPED && CFITSIO_MACHINE != VAXVMS && CFITSIO_MACHINE != ALPHAVMS {
            sptr += 1; /* point to MSBs */
        }

        let shortBuff: &[c_short] = cast_slice(input);

        if scale == 1.0 && zero == 0.0 {
            /* no scaling */
            for ii in 0..(ntodo as usize) {
                let iret = fnan(shortBuff[sptr]);
                if iret != 0 {
                    /* test for NaN or underflow */
                    if iret == 1 {
                        /* is it a NaN? */
                        *anynull = 1;
                        if nullcheck == NullCheckType::SetPixel {
                            output[ii] = nullval;
                        } else {
                            nullarray[ii] = 1;
                        }
                    } else {
                        /* it's an underflow */
                        output[ii] = 0;
                    }
                } else if (input[ii] as f64) < DLONG_MIN {
                    *status = OVERFLOW_ERR;
                    output[ii] = LONG_MIN;
                } else if (input[ii] as f64) > DLONG_MAX {
                    *status = OVERFLOW_ERR;
                    output[ii] = LONG_MAX;
                } else {
                    output[ii] = input[ii] as c_long;
                }
                sptr += 2;
            }
        } else {
            /* must scale the data */
            for ii in 0..(ntodo as usize) {
                let iret = fnan(shortBuff[sptr]);
                if iret != 0 {
                    /* test for NaN or underflow */
                    if iret == 1 {
                        /* is it a NaN? */
                        *anynull = 1;
                        if nullcheck == NullCheckType::SetPixel {
                            output[ii] = nullval;
                        } else {
                            nullarray[ii] = 1;
                        }
                    } else {
                        /* it's an underflow */

                        if zero < DLONG_MIN {
                            *status = OVERFLOW_ERR;
                            output[ii] = LONG_MIN;
                        } else if zero > DLONG_MAX {
                            *status = OVERFLOW_ERR;
                            output[ii] = LONG_MAX;
                        } else {
                            output[ii] = zero as c_long;
                        }
                    }
                } else {
                    dvalue = input[ii] as f64 * scale + zero;

                    if dvalue < DLONG_MIN {
                        *status = OVERFLOW_ERR;
                        output[ii] = LONG_MIN;
                    } else if dvalue > DLONG_MAX {
                        *status = OVERFLOW_ERR;
                        output[ii] = LONG_MAX;
                    } else {
                        output[ii] = dvalue as c_long;
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
/// will = 0.  The anynull parameter will be set = 1 if any of the returned
/// pixels are null, otherwise anynull will be returned with a value = 0;
pub(crate) fn fffr8i4(
    input: &[f64],            /* I - array of values to be converted     */
    ntodo: c_long,            /* I - number of elements in the array     */
    scale: f64,               /* I - FITS TSCALn or BSCALE value         */
    zero: f64,                /* I - FITS TZEROn or BZERO  value         */
    nullcheck: NullCheckType, /* I - null checking code; 0 = don't check */
    /*     1:set null pixels = nullval         */
    /*     2: if null pixel, set nullarray = 1 */
    nullval: c_long,             /* I - set null pixels, if nullcheck = 1   */
    nullarray: &mut [c_char],    /* I - bad pixel array, if nullcheck = 2   */
    anynull: Option<&mut c_int>, /* O - set to 1 if any pixels are null     */
    output: &mut [c_long],       /* O - array of converted pixels           */
    status: &mut c_int,          /* IO - error status                       */
) -> c_int {
    let mut dvalue: f64 = 0.0;
    let mut sptr = 0;
    let iret = 0;

    if nullcheck == NullCheckType::None {
        /* no null checking required */
        if scale == 1.0 && zero == 0.0 {
            /* no scaling */
            for ii in 0..(ntodo as usize) {
                if input[ii] < DLONG_MIN {
                    *status = OVERFLOW_ERR;
                    output[ii] = LONG_MIN;
                } else if input[ii] > DLONG_MAX {
                    *status = OVERFLOW_ERR;
                    output[ii] = LONG_MAX;
                } else {
                    output[ii] = input[ii] as c_long;
                }
            }
        } else {
            /* must scale the data */
            for ii in 0..(ntodo as usize) {
                dvalue = input[ii] * scale + zero;

                if dvalue < DLONG_MIN {
                    *status = OVERFLOW_ERR;
                    output[ii] = LONG_MIN;
                } else if dvalue > DLONG_MAX {
                    *status = OVERFLOW_ERR;
                    output[ii] = LONG_MAX;
                } else {
                    output[ii] = dvalue as c_long;
                }
            }
        }
    } else {
        /* must check for null values */
        let anynull = anynull.unwrap();

        //sptr = (short *) input;

        if BYTESWAPPED && CFITSIO_MACHINE != VAXVMS && CFITSIO_MACHINE != ALPHAVMS {
            sptr += 3; /* point to MSBs */
        }

        let shortBuff: &[c_short] = cast_slice(input);

        if scale == 1.0 && zero == 0.0 {
            /* no scaling */
            for ii in 0..(ntodo as usize) {
                let iret = dnan(shortBuff[sptr]);
                if iret != 0 {
                    /* test for NaN or underflow */
                    if iret == 1 {
                        /* is it a NaN? */
                        *anynull = 1;
                        if nullcheck == NullCheckType::SetPixel {
                            output[ii] = nullval;
                        } else {
                            nullarray[ii] = 1;
                        }
                    } else {
                        /* it's an underflow */
                        output[ii] = 0;
                    }
                } else if input[ii] < DLONG_MIN {
                    *status = OVERFLOW_ERR;
                    output[ii] = LONG_MIN;
                } else if input[ii] > DLONG_MAX {
                    *status = OVERFLOW_ERR;
                    output[ii] = LONG_MAX;
                } else {
                    output[ii] = input[ii] as c_long;
                }
                sptr += 4;
            }
        } else {
            /* must scale the data */
            for ii in 0..(ntodo as usize) {
                let iret = dnan(shortBuff[sptr]);
                if iret != 0 {
                    /* test for NaN or underflow */
                    if iret == 1 {
                        /* is it a NaN? */
                        *anynull = 1;
                        if nullcheck == NullCheckType::SetPixel {
                            output[ii] = nullval;
                        } else {
                            nullarray[ii] = 1;
                        }
                    } else {
                        /* it's an underflow */

                        if zero < DLONG_MIN {
                            *status = OVERFLOW_ERR;
                            output[ii] = LONG_MIN;
                        } else if zero > DLONG_MAX {
                            *status = OVERFLOW_ERR;
                            output[ii] = LONG_MAX;
                        } else {
                            output[ii] = zero as c_long;
                        }
                    }
                } else {
                    dvalue = input[ii] * scale + zero;

                    if dvalue < DLONG_MIN {
                        *status = OVERFLOW_ERR;
                        output[ii] = LONG_MIN;
                    } else if dvalue > DLONG_MAX {
                        *status = OVERFLOW_ERR;
                        output[ii] = LONG_MAX;
                    } else {
                        output[ii] = dvalue as c_long;
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
/// will = 0.  The anynull parameter will be set = 1 if any of the returned
/// pixels are null, otherwise anynull will be returned with a value = 0;
pub(crate) fn fffstri4(
    input: &mut [c_char],     /* I - array of values to be converted     */
    ntodo: c_long,            /* I - number of elements in the array     */
    scale: f64,               /* I - FITS TSCALn or BSCALE value         */
    zero: f64,                /* I - FITS TZEROn or BZERO  value         */
    twidth: c_long,           /* I - width of each substring of chars    */
    implipower: f64,          /* I - power of 10 of implied decimal      */
    nullcheck: NullCheckType, /* I - null checking code; 0 = don't check */
    /*     1:set null pixels = nullval         */
    /*     2: if null pixel, set nullarray = 1 */
    snull: &[c_char],                /* I - value of FITS null string, if any   */
    nullval: c_long,                 /* I - set null pixels, if nullcheck = 1   */
    nullarray: &mut [c_char],        /* I - bad pixel array, if nullcheck = 2   */
    mut anynull: Option<&mut c_int>, /* O - set to 1 if any pixels are null     */
    output: &mut [c_long],           /* O - array of converted pixels           */
    status: &mut c_int,              /* IO - error status                       */
) -> c_int {
    let mut nullen = 0;

    let mut dvalue: f64 = 0.0;

    let mut message: [c_char; FLEN_ERRMSG] = [0; FLEN_ERRMSG];
    let mut tempstore: c_char = 0;
    let chrzero: c_char = bb(b'0'); // 49
    let mut val: f64 = 0.0;
    let mut power: f64 = 0.0;
    let mut exponent: c_int = 0;
    let mut sign: c_int = 0;
    let mut esign: c_int = 0;
    let mut decpt: c_int = 0;

    nullen = strlen_safe(snull);
    let mut cptr = 0; /* pointer to start of input string */
    let mut tpos = 0;
    for ii in 0..(ntodo as usize) {
        let originalCptr = cptr;
        /* temporarily insert a null terminator at end of the string */
        tpos = cptr + twidth as usize;
        tempstore = input[tpos];
        input[tpos] = 0;

        /* check if null value is defined, and if the    */
        /* column string is identical to the null string */
        if snull[0] != ASCII_NULL_UNDEFINED && strncmp_safe(snull, &input[cptr..], nullen) == 0 {
            if nullcheck != NullCheckType::None {
                if let Some(anynull) = anynull.as_deref_mut() {
                    *anynull = 1;
                }
                if nullcheck == NullCheckType::SetPixel {
                    output[ii] = nullval;
                } else {
                    nullarray[ii] = 1;
                }
            }
            cptr += twidth as usize;
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
                val = val * 10.0 + input[cptr] as f64 - chrzero as f64; /* accumulate the value */
                cptr += 1;

                while input[cptr] == bb(b' ') {
                    /* skip embedded blanks in the value */
                    cptr += 1;
                }
            }

            if input[cptr] == bb(b'.') || input[cptr] == bb(b',') {
                /* check for decimal point */
                decpt = 1; /* set flag to show there was a decimal point */
                cptr += 1;
                while input[cptr] == bb(b' ') {
                    /* skip any blanks */
                    cptr += 1;
                }
                while input[cptr] >= bb(b'0') && input[cptr] <= bb(b'9') {
                    val = val * 10.0 + input[cptr] as f64 - chrzero as f64; /* accumulate the value */
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
                    exponent = exponent * 10 + input[cptr] as c_int - chrzero as c_int; /* accumulate exp */
                    cptr += 1;

                    while input[cptr] == bb(b' ') {
                        /* skip embedded blanks */
                        cptr += 1;
                    }
                }
            }

            if input[cptr] != 0 {
                /* should end up at the null terminator */
                let cstring = &input[originalCptr..];
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
                    slice_to_str!(cstring),
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

            if dvalue < DLONG_MIN {
                *status = OVERFLOW_ERR;
                output[ii] = LONG_MIN;
            } else if dvalue > DLONG_MAX {
                *status = OVERFLOW_ERR;
                output[ii] = LONG_MAX;
            } else {
                output[ii] = dvalue as c_long;
            }
        }
        /* restore the char that was overwritten by the null */
        input[tpos] = tempstore;
    }
    *status
}

/* ======================================================================== */
/*      the following routines support the 'long long' data type            */
/* ======================================================================== */

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
pub unsafe extern "C" fn ffgpvjj(
    fptr: *mut fitsfile,  /* I - FITS file pointer                       */
    group: c_long,        /* I - group to read (1 = 1st group)           */
    firstelem: LONGLONG,  /* I - first vector element to read (1 = 1st)  */
    nelem: LONGLONG,      /* I - number of values to read                */
    nulval: LONGLONG,     /* I - value for undefined pixels              */
    array: *mut LONGLONG, /* O - array of values that are returned       */
    anynul: *mut c_int,   /* O - set to 1 if any values are null; else 0 */
    status: *mut c_int,   /* IO - error status                           */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let anynul = anynul.as_mut();

        let array = slice::from_raw_parts_mut(array, nelem as usize);

        ffgpvjj_safe(fptr, group, firstelem, nelem, nulval, array, anynul, status)
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
pub fn ffgpvjj_safe(
    fptr: &mut fitsfile,        /* I - FITS file pointer                       */
    group: c_long,              /* I - group to read (1 = 1st group)           */
    firstelem: LONGLONG,        /* I - first vector element to read (1 = 1st)  */
    nelem: LONGLONG,            /* I - number of values to read                */
    nulval: LONGLONG,           /* I - value for undefined pixels              */
    array: &mut [LONGLONG],     /* O - array of values that are returned       */
    anynul: Option<&mut c_int>, /* O - set to 1 if any values are null; else 0 */
    status: &mut c_int,         /* IO - error status                           */
) -> c_int {
    let mut row = 0;
    let nullcheck = NullCheckType::SetPixel;
    let mut nullvalue: LONGLONG = 0;

    let mut dummy_nularray = vec![0; (nelem) as usize];

    if fits_is_compressed_image_safe(fptr, status) > 0 {
        /* this is a compressed image in a binary table */
        nullvalue = nulval; /* set local variable */

        todo!();
        /*
        fits_read_compressed_pixels(
            fptr, TLONGLONG, firstelem, nelem, nullcheck, &nullvalue, cast_slice_mut(array), None, anynul, status,
        );
        */
        return *status;
    }

    row = cmp::max(1, group);

    ffgcljj(
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
pub unsafe extern "C" fn ffgpfjj(
    fptr: *mut fitsfile,   /* I - FITS file pointer                       */
    group: c_long,         /* I - group to read (1 = 1st group)           */
    firstelem: LONGLONG,   /* I - first vector element to read (1 = 1st)  */
    nelem: LONGLONG,       /* I - number of values to read                */
    array: *mut LONGLONG,  /* O - array of values that are returned       */
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

        ffgpfjj_safe(
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
pub fn ffgpfjj_safe(
    fptr: &mut fitsfile,        /* I - FITS file pointer                       */
    group: c_long,              /* I - group to read (1 = 1st group)           */
    firstelem: LONGLONG,        /* I - first vector element to read (1 = 1st)  */
    nelem: LONGLONG,            /* I - number of values to read                */
    array: &mut [LONGLONG],     /* O - array of values that are returned       */
    nularray: &mut [c_char],    /* O - array of null pixel flags               */
    anynul: Option<&mut c_int>, /* O - set to 1 if any values are null; else 0 */
    status: &mut c_int,         /* IO - error status                           */
) -> c_int {
    let mut row = 0;
    let nullcheck = NullCheckType::SetNullArray;
    let dummy: LONGLONG = 0;

    if fits_is_compressed_image_safe(fptr, status) > 0 {
        /* this is a compressed image in a binary table */
        todo!();
        /*
        fits_read_compressed_pixels(
            fptr, TLONGLONG, firstelem, nelem, nullcheck, None, cast_slice_mut(array), nularray, anynul, status,
        );
        */
        return *status;
    }

    row = cmp::max(1, group);

    ffgcljj(
        fptr,
        2,
        row as LONGLONG,
        firstelem,
        nelem,
        1,
        NullCheckType::SetNullArray,
        dummy,
        array,
        nularray,
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
pub unsafe extern "C" fn ffg2djj(
    fptr: *mut fitsfile,  /* I - FITS file pointer                       */
    group: c_long,        /* I - group to read (1 = 1st group)           */
    nulval: LONGLONG,     /* set undefined pixels equal to this          */
    ncols: LONGLONG,      /* I - number of pixels in each row of array   */
    naxis1: LONGLONG,     /* I - FITS image NAXIS1 value                 */
    naxis2: LONGLONG,     /* I - FITS image NAXIS2 value                 */
    array: *mut LONGLONG, /* O - array to be filled and returned         */
    anynul: *mut c_int,   /* O - set to 1 if any values are null; else 0 */
    status: *mut c_int,   /* IO - error status                           */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let anynul = anynul.as_mut();

        let array = slice::from_raw_parts_mut(array, (naxis1 * naxis2) as usize);

        ffg2djj_safe(
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
pub fn ffg2djj_safe(
    fptr: &mut fitsfile,        /* I - FITS file pointer                       */
    group: c_long,              /* I - group to read (1 = 1st group)           */
    nulval: LONGLONG,           /* set undefined pixels equal to this     */
    ncols: LONGLONG,            /* I - number of pixels in each row of array   */
    naxis1: LONGLONG,           /* I - FITS image NAXIS1 value                 */
    naxis2: LONGLONG,           /* I - FITS image NAXIS2 value                 */
    array: &mut [LONGLONG],     /* O - array to be filled and returned    */
    anynul: Option<&mut c_int>, /* O - set to 1 if any values are null; else 0 */
    status: &mut c_int,         /* IO - error status                           */
) -> c_int {
    /* call the 3D reading routine, with the 3rd dimension = 1 */
    ffg3djj_safe(
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
pub unsafe extern "C" fn ffg3djj(
    fptr: *mut fitsfile,  /* I - FITS file pointer                       */
    group: c_long,        /* I - group to read (1 = 1st group)           */
    nulval: LONGLONG,     /* set undefined pixels equal to this          */
    ncols: LONGLONG,      /* I - number of pixels in each row of array   */
    nrows: LONGLONG,      /* I - number of rows in each plane of array   */
    naxis1: LONGLONG,     /* I - FITS image NAXIS1 value                 */
    naxis2: LONGLONG,     /* I - FITS image NAXIS2 value                 */
    naxis3: LONGLONG,     /* I - FITS image NAXIS3 value                 */
    array: *mut LONGLONG, /* O - array to be filled and returned         */
    anynul: *mut c_int,   /* O - set to 1 if any values are null; else 0 */
    status: *mut c_int,   /* IO - error status                           */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let anynul = anynul.as_mut();

        let array = slice::from_raw_parts_mut(array, (ncols * nrows) as usize);

        ffg3djj_safe(
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
pub fn ffg3djj_safe(
    fptr: &mut fitsfile,            /* I - FITS file pointer                       */
    group: c_long,                  /* I - group to read (1 = 1st group)           */
    nulval: LONGLONG,               /* set undefined pixels equal to this          */
    ncols: LONGLONG,                /* I - number of pixels in each row of array   */
    nrows: LONGLONG,                /* I - number of rows in each plane of array   */
    naxis1: LONGLONG,               /* I - FITS image NAXIS1 value                 */
    naxis2: LONGLONG,               /* I - FITS image NAXIS2 value                 */
    naxis3: LONGLONG,               /* I - FITS image NAXIS3 value                 */
    array: &mut [LONGLONG],         /* O - array to be filled and returned         */
    mut anynul: Option<&mut c_int>, /* O - set to 1 if any values are null; else 0 */
    status: &mut c_int,             /* IO - error status                           */
) -> c_int {
    let mut tablerow: c_long = 0;
    let ii: c_long = 0;
    let jj: c_long = 0;
    let nullcheck = NullCheckType::SetPixel;
    let inc: [c_long; 3] = [1, 1, 1];

    let fpixel: [LONGLONG; 3] = [1, 1, 1];
    let mut nfits: LONGLONG = 0;
    let mut narray: LONGLONG = 0;
    let mut lpixel: [LONGLONG; 3] = [0; 3];
    let mut nullvalue: LONGLONG = 0;

    let mut dummy_nularray = vec![0; (ncols * nrows) as usize];

    if fits_is_compressed_image_safe(fptr, status) > 0 {
        /* this is a compressed image in a binary table */

        lpixel[0] = ncols as LONGLONG;
        lpixel[1] = nrows as LONGLONG;
        lpixel[2] = naxis3 as LONGLONG;
        nullvalue = nulval; /* set local variable */

        todo!();
        /*
        fits_read_compressed_img(
            fptr, TLONGLONG, fpixel, lpixel, inc, nullcheck, &nullvalue, cast_slice_mut(array), None, anynul,
            status,
        );
        */
        return *status;
    }

    tablerow = cmp::max(1, group);

    if ncols == naxis1 && nrows == naxis2 {
        /* arrays have same size? */
        /* all the image pixels are contiguous, so read all at once */
        ffgcljj(
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
            if ffgcljj(
                fptr,
                2,
                tablerow as LONGLONG,
                nfits,
                naxis1,
                1,
                NullCheckType::SetPixel,
                nulval,
                &mut array[(narray as usize)..],
                &mut dummy_nularray,
                anynul.as_deref_mut(),
                status,
            ) > 0
            {
                return *status;
            }

            nfits += naxis1;
            narray += ncols;
        }
        narray += (nrows - naxis2) * ncols;
    }

    *status
}

/*--------------------------------------------------------------------------*/
/// Read a subsection of data values from an image or a table column.
/// This routine is set up to handle a maximum of nine dimensions.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgsvjj(
    fptr: *mut fitsfile,  /* I - FITS file pointer                         */
    colnum: c_int,        /* I - number of the column to read (1 = 1st)    */
    naxis: c_int,         /* I - number of dimensions in the FITS array    */
    naxes: *const c_long, /* I - size of each dimension                    */
    blc: *const c_long,   /* I - 'bottom left corner' of the subsection    */
    trc: *const c_long,   /* I - 'top right corner' of the subsection      */
    inc: *const c_long,   /* I - increment to be applied in each dimension */
    nulval: LONGLONG,     /* I - value to set undefined pixels             */
    array: *mut LONGLONG, /* O - array to be filled and returned           */
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

        let nelem = calculate_subsection_length(blc, trc, inc);

        let array = slice::from_raw_parts_mut(array, nelem as usize);

        ffgsvjj_safe(
            fptr, colnum, naxis, naxes, blc, trc, inc, nulval, array, anynul, status,
        )
    }
}

/*--------------------------------------------------------------------------*/
/// Read a subsection of data values from an image or a table column.
/// This routine is set up to handle a maximum of nine dimensions.
pub fn ffgsvjj_safe(
    fptr: &mut fitsfile,    /* I - FITS file pointer                         */
    colnum: c_int,          /* I - number of the column to read (1 = 1st)    */
    naxis: c_int,           /* I - number of dimensions in the FITS array    */
    naxes: &[c_long],       /* I - size of each dimension                    */
    blc: &[c_long],         /* I - 'bottom left corner' of the subsection    */
    trc: &[c_long],         /* I - 'top right corner' of the subsection      */
    inc: &[c_long],         /* I - increment to be applied in each dimension */
    nulval: LONGLONG,       /* I - value to set undefined pixels             */
    array: &mut [LONGLONG], /* O - array to be filled and returned           */
    mut anynul: Option<&mut c_int>, /* O - set to 1 if any values are null; else 0   */
    status: &mut c_int,     /* IO - error status                             */
) -> c_int {
    let ii: c_long = 0;
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
    let ldummy: c_char = 0;
    let mut msg: [c_char; FLEN_ERRMSG] = [0; FLEN_ERRMSG];
    let nullcheck = NullCheckType::SetPixel;
    let mut nullvalue: LONGLONG = 0;

    let naxis: usize = naxis as usize;

    if naxis < 1 || naxis > 9 {
        int_snprintf!(
            &mut msg,
            FLEN_ERRMSG,
            "NAXIS = {} in call to ffgsvj is out of range",
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
        /*
        fits_read_compressed_img(
            fptr, TLONGLONG, blcll, trcll, inc, nullcheck, &nullvalue, cast_slice_mut(array), None, anynul, status,
        );
         */
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
                    "ffgsvj: illegal range specified for axis {}",
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

                                        if ffgcljj(
                                            fptr,
                                            numcol as c_int,
                                            row as LONGLONG,
                                            felem,
                                            nelem as LONGLONG,
                                            ninc,
                                            nultyp,
                                            nulval,
                                            cast_slice_mut(&mut array[(i0 as usize)..]),
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
                                        i0 += nelem;
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

/*--------------------------------------------------------------------------*/
/// Read a subsection of data values from an image or a table column.
/// This routine is set up to handle a maximum of nine dimensions.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgsfjj(
    fptr: *mut fitsfile,  /* I - FITS file pointer                         */
    colnum: c_int,        /* I - number of the column to read (1 = 1st)    */
    naxis: c_int,         /* I - number of dimensions in the FITS array    */
    naxes: *const c_long, /* I - size of each dimension                    */
    blc: *const c_long,   /* I - 'bottom left corner' of the subsection    */
    trc: *const c_long,   /* I - 'top right corner' of the subsection      */
    inc: *const c_long,   /* I - increment to be applied in each dimension */
    array: *mut LONGLONG, /* O - array to be filled and returned           */
    flagval: *mut c_char, /* O - set to 1 if corresponding value is null   */
    anynul: *mut c_int,   /* O - set to 1 if any values are null; else 0   */
    status: *mut c_int,   /* IO - error status                             */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let naxes = slice::from_raw_parts(naxes, naxis as usize);
        let blc = slice::from_raw_parts(blc, naxis as usize);
        let trc = slice::from_raw_parts(trc, naxis as usize);
        let inc = slice::from_raw_parts(inc, naxis as usize);

        let total_nelem = calculate_subsection_length(blc, trc, inc);
        let array = slice::from_raw_parts_mut(array, total_nelem);
        let flagval = slice::from_raw_parts_mut(flagval, total_nelem);

        let anynul = anynul.as_mut();
        ffgsfjj_safe(
            fptr, colnum, naxis, naxes, blc, trc, inc, array, flagval, anynul, status,
        )
    }
}

/*--------------------------------------------------------------------------*/
/// Read a subsection of data values from an image or a table column.
/// This routine is set up to handle a maximum of nine dimensions.
pub fn ffgsfjj_safe(
    fptr: &mut fitsfile,    /* I - FITS file pointer                         */
    colnum: c_int,          /* I - number of the column to read (1 = 1st)    */
    naxis: c_int,           /* I - number of dimensions in the FITS array    */
    naxes: &[c_long],       /* I - size of each dimension                    */
    blc: &[c_long],         /* I - 'bottom left corner' of the subsection    */
    trc: &[c_long],         /* I - 'top right corner' of the subsection      */
    inc: &[c_long],         /* I - increment to be applied in each dimension */
    array: &mut [LONGLONG], /* O - array to be filled and returned           */
    flagval: &mut [c_char], /* O - set to 1 if corresponding value is null   */
    mut anynul: Option<&mut c_int>, /* O - set to 1 if any values are null; else 0   */
    status: &mut c_int,     /* IO - error status                             */
) -> c_int {
    let ii: c_long = 0;
    let i0: c_long = 0;
    let row: c_long = 0;
    let mut rstr: c_long = 0;
    let mut rstp: c_long = 0;
    let mut rinc: c_long = 0;
    let mut str: [c_long; 9] = [0; 9];
    let mut stp: [c_long; 9] = [0; 9];
    let mut incr: [c_long; 9] = [0; 9];
    let mut dsize: [c_long; 10] = [0; 10];
    let mut blcll: [LONGLONG; 9] = [0; 9];
    let mut trcll: [LONGLONG; 9] = [0; 9];
    let mut felem: c_long = 0;
    let mut nelem: c_long = 0;
    let mut nultyp = NullCheckType::None;
    let mut ninc: c_long = 0;
    let mut numcol: c_long = 0;
    let nulval: LONGLONG = 0;
    let mut hdutype: c_int = 0;
    let mut anyf: c_int = 0;
    let mut msg: [c_char; FLEN_ERRMSG] = [0; FLEN_ERRMSG];
    let nullcheck = NullCheckType::SetNullArray;

    let naxis: usize = naxis as usize;

    if naxis < 1 || naxis > 9 {
        int_snprintf!(
            &mut msg,
            FLEN_ERRMSG,
            "NAXIS = {} in call to ffgsvj is out of range",
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
        /*
        fits_read_compressed_img(
            fptr, TLONGLONG, blcll, trcll, inc, nullcheck, None, cast_slice_mut(array), flagval, anynul, status,
        );
        */
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
            int_snprintf!(
                &mut msg,
                FLEN_ERRMSG,
                "ffgsvj: illegal range specified for axis {}",
                ii + 1,
            );
            ffpmsg_slice(&msg);
            *status = BAD_PIX_NUM;
            return *status;
        }

        str[ii] = blc[ii];
        stp[ii] = trc[ii];
        incr[ii] = inc[ii];
        dsize[ii + 1] = dsize[ii] * naxes[ii] as c_long;
    }

    if naxis == 1 && naxes[0] == 1 {
        /* This is not a vector column, so read all the rows at once */
        nelem = (rstp - rstr) / rinc + 1;
        ninc = rinc;
        rstp = rstr;
    } else {
        /* have to read each row individually, in all dimensions */
        nelem = (stp[0] - str[0]) / inc[0] + 1;
        ninc = incr[0];
    }

    for row in (rstr..=rstp).step_by(rinc as usize) {
        for i8 in ((str[8])..=(stp[8])).step_by(incr[8] as usize) {
            for i7 in ((str[7])..=(stp[7])).step_by(incr[7] as usize) {
                for i6 in ((str[6])..=(stp[6])).step_by(incr[6] as usize) {
                    for i5 in ((str[5])..=(stp[5])).step_by(incr[5] as usize) {
                        for i4 in ((str[4])..=(stp[4])).step_by(incr[4] as usize) {
                            for i3 in ((str[3])..=(stp[3])).step_by(incr[3] as usize) {
                                for i2 in ((str[2])..=(stp[2])).step_by(incr[2] as usize) {
                                    for i1 in ((str[1])..=(stp[1])).step_by(incr[1] as usize) {
                                        felem = (str[0] as c_long)
                                            + (i1 as c_long - 1) * dsize[1]
                                            + (i2 as c_long - 1) * dsize[2]
                                            + (i3 as c_long - 1) * dsize[3]
                                            + (i4 as c_long - 1) * dsize[4]
                                            + (i5 as c_long - 1) * dsize[5]
                                            + (i6 as c_long - 1) * dsize[6]
                                            + (i7 as c_long - 1) * dsize[7]
                                            + (i8 as c_long - 1) * dsize[8];

                                        if ffgcljj(
                                            fptr,
                                            numcol as c_int,
                                            row as LONGLONG,
                                            felem as LONGLONG,
                                            nelem as LONGLONG,
                                            ninc,
                                            nultyp,
                                            nulval,
                                            &mut array[(i0 as usize)..],
                                            &mut flagval[(i0 as usize)..],
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
                                        i0 += nelem;
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

/*--------------------------------------------------------------------------*/
/// Read an array of group parameters from the primary array. Data conversion
/// and scaling will be performed if necessary (e.g, if the datatype of
/// the FITS array is not the same as the array being read).
///
/// The primary array is represented as a binary table:
/// each group of the primary array is a row in the table,
/// where the first column contains the group parameters
/// and the second column contains the image itself.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffggpjj(
    fptr: *mut fitsfile,  /* I - FITS file pointer                       */
    group: c_long,        /* I - group to read (1 = 1st group)           */
    firstelem: c_long,    /* I - first vector element to read (1 = 1st)  */
    nelem: c_long,        /* I - number of values to read                */
    array: *mut LONGLONG, /* O - array of values that are returned       */
    status: *mut c_int,   /* IO - error status                           */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let array = slice::from_raw_parts_mut(array, nelem as usize);

        ffggpjj_safe(fptr, group, firstelem, nelem, array, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Read an array of group parameters from the primary array. Data conversion
/// and scaling will be performed if necessary (e.g, if the datatype of
/// the FITS array is not the same as the array being read).
///
/// The primary array is represented as a binary table:
/// each group of the primary array is a row in the table,
/// where the first column contains the group parameters
/// and the second column contains the image itself.
pub fn ffggpjj_safe(
    fptr: &mut fitsfile,    /* I - FITS file pointer                       */
    group: c_long,          /* I - group to read (1 = 1st group)           */
    firstelem: c_long,      /* I - first vector element to read (1 = 1st)  */
    nelem: c_long,          /* I - number of values to read                */
    array: &mut [LONGLONG], /* O - array of values that are returned       */
    status: &mut c_int,     /* IO - error status                           */
) -> c_int {
    let mut dummy_nularray = vec![0; nelem as usize];
    let mut idummy = 0;
    let dummy: LONGLONG = 0;

    let row = cmp::max(1, group);

    ffgcljj(
        fptr,
        2,
        row as LONGLONG,
        firstelem as LONGLONG,
        nelem as LONGLONG,
        1,
        NullCheckType::SetPixel,
        dummy,
        array,
        &mut dummy_nularray,
        Some(&mut idummy),
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
pub unsafe extern "C" fn ffgcvjj(
    fptr: *mut fitsfile,  /* I - FITS file pointer                       */
    colnum: c_int,        /* I - number of column to read (1 = 1st col)  */
    firstrow: LONGLONG,   /* I - first row to read (1 = 1st row)         */
    firstelem: LONGLONG,  /* I - first vector element to read (1 = 1st)  */
    nelem: LONGLONG,      /* I - number of values to read                */
    nulval: LONGLONG,     /* I - value for null pixels                   */
    array: *mut LONGLONG, /* O - array of values that are read           */
    anynul: *mut c_int,   /* O - set to 1 if any values are null; else 0 */
    status: *mut c_int,   /* IO - error status                           */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let anynul = anynul.as_mut();

        let array = slice::from_raw_parts_mut(array, nelem as usize);

        ffgcvjj_safe(
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
pub fn ffgcvjj_safe(
    fptr: &mut fitsfile,        /* I - FITS file pointer                       */
    colnum: c_int,              /* I - number of column to read (1 = 1st col)  */
    firstrow: LONGLONG,         /* I - first row to read (1 = 1st row)         */
    firstelem: LONGLONG,        /* I - first vector element to read (1 = 1st)  */
    nelem: LONGLONG,            /* I - number of values to read                */
    nulval: LONGLONG,           /* I - value for null pixels                   */
    array: &mut [LONGLONG],     /* O - array of values that are read           */
    anynul: Option<&mut c_int>, /* O - set to 1 if any values are null; else 0 */
    status: &mut c_int,         /* IO - error status                           */
) -> c_int {
    let mut dummy_nularray = vec![0; (nelem) as usize];

    ffgcljj(
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
pub unsafe extern "C" fn ffgcfjj(
    fptr: *mut fitsfile,   /* I - FITS file pointer                       */
    colnum: c_int,         /* I - number of column to read (1 = 1st col)  */
    firstrow: LONGLONG,    /* I - first row to read (1 = 1st row)         */
    firstelem: LONGLONG,   /* I - first vector element to read (1 = 1st)  */
    nelem: LONGLONG,       /* I - number of values to read                */
    array: *mut LONGLONG,  /* O - array of values that are read           */
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

        ffgcfjj_safe(
            fptr,
            colnum,
            firstrow as LONGLONG,
            firstelem as LONGLONG,
            nelem as LONGLONG,
            array,
            nularray,
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
/// Nularray will be set = 1 if the corresponding array pixel is undefined,
/// otherwise nularray will = 0.
pub fn ffgcfjj_safe(
    fptr: &mut fitsfile,        /* I - FITS file pointer                       */
    colnum: c_int,              /* I - number of column to read (1 = 1st col)  */
    firstrow: LONGLONG,         /* I - first row to read (1 = 1st row)         */
    firstelem: LONGLONG,        /* I - first vector element to read (1 = 1st)  */
    nelem: LONGLONG,            /* I - number of values to read                */
    array: &mut [LONGLONG],     /* O - array of values that are read           */
    nularray: &mut [c_char],    /* O - array of flags: 1 if null pixel; else 0 */
    anynul: Option<&mut c_int>, /* O - set to 1 if any values are null; else 0 */
    status: &mut c_int,         /* IO - error status                           */
) -> c_int {
    let dummy: LONGLONG = 0;

    ffgcljj(
        fptr,
        colnum,
        firstrow,
        firstelem,
        nelem,
        1,
        NullCheckType::SetNullArray,
        dummy,
        array,
        nularray,
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
pub(crate) fn ffgcljj(
    fptr: &mut fitsfile,   /* I - FITS file pointer                       */
    colnum: c_int,         /* I - number of column to read (1 = 1st col)  */
    firstrow: LONGLONG,    /* I - first row to read (1 = 1st row)         */
    firstelem: LONGLONG,   /* I - first vector element to read (1 = 1st)  */
    nelem: LONGLONG,       /* I - number of values to read                */
    elemincre: c_long,     /* I - pixel increment; e.g., 2 = every other  */
    nultyp: NullCheckType, /* I - null value handling code:               */
    /*     1: set undefined pixels = nulval        */
    /*     2: set nularray=1 for undefined pixels  */
    nulval: LONGLONG,               /* I - value for null pixels if nultyp = 1     */
    array: &mut [LONGLONG],         /* O - array of values that are read           */
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
    let ii: c_long = 0;
    let mut xwidth: c_long = 0;
    let mut ntodo: c_long = 0;
    let mut convert: bool = false;
    let mut nulcheck = NullCheckType::None;
    let mut readcheck: c_int = 0;
    let mut repeat: LONGLONG = 0;
    let mut startpos: LONGLONG = 0;
    let mut elemnum: LONGLONG = 0;
    let mut readptr: LONGLONG = 0;
    let mut tnull: LONGLONG = 0;
    let mut rowlen: LONGLONG = 0;
    let mut rownum: LONGLONG = 0;
    let mut remain: LONGLONG = 0;
    let mut next: LONGLONG = 0;
    let mut rowincre: LONGLONG = 0;
    let mut maxelem: LONGLONG = 0;
    let mut tform: [c_char; 20] = [0; 20];
    let mut message: [c_char; FLEN_ERRMSG] = [0; FLEN_ERRMSG];
    let mut snull: [c_char; 20] = [0; 20]; /*  the FITS null value if reading from ASCII table  */
    let mut cbuff: [f64; DBUFFSIZE as usize / mem::size_of::<f64>()] =
        [0.0; DBUFFSIZE as usize / mem::size_of::<f64>()]; /* align cbuff on word boundary */

    if *status > 0 || nelem == 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    let buffer: &mut [u8] = cast_slice_mut(&mut cbuff);

    if let Some(anynul) = anynul.as_deref_mut() {
        *anynul = 0;
    }

    if nultyp == NullCheckType::SetNullArray {
        nularray.fill(0); /* initialize nullarray */
    }

    /*---------------------------------------------------*/
    /*  Check input and get parameters about the column: */
    /*---------------------------------------------------*/
    if elemincre < 0 {
        readcheck = -1; /* don't do range checking in this case */
    }
    if ffgcprll(
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
    ) > 0
    {
        return *status;
    }
    maxelem = maxelem2 as LONGLONG;

    incre *= elemincre; /* multiply incre to just get every nth pixel */

    if tcode == TSTRING {
        /* setup for ASCII tables */
        /* get the number of implied decimal places if no explicit decmal point */
        ffasfm_safe(
            &tform,
            Some(&mut xcode),
            Some(&mut xwidth),
            Some(&mut decimals),
            status,
        );

        for ii in 0..(decimals as usize) {
            power *= 10.;
        }
    }

    /*------------------------------------------------------------------*/
    /*  Decide whether to check for null values in the input FITS file: */
    /*------------------------------------------------------------------*/
    nulcheck = nultyp; /* by default check for null values in the FITS file */

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
    if tcode == TLONGLONG {
        /* Special Case:                        */
        /* no type convertion required, so read */
        /* data directly into output buffer.    */

        if nelem < INT32_MAX as LONGLONG / 8 {
            maxelem = nelem;
        } else {
            maxelem = INT32_MAX as LONGLONG / 8;
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

    while remain > 0 {
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

        readptr =
            startpos + (rownum as LONGLONG * rowlen) + (elemnum * (incre / elemincre) as LONGLONG);

        match tcode {
            TLONGLONG => {
                if convert {
                    ffgi8b(fptr, readptr, ntodo, incre, cast_slice_mut(buffer), status);

                    fffi8i8(
                        cast_slice_mut(buffer),
                        ntodo,
                        scale,
                        zero,
                        nulcheck,
                        tnull,
                        nulval,
                        &mut nularray[(next as usize)..],
                        anynul.as_deref_mut(),
                        cast_slice_mut(&mut array[(next as usize)..]),
                        status,
                    );
                } else {
                    ffgi8b(
                        fptr,
                        readptr,
                        ntodo,
                        incre,
                        cast_slice_mut(&mut array[(next as usize)..]),
                        status,
                    );
                }
            }
            TLONG => {
                ffgi4b(fptr, readptr, ntodo, incre, cast_slice_mut(buffer), status);

                fffi4i8(
                    cast_slice(buffer),
                    ntodo,
                    scale,
                    zero,
                    nulcheck,
                    tnull as INT32BIT,
                    nulval,
                    &mut nularray[(next as usize)..],
                    anynul.as_deref_mut(),
                    cast_slice_mut(&mut array[(next as usize)..]),
                    status,
                );
            }
            TBYTE => {
                ffgi1b(fptr, readptr, ntodo, incre, buffer, status);

                fffi1i8(
                    cast_slice(buffer),
                    ntodo,
                    scale,
                    zero,
                    nulcheck,
                    tnull as u8,
                    nulval,
                    cast_slice_mut(&mut nularray[(next as usize)..]),
                    anynul.as_deref_mut(),
                    cast_slice_mut(&mut array[(next as usize)..]),
                    status,
                );
            }
            TSHORT => {
                ffgi2b(fptr, readptr, ntodo, incre, cast_slice_mut(buffer), status);

                fffi2i8(
                    cast_slice(buffer),
                    ntodo,
                    scale,
                    zero,
                    nulcheck,
                    tnull as c_short,
                    nulval,
                    cast_slice_mut(&mut nularray[(next as usize)..]),
                    anynul.as_deref_mut(),
                    cast_slice_mut(&mut array[(next as usize)..]),
                    status,
                );
            }
            TFLOAT => {
                ffgr4b(fptr, readptr, ntodo, incre, cast_slice_mut(buffer), status);

                fffr4i8(
                    cast_slice(buffer),
                    ntodo,
                    scale,
                    zero,
                    nulcheck,
                    nulval,
                    cast_slice_mut(&mut nularray[(next as usize)..]),
                    anynul.as_deref_mut(),
                    cast_slice_mut(&mut array[(next as usize)..]),
                    status,
                );
            }
            TDOUBLE => {
                ffgr8b(fptr, readptr, ntodo, incre, cast_slice_mut(buffer), status);

                fffr8i8(
                    cast_slice(buffer),
                    ntodo,
                    scale,
                    zero,
                    nulcheck,
                    nulval,
                    cast_slice_mut(&mut nularray[(next as usize)..]),
                    anynul.as_deref_mut(),
                    cast_slice_mut(&mut array[(next as usize)..]),
                    status,
                );
            }
            TSTRING => {
                ffmbyt_safe(fptr, readptr, REPORT_EOF, status);

                if incre == twidth {
                    /* contiguous bytes */
                    ffgbyt(fptr, (ntodo * twidth) as LONGLONG, buffer, status);
                } else {
                    ffgbytoff(fptr, twidth, ntodo, incre - twidth, buffer, status);
                }

                fffstri8(
                    cast_slice_mut(buffer),
                    ntodo,
                    scale,
                    zero,
                    twidth,
                    power,
                    nulcheck,
                    &snull,
                    nulval,
                    cast_slice_mut(&mut nularray[(next as usize)..]),
                    anynul.as_deref_mut(),
                    cast_slice_mut(&mut array[(next as usize)..]),
                    status,
                );
            }
            _ => {
                /*  error trap for invalid column format */
                int_snprintf!(
                    &mut message,
                    FLEN_ERRMSG,
                    "Cannot read numbers from column {} which has format {}",
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
                    "Error reading elements {:.0} thru {:.0} from column {} (ffgclj).",
                    dtemp + 1.0,
                    dtemp + ntodo as f64,
                    colnum,
                );
            } else {
                int_snprintf!(
                    &mut message,
                    FLEN_ERRMSG,
                    "Error reading elements {:.0} thru {:.0} from image (ffgclj).",
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
        if remain > 0 {
            next += ntodo as LONGLONG;
            elemnum += (ntodo * elemincre) as LONGLONG;

            if elemnum >= repeat {
                /* completed a row; start on later row */
                rowincre = elemnum / repeat;
                rownum += rowincre;
                elemnum -= rowincre * repeat;
            } else if elemnum < 0 {
                /* completed a row; start on a previous row */
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
/// Copy input to output following reading of the input from a FITS file.
/// Check for null values and do datatype conversion and scaling if required.
/// The nullcheck code value determines how any null values in the input array
/// are treated.  A null value is an input pixel that is equal to tnull.  If
/// nullcheck = 0, then no checking for nulls is performed and any null values
/// will be transformed just like any other pixel.  If nullcheck = 1, then the
/// output pixel will be set = nullval if the corresponding input pixel is null.
/// If nullcheck = 2, then if the pixel is null then the corresponding value of
/// nullarray will be set to 1; the value of nullarray for non-null pixels
/// will = 0.  The anynull parameter will be set = 1 if any of the returned
/// pixels are null, otherwise anynull will be returned with a value = 0;
pub(crate) fn fffi1i8(
    input: &[u8],             /* I - array of values to be converted     */
    ntodo: c_long,            /* I - number of elements in the array     */
    scale: f64,               /* I - FITS TSCALn or BSCALE value         */
    zero: f64,                /* I - FITS TZEROn or BZERO  value         */
    nullcheck: NullCheckType, /* I - null checking code; 0 = don't check */
    /*     1:set null pixels = nullval         */
    /*     2: if null pixel, set nullarray = 1 */
    tnull: u8,                   /* I - value of FITS TNULLn keyword if any */
    nullval: LONGLONG,           /* I - set null pixels, if nullcheck = 1   */
    nullarray: &mut [c_char],    /* I - bad pixel array, if nullcheck = 2   */
    anynull: Option<&mut c_int>, /* O - set to 1 if any pixels are null     */
    output: &mut [LONGLONG],     /* O - array of converted pixels           */
    status: &mut c_int,          /* IO - error status                       */
) -> c_int {
    let mut dvalue: f64 = 0.0;

    if nullcheck == NullCheckType::None {
        /* no null checking required */
        if scale == 1.0 && zero == 0.0 {
            /* no scaling */
            for ii in 0..(ntodo as usize) {
                output[ii] = input[ii] as LONGLONG; /* copy input to output */
            }
        } else {
            /* must scale the data */
            for ii in 0..(ntodo as usize) {
                dvalue = input[ii] as f64 * scale + zero;

                if dvalue < DLONGLONG_MIN {
                    *status = OVERFLOW_ERR;
                    output[ii] = LONGLONG_MIN;
                } else if dvalue > DLONGLONG_MAX {
                    *status = OVERFLOW_ERR;
                    output[ii] = LONGLONG_MAX;
                } else {
                    output[ii] = dvalue as LONGLONG;
                }
            }
        }
    } else {
        /* must check for null values */
        let anynull = anynull.unwrap();

        if scale == 1.0 && zero == 0.0 {
            /* no scaling */
            for ii in 0..(ntodo as usize) {
                if input[ii] == tnull {
                    *anynull = 1;
                    if nullcheck == NullCheckType::SetPixel {
                        output[ii] = nullval;
                    } else {
                        nullarray[ii] = 1;
                    }
                } else {
                    output[ii] = input[ii] as LONGLONG;
                }
            }
        } else {
            /* must scale the data */
            for ii in 0..(ntodo as usize) {
                if input[ii] == tnull {
                    *anynull = 1;
                    if nullcheck == NullCheckType::SetPixel {
                        output[ii] = nullval;
                    } else {
                        nullarray[ii] = 1;
                    }
                } else {
                    dvalue = input[ii] as f64 * scale + zero;

                    if dvalue < DLONGLONG_MIN {
                        *status = OVERFLOW_ERR;
                        output[ii] = LONGLONG_MIN;
                    } else if dvalue > DLONGLONG_MAX {
                        *status = OVERFLOW_ERR;
                        output[ii] = LONGLONG_MAX;
                    } else {
                        output[ii] = dvalue as LONGLONG;
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
/// will = 0.  The anynull parameter will be set = 1 if any of the returned
/// pixels are null, otherwise anynull will be returned with a value = 0;
pub(crate) fn fffi2i8(
    input: &[c_short],        /* I - array of values to be converted     */
    ntodo: c_long,            /* I - number of elements in the array     */
    scale: f64,               /* I - FITS TSCALn or BSCALE value         */
    zero: f64,                /* I - FITS TZEROn or BZERO  value         */
    nullcheck: NullCheckType, /* I - null checking code; 0 = don't check */
    /*     1:set null pixels = nullval         */
    /*     2: if null pixel, set nullarray = 1 */
    tnull: c_short,              /* I - value of FITS TNULLn keyword if any */
    nullval: LONGLONG,           /* I - set null pixels, if nullcheck = 1   */
    nullarray: &mut [c_char],    /* I - bad pixel array, if nullcheck = 2   */
    anynull: Option<&mut c_int>, /* O - set to 1 if any pixels are null     */
    output: &mut [LONGLONG],     /* O - array of converted pixels           */
    status: &mut c_int,          /* IO - error status                       */
) -> c_int {
    let mut dvalue: f64 = 0.0;

    if nullcheck == NullCheckType::None {
        /* no null checking required */
        if scale == 1.0 && zero == 0.0 {
            /* no scaling */
            for ii in 0..(ntodo as usize) {
                output[ii] = input[ii] as LONGLONG; /* copy input to output */
            }
        } else {
            /* must scale the data */
            for ii in 0..(ntodo as usize) {
                dvalue = input[ii] as f64 * scale + zero;

                if dvalue < DLONGLONG_MIN {
                    *status = OVERFLOW_ERR;
                    output[ii] = LONGLONG_MIN;
                } else if dvalue > DLONGLONG_MAX {
                    *status = OVERFLOW_ERR;
                    output[ii] = LONGLONG_MAX;
                } else {
                    output[ii] = dvalue as LONGLONG;
                }
            }
        }
    } else {
        /* must check for null values */
        let anynull = anynull.unwrap();

        if scale == 1.0 && zero == 0.0 {
            /* no scaling */
            for ii in 0..(ntodo as usize) {
                if input[ii] == tnull {
                    *anynull = 1;
                    if nullcheck == NullCheckType::SetPixel {
                        output[ii] = nullval;
                    } else {
                        nullarray[ii] = 1;
                    }
                } else {
                    output[ii] = input[ii] as LONGLONG;
                }
            }
        } else {
            /* must scale the data */
            for ii in 0..(ntodo as usize) {
                if input[ii] == tnull {
                    *anynull = 1;
                    if nullcheck == NullCheckType::SetPixel {
                        output[ii] = nullval;
                    } else {
                        nullarray[ii] = 1;
                    }
                } else {
                    dvalue = input[ii] as f64 * scale + zero;

                    if dvalue < DLONGLONG_MIN {
                        *status = OVERFLOW_ERR;
                        output[ii] = LONGLONG_MIN;
                    } else if dvalue > DLONGLONG_MAX {
                        *status = OVERFLOW_ERR;
                        output[ii] = LONGLONG_MAX;
                    } else {
                        output[ii] = dvalue as LONGLONG;
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
/// will = 0.  The anynull parameter will be set = 1 if any of the returned
/// pixels are null, otherwise anynull will be returned with a value = 0;
pub(crate) fn fffi4i8(
    input: &[INT32BIT],       /* I - array of values to be converted     */
    ntodo: c_long,            /* I - number of elements in the array     */
    scale: f64,               /* I - FITS TSCALn or BSCALE value         */
    zero: f64,                /* I - FITS TZEROn or BZERO  value         */
    nullcheck: NullCheckType, /* I - null checking code; 0 = don't check */
    /*     1:set null pixels = nullval         */
    /*     2: if null pixel, set nullarray = 1 */
    tnull: INT32BIT,             /* I - value of FITS TNULLn keyword if any */
    nullval: LONGLONG,           /* I - set null pixels, if nullcheck = 1   */
    nullarray: &mut [c_char],    /* I - bad pixel array, if nullcheck = 2   */
    anynull: Option<&mut c_int>, /* O - set to 1 if any pixels are null     */
    output: &mut [LONGLONG],     /* O - array of converted pixels           */
    status: &mut c_int,          /* IO - error status                       */
) -> c_int {
    let mut dvalue: f64 = 0.0;

    if nullcheck == NullCheckType::None {
        /* no null checking required */
        if scale == 1.0 && zero == 0.0 {
            /* no scaling */
            for ii in 0..(ntodo as usize) {
                output[ii] = input[ii] as LONGLONG; /* copy input to output */
            }
        } else {
            /* must scale the data */
            for ii in 0..(ntodo as usize) {
                dvalue = input[ii] as f64 * scale + zero;

                if dvalue < DLONGLONG_MIN {
                    *status = OVERFLOW_ERR;
                    output[ii] = LONGLONG_MIN;
                } else if dvalue > DLONGLONG_MAX {
                    *status = OVERFLOW_ERR;
                    output[ii] = LONGLONG_MAX;
                } else {
                    output[ii] = dvalue as LONGLONG;
                }
            }
        }
    } else {
        /* must check for null values */
        let anynull = anynull.unwrap();

        if scale == 1.0 && zero == 0.0 {
            /* no scaling */
            for ii in 0..(ntodo as usize) {
                if input[ii] == tnull {
                    *anynull = 1;
                    if nullcheck == NullCheckType::SetPixel {
                        output[ii] = nullval;
                    } else {
                        nullarray[ii] = 1;
                    }
                } else {
                    output[ii] = input[ii] as LONGLONG;
                }
            }
        } else {
            /* must scale the data */
            for ii in 0..(ntodo as usize) {
                if input[ii] == tnull {
                    *anynull = 1;
                    if nullcheck == NullCheckType::SetPixel {
                        output[ii] = nullval;
                    } else {
                        nullarray[ii] = 1;
                    }
                } else {
                    dvalue = input[ii] as f64 * scale + zero;

                    if dvalue < DLONGLONG_MIN {
                        *status = OVERFLOW_ERR;
                        output[ii] = LONGLONG_MIN;
                    } else if dvalue > DLONGLONG_MAX {
                        *status = OVERFLOW_ERR;
                        output[ii] = LONGLONG_MAX;
                    } else {
                        output[ii] = dvalue as LONGLONG;
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
/// will = 0.  The anynull parameter will be set = 1 if any of the returned
/// pixels are null, otherwise anynull will be returned with a value = 0;
pub(crate) fn fffi8i8(
    input: &[LONGLONG],       /* I - array of values to be converted     */
    ntodo: c_long,            /* I - number of elements in the array     */
    scale: f64,               /* I - FITS TSCALn or BSCALE value         */
    zero: f64,                /* I - FITS TZEROn or BZERO  value         */
    nullcheck: NullCheckType, /* I - null checking code; 0 = don't check */
    /*     1:set null pixels = nullval         */
    /*     2: if null pixel, set nullarray = 1 */
    tnull: LONGLONG,             /* I - value of FITS TNULLn keyword if any */
    nullval: LONGLONG,           /* I - set null pixels, if nullcheck = 1   */
    nullarray: &mut [c_char],    /* I - bad pixel array, if nullcheck = 2   */
    anynull: Option<&mut c_int>, /* O - set to 1 if any pixels are null     */
    output: &mut [LONGLONG],     /* O - array of converted pixels           */
    status: &mut c_int,          /* IO - error status                       */
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
                ulltemp = ((input[ii] as ULONGLONG) ^ 0x8000000000000000) as ULONGLONG;

                if ulltemp > LONGLONG_MAX as ULONGLONG {
                    *status = OVERFLOW_ERR;
                    output[ii] = LONGLONG_MAX;
                } else {
                    output[ii] = ulltemp as LONGLONG;
                }
            }
        } else if scale == 1.0 && zero == 0.0 {
            /* no scaling */
            output[..(ntodo as usize)].copy_from_slice(&input[..(ntodo as usize)]);
        /* copy input to output */
        } else {
            /* must scale the data */
            for ii in 0..(ntodo as usize) {
                dvalue = input[ii] as f64 * scale + zero;

                if dvalue < DLONGLONG_MIN {
                    *status = OVERFLOW_ERR;
                    output[ii] = LONGLONG_MIN;
                } else if dvalue > DLONGLONG_MAX {
                    *status = OVERFLOW_ERR;
                    output[ii] = LONGLONG_MAX;
                } else {
                    output[ii] = dvalue as LONGLONG;
                }
            }
        }
    } else {
        /* must check for null values */
        let anynull = anynull.unwrap();

        if scale == 1.0 && zero == 9223372036854775808.0 {
            /* The column we read contains unsigned long long values. */
            /* Instead of subtracting 9223372036854775808, it is more efficient */
            /* and more precise to just flip the sign bit with the XOR operator */

            for ii in 0..(ntodo as usize) {
                if input[ii] == tnull {
                    *anynull = 1;
                    if nullcheck == NullCheckType::SetPixel {
                        output[ii] = nullval;
                    } else {
                        nullarray[ii] = 1;
                    }
                } else {
                    ulltemp = ((input[ii] as ULONGLONG) ^ 0x8000000000000000) as ULONGLONG;

                    if ulltemp > LONGLONG_MAX as ULONGLONG {
                        *status = OVERFLOW_ERR;
                        output[ii] = LONGLONG_MAX;
                    } else {
                        output[ii] = ulltemp as LONGLONG;
                    }
                }
            }
        } else if scale == 1.0 && zero == 0.0 {
            /* no scaling */
            for ii in 0..(ntodo as usize) {
                if input[ii] == tnull {
                    *anynull = 1;
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
                    *anynull = 1;
                    if nullcheck == NullCheckType::SetPixel {
                        output[ii] = nullval;
                    } else {
                        nullarray[ii] = 1;
                    }
                } else {
                    dvalue = input[ii] as f64 * scale + zero;

                    if dvalue < DLONGLONG_MIN {
                        *status = OVERFLOW_ERR;
                        output[ii] = LONGLONG_MIN;
                    } else if dvalue > DLONGLONG_MAX {
                        *status = OVERFLOW_ERR;
                        output[ii] = LONGLONG_MAX;
                    } else {
                        output[ii] = dvalue as LONGLONG;
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
/// will = 0.  The anynull parameter will be set = 1 if any of the returned
/// pixels are null, otherwise anynull will be returned with a value = 0;
pub(crate) fn fffr4i8(
    input: &[f32],            /* I - array of values to be converted     */
    ntodo: c_long,            /* I - number of elements in the array     */
    scale: f64,               /* I - FITS TSCALn or BSCALE value         */
    zero: f64,                /* I - FITS TZEROn or BZERO  value         */
    nullcheck: NullCheckType, /* I - null checking code; 0 = don't check */
    /*     1:set null pixels = nullval         */
    /*     2: if null pixel, set nullarray = 1 */
    nullval: LONGLONG,           /* I - set null pixels, if nullcheck = 1   */
    nullarray: &mut [c_char],    /* I - bad pixel array, if nullcheck = 2   */
    anynull: Option<&mut c_int>, /* O - set to 1 if any pixels are null     */
    output: &mut [LONGLONG],     /* O - array of converted pixels           */
    status: &mut c_int,          /* IO - error status                       */
) -> c_int {
    let mut dvalue: f64 = 0.0;
    let mut sptr = 0;
    let iret = 0;

    if nullcheck == NullCheckType::None {
        /* no null checking required */
        if scale == 1.0 && zero == 0.0 {
            /* no scaling */
            for ii in 0..(ntodo as usize) {
                if (input[ii] as f64) < DLONGLONG_MIN {
                    *status = OVERFLOW_ERR;
                    output[ii] = LONGLONG_MIN;
                } else if (input[ii] as f64) > DLONGLONG_MAX {
                    *status = OVERFLOW_ERR;
                    output[ii] = LONGLONG_MAX;
                } else {
                    output[ii] = input[ii] as LONGLONG;
                }
            }
        } else {
            /* must scale the data */
            for ii in 0..(ntodo as usize) {
                dvalue = input[ii] as f64 * scale + zero;

                if dvalue < DLONGLONG_MIN {
                    *status = OVERFLOW_ERR;
                    output[ii] = LONGLONG_MIN;
                } else if dvalue > DLONGLONG_MAX {
                    *status = OVERFLOW_ERR;
                    output[ii] = LONGLONG_MAX;
                } else {
                    output[ii] = dvalue as LONGLONG;
                }
            }
        }
    } else {
        /* must check for null values */
        let anynull = anynull.unwrap();

        //sptr = (short *) input;

        if BYTESWAPPED && CFITSIO_MACHINE != VAXVMS && CFITSIO_MACHINE != ALPHAVMS {
            sptr += 1; /* point to MSBs */
        }

        let shortBuff: &[c_short] = cast_slice(input);

        if scale == 1.0 && zero == 0.0 {
            /* no scaling */
            for ii in 0..(ntodo as usize) {
                let iret = fnan(shortBuff[sptr]);
                if iret != 0 {
                    /* test for NaN or underflow */
                    if iret == 1 {
                        /* is it a NaN? */
                        *anynull = 1;
                        if nullcheck == NullCheckType::SetPixel {
                            output[ii] = nullval;
                        } else {
                            nullarray[ii] = 1;
                        }
                    } else {
                        /* it's an underflow */
                        output[ii] = 0;
                    }
                } else if (input[ii] as f64) < DLONGLONG_MIN {
                    *status = OVERFLOW_ERR;
                    output[ii] = LONGLONG_MIN;
                } else if (input[ii] as f64) > DLONGLONG_MAX {
                    *status = OVERFLOW_ERR;
                    output[ii] = LONGLONG_MAX;
                } else {
                    output[ii] = input[ii] as LONGLONG;
                }
                sptr += 2;
            }
        } else {
            /* must scale the data */
            for ii in 0..(ntodo as usize) {
                let iret = fnan(shortBuff[sptr]);
                if iret != 0 {
                    /* test for NaN or underflow */
                    if iret == 1 {
                        /* is it a NaN? */
                        *anynull = 1;
                        if nullcheck == NullCheckType::SetPixel {
                            output[ii] = nullval;
                        } else {
                            nullarray[ii] = 1;
                        }
                    } else {
                        /* it's an underflow */

                        if zero < DLONGLONG_MIN {
                            *status = OVERFLOW_ERR;
                            output[ii] = LONGLONG_MIN;
                        } else if zero > DLONGLONG_MAX {
                            *status = OVERFLOW_ERR;
                            output[ii] = LONGLONG_MAX;
                        } else {
                            output[ii] = zero as LONGLONG;
                        }
                    }
                } else {
                    dvalue = input[ii] as f64 * scale + zero;

                    if dvalue < DLONGLONG_MIN {
                        *status = OVERFLOW_ERR;
                        output[ii] = LONGLONG_MIN;
                    } else if dvalue > DLONGLONG_MAX {
                        *status = OVERFLOW_ERR;
                        output[ii] = LONGLONG_MAX;
                    } else {
                        output[ii] = dvalue as LONGLONG;
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
/// will = 0.  The anynull parameter will be set = 1 if any of the returned
/// pixels are null, otherwise anynull will be returned with a value = 0;
pub(crate) fn fffr8i8(
    input: &[f64],            /* I - array of values to be converted     */
    ntodo: c_long,            /* I - number of elements in the array     */
    scale: f64,               /* I - FITS TSCALn or BSCALE value         */
    zero: f64,                /* I - FITS TZEROn or BZERO  value         */
    nullcheck: NullCheckType, /* I - null checking code; 0 = don't check */
    /*     1:set null pixels = nullval         */
    /*     2: if null pixel, set nullarray = 1 */
    nullval: LONGLONG,           /* I - set null pixels, if nullcheck = 1   */
    nullarray: &mut [c_char],    /* I - bad pixel array, if nullcheck = 2   */
    anynull: Option<&mut c_int>, /* O - set to 1 if any pixels are null     */
    output: &mut [LONGLONG],     /* O - array of converted pixels           */
    status: &mut c_int,          /* IO - error status                       */
) -> c_int {
    let mut dvalue: f64 = 0.0;
    let mut sptr = 0;
    let iret = 0;

    if nullcheck == NullCheckType::None {
        /* no null checking required */
        if scale == 1.0 && zero == 0.0 {
            /* no scaling */
            for ii in 0..(ntodo as usize) {
                if input[ii] < DLONGLONG_MIN {
                    *status = OVERFLOW_ERR;
                    output[ii] = LONGLONG_MIN;
                } else if input[ii] > DLONGLONG_MAX {
                    *status = OVERFLOW_ERR;
                    output[ii] = LONGLONG_MAX;
                } else {
                    output[ii] = input[ii] as LONGLONG;
                }
            }
        } else {
            /* must scale the data */
            for ii in 0..(ntodo as usize) {
                dvalue = input[ii] * scale + zero;

                if dvalue < DLONGLONG_MIN {
                    *status = OVERFLOW_ERR;
                    output[ii] = LONGLONG_MIN;
                } else if dvalue > DLONGLONG_MAX {
                    *status = OVERFLOW_ERR;
                    output[ii] = LONGLONG_MAX;
                } else {
                    output[ii] = dvalue as LONGLONG;
                }
            }
        }
    } else {
        /* must check for null values */
        let anynull = anynull.unwrap();

        //sptr = (short *) input;

        if BYTESWAPPED && CFITSIO_MACHINE != VAXVMS && CFITSIO_MACHINE != ALPHAVMS {
            sptr += 3; /* point to MSBs */
        }

        let shortBuff: &[c_short] = cast_slice(input);

        if scale == 1.0 && zero == 0.0 {
            /* no scaling */
            for ii in 0..(ntodo as usize) {
                let iret = dnan(shortBuff[sptr]);
                if iret != 0 {
                    /* test for NaN or underflow */
                    if iret == 1 {
                        /* is it a NaN? */
                        *anynull = 1;
                        if nullcheck == NullCheckType::SetPixel {
                            output[ii] = nullval;
                        } else {
                            nullarray[ii] = 1;
                        }
                    } else {
                        /* it's an underflow */
                        output[ii] = 0;
                    }
                } else if input[ii] < DLONGLONG_MIN {
                    *status = OVERFLOW_ERR;
                    output[ii] = LONGLONG_MIN;
                } else if input[ii] > DLONGLONG_MAX {
                    *status = OVERFLOW_ERR;
                    output[ii] = LONGLONG_MAX;
                } else {
                    output[ii] = input[ii] as LONGLONG;
                }
                sptr += 4;
            }
        } else {
            /* must scale the data */
            for ii in 0..(ntodo as usize) {
                let iret = dnan(shortBuff[sptr]);
                if iret != 0 {
                    /* test for NaN or underflow */
                    if iret == 1 {
                        /* is it a NaN? */
                        *anynull = 1;
                        if nullcheck == NullCheckType::SetPixel {
                            output[ii] = nullval;
                        } else {
                            nullarray[ii] = 1;
                        }
                    } else {
                        /* it's an underflow */

                        if zero < DLONGLONG_MIN {
                            *status = OVERFLOW_ERR;
                            output[ii] = LONGLONG_MIN;
                        } else if zero > DLONGLONG_MAX {
                            *status = OVERFLOW_ERR;
                            output[ii] = LONGLONG_MAX;
                        } else {
                            output[ii] = zero as LONGLONG;
                        }
                    }
                } else {
                    dvalue = input[ii] * scale + zero;

                    if dvalue < DLONGLONG_MIN {
                        *status = OVERFLOW_ERR;
                        output[ii] = LONGLONG_MIN;
                    } else if dvalue > DLONGLONG_MAX {
                        *status = OVERFLOW_ERR;
                        output[ii] = LONGLONG_MAX;
                    } else {
                        output[ii] = dvalue as LONGLONG;
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
/// will = 0.  The anynull parameter will be set = 1 if any of the returned
/// pixels are null, otherwise anynull will be returned with a value = 0;
pub(crate) fn fffstri8(
    input: &mut [c_char],     /* I - array of values to be converted     */
    ntodo: c_long,            /* I - number of elements in the array     */
    scale: f64,               /* I - FITS TSCALn or BSCALE value         */
    zero: f64,                /* I - FITS TZEROn or BZERO  value         */
    twidth: c_long,           /* I - width of each substring of chars    */
    implipower: f64,          /* I - power of 10 of implied decimal      */
    nullcheck: NullCheckType, /* I - null checking code; 0 = don't check */
    /*     1:set null pixels = nullval         */
    /*     2: if null pixel, set nullarray = 1 */
    snull: &[c_char],                /* I - value of FITS null string, if any   */
    nullval: LONGLONG,               /* I - set null pixels, if nullcheck = 1   */
    nullarray: &mut [c_char],        /* I - bad pixel array, if nullcheck = 2   */
    mut anynull: Option<&mut c_int>, /* O - set to 1 if any pixels are null     */
    output: &mut [LONGLONG],         /* O - array of converted pixels           */
    status: &mut c_int,              /* IO - error status                       */
) -> c_int {
    let mut nullen: c_int = 0;
    let ii: c_long = 0;
    let mut dvalue: f64 = 0.0;

    let mut message: [c_char; FLEN_ERRMSG] = [0; FLEN_ERRMSG];
    let mut tempstore: c_char = 0;
    let chrzero: c_char = bb(b'0'); // 49
    let mut val: f64 = 0.0;
    let mut power: f64 = 0.0;
    let mut exponent: c_int = 0;
    let mut sign: c_int = 0;
    let mut esign: c_int = 0;
    let mut decpt: c_int = 0;

    nullen = strlen_safe(snull) as c_int;
    let mut cptr = 0; /* pointer to start of input string */
    let mut tpos = 0;
    for ii in 0..(ntodo as usize) {
        let originCptr = cptr;
        /* temporarily insert a null terminator at end of the string */
        tpos = cptr + twidth as usize;
        tempstore = input[tpos];
        input[tpos] = 0;

        /* check if null value is defined, and if the    */
        /* column string is identical to the null string */
        if snull[0] != ASCII_NULL_UNDEFINED
            && strncmp_safe(snull, &input[cptr..], nullen as usize) == 0
        {
            if nullcheck != NullCheckType::None {
                if let Some(anynull) = anynull.as_deref_mut() {
                    *anynull = 1;
                }
                if nullcheck == NullCheckType::SetPixel {
                    output[ii] = nullval;
                } else {
                    nullarray[ii] = 1;
                }
            }
            cptr += twidth as usize;
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
                val = val * 10.0 + input[cptr] as f64 - chrzero as f64; /* accumulate the value */
                cptr += 1;

                while input[cptr] == bb(b' ') {
                    /* skip embedded blanks in the value */
                    cptr += 1;
                }
            }

            if input[cptr] == bb(b'.') || input[cptr] == bb(b',') {
                /* check for decimal point */
                decpt = 1; /* set flag to show there was a decimal point */
                cptr += 1;
                while input[cptr] == bb(b' ') {
                    /* skip any blanks */
                    cptr += 1;
                }

                while input[cptr] >= bb(b'0') && input[cptr] <= bb(b'9') {
                    val = val * 10.0 + input[cptr] as f64 - chrzero as f64; /* accumulate the value */
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
                    exponent = exponent * 10 + input[cptr] as c_int - chrzero as c_int; /* accumulate exp */
                    cptr += 1;

                    while input[cptr] == bb(b' ') {
                        /* skip embedded blanks */
                        cptr += 1;
                    }
                }
            }

            if input[cptr] != 0 {
                /* should end up at the null terminator */
                let cstring = &input[originCptr..];

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
                    slice_to_str!(cstring),
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
            dvalue = (sign as f64 * val / power) * (10.0_f64).powf((esign * exponent) as f64);

            dvalue = dvalue * scale + zero; /* apply the scaling */

            if dvalue < DLONGLONG_MIN {
                *status = OVERFLOW_ERR;
                output[ii] = LONGLONG_MIN;
            } else if dvalue > DLONGLONG_MAX {
                *status = OVERFLOW_ERR;
                output[ii] = LONGLONG_MAX;
            } else {
                output[ii] = dvalue as LONGLONG;
            }
        }
        /* restore the char that was overwritten by the null */
        input[tpos] = tempstore;
    }
    *status
}
