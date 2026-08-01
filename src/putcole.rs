/*  This file, putcole.rs, contains routines that write data elements to    */
/*  a FITS image or table, with float datatype.                            */

/*  The FITSIO software was written by William Pence at the High Energy    */
/*  Astrophysic Science Archive Research Center (HEASARC) at the NASA      */
/*  Goddard Space Flight Center.                                           */

use core::ffi::CStr;
use core::slice;
use core::{cmp, mem};

use crate::imcompress::{fits_write_compressed_img, fits_write_compressed_pixels};
use crate::{NullCheckType, NullValue, c_types::*};

use bytemuck::{cast_slice, cast_slice_mut};

use crate::bb;
use crate::fitscore::{
    ffcfmt, ffgcprll, ffmahd_safe, ffpmsg_slice, ffpmsg_str, ffrdef_safe,
    fits_is_compressed_image_safe,
};
use crate::fitsio::*;
use crate::fitsio2::*;
use crate::putcolu::ffpcluc;
use crate::relibc::header::stdio::snprintf_f64;
use crate::wrappers::strlen_safe;
use crate::{buffers::*, int_snprintf, slice_to_str};

/*--------------------------------------------------------------------------*/
/// Write an array of values to the primary array. Data conversion
/// and scaling will be performed if necessary (e.g, if the datatype of
/// the FITS array is not the same as the array being written).
///
/// This routine cannot be called directly by users to write to large
/// arrays with > 2**31 pixels (although CFITSIO can do so by passing
/// the firstelem thru a LONGLONG sized global variable)
///
/// The primary array is represented as a binary table:
/// each group of the primary array is a row in the table,
/// where the first column contains the group parameters
/// and the second column contains the image itself.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffppre(
    fptr: *mut fitsfile, /* I - FITS file pointer                       */
    group: c_long,       /* I - group to write(1 = 1st group)           */
    firstelem: LONGLONG, /* I - first vector element to write(1 = 1st)  */
    nelem: LONGLONG,     /* I - number of values to write               */
    array: *const f32,   /* I - array of values that are written        */
    status: *mut c_int,  /* IO - error status                           */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let array = slice::from_raw_parts(array, nelem as usize);

        ffppre_safe(fptr, group, firstelem, nelem, array, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Write an array of values to the primary array. Data conversion
/// and scaling will be performed if necessary (e.g, if the datatype of
/// the FITS array is not the same as the array being written).
///
/// This routine cannot be called directly by users to write to large
/// arrays with > 2**31 pixels (although CFITSIO can do so by passing
/// the firstelem thru a LONGLONG sized global variable)
///
/// The primary array is represented as a binary table:
/// each group of the primary array is a row in the table,
/// where the first column contains the group parameters
/// and the second column contains the image itself.
pub fn ffppre_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer                       */
    group: c_long,       /* I - group to write(1 = 1st group)           */
    firstelem: LONGLONG, /* I - first vector element to write(1 = 1st)  */
    nelem: LONGLONG,     /* I - number of values to write               */
    array: &[f32],       /* I - array of values that are written        */
    status: &mut c_int,  /* IO - error status                           */
) -> c_int {
    let nullvalue: f32 = 0.0;

    if fits_is_compressed_image_safe(fptr, status) > 0 {
        /* this is a compressed image in a binary table */
        fits_write_compressed_pixels(
            fptr,
            TFLOAT,
            firstelem,
            nelem,
            NullCheckType::None,
            cast_slice(array),
            &Some(NullValue::Float(nullvalue)),
            status,
        );
        return *status;
    }

    let row = cmp::max(1, group);

    ffpcle_safe(
        fptr,
        2,
        row as LONGLONG,
        firstelem as LONGLONG,
        nelem as LONGLONG,
        array,
        status,
    );
    *status
}

/*--------------------------------------------------------------------------*/
/// Write an array of values to the primary array. Data conversion
/// and scaling will be performed if necessary (e.g, if the datatype of the
/// FITS array is not the same as the array being written).  Any array values
/// that are equal to the value of nulval will be replaced with the null
/// pixel value that is appropriate for this column.
///
/// This routine cannot be called directly by users to write to large
/// arrays with > 2**31 pixels (although CFITSIO can do so by passing
/// the firstelem thru a LONGLONG sized global variable)
///
/// The primary array is represented as a binary table:
/// each group of the primary array is a row in the table,
/// where the first column contains the group parameters
/// and the second column contains the image itself.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffppne(
    fptr: *mut fitsfile, /* I - FITS file pointer                       */
    group: c_long,       /* I - group to write(1 = 1st group)           */
    firstelem: LONGLONG, /* I - first vector element to write(1 = 1st)  */
    nelem: LONGLONG,     /* I - number of values to write               */
    array: *const f32,   /* I - array of values that are written        */
    nulval: f32,         /* I - undefined pixel value                   */
    status: *mut c_int,  /* IO - error status                           */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let array = slice::from_raw_parts(array, nelem as usize);

        ffppne_safe(fptr, group, firstelem, nelem, array, nulval, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Write an array of values to the primary array. Data conversion
/// and scaling will be performed if necessary (e.g, if the datatype of the
/// FITS array is not the same as the array being written).  Any array values
/// that are equal to the value of nulval will be replaced with the null
/// pixel value that is appropriate for this column.
///
/// This routine cannot be called directly by users to write to large
/// arrays with > 2**31 pixels (although CFITSIO can do so by passing
/// the firstelem thru a LONGLONG sized global variable)
///
/// The primary array is represented as a binary table:
/// each group of the primary array is a row in the table,
/// where the first column contains the group parameters
/// and the second column contains the image itself.
pub fn ffppne_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer                       */
    group: c_long,       /* I - group to write(1 = 1st group)           */
    firstelem: LONGLONG, /* I - first vector element to write(1 = 1st)  */
    nelem: LONGLONG,     /* I - number of values to write               */
    array: &[f32],       /* I - array of values that are written        */
    nulval: f32,         /* I - undefined pixel value                   */
    status: &mut c_int,  /* IO - error status                           */
) -> c_int {
    let nullvalue: f32;

    if fits_is_compressed_image_safe(fptr, status) > 0 {
        /* this is a compressed image in a binary table */

        nullvalue = nulval; /* set local variable */
        fits_write_compressed_pixels(
            fptr,
            TFLOAT,
            firstelem,
            nelem,
            NullCheckType::SetPixel,
            cast_slice(array),
            &Some(NullValue::Float(nullvalue)),
            status,
        );
        return *status;
    }

    let row = cmp::max(1, group);

    ffpcne_safe(
        fptr,
        2,
        row as LONGLONG,
        firstelem as LONGLONG,
        nelem as LONGLONG,
        array,
        nulval,
        status,
    );
    *status
}

/*--------------------------------------------------------------------------*/
/// Write an entire 2-D array of values to the primary array. Data conversion
/// and scaling will be performed if necessary (e.g, if the datatype of the
/// FITS array is not the same as the array being written).
///
/// This routine does not support writing to large images with
/// more than 2**31 pixels.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffp2de(
    fptr: *mut fitsfile, /* I - FITS file pointer                     */
    group: c_long,       /* I - group to write(1 = 1st group)         */
    ncols: LONGLONG,     /* I - number of pixels in each row of array */
    naxis1: LONGLONG,    /* I - FITS image NAXIS1 value               */
    naxis2: LONGLONG,    /* I - FITS image NAXIS2 value               */
    array: *const f32,   /* I - array to be written                   */
    status: *mut c_int,  /* IO - error status                         */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        let array = slice::from_raw_parts(array, (ncols * naxis2) as usize);

        ffp2de_safe(fptr, group, ncols, naxis1, naxis2, array, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Write an entire 2-D array of values to the primary array. Data conversion
/// and scaling will be performed if necessary (e.g, if the datatype of the
/// FITS array is not the same as the array being written).
///
/// This routine does not support writing to large images with
/// more than 2**31 pixels.
pub fn ffp2de_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer                     */
    group: c_long,       /* I - group to write(1 = 1st group)         */
    ncols: LONGLONG,     /* I - number of pixels in each row of array */
    naxis1: LONGLONG,    /* I - FITS image NAXIS1 value               */
    naxis2: LONGLONG,    /* I - FITS image NAXIS2 value               */
    array: &[f32],       /* I - array to be written                   */
    status: &mut c_int,  /* IO - error status                         */
) -> c_int {
    /* call the 3D writing routine, with the 3rd dimension = 1 */
    ffp3de_safe(fptr, group, ncols, naxis2, naxis1, naxis2, 1, array, status);

    *status
}

/*--------------------------------------------------------------------------*/
/// Write an entire 3-D cube of values to the primary array. Data conversion
/// and scaling will be performed if necessary (e.g, if the datatype of the
/// FITS array is not the same as the array being written).
///
/// This routine does not support writing to large images with
/// more than 2**31 pixels.
///
/// The primary array is represented as a binary table:
/// each group of the primary array is a row in the table,
/// where the first column contains the group parameters
/// and the second column contains the image itself.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffp3de(
    fptr: *mut fitsfile, /* I - FITS file pointer                     */
    group: c_long,       /* I - group to write(1 = 1st group)         */
    ncols: LONGLONG,     /* I - number of pixels in each row of array */
    nrows: LONGLONG,     /* I - number of rows in each plane of array */
    naxis1: LONGLONG,    /* I - FITS image NAXIS1 value               */
    naxis2: LONGLONG,    /* I - FITS image NAXIS2 value               */
    naxis3: LONGLONG,    /* I - FITS image NAXIS3 value               */
    array: *const f32,   /* I - array to be written                   */
    status: *mut c_int,  /* IO - error status                         */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        let array = slice::from_raw_parts(array, (ncols * naxis2 * naxis3) as usize);

        ffp3de_safe(
            fptr, group, ncols, nrows, naxis1, naxis2, naxis3, array, status,
        )
    }
}

/*--------------------------------------------------------------------------*/
/// Write an entire 3-D cube of values to the primary array. Data conversion
/// and scaling will be performed if necessary (e.g, if the datatype of the
/// FITS array is not the same as the array being written).
///
/// This routine does not support writing to large images with
/// more than 2**31 pixels.
///
/// The primary array is represented as a binary table:
/// each group of the primary array is a row in the table,
/// where the first column contains the group parameters
/// and the second column contains the image itself.
pub fn ffp3de_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer                     */
    group: c_long,       /* I - group to write(1 = 1st group)         */
    ncols: LONGLONG,     /* I - number of pixels in each row of array */
    nrows: LONGLONG,     /* I - number of rows in each plane of array */
    naxis1: LONGLONG,    /* I - FITS image NAXIS1 value               */
    naxis2: LONGLONG,    /* I - FITS image NAXIS2 value               */
    naxis3: LONGLONG,    /* I - FITS image NAXIS3 value               */
    array: &[f32],       /* I - array to be written                   */
    status: &mut c_int,  /* IO - error status                         */
) -> c_int {
    let fpixel: [c_long; 3] = [1; 3];
    let mut lpixel: [c_long; 3] = [0; 3];

    if fits_is_compressed_image_safe(fptr, status) > 0 {
        /* this is a compressed image in a binary table */
        lpixel[0] = ncols as c_long;
        lpixel[1] = nrows as c_long;
        lpixel[2] = naxis3 as c_long;

        fits_write_compressed_img(
            fptr,
            TFLOAT,
            &fpixel,
            &lpixel,
            NullCheckType::None,
            cast_slice(array),
            &None,
            status,
        );

        return *status;
    }

    let tablerow = cmp::max(1, group);

    if ncols == naxis1 && nrows == naxis2 {
        /* arrays have same size? */

        /* all the image pixels are contiguous, so write all at once */
        ffpcle_safe(
            fptr,
            2,
            tablerow as LONGLONG,
            1,
            naxis1 * naxis2 * naxis3,
            array,
            status,
        );
        return *status;
    }

    if ncols < naxis1 || nrows < naxis2 {
        *status = BAD_DIMEN;
        return *status;
    }

    let mut nfits = 1; /* next pixel in FITS image to write to */
    let mut narray = 0; /* next pixel in input array to be written */

    /* loop over naxis3 planes in the data cube */
    for _jj in 0..(naxis3 as usize) {
        /* loop over the naxis2 rows in the FITS image, */
        /* writing naxis1 pixels to each row            */

        for _ii in 0..(naxis2 as usize) {
            if ffpcle_safe(
                fptr,
                2,
                tablerow as LONGLONG,
                nfits,
                naxis1,
                &array[narray..],
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
/// Write a subsection of pixels to the primary array or image.
///
/// A subsection is defined to be any contiguous rectangular
/// array of pixels within the n-dimensional FITS data file.
/// Data conversion and scaling will be performed if necessary
/// (e.g, if the datatype of the FITS array is not the same as
/// the array being written).
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffpsse(
    fptr: *mut fitsfile,   /* I - FITS file pointer                       */
    group: c_long,         /* I - group to write(1 = 1st group)           */
    naxis: c_long,         /* I - number of data axes in array            */
    naxes: *const c_long,  /* I - size of each FITS axis                  */
    fpixel: *const c_long, /* I - 1st pixel in each axis to write (1=1st) */
    lpixel: *const c_long, /* I - last pixel in each axis to write        */
    array: *const f32,     /* I - array to be written                     */
    status: *mut c_int,    /* IO - error status                           */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        let fpixel = slice::from_raw_parts(fpixel, naxis as usize);
        let lpixel = slice::from_raw_parts(lpixel, naxis as usize);
        let naxes = slice::from_raw_parts(naxes, naxis as usize);

        let mut nelem = 1;
        for ii in 0..naxis as usize {
            nelem *= (lpixel[ii] - fpixel[ii] + 1) as usize;
        }

        let array = slice::from_raw_parts(array, nelem);

        ffpsse_safe(fptr, group, naxis, naxes, fpixel, lpixel, array, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Write a subsection of pixels to the primary array or image.
///
/// A subsection is defined to be any contiguous rectangular
/// array of pixels within the n-dimensional FITS data file.
/// Data conversion and scaling will be performed if necessary
/// (e.g, if the datatype of the FITS array is not the same as
/// the array being written).
pub fn ffpsse_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer                       */
    group: c_long,       /* I - group to write(1 = 1st group)           */
    naxis: c_long,       /* I - number of data axes in array            */
    naxes: &[c_long],    /* I - size of each FITS axis                  */
    fpixel: &[c_long],   /* I - 1st pixel in each axis to write (1=1st) */
    lpixel: &[c_long],   /* I - last pixel in each axis to write        */
    array: &[f32],       /* I - array to be written                     */
    status: &mut c_int,  /* IO - error status                           */
) -> c_int {
    let mut fpix: [LONGLONG; 7] = [0; 7];
    let mut dimen: [LONGLONG; 7] = [0; 7];
    let mut astart: LONGLONG;
    let mut pstart: LONGLONG;

    let mut st2: LONGLONG;
    let mut st3: LONGLONG;
    let mut st4: LONGLONG;
    let mut st5: LONGLONG;
    let mut st6: LONGLONG;
    let mut st7: LONGLONG;

    let mut irange: [c_long; 7] = [0; 7];

    if *status > 0 {
        return *status;
    }

    if fits_is_compressed_image_safe(fptr, status) > 0 {
        /* this is a compressed image in a binary table */
        fits_write_compressed_img(
            fptr,
            TFLOAT,
            fpixel,
            lpixel,
            NullCheckType::None,
            cast_slice(array),
            &None,
            status,
        );

        return *status;
    }

    if naxis < 1 || naxis > 7 {
        *status = BAD_DIMEN;
        return *status;
    }

    let tablerow = cmp::max(1, group);

    /* calculate the size and number of loops to perform in each dimension */
    for ii in 0..7 {
        fpix[ii] = 1;
        irange[ii] = 1;
        dimen[ii] = 1;
    }

    for ii in 0..(naxis as usize) {
        fpix[ii] = fpixel[ii] as LONGLONG;
        irange[ii] = lpixel[ii] - fpixel[ii] + 1;
        dimen[ii] = naxes[ii] as LONGLONG;
    }

    let i1 = irange[0];

    /* compute the pixel offset between each dimension */
    let off2: LONGLONG = dimen[0];
    let off3: LONGLONG = off2 * dimen[1];
    let off4: LONGLONG = off3 * dimen[2];
    let off5: LONGLONG = off4 * dimen[3];
    let off6: LONGLONG = off5 * dimen[4];
    let off7: LONGLONG = off6 * dimen[5];

    let st10: LONGLONG = fpix[0];
    let st20: LONGLONG = (fpix[1] - 1) * off2;
    let st30: LONGLONG = (fpix[2] - 1) * off3;
    let st40: LONGLONG = (fpix[3] - 1) * off4;
    let st50: LONGLONG = (fpix[4] - 1) * off5;
    let st60: LONGLONG = (fpix[5] - 1) * off6;
    let st70: LONGLONG = (fpix[6] - 1) * off7;

    /* store the initial offset in each dimension */
    let st1: LONGLONG = st10;
    st2 = st20;
    st3 = st30;
    st4 = st40;
    st5 = st50;
    st6 = st60;
    st7 = st70;

    astart = 0;

    for _i7 in 0..irange[6] {
        for _i6 in 0..irange[5] {
            for _i5 in 0..irange[4] {
                for _i4 in 0..irange[3] {
                    for _i3 in 0..irange[2] {
                        pstart = st1 + st2 + st3 + st4 + st5 + st6 + st7;

                        for _i2 in 0..irange[1] {
                            if ffpcle_safe(
                                fptr,
                                2,
                                tablerow as LONGLONG,
                                pstart,
                                i1 as LONGLONG,
                                &array[(astart as usize)..],
                                status,
                            ) > 0
                            {
                                return *status;
                            }

                            astart += i1 as LONGLONG;
                            pstart += off2;
                        }
                        st2 = st20;
                        st3 += off3;
                    }
                    st3 = st30;
                    st4 += off4;
                }
                st4 = st40;
                st5 += off5;
            }
            st5 = st50;
            st6 += off6;
        }
        st6 = st60;
        st7 += off7;
    }
    *status
}

/*--------------------------------------------------------------------------*/
/// Write an array of group parameters to the primary array. Data conversion
/// and scaling will be performed if necessary (e.g, if the datatype of
/// the FITS array is not the same as the array being written).
///
/// The primary array is represented as a binary table:
/// each group of the primary array is a row in the table,
/// where the first column contains the group parameters
/// and the second column contains the image itself.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffpgpe(
    fptr: *mut fitsfile, /* I - FITS file pointer                      */
    group: c_long,       /* I - group to write(1 = 1st group)          */
    firstelem: c_long,   /* I - first vector element to write(1 = 1st) */
    nelem: c_long,       /* I - number of values to write              */
    array: *const f32,   /* I - array of values that are written       */
    status: *mut c_int,  /* IO - error status                          */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        let array = slice::from_raw_parts(array, nelem as usize);

        ffpgpe_safe(fptr, group, firstelem, nelem, array, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Write an array of group parameters to the primary array. Data conversion
/// and scaling will be performed if necessary (e.g, if the datatype of
/// the FITS array is not the same as the array being written).
///
/// The primary array is represented as a binary table:
/// each group of the primary array is a row in the table,
/// where the first column contains the group parameters
/// and the second column contains the image itself.
pub fn ffpgpe_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer                      */
    group: c_long,       /* I - group to write(1 = 1st group)          */
    firstelem: c_long,   /* I - first vector element to write(1 = 1st) */
    nelem: c_long,       /* I - number of values to write              */
    array: &[f32],       /* I - array of values that are written       */
    status: &mut c_int,  /* IO - error status                          */
) -> c_int {
    let row = cmp::max(1, group);

    ffpcle_safe(
        fptr,
        1,
        row as LONGLONG,
        firstelem as LONGLONG,
        nelem as LONGLONG,
        array,
        status,
    );
    *status
}

/*--------------------------------------------------------------------------*/
/// Write an array of values to a column in the current FITS HDU.
///
/// The column number may refer to a real column in an ASCII or binary table,
/// or it may refer to a virtual column in a 1 or more grouped FITS primary
/// array.  FITSIO treats a primary array as a binary table
/// with 2 vector columns: the first column contains the group parameters (often
/// with length = 0) and the second column contains the array of image pixels.
/// Each row of the table represents a group in the case of multigroup FITS
/// images.
///
/// The input array of values will be converted to the datatype of the column
/// and will be inverse-scaled by the FITS TSCALn and TZEROn values if necessary.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffpcle(
    fptr: *mut fitsfile, /* I - FITS file pointer                       */
    colnum: c_int,       /* I - number of column to write (1 = 1st col) */
    firstrow: LONGLONG,  /* I - first row to write (1 = 1st row)        */
    firstelem: LONGLONG, /* I - first vector element to write (1 = 1st) */
    nelem: LONGLONG,     /* I - number of values to write               */
    array: *const f32,   /* I - array of values to write                */
    status: *mut c_int,  /* IO - error status                           */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        let array = slice::from_raw_parts(array, nelem as usize);

        ffpcle_safe(
            fptr,
            colnum,
            firstrow as LONGLONG,
            firstelem as LONGLONG,
            nelem as LONGLONG,
            array,
            status,
        )
    }
}

/*--------------------------------------------------------------------------*/
/// Write an array of values to a column in the current FITS HDU.
///
/// The column number may refer to a real column in an ASCII or binary table,
/// or it may refer to a virtual column in a 1 or more grouped FITS primary
/// array.  FITSIO treats a primary array as a binary table
/// with 2 vector columns: the first column contains the group parameters (often
/// with length = 0) and the second column contains the array of image pixels.
/// Each row of the table represents a group in the case of multigroup FITS
/// images.
///
/// The input array of values will be converted to the datatype of the column
/// and will be inverse-scaled by the FITS TSCALn and TZEROn values if necessary.
pub fn ffpcle_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer                       */
    colnum: c_int,       /* I - number of column to write (1 = 1st col) */
    firstrow: LONGLONG,  /* I - first row to write (1 = 1st row)        */
    firstelem: LONGLONG, /* I - first vector element to write (1 = 1st) */
    nelem: LONGLONG,     /* I - number of values to write               */
    array: &[f32],       /* I - array of values to write                */
    status: &mut c_int,  /* IO - error status                           */
) -> c_int {
    let mut tcode: c_int = 0;
    let mut maxelem2: c_int = 0;
    let mut hdutype: c_int = 0;
    let mut writeraw = false;
    let mut twidth: c_long = 0;
    let mut incre: c_long = 0;
    let mut ntodo: c_long = 0;
    let mut repeat: LONGLONG = 0;
    let mut startpos: LONGLONG = 0;
    let mut elemnum: LONGLONG = 0;
    let mut wrtptr: LONGLONG;
    let mut rowlen: LONGLONG = 0;
    let mut rownum: LONGLONG;
    let mut remain: LONGLONG;
    let mut next: LONGLONG;
    let mut tnull: LONGLONG = 0;
    let mut maxelem: LONGLONG;
    let mut scale: f64 = 0.0;
    let mut zero: f64 = 0.0;
    let mut tform: [c_char; 20] = [0; 20];
    let mut cform: [c_char; 20] = [0; 20];
    let mut message: [c_char; FLEN_ERRMSG] = [0; FLEN_ERRMSG];
    let mut snull: [c_char; 20] = [0; 20]; /*  the FITS null value  */
    let mut buffer: [f64; DBUFFSIZE as usize / mem::size_of::<f64>()] =
        [0.0; DBUFFSIZE as usize / mem::size_of::<f64>()]; /* align cbuff on word boundary */

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    /*---------------------------------------------------*/
    /*  Check input and get parameters about the column: */
    /*---------------------------------------------------*/
    if ffgcprll(
        fptr,
        colnum,
        firstrow,
        firstelem,
        nelem,
        1,
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
    maxelem = LONGLONG::from(maxelem2);

    if tcode == TSTRING {
        ffcfmt(&tform, &mut cform); /* derive C format for writing strings */
    }

    /*
       if there is no scaling and the native machine format is not byteswapped
       then we can simply write the raw data bytes into the FITS file if the
       datatype of the FITS column is the same as the input values.  Otherwise,
       we must convert the raw values into the scaled and/or machine dependent
       format in a temporary buffer that has been allocated for this purpose.
    */
    if scale == 1.0 && zero == 0.0 && CFITSIO_MACHINE == NATIVE && tcode == TFLOAT {
        writeraw = true;
        if nelem < LONGLONG::from(INT32_MAX) {
            maxelem = nelem;
        } else {
            maxelem = LONGLONG::from(INT32_MAX / 4);
        }
    } else {
        writeraw = false;
    }

    /*---------------------------------------------------------------------*/
    /*  Now write the pixels to the FITS column.                           */
    /*  First call the ffXXfYY routine to  (1) convert the datatype        */
    /*  if necessary, and (2) scale the values by the FITS TSCALn and      */
    /*  TZEROn linear scaling parameters into a temporary buffer.          */
    /*---------------------------------------------------------------------*/
    remain = nelem; /* remaining number of values to write  */
    next = 0; /* next element in array to be written  */
    rownum = 0; /* row number, relative to firstrow     */

    while remain > 0 {
        /* limit the number of pixels to process a one time to the number that
           will fit in the buffer space or to the number of pixels that remain
           in the current vector, which ever is smaller.
        */
        ntodo = cmp::min(remain, maxelem) as c_long;
        ntodo = cmp::min(ntodo as LONGLONG, repeat - elemnum) as c_long;

        wrtptr = startpos + ((rownum as LONGLONG) * rowlen) + (elemnum * incre as LONGLONG);

        ffmbyt_safe(fptr, wrtptr, IGNORE_EOF, status); /* move to write position */

        match tcode {
            TFLOAT => {
                if writeraw {
                    /* write raw input bytes without conversion */
                    ffpr4b(fptr, ntodo, incre, &array[(next as usize)..], status);
                } else {
                    /* convert the raw data before writing to FITS file */
                    ffr4fr4(
                        &array[(next as usize)..],
                        ntodo,
                        scale,
                        zero,
                        cast_slice_mut(&mut buffer),
                        status,
                    );
                    ffpr4b(fptr, ntodo, incre, cast_slice_mut(&mut buffer), status);
                }
            }
            TLONGLONG => {
                ffr4fi8(
                    &array[(next as usize)..],
                    ntodo,
                    scale,
                    zero,
                    cast_slice_mut(&mut buffer),
                    status,
                );
                ffpi8b(fptr, ntodo, incre, cast_slice_mut(&mut buffer), status);
            }
            TBYTE => {
                ffr4fi1(
                    &array[(next as usize)..],
                    ntodo,
                    scale,
                    zero,
                    cast_slice_mut(&mut buffer),
                    status,
                );
                ffpi1b(fptr, ntodo, incre, cast_slice_mut(&mut buffer), status);
            }
            TSHORT => {
                ffr4fi2(
                    &array[(next as usize)..],
                    ntodo,
                    scale,
                    zero,
                    cast_slice_mut(&mut buffer),
                    status,
                );
                ffpi2b(fptr, ntodo, incre, cast_slice_mut(&mut buffer), status);
            }
            TLONG => {
                ffr4fi4(
                    &array[(next as usize)..],
                    ntodo,
                    scale,
                    zero,
                    cast_slice_mut(&mut buffer),
                    status,
                );
                ffpi4b(fptr, ntodo, incre, cast_slice_mut(&mut buffer), status);
            }
            TDOUBLE => {
                ffr4fr8(
                    &array[(next as usize)..],
                    ntodo,
                    scale,
                    zero,
                    cast_slice_mut(&mut buffer),
                    status,
                );
                ffpr8b(fptr, ntodo, incre, cast_slice_mut(&mut buffer), status);
            }
            TSTRING => {
                /* numerical column in an ASCII table */

                let formlen = strlen_safe(&cform);

                if hdutype == ASCII_TBL && formlen > 1 {
                    if cform[formlen - 1] == bb(b'f') || cform[formlen - 1] == bb(b'E') {
                        ffr4fstr(
                            &array[(next as usize)..],
                            ntodo,
                            scale,
                            zero,
                            &cform,
                            twidth,
                            cast_slice_mut(&mut buffer),
                            status,
                        );

                        if incre == twidth {
                            /* contiguous bytes */
                            ffpbyt(
                                fptr,
                                (ntodo * twidth) as LONGLONG,
                                cast_slice_mut(&mut buffer),
                                status,
                            );
                        } else {
                            ffpbytoff(
                                fptr,
                                twidth,
                                ntodo,
                                incre - twidth,
                                cast_slice_mut(&mut buffer),
                                status,
                            );
                        }
                    }
                } else {
                    /* can't write to string column, so fall thru to default: */
                    /*  error trap  */
                    int_snprintf!(
                        &mut message,
                        FLEN_ERRMSG,
                        "Cannot write numbers to column {} which has format {}",
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
            }
            _ => {
                /*  error trap  */
                int_snprintf!(
                    &mut message,
                    FLEN_ERRMSG,
                    "Cannot write numbers to column {} which has format {}",
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
            /* test for error during previous write operation */

            int_snprintf!(
                &mut message,
                FLEN_ERRMSG,
                "Error writing elements {:.0} thru {:.0} of input data array (ffpcle).",
                (next + 1) as f64,
                (next + ntodo as LONGLONG) as f64,
            );
            ffpmsg_slice(&message);
            return *status;
        }

        /*--------------------------------------------*/
        /*  increment the counters for the next loop  */
        /*--------------------------------------------*/
        remain -= ntodo as LONGLONG;
        if remain > 0 {
            next += ntodo as LONGLONG;
            elemnum += ntodo as LONGLONG;
            if elemnum == repeat {
                /* completed a row; start on next row */

                elemnum = 0;
                rownum += 1;
            }
        }
    } /*  End of main while Loop  */

    /*--------------------------------*/
    /*  check for numerical overflow  */
    /*--------------------------------*/
    if *status == OVERFLOW_ERR {
        ffpmsg_str("Numerical overflow during type conversion while writing FITS data.");
        *status = NUM_OVERFLOW;
    }

    *status
}

/*--------------------------------------------------------------------------*/
/// Write an array of complex values to a column in the current FITS HDU.
///
/// Each complex number if interpreted as a pair of float values.
/// The column number may refer to a real column in an ASCII or binary table,
/// or it may refer to a virtual column in a 1 or more grouped FITS primary
/// array.  FITSIO treats a primary array as a binary table
/// with 2 vector columns: the first column contains the group parameters (often
/// with length = 0) and the second column contains the array of image pixels.
/// Each row of the table represents a group in the case of multigroup FITS
/// images.
///
/// The input array of values will be converted to the datatype of the column
/// if necessary, but normally complex values should only be written to a binary
/// table with TFORMn = 'rC' where r is an optional repeat count. The TSCALn and
/// TZERO keywords should not be used with complex numbers because mathmatically
/// the scaling should only be applied to the real (first) component of the
/// complex value.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffpclc(
    fptr: *mut fitsfile, /* I - FITS file pointer                       */
    colnum: c_int,       /* I - number of column to write (1 = 1st col) */
    firstrow: LONGLONG,  /* I - first row to write (1 = 1st row)        */
    firstelem: LONGLONG, /* I - first vector element to write (1 = 1st) */
    nelem: LONGLONG,     /* I - number of values to write               */
    array: *const f32,   /* I - array of values to write                */
    status: *mut c_int,  /* IO - error status                           */
) -> c_int {
    unsafe {
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);

        let array = slice::from_raw_parts(array, 2 * nelem as usize);

        ffpclc_safe(
            fptr,
            colnum,
            firstrow as LONGLONG,
            firstelem as LONGLONG,
            nelem as LONGLONG,
            array,
            status,
        )
    }
}

/*--------------------------------------------------------------------------*/
/// Write an array of complex values to a column in the current FITS HDU.
/// Each complex number if interpreted as a pair of float values.
/// The column number may refer to a real column in an ASCII or binary table,
/// or it may refer to a virtual column in a 1 or more grouped FITS primary
/// array.  FITSIO treats a primary array as a binary table
/// with 2 vector columns: the first column contains the group parameters (often
/// with length = 0) and the second column contains the array of image pixels.
/// Each row of the table represents a group in the case of multigroup FITS
/// images.
///
/// The input array of values will be converted to the datatype of the column
/// if necessary, but normally complex values should only be written to a binary
/// table with TFORMn = 'rC' where r is an optional repeat count. The TSCALn and
/// TZERO keywords should not be used with complex numbers because mathmatically
/// the scaling should only be applied to the real (first) component of the
/// complex value.
pub fn ffpclc_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer                       */
    colnum: c_int,       /* I - number of column to write (1 = 1st col) */
    firstrow: LONGLONG,  /* I - first row to write (1 = 1st row)        */
    firstelem: LONGLONG, /* I - first vector element to write (1 = 1st) */
    nelem: LONGLONG,     /* I - number of values to write               */
    array: &[f32],       /* I - array of values to write                */
    status: &mut c_int,  /* IO - error status                           */
) -> c_int {
    /* simply multiply the number of elements by 2, and call ffpcle */
    ffpcle_safe(
        fptr,
        colnum,
        firstrow,
        (firstelem - 1) * 2 + 1,
        nelem * 2,
        array,
        status,
    );
    *status
}

/*--------------------------------------------------------------------------*/
/// Write an array of elements to the specified column of a table.  Any input
/// pixels equal to the value of nulvalue will be replaced by the appropriate
/// null value in the output FITS file.
///
/// The input array of values will be converted to the datatype of the column
/// and will be inverse-scaled by the FITS TSCALn and TZEROn values if necessary
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffpcne(
    fptr: *mut fitsfile, /* I - FITS file pointer                       */
    colnum: c_int,       /* I - number of column to write (1 = 1st col) */
    firstrow: LONGLONG,  /* I - first row to write (1 = 1st row)        */
    firstelem: LONGLONG, /* I - first vector element to write (1 = 1st) */
    nelem: LONGLONG,     /* I - number of values to write               */
    array: *const f32,   /* I - array of values to write                */
    nulvalue: f32,       /* I - value used to flag undefined pixels     */
    status: *mut c_int,  /* IO - error status                           */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        let array = slice::from_raw_parts(array, nelem as usize);

        ffpcne_safe(
            fptr,
            colnum,
            firstrow as LONGLONG,
            firstelem as LONGLONG,
            nelem as LONGLONG,
            array,
            nulvalue,
            status,
        )
    }
}

/*--------------------------------------------------------------------------*/
/// Write an array of elements to the specified column of a table.  Any input
/// pixels equal to the value of nulvalue will be replaced by the appropriate
/// null value in the output FITS file.
///
/// The input array of values will be converted to the datatype of the column
/// and will be inverse-scaled by the FITS TSCALn and TZEROn values if necessary
pub fn ffpcne_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer                       */
    colnum: c_int,       /* I - number of column to write (1 = 1st col) */
    firstrow: LONGLONG,  /* I - first row to write (1 = 1st row)        */
    firstelem: LONGLONG, /* I - first vector element to write (1 = 1st) */
    nelem: LONGLONG,     /* I - number of values to write               */
    array: &[f32],       /* I - array of values to write                */
    nulvalue: f32,       /* I - value used to flag undefined pixels     */
    status: &mut c_int,  /* IO - error status                           */
) -> c_int {
    let mut ngood: LONGLONG = 0;
    let mut nbad: LONGLONG = 0;
    let mut repeat: LONGLONG = 0;

    let mut fstelm: LONGLONG;
    let mut fstrow: LONGLONG;

    let mut tcode = 0;
    let mut overflow = 0;

    if *status > 0 {
        return *status;
    }

    /* reset position to the correct HDU if necessary */
    if fptr.HDUposition != fptr.Fptr.curhdu {
        ffmahd_safe(fptr, (fptr.HDUposition) + 1, None, status);
    } else if fptr.Fptr.datastart == DATA_UNDEFINED as LONGLONG && ffrdef_safe(fptr, status) > 0 {
        /* rescan header */
        return *status;
    }

    /* set pointer to first column */
    let c = fptr.Fptr.get_tableptr_as_slice();
    let ci = colnum as usize - 1; /* offset to correct column structure */

    tcode = c[ci].tdatatype;

    if tcode > 0 {
        repeat = c[ci].trepeat; /* repeat count for this column */
    } else {
        repeat = firstelem - 1 + nelem; /* variable length arrays */
    }

    if tcode.abs() >= TCOMPLEX {
        /* treat complex columns as pairs of numbers */
        repeat *= 2;
    }

    /* if variable length array, first write the whole input vector,
    then go back and fill in the nulls */
    if tcode < 0
        && ffpcle_safe(
            fptr,
            colnum,
            firstrow as LONGLONG,
            firstelem as LONGLONG,
            nelem as LONGLONG,
            array,
            status,
        ) > 0
    {
        if *status == NUM_OVERFLOW {
            /* ignore overflows, which are possibly the null pixel values */
            /*  overflow = 1;   */
            *status = 0;
        } else {
            return *status;
        }
    }

    /* absolute element number in the column */
    let first: LONGLONG = (firstrow - 1) * repeat + firstelem;

    let mut ii: usize = 0;
    while ii < nelem as usize {
        if array[ii] != nulvalue {
            /* is this a good pixel? */

            if nbad > 0 {
                /* write previous string of bad pixels */

                fstelm = (ii as LONGLONG) - nbad + first; /* absolute element number */
                fstrow = (fstelm - 1) / repeat + 1; /* starting row number */
                fstelm -= (fstrow - 1) * repeat; /* relative number */

                /* call ffpcluc, not ffpclu, in case we are writing to a complex ('C') binary table column */
                if ffpcluc(fptr, colnum, fstrow, fstelm, nbad, status) > 0 {
                    return *status;
                }

                nbad = 0;
            }

            ngood += 1; /* the consecutive number of good pixels */
        } else {
            if ngood > 0 {
                /* write previous string of good pixels */

                fstelm = (ii as LONGLONG) - ngood + first; /* absolute element number */
                fstrow = (fstelm - 1) / repeat + 1; /* starting row number */
                fstelm -= (fstrow - 1) * repeat; /* relative number */

                if tcode > 0 {
                    /* variable length arrays have already been written */
                    if ffpcle_safe(
                        fptr,
                        colnum,
                        fstrow,
                        fstelm,
                        ngood,
                        &array[(ii - ngood as usize)..],
                        status,
                    ) > 0
                    {
                        if *status == NUM_OVERFLOW {
                            overflow = 1;
                            *status = 0;
                        } else {
                            return *status;
                        }
                    }
                }
                ngood = 0;
            }

            nbad += 1; /* the consecutive number of bad pixels */
        }
        ii += 1;
    }

    /* finished loop;  now just write the last set of pixels */

    if ngood > 0 {
        /* write last string of good pixels */
        fstelm = (ii as LONGLONG) - ngood + first; /* absolute element number */
        fstrow = (fstelm - 1) / repeat + 1; /* starting row number */
        fstelm -= (fstrow - 1) * repeat; /* relative number */

        if tcode > 0 {
            /* variable length arrays have already been written */
            ffpcle_safe(
                fptr,
                colnum,
                fstrow,
                fstelm,
                ngood,
                &array[(ii - ngood as usize)..],
                status,
            );
        }
    } else if nbad > 0 {
        /* write last string of bad pixels */

        fstelm = (ii as LONGLONG) - nbad + first; /* absolute element number */
        fstrow = (fstelm - 1) / repeat + 1; /* starting row number */
        fstelm -= (fstrow - 1) * repeat; /* relative number */
        ffpcluc(fptr, colnum, fstrow, fstelm, nbad, status);
    }

    if *status <= 0 && overflow > 0 {
        *status = NUM_OVERFLOW;
    }

    *status
}

/*--------------------------------------------------------------------------*/
/// Copy input to output prior to writing output to a FITS file.
/// Do datatype conversion and scaling if required.
pub(crate) fn ffr4fi1(
    input: &[f32],      /* I - array of values to be converted  */
    ntodo: c_long,      /* I - number of elements in the array  */
    scale: f64,         /* I - FITS TSCALn or BSCALE value      */
    zero: f64,          /* I - FITS TZEROn or BZERO  value      */
    output: &mut [u8],  /* O - output array of converted values */
    status: &mut c_int, /* IO - error status                    */
) -> c_int {
    if scale == 1.0 && zero == 0.0 {
        for ii in 0..(ntodo as usize) {
            if f64::from(input[ii]) < DUCHAR_MIN {
                *status = OVERFLOW_ERR;
                output[ii] = 0;
            } else if f64::from(input[ii]) > DUCHAR_MAX {
                *status = OVERFLOW_ERR;
                output[ii] = u8::MAX;
            } else {
                output[ii] = input[ii] as u8;
            }
        }
    } else {
        for ii in 0..(ntodo as usize) {
            let dvalue: f64 = (f64::from(input[ii]) - zero) / scale;

            if dvalue < DUCHAR_MIN {
                *status = OVERFLOW_ERR;
                output[ii] = 0;
            } else if dvalue > DUCHAR_MAX {
                *status = OVERFLOW_ERR;
                output[ii] = u8::MAX;
            } else {
                output[ii] = (dvalue + 0.5) as u8;
            }
        }
    }
    *status
}

/*--------------------------------------------------------------------------*/
/// Copy input to output prior to writing output to a FITS file.
/// Do datatype conversion and scaling if required.
pub(crate) fn ffr4fi2(
    input: &[f32],          /* I - array of values to be converted  */
    ntodo: c_long,          /* I - number of elements in the array  */
    scale: f64,             /* I - FITS TSCALn or BSCALE value      */
    zero: f64,              /* I - FITS TZEROn or BZERO  value      */
    output: &mut [c_short], /* O - output array of converted values */
    status: &mut c_int,     /* IO - error status                    */
) -> c_int {
    if scale == 1.0 && zero == 0.0 {
        for ii in 0..(ntodo as usize) {
            if f64::from(input[ii]) < DSHRT_MIN {
                *status = OVERFLOW_ERR;
                output[ii] = c_short::MIN;
            } else if f64::from(input[ii]) > DSHRT_MAX {
                *status = OVERFLOW_ERR;
                output[ii] = c_short::MAX;
            } else {
                output[ii] = input[ii] as c_short;
            }
        }
    } else {
        for ii in 0..(ntodo as usize) {
            let dvalue: f64 = (f64::from(input[ii]) - zero) / scale;

            if dvalue < DSHRT_MIN {
                *status = OVERFLOW_ERR;
                output[ii] = c_short::MIN;
            } else if dvalue > DSHRT_MAX {
                *status = OVERFLOW_ERR;
                output[ii] = c_short::MAX;
            } else if dvalue >= 0.0 {
                output[ii] = (dvalue + 0.5) as c_short;
            } else {
                output[ii] = (dvalue - 0.5) as c_short;
            }
        }
    }
    *status
}

/*--------------------------------------------------------------------------*/
/// Copy input to output prior to writing output to a FITS file.
/// Do datatype conversion and scaling if required.
pub(crate) fn ffr4fi4(
    input: &[f32],           /* I - array of values to be converted  */
    ntodo: c_long,           /* I - number of elements in the array  */
    scale: f64,              /* I - FITS TSCALn or BSCALE value      */
    zero: f64,               /* I - FITS TZEROn or BZERO  value      */
    output: &mut [INT32BIT], /* O - output array of converted values */
    status: &mut c_int,      /* IO - error status                    */
) -> c_int {
    if scale == 1.0 && zero == 0.0 {
        for ii in 0..(ntodo as usize) {
            if f64::from(input[ii]) < DINT_MIN {
                *status = OVERFLOW_ERR;
                output[ii] = INT32_MIN;
            } else if f64::from(input[ii]) > DINT_MAX {
                *status = OVERFLOW_ERR;
                output[ii] = INT32_MAX;
            } else {
                output[ii] = input[ii] as INT32BIT;
            }
        }
    } else {
        for ii in 0..(ntodo as usize) {
            let dvalue: f64 = (f64::from(input[ii]) - zero) / scale;

            if dvalue < DINT_MIN {
                *status = OVERFLOW_ERR;
                output[ii] = INT32_MIN;
            } else if dvalue > DINT_MAX {
                *status = OVERFLOW_ERR;
                output[ii] = INT32_MAX;
            } else if dvalue >= 0.0 {
                output[ii] = (dvalue + 0.5) as INT32BIT;
            } else {
                output[ii] = (dvalue - 0.5) as INT32BIT;
            }
        }
    }
    *status
}

/*--------------------------------------------------------------------------*/
/// Copy input to output prior to writing output to a FITS file.
/// Do datatype conversion and scaling if required.
pub(crate) fn ffr4fi8(
    input: &[f32],           /* I - array of values to be converted  */
    ntodo: c_long,           /* I - number of elements in the array  */
    scale: f64,              /* I - FITS TSCALn or BSCALE value      */
    zero: f64,               /* I - FITS TZEROn or BZERO  value      */
    output: &mut [LONGLONG], /* O - output array of converted values */
    status: &mut c_int,      /* IO - error status                    */
) -> c_int {
    if scale == 1.0 && zero == 9223372036854775808. {
        /* Writing to unsigned long long column. Input values must not be negative */
        /* Instead of subtracting 9223372036854775808, it is more efficient */
        /* and more precise to just flip the sign bit with the XOR operator */

        for ii in 0..(ntodo as usize) {
            if input[ii] < -0.49 {
                *status = OVERFLOW_ERR;
                output[ii] = LONGLONG_MIN;
            } else if f64::from(input[ii]) > 2.0 * DLONGLONG_MAX {
                *status = OVERFLOW_ERR;
                output[ii] = LONGLONG_MAX;
            } else {
                output[ii] = ((input[ii] as ULONGLONG) ^ 0x8000000000000000) as LONGLONG;
            }
        }
    } else if scale == 1.0 && zero == 0.0 {
        for ii in 0..(ntodo as usize) {
            if f64::from(input[ii]) < DLONGLONG_MIN {
                *status = OVERFLOW_ERR;
                output[ii] = LONGLONG_MIN;
            } else if f64::from(input[ii]) > DLONGLONG_MAX {
                *status = OVERFLOW_ERR;
                output[ii] = LONGLONG_MAX;
            } else {
                output[ii] = input[ii] as LONGLONG;
            }
        }
    } else {
        for ii in 0..(ntodo as usize) {
            let dvalue: f64 = (f64::from(input[ii]) - zero) / scale;

            if dvalue < DLONGLONG_MIN {
                *status = OVERFLOW_ERR;
                output[ii] = LONGLONG_MIN;
            } else if dvalue > DLONGLONG_MAX {
                *status = OVERFLOW_ERR;
                output[ii] = LONGLONG_MAX;
            } else if dvalue >= 0.0 {
                output[ii] = (dvalue + 0.5) as LONGLONG;
            } else {
                output[ii] = (dvalue - 0.5) as LONGLONG;
            }
        }
    }
    *status
}

/*--------------------------------------------------------------------------*/
/// Copy input to output prior to writing output to a FITS file.
/// Do datatype conversion and scaling if required.
pub(crate) fn ffr4fr4(
    input: &[f32],      /* I - array of values to be converted  */
    ntodo: c_long,      /* I - number of elements in the array  */
    scale: f64,         /* I - FITS TSCALn or BSCALE value      */
    zero: f64,          /* I - FITS TZEROn or BZERO  value      */
    output: &mut [f32], /* O - output array of converted values */
    status: &mut c_int, /* IO - error status                    */
) -> c_int {
    let ntodo = ntodo as usize;
    if scale == 1.0 && zero == 0.0 {
        output[..ntodo].copy_from_slice(&input[..ntodo]); /* copy input to output */
    } else {
        for ii in 0..(ntodo) {
            output[ii] = ((f64::from(input[ii]) - zero) / scale) as f32;
        }
    }
    *status
}

/*--------------------------------------------------------------------------*/
/// Copy input to output prior to writing output to a FITS file.
/// Do datatype conversion and scaling if required.
pub(crate) fn ffr4fr8(
    input: &[f32],      /* I - array of values to be converted  */
    ntodo: c_long,      /* I - number of elements in the array  */
    scale: f64,         /* I - FITS TSCALn or BSCALE value      */
    zero: f64,          /* I - FITS TZEROn or BZERO  value      */
    output: &mut [f64], /* O - output array of converted values */
    status: &mut c_int, /* IO - error status                    */
) -> c_int {
    if scale == 1.0 && zero == 0.0 {
        for ii in 0..(ntodo as usize) {
            output[ii] = f64::from(input[ii]);
        }
    } else {
        for ii in 0..(ntodo as usize) {
            output[ii] = (f64::from(input[ii]) - zero) / scale;
        }
    }
    *status
}

/*--------------------------------------------------------------------------*/
/// Copy input to output prior to writing output to a FITS file.
/// Do scaling if required.
pub(crate) fn ffr4fstr(
    input: &[f32],         /* I - array of values to be converted  */
    ntodo: c_long,         /* I - number of elements in the array  */
    scale: f64,            /* I - FITS TSCALn or BSCALE value      */
    zero: f64,             /* I - FITS TZEROn or BZERO  value      */
    cform: &[c_char],      /* I - format for output string values  */
    twidth: c_long,        /* I - width of each field, in chars    */
    output: &mut [c_char], /* O - output array of converted values */
    status: &mut c_int,    /* IO - error status                    */
) -> c_int {
    let mut ci = 0;

    if scale == 1.0 && zero == 0.0 {
        for ii in 0..(ntodo as usize) {
            snprintf_f64(
                &mut output[ci..],
                DBUFFSIZE as usize,
                cform,
                f64::from(input[ii]),
            );
            ci += twidth as usize;

            if output[ci] != 0 {
                /* if this char != \0, then overflow occurred */
                *status = OVERFLOW_ERR;
            }
        }
    } else {
        for ii in 0..(ntodo as usize) {
            let dvalue: f64 = (f64::from(input[ii]) - zero) / scale;
            snprintf_f64(&mut output[ci..], DBUFFSIZE as usize, cform, dvalue);
            ci += twidth as usize;

            if output[ci] != 0 {
                /* if this char != \0, then overflow occurred */
                *status = OVERFLOW_ERR;
            }
        }
    }

    /* replace any commas with periods (e.g., in French locale) */
    for i in output.iter_mut() {
        if *i == bb(b',') {
            *i = bb(b'.');
        }
    }

    *status
}

// Ported from test_putcole.c - float (f32) precision write functions.

#[cfg(test)]
mod tests {
    use crate::aliases::rust_api::*;
    use crate::fitsio::{
        ASCII_TBL, BAD_BTABLE_FORMAT, BAD_DIMEN, BINARY_TBL, BYTE_IMG, FLOAT_IMG, LONGLONG,
        NUM_OVERFLOW, READONLY, TRUE, fitsfile,
    };
    use crate::helpers::testhelpers::{to_buf, with_temp_file};
    use libc::{c_char, c_int, c_long, c_uint};

    /// Make a NUL-terminated `Vec<c_char>` from a `&str`.
    fn cc(s: &str) -> Vec<c_char> {
        let mut v: Vec<c_char> = s.bytes().map(|b| b as c_char).collect();
        v.push(0);
        v
    }

    /// Create a single-column table (binary or ASCII) with an empty BYTE_IMG
    /// primary HDU, returning the open file positioned on the table HDU.
    fn make_table(
        name: &[c_char],
        tbltype: c_int,
        ttype: &str,
        tform: &str,
        nrows: LONGLONG,
        status: &mut c_int,
    ) -> Option<Box<fitsfile>> {
        let mut f: Option<Box<fitsfile>> = None;
        fits_create_file(&mut f, name, status);
        fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], status);
        let ttype_v = [Some(cc(ttype))];
        let ttype_ref: Vec<Option<&[c_char]>> = ttype_v.iter().map(|o| o.as_deref()).collect();
        let tform_v = [cc(tform)];
        let tform_ref: Vec<&[c_char]> = tform_v.iter().map(|v| v.as_slice()).collect();
        fits_create_tbl(
            f.as_deref_mut().unwrap(),
            tbltype,
            nrows,
            1,
            &ttype_ref,
            &tform_ref,
            None,
            None,
            status,
        );
        f
    }

    #[test]
    fn test_write_primary_array() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 1] = [5];
            let data: [f32; 5] = [-1.5e30, -1.0, 0.0, 1.0, 1.5e30];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), FLOAT_IMG, 1, &naxes, &mut status);
            fits_write_img_flt(f.as_deref_mut().unwrap(), 1, 1, 5, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0, "setup");

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut result = [0f32; 5];
            let mut anynull = -1;
            fits_read_img_flt(
                f.as_deref_mut().unwrap(),
                1,
                1,
                5,
                0.0,
                &mut result,
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[0], data[0]);
            assert_eq!(result[2], 0.0);
            assert_eq!(result[4], data[4]);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_primary_with_null() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 1] = [5];
            let data: [f32; 5] = [1.0, 2.0, -999.0, 4.0, 5.0];
            let nulval: f32 = -999.0;

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), FLOAT_IMG, 1, &naxes, &mut status);
            fits_write_imgnull_flt(
                f.as_deref_mut().unwrap(),
                1,
                1,
                5,
                &data,
                nulval,
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0, "setup");

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut result = [0f32; 5];
            let mut nularray = [0 as c_char; 5];
            let mut anynull = -1;
            fits_read_imgnull_flt(
                f.as_deref_mut().unwrap(),
                1,
                1,
                5,
                &mut result,
                &mut nularray,
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[0], 1.0);
            assert_eq!(result[1], 2.0);
            assert_eq!(nularray[2], 1); // Should be null.
            assert_eq!(result[3], 4.0);
            assert_eq!(anynull, 1);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_2d_array() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 2] = [3, 3];
            let data: [f32; 9] = [1.1, 2.2, 3.3, 4.4, 5.5, 6.6, 7.7, 8.8, 9.9];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), FLOAT_IMG, 2, &naxes, &mut status);
            fits_write_2d_flt(f.as_deref_mut().unwrap(), 1, 3, 3, 3, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut result = [0f32; 9];
            let mut anynull = -1;
            fits_read_img_flt(
                f.as_deref_mut().unwrap(),
                1,
                1,
                9,
                0.0,
                &mut result,
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert!((result[0] - 1.1).abs() <= 1e-5);
            assert!((result[8] - 9.9).abs() <= 1e-5);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_3d_array() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 3] = [2, 2, 2];
            let data: [f32; 8] = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), FLOAT_IMG, 3, &naxes, &mut status);
            fits_write_3d_flt(
                f.as_deref_mut().unwrap(),
                1,
                2,
                2,
                2,
                2,
                2,
                &data,
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut result = [0f32; 8];
            let mut anynull = -1;
            fits_read_img_flt(
                f.as_deref_mut().unwrap(),
                1,
                1,
                8,
                0.0,
                &mut result,
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[0], 1.0);
            assert_eq!(result[7], 8.0);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_subsection() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 2] = [4, 4];
            let fpixel: [c_long; 2] = [2, 2];
            let lpixel: [c_long; 2] = [3, 3];
            let data: [f32; 4] = [10.5, 20.5, 30.5, 40.5];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), FLOAT_IMG, 2, &naxes, &mut status);
            fits_write_subset_flt(
                f.as_deref_mut().unwrap(),
                1,
                2,
                &naxes,
                &fpixel,
                &lpixel,
                &data,
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut result = [0f32; 16];
            let mut anynull = -1;
            fits_read_img_flt(
                f.as_deref_mut().unwrap(),
                1,
                1,
                16,
                0.0,
                &mut result,
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert!((result[5] - 10.5).abs() <= 1e-5);
            assert!((result[6] - 20.5).abs() <= 1e-5);
            assert!((result[9] - 30.5).abs() <= 1e-5);
            assert!((result[10] - 40.5).abs() <= 1e-5);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_group_parameters() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 1] = [10];
            let params: [f32; 2] = [3.14159, 2.71828];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_grphdr(
                f.as_deref_mut().unwrap(),
                TRUE as c_int,
                FLOAT_IMG,
                1,
                &naxes,
                2,
                1,
                TRUE as c_int,
                &mut status,
            );
            fits_write_grppar_flt(f.as_deref_mut().unwrap(), 1, 1, 2, &params, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut result = [0f32; 2];
            fits_read_grppar_flt(f.as_deref_mut().unwrap(), 1, 1, 2, &mut result, &mut status);
            assert_eq!(status, 0);
            assert!((result[0] - 3.14159).abs() <= 1e-5);
            assert!((result[1] - 2.71828).abs() <= 1e-5);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_binary_table_column() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [f32; 3] = [-1e30, 0.0, 1e30];

            let mut f = make_table(&name, BINARY_TBL, "FLOATCOL", "1E", 3, &mut status);
            fits_write_col_flt(f.as_deref_mut().unwrap(), 1, 1, 1, 3, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0f32; 3];
            let mut anynull = -1;
            fits_read_col_flt(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                3,
                0.0,
                &mut result,
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[0], data[0]);
            assert_eq!(result[1], 0.0);
            assert_eq!(result[2], data[2]);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_column_with_nulls() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [f32; 5] = [1.0, -999.0, 3.0, -999.0, 5.0];
            let nulval: f32 = -999.0;

            let mut f = make_table(&name, BINARY_TBL, "FLOATCOL", "1E", 5, &mut status);
            fits_write_colnull_flt(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                5,
                &data,
                nulval,
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0f32; 5];
            let mut nularray = [0 as c_char; 5];
            let mut anynull = -1;
            fits_read_colnull_flt(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                5,
                &mut result,
                &mut nularray,
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[0], 1.0);
            assert_eq!(nularray[1], 1);
            assert_eq!(result[2], 3.0);
            assert_eq!(nularray[3], 1);
            assert_eq!(result[4], 5.0);
            assert_eq!(anynull, 1);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_to_byte_column() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [f32; 3] = [0.0, 127.0, 255.0];

            let mut f = make_table(&name, BINARY_TBL, "BYTECOL", "1B", 3, &mut status);
            fits_write_col_flt(f.as_deref_mut().unwrap(), 1, 1, 1, 3, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0u8; 3];
            let mut anynull = -1;
            fits_read_col_byt(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                3,
                0,
                &mut result,
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[0], 0);
            assert_eq!(result[1], 127);
            assert_eq!(result[2], 255);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_to_short_column() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [f32; 3] = [-32768.0, 0.0, 32767.0];

            let mut f = make_table(&name, BINARY_TBL, "SHORTCOL", "1I", 3, &mut status);
            fits_write_col_flt(f.as_deref_mut().unwrap(), 1, 1, 1, 3, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0i16; 3];
            let mut anynull = -1;
            fits_read_col_sht(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                3,
                0,
                &mut result,
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[0], -32768);
            assert_eq!(result[1], 0);
            assert_eq!(result[2], 32767);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_to_long_column() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [f32; 3] = [-1000000.0, 0.0, 1000000.0];

            let mut f = make_table(&name, BINARY_TBL, "LONGCOL", "1J", 3, &mut status);
            fits_write_col_flt(f.as_deref_mut().unwrap(), 1, 1, 1, 3, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0 as c_long; 3];
            let mut anynull = -1;
            fits_read_col_lng(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                3,
                0,
                &mut result,
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[0], -1000000);
            assert_eq!(result[1], 0);
            assert_eq!(result[2], 1000000);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_to_double_column() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [f32; 3] = [-1.5, 0.0, 1.5];

            let mut f = make_table(&name, BINARY_TBL, "DBLCOL", "1D", 3, &mut status);
            fits_write_col_flt(f.as_deref_mut().unwrap(), 1, 1, 1, 3, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0f64; 3];
            let mut anynull = -1;
            fits_read_col_dbl(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                3,
                0.0,
                &mut result,
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert!((result[0] - (-1.5)).abs() <= 1e-6);
            assert_eq!(result[1], 0.0);
            assert!((result[2] - 1.5).abs() <= 1e-6);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_to_longlong_column() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [f32; 3] = [-1e10, 0.0, 1e10];

            let mut f = make_table(&name, BINARY_TBL, "LLCOL", "1K", 3, &mut status);
            fits_write_col_flt(f.as_deref_mut().unwrap(), 1, 1, 1, 3, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0 as LONGLONG; 3];
            let mut anynull = -1;
            fits_read_col_lnglng(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                3,
                0,
                &mut result,
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[1], 0);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_to_ascii_table() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [f32; 3] = [-1.23456e10, 0.0, 9.87654e-5];

            let mut f = make_table(&name, ASCII_TBL, "NUMCOL", "E15.7", 3, &mut status);
            fits_write_col_flt(f.as_deref_mut().unwrap(), 1, 1, 1, 3, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0f32; 3];
            let mut anynull = -1;
            fits_read_col_flt(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                3,
                0.0,
                &mut result,
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[1], 0.0);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_with_scaling() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [f32; 1] = [100.0]; // Stored as (100 - 10) / 2 = 45.

            let mut f = make_table(&name, BINARY_TBL, "SCALED", "1E", 1, &mut status);
            fits_write_key_dbl(
                f.as_deref_mut().unwrap(),
                &cc("TSCAL1"),
                2.0,
                15,
                None,
                &mut status,
            );
            fits_write_key_dbl(
                f.as_deref_mut().unwrap(),
                &cc("TZERO1"),
                10.0,
                15,
                None,
                &mut status,
            );
            fits_write_col_flt(f.as_deref_mut().unwrap(), 1, 1, 1, 1, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0f32; 1];
            let mut anynull = -1;
            fits_read_col_flt(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                1,
                0.0,
                &mut result,
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[0], 100.0);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_complex_column() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [f32; 4] = [1.5, 2.5, 3.5, 4.5]; // Two complex numbers.

            let mut f = make_table(&name, BINARY_TBL, "CPLXCOL", "1C", 2, &mut status);
            fits_write_col_cmp(f.as_deref_mut().unwrap(), 1, 1, 1, 2, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0f32; 4];
            let mut anynull = -1;
            fits_read_col_flt(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                4,
                0.0,
                &mut result,
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[0], 1.5);
            assert_eq!(result[1], 2.5);
            assert_eq!(result[2], 3.5);
            assert_eq!(result[3], 4.5);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_multirow() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [f32; 9] = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0];

            let mut f = make_table(&name, BINARY_TBL, "FLTVEC", "3E", 3, &mut status);
            fits_write_col_flt(f.as_deref_mut().unwrap(), 1, 1, 1, 9, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0f32; 9];
            let mut anynull = -1;
            fits_read_col_flt(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                9,
                0.0,
                &mut result,
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[0], 1.0);
            assert_eq!(result[4], 5.0);
            assert_eq!(result[8], 9.0);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_overflow_to_byte() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [f32; 2] = [-1.0, 256.0]; // Both overflow.

            let mut f = make_table(&name, BINARY_TBL, "BYTECOL", "1B", 2, &mut status);
            fits_write_col_flt(f.as_deref_mut().unwrap(), 1, 1, 1, 2, &data, &mut status);
            assert_eq!(status, NUM_OVERFLOW);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_overflow_to_short() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [f32; 2] = [-40000.0, 40000.0]; // Both overflow.

            let mut f = make_table(&name, BINARY_TBL, "SHORTCOL", "1I", 2, &mut status);
            fits_write_col_flt(f.as_deref_mut().unwrap(), 1, 1, 1, 2, &data, &mut status);
            assert_eq!(status, NUM_OVERFLOW);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_overflow_to_long() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [f32; 2] = [-3e9, 3e9]; // Both overflow 32-bit.

            let mut f = make_table(&name, BINARY_TBL, "LONGCOL", "1J", 2, &mut status);
            fits_write_col_flt(f.as_deref_mut().unwrap(), 1, 1, 1, 2, &data, &mut status);
            assert_eq!(status, NUM_OVERFLOW);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_overflow_to_longlong() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [f32; 1] = [1e19]; // Overflows LONGLONG_MAX.

            let mut f = make_table(&name, BINARY_TBL, "LLCOL", "1K", 1, &mut status);
            fits_write_col_flt(f.as_deref_mut().unwrap(), 1, 1, 1, 1, &data, &mut status);
            assert_eq!(status, NUM_OVERFLOW);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_scaled_overflow_to_byte() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [f32; 1] = [100.0]; // With scale=0.1, 100/0.1 = 1000 > 255.

            let mut f = make_table(&name, BINARY_TBL, "BYTECOL", "1B", 1, &mut status);
            fits_set_tscale(f.as_deref_mut().unwrap(), 1, 0.1, 0.0, &mut status);
            fits_write_col_flt(f.as_deref_mut().unwrap(), 1, 1, 1, 1, &data, &mut status);
            assert_eq!(status, NUM_OVERFLOW);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_scaled_overflow_to_short() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [f32; 1] = [100.0]; // With scale=0.001, 100/0.001 = 100000 > SHRT_MAX.

            let mut f = make_table(&name, BINARY_TBL, "SHORTCOL", "1I", 1, &mut status);
            fits_set_tscale(f.as_deref_mut().unwrap(), 1, 0.001, 0.0, &mut status);
            fits_write_col_flt(f.as_deref_mut().unwrap(), 1, 1, 1, 1, &data, &mut status);
            assert_eq!(status, NUM_OVERFLOW);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_scaled_overflow_to_long() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [f32; 1] = [1000.0]; // With scale=1e-7, 1000/1e-7 = 1e10 > INT32_MAX.

            let mut f = make_table(&name, BINARY_TBL, "LONGCOL", "1J", 1, &mut status);
            fits_set_tscale(f.as_deref_mut().unwrap(), 1, 0.0000001, 0.0, &mut status);
            fits_write_col_flt(f.as_deref_mut().unwrap(), 1, 1, 1, 1, &data, &mut status);
            assert_eq!(status, NUM_OVERFLOW);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_scaled_overflow_to_longlong() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [f32; 1] = [1e10]; // With scale=1e-18, 1e10/1e-18 = 1e28 > LONGLONG_MAX.

            let mut f = make_table(&name, BINARY_TBL, "LLCOL", "1K", 1, &mut status);
            fits_set_tscale(f.as_deref_mut().unwrap(), 1, 1e-18, 0.0, &mut status);
            fits_write_col_flt(f.as_deref_mut().unwrap(), 1, 1, 1, 1, &data, &mut status);
            assert_eq!(status, NUM_OVERFLOW);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_bad_dimen_3d() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 3] = [4, 4, 2];
            let data = [0f32; 32];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), FLOAT_IMG, 3, &naxes, &mut status);
            // ncols=2 < naxis1=4, should fail with BAD_DIMEN.
            fits_write_3d_flt(
                f.as_deref_mut().unwrap(),
                1,
                2,
                4,
                4,
                4,
                2,
                &data,
                &mut status,
            );
            assert_eq!(status, BAD_DIMEN);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_bad_dimen_subsection() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 2] = [2, 2];
            let fpixel: [c_long; 2] = [1, 1];
            let lpixel: [c_long; 2] = [2, 2];
            let data = [0f32; 4];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), FLOAT_IMG, 2, &naxes, &mut status);
            // naxis=8 > 7, should fail with BAD_DIMEN.
            fits_write_subset_flt(
                f.as_deref_mut().unwrap(),
                1,
                8,
                &naxes,
                &fpixel,
                &lpixel,
                &data,
                &mut status,
            );
            assert_eq!(status, BAD_DIMEN);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_variable_length_array() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data1: [f32; 3] = [1.1, 2.2, 3.3];
            let data2: [f32; 2] = [10.0, 20.0];

            let mut f = make_table(&name, BINARY_TBL, "VARCOL", "1PE", 2, &mut status);
            fits_write_col_flt(f.as_deref_mut().unwrap(), 1, 1, 1, 3, &data1, &mut status);
            fits_write_col_flt(f.as_deref_mut().unwrap(), 1, 2, 1, 2, &data2, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut nelem: c_long = 0;
            fits_read_descript(
                f.as_deref_mut().unwrap(),
                1,
                1,
                Some(&mut nelem),
                None,
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(nelem, 3);
            let mut result = [0f32; 3];
            let mut anynull = -1;
            fits_read_col_flt(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                3,
                0.0,
                &mut result,
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert!((result[0] - 1.1).abs() <= 1e-5);
            assert!((result[2] - 3.3).abs() <= 1e-5);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_varlen_with_nulls() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [f32; 3] = [1.0, -999.0, 3.0];
            let nulval: f32 = -999.0;

            let mut f = make_table(&name, BINARY_TBL, "VARCOL", "1PE", 1, &mut status);
            fits_write_colnull_flt(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                3,
                &data,
                nulval,
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0f32; 3];
            let mut nularray = [0 as c_char; 3];
            let mut anynull = -1;
            fits_read_colnull_flt(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                3,
                &mut result,
                &mut nularray,
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[0], 1.0);
            assert_eq!(nularray[1], 1);
            assert_eq!(result[2], 3.0);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_nan_values() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let nan_val = f32::NAN;
            let data: [f32; 3] = [1.0, nan_val, 3.0];

            let mut f = make_table(&name, BINARY_TBL, "FLOATCOL", "1E", 3, &mut status);
            fits_write_col_flt(f.as_deref_mut().unwrap(), 1, 1, 1, 3, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0f32; 3];
            let mut nularray = [0 as c_char; 3];
            let mut anynull = -1;
            fits_read_colnull_flt(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                3,
                &mut result,
                &mut nularray,
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[0], 1.0);
            assert_eq!(nularray[1], 1); // NaN should be flagged as null.
            assert_eq!(result[2], 3.0);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_2d_noncontiguous() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 2] = [2, 3];
            // Input array is 4 columns wide but we only write 2 columns.
            let data: [f32; 12] = [
                1.0, 2.0, 99.0, 99.0, 3.0, 4.0, 99.0, 99.0, 5.0, 6.0, 99.0, 99.0,
            ];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), FLOAT_IMG, 2, &naxes, &mut status);
            // ncols=4, naxis1=2, naxis2=3.
            fits_write_2d_flt(f.as_deref_mut().unwrap(), 1, 4, 2, 3, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut result = [0f32; 6];
            let mut anynull = -1;
            fits_read_img_flt(
                f.as_deref_mut().unwrap(),
                1,
                1,
                6,
                0.0,
                &mut result,
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[0], 1.0);
            assert_eq!(result[1], 2.0);
            assert_eq!(result[2], 3.0);
            assert_eq!(result[3], 4.0);
            assert_eq!(result[4], 5.0);
            assert_eq!(result[5], 6.0);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_3d_noncontiguous() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 3] = [2, 2, 2];
            // nrows=3 > naxis2=2, so there's padding between planes.
            let data: [f32; 12] = [
                1.0, 2.0, 3.0, 4.0, 99.0, 99.0, // Padding row.
                5.0, 6.0, 7.0, 8.0, 99.0, 99.0, // Padding row.
            ];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), FLOAT_IMG, 3, &naxes, &mut status);
            // ncols=2, nrows=3, naxis1=2, naxis2=2, naxis3=2.
            fits_write_3d_flt(
                f.as_deref_mut().unwrap(),
                1,
                2,
                3,
                2,
                2,
                2,
                &data,
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut result = [0f32; 8];
            let mut anynull = -1;
            fits_read_img_flt(
                f.as_deref_mut().unwrap(),
                1,
                1,
                8,
                0.0,
                &mut result,
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[0], 1.0);
            assert_eq!(result[3], 4.0);
            assert_eq!(result[4], 5.0);
            assert_eq!(result[7], 8.0);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_7d_subsection() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 7] = [2, 2, 2, 1, 1, 1, 1];
            let fpixel: [c_long; 7] = [1, 1, 1, 1, 1, 1, 1];
            let lpixel: [c_long; 7] = [2, 2, 2, 1, 1, 1, 1];
            let data: [f32; 8] = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), FLOAT_IMG, 7, &naxes, &mut status);
            fits_write_subset_flt(
                f.as_deref_mut().unwrap(),
                1,
                7,
                &naxes,
                &fpixel,
                &lpixel,
                &data,
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut result = [0f32; 8];
            let mut anynull = -1;
            fits_read_img_flt(
                f.as_deref_mut().unwrap(),
                1,
                1,
                8,
                0.0,
                &mut result,
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            for i in 0..8 {
                assert_eq!(result[i], data[i]);
            }
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_to_unsigned_longlong_column() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [f32; 3] = [0.0, 1e10, 1e15];

            let mut f = make_table(&name, BINARY_TBL, "ULLCOL", "1K", 3, &mut status);
            // Set TZERO to 2^63 to make it unsigned longlong.
            fits_set_tscale(
                f.as_deref_mut().unwrap(),
                1,
                1.0,
                9223372036854775808.0,
                &mut status,
            );
            fits_write_col_flt(f.as_deref_mut().unwrap(), 1, 1, 1, 3, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0 as LONGLONG; 3];
            let mut anynull = -1;
            fits_read_col_lnglng(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                3,
                0,
                &mut result,
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_unsigned_longlong_negative_overflow() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [f32; 1] = [-1.0]; // Negative overflows unsigned.

            let mut f = make_table(&name, BINARY_TBL, "ULLCOL", "1K", 1, &mut status);
            fits_set_tscale(
                f.as_deref_mut().unwrap(),
                1,
                1.0,
                9223372036854775808.0,
                &mut status,
            );
            fits_write_col_flt(f.as_deref_mut().unwrap(), 1, 1, 1, 1, &data, &mut status);
            assert_eq!(status, NUM_OVERFLOW);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_ascii_table_overflow() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            // Value 12345.67 needs more than 6 characters.
            let data: [f32; 1] = [12345.67];

            let mut f = make_table(&name, ASCII_TBL, "NUMCOL", "F6.2", 1, &mut status);
            fits_write_col_flt(f.as_deref_mut().unwrap(), 1, 1, 1, 1, &data, &mut status);
            assert_eq!(status, NUM_OVERFLOW);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_bad_btable_format() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [f32; 1] = [1.0];

            let mut f = make_table(&name, BINARY_TBL, "LOGCOL", "1L", 1, &mut status);
            fits_write_col_flt(f.as_deref_mut().unwrap(), 1, 1, 1, 1, &data, &mut status);
            assert_eq!(status, BAD_BTABLE_FORMAT);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_nulls_at_end() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [f32; 4] = [1.0, 2.0, -999.0, -999.0];
            let nulval: f32 = -999.0;

            let mut f = make_table(&name, BINARY_TBL, "FLOATCOL", "1E", 4, &mut status);
            fits_write_colnull_flt(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                4,
                &data,
                nulval,
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0f32; 4];
            let mut nularray = [0 as c_char; 4];
            let mut anynull = -1;
            fits_read_colnull_flt(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                4,
                &mut result,
                &mut nularray,
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[0], 1.0);
            assert_eq!(result[1], 2.0);
            assert_eq!(nularray[2], 1);
            assert_eq!(nularray[3], 1);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_all_nulls() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [f32; 3] = [-999.0, -999.0, -999.0];
            let nulval: f32 = -999.0;

            let mut f = make_table(&name, BINARY_TBL, "FLOATCOL", "1E", 3, &mut status);
            fits_write_colnull_flt(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                3,
                &data,
                nulval,
                &mut status,
            );
            assert_eq!(status, 0);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_large_vector() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let mut data = vec![0f32; 10000];
            for (i, d) in data.iter_mut().enumerate() {
                *d = i as f32 * 0.001;
            }

            let mut f = make_table(&name, BINARY_TBL, "BIGCOL", "10000E", 1, &mut status);
            fits_write_col_flt(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                10000,
                &data,
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = vec![0f32; 10000];
            let mut anynull = -1;
            fits_read_col_flt(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                10000,
                0.0,
                &mut result,
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            for i in 0..10000 {
                assert!((result[i] - data[i]).abs() <= 1e-6);
            }
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_pcne_complex_column() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [f32; 6] = [1.0, 2.0, -999.0, -999.0, 5.0, 6.0]; // 3 float complex.
            let nulval: f32 = -999.0;

            let mut f = make_table(&name, BINARY_TBL, "CPLXCOL", "3C", 1, &mut status);
            fits_write_colnull_flt(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                6,
                &data,
                nulval,
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0f32; 6];
            let mut nularray = [0 as c_char; 6];
            let mut anynull = -1;
            fits_read_colnull_flt(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                6,
                &mut result,
                &mut nularray,
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[0], 1.0);
            assert_eq!(nularray[2], 1);
            assert_eq!(result[4], 5.0);
            assert_eq!(anynull, 1);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_scaled_underflow_to_byte() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [f32; 1] = [-10.0];

            let mut f = make_table(&name, BINARY_TBL, "BYTECOL", "1B", 1, &mut status);
            fits_set_tscale(f.as_deref_mut().unwrap(), 1, 0.1, 10.0, &mut status);
            fits_write_col_flt(f.as_deref_mut().unwrap(), 1, 1, 1, 1, &data, &mut status);
            assert_eq!(status, NUM_OVERFLOW);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_scaled_normal_to_byte() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [f32; 1] = [50.0];

            let mut f = make_table(&name, BINARY_TBL, "BYTECOL", "1B", 1, &mut status);
            fits_set_tscale(f.as_deref_mut().unwrap(), 1, 0.5, 0.0, &mut status);
            fits_write_col_flt(f.as_deref_mut().unwrap(), 1, 1, 1, 1, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0 as c_uint; 1];
            let mut anynull = -1;
            fits_read_col_uint(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                1,
                0,
                &mut result,
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[0], 100);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_scaled_underflow_to_short() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [f32; 1] = [-100.0];

            let mut f = make_table(&name, BINARY_TBL, "SHORTCOL", "1I", 1, &mut status);
            fits_set_tscale(f.as_deref_mut().unwrap(), 1, 0.001, 0.0, &mut status);
            fits_write_col_flt(f.as_deref_mut().unwrap(), 1, 1, 1, 1, &data, &mut status);
            assert_eq!(status, NUM_OVERFLOW);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_scaled_normal_to_short() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [f32; 2] = [50.0, -50.0];

            let mut f = make_table(&name, BINARY_TBL, "SHORTCOL", "1I", 2, &mut status);
            fits_set_tscale(f.as_deref_mut().unwrap(), 1, 0.1, 0.0, &mut status);
            fits_write_col_flt(f.as_deref_mut().unwrap(), 1, 1, 1, 2, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0i16; 2];
            let mut anynull = -1;
            fits_read_col_sht(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                2,
                0,
                &mut result,
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[0], 500);
            assert_eq!(result[1], -500);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_scaled_underflow_to_long() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [f32; 1] = [-1000.0];

            let mut f = make_table(&name, BINARY_TBL, "LONGCOL", "1J", 1, &mut status);
            fits_set_tscale(f.as_deref_mut().unwrap(), 1, 0.0000001, 0.0, &mut status);
            fits_write_col_flt(f.as_deref_mut().unwrap(), 1, 1, 1, 1, &data, &mut status);
            assert_eq!(status, NUM_OVERFLOW);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_scaled_normal_to_long() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [f32; 2] = [100.0, -100.0];

            let mut f = make_table(&name, BINARY_TBL, "LONGCOL", "1J", 2, &mut status);
            fits_set_tscale(f.as_deref_mut().unwrap(), 1, 0.01, 0.0, &mut status);
            fits_write_col_flt(f.as_deref_mut().unwrap(), 1, 1, 1, 2, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0 as c_long; 2];
            let mut anynull = -1;
            fits_read_col_lng(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                2,
                0,
                &mut result,
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[0], 10000);
            assert_eq!(result[1], -10000);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_scaled_underflow_to_longlong() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [f32; 1] = [-100.0];

            let mut f = make_table(&name, BINARY_TBL, "LLCOL", "1K", 1, &mut status);
            fits_set_tscale(f.as_deref_mut().unwrap(), 1, 1e-18, 0.0, &mut status);
            fits_write_col_flt(f.as_deref_mut().unwrap(), 1, 1, 1, 1, &data, &mut status);
            assert_eq!(status, NUM_OVERFLOW);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_scaled_normal_to_longlong() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [f32; 2] = [100.0, -100.0];

            let mut f = make_table(&name, BINARY_TBL, "LLCOL", "1K", 2, &mut status);
            fits_set_tscale(f.as_deref_mut().unwrap(), 1, 0.001, 0.0, &mut status);
            fits_write_col_flt(f.as_deref_mut().unwrap(), 1, 1, 1, 2, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0 as LONGLONG; 2];
            let mut anynull = -1;
            fits_read_col_lnglng(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                2,
                0,
                &mut result,
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[0], 100000);
            assert_eq!(result[1], -100000);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_scaled_write_to_double_column() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [f32; 1] = [100.0];

            let mut f = make_table(&name, BINARY_TBL, "DBLCOL", "1D", 1, &mut status);
            fits_set_tscale(f.as_deref_mut().unwrap(), 1, 2.0, 10.0, &mut status);
            fits_write_col_flt(f.as_deref_mut().unwrap(), 1, 1, 1, 1, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0f64; 1];
            let mut anynull = -1;
            fits_read_col_dbl(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                1,
                0.0,
                &mut result,
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert!((result[0] - 45.0).abs() <= 0.01);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_scaled_ascii_column() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [f32; 1] = [100.0];

            let mut f = make_table(&name, ASCII_TBL, "VALUE", "E15.6", 1, &mut status);
            fits_set_tscale(f.as_deref_mut().unwrap(), 1, 2.0, 10.0, &mut status);
            fits_write_col_flt(f.as_deref_mut().unwrap(), 1, 1, 1, 1, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0f32; 1];
            let mut anynull = -1;
            fits_read_col_flt(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                1,
                0.0,
                &mut result,
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(anynull, 0);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_unsigned_longlong_positive_overflow() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [f32; 1] = [2e19]; // > 2 * 9.2e18.

            let mut f = make_table(&name, BINARY_TBL, "ULLCOL", "1K", 1, &mut status);
            fits_set_tscale(
                f.as_deref_mut().unwrap(),
                1,
                1.0,
                9223372036854775808.0,
                &mut status,
            );
            fits_write_col_flt(f.as_deref_mut().unwrap(), 1, 1, 1, 1, &data, &mut status);
            assert_eq!(status, NUM_OVERFLOW);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_scaled_ascii_overflow() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            // With scale=0.01, stored = 100/0.01 = 10000, which needs >6 chars.
            let data: [f32; 1] = [100.0];

            let mut f = make_table(&name, ASCII_TBL, "NUMCOL", "F6.2", 1, &mut status);
            fits_set_tscale(f.as_deref_mut().unwrap(), 1, 0.01, 0.0, &mut status);
            fits_write_col_flt(f.as_deref_mut().unwrap(), 1, 1, 1, 1, &data, &mut status);
            assert_eq!(status, NUM_OVERFLOW);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_unscaled_longlong_underflow() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [f32; 1] = [-1e19]; // < LONGLONG_MIN.

            let mut f = make_table(&name, BINARY_TBL, "LLCOL", "1K", 1, &mut status);
            // No scaling - scale=1, zero=0.
            fits_write_col_flt(f.as_deref_mut().unwrap(), 1, 1, 1, 1, &data, &mut status);
            assert_eq!(status, NUM_OVERFLOW);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }
}
