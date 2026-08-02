/*  This file, putcolk.rs, contains routines that write data elements to    */
/*  a FITS image or table, with 'int' datatype.                            */

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
use crate::putcoli::ffpcli_safe;
use crate::putcolj::ffpclj_safe;
use crate::putcolu::ffpclu_safe;
use crate::relibc::header::stdio::snprintf_f64;
use crate::wrappers::strlen_safe;
use crate::{buffers::*, int_snprintf, slice_to_str};

/*--------------------------------------------------------------------------*/
/// Write an array of values to the primary array. Data conversion
/// and scaling will be performed if necessary (e.g, if the datatype of
/// the FITS array is not the same as the array being written).
///
/// The primary array is represented as a binary table:
/// each group of the primary array is a row in the table,
/// where the first column contains the group parameters
/// and the second column contains the image itself.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffpprk(
    fptr: *mut fitsfile, /* I - FITS file pointer                       */
    group: c_long,       /* I - group to write(1 = 1st group)           */
    firstelem: LONGLONG, /* I - first vector element to write(1 = 1st)  */
    nelem: LONGLONG,     /* I - number of values to write               */
    array: *const c_int, /* I - array of values that are written        */
    status: *mut c_int,  /* IO - error status                           */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        let array = slice::from_raw_parts(array, nelem as usize);

        ffpprk_safe(fptr, group, firstelem, nelem, array, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Write an array of values to the primary array. Data conversion
/// and scaling will be performed if necessary (e.g, if the datatype of
/// the FITS array is not the same as the array being written).
///
/// The primary array is represented as a binary table:
/// each group of the primary array is a row in the table,
/// where the first column contains the group parameters
/// and the second column contains the image itself.
pub fn ffpprk_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer                       */
    group: c_long,       /* I - group to write(1 = 1st group)           */
    firstelem: LONGLONG, /* I - first vector element to write(1 = 1st)  */
    nelem: LONGLONG,     /* I - number of values to write               */
    array: &[c_int],     /* I - array of values that are written        */
    status: &mut c_int,  /* IO - error status                           */
) -> c_int {
    let nullvalue: c_int = 0;

    if fits_is_compressed_image_safe(fptr, status) > 0 {
        /* this is a compressed image in a binary table */
        fits_write_compressed_pixels(
            fptr,
            TINT,
            firstelem,
            nelem,
            NullCheckType::None,
            cast_slice(array),
            &Some(NullValue::Int(nullvalue)),
            status,
        );
        return *status;
    }

    let row = cmp::max(1, group);

    ffpclk_safe(
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
/// The primary array is represented as a binary table:
/// each group of the primary array is a row in the table,
/// where the first column contains the group parameters
/// and the second column contains the image itself.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffppnk(
    fptr: *mut fitsfile, /* I - FITS file pointer                       */
    group: c_long,       /* I - group to write(1 = 1st group)           */
    firstelem: LONGLONG, /* I - first vector element to write(1 = 1st)  */
    nelem: LONGLONG,     /* I - number of values to write               */
    array: *const c_int, /* I - array of values that are written        */
    nulval: c_int,       /* I - undefined pixel value                   */
    status: *mut c_int,  /* IO - error status                           */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        let array = slice::from_raw_parts(array, nelem as usize);

        ffppnk_safe(fptr, group, firstelem, nelem, array, nulval, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Write an array of values to the primary array. Data conversion
/// and scaling will be performed if necessary (e.g, if the datatype of the
/// FITS array is not the same as the array being written).  Any array values
/// that are equal to the value of nulval will be replaced with the null
/// pixel value that is appropriate for this column.
///
/// The primary array is represented as a binary table:
/// each group of the primary array is a row in the table,
/// where the first column contains the group parameters
/// and the second column contains the image itself.
pub fn ffppnk_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer                       */
    group: c_long,       /* I - group to write(1 = 1st group)           */
    firstelem: LONGLONG, /* I - first vector element to write(1 = 1st)  */
    nelem: LONGLONG,     /* I - number of values to write               */
    array: &[c_int],     /* I - array of values that are written        */
    nulval: c_int,       /* I - undefined pixel value                   */
    status: &mut c_int,  /* IO - error status                           */
) -> c_int {
    let mut nullvalue: c_int = 0;

    if fits_is_compressed_image_safe(fptr, status) > 0 {
        /* this is a compressed image in a binary table */

        nullvalue = nulval; /* set local variable */
        fits_write_compressed_pixels(
            fptr,
            TINT,
            firstelem,
            nelem,
            NullCheckType::SetPixel,
            cast_slice(array),
            &Some(NullValue::Int(nullvalue)),
            status,
        );
        return *status;
    }

    let row = cmp::max(1, group);

    ffpcnk_safe(
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
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffp2dk(
    fptr: *mut fitsfile, /* I - FITS file pointer                     */
    group: c_long,       /* I - group to write(1 = 1st group)         */
    ncols: LONGLONG,     /* I - number of pixels in each row of array */
    naxis1: LONGLONG,    /* I - FITS image NAXIS1 value               */
    naxis2: LONGLONG,    /* I - FITS image NAXIS2 value               */
    array: *const c_int, /* I - array to be written                   */
    status: *mut c_int,  /* IO - error status                         */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        let array = slice::from_raw_parts(array, (ncols * naxis2) as usize);

        ffp2dk_safe(fptr, group, ncols, naxis1, naxis2, array, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Write an entire 2-D array of values to the primary array. Data conversion
/// and scaling will be performed if necessary (e.g, if the datatype of the
/// FITS array is not the same as the array being written).
pub fn ffp2dk_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer                     */
    group: c_long,       /* I - group to write(1 = 1st group)         */
    ncols: LONGLONG,     /* I - number of pixels in each row of array */
    naxis1: LONGLONG,    /* I - FITS image NAXIS1 value               */
    naxis2: LONGLONG,    /* I - FITS image NAXIS2 value               */
    array: &[c_int],     /* I - array to be written                   */
    status: &mut c_int,  /* IO - error status                         */
) -> c_int {
    /* call the 3D writing routine, with the 3rd dimension = 1 */
    ffp3dk_safe(fptr, group, ncols, naxis2, naxis1, naxis2, 1, array, status);

    *status
}

/*--------------------------------------------------------------------------*/
/// Write an entire 3-D cube of values to the primary array. Data conversion
/// and scaling will be performed if necessary (e.g, if the datatype of the
/// FITS array is not the same as the array being written).
///
/// The primary array is represented as a binary table:
/// each group of the primary array is a row in the table,
/// where the first column contains the group parameters
/// and the second column contains the image itself.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffp3dk(
    fptr: *mut fitsfile, /* I - FITS file pointer                     */
    group: c_long,       /* I - group to write(1 = 1st group)         */
    ncols: LONGLONG,     /* I - number of pixels in each row of array */
    nrows: LONGLONG,     /* I - number of rows in each plane of array */
    naxis1: LONGLONG,    /* I - FITS image NAXIS1 value               */
    naxis2: LONGLONG,    /* I - FITS image NAXIS2 value               */
    naxis3: LONGLONG,    /* I - FITS image NAXIS3 value               */
    array: *const c_int, /* I - array to be written                   */
    status: *mut c_int,  /* IO - error status                         */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        let array = slice::from_raw_parts(array, (ncols * naxis2 * naxis3) as usize);

        ffp3dk_safe(
            fptr, group, ncols, nrows, naxis1, naxis2, naxis3, array, status,
        )
    }
}

/*--------------------------------------------------------------------------*/
/// Write an entire 3-D cube of values to the primary array. Data conversion
/// and scaling will be performed if necessary (e.g, if the datatype of the
/// FITS array is not the same as the array being written).
///
/// The primary array is represented as a binary table:
/// each group of the primary array is a row in the table,
/// where the first column contains the group parameters
/// and the second column contains the image itself.
pub fn ffp3dk_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer                     */
    group: c_long,       /* I - group to write(1 = 1st group)         */
    ncols: LONGLONG,     /* I - number of pixels in each row of array */
    nrows: LONGLONG,     /* I - number of rows in each plane of array */
    naxis1: LONGLONG,    /* I - FITS image NAXIS1 value               */
    naxis2: LONGLONG,    /* I - FITS image NAXIS2 value               */
    naxis3: LONGLONG,    /* I - FITS image NAXIS3 value               */
    array: &[c_int],     /* I - array to be written                   */
    status: &mut c_int,  /* IO - error status                         */
) -> c_int {
    let mut tablerow: c_long = 0;
    let fpixel: [c_long; 3] = [1, 1, 1];
    let mut lpixel: [c_long; 3] = [0, 0, 0];
    let mut nfits: LONGLONG = 0;
    let mut narray: LONGLONG = 0;

    if fits_is_compressed_image_safe(fptr, status) > 0 {
        /* this is a compressed image in a binary table */
        lpixel[0] = ncols as c_long;
        lpixel[1] = nrows as c_long;
        lpixel[2] = naxis3 as c_long;

        fits_write_compressed_img(
            fptr,
            TUINT,
            &fpixel,
            &lpixel,
            NullCheckType::None,
            cast_slice(array),
            &None,
            status,
        );

        return *status;
    }

    tablerow = cmp::max(1, group);

    if ncols == naxis1 && nrows == naxis2 {
        /* arrays have same size? */

        /* all the image pixels are contiguous, so write all at once */
        ffpclk_safe(
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

    nfits = 1; /* next pixel in FITS image to write to */
    narray = 0; /* next pixel in input array to be written */

    /* loop over naxis3 planes in the data cube */
    for _jj in 0..(naxis3 as usize) {
        /* loop over the naxis2 rows in the FITS image, */
        /* writing naxis1 pixels to each row            */
        for _ii in 0..(naxis2 as usize) {
            if ffpclk_safe(
                fptr,
                2,
                tablerow as LONGLONG,
                nfits,
                naxis1,
                &array[(narray as usize)..],
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
/// Write a subsection of pixels to the primary array or image.
///
/// A subsection is defined to be any contiguous rectangular
/// array of pixels within the n-dimensional FITS data file.
/// Data conversion and scaling will be performed if necessary
/// (e.g, if the datatype of the FITS array is not the same as
/// the array being written).
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffpssk(
    fptr: *mut fitsfile,   /* I - FITS file pointer                       */
    group: c_long,         /* I - group to write(1 = 1st group)           */
    naxis: c_long,         /* I - number of data axes in array            */
    naxes: *const c_long,  /* I - size of each FITS axis                  */
    fpixel: *const c_long, /* I - 1st pixel in each axis to write (1=1st) */
    lpixel: *const c_long, /* I - last pixel in each axis to write        */
    array: *const c_int,   /* I - array to be written                     */
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

        ffpssk_safe(fptr, group, naxis, naxes, fpixel, lpixel, array, status)
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
pub fn ffpssk_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer                       */
    group: c_long,       /* I - group to write(1 = 1st group)           */
    naxis: c_long,       /* I - number of data axes in array            */
    naxes: &[c_long],    /* I - size of each FITS axis                  */
    fpixel: &[c_long],   /* I - 1st pixel in each axis to write (1=1st) */
    lpixel: &[c_long],   /* I - last pixel in each axis to write        */
    array: &[c_int],     /* I - array to be written                     */
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
            TUINT,
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
                            if ffpclk_safe(
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
pub unsafe extern "C" fn ffpgpk(
    fptr: *mut fitsfile, /* I - FITS file pointer                      */
    group: c_long,       /* I - group to write(1 = 1st group)          */
    firstelem: c_long,   /* I - first vector element to write(1 = 1st) */
    nelem: c_long,       /* I - number of values to write              */
    array: *const c_int, /* I - array of values that are written       */
    status: *mut c_int,  /* IO - error status                          */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let array = slice::from_raw_parts(array, nelem as usize);

        ffpgpk_safe(fptr, group, firstelem, nelem, array, status)
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
pub fn ffpgpk_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer                      */
    group: c_long,       /* I - group to write(1 = 1st group)          */
    firstelem: c_long,   /* I - first vector element to write(1 = 1st) */
    nelem: c_long,       /* I - number of values to write              */
    array: &[c_int],     /* I - array of values that are written       */
    status: &mut c_int,  /* IO - error status                          */
) -> c_int {
    let row = cmp::max(1, group);

    ffpclk_safe(
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
pub unsafe extern "C" fn ffpclk(
    fptr: *mut fitsfile, /* I - FITS file pointer                       */
    colnum: c_int,       /* I - number of column to write (1 = 1st col) */
    firstrow: LONGLONG,  /* I - first row to write (1 = 1st row)        */
    firstelem: LONGLONG, /* I - first vector element to write (1 = 1st) */
    nelem: LONGLONG,     /* I - number of values to write               */
    array: *const c_int, /* I - array of values to write                */
    status: *mut c_int,  /* IO - error status                           */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        let array = slice::from_raw_parts(array, nelem as usize);

        ffpclk_safe(
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
pub fn ffpclk_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer                       */
    colnum: c_int,       /* I - number of column to write (1 = 1st col) */
    firstrow: LONGLONG,  /* I - first row to write (1 = 1st row)        */
    firstelem: LONGLONG, /* I - first vector element to write (1 = 1st) */
    nelem: LONGLONG,     /* I - number of values to write               */
    array: &[c_int],     /* I - array of values to write                */
    status: &mut c_int,  /* IO - error status                           */
) -> c_int {
    let mut tcode: c_int = 0;
    let mut maxelem2: c_int = 0;
    let mut hdutype: c_int = 0;
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
    let mut scale: f64 = 0.0;
    let mut zero: f64 = 0.0;
    let writeraw;
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

    /* call the 'short' or 'long' version of this routine, if possible */
    if c_int::BITS == c_short::BITS {
        ffpcli_safe(
            fptr,
            colnum,
            firstrow as LONGLONG,
            firstelem as LONGLONG,
            nelem as LONGLONG,
            cast_slice(array),
            status,
        );
    } else if c_int::BITS == c_long::BITS {
        ffpclj_safe(
            fptr,
            colnum,
            firstrow as LONGLONG,
            firstelem as LONGLONG,
            nelem as LONGLONG,
            cast_slice(array),
            status,
        );
    } else {
        /*
        This is a special case: sizeof(int) is not equal to sizeof(short) or
        sizeof(long).  This occurs on Alpha OSF systems where short = 2 bytes,
        int = 4 bytes, and long = 8 bytes.
        */

        /*---------------------------------------------------*/
        /*  Check input and get parameters about the column: */
        /*---------------------------------------------------*/
        if ffgcprll(
            fptr,
            colnum,
            firstrow as LONGLONG,
            firstelem as LONGLONG,
            nelem as LONGLONG,
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

        let mut maxelem = LONGLONG::from(maxelem2);

        if tcode == TSTRING {
            ffcfmt(&tform, &mut cform); /* derive C format for writing strings */
        }

        /*
           if there is no scaling and the native machine format is not byteswapped
           then we can simply write the raw data bytes into the FITS file if the
           datatype of the FITS column is the same as the input values.  Otherwise
           we must convert the raw values into the scaled and/or machine dependent
           format in a temporary buffer that has been allocated for this purpose.
        */

        if scale == 1.0 && zero == 0.0 && CFITSIO_MACHINE == NATIVE && tcode == TLONG {
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
            ntodo = cmp::min(remain, maxelem as LONGLONG) as c_long;
            ntodo = cmp::min(ntodo as LONGLONG, repeat - elemnum) as c_long;

            wrtptr = startpos + ((rownum as LONGLONG) * rowlen) + (elemnum * incre as LONGLONG);

            ffmbyt_safe(fptr, wrtptr, IGNORE_EOF, status); /* move to write position */

            match tcode {
                TLONG => {
                    if writeraw {
                        ffpi4b(
                            fptr,
                            ntodo,
                            incre,
                            cast_slice(&array[(next as usize)..]),
                            status,
                        );
                    } else {
                        /* convert the raw data before writing to FITS file */
                        ffintfi4(
                            &array[(next as usize)..],
                            ntodo,
                            scale,
                            zero,
                            &mut (cast_slice_mut(&mut buffer))[..(ntodo as usize)], //Or more correctly: &mut (cast_slice_mut(&mut buffer))[..(ntodo as usize)]
                            status,
                        );
                        ffpi4b(fptr, ntodo, incre, cast_slice_mut(&mut buffer), status);
                    }
                }
                TLONGLONG => {
                    ffintfi8(
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
                    ffintfi1(
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
                    ffintfi2(
                        &array[(next as usize)..],
                        ntodo,
                        scale,
                        zero,
                        cast_slice_mut(&mut buffer),
                        status,
                    );
                    ffpi2b(fptr, ntodo, incre, cast_slice_mut(&mut buffer), status);
                }
                TFLOAT => {
                    ffintfr4(
                        &array[(next as usize)..],
                        ntodo,
                        scale,
                        zero,
                        cast_slice_mut(&mut buffer),
                        status,
                    );
                    ffpr4b(fptr, ntodo, incre, cast_slice_mut(&mut buffer), status);
                }
                TDOUBLE => {
                    ffintfr8(
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

                    let mut handled = false;
                    if hdutype == ASCII_TBL
                        && formlen > 1
                        && (cform[formlen - 1] == bb(b'f') || cform[formlen - 1] == bb(b'E'))
                    {
                        ffintfstr(
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
                                cast_slice(&buffer),
                                status,
                            );
                        } else {
                            ffpbytoff(
                                fptr,
                                twidth,
                                ntodo,
                                incre - twidth,
                                cast_slice(&buffer),
                                status,
                            );
                        }
                        handled = true;
                    }

                    if !handled {
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
                    "Error writing elements {:.0} thru {:.0} of input data array (ffpclk).",
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
            if remain != 0 {
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
    } /* end of Dec ALPHA special case */

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
pub unsafe extern "C" fn ffpcnk(
    fptr: *mut fitsfile, /* I - FITS file pointer                       */
    colnum: c_int,       /* I - number of column to write (1 = 1st col) */
    firstrow: LONGLONG,  /* I - first row to write (1 = 1st row)        */
    firstelem: LONGLONG, /* I - first vector element to write (1 = 1st) */
    nelem: LONGLONG,     /* I - number of values to write               */
    array: *const c_int, /* I - array of values to write                */
    nulvalue: c_int,     /* I - value used to flag undefined pixels     */
    status: *mut c_int,  /* IO - error status                           */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        let array = slice::from_raw_parts(array, nelem as usize);

        ffpcnk_safe(
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
pub fn ffpcnk_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer                       */
    colnum: c_int,       /* I - number of column to write (1 = 1st col) */
    firstrow: LONGLONG,  /* I - first row to write (1 = 1st row)        */
    firstelem: LONGLONG, /* I - first vector element to write (1 = 1st) */
    nelem: LONGLONG,     /* I - number of values to write               */
    array: &[c_int],     /* I - array of values to write                */
    nulvalue: c_int,     /* I - value used to flag undefined pixels     */
    status: &mut c_int,  /* IO - error status                           */
) -> c_int {
    let mut ngood: LONGLONG = 0;
    let mut nbad: LONGLONG = 0;
    let repeat: LONGLONG;

    let mut fstelm: LONGLONG;
    let mut fstrow: LONGLONG;

    let mut tcode = 0;
    let mut overflow = 0;

    if *status > 0 {
        return *status;
    }

    /* reset position to the correct HDU if necessary */
    if fptr.HDUposition != (fptr.Fptr.curhdu) {
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

    /* if variable length array, first write the whole input vector,
    then go back and fill in the nulls */
    if tcode < 0
        && ffpclk_safe(
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

                if ffpclu_safe(fptr, colnum, fstrow, fstelm, nbad, status) > 0 {
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
                    if ffpclk_safe(
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
            ffpclk_safe(
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

        ffpclu_safe(fptr, colnum, fstrow, fstelm, nbad, status);
    }

    if *status <= 0 && overflow > 0 {
        *status = NUM_OVERFLOW;
    }

    *status
}

/*--------------------------------------------------------------------------*/
/// Copy input to output prior to writing output to a FITS file.
/// Do datatype conversion and scaling if required.
pub(crate) fn ffintfi1(
    input: &[c_int],    /* I - array of values to be converted  */
    ntodo: c_long,      /* I - number of elements in the array  */
    scale: f64,         /* I - FITS TSCALn or BSCALE value      */
    zero: f64,          /* I - FITS TZEROn or BZERO  value      */
    output: &mut [u8],  /* O - output array of converted values */
    status: &mut c_int, /* IO - error status                    */
) -> c_int {
    if scale == 1.0 && zero == 0.0 {
        for ii in 0..(ntodo as usize) {
            if input[ii] < 0 {
                *status = OVERFLOW_ERR;
                output[ii] = 0;
            } else if input[ii] > c_int::from(u8::MAX) {
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
pub(crate) fn ffintfi2(
    input: &[c_int],        /* I - array of values to be converted  */
    ntodo: c_long,          /* I - number of elements in the array  */
    scale: f64,             /* I - FITS TSCALn or BSCALE value      */
    zero: f64,              /* I - FITS TZEROn or BZERO  value      */
    output: &mut [c_short], /* O - output array of converted values */
    status: &mut c_int,     /* IO - error status                    */
) -> c_int {
    if scale == 1.0 && zero == 0.0 {
        for ii in 0..(ntodo as usize) {
            if input[ii] < c_int::from(c_short::MIN) {
                *status = OVERFLOW_ERR;
                output[ii] = c_short::MIN;
            } else if input[ii] > c_int::from(c_short::MAX) {
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
/// Do datatype conversion and scaling if required
pub(crate) fn ffintfi4(
    input: &[c_int],         /* I - array of values to be converted  */
    ntodo: c_long,           /* I - number of elements in the array  */
    scale: f64,              /* I - FITS TSCALn or BSCALE value      */
    zero: f64,               /* I - FITS TZEROn or BZERO  value      */
    output: &mut [INT32BIT], /* O - output array of converted values */
    status: &mut c_int,      /* IO - error status                    */
) -> c_int {
    assert!(output.len() >= ntodo as usize);
    assert!(input.len() >= ntodo as usize);

    if scale == 1.0 && zero == 0.0 {
        output.copy_from_slice(&input[..(ntodo as usize)]);
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
/// Do datatype conversion and scaling if required
pub(crate) fn ffintfi8(
    input: &[c_int],         /* I - array of values to be converted  */
    ntodo: c_long,           /* I - number of elements in the array  */
    scale: f64,              /* I - FITS TSCALn or BSCALE value      */
    zero: f64,               /* I - FITS TZEROn or BZERO  value      */
    output: &mut [LONGLONG], /* O - output array of converted values */
    status: &mut c_int,      /* IO - error status                    */
) -> c_int {
    if scale == 1.0 && zero == 9223372036854775808.0 {
        /* Writing to unsigned long long column. Input values must not be negative */
        /* Instead of subtracting 9223372036854775808, it is more efficient */
        /* and more precise to just flip the sign bit with the XOR operator */

        for ii in 0..(ntodo as usize) {
            if input[ii] < 0 {
                *status = OVERFLOW_ERR;
                output[ii] = LONGLONG_MIN;
            } else {
                output[ii] = (input[ii] as ULONGLONG ^ 0x8000000000000000) as LONGLONG;
            }
        }
    } else if scale == 1.0 && zero == 0.0 {
        for ii in 0..(ntodo as usize) {
            output[ii] = LONGLONG::from(input[ii]);
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
pub(crate) fn ffintfr4(
    input: &[c_int],    /* I - array of values to be converted  */
    ntodo: c_long,      /* I - number of elements in the array  */
    scale: f64,         /* I - FITS TSCALn or BSCALE value      */
    zero: f64,          /* I - FITS TZEROn or BZERO  value      */
    output: &mut [f32], /* O - output array of converted values */
    status: &mut c_int, /* IO - error status                    */
) -> c_int {
    if scale == 1.0 && zero == 0.0 {
        for ii in 0..(ntodo as usize) {
            output[ii] = input[ii] as f32;
        }
    } else {
        for ii in 0..(ntodo as usize) {
            output[ii] = ((f64::from(input[ii]) - zero) / scale) as f32;
        }
    }
    *status
}

/*--------------------------------------------------------------------------*/
/// Copy input to output prior to writing output to a FITS file.
/// Do datatype conversion and scaling if required.
pub(crate) fn ffintfr8(
    input: &[c_int],    /* I - array of values to be converted  */
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
pub(crate) fn ffintfstr(
    input: &[c_int],       /* I - array of values to be converted  */
    ntodo: c_long,         /* I - number of elements in the array  */
    scale: f64,            /* I - FITS TSCALn or BSCALE value      */
    zero: f64,             /* I - FITS TZEROn or BZERO  value      */
    cform: &[c_char],      /* I - format for output string values  */
    twidth: c_long,        /* I - width of each field, in chars    */
    output: &mut [c_char], /* O - output array of converted values */
    status: &mut c_int,    /* IO - error status                    */
) -> c_int {
    let mut oi = 0;

    if scale == 1.0 && zero == 0.0 {
        for ii in 0..(ntodo as usize) {
            snprintf_f64(
                &mut output[oi..],
                DBUFFSIZE as usize,
                cform,
                f64::from(input[ii]),
            );
            oi += twidth as usize;

            if output[oi] != 0 {
                /* if this char != \0, then overflow occurred */
                *status = OVERFLOW_ERR;
            }
        }
    } else {
        for ii in 0..(ntodo as usize) {
            let dvalue: f64 = (f64::from(input[ii]) - zero) / scale;
            snprintf_f64(&mut output[oi..], DBUFFSIZE as usize, cform, dvalue);
            oi += twidth as usize;

            if output[oi] != 0 {
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

// Ported from test_putcolk.c - signed int (type `k` = c_int) write functions.

#[cfg(test)]
mod tests {
    use crate::aliases::rust_api::*;
    use crate::fitsio::{
        ASCII_TBL, BAD_ATABLE_FORMAT, BAD_BTABLE_FORMAT, BAD_DIMEN, BINARY_TBL, BYTE_IMG, LONG_IMG,
        LONGLONG, NO_NULL, NUM_OVERFLOW, READONLY, TRUE, fitsfile,
    };
    use crate::helpers::testhelpers::{to_buf, with_temp_file};
    use libc::{c_char, c_int, c_long};

    /// Make a NUL-terminated `Vec<c_char>` from a `&str`.
    fn cc(s: &str) -> Vec<c_char> {
        let mut v: Vec<c_char> = s.bytes().map(|b| b as c_char).collect();
        v.push(0);
        v
    }

    /// Create a single-column table (binary or ASCII) and return the open file.
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
            let data: [c_int; 5] = [c_int::MIN, -1000, 0, 1000, c_int::MAX];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), LONG_IMG, 1, &naxes, &mut status);
            fits_write_img_int(f.as_deref_mut().unwrap(), 1, 1, 5, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0, "setup");

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut result = [0 as c_int; 5];
            let mut anynull = -1;
            fits_read_img_int(
                f.as_deref_mut().unwrap(),
                1,
                1,
                5,
                0,
                &mut result,
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[0], c_int::MIN);
            assert_eq!(result[1], -1000);
            assert_eq!(result[2], 0);
            assert_eq!(result[3], 1000);
            assert_eq!(result[4], c_int::MAX);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_primary_with_null() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 1] = [5];
            let data: [c_int; 5] = [10, 20, -9999, 40, 50];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), LONG_IMG, 1, &naxes, &mut status);
            fits_write_key_lng(
                f.as_deref_mut().unwrap(),
                &cc("BLANK"),
                -9999,
                None,
                &mut status,
            );
            fits_write_imgnull_int(
                f.as_deref_mut().unwrap(),
                1,
                1,
                5,
                &data,
                -9999,
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0, "setup");

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut result = [0 as c_int; 5];
            let mut nularray = [0 as c_char; 5];
            let mut anynull = -1;
            fits_read_imgnull_int(
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
            assert_eq!(result[0], 10);
            assert_eq!(result[1], 20);
            assert_eq!(nularray[2], 1);
            assert_eq!(result[3], 40);
            assert_eq!(result[4], 50);
            assert_eq!(anynull, 1);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_2d_array() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 2] = [3, 2];
            let data: [c_int; 6] = [-3, -2, -1, 1, 2, 3];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), LONG_IMG, 2, &naxes, &mut status);
            fits_write_2d_int(f.as_deref_mut().unwrap(), 1, 3, 3, 2, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut result = [0 as c_int; 6];
            let mut anynull = -1;
            fits_read_2d_int(
                f.as_deref_mut().unwrap(),
                1,
                0,
                3,
                3,
                2,
                &mut result,
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result, [-3, -2, -1, 1, 2, 3]);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_3d_array() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 3] = [2, 2, 2];
            let data: [c_int; 8] = [-4, -3, -2, -1, 1, 2, 3, 4];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), LONG_IMG, 3, &naxes, &mut status);
            fits_write_3d_int(
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
            let mut result = [0 as c_int; 8];
            let mut anynull = -1;
            fits_read_3d_int(
                f.as_deref_mut().unwrap(),
                1,
                0,
                2,
                2,
                2,
                2,
                2,
                &mut result,
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[0], -4);
            assert_eq!(result[7], 4);
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
            let data: [c_int; 4] = [-100, -200, 100, 200];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), LONG_IMG, 2, &naxes, &mut status);
            fits_write_subset_int(
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
            let mut result = [0 as c_int; 16];
            let mut anynull = -1;
            fits_read_img_int(
                f.as_deref_mut().unwrap(),
                1,
                1,
                16,
                0,
                &mut result,
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[5], -100);
            assert_eq!(result[6], -200);
            assert_eq!(result[9], 100);
            assert_eq!(result[10], 200);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_group_parameters() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 1] = [4];
            let pdata: [c_int; 2] = [-10, 20];
            let idata: [c_int; 4] = [1, 2, 3, 4];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_grphdr(
                f.as_deref_mut().unwrap(),
                TRUE as c_int,
                LONG_IMG,
                1,
                &naxes,
                2,
                1,
                TRUE as c_int,
                &mut status,
            );
            fits_write_grppar_int(f.as_deref_mut().unwrap(), 1, 1, 2, &pdata, &mut status);
            fits_write_img_int(f.as_deref_mut().unwrap(), 1, 1, 4, &idata, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut presult = [0 as c_int; 2];
            fits_read_grppar_int(
                f.as_deref_mut().unwrap(),
                1,
                1,
                2,
                &mut presult,
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(presult[0], -10);
            assert_eq!(presult[1], 20);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_binary_table_column() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [c_int; 3] = [c_int::MIN, 0, c_int::MAX];

            let mut f = make_table(&name, BINARY_TBL, "KCOL", "1J", 3, &mut status);
            fits_write_col_int(f.as_deref_mut().unwrap(), 1, 1, 1, 3, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0 as c_int; 3];
            let mut anynull = -1;
            fits_read_col_int(
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
            assert_eq!(result[0], c_int::MIN);
            assert_eq!(result[1], 0);
            assert_eq!(result[2], c_int::MAX);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_column_with_null() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [c_int; 5] = [10, 20, -9999, 40, 50];

            let mut f = make_table(&name, BINARY_TBL, "KCOL", "5J", 1, &mut status);
            fits_set_btblnull(f.as_deref_mut().unwrap(), 1, -9999, &mut status);
            fits_write_colnull_int(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                5,
                &data,
                -9999,
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0 as c_int; 5];
            let mut anynull = -1;
            fits_read_col_int(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                5,
                0,
                &mut result,
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[0], 10);
            assert_eq!(result[1], 20);
            assert_eq!(result[3], 40);
            assert_eq!(result[4], 50);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_to_byte_column() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [c_int; 3] = [0, 128, 255];

            let mut f = make_table(&name, BINARY_TBL, "BCOL", "1B", 3, &mut status);
            fits_write_col_int(f.as_deref_mut().unwrap(), 1, 1, 1, 3, &data, &mut status);
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
            assert_eq!(result[1], 128);
            assert_eq!(result[2], 255);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_to_short_column() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [c_int; 3] = [i16::MIN as c_int, 0, i16::MAX as c_int];

            let mut f = make_table(&name, BINARY_TBL, "ICOL", "1I", 3, &mut status);
            fits_write_col_int(f.as_deref_mut().unwrap(), 1, 1, 1, 3, &data, &mut status);
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
            assert_eq!(result[0], i16::MIN);
            assert_eq!(result[1], 0);
            assert_eq!(result[2], i16::MAX);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_to_longlong_column() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [c_int; 3] = [c_int::MIN, 0, c_int::MAX];

            let mut f = make_table(&name, BINARY_TBL, "KCOL", "1K", 3, &mut status);
            fits_write_col_int(f.as_deref_mut().unwrap(), 1, 1, 1, 3, &data, &mut status);
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
            assert_eq!(result[0], c_int::MIN as LONGLONG);
            assert_eq!(result[1], 0);
            assert_eq!(result[2], c_int::MAX as LONGLONG);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_to_float_column() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [c_int; 3] = [-1000, 0, 1000];

            let mut f = make_table(&name, BINARY_TBL, "ECOL", "1E", 3, &mut status);
            fits_write_col_int(f.as_deref_mut().unwrap(), 1, 1, 1, 3, &data, &mut status);
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
            assert_eq!(result[0], -1000.0);
            assert_eq!(result[1], 0.0);
            assert_eq!(result[2], 1000.0);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_to_double_column() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [c_int; 3] = [-1000000, 0, 1000000];

            let mut f = make_table(&name, BINARY_TBL, "DCOL", "1D", 3, &mut status);
            fits_write_col_int(f.as_deref_mut().unwrap(), 1, 1, 1, 3, &data, &mut status);
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
            assert_eq!(result[0], -1000000.0);
            assert_eq!(result[1], 0.0);
            assert_eq!(result[2], 1000000.0);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_to_ascii_table() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [c_int; 3] = [-1000000, 0, 1000000];

            let mut f = make_table(&name, ASCII_TBL, "NUMCOL", "I15", 3, &mut status);
            fits_write_col_int(f.as_deref_mut().unwrap(), 1, 1, 1, 3, &data, &mut status);
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
    fn test_write_multirow() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [c_int; 9] = [-4, -3, -2, -1, 0, 1, 2, 3, 4];

            let mut f = make_table(&name, BINARY_TBL, "KCOL", "3J", 3, &mut status);
            fits_write_col_int(f.as_deref_mut().unwrap(), 1, 1, 1, 9, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0 as c_int; 9];
            let mut anynull = -1;
            fits_read_col_int(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                9,
                0,
                &mut result,
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[0], -4);
            assert_eq!(result[4], 0);
            assert_eq!(result[8], 4);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_variable_length_array() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data1: [c_int; 3] = [-30, -20, -10];
            let data2: [c_int; 2] = [10, 20];

            let mut f = make_table(&name, BINARY_TBL, "VARK", "1PJ", 2, &mut status);
            fits_write_col_int(f.as_deref_mut().unwrap(), 1, 1, 1, 3, &data1, &mut status);
            fits_write_col_int(f.as_deref_mut().unwrap(), 1, 2, 1, 2, &data2, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0 as c_int; 5];
            let mut anynull = -1;
            let mut nelem: c_long = 0;
            fits_read_descript(
                f.as_deref_mut().unwrap(),
                1,
                1,
                Some(&mut nelem),
                None,
                &mut status,
            );
            assert_eq!(nelem, 3);
            fits_read_col_int(
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
            assert_eq!(result[0], -30);
            assert_eq!(result[1], -20);
            assert_eq!(result[2], -10);
            fits_read_descript(
                f.as_deref_mut().unwrap(),
                1,
                2,
                Some(&mut nelem),
                None,
                &mut status,
            );
            assert_eq!(nelem, 2);
            fits_read_col_int(
                f.as_deref_mut().unwrap(),
                1,
                2,
                1,
                2,
                0,
                &mut result[3..],
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(result[3], 10);
            assert_eq!(result[4], 20);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_byte_overflow() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [c_int; 1] = [-1];

            let mut f = make_table(&name, BINARY_TBL, "BCOL", "1B", 1, &mut status);
            fits_write_col_int(f.as_deref_mut().unwrap(), 1, 1, 1, 1, &data, &mut status);
            assert_eq!(status, NUM_OVERFLOW);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_short_overflow() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [c_int; 1] = [i16::MIN as c_int - 1];

            let mut f = make_table(&name, BINARY_TBL, "ICOL", "1I", 1, &mut status);
            fits_write_col_int(f.as_deref_mut().unwrap(), 1, 1, 1, 1, &data, &mut status);
            assert_eq!(status, NUM_OVERFLOW);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_bad_status() {
        with_temp_file(|filename| {
            let mut status = 1;
            let name = to_buf(filename);
            let data: [c_int; 1] = [1];

            // Need a valid file pointer; functions return early on status > 0.
            let mut f: Option<Box<fitsfile>> = None;
            let mut st0 = 0;
            fits_create_file(&mut f, &name, &mut st0);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut st0);
            let fp = f.as_deref_mut().unwrap();

            fits_write_img_int(fp, 1, 1, 1, &data, &mut status);
            assert_eq!(status, 1);
            fits_write_imgnull_int(fp, 1, 1, 1, &data, 0, &mut status);
            assert_eq!(status, 1);
            fits_write_col_int(fp, 1, 1, 1, 1, &data, &mut status);
            assert_eq!(status, 1);
            fits_write_colnull_int(fp, 1, 1, 1, 1, &data, 0, &mut status);
            assert_eq!(status, 1);

            fits_close_file(f.take().unwrap(), &mut st0);
        });
    }

    #[test]
    fn test_write_null_column_no_tnull() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [c_int; 5] = [10, 20, -9999, 40, 50];

            let mut f = make_table(&name, BINARY_TBL, "KCOL", "5J", 1, &mut status);
            fits_write_colnull_int(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                5,
                &data,
                -9999,
                &mut status,
            );
            assert_eq!(status, NO_NULL);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_3d_noncontiguous() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 3] = [3, 3, 2];
            let data: [c_int; 32] = core::array::from_fn(|i| (i + 1) as c_int);

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), LONG_IMG, 3, &naxes, &mut status);
            // ncols=4, nrows=4 but naxis1=3, naxis2=3
            fits_write_3d_int(
                f.as_deref_mut().unwrap(),
                1,
                4,
                4,
                3,
                3,
                2,
                &data,
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut result = [0 as c_int; 18];
            let mut anynull = -1;
            fits_read_img_int(
                f.as_deref_mut().unwrap(),
                1,
                1,
                18,
                0,
                &mut result,
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[0], 1);
            assert_eq!(result[3], 5);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_3d_bad_dimen() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 3] = [5, 5, 2];
            let data: [c_int; 50] = [0; 50];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), LONG_IMG, 3, &naxes, &mut status);
            // ncols < naxis1 should fail
            fits_write_3d_int(
                f.as_deref_mut().unwrap(),
                1,
                3,
                3,
                5,
                5,
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
    fn test_subsection_bad_naxis() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 1] = [4];
            let fpixel: [c_long; 1] = [1];
            let lpixel: [c_long; 1] = [4];
            let data: [c_int; 4] = [1, 2, 3, 4];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), LONG_IMG, 1, &naxes, &mut status);
            // naxis = 0 should fail
            fits_write_subset_int(
                f.as_deref_mut().unwrap(),
                1,
                0,
                &naxes,
                &fpixel,
                &lpixel,
                &data,
                &mut status,
            );
            assert_eq!(status, BAD_DIMEN);
            status = 0;
            // naxis = 8 should fail
            fits_write_subset_int(
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
    fn test_write_with_scaling() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [c_int; 3] = [100, 200, 300];

            let mut f = make_table(&name, BINARY_TBL, "ICOL", "1I", 3, &mut status);
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
            fits_write_col_int(f.as_deref_mut().unwrap(), 1, 1, 1, 3, &data, &mut status);
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
            assert_eq!(result, [100.0, 200.0, 300.0]);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    /// Helper for scaled-write-then-read-as-double tests.
    fn scaled_col_roundtrip(tform: &str, tscal: Option<f64>, tzero: Option<f64>) {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [c_int; 3] = [100, 200, 300];

            let mut f = make_table(&name, BINARY_TBL, "COL", tform, 3, &mut status);
            if let Some(s) = tscal {
                fits_write_key_dbl(
                    f.as_deref_mut().unwrap(),
                    &cc("TSCAL1"),
                    s,
                    15,
                    None,
                    &mut status,
                );
            }
            if let Some(z) = tzero {
                fits_write_key_dbl(
                    f.as_deref_mut().unwrap(),
                    &cc("TZERO1"),
                    z,
                    15,
                    None,
                    &mut status,
                );
            }
            fits_write_col_int(f.as_deref_mut().unwrap(), 1, 1, 1, 3, &data, &mut status);
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
            assert_eq!(result, [100.0, 200.0, 300.0]);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_to_byte_column_scaled() {
        scaled_col_roundtrip("1B", Some(2.0), Some(50.0));
    }

    #[test]
    fn test_write_to_long_column_scaled() {
        scaled_col_roundtrip("1J", Some(2.0), Some(50.0));
    }

    #[test]
    fn test_write_to_longlong_column_scaled() {
        scaled_col_roundtrip("1K", Some(2.0), Some(50.0));
    }

    #[test]
    fn test_write_to_float_column_scaled() {
        scaled_col_roundtrip("1E", Some(2.0), Some(50.0));
    }

    #[test]
    fn test_write_to_double_column_scaled() {
        scaled_col_roundtrip("1D", Some(2.0), Some(50.0));
    }

    #[test]
    fn test_write_to_ascii_table_scaled() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [c_int; 3] = [100, 200, 300];

            let mut f = make_table(&name, ASCII_TBL, "NUMCOL", "F15.2", 3, &mut status);
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
                50.0,
                15,
                None,
                &mut status,
            );
            fits_write_col_int(f.as_deref_mut().unwrap(), 1, 1, 1, 3, &data, &mut status);
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
            assert_eq!(result, [100.0, 200.0, 300.0]);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    /// Helper for overflow-on-scaled-write tests.
    fn scaled_write_overflow(tform: &str, tscal: Option<f64>, tzero: Option<f64>, data: &[c_int]) {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);

            let mut f = make_table(
                &name,
                BINARY_TBL,
                "COL",
                tform,
                data.len() as LONGLONG,
                &mut status,
            );
            if let Some(s) = tscal {
                fits_write_key_dbl(
                    f.as_deref_mut().unwrap(),
                    &cc("TSCAL1"),
                    s,
                    15,
                    None,
                    &mut status,
                );
            }
            if let Some(z) = tzero {
                fits_write_key_dbl(
                    f.as_deref_mut().unwrap(),
                    &cc("TZERO1"),
                    z,
                    15,
                    None,
                    &mut status,
                );
            }
            fits_write_col_int(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                data.len() as LONGLONG,
                data,
                &mut status,
            );
            assert_eq!(status, NUM_OVERFLOW);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_byte_overflow_with_scaling() {
        // 200 / 0.5 = 400 > 255, should overflow
        scaled_write_overflow("1B", Some(0.5), None, &[50, 200, 100]);
    }

    #[test]
    fn test_byte_underflow_with_scaling() {
        // -100 - 50 = -150 < 0, should underflow
        scaled_write_overflow("1B", None, Some(50.0), &[-100, 100, 200]);
    }

    #[test]
    fn test_short_underflow_with_scaling() {
        // INT_MIN - 50000 < -32768, should underflow
        scaled_write_overflow("1I", None, Some(50000.0), &[c_int::MIN, 100, 200]);
    }

    #[test]
    fn test_short_overflow_with_scaling() {
        // 20000 / 0.5 = 40000 > 32767, should overflow
        scaled_write_overflow("1I", Some(0.5), None, &[1000, 20000, 5000]);
    }

    #[test]
    fn test_long_underflow_with_scaling() {
        // INT_MIN - 3e9 < -2147483648, should underflow
        scaled_write_overflow("1J", None, Some(3.0e9), &[c_int::MIN, 100, 200]);
    }

    #[test]
    fn test_long_overflow_with_scaling() {
        // INT_MAX / 0.5 > 2147483647, should overflow
        scaled_write_overflow("1J", Some(0.5), None, &[100, c_int::MAX, 200]);
    }

    #[test]
    fn test_longlong_underflow_with_scaling() {
        // INT_MIN - 1e19 < LONGLONG_MIN, should underflow
        scaled_write_overflow("1K", None, Some(1.0e19), &[c_int::MIN, 100, 200]);
    }

    #[test]
    fn test_longlong_overflow_with_scaling() {
        // INT_MAX / 1e-10 > LONGLONG_MAX, should overflow
        scaled_write_overflow("1K", Some(1.0e-10), None, &[100, c_int::MAX, 200]);
    }

    #[test]
    fn test_subsection_3d() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 3] = [4, 4, 2];
            let fpixel: [c_long; 3] = [2, 2, 1];
            let lpixel: [c_long; 3] = [3, 3, 2];
            let data: [c_int; 8] = core::array::from_fn(|i| (i + 100) as c_int);

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), LONG_IMG, 3, &naxes, &mut status);
            fits_write_subset_int(
                f.as_deref_mut().unwrap(),
                1,
                3,
                &naxes,
                &fpixel,
                &lpixel,
                &data,
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut result = [0 as c_int; 32];
            let mut anynull = -1;
            fits_read_img_int(
                f.as_deref_mut().unwrap(),
                1,
                1,
                32,
                0,
                &mut result,
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[5], 100);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_4d_subsection() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 4] = [3, 3, 2, 2];
            let fpixel: [c_long; 4] = [1, 1, 1, 1];
            let lpixel: [c_long; 4] = [2, 2, 2, 2];
            let data: [c_int; 16] = core::array::from_fn(|i| (i + 1) as c_int);

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), LONG_IMG, 4, &naxes, &mut status);
            fits_write_subset_int(
                f.as_deref_mut().unwrap(),
                1,
                4,
                &naxes,
                &fpixel,
                &lpixel,
                &data,
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut result = [0 as c_int; 36];
            let mut anynull = -1;
            fits_read_img_int(
                f.as_deref_mut().unwrap(),
                1,
                1,
                36,
                0,
                &mut result,
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[0], 1);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_5d_subsection() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 5] = [2, 2, 2, 2, 2];
            let fpixel: [c_long; 5] = [1, 1, 1, 1, 1];
            let lpixel: [c_long; 5] = [2, 2, 2, 2, 2];
            let data: [c_int; 32] = core::array::from_fn(|i| (i + 1) as c_int);

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), LONG_IMG, 5, &naxes, &mut status);
            fits_write_subset_int(
                f.as_deref_mut().unwrap(),
                1,
                5,
                &naxes,
                &fpixel,
                &lpixel,
                &data,
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut result = [0 as c_int; 32];
            let mut anynull = -1;
            fits_read_img_int(
                f.as_deref_mut().unwrap(),
                1,
                1,
                32,
                0,
                &mut result,
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[0], 1);
            assert_eq!(result[31], 32);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_nulls_at_end() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [c_int; 4] = [10, 20, -9999, -9999];

            let mut f = make_table(&name, BINARY_TBL, "VAL", "1J", 4, &mut status);
            fits_set_btblnull(f.as_deref_mut().unwrap(), 1, -9999, &mut status);
            fits_write_colnull_int(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                4,
                &data,
                -9999,
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0 as c_long; 4];
            let mut anynull = -1;
            fits_read_col_lng(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                4,
                0,
                &mut result,
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[0], 10);
            assert_eq!(result[1], 20);
            assert_eq!(result[2], -9999);
            assert_eq!(result[3], -9999);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_overflow_in_null_write() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [c_int; 3] = [100, 300, -999];

            let mut f = make_table(&name, BINARY_TBL, "BCOL", "1B", 3, &mut status);
            fits_set_btblnull(f.as_deref_mut().unwrap(), 1, 123, &mut status);
            // 300 overflows byte - tests overflow path in ffpcnk
            fits_write_colnull_int(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                3,
                &data,
                -999,
                &mut status,
            );
            assert_eq!(status, NUM_OVERFLOW);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_short_positive_overflow() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [c_int; 1] = [i16::MAX as c_int + 1];

            let mut f = make_table(&name, BINARY_TBL, "ICOL", "1I", 1, &mut status);
            fits_write_col_int(f.as_deref_mut().unwrap(), 1, 1, 1, 1, &data, &mut status);
            assert_eq!(status, NUM_OVERFLOW);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_negative_dvalue_short() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [c_int; 1] = [0];

            let mut f = make_table(&name, BINARY_TBL, "ICOL", "1I", 1, &mut status);
            // TSCAL=1, TZERO=100: stored = (0 - 100) / 1 = -100
            fits_write_key_dbl(
                f.as_deref_mut().unwrap(),
                &cc("TZERO1"),
                100.0,
                15,
                None,
                &mut status,
            );
            fits_write_col_int(f.as_deref_mut().unwrap(), 1, 1, 1, 1, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            // Read raw value without scaling
            fits_set_tscale(f.as_deref_mut().unwrap(), 1, 1.0, 0.0, &mut status);
            let mut result = [0i16; 1];
            let mut anynull = -1;
            fits_read_col_sht(
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
            assert_eq!(result[0], -100);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_negative_dvalue_int32() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [c_int; 1] = [0];

            let mut f = make_table(&name, BINARY_TBL, "JCOL", "1J", 1, &mut status);
            fits_write_key_dbl(
                f.as_deref_mut().unwrap(),
                &cc("TZERO1"),
                100.0,
                15,
                None,
                &mut status,
            );
            fits_write_col_int(f.as_deref_mut().unwrap(), 1, 1, 1, 1, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            fits_set_tscale(f.as_deref_mut().unwrap(), 1, 1.0, 0.0, &mut status);
            let mut result = [0 as c_long; 1];
            let mut anynull = -1;
            fits_read_col_lng(
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
            assert_eq!(result[0], -100);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_to_unsigned_longlong() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [c_int; 2] = [100, -1];

            let mut f = make_table(&name, BINARY_TBL, "KCOL", "1K", 2, &mut status);
            // Set TZERO to make this an unsigned longlong column
            fits_write_key_dbl(
                f.as_deref_mut().unwrap(),
                &cc("TZERO1"),
                9223372036854775808.0,
                20,
                None,
                &mut status,
            );
            // Negative value should cause overflow
            fits_write_col_int(f.as_deref_mut().unwrap(), 1, 1, 1, 2, &data, &mut status);
            assert_eq!(status, NUM_OVERFLOW);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_negative_dvalue_longlong() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [c_int; 1] = [0];

            let mut f = make_table(&name, BINARY_TBL, "KCOL", "1K", 1, &mut status);
            fits_write_key_dbl(
                f.as_deref_mut().unwrap(),
                &cc("TZERO1"),
                100.0,
                15,
                None,
                &mut status,
            );
            fits_write_col_int(f.as_deref_mut().unwrap(), 1, 1, 1, 1, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            fits_set_tscale(f.as_deref_mut().unwrap(), 1, 1.0, 0.0, &mut status);
            let mut result = [0 as LONGLONG; 1];
            let mut anynull = -1;
            fits_read_col_lnglng(
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
            assert_eq!(result[0], -100);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_ascii_table_overflow() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [c_int; 1] = [10000]; // 5 digits won't fit in 3 chars

            let mut f = make_table(&name, ASCII_TBL, "NUM", "I3", 1, &mut status);
            fits_write_col_int(f.as_deref_mut().unwrap(), 1, 1, 1, 1, &data, &mut status);
            assert_eq!(status, NUM_OVERFLOW);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_ascii_table_overflow_scaled() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [c_int; 1] = [1];

            let mut f = make_table(&name, ASCII_TBL, "NUM", "I3", 1, &mut status);
            // Scaling makes the stored value huge: (1-0)/0.0001 = 10000
            fits_write_key_dbl(
                f.as_deref_mut().unwrap(),
                &cc("TSCAL1"),
                0.0001,
                15,
                None,
                &mut status,
            );
            fits_write_col_int(f.as_deref_mut().unwrap(), 1, 1, 1, 1, &data, &mut status);
            assert_eq!(status, NUM_OVERFLOW);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_vla_with_nulls() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [c_int; 3] = [10, 20, 30];

            let mut f = make_table(&name, BINARY_TBL, "VLA", "1PJ", 1, &mut status);
            fits_set_btblnull(f.as_deref_mut().unwrap(), 1, -999, &mut status);
            // Write VLA with null value handling
            fits_write_colnull_int(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                3,
                &data,
                -999,
                &mut status,
            );
            assert_eq!(status, 0);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_to_ascii_string_column() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [c_int; 1] = [42];

            let mut f = make_table(&name, ASCII_TBL, "STR", "A10", 1, &mut status);
            // Try to write int to string column - should fail
            fits_write_col_int(f.as_deref_mut().unwrap(), 1, 1, 1, 1, &data, &mut status);
            assert_eq!(status, BAD_ATABLE_FORMAT);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_to_binary_string_column() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [c_int; 1] = [42];

            let mut f = make_table(&name, BINARY_TBL, "STR", "10A", 1, &mut status);
            // Try to write int to string column - should fail
            fits_write_col_int(f.as_deref_mut().unwrap(), 1, 1, 1, 1, &data, &mut status);
            assert_eq!(status, BAD_BTABLE_FORMAT);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }
}
