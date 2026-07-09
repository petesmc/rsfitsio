/*  This file, putcolb.rs, contains routines that write data elements to    */
/*  a FITS image or table with char (byte) datatype.                       */

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
use crate::putcolu::ffpclu_safe;
use crate::relibc::header::stdio::snprintf_f64;
use crate::wrappers::*;
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
pub unsafe extern "C" fn ffpprb(
    fptr: *mut fitsfile, /* I - FITS file pointer                       */
    group: c_long,       /* I - group to write(1 = 1st group)           */
    firstelem: LONGLONG, /* I - first vector element to write(1 = 1st)  */
    nelem: LONGLONG,     /* I - number of values to write               */
    array: *const u8,    /* I - array of values that are written   */
    status: *mut c_int,  /* IO - error status                           */
) -> c_int {
    unsafe {
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);

        let array = slice::from_raw_parts(array, nelem as usize);

        ffpprb_safe(fptr, group, firstelem, nelem, array, status)
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
pub fn ffpprb_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer                       */
    group: c_long,       /* I - group to write(1 = 1st group)           */
    firstelem: LONGLONG, /* I - first vector element to write(1 = 1st)  */
    nelem: LONGLONG,     /* I - number of values to write               */
    array: &[u8],        /* I - array of values that are written   */
    status: &mut c_int,  /* IO - error status                           */
) -> c_int {
    let nullvalue: u8 = 0;

    if fits_is_compressed_image_safe(fptr, status) > 0 {
        /* this is a compressed image in a binary table */

        fits_write_compressed_pixels(
            fptr,
            TBYTE,
            firstelem,
            nelem,
            NullCheckType::None,
            cast_slice(array),
            &Some(NullValue::UByte(nullvalue)),
            status,
        );
        return *status;
    }

    let row = cmp::max(1, group);

    ffpclb_safe(
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
pub unsafe extern "C" fn ffppnb(
    fptr: *mut fitsfile, /* I - FITS file pointer                       */
    group: c_long,       /* I - group to write(1 = 1st group)           */
    firstelem: LONGLONG, /* I - first vector element to write(1 = 1st)  */
    nelem: LONGLONG,     /* I - number of values to write               */
    array: *const u8,    /* I - array of values that are written   */
    nulval: u8,          /* I - undefined pixel value              */
    status: *mut c_int,  /* IO - error status                           */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        let array = slice::from_raw_parts(array, nelem as usize);

        ffppnb_safe(fptr, group, firstelem, nelem, array, nulval, status)
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
pub fn ffppnb_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer                       */
    group: c_long,       /* I - group to write(1 = 1st group)           */
    firstelem: LONGLONG, /* I - first vector element to write(1 = 1st)  */
    nelem: LONGLONG,     /* I - number of values to write               */
    array: &[u8],        /* I - array of values that are written   */
    nulval: u8,          /* I - undefined pixel value              */
    status: &mut c_int,  /* IO - error status                           */
) -> c_int {
    let mut nullvalue: u8 = 0;

    if fits_is_compressed_image_safe(fptr, status) > 0 {
        /* this is a compressed image in a binary table */

        nullvalue = nulval; /* set local variable */
        fits_write_compressed_pixels(
            fptr,
            TBYTE,
            firstelem,
            nelem,
            NullCheckType::SetPixel,
            cast_slice(array),
            &Some(NullValue::UByte(nullvalue)),
            status,
        );
        return *status;
    }

    let row = cmp::max(1, group);

    ffpcnb_safe(
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
pub unsafe extern "C" fn ffp2db(
    fptr: *mut fitsfile, /* I - FITS file pointer                     */
    group: c_long,       /* I - group to write(1 = 1st group)         */
    ncols: LONGLONG,     /* I - number of pixels in each row of array */
    naxis1: LONGLONG,    /* I - FITS image NAXIS1 value               */
    naxis2: LONGLONG,    /* I - FITS image NAXIS2 value               */
    array: *const u8,    /* I - array to be written               */
    status: *mut c_int,  /* IO - error status                         */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        let array = slice::from_raw_parts(array, (naxis1 * naxis2) as usize);

        ffp3db_safe(fptr, group, ncols, naxis2, naxis1, naxis2, 1, array, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Write an entire 2-D array of values to the primary array. Data conversion
/// and scaling will be performed if necessary (e.g, if the datatype of the
/// FITS array is not the same as the array being written).
pub fn ffp2db_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer                     */
    group: c_long,       /* I - group to write(1 = 1st group)         */
    ncols: LONGLONG,     /* I - number of pixels in each row of array */
    naxis1: LONGLONG,    /* I - FITS image NAXIS1 value               */
    naxis2: LONGLONG,    /* I - FITS image NAXIS2 value               */
    array: &[u8],        /* I - array to be written               */
    status: &mut c_int,  /* IO - error status                         */
) -> c_int {
    /* call the 3D writing routine, with the 3rd dimension = 1 */

    ffp3db_safe(fptr, group, ncols, naxis2, naxis1, naxis2, 1, array, status);

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
pub unsafe extern "C" fn ffp3db(
    fptr: *mut fitsfile, /* I - FITS file pointer                     */
    group: c_long,       /* I - group to write(1 = 1st group)         */
    ncols: LONGLONG,     /* I - number of pixels in each row of array */
    nrows: LONGLONG,     /* I - number of rows in each plane of array */
    naxis1: LONGLONG,    /* I - FITS image NAXIS1 value               */
    naxis2: LONGLONG,    /* I - FITS image NAXIS2 value               */
    naxis3: LONGLONG,    /* I - FITS image NAXIS3 value               */
    array: *const u8,    /* I - array to be written               */
    status: *mut c_int,  /* IO - error status                         */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        let array = slice::from_raw_parts(array, (ncols * naxis2 * naxis3) as usize);
        ffp3db_safe(
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
pub fn ffp3db_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer                     */
    group: c_long,       /* I - group to write(1 = 1st group)         */
    ncols: LONGLONG,     /* I - number of pixels in each row of array */
    nrows: LONGLONG,     /* I - number of rows in each plane of array */
    naxis1: LONGLONG,    /* I - FITS image NAXIS1 value               */
    naxis2: LONGLONG,    /* I - FITS image NAXIS2 value               */
    naxis3: LONGLONG,    /* I - FITS image NAXIS3 value               */
    array: &[u8],        /* I - array to be written               */
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
            TBYTE,
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
        ffpclb_safe(
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
            if ffpclb_safe(
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
pub unsafe extern "C" fn ffpssb(
    fptr: *mut fitsfile,   /* I - FITS file pointer                       */
    group: c_long,         /* I - group to write(1 = 1st group)           */
    naxis: c_long,         /* I - number of data axes in array            */
    naxes: *const c_long,  /* I - size of each FITS axis                  */
    fpixel: *const c_long, /* I - 1st pixel in each axis to write (1=1st) */
    lpixel: *const c_long, /* I - last pixel in each axis to write        */
    array: *const u8,      /* I - array to be written                 */
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

        ffpssb_safe(fptr, group, naxis, naxes, fpixel, lpixel, array, status)
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
pub fn ffpssb_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer                       */
    group: c_long,       /* I - group to write(1 = 1st group)           */
    naxis: c_long,       /* I - number of data axes in array            */
    naxes: &[c_long],    /* I - size of each FITS axis                  */
    fpixel: &[c_long],   /* I - 1st pixel in each axis to write (1=1st) */
    lpixel: &[c_long],   /* I - last pixel in each axis to write        */
    array: &[u8],        /* I - array to be written                 */
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
            TBYTE,
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
                            if ffpclb_safe(
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
pub unsafe extern "C" fn ffpgpb(
    fptr: *mut fitsfile, /* I - FITS file pointer                      */
    group: c_long,       /* I - group to write(1 = 1st group)          */
    firstelem: c_long,   /* I - first vector element to write(1 = 1st) */
    nelem: c_long,       /* I - number of values to write              */
    array: *const u8,    /* I - array of values that are written   */
    status: *mut c_int,  /* IO - error status                          */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        let array = slice::from_raw_parts(array, nelem as usize);

        ffpgpb_safe(fptr, group, firstelem, nelem, array, status)
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
pub fn ffpgpb_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer                      */
    group: c_long,       /* I - group to write(1 = 1st group)          */
    firstelem: c_long,   /* I - first vector element to write(1 = 1st) */
    nelem: c_long,       /* I - number of values to write              */
    array: &[u8],        /* I - array of values that are written   */
    status: &mut c_int,  /* IO - error status                          */
) -> c_int {
    let row = cmp::max(1, group);

    ffpclb_safe(
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
/// array.  FITSIO treats a primary array as a binary table with
/// 2 vector columns: the first column contains the group parameters (often
/// with length = 0) and the second column contains the array of image pixels.
/// Each row of the table represents a group in the case of multigroup FITS
/// images.
///
/// The input array of values will be converted to the datatype of the column
/// and will be inverse-scaled by the FITS TSCALn and TZEROn values if necessary.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffpclb(
    fptr: *mut fitsfile, /* I - FITS file pointer                       */
    colnum: c_int,       /* I - number of column to write (1 = 1st col) */
    firstrow: LONGLONG,  /* I - first row to write (1 = 1st row)        */
    firstelem: LONGLONG, /* I - first vector element to write (1 = 1st) */
    nelem: LONGLONG,     /* I - number of values to write               */
    array: *const u8,    /* I - array of values to write           */
    status: *mut c_int,  /* IO - error status                           */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        let array = slice::from_raw_parts(array, nelem as usize);

        ffpclb_safe(
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
/// array.  FITSIO treats a primary array as a binary table with
/// 2 vector columns: the first column contains the group parameters (often
/// with length = 0) and the second column contains the array of image pixels.
/// Each row of the table represents a group in the case of multigroup FITS
/// images.
///
/// The input array of values will be converted to the datatype of the column
/// and will be inverse-scaled by the FITS TSCALn and TZEROn values if necessary.
pub fn ffpclb_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer                       */
    colnum: c_int,       /* I - number of column to write (1 = 1st col) */
    firstrow: LONGLONG,  /* I - first row to write (1 = 1st row)        */
    firstelem: LONGLONG, /* I - first vector element to write (1 = 1st) */
    nelem: LONGLONG,     /* I - number of values to write               */
    array: &[u8],        /* I - array of values to write           */
    status: &mut c_int,  /* IO - error status                           */
) -> c_int {
    let mut writemode: c_int = 0;
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
    let mut wrtptr: LONGLONG = 0;
    let mut rowlen: LONGLONG = 0;
    let mut rownum: LONGLONG = 0;
    let mut remain: LONGLONG = 0;
    let mut next: LONGLONG = 0;
    let mut tnull: LONGLONG = 0;
    let mut maxelem: LONGLONG = 0;
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

    /* IMPORTANT NOTE: that the special case of using this subroutine
    to write bytes to a character column are handled internally
    by the call to ffgcprll() below.  It will adjust the effective
    *tcode, repeats, etc, to appear as a TBYTE column. */

    writemode = 17; /* Equivalent to writemode = 1 but allow TSTRING -> TBYTE */

    if ffgcprll(
        fptr,
        colnum,
        firstrow,
        firstelem,
        nelem,
        writemode,
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
      if there is no scaling
      then we can simply write the raw data bytes into the FITS file if the
      datatype of the FITS column is the same as the input values.  Otherwise,
      we must convert the raw values into the scaled and/or machine dependent
      format in a temporary buffer that has been allocated for this purpose.
    */
    if scale == 1.0 && zero == 0.0 && tcode == TBYTE {
        writeraw = true;
        if nelem < LONGLONG::from(INT32_MAX) {
            maxelem = nelem;
        } else {
            maxelem = LONGLONG::from(INT32_MAX);
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

        wrtptr = startpos + (rownum as LONGLONG * rowlen) + (elemnum * incre as LONGLONG);
        ffmbyt_safe(fptr, wrtptr, IGNORE_EOF, status); /* move to write position */

        match tcode {
            TBYTE => {
                if writeraw {
                    /* write raw input bytes without conversion */
                    ffpi1b(fptr, ntodo, incre, &array[(next as usize)..], status);
                } else {
                    /* convert the raw data before writing to FITS file */
                    ffi1fi1(
                        &array[(next as usize)..],
                        ntodo,
                        scale,
                        zero,
                        cast_slice_mut(&mut buffer),
                        status,
                    );
                    ffpi1b(fptr, ntodo, incre, cast_slice(&buffer), status);
                }
            }
            TLONGLONG => {
                ffi1fi8(
                    &array[(next as usize)..],
                    ntodo,
                    scale,
                    zero,
                    cast_slice_mut(&mut buffer),
                    status,
                );
                ffpi8b(fptr, ntodo, incre, cast_slice_mut(&mut buffer), status);
            }
            TSHORT => {
                ffi1fi2(
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
                ffi1fi4(
                    &array[(next as usize)..],
                    ntodo,
                    scale,
                    zero,
                    cast_slice_mut(&mut buffer),
                    status,
                );
                ffpi4b(fptr, ntodo, incre, cast_slice_mut(&mut buffer), status);
            }
            TFLOAT => {
                ffi1fr4(
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
                ffi1fr8(
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
                if strchr_safe(&tform, bb(b'A')).is_some() {
                    /* write raw input bytes without conversion        */
                    /* This case is a hack to let users write a stream */
                    /* of bytes directly to the 'A' format column      */

                    if incre == twidth {
                        ffpbyt(fptr, ntodo as LONGLONG, &array[(next as usize)..], status);
                    } else {
                        ffpbytoff(
                            fptr,
                            twidth,
                            ntodo / twidth,
                            incre - twidth,
                            &array[(next as usize)..],
                            status,
                        );
                    }
                } else if hdutype == ASCII_TBL && formlen > 1 {
                    if cform[formlen - 1] == bb(b'f') || cform[formlen - 1] == bb(b'E') {
                        ffi1fstr(
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
                "Error writing elements {:.0} thru {:.0} of input data array (ffpclb).",
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
/// Write an array of elements to the specified column of a table.  Any input
/// pixels equal to the value of nulvalue will be replaced by the appropriate
/// null value in the output FITS file.
///
/// The input array of values will be converted to the datatype of the column
/// and will be inverse-scaled by the FITS TSCALn and TZEROn values if necessary
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffpcnb(
    fptr: *mut fitsfile, /* I - FITS file pointer                       */
    colnum: c_int,       /* I - number of column to write (1 = 1st col) */
    firstrow: LONGLONG,  /* I - first row to write (1 = 1st row)        */
    firstelem: LONGLONG, /* I - first vector element to write (1 = 1st) */
    nelem: LONGLONG,     /* I - number of values to write               */
    array: *const u8,    /* I - array of values to write         */
    nulvalue: u8,        /* I - flag for undefined pixels        */
    status: *mut c_int,  /* IO - error status                           */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        let array = slice::from_raw_parts(array, nelem as usize);
        ffpcnb_safe(
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
pub fn ffpcnb_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer                       */
    colnum: c_int,       /* I - number of column to write (1 = 1st col) */
    firstrow: LONGLONG,  /* I - first row to write (1 = 1st row)        */
    firstelem: LONGLONG, /* I - first vector element to write (1 = 1st) */
    nelem: LONGLONG,     /* I - number of values to write               */
    array: &[u8],        /* I - array of values to write         */
    nulvalue: u8,        /* I - flag for undefined pixels        */
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

    /* if variable length array, first write the whole input vector,
    then go back and fill in the nulls */
    if tcode < 0
        && ffpclb_safe(
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
                    if ffpclb_safe(
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
            ffpclb_safe(
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
/// Write a stream of bytes to the current FITS HDU.  This primative routine is mainly
/// for writing non-standard "conforming" extensions and should not be used
/// for standard IMAGE, TABLE or BINTABLE extensions.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffpextn(
    fptr: *mut fitsfile,   /* I - FITS file pointer                        */
    offset: LONGLONG,      /* I - byte offset from start of extension data */
    nelem: LONGLONG,       /* I - number of elements to write              */
    buffer: *const c_void, /* I - stream of bytes to write                 */
    status: *mut c_int,    /* IO - error status                            */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        let buffer = slice::from_raw_parts(buffer.cast::<u8>(), nelem as usize);
        ffpextn_safe(fptr, offset, nelem, buffer, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Write a stream of bytes to the current FITS HDU.  This primative routine is mainly
/// for writing non-standard "conforming" extensions and should not be used
/// for standard IMAGE, TABLE or BINTABLE extensions.
pub fn ffpextn_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer                        */
    offset: LONGLONG,    /* I - byte offset from start of extension data */
    nelem: LONGLONG,     /* I - number of elements to write              */
    buffer: &[u8],       /* I - stream of bytes to write                 */
    status: &mut c_int,  /* IO - error status                            */
) -> c_int {
    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    /* reset position to the correct HDU if necessary */
    if fptr.HDUposition != fptr.Fptr.curhdu {
        ffmahd_safe(fptr, (fptr.HDUposition) + 1, None, status);

    /* rescan header if data structure is undefined */
    } else if fptr.Fptr.datastart == DATA_UNDEFINED as LONGLONG && ffrdef_safe(fptr, status) > 0 {
        return *status;
    }

    /* move to write position */
    ffmbyt_safe(fptr, fptr.Fptr.datastart + offset, IGNORE_EOF, status);

    /* write the buffer */
    ffpbyt(fptr, nelem, buffer, status);

    *status
}

/*--------------------------------------------------------------------------*/
/// Copy input to output prior to writing output to a FITS file.
/// Do datatype conversion and scaling if required
pub(crate) fn ffi1fi1(
    input: &[u8],       /* I - array of values to be converted  */
    ntodo: c_long,      /* I - number of elements in the array  */
    scale: f64,         /* I - FITS TSCALn or BSCALE value      */
    zero: f64,          /* I - FITS TZEROn or BZERO  value      */
    output: &mut [u8],  /* O - output array of converted values */
    status: &mut c_int, /* IO - error status                    */
) -> c_int {
    if scale == 1.0 && zero == 0.0 {
        output.copy_from_slice(input); /* just copy input to output */
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
/// Do datatype conversion and scaling if required
pub(crate) fn ffi1fi2(
    input: &[u8],           /* I - array of values to be converted  */
    ntodo: c_long,          /* I - number of elements in the array  */
    scale: f64,             /* I - FITS TSCALn or BSCALE value      */
    zero: f64,              /* I - FITS TZEROn or BZERO  value      */
    output: &mut [c_short], /* O - output array of converted values */
    status: &mut c_int,     /* IO - error status                    */
) -> c_int {
    if scale == 1.0 && zero == 0.0 {
        for ii in 0..(ntodo as usize) {
            output[ii] = c_short::from(input[ii]); /* just copy input to output */
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
pub(crate) fn ffi1fi4(
    input: &[u8],            /* I - array of values to be converted  */
    ntodo: c_long,           /* I - number of elements in the array  */
    scale: f64,              /* I - FITS TSCALn or BSCALE value      */
    zero: f64,               /* I - FITS TZEROn or BZERO  value      */
    output: &mut [INT32BIT], /* O - output array of converted values */
    status: &mut c_int,      /* IO - error status                    */
) -> c_int {
    if scale == 1.0 && zero == 0.0 {
        for ii in 0..(ntodo as usize) {
            output[ii] = INT32BIT::from(input[ii]); /* copy input to output */
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
/// Do datatype conversion and scaling if required
pub(crate) fn ffi1fi8(
    input: &[u8],            /* I - array of values to be converted  */
    ntodo: c_long,           /* I - number of elements in the array  */
    scale: f64,              /* I - FITS TSCALn or BSCALE value      */
    zero: f64,               /* I - FITS TZEROn or BZERO  value      */
    output: &mut [LONGLONG], /* O - output array of converted values */
    status: &mut c_int,      /* IO - error status                    */
) -> c_int {
    if scale == 1.0 && zero == 9223372036854775808. {
        /* Writing to unsigned long long column. */
        /* Instead of subtracting 9223372036854775808, it is more efficient */
        /* and more precise to just flip the sign bit with the XOR operator */

        /* no need to check range limits because all unsigned char values */
        /* are valid ULONGLONG values. */

        for ii in 0..(ntodo as usize) {
            output[ii] = (ULONGLONG::from(input[ii]) ^ 0x8000000000000000) as LONGLONG;
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
pub(crate) fn ffi1fr4(
    input: &[u8],       /* I - array of values to be converted  */
    ntodo: c_long,      /* I - number of elements in the array  */
    scale: f64,         /* I - FITS TSCALn or BSCALE value      */
    zero: f64,          /* I - FITS TZEROn or BZERO  value      */
    output: &mut [f32], /* O - output array of converted values */
    status: &mut c_int, /* IO - error status                    */
) -> c_int {
    if scale == 1.0 && zero == 0.0 {
        for ii in 0..(ntodo as usize) {
            output[ii] = f32::from(input[ii]);
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
pub(crate) fn ffi1fr8(
    input: &[u8],       /* I - array of values to be converted  */
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
pub(crate) fn ffi1fstr(
    input: &[u8],          /* I - array of values to be converted  */
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
            //unsafe { sprintf(output[ci..].as_mut_ptr(), cform.as_ptr(), input[ii] as f64) };
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
            // unsafe { sprintf(output[ci..].as_mut_ptr(), cform.as_ptr(), dvalue) };
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

// Ported from test_putcolb.c - unsigned char (byte) write functions.

#[cfg(test)]
mod tests {
    use super::*;
    use crate::aliases::rust_api::*;
    use crate::fitsio::{
        ASCII_TBL, BAD_BTABLE_FORMAT, BAD_COL_NUM, BAD_DIMEN, BINARY_TBL, BYTE_IMG, LONGLONG,
        NUM_OVERFLOW, READONLY, READWRITE, fitsfile,
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
            let data: [u8; 5] = [0, 50, 100, 200, 255];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 1, &naxes, &mut status);
            fits_write_img_byt(f.as_deref_mut().unwrap(), 1, 1, 5, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0, "setup");

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut result = [0u8; 5];
            let mut anynull = 0;
            fits_read_img_byt(
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
            assert_eq!(result, [0, 50, 100, 200, 255]);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_primary_with_null() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 1] = [5];
            let data: [u8; 5] = [10, 20, 30, 40, 50];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 1, &naxes, &mut status);
            fits_write_key_lng(
                f.as_deref_mut().unwrap(),
                &cc("BLANK"),
                30,
                None,
                &mut status,
            ); // Define null value first
            fits_write_imgnull_byt(f.as_deref_mut().unwrap(), 1, 1, 5, &data, 30, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0, "setup");

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut result = [0u8; 5];
            let mut nularray = [0 as c_char; 5];
            let mut anynull = 0;
            fits_read_imgnull_byt(
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
            assert_eq!(nularray[2], 1); // 30 = null
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
            let data: [u8; 6] = [1, 2, 3, 4, 5, 6];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 2, &naxes, &mut status);
            fits_write_2d_byt(f.as_deref_mut().unwrap(), 1, 3, 3, 2, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut result = [0u8; 6];
            let mut anynull = 0;
            fits_read_2d_byt(
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
            assert_eq!(result[0], 1);
            assert_eq!(result[5], 6);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_3d_array() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 3] = [2, 2, 2];
            let data: [u8; 8] = [1, 2, 3, 4, 5, 6, 7, 8];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 3, &naxes, &mut status);
            fits_write_3d_byt(
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
            let mut result = [0u8; 8];
            let mut anynull = 0;
            fits_read_3d_byt(
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
            assert_eq!(result[0], 1);
            assert_eq!(result[7], 8);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_subset() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 2] = [4, 4];
            let data: [u8; 4] = [99, 99, 99, 99];
            let fpixel: [c_long; 2] = [2, 2];
            let lpixel: [c_long; 2] = [3, 3];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 2, &naxes, &mut status);
            // Initialize with zeros
            let zeros = [0u8; 16];
            fits_write_img_byt(f.as_deref_mut().unwrap(), 1, 1, 16, &zeros, &mut status);
            // Write subset
            fits_write_subset_byt(
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
            let mut result = [0u8; 16];
            let mut anynull = 0;
            fits_read_img_byt(
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
            assert_eq!(result[0], 0);
            assert_eq!(result[5], 99); // (2,2)
            assert_eq!(result[6], 99); // (3,2)
            assert_eq!(result[9], 99); // (2,3)
            assert_eq!(result[10], 99); // (3,3)
            assert_eq!(result[15], 0);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_group_parameters() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 1] = [4];
            let pdata: [u8; 2] = [10, 20];
            let idata: [u8; 4] = [1, 2, 3, 4];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_grphdr(
                f.as_deref_mut().unwrap(),
                1,
                BYTE_IMG,
                1,
                &naxes,
                2,
                1,
                1,
                &mut status,
            );
            fits_write_grppar_byt(f.as_deref_mut().unwrap(), 1, 1, 2, &pdata, &mut status);
            fits_write_img_byt(f.as_deref_mut().unwrap(), 1, 1, 4, &idata, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut presult = [0u8; 2];
            fits_read_grppar_byt(
                f.as_deref_mut().unwrap(),
                1,
                1,
                2,
                &mut presult,
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(presult[0], 10);
            assert_eq!(presult[1], 20);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_column() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [u8; 3] = [0, 127, 255];

            let mut f = make_table(&name, BINARY_TBL, "BYTECOL", "1B", 3, &mut status);
            fits_write_col_byt(f.as_deref_mut().unwrap(), 1, 1, 1, 3, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0u8; 3];
            let mut anynull = 0;
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
    fn test_write_column_with_null() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [u8; 3] = [10, 20, 30];

            let mut f = make_table(&name, BINARY_TBL, "BYTECOL", "1B", 3, &mut status);
            fits_write_key_lng(
                f.as_deref_mut().unwrap(),
                &cc("TNULL1"),
                20,
                None,
                &mut status,
            );
            fits_write_colnull_byt(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                3,
                &data,
                20,
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0u8; 3];
            let mut nularray = [0 as c_char; 3];
            let mut anynull = 0;
            fits_read_colnull_byt(
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
            assert_eq!(result[0], 10);
            assert_eq!(nularray[1], 1); // 20 = null
            assert_eq!(result[2], 30);
            assert_eq!(anynull, 1);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_vector_column() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [u8; 5] = [1, 2, 3, 4, 5];

            let mut f = make_table(&name, BINARY_TBL, "BYTEVEC", "5B", 1, &mut status);
            fits_write_col_byt(f.as_deref_mut().unwrap(), 1, 1, 1, 5, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0u8; 5];
            let mut anynull = 0;
            fits_read_col_byt(
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
            assert_eq!(result[0], 1);
            assert_eq!(result[4], 5);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_to_short_column() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [u8; 3] = [0, 127, 255];

            let mut f = make_table(&name, BINARY_TBL, "SHORTCOL", "1I", 3, &mut status);
            fits_write_col_byt(f.as_deref_mut().unwrap(), 1, 1, 1, 3, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0i16; 3];
            let mut anynull = 0;
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
            assert_eq!(result[0], 0);
            assert_eq!(result[1], 127);
            assert_eq!(result[2], 255);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_to_long_column() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [u8; 3] = [0, 100, 255];

            let mut f = make_table(&name, BINARY_TBL, "LONGCOL", "1J", 3, &mut status);
            fits_write_col_byt(f.as_deref_mut().unwrap(), 1, 1, 1, 3, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0 as c_long; 3];
            let mut anynull = 0;
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
            assert_eq!(result[0], 0);
            assert_eq!(result[1], 100);
            assert_eq!(result[2], 255);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_to_float_column() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [u8; 3] = [0, 100, 255];

            let mut f = make_table(&name, BINARY_TBL, "FLOATCOL", "1E", 3, &mut status);
            fits_write_col_byt(f.as_deref_mut().unwrap(), 1, 1, 1, 3, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0f32; 3];
            let mut anynull = 0;
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
            assert_eq!(result[0], 0.0);
            assert_eq!(result[1], 100.0);
            assert_eq!(result[2], 255.0);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_to_double_column() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [u8; 3] = [0, 100, 255];

            let mut f = make_table(&name, BINARY_TBL, "DBLCOL", "1D", 3, &mut status);
            fits_write_col_byt(f.as_deref_mut().unwrap(), 1, 1, 1, 3, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0f64; 3];
            let mut anynull = 0;
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
            assert_eq!(result[0], 0.0);
            assert_eq!(result[1], 100.0);
            assert_eq!(result[2], 255.0);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_ascii_table() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [u8; 3] = [0, 100, 255];

            let mut f = make_table(&name, ASCII_TBL, "BYTECOL", "I4", 3, &mut status);
            fits_write_col_byt(f.as_deref_mut().unwrap(), 1, 1, 1, 3, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0u8; 3];
            let mut anynull = 0;
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
            assert_eq!(result[1], 100);
            assert_eq!(result[2], 255);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_with_scaling() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [u8; 1] = [110]; // Will be stored as (110 - 10) / 2 = 50

            let mut f = make_table(&name, BINARY_TBL, "SCALED", "1B", 1, &mut status);
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
            fits_write_col_byt(f.as_deref_mut().unwrap(), 1, 1, 1, 1, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0u8; 1];
            let mut anynull = 0;
            fits_read_col_byt(
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
            assert_eq!(result[0], 110);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_to_longlong_column() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [u8; 3] = [0, 128, 255];

            let mut f = make_table(&name, BINARY_TBL, "LLCOL", "1K", 3, &mut status);
            fits_write_col_byt(f.as_deref_mut().unwrap(), 1, 1, 1, 3, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0 as LONGLONG; 3];
            let mut anynull = 0;
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
            assert_eq!(result[0], 0);
            assert_eq!(result[1], 128);
            assert_eq!(result[2], 255);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_bad_col_num() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [u8; 1] = [1];

            let mut f = make_table(&name, BINARY_TBL, "BYTECOL", "1B", 1, &mut status);
            fits_write_col_byt(f.as_deref_mut().unwrap(), 0, 1, 1, 1, &data, &mut status);
            assert_eq!(status, BAD_COL_NUM);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_multirow() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [u8; 9] = [1, 2, 3, 4, 5, 6, 7, 8, 9];

            let mut f = make_table(&name, BINARY_TBL, "BYTEVEC", "3B", 3, &mut status);
            fits_write_col_byt(f.as_deref_mut().unwrap(), 1, 1, 1, 9, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0u8; 9];
            let mut anynull = 0;
            fits_read_col_byt(
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
            assert_eq!(result[0], 1);
            assert_eq!(result[4], 5);
            assert_eq!(result[8], 9);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_3d_noncontiguous() {
        // Test 3D non-contiguous write path
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 3] = [2, 2, 2];
            // Array is 3x3x2 but we only write 2x2x2
            let data: [u8; 18] = [1, 2, 0, 3, 4, 0, 0, 0, 0, 5, 6, 0, 7, 8, 0, 0, 0, 0];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 3, &naxes, &mut status);
            fits_write_3d_byt(
                f.as_deref_mut().unwrap(),
                1,
                3,
                3,
                2,
                2,
                2,
                &data,
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut result = [0u8; 8];
            let mut anynull = 0;
            fits_read_img_byt(
                f.as_deref_mut().unwrap(),
                1,
                1,
                8,
                0,
                &mut result,
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[0], 1);
            assert_eq!(result[1], 2);
            assert_eq!(result[2], 3);
            assert_eq!(result[3], 4);
            assert_eq!(result[4], 5);
            assert_eq!(result[5], 6);
            assert_eq!(result[6], 7);
            assert_eq!(result[7], 8);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_3d_bad_dimen() {
        // Test BAD_DIMEN error for 3D write with too small dimensions
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 3] = [4, 4, 2];
            let data: [u8; 8] = [1, 2, 3, 4, 5, 6, 7, 8];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 3, &naxes, &mut status);
            // ncols=2 < naxis1=4, should fail
            fits_write_3d_byt(
                f.as_deref_mut().unwrap(),
                1,
                2,
                2,
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
    fn test_subsection_bad_naxis() {
        // Test bad naxis in subset write
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 1] = [10];
            let fpix: [c_long; 1] = [1];
            let lpix: [c_long; 1] = [5];
            let data: [u8; 5] = [1, 2, 3, 4, 5];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 1, &naxes, &mut status);
            // naxis=0 should fail
            fits_write_subset_byt(
                f.as_deref_mut().unwrap(),
                1,
                0,
                &naxes,
                &fpix,
                &lpix,
                &data,
                &mut status,
            );
            assert_eq!(status, BAD_DIMEN);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_short_overflow_with_scaling() {
        // Test writing byte that overflows short when scaled
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            // With TSCAL=0.001, dvalue = input/0.001 = input*1000
            // 255 * 1000 = 255000 > 32767
            let data: [u8; 3] = [1, 255, 100];

            let mut f = make_table(&name, BINARY_TBL, "ICOL", "1I", 3, &mut status);
            fits_write_key_dbl(
                f.as_deref_mut().unwrap(),
                &cc("TSCAL1"),
                0.001,
                15,
                None,
                &mut status,
            );
            fits_write_col_byt(f.as_deref_mut().unwrap(), 1, 1, 1, 3, &data, &mut status);
            assert_eq!(status, NUM_OVERFLOW);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_short_underflow_with_scaling() {
        // Test writing byte that underflows short when scaled
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            // With TZERO=100000, dvalue = (input - 100000)/1 < -32768
            let data: [u8; 3] = [1, 2, 3];

            let mut f = make_table(&name, BINARY_TBL, "ICOL", "1I", 3, &mut status);
            fits_write_key_dbl(
                f.as_deref_mut().unwrap(),
                &cc("TZERO1"),
                100000.0,
                15,
                None,
                &mut status,
            );
            fits_write_col_byt(f.as_deref_mut().unwrap(), 1, 1, 1, 3, &data, &mut status);
            assert_eq!(status, NUM_OVERFLOW);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_int32_overflow_with_scaling() {
        // Test writing byte that overflows int32 when scaled
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            // With TSCAL=0.00000001, dvalue = input/0.00000001 = input*1e8
            // 255 * 1e8 > INT32_MAX
            let data: [u8; 3] = [1, 255, 100];

            let mut f = make_table(&name, BINARY_TBL, "JCOL", "1J", 3, &mut status);
            fits_write_key_dbl(
                f.as_deref_mut().unwrap(),
                &cc("TSCAL1"),
                0.00000001,
                15,
                None,
                &mut status,
            );
            fits_write_col_byt(f.as_deref_mut().unwrap(), 1, 1, 1, 3, &data, &mut status);
            assert_eq!(status, NUM_OVERFLOW);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_int32_underflow_with_scaling() {
        // Test writing byte that underflows int32 when scaled
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            // With TZERO=1e10, dvalue = (input - 1e10)/1 < INT32_MIN
            let data: [u8; 3] = [1, 2, 3];

            let mut f = make_table(&name, BINARY_TBL, "JCOL", "1J", 3, &mut status);
            fits_write_key_dbl(
                f.as_deref_mut().unwrap(),
                &cc("TZERO1"),
                1e10,
                15,
                None,
                &mut status,
            );
            fits_write_col_byt(f.as_deref_mut().unwrap(), 1, 1, 1, 3, &data, &mut status);
            assert_eq!(status, NUM_OVERFLOW);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_to_unsigned_longlong() {
        // Test writing byte to unsigned longlong column (TZERO = 2^63)
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [u8; 3] = [0, 128, 255];

            let mut f = make_table(&name, BINARY_TBL, "UKCOL", "1K", 3, &mut status);
            fits_write_key_dbl(
                f.as_deref_mut().unwrap(),
                &cc("TZERO1"),
                9223372036854775808.0,
                20,
                None,
                &mut status,
            );
            fits_write_col_byt(f.as_deref_mut().unwrap(), 1, 1, 1, 3, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            fits_set_tscale(f.as_deref_mut().unwrap(), 1, 1.0, 0.0, &mut status);
            let mut result = [0 as LONGLONG; 3];
            let mut anynull = 0;
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
    fn test_longlong_overflow_with_scaling() {
        // Test writing byte that overflows longlong when scaled
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            // With TSCAL=1e-19, dvalue = input/1e-19 = input*1e19 > LLONG_MAX
            let data: [u8; 3] = [1, 255, 100];

            let mut f = make_table(&name, BINARY_TBL, "KCOL", "1K", 3, &mut status);
            fits_write_key_dbl(
                f.as_deref_mut().unwrap(),
                &cc("TSCAL1"),
                1e-19,
                15,
                None,
                &mut status,
            );
            fits_write_col_byt(f.as_deref_mut().unwrap(), 1, 1, 1, 3, &data, &mut status);
            assert_eq!(status, NUM_OVERFLOW);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_longlong_underflow_with_scaling() {
        // Test writing byte that underflows longlong when scaled
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            // With TZERO=1e19, dvalue = (input - 1e19)/1 < LLONG_MIN
            let data: [u8; 3] = [1, 2, 3];

            let mut f = make_table(&name, BINARY_TBL, "KCOL", "1K", 3, &mut status);
            fits_write_key_dbl(
                f.as_deref_mut().unwrap(),
                &cc("TZERO1"),
                1e19,
                15,
                None,
                &mut status,
            );
            fits_write_col_byt(f.as_deref_mut().unwrap(), 1, 1, 1, 3, &data, &mut status);
            assert_eq!(status, NUM_OVERFLOW);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_to_float_column_scaled() {
        // Test writing byte to float column with scaling
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [u8; 3] = [100, 200, 255];

            let mut f = make_table(&name, BINARY_TBL, "ECOL", "1E", 3, &mut status);
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
            fits_write_col_byt(f.as_deref_mut().unwrap(), 1, 1, 1, 3, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0f32; 3];
            let mut anynull = 0;
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
            assert_eq!(result[0], 100.0);
            assert_eq!(result[1], 200.0);
            assert_eq!(result[2], 255.0);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_to_double_column_scaled() {
        // Test writing byte to double column with scaling
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [u8; 3] = [100, 200, 255];

            let mut f = make_table(&name, BINARY_TBL, "DCOL", "1D", 3, &mut status);
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
            fits_write_col_byt(f.as_deref_mut().unwrap(), 1, 1, 1, 3, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0f64; 3];
            let mut anynull = 0;
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
            assert_eq!(result[0], 100.0);
            assert_eq!(result[1], 200.0);
            assert_eq!(result[2], 255.0);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_ascii_table_scaled() {
        // Test writing byte to ASCII table with scaling
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [u8; 3] = [100, 200, 255];

            let mut f = make_table(&name, ASCII_TBL, "NUM", "F10.2", 3, &mut status);
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
            fits_write_col_byt(f.as_deref_mut().unwrap(), 1, 1, 1, 3, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0f64; 3];
            let mut anynull = 0;
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
            assert_eq!(result[0], 100.0);
            assert_eq!(result[1], 200.0);
            assert_eq!(result[2], 255.0);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_ascii_table_overflow_scaled() {
        // Test ASCII table overflow with scaling
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [u8; 1] = [1];

            // tform "I3" - only 3 chars width
            let mut f = make_table(&name, ASCII_TBL, "NUM", "I3", 1, &mut status);
            // With TSCAL=0.001, dvalue = 1/0.001 = 1000, needs 4 digits
            fits_write_key_dbl(
                f.as_deref_mut().unwrap(),
                &cc("TSCAL1"),
                0.001,
                15,
                None,
                &mut status,
            );
            fits_write_col_byt(f.as_deref_mut().unwrap(), 1, 1, 1, 1, &data, &mut status);
            assert_eq!(status, NUM_OVERFLOW);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_negative_dvalue_short() {
        // Test negative dvalue rounding for byte to short
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [u8; 1] = [0];

            let mut f = make_table(&name, BINARY_TBL, "ICOL", "1I", 1, &mut status);
            // TZERO=100 makes dvalue = (0-100)/1 = -100
            fits_write_key_dbl(
                f.as_deref_mut().unwrap(),
                &cc("TZERO1"),
                100.0,
                15,
                None,
                &mut status,
            );
            fits_write_col_byt(f.as_deref_mut().unwrap(), 1, 1, 1, 1, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            fits_set_tscale(f.as_deref_mut().unwrap(), 1, 1.0, 0.0, &mut status);
            let mut result = [0i16; 1];
            let mut anynull = 0;
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
    fn test_positive_dvalue_short() {
        // Test positive dvalue rounding for byte to short with scaling
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [u8; 1] = [100];

            let mut f = make_table(&name, BINARY_TBL, "ICOL", "1I", 1, &mut status);
            // TSCAL=0.5 makes dvalue = 100/0.5 = 200 (positive)
            fits_write_key_dbl(
                f.as_deref_mut().unwrap(),
                &cc("TSCAL1"),
                0.5,
                15,
                None,
                &mut status,
            );
            fits_write_col_byt(f.as_deref_mut().unwrap(), 1, 1, 1, 1, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            fits_set_tscale(f.as_deref_mut().unwrap(), 1, 1.0, 0.0, &mut status);
            let mut result = [0i16; 1];
            let mut anynull = 0;
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
            assert_eq!(result[0], 200);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_negative_dvalue_int32() {
        // Test negative dvalue rounding for byte to int32
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [u8; 1] = [0];

            let mut f = make_table(&name, BINARY_TBL, "JCOL", "1J", 1, &mut status);
            // TZERO=100 makes dvalue = (0-100)/1 = -100
            fits_write_key_dbl(
                f.as_deref_mut().unwrap(),
                &cc("TZERO1"),
                100.0,
                15,
                None,
                &mut status,
            );
            fits_write_col_byt(f.as_deref_mut().unwrap(), 1, 1, 1, 1, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            fits_set_tscale(f.as_deref_mut().unwrap(), 1, 1.0, 0.0, &mut status);
            let mut result = [0 as c_long; 1];
            let mut anynull = 0;
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
    fn test_positive_dvalue_int32() {
        // Test positive dvalue rounding for byte to int32 with scaling
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [u8; 1] = [100];

            let mut f = make_table(&name, BINARY_TBL, "JCOL", "1J", 1, &mut status);
            // TSCAL=0.5 makes dvalue = 100/0.5 = 200 (positive)
            fits_write_key_dbl(
                f.as_deref_mut().unwrap(),
                &cc("TSCAL1"),
                0.5,
                15,
                None,
                &mut status,
            );
            fits_write_col_byt(f.as_deref_mut().unwrap(), 1, 1, 1, 1, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            fits_set_tscale(f.as_deref_mut().unwrap(), 1, 1.0, 0.0, &mut status);
            let mut result = [0 as c_long; 1];
            let mut anynull = 0;
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
            assert_eq!(result[0], 200);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_negative_dvalue_longlong() {
        // Test negative dvalue rounding for byte to longlong
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [u8; 1] = [0];

            let mut f = make_table(&name, BINARY_TBL, "KCOL", "1K", 1, &mut status);
            // TZERO=100 makes dvalue = (0-100)/1 = -100
            fits_write_key_dbl(
                f.as_deref_mut().unwrap(),
                &cc("TZERO1"),
                100.0,
                15,
                None,
                &mut status,
            );
            fits_write_col_byt(f.as_deref_mut().unwrap(), 1, 1, 1, 1, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            fits_set_tscale(f.as_deref_mut().unwrap(), 1, 1.0, 0.0, &mut status);
            let mut result = [0 as LONGLONG; 1];
            let mut anynull = 0;
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
    fn test_positive_dvalue_longlong() {
        // Test positive dvalue rounding for byte to longlong with scaling
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [u8; 1] = [100];

            let mut f = make_table(&name, BINARY_TBL, "KCOL", "1K", 1, &mut status);
            // TSCAL=0.5 makes dvalue = 100/0.5 = 200 (positive, triggers +0.5 rounding)
            fits_write_key_dbl(
                f.as_deref_mut().unwrap(),
                &cc("TSCAL1"),
                0.5,
                15,
                None,
                &mut status,
            );
            fits_write_col_byt(f.as_deref_mut().unwrap(), 1, 1, 1, 1, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            fits_set_tscale(f.as_deref_mut().unwrap(), 1, 1.0, 0.0, &mut status);
            let mut result = [0 as LONGLONG; 1];
            let mut anynull = 0;
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
            assert_eq!(result[0], 200);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_vla_with_nulls() {
        // Test VLA write with nulls
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [u8; 3] = [10, 20, 30];

            let mut f = make_table(&name, BINARY_TBL, "VLA", "1PB", 1, &mut status);
            fits_set_btblnull(f.as_deref_mut().unwrap(), 1, 255, &mut status);
            fits_write_colnull_byt(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                3,
                &data,
                255,
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0);
        });
    }

    #[test]
    fn test_vla_with_overflow_null() {
        // Test VLA null handling where ffpclb returns overflow
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            // tform "1PI" - VLA of short
            // 40000 overflows short, triggers overflow in VLA path
            let data: [u8; 3] = [100, 200, 255];

            let mut f = make_table(&name, BINARY_TBL, "VLA", "1PI", 1, &mut status);
            fits_set_btblnull(f.as_deref_mut().unwrap(), 1, 255, &mut status);
            fits_write_key_dbl(
                f.as_deref_mut().unwrap(),
                &cc("TSCAL1"),
                0.001,
                15,
                None,
                &mut status,
            ); // 200/0.001 = 200000 > 32767
            fits_write_colnull_byt(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                3,
                &data,
                255,
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_vla_good_overflow() {
        // Test VLA write where writing good pixels causes overflow
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            // tform "1PI" - VLA of short
            // Pattern: null, overflow-good, null
            let data: [u8; 3] = [255, 200, 255];

            let mut f = make_table(&name, BINARY_TBL, "VLA", "1PI", 1, &mut status);
            fits_set_btblnull(f.as_deref_mut().unwrap(), 1, 255, &mut status);
            fits_write_key_dbl(
                f.as_deref_mut().unwrap(),
                &cc("TSCAL1"),
                0.001,
                15,
                None,
                &mut status,
            ); // 200/0.001 = 200000 > 32767
            fits_write_colnull_byt(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                3,
                &data,
                255,
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_byte_to_byte_overflow_with_scaling() {
        // Test writing byte to byte column with scaling that causes overflow
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            // With TSCAL=0.5, dvalue = input/0.5 = input*2
            // 255 * 2 = 510 > 255
            let data: [u8; 3] = [100, 255, 200];

            let mut f = make_table(&name, BINARY_TBL, "BCOL", "1B", 3, &mut status);
            fits_write_key_dbl(
                f.as_deref_mut().unwrap(),
                &cc("TSCAL1"),
                0.5,
                15,
                None,
                &mut status,
            );
            fits_write_col_byt(f.as_deref_mut().unwrap(), 1, 1, 1, 3, &data, &mut status);
            assert_eq!(status, NUM_OVERFLOW);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_byte_to_byte_underflow_with_scaling() {
        // Test writing byte to byte column with scaling that causes underflow
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            // With TZERO=1000, dvalue = (input-1000)/1 < 0
            let data: [u8; 3] = [0, 1, 2];

            let mut f = make_table(&name, BINARY_TBL, "BCOL", "1B", 3, &mut status);
            fits_write_key_dbl(
                f.as_deref_mut().unwrap(),
                &cc("TZERO1"),
                1000.0,
                15,
                None,
                &mut status,
            );
            fits_write_col_byt(f.as_deref_mut().unwrap(), 1, 1, 1, 3, &data, &mut status);
            assert_eq!(status, NUM_OVERFLOW);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_subsection_3d() {
        // Test 3D subsection write
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 3] = [10, 10, 5];
            let fpixel: [c_long; 3] = [2, 3, 1];
            let lpixel: [c_long; 3] = [4, 5, 2];
            let data: [u8; 18] = [
                1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18,
            ];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 3, &naxes, &mut status);
            fits_write_subset_byt(
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
            let mut result = [0u8; 3];
            let mut anynull = 0;
            fits_read_img_byt(
                f.as_deref_mut().unwrap(),
                1,
                22,
                3,
                0,
                &mut result,
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[0], 1);
            assert_eq!(result[1], 2);
            assert_eq!(result[2], 3);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_4d_subsection() {
        // Test 4D subsection write
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 4] = [4, 4, 4, 2];
            let fpixel: [c_long; 4] = [1, 1, 1, 1];
            let lpixel: [c_long; 4] = [2, 2, 2, 2];
            let data: [u8; 16] = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 4, &naxes, &mut status);
            fits_write_subset_byt(
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
            assert_eq!(status, 0);
        });
    }

    #[test]
    fn test_5d_subsection() {
        // Test 5D subsection write
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 5] = [3, 3, 3, 3, 2];
            let fpixel: [c_long; 5] = [1, 1, 1, 1, 1];
            let lpixel: [c_long; 5] = [2, 2, 2, 2, 2];
            let data: [u8; 32] = core::array::from_fn(|i| (i + 1) as u8);

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 5, &naxes, &mut status);
            fits_write_subset_byt(
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
            assert_eq!(status, 0);
        });
    }

    #[test]
    fn test_write_to_logical_column() {
        // Test BAD_BTABLE_FORMAT error - logical columns not supported
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [u8; 1] = [1];

            let mut f = make_table(&name, BINARY_TBL, "FLAG", "1L", 1, &mut status);
            fits_write_col_byt(f.as_deref_mut().unwrap(), 1, 1, 1, 1, &data, &mut status);
            assert_eq!(status, BAD_BTABLE_FORMAT);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_float_to_byte_unscaled() {
        // Test writing float data to unscaled byte column (scale=1, zero=0).
        // This hits the memcpy path in ffr4fi1.
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [f32; 3] = [10.0, 20.0, 100.0];

            let mut f = make_table(&name, BINARY_TBL, "BCOL", "1B", 3, &mut status);
            fits_write_col_flt(f.as_deref_mut().unwrap(), 1, 1, 1, 3, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0u8; 3];
            let mut anynull = 0;
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
            assert_eq!(result[0], 10);
            assert_eq!(result[1], 20);
            assert_eq!(result[2], 100);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_ascii_table_overflow() {
        // Test ASCII table overflow - large byte values that exceed field width.
        // This hits ffi1fstr overflow.
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            // tform "I2" - 2-character field - 255 needs 3 chars
            let data: [u8; 1] = [255];

            let mut f = make_table(&name, ASCII_TBL, "BCOL", "I2", 1, &mut status);
            fits_write_col_byt(f.as_deref_mut().unwrap(), 1, 1, 1, 1, &data, &mut status);
            assert_eq!(status, NUM_OVERFLOW);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_byte_to_string_column_raw() {
        // Test writing bytes to fixed-width A (string) column.
        // This hits the raw byte write path.
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [u8; 5] = [b'H', b'E', b'L', b'L', b'O'];

            let mut f = make_table(&name, BINARY_TBL, "STR", "5A", 1, &mut status);
            fits_write_col_byt(f.as_deref_mut().unwrap(), 1, 1, 1, 5, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0 as c_char; 6];
            let mut results: [&mut [c_char]; 1] = [&mut result];
            let mut anynull = 0;
            let empty = cc("");
            fits_read_col_str(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                1,
                Some(&empty),
                &mut results,
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(&result[..5], &cc("HELLO")[..5]);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_byte_to_ascii_a_column() {
        // Test writing bytes to ASCII table A (string) column.
        // This hits the TSTRING case with strchr(tform,'A').
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [u8; 5] = [b'H', b'E', b'L', b'L', b'O'];

            let mut f = make_table(&name, ASCII_TBL, "STR", "A5", 1, &mut status);
            fits_write_col_byt(f.as_deref_mut().unwrap(), 1, 1, 1, 5, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0 as c_char; 6];
            let mut results: [&mut [c_char]; 1] = [&mut result];
            let mut anynull = 0;
            let empty = cc("");
            fits_read_col_str(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                1,
                Some(&empty),
                &mut results,
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(&result[..5], &cc("HELLO")[..5]);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_vla_byte_overflow_cleared() {
        // Test VLA write where overflow is caught and cleared in null-handling path.
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            // tform "1PB" - VLA of byte
            // With TSCAL=0.5, 200/0.5 = 400 > 255, causes overflow
            let data: [u8; 3] = [50, 200, 100];

            let mut f = make_table(&name, BINARY_TBL, "VLA", "1PB", 1, &mut status);
            fits_set_btblnull(f.as_deref_mut().unwrap(), 1, 100, &mut status); // 100 is the null value
            fits_write_key_dbl(
                f.as_deref_mut().unwrap(),
                &cc("TSCAL1"),
                0.5,
                15,
                None,
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);

            // Reopen to pick up TSCAL
            fits_open_file(&mut f, &name, READWRITE, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            // This should overflow on 200 but continue for null handling
            fits_write_colnull_byt(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                3,
                &data,
                100,
                &mut status,
            );
            // Overflow is caught and cleared in VLA path, so final status is 0
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_fixed_array_null_overflow() {
        // Test fixed-length array null handling where writing good pixels overflows.
        // When good pixels precede null, the good pixel write can overflow.
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            // With TSCAL=0.5: 200/0.5=400 (overflow!), 100=null, 50/0.5=100 OK
            // Data order matters: overflow value BEFORE null triggers the path
            let data: [u8; 3] = [200, 100, 50];

            let mut f = make_table(&name, BINARY_TBL, "BCOL", "1B", 3, &mut status);
            fits_write_key_dbl(
                f.as_deref_mut().unwrap(),
                &cc("TSCAL1"),
                0.5,
                15,
                None,
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);

            // Reopen to pick up TSCAL
            fits_open_file(&mut f, &name, READWRITE, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            // Write with 100 as null - 200 before null will overflow on write
            fits_write_colnull_byt(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                3,
                &data,
                100,
                &mut status,
            );
            // Overflow should be caught and cleared
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_ffpextn() {
        // Test ffpextn - write raw bytes to extension data area.
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 1] = [10];
            let data: [u8; 10] = [0xAA, 0xBB, 0xCC, 0xDD, 0, 0, 0, 0, 0, 0];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 1, &naxes, &mut status);
            fits_write_img_byt(f.as_deref_mut().unwrap(), 1, 1, 10, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READWRITE, &mut status);
            fits_write_ext(f.as_deref_mut().unwrap(), 0, 4, &data, &mut status);
            assert_eq!(status, 0);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut result = [0u8; 4];
            let mut anynull = 0;
            fits_read_img_byt(
                f.as_deref_mut().unwrap(),
                1,
                1,
                4,
                0,
                &mut result,
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[0], 0xAA);
            assert_eq!(result[1], 0xBB);
            assert_eq!(result[2], 0xCC);
            assert_eq!(result[3], 0xDD);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }
}
