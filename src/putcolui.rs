/*  This file, putcolui.rs, contains routines that write data elements to   */
/*  a FITS image or table, with unsigned short datatype.                   */

/*  The FITSIO software was written by William Pence at the High Energy    */
/*  Astrophysic Science Archive Research Center (HEASARC) at the NASA      */
/*  Goddard Space Flight Center.                                           */

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
pub unsafe extern "C" fn ffpprui(
    fptr: *mut fitsfile,    /* I - FITS file pointer                       */
    group: c_long,          /* I - group to write (1 = 1st group)          */
    firstelem: LONGLONG,    /* I - first vector element to write (1 = 1st) */
    nelem: LONGLONG,        /* I - number of values to write               */
    array: *const c_ushort, /* I - array of values that are written        */
    status: *mut c_int,     /* IO - error status                           */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let array = slice::from_raw_parts(array, nelem as usize);

        ffpprui_safe(fptr, group, firstelem, nelem, array, status)
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
pub fn ffpprui_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer                       */
    group: c_long,       /* I - group to write (1 = 1st group)          */
    firstelem: LONGLONG, /* I - first vector element to write (1 = 1st) */
    nelem: LONGLONG,     /* I - number of values to write               */
    array: &[c_ushort],  /* I - array of values that are written        */
    status: &mut c_int,  /* IO - error status                           */
) -> c_int {
    let nullvalue: c_ushort = 0;

    if fits_is_compressed_image_safe(fptr, status) > 0 {
        /* this is a compressed image in a binary table */

        fits_write_compressed_pixels(
            fptr,
            TUSHORT,
            firstelem,
            nelem,
            NullCheckType::None,
            cast_slice(array),
            &Some(NullValue::UShort(nullvalue)),
            status,
        );
        return *status;
    }

    let row = cmp::max(1, group);

    ffpclui_safe(
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
pub unsafe extern "C" fn ffppnui(
    fptr: *mut fitsfile,    /* I - FITS file pointer                       */
    group: c_long,          /* I - group to write(1 = 1st group)           */
    firstelem: LONGLONG,    /* I - first vector element to write(1 = 1st)  */
    nelem: LONGLONG,        /* I - number of values to write               */
    array: *const c_ushort, /* I - array of values that are written        */
    nulval: c_ushort,       /* I - undefined pixel value                   */
    status: *mut c_int,     /* IO - error status                           */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        let array = slice::from_raw_parts(array, nelem as usize);

        ffppnui_safe(fptr, group, firstelem, nelem, array, nulval, status)
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
pub fn ffppnui_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer                       */
    group: c_long,       /* I - group to write(1 = 1st group)           */
    firstelem: LONGLONG, /* I - first vector element to write(1 = 1st)  */
    nelem: LONGLONG,     /* I - number of values to write               */
    array: &[c_ushort],  /* I - array of values that are written        */
    nulval: c_ushort,    /* I - undefined pixel value                   */
    status: &mut c_int,  /* IO - error status                           */
) -> c_int {
    let mut nullvalue: c_ushort = 0;

    if fits_is_compressed_image_safe(fptr, status) > 0 {
        /* this is a compressed image in a binary table */

        nullvalue = nulval; /* set local variable */
        fits_write_compressed_pixels(
            fptr,
            TUSHORT,
            firstelem,
            nelem,
            NullCheckType::SetPixel,
            cast_slice(array),
            &Some(NullValue::UShort(nullvalue)),
            status,
        );
        return *status;
    }

    let row = cmp::max(1, group);

    ffpcnui_safe(
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
pub unsafe extern "C" fn ffp2dui(
    fptr: *mut fitsfile,    /* I - FITS file pointer                     */
    group: c_long,          /* I - group to write(1 = 1st group)         */
    ncols: LONGLONG,        /* I - number of pixels in each row of array */
    naxis1: LONGLONG,       /* I - FITS image NAXIS1 value               */
    naxis2: LONGLONG,       /* I - FITS image NAXIS2 value               */
    array: *const c_ushort, /* I - array to be written                   */
    status: *mut c_int,     /* IO - error status                         */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        let array = slice::from_raw_parts(array, (ncols * naxis2) as usize);

        ffp2dui_safe(fptr, group, ncols, naxis1, naxis2, array, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Write an entire 2-D array of values to the primary array. Data conversion
/// and scaling will be performed if necessary (e.g, if the datatype of the
/// FITS array is not the same as the array being written).
pub fn ffp2dui_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer                     */
    group: c_long,       /* I - group to write(1 = 1st group)         */
    ncols: LONGLONG,     /* I - number of pixels in each row of array */
    naxis1: LONGLONG,    /* I - FITS image NAXIS1 value               */
    naxis2: LONGLONG,    /* I - FITS image NAXIS2 value               */
    array: &[c_ushort],  /* I - array to be written                   */
    status: &mut c_int,  /* IO - error status                         */
) -> c_int {
    /* call the 3D writing routine, with the 3rd dimension = 1 */
    ffp3dui_safe(fptr, group, ncols, naxis2, naxis1, naxis2, 1, array, status);

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
pub unsafe extern "C" fn ffp3dui(
    fptr: *mut fitsfile,    /* I - FITS file pointer                     */
    group: c_long,          /* I - group to write(1 = 1st group)         */
    ncols: LONGLONG,        /* I - number of pixels in each row of array */
    nrows: LONGLONG,        /* I - number of rows in each plane of array */
    naxis1: LONGLONG,       /* I - FITS image NAXIS1 value               */
    naxis2: LONGLONG,       /* I - FITS image NAXIS2 value               */
    naxis3: LONGLONG,       /* I - FITS image NAXIS3 value               */
    array: *const c_ushort, /* I - array to be written                   */
    status: *mut c_int,     /* IO - error status                         */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        let array = slice::from_raw_parts(array, (ncols * naxis2 * naxis3) as usize);

        ffp3dui_safe(
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
pub fn ffp3dui_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer                     */
    group: c_long,       /* I - group to write(1 = 1st group)         */
    ncols: LONGLONG,     /* I - number of pixels in each row of array */
    nrows: LONGLONG,     /* I - number of rows in each plane of array */
    naxis1: LONGLONG,    /* I - FITS image NAXIS1 value               */
    naxis2: LONGLONG,    /* I - FITS image NAXIS2 value               */
    naxis3: LONGLONG,    /* I - FITS image NAXIS3 value               */
    array: &[c_ushort],  /* I - array to be written                   */
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
            TUSHORT,
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

    /* arrays have same size? */
    if ncols == naxis1 && nrows == naxis2 {
        /* all the image pixels are contiguous, so write all at once */
        ffpclui_safe(
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
            if ffpclui_safe(
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
pub unsafe extern "C" fn ffpssui(
    fptr: *mut fitsfile,    /* I - FITS file pointer                       */
    group: c_long,          /* I - group to write(1 = 1st group)           */
    naxis: c_long,          /* I - number of data axes in array            */
    naxes: *const c_long,   /* I - size of each FITS axis                  */
    fpixel: *const c_long,  /* I - 1st pixel in each axis to write (1=1st) */
    lpixel: *const c_long,  /* I - last pixel in each axis to write        */
    array: *const c_ushort, /* I - array to be written                     */
    status: *mut c_int,     /* IO - error status                           */
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

        ffpssui_safe(fptr, group, naxis, naxes, fpixel, lpixel, array, status)
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
pub fn ffpssui_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer                       */
    group: c_long,       /* I - group to write(1 = 1st group)           */
    naxis: c_long,       /* I - number of data axes in array            */
    naxes: &[c_long],    /* I - size of each FITS axis                  */
    fpixel: &[c_long],   /* I - 1st pixel in each axis to write (1=1st) */
    lpixel: &[c_long],   /* I - last pixel in each axis to write        */
    array: &[c_ushort],  /* I - array to be written                     */
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
            TUSHORT,
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
                            if ffpclui_safe(
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
pub unsafe extern "C" fn ffpgpui(
    fptr: *mut fitsfile,    /* I - FITS file pointer                      */
    group: c_long,          /* I - group to write(1 = 1st group)          */
    firstelem: c_long,      /* I - first vector element to write(1 = 1st) */
    nelem: c_long,          /* I - number of values to write              */
    array: *const c_ushort, /* I - array of values that are written       */
    status: *mut c_int,     /* IO - error status                          */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        let array = slice::from_raw_parts(array, nelem as usize);

        ffpgpui_safe(fptr, group, firstelem, nelem, array, status)
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
pub fn ffpgpui_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer                      */
    group: c_long,       /* I - group to write(1 = 1st group)          */
    firstelem: c_long,   /* I - first vector element to write(1 = 1st) */
    nelem: c_long,       /* I - number of values to write              */
    array: &[c_ushort],  /* I - array of values that are written       */
    status: &mut c_int,  /* IO - error status                          */
) -> c_int {
    let row = cmp::max(1, group);

    ffpclui_safe(
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
pub unsafe extern "C" fn ffpclui(
    fptr: *mut fitsfile,    /* I - FITS file pointer                       */
    colnum: c_int,          /* I - number of column to write (1 = 1st col) */
    firstrow: LONGLONG,     /* I - first row to write (1 = 1st row)        */
    firstelem: LONGLONG,    /* I - first vector element to write (1 = 1st) */
    nelem: LONGLONG,        /* I - number of values to write               */
    array: *const c_ushort, /* I - array of values to write                */
    status: *mut c_int,     /* IO - error status                           */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        let array = slice::from_raw_parts(array, nelem as usize);

        ffpclui_safe(
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
pub fn ffpclui_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer                       */
    colnum: c_int,       /* I - number of column to write (1 = 1st col) */
    firstrow: LONGLONG,  /* I - first row to write (1 = 1st row)        */
    firstelem: LONGLONG, /* I - first vector element to write (1 = 1st) */
    nelem: LONGLONG,     /* I - number of values to write               */
    array: &[c_ushort],  /* I - array of values to write                */
    status: &mut c_int,  /* IO - error status                           */
) -> c_int {
    let mut tcode: c_int = 0;
    let mut maxelem: c_int = 0;
    let mut hdutype: c_int = 0;
    let mut twidth: c_long = 0;
    let mut incre: c_long = 0;
    let mut ntodo: c_long;
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
        &mut maxelem,
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

    if tcode == TSTRING {
        ffcfmt(&tform, &mut cform); /* derive C format for writing strings */
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
        ntodo = cmp::min(remain, LONGLONG::from(maxelem)) as c_long;
        ntodo = cmp::min(ntodo as LONGLONG, repeat - elemnum) as c_long;

        wrtptr = startpos + (rownum as LONGLONG * rowlen) + (elemnum * incre as LONGLONG);

        ffmbyt_safe(fptr, wrtptr, IGNORE_EOF, status); /* move to write position */

        match tcode {
            TSHORT => {
                /* convert the raw data before writing to FITS file */
                ffu2fi2(
                    &array[(next as usize)..],
                    ntodo,
                    scale,
                    zero,
                    cast_slice_mut(&mut buffer),
                    status,
                );
                ffpi2b(fptr, ntodo, incre, cast_slice_mut(&mut buffer), status);
            }
            TLONGLONG => {
                ffu2fi8(
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
                ffu2fi1(
                    &array[(next as usize)..],
                    ntodo,
                    scale,
                    zero,
                    cast_slice_mut(&mut buffer),
                    status,
                );
                ffpi1b(fptr, ntodo, incre, cast_slice_mut(&mut buffer), status);
            }
            TLONG => {
                ffu2fi4(
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
                ffu2fr4(
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
                ffu2fr8(
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
                    ffu2fstr(
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
                "Error writing elements {:.0} thru {:.0} of input data array (ffpclui).",
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
pub unsafe extern "C" fn ffpcnui(
    fptr: *mut fitsfile,    /* I - FITS file pointer                       */
    colnum: c_int,          /* I - number of column to write (1 = 1st col) */
    firstrow: LONGLONG,     /* I - first row to write (1 = 1st row)        */
    firstelem: LONGLONG,    /* I - first vector element to write (1 = 1st) */
    nelem: LONGLONG,        /* I - number of values to write               */
    array: *const c_ushort, /* I - array of values to write                */
    nulvalue: c_ushort,     /* I - value used to flag undefined pixels     */
    status: *mut c_int,     /* IO - error status                           */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        let array = slice::from_raw_parts(array, nelem as usize);

        ffpcnui_safe(
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
pub fn ffpcnui_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer                       */
    colnum: c_int,       /* I - number of column to write (1 = 1st col) */
    firstrow: LONGLONG,  /* I - first row to write (1 = 1st row)        */
    firstelem: LONGLONG, /* I - first vector element to write (1 = 1st) */
    nelem: LONGLONG,     /* I - number of values to write               */
    array: &[c_ushort],  /* I - array of values to write                */
    nulvalue: c_ushort,  /* I - value used to flag undefined pixels     */
    status: &mut c_int,  /* IO - error status                           */
) -> c_int {
    let mut ngood: LONGLONG = 0;
    let mut nbad: LONGLONG = 0;
    let mut repeat: LONGLONG = 0;
    let mut first: LONGLONG = 0;
    let mut fstelm: LONGLONG = 0;
    let mut fstrow: LONGLONG = 0;

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
        && ffpclui_safe(
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
    first = (firstrow - 1) * repeat + firstelem;

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
                    if ffpclui_safe(
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
            ffpclui_safe(
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
/// Do datatype conversion and scaling if required
pub(crate) fn ffu2fi1(
    input: &[c_ushort], /* I - array of values to be converted  */
    ntodo: c_long,      /* I - number of elements in the array  */
    scale: f64,         /* I - FITS TSCALn or BSCALE value      */
    zero: f64,          /* I - FITS TZEROn or BZERO  value      */
    output: &mut [u8],  /* O - output array of converted values */
    status: &mut c_int, /* IO - error status                    */
) -> c_int {
    if scale == 1.0 && zero == 0.0 {
        for ii in 0..(ntodo as usize) {
            if input[ii] > c_ushort::from(u8::MAX) {
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
/// Do datatype conversion and scaling if required
pub(crate) fn ffu2fi2(
    input: &[c_ushort],     /* I - array of values to be converted  */
    ntodo: c_long,          /* I - number of elements in the array  */
    scale: f64,             /* I - FITS TSCALn or BSCALE value      */
    zero: f64,              /* I - FITS TZEROn or BZERO  value      */
    output: &mut [c_short], /* O - output array of converted values */
    status: &mut c_int,     /* IO - error status                    */
                            /*

                            */
) -> c_int {
    if scale == 1.0 && zero == 32768. {
        /* Instead of subtracting 32768, it is more efficient */
        /* to just flip the sign bit with the XOR operator */

        for ii in 0..(ntodo as usize) {
            output[ii] = (input[ii] ^ 0x8000) as c_short;
        }
    } else if scale == 1.0 && zero == 0.0 {
        for ii in 0..(ntodo as usize) {
            if input[ii] > c_short::MAX as c_ushort {
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
pub(crate) fn ffu2fi4(
    input: &[c_ushort],      /* I - array of values to be converted  */
    ntodo: c_long,           /* I - number of elements in the array  */
    scale: f64,              /* I - FITS TSCALn or BSCALE value      */
    zero: f64,               /* I - FITS TZEROn or BZERO  value      */
    output: &mut [INT32BIT], /* O - output array of converted values */
    status: &mut c_int,      /* IO - error status                    */
                             /*
                             Copy input to output prior to writing output to a FITS file.
                             Do datatype conversion and scaling if required
                             */
) -> c_int {
    if scale == 1.0 && zero == 0.0 {
        for ii in 0..(ntodo as usize) {
            output[ii] = INT32BIT::from(input[ii]); /* just copy input to output */
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
pub(crate) fn ffu2fi8(
    input: &[c_ushort],      /* I - array of values to be converted  */
    ntodo: c_long,           /* I - number of elements in the array  */
    scale: f64,              /* I - FITS TSCALn or BSCALE value      */
    zero: f64,               /* I - FITS TZEROn or BZERO  value      */
    output: &mut [LONGLONG], /* O - output array of converted values */
    status: &mut c_int,      /* IO - error status                    */
                             /*
                             Copy input to output prior to writing output to a FITS file.
                             Do datatype conversion and scaling if required
                             */
) -> c_int {
    if scale == 1.0 && zero == 9223372036854775808. {
        /* Writing to unsigned long long column. Input values must not be negative */
        /* Instead of subtracting 9223372036854775808, it is more efficient */
        /* and more precise to just flip the sign bit with the XOR operator */

        /* no need to check range limits because all unsigned short values */
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
pub(crate) fn ffu2fr4(
    input: &[c_ushort], /* I - array of values to be converted  */
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
pub(crate) fn ffu2fr8(
    input: &[c_ushort], /* I - array of values to be converted  */
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
pub(crate) fn ffu2fstr(
    input: &[c_ushort],    /* I - array of values to be converted  */
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

#[cfg(test)]
mod tests {
    use crate::aliases::rust_api::*;
    use crate::fitsio::{
        ASCII_TBL, BAD_ATABLE_FORMAT, BAD_BTABLE_FORMAT, BAD_DIMEN, BINARY_TBL, BYTE_IMG, LONGLONG,
        NUM_OVERFLOW, READONLY, USHORT_IMG, fitsfile,
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
            let naxes: [c_long; 1] = [6];
            let data: [u16; 6] = [0, 100, 1000, 10000, 50000, 65535];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(
                f.as_deref_mut().unwrap(),
                USHORT_IMG,
                1,
                &naxes,
                &mut status,
            );
            fits_write_img_usht(f.as_deref_mut().unwrap(), 1, 1, 6, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0, "setup");

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut result = [0u16; 6];
            let mut anynull = 0;
            fits_read_img_usht(
                f.as_deref_mut().unwrap(),
                1,
                1,
                6,
                0,
                &mut result,
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result, [0, 100, 1000, 10000, 50000, 65535]);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_primary_with_nulls() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 1] = [5];
            let data: [u16; 5] = [100, 65535, 200, 65535, 300];
            let nulval: u16 = 65535;

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(
                f.as_deref_mut().unwrap(),
                USHORT_IMG,
                1,
                &naxes,
                &mut status,
            );
            fits_set_imgnull(f.as_deref_mut().unwrap(), nulval as LONGLONG, &mut status);
            fits_write_imgnull_usht(
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
            let mut result = [0u16; 5];
            let mut anynull = 0;
            fits_read_img_usht(
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
            assert_eq!(result[0], 100);
            assert_eq!(result[2], 200);
            assert_eq!(result[4], 300);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_2d_array() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 2] = [3, 2];
            let data: [u16; 6] = [1, 2, 3, 4, 5, 6];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(
                f.as_deref_mut().unwrap(),
                USHORT_IMG,
                2,
                &naxes,
                &mut status,
            );
            fits_write_2d_usht(f.as_deref_mut().unwrap(), 1, 3, 3, 2, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut result = [0u16; 6];
            let mut anynull = 0;
            fits_read_img_usht(
                f.as_deref_mut().unwrap(),
                1,
                1,
                6,
                0,
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
            let data: [u16; 8] = [1, 2, 3, 4, 5, 6, 7, 8];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(
                f.as_deref_mut().unwrap(),
                USHORT_IMG,
                3,
                &naxes,
                &mut status,
            );
            fits_write_3d_usht(
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
            let mut result = [0u16; 8];
            let mut anynull = 0;
            fits_read_img_usht(
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
            assert_eq!(result[7], 8);
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
            let data: [u16; 4] = [10, 20, 30, 40];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(
                f.as_deref_mut().unwrap(),
                USHORT_IMG,
                2,
                &naxes,
                &mut status,
            );
            fits_write_subset_usht(
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
            let mut result = [0u16; 16];
            let mut anynull = 0;
            fits_read_img_usht(
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
            // Check that the 2x2 subsection was written at positions [1][1] to [2][2]
            assert_eq!(result[5], 10); // [1][1] = index 4+1 = 5
            assert_eq!(result[6], 20); // [2][1] = index 4+2 = 6
            assert_eq!(result[9], 30); // [1][2] = index 8+1 = 9
            assert_eq!(result[10], 40); // [2][2] = index 8+2 = 10
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_binary_table_column() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [u16; 3] = [0, 32768, 65535];

            let mut f = make_table(&name, BINARY_TBL, "UICOL", "1U", 3, &mut status);
            fits_write_col_usht(f.as_deref_mut().unwrap(), 1, 1, 1, 3, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0u16; 3];
            let mut anynull = 0;
            fits_read_col_usht(
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
            assert_eq!(result[1], 32768);
            assert_eq!(result[2], 65535);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_column_with_nulls() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [u16; 5] = [100, 0, 200, 0, 300];
            let nulval: u16 = 0;

            let mut f = make_table(&name, BINARY_TBL, "UICOL", "5U", 1, &mut status);
            fits_set_btblnull(f.as_deref_mut().unwrap(), 1, 1, &mut status);
            fits_write_colnull_usht(
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
            let mut result = [0u16; 5];
            let mut anynull = 0;
            fits_read_col_usht(
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
            assert_eq!(result[0], 100);
            assert_eq!(result[2], 200);
            assert_eq!(result[4], 300);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_to_short_column() {
        // Write unsigned short to signed short column with TZERO=32768.
        // This exercises the XOR path in ffu2fi2.
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [u16; 3] = [0, 32768, 65535];

            let mut f = make_table(&name, BINARY_TBL, "ICOL", "1I", 3, &mut status);
            // Set TZERO=32768 to store unsigned values in signed column
            fits_set_tscale(f.as_deref_mut().unwrap(), 1, 1.0, 32768.0, &mut status);
            fits_write_col_usht(f.as_deref_mut().unwrap(), 1, 1, 1, 3, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0);
        });
    }

    #[test]
    fn test_write_to_long_column() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [u16; 3] = [0, 32768, 65535];

            let mut f = make_table(&name, BINARY_TBL, "JCOL", "1J", 3, &mut status);
            fits_write_col_usht(f.as_deref_mut().unwrap(), 1, 1, 1, 3, &data, &mut status);
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
            assert_eq!(result[1], 32768);
            assert_eq!(result[2], 65535);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_to_float_column() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [u16; 3] = [0, 32768, 65535];

            let mut f = make_table(&name, BINARY_TBL, "ECOL", "1E", 3, &mut status);
            fits_write_col_usht(f.as_deref_mut().unwrap(), 1, 1, 1, 3, &data, &mut status);
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
            assert_eq!(result[1], 32768.0);
            assert_eq!(result[2], 65535.0);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_to_double_column() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [u16; 3] = [0, 32768, 65535];

            let mut f = make_table(&name, BINARY_TBL, "DCOL", "1D", 3, &mut status);
            fits_write_col_usht(f.as_deref_mut().unwrap(), 1, 1, 1, 3, &data, &mut status);
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
            assert_eq!(result[1], 32768.0);
            assert_eq!(result[2], 65535.0);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_to_longlong_column() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [u16; 3] = [0, 32768, 65535];

            let mut f = make_table(&name, BINARY_TBL, "KCOL", "1K", 3, &mut status);
            fits_write_col_usht(f.as_deref_mut().unwrap(), 1, 1, 1, 3, &data, &mut status);
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
            assert_eq!(result[1], 32768);
            assert_eq!(result[2], 65535);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_to_byte_column() {
        // Write unsigned short to byte column - values <= 255 fit.
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [u16; 3] = [0, 127, 255];

            let mut f = make_table(&name, BINARY_TBL, "BCOL", "1B", 3, &mut status);
            fits_write_col_usht(f.as_deref_mut().unwrap(), 1, 1, 1, 3, &data, &mut status);
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
    fn test_byte_overflow() {
        // Write unsigned short > 255 to byte column - should overflow.
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [u16; 1] = [256];

            let mut f = make_table(&name, BINARY_TBL, "BCOL", "1B", 1, &mut status);
            fits_write_col_usht(f.as_deref_mut().unwrap(), 1, 1, 1, 1, &data, &mut status);
            assert_eq!(status, NUM_OVERFLOW);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_to_ascii_table() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [u16; 3] = [0, 32768, 65535];

            let mut f = make_table(&name, ASCII_TBL, "NUMCOL", "F10.1", 3, &mut status);
            fits_write_col_usht(f.as_deref_mut().unwrap(), 1, 1, 1, 3, &data, &mut status);
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
            assert_eq!(result[1], 32768.0);
            assert_eq!(result[2], 65535.0);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_group_parameters() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 1] = [4];
            let pdata: [u16; 2] = [100, 200];
            let idata: [u16; 4] = [1, 2, 3, 4];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_grphdr(
                f.as_deref_mut().unwrap(),
                1,
                USHORT_IMG,
                1,
                &naxes,
                2,
                1,
                1,
                &mut status,
            );
            fits_write_grppar_usht(f.as_deref_mut().unwrap(), 1, 1, 2, &pdata, &mut status);
            fits_write_img_usht(f.as_deref_mut().unwrap(), 1, 1, 4, &idata, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut presult = [0f32; 2];
            fits_read_grppar_flt(
                f.as_deref_mut().unwrap(),
                1,
                1,
                2,
                &mut presult,
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(presult[0], 100.0);
            assert_eq!(presult[1], 200.0);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_with_scaling() {
        // Test scaled write path - verify write succeeds.
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [u16; 3] = [100, 200, 300];

            let mut f = make_table(&name, BINARY_TBL, "ICOL", "1I", 3, &mut status);
            // scale=2.0, zero=50 => stored = (input - 50) / 2
            fits_set_tscale(f.as_deref_mut().unwrap(), 1, 2.0, 50.0, &mut status);
            fits_write_col_usht(f.as_deref_mut().unwrap(), 1, 1, 1, 3, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0);
        });
    }

    #[test]
    fn test_write_noncontiguous_3d() {
        // Test 3D write where array dimensions differ from FITS dimensions.
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 3] = [2, 2, 2];
            // Input array is 3x3x2 but we only write 2x2x2
            let data: [u16; 18] = [
                1, 2, 99, 3, 4, 99, 99, 99, 99, 5, 6, 99, 7, 8, 99, 99, 99, 99,
            ];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(
                f.as_deref_mut().unwrap(),
                USHORT_IMG,
                3,
                &naxes,
                &mut status,
            );
            fits_write_3d_usht(
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
            let mut result = [0u16; 8];
            let mut anynull = 0;
            fits_read_img_usht(
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
    fn test_write_multirow() {
        // Write data that spans multiple rows.
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [u16; 9] = [1, 2, 3, 4, 5, 6, 7, 8, 9];

            let mut f = make_table(&name, BINARY_TBL, "UICOL", "3U", 3, &mut status);
            fits_write_col_usht(f.as_deref_mut().unwrap(), 1, 1, 1, 9, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0u16; 9];
            let mut anynull = 0;
            fits_read_col_usht(
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
    fn test_write_variable_length_array() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data1: [u16; 3] = [100, 200, 300];
            let data2: [u16; 2] = [400, 500];

            let mut f = make_table(&name, BINARY_TBL, "VARUI", "1PU", 2, &mut status);
            fits_write_col_usht(f.as_deref_mut().unwrap(), 1, 1, 1, 3, &data1, &mut status);
            fits_write_col_usht(f.as_deref_mut().unwrap(), 1, 2, 1, 2, &data2, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0u16; 5];
            let mut anynull = 0;
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
            fits_read_col_usht(
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
            assert_eq!(result[0], 100);
            assert_eq!(result[1], 200);
            assert_eq!(result[2], 300);
            fits_read_descript(
                f.as_deref_mut().unwrap(),
                1,
                2,
                Some(&mut nelem),
                None,
                &mut status,
            );
            assert_eq!(nelem, 2);
            fits_read_col_usht(
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
            assert_eq!(result[3], 400);
            assert_eq!(result[4], 500);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_short_overflow() {
        // Write to signed short without TZERO - values > SHRT_MAX overflow.
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [u16; 1] = [32768]; // > SHRT_MAX

            let mut f = make_table(&name, BINARY_TBL, "ICOL", "1I", 1, &mut status);
            // No TZERO, so 32768 > SHRT_MAX causes overflow
            fits_write_col_usht(f.as_deref_mut().unwrap(), 1, 1, 1, 1, &data, &mut status);
            assert_eq!(status, NUM_OVERFLOW);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_to_unsigned_longlong_column() {
        // Write to unsigned longlong column with TZERO=2^63.
        // This exercises the XOR path in ffu2fi8.
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [u16; 3] = [0, 32768, 65535];

            let mut f = make_table(&name, BINARY_TBL, "ULLCOL", "1K", 3, &mut status);
            fits_set_tscale(
                f.as_deref_mut().unwrap(),
                1,
                1.0,
                9223372036854775808.0,
                &mut status,
            );
            fits_write_col_usht(f.as_deref_mut().unwrap(), 1, 1, 1, 3, &data, &mut status);
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
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_scaled_byte_underflow() {
        // With scale and zero, value becomes < 0.
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            // With scale=1 and zero=100, stored = (0 - 100) / 1 = -100 < 0
            let data: [u16; 1] = [0];

            let mut f = make_table(&name, BINARY_TBL, "BCOL", "1B", 1, &mut status);
            fits_set_tscale(f.as_deref_mut().unwrap(), 1, 1.0, 100.0, &mut status);
            fits_write_col_usht(f.as_deref_mut().unwrap(), 1, 1, 1, 1, &data, &mut status);
            assert_eq!(status, NUM_OVERFLOW);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_scaled_byte_overflow() {
        // With scale and zero, value becomes > 255.
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            // With scale=0.1 and zero=0, stored = 65535/0.1 = 655350 > 255
            let data: [u16; 1] = [65535];

            let mut f = make_table(&name, BINARY_TBL, "BCOL", "1B", 1, &mut status);
            fits_set_tscale(f.as_deref_mut().unwrap(), 1, 0.1, 0.0, &mut status);
            fits_write_col_usht(f.as_deref_mut().unwrap(), 1, 1, 1, 1, &data, &mut status);
            assert_eq!(status, NUM_OVERFLOW);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_scaled_short_underflow() {
        // Scaled value underflows signed short.
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            // With scale=1 and zero=100000, stored = (0 - 100000) = -100000 < SHRT_MIN
            let data: [u16; 1] = [0];

            let mut f = make_table(&name, BINARY_TBL, "ICOL", "1I", 1, &mut status);
            fits_set_tscale(f.as_deref_mut().unwrap(), 1, 1.0, 100000.0, &mut status);
            fits_write_col_usht(f.as_deref_mut().unwrap(), 1, 1, 1, 1, &data, &mut status);
            assert_eq!(status, NUM_OVERFLOW);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_scaled_short_overflow() {
        // Scaled value overflows signed short.
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            // With scale=0.5 and zero=0, stored = 65535/0.5 = 131070 > SHRT_MAX
            let data: [u16; 1] = [65535];

            let mut f = make_table(&name, BINARY_TBL, "ICOL", "1I", 1, &mut status);
            fits_set_tscale(f.as_deref_mut().unwrap(), 1, 0.5, 0.0, &mut status);
            fits_write_col_usht(f.as_deref_mut().unwrap(), 1, 1, 1, 1, &data, &mut status);
            assert_eq!(status, NUM_OVERFLOW);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_scaled_long_underflow() {
        // Scaled value underflows INT32.
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            // With scale=1 and zero=3000000000, stored = (0 - 3e9) < INT_MIN
            let data: [u16; 1] = [0];

            let mut f = make_table(&name, BINARY_TBL, "JCOL", "1J", 1, &mut status);
            fits_set_tscale(f.as_deref_mut().unwrap(), 1, 1.0, 3000000000.0, &mut status);
            fits_write_col_usht(f.as_deref_mut().unwrap(), 1, 1, 1, 1, &data, &mut status);
            assert_eq!(status, NUM_OVERFLOW);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_scaled_long_overflow() {
        // Scaled value overflows INT32.
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            // With scale=0.00001 and zero=0, stored = 65535/0.00001 = 6.5e9 > INT_MAX
            let data: [u16; 1] = [65535];

            let mut f = make_table(&name, BINARY_TBL, "JCOL", "1J", 1, &mut status);
            fits_set_tscale(f.as_deref_mut().unwrap(), 1, 0.00001, 0.0, &mut status);
            fits_write_col_usht(f.as_deref_mut().unwrap(), 1, 1, 1, 1, &data, &mut status);
            assert_eq!(status, NUM_OVERFLOW);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_scaled_longlong_underflow() {
        // Scaled value underflows LONGLONG.
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            // With scale=1e-17 and zero=1e19, stored underflows
            let data: [u16; 1] = [0];

            let mut f = make_table(&name, BINARY_TBL, "KCOL", "1K", 1, &mut status);
            fits_set_tscale(f.as_deref_mut().unwrap(), 1, 1e-17, 1e19, &mut status);
            fits_write_col_usht(f.as_deref_mut().unwrap(), 1, 1, 1, 1, &data, &mut status);
            assert_eq!(status, NUM_OVERFLOW);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_scaled_longlong_overflow() {
        // Scaled value overflows LONGLONG.
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            // With scale=1e-17 and zero=0, stored = 65535/1e-17 = 6.5e21 > LONGLONG_MAX
            let data: [u16; 1] = [65535];

            let mut f = make_table(&name, BINARY_TBL, "KCOL", "1K", 1, &mut status);
            fits_set_tscale(f.as_deref_mut().unwrap(), 1, 1e-17, 0.0, &mut status);
            fits_write_col_usht(f.as_deref_mut().unwrap(), 1, 1, 1, 1, &data, &mut status);
            assert_eq!(status, NUM_OVERFLOW);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_scaled_ascii_overflow() {
        // Scaled ASCII value overflows column width.
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            // With scale=0.001, 65535/0.001 = 65535000 needs 10 chars
            let data: [u16; 1] = [65535];

            let mut f = make_table(&name, ASCII_TBL, "NUMCOL", "F6.1", 1, &mut status);
            fits_set_tscale(f.as_deref_mut().unwrap(), 1, 0.001, 0.0, &mut status);
            fits_write_col_usht(f.as_deref_mut().unwrap(), 1, 1, 1, 1, &data, &mut status);
            assert_eq!(status, NUM_OVERFLOW);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_bad_dimen_3d() {
        // Test BAD_DIMEN error when ncols < naxis1 in ffp3dui.
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 3] = [4, 4, 2];
            let data = [0u16; 32];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(
                f.as_deref_mut().unwrap(),
                USHORT_IMG,
                3,
                &naxes,
                &mut status,
            );
            // ncols=2 < naxis1=4, should fail
            fits_write_3d_usht(
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
        // Test BAD_DIMEN error when naxis > 7 in ffpssui.
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 2] = [2, 2];
            let fpixel: [c_long; 2] = [1, 1];
            let lpixel: [c_long; 2] = [2, 2];
            let data = [0u16; 4];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(
                f.as_deref_mut().unwrap(),
                USHORT_IMG,
                2,
                &naxes,
                &mut status,
            );
            // naxis=8 > 7, should fail
            fits_write_subset_usht(
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
    fn test_bad_btable_format() {
        // Test writing to binary table column with unsupported format.
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [u16; 1] = [1];

            let mut f = make_table(&name, BINARY_TBL, "LOGCOL", "1L", 1, &mut status);
            fits_write_col_usht(f.as_deref_mut().unwrap(), 1, 1, 1, 1, &data, &mut status);
            assert_eq!(status, BAD_BTABLE_FORMAT);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_scaled_float() {
        // Test scaled float conversion - exercises ffu2fr4 with scaling.
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [u16; 3] = [100, 200, 300];

            let mut f = make_table(&name, BINARY_TBL, "ECOL", "1E", 3, &mut status);
            fits_set_tscale(f.as_deref_mut().unwrap(), 1, 2.0, 50.0, &mut status);
            fits_write_col_usht(f.as_deref_mut().unwrap(), 1, 1, 1, 3, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0);
        });
    }

    #[test]
    fn test_scaled_double() {
        // Test scaled double conversion - exercises ffu2fr8 with scaling.
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [u16; 3] = [100, 200, 300];

            let mut f = make_table(&name, BINARY_TBL, "DCOL", "1D", 3, &mut status);
            fits_set_tscale(f.as_deref_mut().unwrap(), 1, 2.0, 50.0, &mut status);
            fits_write_col_usht(f.as_deref_mut().unwrap(), 1, 1, 1, 3, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0);
        });
    }

    #[test]
    fn test_scaled_ascii() {
        // Test scaled ASCII conversion - exercises ffu2fstr with scaling.
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [u16; 3] = [100, 200, 300];

            let mut f = make_table(&name, ASCII_TBL, "NUMCOL", "F10.2", 3, &mut status);
            fits_set_tscale(f.as_deref_mut().unwrap(), 1, 2.0, 50.0, &mut status);
            fits_write_col_usht(f.as_deref_mut().unwrap(), 1, 1, 1, 3, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0);
        });
    }

    #[test]
    fn test_overflow_with_nulls() {
        // Test overflow tracking in ffpcnui with fixed-length column. Good
        // values that overflow should set the overflow flag.
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            // 65535 > SHRT_MAX when no TZERO, so it overflows.
            // Sequence: 65535 (overflow), 0 (null), 65535 (overflow), 0 (null)
            let data: [u16; 4] = [65535, 0, 65535, 0];
            let nulval: u16 = 0;

            let mut f = make_table(&name, BINARY_TBL, "ICOL", "5I", 1, &mut status);
            fits_set_btblnull(
                f.as_deref_mut().unwrap(),
                1,
                nulval as LONGLONG,
                &mut status,
            );
            fits_write_colnull_usht(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                4,
                &data,
                nulval,
                &mut status,
            );
            assert_eq!(status, NUM_OVERFLOW);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_scaled_byte_valid() {
        // Test successful scaled conversion to byte. Covers ffu2fi1 path
        // where dvalue is within valid range.
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            // With scale=0.5 and zero=0, stored = (50 - 0) / 0.5 = 100, valid byte
            let data: [u16; 1] = [50];

            let mut f = make_table(&name, BINARY_TBL, "BCOL", "1B", 1, &mut status);
            fits_set_tscale(f.as_deref_mut().unwrap(), 1, 0.5, 0.0, &mut status);
            fits_write_col_usht(f.as_deref_mut().unwrap(), 1, 1, 1, 1, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0);
        });
    }

    #[test]
    fn test_scaled_short_negative_valid() {
        // Test scaled conversion to short with negative dvalue in valid
        // range. Covers ffu2fi2 negative rounding path.
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            // With scale=1 and zero=200, stored = (100 - 200) / 1 = -100, valid
            let data: [u16; 1] = [100];

            let mut f = make_table(&name, BINARY_TBL, "ICOL", "1I", 1, &mut status);
            fits_set_tscale(f.as_deref_mut().unwrap(), 1, 1.0, 200.0, &mut status);
            fits_write_col_usht(f.as_deref_mut().unwrap(), 1, 1, 1, 1, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0);
        });
    }

    #[test]
    fn test_scaled_long_valid() {
        // Test successful scaled conversion to long (int32). Covers ffu2fi4
        // path where dvalue is within valid range.
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            // With scale=0.5 and zero=0, stored = (1000 - 0) / 0.5 = 2000, valid
            let data: [u16; 1] = [1000];

            let mut f = make_table(&name, BINARY_TBL, "JCOL", "1J", 1, &mut status);
            fits_set_tscale(f.as_deref_mut().unwrap(), 1, 0.5, 0.0, &mut status);
            fits_write_col_usht(f.as_deref_mut().unwrap(), 1, 1, 1, 1, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0);
        });
    }

    #[test]
    fn test_scaled_longlong_valid() {
        // Test successful scaled conversion to longlong (int64). Covers
        // ffu2fi8 path where dvalue is within valid range.
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            // With scale=0.5 and zero=0, stored = (1000 - 0) / 0.5 = 2000, valid
            let data: [u16; 1] = [1000];

            let mut f = make_table(&name, BINARY_TBL, "KCOL", "1K", 1, &mut status);
            fits_set_tscale(f.as_deref_mut().unwrap(), 1, 0.5, 0.0, &mut status);
            fits_write_col_usht(f.as_deref_mut().unwrap(), 1, 1, 1, 1, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0);
        });
    }

    #[test]
    fn test_ascii_unscaled_overflow() {
        // Test ASCII table overflow with unscaled value. Covers ffu2fstr
        // overflow check in the scale=1, zero=0 path.
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            // 65535 won't fit in 4 chars with decimal
            let data: [u16; 1] = [65535];

            let mut f = make_table(&name, ASCII_TBL, "FCOL", "F4.1", 1, &mut status);
            fits_write_col_usht(f.as_deref_mut().unwrap(), 1, 1, 1, 1, &data, &mut status);
            assert_eq!(status, NUM_OVERFLOW);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_scaled_long_negative_valid() {
        // Test scaled conversion to long with negative dvalue in valid range.
        // Covers ffu2fi4 negative rounding path.
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            // With scale=1 and zero=200, stored = (100 - 200) / 1 = -100, valid
            let data: [u16; 1] = [100];

            let mut f = make_table(&name, BINARY_TBL, "JCOL", "1J", 1, &mut status);
            fits_set_tscale(f.as_deref_mut().unwrap(), 1, 1.0, 200.0, &mut status);
            fits_write_col_usht(f.as_deref_mut().unwrap(), 1, 1, 1, 1, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0);
        });
    }

    #[test]
    fn test_scaled_longlong_negative_valid() {
        // Test scaled conversion to longlong with negative dvalue in valid
        // range. Covers ffu2fi8 negative rounding path.
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            // With scale=1 and zero=200, stored = (100 - 200) / 1 = -100, valid
            let data: [u16; 1] = [100];

            let mut f = make_table(&name, BINARY_TBL, "KCOL", "1K", 1, &mut status);
            fits_set_tscale(f.as_deref_mut().unwrap(), 1, 1.0, 200.0, &mut status);
            fits_write_col_usht(f.as_deref_mut().unwrap(), 1, 1, 1, 1, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0);
        });
    }

    #[test]
    fn test_variable_length_array_with_nulls() {
        // Test ffpcnui with a variable length array column. Covers tcode < 0
        // path in ffpcnui.
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            // Mix of good values and null marker
            let data: [u16; 5] = [100, 0, 200, 0, 300];
            let nulval: u16 = 0;

            let mut f = make_table(&name, BINARY_TBL, "VARUI", "1PU", 1, &mut status);
            fits_set_btblnull(
                f.as_deref_mut().unwrap(),
                1,
                nulval as LONGLONG,
                &mut status,
            );
            fits_write_colnull_usht(
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
            assert_eq!(status, 0);
        });
    }

    #[test]
    fn test_bad_atable_format() {
        // Test writing unsigned short to ASCII table string column. Covers
        // BAD_ATABLE_FORMAT return in ffpclui.
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [u16; 1] = [100];

            let mut f = make_table(&name, ASCII_TBL, "ACOL", "A5", 1, &mut status);
            fits_write_col_usht(f.as_deref_mut().unwrap(), 1, 1, 1, 1, &data, &mut status);
            assert_eq!(status, BAD_ATABLE_FORMAT);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }
}
