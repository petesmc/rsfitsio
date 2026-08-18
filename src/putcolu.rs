/*  This file, putcolu.rs, contains routines that write data elements to    */
/*  a FITS image or table.  Writes null values.                            */

/*  The FITSIO software was written by William Pence at the High Energy    */
/*  Astrophysic Science Archive Research Center (HEASARC) at the NASA      */
/*  Goddard Space Flight Center.                                           */

use core::cmp;

use crate::c_types::*;

use bytemuck::{cast_slice, cast_slice_mut};

use crate::fitscore::{
    ffgcprll, ffgncl_safe, ffgnrwll_safe, ffgtcl_safe, ffgtclll_safe, ffpmsg_slice, ffpmsg_str,
    fits_is_compressed_image_safe,
};
use crate::fitsio::*;
use crate::fitsio2::*;
use crate::wrappers::*;
use crate::{buffers::*, int_snprintf, slice_to_str};

/*--------------------------------------------------------------------------*/
/// Write null values to the primary array.
///
/// The primary array is represented as a binary table:
/// each group of the primary array is a row in the table,
/// where the first column contains the group parameters
/// and the second column contains the image itself.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffppru(
    fptr: *mut fitsfile, /* I - FITS file pointer                       */
    group: c_long,       /* I - group to write(1 = 1st group)          */
    firstelem: LONGLONG, /* I - first vector element to write(1 = 1st) */
    nelem: LONGLONG,     /* I - number of values to write              */
    status: *mut c_int,  /* IO - error status                          */
) -> c_int {
    // FFI WRAPPER
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        ffppru_safe(fptr, group, firstelem, nelem, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Write null values to the primary array.
///
/// The primary array is represented as a binary table:
/// each group of the primary array is a row in the table,
/// where the first column contains the group parameters
/// and the second column contains the image itself.
pub fn ffppru_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer                       */
    group: c_long,       /* I - group to write(1 = 1st group)          */
    firstelem: LONGLONG, /* I - first vector element to write(1 = 1st) */
    nelem: LONGLONG,     /* I - number of values to write              */
    status: &mut c_int,  /* IO - error status                          */
) -> c_int {
    if fits_is_compressed_image_safe(fptr, status) > 0 {
        /* this is a compressed image in a binary table */
        ffpmsg_str("writing to compressed image is not supported");

        *status = DATA_COMPRESSION_ERR;
        return *status;
    }

    let row = cmp::max(1, group);

    ffpclu_safe(
        fptr,
        2,
        row as LONGLONG,
        firstelem as LONGLONG,
        nelem as LONGLONG,
        status,
    );
    *status
}

/*--------------------------------------------------------------------------*/
/// Write null values to the primary array. (Doesn't support groups).
///
/// The primary array is represented as a binary table:
/// each group of the primary array is a row in the table,
/// where the first column contains the group parameters
/// and the second column contains the image itself.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffpprn(
    fptr: *mut fitsfile, /* I - FITS file pointer                       */
    firstelem: LONGLONG, /* I - first vector element to write(1 = 1st) */
    nelem: LONGLONG,     /* I - number of values to write              */
    status: *mut c_int,  /* IO - error status                          */
) -> c_int {
    // FFI WRAPPER
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        ffpprn_safe(fptr, firstelem, nelem, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Write null values to the primary array. (Doesn't support groups).
///
/// The primary array is represented as a binary table:
/// each group of the primary array is a row in the table,
/// where the first column contains the group parameters
/// and the second column contains the image itself.
pub fn ffpprn_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer                       */
    firstelem: LONGLONG, /* I - first vector element to write(1 = 1st) */
    nelem: LONGLONG,     /* I - number of values to write              */
    status: &mut c_int,  /* IO - error status                          */
) -> c_int {
    let row: c_long = 1;

    if fits_is_compressed_image_safe(fptr, status) > 0 {
        /* this is a compressed image in a binary table */

        ffpmsg_str("writing to compressed image is not supported");

        *status = DATA_COMPRESSION_ERR;
        return *status;
    }

    ffpclu_safe(
        fptr,
        2,
        row as LONGLONG,
        firstelem as LONGLONG,
        nelem as LONGLONG,
        status,
    );
    *status
}

/*--------------------------------------------------------------------------*/
/// Set elements of a table column to the appropriate null value for the column
/// The column number may refer to a real column in an ASCII or binary table,
/// or it may refer to a virtual column in a 1 or more grouped FITS primary
/// array.  FITSIO treats a primary array as a binary table
/// with 2 vector columns: the first column contains the group parameters (often
/// with length = 0) and the second column contains the array of image pixels.
/// Each row of the table represents a group in the case of multigroup FITS
/// images.
///
/// This routine support COMPLEX and DOUBLE COMPLEX binary table columns, and
/// sets both the real and imaginary components of the element to a NaN.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffpclu(
    fptr: *mut fitsfile, /* I - FITS file pointer                       */
    colnum: c_int,       /* I - number of column to write (1 = 1st col) */
    firstrow: LONGLONG,  /* I - first row to write (1 = 1st row)        */
    firstelem: LONGLONG, /* I - first vector element to write (1 = 1st) */
    nelempar: LONGLONG,  /* I - number of values to write               */
    status: *mut c_int,  /* IO - error status                           */
) -> c_int {
    // FFI WRAPPER
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        ffpclu_safe(fptr, colnum, firstrow, firstelem, nelempar, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Set elements of a table column to the appropriate null value for the column
/// The column number may refer to a real column in an ASCII or binary table,
/// or it may refer to a virtual column in a 1 or more grouped FITS primary
/// array.  FITSIO treats a primary array as a binary table
/// with 2 vector columns: the first column contains the group parameters (often
/// with length = 0) and the second column contains the array of image pixels.
/// Each row of the table represents a group in the case of multigroup FITS
/// images.
///
/// This routine support COMPLEX and DOUBLE COMPLEX binary table columns, and
/// sets both the real and imaginary components of the element to a NaN.
pub fn ffpclu_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer                       */
    colnum: c_int,       /* I - number of column to write (1 = 1st col) */
    firstrow: LONGLONG,  /* I - first row to write (1 = 1st row)        */
    firstelem: LONGLONG, /* I - first vector element to write (1 = 1st) */
    nelempar: LONGLONG,  /* I - number of values to write               */
    status: &mut c_int,  /* IO - error status                           */
) -> c_int {
    let mut tcode: c_int = 0;
    let mut maxelem: c_int = 0;
    let mut hdutype: c_int = 0;
    let mut writemode: c_int = 2;
    let mut leng: usize = 0;
    let mut i2null: c_short = 0;
    let mut i4null: c_int = 0;
    let mut twidth: c_long = 0;
    let mut incre: c_long = 0;
    let mut largeelem: LONGLONG = 0;
    let mut nelem: LONGLONG;
    let mut tnull: LONGLONG = 0;
    let mut i8null: LONGLONG = 0;
    let mut repeat: LONGLONG = 0;
    let mut startpos: LONGLONG = 0;
    let mut elemnum: LONGLONG = 0;
    let mut wrtptr: LONGLONG;
    let mut rowlen: LONGLONG = 0;
    let mut rownum: LONGLONG;
    let mut remain: LONGLONG;
    let mut next: LONGLONG;
    let mut ntodo: LONGLONG;
    let mut scale: f64 = 0.0;
    let mut zero: f64 = 0.0;
    let mut i1null: c_uchar = 0;
    let lognul: c_uchar = 0;
    let mut tform: [c_char; 20] = [0; 20];
    let mut message: [c_char; FLEN_ERRMSG] = [0; FLEN_ERRMSG];
    let mut snull: [c_char; 20] = [0; 20]; /*  the FITS null value  */
    let jbuff: [c_long; 2] = [-1, -1]; /* all bits set is equivalent to a NaN */
    let mut buffsize: usize = 0;

    let mut cstring: Vec<c_char> = Vec::new();

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    nelem = nelempar;

    largeelem = firstelem;

    /*---------------------------------------------------*/
    /*  Check input and get parameters about the column: */
    /*---------------------------------------------------*/

    /* note that writemode = 2 by default (not 1), so that the returned */
    /* repeat and incre values will be the actual values for this column. */

    /* If writing nulls to a variable length column then dummy data values  */
    /* must have already been written to the heap. */
    /* We just have to overwrite the previous values with null values. */
    /* Set writemode = 0 in this case, to test that values have been written */

    ffgtcl_safe(fptr, colnum, Some(&mut tcode), None, None, status);
    if tcode < 0 {
        writemode = 0; /* this is a variable length column */
    }

    if tcode.abs() >= TCOMPLEX {
        /* treat complex columns as pairs of numbers */
        largeelem = (largeelem - 1) * 2 + 1;
        nelem *= 2;
    }

    if ffgcprll(
        fptr,
        colnum,
        firstrow,
        largeelem,
        nelem,
        writemode,
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
        if snull[0] == ASCII_NULL_UNDEFINED {
            ffpmsg_str("Null value string for ASCII table column is not defined (FTPCLU).");
            *status = NO_NULL;
            return *status;
        }

        /* allocate buffer to hold the null string.  Must write the entire */
        /* width of the column (twidth bytes) to avoid possible problems */
        /* with uninitialized FITS blocks, in case the field spans blocks */

        buffsize = cmp::max(20, twidth) as usize;
        cstring = vec![b' ' as c_char; buffsize]; /* initialize  with blanks */

        leng = strlen_safe(&snull);
        if hdutype == BINARY_TBL {
            leng += 1; /* copy the terminator too in binary tables */
        }
        strncpy_safe(&mut cstring, &snull, leng); /* copy null string to temp buffer */
    } else if tcode == TBYTE || tcode == TSHORT || tcode == TLONG || tcode == TLONGLONG {
        if tnull == LONGLONG::from(NULL_UNDEFINED) {
            ffpmsg_str("Null value for integer table column is not defined (FTPCLU).");
            *status = NO_NULL;
            return *status;
        }

        if tcode == TBYTE {
            i1null = tnull as u8;
        } else if tcode == TSHORT {
            i2null = tnull as c_short;
            if BYTESWAPPED {
                i2null = i2null.swap_bytes(); /* reverse order of bytes */
            }
        } else if tcode == TLONG {
            i4null = tnull as INT32BIT;
            if BYTESWAPPED {
                i4null = i4null.swap_bytes(); /* reverse order of bytes */
            }
        } else {
            i8null = tnull;
            if BYTESWAPPED {
                i8null = i8null.swap_bytes(); /* reverse order of bytes */
            }
        }
    }

    /*---------------------------------------------------------------------*/
    /*  Now write the pixels to the FITS column.                           */
    /*---------------------------------------------------------------------*/
    remain = nelem; /* remaining number of values to write  */
    next = 0; /* next element in array to be written  */
    rownum = 0; /* row number, relative to firstrow     */
    ntodo = remain; /* number of elements to write at one time */

    while ntodo > 0 {
        /* limit the number of pixels to process at one time to the number that
           will fit in the buffer space or to the number of pixels that remain
           in the current vector, which ever is smaller.
        */
        ntodo = cmp::min(ntodo, repeat - elemnum);
        wrtptr = startpos + ((rownum as LONGLONG) * rowlen) + (elemnum * incre as LONGLONG);

        ffmbyt_safe(fptr, wrtptr, IGNORE_EOF, status); /* move to write position */

        match tcode {
            TBYTE => {
                for _ii in 0..(ntodo as usize) {
                    ffpbyt(fptr, 1, &[i1null], status);
                }
            }
            TSHORT => {
                for _ii in 0..(ntodo as usize) {
                    ffpbyt(fptr, 2, cast_slice(&[i2null]), status);
                }
            }
            TLONG => {
                for _ii in 0..(ntodo as usize) {
                    ffpbyt(fptr, 4, cast_slice(&[i4null]), status);
                }
            }
            TLONGLONG => {
                for _ii in 0..(ntodo as usize) {
                    ffpbyt(fptr, 8, cast_slice(&[i8null]), status);
                }
            }
            TFLOAT => {
                for _ii in 0..(ntodo as usize) {
                    ffpbyt(fptr, 4, cast_slice(&jbuff), status);
                }
            }
            TDOUBLE => {
                for _ii in 0..(ntodo as usize) {
                    ffpbyt(fptr, 8, cast_slice(&jbuff), status);
                }
            }
            TLOGICAL => {
                for _ii in 0..(ntodo as usize) {
                    ffpbyt(fptr, 1, cast_slice(&[lognul]), status);
                }
            }
            TSTRING => {
                /* an ASCII table column */
                /* repeat always = 1, so ntodo is also guaranteed to = 1 */
                ffpbyt(fptr, twidth as LONGLONG, cast_slice(&cstring), status);
            }
            _ => {
                /*  error trap  */
                int_snprintf!(
                    &mut message,
                    FLEN_ERRMSG,
                    "Cannot write null value to column {} which has format {}",
                    colnum,
                    slice_to_str!(&tform),
                );
                ffpmsg_slice(&message);
                return *status;
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
                "Error writing {:.0} thru {:.0} of null values (ffpclu).",
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
        ntodo = remain; /* this is the maximum number to do in next loop */
    } /*  End of main while Loop  */

    *status
}

/*--------------------------------------------------------------------------*/
/// Set elements of a table column to the appropriate null value for the column
/// The column number may refer to a real column in an ASCII or binary table,
/// or it may refer to a virtual column in a 1 or more grouped FITS primary
/// array.  FITSIO treats a primary array as a binary table
/// with 2 vector columns: the first column contains the group parameters (often
/// with length = 0) and the second column contains the array of image pixels.
/// Each row of the table represents a group in the case of multigroup FITS
/// images.
///
/// This routine does not do anything special in the case of COMPLEX table columns
/// (unlike the similar ffpclu routine).  This routine is mainly for use by
/// ffpcne which already compensates for the effective doubling of the number of
/// elements in a complex column.
pub(crate) fn ffpcluc(
    fptr: &mut fitsfile, /* I - FITS file pointer                       */
    colnum: c_int,       /* I - number of column to write (1 = 1st col) */
    firstrow: LONGLONG,  /* I - first row to write (1 = 1st row)        */
    firstelem: LONGLONG, /* I - first vector element to write (1 = 1st) */
    nelem: LONGLONG,     /* I - number of values to write               */
    status: &mut c_int,  /* IO - error status                           */
) -> c_int {
    let mut tcode: c_int = 0;
    let mut maxelem: c_int = 0;
    let mut hdutype: c_int = 0;
    let mut writemode: c_int = 2;
    let mut leng: usize = 0;
    let mut i2null: c_short = 0;
    let mut i4null: c_int = 0;
    let mut twidth: c_long = 0;
    let mut incre: c_long = 0;
    let mut tnull: LONGLONG = 0;
    let mut i8null: LONGLONG = 0;
    let mut repeat: LONGLONG = 0;
    let mut startpos: LONGLONG = 0;
    let mut elemnum: LONGLONG = 0;
    let mut wrtptr: LONGLONG = 0;
    let mut rowlen: LONGLONG = 0;
    let mut rownum: LONGLONG;
    let mut remain: LONGLONG;
    let mut next: LONGLONG;
    let mut ntodo: LONGLONG;
    let mut scale: f64 = 0.0;
    let mut zero: f64 = 0.0;
    let mut i1null: c_uchar = 0;
    let lognul: c_uchar = 0;
    let mut tform: [c_char; 20] = [0; 20];
    let mut message: [c_char; FLEN_ERRMSG] = [0; FLEN_ERRMSG];
    let mut snull: [c_char; 20] = [0; 20]; /*  the FITS null value  */
    let jbuff: [c_long; 2] = [-1, -1]; /* all bits set is equivalent to a NaN */

    let mut cstring: Vec<c_char> = Vec::new();

    let mut buffsize: usize = 0;

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    /*---------------------------------------------------*/
    /*  Check input and get parameters about the column: */
    /*---------------------------------------------------*/

    /* note that writemode = 2 by default (not 1), so that the returned */
    /* repeat and incre values will be the actual values for this column. */

    /* If writing nulls to a variable length column then dummy data values  */
    /* must have already been written to the heap. */
    /* We just have to overwrite the previous values with null values. */
    /* Set writemode = 0 in this case, to test that values have been written */

    ffgtcl_safe(fptr, colnum, Some(&mut tcode), None, None, status);

    if tcode < 0 {
        writemode = 0; /* this is a variable length column */
    }

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
        if snull[0] == ASCII_NULL_UNDEFINED {
            ffpmsg_str("Null value string for ASCII table column is not defined (FTPCLU).");
            *status = NO_NULL;
            return *status;
        }

        /* allocate buffer to hold the null string.  Must write the entire */
        /* width of the column (twidth bytes) to avoid possible problems */
        /* with uninitialized FITS blocks, in case the field spans blocks */

        buffsize = cmp::max(20, twidth as usize);

        cstring = vec![b' ' as c_char; buffsize]; /* initialize  with blanks */

        leng = strlen_safe(&snull);
        if hdutype == BINARY_TBL {
            leng += 1; /* copy the termsnator too in binary tables */
        }
        strncpy_safe(&mut cstring, &snull, leng); /* copy null string to temp buffer */
    } else if tcode == TBYTE || tcode == TSHORT || tcode == TLONG || tcode == TLONGLONG {
        if tnull == LONGLONG::from(NULL_UNDEFINED) {
            ffpmsg_str("Null value for integer table column is not defined (FTPCLU).");
            *status = NO_NULL;
            return *status;
        }

        if tcode == TBYTE {
            i1null = tnull as u8;
        } else if tcode == TSHORT {
            i2null = tnull as c_short;
            if BYTESWAPPED {
                i2null = i2null.swap_bytes(); /* reverse order of bytes */
            }
        } else if tcode == TLONG {
            i4null = tnull as INT32BIT;
            if BYTESWAPPED {
                i4null = i4null.swap_bytes(); /* reverse order of bytes */
            }
        } else {
            i8null = tnull;
            if BYTESWAPPED {
                i8null = i8null.swap_bytes(); /* reverse order of bytes */
            }
        }
    }

    /*---------------------------------------------------------------------*/
    /*  Now write the pixels to the FITS column.                           */
    /*---------------------------------------------------------------------*/
    remain = nelem; /* remaining number of values to write  */
    next = 0; /* next element in array to be written  */
    rownum = 0; /* row number, relative to firstrow     */
    ntodo = remain; /* number of elements to write at one time */

    while ntodo > 0 {
        /* limit the number of pixels to process at one time to the number that
           will fit in the buffer space or to the number of pixels that remain
           in the current vector, which ever is smaller.
        */
        ntodo = cmp::min(ntodo, repeat - elemnum);
        wrtptr = startpos + ((rownum as LONGLONG) * rowlen) + (elemnum * incre as LONGLONG);

        ffmbyt_safe(fptr, wrtptr, IGNORE_EOF, status); /* move to write position */

        match tcode {
            TBYTE => {
                for _ii in 0..(ntodo as usize) {
                    ffpbyt(fptr, 1, &[i1null], status);
                }
            }
            TSHORT => {
                for _ii in 0..(ntodo as usize) {
                    ffpbyt(fptr, 2, cast_slice(&[i2null]), status);
                }
            }
            TLONG => {
                for _ii in 0..(ntodo as usize) {
                    ffpbyt(fptr, 4, cast_slice(&[i4null]), status);
                }
            }
            TLONGLONG => {
                for _ii in 0..(ntodo as usize) {
                    ffpbyt(fptr, 8, cast_slice(&[i8null]), status);
                }
            }
            TFLOAT => {
                for _ii in 0..(ntodo as usize) {
                    ffpbyt(fptr, 4, cast_slice(&jbuff), status);
                }
            }
            TDOUBLE => {
                for _ii in 0..(ntodo as usize) {
                    ffpbyt(fptr, 8, cast_slice(&jbuff), status);
                }
            }
            TLOGICAL => {
                for _ii in 0..(ntodo as usize) {
                    ffpbyt(fptr, 1, cast_slice(&[lognul]), status);
                }
            }
            TSTRING => {
                /* an ASCII table column */
                /* repeat always = 1, so ntodo is also guaranteed to = 1 */
                ffpbyt(
                    fptr,
                    twidth as LONGLONG,
                    cast_slice(cstring.as_slice()),
                    status,
                );
            }
            _ => {
                /*  error trap  */
                int_snprintf!(
                    &mut message,
                    FLEN_ERRMSG,
                    "Cannot write null value to column {} which has format {}",
                    colnum,
                    slice_to_str!(&tform),
                );
                ffpmsg_slice(&message);
                return *status;
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
                "Error writing {:.0} thru {:.0} of null values (ffpclu).",
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
        ntodo = remain; /* this is the maximum number to do in next loop */
    } /*  End of main while Loop  */

    *status
}

/*--------------------------------------------------------------------------*/
/// fits_write_nullrows / ffprwu - write TNULLs to all columns in one or more rows
///
/// fitsfile *fptr - pointer to FITS HDU opened for read/write
/// long int firstrow - first table row to set to null. (firstrow >= 1)
/// long int nrows - total number or rows to set to null. (nrows >= 1)
/// int *status - upon return, *status contains CFITSIO status code
///
/// RETURNS: CFITSIO status code
///
/// written by Craig Markwardt, GSFC
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffprwu(
    fptr: *mut fitsfile,
    firstrow: LONGLONG,
    nrows: LONGLONG,
    status: *mut c_int,
) -> c_int {
    // FFI WRAPPER
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        ffprwu_safe(fptr, firstrow, nrows, status)
    }
}

/*--------------------------------------------------------------------------*/
/// fits_write_nullrows / ffprwu - write TNULLs to all columns in one or more rows
///
/// fitsfile *fptr - pointer to FITS HDU opened for read/write
/// long int firstrow - first table row to set to null. (firstrow >= 1)
/// long int nrows - total number or rows to set to null. (nrows >= 1)
/// int *status - upon return, *status contains CFITSIO status code
///
/// RETURNS: CFITSIO status code
///
/// written by Craig Markwardt, GSFC
pub fn ffprwu_safe(
    fptr: &mut fitsfile,
    firstrow: LONGLONG,
    nrows: LONGLONG,
    status: &mut c_int,
) -> c_int {
    let mut ntotrows: LONGLONG = 0;
    let mut ncols: c_int = 0;
    let mut typecode: c_int;
    let mut repeat: LONGLONG;
    let mut width;
    let mut nullstatus: c_int;

    if *status > 0 {
        return *status;
    }

    if (firstrow <= 0) || (nrows <= 0) {
        *status = BAD_ROW_NUM;
        return *status;
    }

    ffgnrwll_safe(fptr, &mut ntotrows, status);

    if firstrow + nrows - 1 > ntotrows {
        *status = BAD_ROW_NUM;
        return *status;
    }

    ffgncl_safe(fptr, &mut ncols, status);
    if *status > 0 {
        return *status;
    }

    /* Loop through each column and write nulls */
    for i in 1..=(ncols as usize) {
        repeat = 0;
        typecode = 0;
        width = 0;
        ffgtclll_safe(
            fptr,
            i as c_int,
            Some(&mut typecode),
            Some(&mut repeat),
            Some(&mut width),
            status,
        );
        if *status > 0 {
            break;
        }

        /* NOTE: data of TSTRING type must not write the total repeat
        count, since the repeat count is the *character* count, not the
        nstring count.  Divide by string width to get number of
        strings. */

        if typecode == TSTRING {
            repeat /= width;
        }

        /* Write NULLs */
        nullstatus = 0;
        ffpclu_safe(
            fptr,
            i as c_int,
            firstrow,
            1,
            repeat * nrows,
            &mut nullstatus,
        );

        /* ignore error if no null value is defined for the column */
        if nullstatus > 0 && nullstatus != NO_NULL {
            *status = nullstatus;
            return *status;
        }
    }

    *status
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::aliases::rust_api::*;
    use crate::fitsio::{
        ASCII_TBL, BAD_ROW_NUM, BINARY_TBL, BYTE_IMG, DOUBLE_IMG, FLOAT_IMG, LONGLONG, NO_NULL,
        TDOUBLE, TFLOAT, fitsfile,
    };
    use crate::helpers::testhelpers::{to_buf, with_temp_file};
    use bytemuck::cast_slice;
    use libc::{c_char, c_int, c_long};

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
    fn test_write_null_primary_float() {
        // Float/double images use NaN for null, no TNULL needed
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 1] = [10];
            let data: [f32; 10] = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), FLOAT_IMG, 1, &naxes, &mut status);
            fits_write_img(
                f.as_deref_mut().unwrap(),
                TFLOAT,
                1,
                10,
                cast_slice(&data),
                &mut status,
            );
            fits_write_img_null(f.as_deref_mut().unwrap(), 1, 1, 5, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0);
        });
    }

    #[test]
    fn test_write_null_primary_double() {
        // Double image uses NaN for null
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 1] = [10];
            let data: [f64; 10] = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(
                f.as_deref_mut().unwrap(),
                DOUBLE_IMG,
                1,
                &naxes,
                &mut status,
            );
            fits_write_img(
                f.as_deref_mut().unwrap(),
                TDOUBLE,
                1,
                10,
                cast_slice(&data),
                &mut status,
            );
            fits_write_null_img(f.as_deref_mut().unwrap(), 1, 5, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0);
        });
    }

    #[test]
    fn test_write_null_byte_column() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [u8; 10] = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10];

            let mut f = make_table(&name, BINARY_TBL, "COL1", "10B", 1, &mut status);
            fits_set_btblnull(f.as_deref_mut().unwrap(), 1, 255, &mut status);
            fits_write_col_byt(f.as_deref_mut().unwrap(), 1, 1, 1, 10, &data, &mut status);
            fits_write_col_null(f.as_deref_mut().unwrap(), 1, 1, 1, 3, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0);
        });
    }

    #[test]
    fn test_write_null_short_column() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [i16; 5] = [100, 200, 300, 400, 500];

            let mut f = make_table(&name, BINARY_TBL, "COL1", "5I", 1, &mut status);
            fits_set_btblnull(f.as_deref_mut().unwrap(), 1, -32768, &mut status);
            fits_write_col_sht(f.as_deref_mut().unwrap(), 1, 1, 1, 5, &data, &mut status);
            fits_write_col_null(f.as_deref_mut().unwrap(), 1, 1, 1, 2, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0);
        });
    }

    #[test]
    fn test_write_null_long_column() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [c_long; 3] = [1000, 2000, 3000];

            let mut f = make_table(&name, BINARY_TBL, "COL1", "3J", 1, &mut status);
            fits_set_btblnull(f.as_deref_mut().unwrap(), 1, -2147483647, &mut status);
            fits_write_col_lng(f.as_deref_mut().unwrap(), 1, 1, 1, 3, &data, &mut status);
            fits_write_col_null(f.as_deref_mut().unwrap(), 1, 1, 1, 1, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0);
        });
    }

    #[test]
    fn test_write_null_longlong_column() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [LONGLONG; 2] = [1000000000, 2000000000];

            let mut f = make_table(&name, BINARY_TBL, "COL1", "2K", 1, &mut status);
            fits_set_btblnull(
                f.as_deref_mut().unwrap(),
                1,
                -9223372036854775807,
                &mut status,
            );
            fits_write_col_lnglng(f.as_deref_mut().unwrap(), 1, 1, 1, 2, &data, &mut status);
            fits_write_col_null(f.as_deref_mut().unwrap(), 1, 1, 1, 1, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0);
        });
    }

    #[test]
    fn test_write_null_float_column() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [f32; 4] = [1.1, 2.2, 3.3, 4.4];

            let mut f = make_table(&name, BINARY_TBL, "COL1", "4E", 1, &mut status);
            fits_write_col_flt(f.as_deref_mut().unwrap(), 1, 1, 1, 4, &data, &mut status);
            fits_write_col_null(f.as_deref_mut().unwrap(), 1, 1, 1, 2, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0);
        });
    }

    #[test]
    fn test_write_null_double_column() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [f64; 3] = [1.11, 2.22, 3.33];

            let mut f = make_table(&name, BINARY_TBL, "COL1", "3D", 1, &mut status);
            fits_write_col_dbl(f.as_deref_mut().unwrap(), 1, 1, 1, 3, &data, &mut status);
            fits_write_col_null(f.as_deref_mut().unwrap(), 1, 1, 1, 1, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0);
        });
    }

    #[test]
    fn test_write_null_logical_column() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [c_char; 5] = [
                b'T' as c_char,
                b'F' as c_char,
                b'T' as c_char,
                b'F' as c_char,
                b'T' as c_char,
            ];

            let mut f = make_table(&name, BINARY_TBL, "COL1", "5L", 1, &mut status);
            fits_write_col_log(f.as_deref_mut().unwrap(), 1, 1, 1, 5, &data, &mut status);
            fits_write_col_null(f.as_deref_mut().unwrap(), 1, 1, 2, 2, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0);
        });
    }

    #[test]
    fn test_write_null_string_column() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let d0 = cc("hello");
            let d1 = cc("world");
            let data: [&[c_char]; 2] = [d0.as_slice(), d1.as_slice()];

            let mut f = make_table(&name, BINARY_TBL, "COL1", "10A", 2, &mut status);
            fits_write_key_str(
                f.as_deref_mut().unwrap(),
                &cc("TNULL1"),
                &cc("*"),
                None,
                &mut status,
            );
            fits_write_col_str(f.as_deref_mut().unwrap(), 1, 1, 1, 2, &data, &mut status);
            fits_write_col_null(f.as_deref_mut().unwrap(), 1, 1, 1, 1, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0);
        });
    }

    #[test]
    fn test_write_null_rows() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let fdata: [f32; 5] = [1.0, 2.0, 3.0, 4.0, 5.0];
            let ddata: [f64; 3] = [1.1, 2.2, 3.3];

            // 2-column table: COL1 = 5E, COL2 = 3D, with 3 rows.
            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            let ttype_v = [Some(cc("COL1")), Some(cc("COL2"))];
            let ttype_ref: Vec<Option<&[c_char]>> = ttype_v.iter().map(|o| o.as_deref()).collect();
            let tform_v = [cc("5E"), cc("3D")];
            let tform_ref: Vec<&[c_char]> = tform_v.iter().map(|v| v.as_slice()).collect();
            fits_create_tbl(
                f.as_deref_mut().unwrap(),
                BINARY_TBL,
                3,
                2,
                &ttype_ref,
                &tform_ref,
                None,
                None,
                &mut status,
            );
            fits_write_col_flt(f.as_deref_mut().unwrap(), 1, 1, 1, 5, &fdata, &mut status);
            fits_write_col_flt(f.as_deref_mut().unwrap(), 1, 2, 1, 5, &fdata, &mut status);
            fits_write_col_flt(f.as_deref_mut().unwrap(), 1, 3, 1, 5, &fdata, &mut status);
            fits_write_col_dbl(f.as_deref_mut().unwrap(), 2, 1, 1, 3, &ddata, &mut status);
            fits_write_col_dbl(f.as_deref_mut().unwrap(), 2, 2, 1, 3, &ddata, &mut status);
            fits_write_col_dbl(f.as_deref_mut().unwrap(), 2, 3, 1, 3, &ddata, &mut status);
            fits_write_nullrows(f.as_deref_mut().unwrap(), 2, 1, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0);
        });
    }

    #[test]
    fn test_write_null_complex_column() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [f32; 6] = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];

            let mut f = make_table(&name, BINARY_TBL, "COL1", "3C", 1, &mut status);
            fits_write_col_flt(f.as_deref_mut().unwrap(), 1, 1, 1, 6, &data, &mut status);
            fits_write_col_null(f.as_deref_mut().unwrap(), 1, 1, 1, 2, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0);
        });
    }

    #[test]
    fn test_pclu_bad_status() {
        // Pre-existing error status must be inherited unchanged.
        with_temp_file(|filename| {
            let name = to_buf(filename);
            let mut status = 0;
            let mut f = make_table(&name, BINARY_TBL, "COL1", "1B", 1, &mut status);
            assert_eq!(status, 0);
            status = 1;
            fits_write_col_null(f.as_deref_mut().unwrap(), 1, 1, 1, 1, &mut status);
            assert!(status == 1);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_prwu_bad_row() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [f32; 5] = [1.0, 2.0, 3.0, 4.0, 5.0];

            let mut f = make_table(&name, BINARY_TBL, "COL1", "5E", 1, &mut status);
            fits_write_col_flt(f.as_deref_mut().unwrap(), 1, 1, 1, 5, &data, &mut status);
            fits_write_col_flt(f.as_deref_mut().unwrap(), 1, 2, 1, 5, &data, &mut status);
            status = 0;
            fits_write_nullrows(f.as_deref_mut().unwrap(), 0, 1, &mut status);
            assert!(status == BAD_ROW_NUM);
            status = 0;
            fits_write_nullrows(f.as_deref_mut().unwrap(), 1, 10, &mut status);
            assert!(status == BAD_ROW_NUM);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0);
        });
    }

    #[test]
    fn test_pclu_no_null_defined() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [u8; 10] = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10];

            let mut f = make_table(&name, BINARY_TBL, "COL1", "10B", 1, &mut status);
            fits_write_col_byt(f.as_deref_mut().unwrap(), 1, 1, 1, 10, &data, &mut status);
            status = 0;
            fits_write_col_null(f.as_deref_mut().unwrap(), 1, 1, 1, 3, &mut status);
            assert!(status == NO_NULL);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0);
        });
    }

    #[test]
    fn test_ascii_table_null() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [c_long; 3] = [100, 200, 300];

            // ASCII table, single I10 column with explicit empty unit.
            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            let ttype_v = [Some(cc("COL1"))];
            let ttype_ref: Vec<Option<&[c_char]>> = ttype_v.iter().map(|o| o.as_deref()).collect();
            let tform_v = [cc("I10")];
            let tform_ref: Vec<&[c_char]> = tform_v.iter().map(|v| v.as_slice()).collect();
            let tunit_v = [Some(cc(""))];
            let tunit_ref: Vec<Option<&[c_char]>> = tunit_v.iter().map(|o| o.as_deref()).collect();
            fits_create_tbl(
                f.as_deref_mut().unwrap(),
                ASCII_TBL,
                3,
                1,
                &ttype_ref,
                &tform_ref,
                Some(&tunit_ref),
                None,
                &mut status,
            );
            fits_write_key_str(
                f.as_deref_mut().unwrap(),
                &cc("TNULL1"),
                &cc("-999"),
                None,
                &mut status,
            );
            fits_write_col_lng(f.as_deref_mut().unwrap(), 1, 1, 1, 3, &data, &mut status);
            fits_write_col_null(f.as_deref_mut().unwrap(), 1, 2, 1, 1, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0);
        });
    }

    #[test]
    fn test_pcluc_byte() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [u8; 10] = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10];

            let mut f = make_table(&name, BINARY_TBL, "COL1", "10B", 1, &mut status);
            fits_set_btblnull(f.as_deref_mut().unwrap(), 1, 255, &mut status);
            fits_write_col_byt(f.as_deref_mut().unwrap(), 1, 1, 1, 10, &data, &mut status);
            ffpcluc(f.as_deref_mut().unwrap(), 1, 1, 1, 3, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0);
        });
    }

    #[test]
    fn test_pcluc_longlong() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [LONGLONG; 5] = [100, 200, 300, 400, 500];

            let mut f = make_table(&name, BINARY_TBL, "COL1", "5K", 1, &mut status);
            fits_set_btblnull(
                f.as_deref_mut().unwrap(),
                1,
                -9223372036854775807,
                &mut status,
            );
            fits_write_col_lnglng(f.as_deref_mut().unwrap(), 1, 1, 1, 5, &data, &mut status);
            ffpcluc(f.as_deref_mut().unwrap(), 1, 1, 1, 2, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0);
        });
    }

    #[test]
    fn test_pcluc_logical() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [c_char; 10] = [
                b'T' as c_char,
                b'F' as c_char,
                b'T' as c_char,
                b'F' as c_char,
                b'T' as c_char,
                b'T' as c_char,
                b'F' as c_char,
                b'T' as c_char,
                b'F' as c_char,
                b'T' as c_char,
            ];

            let mut f = make_table(&name, BINARY_TBL, "COL1", "10L", 1, &mut status);
            fits_write_col_log(f.as_deref_mut().unwrap(), 1, 1, 1, 10, &data, &mut status);
            ffpcluc(f.as_deref_mut().unwrap(), 1, 1, 3, 4, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0);
        });
    }

    #[test]
    fn test_pcluc_string() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let d0 = cc("hello");
            let d1 = cc("world");
            let d2 = cc("test");
            let data: [&[c_char]; 3] = [d0.as_slice(), d1.as_slice(), d2.as_slice()];

            let mut f = make_table(&name, BINARY_TBL, "COL1", "10A", 3, &mut status);
            fits_write_key_str(
                f.as_deref_mut().unwrap(),
                &cc("TNULL1"),
                &cc("*"),
                None,
                &mut status,
            );
            fits_write_col_str(f.as_deref_mut().unwrap(), 1, 1, 1, 3, &data, &mut status);
            ffpcluc(f.as_deref_mut().unwrap(), 1, 2, 1, 1, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0);
        });
    }

    #[test]
    fn test_pcluc_float() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [f32; 5] = [1.1, 2.2, 3.3, 4.4, 5.5];

            let mut f = make_table(&name, BINARY_TBL, "COL1", "5E", 1, &mut status);
            fits_write_col_flt(f.as_deref_mut().unwrap(), 1, 1, 1, 5, &data, &mut status);
            ffpcluc(f.as_deref_mut().unwrap(), 1, 1, 1, 2, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0);
        });
    }

    #[test]
    fn test_pcluc_double() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [f64; 4] = [1.11, 2.22, 3.33, 4.44];

            let mut f = make_table(&name, BINARY_TBL, "COL1", "4D", 1, &mut status);
            fits_write_col_dbl(f.as_deref_mut().unwrap(), 1, 1, 1, 4, &data, &mut status);
            ffpcluc(f.as_deref_mut().unwrap(), 1, 1, 1, 2, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0);
        });
    }

    #[test]
    fn test_pcluc_bad_status() {
        // Pre-existing error status must be inherited unchanged.
        with_temp_file(|filename| {
            let name = to_buf(filename);
            let mut status = 0;
            let mut f = make_table(&name, BINARY_TBL, "COL1", "1B", 1, &mut status);
            assert_eq!(status, 0);
            status = 1;
            ffpcluc(f.as_deref_mut().unwrap(), 1, 1, 1, 1, &mut status);
            assert!(status == 1);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    // Test writing nulls that span multiple rows.
    #[test]
    fn test_pcluc_multirow() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [f32; 5] = [1.1, 2.2, 3.3, 4.4, 5.5];

            let mut f = make_table(&name, BINARY_TBL, "COL1", "5E", 3, &mut status);
            fits_write_col_flt(f.as_deref_mut().unwrap(), 1, 1, 1, 5, &data, &mut status);
            fits_write_col_flt(f.as_deref_mut().unwrap(), 1, 2, 1, 5, &data, &mut status);
            fits_write_col_flt(f.as_deref_mut().unwrap(), 1, 3, 1, 5, &data, &mut status);
            // Write nulls across rows - start at elem 4 of row 1, continue to row 2
            ffpcluc(f.as_deref_mut().unwrap(), 1, 1, 4, 7, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0);
        });
    }

    // Test writing nulls that span multiple rows with ffpclu.
    #[test]
    fn test_pclu_multirow() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [f64; 4] = [1.1, 2.2, 3.3, 4.4];

            let mut f = make_table(&name, BINARY_TBL, "COL1", "4D", 3, &mut status);
            fits_write_col_dbl(f.as_deref_mut().unwrap(), 1, 1, 1, 4, &data, &mut status);
            fits_write_col_dbl(f.as_deref_mut().unwrap(), 1, 2, 1, 4, &data, &mut status);
            fits_write_col_dbl(f.as_deref_mut().unwrap(), 1, 3, 1, 4, &data, &mut status);
            // Write nulls across rows - start at elem 3 of row 1, continue to row 3
            fits_write_col_null(f.as_deref_mut().unwrap(), 1, 1, 3, 8, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0);
        });
    }

    // Test NO_NULL error for ASCII table string column without null defined.
    // This exercises the TSTRING path in ffpclu.
    #[test]
    fn test_pclu_ascii_string_no_null() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);

            // ASCII table, single A8 string column, no null string defined.
            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            let ttype_v = [Some(cc("STRCOL"))];
            let ttype_ref: Vec<Option<&[c_char]>> = ttype_v.iter().map(|o| o.as_deref()).collect();
            let tform_v = [cc("A8")];
            let tform_ref: Vec<&[c_char]> = tform_v.iter().map(|v| v.as_slice()).collect();
            fits_create_tbl(
                f.as_deref_mut().unwrap(),
                ASCII_TBL,
                1,
                1,
                &ttype_ref,
                &tform_ref,
                None,
                None,
                &mut status,
            );
            // Try to write nulls without defining null string - should fail.
            fits_write_col_null(f.as_deref_mut().unwrap(), 1, 1, 1, 1, &mut status);
            assert!(status == NO_NULL);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0);
        });
    }

    // Test NO_NULL error for ASCII table string column without null defined.
    // This exercises the TSTRING path in ffpcluc.
    #[test]
    fn test_pcluc_ascii_string_no_null() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            let ttype_v = [Some(cc("STRCOL"))];
            let ttype_ref: Vec<Option<&[c_char]>> = ttype_v.iter().map(|o| o.as_deref()).collect();
            let tform_v = [cc("A8")];
            let tform_ref: Vec<&[c_char]> = tform_v.iter().map(|v| v.as_slice()).collect();
            fits_create_tbl(
                f.as_deref_mut().unwrap(),
                ASCII_TBL,
                1,
                1,
                &ttype_ref,
                &tform_ref,
                None,
                None,
                &mut status,
            );
            // Try to write nulls without defining null string - should fail.
            ffpcluc(f.as_deref_mut().unwrap(), 1, 1, 1, 1, &mut status);
            assert!(status == NO_NULL);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0);
        });
    }

    // Test NO_NULL error for integer column without null defined in ffpcluc.
    // This exercises the integer NO_NULL path in ffpcluc.
    #[test]
    fn test_pcluc_no_null_int() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [c_long; 3] = [1, 2, 3];

            let mut f = make_table(&name, BINARY_TBL, "INTCOL", "1J", 3, &mut status);
            fits_write_col_lng(f.as_deref_mut().unwrap(), 1, 1, 1, 3, &data, &mut status);
            // Try to write nulls without defining TNULL - should fail.
            ffpcluc(f.as_deref_mut().unwrap(), 1, 1, 1, 1, &mut status);
            assert!(status == NO_NULL);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0);
        });
    }
}
