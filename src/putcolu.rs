/*  This file, putcolu.c, contains routines that write data elements to    */
/*  a FITS image or table.  Writes null values.                            */

/*  The FITSIO software was written by William Pence at the High Energy    */
/*  Astrophysic Science Archive Research Center (HEASARC) at the NASA      */
/*  Goddard Space Flight Center.                                           */

use std::cmp;
use std::ffi::CStr;

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
    let ii: LONGLONG = 0;
    let mut largeelem: LONGLONG = 0;
    let mut nelem: LONGLONG = 0;
    let mut tnull: LONGLONG = 0;
    let mut i8null: LONGLONG = 0;
    let mut repeat: LONGLONG = 0;
    let mut startpos: LONGLONG = 0;
    let mut elemnum: LONGLONG = 0;
    let mut wrtptr: LONGLONG = 0;
    let mut rowlen: LONGLONG = 0;
    let mut rownum: LONGLONG = 0;
    let mut remain: LONGLONG = 0;
    let mut next: LONGLONG = 0;
    let mut ntodo: LONGLONG = 0;
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
        if tnull == NULL_UNDEFINED as LONGLONG {
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
                for ii in 0..(ntodo as usize) {
                    ffpbyt(fptr, 1, &[i1null], status);
                }
            }
            TSHORT => {
                for ii in 0..(ntodo as usize) {
                    ffpbyt(fptr, 2, cast_slice(&[i2null]), status);
                }
            }
            TLONG => {
                for ii in 0..(ntodo as usize) {
                    ffpbyt(fptr, 4, cast_slice(&[i4null]), status);
                }
            }
            TLONGLONG => {
                for ii in 0..(ntodo as usize) {
                    ffpbyt(fptr, 8, cast_slice(&[i8null]), status);
                }
            }
            TFLOAT => {
                for ii in 0..(ntodo as usize) {
                    ffpbyt(fptr, 4, cast_slice(&jbuff), status);
                }
            }
            TDOUBLE => {
                for ii in 0..(ntodo as usize) {
                    ffpbyt(fptr, 8, cast_slice(&jbuff), status);
                }
            }
            TLOGICAL => {
                for ii in 0..(ntodo as usize) {
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
    let ii: LONGLONG = 0;
    let mut tnull: LONGLONG = 0;
    let mut i8null: LONGLONG = 0;
    let mut repeat: LONGLONG = 0;
    let mut startpos: LONGLONG = 0;
    let mut elemnum: LONGLONG = 0;
    let mut wrtptr: LONGLONG = 0;
    let mut rowlen: LONGLONG = 0;
    let mut rownum: LONGLONG = 0;
    let mut remain: LONGLONG = 0;
    let mut next: LONGLONG = 0;
    let mut ntodo: LONGLONG = 0;
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
        if tnull == NULL_UNDEFINED as LONGLONG {
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
                for ii in 0..(ntodo as usize) {
                    ffpbyt(fptr, 1, &[i1null], status);
                }
            }
            TSHORT => {
                for ii in 0..(ntodo as usize) {
                    ffpbyt(fptr, 2, cast_slice(&[i2null]), status);
                }
            }
            TLONG => {
                for ii in 0..(ntodo as usize) {
                    ffpbyt(fptr, 4, cast_slice(&[i4null]), status);
                }
            }
            TLONGLONG => {
                for ii in 0..(ntodo as usize) {
                    ffpbyt(fptr, 8, cast_slice(&[i8null]), status);
                }
            }
            TFLOAT => {
                for ii in 0..(ntodo as usize) {
                    ffpbyt(fptr, 4, cast_slice(&jbuff), status);
                }
            }
            TDOUBLE => {
                for ii in 0..(ntodo as usize) {
                    ffpbyt(fptr, 8, cast_slice(&jbuff), status);
                }
            }
            TLOGICAL => {
                for ii in 0..(ntodo as usize) {
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
    let mut typecode: c_int = 0;
    let mut repeat: LONGLONG = 0;
    let mut width = 0;
    let mut nullstatus: c_int = 0;

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
