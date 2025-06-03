/*  This file, putcols.c, contains routines that write data elements to    */
/*  a FITS image or table, of type character string.                       */

/*  The FITSIO software was written by William Pence at the High Energy    */
/*  Astrophysic Science Archive Research Center (HEASARC) at the NASA      */
/*  Goddard Space Flight Center.                                           */
use core::slice;
use std::ffi::CStr;
use std::{cmp, mem, ptr};

use crate::c_types::*;

use bytemuck::{cast_slice, cast_slice_mut};

use crate::bb;
use crate::fitscore::{ffgcprll, ffgtcl_safe, ffmahd_safe, ffpmsg_slice, ffrdef_safe};
use crate::fitsio::*;
use crate::fitsio2::*;
use crate::putcolu::ffpclu_safe;
use crate::wrappers::*;
use crate::{buffers::*, int_snprintf};

/*--------------------------------------------------------------------------*/
/// Write an array of string values to a column in the current FITS HDU.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffpcls(
    fptr: *mut fitsfile,         /* I - FITS file pointer                       */
    colnum: c_int,               /* I - number of column to write (1 = 1st col) */
    firstrow: LONGLONG,          /* I - first row to write (1 = 1st row)        */
    firstelem: LONGLONG,         /* I - first vector element to write (1 = 1st) */
    nelem: LONGLONG,             /* I - number of strings to write              */
    array: *const *const c_char, /* I - array of pointers to strings            */
    status: *mut c_int,          /* IO - error status                           */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        let array = slice::from_raw_parts(array, nelem as usize);
        let mut v_array = Vec::new();

        for item in array {
            let tform_item = cast_slice(CStr::from_ptr(*item).to_bytes_with_nul());
            v_array.push(tform_item);
        }

        ffpcls_safe(
            fptr,
            colnum,
            firstrow as LONGLONG,
            firstelem as LONGLONG,
            nelem as LONGLONG,
            &v_array,
            status,
        )
    }
}

/*--------------------------------------------------------------------------*/
/// Write an array of string values to a column in the current FITS HDU.
pub fn ffpcls_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer                       */
    colnum: c_int,       /* I - number of column to write (1 = 1st col) */
    firstrow: LONGLONG,  /* I - first row to write (1 = 1st row)        */
    firstelem: LONGLONG, /* I - first vector element to write (1 = 1st) */
    nelem: LONGLONG,     /* I - number of strings to write              */
    array: &[&[c_char]], /* I - array of pointers to strings            */
    status: &mut c_int,  /* IO - error status                           */
) -> c_int {
    let mut tcode: c_int = 0;
    let mut maxelem: c_int = 0;
    let mut hdutype: c_int = 0;
    let mut nchar: c_int = 0;
    let mut twidth: c_long = 0;
    let mut incre: c_long = 0;
    let ii: c_long = 0;
    let jj: c_long = 0;
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
    let mut scale: f64 = 0.0;
    let mut zero: f64 = 0.0;
    let mut tform: [c_char; 20] = [0; 20];
    let mut blanks: &[c_char];
    let mut message: [c_char; FLEN_ERRMSG] = [0; FLEN_ERRMSG];
    let mut snull: [c_char; 20] = [0; 20]; /*  the FITS null value  */
    let colptr: *mut tcolumn = ptr::null_mut();
    let mut cbuff: [f64; DBUFFSIZE as usize / mem::size_of::<f64>()] =
        [0.0; DBUFFSIZE as usize / mem::size_of::<f64>()]; /* align cbuff on word boundary */

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    /* reset position to the correct HDU if necessary */
    if fptr.HDUposition != fptr.Fptr.curhdu {
        ffmahd_safe(fptr, (fptr.HDUposition) + 1, None, status);
    } else if fptr.Fptr.datastart == DATA_UNDEFINED as LONGLONG && ffrdef_safe(fptr, status) > 0 {
        /* rescan header */
        return *status;
    }

    /*---------------------------------------------------*/
    /*  Check input and get parameters about the column: */
    /*---------------------------------------------------*/
    if colnum < 1 || colnum > fptr.Fptr.tfield {
        int_snprintf!(
            &mut message,
            FLEN_ERRMSG,
            "Specified column number is out of range: {}",
            colnum,
        );
        ffpmsg_slice(&message);
        *status = BAD_COL_NUM;
        return *status;
    }

    let colptr = fptr.Fptr.tableptr; /* point to first column structure */
    let c = fptr.Fptr.get_tableptr_as_slice();
    let ci = (colnum - 1) as usize; /* offset to the correct column */

    tcode = c[ci].tdatatype;

    if tcode == -TSTRING {
        /* variable length column in a binary table? */

        /* only write a single string; ignore value of firstelem */
        nchar = cmp::max(1, strlen_safe(array[0])) as c_int; /* will write at least 1 char */
        /* even if input string is null */

        if ffgcprll(
            fptr,
            colnum,
            firstrow,
            1,
            nchar as LONGLONG,
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

        /* simply move to write position, then write the string */
        let arr = array[0];
        ffmbyt_safe(fptr, startpos, IGNORE_EOF, status);
        ffpbyt(fptr, nchar as LONGLONG, cast_slice(arr), status);

        if *status > 0 {
            /* test for error during previous write operation */

            int_snprintf!(
                &mut message,
                FLEN_ERRMSG,
                "Error writing to variable length string column (ffpcls).",
            );
            ffpmsg_slice(&message);
        }

        return *status;
    } else if tcode == TSTRING {
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
        /* if string length is greater than a FITS block (2880 char) then must */
        /* only write 1 string at a time, to force writein by ffpbyt instead of */
        /* ffpbytoff (ffpbytoff can't handle this case) */
        if (twidth as LONGLONG) > IOBUFLEN {
            maxelem = 1;
            incre = twidth;
            repeat = 1;
        }

        let mut v = vec![0; twidth as usize];
        let blanks = &mut v; /* string for blank fill values */

        /* if (!blanks) {
            ffpmsg_str("Could not allocate memory for string (ffpcls)");
            *status = ARRAY_TOO_BIG;
            return *status;
        } */

        for ii in 0..(twidth as usize) {
            blanks[ii] = bb(b' '); /* fill string with blanks */
        }
        remain = nelem; /* remaining number of values to write  */
    } else {
        *status = NOT_ASCII_COL;
        return *status;
    }
    /*-------------------------------------------------------*/
    /*  Now write the strings to the FITS column.            */
    /*-------------------------------------------------------*/

    next = 0; /* next element in array to be written  */
    rownum = 0; /* row number, relative to firstrow     */

    while remain > 0 {
        /* limit the number of pixels to process at one time to the number that
           will fit in the buffer space or to the number of pixels that remain
           in the current vector, which ever is smaller.
        */
        ntodo = cmp::min(remain, maxelem as LONGLONG) as c_long;
        ntodo = cmp::min(ntodo as LONGLONG, repeat - elemnum) as c_long;

        wrtptr = startpos + (rownum * rowlen) + (elemnum * incre as LONGLONG);
        ffmbyt_safe(fptr, wrtptr, IGNORE_EOF, status); /* move to write position */

        let buffer: &mut [c_char] = cast_slice_mut(&mut cbuff);
        let mut bi = 0;
        /* copy the user's strings into the buffer */
        for ii in 0..(ntodo as usize) {
            let arrayptr = array[next as usize];
            let mut ai = 0;

            let mut jj = 0;
            while jj < twidth {
                /*  copy the string, char by char */

                if arrayptr[ai] > 0 {
                    buffer[bi] = arrayptr[ai];
                    bi += 1;
                    ai += 1;
                } else {
                    break;
                }
                jj += 1;
            }

            while jj < twidth {
                /* fill field with blanks, if needed */

                buffer[bi] = bb(b' ');
                bi += 1;
                jj += 1;
            }

            next += 1;
        }

        /* write the buffer full of strings to the FITS file */
        if incre == twidth {
            ffpbyt(
                fptr,
                (ntodo * twidth) as LONGLONG,
                cast_slice_mut(&mut cbuff),
                status,
            );
        } else {
            ffpbytoff(
                fptr,
                twidth,
                ntodo,
                incre - twidth,
                cast_slice(&cbuff),
                status,
            );
        }
        if *status > 0 {
            /* test for error during previous write operation */

            int_snprintf!(
                &mut message,
                FLEN_ERRMSG,
                "Error writing elements {:.0} thru {:.0} of input data array (ffpcls).",
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
            elemnum += ntodo as LONGLONG;
            if elemnum == repeat {
                /* completed a row; start on next row */

                elemnum = 0;
                rownum += 1;
            }
        }
    } /*  End of main while Loop  */

    *status
}

/*--------------------------------------------------------------------------*/
///  Write an array of elements to the specified column of a table.  Any input
///  pixels flagged as null will be replaced by the appropriate
///  null value in the output FITS file.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffpcns(
    fptr: *mut fitsfile,         /* I - FITS file pointer                       */
    colnum: c_int,               /* I - number of column to write (1 = 1st col) */
    firstrow: LONGLONG,          /* I - first row to write (1 = 1st row)        */
    firstelem: LONGLONG,         /* I - first vector element to write (1 = 1st) */
    nelem: LONGLONG,             /* I - number of values to write               */
    array: *const *const c_char, /* I - array of values to write                */
    nulvalue: *const c_char,     /* I - string representing a null value        */
    status: *mut c_int,          /* IO - error status                           */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        let array = slice::from_raw_parts(array, nelem as usize);
        let mut v_array = Vec::new();

        for item in array {
            let tform_item = cast_slice(CStr::from_ptr(*item).to_bytes_with_nul());
            v_array.push(tform_item);
        }

        let nulvalue: &[c_char] = cast_slice(CStr::from_ptr(nulvalue).to_bytes_with_nul());

        ffpcns_safe(
            fptr,
            colnum,
            firstrow as LONGLONG,
            firstelem as LONGLONG,
            nelem as LONGLONG,
            &v_array,
            nulvalue,
            status,
        )
    }
}

/*--------------------------------------------------------------------------*/
///  Write an array of elements to the specified column of a table.  Any input
///  pixels flagged as null will be replaced by the appropriate
///  null value in the output FITS file.
pub fn ffpcns_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer                       */
    colnum: c_int,       /* I - number of column to write (1 = 1st col) */
    firstrow: LONGLONG,  /* I - first row to write (1 = 1st row)        */
    firstelem: LONGLONG, /* I - first vector element to write (1 = 1st) */
    nelem: LONGLONG,     /* I - number of values to write               */
    array: &[&[c_char]], /* I - array of values to write                */
    nulvalue: &[c_char], /* I - string representing a null value        */
    status: &mut c_int,  /* IO - error status                           */
) -> c_int {
    let mut repeat: c_long = 0;
    let mut width: c_long = 0;
    let mut ngood: LONGLONG = 0;
    let mut nbad: LONGLONG = 0;
    let mut first: LONGLONG = 0;
    let mut fstelm: LONGLONG = 0;
    let mut fstrow: LONGLONG = 0;

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

    /* get the vector repeat length of the column */
    ffgtcl_safe(
        fptr,
        colnum,
        None,
        Some(&mut repeat),
        Some(&mut width),
        status,
    );

    if fptr.Fptr.hdutype == BINARY_TBL {
        repeat /= width; /* convert from chars to unit strings */
    }

    /* absolute element number in the column */
    first = (firstrow - 1) * repeat as LONGLONG + firstelem;

    let mut ii: LONGLONG = 0;
    while ii < nelem {
        let arr = array[ii as usize];
        if strcmp_safe(nulvalue, arr) > 0 {
            /* is this a good pixel? */

            if nbad != 0 {
                /* write previous string of bad pixels */

                fstelm = ii - nbad + first; /* absolute element number */
                fstrow = (fstelm - 1) / repeat as LONGLONG + 1; /* starting row number */
                fstelm -= (fstrow - 1) * repeat as LONGLONG; /* relative number */

                if ffpclu_safe(fptr, colnum, fstrow, fstelm, nbad, status) > 0 {
                    return *status;
                }
                nbad = 0;
            }

            ngood += 1; /* the consecutive number of good pixels */
        } else {
            if ngood > 0 {
                /* write previous string of good pixels */

                fstelm = ii - ngood + first; /* absolute element number */
                fstrow = (fstelm - 1) / repeat as LONGLONG + 1; /* starting row number */
                fstelm -= (fstrow - 1) * repeat as LONGLONG; /* relative number */

                if ffpcls_safe(
                    fptr,
                    colnum,
                    fstrow,
                    fstelm,
                    ngood,
                    &array[(ii - ngood) as usize..],
                    status,
                ) > 0
                {
                    return *status;
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
        fstelm = ii as LONGLONG - ngood + first; /* absolute element number */
        fstrow = (fstelm - 1) / repeat as LONGLONG + 1; /* starting row number */
        fstelm -= (fstrow - 1) * repeat as LONGLONG; /* relative number */

        ffpcls_safe(
            fptr,
            colnum,
            fstrow,
            fstelm,
            ngood,
            &array[(ii - ngood) as usize..],
            status,
        );
    } else if nbad != 0 {
        /* write last string of bad pixels */

        fstelm = ii - nbad + first; /* absolute element number */
        fstrow = (fstelm - 1) / repeat as LONGLONG + 1; /* starting row number */
        fstelm -= (fstrow - 1) * repeat as LONGLONG; /* relative number */

        ffpclu_safe(fptr, colnum, fstrow, fstelm, nbad, status);
    }

    *status
}
