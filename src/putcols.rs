/*  This file, putcols.rs, contains routines that write data elements to    */
/*  a FITS image or table, of type character string.                       */

/*  The FITSIO software was written by William Pence at the High Energy    */
/*  Astrophysic Science Archive Research Center (HEASARC) at the NASA      */
/*  Goddard Space Flight Center.                                           */
use core::ffi::CStr;
use core::slice;
use core::{cmp, mem};

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
    let mut tcode: c_int;
    let mut maxelem: c_int = 0;
    let mut hdutype: c_int = 0;
    let nchar: c_int;
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
    let mut message: [c_char; FLEN_ERRMSG] = [0; FLEN_ERRMSG];
    let mut snull: [c_char; 20] = [0; 20]; /*  the FITS null value  */
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

    /* point to first column structure */
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
            LONGLONG::from(nchar),
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
        ffpbyt(fptr, LONGLONG::from(nchar), cast_slice(arr), status);

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
        ntodo = cmp::min(remain, LONGLONG::from(maxelem)) as c_long;
        ntodo = cmp::min(ntodo as LONGLONG, repeat - elemnum) as c_long;

        wrtptr = startpos + (rownum * rowlen) + (elemnum * incre as LONGLONG);
        ffmbyt_safe(fptr, wrtptr, IGNORE_EOF, status); /* move to write position */

        let buffer: &mut [c_char] = cast_slice_mut(&mut cbuff);
        let mut bi = 0;
        /* copy the user's strings into the buffer */
        for _ii in 0..(ntodo as usize) {
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
        if strcmp_safe(nulvalue, arr) != 0 {
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

#[cfg(test)]
mod tests {
    use crate::aliases::rust_api::*;
    use crate::fitsio::{
        ASCII_TBL, BAD_COL_NUM, BAD_ROW_NUM, BINARY_TBL, BYTE_IMG, LONGLONG, NOT_ASCII_COL,
        READONLY, READWRITE, fitsfile,
    };
    use crate::helpers::testhelpers::{from_buf, to_buf, with_temp_file};
    use libc::{c_char, c_int};

    /// Make a NUL-terminated `Vec<c_char>` from a `&str`.
    fn cc(s: &str) -> Vec<c_char> {
        let mut v: Vec<c_char> = s.bytes().map(|b| b as c_char).collect();
        v.push(0);
        v
    }

    /// Create a single-column table (binary or ascii) and return the open file.
    fn make_table(
        name: &[c_char],
        hdutype: c_int,
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
            hdutype,
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
    fn test_write_string_column() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);

            let mut f = make_table(&name, BINARY_TBL, "NAME", "10A", 3, &mut status);
            let data = [cc("Alice"), cc("Bob"), cc("Charlie")];
            let data_ref: Vec<&[c_char]> = data.iter().map(|v| v.as_slice()).collect();
            fits_write_col_str(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                3,
                &data_ref,
                &mut status,
            );
            assert_eq!(status, 0);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            assert_eq!(status, 0);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            assert_eq!(status, 0);

            let mut r1 = [0 as c_char; 11];
            let mut results: [&mut [c_char]; 1] = [&mut r1];
            let mut anynull = -1;
            fits_read_col_str(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                1,
                None,
                &mut results,
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(from_buf(&r1), "Alice");
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_ascii_table_string() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);

            let mut f = make_table(&name, ASCII_TBL, "NAME", "A10", 2, &mut status);
            let data = [cc("Test1"), cc("Test2")];
            let data_ref: Vec<&[c_char]> = data.iter().map(|v| v.as_slice()).collect();
            fits_write_col_str(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                2,
                &data_ref,
                &mut status,
            );
            assert_eq!(status, 0);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            assert_eq!(status, 0);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            assert_eq!(status, 0);

            let mut r1 = [0 as c_char; 11];
            let mut results: [&mut [c_char]; 1] = [&mut r1];
            let mut anynull = -1;
            fits_read_col_str(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                1,
                None,
                &mut results,
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(from_buf(&r1), "Test1");
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_variable_length_string() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);

            let mut f = make_table(&name, BINARY_TBL, "VARSTR", "1PA", 1, &mut status);
            let data = [cc("Variable length string")];
            let data_ref: Vec<&[c_char]> = data.iter().map(|v| v.as_slice()).collect();
            fits_write_col_str(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                1,
                &data_ref,
                &mut status,
            );
            assert_eq!(status, 0);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            assert_eq!(status, 0);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            assert_eq!(status, 0);

            let mut r1 = [0 as c_char; 64];
            let mut results: [&mut [c_char]; 1] = [&mut r1];
            let mut anynull = -1;
            fits_read_col_str(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                1,
                None,
                &mut results,
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(from_buf(&r1), "Variable length string");
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_bad_col_num() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);

            let mut f = make_table(&name, BINARY_TBL, "COL", "10A", 1, &mut status);
            let data = [cc("test")];
            let data_ref: Vec<&[c_char]> = data.iter().map(|v| v.as_slice()).collect();
            fits_write_col_str(
                f.as_deref_mut().unwrap(),
                0,
                1,
                1,
                1,
                &data_ref,
                &mut status,
            );
            assert_eq!(status, BAD_COL_NUM);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_not_ascii_col() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);

            let mut f = make_table(&name, BINARY_TBL, "INTCOL", "1J", 1, &mut status);
            let data = [cc("test")];
            let data_ref: Vec<&[c_char]> = data.iter().map(|v| v.as_slice()).collect();
            fits_write_col_str(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                1,
                &data_ref,
                &mut status,
            );
            assert_eq!(status, NOT_ASCII_COL);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_strings_with_nulls() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);

            let mut f = make_table(&name, BINARY_TBL, "NAME", "10A", 4, &mut status);
            let data = [cc("Good1"), cc("NULL"), cc("Good2"), cc("NULL")];
            let data_ref: Vec<&[c_char]> = data.iter().map(|v| v.as_slice()).collect();
            let nulval = cc("NULL");
            fits_write_colnull_str(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                4,
                &data_ref,
                &nulval,
                &mut status,
            );
            assert_eq!(status, 0);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            assert_eq!(status, 0);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            assert_eq!(status, 0);

            let mut r1 = [0 as c_char; 11];
            let mut results: [&mut [c_char]; 1] = [&mut r1];
            let mut anynull = -1;
            fits_read_col_str(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                1,
                None,
                &mut results,
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(from_buf(&r1), "Good1");
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_nulls_at_end() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);

            let mut f = make_table(&name, BINARY_TBL, "NAME", "10A", 4, &mut status);
            let data = [cc("Good1"), cc("Good2"), cc("NULL"), cc("NULL")];
            let data_ref: Vec<&[c_char]> = data.iter().map(|v| v.as_slice()).collect();
            let nulval = cc("NULL");
            fits_write_colnull_str(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                4,
                &data_ref,
                &nulval,
                &mut status,
            );
            assert_eq!(status, 0);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_all_nulls() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);

            let mut f = make_table(&name, BINARY_TBL, "NAME", "10A", 3, &mut status);
            let data = [cc("NULL"), cc("NULL"), cc("NULL")];
            let data_ref: Vec<&[c_char]> = data.iter().map(|v| v.as_slice()).collect();
            let nulval = cc("NULL");
            fits_write_colnull_str(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                3,
                &data_ref,
                &nulval,
                &mut status,
            );
            assert_eq!(status, 0);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_string_vector() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);

            let mut f = make_table(&name, BINARY_TBL, "NAMES", "30A", 1, &mut status);
            let data = [cc("First"), cc("Second"), cc("Third")];
            let data_ref: Vec<&[c_char]> = data.iter().map(|v| v.as_slice()).collect();
            fits_write_col_str(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                3,
                &data_ref,
                &mut status,
            );
            assert_eq!(status, 0);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            assert_eq!(status, 0);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            assert_eq!(status, 0);

            let mut r1 = [0 as c_char; 11];
            let mut results: [&mut [c_char]; 1] = [&mut r1];
            let mut anynull = -1;
            fits_read_col_str(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                1,
                None,
                &mut results,
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(from_buf(&r1), "First");
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_multirow_strings() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);

            let mut f = make_table(&name, BINARY_TBL, "NAME", "20A", 2, &mut status);
            let data = [
                cc("Row1Col1"),
                cc("Row1Col2"),
                cc("Row2Col1"),
                cc("Row2Col2"),
            ];
            let data_ref: Vec<&[c_char]> = data.iter().map(|v| v.as_slice()).collect();
            fits_write_col_str(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                4,
                &data_ref,
                &mut status,
            );
            assert_eq!(status, 0);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            assert_eq!(status, 0);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            assert_eq!(status, 0);

            let mut r1 = [0 as c_char; 21];
            let mut results: [&mut [c_char]; 1] = [&mut r1];
            let mut anynull = -1;
            fits_read_col_str(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                1,
                None,
                &mut results,
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(from_buf(&r1), "Row1Col1");
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_empty_string() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);

            let mut f = make_table(&name, BINARY_TBL, "NAME", "10A", 3, &mut status);
            let data = [cc(""), cc("NotEmpty"), cc("")];
            let data_ref: Vec<&[c_char]> = data.iter().map(|v| v.as_slice()).collect();
            fits_write_col_str(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                3,
                &data_ref,
                &mut status,
            );
            assert_eq!(status, 0);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            assert_eq!(status, 0);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            assert_eq!(status, 0);

            let mut r1 = [0 as c_char; 11];
            let mut results: [&mut [c_char]; 1] = [&mut r1];
            let mut anynull = -1;
            fits_read_col_str(
                f.as_deref_mut().unwrap(),
                1,
                2,
                1,
                1,
                None,
                &mut results,
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(from_buf(&r1), "NotEmpty");
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    /// Test string columns wider than IOBUFLEN (2880)
    #[test]
    fn test_write_very_wide_string() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);

            // Fill with pattern
            let mut data0 = [0 as c_char; 3001];
            let mut data1 = [0 as c_char; 3001];
            for i in 0..3000 {
                data0[i] = (b'A' + (i % 26) as u8) as c_char;
                data1[i] = (b'a' + (i % 26) as u8) as c_char;
            }
            data0[3000] = 0;
            data1[3000] = 0;

            let mut f = make_table(&name, BINARY_TBL, "WIDE", "3000A", 2, &mut status);
            let data_ref: Vec<&[c_char]> = vec![&data0, &data1];
            fits_write_col_str(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                2,
                &data_ref,
                &mut status,
            );
            assert_eq!(status, 0);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            assert_eq!(status, 0);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            assert_eq!(status, 0);

            let mut r1 = [0 as c_char; 3001];
            let mut results: [&mut [c_char]; 1] = [&mut r1];
            let mut anynull = -1;
            fits_read_col_str(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                1,
                None,
                &mut results,
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            // Just verify first few chars
            assert_eq!(r1[0], b'A' as c_char);
            assert_eq!(r1[1], b'B' as c_char);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    /// Test that functions return early when status > 0
    #[test]
    fn test_calls_with_error_status() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);

            let mut f = make_table(&name, BINARY_TBL, "NAME", "10A", 1, &mut status);
            let data = [cc("test")];
            let data_ref: Vec<&[c_char]> = data.iter().map(|v| v.as_slice()).collect();
            let nulval = cc("");

            // Call ffpcls with pre-existing error
            status = 1;
            fits_write_col_str(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                1,
                &data_ref,
                &mut status,
            );
            assert_eq!(status, 1);

            // Call ffpcns with pre-existing error
            status = 1;
            fits_write_colnull_str(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                1,
                &data_ref,
                &nulval,
                &mut status,
            );
            assert_eq!(status, 1);

            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    /// Test ffgcprll error paths in ffpcls. Tests variable length string
    /// and fixed string ffgcprll error paths
    #[test]
    fn test_ffgcprll_errors() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            let ttype_v = [Some(cc("VARSTR")), Some(cc("FIXSTR"))];
            let ttype_ref: Vec<Option<&[c_char]>> = ttype_v.iter().map(|o| o.as_deref()).collect();
            let tform_v = [cc("1PA"), cc("10A")];
            let tform_ref: Vec<&[c_char]> = tform_v.iter().map(|v| v.as_slice()).collect();
            fits_create_tbl(
                f.as_deref_mut().unwrap(),
                BINARY_TBL,
                1,
                2,
                &ttype_ref,
                &tform_ref,
                None,
                None,
                &mut status,
            );
            assert_eq!(status, 0);

            let data = [cc("test")];
            let data_ref: Vec<&[c_char]> = data.iter().map(|v| v.as_slice()).collect();

            // Bad row number for variable length column
            fits_write_col_str(
                f.as_deref_mut().unwrap(),
                1,
                0,
                1,
                1,
                &data_ref,
                &mut status,
            );
            assert_eq!(status, BAD_ROW_NUM);
            status = 0;

            // Bad row number for fixed length column
            fits_write_col_str(
                f.as_deref_mut().unwrap(),
                2,
                0,
                1,
                1,
                &data_ref,
                &mut status,
            );
            assert_eq!(status, BAD_ROW_NUM);
            status = 0;

            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    /// Test error paths within ffpcns when ffpclu/ffpcls fail
    #[test]
    fn test_ffpcns_internal_errors() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);

            let mut f = make_table(&name, BINARY_TBL, "NAME", "10A", 2, &mut status);
            let data1 = [cc("NULL"), cc("good")];
            let data1_ref: Vec<&[c_char]> = data1.iter().map(|v| v.as_slice()).collect();
            let data2 = [cc("good"), cc("NULL")];
            let data2_ref: Vec<&[c_char]> = data2.iter().map(|v| v.as_slice()).collect();
            let nulval = cc("NULL");

            // "NULL" then "good" - should fail when writing nulls
            fits_write_colnull_str(
                f.as_deref_mut().unwrap(),
                1,
                0,
                1,
                2,
                &data1_ref,
                &nulval,
                &mut status,
            );
            assert_ne!(status, 0);
            status = 0;

            // "good" then "NULL" - should fail when writing good values
            fits_write_colnull_str(
                f.as_deref_mut().unwrap(),
                1,
                0,
                1,
                2,
                &data2_ref,
                &nulval,
                &mut status,
            );
            assert_ne!(status, 0);
            status = 0;

            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    /// Test HDU position mismatch recovery in ffpcls and ffpcns. Create a
    /// file with two tables, then manipulate HDUposition to trigger the
    /// mismatch correction code
    #[test]
    fn test_hdu_position_mismatch() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);

            let data = [cc("test")];
            let data_ref: Vec<&[c_char]> = data.iter().map(|v| v.as_slice()).collect();
            let nulval = cc("X");

            // Create file with primary + 2 tables
            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            let ttype_v = [Some(cc("NAME"))];
            let ttype_ref: Vec<Option<&[c_char]>> = ttype_v.iter().map(|o| o.as_deref()).collect();
            let tform_v = [cc("10A")];
            let tform_ref: Vec<&[c_char]> = tform_v.iter().map(|v| v.as_slice()).collect();
            fits_create_tbl(
                f.as_deref_mut().unwrap(),
                BINARY_TBL,
                1,
                1,
                &ttype_ref,
                &tform_ref,
                None,
                Some(&cc("TABLE1")),
                &mut status,
            );
            fits_create_tbl(
                f.as_deref_mut().unwrap(),
                BINARY_TBL,
                1,
                1,
                &ttype_ref,
                &tform_ref,
                None,
                Some(&cc("TABLE2")),
                &mut status,
            );
            assert_eq!(status, 0);
            fits_close_file(f.take().unwrap(), &mut status);

            // Reopen and move to HDU 3 (TABLE2)
            fits_open_file(&mut f, &name, READWRITE, &mut status);
            assert_eq!(status, 0);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 3, None, &mut status);
            assert_eq!(status, 0);

            // Set HDUposition to 1 (TABLE1) to trigger mismatch
            // curhdu=2 (0-indexed, TABLE2), HDUposition=1 (TABLE1)
            f.as_deref_mut().unwrap().HDUposition = 1;

            // ffpcls should call ffmahd to move to HDU 2 (TABLE1)
            fits_write_col_str(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                1,
                &data_ref,
                &mut status,
            );
            assert_eq!(status, 0);

            // Now move back to HDU 3 and set mismatch for ffpcns
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 3, None, &mut status);
            assert_eq!(status, 0);
            f.as_deref_mut().unwrap().HDUposition = 1;
            fits_write_colnull_str(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                1,
                &data_ref,
                &nulval,
                &mut status,
            );
            assert_eq!(status, 0);

            fits_close_file(f.take().unwrap(), &mut status);
        });
    }
}
