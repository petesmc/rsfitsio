//! Routines that write data elements to a FITS image or table, with logical
//! datatype.
//!
//! The [`TLOGICAL`] arm of the typed column I/O family. [`crate::putcol`]
//! dispatches here when a caller asks for this datatype at run time; the
//! write-side counterpart is [`crate::getcoll`].
//!
//! Ported from CFITSIO's `putcoll.c`, written by William Pence at the High
//! Energy Astrophysics Science Archive Research Center (HEASARC), NASA Goddard
//! Space Flight Center.
#![warn(missing_docs)]

use core::slice;

use crate::c_types::*;

use bytemuck::cast_slice_mut;

use crate::bb;
use crate::fitscore::{
    ffgcprll, ffgdesll_safe, ffmahd_safe, ffpdes_safe, ffpmsg_slice, ffrdef_safe,
};
use crate::fitsio::*;
use crate::fitsio2::*;
use crate::putcolu::ffpclu_safe;
use crate::{buffers::*, int_snprintf};

/// Write an array of logical values to a column in the current FITS HDU.
///
/// # Parameters
///
/// * `fptr`      — (I) FITS file pointer
/// * `colnum`    — (I) number of column to write (1 = 1st col)
/// * `firstrow`  — (I) first row to write (1 = 1st row)
/// * `firstelem` — (I) first vector element to write (1 = 1st)
/// * `nelem`     — (I) number of values to write
/// * `array`     — (I) array of values to write
/// * `status`    — (IO) error status
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffpcll(
    fptr: *mut fitsfile,
    colnum: c_int,
    firstrow: LONGLONG,
    firstelem: LONGLONG,
    nelem: LONGLONG,
    array: *const c_char,
    status: *mut c_int,
) -> c_int {
    // FFI WRAPPER
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        let array = slice::from_raw_parts(array, nelem as usize);

        ffpcll_safe(
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

/// Write an array of logical values to a column in the current FITS HDU.
///
/// # Parameters
///
/// * `fptr`      — (I) FITS file pointer
/// * `colnum`    — (I) number of column to write (1 = 1st col)
/// * `firstrow`  — (I) first row to write (1 = 1st row)
/// * `firstelem` — (I) first vector element to write (1 = 1st)
/// * `nelem`     — (I) number of values to write
/// * `array`     — (I) array of values to write
/// * `status`    — (IO) error status
pub fn ffpcll_safe(
    fptr: &mut fitsfile,
    colnum: c_int,
    firstrow: LONGLONG,
    firstelem: LONGLONG,
    nelem: LONGLONG,
    array: &[c_char],
    status: &mut c_int,
) -> c_int {
    let mut tcode: c_int = 0;
    let mut maxelem: c_int = 0;
    let mut hdutype: c_int = 0;
    let mut twidth: c_long = 0;
    let mut incre: c_long = 0;
    let mut repeat: LONGLONG = 0;
    let mut startpos: LONGLONG = 0;
    let mut elemnum: LONGLONG = 0;
    let mut wrtptr: LONGLONG = 0;
    let mut rowlen: LONGLONG = 0;
    let mut rownum: LONGLONG = 0;
    let mut remain: LONGLONG = 0;
    let mut next = 0;
    let mut tnull: LONGLONG = 0;
    let mut scale: f64 = 0.0;
    let mut zero: f64 = 0.0;
    let mut tform: [c_char; 20] = [0; 20];
    let ctrue: c_char = bb(b'T');
    let cfalse: c_char = bb(b'F');
    let mut message: [c_char; FLEN_ERRMSG] = [0; FLEN_ERRMSG];
    let mut snull: [c_char; 20] = [0; 20]; /*  the FITS null value  */

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    /*  Check input and get parameters about the column: */
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

    if tcode != TLOGICAL {
        *status = NOT_LOGICAL_COL;
        return *status;
    }

    /*  Now write the logical values one at a time to the FITS column.     */
    remain = nelem; /* remaining number of values to write  */
    next = 0; /* next element in array to be written  */
    rownum = 0; /* row number, relative to firstrow     */

    while remain != 0 {
        wrtptr = startpos + (rowlen * rownum) + (elemnum * incre as LONGLONG);

        ffmbyt_safe(fptr, wrtptr, IGNORE_EOF, status); /* move to write position */

        if array[next] != 0 {
            ffpbyt(fptr, 1, cast_slice_mut(&mut [ctrue]), status);
        } else {
            ffpbyt(fptr, 1, cast_slice_mut(&mut [cfalse]), status);
        }

        if *status > 0 {
            /* test for error during previous write operation */

            int_snprintf!(
                &mut message,
                FLEN_ERRMSG,
                "Error writing element {:.0} of input array of logicals (ffpcll).",
                (next + 1) as f64,
            );
            ffpmsg_slice(&message);
            return *status;
        }

        /*  increment the counters for the next loop  */
        remain -= 1;
        if remain != 0 {
            next += 1;
            elemnum += 1;
            if elemnum == repeat {
                /* completed a row; start on next row */

                elemnum = 0;
                rownum += 1;
            }
        }
    } /*  End of main while Loop  */

    *status
}

/// Write an array of elements to the specified column of a table.  Any input
/// pixels flagged as null will be replaced by the appropriate
/// null value in the output FITS file.
///
/// # Parameters
///
/// * `fptr`      — (I) FITS file pointer
/// * `colnum`    — (I) number of column to write (1 = 1st col)
/// * `firstrow`  — (I) first row to write (1 = 1st row)
/// * `firstelem` — (I) first vector element to write (1 = 1st)
/// * `nelem`     — (I) number of values to write
/// * `array`     — (I) array of values to write
/// * `nulvalue`  — (I) array flagging undefined pixels if true
/// * `status`    — (IO) error status
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffpcnl(
    fptr: *mut fitsfile,
    colnum: c_int,
    firstrow: LONGLONG,
    firstelem: LONGLONG,
    nelem: LONGLONG,
    array: *const c_char,
    nulvalue: c_char,
    status: *mut c_int,
) -> c_int {
    // FFI WRAPPER
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        let array = slice::from_raw_parts(array, nelem as usize);

        ffpcnl_safe(
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

/// Write an array of elements to the specified column of a table.  Any input
/// pixels flagged as null will be replaced by the appropriate
/// null value in the output FITS file.
///
/// # Parameters
///
/// * `fptr`      — (I) FITS file pointer
/// * `colnum`    — (I) number of column to write (1 = 1st col)
/// * `firstrow`  — (I) first row to write (1 = 1st row)
/// * `firstelem` — (I) first vector element to write (1 = 1st)
/// * `nelem`     — (I) number of values to write
/// * `array`     — (I) array of values to write
/// * `nulvalue`  — (I) array flagging undefined pixels if true
/// * `status`    — (IO) error status
pub fn ffpcnl_safe(
    fptr: &mut fitsfile,
    colnum: c_int,
    firstrow: LONGLONG,
    firstelem: LONGLONG,
    nelem: LONGLONG,
    array: &[c_char],
    nulvalue: c_char,
    status: &mut c_int,
) -> c_int {
    let mut ngood: LONGLONG = 0;
    let mut nbad: LONGLONG = 0;
    let mut ii: LONGLONG = 0;
    let mut repeat: LONGLONG = 0;
    let mut first: LONGLONG = 0;
    let mut fstelm: LONGLONG = 0;
    let mut fstrow: LONGLONG = 0;
    let mut tcode: c_int = 0;

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
    let ci = colnum as usize - 1; /* increment to the correct column */

    tcode = c[ci].tdatatype;

    if tcode > 0 {
        repeat = c[ci].trepeat; /* repeat count for this column */
    } else {
        repeat = firstelem - 1 + nelem; /* variable length arrays */
    }

    /* first write the whole input vector, then go back and fill in the nulls */
    if ffpcll_safe(
        fptr,
        colnum,
        firstrow as LONGLONG,
        firstelem as LONGLONG,
        nelem as LONGLONG,
        array,
        status,
    ) > 0
    {
        return *status;
    }

    /* absolute element number in the column */
    first = (firstrow - 1) * repeat + firstelem;

    ii = 0;
    while ii < nelem {
        if array[ii as usize] != nulvalue {
            /* is this a good pixel? */

            if nbad != 0 {
                /* write previous string of bad pixels */

                fstelm = ii - nbad + first; /* absolute element number */
                fstrow = (fstelm - 1) / repeat + 1; /* starting row number */
                fstelm -= (fstrow - 1) * repeat; /* relative number */

                if ffpclu_safe(fptr, colnum, fstrow, fstelm, nbad, status) > 0 {
                    return *status;
                }

                nbad = 0;
            }

            ngood += 1; /* the consecutive number of good pixels */
        } else {
            if ngood != 0
            /* write previous string of good pixels */
            {
                fstelm = ii - ngood + first; /* absolute element number */
                fstrow = (fstelm - 1) / repeat + 1; /* starting row number */
                fstelm -= (fstrow - 1) * repeat; /* relative number */

                /*  good values have already been written
                            if (ffpcll(fptr, colnum, fstrow, fstelm, ngood, &array[ii-ngood],
                                status) > 0)
                                return(*status);
                */
                ngood = 0;
            }

            nbad += 1; /* the consecutive number of bad pixels */
        }
        ii += 1;
    }

    /* finished loop;  now just write the last set of pixels */

    if ngood != 0 {
        /* write last string of good pixels */

        fstelm = ii - ngood + first; /* absolute element number */
        fstrow = (fstelm - 1) / repeat + 1; /* starting row number */
        fstelm -= (fstrow - 1) * repeat; /* relative number */

    /*  these have already been written
          ffpcll(fptr, colnum, fstrow, fstelm, ngood, &array[ii-ngood], status);
    */
    } else if nbad != 0 {
        /* write last string of bad pixels */

        fstelm = ii - nbad + first; /* absolute element number */
        fstrow = (fstelm - 1) / repeat + 1; /* starting row number */
        fstelm -= (fstrow - 1) * repeat; /* relative number */

        ffpclu_safe(fptr, colnum, fstrow, fstelm, nbad, status);
    }

    *status
}

/// write an array of logical values to a specified bit or byte
/// column of the binary table.   If larray is TRUE, then the corresponding
/// bit is set to 1, otherwise it is set to 0.
/// The binary table column being written to must have datatype 'B' or 'X'.
///
/// # Parameters
///
/// * `fptr`   — (I) FITS file pointer
/// * `colnum` — (I) number of column to write (1 = 1st col)
/// * `frow`   — (I) first row to write (1 = 1st row)
/// * `fbit`   — (I) first bit to write (1 = 1st)
/// * `nbit`   — (I) number of bits to write
/// * `larray` — (I) array of logicals corresponding to bits
/// * `status` — (IO) error status
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffpclx(
    fptr: *mut fitsfile,
    colnum: c_int,
    frow: LONGLONG,
    fbit: c_long,
    nbit: c_long,
    larray: *const c_char,
    status: *mut c_int,
) -> c_int {
    // FFI WRAPPER
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        let larray = slice::from_raw_parts(larray, nbit as usize);

        ffpclx_safe(fptr, colnum, frow, fbit, nbit, larray, status)
    }
}

/// write an array of logical values to a specified bit or byte
/// column of the binary table.   If larray is TRUE, then the corresponding
/// bit is set to 1, otherwise it is set to 0.
/// The binary table column being written to must have datatype 'B' or 'X'.
///
/// # Parameters
///
/// * `fptr`   — (I) FITS file pointer
/// * `colnum` — (I) number of column to write (1 = 1st col)
/// * `frow`   — (I) first row to write (1 = 1st row)
/// * `fbit`   — (I) first bit to write (1 = 1st)
/// * `nbit`   — (I) number of bits to write
/// * `larray` — (I) array of logicals corresponding to bits
/// * `status` — (IO) error status
pub fn ffpclx_safe(
    fptr: &mut fitsfile,
    colnum: c_int,
    frow: LONGLONG,
    fbit: c_long,
    nbit: c_long,
    larray: &[c_char],
    status: &mut c_int,
) -> c_int {
    let mut bstart: LONGLONG = 0;
    let mut repeat: LONGLONG = 0;
    let mut rowlen: LONGLONG = 0;
    let mut elemnum: LONGLONG = 0;
    let mut rstart: LONGLONG = 0;
    let mut estart: LONGLONG = 0;
    let mut tnull: LONGLONG = 0;
    let mut fbyte: c_long = 0;
    let mut lbyte: c_long = 0;
    let mut nbyte: c_long = 0;
    let mut bitloc: c_long = 0;
    let mut ndone: c_long = 0;
    let mut twidth: c_long = 0;
    let mut incre: c_long = 0;
    let mut tcode: c_int = 0;
    let mut descrp: c_int = 0;
    let mut maxelem: c_int = 0;
    let mut hdutype: c_int = 0;
    let mut dummyd: f64 = 0.0;
    let mut dummyd2: f64 = 0.0;
    let mut tform: [c_char; 12] = [0; 12];
    let mut snull: [c_char; 12] = [0; 12];
    let mut cbuff: [u8; 1] = [0];
    static ONBIT: [u8; 8] = [128, 64, 32, 16, 8, 4, 2, 1];
    static OFFBIT: [u8; 8] = [127, 191, 223, 239, 247, 251, 253, 254];
    let mut heapoffset: LONGLONG = 0;
    let mut lrepeat: LONGLONG = 0;

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    /*  check input parameters */
    if nbit < 1 {
        return *status;
    } else if frow < 1 {
        *status = BAD_ROW_NUM;
        return *status;
    } else if fbit < 1 {
        *status = BAD_ELEM_NUM;
        return *status;
    }

    /* reset position to the correct HDU if necessary */
    if fptr.HDUposition != fptr.Fptr.curhdu {
        ffmahd_safe(fptr, (fptr.HDUposition) + 1, None, status);

    /* rescan header if data structure is undefined */
    } else if fptr.Fptr.datastart == DATA_UNDEFINED as LONGLONG && ffrdef_safe(fptr, status) > 0 {
        return *status;
    }

    fbyte = (fbit + 7) / 8;
    lbyte = (fbit + nbit + 6) / 8;
    nbyte = lbyte - fbyte + 1;

    /* Save the current heapsize; ffgcprll will increment the value if */
    /* we are writing to a variable length column. */
    // offset = fptr.Fptr.heapsize;

    /* call ffgcprll in case we are writing beyond the current end of   */
    /* the table; it will allocate more space and shift any following */
    /* HDU's.  Otherwise, we have little use for most of the returned */
    /* parameters, therefore just use dummy parameters.               */

    if ffgcprll(
        fptr,
        colnum,
        frow,
        fbyte as LONGLONG,
        nbyte as LONGLONG,
        1,
        &mut dummyd,
        &mut dummyd2,
        &mut tform,
        &mut twidth,
        &mut tcode,
        &mut maxelem,
        &mut bstart,
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

    bitloc = fbit - 1 - ((fbit - 1) / 8 * 8);
    ndone = 0;
    rstart = frow - 1;
    estart = fbyte as LONGLONG - 1;

    /* set pointer to first column */
    let c = fptr.Fptr.get_tableptr_as_slice();
    let ci = colnum as usize - 1; /* increment to the correct column */

    tcode = c[ci].tdatatype;

    if (tcode.abs()) > TBYTE {
        *status = NOT_LOGICAL_COL;
        return *status;
    } /* not correct datatype column */

    if tcode > 0 {
        descrp = FALSE as c_int; /* not a variable length descriptor column */
        repeat = c[ci].trepeat;

        if tcode == TBIT {
            repeat = (repeat + 7) / 8; /* convert from bits to bytes */
        }

        if fbyte as LONGLONG > repeat {
            *status = BAD_ELEM_NUM;
            return *status;
        }

        /* calc the i/o pointer location to start of sequence of pixels */
        bstart = fptr.Fptr.datastart + (fptr.Fptr.rowlength * rstart) + c[ci].tbcol + estart;
    } else {
        descrp = TRUE as c_int; /* a variable length descriptor column */
        /* only bit arrays (tform = 'X') are supported for variable */
        /* length arrays.  REPEAT is the number of BITS in the array. */

        repeat = (fbit + nbit - 1) as LONGLONG;

        /* write the number of elements and the starting offset.    */
        /* Note: ffgcprll previous wrote the descripter, but with the */
        /* wrong repeat value  (gave bytes instead of bits).        */
        /* Make sure to not change the current heap offset value!  */

        if tcode == -TBIT {
            ffgdesll_safe(
                fptr,
                colnum,
                frow,
                Some(&mut lrepeat),
                Some(&mut heapoffset),
                status,
            );
            ffpdes_safe(fptr, colnum, frow, repeat as LONGLONG, heapoffset, status);
        }

        /* Calc the i/o pointer location to start of sequence of pixels.   */
        /* ffgcprll has already calculated a value for bstart that         */
        /* points to the first element of the vector; we just have to      */
        /* increment it to point to the first element we want to write to. */
        /* Note: ffgcprll also already updated the size of the heap, so we */
        /* don't have to do that again here.                               */

        bstart += estart;
    }

    /* move the i/o pointer to the start of the pixel sequence */
    ffmbyt_safe(fptr, bstart, IGNORE_EOF, status);

    /* read the next byte (we may only be modifying some of the bits) */
    loop {
        if ffgbyt(fptr, 1, cast_slice_mut(&mut cbuff), status) == END_OF_FILE {
            /* hit end of file trying to read the byte, so just set byte = 0 */
            *status = 0;
            cbuff[0] = 0;
        }

        /* move back, to be able to overwrite the byte */
        ffmbyt_safe(fptr, bstart, IGNORE_EOF, status);

        let mut ii = bitloc as usize;
        while (ii < 8) && (ndone < nbit) {
            if larray[ndone as usize] != 0 {
                cbuff[0] |= ONBIT[ii];
            } else {
                cbuff[0] &= OFFBIT[ii];
            }
            ii += 1;
            ndone += 1;
        }

        ffpbyt(fptr, 1, cast_slice_mut(&mut cbuff), status); /* write the modified byte */
        if ndone == nbit {
            /* finished all the bits */
            return *status;
        }

        /* not done, so get the next byte */
        bstart += 1;
        if descrp == 0 {
            estart += 1;
            if estart == repeat {
                let c = fptr.Fptr.get_tableptr_as_slice();

                /* move the i/o pointer to the next row of pixels */
                estart = 0;
                rstart += 1;
                bstart = fptr.Fptr.datastart + (fptr.Fptr.rowlength * rstart) + c[ci].tbcol;

                ffmbyt_safe(fptr, bstart, IGNORE_EOF, status);
            }
        }
        bitloc = 0;
    }
}

#[cfg(test)]
mod tests {
    use crate::aliases::rust_api::*;
    use crate::fitsio::{
        BAD_ELEM_NUM, BAD_ROW_NUM, BINARY_TBL, BYTE_IMG, LONGLONG, NOT_LOGICAL_COL, READONLY,
        fitsfile,
    };
    use crate::helpers::testhelpers::{to_buf, with_temp_file};
    use libc::{c_char, c_int};

    /// Make a NUL-terminated `Vec<c_char>` from a `&str`.
    fn cc(s: &str) -> Vec<c_char> {
        let mut v: Vec<c_char> = s.bytes().map(|b| b as c_char).collect();
        v.push(0);
        v
    }

    /// Create a single-column binary table and return the open file.
    fn make_table(
        name: &[c_char],
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
            BINARY_TBL,
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
    fn test_write_logical_column() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [c_char; 5] = [1, 0, 1, 1, 0];

            let mut f = make_table(&name, "LOGCOL", "1L", 5, &mut status);
            fits_write_col_log(f.as_deref_mut().unwrap(), 1, 1, 1, 5, &data, &mut status);
            assert_eq!(status, 0);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            assert_eq!(status, 0);
            let mut result = [0 as c_char; 5];
            let mut anynull = 0;
            fits_read_col_log(
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
            assert!(!(result[0] != 1));
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_not_logical_column() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [c_char; 1] = [1];

            let mut f = make_table(&name, "INTCOL", "1J", 1, &mut status);
            assert_eq!(status, 0);
            fits_write_col_log(f.as_deref_mut().unwrap(), 1, 1, 1, 1, &data, &mut status);
            assert!(!(status != NOT_LOGICAL_COL));
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_logical_with_nulls() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [c_char; 5] = [1, 99, 0, 99, 1];
            let nulval: c_char = 99;

            let mut f = make_table(&name, "LOGCOL", "1L", 5, &mut status);
            fits_write_colnull_log(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                5,
                &data,
                nulval,
                &mut status,
            );
            assert_eq!(status, 0);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_nulls_at_end() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [c_char; 4] = [1, 0, 99, 99];
            let nulval: c_char = 99;

            let mut f = make_table(&name, "LOGCOL", "1L", 4, &mut status);
            fits_write_colnull_log(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                4,
                &data,
                nulval,
                &mut status,
            );
            assert_eq!(status, 0);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_bit_column() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let bits: [c_char; 8] = [1, 0, 1, 0, 1, 0, 1, 0];

            let mut f = make_table(&name, "BITCOL", "8X", 1, &mut status);
            fits_write_col_bit(f.as_deref_mut().unwrap(), 1, 1, 1, 8, &bits, &mut status);
            assert_eq!(status, 0);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            assert_eq!(status, 0);
            let mut result = [0 as c_char; 8];
            fits_read_col_bit(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                8,
                &mut result,
                &mut status,
            );
            assert_eq!(status, 0);
            assert!(!(result[0] != 1));
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_bit_bad_row() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let bits: [c_char; 1] = [1];

            let mut f = make_table(&name, "BITCOL", "8X", 1, &mut status);
            assert_eq!(status, 0);
            fits_write_col_bit(f.as_deref_mut().unwrap(), 1, 0, 1, 1, &bits, &mut status);
            assert!(!(status != BAD_ROW_NUM));
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_bit_bad_elem() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let bits: [c_char; 1] = [1];

            let mut f = make_table(&name, "BITCOL", "8X", 1, &mut status);
            assert_eq!(status, 0);
            fits_write_col_bit(f.as_deref_mut().unwrap(), 1, 1, 0, 1, &bits, &mut status);
            assert!(!(status != BAD_ELEM_NUM));
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_bit_not_logical() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let bits: [c_char; 1] = [1];

            let mut f = make_table(&name, "INTCOL", "1J", 1, &mut status);
            assert_eq!(status, 0);
            fits_write_col_bit(f.as_deref_mut().unwrap(), 1, 1, 1, 1, &bits, &mut status);
            assert!(!(status != NOT_LOGICAL_COL));
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_bit_beyond_repeat() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let bits: [c_char; 1] = [1];

            let mut f = make_table(&name, "BITCOL", "8X", 1, &mut status);
            assert_eq!(status, 0);
            fits_write_col_bit(f.as_deref_mut().unwrap(), 1, 1, 9, 1, &bits, &mut status);
            assert!(!(status != BAD_ELEM_NUM));
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_varlen_bit() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let bits: [c_char; 6] = [1, 0, 1, 1, 0, 1];

            let mut f = make_table(&name, "VARBIT", "1PX", 1, &mut status);
            fits_write_col_bit(f.as_deref_mut().unwrap(), 1, 1, 1, 6, &bits, &mut status);
            assert_eq!(status, 0);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            assert_eq!(status, 0);
            let mut result = [0 as c_char; 6];
            fits_read_col_bit(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                6,
                &mut result,
                &mut status,
            );
            assert_eq!(status, 0);
            assert!(!(result[0] != 1));
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_logical_spanning_rows() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [c_char; 9] = [1, 0, 1, 0, 1, 0, 1, 0, 1];

            let mut f = make_table(&name, "LOGCOL", "3L", 3, &mut status);
            fits_write_col_log(f.as_deref_mut().unwrap(), 1, 1, 1, 9, &data, &mut status);
            assert_eq!(status, 0);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            assert_eq!(status, 0);
            let mut result = [0 as c_char; 9];
            let mut anynull = 0;
            fits_read_col_log(
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
            for i in 0..9 {
                assert!(!(result[i] != data[i]));
            }
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_varlen_logical_with_nulls() {
        // Test writing variable length logical array with nulls.
        // This exercises the tcode < 0 path in ffpcnl.
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [c_char; 5] = [1, 0, 1, 0, 1];
            let nulval: c_char = 0;

            let mut f = make_table(&name, "VARCOL", "1PL", 1, &mut status);
            fits_write_colnull_log(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                5,
                &data,
                nulval,
                &mut status,
            );
            assert_eq!(status, 0);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            assert_eq!(status, 0);
            let mut result = [0 as c_char; 5];
            let mut anynull = 0;
            fits_read_col_log(
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
            assert!(!(result[0] != 1));
            assert!(!(result[2] != 1));
            assert!(!(result[4] != 1));
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_logical_nulls_mixed() {
        // Test writing logical values with mixed null/non-null values.
        // Exercise the nbad path with nulls at different positions.
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [c_char; 5] = [1, 99, 99, 1, 99];
            let nulval: c_char = 99;

            let mut f = make_table(&name, "LOGCOL", "5L", 1, &mut status);
            fits_write_colnull_log(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                5,
                &data,
                nulval,
                &mut status,
            );
            assert_eq!(status, 0);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            assert_eq!(status, 0);
            let mut result = [0 as c_char; 5];
            let mut anynull = 0;
            fits_read_col_log(
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
            // Non-null values should be 1
            assert!(!(result[0] != 1));
            assert!(!(result[3] != 1));
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_bit_bad_elem_beyond() {
        // Test BAD_ELEM_NUM error when fbit exceeds repeat count.
        // This exercises the fbyte > repeat check in ffpclx.
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let bits: [c_char; 1] = [1];

            let mut f = make_table(&name, "BITCOL", "8X", 1, &mut status);
            assert_eq!(status, 0);
            // Try to write at bit position 100 which exceeds 8-bit column
            fits_write_col_bit(f.as_deref_mut().unwrap(), 1, 1, 100, 1, &bits, &mut status);
            assert!(!(status != BAD_ELEM_NUM));
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_bit_byte_col_overflow() {
        // Test BAD_ELEM_NUM when writing bits beyond byte column width.
        // Uses single byte column to trigger the fbyte > repeat path.
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let bits: [c_char; 1] = [1];

            let mut f = make_table(&name, "BYTECOL", "1B", 1, &mut status); // 1 byte = 8 bits
            assert_eq!(status, 0);
            // Try to write at bit 9 which exceeds 8 bits (1 byte)
            fits_write_col_bit(f.as_deref_mut().unwrap(), 1, 1, 9, 1, &bits, &mut status);
            assert!(!(status != BAD_ELEM_NUM));
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_bits_spanning_rows() {
        // Test writing bits that span multiple rows.
        // Exercises the row boundary handling in ffpclx.
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let bits: [c_char; 12] = [1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0];

            let mut f = make_table(&name, "BITCOL", "4X", 3, &mut status); // 4 bits per row
            // Write 12 bits starting at row 1, which spans all 3 rows
            fits_write_col_bit(f.as_deref_mut().unwrap(), 1, 1, 1, 12, &bits, &mut status);
            assert_eq!(status, 0);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            assert_eq!(status, 0);
            let mut result = [0 as c_char; 12];
            fits_read_col_bit(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                4,
                &mut result[0..4],
                &mut status,
            );
            assert_eq!(status, 0);
            assert!(!(result[0] != 1));
            fits_read_col_bit(
                f.as_deref_mut().unwrap(),
                1,
                2,
                1,
                4,
                &mut result[4..8],
                &mut status,
            );
            assert_eq!(status, 0);
            assert!(!(result[4] != 1));
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_bits_to_new_row() {
        // Test writing bits to a row beyond current table size. When writing
        // past EOF, ffpclu reads cbuff which triggers the END_OF_FILE case and
        // sets status=0, cbuff=0
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let bits: [c_char; 4] = [1, 0, 1, 0];

            // Create table with 0 rows initially
            let mut f = make_table(&name, "BITCOL", "8X", 0, &mut status);
            // Write bits to row 1 which extends the table
            fits_write_col_bit(f.as_deref_mut().unwrap(), 1, 1, 1, 4, &bits, &mut status);
            assert_eq!(status, 0);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            assert_eq!(status, 0);
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
            // Bits 1,0,1,0 starting at bit 1 = binary 10100000 = 0xA0 = 160
            assert!(!(result[0] != 0xA0));
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_write_zero_bits() {
        // Test ffpclx with nbit=0 for early return
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let bits: [c_char; 1] = [1];

            let mut f = make_table(&name, "BITCOL", "8X", 1, &mut status);
            assert_eq!(status, 0);
            // nbit=0 should return early without error
            fits_write_col_bit(f.as_deref_mut().unwrap(), 1, 1, 1, 0, &bits, &mut status);
            assert!(!(status != 0));
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_calls_with_error_status() {
        // Test that functions return early when status > 0
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [c_char; 1] = [1];

            let mut f = make_table(&name, "LOGCOL", "1L", 1, &mut status);
            assert_eq!(status, 0);

            // Call ffpcll with pre-existing error
            status = 1;
            fits_write_col_log(f.as_deref_mut().unwrap(), 1, 1, 1, 1, &data, &mut status);
            assert!(!(status != 1));

            // Call ffpcnl with pre-existing error
            status = 1;
            fits_write_colnull_log(f.as_deref_mut().unwrap(), 1, 1, 1, 1, &data, 0, &mut status);
            assert!(!(status != 1));

            // Call ffpclx with pre-existing error
            status = 1;
            fits_write_col_bit(f.as_deref_mut().unwrap(), 1, 1, 1, 1, &data, &mut status);
            assert!(!(status != 1));

            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }
}
