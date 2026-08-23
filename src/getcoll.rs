//! Routines that read data elements from a FITS image or table, with
//! logical datatype.
//!
//! The [`TLOGICAL`] arm of the typed column I/O family. [`crate::getcol`]
//! dispatches here when a caller asks for this datatype at run time; the
//! write-side counterpart is [`crate::putcoll`].
//!
//! Ported from CFITSIO's `getcoll.c`, written by William Pence at the High
//! Energy Astrophysics Science Archive Research Center (HEASARC), NASA Goddard
//! Space Flight Center.
#![warn(missing_docs)]

use core::cmp;
use core::slice;

use crate::c_types::{c_char, c_int, c_long, c_uchar, c_uint, c_ushort};

use bytemuck::cast_slice_mut;

use crate::NullCheckType;
use crate::fitscore::{ffgcprll, ffgdes_safe, ffmahd_safe, ffpmsg_slice, ffpmsg_str, ffrdef_safe};
use crate::fitsio::*;
use crate::fitsio2::*;
use crate::getcolui::ffgcvui_safe;
use crate::getcoluk::ffgcvuk_safe;
use crate::{buffers::*, int_snprintf};

/// Read an array of logical values from a column in the current FITS HDU.
///
/// Any undefined pixels will be set equal to the value of 'nulval' unless
/// nulval = 0 in which case no checks for undefined pixels will be made.
///
/// # Parameters
///
/// * `fptr`      — (I) FITS file pointer
/// * `colnum`    — (I) number of column to read (1 = 1st col)
/// * `firstrow`  — (I) first row to read (1 = 1st row)
/// * `firstelem` — (I) first vector element to read (1 = 1st)
/// * `nelem`     — (I) number of values to read
/// * `nulval`    — (I) value for null pixels
/// * `array`     — (O) array of values
/// * `anynul`    — (O) set to 1 if any values are null; else 0
/// * `status`    — (IO) error status
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgcvl(
    fptr: *mut fitsfile,
    colnum: c_int,
    firstrow: LONGLONG,
    firstelem: LONGLONG,
    nelem: LONGLONG,
    nulval: c_char,
    array: *mut c_char,
    anynul: *mut c_int,
    status: *mut c_int,
) -> c_int {
    // FFI WRAPPER
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let anynul = anynul.as_mut();

        let array = slice::from_raw_parts_mut(array, nelem as usize);

        ffgcvl_safe(
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

/// Read an array of logical values from a column in the current FITS HDU.
///
/// Any undefined pixels will be set equal to the value of 'nulval' unless
/// nulval = 0 in which case no checks for undefined pixels will be made.
///
/// # Parameters
///
/// * `fptr`      — (I) FITS file pointer
/// * `colnum`    — (I) number of column to read (1 = 1st col)
/// * `firstrow`  — (I) first row to read (1 = 1st row)
/// * `firstelem` — (I) first vector element to read (1 = 1st)
/// * `nelem`     — (I) number of values to read
/// * `nulval`    — (I) value for null pixels
/// * `array`     — (O) array of values
/// * `anynul`    — (O) set to 1 if any values are null; else 0
/// * `status`    — (IO) error status
pub fn ffgcvl_safe(
    fptr: &mut fitsfile,
    colnum: c_int,
    firstrow: LONGLONG,
    firstelem: LONGLONG,
    nelem: LONGLONG,
    nulval: c_char,
    array: &mut [c_char],
    anynul: Option<&mut c_int>,
    status: &mut c_int,
) -> c_int {
    let mut cdummy: [c_char; 1] = [0];

    ffgcll(
        fptr,
        colnum,
        firstrow,
        firstelem,
        nelem,
        NullCheckType::SetPixel,
        nulval,
        array,
        &mut cdummy,
        anynul,
        status,
    );

    *status
}

/// Read an array of logical values from a column in the current FITS HDU.
///
/// !!!! THIS ROUTINE IS DEPRECATED AND SHOULD NOT BE USED !!!!!!
///           !!!! USE ffgcvl INSTEAD  !!!!!!
/// No checking for null values will be performed.
// #[deprecated]
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
/// # Parameters
///
/// * `fptr`      — (I) FITS file pointer
/// * `colnum`    — (I) number of column to read (1 = 1st col)
/// * `firstrow`  — (I) first row to read (1 = 1st row)
/// * `firstelem` — (I) first vector element to read (1 = 1st)
/// * `nelem`     — (I) number of values to read
/// * `array`     — (O) array of values
/// * `status`    — (IO) error status
pub unsafe extern "C" fn ffgcl(
    fptr: *mut fitsfile,
    colnum: c_int,
    firstrow: LONGLONG,
    firstelem: LONGLONG,
    nelem: LONGLONG,
    array: *mut c_char,
    status: *mut c_int,
) -> c_int {
    // FFI WRAPPER
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        let array = slice::from_raw_parts_mut(array, nelem as usize);

        #[allow(deprecated)]
        ffgcl_safe(
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

/// Read an array of logical values from a column in the current FITS HDU.
///
/// !!!! THIS ROUTINE IS DEPRECATED AND SHOULD NOT BE USED !!!!!!
///           !!!! USE ffgcvl INSTEAD  !!!!!!
/// No checking for null values will be performed.
///
/// # Parameters
///
/// * `fptr`      — (I) FITS file pointer
/// * `colnum`    — (I) number of column to read (1 = 1st col)
/// * `firstrow`  — (I) first row to read (1 = 1st row)
/// * `firstelem` — (I) first vector element to read (1 = 1st)
/// * `nelem`     — (I) number of values to read
/// * `array`     — (O) array of values
/// * `status`    — (IO) error status
#[deprecated]
pub fn ffgcl_safe(
    fptr: &mut fitsfile,
    colnum: c_int,
    firstrow: LONGLONG,
    firstelem: LONGLONG,
    nelem: LONGLONG,
    array: &mut [c_char],
    status: &mut c_int,
) -> c_int {
    let mut anynul: c_int = 0;
    let nulval: c_char = 0;

    ffgcvl_safe(
        fptr,
        colnum,
        firstrow as LONGLONG,
        firstelem as LONGLONG,
        nelem as LONGLONG,
        nulval,
        array,
        Some(&mut anynul),
        status,
    )
}

/// Read an array of logical values from a column in the current FITS HDU.
///
/// # Parameters
///
/// * `fptr`      — (I) FITS file pointer
/// * `colnum`    — (I) number of column to read (1 = 1st col)
/// * `firstrow`  — (I) first row to read (1 = 1st row)
/// * `firstelem` — (I) first vector element to read (1 = 1st)
/// * `nelem`     — (I) number of values to read
/// * `array`     — (O) array of values
/// * `nularray`  — (O) array of flags = 1 if nultyp = 2
/// * `anynul`    — (O) set to 1 if any values are null; else 0
/// * `status`    — (IO) error status
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgcfl(
    fptr: *mut fitsfile,
    colnum: c_int,
    firstrow: LONGLONG,
    firstelem: LONGLONG,
    nelem: LONGLONG,
    array: *mut c_char,
    nularray: *mut c_char,
    anynul: *mut c_int,
    status: *mut c_int,
) -> c_int {
    // FFI WRAPPER
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let anynul = anynul.as_mut();

        let array = slice::from_raw_parts_mut(array, nelem as usize);
        let nularray = slice::from_raw_parts_mut(nularray, nelem as usize);

        ffgcfl_safe(
            fptr, colnum, firstrow, firstelem, nelem, array, nularray, anynul, status,
        )
    }
}

/// Read an array of logical values from a column in the current FITS HDU.
///
/// # Parameters
///
/// * `fptr`      — (I) FITS file pointer
/// * `colnum`    — (I) number of column to read (1 = 1st col)
/// * `firstrow`  — (I) first row to read (1 = 1st row)
/// * `firstelem` — (I) first vector element to read (1 = 1st)
/// * `nelem`     — (I) number of values to read
/// * `array`     — (O) array of values
/// * `nularray`  — (O) array of flags = 1 if nultyp = 2
/// * `anynul`    — (O) set to 1 if any values are null; else 0
/// * `status`    — (IO) error status
pub fn ffgcfl_safe(
    fptr: &mut fitsfile,
    colnum: c_int,
    firstrow: LONGLONG,
    firstelem: LONGLONG,
    nelem: LONGLONG,
    array: &mut [c_char],
    nularray: &mut [c_char],
    anynul: Option<&mut c_int>,
    status: &mut c_int,
) -> c_int {
    let nulval: c_char = 0;

    ffgcll(
        fptr,
        colnum,
        firstrow,
        firstelem,
        nelem,
        NullCheckType::SetNullArray,
        nulval,
        array,
        nularray,
        anynul,
        status,
    );
    *status
}

/// Read an array of logical values from a column in the current FITS HDU.
///
/// # Parameters
///
/// * `fptr`      — (I) FITS file pointer
/// * `colnum`    — (I) number of column to read (1 = 1st col)
/// * `firstrow`  — (I) first row to read (1 = 1st row)
/// * `firstelem` — (I) first vector element to read (1 = 1st)
/// * `nelem`     — (I) number of values to read
/// * `nultyp`    — (I) null value handling code:
/// * `nulval`    — (I) value for null pixels if nultyp = 1
/// * `array`     — (O) array of values
/// * `nularray`  — (O) array of flags = 1 if nultyp = 2
/// * `anynul`    — (O) set to 1 if any values are null; else 0
/// * `status`    — (IO) error status
pub(crate) fn ffgcll(
    fptr: &mut fitsfile,
    colnum: c_int,
    firstrow: LONGLONG,
    firstelem: LONGLONG,
    nelem: LONGLONG,
    nultyp: NullCheckType,
    /*     1: set undefined pixels = nulval        */
    /*     2: set nularray=1 for undefined pixels  */
    nulval: c_char,
    array: &mut [c_char],
    nularray: &mut [c_char],
    mut anynul: Option<&mut c_int>,
    status: &mut c_int,
) -> c_int {
    let mut dtemp: f64 = 0.0;
    let mut tcode: c_int = 0;
    let mut maxelem: c_int = 0;
    let mut hdutype: c_int = 0;
    let mut nulcheck = NullCheckType::None;
    let mut twidth: c_long = 0;
    let mut incre: c_long = 0;
    let mut ntodo: c_long = 0;
    let mut repeat: LONGLONG = 0;
    let mut startpos: LONGLONG = 0;
    let mut elemnum: LONGLONG = 0;
    let mut readptr: LONGLONG = 0;
    let mut tnull: LONGLONG = 0;
    let mut rowlen: LONGLONG = 0;
    let mut rownum: LONGLONG = 0;
    let mut remain: LONGLONG = 0;
    let mut next: usize = 0;
    let mut scale: f64 = 0.0;
    let mut zero: f64 = 0.0;
    let mut tform: [c_char; 20] = [0; 20];
    let mut message: [c_char; FLEN_ERRMSG] = [0; FLEN_ERRMSG];
    let mut snull: [c_char; 20] = [0; 20]; /*  the FITS null value  */
    let mut buffer: [c_uchar; DBUFFSIZE as usize] = [0; DBUFFSIZE as usize];

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
    /*  Check input and get parameters about the column: */
    if ffgcprll(
        fptr,
        colnum,
        firstrow,
        firstelem,
        nelem,
        0,
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
    /*  Decide whether to check for null values in the input FITS file: */
    nulcheck = nultyp; /* by default, check for null values in the FITS file */

    if nultyp == NullCheckType::SetPixel && nulval == 0 {
        nulcheck = NullCheckType::None; /* calling routine does not want to check for nulls */
    }
    /*  Now read the logical values from the FITS column.                  */

    remain = nelem; /* remaining number of values to read */
    next = 0; /* next element in array to be read   */
    rownum = 0; /* row number, relative to firstrow   */
    ntodo = remain as c_long; /* max number of elements to read at one time */

    while ntodo > 0 {
        /*
         limit the number of pixels to read at one time to the number that
         remain in the current vector.
        */
        ntodo = cmp::min(ntodo as LONGLONG, LONGLONG::from(maxelem)) as c_long;
        ntodo = cmp::min(ntodo as LONGLONG, repeat - elemnum) as c_long;

        readptr = startpos + (rowlen * rownum) + (elemnum * incre as LONGLONG);

        ffgi1b(fptr, readptr, ntodo, incre, &mut buffer, status);

        /* convert from T or F to 1 or 0 */

        for bi in 0..(ntodo as usize) {
            if buffer[bi] == b'T' {
                array[next] = 1;
            } else if buffer[bi] == b'F' {
                array[next] = 0;
            } else if buffer[bi] == 0 {
                array[next] = nulval; /* set null values to input nulval */
                if let Some(anynul) = anynul.as_deref_mut() {
                    *anynul = 1;
                }
                if nulcheck == NullCheckType::SetNullArray {
                    nularray[next] = 1; /* set null flags */
                }
            } else {
                /* some other illegal character; return the char value */

                if buffer[bi] == 1 {
                    /* this is an unfortunate case where the illegal value is the same
                    as what we set True values to, so set the value to the character '1'
                    instead, which has ASCII value 49.  */
                    array[next] = 49;
                } else {
                    array[next] = buffer[bi] as c_char;
                }
            }
            next += 1;
        }

        if *status > 0 {
            /* test for error during previous read operation */

            dtemp = next as f64;
            int_snprintf!(
                &mut message,
                FLEN_ERRMSG,
                "Error reading elements {:.0} thruough {:.0} of logical array (ffgcl).",
                dtemp + 1.0,
                dtemp + ntodo as f64,
            );
            ffpmsg_slice(&message);
            return *status;
        }

        /*  increment the counters for the next loop  */
        remain -= ntodo as LONGLONG;
        if remain > 0 {
            elemnum += ntodo as LONGLONG;

            if elemnum == repeat
            /* completed a row; start on later row */
            {
                elemnum = 0;
                rownum += 1;
            }
        }
        ntodo = remain as c_long; /* this is the maximum number to do in next loop */
    } /*  End of main while Loop  */

    *status
}

/// read an array of logical values from a specified bit or byte
/// column of the binary table.    larray is set = TRUE, if the corresponding
/// bit = 1, otherwise it is set to FALSE.
/// The binary table column being read from must have datatype 'B' or 'X'.
///
/// # Parameters
///
/// * `fptr`   — (I) FITS file pointer
/// * `colnum` — (I) number of column to write (1 = 1st col)
/// * `frow`   — (I) first row to write (1 = 1st row)
/// * `fbit`   — (I) first bit to write (1 = 1st)
/// * `nbit`   — (I) number of bits to write
/// * `larray` — (O) array of logicals corresponding to bits
/// * `status` — (IO) error status
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgcx(
    fptr: *mut fitsfile,
    colnum: c_int,
    frow: LONGLONG,
    fbit: LONGLONG,
    nbit: LONGLONG,
    larray: *mut c_char,
    status: *mut c_int,
) -> c_int {
    // FFI WRAPPER
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        let larray = slice::from_raw_parts_mut(larray, nbit as usize);

        ffgcx_safe(fptr, colnum, frow, fbit, nbit, larray, status)
    }
}

/// read an array of logical values from a specified bit or byte
/// column of the binary table.    larray is set = TRUE, if the corresponding
/// bit = 1, otherwise it is set to FALSE.
/// The binary table column being read from must have datatype 'B' or 'X'.
///
/// # Parameters
///
/// * `fptr`   — (I) FITS file pointer
/// * `colnum` — (I) number of column to write (1 = 1st col)
/// * `frow`   — (I) first row to write (1 = 1st row)
/// * `fbit`   — (I) first bit to write (1 = 1st)
/// * `nbit`   — (I) number of bits to write
/// * `larray` — (O) array of logicals corresponding to bits
/// * `status` — (IO) error status
pub fn ffgcx_safe(
    fptr: &mut fitsfile,
    colnum: c_int,
    frow: LONGLONG,
    fbit: LONGLONG,
    nbit: LONGLONG,
    larray: &mut [c_char],
    status: &mut c_int,
) -> c_int {
    let mut bstart: LONGLONG = 0;
    let mut offset: c_long = 0;
    let mut ndone: c_long = 0;
    let mut repeat: c_long = 0;
    let mut bitloc: c_long = 0;
    let mut fbyte: c_long = 0;
    let mut rstart: LONGLONG = 0;
    let mut estart: LONGLONG = 0;
    let mut tcode: c_int = 0;
    let mut descrp: c_int = 0;
    let mut cbuff: [u8; 1] = [0];
    static ONBIT: [u8; 8] = [128, 64, 32, 16, 8, 4, 2, 1];

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

    /* position to the correct HDU */
    if fptr.HDUposition != fptr.Fptr.curhdu {
        ffmahd_safe(fptr, (fptr.HDUposition) + 1, None, status);

    /* rescan header if data structure is undefined */
    } else if fptr.Fptr.datastart == DATA_UNDEFINED as LONGLONG && ffrdef_safe(fptr, status) > 0 {
        return *status;
    }

    fbyte = ((fbit + 7) / 8) as c_long;
    bitloc = (fbit - 1 - ((fbit - 1) / 8 * 8)) as c_long;
    ndone = 0;
    rstart = frow - 1;
    estart = fbyte as LONGLONG - 1;

    let c = fptr.Fptr.get_tableptr_as_slice(); /* point to first column structure */
    let ci = (colnum - 1) as usize; /* offset to the correct column */
    let tcol = c[ci];

    tcode = tcol.tdatatype;

    if (tcode.abs()) > TBYTE {
        *status = NOT_LOGICAL_COL;
        return *status;
    } /* not correct datatype column */

    if tcode > 0 {
        descrp = FALSE as c_int; /* not a variable length descriptor column */
        /* N.B: REPEAT is the number of bytes, not number of bits */
        repeat = tcol.trepeat as c_long;

        if tcode == TBIT {
            repeat = (repeat + 7) / 8; /* convert from bits to bytes */
        }

        if fbyte > repeat {
            *status = BAD_ELEM_NUM;
            return *status;
        }

        /* calc the i/o pointer location to start of sequence of pixels */
        bstart = fptr.Fptr.datastart + (fptr.Fptr.rowlength * rstart) + tcol.tbcol + estart;
    } else {
        descrp = TRUE as c_int; /* a variable length descriptor column */
        /* only bit arrays (tform = 'X') are supported for variable */
        /* length arrays.  REPEAT is the number of BITS in the array. */

        ffgdes_safe(
            fptr,
            colnum,
            frow,
            Some(&mut repeat),
            Some(&mut offset),
            status,
        );

        if tcode == -TBIT {
            repeat = (repeat + 7) / 8;
        }

        if (fbit + nbit + 6) / 8 > repeat as LONGLONG {
            *status = BAD_ELEM_NUM;
            return *status;
        }

        /* calc the i/o pointer location to start of sequence of pixels */
        bstart = fptr.Fptr.datastart + offset as LONGLONG + fptr.Fptr.heapstart + estart;
    }

    /* move the i/o pointer to the start of the pixel sequence */
    if ffmbyt_safe(fptr, bstart, REPORT_EOF, status) > 0 {
        return *status;
    }

    /* read the next byte */
    loop {
        if ffgbyt(fptr, 1, &mut cbuff, status) > 0 {
            return *status;
        }

        let mut ii = bitloc as usize;
        while (ii < 8) && ((ndone as LONGLONG) < nbit) {
            if (cbuff[0] & ONBIT[ii]) != 0 {
                /* test if bit is set */
                larray[ndone as usize] = TRUE as c_char;
            } else {
                larray[ndone as usize] = FALSE as c_char;
            }
            ii += 1;
            ndone += 1;
        }

        if ndone as LONGLONG == nbit {
            /* finished all the bits */
            return *status;
        }

        /* not done, so get the next byte */
        if descrp == 0 {
            estart += 1;
            if estart == repeat as LONGLONG {
                /* move the i/o pointer to the next row of pixels */
                estart = 0;
                rstart += 1;
                bstart = fptr.Fptr.datastart + (fptr.Fptr.rowlength * rstart) + tcol.tbcol;

                ffmbyt_safe(fptr, bstart, REPORT_EOF, status);
            }
        }
        bitloc = 0;
    }
}

/// Read a consecutive string of bits from an 'X' or 'B' column and
/// interpret them as an unsigned integer.  The number of bits must be
/// less than or equal to 16 or the total number of bits in the column,
/// which ever is less.
///
/// # Parameters
///
/// * `fptr`            — (I) FITS file pointer
/// * `colnum`          — (I) number of column to read (1 = 1st col)
/// * `firstrow`        — (I) first row to read (1 = 1st row)
/// * `nrows`           — (I) no. of rows to read
/// * `input_first_bit` — (I) first bit to read (1 = 1st)
/// * `input_nbits`     — (I) number of bits to read (<= 32)
/// * `array`           — (O) array of integer values
/// * `status`          — (IO) error status
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgcxui(
    fptr: *mut fitsfile,
    colnum: c_int,
    firstrow: LONGLONG,
    nrows: LONGLONG,
    input_first_bit: c_long,
    input_nbits: c_int,
    array: *mut c_ushort,
    status: *mut c_int,
) -> c_int {
    // FFI WRAPPER
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        let array = slice::from_raw_parts_mut(array, nrows as usize);

        ffgcxui_safe(
            fptr,
            colnum,
            firstrow,
            nrows,
            input_first_bit,
            input_nbits,
            array,
            status,
        )
    }
}

/// Read a consecutive string of bits from an 'X' or 'B' column and
/// interpret them as an unsigned integer.  The number of bits must be
/// less than or equal to 16 or the total number of bits in the column,
/// which ever is less.
#[allow(clippy::if_same_then_else)]
// C dispatch chain: distinct conditions deliberately share an action.
/// # Parameters
///
/// * `fptr`            — (I) FITS file pointer
/// * `colnum`          — (I) number of column to read (1 = 1st col)
/// * `firstrow`        — (I) first row to read (1 = 1st row)
/// * `nrows`           — (I) no. of rows to read
/// * `input_first_bit` — (I) first bit to read (1 = 1st)
/// * `input_nbits`     — (I) number of bits to read (<= 32)
/// * `array`           — (O) array of integer values
/// * `status`          — (IO) error status
pub fn ffgcxui_safe(
    fptr: &mut fitsfile,
    colnum: c_int,
    firstrow: LONGLONG,
    nrows: LONGLONG,
    input_first_bit: c_long,
    input_nbits: c_int,
    array: &mut [c_ushort],
    status: &mut c_int,
) -> c_int {
    let mut colbyte: [c_ushort; 5] = [0; 5];
    let mut message: [c_char; FLEN_ERRMSG] = [0; FLEN_ERRMSG];

    if *status > 0 || nrows == 0 {
        return *status;
    }

    /*  check input parameters */
    if firstrow < 1 {
        int_snprintf!(
            &mut message,
            FLEN_ERRMSG,
            "Starting row number is less than 1: {} (ffgcxui)",
            firstrow
        );
        ffpmsg_slice(&message);
        *status = BAD_ROW_NUM;
        return *status;
    } else if input_first_bit < 1 {
        int_snprintf!(
            &mut message,
            FLEN_ERRMSG,
            "Starting bit number is less than 1: {} (ffgcxui)",
            input_first_bit
        );
        ffpmsg_slice(&message);
        *status = BAD_ELEM_NUM;
        return *status;
    } else if input_nbits > 16 {
        int_snprintf!(
            &mut message,
            FLEN_ERRMSG,
            "Number of bits to read is > 16: {} (ffgcxui)",
            input_nbits
        );
        ffpmsg_slice(&message);
        *status = BAD_ELEM_NUM;
        return *status;
    }

    /* position to the correct HDU */
    if fptr.HDUposition != fptr.Fptr.curhdu {
        ffmahd_safe(fptr, (fptr.HDUposition) + 1, None, status);

    /* rescan header if data structure is undefined */
    } else if fptr.Fptr.datastart == DATA_UNDEFINED as LONGLONG && ffrdef_safe(fptr, status) > 0 {
        return *status;
    }

    if fptr.Fptr.hdutype != BINARY_TBL {
        ffpmsg_str("This is not a binary table extension (ffgcxui)");
        *status = NOT_BTABLE;
        return *status;
    }

    if colnum > fptr.Fptr.tfield {
        int_snprintf!(
            &mut message,
            FLEN_ERRMSG,
            "Specified column number is out of range: {} (ffgcxui)",
            colnum
        );
        ffpmsg_slice(&message);
        int_snprintf!(
            &mut message,
            FLEN_ERRMSG,
            "  There are {} columns in this table.",
            fptr.Fptr.tfield
        );
        ffpmsg_slice(&message);

        *status = BAD_COL_NUM;
        return *status;
    }

    let c = fptr.Fptr.get_tableptr_as_slice(); /* point to first column */
    let ci = (colnum - 1) as usize; /* offset to correct column structure */
    let tcol = c[ci];

    if tcol.tdatatype.abs() > TBYTE {
        ffpmsg_str("Can only read bits from X or B type columns. (ffgcxui)");
        *status = NOT_LOGICAL_COL; /* not correct datatype column */
        return *status;
    }

    let firstbyte: c_int = ((input_first_bit - 1) / 8 + 1) as c_int;
    let lastbyte: c_int = ((input_first_bit + input_nbits as c_long - 2) / 8 + 1) as c_int;
    let nbytes: c_int = lastbyte - firstbyte + 1;

    if tcol.tdatatype == TBIT
        && input_first_bit + input_nbits as c_long - 1 > tcol.trepeat as c_long
    {
        ffpmsg_str("Too many bits. Tried to read past width of column (ffgcxui)");
        *status = BAD_ELEM_NUM;
        return *status;
    } else if tcol.tdatatype == TBYTE && lastbyte as c_long > tcol.trepeat as c_long {
        ffpmsg_str("Too many bits. Tried to read past width of column (ffgcxui)");
        *status = BAD_ELEM_NUM;
        return *status;
    }

    for ii in 0..nrows {
        /* read the relevant bytes from the row */
        if ffgcvui_safe(
            fptr,
            colnum,
            firstrow + ii,
            firstbyte as LONGLONG,
            nbytes as LONGLONG,
            0,
            &mut colbyte,
            None,
            status,
        ) > 0
        {
            ffpmsg_str("Error reading bytes from column (ffgcxui)");
            return *status;
        }

        let mut firstbit: c_int = ((input_first_bit - 1) % 8) as c_int; /* modulus operator */
        let mut nbits: c_int = input_nbits;

        array[ii as usize] = 0;

        /* select and shift the bits from each byte into the output word */
        while nbits != 0 {
            let bytenum: c_int = firstbit / 8;

            let startbit: c_int = firstbit % 8;
            let numbits: c_int = cmp::min(nbits, 8 - startbit);
            let endbit: c_int = startbit + numbits - 1;

            let rshift: c_int = 7 - endbit;
            let lshift: c_int = nbits - numbits;

            array[ii as usize] |= (colbyte[bytenum as usize] >> rshift) << lshift;

            nbits -= numbits;
            firstbit += numbits;
        }
    }

    *status
}

/// Read a consecutive string of bits from an 'X' or 'B' column and
/// interpret them as an unsigned integer.  The number of bits must be
/// less than or equal to 32 or the total number of bits in the column,
/// which ever is less.
///
/// # Parameters
///
/// * `fptr`            — (I) FITS file pointer
/// * `colnum`          — (I) number of column to read (1 = 1st col)
/// * `firstrow`        — (I) first row to read (1 = 1st row)
/// * `nrows`           — (I) no. of rows to read
/// * `input_first_bit` — (I) first bit to read (1 = 1st)
/// * `input_nbits`     — (I) number of bits to read (<= 32)
/// * `array`           — (O) array of integer values
/// * `status`          — (IO) error status
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgcxuk(
    fptr: *mut fitsfile,
    colnum: c_int,
    firstrow: LONGLONG,
    nrows: LONGLONG,
    input_first_bit: c_long,
    input_nbits: c_int,
    array: *mut c_uint,
    status: *mut c_int,
) -> c_int {
    // FFI WRAPPER
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        let array = slice::from_raw_parts_mut(array, nrows as usize);

        ffgcxuk_safe(
            fptr,
            colnum,
            firstrow,
            nrows,
            input_first_bit,
            input_nbits,
            array,
            status,
        )
    }
}

/// Read a consecutive string of bits from an 'X' or 'B' column and
/// interpret them as an unsigned integer.  The number of bits must be
/// less than or equal to 32 or the total number of bits in the column,
/// which ever is less.
#[allow(clippy::if_same_then_else)]
// C dispatch chain: distinct conditions deliberately share an action.
/// # Parameters
///
/// * `fptr`            — (I) FITS file pointer
/// * `colnum`          — (I) number of column to read (1 = 1st col)
/// * `firstrow`        — (I) first row to read (1 = 1st row)
/// * `nrows`           — (I) no. of rows to read
/// * `input_first_bit` — (I) first bit to read (1 = 1st)
/// * `input_nbits`     — (I) number of bits to read (<= 32)
/// * `array`           — (O) array of integer values
/// * `status`          — (IO) error status
pub fn ffgcxuk_safe(
    fptr: &mut fitsfile,
    colnum: c_int,
    firstrow: LONGLONG,
    nrows: LONGLONG,
    input_first_bit: c_long,
    input_nbits: c_int,
    array: &mut [c_uint],
    status: &mut c_int,
) -> c_int {
    let mut colbyte: [c_uint; 5] = [0; 5];
    let mut message: [c_char; FLEN_ERRMSG] = [0; FLEN_ERRMSG];

    if *status > 0 || nrows == 0 {
        return *status;
    }

    /*  check input parameters */
    if firstrow < 1 {
        int_snprintf!(
            &mut message,
            FLEN_ERRMSG,
            "Starting row number is less than 1: {} (ffgcxuk)",
            firstrow
        );
        ffpmsg_slice(&message);
        *status = BAD_ROW_NUM;
        return *status;
    } else if input_first_bit < 1 {
        int_snprintf!(
            &mut message,
            FLEN_ERRMSG,
            "Starting bit number is less than 1: {} (ffgcxuk)",
            input_first_bit
        );
        ffpmsg_slice(&message);
        *status = BAD_ELEM_NUM;
        return *status;
    } else if input_nbits > 32 {
        int_snprintf!(
            &mut message,
            FLEN_ERRMSG,
            "Number of bits to read is > 32: {} (ffgcxuk)",
            input_nbits
        );
        ffpmsg_slice(&message);
        *status = BAD_ELEM_NUM;
        return *status;
    }

    /* position to the correct HDU */
    if fptr.HDUposition != fptr.Fptr.curhdu {
        ffmahd_safe(fptr, (fptr.HDUposition) + 1, None, status);

    /* rescan header if data structure is undefined */
    } else if fptr.Fptr.datastart == DATA_UNDEFINED as LONGLONG && ffrdef_safe(fptr, status) > 0 {
        return *status;
    }

    if fptr.Fptr.hdutype != BINARY_TBL {
        ffpmsg_str("This is not a binary table extension (ffgcxuk)");
        *status = NOT_BTABLE;
        return *status;
    }

    if colnum > fptr.Fptr.tfield {
        int_snprintf!(
            &mut message,
            FLEN_ERRMSG,
            "Specified column number is out of range: {} (ffgcxuk)",
            colnum
        );
        ffpmsg_slice(&message);
        int_snprintf!(
            &mut message,
            FLEN_ERRMSG,
            "  There are {} columns in this table.",
            fptr.Fptr.tfield
        );
        ffpmsg_slice(&message);

        *status = BAD_COL_NUM;
        return *status;
    }

    let c = fptr.Fptr.get_tableptr_as_slice(); /* point to first column */
    let ci = (colnum - 1) as usize; /* offset to correct column structure */
    let tcol = c[ci];

    if tcol.tdatatype.abs() > TBYTE {
        ffpmsg_str("Can only read bits from X or B type columns. (ffgcxuk)");
        *status = NOT_LOGICAL_COL; /* not correct datatype column */
        return *status;
    }

    let firstbyte: c_int = ((input_first_bit - 1) / 8 + 1) as c_int;
    let lastbyte: c_int = ((input_first_bit + input_nbits as c_long - 2) / 8 + 1) as c_int;
    let nbytes: c_int = lastbyte - firstbyte + 1;

    if tcol.tdatatype == TBIT
        && input_first_bit + input_nbits as c_long - 1 > tcol.trepeat as c_long
    {
        ffpmsg_str("Too many bits. Tried to read past width of column (ffgcxuk)");
        *status = BAD_ELEM_NUM;
        return *status;
    } else if tcol.tdatatype == TBYTE && lastbyte as c_long > tcol.trepeat as c_long {
        ffpmsg_str("Too many bits. Tried to read past width of column (ffgcxuk)");
        *status = BAD_ELEM_NUM;
        return *status;
    }

    for ii in 0..nrows {
        /* read the relevant bytes from the row */
        if ffgcvuk_safe(
            fptr,
            colnum,
            firstrow + ii,
            firstbyte as LONGLONG,
            nbytes as LONGLONG,
            0,
            &mut colbyte,
            None,
            status,
        ) > 0
        {
            ffpmsg_str("Error reading bytes from column (ffgcxuk)");
            return *status;
        }

        let mut firstbit: c_int = ((input_first_bit - 1) % 8) as c_int; /* modulus operator */
        let mut nbits: c_int = input_nbits;

        array[ii as usize] = 0;

        /* select and shift the bits from each byte into the output word */
        while nbits != 0 {
            let bytenum: c_int = firstbit / 8;

            let startbit: c_int = firstbit % 8;
            let numbits: c_int = cmp::min(nbits, 8 - startbit);
            let endbit: c_int = startbit + numbits - 1;

            let rshift: c_int = 7 - endbit;
            let lshift: c_int = nbits - numbits;

            array[ii as usize] |= (colbyte[bytenum as usize] >> rshift) << lshift;

            nbits -= numbits;
            firstbit += numbits;
        }
    }

    *status
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::aliases::rust_api::*;
    use crate::fitsio::{
        BAD_COL_NUM, BAD_ELEM_NUM, BAD_ROW_NUM, BINARY_TBL, BYTE_IMG, LONGLONG, NOT_BTABLE,
        NOT_LOGICAL_COL, READONLY, TLOGICAL, fitsfile,
    };
    use crate::helpers::testhelpers::{to_buf, with_temp_file};
    use libc::{c_char, c_int, c_long};

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
    fn test_read_logical_column() {
        // Test ffgcvl - read logical values with null substitution.
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [c_char; 5] = [1, 0, 1, 0, 1];

            let mut f = make_table(&name, "LOGCOL", "1L", 5, &mut status);
            fits_write_col_log(f.as_deref_mut().unwrap(), 1, 1, 1, 5, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0, "setup");

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0 as c_char; 5];
            let mut anynull = -1;
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
            assert_eq!(result, [1, 0, 1, 0, 1]);
            assert_eq!(anynull, 0);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_read_logical_deprecated() {
        // Test ffgcl - deprecated logical read function.
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [c_char; 3] = [1, 0, 1];

            let mut f = make_table(&name, "LOGCOL", "1L", 3, &mut status);
            fits_write_col_log(f.as_deref_mut().unwrap(), 1, 1, 1, 3, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0, "setup");

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0 as c_char; 3];
            #[allow(deprecated)]
            ffgcl_safe(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                3,
                &mut result,
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result, [1, 0, 1]);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_read_logical_with_null_flags() {
        // Test ffgcfl - read logical with null flag array.
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [c_char; 2] = [1, 0];

            let mut f = make_table(&name, "LOGCOL", "1L", 3, &mut status);
            // Write T, F, then leave row 3 as undefined.
            fits_write_col_log(f.as_deref_mut().unwrap(), 1, 1, 1, 2, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0, "setup");

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0 as c_char; 3];
            let mut nularray = [0 as c_char; 3];
            let mut anynull = -1;
            fits_read_colnull_log(
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
            assert_eq!(result[0], 1);
            assert_eq!(result[1], 0);
            assert_eq!(nularray[2], 1); // Row 3 is undefined.
            assert_eq!(anynull, 1);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_read_logical_vector() {
        // Test reading logical vector column.
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [c_char; 4] = [1, 0, 1, 0];

            let mut f = make_table(&name, "LOGVEC", "4L", 1, &mut status);
            fits_write_col_log(f.as_deref_mut().unwrap(), 1, 1, 1, 4, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0, "setup");

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0 as c_char; 4];
            let mut anynull = -1;
            fits_read_col_log(
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
            assert_eq!(result, [1, 0, 1, 0]);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_read_logical_multiple_rows() {
        // Test reading logical values across multiple rows.
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [c_char; 6] = [1, 0, 0, 1, 1, 1];

            let mut f = make_table(&name, "LOGCOL", "2L", 3, &mut status);
            fits_write_col_log(f.as_deref_mut().unwrap(), 1, 1, 1, 6, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0, "setup");

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0 as c_char; 6];
            let mut anynull = -1;
            fits_read_col_log(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                6,
                0,
                &mut result,
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            // Row 1: T, F; Row 2: F, T; Row 3: T, T.
            assert_eq!(result, [1, 0, 0, 1, 1, 1]);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_read_logical_with_nulval() {
        // Test null value substitution with ffgcvl.
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [c_char; 1] = [1];

            let mut f = make_table(&name, "LOGCOL", "1L", 2, &mut status);
            fits_write_col_log(f.as_deref_mut().unwrap(), 1, 1, 1, 1, &data, &mut status);
            // Row 2 is undefined.
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0, "setup");

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            // Use 99 as the null substitute value.
            let mut result = [0 as c_char; 2];
            let mut anynull = -1;
            fits_read_col_log(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                2,
                99,
                &mut result,
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[0], 1);
            assert_eq!(result[1], 99); // Null replaced with 99.
            assert_eq!(anynull, 1);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_read_bits_as_logical() {
        // Test ffgcx - read bits from X column as logical array.
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let bits: [c_char; 8] = [1, 0, 1, 0, 1, 0, 1, 0]; // 8 bits per row.

            let mut f = make_table(&name, "BITCOL", "8X", 1, &mut status);
            fits_write_col_bit(f.as_deref_mut().unwrap(), 1, 1, 1, 8, &bits, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0, "setup");

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
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
            assert_eq!(result[0], 1);
            assert_eq!(result[1], 0);
            assert_eq!(result[2], 1);
            assert_eq!(result[3], 0);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_read_bits_from_byte_column() {
        // Test ffgcx with B column (byte).
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [u8; 1] = [0xAA]; // 10101010 in binary.

            let mut f = make_table(&name, "BYTECOL", "1B", 1, &mut status);
            fits_write_col_byt(f.as_deref_mut().unwrap(), 1, 1, 1, 1, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0, "setup");

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
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
            // 0xAA = 10101010: bit 1=1, bit 2=0, bit 3=1, bit 4=0, etc.
            assert_eq!(result, [1, 0, 1, 0, 1, 0, 1, 0]);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_read_bits_spanning_rows() {
        // Test ffgcx reading bits that span multiple rows.
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let bits_row1: [c_char; 4] = [1, 1, 0, 0]; // 4 bits per row.
            let bits_row2: [c_char; 4] = [0, 0, 1, 1];

            let mut f = make_table(&name, "BITCOL", "4X", 2, &mut status);
            fits_write_col_bit(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                4,
                &bits_row1,
                &mut status,
            );
            fits_write_col_bit(
                f.as_deref_mut().unwrap(),
                1,
                2,
                1,
                4,
                &bits_row2,
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0, "setup");

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0 as c_char; 8];
            // Read 4 bits from row 1, then 4 bits from row 2.
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
            assert_eq!(result[0], 1);
            assert_eq!(result[1], 1);
            assert_eq!(result[2], 0);
            assert_eq!(result[3], 0);
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
            assert_eq!(result[4], 0);
            assert_eq!(result[5], 0);
            assert_eq!(result[6], 1);
            assert_eq!(result[7], 1);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_read_bits_as_ushort() {
        // Test ffgcxui - read bits as unsigned short.
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);

            // Set up bit pattern: 0xFFFF = all 1s, 0x0000 = all 0s.
            let bits_all_one = [1 as c_char; 16];
            let bits_all_zero = [0 as c_char; 16];

            let mut f = make_table(&name, "BITCOL", "16X", 2, &mut status);
            fits_write_col_bit(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                16,
                &bits_all_one,
                &mut status,
            );
            fits_write_col_bit(
                f.as_deref_mut().unwrap(),
                1,
                2,
                1,
                16,
                &bits_all_zero,
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0, "setup");

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0u16; 2];
            fits_read_col_bit_usht(
                f.as_deref_mut().unwrap(),
                1,
                1,
                2,
                1,
                16,
                &mut result,
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[0], 0xFFFF);
            assert_eq!(result[1], 0x0000);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_read_bits_as_ushort_partial() {
        // Test ffgcxui - read partial bits as unsigned short.
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let bits: [c_char; 16] = [1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];

            let mut f = make_table(&name, "BITCOL", "16X", 1, &mut status);
            fits_write_col_bit(f.as_deref_mut().unwrap(), 1, 1, 1, 16, &bits, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0, "setup");

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            // Read only first 4 bits: 1010 = 10 in decimal.
            let mut result = [0u16; 1];
            fits_read_col_bit_usht(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                1,
                4,
                &mut result,
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[0], 10);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_read_bits_as_uint() {
        // Test ffgcxuk - read bits as unsigned int.
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);

            // Set up bit pattern: high byte = 0xFF, rest = 0.
            let mut bits = [0 as c_char; 32];
            for b in bits.iter_mut().take(8) {
                *b = 1;
            }

            let mut f = make_table(&name, "BITCOL", "32X", 1, &mut status);
            fits_write_col_bit(f.as_deref_mut().unwrap(), 1, 1, 1, 32, &bits, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0, "setup");

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0u32; 1];
            fits_read_col_bit_uint(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                1,
                32,
                &mut result,
                &mut status,
            );
            assert_eq!(status, 0);
            // 11111111 00000000 00000000 00000000 = 0xFF000000.
            assert_eq!(result[0], 0xFF000000);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_read_bits_as_uint_partial() {
        // Test ffgcxuk - read partial bits as unsigned int.
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let bits: [c_char; 32] = [
                1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0,
            ];

            let mut f = make_table(&name, "BITCOL", "32X", 1, &mut status);
            fits_write_col_bit(f.as_deref_mut().unwrap(), 1, 1, 1, 32, &bits, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0, "setup");

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            // Read only first 8 bits: 11110000 = 240 in decimal.
            let mut result = [0u32; 1];
            fits_read_col_bit_uint(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                1,
                8,
                &mut result,
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[0], 240);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_read_bits_bad_row() {
        // Test ffgcx with bad row number.
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);

            let mut f = make_table(&name, "BITCOL", "8X", 1, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0, "setup");

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0 as c_char; 8];
            fits_read_col_bit(
                f.as_deref_mut().unwrap(),
                1,
                0,
                1,
                8,
                &mut result,
                &mut status,
            ); // Row 0 is invalid.
            assert_eq!(status, BAD_ROW_NUM);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_read_bits_bad_elem() {
        // Test ffgcx with bad element number.
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);

            let mut f = make_table(&name, "BITCOL", "8X", 1, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0, "setup");

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0 as c_char; 8];
            fits_read_col_bit(
                f.as_deref_mut().unwrap(),
                1,
                1,
                0,
                8,
                &mut result,
                &mut status,
            ); // Bit 0 is invalid.
            assert_eq!(status, BAD_ELEM_NUM);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_read_bits_ushort_bad_row() {
        // Test ffgcxui with bad row number.
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);

            let mut f = make_table(&name, "BITCOL", "16X", 1, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0, "setup");

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0u16; 1];
            fits_read_col_bit_usht(
                f.as_deref_mut().unwrap(),
                1,
                0,
                1,
                1,
                16,
                &mut result,
                &mut status,
            ); // Row 0 is invalid.
            assert_eq!(status, BAD_ROW_NUM);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_read_bits_ushort_bad_bit() {
        // Test ffgcxui with bad first bit number.
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);

            let mut f = make_table(&name, "BITCOL", "16X", 1, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0, "setup");

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0u16; 1];
            fits_read_col_bit_usht(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                0,
                16,
                &mut result,
                &mut status,
            ); // Bit 0 is invalid.
            assert_eq!(status, BAD_ELEM_NUM);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_read_bits_ushort_too_many() {
        // Test ffgcxui with too many bits (>16).
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);

            let mut f = make_table(&name, "BITCOL", "32X", 1, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0, "setup");

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0u16; 1];
            fits_read_col_bit_usht(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                1,
                17,
                &mut result,
                &mut status,
            ); // > 16 bits.
            assert_eq!(status, BAD_ELEM_NUM);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_read_bits_uint_bad_row() {
        // Test ffgcxuk with bad row number.
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);

            let mut f = make_table(&name, "BITCOL", "32X", 1, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0, "setup");

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0u32; 1];
            fits_read_col_bit_uint(
                f.as_deref_mut().unwrap(),
                1,
                0,
                1,
                1,
                32,
                &mut result,
                &mut status,
            ); // Row 0 is invalid.
            assert_eq!(status, BAD_ROW_NUM);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_read_bits_uint_bad_bit() {
        // Test ffgcxuk with bad first bit number.
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);

            let mut f = make_table(&name, "BITCOL", "32X", 1, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0, "setup");

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0u32; 1];
            fits_read_col_bit_uint(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                0,
                32,
                &mut result,
                &mut status,
            ); // Bit 0 is invalid.
            assert_eq!(status, BAD_ELEM_NUM);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_read_bits_uint_too_many() {
        // Test ffgcxuk with too many bits (>32).
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);

            let mut f = make_table(&name, "BITCOL", "64X", 1, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0, "setup");

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0u32; 1];
            fits_read_col_bit_uint(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                1,
                33,
                &mut result,
                &mut status,
            ); // > 32 bits.
            assert_eq!(status, BAD_ELEM_NUM);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_read_logical_not_logical_col() {
        // Test ffgcvl with non-logical column returns error.
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [c_long; 1] = [1];

            let mut f = make_table(&name, "INTCOL", "1J", 1, &mut status);
            fits_write_col_lng(f.as_deref_mut().unwrap(), 1, 1, 1, 1, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0, "setup");

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0 as c_char; 1];
            let mut anynull = -1;
            fits_read_col_log(
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
            assert_eq!(status, NOT_LOGICAL_COL);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_read_bits_not_bit_col() {
        // Test ffgcx with non-bit/byte column returns error.
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [c_long; 1] = [1];

            let mut f = make_table(&name, "INTCOL", "1J", 1, &mut status);
            fits_write_col_lng(f.as_deref_mut().unwrap(), 1, 1, 1, 1, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0, "setup");

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
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
            assert_eq!(status, NOT_LOGICAL_COL);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_read_bits_ushort_not_bit_col() {
        // Test ffgcxui with non-bit/byte column returns error.
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [c_long; 1] = [1];

            let mut f = make_table(&name, "INTCOL", "1J", 1, &mut status);
            fits_write_col_lng(f.as_deref_mut().unwrap(), 1, 1, 1, 1, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0, "setup");

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0u16; 1];
            fits_read_col_bit_usht(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                1,
                16,
                &mut result,
                &mut status,
            );
            assert_eq!(status, NOT_LOGICAL_COL);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_read_bits_uint_not_bit_col() {
        // Test ffgcxuk with non-bit/byte column returns error.
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [c_long; 1] = [1];

            let mut f = make_table(&name, "INTCOL", "1J", 1, &mut status);
            fits_write_col_lng(f.as_deref_mut().unwrap(), 1, 1, 1, 1, &data, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0, "setup");

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0u32; 1];
            fits_read_col_bit_uint(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                1,
                32,
                &mut result,
                &mut status,
            );
            assert_eq!(status, NOT_LOGICAL_COL);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_read_bits_ushort_bad_col() {
        // Test ffgcxui with bad column number.
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);

            let mut f = make_table(&name, "BITCOL", "16X", 1, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0, "setup");

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0u16; 1];
            fits_read_col_bit_usht(
                f.as_deref_mut().unwrap(),
                99,
                1,
                1,
                1,
                16,
                &mut result,
                &mut status,
            ); // Column 99 invalid.
            assert_eq!(status, BAD_COL_NUM);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_read_bits_uint_bad_col() {
        // Test ffgcxuk with bad column number.
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);

            let mut f = make_table(&name, "BITCOL", "32X", 1, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0, "setup");

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0u32; 1];
            fits_read_col_bit_uint(
                f.as_deref_mut().unwrap(),
                99,
                1,
                1,
                1,
                32,
                &mut result,
                &mut status,
            ); // Column 99 invalid.
            assert_eq!(status, BAD_COL_NUM);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_read_bits_ushort_not_bintable() {
        // Test ffgcxui when not in binary table.
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 1] = [10];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 1, &naxes, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0, "setup");

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut result = [0u16; 1];
            fits_read_col_bit_usht(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                1,
                8,
                &mut result,
                &mut status,
            );
            assert_eq!(status, NOT_BTABLE);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_read_bits_uint_not_bintable() {
        // Test ffgcxuk when not in binary table.
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 1] = [10];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 1, &naxes, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0, "setup");

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut result = [0u32; 1];
            fits_read_col_bit_uint(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                1,
                8,
                &mut result,
                &mut status,
            );
            assert_eq!(status, NOT_BTABLE);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_read_bits_overflow_x_col() {
        // Test ffgcxui reading past end of X column.
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);

            let mut f = make_table(&name, "BITCOL", "8X", 1, &mut status); // Only 8 bits.
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0, "setup");

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            // Try to read 16 bits from an 8-bit column.
            let mut result = [0u16; 1];
            fits_read_col_bit_usht(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                1,
                16,
                &mut result,
                &mut status,
            );
            assert_eq!(status, BAD_ELEM_NUM);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_read_bits_overflow_b_col() {
        // Test ffgcxui reading past end of B column.
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);

            let mut f = make_table(&name, "BYTECOL", "1B", 1, &mut status); // Only 8 bits (1 byte).
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0, "setup");

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            // Try to read 16 bits starting at bit 1 from a 1-byte column.
            let mut result = [0u16; 1];
            fits_read_col_bit_usht(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                1,
                16,
                &mut result,
                &mut status,
            );
            assert_eq!(status, BAD_ELEM_NUM);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_read_illegal_logical_char() {
        // Test reading a logical column that contains illegal characters.
        // Write raw bytes to a logical column to simulate corruption.
        // This exercises the "other illegal character" path in ffgcll.
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let illegal_data: [u8; 4] = [b'T', b'X', b'F', 1]; // 'X' is illegal, 1 is special case

            let mut f = make_table(&name, "LOGCOL", "4L", 1, &mut status);
            // Write raw bytes directly to bypass validation
            fits_write_tblbytes(
                f.as_deref_mut().unwrap(),
                1,
                1,
                4,
                &illegal_data,
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0, "setup");

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0u8; 4];
            let mut anynull = -1;
            fits_read_col(
                f.as_deref_mut().unwrap(),
                TLOGICAL,
                1,
                1,
                1,
                4,
                None,
                &mut result,
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            // 'T' should become 1, 'X' should return as char value, 'F' should become 0
            // The value 1 is a special case - returns 49 (ASCII '1')
            assert_eq!(result[0], 1); // 'T'
            assert_eq!(result[1], b'X'); // illegal char returned as-is
            assert_eq!(result[2], 0); // 'F'
            assert_eq!(result[3], 49); // byte value 1 becomes ASCII '1' (49)
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_read_bits_varlen_overflow() {
        // Test BAD_ELEM_NUM error for variable length bit column.
        // Try to read more bits than exist in a variable length descriptor.
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let bits: [c_char; 4] = [1, 0, 1, 0];

            let mut f = make_table(&name, "VARBITS", "1PX", 1, &mut status);
            fits_write_col_bit(f.as_deref_mut().unwrap(), 1, 1, 1, 4, &bits, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0, "setup");

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            // Try to read 100 bits but we only wrote 4
            let mut result = [0 as c_char; 10];
            fits_read_col_bit(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                100,
                &mut result,
                &mut status,
            );
            assert_eq!(status, BAD_ELEM_NUM);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_read_bits_uint_overflow_x_col() {
        // Test ffgcxuk reading past end of X column (TBIT).
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);

            let mut f = make_table(&name, "BITCOL", "8X", 1, &mut status); // Only 8 bits.
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0, "setup");

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            // Try to read 32 bits from an 8-bit column.
            let mut result = [0u32; 1];
            fits_read_col_bit_uint(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                1,
                32,
                &mut result,
                &mut status,
            );
            assert_eq!(status, BAD_ELEM_NUM);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_read_bits_uint_overflow_b_col() {
        // Test ffgcxuk reading past end of B column (TBYTE).
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);

            let mut f = make_table(&name, "BYTECOL", "1B", 1, &mut status); // Only 8 bits (1 byte).
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0, "setup");

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            // Try to read 32 bits starting at bit 1 from a 1-byte column.
            let mut result = [0u32; 1];
            fits_read_col_bit_uint(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                1,
                32,
                &mut result,
                &mut status,
            );
            assert_eq!(status, BAD_ELEM_NUM);
            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }
}
