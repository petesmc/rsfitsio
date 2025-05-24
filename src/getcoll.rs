/*  This file, getcoll.c, contains routines that read data elements from   */
/*  a FITS image or table, with logical datatype.                          */

/*  The FITSIO software was written by William Pence at the High Energy    */
/*  Astrophysic Science Archive Research Center (HEASARC) at the NASA      */
/*  Goddard Space Flight Center.                                           */

use core::slice;
use std::cmp;

use crate::c_types::{c_char, c_int, c_long, c_uchar};

use bytemuck::cast_slice_mut;

use crate::NullCheckType;
use crate::fitscore::{ffgcprll, ffgdes_safe, ffmahd_safe, ffpmsg_slice, ffrdef_safe};
use crate::fitsio::*;
use crate::fitsio2::*;
use crate::{buffers::*, int_snprintf};

/*--------------------------------------------------------------------------*/
/// Read an array of logical values from a column in the current FITS HDU.
///
/// Any undefined pixels will be set equal to the value of 'nulval' unless
/// nulval = 0 in which case no checks for undefined pixels will be made.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgcvl(
    fptr: *mut fitsfile, /* I - FITS file pointer                       */
    colnum: c_int,       /* I - number of column to read (1 = 1st col)  */
    firstrow: LONGLONG,  /* I - first row to read (1 = 1st row)         */
    firstelem: LONGLONG, /* I - first vector element to read (1 = 1st)  */
    nelem: LONGLONG,     /* I - number of values to read                */
    nulval: c_char,      /* I - value for null pixels                   */
    array: *mut c_char,  /* O - array of values                         */
    anynul: *mut c_int,  /* O - set to 1 if any values are null; else 0 */
    status: *mut c_int,  /* IO - error status                           */
) -> c_int {
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

/*--------------------------------------------------------------------------*/
/// Read an array of logical values from a column in the current FITS HDU.
///
/// Any undefined pixels will be set equal to the value of 'nulval' unless
/// nulval = 0 in which case no checks for undefined pixels will be made.
pub fn ffgcvl_safe(
    fptr: &mut fitsfile,        /* I - FITS file pointer                       */
    colnum: c_int,              /* I - number of column to read (1 = 1st col)  */
    firstrow: LONGLONG,         /* I - first row to read (1 = 1st row)         */
    firstelem: LONGLONG,        /* I - first vector element to read (1 = 1st)  */
    nelem: LONGLONG,            /* I - number of values to read                */
    nulval: c_char,             /* I - value for null pixels                   */
    array: &mut [c_char],       /* O - array of values                         */
    anynul: Option<&mut c_int>, /* O - set to 1 if any values are null; else 0 */
    status: &mut c_int,         /* IO - error status                           */
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

/*--------------------------------------------------------------------------*/
/// Read an array of logical values from a column in the current FITS HDU.
///
/// !!!! THIS ROUTINE IS DEPRECATED AND SHOULD NOT BE USED !!!!!!
///           !!!! USE ffgcvl INSTEAD  !!!!!!
/// No checking for null values will be performed.
// #[deprecated]
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgcl(
    fptr: *mut fitsfile, /* I - FITS file pointer                       */
    colnum: c_int,       /* I - number of column to read (1 = 1st col)  */
    firstrow: LONGLONG,  /* I - first row to read (1 = 1st row)         */
    firstelem: LONGLONG, /* I - first vector element to read (1 = 1st)  */
    nelem: LONGLONG,     /* I - number of values to read                */
    array: *mut c_char,  /* O - array of values                         */
    status: *mut c_int,  /* IO - error status                           */
) -> c_int {
    unsafe {
        let nulval: c_char = 0;
        let mut anynul: c_int = 0;

        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        let array = slice::from_raw_parts_mut(array, nelem as usize);

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
        );

        *status
    }
}

/*--------------------------------------------------------------------------*/
/// Read an array of logical values from a column in the current FITS HDU.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgcfl(
    fptr: *mut fitsfile,   /* I - FITS file pointer                       */
    colnum: c_int,         /* I - number of column to read (1 = 1st col)  */
    firstrow: LONGLONG,    /* I - first row to read (1 = 1st row)         */
    firstelem: LONGLONG,   /* I - first vector element to read (1 = 1st)  */
    nelem: LONGLONG,       /* I - number of values to read                */
    array: *mut c_char,    /* O - array of values                         */
    nularray: *mut c_char, /* O - array of flags = 1 if nultyp = 2        */
    anynul: *mut c_int,    /* O - set to 1 if any values are null; else 0 */
    status: *mut c_int,    /* IO - error status                           */
) -> c_int {
    unsafe {
        let nulval: c_char = 0;

        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        let array = slice::from_raw_parts_mut(array, nelem as usize);

        let anynul = anynul.as_mut();

        let nularray = slice::from_raw_parts_mut(nularray, nelem as usize);

        ffgcll(
            fptr,
            colnum,
            firstrow as LONGLONG,
            firstelem as LONGLONG,
            nelem as LONGLONG,
            NullCheckType::SetNullArray,
            nulval,
            array,
            nularray,
            anynul,
            status,
        );

        *status
    }
}

/*--------------------------------------------------------------------------*/
/// Read an array of logical values from a column in the current FITS HDU.
pub(crate) fn ffgcll(
    fptr: &mut fitsfile,   /* I - FITS file pointer                       */
    colnum: c_int,         /* I - number of column to read (1 = 1st col)  */
    firstrow: LONGLONG,    /* I - first row to read (1 = 1st row)         */
    firstelem: LONGLONG,   /* I - first vector element to read (1 = 1st)  */
    nelem: LONGLONG,       /* I - number of values to read                */
    nultyp: NullCheckType, /* I - null value handling code:               */
    /*     1: set undefined pixels = nulval        */
    /*     2: set nularray=1 for undefined pixels  */
    nulval: c_char,                 /* I - value for null pixels if nultyp = 1     */
    array: &mut [c_char],           /* O - array of values                         */
    nularray: &mut [c_char],        /* O - array of flags = 1 if nultyp = 2        */
    mut anynul: Option<&mut c_int>, /* O - set to 1 if any values are null; else 0 */
    status: &mut c_int,             /* IO - error status                           */
) -> c_int {
    let mut dtemp: f64 = 0.0;
    let mut tcode: c_int = 0;
    let mut maxelem: c_int = 0;
    let mut hdutype: c_int = 0;
    let ii: c_int = 0;
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
    /*---------------------------------------------------*/
    /*  Check input and get parameters about the column: */
    /*---------------------------------------------------*/
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
    /*------------------------------------------------------------------*/
    /*  Decide whether to check for null values in the input FITS file: */
    /*------------------------------------------------------------------*/
    nulcheck = nultyp; /* by default, check for null values in the FITS file */

    if nultyp == NullCheckType::SetPixel && nulval == 0 {
        nulcheck = NullCheckType::None; /* calling routine does not want to check for nulls */
    }
    /*---------------------------------------------------------------------*/
    /*  Now read the logical values from the FITS column.                  */
    /*---------------------------------------------------------------------*/

    remain = nelem; /* remaining number of values to read */
    next = 0; /* next element in array to be read   */
    rownum = 0; /* row number, relative to firstrow   */
    ntodo = remain as c_long; /* max number of elements to read at one time */

    while ntodo > 0 {
        /*
         limit the number of pixels to read at one time to the number that
         remain in the current vector.
        */
        ntodo = cmp::min(ntodo as LONGLONG, maxelem as LONGLONG) as c_long;
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

        /*--------------------------------------------*/
        /*  increment the counters for the next loop  */
        /*--------------------------------------------*/
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

/*--------------------------------------------------------------------------*/
/// read an array of logical values from a specified bit or byte
/// column of the binary table.    larray is set = TRUE, if the corresponding
/// bit = 1, otherwise it is set to FALSE.
/// The binary table column being read from must have datatype 'B' or 'X'.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffgcx(
    fptr: *mut fitsfile, /* I - FITS file pointer                       */
    colnum: c_int,       /* I - number of column to write (1 = 1st col) */
    frow: LONGLONG,      /* I - first row to write (1 = 1st row)        */
    fbit: LONGLONG,      /* I - first bit to write (1 = 1st)            */
    nbit: LONGLONG,      /* I - number of bits to write                 */
    larray: *mut c_char, /* O - array of logicals corresponding to bits */
    status: *mut c_int,  /* IO - error status                           */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        let larray = slice::from_raw_parts_mut(larray, nbit as usize);

        ffgcx_safe(fptr, colnum, frow, fbit, nbit, larray, status)
    }
}

/*--------------------------------------------------------------------------*/
/// read an array of logical values from a specified bit or byte
/// column of the binary table.    larray is set = TRUE, if the corresponding
/// bit = 1, otherwise it is set to FALSE.
/// The binary table column being read from must have datatype 'B' or 'X'.
pub fn ffgcx_safe(
    fptr: &mut fitsfile,   /* I - FITS file pointer                       */
    colnum: c_int,         /* I - number of column to write (1 = 1st col) */
    frow: LONGLONG,        /* I - first row to write (1 = 1st row)        */
    fbit: LONGLONG,        /* I - first bit to write (1 = 1st)            */
    nbit: LONGLONG,        /* I - number of bits to write                 */
    larray: &mut [c_char], /* O - array of logicals corresponding to bits */
    status: &mut c_int,    /* IO - error status                           */
) -> c_int {
    let mut bstart: LONGLONG = 0;
    let mut offset: c_long = 0;
    let mut ndone: c_long = 0;
    let ii: c_long = 0;
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
