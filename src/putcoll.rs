/*  This file, putcoll.c, contains routines that write data elements to    */
/*  a FITS image or table, with logical datatype.                          */

/*  The FITSIO software was written by William Pence at the High Energy    */
/*  Astrophysic Science Archive Research Center (HEASARC) at the NASA      */
/*  Goddard Space Flight Center.                                           */

use core::slice;
use std::cmp;

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

/*--------------------------------------------------------------------------*/
/// Write an array of logical values to a column in the current FITS HDU.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffpcll(
    fptr: *mut fitsfile,  /* I - FITS file pointer                       */
    colnum: c_int,        /* I - number of column to write (1 = 1st col) */
    firstrow: LONGLONG,   /* I - first row to write (1 = 1st row)        */
    firstelem: LONGLONG,  /* I - first vector element to write (1 = 1st) */
    nelem: LONGLONG,      /* I - number of values to write               */
    array: *const c_char, /* I - array of values to write                */
    status: *mut c_int,   /* IO - error status                           */
) -> c_int {
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

/*--------------------------------------------------------------------------*/
/// Write an array of logical values to a column in the current FITS HDU.
pub fn ffpcll_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer                       */
    colnum: c_int,       /* I - number of column to write (1 = 1st col) */
    firstrow: LONGLONG,  /* I - first row to write (1 = 1st row)        */
    firstelem: LONGLONG, /* I - first vector element to write (1 = 1st) */
    nelem: LONGLONG,     /* I - number of values to write               */
    array: &[c_char],    /* I - array of values to write                */
    status: &mut c_int,  /* IO - error status                           */
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

    if tcode != TLOGICAL {
        *status = NOT_LOGICAL_COL;
        return *status;
    }

    /*---------------------------------------------------------------------*/
    /*  Now write the logical values one at a time to the FITS column.     */
    /*---------------------------------------------------------------------*/
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

        /*--------------------------------------------*/
        /*  increment the counters for the next loop  */
        /*--------------------------------------------*/
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

/*--------------------------------------------------------------------------*/
/// Write an array of elements to the specified column of a table.  Any input
/// pixels flagged as null will be replaced by the appropriate
/// null value in the output FITS file.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffpcnl(
    fptr: *mut fitsfile,  /* I - FITS file pointer                       */
    colnum: c_int,        /* I - number of column to write (1 = 1st col) */
    firstrow: LONGLONG,   /* I - first row to write (1 = 1st row)        */
    firstelem: LONGLONG,  /* I - first vector element to write (1 = 1st) */
    nelem: LONGLONG,      /* I - number of values to write               */
    array: *const c_char, /* I - array of values to write                */
    nulvalue: c_char,     /* I - array flagging undefined pixels if true */
    status: *mut c_int,   /* IO - error status                           */
) -> c_int {
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

/*--------------------------------------------------------------------------*/
/// Write an array of elements to the specified column of a table.  Any input
/// pixels flagged as null will be replaced by the appropriate
/// null value in the output FITS file.
pub fn ffpcnl_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer                       */
    colnum: c_int,       /* I - number of column to write (1 = 1st col) */
    firstrow: LONGLONG,  /* I - first row to write (1 = 1st row)        */
    firstelem: LONGLONG, /* I - first vector element to write (1 = 1st) */
    nelem: LONGLONG,     /* I - number of values to write               */
    array: &[c_char],    /* I - array of values to write                */
    nulvalue: c_char,    /* I - array flagging undefined pixels if true */
    status: &mut c_int,  /* IO - error status                           */
) -> c_int {
    let mut ngood: LONGLONG = 0;
    let mut nbad: LONGLONG = 0;
    let ii: LONGLONG = 0;
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

    let colptr = fptr.Fptr.tableptr; /* set pointer to first column */
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

    for ii in 0..(nelem as usize) {
        if array[ii] != nulvalue {
            /* is this a good pixel? */

            if nbad != 0 {
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
            if ngood != 0
            /* write previous string of good pixels */
            {
                fstelm = (ii as LONGLONG) - ngood + first; /* absolute element number */
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

/*--------------------------------------------------------------------------*/
/// write an array of logical values to a specified bit or byte
/// column of the binary table.   If larray is TRUE, then the corresponding
/// bit is set to 1, otherwise it is set to 0.
/// The binary table column being written to must have datatype 'B' or 'X'.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffpclx(
    fptr: *mut fitsfile,   /* I - FITS file pointer                       */
    colnum: c_int,         /* I - number of column to write (1 = 1st col) */
    frow: LONGLONG,        /* I - first row to write (1 = 1st row)        */
    fbit: c_long,          /* I - first bit to write (1 = 1st)            */
    nbit: c_long,          /* I - number of bits to write                 */
    larray: *const c_char, /* I - array of logicals corresponding to bits */
    status: *mut c_int,    /* IO - error status                           */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        let larray = slice::from_raw_parts(larray, nbit as usize);

        ffpclx_safe(fptr, colnum, frow, fbit, nbit, larray, status)
    }
}

/*--------------------------------------------------------------------------*/
/// write an array of logical values to a specified bit or byte
/// column of the binary table.   If larray is TRUE, then the corresponding
/// bit is set to 1, otherwise it is set to 0.
/// The binary table column being written to must have datatype 'B' or 'X'.
pub fn ffpclx_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer                       */
    colnum: c_int,       /* I - number of column to write (1 = 1st col) */
    frow: LONGLONG,      /* I - first row to write (1 = 1st row)        */
    fbit: c_long,        /* I - first bit to write (1 = 1st)            */
    nbit: c_long,        /* I - number of bits to write                 */
    larray: &[c_char],   /* I - array of logicals corresponding to bits */
    status: &mut c_int,  /* IO - error status                           */
) -> c_int {
    let mut offset: LONGLONG = 0;
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
    let ii: c_long = 0;
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
    offset = fptr.Fptr.heapsize;

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

    let colptr = fptr.Fptr.tableptr; /* set pointer to first column */
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
