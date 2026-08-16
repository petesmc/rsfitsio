/* Transpiled from cfitsio/utilities/fvrf_data.c

The CFITSIO iterator hands its work function a `void *userPointer' and an
array of iteratorCol whose `array' member is a raw buffer, so test_data /
iterdata keep the C's pointer arithmetic inside small unsafe blocks.  All
the other column reads use the typed safe wrappers.  */

use std::cell::{Cell, RefCell};

use bytemuck::cast_slice;
use rsfitsio::aliases::rust_api::{
    fits_get_coltype, fits_get_num_rows, fits_get_num_rowsll, fits_get_rowsize, fits_iterate_data,
    fits_read_col_log, fits_read_col_str, fits_read_descriptll, fits_read_key_lng,
    fits_read_key_lnglng, fits_read_key_str, fits_read_tblbytes, fits_verify_chksum,
};
use rsfitsio::buffers::ffgtbb_safe;
use std::ffi::CStr;

use rsfitsio::c_types::{c_char, c_int, c_long, c_uchar, c_void};
use rsfitsio::cs;
use rsfitsio::fitscore::{ffasfm_safe, ffcdfl_safe};
use rsfitsio::fitsio::{
    ASCII_TBL, BINARY_TBL, FLEN_VALUE, LONGLONG, TBIT, TBYTE, TCOMPLEX, TDBLCOMPLEX, TDOUBLE,
    TFLOAT, TLOGICAL, TLONG, TSHORT, TSTRING, fitsfile, iteratorCol,
};

use crate::common::*;
use crate::fvrf_head::parse_vtform;
use crate::fvrf_misc::*;
use crate::{scat, spf};

/* flag for input only iterator column (fitsio.h InputCol) */
const InputCol: c_int = 0;

struct UserIter {
    nnum: c_int,
    ncmp: c_int,
    nfloat: c_int,
    indatatyp: Vec<c_int>,
    datamax: Vec<f64>,
    datamin: Vec<f64>,
    tnull: Vec<f64>,
    mask: Vec<c_uchar>, /* for bit X column only */
    ntxt: c_int,
    out: Out,
}

/*************************************************************
*
*      test_data
*
*   Test the HDU data
*
*  This routine reads every row and column of ASCII tables to
*  verify that the values have the correct format.
*
*  This routine checks the following types of columns in binary tables:
*
*    Logical L - value must be T, F or zero
*    Bit    nX - if n != a multiple of 8, then check that fill bits = 0
*    String  A - must contain ascii text, or zero
*
*   It is impossible to write an invalid value to the other types of
*   columns in binary tables (B, I, J, K, E, D, C and M) so these
*   columns are not read.
*
*  Since it is impossible to write an invalid value in a FITS image,
*  this routine does not read the image pixels.
*
*************************************************************/
pub(crate) fn test_data(
    infits: &mut fitsfile, /* input fits file   */
    out: Out,              /* output ascii file */
    hduptr: &mut FitsHdu,  /* fits hdu pointer  */
) {
    let mut iter_col: Vec<iteratorCol> = Vec::new();
    let ncols: c_int;

    let mut nnum = 0; /* the list of the column  whose data
    type is numerical(scalar and complex) */
    let mut numlist: Vec<c_int>;
    let mut nfloat = 0; /* the list of the floating point columns in ASCII table */
    let mut floatlist: Vec<c_int>;

    let ncmp = 0; /* the list of the column  whose data
    type is numerical(scalar and complex) */
    let cmplist: Vec<c_int>;
    let mut ntxt = 0; /* the list of column  whose data type is
    string, logical, bit or complex */
    let mut txtlist: Vec<c_int>;
    let niter: c_int; /* total columns read into  the literator function */

    let mut ndesc = 0; /* the list of column which is the descriptor of
    the variable length array. */
    let mut desclist: Vec<c_int>;
    let mut isVarQFormat: Vec<c_int>; /* Format type for each of the ndesc variable-length
    columns: 0 = type 'P', 1 = type 'Q' */

    let mut rows_per_loop: c_long = 0;
    let offset: c_long;

    let mut datatype: c_int = 0;
    let mut repeat: c_long = 0;

    let mut totalrows: c_long = 0;
    let mut length: LONGLONG = 0;
    let mut toffset: LONGLONG = 0;
    let mut maxlen: Vec<c_long>;
    let mut icol: c_int;
    let mut cdata: Vec<c_char>;
    let mut dflag: Vec<c_int>;
    let lnull: c_char = 2;
    let mut anynul: c_int = 0;

    let mut rlength: c_long;
    let mut bytelength: c_long;
    let mut maxmax: c_long;

    let mut i: usize;
    let mut j: usize;
    let mut jl: c_long;
    let mut k: c_long;
    let mut status: c_int = 0;
    let mut errtmp: [c_char; 80] = [0; 80];
    let mut perbyte: Vec<c_int>;
    let mut errmes: [c_char; ERRMES_LEN] = [0; ERRMES_LEN];

    let mut naxis2: LONGLONG = 0;

    let mut largeVarLengthWarned = 0;
    let mut largeVarOffsetWarned = 0;

    if testcsum() != 0 {
        test_checksum(infits, out);
    }

    if testfill() != 0 {
        test_agap(infits, out, hduptr); /* test the bytes between the
        ascii table columns. */
        if ffcdfl_safe(infits, &mut status) != 0 {
            wrtferr_str(out, "checking data fill: ", &mut status, 1);
            status = 0;
        }
    }

    if hduptr.hdutype != ASCII_TBL && hduptr.hdutype != BINARY_TBL {
        return;
    }

    ncols = hduptr.ncols;
    if ncols <= 0 {
        return;
    }

    fits_read_key_lnglng(infits, cs!(c"NAXIS2"), &mut naxis2, None, &mut status);

    if naxis2 > 2147483647 {
        wrtout_str(
            out,
            "Cannot test data in tables with more than 2**31 (2147483647) rows.",
        );
        return;
    }

    /* separate the numerical, complex, text and
    the variable length vector columns */
    numlist = vec![0; ncols as usize];
    floatlist = vec![0; ncols as usize];
    cmplist = vec![0; ncols as usize];
    txtlist = vec![0; ncols as usize];
    desclist = vec![0; ncols as usize];
    let _ = &cmplist;

    if hduptr.hdutype == ASCII_TBL {
        /*read every column of an ASCII table */
        rows_per_loop = 0;
        for i in 0..ncols {
            if fits_get_coltype(infits, i + 1, Some(&mut datatype), None, None, &mut status) != 0 {
                spf!(errmes; "Column #", i, ": ");
                wrtferr(out, &errmes, &mut status, 2);
            }
            if datatype != TSTRING {
                numlist[nnum as usize] = i + 1;
                nnum += 1;

                if datatype > TLONG {
                    /* floating point number column */
                    floatlist[nfloat as usize] = i + 1;
                    nfloat += 1;
                }
            } else {
                txtlist[ntxt as usize] = i + 1;
                ntxt += 1;
            }
        }
    } else if hduptr.hdutype == BINARY_TBL {
        /* only check Bit, Logical and String columns in Binary tables */
        rows_per_loop = 0;
        for i in 0..ncols {
            if fits_get_coltype(
                infits,
                i + 1,
                Some(&mut datatype),
                Some(&mut repeat),
                None,
                &mut status,
            ) != 0
            {
                spf!(errmes; "Column #", i, ": ");
                wrtferr(out, &errmes, &mut status, 2);
            }

            if datatype < 0 {
                /* variable length column */
                desclist[ndesc as usize] = i + 1;
                ndesc += 1;
            } else if datatype == TBIT && (repeat % 8) != 0 {
                /* bit column that does not have a multiple of 8 bits */
                numlist[nnum as usize] = i + 1;
                nnum += 1;
            } else if datatype == TLOGICAL || datatype == TSTRING {
                txtlist[ntxt as usize] = i + 1;
                ntxt += 1;
            }
            /* ignore all other types of columns (B I J K E D C and M ) */
        }
    }

    /*  Use Iterator to read the columns that are not variable length arrays */
    /* columns from  1 to nnum are scalar numerical columns.
    columns from  nnum+1 to  nnum+ncmp are complex columns.
    columns from  nnum+ncmp are text columns */
    niter = nnum + ncmp + ntxt + nfloat;

    if niter != 0 {
        iter_col = vec![iteratorCol::default(); niter as usize];
    }

    /* fits_iter_set_by_num() just fills in these four members */
    let fptr_raw: *mut fitsfile = infits;
    for i in 0..nnum as usize {
        iter_col[i].fptr = fptr_raw;
        iter_col[i].colnum = numlist[i];
        iter_col[i].datatype = TDOUBLE;
        iter_col[i].iotype = InputCol;
    }
    for i in 0..ncmp as usize {
        let j = nnum as usize + i;
        iter_col[j].fptr = fptr_raw;
        iter_col[j].colnum = cmplist[i];
        iter_col[j].datatype = TDBLCOMPLEX;
        iter_col[j].iotype = InputCol;
    }
    for i in 0..ntxt as usize {
        let j = (nnum + ncmp) as usize + i;
        iter_col[j].fptr = fptr_raw;
        iter_col[j].colnum = txtlist[i];
        iter_col[j].datatype = 0;
        iter_col[j].iotype = InputCol;
    }
    for i in 0..nfloat as usize {
        let j = (nnum + ncmp + ntxt) as usize + i;
        iter_col[j].fptr = fptr_raw;
        iter_col[j].colnum = floatlist[i];
        iter_col[j].datatype = TSTRING;
        iter_col[j].iotype = InputCol;
    }

    offset = 0;
    let mut usrdata = UserIter {
        nnum,
        ncmp,
        nfloat,
        indatatyp: Vec::new(),
        datamax: Vec::new(),
        datamin: Vec::new(),
        tnull: Vec::new(),
        mask: Vec::new(),
        ntxt,
        out,
    };
    if nnum > 0 || ncmp > 0 {
        usrdata.datamax = vec![0.0; (nnum + ncmp) as usize];
        usrdata.datamin = vec![0.0; (nnum + ncmp) as usize];
    }
    usrdata.tnull = vec![0.0; ncols as usize];

    /* get the mask for the bit X column
        for column other than the X, it always 255
        for Column nX, it will be 000...111, where # of 0 is n%8,
        # of 1 is 8 - n%8.
    */

    if nnum > 0 {
        usrdata.mask = vec![0; nnum as usize];
        usrdata.indatatyp = vec![0; nnum as usize];
    }
    for i in 0..nnum as usize {
        let jc = iter_col[i].colnum; /* fits_iter_get_colnum() */
        if fits_get_coltype(
            infits,
            jc,
            Some(&mut datatype),
            Some(&mut repeat),
            None,
            &mut status,
        ) != 0
        {
            spf!(errmes; "Column #", i, ": ");
            wrtferr(out, &errmes, &mut status, 2);
        }
        usrdata.indatatyp[i] = datatype;
        usrdata.mask[i] = 255;
        if datatype == TBIT {
            repeat %= 8;
            usrdata.mask[i] >>= repeat;
            if repeat == 0 {
                usrdata.mask[i] = 0;
            }
        }
    }

    if niter > 0 {
        iterdata_reset();
        let up: *mut c_void = core::ptr::from_mut::<UserIter>(&mut usrdata).cast();
        if fits_iterate_data(
            niter,
            &mut iter_col,
            offset,
            rows_per_loop,
            iterdata,
            up,
            &mut status,
        ) != 0
        {
            wrtserr_str(out, "When Reading data, ", &mut status, 2);
        }
    }

    /* the C free()s iter_col, numlist, floatlist, cmplist, txtlist and the
    usrdata arrays here; RAII does that. */

    if ndesc != 0 {
        /* ------------read the variable length vectors -------------------*/
        usrdata.datamax = vec![0.0; ndesc as usize];
        usrdata.datamin = vec![0.0; ndesc as usize];
        usrdata.tnull = vec![0.0; ndesc as usize];
        maxlen = vec![0; ndesc as usize];
        dflag = vec![0; ndesc as usize];
        perbyte = vec![0; ndesc as usize];
        isVarQFormat = vec![0; ndesc as usize];
        fits_get_num_rows(infits, &mut totalrows, &mut status);
        status = 0;

        /* this routine now only reads and test BIT, LOGICAL, and STRING columns */
        /* There is no point in reading the other columns because the other datatypes */
        /* have no possible invalid values.  */

        for i in 0..ndesc as usize {
            icol = desclist[i];
            parse_vtform(
                infits,
                out,
                hduptr,
                icol,
                &mut datatype,
                &mut maxlen[i],
                &mut isVarQFormat[i],
            );
            dflag[i] = 4;
            match datatype {
                d if d == -TBIT => {
                    dflag[i] = 1;
                    perbyte[i] = -8;
                }
                d if d == -TBYTE => {
                    perbyte[i] = 1;
                }
                d if d == -TLOGICAL => {
                    dflag[i] = 3;
                    perbyte[i] = 1;
                }
                d if d == -TSTRING => {
                    dflag[i] = 0;
                    perbyte[i] = 1;
                }
                d if d == -TSHORT => {
                    perbyte[i] = 2;
                }
                d if d == -TLONG => {
                    perbyte[i] = 4;
                }
                d if d == -TFLOAT => {
                    perbyte[i] = 4;
                }
                d if d == -TDOUBLE => {
                    perbyte[i] = 8;
                }
                d if d == -TCOMPLEX => {
                    dflag[i] = 2;
                    perbyte[i] = 8;
                }
                d if d == -TDBLCOMPLEX => {
                    dflag[i] = 2;
                    perbyte[i] = 16;
                }
                _ => {}
            }
        }

        maxmax = maxlen[0];
        for i in 1..ndesc as usize {
            if maxmax < maxlen[i] {
                maxmax = maxlen[i];
            }
        }
        if maxmax < 0 {
            maxmax = 100;
        }
        cdata = vec![0; (maxmax + 1) as usize];

        jl = 1;
        while jl <= totalrows {
            for i in 0..ndesc as usize {
                icol = desclist[i];

                /* read and check the descriptor length and offset values */
                if fits_read_descriptll(
                    infits,
                    icol,
                    jl as LONGLONG,
                    Some(&mut length),
                    Some(&mut toffset),
                    &mut status,
                ) != 0
                {
                    spf!(errtmp; "Row #", jl as i64, " Col.#", icol, ": ");
                    wrtferr(out, &errtmp, &mut status, 2);
                }
                if isVarQFormat[i] == 0 {
                    if largeVarLengthWarned == 0 && length > 2147483647 {
                        spf!(errmes;
                            "Var row length exceeds maximum 32-bit signed int.  ");
                        spf!(errtmp;
                            "First detected for Row #", jl as i64, " Column #", icol);
                        scat!(errmes; CS(&errtmp));
                        wrtwrn(out, &errmes, 0);
                        largeVarLengthWarned = 1;
                    }
                    if largeVarOffsetWarned == 0 && toffset > 2147483647 {
                        spf!(errmes;
                            "Heap offset for var length row exceeds maximum 32-bit signed int.  ");
                        spf!(errtmp;
                            "First detected for Row #", jl as i64, " Column #", icol);
                        scat!(errmes; CS(&errtmp));
                        wrtwrn(out, &errmes, 0);
                        largeVarOffsetWarned = 1;
                    }
                }

                if length > maxlen[i] as LONGLONG && maxlen[i] > -1 {
                    spf!(errmes;
                        "Descriptor of Column #", icol, " at Row ", jl as i64, ": ");
                    spf!(errtmp;
                        "nelem(", length as i64, ") > maxlen(", maxlen[i] as i64,
                        ") given by TFORM", icol, ".");
                    scat!(errmes; CS(&errtmp));
                    wrterr(out, &errmes, 1);
                }

                if perbyte[i] < 0 {
                    bytelength = (length / 8) as c_long;
                } else {
                    bytelength = (length * perbyte[i] as LONGLONG) as c_long;
                }

                if toffset + bytelength as LONGLONG > hduptr.pcount {
                    spf!(errmes;
                        "Descriptor of Column #", icol, " at Row ", jl as i64, ": ");
                    spf!(errtmp;
                        " offset of first element(", toffset as i64, ") + nelem(",
                        length as i64, ")");
                    scat!(errmes; CS(&errtmp));
                    if perbyte[i] < 0 {
                        spf!(errtmp;
                            "/8 >  total heap area  = ", hduptr.pcount as i64, ".");
                    } else {
                        spf!(errtmp;
                            "*", perbyte[i], " >  total heap area  = ", hduptr.pcount as i64,
                            ".");
                    }
                    scat!(errmes; CS(&errtmp));
                    wrterr(out, &errmes, 2);
                }

                if length == 0 {
                    continue;
                } /* skip the 0 length array */

                /* now check the values in BIT, LOGICAL, and String columns */
                rlength = length as c_long;
                if length > maxmax as LONGLONG {
                    rlength = maxmax;
                }

                if dflag[i] == 1 { /* read BIT column */

                    /*  NOT YET IMPLEMENTED:  This code should test that the fill bits that
                    pad out the last byte are all zero. */
                } else if dflag[i] == 0 {
                    /* read String column */
                    anynul = 0;
                    let mut sbuf: Vec<c_char> = vec![0; (maxmax + 1) as usize];
                    {
                        let mut arr: [&mut [c_char]; 1] = [&mut sbuf];
                        if fits_read_col_str(
                            infits,
                            icol,
                            jl as LONGLONG,
                            1,
                            rlength as LONGLONG,
                            None,
                            &mut arr,
                            Some(&mut anynul),
                            &mut status,
                        ) != 0
                        {
                            spf!(errtmp; "Row #", jl as i64, " Col.#", icol, ": ");
                            wrtferr(out, &errtmp, &mut status, 2);
                        } else {
                            j = 0;
                            while sbuf[j] != 0 {
                                if (sbuf[j] as u8) > 126 || (sbuf[j] as u8) < 32 {
                                    spf!(errmes;
                                        "String in row #", jl as i64, ", and column #", icol,
                                        " contains non-ASCII text.");
                                    wrterr(out, &errmes, 1);
                                    spf!(errmes;
                                        "             (This error is reported only once; other rows may have errors).");
                                    print_fmt(out, &errmes, 13);
                                    break;
                                }
                                j += 1;
                            }
                        }
                    }
                    cdata = sbuf;
                } else if dflag[i] == 3 {
                    /* read Logical column */
                    anynul = 0;
                    if fits_read_col_log(
                        infits,
                        icol,
                        jl as LONGLONG,
                        1,
                        rlength as LONGLONG,
                        lnull,
                        &mut cdata,
                        Some(&mut anynul),
                        &mut status,
                    ) != 0
                    {
                        spf!(errtmp; "Row #", jl as i64, " Col.#", icol, ": ");
                        wrtferr(out, &errtmp, &mut status, 2);
                    } else {
                        k = 0;
                        while k < rlength {
                            if (cdata[k as usize] as i8) > 2 {
                                spf!(errmes;
                                    "Logical value in row #", jl as i64, ", column #", icol,
                                    " not equal to 'T', 'F', or 0");
                                wrterr(out, &errmes, 1);
                                spf!(errmes;
                                    "             (This error is reported only once; other rows may have errors).");
                                print_fmt(out, &errmes, 13);
                                break;
                            }
                            k += 1;
                        }
                    }
                }
            }
            jl += 1;
        }
    }

    /* data_end: */
    i = 0;
    while i < ncols as usize {
        hduptr.datamax[i][12] = 0;
        hduptr.datamin[i][12] = 0;
        hduptr.tnull[i][11] = 0;
        i += 1;
    }
}

/***********************************************************************/
/* iterator work function */

/* The C keeps its per-table state in function statics; the same state lives
here so that it survives across the iterator's repeated calls. */
thread_local! {
    static REPEAT_V: RefCell<Vec<c_long>> = const { RefCell::new(Vec::new()) };
    static DATATYPE_V: RefCell<Vec<c_int>> = const { RefCell::new(Vec::new()) };
}
thread_local! { static IT_NNUM: Cell<c_int> = const { Cell::new(0) }; }
thread_local! { static IT_NTXT: Cell<c_int> = const { Cell::new(0) }; }
thread_local! { static IT_NCMP: Cell<c_int> = const { Cell::new(0) }; }
thread_local! { static IT_NFLOAT: Cell<c_int> = const { Cell::new(0) }; }
thread_local! { static FIND_BADBIT: Cell<c_int> = const { Cell::new(0) }; }
thread_local! { static FIND_BADDOT: Cell<c_int> = const { Cell::new(0) }; }
thread_local! { static FIND_BADSPACE: Cell<c_int> = const { Cell::new(0) }; }
thread_local! { static FIND_BADCHAR: Cell<c_int> = const { Cell::new(0) }; }
thread_local! { static FIND_BADLOG: Cell<c_int> = const { Cell::new(0) }; }

fn iterdata_reset() {
    FIND_BADBIT.with(|c| c.set(0));
    FIND_BADDOT.with(|c| c.set(0));
    FIND_BADSPACE.with(|c| c.set(0));
    FIND_BADCHAR.with(|c| c.set(0));
    FIND_BADLOG.with(|c| c.set(0));
}

pub extern "C" fn iterdata(
    totaln: c_long,
    _offset: c_long,
    firstn: c_long,
    nrows: c_long,
    narray: c_int,
    iter_col: *mut iteratorCol,
    usrdata: *mut c_void,
) -> c_int {
    let mut i: c_int;
    let mut j: c_long;
    let mut k: c_long;
    let mut l: c_long;
    let mut nelem: c_long;
    let mut errmes: [c_char; ERRMES_LEN] = [0; ERRMES_LEN];
    let mut comm: [c_char; COMM_LEN] = [0; COMM_LEN];

    // SAFETY: `usrdata` is the &mut UserIter that test_data handed to
    // fits_iterate_data, and `iter_col` is the narray-long column array the
    // iterator owns for the duration of this call.
    let usrpt: &UserIter = unsafe { &*(usrdata as *const UserIter) };
    let cols: &[iteratorCol] = unsafe { std::slice::from_raw_parts(iter_col, narray as usize) };

    if firstn == 1 {
        /* first time for this table, so initialize */
        IT_NNUM.with(|c| c.set(usrpt.nnum));
        IT_NCMP.with(|c| c.set(usrpt.ncmp));
        IT_NTXT.with(|c| c.set(usrpt.ntxt));
        IT_NFLOAT.with(|c| c.set(usrpt.nfloat));
        REPEAT_V.with(|r| {
            let mut r = r.borrow_mut();
            r.clear();
            for c in cols.iter() {
                r.push(c.repeat); /* fits_iter_get_repeat() */
            }
        });
        DATATYPE_V.with(|d| {
            let mut d = d.borrow_mut();
            d.clear();
            for c in cols.iter() {
                d.push(c.datatype); /* fits_iter_get_datatype() */
            }
        });
        iterdata_reset();
    }

    let nnum = IT_NNUM.with(Cell::get);
    let ncmp = IT_NCMP.with(Cell::get);
    let ntxt = IT_NTXT.with(Cell::get);
    let nfloat = IT_NFLOAT.with(Cell::get);
    let repeat = REPEAT_V.with(|r| r.borrow().clone());
    let datatype = DATATYPE_V.with(|d| d.borrow().clone());

    /* columns from  1 to nnum are scalar numerical columns.
    columns from  nnum+1 to  nnum+ncmp are complex columns. (not used any more)
    columns from  nnum+ncmp are text columns */

    /* deal with the numerical column */
    i = 0;
    while i < nnum + ncmp {
        let iu = i as usize;
        // SAFETY: for a numeric column the iterator's array is repeat*nrows+1 doubles.
        let data: *const f64 = cols[iu].array.cast::<f64>();
        nelem = nrows * repeat[iu];
        if i >= nnum {
            nelem = 2 * nrows * repeat[iu];
        }
        if nelem == 0 {
            i += 1;
            continue;
        }
        FIND_BADBIT.with(|c| c.set(0));

        /* check for the bit jurisfication  */
        if FIND_BADBIT.with(Cell::get) == 0 && usrpt.indatatyp[iu] == TBIT {
            k = 0;
            while k < nrows {
                j = (k + 1) * repeat[iu];
                let bdata = unsafe { *data.offset(j as isize) } as c_uchar;
                if bdata & usrpt.mask[iu] != 0 {
                    spf!(errmes;
                        "Row #", (firstn + k) as i64, ", and Column #", cols[iu].colnum,
                        ": X vector ");
                    l = 1;
                    while l <= repeat[iu] {
                        let v = unsafe { *data.offset((k * repeat[iu] + l) as isize) } as c_uchar;
                        spf!(comm; "0x", format!("{v:02x}"), " ");
                        scat!(errmes; CS(&comm));
                        l += 1;
                    }
                    scat!(errmes; "is not left justified.");
                    wrterr(usrpt.out, &errmes, 2);
                    spf!(errmes; "             (Other rows may have errors).");
                    print_fmt(usrpt.out, &errmes, 13);
                    FIND_BADBIT.with(|c| c.set(1));
                    break;
                }
                k += 1;
            }
        }
        i += 1;
    }

    /* deal with character and logical columns */
    i = nnum + ncmp;
    while i < nnum + ncmp + ntxt {
        let iu = i as usize;
        if datatype[iu] == TSTRING {
            /* character */
            nelem = nrows;
            if nelem == 0 {
                i += 1;
                continue;
            }
            // SAFETY: a TSTRING iterator column's array is a char*[nrows+1];
            // slot 0 is the null-value string, slots 1..=nrows the row values,
            // all NUL-terminated and owned by the iterator for this call.
            let cdata: *const *const c_char = cols[iu].array.cast::<*const c_char>();
            let row = |n: c_long| -> &[u8] {
                unsafe { CStr::from_ptr(*cdata.offset(n as isize)) }.to_bytes()
            };
            FIND_BADCHAR.with(|c| c.set(0));

            /* test for illegal ASCII text characters > 126  or < 32 */
            if FIND_BADCHAR.with(Cell::get) == 0 {
                k = 0;
                while k < nrows {
                    /* NB: the C breaks out of the inner scan only, so the row
                    loop carries on and every offending row is reported. */
                    if row(k + 1).iter().any(|&c| !(32..=126).contains(&c)) {
                        spf!(errmes;
                            "String in row #", (firstn + k) as i64, ", column #",
                            cols[iu].colnum, " contains non-ASCII text.");
                        wrterr(usrpt.out, &errmes, 1);
                        spf!(errmes; "             (Other rows may have errors).");
                        print_fmt(usrpt.out, &errmes, 13);
                        FIND_BADCHAR.with(|c| c.set(1));
                    }
                    k += 1;
                }
            }
        } else {
            /* logical value */
            nelem = nrows * repeat[iu];
            if nelem == 0 {
                i += 1;
                continue;
            }
            // SAFETY: for a logical column the array is repeat*nrows+1 bytes.
            let ldata: *const c_uchar = cols[iu].array.cast::<c_uchar>();

            /* test for illegal logical column values */
            /* The first element in the array gives the value that is used to
            represent nulls */
            if FIND_BADLOG.with(Cell::get) == 0 {
                j = 1;
                while j <= nrows * repeat[iu] {
                    if unsafe { *ldata.offset(j as isize) } > 2 {
                        spf!(errmes;
                            "Logical value in row #",
                            ((firstn + j - 2) / repeat[iu] + 1) as i64, ", column #",
                            cols[iu].colnum, " not equal to 'T', 'F', or 0");
                        wrterr(usrpt.out, &errmes, 1);
                        spf!(errmes; "             (Other rows may have similar errors).");
                        print_fmt(usrpt.out, &errmes, 13);
                        FIND_BADLOG.with(|c| c.set(1));
                        break;
                    }
                    j += 1;
                }
            }
        }
        i += 1;
    }

    i = nnum + ncmp + ntxt;
    while i < nnum + ncmp + ntxt + nfloat {
        let iu = i as usize;
        nelem = nrows;
        if nelem == 0 {
            i += 1;
            continue;
        }
        // SAFETY: TSTRING iterator column, char*[nrows+1] as above.
        // SAFETY: as above -- a TSTRING iterator column.
        let cdata: *const *const c_char = cols[iu].array.cast::<*const c_char>();
        let row = |n: c_long| -> &[u8] {
            unsafe { CStr::from_ptr(*cdata.offset(n as isize)) }.to_bytes()
        };
        FIND_BADDOT.with(|c| c.set(0));
        FIND_BADSPACE.with(|c| c.set(0));

        /* test for missing (implicit) decimal point in floating point numbers */
        if FIND_BADDOT.with(Cell::get) == 0 {
            k = 0;
            while k < nrows {
                let floatvalue = row(k + 1);
                if floatvalue != row(0) && !floatvalue.contains(&b'.') {
                    let floatvalue = skip_spaces(floatvalue); /* skip leading spaces */

                    if !floatvalue.is_empty() {
                        /* ignore completely blank fields */

                        spf!(errmes;
                            "Number in row #", (firstn + k) as i64, ", column #",
                            cols[iu].colnum, " has no decimal point:");
                        wrterr(usrpt.out, &errmes, 1);
                        spf!(errmes; BS(floatvalue));
                        scat!(errmes; "  (Other rows may have similar errors).");
                        print_fmt(usrpt.out, &errmes, 13);
                        FIND_BADDOT.with(|c| c.set(1));
                        break;
                    }
                }
                k += 1;
            }
        }

        if FIND_BADSPACE.with(Cell::get) == 0 {
            k = 0;
            while k < nrows {
                let floatvalue = row(k + 1);

                if floatvalue != row(0) {
                    /* not a null value? */
                    /* The C strips the trailing blanks in place; nothing reads
                    the row again afterwards, so a trimmed view is equivalent. */
                    let floatvalue = trim_end_spaces(skip_spaces(floatvalue));

                    if floatvalue.contains(&b' ') {
                        spf!(errmes;
                            "Number in row #", (firstn + k) as i64, ", column #",
                            cols[iu].colnum, " has embedded space:");
                        wrterr(usrpt.out, &errmes, 1);
                        spf!(errmes; BS(floatvalue));
                        scat!(errmes; "  (Other rows may have similar errors).");
                        print_fmt(usrpt.out, &errmes, 13);
                        FIND_BADSPACE.with(|c| c.set(1));
                        break;
                    }
                }
                k += 1;
            }
        }
        i += 1;
    }

    if firstn + nrows - 1 == totaln {
        /* the C free()s flag_minmax, datatype and repeat here */
        REPEAT_V.with(|r| r.borrow_mut().clear());
        DATATYPE_V.with(|d| d.borrow_mut().clear());
    }
    0
}

/*************************************************************
*
*      test_agap
*
*   Test the bytes between the ASCII table column.
*
*************************************************************/
fn test_agap(
    infits: &mut fitsfile, /* input fits file   */
    out: Out,              /* output ascii file */
    hduptr: &mut FitsHdu,  /* fits hdu pointer  */
) {
    let ncols: c_int;
    let mut nrows: LONGLONG = 0;
    let mut irows: c_long = 0;
    let rowlen: LONGLONG;
    let mut data: Vec<c_uchar>;
    let mut temp: Vec<c_int>;
    let mut p: usize;
    let mut i: LONGLONG;
    let mut j: LONGLONG;
    let mut ntodo: c_long;
    let mut nerr: c_long = 0;
    let mut status: c_int = 0;
    let mut keyname: [c_char; 9] = [0; 9];
    let mut tform: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut comment: [c_char; 256] = [0; 256];
    let mut typecode: c_int = 0;
    let mut decimals: c_int = 0;
    let mut width: c_long = 0;
    let mut tbcol: c_long = 0;
    let mut firstrow: c_long = 1;
    let mut errmes: [c_char; ERRMES_LEN] = [0; ERRMES_LEN];
    nerr = 0;

    if hduptr.hdutype != ASCII_TBL {
        return;
    }
    ncols = hduptr.ncols;
    fits_get_num_rowsll(infits, &mut nrows, &mut status);
    status = 0;

    fits_get_rowsize(infits, &mut irows, &mut status);
    status = 0;
    rowlen = hduptr.naxes[0];
    data = vec![0; (rowlen * irows as LONGLONG) as usize];

    /* Create a template row with data fields filled with 1s.
    Used below - different ASCII rules apply within data columns
    vs. between data columns. */

    temp = vec![0; rowlen as usize];
    for k in 1..=ncols {
        spf!(keyname; "TFORM", k);
        let mut cbuf: [c_char; rsfitsio::fitsio::FLEN_COMMENT] =
            [0; rsfitsio::fitsio::FLEN_COMMENT];
        fits_read_key_str(infits, &keyname, &mut tform, Some(&mut cbuf), &mut status);
        let _ = &mut comment;
        if ffasfm_safe(
            &tform,
            Some(&mut typecode),
            Some(&mut width),
            Some(&mut decimals),
            &mut status,
        ) != 0
        {
            wrtferr_str(out, "", &mut status, 1);
        }
        spf!(keyname; "TBCOL", k);
        fits_read_key_lng(infits, &keyname, &mut tbcol, Some(&mut cbuf), &mut status);
        let mut t = tbcol;
        while t < tbcol + width {
            if t >= 1 && (t as LONGLONG) <= rowlen {
                temp[(t - 1) as usize] = 1;
            }
            t += 1;
        }
    }

    i = nrows;
    while i > 0 {
        if i > irows as LONGLONG {
            ntodo = irows;
        } else {
            ntodo = i as c_long;
        }

        p = 0;
        if ffgtbb_safe(
            infits,
            firstrow as LONGLONG,
            1,
            rowlen * ntodo as LONGLONG,
            &mut data,
            &mut status,
        ) != 0
        {
            wrtferr_str(out, "", &mut status, 1);
        }
        j = 0;
        while j < rowlen * ntodo as LONGLONG {
            let c = data[p];
            if !isascii_c(c as c_char) {
                if nerr == 0 {
                    spf!(errmes; "row ", (j / rowlen + 1) as i64,
                        " contains non-ASCII characters.");
                    wrterr(out, &errmes, 1);
                }
                nerr += 1;
            } else if isascii_c(c as c_char)
                && !isprint_c(c as c_char)
                && temp[(j % rowlen) as usize] != 0
            {
                if nerr == 0 {
                    spf!(errmes; "row ", (j / rowlen + 1) as i64,
                        " data contains non-ASCII-text characters.");
                    wrterr(out, &errmes, 1);
                }
                nerr += 1;
            }
            p += 1;
            j += 1;
        }
        firstrow += ntodo;
        i -= ntodo as LONGLONG;
    }
    if nerr != 0 {
        spf!(errmes;
            "This ASCII table contains ", nerr as i64, " non-ASCII-text characters");
        wrterr(out, &errmes, 1);
    }
}

/*************************************************************
*
*      test_checksum
*
*   Test the checksum of the hdu
*
*************************************************************/
fn test_checksum(
    infits: &mut fitsfile, /* input fits file   */
    out: Out,              /* output ascii file */
) {
    let mut status: c_int = 0;
    let mut dataok: c_int = 0;
    let mut hduok: c_int = 0;

    if fits_verify_chksum(infits, &mut dataok, &mut hduok, &mut status) != 0 {
        wrtferr_str(out, "verifying checksums: ", &mut status, 2);
        return;
    }

    if dataok == -1 {
        wrtwrn_str(
            out,
            "Data checksum is not consistent with  the DATASUM keyword",
            0,
        );
    }

    if hduok == -1 {
        if dataok == 1 {
            wrtwrn_str(
                out,
                "Invalid CHECKSUM means header has been modified. (DATASUM is OK) ",
                0,
            );
        } else {
            wrtwrn_str(out, "HDU checksum is not in agreement with CHECKSUM.", 0);
        }
    }
}

#[allow(dead_code)]
fn _unused() {
    let _ = fits_read_tblbytes;
    let _ = fits_iterate_data;
}
