/************************************************************************/
/*                                                                      */
/*                       CFITSIO Lexical Parser                         */
/*                                                                      */
/* This file is one of 3 files containing code which parses an          */
/* arithmetic expression and evaluates it in the context of an input    */
/* FITS file table extension.  The CFITSIO lexical parser is divided    */
/* into the following 3 parts/files: the CFITSIO "front-end",           */
/* eval_f.c, contains the interface between the user/CFITSIO and the    */
/* real core of the parser; the FLEX interpreter, eval_l.c, takes the   */
/* input string and parses it into tokens and identifies the FITS       */
/* information required to evaluate the expression (ie, keywords and    */
/* columns); and, the BISON grammar and evaluation routines, eval_y.c,  */
/* receives the FLEX output and determines and performs the actual      */
/* operations.  The files eval_l.c and eval_y.c are produced from       */
/* running flex and bison on the files eval.l and eval.y, respectively. */
/* (flex and bison are available from any GNU archive: see www.gnu.org) */
/*                                                                      */
/* The grammar rules, rather than evaluating the expression in situ,    */
/* builds a tree, or Nodal, structure mapping out the order of          */
/* operations and expression dependencies.  This "compilation" process  */
/* allows for much faster processing of multiple rows.  This technique  */
/* was developed by Uwe Lammers of the XMM Science Analysis System,     */
/* although the CFITSIO implementation is entirely code original.       */
/*                                                                      */
/*                                                                      */
/* Modification History:                                                */
/*                                                                      */
/*   Kent Blackburn      c1992  Original parser code developed for the  */
/*                              FTOOLS software package, in particular, */
/*                              the fselect task.                       */
/*   Kent Blackburn      c1995  BIT column support added                */
/*   Peter D Wilson   Feb 1998  Vector column support added             */
/*   Peter D Wilson   May 1998  Ported to CFITSIO library.  User        */
/*                              interface routines written, in essence  */
/*                              making fselect, fcalc, and maketime     */
/*                              capabilities available to all tools     */
/*                              via single function calls.              */
/*   Peter D Wilson   Jun 1998  Major rewrite of parser core, so as to  */
/*                              create a run-time evaluation tree,      */
/*                              inspired by the work of Uwe Lammers,    */
/*                              resulting in a speed increase of        */
/*                              10-100 times.                           */
/*   Peter D Wilson   Jul 1998  gtifilter(a,b,c,d) function added       */
/*   Peter D Wilson   Aug 1998  regfilter(a,b,c,d) function added       */
/*   Peter D Wilson   Jul 1999  Make parser fitsfile-independent,       */
/*                              allowing a purely vector-based usage    */
/*   Peter D Wilson   Aug 1999  Add row-offset capability               */
/*   Peter D Wilson   Sep 1999  Add row-range capability to ffcalc_rng  */
/*                                                                      */
/************************************************************************/

use libc::{c_float, memcpy, memset};
use std::{cmp, ptr};

use crate::buffers::{ffgbyt, ffgtbb_safe, ffmbyt_safe, ffpbyt, ffptbb_safe};
use crate::c_types::{c_char, c_double, c_int, c_long, c_short, c_uchar, c_void};

use crate::aliases::rust_api::{fits_binary_tform, fits_set_atblnull, fits_set_btblnull};
use crate::aliases::rust_api::{
    fits_create_img, fits_get_colnum, fits_get_coltype, fits_get_hdrspace, fits_get_hdu_type,
    fits_get_img_param, fits_get_keyclass, fits_get_keytype, fits_read_imgnull, fits_read_key_dbl,
    fits_read_key_lng, fits_read_key_log, fits_read_key_str, fits_read_keyword, fits_read_record,
    fits_read_tdim, fits_set_imgnull, fits_update_key_lng, fits_write_record,
};
use crate::cfileio::ffimport_file_safer;
use crate::editcol::{ffdrow_safe, fficol_safe, ffirow_safe};
use crate::eval_defs::{
    CONST_OP, DataInfo, MAX_STRLEN, MAXDIMS, MAXVARNAME, Node, P_ERROR, ParseData,
    ParseStatusVariables, parseInfo,
};
use crate::eval_l::{
    fits_parser_yylex_destroy, fits_parser_yylex_init_extra, fits_parser_yyrestart, yyguts_t,
};
use crate::eval_tab::{FITS_PARSER_YYSTYPE, fits_parser_yytokentype};
use crate::eval_y::{Evaluate_Parser, fits_parser_yyparse, gtifilt_fct, regfilt_fct};
use crate::fitscore::{
    ffcmph_safer, ffcmsg_safe, ffgcno_safe, ffgdesll_safe, ffgncl_safe, ffgnrw_safe, ffiblk,
    ffkeyn_safe, ffmahd_safe, ffpdes_safe, ffpmrk_safe, fits_strcasecmp,
};
use crate::fitscore::{ffpmsg_slice, ffpmsg_str, ffrdef_safe};
use crate::fitsio2::{FSTRCMP, IGNORE_EOF, REPORT_EOF};
use crate::getcolb::{fffi2i1, fffr4i1, fffr8i1, ffgcvb_safe};
use crate::getcold::{ffgcfd_safe, ffgcvd_safe};
use crate::getcole::fffr8r4;
use crate::getcoli::{fffr4i2, fffr8i2};
use crate::getcolj::{fffr4i4, fffr4i8, fffr8i4, fffr8i8, ffgcfj_safe, ffgcvj_safe};
use crate::getcolk::{fffr4int, fffr8int};
use crate::getcoll::ffgcfl_safe;
use crate::getcols::{ffgcfs_safe, ffgcvs_safe};
use crate::getkey::ffgkyj_safe;
use crate::getkey::{ffgcrd_safe, ffgknjj_safe};
use crate::modkey::{ffdkey_safe, ffukyd_safe, ffukyj_safe, ffukyl_safe, ffukys_safe};
use crate::putcol::{ffiter_safe, fits_iter_set_by_num_safe};
use crate::putkey::{ffpcom_safe, ffphis_safe, ffpkyj_safe, ffpkys_safe, ffptdm_safe};
use crate::region::{SAORegion, fits_free_region};
use crate::wrappers::{strcat, strcpy, strlen, strlen_safe, strncpy_safe};
use crate::{BL, NullCheckType, cs, fitsio::*, int_snprintf, raw_to_slice};
use bytemuck::{cast_slice, cast_slice_mut};
use core::ffi::CStr;
use libc::{free, malloc};

// Constants needed for the function
const UCHAR_MAX: LONGLONG = 255;
const SHRT_MIN: LONGLONG = c_short::MIN as LONGLONG;
const SHRT_MAX: LONGLONG = c_short::MAX as LONGLONG;
const INT_MIN: LONGLONG = c_int::MIN as LONGLONG;
const LONG_MIN: LONGLONG = c_long::MIN as LONGLONG;
use std::mem::size_of;

pub(crate) struct ffffrw_workdata {
    prownum: *mut c_long,
    lParse: *mut ParseData,
}

static DEBUG_PIXFILTER: c_int = 0;

// Free macro for raw pointers
macro_rules! FREE {
    ($x:expr) => {
        if !$x.is_null() {
            unsafe {
                libc::free($x as *mut libc::c_void);
            }
            $x = std::ptr::null_mut();
        } else {
            println!("invalid free at {}:{}", file!(), line!());
        }
    };
}

/*---------------------------------------------------------------------------*/
/// Evaluate a boolean expression using the indicated rows, returning an
/// array of flags indicating which rows evaluated to TRUE/FALSE
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn fffrow(
    fptr: *mut fitsfile,      /* I - Input FITS file                   */
    expr: *mut c_char,        /* I - Boolean expression                */
    firstrow: c_long,         /* I - First row of table to eval        */
    nrows: c_long,            /* I - Number of rows to evaluate        */
    n_good_rows: *mut c_long, /* O - Number of rows eval to True       */
    row_status: *mut c_char,  /* O - Array of boolean results          */
    status: *mut c_int,       /* O - Error status                      */
) -> c_int {
    unsafe {
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);
        let n_good_rows = n_good_rows.as_mut().expect(NULL_MSG);

        let row_status = std::slice::from_raw_parts_mut(row_status, nrows as usize);

        raw_to_slice!(expr);

        fffrow_safe(fptr, expr, firstrow, nrows, n_good_rows, row_status, status)
    }
}

/*---------------------------------------------------------------------------*/
/// Evaluate a boolean expression using the indicated rows, returning an
/// array of flags indicating which rows evaluated to TRUE/FALSE
pub fn fffrow_safe(
    fptr: &mut fitsfile,       /* I - Input FITS file                   */
    expr: &[c_char],           /* I - Boolean expression                */
    firstrow: c_long,          /* I - First row of table to eval        */
    nrows: c_long,             /* I - Number of rows to evaluate        */
    n_good_rows: &mut c_long,  /* O - Number of rows eval to True       */
    row_status: &mut [c_char], /* O - Array of boolean results          */
    status: &mut c_int,        /* O - Error status                      */
) -> c_int {
    let mut naxis: c_int = 0;
    let constant: bool;
    let mut nelem: c_long = 0;
    let mut naxes: [c_long; MAXDIMS as usize] = [0; MAXDIMS as usize];
    let mut elem: c_long;
    let result: c_char;
    let mut lParse: ParseData = ParseData::default();

    if *status != 0 {
        return *status;
    }

    let mut Info: parseInfo = Default::default();

    if ffiprs(
        fptr,
        0,
        expr,
        MAXDIMS,
        &mut Info.datatype,
        &mut nelem,
        &mut naxis,
        &mut naxes,
        &mut lParse,
        status,
    ) != 0
    {
        ffcprs(&mut lParse);
        return *status;
    }

    if nelem < 0 {
        constant = true;
        nelem = -nelem;
    } else {
        constant = false;
    }

    if Info.datatype != TLOGICAL || nelem != 1 {
        ffcprs(&mut lParse);
        ffpmsg_str("Expression does not evaluate to a logical scalar.");
        *status = PARSE_BAD_TYPE;
        return *status;
    }

    if constant {
        /* No need to call parser... have result from ffiprs */
        result = unsafe { (lParse.Nodes[lParse.resultNode as usize]).value.data.log };
        *n_good_rows = nrows;

        for elem in 0..nrows {
            row_status[elem as usize] = result;
        }
    } else {
        let firstrow = if firstrow > 1 { firstrow } else { 1 };
        Info.dataPtr = row_status.as_mut_ptr() as *mut c_void;
        Info.nullPtr = ptr::null_mut();
        Info.maxRows = nrows;
        Info.parseData = &mut lParse;

        let colData_slice = &mut lParse.colData[..];
        if ffiter_safe(
            lParse.nCols,
            colData_slice,
            firstrow - 1,
            0,
            fits_parser_workfn,
            &Info as *const _ as *mut c_void,
            status,
        ) == -1
        {
            *status = 0; /* -1 indicates exited without error before end... OK */
        }

        if *status != 0 {

            /***********************/
            /* Error... Do nothing */
            /***********************/
        } else {
            /***********************************/
            /* Count number of good rows found */
            /***********************************/

            *n_good_rows = 0;

            for elem in 0..Info.maxRows {
                if row_status[elem as usize] == 1 {
                    *n_good_rows += 1;
                }
            }
        }
    }

    ffcprs(&mut lParse);
    *status
}

/*--------------------------------------------------------------------------*/
/// Evaluate an expression on all rows of a table.  If the input and output
/// files are not the same, copy the TRUE rows to the output file.  If the
/// files are the same, delete the FALSE rows (preserve the TRUE rows).
/// Can copy rows between extensions of the same file, *BUT* if output
/// extension is before the input extension, the second extension *MUST* be
/// opened using ffreopen, so that CFITSIO can handle changing file lengths
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffsrow(
    infptr: *mut fitsfile,  /* I - Input FITS file                      */
    outfptr: *mut fitsfile, /* I - Output FITS file                     */
    expr: *const c_char,    /* I - Boolean expression                   */
    status: *mut c_int,     /* O - Error status                         */
) -> c_int {
    unsafe {
        let infptr = infptr.as_mut().expect(NULL_MSG);
        let outfptr = outfptr.as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);

        raw_to_slice!(expr);

        ffsrow_safe(infptr, outfptr, expr, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Evaluate an expression on all rows of a table.  If the input and output
/// files are not the same, copy the TRUE rows to the output file.  If the
/// files are the same, delete the FALSE rows (preserve the TRUE rows).
/// Can copy rows between extensions of the same file, *BUT* if output
/// extension is before the input extension, the second extension *MUST* be
/// opened using ffreopen, so that CFITSIO can handle changing file lengths
pub fn ffsrow_safe(
    infptr: &mut fitsfile,  /* I - Input FITS file                      */
    outfptr: &mut fitsfile, /* I - Output FITS file                     */
    expr: &[c_char],        /* I - Boolean expression                   */
    status: &mut c_int,     /* O - Error status                         */
) -> c_int {
    let mut Info: parseInfo = Default::default();
    let mut naxis: c_int = 0;
    let mut constant: c_int = 0;
    let mut nelem: c_long = 0;
    let mut rdlen: c_long = 0;
    let mut naxes: [c_long; MAXDIMS as usize] = [0; MAXDIMS as usize];
    let mut maxrows: c_long = 0;
    let mut nbuff: c_long = 0;
    let mut nGood: c_long = 0;
    let mut inloc: c_long = 0;
    let mut outloc: c_long = 0;
    let mut ntodo: LONGLONG = 0;
    let mut inbyteloc: LONGLONG = 0;
    let mut outbyteloc: LONGLONG = 0;
    let mut hsize: LONGLONG = 0;
    let mut freespace: c_long = 0;
    let mut buffer: Vec<u8> = Vec::new();
    let mut result: c_char = 0;

    #[derive(Default)]
    struct in_out_ext {
        rowLength: LONGLONG,
        numRows: LONGLONG,
        heapSize: LONGLONG,
        dataStart: LONGLONG,
        heapStart: LONGLONG,
    }

    let mut inExt: in_out_ext = in_out_ext::default();
    let mut outExt: in_out_ext = in_out_ext::default();

    let mut lParse: ParseData = ParseData::default();

    unsafe {
        if *status != 0 {
            return *status;
        }

        if ffiprs(
            infptr,
            0,
            expr,
            MAXDIMS,
            &mut Info.datatype,
            &mut nelem,
            &mut naxis,
            &mut naxes,
            &mut lParse,
            status,
        ) != 0
        {
            ffcprs(&mut lParse);
            return *status;
        }

        if nelem < 0 {
            constant = 1;
            nelem = -nelem;
        } else {
            constant = 0;
        }

        /**********************************************************************/
        /* Make sure expression evaluates to the right type... logical scalar */
        /**********************************************************************/

        if Info.datatype != TLOGICAL || nelem != 1 {
            ffcprs(&mut lParse);
            ffpmsg_str("Expression does not evaluate to a logical scalar.");
            *status = PARSE_BAD_TYPE;
            return *status;
        }

        /***********************************************************/
        /*  Extract various table information from each extension  */
        /***********************************************************/

        if infptr.HDUposition != (infptr.Fptr).curhdu {
            ffmahd_safe(infptr, (infptr.HDUposition) + 1, None, status);
        }
        if *status != 0 {
            ffcprs(&mut lParse);
            return *status;
        }
        inExt.rowLength = (infptr.Fptr).rowlength as c_long;
        inExt.numRows = (infptr.Fptr).numrows;
        inExt.heapSize = (infptr.Fptr).heapsize;
        if inExt.numRows == 0 {
            /* Nothing to copy */
            ffcprs(&mut lParse);
            return *status;
        }

        if outfptr.HDUposition != (outfptr.Fptr).curhdu {
            ffmahd_safe(outfptr, (outfptr.HDUposition) + 1, None, status);
        }
        if (outfptr.Fptr).datastart < 0 {
            ffrdef_safe(outfptr, status);
        }
        if *status != 0 {
            ffcprs(&mut lParse);
            return *status;
        }
        outExt.rowLength = (outfptr.Fptr).rowlength as c_long;
        outExt.numRows = (outfptr.Fptr).numrows;
        if outExt.numRows == 0 {
            (outfptr.Fptr).heapsize = 0;
        }
        outExt.heapSize = (outfptr.Fptr).heapsize;

        if inExt.rowLength != outExt.rowLength {
            ffpmsg_str("Output table has different row length from input");
            ffcprs(&mut lParse);
            *status = PARSE_BAD_OUTPUT;
            return *status;
        }

        /***********************************/
        /*  Fill out Info data for parser  */
        /***********************************/

        Info.dataPtr = malloc((inExt.numRows + 1) as usize * size_of::<c_char>()) as *mut c_void;
        Info.nullPtr = ptr::null_mut();
        Info.maxRows = inExt.numRows as c_long;
        Info.parseData = &mut lParse as *mut ParseData;
        if Info.dataPtr.is_null() {
            ffpmsg_str("Unable to allocate memory for row selection");
            ffcprs(&mut lParse);
            *status = MEMORY_ALLOCATION;
            return *status;
        }

        /* make sure array is zero terminated */
        *(Info.dataPtr as *mut c_char).add(inExt.numRows as usize) = 0;

        if constant != 0 {
            /*  Set all rows to the same value from constant result  */

            result = (lParse.Nodes[lParse.resultNode as usize]).value.data.log;

            for ntodo in 0..inExt.numRows {
                *(Info.dataPtr as *mut c_char).add(ntodo as usize) = result;
            }
            nGood = if result != 0 { inExt.numRows } else { 0 } as c_long;
        } else {
            let col_slice = &mut lParse.colData[..];
            ffiter_safe(
                lParse.nCols,
                col_slice,
                0,
                0,
                fits_parser_workfn,
                &mut Info as *mut parseInfo as *mut c_void,
                status,
            );

            nGood = 0;

            for ntodo in 0..inExt.numRows {
                if (*(Info.dataPtr as *mut c_char).add(ntodo as usize)) != 0 {
                    nGood += 1;
                }
            }
        }

        if *status != 0 {
            /* Error... Do nothing */
        } else {
            rdlen = inExt.rowLength as c_long;

            if buffer
                .try_reserve_exact(cmp::max(500000, rdlen) as usize)
                .is_err()
            {
                FREE!(Info.dataPtr);
                ffcprs(&mut lParse);
                *status = MEMORY_ALLOCATION;
                return *status;
            } else {
                buffer.resize(cmp::max(500000, rdlen) as usize, 0);
            }

            maxrows = cmp::max(500000 / rdlen, 1);
            nbuff = 0;
            inloc = 1;
            if std::ptr::eq(infptr, outfptr) {
                /* Skip initial good rows if input==output file */
                while *(Info.dataPtr as *mut c_char).add((inloc - 1) as usize) != 0 {
                    inloc += 1;
                }
                outloc = inloc;
            } else {
                outloc = (outExt.numRows + 1) as c_long;
                if outloc > 1 {
                    ffirow_safe(outfptr, outExt.numRows, nGood, status);
                }
            }

            loop {
                if *(Info.dataPtr as *mut c_char).add((inloc - 1) as usize) != 0 {
                    let buffer_part = &mut buffer[((rdlen * nbuff) as usize)..];

                    ffgtbb_safe(infptr, inloc, 1, rdlen, buffer_part, status);
                    nbuff += 1;
                    if nbuff == maxrows {
                        ffptbb_safe(outfptr, outloc, 1, rdlen * nbuff, &buffer, status);
                        outloc += nbuff;
                        nbuff = 0;
                    }
                }
                inloc += 1;
                if *status != 0 || inloc > inExt.numRows {
                    break;
                }
            }

            if nbuff != 0 {
                let buffer_part = &buffer[((rdlen * nbuff) as usize)..];

                ffptbb_safe(outfptr, outloc, 1, rdlen * nbuff, buffer_part, status);
                outloc += nbuff;
            }

            if std::ptr::eq(infptr, outfptr) {
                if outloc <= inExt.numRows {
                    ffdrow_safe(infptr, outloc, inExt.numRows - outloc + 1, status);
                }
            } else if inExt.heapSize != 0 && nGood != 0 {
                /* Copy heap, if it exists and at least one row copied */

                /********************************************************/
                /*  Get location information from the output extension  */
                /********************************************************/

                if outfptr.HDUposition != (outfptr.Fptr).curhdu {
                    ffmahd_safe(outfptr, (outfptr.HDUposition) + 1, None, status);
                }
                outExt.dataStart = (outfptr.Fptr).datastart;
                outExt.heapStart = (outfptr.Fptr).heapstart;

                /*************************************************/
                /*  Insert more space into outfptr if necessary  */
                /*************************************************/

                hsize = outExt.heapStart + outExt.heapSize;
                freespace = ((((hsize + (BL!() - 1)) / BL!()) * BL!()) - hsize) as c_long;
                ntodo = inExt.heapSize;

                if (freespace - ntodo) < 0 {
                    /* not enough existing space? */
                    ntodo = (ntodo - freespace + (BL!() - 1)) / BL!(); /* number of blocks  */
                    ffiblk(outfptr, ntodo as c_long, 1, status); /* insert the blocks */
                }
                ffukyj_safe(
                    outfptr,
                    cs!(c"PCOUNT"),
                    inExt.heapSize + outExt.heapSize,
                    None,
                    status,
                );

                /*******************************************************/
                /*  Get location information from the input extension  */
                /*******************************************************/

                if infptr.HDUposition != (infptr.Fptr).curhdu {
                    ffmahd_safe(infptr, (infptr.HDUposition) + 1, None, status);
                }
                inExt.dataStart = (infptr.Fptr).datastart;
                inExt.heapStart = (infptr.Fptr).heapstart;

                /**********************************/
                /*  Finally copy heap to outfptr  */
                /**********************************/

                ntodo = inExt.heapSize;
                inbyteloc = inExt.heapStart + inExt.dataStart;
                outbyteloc = outExt.heapStart + outExt.dataStart + outExt.heapSize;

                while ntodo != 0 && *status == 0 {
                    rdlen = cmp::min(ntodo, 500000) as c_long;
                    ffmbyt_safe(infptr, inbyteloc, REPORT_EOF, status);
                    ffgbyt(infptr, rdlen, &mut buffer[..rdlen as usize], status);
                    ffmbyt_safe(outfptr, outbyteloc, IGNORE_EOF, status);
                    ffpbyt(outfptr, rdlen, &buffer[..rdlen as usize], status);
                    inbyteloc += rdlen;
                    outbyteloc += rdlen;
                    ntodo -= rdlen;
                }

                /***********************************************************/
                /*  But must update DES if data is being appended to a     */
                /*  pre-existing heap space.  Edit each new entry in file  */
                /***********************************************************/

                if outExt.heapSize != 0 {
                    let mut repeat: LONGLONG = 0;
                    let mut offset: LONGLONG = 0;
                    for i in 1..=(outfptr.Fptr).tfield {
                        if (*(outfptr.Fptr).tableptr.add((i - 1) as usize)).tdatatype < 0 {
                            for j in (outExt.numRows + 1)..=(outExt.numRows + nGood) {
                                ffgdesll_safe(
                                    outfptr,
                                    i,
                                    j,
                                    Some(&mut repeat),
                                    Some(&mut offset),
                                    status,
                                );
                                offset += outExt.heapSize;
                                ffpdes_safe(outfptr, i, j, repeat, offset, status);
                            }
                        }
                    }
                }
            } /*  End of HEAP copy  */
        }

        FREE!(Info.dataPtr);
        ffcprs(&mut lParse);

        ffcmph_safer(outfptr, status); /* compress heap, deleting any orphaned data */
        *status
    }
}

/*---------------------------------------------------------------------------*/
/// Calculate an expression for the indicated rows of a table, returning
/// the results, cast as datatype (TSHORT, TDOUBLE, etc), in array.  If
/// nulval==NULL, UNDEFs will be zeroed out.  For vector results, the number
/// of elements returned may be less than nelements if nelements is not an
/// even multiple of the result dimension.  Call fftexp to obtain the
/// dimensions of the results.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffcrow(
    fptr: *mut fitsfile,   /* I - Input FITS file                      */
    datatype: c_int,       /* I - Datatype to return results as        */
    expr: *const c_char,   /* I - Arithmetic expression                */
    firstrow: c_long,      /* I - First row to evaluate                */
    nelements: c_long,     /* I - Number of elements to return         */
    nulval: *const c_void, /* I - Ptr to value to use as UNDEF         */
    array: *mut c_void,    /* O - Array of results                     */
    anynul: *mut c_int,    /* O - Were any UNDEFs encountered?         */
    status: *mut c_int,    /* O - Error status                         */
) -> c_int {
    unsafe {
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);
        let anynul = anynul.as_mut().expect(NULL_MSG);

        let array = std::slice::from_raw_parts_mut(array as *mut u8, nelements as usize);

        raw_to_slice!(expr);

        ffcrow_safe(
            fptr, datatype, expr, firstrow, nelements, nulval, array, anynul, status,
        )
    }
}

/*---------------------------------------------------------------------------*/
/// Calculate an expression for the indicated rows of a table, returning
/// the results, cast as datatype (TSHORT, TDOUBLE, etc), in array.  If
/// nulval==NULL, UNDEFs will be zeroed out.  For vector results, the number
/// of elements returned may be less than nelements if nelements is not an
/// even multiple of the result dimension.  Call fftexp to obtain the
/// dimensions of the results.
pub fn ffcrow_safe(
    fptr: &mut fitsfile,   /* I - Input FITS file                      */
    datatype: c_int,       /* I - Datatype to return results as        */
    expr: &[c_char],       /* I - Arithmetic expression                */
    firstrow: c_long,      /* I - First row to evaluate                */
    nelements: c_long,     /* I - Number of elements to return         */
    nulval: *const c_void, /* I - Ptr to value to use as UNDEF         */
    array: &mut [u8],      /* O - Array of results                     */
    anynul: &mut c_int,    /* O - Were any UNDEFs encountered?         */
    status: &mut c_int,    /* O - Error status                         */
) -> c_int {
    let mut Info: parseInfo = Default::default();
    let mut naxis: c_int = 0;
    let mut nelem1: c_long = 0;
    let mut naxes: [c_long; MAXDIMS as usize] = [0; MAXDIMS as usize];
    let mut lParse: ParseData = ParseData::default();

    if *status != 0 {
        return *status;
    }

    if ffiprs(
        fptr,
        0,
        expr,
        MAXDIMS,
        &mut Info.datatype,
        &mut nelem1,
        &mut naxis,
        &mut naxes,
        &mut lParse,
        status,
    ) != 0
    {
        ffcprs(&mut lParse);
        return *status;
    }
    if nelem1 < 0 {
        nelem1 = -nelem1;
    }

    if nelements < nelem1 {
        ffcprs(&mut lParse);
        ffpmsg_str("Array not large enough to hold at least one row of data.");
        *status = PARSE_LRG_VECTOR;
        return *status;
    }

    let firstrow = if firstrow > 1 { firstrow } else { 1 };

    if datatype != 0 {
        Info.datatype = datatype;
    }

    Info.dataPtr = array.as_mut_ptr() as *mut c_void;
    Info.nullPtr = nulval as *mut c_void;
    Info.maxRows = nelements / nelem1;
    Info.parseData = &mut lParse;

    let colData_slice = &mut lParse.colData[..];
    if ffiter_safe(
        lParse.nCols,
        colData_slice,
        firstrow - 1,
        0,
        fits_parser_workfn,
        (&mut Info) as *mut parseInfo as *mut c_void,
        status,
    ) == -1
    {
        *status = 0; /* -1 indicates exitted without error before end... OK */
    }

    *anynul = Info.anyNull;
    ffcprs(&mut lParse);
    *status
}

/*--------------------------------------------------------------------------*/
/// Evaluate an expression for all rows of a table.  Call ffcalc_rng with
/// a row range of 1-MAX.              
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffcalc(
    infptr: *mut fitsfile,  /* I - Input FITS file                      */
    expr: *const c_char,    /* I - Arithmetic expression                */
    outfptr: *mut fitsfile, /* I - Output fits file                     */
    parName: *const c_char, /* I - Name of output parameter             */
    parInfo: *const c_char, /* I - Extra information on parameter       */
    status: *mut c_int,     /* O - Error status                         */
) -> c_int {
    unsafe {
        let infptr = infptr.as_mut().expect(NULL_MSG);
        let outfptr = outfptr.as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);

        raw_to_slice!(expr);
        raw_to_slice!(parName);
        raw_to_slice!(parInfo);

        ffcalc_safe(infptr, expr, outfptr, parName, parInfo, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Evaluate an expression for all rows of a table.  Call ffcalc_rng with
/// a row range of 1-MAX.              
pub fn ffcalc_safe(
    infptr: &mut fitsfile,  /* I - Input FITS file                      */
    expr: &[c_char],        /* I - Arithmetic expression                */
    outfptr: &mut fitsfile, /* I - Output fits file                     */
    parName: &[c_char],     /* I - Name of output parameter             */
    parInfo: &[c_char],     /* I - Extra information on parameter       */
    status: &mut c_int,     /* O - Error status                         */
) -> c_int {
    let mut start: c_long = 1;
    let mut end: c_long = LONG_MAX;

    ffcalc_rng_safe(
        infptr, expr, outfptr, parName, parInfo, 1, &mut start, &mut end, status,
    )
}

/*--------------------------------------------------------------------------*/
/// Evaluate an expression using the data in the input FITS file and place  
/// the results into either a column or keyword in the output fits file,    
/// depending on the value of parName (keywords normally prefixed with '#')
/// and whether the expression evaluates to a constant or a table column.   
/// The logic is as follows:                                                
///    (1) If a column exists with name, parName, put results there.        
///    (2) If parName starts with '#', as in #NAXIS, put result there,      
///        with parInfo used as the comment. If expression does not evaluate
///        to a constant, flag an error.                                    
///    (3) If a keyword exists with name, parName, and expression is a      
///        constant, put result there, using parInfo as the new comment.    
///    (4) Else, create a new column with name parName and TFORM parInfo.   
///        If parInfo is NULL, use a default data type for the column.      
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffcalc_rng(
    infptr: *mut fitsfile,  /* I - Input FITS file                  */
    expr: *const c_char,    /* I - Arithmetic expression            */
    outfptr: *mut fitsfile, /* I - Output fits file                 */
    parName: *const c_char, /* I - Name of output parameter         */
    parInfo: *const c_char, /* I - Extra information on parameter   */
    nRngs: c_int,           /* I - Row range info                   */
    start: *mut c_long,     /* I - Row range info                   */
    end: *mut c_long,       /* I - Row range info                   */
    status: *mut c_int,     /* O - Error status                     */
) -> c_int {
    unsafe {
        raw_to_slice!(expr);
        raw_to_slice!(parName);
        raw_to_slice!(parInfo);
        let infptr = infptr.as_mut().expect(NULL_MSG);
        let outfptr = outfptr.as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);
        let start = start.as_mut().expect(NULL_MSG);
        let end = end.as_mut().expect(NULL_MSG);

        ffcalc_rng_safe(
            infptr, expr, outfptr, parName, parInfo, nRngs, start, end, status,
        )
    }
}

/*--------------------------------------------------------------------------*/
/// Evaluate an expression using the data in the input FITS file and place  
/// the results into either a column or keyword in the output fits file,    
/// depending on the value of parName (keywords normally prefixed with '#')
/// and whether the expression evaluates to a constant or a table column.   
/// The logic is as follows:                                                
///    (1) If a column exists with name, parName, put results there.        
///    (2) If parName starts with '#', as in #NAXIS, put result there,      
///        with parInfo used as the comment. If expression does not evaluate
///        to a constant, flag an error.                                    
///    (3) If a keyword exists with name, parName, and expression is a      
///        constant, put result there, using parInfo as the new comment.    
///    (4) Else, create a new column with name parName and TFORM parInfo.   
///        If parInfo is NULL, use a default data type for the column.      
pub fn ffcalc_rng_safe(
    infptr: &mut fitsfile,  /* I - Input FITS file                  */
    expr: &[c_char],        /* I - Arithmetic expression            */
    outfptr: &mut fitsfile, /* I - Output fits file                 */
    parName: &[c_char],     /* I - Name of output parameter         */
    parInfo: &[c_char],     /* I - Extra information on parameter   */
    nRngs: c_int,           /* I - Row range info                   */
    start: &mut c_long,     /* I - Row range info                   */
    end: &mut c_long,       /* I - Row range info                   */
    status: &mut c_int,     /* O - Error status                     */
) -> c_int {
    unsafe {
        let mut Info: parseInfo = Default::default();
        let mut naxis: c_int = 0;
        let constant: c_int;
        let mut typecode: c_int = 0;
        let mut newNullKwd: c_int = 0;
        let mut nelem: c_long = 0;
        let mut naxes: [c_long; MAXDIMS as usize] = [0; MAXDIMS as usize];
        let mut repeat: c_long = 0;
        let mut width: c_long = 0;
        let col_cnt: c_int;
        let mut colNo: c_int;
        let result: &mut Node;
        let mut card: [c_char; FLEN_CARD] = [0; FLEN_CARD];
        let mut tform: [c_char; 16] = [0; 16];
        let mut nullKwd: [c_char; 9] = [0; 9];
        let mut tdimKwd: [c_char; 9] = [0; 9];
        let mut lParse: ParseData = ParseData::default();

        let mut parInfo = parInfo;

        if *status != 0 {
            return *status;
        }

        if ffiprs(
            infptr,
            0,
            expr,
            MAXDIMS,
            &mut Info.datatype,
            &mut nelem,
            &mut naxis,
            &mut naxes,
            &mut lParse,
            status,
        ) != 0
        {
            ffcprs(&mut lParse);
            return *status;
        }
        if nelem < 0 {
            constant = 1;
            nelem = -nelem;
        } else {
            constant = 0;
        }

        Info.parseData = &mut lParse;
        /*  Case (1): If column exists put it there  */

        colNo = 0;
        ffpmrk_safe(); /* prevent lack of column name from sullying the stack */
        ffgcno_safe(outfptr, CASEINSEN as c_int, parName, &mut colNo, status);
        ffcmsg_safe();
        if *status == 0 {
            /*  Output column doesn't exist.  Test for keyword. */

            /* Case (2): Does parName indicate result should be put into keyword */

            *status = 0;
            if parName[0] == b'#' as c_char {
                if constant == 0 {
                    ffcprs(&mut lParse);
                    ffpmsg_str("Cannot put tabular result into keyword (ffcalc)");
                    *status = PARSE_BAD_TYPE;
                    return *status;
                }
                let parName = &parName[1..]; /* Advance past '#' */
                if (fits_strcasecmp(parName, cs!(c"HISTORY")) == 0
                    || fits_strcasecmp(parName, cs!(c"COMMENT")) == 0)
                    && Info.datatype != TSTRING
                {
                    ffcprs(&mut lParse);
                    ffpmsg_str("HISTORY and COMMENT values must be strings (ffcalc)");
                    *status = PARSE_BAD_TYPE;
                    return *status;
                }
            } else if constant != 0 {
                /* Case (3): Does a keyword named parName already exist */

                if ffgcrd_safe(outfptr, parName, &mut card, status) == KEY_NO_EXIST {
                    colNo = -1;
                } else if *status != 0 {
                    ffcprs(&mut lParse);
                    return *status;
                }
            } else {
                colNo = -1;
            }

            if colNo < 0 {
                /* Case (4): Create new column */

                *status = 0;
                ffgncl_safe(outfptr, &mut colNo, status);
                colNo += 1;
                if parInfo.is_empty() || parInfo[0] == 0 {
                    /*  Figure out best default column type  */
                    if lParse.hdutype == BINARY_TBL {
                        int_snprintf!(&mut tform, 15, "{}", nelem);
                        match Info.datatype {
                            TLOGICAL => {
                                strcat(tform.as_mut_ptr(), (c"L").as_ptr());
                            }
                            TLONG => {
                                strcat(tform.as_mut_ptr(), (c"J").as_ptr());
                            }
                            TDOUBLE => {
                                strcat(tform.as_mut_ptr(), (c"D").as_ptr());
                            }
                            TSTRING => {
                                strcat(tform.as_mut_ptr(), (c"A").as_ptr());
                            }
                            TBIT => {
                                strcat(tform.as_mut_ptr(), (c"X").as_ptr());
                            }
                            TLONGLONG => {
                                strcat(tform.as_mut_ptr(), (c"K").as_ptr());
                            }
                            _ => {}
                        }
                    } else {
                        match Info.datatype {
                            TLOGICAL => {
                                ffcprs(&mut lParse);
                                ffpmsg_str("Cannot create LOGICAL column in ASCII table");
                                *status = NOT_BTABLE;
                                return *status;
                            }
                            TLONG => {
                                strcpy(tform.as_mut_ptr(), (c"I11").as_ptr());
                            }
                            TDOUBLE => {
                                strcpy(tform.as_mut_ptr(), (c"D23.15").as_ptr());
                            }
                            TSTRING | TBIT => {
                                int_snprintf!(&mut tform, 16, "A{}", nelem);
                            }
                            _ => {}
                        }
                    }
                    parInfo = &tform;
                } else if !((parInfo[0] as u8) as char).is_ascii_digit()
                    && lParse.hdutype == BINARY_TBL
                {
                    if Info.datatype == TBIT && parInfo[0] == b'B' as c_char {
                        nelem = (nelem + 7) / 8;
                    }
                    int_snprintf!(
                        &mut tform,
                        16,
                        "{}{}",
                        nelem,
                        std::str::from_utf8(cast_slice(parInfo)).unwrap_or("")
                    );
                    parInfo = &tform;
                }
                fficol_safe(outfptr, colNo, parName, parInfo, status);

                if naxis > 1 {
                    ffptdm_safe(outfptr, colNo, naxis, &naxes[..naxis as usize], status);
                }

                /*  Setup TNULLn keyword in case NULLs are encountered  */

                ffkeyn_safe(cs!(c"TNULL"), colNo, &mut nullKwd, status);
                if ffgcrd_safe(outfptr, &nullKwd, &mut card, status) == KEY_NO_EXIST {
                    *status = 0;
                    if lParse.hdutype == BINARY_TBL {
                        let mut nullVal: LONGLONG = 0;
                        fits_binary_tform(
                            parInfo,
                            Some(&mut typecode),
                            Some(&mut repeat),
                            Some(&mut width),
                            status,
                        );
                        if typecode == TBYTE {
                            nullVal = UCHAR_MAX;
                        } else if typecode == TSHORT {
                            nullVal = SHRT_MIN;
                        } else if typecode == TINT {
                            nullVal = INT_MIN;
                        } else if typecode == TLONG {
                            if std::mem::size_of::<c_long>() == 8
                                && std::mem::size_of::<c_int>() == 4
                            {
                                nullVal = INT_MIN;
                            } else {
                                nullVal = LONG_MIN;
                            }
                        } else if typecode == TLONGLONG {
                            nullVal = LONGLONG_MIN;
                        }

                        if nullVal != 0 {
                            ffpkyj_safe(
                                outfptr,
                                &nullKwd,
                                nullVal,
                                Some(cs!(c"Null value")),
                                status,
                            );
                            fits_set_btblnull(outfptr, colNo, nullVal, status);
                            newNullKwd = 1;
                        }
                    } else if lParse.hdutype == ASCII_TBL {
                        ffpkys_safe(
                            outfptr,
                            &nullKwd,
                            cs!(c"NULL"),
                            Some(cs!(c"Null value string")),
                            status,
                        );
                        fits_set_atblnull(outfptr, colNo, cs!(c"NULL"), status);
                        newNullKwd = 1;
                    }
                }
            }
        } else if *status != 0 {
            ffcprs(&mut lParse);
            return *status;
        } else {
            /********************************************************/
            /*  Check if a TDIM keyword should be written/updated.  */
            /********************************************************/

            ffkeyn_safe(cs!(c"TDIM"), colNo, &mut tdimKwd, status);
            ffgcrd_safe(outfptr, &tdimKwd, &mut card, status);
            if *status == 0 {
                /*  TDIM exists, so update it with result's dimension  */
                ffptdm_safe(outfptr, colNo, naxis, &naxes[..naxis as usize], status);
            } else if *status == KEY_NO_EXIST {
                /*  TDIM does not exist, so clear error stack and     */
                /*  write a TDIM only if result is multi-dimensional  */
                *status = 0;
                ffcmsg_safe();
                if naxis > 1 {
                    ffptdm_safe(outfptr, colNo, naxis, &naxes[..naxis as usize], status);
                }
            }
            if *status != 0 {
                /*  Either some other error happened in ffgcrd   */
                /*  or one happened in ffptdm                    */
                ffcprs(&mut lParse);
                return *status;
            }
        }

        if colNo > 0 {
            /*  Output column exists (now)... put results into it  */

            let mut anyNull: c_int = 0;
            let mut nPerLp: c_int;
            let mut i: c_int;
            let mut totaln: c_long = 0;

            ffgkyj_safe(infptr, cs!(c"NAXIS2"), &mut totaln, None, status);

            /*************************************/
            /* Create new iterator Output Column */
            /*************************************/

            col_cnt = lParse.nCols;
            if fits_parser_allocateCol(&mut lParse, col_cnt, status) != 0 {
                ffcprs(&mut lParse);
                return *status;
            }

            fits_iter_set_by_num_safe(
                &mut lParse.colData[col_cnt as usize],
                outfptr,
                colNo,
                0,
                OutputCol,
            );

            lParse.nCols += 1;

            for i in 0..nRngs {
                Info.dataPtr = std::ptr::null_mut();
                let start_slice =
                    unsafe { std::slice::from_raw_parts(start as *const c_long, nRngs as usize) };
                let end_slice =
                    unsafe { std::slice::from_raw_parts(end as *const c_long, nRngs as usize) };
                Info.maxRows = end_slice[i as usize] - start_slice[i as usize] + 1;

                /*
                   If there is only 1 range, and it includes all the rows,
                   and there are 10 or more rows, then set nPerLp = 0 so
                   that the iterator function will dynamically choose the
                   most efficient number of rows to process in each loop.
                   Otherwise, set nPerLp to the number of rows in this range.
                */

                if (Info.maxRows >= 10)
                    && (nRngs == 1)
                    && (start_slice[0] == 1)
                    && (end_slice[0] == totaln)
                {
                    nPerLp = 0;
                } else {
                    nPerLp = Info.maxRows as c_int;
                }

                let colData_slice = &mut lParse.colData[..];
                if ffiter_safe(
                    lParse.nCols,
                    colData_slice,
                    start_slice[i as usize] - 1,
                    nPerLp as c_long,
                    fits_parser_workfn,
                    &Info as *const _ as *mut c_void,
                    status,
                ) == -1
                {
                    *status = 0;
                } else if *status != 0 {
                    ffcprs(&mut lParse);
                    return *status;
                }
                if Info.anyNull != 0 {
                    anyNull = 1;
                }
            }

            if newNullKwd != 0 && anyNull == 0 {
                ffdkey_safe(outfptr, &nullKwd, status);
            }
        } else {
            /* Put constant result into keyword */

            result = &mut lParse.Nodes[lParse.resultNode as usize];
            match Info.datatype {
                TDOUBLE => {
                    ffukyd_safe(
                        outfptr,
                        parName,
                        unsafe { result.value.data.dbl },
                        15,
                        Some(parInfo),
                        status,
                    );
                }
                TLONG => {
                    ffukyj_safe(
                        outfptr,
                        parName,
                        unsafe { result.value.data.lng },
                        Some(parInfo),
                        status,
                    );
                }
                TLOGICAL => {
                    ffukyl_safe(
                        outfptr,
                        parName,
                        unsafe { result.value.data.log } as i32,
                        Some(parInfo),
                        status,
                    );
                }
                TBIT | TSTRING => {
                    if fits_strcasecmp(parName, cs!(c"HISTORY")) == 0 {
                        ffphis_safe(outfptr, unsafe { &result.value.data.astr }, status);
                    } else if fits_strcasecmp(parName, cs!(c"COMMENT")) == 0 {
                        ffpcom_safe(outfptr, unsafe { &result.value.data.astr }, status);
                    } else {
                        ffukys_safe(
                            outfptr,
                            parName,
                            unsafe { &result.value.data.astr },
                            Some(parInfo),
                            status,
                        );
                    }
                }
                _ => {}
            }
        }

        ffcprs(&mut lParse);
        *status
    }
}

/*--------------------------------------------------------------------------*/
/// Evaluate the given expression and return information on the (*result).
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn fftexp(
    fptr: *mut fitsfile,  /* I - Input FITS file                     */
    expr: *const c_char,  /* I - Arithmetic expression               */
    maxdim: c_int,        /* I - Max Dimension of naxes              */
    datatype: *mut c_int, /* O - Data type of result                 */
    nelem: *mut c_long,   /* O - Vector length of result             */
    naxis: *mut c_int,    /* O - # of dimensions of result           */
    naxes: *mut c_long,   /* O - Size of each dimension              */
    status: *mut c_int,   /* O - Error status                        */
) -> c_int {
    unsafe {
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);
        let datatype = datatype.as_mut().expect(NULL_MSG);
        let nelem = nelem.as_mut().expect(NULL_MSG);
        let naxis = naxis.as_mut().expect(NULL_MSG);
        let naxes = std::slice::from_raw_parts_mut(naxes, maxdim as usize);

        raw_to_slice!(expr);

        fftexp_safe(fptr, expr, maxdim, datatype, nelem, naxis, naxes, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Evaluate the given expression and return information on the (*result).
pub fn fftexp_safe(
    fptr: &mut fitsfile,  /* I - Input FITS file                     */
    expr: &[c_char],      /* I - Arithmetic expression               */
    maxdim: c_int,        /* I - Max Dimension of naxes              */
    datatype: &mut c_int, /* O - Data type of result                 */
    nelem: &mut c_long,   /* O - Vector length of result             */
    naxis: &mut c_int,    /* O - # of dimensions of result           */
    naxes: &mut [c_long], /* O - Size of each dimension              */
    status: &mut c_int,   /* O - Error status                        */
) -> c_int {
    let mut lParse: ParseData = ParseData::default();

    ffiprs(
        fptr,
        0,
        expr,
        maxdim,
        datatype,
        nelem,
        naxis,
        naxes,
        &mut lParse,
        status,
    );
    ffcprs(&mut lParse);
    *status
}

/*--------------------------------------------------------------------------*/
/// Initialize the parser and determine what type of result the expression
/// produces.
pub(crate) fn ffiprs(
    fptr: &mut fitsfile,    /* I - Input FITS file                     */
    compressed: c_int,      /* I - Is FITS file hkunexpanded?          */
    expr: &[c_char],        /* I - Arithmetic expression               */
    maxdim: c_int,          /* I - Max Dimension of naxes              */
    datatype: &mut c_int,   /* O - Data type of result                 */
    nelem: &mut c_long,     /* O - Vector length of result             */
    naxis: &mut c_int,      /* O - # of dimensions of result           */
    naxes: &mut [c_long],   /* O - Size of each dimension              */
    lParse: &mut ParseData, /* O - parser status                       */
    status: &mut c_int,     /* O - Error status                        */
) -> c_int {
    unsafe {
        let i: c_int = 0;
        let mut lexpr: c_int = 0;
        let mut tstatus: c_int = 0;
        let mut xaxis: c_int = 0;
        let mut bitpix: c_int = 0;
        let mut xaxes: [c_long; 9] = [0; 9];

        if *status != 0 {
            return *status;
        }

        /* make sure all internal structures for this HDU are current */
        if ffrdef_safe(fptr, status) != 0 {
            return *status;
        }

        /*  Initialize the Parser structure  */

        /* Unfortunately we need to preserve the pixFilter value since it is pre-set when ffiprs() is called */
        let pixFilter: *mut PixelFilter = lParse.pixFilter;

        lParse.reset();

        lParse.pixFilter = pixFilter;

        lParse.def_fptr = fptr;
        lParse.compressed = compressed;
        lParse.nCols = 0;
        lParse.colData = Vec::new();
        lParse.varData = Vec::new();
        lParse.getData = Some(find_column);
        lParse.loadData = Some(load_column);
        lParse.Nodes = Vec::new();
        lParse.nNodesAlloc = 0;
        lParse.nNodes = 0;
        lParse.hdutype = 0;
        lParse.status = 0;

        fits_get_hdu_type(fptr, &mut (lParse.hdutype), status);

        if lParse.hdutype == IMAGE_HDU {
            fits_get_img_param(
                fptr,
                9,
                Some(&mut bitpix),
                Some(&mut xaxis),
                Some(&mut xaxes[..]),
                status,
            );
            if *status != 0 {
                ffpmsg_str("ffiprs: unable to get image dimensions");
                return *status;
            }

            lParse.totalRows = if xaxis > 0 { 1 } else { 0 };

            for i in 0..(xaxis as usize) {
                lParse.totalRows *= xaxes[i];
            }

            if DEBUG_PIXFILTER != 0 {
                println!("naxis={}, (*lParse).totalRows={}", xaxis, lParse.totalRows,);
            }
        } else if ffgkyj_safe(
            fptr,
            cs!(c"NAXIS2"),
            &mut lParse.totalRows,
            None,
            &mut tstatus,
        ) != 0
        {
            /* this might be a 1D or null image with no NAXIS2 keyword */
            lParse.totalRows = 0;
        }

        /*  Copy expression into parser... read from file if necessary  */

        if expr[0] == b'@' as c_char {
            if ffimport_file_safer(&expr[1..], &mut lParse.expr, status) != 0 {
                return *status;
            }
            lexpr = strlen(lParse.expr) as c_int;
        } else {
            lexpr = strlen_safe(expr) as c_int;
            lParse.expr =
                malloc(((2 + lexpr as usize) * size_of::<c_char>()) as usize) as *mut c_char;
            strcpy(lParse.expr, expr.as_ptr());
        }
        strcat(
            (lParse.expr as *mut c_char).wrapping_add(lexpr as usize),
            c"\n".as_ptr(),
        );
        lParse.index = 0;
        lParse.is_eobuf = 0;

        /*  Parse the expression, building the Nodes and determing  */
        /*  which columns are needed and what data type is returned  */

        let mut yylex_scanner: Option<Box<yyguts_t>> = None; /* Used internally by FLEX lexer */

        fits_parser_yylex_init_extra(lParse, &mut yylex_scanner);
        fits_parser_yyrestart(ptr::null_mut(), yylex_scanner.as_deref_mut().unwrap());
        *status = fits_parser_yyparse(yylex_scanner.as_deref_mut().unwrap(), lParse);
        fits_parser_yylex_destroy(yylex_scanner.unwrap());

        if *status != 0 {
            *status = PARSE_SYNTAX_ERR;
            return *status;
        }

        /*  Check results  */
        *status = lParse.status;
        if *status != 0 {
            return *status;
        }

        if lParse.nNodes == 0 {
            ffpmsg_str("Blank expression");
            *status = PARSE_SYNTAX_ERR;
            return *status;
        }
        if lParse.nCols == 0 {
            if lParse.colData.try_reserve_exact(1).is_err() {
                ffpmsg_str("memory allocation failed (ffiprs)");
                *status = MEMORY_ALLOCATION;
                return *status;
            } else {
                lParse.colData.resize(1, iteratorCol::default());
            }

            /* This allows iterator to know value of fptr when no columns are referenced   */
            (lParse.colData[0]).fptr = fptr;
        }

        let result: &mut Node = &mut lParse.Nodes[lParse.resultNode as usize];

        lParse.nAxis = result.value.naxis;
        *naxis = lParse.nAxis;
        lParse.nElements = result.value.nelem;
        *nelem = lParse.nElements;

        for i in 0..cmp::max(*naxis, maxdim) {
            lParse.nAxes[i as usize] = result.value.naxes[i as usize];
            naxes[i as usize] = lParse.nAxes[i as usize];
        }

        match result.ntype {
            val if val == fits_parser_yytokentype::BOOLEAN as c_int => *datatype = TLOGICAL,

            val if val == fits_parser_yytokentype::LONG as c_int => *datatype = TLONG,

            val if val == fits_parser_yytokentype::DOUBLE as c_int => *datatype = TDOUBLE,

            val if val == fits_parser_yytokentype::BITSTR as c_int => *datatype = TBIT,

            val if val == fits_parser_yytokentype::STRING as c_int => *datatype = TSTRING,

            _ => {
                *datatype = 0;
                ffpmsg_str("Bad return data type");
                lParse.status = PARSE_BAD_TYPE;
                *status = lParse.status;
            }
        }
        lParse.datatype = *datatype;
        FREE!(lParse.expr);

        if result.operation == CONST_OP {
            *nelem = -*nelem;
        }
        *status
    }
}

/*---------------------------------------------------------------------------*/
/// Clear the parser, making it ready to accept a new expression.
fn ffcprs(lParse: &mut ParseData) {
    unsafe {
        let col: c_int = 0;
        let mut node: c_int = 0;
        let mut i: usize = 0;

        if lParse.nCols > 0 {
            lParse.colData.clear();

            for col in 0..lParse.nCols {
                if lParse.varData[col as usize].undef.is_null() {
                    continue;
                }
                if lParse.varData[col as usize].dtype == fits_parser_yytokentype::BITSTR as c_int {
                    let data_ptr = lParse.varData[col as usize].data as *mut *mut c_char;
                    let mut first_ptr = *data_ptr;
                    FREE!(first_ptr);
                }
                free(lParse.varData[col as usize].undef as *mut c_void);
            }
            lParse.varData.clear();
            lParse.nCols = 0;
        } else if !lParse.colData.is_empty() {
            /* Special case if colData needed to be created with no columns */
            lParse.colData.clear();
        }

        if lParse.nNodes > 0 {
            node = lParse.nNodes;
            while node != 0 {
                node -= 1;
                if (lParse.Nodes[node as usize]).operation == gtifilt_fct as c_int {
                    i = lParse.Nodes[node as usize].SubNodes[0];
                    if !(lParse.Nodes[i]).value.data.ptr.is_null() {
                        let mut data_ptr = (lParse.Nodes[i]).value.data.ptr;
                        FREE!(data_ptr);
                    }
                } else if (lParse.Nodes[node as usize]).operation == regfilt_fct as c_int {
                    i = (lParse.Nodes[node as usize]).SubNodes[0];
                    fits_free_region(Box::from_raw(
                        (lParse.Nodes[i]).value.data.ptr as *mut SAORegion,
                    ));
                }
            }
            lParse.nNodes = 0;
        }

        lParse.Nodes = Vec::new(); // Frees

        lParse.hdutype = ANY_HDU;
        lParse.pixFilter = std::ptr::null_mut();
        lParse.nDataRows = 0;
        lParse.nPrevDataRows = 0;
    }
}

/*---------------------------------------------------------------------------*/
/// Iterator work function which calls the parser and copies the results
/// into either an OutputCol or a data pointer supplied in the userPtr
/// structure.
extern "C" fn fits_parser_workfn(
    totalrows: c_long,         /* I - Total rows to be processed     */
    offset: c_long,            /* I - Number of rows skipped at start*/
    firstrow: c_long,          /* I - First row of this iteration    */
    nrows: c_long,             /* I - Number of rows in this iter    */
    nCols: c_int,              /* I - Number of columns in use       */
    colData: *mut iteratorCol, /* IO- Column information/data        */
    userPtr: *mut c_void,      /* I - Data handling instructions     */
) -> c_int {
    let colData = unsafe { std::slice::from_raw_parts_mut(colData, nCols as usize) };

    fits_parser_workfn_safe(totalrows, offset, firstrow, nrows, nCols, colData, userPtr)
}

/*---------------------------------------------------------------------------*/
/// Iterator work function which calls the parser and copies the results
/// into either an OutputCol or a data pointer supplied in the userPtr
/// structure.
fn fits_parser_workfn_safe(
    totalrows: c_long,           /* I - Total rows to be processed     */
    offset: c_long,              /* I - Number of rows skipped at start*/
    mut firstrow: c_long,        /* I - First row of this iteration    */
    mut nrows: c_long,           /* I - Number of rows in this iter    */
    nCols: c_int,                /* I - Number of columns in use       */
    colData: &mut [iteratorCol], /* IO- Column information/data        */
    userPtr: *mut c_void,        /* I - Data handling instructions     */
) -> c_int {
    unsafe {
        let mut status: c_int = 0;
        let mut constant: c_int = 0;
        let mut anyNullThisTime: c_int = 0;
        let mut jj: c_long = 0;
        let kk: c_long = 0;
        let mut idx: c_long = 0;
        let mut remain: c_long = 0;
        let mut ntodo: c_long = 0;
        let mut result: &mut Node;
        let mut outcol: &mut iteratorCol;
        let lParse: &mut ParseData = unsafe {
            (userPtr as *mut parseInfo)
                .as_mut()
                .unwrap()
                .parseData
                .as_mut()
                .unwrap()
        };
        let pv: &mut ParseStatusVariables =
            unsafe { &mut (userPtr as *mut parseInfo).as_mut().unwrap().parseVariables };

        /* declare variables static to preserve their values between calls */
        let mut zeros: [c_long; 4] = [0, 0, 0, 0];

        if DEBUG_PIXFILTER != 0 {
            println!(
                "fits_parser_workfn(total={}, offset={}, first={}, rows={}, cols={})",
                totalrows, offset, firstrow, nrows, nCols
            );
        }
        /*--------------------------------------------------------*/
        /*  Initialization procedures: execute on the first call  */
        /*--------------------------------------------------------*/

        if firstrow == offset + 1 {
            (pv.userInfo) = userPtr as *mut parseInfo;
            (*(pv.userInfo)).anyNull = 0;

            /* Unfortunately there are two copies of the iterator columns,
            one inside the parser and one outside maintained by the
            higher level.  (This could happen if the histogramming
            routines are binning multiple columns, and so there are
            multiple parsers being managed at one time.) Upon the first
            call we make sure they match */
            for jj in 0..(nCols as usize) {
                lParse.colData[jj].repeat = colData[jj].repeat;
            }

            if (*(pv.userInfo)).maxRows > 0 {
                (*(pv.userInfo)).maxRows = cmp::min(totalrows, (*(pv.userInfo)).maxRows);
            } else if (*(pv.userInfo)).maxRows < 0 {
                (*(pv.userInfo)).maxRows = totalrows;
            } else {
                (*(pv.userInfo)).maxRows = nrows;
            }

            (pv.lastRow) = firstrow + (*(pv.userInfo)).maxRows - 1;

            /* dataPtr == NULL indicates an iterator-derived column, which
            means that the first value will be a null value and the remaining
            values will be the where the outputs are placed */
            if (*(pv.userInfo)).dataPtr.is_null() {
                outcol = &mut colData[(nCols - 1) as usize]; // Re-bind
                if outcol.iotype == InputCol {
                    ffpmsg_str("Output column for parser results not found!");
                    return PARSE_NO_OUTPUT;
                }
                /* Data gets set later */
                (pv.Null) = outcol.array;
                (*(pv.userInfo)).datatype = outcol.datatype;

                /* Check for a TNULL/BLANK keyword for output column/image */

                status = 0;
                (pv.jnull) = 0;
                if lParse.hdutype == IMAGE_HDU {
                    if (*(lParse.pixFilter)).blank != 0 {
                        (pv.jnull) = (*(lParse.pixFilter)).blank as LONGLONG;
                    }
                } else {
                    if outcol.iotype != TemporaryCol {
                        let mut jj_tmp: c_int = 0;
                        ffgknjj_safe(
                            outcol.fptr.as_mut().unwrap(),
                            cs!(c"TNULL"),
                            outcol.colnum,
                            1,
                            &mut [(pv.jnull)],
                            &mut jj_tmp,
                            &mut status,
                        );
                        jj = jj_tmp as c_long;
                    }

                    if status == BAD_INTKEY || outcol.iotype == TemporaryCol {
                        /*  Probably ASCII table with text TNULL keyword  */
                        match (*(pv.userInfo)).datatype {
                            TSHORT => {
                                (pv.jnull) = SHRT_MIN as LONGLONG;
                            }
                            TINT => {
                                (pv.jnull) = INT_MIN as LONGLONG;
                            }
                            TLONG => {
                                (pv.jnull) = LONG_MIN as LONGLONG;
                            }
                            _ => {}
                        }
                    }
                }
                (pv.repeat) = outcol.repeat;
                /*
                          if (DEBUG_PIXFILTER)
                            printf("fits_parser_workfn: using null value %ld\n", (pv.jnull));
                */
            } else {
                /* This clause applies if the user is passing user-allocated
                data arrays, which is where the data will be placed.  This
                means they should also be passing null values */
                (pv.Data) = (*(pv.userInfo)).dataPtr;
                (pv.Null) = if !(*(pv.userInfo)).nullPtr.is_null() {
                    (*(pv.userInfo)).nullPtr
                } else {
                    zeros.as_ptr() as *mut c_void
                };
                (pv.repeat) = lParse.Nodes[lParse.resultNode as usize].value.nelem;
            }

            /* Determine the size of each element of the returned result */

            match (*(pv.userInfo)).datatype {
                TBIT | TLOGICAL | TBYTE => {
                    (pv.datasize) = size_of::<c_char>() as c_int;
                }
                TSHORT => {
                    (pv.datasize) = size_of::<c_short>() as c_int;
                }
                TINT => {
                    (pv.datasize) = size_of::<c_int>() as c_int;
                }
                TLONG => {
                    (pv.datasize) = size_of::<c_long>() as c_int;
                }
                TLONGLONG => {
                    (pv.datasize) = size_of::<LONGLONG>() as c_int;
                }
                TFLOAT => {
                    (pv.datasize) = size_of::<f32>() as c_int;
                }
                TDOUBLE => {
                    (pv.datasize) = size_of::<f64>() as c_int;
                }
                TSTRING => {
                    (pv.datasize) = size_of::<*mut c_char>() as c_int;
                }
                _ => {}
            }

            /* Determine the size of each element of the calculated result */
            /*   (only matters for numeric/logical data)                   */

            match lParse.Nodes[lParse.resultNode as usize].ntype {
                BOOLEAN => {
                    (pv.resDataSize) = size_of::<c_char>() as c_long;
                }
                LONG => {
                    (pv.resDataSize) = size_of::<c_long>() as c_long;
                }
                DOUBLE => {
                    (pv.resDataSize) = size_of::<f64>() as c_long;
                }
            }
        }

        /*-------------------------------------------*/
        /*  Main loop: process all the rows of data  */
        /*-------------------------------------------*/

        /*  If writing to output column, set first element to appropriate  */
        /*  null value.  If no NULLs encounter, zero out before returning. */
        /*
                  if (DEBUG_PIXFILTER)
                    printf("fits_parser_workfn: using null value %ld\n", (pv.jnull));
        */

        if (*(pv.userInfo)).dataPtr.is_null() {
            outcol = &mut colData[(nCols - 1) as usize]; // Re-bind

            /* First, reset Data pointer to start of output array, plus 1
            because the 0th element is the null value (cute undocumented
            feature of the iterator!) */
            (pv.Data) =
                (outcol.array as *mut c_char).add(pv.datasize.try_into().unwrap()) as *mut c_void;

            /* A TemporaryCol with null value specified explicitly */
            if outcol.iotype == TemporaryCol && !(*(pv.userInfo)).nullPtr.is_null() {
                pv.Null = (*(pv.userInfo)).nullPtr;
            } else {
                /* ... or an OutputCol or TemporaryCol with no explicit null */
                match (*(pv.userInfo)).datatype {
                    TLOGICAL => {
                        *(pv.Null as *mut c_char) = b'U' as c_char;
                    }
                    TBYTE => {
                        *(pv.Null as *mut c_char) = (pv.jnull) as c_char;
                    }
                    TSHORT => {
                        *(pv.Null as *mut c_short) = (pv.jnull) as c_short;
                    }
                    TINT => {
                        *(pv.Null as *mut c_int) = (pv.jnull) as c_int;
                    }
                    TLONG => {
                        *(pv.Null as *mut c_long) = (pv.jnull) as c_long;
                    }
                    TLONGLONG => {
                        *(pv.Null as *mut LONGLONG) = (pv.jnull) as LONGLONG;
                    }
                    TFLOAT => {
                        *(pv.Null as *mut f32) = FLOATNULLVALUE;
                    }
                    TDOUBLE => {
                        *(pv.Null as *mut f64) = DOUBLENULLVALUE;
                    }
                    TSTRING => {
                        (**(pv.Null as *mut *mut c_char)) = 1;
                        *(*(pv.Null as *mut *mut c_char)).add(1) = 0;
                    }
                    _ => {}
                }
            }
        }

        /* Alter nrows in case calling routine didn't want to do all rows */

        let Data0: *mut c_void = pv.Data; /* Record starting point */
        nrows = cmp::min(nrows, (pv.lastRow) - firstrow + 1);

        Setup_DataArrays(lParse, nCols, colData, firstrow, nrows);

        /* Parser allocates arrays for each column and calculation it performs. */
        /* Limit number of rows processed during each pass to reduce memory     */
        /* requirements... In most cases, iterator will limit rows to less      */
        /* than 10000 rows per iteration, so this is really only relevant for    */
        /* hk-compressed files which must be decompressed in memory and sent    */
        /* whole to fits_parser_workfn in a single iteration.                           */

        remain = nrows;
        while remain != 0 {
            ntodo = cmp::min(remain, 10000);
            Evaluate_Parser(lParse, firstrow, ntodo);
            if lParse.status != 0 {
                break;
            }

            firstrow += ntodo;
            remain -= ntodo;

            /*  Copy results into data array  */

            result = &mut lParse.Nodes[lParse.resultNode as usize];
            if result.operation == CONST_OP {
                constant = 1;
            }

            match result.ntype.into() {
                fits_parser_yytokentype::BOOLEAN
                | fits_parser_yytokentype::LONG
                | fits_parser_yytokentype::DOUBLE => {
                    if constant != 0 {
                        let undef: c_char = 0;
                        for kk in 0..ntodo {
                            for jj in 0..pv.repeat {
                                ffcvtn(
                                    lParse.datatype,
                                    &(result.value.data) as *const _ as *const c_void,
                                    &undef,
                                    result.value.nelem, /* 1 */
                                    (*(pv.userInfo)).datatype,
                                    pv.Null,
                                    (pv.Data as *mut c_char).add(
                                        ((kk * (pv.repeat) + jj) * (pv.datasize as c_long))
                                            .try_into()
                                            .unwrap(),
                                    ) as *mut c_void,
                                    &mut anyNullThisTime,
                                    &mut lParse.status,
                                );
                            }
                        }
                    } else {
                        if (pv.repeat) == result.value.nelem {
                            ffcvtn(
                                lParse.datatype,
                                result.value.data.ptr,
                                result.value.undef,
                                result.value.nelem * ntodo,
                                (*(pv.userInfo)).datatype,
                                pv.Null,
                                pv.Data,
                                &mut anyNullThisTime,
                                &mut lParse.status,
                            );
                        } else if result.value.nelem == 1 {
                            for kk in 0..ntodo {
                                for jj in 0..pv.repeat {
                                    ffcvtn(
                                        lParse.datatype,
                                        (result.value.data.ptr as *mut c_char)
                                            .add((kk * (pv.resDataSize)).try_into().unwrap())
                                            as *const c_void,
                                        (result.value.undef as *mut c_char)
                                            .add(kk.try_into().unwrap()),
                                        1,
                                        (*(pv.userInfo)).datatype,
                                        pv.Null,
                                        (pv.Data as *mut c_char).add(
                                            ((kk * (pv.repeat) + jj) * (pv.datasize as c_long))
                                                .try_into()
                                                .unwrap(),
                                        ) as *mut c_void,
                                        &mut anyNullThisTime,
                                        &mut lParse.status,
                                    );
                                }
                            }
                        } else {
                            let nCopy: c_int = cmp::min(pv.repeat, result.value.nelem) as c_int;
                            for kk in 0..ntodo {
                                ffcvtn(
                                    lParse.datatype,
                                    (result.value.data.ptr as *const c_char).add(
                                        (kk * result.value.nelem * (pv.resDataSize))
                                            .try_into()
                                            .unwrap(),
                                    ) as *const c_void,
                                    (result.value.undef as *mut c_char)
                                        .add((kk * result.value.nelem).try_into().unwrap()),
                                    nCopy as c_long,
                                    (*(pv.userInfo)).datatype,
                                    pv.Null,
                                    (pv.Data as *mut c_char).add(
                                        ((kk * (pv.repeat)) * (pv.datasize as c_long))
                                            .try_into()
                                            .unwrap(),
                                    ) as *mut c_void,
                                    &mut anyNullThisTime,
                                    &mut lParse.status,
                                );
                                if (nCopy as c_long) < (pv.repeat) {
                                    memset(
                                        (pv.Data as *mut c_char).add(
                                            ((kk * (pv.repeat) + nCopy as c_long)
                                                * (pv.datasize as c_long))
                                                .try_into()
                                                .unwrap(),
                                        ) as *mut c_void,
                                        0,
                                        (((pv.repeat) - nCopy as c_long) * (pv.datasize as c_long))
                                            as usize,
                                    );
                                }
                            }
                        }
                        if result.operation > 0 {
                            FREE!(result.value.data.ptr);
                        }
                    }
                    if lParse.status == OVERFLOW_ERR {
                        lParse.status = NUM_OVERFLOW;
                        ffpmsg_str(
                            "Numerical overflow while converting expression to necessary datatype",
                        );
                    }
                }

                fits_parser_yytokentype::BITSTR => {
                    match (*(pv.userInfo)).datatype {
                        TBYTE => {
                            idx = -1;
                            for kk in 0..(ntodo) {
                                for jj in 0..result.value.nelem {
                                    if jj % 8 == 0 {
                                        idx += 1;
                                        *((pv.Data as *mut c_char).add(idx.try_into().unwrap())) =
                                            0;
                                    }
                                    if constant != 0 {
                                        if result.value.data.astr[jj as usize] == b'1' as c_char {
                                            *((pv.Data as *mut c_uchar)
                                                .add(idx.try_into().unwrap())) |= 128 >> (jj % 8);
                                        }
                                    } else if *(*(result
                                        .value
                                        .data
                                        .strptr
                                        .add(kk.try_into().unwrap())))
                                    .add(jj.try_into().unwrap())
                                        == b'1' as c_char
                                    {
                                        *((pv.Data as *mut c_uchar)
                                            .add(idx.try_into().unwrap())) |= 128 >> (jj % 8);
                                    }
                                }
                            }
                        }
                        TBIT | TLOGICAL => {
                            if constant != 0 {
                                for kk in 0..ntodo {
                                    for jj in 0..result.value.nelem {
                                        let r = if (result.value.data.astr
                                            [<i64 as TryInto<usize>>::try_into(jj).unwrap()])
                                            == b'1' as c_char
                                        {
                                            1
                                        } else {
                                            0
                                        };
                                        *((pv.Data as *mut c_char).add(
                                            (jj + kk * result.value.nelem).try_into().unwrap(),
                                        )) = r;
                                    }
                                }
                            } else {
                                for kk in 0..ntodo {
                                    for jj in 0..result.value.nelem {
                                        let r = if (*(*(result
                                            .value
                                            .data
                                            .strptr
                                            .add(kk.try_into().unwrap())))
                                        .add(jj.try_into().unwrap()))
                                            == b'1' as c_char
                                        {
                                            1
                                        } else {
                                            0
                                        };
                                        *((pv.Data as *mut c_char).add(
                                            (jj + kk * result.value.nelem).try_into().unwrap(),
                                        )) = r;
                                    }
                                }
                            }
                        }
                        TSTRING => {
                            if constant != 0 {
                                for jj in 0..ntodo {
                                    strcpy(
                                        *((pv.Data as *mut *mut c_char)
                                            .add(jj.try_into().unwrap())),
                                        result.value.data.astr.as_ptr(),
                                    );
                                }
                            } else {
                                for jj in 0..ntodo {
                                    strcpy(
                                        *((pv.Data as *mut *mut c_char)
                                            .add(jj.try_into().unwrap())),
                                        *(result.value.data.strptr.add(jj.try_into().unwrap())),
                                    );
                                }
                            }
                        }
                        _ => {
                            ffpmsg_str("Cannot convert bit expression to desired type.");
                            lParse.status = PARSE_BAD_TYPE;
                        }
                    }
                    if result.operation > 0 {
                        FREE!(*(result.value.data.strptr));
                        FREE!(result.value.data.strptr);
                    }
                }
                fits_parser_yytokentype::STRING => {
                    if (*(pv.userInfo)).datatype == TSTRING {
                        if constant != 0 {
                            for jj in 0..ntodo {
                                strcpy(
                                    *((pv.Data as *mut *mut c_char).add(jj.try_into().unwrap())),
                                    result.value.data.astr.as_ptr(),
                                );
                            }
                        } else {
                            for jj in 0..ntodo {
                                if *(result.value.undef.add(jj.try_into().unwrap())) != 0 {
                                    anyNullThisTime = 1;
                                    strcpy(
                                        *((pv.Data as *mut *mut c_char)
                                            .add(jj.try_into().unwrap())),
                                        *(pv.Null as *mut *mut c_char),
                                    );
                                } else {
                                    strcpy(
                                        *((pv.Data as *mut *mut c_char)
                                            .add(jj.try_into().unwrap())),
                                        *(result.value.data.strptr.add(jj.try_into().unwrap())),
                                    );
                                }
                            }
                        }
                    } else {
                        ffpmsg_str("Cannot convert string expression to desired type.");
                        lParse.status = PARSE_BAD_TYPE;
                    }
                    if result.operation > 0 {
                        FREE!(*(result.value.data.strptr));
                        FREE!(result.value.data.strptr);
                    }
                }
                _ => {}
            }

            if lParse.status != 0 {
                break;
            }

            /*  Increment Data to point to where the next block should go  */

            if result.ntype == fits_parser_yytokentype::BITSTR as c_int
                && (*(pv.userInfo)).datatype == TBYTE
            {
                (pv.Data) = (pv.Data as *mut c_char).add(
                    ((pv.datasize as c_long) * ((result.value.nelem + 7) / 8) * ntodo)
                        .try_into()
                        .unwrap(),
                ) as *mut c_void;
            } else if result.ntype == fits_parser_yytokentype::STRING as c_int {
                (pv.Data) = (pv.Data as *mut c_char)
                    .add(((pv.datasize as c_long) * ntodo).try_into().unwrap())
                    as *mut c_void;
            } else {
                (pv.Data) = (pv.Data as *mut c_char).add(
                    ((pv.datasize as c_long) * ntodo * (pv.repeat))
                        .try_into()
                        .unwrap(),
                ) as *mut c_void;
            }
        }

        /* If a TemporaryCol output is used, we want to inform the caller
        what the null value is expected to be */
        
        // WARNING - THIS IS DANGEROUS. If nCols = 0, points to invalid memory.
        // In the case of the where expr = "#ROW > 2" there are no columns.
        outcol = colData.as_mut_ptr().offset((nCols - 1) as isize).as_mut().unwrap(); 
        // outcol = &mut colData[(nCols - 1) as usize]; // Re-bind 

        if pv.Null != outcol.array
            && (Data0)
                == (outcol.array as *mut c_char).add((pv.datasize).try_into().unwrap())
                    as *mut c_void
        {
            if (*(pv.userInfo)).datatype == TSTRING {
                memcpy(
                    outcol.array,
                    (*(pv.Null as *mut *mut c_char)) as *mut c_void,
                    2,
                );
            } else {
                memcpy(outcol.array, pv.Null, pv.datasize.try_into().unwrap());
            }
        }

        /* If no NULLs encountered during this pass, set Null value to */
        /* zero to make the writing of the output column data faster   */

        if anyNullThisTime != 0 {
            (*(pv.userInfo)).anyNull = 1;
        } else if pv.Null == outcol.array {
            if (*(pv.userInfo)).datatype == TSTRING {
                memcpy(
                    (*(pv.Null as *mut *mut c_char)) as *mut c_void,
                    zeros.as_mut_ptr() as *mut c_void,
                    2,
                );
            } else {
                memcpy(
                    pv.Null,
                    zeros.as_mut_ptr() as *mut c_void,
                    (pv.datasize).try_into().unwrap(),
                );
            }
        }

        /*-------------------------------------------------------*/
        /*  Clean up procedures:  after processing all the rows  */
        /*-------------------------------------------------------*/

        /*  if the calling routine specified that only a limited number    */
        /*  of rows in the table should be processed, return a value of -1 */
        /*  once all the rows have been done, if no other error occurred.  */

        if lParse.hdutype != IMAGE_HDU
            && firstrow - 1 == (pv.lastRow)
            && lParse.status == 0
            && (*(pv.userInfo)).maxRows < totalrows
        {
            return -1;
        }

        lParse.status /* return successful status */
    }
}

/***********************************************************************/
/// Setup the varData array in gParse to contain the fits column data.
/// Then, allocate and initialize the necessary UNDEF arrays for each
/// column used by the parser.
fn Setup_DataArrays(
    lParse: &mut ParseData,
    nCols: c_int,
    cols: &mut [iteratorCol],
    fRow: c_long,
    nRows: c_long,
) {
    unsafe {
        let i: c_int = 0;
        let mut nelem: c_long = 0;
        let mut len: c_long = 0;
        let mut row: c_long = 0;
        let mut idx: c_long = 0;

        let mut bitStrs: *mut *mut c_char = std::ptr::null_mut();
        let mut sptr: *mut *mut c_char = std::ptr::null_mut();
        let mut barray: *mut c_char = std::ptr::null_mut();
        let mut iarray: *mut c_long = std::ptr::null_mut();
        let mut rarray: *mut c_double = std::ptr::null_mut();
        let mut msg: [c_char; 80] = [0; 80];
        let mut do_realloc: c_int = 0;

        lParse.firstDataRow = fRow;
        lParse.nDataRows = nRows;
        /* Only perform reallocations if the number of rows changed */
        if lParse.nPrevDataRows != nRows {
            do_realloc = 1;
        }

        /*  Resize and fill in UNDEF arrays for each column  */
        for i in 0..nCols as usize {
            let icol: &mut iteratorCol = &mut cols[i];
            let varData = &mut lParse.varData[i];

            if icol.iotype == OutputCol || icol.iotype == TemporaryCol {
                continue;
            }

            nelem = varData.nelem;
            len = nelem * nRows;

            match varData.dtype.into() {
                fits_parser_yytokentype::BITSTR => {
                    /* No need for UNDEF array, but must make string DATA array */
                    len = (nelem + 1) * nRows; /* Count '\0' */
                    bitStrs = varData.data as *mut *mut c_char;
                    if do_realloc != 0 {
                        if !bitStrs.is_null() {
                            FREE!(*bitStrs);
                        }
                        free(bitStrs as *mut c_void);
                        bitStrs =
                            malloc(nRows as usize * size_of::<*mut c_char>()) as *mut *mut c_char;
                        if bitStrs.is_null() {
                            varData.undef = ptr::null_mut();
                            varData.data = ptr::null_mut();
                            lParse.status = MEMORY_ALLOCATION;
                            break;
                        }
                        *bitStrs = malloc(len as usize * size_of::<c_char>()) as *mut c_char;
                        if (*bitStrs).is_null() {
                            free(bitStrs as *mut c_void);
                            varData.undef = ptr::null_mut();
                            varData.data = ptr::null_mut();
                            lParse.status = MEMORY_ALLOCATION;
                            break;
                        }
                    }

                    for row in 0..nRows as usize {
                        *bitStrs.add(row) =
                            (*bitStrs.add(0)).add((row as c_long * (nelem + 1)) as usize);
                        idx = (row as c_long) * ((nelem + 7) / 8) + 1;
                        for len in 0..nelem as usize {
                            if (*(icol.array as *mut c_char).add(idx as usize)
                                & (1 << (7 - (len % 8))))
                                != 0
                            {
                                *(*bitStrs.add(row)).add(len) = b'1' as c_char;
                            } else {
                                *(*bitStrs.add(row)).add(len) = b'0' as c_char;
                            }
                            if len % 8 == 7 {
                                idx += 1;
                            }
                        }
                        *(*bitStrs.add(row)).add(len as usize) = 0;
                    }
                    varData.undef = bitStrs as *mut c_char;
                    varData.data = bitStrs as *mut c_void;
                }

                fits_parser_yytokentype::STRING => {
                    sptr = icol.array as *mut *mut c_char;
                    if do_realloc != 0 {
                        if !varData.undef.is_null() {
                            free(varData.undef as *mut c_void);
                        }
                        varData.undef = malloc(nRows as usize * size_of::<c_char>()) as *mut c_char;
                        if varData.undef.is_null() {
                            lParse.status = MEMORY_ALLOCATION;
                            break;
                        }
                    }
                    row = nRows;
                    while row > 0 {
                        row -= 1;
                        *varData.undef.add(row as usize) = (**sptr != 0
                            && FSTRCMP(
                                cast_slice(CStr::from_ptr(*sptr.add(0)).to_bytes_with_nul()),
                                cast_slice(
                                    CStr::from_ptr(*sptr.add((row + 1) as usize))
                                        .to_bytes_with_nul(),
                                ),
                            ) == 0)
                            as c_char;
                    }
                    varData.data = sptr.add(1) as *mut c_void;
                }

                fits_parser_yytokentype::BOOLEAN => {
                    barray = icol.array as *mut c_char;
                    if do_realloc != 0 {
                        if !varData.undef.is_null() {
                            free(varData.undef as *mut c_void);
                        }
                        varData.undef = malloc(len as usize * size_of::<c_char>()) as *mut c_char;
                        if varData.undef.is_null() {
                            lParse.status = MEMORY_ALLOCATION;
                            break;
                        }
                    }
                    while len > 0 {
                        len -= 1;
                        *varData.undef.add(len as usize) = (*barray.add(0) != 0
                            && *barray.add(0) == *barray.add((len + 1) as usize))
                            as c_char;
                    }
                    varData.data = barray.add(1) as *mut c_void;
                }

                fits_parser_yytokentype::LONG => {
                    iarray = icol.array as *mut c_long;
                    if do_realloc != 0 {
                        if !varData.undef.is_null() {
                            free(varData.undef as *mut c_void);
                        }
                        varData.undef = malloc(len as usize * size_of::<c_char>()) as *mut c_char;
                        if varData.undef.is_null() {
                            lParse.status = MEMORY_ALLOCATION;
                            break;
                        }
                    }
                    while len > 0 {
                        len -= 1;
                        *varData.undef.add(len as usize) = (*iarray.add(0) != 0
                            && *iarray.add(0) == *iarray.add((len + 1) as usize))
                            as c_char;
                    }
                    varData.data = iarray.add(1) as *mut c_void;
                }

                fits_parser_yytokentype::DOUBLE => {
                    rarray = icol.array as *mut c_double;
                    if do_realloc != 0 {
                        if !varData.undef.is_null() {
                            free(varData.undef as *mut c_void);
                        }
                        varData.undef = malloc(len as usize * size_of::<c_char>()) as *mut c_char;
                        if varData.undef.is_null() {
                            lParse.status = MEMORY_ALLOCATION;
                            break;
                        }
                    }
                    while len > 0 {
                        len -= 1;
                        *varData.undef.add(len as usize) = (*rarray.add(0) != 0.0
                            && *rarray.add(0) == *rarray.add((len + 1) as usize))
                            as c_char;
                    }
                    varData.data = rarray.add(1) as *mut c_void;
                }

                _ => {
                    int_snprintf!(
                        &mut msg,
                        80,
                        "SetupDataArrays, unhandled type {}\n",
                        varData.dtype
                    );
                    ffpmsg_slice(&msg);
                }
            }

            if lParse.status != 0 {
                /*  Deallocate NULL arrays of previous columns */
                let mut j = i;
                while j > 0 {
                    j -= 1;
                    let varData = &mut lParse.varData[j];
                    if varData.dtype == fits_parser_yytokentype::BITSTR as c_int {
                        FREE!(*(varData.data as *mut *mut c_char).add(0));
                    }
                    FREE!(varData.undef);
                    varData.undef = ptr::null_mut();
                }
                lParse.nPrevDataRows = 0;
                return;
            }
        }

        lParse.nPrevDataRows = nRows;
    }
}

/*--------------------------------------------------------------------------*/
/// Convert an array of any input data type to an array of any output
/// data type, using an array of UNDEF flags to assign nulvals to
fn ffcvtn(
    inputType: c_int,      /* I - Data type of input array               */
    input: *const c_void,  /* I - Input array of type inputType          */
    undef: *const c_char,  /* I - Array of flags indicating UNDEF elems  */
    ntodo: c_long,         /* I - Number of elements to process          */
    outputType: c_int,     /* I - Data type of output array              */
    nulval: *const c_void, /* I - Ptr to value to use for UNDEF elements */
    output: *mut c_void,   /* O - Output array of type outputType        */
    anynull: &mut c_int,   /* O - Any nulls flagged?                     */
    status: &mut c_int,    /* O - Error status                           */
) -> c_int {
    unsafe {
        let i: c_long = 0;
        let mut dummy_nullarray: [c_char; 0] = [0; 0];

        match outputType {
            TLOGICAL => {
                match inputType {
                    TLOGICAL | TBYTE => {
                        for i in 0..ntodo {
                            if (*((input as *const c_uchar).add(i.try_into().unwrap()))) != 0 {
                                *((output as *mut c_uchar).add(i.try_into().unwrap())) = 1;
                            } else {
                                *((output as *mut c_uchar).add(i.try_into().unwrap())) = 0;
                            }
                        }
                    }
                    TSHORT => {
                        for i in 0..ntodo {
                            if (*((input as *const c_short).add(i.try_into().unwrap()))) != 0 {
                                *((output as *mut c_uchar).add(i.try_into().unwrap())) = 1;
                            } else {
                                *((output as *mut c_uchar).add(i.try_into().unwrap())) = 0;
                            }
                        }
                    }
                    TLONG => {
                        for i in 0..ntodo {
                            if (*((input as *const c_long).add(i.try_into().unwrap()))) != 0 {
                                *((output as *mut c_uchar).add(i.try_into().unwrap())) = 1;
                            } else {
                                *((output as *mut c_uchar).add(i.try_into().unwrap())) = 0;
                            }
                        }
                    }
                    TFLOAT => {
                        for i in 0..ntodo {
                            if (*((input as *const c_float).add(i.try_into().unwrap()))) != 0.0 {
                                *((output as *mut c_uchar).add(i.try_into().unwrap())) = 1;
                            } else {
                                *((output as *mut c_uchar).add(i.try_into().unwrap())) = 0;
                            }
                        }
                    }
                    TDOUBLE => {
                        for i in 0..ntodo {
                            if (*((input as *const c_double).add(i.try_into().unwrap()))) != 0.0 {
                                *((output as *mut c_uchar).add(i.try_into().unwrap())) = 1;
                            } else {
                                *((output as *mut c_uchar).add(i.try_into().unwrap())) = 0;
                            }
                        }
                    }
                    _ => {
                        *status = BAD_DATATYPE;
                    }
                }
                for i in 0..ntodo {
                    if *((undef).add(i.try_into().unwrap())) != 0 {
                        *((output as *mut c_uchar).add(i.try_into().unwrap())) =
                            *(nulval as *const c_uchar);
                        *anynull = 1;
                    }
                }
            }
            TBYTE => {
                match inputType {
                    TLOGICAL | TBYTE => {
                        for i in 0..ntodo {
                            *((output as *mut c_uchar).add(i.try_into().unwrap())) =
                                *((input as *const c_uchar).add(i.try_into().unwrap()));
                        }
                    }
                    TSHORT => {
                        let input_slice = unsafe {
                            std::slice::from_raw_parts(
                                input as *const c_short,
                                ntodo.try_into().unwrap(),
                            )
                        };

                        let output_slice = unsafe {
                            std::slice::from_raw_parts_mut(
                                output as *mut c_uchar,
                                ntodo.try_into().unwrap(),
                            )
                        };

                        fffi2i1(
                            input_slice,
                            ntodo,
                            1.,
                            0.,
                            NullCheckType::None,
                            0,
                            0,
                            &mut dummy_nullarray,
                            None,
                            output_slice,
                            status,
                        );
                    }
                    TLONG => {
                        for i in 0..ntodo {
                            if *((undef).add(i.try_into().unwrap())) != 0 {
                                *((output as *mut c_uchar).add(i.try_into().unwrap())) =
                                    *(nulval as *const c_uchar);
                                *anynull = 1;
                            } else if *((input as *const c_long).add(i.try_into().unwrap())) < 0 {
                                *status = OVERFLOW_ERR;
                                *((output as *mut c_uchar).add(i.try_into().unwrap())) = 0;
                            } else if *((input as *const c_long).add(i.try_into().unwrap()))
                                > UCHAR_MAX
                            {
                                *status = OVERFLOW_ERR;
                                *((output as *mut c_uchar).add(i.try_into().unwrap())) =
                                    UCHAR_MAX as c_uchar;
                            } else {
                                *((output as *mut c_uchar).add(i.try_into().unwrap())) =
                                    *((input as *const c_long).add(i.try_into().unwrap()))
                                        as c_uchar;
                            }
                        }
                        return *status;
                    }
                    TFLOAT => {
                        let input_slice = unsafe {
                            std::slice::from_raw_parts(
                                input as *const f32,
                                ntodo.try_into().unwrap(),
                            )
                        };

                        let output_slice = unsafe {
                            std::slice::from_raw_parts_mut(
                                output as *mut c_uchar,
                                ntodo.try_into().unwrap(),
                            )
                        };

                        fffr4i1(
                            input_slice,
                            ntodo,
                            1.,
                            0.,
                            NullCheckType::None,
                            0,
                            &mut dummy_nullarray,
                            None,
                            output_slice,
                            status,
                        );
                    }
                    TDOUBLE => {
                        let input_slice = unsafe {
                            std::slice::from_raw_parts(
                                input as *const f64,
                                ntodo.try_into().unwrap(),
                            )
                        };

                        let output_slice = unsafe {
                            std::slice::from_raw_parts_mut(
                                output as *mut c_uchar,
                                ntodo.try_into().unwrap(),
                            )
                        };

                        fffr8i1(
                            input_slice,
                            ntodo,
                            1.,
                            0.,
                            NullCheckType::None,
                            0,
                            &mut dummy_nullarray,
                            None,
                            output_slice,
                            status,
                        );
                    }
                    _ => {
                        *status = BAD_DATATYPE;
                    }
                }

                for i in 0..ntodo {
                    if *((undef).add(i.try_into().unwrap())) != 0 {
                        *((output as *mut c_uchar).add(i.try_into().unwrap())) =
                            *(nulval as *const c_uchar);
                        *anynull = 1;
                    }
                }
            }
            TSHORT => {
                match inputType {
                    TLOGICAL | TBYTE => {
                        for i in 0..ntodo {
                            *((output as *mut c_short).add(i.try_into().unwrap())) =
                                *((input as *const c_uchar).add(i.try_into().unwrap())) as c_short;
                        }
                    }
                    TSHORT => {
                        for i in 0..ntodo {
                            *((output as *mut c_short).add(i.try_into().unwrap())) =
                                *((input as *const c_short).add(i.try_into().unwrap()));
                        }
                    }
                    TLONG => {
                        for i in 0..ntodo {
                            if *((undef).add(i.try_into().unwrap())) != 0 {
                                *((output as *mut c_short).add(i.try_into().unwrap())) =
                                    *(nulval as *const c_short);
                                *anynull = 1;
                            } else if *((input as *const c_long).add(i.try_into().unwrap()))
                                < SHRT_MIN
                            {
                                *status = OVERFLOW_ERR;
                                *((output as *mut c_short).add(i.try_into().unwrap())) =
                                    SHRT_MIN as c_short;
                            } else if *((input as *const c_long).add(i.try_into().unwrap()))
                                > SHRT_MAX as c_long
                            {
                                *status = OVERFLOW_ERR;
                                *((output as *mut c_short).add(i.try_into().unwrap())) =
                                    SHRT_MAX as c_short;
                            } else {
                                *((output as *mut c_short).add(i.try_into().unwrap())) =
                                    *((input as *const c_long).add(i.try_into().unwrap()))
                                        as c_short;
                            }
                        }
                        return *status;
                    }
                    TFLOAT => {
                        let input_slice = unsafe {
                            std::slice::from_raw_parts(
                                input as *const f32,
                                ntodo.try_into().unwrap(),
                            )
                        };

                        let output_slice = unsafe {
                            std::slice::from_raw_parts_mut(
                                output as *mut c_short,
                                ntodo.try_into().unwrap(),
                            )
                        };

                        fffr4i2(
                            input_slice,
                            ntodo,
                            1.,
                            0.,
                            NullCheckType::None,
                            0,
                            &mut dummy_nullarray,
                            None,
                            output_slice,
                            status,
                        );
                    }
                    TDOUBLE => {
                        let input_slice = unsafe {
                            std::slice::from_raw_parts(
                                input as *const f64,
                                ntodo.try_into().unwrap(),
                            )
                        };

                        let output_slice = unsafe {
                            std::slice::from_raw_parts_mut(
                                output as *mut c_short,
                                ntodo.try_into().unwrap(),
                            )
                        };

                        fffr8i2(
                            input_slice,
                            ntodo,
                            1.,
                            0.,
                            NullCheckType::None,
                            0,
                            &mut dummy_nullarray,
                            None,
                            output_slice,
                            status,
                        );
                    }
                    _ => {
                        *status = BAD_DATATYPE;
                    }
                }
                for i in 0..ntodo {
                    if *((undef).add(i.try_into().unwrap())) != 0 {
                        *((output as *mut c_short).add(i.try_into().unwrap())) =
                            *(nulval as *const c_short);
                        *anynull = 1;
                    }
                }
            }
            TINT => {
                match inputType {
                    TLOGICAL | TBYTE => {
                        for i in 0..ntodo {
                            *((output as *mut c_int).add(i.try_into().unwrap())) =
                                *((input as *const c_uchar).add(i.try_into().unwrap())) as c_int;
                        }
                    }
                    TSHORT => {
                        for i in 0..ntodo {
                            *((output as *mut c_int).add(i.try_into().unwrap())) =
                                *((input as *const c_short).add(i.try_into().unwrap())) as c_int;
                        }
                    }
                    TLONG => {
                        for i in 0..ntodo {
                            *((output as *mut c_int).add(i.try_into().unwrap())) =
                                *((input as *const c_long).add(i.try_into().unwrap())) as c_int;
                        }
                    }
                    TFLOAT => {
                        let input_slice = unsafe {
                            std::slice::from_raw_parts(
                                input as *const f32,
                                ntodo.try_into().unwrap(),
                            )
                        };

                        let output_slice = unsafe {
                            std::slice::from_raw_parts_mut(
                                output as *mut c_int,
                                ntodo.try_into().unwrap(),
                            )
                        };

                        fffr4int(
                            input_slice,
                            ntodo,
                            1.,
                            0.,
                            NullCheckType::None,
                            0,
                            &mut dummy_nullarray,
                            None,
                            output_slice,
                            status,
                        );
                    }
                    TDOUBLE => {
                        let input_slice = unsafe {
                            std::slice::from_raw_parts(
                                input as *const f64,
                                ntodo.try_into().unwrap(),
                            )
                        };

                        let output_slice = unsafe {
                            std::slice::from_raw_parts_mut(
                                output as *mut c_int,
                                ntodo.try_into().unwrap(),
                            )
                        };

                        fffr8int(
                            input_slice,
                            ntodo,
                            1.,
                            0.,
                            NullCheckType::None,
                            0,
                            &mut dummy_nullarray,
                            None,
                            output_slice,
                            status,
                        );
                    }
                    _ => {
                        *status = BAD_DATATYPE;
                    }
                }
                for i in 0..ntodo {
                    if *((undef).add(i.try_into().unwrap())) != 0 {
                        *((output as *mut c_int).add(i.try_into().unwrap())) =
                            *(nulval as *const c_int);
                        *anynull = 1;
                    }
                }
            }
            TLONG => {
                match inputType {
                    TLOGICAL | TBYTE => {
                        for i in 0..ntodo {
                            *((output as *mut c_long).add(i.try_into().unwrap())) =
                                *((input as *const c_uchar).add(i.try_into().unwrap())) as c_long;
                        }
                    }
                    TSHORT => {
                        for i in 0..ntodo {
                            *((output as *mut c_long).add(i.try_into().unwrap())) =
                                *((input as *const c_short).add(i.try_into().unwrap())) as c_long;
                        }
                    }
                    TLONG => {
                        for i in 0..ntodo {
                            *((output as *mut c_long).add(i.try_into().unwrap())) =
                                *((input as *const c_long).add(i.try_into().unwrap()));
                        }
                    }
                    TFLOAT => {
                        let input_slice = unsafe {
                            std::slice::from_raw_parts(
                                input as *const f32,
                                ntodo.try_into().unwrap(),
                            )
                        };

                        let output_slice = unsafe {
                            std::slice::from_raw_parts_mut(
                                output as *mut c_long,
                                ntodo.try_into().unwrap(),
                            )
                        };

                        fffr4i4(
                            input_slice,
                            ntodo,
                            1.,
                            0.,
                            NullCheckType::None,
                            0,
                            &mut dummy_nullarray,
                            None,
                            output_slice,
                            status,
                        );
                    }
                    TDOUBLE => {
                        let input_slice = unsafe {
                            std::slice::from_raw_parts(
                                input as *const f64,
                                ntodo.try_into().unwrap(),
                            )
                        };

                        let output_slice = unsafe {
                            std::slice::from_raw_parts_mut(
                                output as *mut c_long,
                                ntodo.try_into().unwrap(),
                            )
                        };

                        fffr8i4(
                            input_slice,
                            ntodo,
                            1.,
                            0.,
                            NullCheckType::None,
                            0,
                            &mut dummy_nullarray,
                            None,
                            output_slice,
                            status,
                        );
                    }
                    _ => {
                        *status = BAD_DATATYPE;
                    }
                }
                for i in 0..ntodo {
                    if *((undef).add(i.try_into().unwrap())) != 0 {
                        *((output as *mut c_long).add(i.try_into().unwrap())) =
                            *(nulval as *const c_long);
                        *anynull = 1;
                    }
                }
            }
            TLONGLONG => {
                match inputType {
                    TLOGICAL | TBYTE => {
                        for i in 0..ntodo {
                            *((output as *mut LONGLONG).add(i.try_into().unwrap())) =
                                *((input as *const c_uchar).add(i.try_into().unwrap())) as LONGLONG;
                        }
                    }
                    TSHORT => {
                        for i in 0..ntodo {
                            *((output as *mut LONGLONG).add(i.try_into().unwrap())) =
                                *((input as *const c_short).add(i.try_into().unwrap())) as LONGLONG;
                        }
                    }
                    TLONG => {
                        for i in 0..ntodo {
                            *((output as *mut LONGLONG).add(i.try_into().unwrap())) =
                                *((input as *const c_long).add(i.try_into().unwrap()));
                        }
                    }
                    TFLOAT => {
                        let input_slice = unsafe {
                            std::slice::from_raw_parts(
                                input as *const f32,
                                ntodo.try_into().unwrap(),
                            )
                        };

                        let output_slice = unsafe {
                            std::slice::from_raw_parts_mut(
                                output as *mut LONGLONG,
                                ntodo.try_into().unwrap(),
                            )
                        };

                        fffr4i8(
                            input_slice,
                            ntodo,
                            1.,
                            0.,
                            NullCheckType::None,
                            0,
                            &mut dummy_nullarray,
                            None,
                            output_slice,
                            status,
                        );
                    }
                    TDOUBLE => {
                        let input_slice = unsafe {
                            std::slice::from_raw_parts(
                                input as *const f64,
                                ntodo.try_into().unwrap(),
                            )
                        };

                        let output_slice = unsafe {
                            std::slice::from_raw_parts_mut(
                                output as *mut LONGLONG,
                                ntodo.try_into().unwrap(),
                            )
                        };

                        fffr8i8(
                            input_slice,
                            ntodo,
                            1.,
                            0.,
                            NullCheckType::None,
                            0,
                            &mut dummy_nullarray,
                            None,
                            output_slice,
                            status,
                        );
                    }
                    _ => {
                        *status = BAD_DATATYPE;
                    }
                }
                for i in 0..ntodo {
                    if *((undef).add(i.try_into().unwrap())) != 0 {
                        *((output as *mut LONGLONG).add(i.try_into().unwrap())) =
                            *(nulval as *const LONGLONG);
                        *anynull = 1;
                    }
                }
            }
            TFLOAT => {
                match inputType {
                    TLOGICAL | TBYTE => {
                        for i in 0..ntodo {
                            *((output as *mut c_float).add(i.try_into().unwrap())) =
                                *((input as *const c_uchar).add(i.try_into().unwrap())) as c_float;
                        }
                    }
                    TSHORT => {
                        for i in 0..ntodo {
                            *((output as *mut c_float).add(i.try_into().unwrap())) =
                                *((input as *const c_short).add(i.try_into().unwrap())) as c_float;
                        }
                    }
                    TLONG => {
                        for i in 0..ntodo {
                            *((output as *mut c_float).add(i.try_into().unwrap())) =
                                *((input as *const c_long).add(i.try_into().unwrap())) as c_float;
                        }
                    }
                    TFLOAT => {
                        for i in 0..ntodo {
                            *((output as *mut c_float).add(i.try_into().unwrap())) =
                                *((input as *const c_float).add(i.try_into().unwrap()));
                        }
                    }
                    TDOUBLE => {
                        let input_slice = unsafe {
                            std::slice::from_raw_parts(
                                input as *const f64,
                                ntodo.try_into().unwrap(),
                            )
                        };

                        let output_slice = unsafe {
                            std::slice::from_raw_parts_mut(
                                output as *mut f32,
                                ntodo.try_into().unwrap(),
                            )
                        };

                        fffr8r4(
                            input_slice,
                            ntodo,
                            1.,
                            0.,
                            NullCheckType::None,
                            0.0,
                            &mut dummy_nullarray,
                            None,
                            output_slice,
                            status,
                        );
                    }
                    _ => {
                        *status = BAD_DATATYPE;
                    }
                }
                for i in 0..ntodo {
                    if *((undef).add(i.try_into().unwrap())) != 0 {
                        *((output as *mut c_float).add(i.try_into().unwrap())) =
                            *(nulval as *const c_float);
                        *anynull = 1;
                    }
                }
            }
            TDOUBLE => {
                match inputType {
                    TLOGICAL | TBYTE => {
                        for i in 0..ntodo {
                            *((output as *mut c_double).add(i.try_into().unwrap())) =
                                *((input as *const c_uchar).add(i.try_into().unwrap())) as c_double;
                        }
                    }
                    TSHORT => {
                        for i in 0..ntodo {
                            *((output as *mut c_double).add(i.try_into().unwrap())) =
                                *((input as *const c_short).add(i.try_into().unwrap())) as c_double;
                        }
                    }
                    TLONG => {
                        for i in 0..ntodo {
                            *((output as *mut c_double).add(i.try_into().unwrap())) =
                                *((input as *const c_long).add(i.try_into().unwrap())) as c_double;
                        }
                    }
                    TFLOAT => {
                        for i in 0..ntodo {
                            *((output as *mut c_double).add(i.try_into().unwrap())) =
                                *((input as *const c_float).add(i.try_into().unwrap())) as c_double;
                        }
                    }
                    TDOUBLE => {
                        for i in 0..ntodo {
                            *((output as *mut c_double).add(i.try_into().unwrap())) =
                                *((input as *const c_double).add(i.try_into().unwrap()))
                                    as c_double;
                        }
                    }
                    _ => {
                        *status = BAD_DATATYPE;
                    }
                }
                for i in 0..ntodo {
                    if *((undef).add(i.try_into().unwrap())) != 0 {
                        *((output as *mut c_double).add(i.try_into().unwrap())) =
                            *(nulval as *const c_double);
                        *anynull = 1;
                    }
                }
            }
            _ => {
                *status = BAD_DATATYPE;
            }
        }

        *status
    }
}

/*---------------------------------------------------------------------------*/
/// Evaluate a boolean expression for each time in a compressed file,
/// returning an array of flags indicating which times evaluated to TRUE/FALSE
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn fffrwc(
    fptr: *mut fitsfile,      /* I - Input FITS file                    */
    expr: *const c_char,      /* I - Boolean expression                 */
    timeCol: *const c_char,   /* I - Name of time column                */
    parCol: *const c_char,    /* I - Name of parameter column           */
    valCol: *const c_char,    /* I - Name of value column               */
    ntimes: c_long,           /* I - Number of distinct times in file   */
    times: *mut f64,          /* O - Array of times in file             */
    time_status: *mut c_char, /* O - Array of boolean results           */
    status: *mut c_int,       /* O - Error status                       */
) -> c_int {
    unsafe {
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);

        let times = std::slice::from_raw_parts_mut(times, ntimes as usize);
        let time_status = std::slice::from_raw_parts_mut(time_status, ntimes as usize);

        raw_to_slice!(expr);
        raw_to_slice!(timeCol);
        raw_to_slice!(parCol);
        raw_to_slice!(valCol);

        fffrwc_safe(
            fptr,
            expr,
            timeCol,
            parCol,
            valCol,
            ntimes,
            times,
            time_status,
            status,
        )
    }
}

/*---------------------------------------------------------------------------*/
/// Evaluate a boolean expression for each time in a compressed file,
/// returning an array of flags indicating which times evaluated to TRUE/FALSE
pub fn fffrwc_safe(
    fptr: &mut fitsfile,        /* I - Input FITS file                    */
    expr: &[c_char],            /* I - Boolean expression                 */
    timeCol: &[c_char],         /* I - Name of time column                */
    parCol: &[c_char],          /* I - Name of parameter column           */
    valCol: &[c_char],          /* I - Name of value column               */
    ntimes: c_long,             /* I - Number of distinct times in file   */
    times: &mut [f64],          /* O - Array of times in file             */
    time_status: &mut [c_char], /* O - Array of boolean results           */
    status: &mut c_int,         /* O - Error status                       */
) -> c_int {
    unsafe {
        let mut Info: parseInfo = Default::default();
        let mut alen: c_long = 0;
        let mut width: c_long = 0;
        let mut parNo: c_int = 0;
        let mut typecode: c_int = 0;
        let mut naxis: c_int = 0;
        let mut constant: c_int = 0;
        let mut nCol: c_int = 0;
        let mut nelem: c_long = 0;
        let mut naxes: [c_long; MAXDIMS as usize] = [0; MAXDIMS as usize];
        let mut result: c_char = 0;
        let mut lParse: ParseData = ParseData::default();

        if *status != 0 {
            return *status;
        }

        if ffiprs(
            fptr,
            1,
            expr,
            MAXDIMS,
            &mut Info.datatype,
            &mut nelem,
            &mut naxis,
            &mut naxes,
            &mut lParse,
            status,
        ) != 0
        {
            ffcprs(&mut lParse);
            return *status;
        }

        fits_get_colnum(
            fptr,
            CASEINSEN.try_into().unwrap(),
            timeCol,
            &mut lParse.timeCol,
            status,
        );
        fits_get_colnum(
            fptr,
            CASEINSEN.try_into().unwrap(),
            parCol,
            &mut lParse.parCol,
            status,
        );
        fits_get_colnum(
            fptr,
            CASEINSEN.try_into().unwrap(),
            valCol,
            &mut lParse.valCol,
            status,
        );
        if *status != 0 {
            return *status;
        }

        if nelem < 0 {
            constant = 1;
            nelem = -nelem;
            nCol = lParse.nCols;
            lParse.nCols = 0; /*  Ignore all column references  */
        } else {
            constant = 0;
        }

        if Info.datatype != TLOGICAL || nelem != 1 {
            ffcprs(&mut lParse);
            ffpmsg_str("Expression does not evaluate to a logical scalar.");
            *status = PARSE_BAD_TYPE;
            return *status;
        }

        /*******************************************/
        /* Allocate data arrays for each parameter */
        /*******************************************/

        parNo = lParse.nCols;
        while parNo > 0 {
            parNo -= 1;
            let mut col_data = lParse.colData[parNo as usize];
            match col_data.datatype {
                TLONG => {
                    col_data.array =
                        malloc(((ntimes + 1) as usize) * size_of::<c_long>()) as *mut c_void;
                    if !col_data.array.is_null() {
                        let array_ptr = col_data.array as *mut c_long;
                        *array_ptr.wrapping_add(0) = 1234554321;
                    } else {
                        *status = MEMORY_ALLOCATION;
                    }
                }
                TDOUBLE => {
                    col_data.array =
                        malloc(((ntimes + 1) as usize) * size_of::<c_double>()) as *mut c_void;
                    if !col_data.array.is_null() {
                        let array_ptr = col_data.array as *mut c_double;
                        *array_ptr.wrapping_add(0) = DOUBLENULLVALUE;
                    } else {
                        *status = MEMORY_ALLOCATION;
                    }
                }
                TSTRING => {
                    if fits_get_coltype(
                        fptr,
                        lParse.valCol,
                        Some(&mut typecode),
                        Some(&mut alen),
                        Some(&mut width),
                        status,
                    ) == 0
                    {
                        alen += 1;
                        col_data.array = malloc(((ntimes + 1) as usize) * size_of::<*mut c_char>())
                            as *mut c_void;
                        if !col_data.array.is_null() {
                            let str_array = col_data.array as *mut *mut c_char;
                            let first_str = malloc(((ntimes + 1) * alen) as usize) as *mut c_char;
                            *str_array.wrapping_add(0) = first_str;
                            if !first_str.is_null() {
                                for elem in 1..=(ntimes as usize) {
                                    *str_array.wrapping_add(elem) = (*str_array
                                        .wrapping_add(elem - 1))
                                    .wrapping_add(alen as usize);
                                }
                                *(*str_array.wrapping_add(0)).wrapping_add(0) = 0;
                            } else {
                                free(col_data.array);
                                *status = MEMORY_ALLOCATION;
                            }
                        } else {
                            *status = MEMORY_ALLOCATION;
                        }
                    }
                }
                _ => {}
            }

            if *status != 0 {
                while parNo > 0 {
                    parNo -= 1;
                    let col_data2 = lParse.colData[parNo as usize];
                    if (col_data2).datatype == TSTRING {
                        let str_array = (col_data2).array as *mut *mut c_char;
                        if !str_array.is_null() {
                            let mut first_str = *str_array.wrapping_add(0);
                            FREE!(first_str);
                        }
                    }
                    let mut array_ptr = (col_data2).array;
                    FREE!(array_ptr);
                }
                return *status;
            }
        }

        /**********************************************************************/
        /* Read data from columns needed for the expression and then parse it */
        /**********************************************************************/

        if fits_uncompress_hkdata(&mut lParse, fptr, ntimes, times, status) == 0 {
            if constant != 0 {
                let result_node = lParse.Nodes[lParse.resultNode as usize];
                result = (result_node).value.data.log;
                let mut elem = ntimes;
                while elem > 0 {
                    elem -= 1;
                    time_status[elem as usize] = result;
                }
            } else {
                Info.dataPtr = time_status.as_mut_ptr() as *mut c_void;
                Info.nullPtr = ptr::null_mut();
                Info.maxRows = ntimes;
                *status = fits_parser_workfn_safe(
                    ntimes,
                    0,
                    1,
                    ntimes,
                    lParse.nCols,
                    &mut lParse.colData,
                    (&mut Info) as *mut parseInfo as *mut c_void,
                );
            }
        }

        /************/
        /* Clean up */
        /************/

        parNo = lParse.nCols;
        while parNo > 0 {
            parNo -= 1;
            let col_data = lParse.colData[parNo as usize];
            if (col_data).datatype == TSTRING {
                let str_array = (col_data).array as *mut *mut c_char;
                if !str_array.is_null() {
                    let mut first_str = *str_array.wrapping_add(0);
                    FREE!(first_str);
                }
            }
            let mut array_ptr = (col_data).array;
            FREE!(array_ptr);
        }

        if constant != 0 {
            lParse.nCols = nCol;
        }

        ffcprs(&mut lParse);
        *status
    }
}

/*---------------------------------------------------------------------------*/
/// Evaluate a boolean expression, returning the row number of the first
/// row which evaluates to TRUE
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffffrw(
    fptr: *mut fitsfile, /* I - Input FITS file                   */
    expr: *mut c_char,   /* I - Boolean expression                */
    rownum: *mut c_long, /* O - First row of table to eval to T   */
    status: *mut c_int,  /* O - Error status                      */
) -> c_int {
    unsafe {
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);
        let rownum = rownum.as_mut().expect(NULL_MSG);

        raw_to_slice!(expr);

        ffffrw_safer(fptr, expr, rownum, status)
    }
}

/*---------------------------------------------------------------------------*/
/// Evaluate a boolean expression, returning the row number of the first
/// row which evaluates to TRUE
pub fn ffffrw_safer(
    fptr: &mut fitsfile, /* I - Input FITS file                   */
    expr: &[c_char],     /* I - Boolean expression                */
    rownum: &mut c_long, /* O - First row of table to eval to T   */
    status: &mut c_int,  /* O - Error status                      */
) -> c_int {
    let mut naxis = 0; /* No need to call parser... have result from ffiprs */
    let mut constant = 0; /*  Make sure there is at least 1 row in table  */
    let mut dtype = 0; /* -1 indicates exitted without error before end... OK */
    let mut nelem: c_long = 0;
    let mut naxes: [c_long; 5] = [0; 5];
    let mut result: c_char = 0;

    let mut lParse: ParseData = ParseData::default();

    if *status != 0 {
        return *status;
    }
    if ffiprs(
        fptr,
        0,
        expr,
        MAXDIMS,
        &mut dtype,
        &mut nelem,
        &mut naxis,
        naxes.as_mut_slice(),
        &mut lParse,
        status,
    ) != 0
    {
        ffcprs(&mut lParse);
        return *status;
    }
    if nelem < 0 {
        constant = 1;
        nelem = -nelem;
    } else {
        constant = 0;
    }

    if dtype != TLOGICAL || nelem != 1 {
        ffcprs(&mut lParse);
        ffpmsg_str("Expression does not evaluate to a logical scalar.");
        return {
            *status = PARSE_BAD_TYPE;
            *status
        };
    }
    *rownum = 0;
    {
        if constant != 0 {
            result = unsafe { (lParse.Nodes[lParse.resultNode as usize]).value.data.log };
            if result != 0 {
                ffgnrw_safe(fptr, &mut nelem, status);
                if nelem != 0 {
                    *rownum = 1;
                };
            };
        } else {
            let mut workData = ffffrw_workdata {
                prownum: rownum,
                lParse: &mut lParse as *mut ParseData,
            };
            let colData_slice = &mut lParse.colData[..];
            if ffiter_safe(
                (lParse).nCols,
                colData_slice,
                0,
                0,
                ffffrw_work,
                (&mut workData as *mut ffffrw_workdata) as *mut c_void,
                status,
            ) == -1
            {
                *status = 0;
            };
        }
    }
    ffcprs(&mut lParse);
    *status
}

/*---------------------------------------------------------------------------*/
fn fits_parser_set_temporary_col(
    lParse: &mut ParseData,
    Info: &mut parseInfo,
    nrows: c_long,
    nulval: *mut c_void,
    status: &mut c_int,
) -> c_int {
    unsafe {
        /* Setup iterator column and parser information to be ready to compute
        temporary calculator expression */

        if *status != 0 {
            return *status;
        }

        let col_cnt: c_int = lParse.nCols;

        if fits_parser_allocateCol(lParse, col_cnt, unsafe { &mut *status }) != 0 {
            return *status;
        }

        /* Set important variables for TemporaryCol where calculated results end up */

        fits_iter_set_by_num_safe(
            &mut lParse.colData[col_cnt as usize],
            ptr::null_mut(),
            0,
            TDOUBLE,
            TemporaryCol,
        );

        (lParse.colData[col_cnt as usize]).repeat = lParse.nElements;

        Info.dataPtr = ptr::null_mut();
        Info.nullPtr = nulval;
        Info.maxRows = nrows;
        Info.parseData = lParse;
        lParse.nCols += 1;

        0
    }
}

/*---------------------------------------------------------------------------*/
/// Uncompress housekeeping data from a compressed FITS table
fn fits_uncompress_hkdata(
    lParse: &mut ParseData,
    fptr: &mut fitsfile,
    ntimes: c_long,
    times: &mut [f64],
    status: &mut c_int,
) -> c_int {
    unsafe {
        let mut parName: [c_char; 256] = [0; 256];
        let sPtr: [*mut c_char; 1] = [parName.as_mut_ptr()];
        let mut found: [c_char; 1000] = [0; 1000];
        let mut parNo: c_int;
        let mut anynul: c_int = 0;
        let mut naxis2: c_long = 0;
        let mut row: c_long;
        let mut currelem: c_long;
        let mut currtime: f64;
        let mut newtime: f64 = 0.0;

        currelem = 0;
        currtime = -1e38;

        parNo = lParse.nCols;
        while parNo > 0 {
            parNo -= 1;
            found[parNo as usize] = 0;
        }

        if ffgkyj_safe(fptr, cs!(c"NAXIS2"), &mut naxis2, None, status) != 0 {
            return *status;
        }

        for row in 1..=naxis2 {
            if ffgcvd_safe(
                fptr,
                lParse.timeCol,
                row,
                1,
                1,
                0.0,
                std::slice::from_mut(&mut newtime),
                Some(&mut anynul),
                status,
            ) != 0
            {
                return *status;
            }

            if newtime != currtime {
                /*  New time encountered... propogate parameters to next row  */
                if currelem == ntimes {
                    ffpmsg_str("Found more unique time stamps than caller indicated");
                    *status = PARSE_BAD_COL;
                    return *status;
                }

                times[currelem as usize] = newtime;
                currtime = newtime;
                currelem += 1;

                parNo = lParse.nCols;
                while parNo > 0 {
                    parNo -= 1;
                    let col_data = lParse.colData[parNo as usize];
                    match col_data.datatype {
                        TLONG => {
                            let array = col_data.array as *mut c_long;
                            *array.wrapping_add(currelem as usize) =
                                *array.wrapping_add((currelem - 1) as usize);
                        }
                        TDOUBLE => {
                            let array = col_data.array as *mut f64;
                            *array.wrapping_add(currelem as usize) =
                                *array.wrapping_add((currelem - 1) as usize);
                        }
                        TSTRING => {
                            let str_array = col_data.array as *mut *mut c_char;
                            strcpy(
                                *str_array.wrapping_add(currelem as usize),
                                *str_array.wrapping_add((currelem - 1) as usize),
                            );
                        }
                        _ => {}
                    }
                }
            }

            if ffgcvs_safe(
                fptr,
                lParse.parCol,
                row,
                1,
                1,
                Some(cs!(c"")),
                &mut [std::slice::from_raw_parts_mut(sPtr[0], 256)],
                Some(&mut anynul),
                status,
            ) != 0
            {
                return *status;
            }

            parNo = lParse.nCols;
            while parNo > 0 {
                parNo -= 1;
                let var_data = &lParse.varData[parNo as usize];
                if fits_strcasecmp(
                    std::slice::from_raw_parts(parName.as_ptr(), strlen_safe(&parName)),
                    std::slice::from_raw_parts(
                        (var_data).name.as_ptr(),
                        strlen_safe(&(var_data).name),
                    ),
                ) == 0
                {
                    break;
                }
            }

            if parNo >= 0 {
                found[parNo as usize] = 1; /* Flag this parameter as found */
                let col_data = lParse.colData[parNo as usize];
                match col_data.datatype {
                    TLONG => {
                        let array = col_data.array as *mut c_long;
                        ffgcvj_safe(
                            fptr,
                            lParse.valCol,
                            row,
                            1,
                            1,
                            *array.wrapping_add(0),
                            std::slice::from_raw_parts_mut(
                                array.wrapping_add(currelem as usize),
                                1,
                            ),
                            Some(&mut anynul),
                            status,
                        );
                    }
                    TDOUBLE => {
                        let array = col_data.array as *mut f64;
                        ffgcvd_safe(
                            fptr,
                            lParse.valCol,
                            row,
                            1,
                            1,
                            *array.wrapping_add(0),
                            std::slice::from_raw_parts_mut(
                                array.wrapping_add(currelem as usize),
                                1,
                            ),
                            Some(&mut anynul),
                            status,
                        );
                    }
                    TSTRING => {
                        let str_array = col_data.array as *mut *mut c_char;
                        ffgcvs_safe(
                            fptr,
                            lParse.valCol,
                            row,
                            1,
                            1,
                            Some(std::slice::from_raw_parts(
                                *str_array.wrapping_add(0),
                                strlen(*str_array.wrapping_add(0)) as usize,
                            )),
                            &mut [std::slice::from_raw_parts_mut(
                                *str_array.wrapping_add(currelem as usize),
                                256,
                            )],
                            Some(&mut anynul),
                            status,
                        );
                    }
                    _ => {}
                }
                if *status != 0 {
                    return *status;
                }
            }
        }

        if currelem < ntimes {
            ffpmsg_str("Found fewer unique time stamps than caller indicated");
            *status = PARSE_BAD_COL;
            return *status;
        }

        /*  Check for any parameters which were not located in the table  */
        parNo = lParse.nCols;
        while parNo > 0 {
            parNo -= 1;
            if found[parNo as usize] == 0 {
                let var_data = &lParse.varData[parNo as usize];
                int_snprintf!(
                    &mut parName,
                    256,
                    "Parameter not found: {:<30}",
                    str::from_utf8(cast_slice(&(var_data).name)).unwrap(),
                );
                ffpmsg_slice(&parName);
                *status = PARSE_SYNTAX_ERR;
            }
        }

        *status
    }
}

/*---------------------------------------------------------------------------*/
/// Iterator work function which calls the parser and searches for the
/// first row which evaluates to TRUE.
pub(crate) extern "C" fn ffffrw_work(
    totalrows: c_long,         /* I - Total rows to be processed     */
    offset: c_long,            /* I - Number of rows skipped at start*/
    firstrow: c_long,          /* I - First row of this iteration    */
    nrows: c_long,             /* I - Number of rows in this iter    */
    nCols: c_int,              /* I - Number of columns in use       */
    colData: *mut iteratorCol, /* IO- Column information/data        */
    userPtr: *mut c_void,
) -> c_int {
    ffffrw_work_safe(totalrows, offset, firstrow, nrows, nCols, colData, userPtr)
}

/*---------------------------------------------------------------------------*/
/// Iterator work function which calls the parser and searches for the
/// first row which evaluates to TRUE.
pub(crate) fn ffffrw_work_safe(
    totalrows: c_long,         /* I - Total rows to be processed     */
    offset: c_long,            /* I - Number of rows skipped at start*/
    firstrow: c_long,          /* I - First row of this iteration    */
    nrows: c_long,             /* I - Number of rows in this iter    */
    nCols: c_int,              /* I - Number of columns in use       */
    colData: *mut iteratorCol, /* IO- Column information/data        */
    userPtr: *mut c_void,
) -> c_int {
    unsafe {
        let result: &mut Node;
        let workData: *mut ffffrw_workdata = userPtr as *mut ffffrw_workdata;
        let lParse: &mut ParseData = &mut *(*workData).lParse;

        Evaluate_Parser(lParse, firstrow, nrows);

        if (lParse.status) == 0 {
            result = &mut lParse.Nodes[lParse.resultNode as usize];
            if result.operation == CONST_OP {
                if result.value.data.log != 0 {
                    *((*workData).prownum) = firstrow;
                    return -1;
                }
            } else {
                for idx in 0..(nrows as usize) {
                    if *(result.value.data.logptr.add(idx)) != 0
                        && *(result.value.undef.add(idx)) == 0
                    {
                        *((*workData).prownum) = firstrow + idx as c_long;
                        return -1;
                    }
                }
            }
        }

        lParse.status
    }
}

/*--------------------------------------------------------------------------*/
/// Apply pixel filtering operations
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn fits_pixel_filter(
    filter: *mut PixelFilter, /* I - pixel filter structure */
    status: *mut c_int,       /* IO - error status */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect("Null status pointer");
        let filter = filter.as_mut().expect("Null filter pointer");

        fits_pixel_filter_safer(filter, status)
    }
}

/// Apply pixel filtering operations (safe version)
pub fn fits_pixel_filter_safer(
    filter: &mut PixelFilter, /* I - pixel filter structure */
    status: &mut c_int,       /* IO - error status */
) -> c_int {
    unsafe {
        let mut info: parseInfo = parseInfo::default();
        let mut naxis: c_int = 0;
        let mut bitpix: c_int = 0;
        let mut nelem: c_long = 0;
        let mut naxes: [c_long; MAXDIMS as usize] = [0; MAXDIMS as usize];
        let mut col_cnt: c_int = 0;
        let result: *mut Node = std::ptr::null_mut();
        let mut datatype: c_int = 0;

        let default_tags: [c_char; 2] = [b'X' as c_char, 0];
        let mut msg: [c_char; 256] = [0; 256];
        let mut write_blank_kwd: c_int = 0; /* write BLANK if any output nulls? */
        let mut lParse: ParseData = ParseData::default();

        let debug_pixfilter = if std::env::var("DEBUG_PIXFILTER").is_ok() {
            1
        } else {
            0
        };

        if *status != 0 {
            return *status;
        }

        unsafe {
            if filter.tag.is_null() || (*filter.tag).is_null() || **filter.tag == 0 {
                filter.tag = default_tags.as_ptr() as *mut *mut c_char;
                if debug_pixfilter != 0 {
                    println!("using default tag '{}'", **filter.tag);
                }
            }
        }

        let infptr: *mut fitsfile = unsafe { *filter.ifptr };
        let outfptr: *mut fitsfile = filter.ofptr;
        lParse.pixFilter = filter;

        let filter_expr = cast_slice(CStr::from_ptr(filter.expression).to_bytes_with_nul());
        if ffiprs(
            unsafe { infptr.as_mut().unwrap() },
            0,
            filter_expr,
            MAXDIMS,
            &mut info.datatype,
            &mut nelem,
            &mut naxis,
            &mut naxes,
            &mut lParse,
            status,
        ) != 0
        {
            ffcprs(&mut lParse);
            return *status;
        }

        if nelem < 0 {
            nelem = -nelem;
        }

        {
            /* validate result type */
            let rtype: &str;
            match info.datatype {
                TLOGICAL => {
                    rtype = "LOGICAL";
                }
                TLONG => {
                    rtype = "LONG";
                }
                TDOUBLE => {
                    rtype = "DOUBLE";
                }
                TSTRING => {
                    rtype = "STRING";
                    *status = P_ERROR;
                    ffpmsg_str("pixel_filter: cannot have string image");
                }
                TBIT => {
                    rtype = "BIT";
                    if debug_pixfilter != 0 {
                        println!("hmm, image from bits?");
                    }
                }
                _ => {
                    rtype = "UNKNOWN?!";
                    *status = P_ERROR;
                    ffpmsg_str("pixel_filter: unexpected result datatype");
                }
            }
            if debug_pixfilter != 0 {
                println!("result type is {} [{}]", rtype, info.datatype);
            }
            if *status != 0 {
                ffcprs(&mut lParse);
                return *status;
            }
        }

        if fits_get_img_param(
            unsafe { infptr.as_mut().unwrap() },
            MAXDIMS,
            Some(&mut bitpix),
            Some(&mut naxis),
            Some(&mut naxes),
            status,
        ) != 0
        {
            ffpmsg_str("pixel_filter: unable to read input image parameters");
            ffcprs(&mut lParse);
            return *status;
        }

        if debug_pixfilter != 0 {
            println!("input bitpix {}", bitpix);
        }

        if info.datatype == TDOUBLE {
            // for floating point expressions, set the default output image to
            // bitpix = -32 (float) unless the default is already a double
            if bitpix != DOUBLE_IMG {
                bitpix = FLOAT_IMG;
            }
        }

        // override output image bitpix if specified by caller
        if filter.bitpix != 0 {
            bitpix = filter.bitpix;
        }
        if debug_pixfilter != 0 {
            println!("output bitpix {}", bitpix);
        }

        if unsafe {
            fits_create_img(
                outfptr.as_mut().unwrap(),
                bitpix,
                naxis,
                &naxes[..naxis as usize],
                status,
            )
        } != 0
        {
            ffpmsg_str("pixel_filter: unable to create output image");
            ffcprs(&mut lParse);
            return *status;
        }

        // transfer keycards
        {
            let mut ncards: c_int = 0;
            let mut more: c_int = 0;
            if fits_get_hdrspace(
                unsafe { infptr.as_mut().unwrap() },
                Some(&mut ncards),
                Some(&mut more),
                status,
            ) != 0
            {
                ffpmsg_str("pixel_filter: unable to determine number of keycards");
                ffcprs(&mut lParse);
                return *status;
            }

            for i in 1..=ncards {
                let mut keyclass: c_int = 0;
                let mut card: [c_char; FLEN_CARD] = [0; FLEN_CARD];

                if fits_read_record(
                    unsafe { infptr.as_mut().unwrap() },
                    i,
                    Some(&mut card),
                    status,
                ) != 0
                {
                    unsafe {
                        int_snprintf!(&mut msg, 256, "pixel_filter: unable to read keycard {}", i,);
                    }
                    ffpmsg_slice(&msg);
                    ffcprs(&mut lParse);
                    return *status;
                }

                keyclass = fits_get_keyclass(&card);
                if keyclass == TYP_STRUC_KEY {
                    // output structure defined by fits_create_img
                } else if keyclass == TYP_COMM_KEY && i < 12 {
                    // assume this is one of the FITS standard comments
                } else if keyclass == TYP_NULL_KEY && bitpix < 0 {
                    // do not transfer BLANK to real output image
                } else if keyclass == TYP_SCAL_KEY && bitpix < 0 {
                    // do not transfer BZERO, BSCALE to real output image
                } else if fits_write_record(unsafe { outfptr.as_mut().unwrap() }, &card, status)
                    != 0
                {
                    unsafe {
                        int_snprintf!(
                            &mut msg,
                            256,
                            "pixel_filter: unable to write keycard [{}]",
                            *status,
                        );
                    }
                    ffpmsg_slice(&msg);
                    ffcprs(&mut lParse);
                    return *status;
                }
            }
        }

        match bitpix {
            BYTE_IMG => {
                datatype = TLONG;
                info.datatype = TBYTE;
            }
            SHORT_IMG => {
                datatype = TLONG;
                info.datatype = TSHORT;
            }
            LONG_IMG => {
                datatype = TLONG;
                info.datatype = TLONG;
            }
            FLOAT_IMG => {
                datatype = TDOUBLE;
                info.datatype = TFLOAT;
            }
            DOUBLE_IMG => {
                datatype = TDOUBLE;
                info.datatype = TDOUBLE;
            }

            _ => {
                unsafe {
                    int_snprintf!(
                        &mut msg,
                        256,
                        "pixel_filter: unexpected output bitpix {}",
                        bitpix,
                    );
                }
                ffpmsg_slice(&msg);
                *status = P_ERROR;
                ffcprs(&mut lParse);
                return *status;
            }
        }

        if bitpix > 0 {
            // arrange for NULLs in output
            let mut null_val: c_long = filter.blank;
            if filter.blank == 0 {
                let mut tstatus: c_int = 0;
                if fits_read_key_lng(
                    unsafe { infptr.as_mut().unwrap() },
                    unsafe { cs!(c"BLANK") },
                    &mut null_val,
                    None,
                    &mut tstatus,
                ) != 0
                {
                    write_blank_kwd = 1;

                    if bitpix == BYTE_IMG {
                        null_val = UCHAR_MAX as c_long;
                    } else if bitpix == SHORT_IMG {
                        null_val = SHRT_MIN as c_long;
                    } else if bitpix == LONG_IMG {
                        if std::mem::size_of::<c_long>() == 8 && std::mem::size_of::<c_int>() == 4 {
                            null_val = INT_MIN as c_long;
                        } else {
                            null_val = LONG_MIN;
                        }
                    } else {
                        println!("unhandled positive output BITPIX {}", bitpix);
                    }
                }

                filter.blank = null_val;
            }

            fits_set_imgnull(unsafe { outfptr.as_mut().unwrap() }, filter.blank, status);
            if debug_pixfilter != 0 {
                println!("using blank {}", null_val);
            }
        }

        if unsafe { *filter.keyword.as_ptr() } == 0 {
            /*************************************/
            /* Create new iterator Output Column */
            /*************************************/
            col_cnt = lParse.nCols;
            if fits_parser_allocateCol(&mut lParse, col_cnt, status) != 0 {
                ffcprs(&mut lParse);
                return *status;
            }
            lParse.nCols += 1;

            let colIter: &mut iteratorCol = &mut lParse.colData[col_cnt as usize];
            unsafe { colIter.fptr = filter.ofptr };
            unsafe { colIter.iotype = OutputCol };

            set_image_col_types(
                &mut lParse.status,
                unsafe { colIter.fptr.as_mut().unwrap() },
                cs!(c"CREATED"),
                bitpix,
                &mut lParse.varData[col_cnt as usize],
                colIter,
            );

            info.maxRows = -1;
            info.parseData = &mut lParse;

            if unsafe {
                let colData_slice = &mut lParse.colData[..];
                ffiter_safe(
                    lParse.nCols,
                    colData_slice,
                    0,
                    0,
                    fits_parser_workfn,
                    (&mut info as *mut parseInfo) as *mut c_void,
                    status,
                )
            } == -1
            {
                *status = 0;
            } else if *status != 0 {
                ffcprs(&mut lParse);
                return *status;
            }

            if info.anyNull != 0 {
                if write_blank_kwd != 0 {
                    fits_update_key_lng(
                        unsafe { outfptr.as_mut().unwrap() },
                        unsafe { cs!(c"BLANK") },
                        filter.blank,
                        Some(unsafe { cs!(c"NULL pixel value") }),
                        status,
                    );
                    if *status != 0 {
                        ffpmsg_str("pixel_filter: unable to write BLANK keyword");
                    }
                    if debug_pixfilter != 0 {
                        println!("output has NULLs");
                        println!("wrote blank [{}]", *status);
                    }
                }
            } else if bitpix > 0 {
                // never used a null
                if fits_set_imgnull(outfptr.as_mut().unwrap(), -1234554321, status) != 0 {
                    ffpmsg_str("pixel_filter: unable to reset imgnull");
                }
            }
        } else {
            // Put constant result into keyword
            let par_name = unsafe { &filter.keyword };
            let par_info: Option<&[c_char]> = if filter.comment[0] == 0 {
                None
            } else {
                Some(&filter.comment)
            };

            let result = &mut lParse.Nodes[lParse.resultNode as usize];
            match info.datatype {
                TDOUBLE => {
                    ffukyd_safe(
                        outfptr.as_mut().unwrap(),
                        par_name,
                        result.value.data.dbl,
                        15,
                        par_info,
                        status,
                    );
                }
                TLONG => {
                    ffukyj_safe(
                        outfptr.as_mut().unwrap(),
                        par_name,
                        result.value.data.lng,
                        par_info,
                        status,
                    );
                }
                TLOGICAL => {
                    ffukyl_safe(
                        outfptr.as_mut().unwrap(),
                        par_name,
                        result.value.data.log.into(),
                        par_info,
                        status,
                    );
                }
                TBIT | TSTRING => {
                    let str_val = unsafe {
                        std::slice::from_raw_parts(
                            result.value.data.astr.as_ptr(),
                            strlen_safe(&result.value.data.astr),
                        )
                    };
                    ffukys_safe(
                        outfptr.as_mut().unwrap(),
                        par_name,
                        str_val,
                        par_info,
                        status,
                    );
                }
                _ => {
                    int_snprintf!(
                        &mut msg,
                        256,
                        "pixel_filter: unexpected constant result type [{}]",
                        info.datatype,
                    );

                    ffpmsg_slice(&msg);
                }
            }
        }

        ffcprs(&mut lParse);
        *status
    }
}

fn set_image_col_types(
    lParse_status: &mut c_int,
    fptr: &mut fitsfile,
    name: &[c_char],
    bitpix: c_int,
    varInfo: &mut DataInfo,
    colIter: &mut iteratorCol,
) -> c_int {
    unsafe {
        let mut istatus: c_int = 0;
        let mut tscale: f64 = 0.0;
        let mut tzero: f64 = 0.0;
        let mut temp: [c_char; 80] = [0; 80];

        match bitpix {
            BYTE_IMG | SHORT_IMG | LONG_IMG => {
                istatus = 0;
                if fits_read_key_dbl(fptr, cs!(c"BZERO"), &mut tzero, None, &mut istatus) != 0 {
                    tzero = 0.0;
                }

                istatus = 0;
                if fits_read_key_dbl(fptr, cs!(c"BSCALE"), &mut tscale, None, &mut istatus) != 0 {
                    tscale = 1.0;
                }

                if tscale == 1.0 && (tzero == 0.0 || tzero == 32768.0) {
                    varInfo.dtype = fits_parser_yytokentype::LONG as c_int;
                    colIter.datatype = TLONG;
                } else {
                    varInfo.dtype = fits_parser_yytokentype::DOUBLE as c_int;
                    colIter.datatype = TDOUBLE;
                    if DEBUG_PIXFILTER != 0 {
                        println!(
                            "usefits_parser_yytokentype::DOUBLE as c_int for {} with BSCALE={}/BZERO={}",
                            str::from_utf8(cast_slice(name)).unwrap(),
                            tscale,
                            tzero,
                        );
                    }
                }
            }
            LONGLONG_IMG | FLOAT_IMG | DOUBLE_IMG => {
                varInfo.dtype = fits_parser_yytokentype::DOUBLE as c_int;
                colIter.datatype = TDOUBLE;
            }
            _ => {
                int_snprintf!(
                    &mut temp,
                    80,
                    "set_image_col_types: unrecognized image bitpix [{}]\n",
                    bitpix
                );
                ffpmsg_slice(&temp);
                *lParse_status = PARSE_BAD_TYPE;
                return *lParse_status;
            }
        }
        0
    }
}

/*************************************************************************

       Functions used by the evaluator to access FITS data
           (find_column, find_keywd, fits_parser_allocateCol, load_column)

*************************************************************************/

fn find_column(lParse: &mut ParseData, colName: &[c_char], itslval: *mut c_void) -> c_int {
    unsafe {
        let thelval: *mut FITS_PARSER_YYSTYPE = itslval as *mut FITS_PARSER_YYSTYPE;

        let mut status: c_int;
        let mut colnum: c_int = 0;
        let mut typecode: c_int = 0;
        let mut ktype: c_int = 0;
        let mut repeat: c_long = 0;
        let mut width: c_long = 0;
        let mut fptr: *mut fitsfile;
        let mut temp: [c_char; 80] = [0; 80];
        let mut tzero: f64 = 0.0;
        let mut tscale: f64 = 1.0;
        let mut istatus: c_int;
        let varInfo: &mut DataInfo;
        let colIter: &mut iteratorCol;

        if DEBUG_PIXFILTER != 0 {
            println!(
                "find_column({})",
                str::from_utf8(cast_slice(colName)).unwrap()
            );
        }

        if colName[0] == b'#' as c_char {
            return find_keywd(lParse, &colName[1..], itslval);
        }

        fptr = lParse.def_fptr;

        status = 0;
        let col_cnt: c_int = lParse.nCols;

        if lParse.hdutype == IMAGE_HDU {
            let i: c_int = 0;
            if lParse.pixFilter.is_null() {
                lParse.status = COL_NOT_FOUND;
                ffpmsg_str("find_column: IMAGE_HDU but no PixelFilter");
                return P_ERROR;
            }

            colnum = -1;
            for i in 0..unsafe { (*lParse.pixFilter).count } {
                if fits_strcasecmp(colName, unsafe {
                    std::slice::from_raw_parts(
                        (*lParse.pixFilter).tag.wrapping_add(i as usize).read(),
                        strlen((*lParse.pixFilter).tag.wrapping_add(i as usize).read()),
                    )
                }) == 0
                {
                    colnum = i;
                }
            }
            if colnum < 0 {
                int_snprintf!(
                    &mut temp,
                    80,
                    "find_column: PixelFilter tag {} not found",
                    str::from_utf8(cast_slice(colName)).unwrap()
                );
                ffpmsg_slice(&temp);
                lParse.status = COL_NOT_FOUND;
                return P_ERROR;
            }

            let mut tstatus = 0;
            if fits_parser_allocateCol(lParse, col_cnt, &mut tstatus) != 0 {
                lParse.status = tstatus;
                return P_ERROR;
            }
            lParse.status = tstatus;

            varInfo = &mut lParse.varData[col_cnt as usize];
            colIter = &mut lParse.colData[col_cnt as usize];

            fptr = unsafe {
                (*lParse.pixFilter)
                    .ifptr
                    .wrapping_add(colnum as usize)
                    .read()
            };

            fits_get_img_param(
                unsafe { &mut *fptr },
                MAXDIMS,
                Some(&mut typecode), /* actually bitpix */
                Some(unsafe { &mut varInfo.naxis }),
                Some(unsafe { &mut varInfo.naxes }),
                &mut status,
            );

            varInfo.nelem = 1;
            ktype = fits_parser_yytokentype::COLUMN as c_int;

            if set_image_col_types(
                &mut lParse.status,
                unsafe { &mut *fptr },
                colName,
                typecode,
                &mut *varInfo,
                colIter,
            ) != 0
            {
                return P_ERROR;
            }
            colIter.fptr = fptr;
            colIter.iotype = InputCol;
        } else {
            /* HDU holds a table */
            if lParse.compressed != 0 {
                colnum = lParse.valCol;
            } else if fits_get_colnum(
                unsafe { &mut *fptr },
                CASEINSEN.try_into().unwrap(),
                colName,
                &mut colnum,
                &mut status,
            ) != 0
            {
                if status == COL_NOT_FOUND {
                    ktype = find_keywd(lParse, colName, itslval);
                    if ktype != P_ERROR {
                        ffcmsg_safe();
                    }
                    return ktype;
                }
                lParse.status = status;
                return P_ERROR;
            }

            if fits_get_coltype(
                unsafe { &mut *fptr },
                colnum,
                Some(&mut typecode),
                Some(&mut repeat),
                Some(&mut width),
                &mut status,
            ) != 0
            {
                lParse.status = status;
                return P_ERROR;
            }

            let mut tstatus = 0;
            if fits_parser_allocateCol(lParse, col_cnt, &mut tstatus) != 0 {
                lParse.status = tstatus;
                return P_ERROR;
            }
            lParse.status = tstatus;

            varInfo = &mut lParse.varData[col_cnt as usize];
            colIter = &mut lParse.colData[col_cnt as usize];

            fits_iter_set_by_num_safe(colIter, fptr, colnum, 0, InputCol);
        }

        /*  Make sure we don't overflow variable name array  */
        strncpy_safe(&mut varInfo.name, colName, MAXVARNAME);
        varInfo.name[MAXVARNAME] = 0;

        if lParse.hdutype != IMAGE_HDU {
            match typecode {
                TBIT => {
                    varInfo.dtype = fits_parser_yytokentype::BITSTR as c_int;
                    colIter.datatype = TBYTE;
                    ktype = fits_parser_yytokentype::BITCOL as c_int;
                }
                TBYTE | TSHORT | TLONG => {
                    /* The datatype of column with TZERO and TSCALE keywords might be
                       float or double.
                    */
                    int_snprintf!(&mut temp, 80, "TZERO{}", colnum,);
                    istatus = 0;
                    if fits_read_key_dbl(
                        unsafe { &mut *fptr },
                        &temp,
                        &mut tzero,
                        None,
                        &mut istatus,
                    ) != 0
                    {
                        tzero = 0.0;
                    }
                    int_snprintf!(&mut temp, 80, "TSCAL{}", colnum,);
                    istatus = 0;
                    if fits_read_key_dbl(
                        unsafe { &mut *fptr },
                        &temp,
                        &mut tscale,
                        None,
                        &mut istatus,
                    ) != 0
                    {
                        tscale = 1.0;
                    }
                    if tscale == 1.0 && (tzero == 0.0 || tzero == 32768.0) {
                        varInfo.dtype = fits_parser_yytokentype::LONG as c_int;
                        colIter.datatype = TLONG;
                    /*    Reading an unsigned long column as a long can cause overflow errors.
                         Treat the column as a double instead.
                         } else if (tscale == 1.0 &&  tzero == 2147483648.0 ) {
                             (*varInfo).dtype     =fits_parser_yytokentype::LONG as c_int;
                             (*colIter).datatype = TULONG;
                    */
                    } else {
                        varInfo.dtype = fits_parser_yytokentype::DOUBLE as c_int;
                        colIter.datatype = TDOUBLE;
                    }
                    ktype = fits_parser_yytokentype::COLUMN as c_int;
                }
                /*
                  For now, treat 8-byte integer columns as type double.
                  This can lose precision, so the better long term solution
                  will be to add support for TLONGLONG as a separate datatype.
                */
                TLONGLONG | TFLOAT | TDOUBLE => {
                    varInfo.dtype = fits_parser_yytokentype::DOUBLE as c_int;
                    colIter.datatype = TDOUBLE;
                    ktype = fits_parser_yytokentype::COLUMN as c_int;
                }
                TLOGICAL => {
                    varInfo.dtype = fits_parser_yytokentype::BOOLEAN as c_int;
                    colIter.datatype = TLOGICAL;
                    ktype = fits_parser_yytokentype::BCOLUMN as c_int;
                }
                TSTRING => {
                    varInfo.dtype = fits_parser_yytokentype::STRING as c_int;
                    colIter.datatype = TSTRING;
                    ktype = fits_parser_yytokentype::SCOLUMN as c_int;
                    if width >= MAX_STRLEN.into() {
                        int_snprintf!(
                            &mut temp,
                            80,
                            "column {} is wider than maximum {} characters",
                            colnum,
                            MAX_STRLEN - 1,
                        );
                        ffpmsg_slice(&temp);
                        lParse.status = PARSE_LRG_VECTOR;
                        return P_ERROR;
                    }
                    if lParse.hdutype == ASCII_TBL {
                        repeat = width;
                    }
                }
                _ => {
                    if typecode < 0 {
                        int_snprintf!(
                            &mut temp,
                            80,
                            "variable-length array columns are not supported. typecode = {}",
                            typecode,
                        );
                        ffpmsg_slice(&temp);
                    }
                    lParse.status = PARSE_BAD_TYPE;
                    return P_ERROR;
                }
            }
            varInfo.nelem = repeat;
            colIter.repeat = 0; /* ffiter() will fill in this value */
            if repeat > 1 && typecode != TSTRING {
                if fits_read_tdim(
                    unsafe { &mut *fptr },
                    colnum,
                    MAXDIMS,
                    unsafe { &mut varInfo.naxis },
                    unsafe { &mut varInfo.naxes },
                    &mut status,
                ) != 0
                {
                    lParse.status = status;
                    return P_ERROR;
                }
            } else {
                varInfo.naxis = 1;
                varInfo.naxes[0] = 1;
            }
        }
        lParse.nCols += 1;
        (*thelval).lng = col_cnt as c_long;

        ktype
    }
}

fn find_keywd(lParse: &mut ParseData, keyname: &[c_char], itslval: *mut c_void) -> c_int {
    unsafe {
        let thelval: *mut FITS_PARSER_YYSTYPE = itslval as *mut FITS_PARSER_YYSTYPE;

        let mut status: c_int = 0;
        let mut ktype: c_int = 0;

        let mut keyvalue: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
        let mut dtype: c_char = 0;

        let mut rval: c_double = 0.0;
        let mut bval: c_int = 0;
        let mut ival: c_long = 0;

        status = 0;
        let fptr: *mut fitsfile = lParse.def_fptr;
        if fits_read_keyword(
            unsafe { &mut *fptr },
            keyname,
            &mut keyvalue,
            None,
            &mut status,
        ) != 0
        {
            if status == KEY_NO_EXIST {
                /*  Do this since ffgkey doesn't put an error message on stack  */
                int_snprintf!(
                    &mut keyvalue,
                    FLEN_VALUE,
                    "ffgkey could not find keyword: {}",
                    CStr::from_bytes_until_nul(cast_slice(keyname))
                        .unwrap()
                        .to_str()
                        .unwrap(),
                );
                ffpmsg_slice(&keyvalue);
            }
            lParse.status = status;
            return P_ERROR;
        }

        if fits_get_keytype(&keyvalue, &mut dtype, &mut status) != 0 {
            lParse.status = status;
            return P_ERROR;
        }

        /* Read appropriate value type and set to CONST_OP */
        match dtype as u8 {
            b'C' => {
                // 'C' as c_char
                fits_read_key_str(
                    unsafe { &mut *fptr },
                    keyname,
                    &mut keyvalue,
                    None,
                    &mut status,
                );
                ktype = fits_parser_yytokentype::STRING as c_int;
                strcpy(unsafe { (*thelval).astr.as_mut_ptr() }, keyvalue.as_ptr());
            }
            b'L' => {
                // 'L' as c_char
                fits_read_key_log(unsafe { &mut *fptr }, keyname, &mut bval, None, &mut status);
                ktype = fits_parser_yytokentype::BOOLEAN as c_int;
                unsafe {
                    (*thelval).log = bval as c_char;
                }
            }
            b'I' => {
                // 'I' as c_char
                fits_read_key_lng(unsafe { &mut *fptr }, keyname, &mut ival, None, &mut status);
                ktype = fits_parser_yytokentype::LONG as c_int;
                unsafe {
                    (*thelval).lng = ival;
                }
            }
            b'F' => {
                // 'F' as c_char
                fits_read_key_dbl(unsafe { &mut *fptr }, keyname, &mut rval, None, &mut status);
                ktype = fits_parser_yytokentype::DOUBLE as c_int;
                unsafe {
                    (*thelval).dbl = rval;
                }
            }
            _ => {
                ktype = P_ERROR;
            }
        }

        if status != 0 {
            lParse.status = status;
            return P_ERROR;
        }

        ktype
    }
}

/// Allocates parser iterator column storage for 25 columns *more* than nCols
fn fits_parser_allocateCol(lParse: &mut ParseData, nCol: c_int, status: &mut c_int) -> c_int {
    if (nCol % 25) == 0 {
        let existing_len = lParse.colData.len();
        if lParse.colData.try_reserve_exact(25).is_err() {
            lParse.colData.clear();
            lParse.varData.clear();

            *status = MEMORY_ALLOCATION;
            return *status;
        } else {
            lParse
                .colData
                .resize(existing_len + 25, iteratorCol::default());
        }

        let existing_len = lParse.varData.len();
        if lParse.varData.try_reserve_exact(25).is_err() {
            lParse.colData.clear();
            lParse.varData.clear();

            *status = MEMORY_ALLOCATION;
            return *status;
        } else {
            lParse
                .varData
                .resize(existing_len + 25, DataInfo::default());
        }
    }
    (lParse.varData[nCol as usize]).data = ptr::null_mut();
    (lParse.varData[nCol as usize]).undef = ptr::null_mut();

    0
}

fn load_column(
    lParse: &mut ParseData,
    varNum: c_int,
    fRow: c_long,
    nRows: c_long,
    data: *mut c_void,
    undef: *mut c_char,
) -> c_int {
    unsafe {
        let mut nelem: c_long = 0;
        let mut nbytes: c_long = 0;
        let row: c_long = 0;
        let len: c_long = 0;
        let mut idx: c_long = 0;
        let bitStrs: *mut *mut c_char;
        let mut msg: [c_char; 80] = [0; 80];
        let mut bytes: *mut u8;
        let mut status: c_int = 0;
        let mut anynul: c_int = 0;

        let var: &mut iteratorCol = &mut lParse.colData[varNum as usize];
        if lParse.hdutype == IMAGE_HDU {
            /* This test would need to be on a per varNum basis to support
             * cross HDU operations */
            fits_read_imgnull(
                unsafe { &mut *var.fptr },
                var.datatype,
                fRow,
                nRows,
                unsafe { std::slice::from_raw_parts_mut(data as *mut u8, (nRows * 8) as usize) }, // Assuming 8 bytes per element
                unsafe { std::slice::from_raw_parts_mut(undef, nRows as usize) },
                Some(&mut anynul),
                &mut status,
            );
            if DEBUG_PIXFILTER != 0 {
                println!(
                    "load_column: IMAGE_HDU fRow={}, nRows={} => {}",
                    fRow, nRows, status,
                );
            }
        } else {
            nelem = nRows * var.repeat;

            match var.datatype {
                TBYTE => {
                    nbytes = ((var.repeat + 7) / 8) * nRows;
                    bytes = malloc((nbytes as usize) * size_of::<c_char>()) as *mut c_uchar;

                    ffgcvb_safe(
                        unsafe { &mut *var.fptr },
                        var.colnum,
                        fRow,
                        1,
                        nbytes,
                        0,
                        unsafe { std::slice::from_raw_parts_mut(bytes, nbytes as usize) },
                        Some(&mut anynul),
                        &mut status,
                    );

                    nelem = var.repeat;
                    bitStrs = data as *mut *mut c_char;
                    for row in 0..nRows {
                        idx = (row) * ((nelem + 7) / 8) + 1;
                        for len in 0..nelem {
                            if unsafe { *bytes.wrapping_add(idx as usize) } & (1 << (7 - len % 8))
                                != 0
                            {
                                unsafe {
                                    *(*bitStrs.wrapping_add(row as usize))
                                        .wrapping_add(len as usize) = b'1' as c_char;
                                }
                            } else {
                                unsafe {
                                    *(*bitStrs.wrapping_add(row as usize))
                                        .wrapping_add(len as usize) = b'0' as c_char;
                                }
                            }
                            if len % 8 == 7 {
                                idx += 1;
                            }
                        }
                        unsafe {
                            *(*bitStrs.wrapping_add(row as usize)).wrapping_add(len as usize) = 0;
                        }
                    }

                    FREE!(bytes);
                }
                TSTRING => {
                    // Convert data to Vec<&mut [c_char]> and undef to &mut [c_char]
                    let data_ptr_array = data as *mut *mut c_char;
                    let mut string_vec = Vec::new();
                    for i in 0..nRows {
                        let str_ptr = unsafe { *data_ptr_array.wrapping_add(i as usize) };
                        let str_len = strlen(str_ptr as *const c_char);
                        let str_slice =
                            unsafe { std::slice::from_raw_parts_mut(str_ptr, str_len + 1) };
                        string_vec.push(str_slice);
                    }
                    let undef_slice =
                        unsafe { std::slice::from_raw_parts_mut(undef, nRows as usize) };

                    ffgcfs_safe(
                        unsafe { &mut *var.fptr },
                        var.colnum,
                        fRow,
                        1,
                        nRows,
                        &mut string_vec,
                        undef_slice,
                        Some(&mut anynul),
                        &mut status,
                    );
                }
                TLOGICAL => {
                    let data_slice = unsafe {
                        std::slice::from_raw_parts_mut(data as *mut c_char, nelem as usize)
                    };
                    let undef_slice =
                        unsafe { std::slice::from_raw_parts_mut(undef, nelem as usize) };

                    ffgcfl_safe(
                        unsafe { &mut *var.fptr },
                        var.colnum,
                        fRow,
                        1,
                        nelem,
                        data_slice,
                        undef_slice,
                        Some(&mut anynul),
                        &mut status,
                    );
                }
                TLONG => {
                    ffgcfj_safe(
                        &mut *var.fptr,
                        var.colnum,
                        fRow,
                        1,
                        nelem,
                        unsafe {
                            std::slice::from_raw_parts_mut(data as *mut c_long, nelem as usize)
                        },
                        unsafe { std::slice::from_raw_parts_mut(undef, nelem as usize) },
                        Some(&mut anynul),
                        &mut status,
                    );
                }
                TDOUBLE => {
                    let data_slice =
                        unsafe { std::slice::from_raw_parts_mut(data as *mut f64, nelem as usize) };
                    let undef_slice =
                        unsafe { std::slice::from_raw_parts_mut(undef, nelem as usize) };

                    ffgcfd_safe(
                        unsafe { &mut *var.fptr },
                        var.colnum,
                        fRow,
                        1,
                        nelem,
                        data_slice,
                        undef_slice,
                        Some(&mut anynul),
                        &mut status,
                    );
                }
                _ => {
                    int_snprintf!(
                        &mut msg,
                        80,
                        "load_column: unexpected datatype {}",
                        var.datatype,
                    );
                    ffpmsg_slice(&msg);
                }
            }
        }
        if status != 0 {
            lParse.status = status;
            return P_ERROR;
        }

        0
    }
}
