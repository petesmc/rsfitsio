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

use core::slice;
use core::{cmp, ptr};
use libc::{c_float, memcpy, memset};

use crate::buffers::{ffgbyt, ffgtbb_safe, ffmbyt_safe, ffpbyt, ffptbb_safe};
use crate::c_types::{c_char, c_double, c_int, c_long, c_short, c_uchar, c_void};

use crate::aliases::rust_api::{fits_binary_tform, fits_set_atblnull, fits_set_btblnull};
use crate::aliases::rust_api::{
    fits_create_img, fits_get_colnum, fits_get_coltype, fits_get_hdrspace, fits_get_hdu_type,
    fits_get_img_param, fits_get_keyclass, fits_get_keytype, fits_read_imgnull, fits_read_key_dbl,
    fits_read_key_lng, fits_read_key_log, fits_read_key_str, fits_read_keyword, fits_read_record,
    fits_read_tdim, fits_set_imgnull, fits_update_key_lng, fits_write_record,
};
use crate::cfileio::ffimport_file_safe;
use crate::editcol::{ffdrow_safe, fficol_safe, ffirow_safe};
use crate::eval_defs::{
    DataInfo, MAX_STRLEN, MAXDIMS, MAXVARNAME, Node, P_ERROR, ParseData, ParseStatusVariables,
    ValueSort, parseInfo,
};
use crate::eval_l::{
    fits_parser_yylex_destroy, fits_parser_yylex_init, fits_parser_yyrestart, yyguts_t,
};
use crate::eval_tab::{FITS_PARSER_YYSTYPE, fits_parser_yytokentype};
use crate::eval_y::{Evaluate_Parser, fits_parser_yyparse, funcOp};
use crate::fitscore::{
    ffcmph_safe, ffcmsg_safe, ffgcno_safe, ffgdesll_safe, ffgncl_safe, ffgnrw_safe, ffiblk,
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
use crate::wrappers::{strcat_safe, strcpy_safe, strlen_safe, strncpy_safe};
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
use core::mem::size_of;

pub(crate) struct ffffrw_workdata {
    prownum: *mut c_long,
    lParse: *mut ParseData,
}

static DEBUG_PIXFILTER: c_int = 0;

// Free macro for raw pointers
macro_rules! FREE {
    ($x:expr) => {
        if !$x.is_null() {
            // Every current expansion sits inside an unsafe fn or block, but
            // keep the block so the macro stays usable from safe code.
            #[allow(unused_unsafe)]
            unsafe {
                libc::free($x.cast::<libc::c_void>());
            }
            $x = core::ptr::null_mut();
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
    // FFI WRAPPER
    unsafe {
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);
        let n_good_rows = n_good_rows.as_mut().expect(NULL_MSG);

        let row_status = core::slice::from_raw_parts_mut(row_status, nrows as usize);

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
        result = (lParse.Nodes[lParse.resultNode as usize]).value.data.log();
        *n_good_rows = nrows;

        for elem in 0..nrows {
            row_status[elem as usize] = result;
        }
    } else {
        let firstrow = if firstrow > 1 { firstrow } else { 1 };
        Info.dataPtr = row_status.as_mut_ptr().cast::<c_void>();
        Info.nullPtr = ptr::null_mut();
        Info.maxRows = nrows;
        Info.parseData = &raw mut lParse;

        let colData_slice = &mut lParse.colData[..];
        if ffiter_safe(
            lParse.nCols,
            colData_slice,
            firstrow - 1,
            0,
            fits_parser_workfn,
            core::ptr::from_ref(&Info) as *mut c_void,
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
/// opened using ffreopen, so that CFITSIO can handle changing file lengths.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffsrow(
    infptr: *mut fitsfile,  /* I - Input FITS file                      */
    outfptr: *mut fitsfile, /* I - Output FITS file                     */
    expr: *const c_char,    /* I - Boolean expression                   */
    status: *mut c_int,     /* O - Error status                         */
) -> c_int {
    // FFI WRAPPER
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
        inExt.rowLength = (infptr.Fptr).rowlength as LONGLONG;
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
        outExt.rowLength = (outfptr.Fptr).rowlength as LONGLONG;
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

        /* The row-selection flags are owned here and handed to the parser as a
        raw pointer, so the error returns below cannot leak them. Everything
        after this point reaches them through Info.dataPtr, never through
        `row_flags` directly, so the pointer stays the only live borrow. */
        let mut row_flags: Vec<c_char> = Vec::new();
        if row_flags
            .try_reserve_exact((inExt.numRows + 1) as usize)
            .is_err()
        {
            ffpmsg_str("Unable to allocate memory for row selection");
            ffcprs(&mut lParse);
            *status = MEMORY_ALLOCATION;
            return *status;
        }
        /* resize also zero terminates the array */
        row_flags.resize((inExt.numRows + 1) as usize, 0);

        Info.dataPtr = row_flags.as_mut_ptr().cast::<c_void>();
        Info.nullPtr = ptr::null_mut();
        Info.maxRows = inExt.numRows as c_long;
        Info.parseData = core::ptr::from_mut::<ParseData>(&mut lParse);

        if constant != 0 {
            /*  Set all rows to the same value from constant result  */

            result = (lParse.Nodes[lParse.resultNode as usize]).value.data.log();

            for ntodo in 0..inExt.numRows {
                *Info.dataPtr.cast::<c_char>().add(ntodo as usize) = result;
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
                core::ptr::from_mut::<parseInfo>(&mut Info).cast::<c_void>(),
                status,
            );

            nGood = 0;

            for ntodo in 0..inExt.numRows {
                if (*Info.dataPtr.cast::<c_char>().add(ntodo as usize)) != 0 {
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
                ffcprs(&mut lParse);
                *status = MEMORY_ALLOCATION;
                return *status;
            } else {
                buffer.resize(cmp::max(500000, rdlen) as usize, 0);
            }

            maxrows = cmp::max(500000 / rdlen, 1);
            nbuff = 0;
            inloc = 1;
            if core::ptr::eq(infptr, outfptr) {
                /* Skip initial good rows if input==output file */
                while *Info.dataPtr.cast::<c_char>().add((inloc - 1) as usize) != 0 {
                    inloc += 1;
                }
                outloc = inloc;
            } else {
                outloc = (outExt.numRows + 1) as c_long;
                if outloc > 1 {
                    ffirow_safe(outfptr, outExt.numRows, nGood as LONGLONG, status);
                }
            }

            loop {
                if *Info.dataPtr.cast::<c_char>().add((inloc - 1) as usize) != 0 {
                    let buffer_part = &mut buffer[((rdlen * nbuff) as usize)..];

                    ffgtbb_safe(
                        infptr,
                        inloc as LONGLONG,
                        1,
                        rdlen as LONGLONG,
                        buffer_part,
                        status,
                    );
                    nbuff += 1;
                    if nbuff == maxrows {
                        ffptbb_safe(
                            outfptr,
                            outloc as LONGLONG,
                            1,
                            (rdlen * nbuff) as LONGLONG,
                            &buffer,
                            status,
                        );
                        outloc += nbuff;
                        nbuff = 0;
                    }
                }
                inloc += 1;
                if *status != 0 || (inloc as LONGLONG) > inExt.numRows {
                    break;
                }
            }

            if nbuff != 0 {
                ffptbb_safe(
                    outfptr,
                    outloc as LONGLONG,
                    1,
                    (rdlen * nbuff) as LONGLONG,
                    &buffer,
                    status,
                );
                outloc += nbuff;
            }

            if core::ptr::eq(infptr, outfptr) {
                if (outloc as LONGLONG) <= inExt.numRows {
                    ffdrow_safe(
                        infptr,
                        outloc as LONGLONG,
                        inExt.numRows - outloc as LONGLONG + 1,
                        status,
                    );
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

                if (freespace as LONGLONG - ntodo) < 0 {
                    /* not enough existing space? */
                    ntodo = (ntodo - (freespace as LONGLONG) + (BL!() - 1)) / BL!(); /* number of blocks  */
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
                    ffgbyt(
                        infptr,
                        rdlen as LONGLONG,
                        &mut buffer[..rdlen as usize],
                        status,
                    );
                    ffmbyt_safe(outfptr, outbyteloc, IGNORE_EOF, status);
                    ffpbyt(
                        outfptr,
                        rdlen as LONGLONG,
                        &buffer[..rdlen as usize],
                        status,
                    );
                    inbyteloc += rdlen as LONGLONG;
                    outbyteloc += rdlen as LONGLONG;
                    ntodo -= rdlen as LONGLONG;
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
                            for j in (outExt.numRows + 1)..=(outExt.numRows + nGood as LONGLONG) {
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

        ffcprs(&mut lParse);

        ffcmph_safe(outfptr, status); /* compress heap, deleting any orphaned data */
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
    // FFI WRAPPER
    unsafe {
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);
        let anynul = anynul.as_mut().expect(NULL_MSG);

        let array = core::slice::from_raw_parts_mut(array.cast::<u8>(), nelements as usize);

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

    Info.dataPtr = array.as_mut_ptr().cast::<c_void>();
    Info.nullPtr = nulval.cast_mut();
    Info.maxRows = nelements / nelem1;
    Info.parseData = &raw mut lParse;

    let colData_slice = &mut lParse.colData[..];
    if ffiter_safe(
        lParse.nCols,
        colData_slice,
        firstrow - 1,
        0,
        fits_parser_workfn,
        core::ptr::from_mut::<parseInfo>(&mut Info).cast::<c_void>(),
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
    // FFI WRAPPER
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
    let start: c_long = 1;
    let end: c_long = LONG_MAX;

    ffcalc_rng_safe(
        infptr,
        expr,
        outfptr,
        parName,
        parInfo,
        1,
        &[start],
        &[end],
        status,
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
    start: *const c_long,   /* I - Row range info                   */
    end: *const c_long,     /* I - Row range info                   */
    status: *mut c_int,     /* O - Error status                     */
) -> c_int {
    // FFI WRAPPER
    unsafe {
        raw_to_slice!(expr);
        raw_to_slice!(parName);
        raw_to_slice!(parInfo);
        let infptr = infptr.as_mut().expect(NULL_MSG);
        let outfptr = outfptr.as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);
        let start = slice::from_raw_parts(start, nRngs as usize);
        let end = slice::from_raw_parts(end, nRngs as usize);

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
    start: &[c_long],       /* I - Row range info                   */
    end: &[c_long],         /* I - Row range info                   */
    status: &mut c_int,     /* O - Error status                     */
) -> c_int {
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
    let mut parName = parName;

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

    Info.parseData = &raw mut lParse;
    /*  Case (1): If column exists put it there  */

    colNo = 0;
    ffpmrk_safe(); /* prevent lack of column name from sullying the stack */
    ffgcno_safe(outfptr, CASEINSEN as c_int, parName, &mut colNo, status);
    ffcmsg_safe();
    if *status != 0 {
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
            parName = &parName[1..]; /* Advance past '#' */
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
                            strcat_safe(&mut tform, cs!(c"L"));
                        }
                        TLONG => {
                            strcat_safe(&mut tform, cs!(c"J"));
                        }
                        TDOUBLE => {
                            strcat_safe(&mut tform, cs!(c"D"));
                        }
                        TSTRING => {
                            strcat_safe(&mut tform, cs!(c"A"));
                        }
                        TBIT => {
                            strcat_safe(&mut tform, cs!(c"X"));
                        }
                        TLONGLONG => {
                            strcat_safe(&mut tform, cs!(c"K"));
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
                            strcpy_safe(&mut tform, cs!(c"I11"));
                        }
                        TDOUBLE => {
                            strcpy_safe(&mut tform, cs!(c"D23.15"));
                        }
                        TSTRING | TBIT => {
                            int_snprintf!(&mut tform, 16, "A{}", nelem);
                        }
                        _ => {}
                    }
                }
                parInfo = &tform;
            } else if !((parInfo[0] as u8) as char).is_ascii_digit() && lParse.hdutype == BINARY_TBL
            {
                if Info.datatype == TBIT && parInfo[0] == b'B' as c_char {
                    nelem = (nelem + 7) / 8;
                }
                int_snprintf!(
                    &mut tform,
                    16,
                    "{}{}",
                    nelem,
                    core::str::from_utf8(cast_slice(parInfo)).unwrap_or("")
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
                        if core::mem::size_of::<c_long>() == 8 && core::mem::size_of::<c_int>() == 4
                        {
                            nullVal = INT_MIN;
                        } else {
                            nullVal = LONG_MIN;
                        }
                    } else if typecode == TLONGLONG {
                        nullVal = LONGLONG_MIN;
                    }

                    if nullVal != 0 {
                        ffpkyj_safe(outfptr, &nullKwd, nullVal, Some(cs!(c"Null value")), status);
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
            OUTPUT_COL,
        );

        lParse.nCols += 1;

        for i in 0..nRngs {
            Info.dataPtr = core::ptr::null_mut();
            Info.maxRows = end[i as usize] - start[i as usize] + 1;

            /*
               If there is only 1 range, and it includes all the rows,
               and there are 10 or more rows, then set nPerLp = 0 so
               that the iterator function will dynamically choose the
               most efficient number of rows to process in each loop.
               Otherwise, set nPerLp to the number of rows in this range.
            */

            if (Info.maxRows >= 10) && (nRngs == 1) && (start[0] == 1) && (end[0] == totaln) {
                nPerLp = 0;
            } else {
                nPerLp = Info.maxRows as c_int;
            }

            let colData_slice = &mut lParse.colData[..];
            if ffiter_safe(
                lParse.nCols,
                colData_slice,
                start[i as usize] - 1,
                c_long::from(nPerLp),
                fits_parser_workfn,
                core::ptr::from_ref(&Info) as *mut c_void,
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
                    result.value.data.dbl(),
                    15,
                    Some(parInfo),
                    status,
                );
            }
            TLONG => {
                ffukyj_safe(
                    outfptr,
                    parName,
                    result.value.data.lng() as LONGLONG,
                    Some(parInfo),
                    status,
                );
            }
            TLOGICAL => {
                ffukyl_safe(
                    outfptr,
                    parName,
                    i32::from(result.value.data.log()),
                    Some(parInfo),
                    status,
                );
            }
            TBIT | TSTRING => {
                if fits_strcasecmp(parName, cs!(c"HISTORY")) == 0 {
                    ffphis_safe(outfptr, result.value.data.text(), status);
                } else if fits_strcasecmp(parName, cs!(c"COMMENT")) == 0 {
                    ffpcom_safe(outfptr, result.value.data.text(), status);
                } else {
                    ffukys_safe(
                        outfptr,
                        parName,
                        result.value.data.text(),
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

/*--------------------------------------------------------------------------*/
/// Evaluate the given expression and return information on the result.      
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
    // FFI WRAPPER
    unsafe {
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);
        let datatype = datatype.as_mut().expect(NULL_MSG);
        let nelem = nelem.as_mut().expect(NULL_MSG);
        let naxis = naxis.as_mut().expect(NULL_MSG);
        let naxes = core::slice::from_raw_parts_mut(naxes, maxdim as usize);

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
        // Handle file import case
        let mut temp_ptr = None;
        if ffimport_file_safe(&expr[1..], &mut temp_ptr, status) != 0 {
            return *status;
        }

        if let Some(temp_ptr) = temp_ptr {
            lexpr = strlen_safe(&temp_ptr) as c_int;
            // Convert C string to Vec<u8>, appending the trailing newline that
            // the parser requires (matches C's strcat(lParse->expr+lexpr,"\n")).
            let lexpr_cstr = CStr::from_bytes_until_nul(cast_slice(&temp_ptr)).unwrap();
            let mut vec = lexpr_cstr.to_bytes().to_vec();
            vec.push(b'\n');
            vec.push(0);
            lParse.expr = Some(vec.into_boxed_slice());
        }
    } else {
        lexpr = strlen_safe(expr) as c_int;
        // Create a new boxed slice with the expression
        let mut vec = Vec::new();
        if vec.try_reserve_exact((lexpr + 2) as usize).is_err() {
            ffpmsg_str("memory allocation failed (ffiprs)");
            *status = MEMORY_ALLOCATION;
            return *status;
        }

        for i in 0..lexpr {
            vec.push(expr[i as usize] as u8);
        }
        vec.push(b'\n');
        vec.push(0);
        lParse.expr = Some(vec.into_boxed_slice());
    }
    lParse.index = 0;
    lParse.is_eobuf = 0;

    /*  Parse the expression, building the Nodes and determing  */
    /*  which columns are needed and what data type is returned  */

    let mut yylex_scanner: Option<Box<yyguts_t>> = None; /* Used internally by FLEX lexer */

    fits_parser_yylex_init(&mut yylex_scanner);
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

    let result: &Node = &lParse.Nodes[lParse.resultNode as usize];

    lParse.nAxis = result.value.naxis;
    *naxis = lParse.nAxis;
    lParse.nElements = result.value.nelem;
    *nelem = lParse.nElements;

    for i in 0..cmp::min(*naxis, maxdim) {
        lParse.nAxes[i as usize] = result.value.naxes[i as usize];
        naxes[i as usize] = lParse.nAxes[i as usize];
    }

    /* ValueSort has no other variants, so this arm cannot be reached; it is
    kept because it is the C's error path and the field may widen again. */
    #[allow(unreachable_patterns)]
    match result.ntype {
        ValueSort::Boolean => *datatype = TLOGICAL,

        ValueSort::Long => *datatype = TLONG,

        ValueSort::Double => *datatype = TDOUBLE,

        ValueSort::Bits => *datatype = TBIT,

        ValueSort::String => *datatype = TSTRING,

        _ => {
            *datatype = 0;
            ffpmsg_str("Bad return data type");
            lParse.status = PARSE_BAD_TYPE;
            *status = lParse.status;
        }
    }
    lParse.datatype = *datatype;
    lParse.expr = None; // Clear the Option<Box<[u8]>> instead of using FREE!

    if result.is_const() {
        *nelem = -*nelem;
    }
    *status
}

/*---------------------------------------------------------------------------*/
/// Clear the parser, making it ready to accept a new expression.
pub(crate) fn ffcprs(lParse: &mut ParseData) {
    unsafe {
        let col: c_int = 0;
        let mut node: c_int = 0;
        let mut i: usize = 0;

        if lParse.expr.is_some() {
            lParse.expr = None; // Clear the Option<Box<[u8]>> instead of using FREE!
        }

        if lParse.nCols > 0 {
            lParse.colData.clear();

            for col in 0..lParse.nCols {
                if lParse.varData[col as usize].undef.is_none() {
                    continue;
                }

                if lParse.varData[col as usize].dtype == ValueSort::Bits {
                    let data_ptr = lParse.varData[col as usize].data.cast::<*mut c_char>();
                    let mut first_ptr = *data_ptr;
                    FREE!(first_ptr);
                }

                lParse.varData[col as usize].undef = None;
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
                if (lParse.Nodes[node as usize]).operation == funcOp::GTIFILT_FCT as c_int {
                    i = lParse.Nodes[node as usize].SubNodes[0];
                    if !(lParse.Nodes[i]).value.data.raw().is_null() {
                        (lParse.Nodes[i]).value.data.free_buffer();
                    }
                }
                /* REGFILT nodes used to be freed here by hand. The region now
                lives in lParse.regions behind an Rc, so it goes when the
                ParseData does. */
            }
            lParse.nNodes = 0;
        }

        lParse.Nodes = Vec::new(); // Frees

        lParse.hdutype = ANY_HDU;
        lParse.pixFilter = core::ptr::null_mut();
        lParse.nDataRows = 0;
        lParse.nPrevDataRows = 0;
    }
}

/*---------------------------------------------------------------------------*/
/// Iterator work function which calls the parser and copies the results
/// into either an OutputCol or a data pointer supplied in the userPtr
/// structure.
pub(crate) extern "C" fn fits_parser_workfn(
    totalrows: c_long,         /* I - Total rows to be processed     */
    offset: c_long,            /* I - Number of rows skipped at start*/
    firstrow: c_long,          /* I - First row of this iteration    */
    nrows: c_long,             /* I - Number of rows in this iter    */
    nCols: c_int,              /* I - Number of columns in use       */
    colData: *mut iteratorCol, /* IO- Column information/data        */
    userPtr: *mut c_void,      /* I - Data handling instructions     */
) -> c_int {
    let colData = unsafe { core::slice::from_raw_parts_mut(colData, nCols as usize) };

    fits_parser_workfn_safe(totalrows, offset, firstrow, nrows, nCols, colData, userPtr)
}

/*---------------------------------------------------------------------------*/
/// Iterator work function which calls the parser and copies the results
/// into either an OutputCol or a data pointer supplied in the userPtr
/// structure.
/// `strcpy` into one of the iterator's string elements, bounded by the source.
///
/// The element width of an iterator string column is not reachable from the
/// work function: `ffiter` sets `iteratorCol::repeat` to 1 for string columns
/// and `pv.datasize` is the size of a *pointer*, not of a string. So the
/// destination slice is built to exactly the bytes a C `strcpy` would write,
/// terminator included -- the copy is bounded and its length explicit, but the
/// destination's own capacity still cannot be checked here.
///
/// # Safety
/// `dst` must point at an iterator string element with room for `src` and its
/// NUL, which is the contract the C `strcpy` already relied on.
unsafe fn copy_str_elem(dst: *mut c_char, src: &[c_char]) {
    let n = strlen_safe(src) + 1;
    unsafe { core::slice::from_raw_parts_mut(dst, n) }.copy_from_slice(&src[..n]);
}

pub(crate) fn fits_parser_workfn_safe(
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
        let lParse: &mut ParseData = userPtr
            .cast::<parseInfo>()
            .as_mut()
            .unwrap()
            .parseData
            .as_mut()
            .unwrap();
        let pv: &mut ParseStatusVariables =
            &mut userPtr.cast::<parseInfo>().as_mut().unwrap().parseVariables;

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
            (pv.userInfo) = userPtr.cast::<parseInfo>();
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
                if outcol.iotype == INPUT_COL {
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
                    if outcol.iotype != TEMPORARY_COL {
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
                        jj = c_long::from(jj_tmp);
                    }

                    if status == BAD_INTKEY || outcol.iotype == TEMPORARY_COL {
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
                ValueSort::Boolean => {
                    (pv.resDataSize) = size_of::<c_char>() as c_long;
                }
                ValueSort::Long => {
                    (pv.resDataSize) = size_of::<c_long>() as c_long;
                }
                ValueSort::Double => {
                    (pv.resDataSize) = size_of::<f64>() as c_long;
                }
                _ => {}
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
            (pv.Data) = outcol
                .array
                .cast::<c_char>()
                .add(pv.datasize.try_into().unwrap())
                .cast::<c_void>();

            /* A TemporaryCol with null value specified explicitly */
            if outcol.iotype == TEMPORARY_COL && !(*(pv.userInfo)).nullPtr.is_null() {
                pv.Null = (*(pv.userInfo)).nullPtr;
            } else {
                /* ... or an OutputCol or TemporaryCol with no explicit null */
                match (*(pv.userInfo)).datatype {
                    TLOGICAL => {
                        *pv.Null.cast::<c_char>() = b'U' as c_char;
                    }
                    TBYTE => {
                        *pv.Null.cast::<c_char>() = (pv.jnull) as c_char;
                    }
                    TSHORT => {
                        *pv.Null.cast::<c_short>() = (pv.jnull) as c_short;
                    }
                    TINT => {
                        *pv.Null.cast::<c_int>() = (pv.jnull) as c_int;
                    }
                    TLONG => {
                        *pv.Null.cast::<c_long>() = (pv.jnull) as c_long;
                    }
                    TLONGLONG => {
                        *pv.Null.cast::<LONGLONG>() = (pv.jnull) as LONGLONG;
                    }
                    TFLOAT => {
                        *pv.Null.cast::<f32>() = FLOATNULLVALUE;
                    }
                    TDOUBLE => {
                        *pv.Null.cast::<f64>() = DOUBLENULLVALUE;
                    }
                    TSTRING => {
                        (**pv.Null.cast::<*mut c_char>()) = 1;
                        *(*pv.Null.cast::<*mut c_char>()).add(1) = 0;
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
            if result.is_const() {
                constant = 1;
            }

            match result.ntype {
                ValueSort::Boolean | ValueSort::Long | ValueSort::Double => {
                    if constant != 0 {
                        let undef: c_char = 0;
                        for kk in 0..ntodo {
                            for jj in 0..pv.repeat {
                                ffcvtn(
                                    lParse.datatype,
                                    result.value.data.scalar_ptr(),
                                    &raw const undef,
                                    result.value.nelem, /* 1 */
                                    (*(pv.userInfo)).datatype,
                                    pv.Null,
                                    pv.Data
                                        .cast::<c_char>()
                                        .add(
                                            ((kk * (pv.repeat) + jj) * c_long::from(pv.datasize))
                                                .try_into()
                                                .unwrap(),
                                        )
                                        .cast::<c_void>(),
                                    &mut anyNullThisTime,
                                    &mut lParse.status,
                                );
                            }
                        }
                    } else {
                        if (pv.repeat) == result.value.nelem {
                            ffcvtn(
                                lParse.datatype,
                                result.value.data.raw(),
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
                                        result
                                            .value
                                            .data
                                            .raw()
                                            .cast::<c_char>()
                                            .add((kk * (pv.resDataSize)).try_into().unwrap())
                                            as *const c_void,
                                        result
                                            .value
                                            .undef
                                            .cast::<c_char>()
                                            .add(kk.try_into().unwrap()),
                                        1,
                                        (*(pv.userInfo)).datatype,
                                        pv.Null,
                                        pv.Data
                                            .cast::<c_char>()
                                            .add(
                                                ((kk * (pv.repeat) + jj)
                                                    * c_long::from(pv.datasize))
                                                .try_into()
                                                .unwrap(),
                                            )
                                            .cast::<c_void>(),
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
                                    (result.value.data.raw() as *const c_char)
                                        .add(
                                            (kk * result.value.nelem * (pv.resDataSize))
                                                .try_into()
                                                .unwrap(),
                                        )
                                        .cast::<c_void>(),
                                    result
                                        .value
                                        .undef
                                        .cast::<c_char>()
                                        .add((kk * result.value.nelem).try_into().unwrap()),
                                    c_long::from(nCopy),
                                    (*(pv.userInfo)).datatype,
                                    pv.Null,
                                    pv.Data
                                        .cast::<c_char>()
                                        .add(
                                            ((kk * (pv.repeat)) * c_long::from(pv.datasize))
                                                .try_into()
                                                .unwrap(),
                                        )
                                        .cast::<c_void>(),
                                    &mut anyNullThisTime,
                                    &mut lParse.status,
                                );
                                if c_long::from(nCopy) < (pv.repeat) {
                                    memset(
                                        pv.Data
                                            .cast::<c_char>()
                                            .add(
                                                ((kk * (pv.repeat) + c_long::from(nCopy))
                                                    * c_long::from(pv.datasize))
                                                .try_into()
                                                .unwrap(),
                                            )
                                            .cast::<c_void>(),
                                        0,
                                        (((pv.repeat) - c_long::from(nCopy))
                                            * c_long::from(pv.datasize))
                                            as usize,
                                    );
                                }
                            }
                        }
                        if result.is_computed() {
                            result.value.data.free_buffer();
                        }
                    }
                    if lParse.status == OVERFLOW_ERR {
                        lParse.status = NUM_OVERFLOW;
                        ffpmsg_str(
                            "Numerical overflow while converting expression to necessary datatype",
                        );
                    }
                }

                ValueSort::Bits => {
                    match (*(pv.userInfo)).datatype {
                        TBYTE => {
                            idx = -1;
                            for kk in 0..(ntodo) {
                                for jj in 0..result.value.nelem {
                                    if jj % 8 == 0 {
                                        idx += 1;
                                        *(pv.Data.cast::<c_char>().add(idx.try_into().unwrap())) =
                                            0;
                                    }
                                    if constant != 0 {
                                        if result.value.data.text()[jj as usize] == b'1' as c_char {
                                            *(pv.Data
                                                .cast::<c_uchar>()
                                                .add(idx.try_into().unwrap())) |= 128 >> (jj % 8);
                                        }
                                    } else if *(*(result
                                        .value
                                        .data
                                        .str_buf()
                                        .add(kk.try_into().unwrap())))
                                    .add(jj.try_into().unwrap())
                                        == b'1' as c_char
                                    {
                                        *(pv.Data
                                            .cast::<c_uchar>()
                                            .add(idx.try_into().unwrap())) |= 128 >> (jj % 8);
                                    }
                                }
                            }
                        }
                        TBIT | TLOGICAL => {
                            /*  Fall through to TBYTE  */
                            if constant != 0 {
                                for kk in 0..ntodo {
                                    for jj in 0..result.value.nelem {
                                        let r = if (result.value.data.text()[jj as usize])
                                            == b'1' as c_char
                                        {
                                            1
                                        } else {
                                            0
                                        };
                                        *(pv.Data.cast::<c_char>().add(
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
                                            .str_buf()
                                            .add(kk.try_into().unwrap())))
                                        .add(jj.try_into().unwrap()))
                                            == b'1' as c_char
                                        {
                                            1
                                        } else {
                                            0
                                        };
                                        *(pv.Data.cast::<c_char>().add(
                                            (jj + kk * result.value.nelem).try_into().unwrap(),
                                        )) = r;
                                    }
                                }
                            }
                        }
                        TSTRING => {
                            if constant != 0 {
                                let src = *result.value.data.text();
                                for jj in 0..ntodo {
                                    copy_str_elem(
                                        *(pv.Data
                                            .cast::<*mut c_char>()
                                            .add(jj.try_into().unwrap())),
                                        &src,
                                    );
                                }
                            } else {
                                for jj in 0..ntodo {
                                    copy_str_elem(
                                        *(pv.Data
                                            .cast::<*mut c_char>()
                                            .add(jj.try_into().unwrap())),
                                        result.value.str_row(jj),
                                    );
                                }
                            }
                        }
                        _ => {
                            ffpmsg_str("Cannot convert bit expression to desired type.");
                            lParse.status = PARSE_BAD_TYPE;
                        }
                    }
                    if result.is_computed() {
                        /* the row pointers all index one block, allocated at [0] */
                        let mut block = *(result.value.data.str_buf());
                        FREE!(block);
                        result.value.data.free_buffer();
                    }
                }
                ValueSort::String => {
                    if (*(pv.userInfo)).datatype == TSTRING {
                        if constant != 0 {
                            let src = *result.value.data.text();
                            for jj in 0..ntodo {
                                copy_str_elem(
                                    *(pv.Data.cast::<*mut c_char>().add(jj.try_into().unwrap())),
                                    &src,
                                );
                            }
                        } else {
                            for jj in 0..ntodo {
                                if *(result.value.undef.add(jj.try_into().unwrap())) != 0 {
                                    anyNullThisTime = 1;
                                    let nullstr: &[c_char] = cast_slice(
                                        CStr::from_ptr(*pv.Null.cast::<*mut c_char>())
                                            .to_bytes_with_nul(),
                                    );
                                    copy_str_elem(
                                        *(pv.Data
                                            .cast::<*mut c_char>()
                                            .add(jj.try_into().unwrap())),
                                        nullstr,
                                    );
                                } else {
                                    copy_str_elem(
                                        *(pv.Data
                                            .cast::<*mut c_char>()
                                            .add(jj.try_into().unwrap())),
                                        result.value.str_row(jj),
                                    );
                                }
                            }
                        }
                    } else {
                        ffpmsg_str("Cannot convert string expression to desired type.");
                        lParse.status = PARSE_BAD_TYPE;
                    }
                    if result.is_computed() {
                        /* the row pointers all index one block, allocated at [0] */
                        let mut block = *(result.value.data.str_buf());
                        FREE!(block);
                        result.value.data.free_buffer();
                    }
                }
            }

            if lParse.status != 0 {
                break;
            }

            /*  Increment Data to point to where the next block should go  */

            if result.ntype == ValueSort::Bits && (*(pv.userInfo)).datatype == TBYTE {
                (pv.Data) = pv
                    .Data
                    .cast::<c_char>()
                    .add(
                        (c_long::from(pv.datasize) * ((result.value.nelem + 7) / 8) * ntodo)
                            .try_into()
                            .unwrap(),
                    )
                    .cast::<c_void>();
            } else if result.ntype == ValueSort::String {
                (pv.Data) = pv
                    .Data
                    .cast::<c_char>()
                    .add((c_long::from(pv.datasize) * ntodo).try_into().unwrap())
                    .cast::<c_void>();
            } else {
                (pv.Data) = pv
                    .Data
                    .cast::<c_char>()
                    .add(
                        (c_long::from(pv.datasize) * ntodo * (pv.repeat))
                            .try_into()
                            .unwrap(),
                    )
                    .cast::<c_void>();
            }
        }

        /* If a TemporaryCol output is used, we want to inform the caller
        what the null value is expected to be */

        // WARNING - THIS IS DANGEROUS. If nCols = 0, points to invalid memory.
        // In the case of the where expr = "#ROW > 2" there are no columns.
        outcol = colData
            .as_mut_ptr()
            .offset((nCols - 1) as isize)
            .as_mut()
            .unwrap();
        // outcol = &mut colData[(nCols - 1) as usize]; // Re-bind

        if pv.Null != outcol.array
            && (Data0)
                == outcol
                    .array
                    .cast::<c_char>()
                    .add((pv.datasize).try_into().unwrap())
                    .cast::<c_void>()
        {
            if (*(pv.userInfo)).datatype == TSTRING {
                memcpy(
                    outcol.array,
                    (*pv.Null.cast::<*mut c_char>()).cast::<c_void>(),
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
                    (*pv.Null.cast::<*mut c_char>()).cast::<c_void>(),
                    zeros.as_mut_ptr().cast::<c_void>(),
                    2,
                );
            } else {
                memcpy(
                    pv.Null,
                    zeros.as_mut_ptr().cast::<c_void>(),
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

/*  Internal routines needed to allow the evaluator to operate on FITS data  */

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

        let mut bitStrs: *mut *mut c_char = core::ptr::null_mut();
        let mut sptr: *mut *mut c_char = core::ptr::null_mut();
        let mut barray: *mut c_char = core::ptr::null_mut();
        let mut iarray: *mut c_long = core::ptr::null_mut();
        let mut rarray: *mut c_double = core::ptr::null_mut();
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

            if icol.iotype == OUTPUT_COL || icol.iotype == TEMPORARY_COL {
                continue;
            }

            nelem = varData.nelem;
            len = nelem * nRows;

            /* ValueSort has no other variants, so this arm cannot be reached; it is
            kept because it is the C's error path and the field may widen again. */
            #[allow(unreachable_patterns)]
            match varData.dtype {
                ValueSort::Bits => {
                    /* No need for UNDEF array, but must make string DATA array */
                    len = (nelem + 1) * nRows; /* Count '\0' */
                    bitStrs = varData.data.cast::<*mut c_char>();
                    if do_realloc != 0 {
                        if !bitStrs.is_null() {
                            FREE!(*bitStrs);
                        }
                        free(bitStrs.cast::<c_void>());
                        bitStrs =
                            malloc(nRows as usize * size_of::<*mut c_char>()).cast::<*mut c_char>();
                        if bitStrs.is_null() {
                            varData.undef = None;
                            varData.data = ptr::null_mut();
                            lParse.status = MEMORY_ALLOCATION;
                            break;
                        }
                        *bitStrs = malloc(len as usize * size_of::<c_char>()).cast::<c_char>();
                        if (*bitStrs).is_null() {
                            free(bitStrs.cast::<c_void>());
                            varData.undef = None;
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
                            if (*icol.array.cast::<c_char>().add(idx as usize)
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
                        *(*bitStrs.add(row)).add(nelem as usize) = 0;
                    }

                    // WARNING: THIS SHOULD NOT BE COMMENTED OUT, BUT IT CAUSES A
                    // varData.undef = bitStrs.cast::<c_char>();
                    varData.data = bitStrs.cast::<c_void>();
                }

                ValueSort::String => {
                    sptr = icol.array.cast::<*mut c_char>();
                    if do_realloc != 0 {
                        if varData.undef.is_some() {
                            varData.undef = None;
                        }

                        let mut v = Vec::new();
                        if v.try_reserve_exact(len as usize).is_err() {
                            lParse.status = MEMORY_ALLOCATION;
                            break;
                        } else {
                            v.resize(len as usize, 0);
                        }

                        varData.undef = Some(v.into_boxed_slice());
                    }
                    row = nRows;
                    while row > 0 {
                        row -= 1;
                        varData.undef.as_deref_mut().unwrap()[row as usize] = c_char::from(
                            **sptr != 0
                                && FSTRCMP(
                                    cast_slice(CStr::from_ptr(*sptr.add(0)).to_bytes_with_nul()),
                                    cast_slice(
                                        CStr::from_ptr(*sptr.add((row + 1) as usize))
                                            .to_bytes_with_nul(),
                                    ),
                                ) == 0,
                        );
                    }
                    varData.data = sptr.add(1).cast::<c_void>();
                }

                ValueSort::Boolean => {
                    barray = icol.array.cast::<c_char>();
                    if do_realloc != 0 {
                        if varData.undef.is_some() {
                            varData.undef = None;
                        }
                        let mut v = Vec::new();
                        if v.try_reserve_exact(len as usize).is_err() {
                            lParse.status = MEMORY_ALLOCATION;
                            break;
                        } else {
                            v.resize(len as usize, 0);
                        }

                        varData.undef = Some(v.into_boxed_slice());
                    }
                    while len > 0 {
                        len -= 1;
                        varData.undef.as_deref_mut().unwrap()[len as usize] = c_char::from(
                            *barray.add(0) != 0
                                && *barray.add(0) == *barray.add((len + 1) as usize),
                        );
                    }
                    varData.data = barray.add(1).cast::<c_void>();
                }

                ValueSort::Long => {
                    iarray = icol.array.cast::<c_long>();
                    if do_realloc != 0 {
                        if varData.undef.is_some() {
                            varData.undef = None;
                        }
                        let mut v = Vec::new();
                        if v.try_reserve_exact(len as usize).is_err() {
                            lParse.status = MEMORY_ALLOCATION;
                            break;
                        } else {
                            v.resize(len as usize, 0);
                        }

                        varData.undef = Some(v.into_boxed_slice());
                    }
                    while len > 0 {
                        len -= 1;
                        varData.undef.as_deref_mut().unwrap()[len as usize] = c_char::from(
                            *iarray.add(0) != 0
                                && *iarray.add(0) == *iarray.add((len + 1) as usize),
                        );
                    }
                    varData.data = iarray.add(1).cast::<c_void>();
                }

                ValueSort::Double => {
                    rarray = icol.array.cast::<c_double>();
                    if do_realloc != 0 {
                        if varData.undef.is_some() {
                            varData.undef = None;
                        }

                        let mut v = Vec::new();
                        if v.try_reserve_exact(len as usize).is_err() {
                            lParse.status = MEMORY_ALLOCATION;
                            break;
                        } else {
                            v.resize(len as usize, 0);
                        }

                        varData.undef = Some(v.into_boxed_slice());
                    }
                    while len > 0 {
                        len -= 1;
                        varData.undef.as_deref_mut().unwrap()[len as usize] = c_char::from(
                            *rarray.add(0) != 0.0
                                && *rarray.add(0) == *rarray.add((len + 1) as usize),
                        );
                    }
                    varData.data = rarray.add(1).cast::<c_void>();
                }

                _ => {
                    int_snprintf!(
                        &mut msg,
                        80,
                        "SetupDataArrays, unhandled type {}\n",
                        varData.dtype.code()
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
                    if varData.dtype == ValueSort::Bits {
                        FREE!(*varData.data.cast::<*mut c_char>().add(0));
                    }
                    varData.undef = None;
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
                            if (*(input.cast::<c_uchar>().add(i.try_into().unwrap()))) != 0 {
                                *(output.cast::<c_uchar>().add(i.try_into().unwrap())) = 1;
                            } else {
                                *(output.cast::<c_uchar>().add(i.try_into().unwrap())) = 0;
                            }
                        }
                    }
                    TSHORT => {
                        for i in 0..ntodo {
                            if (*(input.cast::<c_short>().add(i.try_into().unwrap()))) != 0 {
                                *(output.cast::<c_uchar>().add(i.try_into().unwrap())) = 1;
                            } else {
                                *(output.cast::<c_uchar>().add(i.try_into().unwrap())) = 0;
                            }
                        }
                    }
                    TLONG => {
                        for i in 0..ntodo {
                            if (*(input.cast::<c_long>().add(i.try_into().unwrap()))) != 0 {
                                *(output.cast::<c_uchar>().add(i.try_into().unwrap())) = 1;
                            } else {
                                *(output.cast::<c_uchar>().add(i.try_into().unwrap())) = 0;
                            }
                        }
                    }
                    TFLOAT => {
                        for i in 0..ntodo {
                            if (*(input.cast::<c_float>().add(i.try_into().unwrap()))) != 0.0 {
                                *(output.cast::<c_uchar>().add(i.try_into().unwrap())) = 1;
                            } else {
                                *(output.cast::<c_uchar>().add(i.try_into().unwrap())) = 0;
                            }
                        }
                    }
                    TDOUBLE => {
                        for i in 0..ntodo {
                            if (*(input.cast::<c_double>().add(i.try_into().unwrap()))) != 0.0 {
                                *(output.cast::<c_uchar>().add(i.try_into().unwrap())) = 1;
                            } else {
                                *(output.cast::<c_uchar>().add(i.try_into().unwrap())) = 0;
                            }
                        }
                    }
                    _ => {
                        *status = BAD_DATATYPE;
                    }
                }
                for i in 0..ntodo {
                    if *((undef).add(i.try_into().unwrap())) != 0 {
                        *(output.cast::<c_uchar>().add(i.try_into().unwrap())) =
                            *nulval.cast::<c_uchar>();
                        *anynull = 1;
                    }
                }
            }
            TBYTE => {
                match inputType {
                    TLOGICAL | TBYTE => {
                        for i in 0..ntodo {
                            *(output.cast::<c_uchar>().add(i.try_into().unwrap())) =
                                *(input.cast::<c_uchar>().add(i.try_into().unwrap()));
                        }
                    }
                    TSHORT => {
                        let input_slice = core::slice::from_raw_parts(
                            input.cast::<c_short>(),
                            ntodo.try_into().unwrap(),
                        );

                        let output_slice = core::slice::from_raw_parts_mut(
                            output.cast::<c_uchar>(),
                            ntodo.try_into().unwrap(),
                        );

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
                                *(output.cast::<c_uchar>().add(i.try_into().unwrap())) =
                                    *nulval.cast::<c_uchar>();
                                *anynull = 1;
                            } else if *(input.cast::<c_long>().add(i.try_into().unwrap())) < 0 {
                                *status = OVERFLOW_ERR;
                                *(output.cast::<c_uchar>().add(i.try_into().unwrap())) = 0;
                            } else if *(input.cast::<c_long>().add(i.try_into().unwrap()))
                                > UCHAR_MAX as c_long
                            {
                                *status = OVERFLOW_ERR;
                                *(output.cast::<c_uchar>().add(i.try_into().unwrap())) =
                                    UCHAR_MAX as c_uchar;
                            } else {
                                *(output.cast::<c_uchar>().add(i.try_into().unwrap())) =
                                    *(input.cast::<c_long>().add(i.try_into().unwrap())) as c_uchar;
                            }
                        }
                        return *status;
                    }
                    TFLOAT => {
                        let input_slice = core::slice::from_raw_parts(
                            input.cast::<f32>(),
                            ntodo.try_into().unwrap(),
                        );

                        let output_slice = core::slice::from_raw_parts_mut(
                            output.cast::<c_uchar>(),
                            ntodo.try_into().unwrap(),
                        );

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
                        let input_slice = core::slice::from_raw_parts(
                            input.cast::<f64>(),
                            ntodo.try_into().unwrap(),
                        );

                        let output_slice = core::slice::from_raw_parts_mut(
                            output.cast::<c_uchar>(),
                            ntodo.try_into().unwrap(),
                        );

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
                        *(output.cast::<c_uchar>().add(i.try_into().unwrap())) =
                            *nulval.cast::<c_uchar>();
                        *anynull = 1;
                    }
                }
            }
            TSHORT => {
                match inputType {
                    TLOGICAL | TBYTE => {
                        for i in 0..ntodo {
                            *(output.cast::<c_short>().add(i.try_into().unwrap())) = c_short::from(
                                *(input.cast::<c_uchar>().add(i.try_into().unwrap())),
                            );
                        }
                    }
                    TSHORT => {
                        for i in 0..ntodo {
                            *(output.cast::<c_short>().add(i.try_into().unwrap())) =
                                *(input.cast::<c_short>().add(i.try_into().unwrap()));
                        }
                    }
                    TLONG => {
                        for i in 0..ntodo {
                            if *((undef).add(i.try_into().unwrap())) != 0 {
                                *(output.cast::<c_short>().add(i.try_into().unwrap())) =
                                    *nulval.cast::<c_short>();
                                *anynull = 1;
                            } else if *(input.cast::<c_long>().add(i.try_into().unwrap()))
                                < SHRT_MIN as c_long
                            {
                                *status = OVERFLOW_ERR;
                                *(output.cast::<c_short>().add(i.try_into().unwrap())) =
                                    SHRT_MIN as c_short;
                            } else if *(input.cast::<c_long>().add(i.try_into().unwrap()))
                                > SHRT_MAX as c_long
                            {
                                *status = OVERFLOW_ERR;
                                *(output.cast::<c_short>().add(i.try_into().unwrap())) =
                                    SHRT_MAX as c_short;
                            } else {
                                *(output.cast::<c_short>().add(i.try_into().unwrap())) =
                                    *(input.cast::<c_long>().add(i.try_into().unwrap())) as c_short;
                            }
                        }
                        return *status;
                    }
                    TFLOAT => {
                        let input_slice = core::slice::from_raw_parts(
                            input.cast::<f32>(),
                            ntodo.try_into().unwrap(),
                        );

                        let output_slice = core::slice::from_raw_parts_mut(
                            output.cast::<c_short>(),
                            ntodo.try_into().unwrap(),
                        );

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
                        let input_slice = core::slice::from_raw_parts(
                            input.cast::<f64>(),
                            ntodo.try_into().unwrap(),
                        );

                        let output_slice = core::slice::from_raw_parts_mut(
                            output.cast::<c_short>(),
                            ntodo.try_into().unwrap(),
                        );

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
                        *(output.cast::<c_short>().add(i.try_into().unwrap())) =
                            *nulval.cast::<c_short>();
                        *anynull = 1;
                    }
                }
            }
            TINT => {
                match inputType {
                    TLOGICAL | TBYTE => {
                        for i in 0..ntodo {
                            *(output.cast::<c_int>().add(i.try_into().unwrap())) =
                                c_int::from(*(input.cast::<c_uchar>().add(i.try_into().unwrap())));
                        }
                    }
                    TSHORT => {
                        for i in 0..ntodo {
                            *(output.cast::<c_int>().add(i.try_into().unwrap())) =
                                c_int::from(*(input.cast::<c_short>().add(i.try_into().unwrap())));
                        }
                    }
                    TLONG => {
                        for i in 0..ntodo {
                            *(output.cast::<c_int>().add(i.try_into().unwrap())) =
                                *(input.cast::<c_long>().add(i.try_into().unwrap())) as c_int;
                        }
                    }
                    TFLOAT => {
                        let input_slice = core::slice::from_raw_parts(
                            input.cast::<f32>(),
                            ntodo.try_into().unwrap(),
                        );

                        let output_slice = core::slice::from_raw_parts_mut(
                            output.cast::<c_int>(),
                            ntodo.try_into().unwrap(),
                        );

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
                        let input_slice = core::slice::from_raw_parts(
                            input.cast::<f64>(),
                            ntodo.try_into().unwrap(),
                        );

                        let output_slice = core::slice::from_raw_parts_mut(
                            output.cast::<c_int>(),
                            ntodo.try_into().unwrap(),
                        );

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
                        *(output.cast::<c_int>().add(i.try_into().unwrap())) =
                            *nulval.cast::<c_int>();
                        *anynull = 1;
                    }
                }
            }
            TLONG => {
                match inputType {
                    TLOGICAL | TBYTE => {
                        for i in 0..ntodo {
                            *(output.cast::<c_long>().add(i.try_into().unwrap())) =
                                c_long::from(*(input.cast::<c_uchar>().add(i.try_into().unwrap())));
                        }
                    }
                    TSHORT => {
                        for i in 0..ntodo {
                            *(output.cast::<c_long>().add(i.try_into().unwrap())) =
                                c_long::from(*(input.cast::<c_short>().add(i.try_into().unwrap())));
                        }
                    }
                    TLONG => {
                        for i in 0..ntodo {
                            *(output.cast::<c_long>().add(i.try_into().unwrap())) =
                                *(input.cast::<c_long>().add(i.try_into().unwrap()));
                        }
                    }
                    TFLOAT => {
                        let input_slice = core::slice::from_raw_parts(
                            input.cast::<f32>(),
                            ntodo.try_into().unwrap(),
                        );

                        let output_slice = core::slice::from_raw_parts_mut(
                            output.cast::<c_long>(),
                            ntodo.try_into().unwrap(),
                        );

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
                        let input_slice = core::slice::from_raw_parts(
                            input.cast::<f64>(),
                            ntodo.try_into().unwrap(),
                        );

                        let output_slice = core::slice::from_raw_parts_mut(
                            output.cast::<c_long>(),
                            ntodo.try_into().unwrap(),
                        );

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
                        *(output.cast::<c_long>().add(i.try_into().unwrap())) =
                            *nulval.cast::<c_long>();
                        *anynull = 1;
                    }
                }
            }
            TLONGLONG => {
                match inputType {
                    TLOGICAL | TBYTE => {
                        for i in 0..ntodo {
                            *(output.cast::<LONGLONG>().add(i.try_into().unwrap())) =
                                LONGLONG::from(
                                    *(input.cast::<c_uchar>().add(i.try_into().unwrap())),
                                );
                        }
                    }
                    TSHORT => {
                        for i in 0..ntodo {
                            *(output.cast::<LONGLONG>().add(i.try_into().unwrap())) =
                                LONGLONG::from(
                                    *(input.cast::<c_short>().add(i.try_into().unwrap())),
                                );
                        }
                    }
                    TLONG => {
                        for i in 0..ntodo {
                            *(output.cast::<LONGLONG>().add(i.try_into().unwrap())) =
                                *(input.cast::<c_long>().add(i.try_into().unwrap())) as LONGLONG;
                        }
                    }
                    TFLOAT => {
                        let input_slice = core::slice::from_raw_parts(
                            input.cast::<f32>(),
                            ntodo.try_into().unwrap(),
                        );

                        let output_slice = core::slice::from_raw_parts_mut(
                            output.cast::<LONGLONG>(),
                            ntodo.try_into().unwrap(),
                        );

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
                        let input_slice = core::slice::from_raw_parts(
                            input.cast::<f64>(),
                            ntodo.try_into().unwrap(),
                        );

                        let output_slice = core::slice::from_raw_parts_mut(
                            output.cast::<LONGLONG>(),
                            ntodo.try_into().unwrap(),
                        );

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
                        *(output.cast::<LONGLONG>().add(i.try_into().unwrap())) =
                            *nulval.cast::<LONGLONG>();
                        *anynull = 1;
                    }
                }
            }
            TFLOAT => {
                match inputType {
                    TLOGICAL | TBYTE => {
                        for i in 0..ntodo {
                            *(output.cast::<c_float>().add(i.try_into().unwrap())) = c_float::from(
                                *(input.cast::<c_uchar>().add(i.try_into().unwrap())),
                            );
                        }
                    }
                    TSHORT => {
                        for i in 0..ntodo {
                            *(output.cast::<c_float>().add(i.try_into().unwrap())) = c_float::from(
                                *(input.cast::<c_short>().add(i.try_into().unwrap())),
                            );
                        }
                    }
                    TLONG => {
                        for i in 0..ntodo {
                            *(output.cast::<c_float>().add(i.try_into().unwrap())) =
                                *(input.cast::<c_long>().add(i.try_into().unwrap())) as c_float;
                        }
                    }
                    TFLOAT => {
                        for i in 0..ntodo {
                            *(output.cast::<c_float>().add(i.try_into().unwrap())) =
                                *(input.cast::<c_float>().add(i.try_into().unwrap()));
                        }
                    }
                    TDOUBLE => {
                        let input_slice = core::slice::from_raw_parts(
                            input.cast::<f64>(),
                            ntodo.try_into().unwrap(),
                        );

                        let output_slice = core::slice::from_raw_parts_mut(
                            output.cast::<f32>(),
                            ntodo.try_into().unwrap(),
                        );

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
                        *(output.cast::<c_float>().add(i.try_into().unwrap())) =
                            *nulval.cast::<c_float>();
                        *anynull = 1;
                    }
                }
            }
            TDOUBLE => {
                match inputType {
                    TLOGICAL | TBYTE => {
                        for i in 0..ntodo {
                            *(output.cast::<c_double>().add(i.try_into().unwrap())) =
                                c_double::from(
                                    *(input.cast::<c_uchar>().add(i.try_into().unwrap())),
                                );
                        }
                    }
                    TSHORT => {
                        for i in 0..ntodo {
                            *(output.cast::<c_double>().add(i.try_into().unwrap())) =
                                c_double::from(
                                    *(input.cast::<c_short>().add(i.try_into().unwrap())),
                                );
                        }
                    }
                    TLONG => {
                        for i in 0..ntodo {
                            *(output.cast::<c_double>().add(i.try_into().unwrap())) =
                                *(input.cast::<c_long>().add(i.try_into().unwrap())) as c_double;
                        }
                    }
                    TFLOAT => {
                        for i in 0..ntodo {
                            *(output.cast::<c_double>().add(i.try_into().unwrap())) =
                                c_double::from(
                                    *(input.cast::<c_float>().add(i.try_into().unwrap())),
                                );
                        }
                    }
                    TDOUBLE => {
                        for i in 0..ntodo {
                            *(output.cast::<c_double>().add(i.try_into().unwrap())) =
                                *(input.cast::<c_double>().add(i.try_into().unwrap())) as c_double;
                        }
                    }
                    _ => {
                        *status = BAD_DATATYPE;
                    }
                }
                for i in 0..ntodo {
                    if *((undef).add(i.try_into().unwrap())) != 0 {
                        *(output.cast::<c_double>().add(i.try_into().unwrap())) =
                            *nulval.cast::<c_double>();
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
    // FFI WRAPPER
    unsafe {
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);

        let times = core::slice::from_raw_parts_mut(times, ntimes as usize);
        let time_status = core::slice::from_raw_parts_mut(time_status, ntimes as usize);

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
// The TSTRING arm dispatches on fits_get_coltype's result; hoisting that call
// into a match guard would hide it. Left as the C writes it.
#[allow(clippy::collapsible_match, clippy::collapsible_if)]
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

        /* The C wrote each array straight into lParse.colData[parNo].array.
        `iteratorCol` is Copy here, so pulling the descriptor out into a local
        and assigning through that -- as this used to -- set the field on a
        copy that was then dropped: the allocation leaked and the parser saw
        whatever array the descriptor already had. Index the vector directly.

        The arrays are owned for the length of the call instead of malloc'd,
        so the error return below and every path after it release them without
        a cleanup loop. Pushing to these outer vectors moves the inner Vec
        headers but not their heap data, so the pointers handed to colData
        stay valid. */
        let mut lng_arrays: Vec<Vec<c_long>> = Vec::new();
        let mut dbl_arrays: Vec<Vec<c_double>> = Vec::new();
        let mut str_blocks: Vec<Vec<c_char>> = Vec::new();
        let mut str_ptrs: Vec<Vec<*mut c_char>> = Vec::new();

        parNo = lParse.nCols;
        while parNo > 0 {
            parNo -= 1;
            match lParse.colData[parNo as usize].datatype {
                TLONG => {
                    let mut v: Vec<c_long> = Vec::new();
                    if v.try_reserve_exact((ntimes + 1) as usize).is_err() {
                        *status = MEMORY_ALLOCATION;
                    } else {
                        v.resize((ntimes + 1) as usize, 0);
                        v[0] = 1234554321;
                        lParse.colData[parNo as usize].array = v.as_mut_ptr().cast::<c_void>();
                        lng_arrays.push(v);
                    }
                }
                TDOUBLE => {
                    let mut v: Vec<c_double> = Vec::new();
                    if v.try_reserve_exact((ntimes + 1) as usize).is_err() {
                        *status = MEMORY_ALLOCATION;
                    } else {
                        v.resize((ntimes + 1) as usize, 0.0);
                        v[0] = DOUBLENULLVALUE;
                        lParse.colData[parNo as usize].array = v.as_mut_ptr().cast::<c_void>();
                        dbl_arrays.push(v);
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
                        let nelems = (ntimes + 1) as usize;
                        let mut block: Vec<c_char> = Vec::new();
                        let mut ptrs: Vec<*mut c_char> = Vec::new();
                        if block.try_reserve_exact(nelems * alen as usize).is_err()
                            || ptrs.try_reserve_exact(nelems).is_err()
                        {
                            *status = MEMORY_ALLOCATION;
                        } else {
                            /* one block, `alen` characters per element -- the
                            zero fill also gives element 0 its empty string */
                            block.resize(nelems * alen as usize, 0);
                            let base = block.as_mut_ptr();
                            for elem in 0..nelems {
                                ptrs.push(base.add(elem * alen as usize));
                            }
                            lParse.colData[parNo as usize].array =
                                ptrs.as_mut_ptr().cast::<c_void>();
                            str_blocks.push(block);
                            str_ptrs.push(ptrs);
                        }
                    }
                }
                _ => {}
            }

            if *status != 0 {
                return *status;
            }
        }

        /**********************************************************************/
        /* Read data from columns needed for the expression and then parse it */
        /**********************************************************************/

        if fits_uncompress_hkdata(&mut lParse, fptr, ntimes, times, status) == 0 {
            if constant != 0 {
                let result_node = &lParse.Nodes[lParse.resultNode as usize];
                result = (result_node).value.data.log();
                let mut elem = ntimes;
                while elem > 0 {
                    elem -= 1;
                    time_status[elem as usize] = result;
                }
            } else {
                Info.dataPtr = time_status.as_mut_ptr().cast::<c_void>();
                Info.nullPtr = ptr::null_mut();
                Info.maxRows = ntimes;
                *status = fits_parser_workfn_safe(
                    ntimes,
                    0,
                    1,
                    ntimes,
                    lParse.nCols,
                    &mut lParse.colData,
                    core::ptr::from_mut::<parseInfo>(&mut Info).cast::<c_void>(),
                );
            }
        }

        /************/
        /* Clean up */
        /************/

        /* No cleanup loop: lng_arrays/dbl_arrays/str_blocks/str_ptrs own the
        column arrays and drop at the end of this function. */

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
    // FFI WRAPPER
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
            result = (lParse.Nodes[lParse.resultNode as usize]).value.data.log();
            if result != 0 {
                ffgnrw_safe(fptr, &mut nelem, status);
                if nelem != 0 {
                    *rownum = 1;
                };
            };
        } else {
            let mut workData = ffffrw_workdata {
                prownum: rownum,
                lParse: core::ptr::from_mut::<ParseData>(&mut lParse),
            };
            let colData_slice = &mut lParse.colData[..];
            if ffiter_safe(
                (lParse).nCols,
                colData_slice,
                0,
                0,
                ffffrw_work,
                core::ptr::from_mut::<ffffrw_workdata>(&mut workData).cast::<c_void>(),
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
///  Setup iterator column and parser information to be ready to compute temporary calculator expression
pub(crate) fn fits_parser_set_temporary_col(
    lParse: &mut ParseData,
    Info: &mut parseInfo,
    nrows: c_long,
    nulval: *mut c_void,
    status: &mut c_int,
) -> c_int {
    if *status != 0 {
        return *status;
    }

    let col_cnt: c_int = lParse.nCols;

    if fits_parser_allocateCol(lParse, col_cnt, status) != 0 {
        return *status;
    }

    /* Set important variables for TemporaryCol where calculated results end up */

    fits_iter_set_by_num_safe(
        &mut lParse.colData[col_cnt as usize],
        ptr::null_mut(),
        0,
        TDOUBLE,
        TEMPORARY_COL,
    );

    (lParse.colData[col_cnt as usize]).repeat = lParse.nElements;

    Info.dataPtr = ptr::null_mut();
    Info.nullPtr = nulval;
    Info.maxRows = nrows;
    Info.parseData = lParse;
    lParse.nCols += 1;

    0
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
                row as LONGLONG,
                1,
                1,
                0.0,
                core::slice::from_mut(&mut newtime),
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
                            let array = col_data.array.cast::<c_long>();
                            *array.wrapping_add(currelem as usize) =
                                *array.wrapping_add((currelem - 1) as usize);
                        }
                        TDOUBLE => {
                            let array = col_data.array.cast::<f64>();
                            *array.wrapping_add(currelem as usize) =
                                *array.wrapping_add((currelem - 1) as usize);
                        }
                        TSTRING => {
                            let str_array = col_data.array.cast::<*mut c_char>();
                            let prev: &[c_char] = cast_slice(
                                CStr::from_ptr(*str_array.wrapping_add((currelem - 1) as usize))
                                    .to_bytes_with_nul(),
                            );
                            copy_str_elem(*str_array.wrapping_add(currelem as usize), prev);
                        }
                        _ => {}
                    }
                }
            }

            if ffgcvs_safe(
                fptr,
                lParse.parCol,
                row as LONGLONG,
                1,
                1,
                Some(cs!(c"")),
                &mut [core::slice::from_raw_parts_mut(sPtr[0], 256)],
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
                    core::slice::from_raw_parts(parName.as_ptr(), strlen_safe(&parName)),
                    core::slice::from_raw_parts(
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
                        let array = col_data.array.cast::<c_long>();
                        ffgcvj_safe(
                            fptr,
                            lParse.valCol,
                            row as LONGLONG,
                            1,
                            1,
                            *array.wrapping_add(0),
                            core::slice::from_raw_parts_mut(
                                array.wrapping_add(currelem as usize),
                                1,
                            ),
                            Some(&mut anynul),
                            status,
                        );
                    }
                    TDOUBLE => {
                        let array = col_data.array.cast::<f64>();
                        ffgcvd_safe(
                            fptr,
                            lParse.valCol,
                            row as LONGLONG,
                            1,
                            1,
                            *array.wrapping_add(0),
                            core::slice::from_raw_parts_mut(
                                array.wrapping_add(currelem as usize),
                                1,
                            ),
                            Some(&mut anynul),
                            status,
                        );
                    }
                    TSTRING => {
                        let str_array = col_data.array.cast::<*mut c_char>();
                        ffgcvs_safe(
                            fptr,
                            lParse.valCol,
                            row as LONGLONG,
                            1,
                            1,
                            Some(cast_slice::<u8, c_char>(
                                CStr::from_ptr(*str_array).to_bytes(),
                            )),
                            &mut [core::slice::from_raw_parts_mut(
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
    let workData: &mut ffffrw_workdata =
        unsafe { userPtr.cast::<ffffrw_workdata>().as_mut().expect(NULL_MSG) };

    ffffrw_work_safe(totalrows, offset, firstrow, nrows, nCols, workData)
}

/*---------------------------------------------------------------------------*/
/// Iterator work function which calls the parser and searches for the
/// first row which evaluates to TRUE.
pub(crate) fn ffffrw_work_safe(
    totalrows: c_long, /* I - Total rows to be processed     */
    offset: c_long,    /* I - Number of rows skipped at start*/
    firstrow: c_long,  /* I - First row of this iteration    */
    nrows: c_long,     /* I - Number of rows in this iter    */
    nCols: c_int,      /* I - Number of columns in use       */
    workData: &mut ffffrw_workdata,
) -> c_int {
    unsafe {
        let result: &mut Node;

        let lParse: &mut ParseData = &mut *workData.lParse;

        Evaluate_Parser(lParse, firstrow, nrows);

        if (lParse.status) == 0 {
            result = &mut lParse.Nodes[lParse.resultNode as usize];
            if result.is_const() {
                if result.value.data.log() != 0 {
                    *(workData.prownum) = firstrow;
                    return -1;
                }
            } else {
                for idx in 0..(nrows as usize) {
                    if *(result.value.data.log_buf().add(idx)) != 0
                        && *(result.value.undef.add(idx)) == 0
                    {
                        *(workData.prownum) = firstrow + idx as c_long;
                        return -1;
                    }
                }
            }
        }

        lParse.status
    }
}

/// Backing store for the C's `DEFAULT_TAGS`, which fits_pixel_filter hands to
/// the caller through `filter->tag`.
struct DefaultTags([*mut c_char; 1]);

/* SAFETY: the single element points into a 'static string literal and is only
ever read; nothing writes through it or through the array. */
unsafe impl Sync for DefaultTags {}

static DEFAULT_TAGS: DefaultTags = DefaultTags([c"X".as_ptr().cast_mut()]);

/*--------------------------------------------------------------------------*/
/// Apply pixel filtering operations
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn fits_pixel_filter(
    filter: *mut PixelFilter, /* I - pixel filter structure */
    status: *mut c_int,       /* IO - error status */
) -> c_int {
    // FFI WRAPPER
    unsafe {
        let status = status.as_mut().expect("Null status pointer");
        let filter = filter.as_mut().expect("Null filter pointer");

        fits_pixel_filter_safer(filter, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Apply pixel filtering operations (safe version)
/* Evaluate an expression using the data in the input FITS file(s)          */
/*--------------------------------------------------------------------------*/
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
        let result: *mut Node = core::ptr::null_mut();
        let mut datatype: c_int = 0;

        /* C: `char *DEFAULT_TAGS[] = { "X" };` -- an array holding one string
        pointer, not a string.  It is assigned to filter->tag, which the caller
        owns and reads back as tag[0][0], so it must both have pointer type and
        outlive this frame. */
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

        if filter.tag.is_null() || (*filter.tag).is_null() || **filter.tag == 0 {
            filter.tag = DEFAULT_TAGS.0.as_ptr().cast_mut();
            if debug_pixfilter != 0 {
                println!("using default tag '{}'", **filter.tag);
            }
        }

        let infptr: *mut fitsfile = *filter.ifptr;
        let outfptr: *mut fitsfile = filter.ofptr;
        lParse.pixFilter = filter;

        let filter_expr = cast_slice(CStr::from_ptr(filter.expression).to_bytes_with_nul());
        if ffiprs(
            infptr.as_mut().unwrap(),
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
            infptr.as_mut().unwrap(),
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
            /*  for floating point expressions, set the default output image to
            bitpix = -32 (float) unless the default is already a double */
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

        if fits_create_img(
            outfptr.as_mut().unwrap(),
            bitpix,
            naxis,
            &naxes[..naxis as usize],
            status,
        ) != 0
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
                infptr.as_mut().unwrap(),
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

                if fits_read_record(infptr.as_mut().unwrap(), i, Some(&mut card), status) != 0 {
                    int_snprintf!(&mut msg, 256, "pixel_filter: unable to read keycard {}", i,);

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
                } else if fits_write_record(outfptr.as_mut().unwrap(), &card, status) != 0 {
                    int_snprintf!(
                        &mut msg,
                        256,
                        "pixel_filter: unable to write keycard [{}]",
                        *status,
                    );
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
                int_snprintf!(
                    &mut msg,
                    256,
                    "pixel_filter: unexpected output bitpix {}",
                    bitpix,
                );
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
                    infptr.as_mut().unwrap(),
                    cs!(c"BLANK"),
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
                        if core::mem::size_of::<c_long>() == 8 && core::mem::size_of::<c_int>() == 4
                        {
                            null_val = INT_MIN as c_long;
                        } else {
                            null_val = LONG_MIN as c_long;
                        }
                    } else {
                        println!("unhandled positive output BITPIX {}", bitpix);
                    }
                }

                filter.blank = null_val;
            }

            fits_set_imgnull(outfptr.as_mut().unwrap(), filter.blank as LONGLONG, status);
            if debug_pixfilter != 0 {
                println!("using blank {}", null_val);
            }
        }

        if *filter.keyword.as_ptr() == 0 {
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
            colIter.fptr = filter.ofptr;
            colIter.iotype = OUTPUT_COL;

            set_image_col_types(
                &mut lParse.status,
                colIter.fptr.as_mut().unwrap(),
                cs!(c"CREATED"),
                bitpix,
                &mut lParse.varData[col_cnt as usize],
                colIter,
            );

            info.maxRows = -1;
            info.parseData = &raw mut lParse;

            if {
                let colData_slice = &mut lParse.colData[..];
                ffiter_safe(
                    lParse.nCols,
                    colData_slice,
                    0,
                    0,
                    fits_parser_workfn,
                    core::ptr::from_mut::<parseInfo>(&mut info).cast::<c_void>(),
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
                        outfptr.as_mut().unwrap(),
                        cs!(c"BLANK"),
                        filter.blank as LONGLONG,
                        Some(cs!(c"NULL pixel value")),
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
            let par_name = &filter.keyword;
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
                        result.value.data.dbl(),
                        15,
                        par_info,
                        status,
                    );
                }
                TLONG => {
                    ffukyj_safe(
                        outfptr.as_mut().unwrap(),
                        par_name,
                        result.value.data.lng() as LONGLONG,
                        par_info,
                        status,
                    );
                }
                TLOGICAL => {
                    ffukyl_safe(
                        outfptr.as_mut().unwrap(),
                        par_name,
                        result.value.data.log().into(),
                        par_info,
                        status,
                    );
                }
                TBIT | TSTRING => {
                    let str_val = core::slice::from_raw_parts(
                        result.value.data.text().as_ptr(),
                        strlen_safe(result.value.data.text()),
                    );
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
                varInfo.dtype = ValueSort::Long;
                colIter.datatype = TLONG;
            /*    Reading an unsigned long column as a long can cause overflow errors.
            Treat the column as a double instead.
            } else if (tscale == 1.0 &&  tzero == 2147483648.0 ) {
                varInfo->type     = LONG;
                colIter->datatype = TULONG;
             */
            } else {
                varInfo.dtype = ValueSort::Double;
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
            varInfo.dtype = ValueSort::Double;
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

/*************************************************************************

       Functions used by the evaluator to access FITS data
           (find_column, find_keywd, fits_parser_allocateCol, load_column)

*************************************************************************/

fn find_column(
    lParse: &mut ParseData,
    colName: &[c_char],
    itslval: &mut FITS_PARSER_YYSTYPE,
) -> c_int {
    unsafe {
        let thelval = itslval;

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
            return find_keywd(lParse, &colName[1..], thelval);
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
            for i in 0..(*lParse.pixFilter).count {
                if !(*lParse.pixFilter).tag.add(i as usize).is_null()
                    && fits_strcasecmp(
                        colName,
                        cast_slice::<u8, c_char>(
                            CStr::from_ptr((*lParse.pixFilter).tag.wrapping_add(i as usize).read())
                                .to_bytes(),
                        ),
                    ) == 0
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

            /* The C passes &lParse->status straight in, so a status already
            set stays set: allocateCol only ever writes MEMORY_ALLOCATION on
            failure, never clears. Seeding the temporary from lParse.status
            reproduces that aliasing -- writing back a fresh 0 would discard a
            pending error, which is how a syntax error raised earlier in the
            parse used to be lost. */
            let mut tstatus = lParse.status;
            if fits_parser_allocateCol(lParse, col_cnt, &mut tstatus) != 0 {
                lParse.status = tstatus;
                return P_ERROR;
            }
            lParse.status = tstatus;

            varInfo = &mut lParse.varData[col_cnt as usize];
            colIter = &mut lParse.colData[col_cnt as usize];

            if (*lParse.pixFilter).ifptr.add(colnum as usize).is_null() {
                ffpmsg_str("find_column: PixelFilter missing image pointer");
                lParse.status = COL_NOT_FOUND;
                return P_ERROR;
            }

            fptr = (*lParse.pixFilter)
                .ifptr
                .wrapping_add(colnum as usize)
                .read();

            fits_get_img_param(
                &mut *fptr,
                MAXDIMS,
                Some(&mut typecode), /* actually bitpix */
                Some(&mut varInfo.naxis),
                Some(&mut varInfo.naxes),
                &mut status,
            );

            varInfo.nelem = 1;
            ktype = fits_parser_yytokentype::COLUMN as c_int;

            if set_image_col_types(
                &mut lParse.status,
                &mut *fptr,
                colName,
                typecode,
                &mut *varInfo,
                colIter,
            ) != 0
            {
                return P_ERROR;
            }
            colIter.fptr = fptr;
            colIter.iotype = INPUT_COL;
        } else {
            /* HDU holds a table */
            if lParse.compressed != 0 {
                colnum = lParse.valCol;
            } else if fits_get_colnum(
                &mut *fptr,
                CASEINSEN.try_into().unwrap(),
                colName,
                &mut colnum,
                &mut status,
            ) != 0
            {
                if status == COL_NOT_FOUND {
                    ktype = find_keywd(lParse, colName, thelval);
                    if ktype != P_ERROR {
                        ffcmsg_safe();
                    }
                    return ktype;
                }
                lParse.status = status;
                return P_ERROR;
            }

            if fits_get_coltype(
                &mut *fptr,
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

            /* The C passes &lParse->status straight in, so a status already
            set stays set: allocateCol only ever writes MEMORY_ALLOCATION on
            failure, never clears. Seeding the temporary from lParse.status
            reproduces that aliasing -- writing back a fresh 0 would discard a
            pending error, which is how a syntax error raised earlier in the
            parse used to be lost. */
            let mut tstatus = lParse.status;
            if fits_parser_allocateCol(lParse, col_cnt, &mut tstatus) != 0 {
                lParse.status = tstatus;
                return P_ERROR;
            }
            lParse.status = tstatus;

            varInfo = &mut lParse.varData[col_cnt as usize];
            colIter = &mut lParse.colData[col_cnt as usize];

            fits_iter_set_by_num_safe(colIter, fptr, colnum, 0, INPUT_COL);
        }

        /*  Make sure we don't overflow variable name array  */
        strncpy_safe(&mut varInfo.name, colName, MAXVARNAME);
        varInfo.name[MAXVARNAME] = 0;

        if lParse.hdutype != IMAGE_HDU {
            match typecode {
                TBIT => {
                    varInfo.dtype = ValueSort::Bits;
                    colIter.datatype = TBYTE;
                    ktype = fits_parser_yytokentype::BITCOL as c_int;
                }
                TBYTE | TSHORT | TLONG => {
                    /* The datatype of column with TZERO and TSCALE keywords might be
                       float or double.
                    */
                    int_snprintf!(&mut temp, 80, "TZERO{}", colnum,);
                    istatus = 0;
                    if fits_read_key_dbl(&mut *fptr, &temp, &mut tzero, None, &mut istatus) != 0 {
                        tzero = 0.0;
                    }
                    int_snprintf!(&mut temp, 80, "TSCAL{}", colnum,);
                    istatus = 0;
                    if fits_read_key_dbl(&mut *fptr, &temp, &mut tscale, None, &mut istatus) != 0 {
                        tscale = 1.0;
                    }
                    if tscale == 1.0 && (tzero == 0.0 || tzero == 32768.0) {
                        varInfo.dtype = ValueSort::Long;
                        colIter.datatype = TLONG;
                    /*    Reading an unsigned long column as a long can cause overflow errors.
                         Treat the column as a double instead.
                         } else if (tscale == 1.0 &&  tzero == 2147483648.0 ) {
                             (*varInfo).dtype     =ValueSort::Long;
                             (*colIter).datatype = TULONG;
                    */
                    } else {
                        varInfo.dtype = ValueSort::Double;
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
                    varInfo.dtype = ValueSort::Double;
                    colIter.datatype = TDOUBLE;
                    ktype = fits_parser_yytokentype::COLUMN as c_int;
                }
                TLOGICAL => {
                    varInfo.dtype = ValueSort::Boolean;
                    colIter.datatype = TLOGICAL;
                    ktype = fits_parser_yytokentype::BCOLUMN as c_int;
                }
                TSTRING => {
                    varInfo.dtype = ValueSort::String;
                    colIter.datatype = TSTRING;
                    ktype = fits_parser_yytokentype::SCOLUMN as c_int;
                    if width >= c_long::from(MAX_STRLEN) {
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
                    &mut *fptr,
                    colnum,
                    MAXDIMS,
                    &mut varInfo.naxis,
                    &mut varInfo.naxes,
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
        *thelval = FITS_PARSER_YYSTYPE::Long(c_long::from(col_cnt));

        ktype
    }
}

fn find_keywd(
    lParse: &mut ParseData,
    keyname: &[c_char],
    itslval: &mut FITS_PARSER_YYSTYPE,
) -> c_int {
    unsafe {
        let thelval = itslval;

        let mut status: c_int = 0;
        let mut ktype: c_int = 0;

        let mut keyvalue: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
        let mut dtype: c_char = 0;

        let mut rval: c_double = 0.0;
        let mut bval: c_int = 0;
        let mut ival: c_long = 0;

        status = 0;

        if lParse.def_fptr.is_null() {
            ffpmsg_str("find_keywd: no default fitsfile defined");
            lParse.status = P_ERROR;
            return P_ERROR;
        }

        let fptr: &mut fitsfile = lParse.def_fptr.as_mut().expect(NULL_MSG);

        if fits_read_keyword(fptr, keyname, &mut keyvalue, None, &mut status) != 0 {
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
                fits_read_key_str(fptr, keyname, &mut keyvalue, None, &mut status);
                ktype = fits_parser_yytokentype::STRING as c_int;
                strcpy_safe(thelval.text_mut(), &keyvalue);
            }
            b'L' => {
                // 'L' as c_char
                fits_read_key_log(fptr, keyname, &mut bval, None, &mut status);
                ktype = fits_parser_yytokentype::BOOLEAN as c_int;
                *thelval = FITS_PARSER_YYSTYPE::Logical(bval as c_char);
            }
            b'I' => {
                // 'I' as c_char
                fits_read_key_lng(fptr, keyname, &mut ival, None, &mut status);
                ktype = fits_parser_yytokentype::LONG as c_int;
                *thelval = FITS_PARSER_YYSTYPE::Long(ival);
            }
            b'F' => {
                // 'F' as c_char
                fits_read_key_dbl(fptr, keyname, &mut rval, None, &mut status);
                ktype = fits_parser_yytokentype::DOUBLE as c_int;
                *thelval = FITS_PARSER_YYSTYPE::Double(rval);
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
                .resize_with(existing_len + 25, DataInfo::default);
        }
    }
    (lParse.varData[nCol as usize]).data = ptr::null_mut();
    (lParse.varData[nCol as usize]).undef = None;

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
                &mut *var.fptr,
                var.datatype,
                fRow as LONGLONG,
                nRows as LONGLONG,
                core::slice::from_raw_parts_mut(data.cast::<u8>(), (nRows * 8) as usize), // Assuming 8 bytes per element
                core::slice::from_raw_parts_mut(undef, nRows as usize),
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
                    let mut bytes = vec![0; nbytes as usize];

                    ffgcvb_safe(
                        &mut *var.fptr,
                        var.colnum,
                        fRow as LONGLONG,
                        1,
                        nbytes as LONGLONG,
                        0,
                        &mut bytes,
                        Some(&mut anynul),
                        &mut status,
                    );

                    nelem = var.repeat;
                    bitStrs = data.cast::<*mut c_char>();
                    for row in 0..nRows {
                        idx = (row) * ((nelem + 7) / 8) + 1;
                        for len in 0..nelem {
                            if bytes[idx as usize] & (1 << (7 - len % 8)) != 0 {
                                *(*bitStrs.wrapping_add(row as usize)).wrapping_add(len as usize) =
                                    b'1' as c_char;
                            } else {
                                *(*bitStrs.wrapping_add(row as usize)).wrapping_add(len as usize) =
                                    b'0' as c_char;
                            }
                            if len % 8 == 7 {
                                idx += 1;
                            }
                        }
                        *(*bitStrs.wrapping_add(row as usize)).wrapping_add(len as usize) = 0;
                    }
                }
                TSTRING => {
                    // Convert data to Vec<&mut [c_char]> and undef to &mut [c_char]
                    let data_ptr_array = data.cast::<*mut c_char>();
                    let mut string_vec = Vec::new();
                    for i in 0..nRows {
                        let str_ptr = *data_ptr_array.wrapping_add(i as usize);
                        let str_len = CStr::from_ptr(str_ptr).to_bytes().len();
                        let str_slice = core::slice::from_raw_parts_mut(str_ptr, str_len + 1);
                        string_vec.push(str_slice);
                    }
                    let undef_slice = core::slice::from_raw_parts_mut(undef, nRows as usize);

                    ffgcfs_safe(
                        &mut *var.fptr,
                        var.colnum,
                        fRow as LONGLONG,
                        1,
                        nRows as LONGLONG,
                        &mut string_vec,
                        undef_slice,
                        Some(&mut anynul),
                        &mut status,
                    );
                }
                TLOGICAL => {
                    let data_slice =
                        core::slice::from_raw_parts_mut(data.cast::<c_char>(), nelem as usize);
                    let undef_slice = core::slice::from_raw_parts_mut(undef, nelem as usize);

                    ffgcfl_safe(
                        &mut *var.fptr,
                        var.colnum,
                        fRow as LONGLONG,
                        1,
                        nelem as LONGLONG,
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
                        fRow as LONGLONG,
                        1,
                        nelem as LONGLONG,
                        core::slice::from_raw_parts_mut(data.cast::<c_long>(), nelem as usize),
                        core::slice::from_raw_parts_mut(undef, nelem as usize),
                        Some(&mut anynul),
                        &mut status,
                    );
                }
                TDOUBLE => {
                    let data_slice =
                        core::slice::from_raw_parts_mut(data.cast::<f64>(), nelem as usize);
                    let undef_slice = core::slice::from_raw_parts_mut(undef, nelem as usize);

                    ffgcfd_safe(
                        &mut *var.fptr,
                        var.colnum,
                        fRow as LONGLONG,
                        1,
                        nelem as LONGLONG,
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::aliases::rust_api::*;
    use crate::fitsio::{
        ASCII_TBL, BINARY_TBL, BYTE_IMG, CASEINSEN, FLEN_VALUE, LONGLONG, NOT_BTABLE,
        PARSE_BAD_TYPE, PARSE_SYNTAX_ERR, READONLY, READWRITE, TBYTE, TDOUBLE, TFLOAT, TINT,
        TLOGICAL, TLONG, TLONGLONG, TSHORT, fitsfile,
    };
    use crate::helpers::testhelpers::{to_buf, with_temp_file};
    use libc::{c_char, c_int, c_long, c_void};

    /// Make a NUL-terminated `Vec<c_char>` from a `&str`.
    fn cc(s: &str) -> Vec<c_char> {
        let mut v: Vec<c_char> = s.bytes().map(|b| b as c_char).collect();
        v.push(0);
        v
    }

    /// View a typed slice as a mutable byte slice (native layout) - this is how
    /// the C `void *array` argument of ffcrow is interpreted.
    fn as_bytes_mut<T>(s: &mut [T]) -> &mut [u8] {
        unsafe {
            core::slice::from_raw_parts_mut(s.as_mut_ptr().cast::<u8>(), core::mem::size_of_val(s))
        }
    }

    /// Build the ttype/tform argument vectors needed by fits_create_tbl.
    fn make_cols(ttype: &[&str], tform: &[&str]) -> (Vec<Option<Vec<c_char>>>, Vec<Vec<c_char>>) {
        let tt: Vec<Option<Vec<c_char>>> = ttype.iter().map(|s| Some(cc(s))).collect();
        let tf: Vec<Vec<c_char>> = tform.iter().map(|s| cc(s)).collect();
        (tt, tf)
    }

    fn create_tbl(
        f: &mut fitsfile,
        tbltype: c_int,
        nrows: LONGLONG,
        ttype: &[&str],
        tform: &[&str],
        tunit: Option<&[&str]>,
        status: &mut c_int,
    ) {
        let (tt, tf) = make_cols(ttype, tform);
        let tt_ref: Vec<Option<&[c_char]>> = tt.iter().map(|o| o.as_deref()).collect();
        let tf_ref: Vec<&[c_char]> = tf.iter().map(|v| v.as_slice()).collect();
        let tu_owned: Option<Vec<Option<Vec<c_char>>>> =
            tunit.map(|u| u.iter().map(|s| Some(cc(s))).collect());
        let tu_ref: Option<Vec<Option<&[c_char]>>> = tu_owned
            .as_ref()
            .map(|v| v.iter().map(|o| o.as_deref()).collect());
        fits_create_tbl(
            f,
            tbltype,
            nrows,
            ttype.len() as c_int,
            &tt_ref,
            &tf_ref,
            tu_ref.as_deref(),
            None,
            status,
        );
    }

    /// Create a test table (binary) with various column types for expression
    /// testing.  Returns an open file pointing at the table HDU.
    fn create_test_table(name: &[c_char]) -> Box<fitsfile> {
        let mut status = 0;
        let intdata: [c_long; 10] = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
        let floatdata: [f32; 10] = [1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5];
        let logdata: [c_char; 10] = [1, 0, 1, 0, 1, 0, 1, 0, 1, 0];
        let strdata = [
            "alpha", "beta", "gamma", "delta", "epsilon", "zeta", "eta", "theta", "iota", "kappa",
        ];

        let mut f: Option<Box<fitsfile>> = None;
        fits_create_file(&mut f, name, &mut status);
        fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
        create_tbl(
            f.as_deref_mut().unwrap(),
            BINARY_TBL,
            10,
            &["INTCOL", "FLOATCOL", "STRCOL", "BOOLCOL"],
            &["1J", "1E", "10A", "1L"],
            None,
            &mut status,
        );

        fits_write_col_lng(
            f.as_deref_mut().unwrap(),
            1,
            1,
            1,
            10,
            &intdata,
            &mut status,
        );
        fits_write_col_flt(
            f.as_deref_mut().unwrap(),
            2,
            1,
            1,
            10,
            &floatdata,
            &mut status,
        );
        for (i, s) in strdata.iter().enumerate() {
            let sv = cc(s);
            let arr: [&[c_char]; 1] = [sv.as_slice()];
            fits_write_col_str(
                f.as_deref_mut().unwrap(),
                3,
                (i + 1) as LONGLONG,
                1,
                1,
                &arr,
                &mut status,
            );
        }
        fits_write_col_log(
            f.as_deref_mut().unwrap(),
            4,
            1,
            1,
            10,
            &logdata,
            &mut status,
        );
        assert_eq!(status, 0, "create_test_table setup failed");
        f.unwrap()
    }

    /// Create an ASCII table for testing ASCII-specific code paths.
    fn create_ascii_table(name: &[c_char]) -> Box<fitsfile> {
        let mut status = 0;
        let intdata: [c_long; 10] = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
        let floatdata: [f64; 10] = [1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5];
        let strdata = [
            "alpha", "beta", "gamma", "delta", "epsilon", "zeta", "eta", "theta", "iota", "kappa",
        ];

        let mut f: Option<Box<fitsfile>> = None;
        fits_create_file(&mut f, name, &mut status);
        fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
        create_tbl(
            f.as_deref_mut().unwrap(),
            ASCII_TBL,
            10,
            &["INTCOL", "FLOATCOL", "STRCOL"],
            &["I11", "D23.15", "A10"],
            Some(&["", "", ""]),
            &mut status,
        );

        for i in 0..10 {
            fits_write_col_lng(
                f.as_deref_mut().unwrap(),
                1,
                (i + 1) as LONGLONG,
                1,
                1,
                &intdata[i..i + 1],
                &mut status,
            );
        }
        for i in 0..10 {
            fits_write_col_dbl(
                f.as_deref_mut().unwrap(),
                2,
                (i + 1) as LONGLONG,
                1,
                1,
                &floatdata[i..i + 1],
                &mut status,
            );
        }
        for (i, s) in strdata.iter().enumerate() {
            let sv = cc(s);
            let arr: [&[c_char]; 1] = [sv.as_slice()];
            fits_write_col_str(
                f.as_deref_mut().unwrap(),
                3,
                (i + 1) as LONGLONG,
                1,
                1,
                &arr,
                &mut status,
            );
        }
        assert_eq!(status, 0, "create_ascii_table setup failed");
        f.unwrap()
    }

    // ===================== fffrow tests =====================

    #[test]
    fn test_fffrow_basic() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));
            let mut n_good_rows = 0;
            let mut row_status = [0 as c_char; 10];

            /* Test simple comparison: INTCOL > 5 */
            fits_find_rows(
                &mut f,
                &cc("INTCOL > 5"),
                1,
                10,
                &mut n_good_rows,
                &mut row_status,
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(n_good_rows, 5); /* rows 6,7,8,9,10 */
            assert_eq!(row_status[0], 0); /* row 1: 1 > 5 is false */
            assert_eq!(row_status[5], 1); /* row 6: 6 > 5 is true */

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_fffrow_float_comparison() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));
            let mut n_good_rows = 0;
            let mut row_status = [0 as c_char; 10];

            fits_find_rows(
                &mut f,
                &cc("FLOATCOL < 5.0"),
                1,
                10,
                &mut n_good_rows,
                &mut row_status,
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(n_good_rows, 4); /* rows 1,2,3,4 */

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_fffrow_logical_column() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));
            let mut n_good_rows = 0;
            let mut row_status = [0 as c_char; 10];

            fits_find_rows(
                &mut f,
                &cc("BOOLCOL"),
                1,
                10,
                &mut n_good_rows,
                &mut row_status,
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(n_good_rows, 5); /* odd rows are true */

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_fffrow_compound_expression() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));
            let mut n_good_rows = 0;
            let mut row_status = [0 as c_char; 10];

            fits_find_rows(
                &mut f,
                &cc("INTCOL > 3 && INTCOL < 8"),
                1,
                10,
                &mut n_good_rows,
                &mut row_status,
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(n_good_rows, 4); /* rows 4,5,6,7 */

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_fffrow_constant_true() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));
            let mut n_good_rows = 0;
            let mut row_status = [0 as c_char; 10];

            fits_find_rows(
                &mut f,
                &cc("1 == 1"),
                1,
                10,
                &mut n_good_rows,
                &mut row_status,
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(n_good_rows, 10); /* all rows match */
            assert_eq!(row_status[0], 1);

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_fffrow_constant_false() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));
            let mut n_good_rows = 0;
            let mut row_status = [0 as c_char; 10];

            fits_find_rows(
                &mut f,
                &cc("1 == 0"),
                1,
                10,
                &mut n_good_rows,
                &mut row_status,
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(row_status[0], 0);

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_fffrow_or_expression() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));
            let mut n_good_rows = 0;
            let mut row_status = [0 as c_char; 10];

            fits_find_rows(
                &mut f,
                &cc("INTCOL == 1 || INTCOL == 10"),
                1,
                10,
                &mut n_good_rows,
                &mut row_status,
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(n_good_rows, 2);
            assert_eq!(row_status[0], 1);
            assert_eq!(row_status[9], 1);

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_fffrow_not_expression() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));
            let mut n_good_rows = 0;
            let mut row_status = [0 as c_char; 10];

            fits_find_rows(
                &mut f,
                &cc("!(INTCOL > 5)"),
                1,
                10,
                &mut n_good_rows,
                &mut row_status,
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(n_good_rows, 5); /* rows 1-5 */

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_fffrow_modulo() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));
            let mut n_good_rows = 0;
            let mut row_status = [0 as c_char; 10];

            fits_find_rows(
                &mut f,
                &cc("INTCOL % 2 == 0"),
                1,
                10,
                &mut n_good_rows,
                &mut row_status,
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(n_good_rows, 5); /* rows 2,4,6,8,10 */

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_fffrow_between() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));
            let mut n_good_rows = 0;
            let mut row_status = [0 as c_char; 10];

            fits_find_rows(
                &mut f,
                &cc("INTCOL >= 3 && INTCOL <= 7"),
                1,
                10,
                &mut n_good_rows,
                &mut row_status,
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(n_good_rows, 5); /* rows 3,4,5,6,7 */

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_fffrow_not_equal() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));
            let mut n_good_rows = 0;
            let mut row_status = [0 as c_char; 10];

            fits_find_rows(
                &mut f,
                &cc("INTCOL != 5"),
                1,
                10,
                &mut n_good_rows,
                &mut row_status,
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(n_good_rows, 9); /* all except row 5 */
            assert_eq!(row_status[4], 0); /* row 5 should be false */

            fits_close_file(f, &mut status);
        });
    }

    // ===================== fftexp tests =====================

    #[test]
    fn test_fftexp_integer() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));
            let mut datatype = 0;
            let mut naxis = 0;
            let mut nelem = 0;
            let mut naxes = [0 as c_long; 10];

            fits_test_expr(
                &mut f,
                &cc("INTCOL"),
                10,
                &mut datatype,
                &mut nelem,
                &mut naxis,
                &mut naxes,
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(datatype, TLONG);
            assert_eq!(nelem, 1);

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_fftexp_float() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));
            let mut datatype = 0;
            let mut naxis = 0;
            let mut nelem = 0;
            let mut naxes = [0 as c_long; 10];

            fits_test_expr(
                &mut f,
                &cc("FLOATCOL"),
                10,
                &mut datatype,
                &mut nelem,
                &mut naxis,
                &mut naxes,
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(datatype, TDOUBLE); /* floats promoted to double in expressions */
            assert_eq!(nelem, 1);

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_fftexp_arithmetic() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));
            let mut datatype = 0;
            let mut naxis = 0;
            let mut nelem = 0;
            let mut naxes = [0 as c_long; 10];

            fits_test_expr(
                &mut f,
                &cc("INTCOL * 2 + 1"),
                10,
                &mut datatype,
                &mut nelem,
                &mut naxis,
                &mut naxes,
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(datatype, TLONG);
            assert_eq!(nelem, 1);

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_fftexp_boolean() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));
            let mut datatype = 0;
            let mut naxis = 0;
            let mut nelem = 0;
            let mut naxes = [0 as c_long; 10];

            fits_test_expr(
                &mut f,
                &cc("INTCOL > 5"),
                10,
                &mut datatype,
                &mut nelem,
                &mut naxis,
                &mut naxes,
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(datatype, TLOGICAL);
            assert_eq!(nelem, 1);

            fits_close_file(f, &mut status);
        });
    }

    // ===================== ffcrow tests =====================

    #[test]
    fn test_ffcrow_arithmetic() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));
            let mut results = [0.0f64; 10];
            let mut anynul = 0;

            fits_calc_rows(
                &mut f,
                TDOUBLE,
                &cc("INTCOL * 2"),
                1,
                5,
                core::ptr::null(),
                as_bytes_mut(&mut results),
                &mut anynul,
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(results[0], 2.0);
            assert_eq!(results[1], 4.0);
            assert_eq!(results[4], 10.0);

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_ffcrow_add_columns() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));
            let mut results = [0.0f64; 10];
            let mut anynul = 0;

            fits_calc_rows(
                &mut f,
                TDOUBLE,
                &cc("INTCOL + FLOATCOL"),
                1,
                3,
                core::ptr::null(),
                as_bytes_mut(&mut results),
                &mut anynul,
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(results[0], 2.5);
            assert_eq!(results[1], 4.5);
            assert_eq!(results[2], 6.5);

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_ffcrow_math_functions() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));
            let mut results = [0.0f64; 10];
            let mut anynul = 0;

            fits_calc_rows(
                &mut f,
                TDOUBLE,
                &cc("ABS(INTCOL - 5)"),
                1,
                5,
                core::ptr::null(),
                as_bytes_mut(&mut results),
                &mut anynul,
                &mut status,
            );
            assert_eq!(status, 0);
            assert!((results[0] - 4.0).abs() < 0.1); /* |1-5| = 4 */
            assert!(results[4].abs() < 0.1); /* |5-5| = 0 */

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_ffcrow_min_max() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));
            let mut results = [0.0f64; 10];
            let mut anynul = 0;

            fits_calc_rows(
                &mut f,
                TDOUBLE,
                &cc("MIN(INTCOL, 3)"),
                1,
                5,
                core::ptr::null(),
                as_bytes_mut(&mut results),
                &mut anynul,
                &mut status,
            );
            assert_eq!(status, 0);
            assert!((results[0] - 1.0).abs() < 0.1); /* min(1,3) = 1 */
            assert!((results[4] - 3.0).abs() < 0.1); /* min(5,3) = 3 */

            fits_calc_rows(
                &mut f,
                TDOUBLE,
                &cc("MAX(INTCOL, 3)"),
                1,
                5,
                core::ptr::null(),
                as_bytes_mut(&mut results),
                &mut anynul,
                &mut status,
            );
            assert_eq!(status, 0);
            assert!((results[0] - 3.0).abs() < 0.1); /* max(1,3) = 3 */
            assert!((results[4] - 5.0).abs() < 0.1); /* max(5,3) = 5 */

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_ffcrow_sqrt() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));
            let mut results = [0.0f64; 10];
            let mut anynul = 0;

            fits_calc_rows(
                &mut f,
                TDOUBLE,
                &cc("SQRT(FLOATCOL)"),
                1,
                4,
                core::ptr::null(),
                as_bytes_mut(&mut results),
                &mut anynul,
                &mut status,
            );
            assert_eq!(status, 0);
            assert!(results[0] >= 1.2 && results[0] <= 1.3);
            assert!(results[3] >= 2.1 && results[3] <= 2.2); /* sqrt(4.5) ~ 2.12 */

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_ffcrow_power() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));
            let mut results = [0.0f64; 10];
            let mut anynul = 0;

            fits_calc_rows(
                &mut f,
                TDOUBLE,
                &cc("INTCOL ** 2"),
                1,
                5,
                core::ptr::null(),
                as_bytes_mut(&mut results),
                &mut anynul,
                &mut status,
            );
            assert_eq!(status, 0);
            assert!((results[0] - 1.0).abs() < 0.1);
            assert!((results[1] - 4.0).abs() < 0.1);
            assert!((results[2] - 9.0).abs() < 0.1);
            assert!((results[4] - 25.0).abs() < 0.1);

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_ffcrow_log() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));
            let mut results = [0.0f64; 10];
            let mut anynul = 0;

            fits_calc_rows(
                &mut f,
                TDOUBLE,
                &cc("LOG(INTCOL)"),
                1,
                3,
                core::ptr::null(),
                as_bytes_mut(&mut results),
                &mut anynul,
                &mut status,
            );
            assert_eq!(status, 0);
            assert!(results[0].abs() < 0.1); /* ln(1) = 0 */
            assert!(results[1] >= 0.6 && results[1] <= 0.75); /* ln(2) ~ 0.693 */

            fits_calc_rows(
                &mut f,
                TDOUBLE,
                &cc("LOG10(INTCOL)"),
                1,
                3,
                core::ptr::null(),
                as_bytes_mut(&mut results),
                &mut anynul,
                &mut status,
            );
            assert_eq!(status, 0);
            assert!(results[0].abs() < 0.1); /* log10(1) = 0 */
            assert!(results[1] >= 0.29 && results[1] <= 0.32); /* log10(2) ~ 0.301 */

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_ffcrow_trig() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));
            let mut results = [0.0f64; 10];
            let mut anynul = 0;

            fits_calc_rows(
                &mut f,
                TDOUBLE,
                &cc("SIN(FLOATCOL)"),
                1,
                3,
                core::ptr::null(),
                as_bytes_mut(&mut results),
                &mut anynul,
                &mut status,
            );
            assert_eq!(status, 0);
            assert!(results[0] >= -1.1 && results[0] <= 1.1);

            fits_calc_rows(
                &mut f,
                TDOUBLE,
                &cc("COS(FLOATCOL)"),
                1,
                3,
                core::ptr::null(),
                as_bytes_mut(&mut results),
                &mut anynul,
                &mut status,
            );
            assert_eq!(status, 0);
            assert!(results[0] >= -1.1 && results[0] <= 1.1);

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_ffcrow_conditional() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));
            let mut results = [0.0f64; 10];
            let mut anynul = 0;

            fits_calc_rows(
                &mut f,
                TDOUBLE,
                &cc("(INTCOL > 5) ? 100 : 0"),
                1,
                10,
                core::ptr::null(),
                as_bytes_mut(&mut results),
                &mut anynul,
                &mut status,
            );
            assert_eq!(status, 0);
            assert!(results[0] <= 0.1); /* row 1: false -> 0 */
            assert!(results[5] >= 99.9); /* row 6: true -> 100 */
            assert!(results[9] >= 99.9); /* row 10: true */

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_ffcrow_division() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));
            let mut results = [0.0f64; 10];
            let mut anynul = 0;

            fits_calc_rows(
                &mut f,
                TDOUBLE,
                &cc("FLOATCOL / INTCOL"),
                1,
                5,
                core::ptr::null(),
                as_bytes_mut(&mut results),
                &mut anynul,
                &mut status,
            );
            assert_eq!(status, 0);
            assert!(results[0] >= 1.4 && results[0] <= 1.6); /* 1.5 / 1 */
            assert!(results[1] >= 1.2 && results[1] <= 1.3); /* 2.5 / 2 */

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_ffcrow_subtraction() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));
            let mut results = [0.0f64; 10];
            let mut anynul = 0;

            fits_calc_rows(
                &mut f,
                TDOUBLE,
                &cc("FLOATCOL - INTCOL"),
                1,
                5,
                core::ptr::null(),
                as_bytes_mut(&mut results),
                &mut anynul,
                &mut status,
            );
            assert_eq!(status, 0);
            assert!(results[0] >= 0.4 && results[0] <= 0.6); /* 1.5 - 1 */
            assert!(results[4] >= 0.4 && results[4] <= 0.6); /* 5.5 - 5 */

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_ffcrow_row_number() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));
            let mut results = [0.0f64; 10];
            let mut anynul = 0;

            fits_calc_rows(
                &mut f,
                TDOUBLE,
                &cc("#row"),
                1,
                5,
                core::ptr::null(),
                as_bytes_mut(&mut results),
                &mut anynul,
                &mut status,
            );
            assert_eq!(status, 0);
            assert!((results[0] - 1.0).abs() < 0.1);
            assert!((results[4] - 5.0).abs() < 0.1);

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_ffcrow_negation() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));
            let mut results = [0.0f64; 10];
            let mut anynul = 0;

            fits_calc_rows(
                &mut f,
                TDOUBLE,
                &cc("-INTCOL"),
                1,
                5,
                core::ptr::null(),
                as_bytes_mut(&mut results),
                &mut anynul,
                &mut status,
            );
            assert_eq!(status, 0);
            assert!((results[0] + 1.0).abs() < 0.1); /* -1 */
            assert!((results[4] + 5.0).abs() < 0.1); /* -5 */

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_ffcrow_floor_ceil() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));
            let mut results = [0.0f64; 10];
            let mut anynul = 0;

            fits_calc_rows(
                &mut f,
                TDOUBLE,
                &cc("FLOOR(FLOATCOL)"),
                1,
                3,
                core::ptr::null(),
                as_bytes_mut(&mut results),
                &mut anynul,
                &mut status,
            );
            assert_eq!(status, 0);
            assert!((results[0] - 1.0).abs() < 0.1); /* floor(1.5) = 1 */
            assert!((results[1] - 2.0).abs() < 0.1); /* floor(2.5) = 2 */

            fits_calc_rows(
                &mut f,
                TDOUBLE,
                &cc("CEIL(FLOATCOL)"),
                1,
                3,
                core::ptr::null(),
                as_bytes_mut(&mut results),
                &mut anynul,
                &mut status,
            );
            assert_eq!(status, 0);
            assert!((results[0] - 2.0).abs() < 0.1); /* ceil(1.5) = 2 */
            assert!((results[1] - 3.0).abs() < 0.1); /* ceil(2.5) = 3 */

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_ffcrow_round() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));
            let mut results = [0.0f64; 10];
            let mut anynul = 0;

            fits_calc_rows(
                &mut f,
                TDOUBLE,
                &cc("ROUND(FLOATCOL)"),
                1,
                3,
                core::ptr::null(),
                as_bytes_mut(&mut results),
                &mut anynul,
                &mut status,
            );
            assert_eq!(status, 0);
            assert!((results[0] - 2.0).abs() < 0.1); /* round(1.5) = 2 */
            assert!((results[1] - 3.0).abs() < 0.1); /* round(2.5) = 3 */

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_ffcrow_exp() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));
            let mut results = [0.0f64; 10];
            let mut anynul = 0;

            fits_calc_rows(
                &mut f,
                TDOUBLE,
                &cc("EXP(0)"),
                1,
                1,
                core::ptr::null(),
                as_bytes_mut(&mut results),
                &mut anynul,
                &mut status,
            );
            assert_eq!(status, 0);
            assert!((results[0] - 1.0).abs() < 0.01); /* exp(0) = 1 */

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_ffcrow_atan() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));
            let mut results = [0.0f64; 10];
            let mut anynul = 0;

            fits_calc_rows(
                &mut f,
                TDOUBLE,
                &cc("ATAN(0)"),
                1,
                1,
                core::ptr::null(),
                as_bytes_mut(&mut results),
                &mut anynul,
                &mut status,
            );
            assert_eq!(status, 0);
            assert!(results[0].abs() < 0.01); /* atan(0) = 0 */

            fits_calc_rows(
                &mut f,
                TDOUBLE,
                &cc("ARCTAN2(1, 1)"),
                1,
                1,
                core::ptr::null(),
                as_bytes_mut(&mut results),
                &mut anynul,
                &mut status,
            );
            assert_eq!(status, 0);
            assert!(results[0] >= 0.78 && results[0] <= 0.79); /* pi/4 ~ 0.785 */

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_ffcrow_tan() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));
            let mut results = [0.0f64; 10];
            let mut anynul = 0;

            fits_calc_rows(
                &mut f,
                TDOUBLE,
                &cc("TAN(0)"),
                1,
                1,
                core::ptr::null(),
                as_bytes_mut(&mut results),
                &mut anynul,
                &mut status,
            );
            assert_eq!(status, 0);
            assert!(results[0].abs() < 0.01); /* tan(0) = 0 */

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_ffcrow_asin_acos() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));
            let mut results = [0.0f64; 10];
            let mut anynul = 0;

            fits_calc_rows(
                &mut f,
                TDOUBLE,
                &cc("ASIN(0)"),
                1,
                1,
                core::ptr::null(),
                as_bytes_mut(&mut results),
                &mut anynul,
                &mut status,
            );
            assert_eq!(status, 0);
            assert!(results[0].abs() < 0.01); /* asin(0) = 0 */

            fits_calc_rows(
                &mut f,
                TDOUBLE,
                &cc("ACOS(1)"),
                1,
                1,
                core::ptr::null(),
                as_bytes_mut(&mut results),
                &mut anynul,
                &mut status,
            );
            assert_eq!(status, 0);
            assert!(results[0].abs() < 0.01); /* acos(1) = 0 */

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_ffcrow_short_output() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));
            let mut results = [0i16; 10];
            let mut anynul = 0;

            fits_calc_rows(
                &mut f,
                TSHORT,
                &cc("INTCOL * 2"),
                1,
                10,
                core::ptr::null(),
                as_bytes_mut(&mut results),
                &mut anynul,
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(results[0], 2);
            assert_eq!(results[9], 20);

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_ffcrow_float_output() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));
            let mut results = [0.0f32; 10];
            let mut anynul = 0;

            fits_calc_rows(
                &mut f,
                TFLOAT,
                &cc("FLOATCOL + 1.0"),
                1,
                10,
                core::ptr::null(),
                as_bytes_mut(&mut results),
                &mut anynul,
                &mut status,
            );
            assert_eq!(status, 0);
            assert!(results[0] >= 2.4 && results[0] <= 2.6); /* 1.5 + 1.0 = 2.5 */
            assert!(results[9] >= 11.4 && results[9] <= 11.6); /* 10.5 + 1.0 = 11.5 */

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_ffcrow_long_output() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));
            let mut results = [0 as c_long; 10];
            let mut anynul = 0;

            fits_calc_rows(
                &mut f,
                TLONG,
                &cc("INTCOL + 100"),
                1,
                10,
                core::ptr::null(),
                as_bytes_mut(&mut results),
                &mut anynul,
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(results[0], 101);
            assert_eq!(results[9], 110);

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_ffcrow_logical_output() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));
            let mut results = [0 as c_char; 10];
            let mut anynul = 0;

            fits_calc_rows(
                &mut f,
                TLOGICAL,
                &cc("INTCOL > 5"),
                1,
                10,
                core::ptr::null(),
                as_bytes_mut(&mut results),
                &mut anynul,
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(results[0], 0); /* 1 > 5 false */
            assert_eq!(results[5], 1); /* 6 > 5 true */

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_ffcrow_with_nullval() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));
            let mut results = [0 as c_long; 10];
            let nulval: c_long = -999;
            let mut anynul = 0;

            fits_calc_rows(
                &mut f,
                TLONG,
                &cc("INTCOL * 2"),
                1,
                10,
                core::ptr::from_ref::<c_long>(&nulval).cast::<c_void>(),
                as_bytes_mut(&mut results),
                &mut anynul,
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(results[0], 2);
            assert_eq!(results[9], 20);
            assert_eq!(anynul, 0);

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_ffcrow_vector_expression() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let vecdata: [f64; 5] = [1.0, 2.0, 3.0, 4.0, 5.0];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            create_tbl(
                f.as_deref_mut().unwrap(),
                BINARY_TBL,
                1,
                &["VECCOL"],
                &["5D"],
                None,
                &mut status,
            );
            fits_write_col_dbl(f.as_deref_mut().unwrap(), 1, 1, 1, 5, &vecdata, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);

            let mut results = [0.0f64; 5];
            let mut anynul = 0;
            fits_calc_rows(
                f.as_deref_mut().unwrap(),
                TDOUBLE,
                &cc("VECCOL * 2"),
                1,
                5,
                core::ptr::null(),
                as_bytes_mut(&mut results),
                &mut anynul,
                &mut status,
            );
            assert_eq!(status, 0);
            assert!((results[0] - 2.0).abs() < 0.1);
            assert!((results[4] - 10.0).abs() < 0.1);

            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_ffcrow_defnull() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));
            let mut results = [0 as c_long; 10];
            let mut anynul = 0;

            fits_calc_rows(
                &mut f,
                TLONG,
                &cc("DEFNULL(INTCOL, 0)"),
                1,
                10,
                core::ptr::null(),
                as_bytes_mut(&mut results),
                &mut anynul,
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(results[0], 1);
            assert_eq!(results[9], 10);

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_ffcrow_byte_output() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));
            let mut results = [0u8; 10];
            let mut anynul = 0;

            fits_calc_rows(
                &mut f,
                TBYTE,
                &cc("INTCOL"),
                1,
                10,
                core::ptr::null(),
                &mut results,
                &mut anynul,
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(results[0], 1);
            assert_eq!(results[9], 10);

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    #[ignore = "Long long not supported on Windows"]
    fn test_ffcrow_longlong_output() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));
            let mut results = [0 as LONGLONG; 10];
            let mut anynul = 0;

            fits_calc_rows(
                &mut f,
                TLONGLONG,
                &cc("INTCOL * 1000000000"),
                1,
                10,
                core::ptr::null(),
                as_bytes_mut(&mut results),
                &mut anynul,
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(results[0], 1000000000);
            assert_eq!(results[9], 10000000000);

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_ffcrow_int_output() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));
            let mut results = [0 as c_int; 10];
            let mut anynul = 0;

            fits_calc_rows(
                &mut f,
                TINT,
                &cc("INTCOL + 1000"),
                1,
                10,
                core::ptr::null(),
                as_bytes_mut(&mut results),
                &mut anynul,
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(results[0], 1001);
            assert_eq!(results[9], 1010);

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_ffcrow_double_output() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));
            let mut results = [0.0f64; 10];
            let mut anynul = 0;

            fits_calc_rows(
                &mut f,
                TDOUBLE,
                &cc("FLOATCOL * 2.5"),
                1,
                10,
                core::ptr::null(),
                as_bytes_mut(&mut results),
                &mut anynul,
                &mut status,
            );
            assert_eq!(status, 0);
            assert!(results[0] >= 3.7 && results[0] <= 3.8); /* 1.5 * 2.5 = 3.75 */
            assert!(results[9] >= 26.2 && results[9] <= 26.3); /* 10.5 * 2.5 = 26.25 */

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_ffcrow_bitwise() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));
            let mut results = [0 as c_long; 10];
            let mut anynul = 0;

            fits_calc_rows(
                &mut f,
                TLONG,
                &cc("INTCOL & 1"),
                1,
                10,
                core::ptr::null(),
                as_bytes_mut(&mut results),
                &mut anynul,
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(results[0], 1); /* 1 & 1 = 1 */
            assert_eq!(results[1], 0); /* 2 & 1 = 0 */

            fits_calc_rows(
                &mut f,
                TLONG,
                &cc("INTCOL | 16"),
                1,
                10,
                core::ptr::null(),
                as_bytes_mut(&mut results),
                &mut anynul,
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(results[0], 17); /* 1 | 16 = 17 */
            assert_eq!(results[1], 18); /* 2 | 16 = 18 */

            fits_close_file(f, &mut status);
        });
    }

    // ===================== ffcalc tests =====================

    #[test]
    fn test_ffcalc_new_column() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));
            let mut results = [0.0f64; 10];
            let mut anynull = 0;

            let fp_self: *mut fitsfile = &raw mut *f;
            fits_calculator(
                unsafe { &mut *fp_self },
                &cc("INTCOL * 2.5"),
                unsafe { &mut *fp_self },
                &cc("DOUBLECOL"),
                &cc(""),
                &mut status,
            );
            assert_eq!(status, 0);

            fits_read_col_dbl(
                &mut f,
                5,
                1,
                1,
                10,
                0.0,
                &mut results,
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert!(results[0] >= 2.4 && results[0] <= 2.6);
            assert!(results[1] >= 4.9 && results[1] <= 5.1);
            assert!(results[9] >= 24.9 && results[9] <= 25.1);

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_ffcalc_overwrite_column() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));
            let mut results = [0.0f64; 10];
            let mut anynull = 0;

            let fp_self: *mut fitsfile = &raw mut *f;
            fits_calculator(
                unsafe { &mut *fp_self },
                &cc("INTCOL + 100"),
                unsafe { &mut *fp_self },
                &cc("FLOATCOL"),
                &cc(""),
                &mut status,
            );
            assert_eq!(status, 0);

            fits_read_col_dbl(
                &mut f,
                2,
                1,
                1,
                5,
                0.0,
                &mut results,
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert!(results[0] >= 100.9 && results[0] <= 101.1);
            assert!(results[4] >= 104.9 && results[4] <= 105.1);

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_ffcalc_with_tform() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));
            let mut results = [0.0f64; 10];
            let mut anynull = 0;

            let fp_self: *mut fitsfile = &raw mut *f;
            fits_calculator(
                unsafe { &mut *fp_self },
                &cc("INTCOL * 3.14159"),
                unsafe { &mut *fp_self },
                &cc("NEWCOL"),
                &cc("1D"),
                &mut status,
            );
            assert_eq!(status, 0);

            fits_read_col_dbl(
                &mut f,
                5,
                1,
                1,
                10,
                0.0,
                &mut results,
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert!(results[0] >= 3.1 && results[0] <= 3.2);
            assert!(results[9] >= 31.4 && results[9] <= 31.5);

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_ffcalc_keyword_output() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));
            let mut keyval = 0.0;

            let fp_self: *mut fitsfile = &raw mut *f;
            fits_calculator(
                unsafe { &mut *fp_self },
                &cc("3.14159"),
                unsafe { &mut *fp_self },
                &cc("#MYCONST"),
                &cc(""),
                &mut status,
            );
            assert_eq!(status, 0);

            fits_read_key_dbl(&mut f, &cc("MYCONST"), &mut keyval, None, &mut status);
            assert_eq!(status, 0);
            assert!(keyval >= 3.14 && keyval <= 3.15);

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_ffcalc_logical_column() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));
            let mut ncols = 0;
            let mut results = [0 as c_char; 10];
            let mut nullarr = [0 as c_char; 10];
            let mut anynul = 0;

            let fp_self: *mut fitsfile = &raw mut *f;
            fits_calculator(
                unsafe { &mut *fp_self },
                &cc("INTCOL > 5"),
                unsafe { &mut *fp_self },
                &cc("BIGVAL"),
                &cc("1L"),
                &mut status,
            );
            assert_eq!(status, 0);

            fits_get_num_cols(&mut f, &mut ncols, &mut status);
            assert_eq!(ncols, 5);

            fits_read_colnull_log(
                &mut f,
                5,
                1,
                1,
                10,
                &mut results,
                &mut nullarr,
                Some(&mut anynul),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(results[0], 0);
            assert_eq!(results[5], 1);

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_ffcalc_string_column() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));

            let fp_self: *mut fitsfile = &raw mut *f;
            fits_calculator(
                unsafe { &mut *fp_self },
                &cc("'test'"),
                unsafe { &mut *fp_self },
                &cc("STRCOL2"),
                &cc("4A"),
                &mut status,
            );
            assert_eq!(status, 0);

            let mut buffer = [0 as c_char; 80];
            let mut nullarr = [0 as c_char; 1];
            let mut anynul = 0;
            {
                let mut arr: Vec<&mut [c_char]> = vec![&mut buffer[..]];
                fits_read_colnull_str(
                    &mut f,
                    5,
                    1,
                    1,
                    1,
                    &mut arr,
                    &mut nullarr,
                    Some(&mut anynul),
                    &mut status,
                );
            }
            assert_eq!(status, 0);
            assert_eq!(crate::helpers::testhelpers::from_buf(&buffer), "test");

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_ffcalc_keyword_nonconstant_error() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));

            let fp_self: *mut fitsfile = &raw mut *f;
            fits_calculator(
                unsafe { &mut *fp_self },
                &cc("INTCOL * 2"),
                unsafe { &mut *fp_self },
                &cc("#BADKEY"),
                &cc(""),
                &mut status,
            );
            assert_ne!(status, 0); /* Should fail */

            status = 0;
            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_ffcalc_byte_column() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));
            let mut ncols = 0;
            let mut results = [0u8; 10];
            let mut anynul = 0;

            let fp_self: *mut fitsfile = &raw mut *f;
            fits_calculator(
                unsafe { &mut *fp_self },
                &cc("INTCOL"),
                unsafe { &mut *fp_self },
                &cc("BYTECOL"),
                &cc("1B"),
                &mut status,
            );
            assert_eq!(status, 0);

            fits_get_num_cols(&mut f, &mut ncols, &mut status);
            assert_eq!(ncols, 5);

            fits_read_col_byt(
                &mut f,
                5,
                1,
                1,
                10,
                0,
                &mut results,
                Some(&mut anynul),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(results[0], 1);
            assert_eq!(results[9], 10);

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_ffcalc_short_column() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));
            let mut ncols = 0;
            let mut results = [0i16; 10];
            let mut anynul = 0;

            let fp_self: *mut fitsfile = &raw mut *f;
            fits_calculator(
                unsafe { &mut *fp_self },
                &cc("INTCOL * 100"),
                unsafe { &mut *fp_self },
                &cc("SHORTCOL"),
                &cc("1I"),
                &mut status,
            );
            assert_eq!(status, 0);

            fits_get_num_cols(&mut f, &mut ncols, &mut status);
            assert_eq!(ncols, 5);

            fits_read_col_sht(
                &mut f,
                5,
                1,
                1,
                10,
                0,
                &mut results,
                Some(&mut anynul),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(results[0], 100);
            assert_eq!(results[9], 1000);

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    #[ignore = "Long long not supported on Windows"]
    fn test_ffcalc_longlong_column() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));
            let mut ncols = 0;
            let mut results = [0 as LONGLONG; 10];
            let mut anynul = 0;

            let fp_self: *mut fitsfile = &raw mut *f;
            fits_calculator(
                unsafe { &mut *fp_self },
                &cc("INTCOL * 1000000000"),
                unsafe { &mut *fp_self },
                &cc("BIGCOL"),
                &cc("1K"),
                &mut status,
            );
            assert_eq!(status, 0);

            fits_get_num_cols(&mut f, &mut ncols, &mut status);
            assert_eq!(ncols, 5);

            fits_read_col_lnglng(
                &mut f,
                5,
                1,
                1,
                10,
                0,
                &mut results,
                Some(&mut anynul),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(results[0], 1000000000);
            assert_eq!(results[9], 10000000000);

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_ffcalc_int_column() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));
            let mut ncols = 0;
            let mut results = [0 as c_int; 10];
            let mut anynul = 0;

            let fp_self: *mut fitsfile = &raw mut *f;
            fits_calculator(
                unsafe { &mut *fp_self },
                &cc("INTCOL + 500"),
                unsafe { &mut *fp_self },
                &cc("INTCOL2"),
                &cc("J"),
                &mut status,
            );
            assert_eq!(status, 0);

            fits_get_num_cols(&mut f, &mut ncols, &mut status);
            assert_eq!(ncols, 5);

            fits_read_col_int(
                &mut f,
                5,
                1,
                1,
                10,
                0,
                &mut results,
                Some(&mut anynul),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(results[0], 501);
            assert_eq!(results[9], 510);

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_ffcalc_keyword_int() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));
            let mut keyval = 0;

            let fp_self: *mut fitsfile = &raw mut *f;
            fits_calculator(
                unsafe { &mut *fp_self },
                &cc("42"),
                unsafe { &mut *fp_self },
                &cc("#INTKEY"),
                &cc(""),
                &mut status,
            );
            assert_eq!(status, 0);

            fits_read_key_lng(&mut f, &cc("INTKEY"), &mut keyval, None, &mut status);
            assert_eq!(status, 0);
            assert_eq!(keyval, 42);

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_ffcalc_keyword_logical() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));
            let mut keyval = 0;

            let fp_self: *mut fitsfile = &raw mut *f;
            fits_calculator(
                unsafe { &mut *fp_self },
                &cc("1==1"),
                unsafe { &mut *fp_self },
                &cc("#LOGKEY"),
                &cc(""),
                &mut status,
            );
            assert_eq!(status, 0);

            fits_read_key_log(&mut f, &cc("LOGKEY"), &mut keyval, None, &mut status);
            assert_eq!(status, 0);
            assert_eq!(keyval, 1);

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_ffcalc_keyword_string() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));
            let mut keyval = [0 as c_char; FLEN_VALUE];

            let fp_self: *mut fitsfile = &raw mut *f;
            fits_calculator(
                unsafe { &mut *fp_self },
                &cc("'hello'"),
                unsafe { &mut *fp_self },
                &cc("#STRKEY"),
                &cc(""),
                &mut status,
            );
            assert_eq!(status, 0);

            fits_read_key_str(&mut f, &cc("STRKEY"), &mut keyval, None, &mut status);
            assert_eq!(status, 0);
            assert_eq!(crate::helpers::testhelpers::from_buf(&keyval), "hello");

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_ffcalc_keyword_history() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));

            let fp_self: *mut fitsfile = &raw mut *f;
            fits_calculator(
                unsafe { &mut *fp_self },
                &cc("'Test history entry'"),
                unsafe { &mut *fp_self },
                &cc("#HISTORY"),
                &cc(""),
                &mut status,
            );
            assert_eq!(status, 0);

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_ffcalc_keyword_comment() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));

            let fp_self: *mut fitsfile = &raw mut *f;
            fits_calculator(
                unsafe { &mut *fp_self },
                &cc("'Test comment entry'"),
                unsafe { &mut *fp_self },
                &cc("#COMMENT"),
                &cc(""),
                &mut status,
            );
            assert_eq!(status, 0);

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_ffcalc_history_nonstring_error() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));

            let fp_self: *mut fitsfile = &raw mut *f;
            fits_calculator(
                unsafe { &mut *fp_self },
                &cc("42"),
                unsafe { &mut *fp_self },
                &cc("#HISTORY"),
                &cc(""),
                &mut status,
            );
            assert_ne!(status, 0);

            status = 0;
            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_ffcalc_logical_no_tform() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));
            let mut ncols = 0;
            let mut results = [0 as c_char; 10];
            let mut nullarr = [0 as c_char; 10];
            let mut anynul = 0;

            let fp_self: *mut fitsfile = &raw mut *f;
            fits_calculator(
                unsafe { &mut *fp_self },
                &cc("INTCOL > 5"),
                unsafe { &mut *fp_self },
                &cc("LOGCOL2"),
                &cc(""),
                &mut status,
            );
            assert_eq!(status, 0);

            fits_get_num_cols(&mut f, &mut ncols, &mut status);
            assert_eq!(ncols, 5);

            fits_read_colnull_log(
                &mut f,
                5,
                1,
                1,
                10,
                &mut results,
                &mut nullarr,
                Some(&mut anynul),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(results[0], 0);
            assert_eq!(results[5], 1);

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_ffcalc_string_no_tform() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));
            let mut ncols = 0;

            let fp_self: *mut fitsfile = &raw mut *f;
            fits_calculator(
                unsafe { &mut *fp_self },
                &cc("'test'"),
                unsafe { &mut *fp_self },
                &cc("STRCOL3"),
                &cc(""),
                &mut status,
            );
            assert_eq!(status, 0);

            fits_get_num_cols(&mut f, &mut ncols, &mut status);
            assert_eq!(ncols, 5);

            let mut buffer = [0 as c_char; 80];
            let mut nullarr = [0 as c_char; 1];
            let mut anynul = 0;
            {
                let mut arr: Vec<&mut [c_char]> = vec![&mut buffer[..]];
                fits_read_colnull_str(
                    &mut f,
                    5,
                    1,
                    1,
                    1,
                    &mut arr,
                    &mut nullarr,
                    Some(&mut anynul),
                    &mut status,
                );
            }
            assert_eq!(status, 0);
            assert_eq!(crate::helpers::testhelpers::from_buf(&buffer), "test");

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_ffcalc_double_no_tform() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));
            let mut ncols = 0;
            let mut results = [0.0f64; 10];
            let mut anynull = 0;

            let fp_self: *mut fitsfile = &raw mut *f;
            fits_calculator(
                unsafe { &mut *fp_self },
                &cc("FLOATCOL * 2.5"),
                unsafe { &mut *fp_self },
                &cc("DBLCOL"),
                &cc(""),
                &mut status,
            );
            assert_eq!(status, 0);

            fits_get_num_cols(&mut f, &mut ncols, &mut status);
            assert_eq!(ncols, 5);

            fits_read_col_dbl(
                &mut f,
                5,
                1,
                1,
                10,
                0.0,
                &mut results,
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert!(results[0] >= 3.7 && results[0] <= 3.8);

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_ffcalc_bit_column() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));

            let fp_self: *mut fitsfile = &raw mut *f;
            fits_calculator(
                unsafe { &mut *fp_self },
                &cc("b11110000"),
                unsafe { &mut *fp_self },
                &cc("BITCOL"),
                &cc(""),
                &mut status,
            );
            assert_eq!(status, 0);

            let mut colnum = 0;
            fits_get_colnum(
                &mut f,
                CASEINSEN as c_int,
                &cc("BITCOL"),
                &mut colnum,
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(colnum, 5);

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    #[ignore = "Long long not supported on Windows"]
    fn test_ffcalc_longlong_no_tform() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));
            let mut results = [0 as LONGLONG; 10];
            let mut anynul = 0;

            let fp_self: *mut fitsfile = &raw mut *f;
            fits_calculator(
                unsafe { &mut *fp_self },
                &cc("INTCOL * 1000000000"),
                unsafe { &mut *fp_self },
                &cc("BIGNUM"),
                &cc("1K"),
                &mut status,
            );
            assert_eq!(status, 0);

            fits_read_col_lnglng(
                &mut f,
                5,
                1,
                1,
                10,
                0,
                &mut results,
                Some(&mut anynul),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(results[0], 1000000000);
            assert_eq!(results[4], 5000000000);

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_ffcalc_multidim_column() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let matdata: [c_long; 9] = [1, 2, 3, 4, 5, 6, 7, 8, 9];
            let naxes: [c_long; 2] = [3, 3];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            create_tbl(
                f.as_deref_mut().unwrap(),
                BINARY_TBL,
                1,
                &["MATRIX"],
                &["9J"],
                None,
                &mut status,
            );
            fits_write_tdim(f.as_deref_mut().unwrap(), 1, 2, &naxes, &mut status);
            fits_write_col_lng(f.as_deref_mut().unwrap(), 1, 1, 1, 9, &matdata, &mut status);

            let fp: *mut fitsfile = f.as_deref_mut().unwrap();
            fits_calculator(
                unsafe { &mut *fp },
                &cc("MATRIX * 2"),
                unsafe { &mut *fp },
                &cc("DOUBLED"),
                &cc("9J"),
                &mut status,
            );
            assert_eq!(status, 0);

            let mut results = [0 as c_long; 9];
            let mut anynul = 0;
            fits_read_col_lng(
                f.as_deref_mut().unwrap(),
                2,
                1,
                1,
                9,
                0,
                &mut results,
                Some(&mut anynul),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(results[0], 2);
            assert_eq!(results[8], 18);

            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    // ===================== ffcalc_rng tests =====================

    #[test]
    fn test_ffcalc_rng_basic() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));
            let mut results = [0.0f64; 10];
            let mut anynull = 0;
            let start: [c_long; 1] = [3];
            let end: [c_long; 1] = [7];

            let fp_self: *mut fitsfile = &raw mut *f;
            fits_calculator_rng(
                unsafe { &mut *fp_self },
                &cc("INTCOL * 10"),
                unsafe { &mut *fp_self },
                &cc("COMPUTED"),
                &cc(""),
                1,
                &start,
                &end,
                &mut status,
            );
            assert_eq!(status, 0);

            fits_read_col_dbl(
                &mut f,
                5,
                3,
                1,
                5,
                0.0,
                &mut results,
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert!(results[0] >= 29.9 && results[0] <= 30.1); /* 3 * 10 */
            assert!(results[4] >= 69.9 && results[4] <= 70.1); /* 7 * 10 */

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_ffcalc_rng_multiple_ranges() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));
            let mut results = [0.0f64; 10];
            let mut anynull = 0;
            let start: [c_long; 2] = [1, 8];
            let end: [c_long; 2] = [3, 10];

            let fp_self: *mut fitsfile = &raw mut *f;
            fits_calculator_rng(
                unsafe { &mut *fp_self },
                &cc("INTCOL * 5"),
                unsafe { &mut *fp_self },
                &cc("RANGED"),
                &cc(""),
                2,
                &start,
                &end,
                &mut status,
            );
            assert_eq!(status, 0);

            fits_read_col_dbl(
                &mut f,
                5,
                1,
                1,
                3,
                0.0,
                &mut results,
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert!(results[0] >= 4.9 && results[0] <= 5.1); /* 1 * 5 */
            assert!(results[2] >= 14.9 && results[2] <= 15.1); /* 3 * 5 */

            fits_close_file(f, &mut status);
        });
    }

    // ===================== string / error tests =====================

    #[test]
    fn test_ffcrow_string_comparison() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));
            let mut n_good_rows = 0;
            let mut row_status = [0 as c_char; 10];

            fits_find_rows(
                &mut f,
                &cc("STRCOL == 'alpha'"),
                1,
                10,
                &mut n_good_rows,
                &mut row_status,
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(n_good_rows, 1);
            assert_eq!(row_status[0], 1);

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_invalid_expression() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));
            let mut n_good_rows = 0;
            let mut row_status = [0 as c_char; 10];

            fits_find_rows(
                &mut f,
                &cc("NONEXISTENT > 5"),
                1,
                10,
                &mut n_good_rows,
                &mut row_status,
                &mut status,
            );
            assert_ne!(status, 0);

            status = 0;
            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_non_boolean_expression() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));
            let mut n_good_rows = 0;
            let mut row_status = [0 as c_char; 10];

            fits_find_rows(
                &mut f,
                &cc("INTCOL + 1"),
                1,
                10,
                &mut n_good_rows,
                &mut row_status,
                &mut status,
            );
            assert_ne!(status, 0);

            status = 0;
            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_fffrow_expr_from_file() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));
            let mut n_good_rows = 0;
            let mut row_status = [0 as c_char; 10];

            // Create expression file
            let exprpath = format!("{}.expr.txt", filename);
            std::fs::write(&exprpath, "INTCOL > 5\n").unwrap();
            let exprarg = cc(&format!("@{}", exprpath));

            fits_find_rows(
                &mut f,
                &exprarg,
                1,
                10,
                &mut n_good_rows,
                &mut row_status,
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(n_good_rows, 5); /* Rows 6,7,8,9,10 match */

            fits_close_file(f, &mut status);
            let _ = std::fs::remove_file(&exprpath);
        });
    }

    // ===================== ffsrow tests =====================

    #[test]
    fn test_ffsrow_copy_to_different_file() {
        with_temp_file(|inname| {
            with_temp_file(|outname| {
                let mut status = 0;
                let iname = to_buf(inname);
                let oname = to_buf(outname);

                let f = create_test_table(&iname);
                let mut fb = Some(f);
                fits_close_file(fb.take().unwrap(), &mut status);

                fits_open_file(&mut fb, &iname, READONLY, &mut status);
                fits_movabs_hdu(fb.as_deref_mut().unwrap(), 2, None, &mut status);

                let mut out: Option<Box<fitsfile>> = None;
                fits_create_file(&mut out, &oname, &mut status);
                fits_write_imghdr(out.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
                create_tbl(
                    out.as_deref_mut().unwrap(),
                    BINARY_TBL,
                    0,
                    &["INTCOL", "FLOATCOL", "STRCOL", "BOOLCOL"],
                    &["1J", "1E", "10A", "1L"],
                    None,
                    &mut status,
                );
                assert_eq!(status, 0, "out setup");

                fits_select_rows(
                    fb.as_deref_mut().unwrap(),
                    out.as_deref_mut().unwrap(),
                    &cc("INTCOL > 5"),
                    &mut status,
                );
                assert_eq!(status, 0);

                let mut nrows = 0;
                fits_get_num_rows(out.as_deref_mut().unwrap(), &mut nrows, &mut status);
                assert_eq!(nrows, 5);

                let mut intdata = [0 as c_long; 10];
                let mut anynull = 0;
                fits_read_col_lng(
                    out.as_deref_mut().unwrap(),
                    1,
                    1,
                    1,
                    5,
                    0,
                    &mut intdata,
                    Some(&mut anynull),
                    &mut status,
                );
                assert_eq!(status, 0);
                assert_eq!(intdata[0], 6);
                assert_eq!(intdata[4], 10);

                fits_close_file(fb.take().unwrap(), &mut status);
                fits_close_file(out.take().unwrap(), &mut status);
            });
        });
    }

    #[test]
    fn test_ffsrow_same_file_filter() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let coldata: [c_long; 10] = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            create_tbl(
                f.as_deref_mut().unwrap(),
                BINARY_TBL,
                10,
                &["INTCOL"],
                &["1J"],
                None,
                &mut status,
            );
            fits_write_col_lng(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                10,
                &coldata,
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READWRITE, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);

            let fp: *mut fitsfile = f.as_deref_mut().unwrap();
            fits_select_rows(
                unsafe { &mut *fp },
                unsafe { &mut *fp },
                &cc("INTCOL <= 5"),
                &mut status,
            );
            assert_eq!(status, 0);

            let mut nrows = 0;
            fits_get_num_rows(f.as_deref_mut().unwrap(), &mut nrows, &mut status);
            assert_eq!(nrows, 5);

            let mut intdata = [0 as c_long; 5];
            let mut anynull = 0;
            fits_read_col_lng(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                5,
                0,
                &mut intdata,
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(intdata[0], 1);
            assert_eq!(intdata[4], 5);

            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_ffsrow_constant_true() {
        with_temp_file(|inname| {
            with_temp_file(|outname| {
                let mut status = 0;
                let iname = to_buf(inname);
                let oname = to_buf(outname);
                let coldata: [c_long; 10] = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10];

                let mut infile: Option<Box<fitsfile>> = None;
                fits_create_file(&mut infile, &iname, &mut status);
                fits_write_imghdr(
                    infile.as_deref_mut().unwrap(),
                    BYTE_IMG,
                    0,
                    &[],
                    &mut status,
                );
                create_tbl(
                    infile.as_deref_mut().unwrap(),
                    BINARY_TBL,
                    10,
                    &["INTCOL"],
                    &["1J"],
                    None,
                    &mut status,
                );
                fits_write_col_lng(
                    infile.as_deref_mut().unwrap(),
                    1,
                    1,
                    1,
                    10,
                    &coldata,
                    &mut status,
                );

                let mut outfile: Option<Box<fitsfile>> = None;
                fits_create_file(&mut outfile, &oname, &mut status);
                fits_write_imghdr(
                    outfile.as_deref_mut().unwrap(),
                    BYTE_IMG,
                    0,
                    &[],
                    &mut status,
                );
                create_tbl(
                    outfile.as_deref_mut().unwrap(),
                    BINARY_TBL,
                    0,
                    &["INTCOL"],
                    &["1J"],
                    None,
                    &mut status,
                );

                fits_movabs_hdu(infile.as_deref_mut().unwrap(), 2, None, &mut status);
                fits_movabs_hdu(outfile.as_deref_mut().unwrap(), 2, None, &mut status);

                fits_select_rows(
                    infile.as_deref_mut().unwrap(),
                    outfile.as_deref_mut().unwrap(),
                    &cc("1==1"),
                    &mut status,
                );
                assert_eq!(status, 0);

                let mut nrows = 0;
                fits_get_num_rows(outfile.as_deref_mut().unwrap(), &mut nrows, &mut status);
                assert_eq!(nrows, 10);

                fits_close_file(infile.take().unwrap(), &mut status);
                fits_close_file(outfile.take().unwrap(), &mut status);
            });
        });
    }

    #[test]
    fn test_ffsrow_constant_false() {
        with_temp_file(|inname| {
            with_temp_file(|outname| {
                let mut status = 0;
                let iname = to_buf(inname);
                let oname = to_buf(outname);
                let coldata: [c_long; 10] = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10];

                let mut infile: Option<Box<fitsfile>> = None;
                fits_create_file(&mut infile, &iname, &mut status);
                fits_write_imghdr(
                    infile.as_deref_mut().unwrap(),
                    BYTE_IMG,
                    0,
                    &[],
                    &mut status,
                );
                create_tbl(
                    infile.as_deref_mut().unwrap(),
                    BINARY_TBL,
                    10,
                    &["INTCOL"],
                    &["1J"],
                    None,
                    &mut status,
                );
                fits_write_col_lng(
                    infile.as_deref_mut().unwrap(),
                    1,
                    1,
                    1,
                    10,
                    &coldata,
                    &mut status,
                );

                let mut outfile: Option<Box<fitsfile>> = None;
                fits_create_file(&mut outfile, &oname, &mut status);
                fits_write_imghdr(
                    outfile.as_deref_mut().unwrap(),
                    BYTE_IMG,
                    0,
                    &[],
                    &mut status,
                );
                create_tbl(
                    outfile.as_deref_mut().unwrap(),
                    BINARY_TBL,
                    0,
                    &["INTCOL"],
                    &["1J"],
                    None,
                    &mut status,
                );

                fits_movabs_hdu(infile.as_deref_mut().unwrap(), 2, None, &mut status);
                fits_movabs_hdu(outfile.as_deref_mut().unwrap(), 2, None, &mut status);

                fits_select_rows(
                    infile.as_deref_mut().unwrap(),
                    outfile.as_deref_mut().unwrap(),
                    &cc("1==0"),
                    &mut status,
                );
                assert_eq!(status, 0);

                let mut nrows = 0;
                fits_get_num_rows(outfile.as_deref_mut().unwrap(), &mut nrows, &mut status);
                assert_eq!(nrows, 0);

                fits_close_file(infile.take().unwrap(), &mut status);
                fits_close_file(outfile.take().unwrap(), &mut status);
            });
        });
    }

    #[test]
    fn test_ffsrow_varlen_copy() {
        with_temp_file(|inname| {
            with_temp_file(|outname| {
                let mut status = 0;
                let iname = to_buf(inname);
                let oname = to_buf(outname);
                let intdata: [c_long; 5] = [1, 2, 3, 4, 5];
                let vardata: [c_long; 5] = [10, 20, 30, 40, 50];

                let mut infile: Option<Box<fitsfile>> = None;
                fits_create_file(&mut infile, &iname, &mut status);
                fits_write_imghdr(
                    infile.as_deref_mut().unwrap(),
                    BYTE_IMG,
                    0,
                    &[],
                    &mut status,
                );
                create_tbl(
                    infile.as_deref_mut().unwrap(),
                    BINARY_TBL,
                    5,
                    &["INTCOL", "VARCOL"],
                    &["1J", "1PJ"],
                    None,
                    &mut status,
                );

                for i in 0..5 {
                    fits_write_col_lng(
                        infile.as_deref_mut().unwrap(),
                        1,
                        (i + 1) as LONGLONG,
                        1,
                        1,
                        &intdata[i..i + 1],
                        &mut status,
                    );
                    fits_write_col(
                        infile.as_deref_mut().unwrap(),
                        TLONG,
                        2,
                        (i + 1) as LONGLONG,
                        1,
                        1,
                        as_bytes_const(&vardata[i..i + 1]),
                        &mut status,
                    );
                }

                let mut outfile: Option<Box<fitsfile>> = None;
                fits_create_file(&mut outfile, &oname, &mut status);
                fits_write_imghdr(
                    outfile.as_deref_mut().unwrap(),
                    BYTE_IMG,
                    0,
                    &[],
                    &mut status,
                );
                create_tbl(
                    outfile.as_deref_mut().unwrap(),
                    BINARY_TBL,
                    0,
                    &["INTCOL", "VARCOL"],
                    &["1J", "1PJ"],
                    None,
                    &mut status,
                );
                assert_eq!(status, 0, "varlen setup");

                fits_select_rows(
                    infile.as_deref_mut().unwrap(),
                    outfile.as_deref_mut().unwrap(),
                    &cc("INTCOL > 2"),
                    &mut status,
                );
                assert_eq!(status, 0);

                let mut nrows = 0;
                fits_get_num_rows(outfile.as_deref_mut().unwrap(), &mut nrows, &mut status);
                assert_eq!(nrows, 3);

                fits_close_file(infile.take().unwrap(), &mut status);
                fits_close_file(outfile.take().unwrap(), &mut status);
            });
        });
    }

    #[test]
    fn test_ffsrow_varlen_heap_copy() {
        with_temp_file(|inname| {
            with_temp_file(|outname| {
                let mut status = 0;
                let iname = to_buf(inname);
                let oname = to_buf(outname);
                let intdata: [c_long; 5] = [1, 2, 3, 4, 5];
                let mut vardata = [0 as c_long; 100];
                for (j, v) in vardata.iter_mut().enumerate() {
                    *v = (j * 10) as c_long;
                }

                let mut infile: Option<Box<fitsfile>> = None;
                fits_create_file(&mut infile, &iname, &mut status);
                fits_write_imghdr(
                    infile.as_deref_mut().unwrap(),
                    BYTE_IMG,
                    0,
                    &[],
                    &mut status,
                );
                create_tbl(
                    infile.as_deref_mut().unwrap(),
                    BINARY_TBL,
                    5,
                    &["INTCOL", "VARCOL"],
                    &["1J", "PJ(100)"],
                    None,
                    &mut status,
                );

                for i in 0..5 {
                    fits_write_col_lng(
                        infile.as_deref_mut().unwrap(),
                        1,
                        (i + 1) as LONGLONG,
                        1,
                        1,
                        &intdata[i..i + 1],
                        &mut status,
                    );
                    fits_write_col(
                        infile.as_deref_mut().unwrap(),
                        TLONG,
                        2,
                        (i + 1) as LONGLONG,
                        1,
                        100,
                        as_bytes_const(&vardata),
                        &mut status,
                    );
                }

                let mut outfile: Option<Box<fitsfile>> = None;
                fits_create_file(&mut outfile, &oname, &mut status);
                fits_write_imghdr(
                    outfile.as_deref_mut().unwrap(),
                    BYTE_IMG,
                    0,
                    &[],
                    &mut status,
                );
                create_tbl(
                    outfile.as_deref_mut().unwrap(),
                    BINARY_TBL,
                    0,
                    &["INTCOL", "VARCOL"],
                    &["1J", "PJ(100)"],
                    None,
                    &mut status,
                );
                assert_eq!(status, 0, "heap setup");

                fits_select_rows(
                    infile.as_deref_mut().unwrap(),
                    outfile.as_deref_mut().unwrap(),
                    &cc("INTCOL > 2"),
                    &mut status,
                );
                assert_eq!(status, 0);

                let mut nrows = 0;
                fits_get_num_rows(outfile.as_deref_mut().unwrap(), &mut nrows, &mut status);
                assert_eq!(nrows, 3);

                fits_close_file(infile.take().unwrap(), &mut status);
                fits_close_file(outfile.take().unwrap(), &mut status);
            });
        });
    }

    #[test]
    fn test_ffsrow_nonboolean_error() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));

            let fp: *mut fitsfile = &raw mut *f;
            fits_select_rows(
                unsafe { &mut *fp },
                unsafe { &mut *fp },
                &cc("INTCOL + 5"),
                &mut status,
            );
            assert_eq!(status, PARSE_BAD_TYPE);

            status = 0;
            fits_close_file(f, &mut status);
        });
    }

    // ===================== ASCII table tests =====================

    #[test]
    fn test_ffcalc_ascii_int_column() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_ascii_table(&to_buf(filename));
            let mut results = [0 as c_long; 10];
            let mut anynul = 0;

            let fp_self: *mut fitsfile = &raw mut *f;
            fits_calculator(
                unsafe { &mut *fp_self },
                &cc("INTCOL * 2"),
                unsafe { &mut *fp_self },
                &cc("DOUBLED"),
                &cc(""),
                &mut status,
            );
            assert_eq!(status, 0);

            fits_read_col_lng(
                &mut f,
                4,
                1,
                1,
                10,
                0,
                &mut results,
                Some(&mut anynul),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(results[0], 2);
            assert_eq!(results[4], 10);
            assert_eq!(results[9], 20);

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_ffcalc_ascii_double_column() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_ascii_table(&to_buf(filename));
            let mut results = [0.0f64; 10];
            let mut anynul = 0;

            let fp_self: *mut fitsfile = &raw mut *f;
            fits_calculator(
                unsafe { &mut *fp_self },
                &cc("FLOATCOL + 0.5"),
                unsafe { &mut *fp_self },
                &cc("SHIFTED"),
                &cc(""),
                &mut status,
            );
            assert_eq!(status, 0);

            fits_read_col_dbl(
                &mut f,
                4,
                1,
                1,
                10,
                0.0,
                &mut results,
                Some(&mut anynul),
                &mut status,
            );
            assert_eq!(status, 0);
            assert!(results[0] >= 1.9 && results[0] <= 2.1); /* 1.5 + 0.5 */
            assert!(results[4] >= 5.9 && results[4] <= 6.1); /* 5.5 + 0.5 */

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_ffcalc_ascii_string_column() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_ascii_table(&to_buf(filename));

            let fp_self: *mut fitsfile = &raw mut *f;
            fits_calculator(
                unsafe { &mut *fp_self },
                &cc("'test'"),
                unsafe { &mut *fp_self },
                &cc("NEWSTR"),
                &cc(""),
                &mut status,
            );
            assert_eq!(status, 0);

            let mut buffer = [0 as c_char; 20];
            let mut anynull = 0;
            {
                let mut arr: Vec<&mut [c_char]> = vec![&mut buffer[..]];
                fits_read_col_str(
                    &mut f,
                    4,
                    1,
                    1,
                    1,
                    None,
                    &mut arr,
                    Some(&mut anynull),
                    &mut status,
                );
            }
            assert_eq!(status, 0);
            assert_eq!(&crate::helpers::testhelpers::from_buf(&buffer)[..4], "test");

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_ffcalc_ascii_logical_error() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_ascii_table(&to_buf(filename));

            let fp_self: *mut fitsfile = &raw mut *f;
            fits_calculator(
                unsafe { &mut *fp_self },
                &cc("INTCOL > 5"),
                unsafe { &mut *fp_self },
                &cc("LOGCOL"),
                &cc(""),
                &mut status,
            );
            assert_eq!(status, NOT_BTABLE);

            status = 0;
            fits_close_file(f, &mut status);
        });
    }

    // ===================== overflow / div-by-zero / deref tests =====================

    #[test]
    #[ignore = "Long long not supported on Windows"]
    fn test_ffcrow_long_div_overflow() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let minval: c_long = i32::MIN as c_long;
            let negone: c_long = -1;

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            create_tbl(
                f.as_deref_mut().unwrap(),
                BINARY_TBL,
                1,
                &["MINCOL", "NEGCOL"],
                &["1J", "1J"],
                None,
                &mut status,
            );
            fits_write_col_lng(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                1,
                &[minval],
                &mut status,
            );
            fits_write_col_lng(
                f.as_deref_mut().unwrap(),
                2,
                1,
                1,
                1,
                &[negone],
                &mut status,
            );

            let mut div_result = [0 as c_long; 1];
            let mut mod_result = [0 as c_long; 1];
            let mut anynul = 0;

            fits_calc_rows(
                f.as_deref_mut().unwrap(),
                TLONG,
                &cc("MINCOL / NEGCOL"),
                1,
                1,
                core::ptr::null(),
                as_bytes_mut(&mut div_result),
                &mut anynul,
                &mut status,
            );
            assert_eq!(status, 0);
            if core::mem::size_of::<c_long>() == 4 {
                assert_eq!(div_result[0], c_long::MAX);
            } else {
                assert_eq!(div_result[0], -(i32::MIN as i64) as c_long);
            }

            fits_calc_rows(
                f.as_deref_mut().unwrap(),
                TLONG,
                &cc("MINCOL % NEGCOL"),
                1,
                1,
                core::ptr::null(),
                as_bytes_mut(&mut mod_result),
                &mut anynul,
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(mod_result[0], 0);

            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_ffcrow_long_div_by_zero() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            create_tbl(
                f.as_deref_mut().unwrap(),
                BINARY_TBL,
                1,
                &["NUMCOL", "DENCOL"],
                &["1J", "1J"],
                None,
                &mut status,
            );
            fits_write_col_lng(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                1,
                &[42 as c_long],
                &mut status,
            );
            fits_write_col_lng(
                f.as_deref_mut().unwrap(),
                2,
                1,
                1,
                1,
                &[0 as c_long],
                &mut status,
            );

            let mut div_result = [0 as c_long; 1];
            let mut mod_result = [0 as c_long; 1];

            let mut anynul = 0;
            fits_calc_rows(
                f.as_deref_mut().unwrap(),
                TLONG,
                &cc("NUMCOL / DENCOL"),
                1,
                1,
                core::ptr::null(),
                as_bytes_mut(&mut div_result),
                &mut anynul,
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(anynul, 1);

            let mut anynul = 0;
            fits_calc_rows(
                f.as_deref_mut().unwrap(),
                TLONG,
                &cc("NUMCOL % DENCOL"),
                1,
                1,
                core::ptr::null(),
                as_bytes_mut(&mut mod_result),
                &mut anynul,
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(anynul, 1);

            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_ffcrow_vector_index_out_of_range() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let vecdata: [c_long; 5] = [10, 20, 30, 40, 50];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            create_tbl(
                f.as_deref_mut().unwrap(),
                BINARY_TBL,
                1,
                &["VECCOL"],
                &["5J"],
                None,
                &mut status,
            );
            fits_write_col_lng(f.as_deref_mut().unwrap(), 1, 1, 1, 5, &vecdata, &mut status);

            let fp: *mut fitsfile = f.as_deref_mut().unwrap();
            fits_calculator(
                unsafe { &mut *fp },
                &cc("VECCOL[99]"),
                unsafe { &mut *fp },
                &cc("BADCOL"),
                &cc(""),
                &mut status,
            );
            assert_ne!(status, 0);

            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_ffcrow_vector_index_zero() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let vecdata: [c_long; 5] = [10, 20, 30, 40, 50];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            create_tbl(
                f.as_deref_mut().unwrap(),
                BINARY_TBL,
                1,
                &["VECCOL"],
                &["5J"],
                None,
                &mut status,
            );
            fits_write_col_lng(f.as_deref_mut().unwrap(), 1, 1, 1, 5, &vecdata, &mut status);

            let fp: *mut fitsfile = f.as_deref_mut().unwrap();
            fits_calculator(
                unsafe { &mut *fp },
                &cc("VECCOL[0]"),
                unsafe { &mut *fp },
                &cc("BADCOL"),
                &cc(""),
                &mut status,
            );
            assert_ne!(status, 0);

            status = 0;
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_ffcrow_double_column_arithmetic() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let v1: f64 = 1.0e308;
            let v2: f64 = 2.0;

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            create_tbl(
                f.as_deref_mut().unwrap(),
                BINARY_TBL,
                1,
                &["DCOL1", "DCOL2"],
                &["1D", "1D"],
                None,
                &mut status,
            );
            fits_write_col_dbl(f.as_deref_mut().unwrap(), 1, 1, 1, 1, &[v1], &mut status);
            fits_write_col_dbl(f.as_deref_mut().unwrap(), 2, 1, 1, 1, &[v2], &mut status);

            let mut result = [0.0f64; 1];
            let mut anynul = 0;
            fits_calc_rows(
                f.as_deref_mut().unwrap(),
                TDOUBLE,
                &cc("DCOL1 / DCOL2"),
                1,
                1,
                core::ptr::null(),
                as_bytes_mut(&mut result),
                &mut anynul,
                &mut status,
            );
            assert_eq!(status, 0);
            assert!(result[0] >= 4.9e307 && result[0] <= 5.1e307);

            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    // ============ within-range ("in-range") operator: x = a:b ============
    //
    // "In the within-range expression, "x=a:b", the value of x is tested to
    //  be within the range a through b, where all of those values may be
    //  columns or expressions.  The result value is a boolean true or false.
    //  The expression is a shorthand notation equivalent to the expression
    //  "((a.le.x).and.(x.le.b))"."

    #[test]
    fn test_fffrow_within_range() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));
            let mut n_good_rows = 0;
            let mut row_status = [0 as c_char; 10];

            fits_find_rows(
                &mut f,
                &cc("(INTCOL = 3:7)"),
                1,
                10,
                &mut n_good_rows,
                &mut row_status,
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(n_good_rows, 5); /* rows 3,4,5,6,7 */
            assert_eq!(row_status[1], 0); /* 2 is below the range */
            assert_eq!(row_status[2], 1); /* lower bound is inclusive */
            assert_eq!(row_status[6], 1); /* upper bound is inclusive */
            assert_eq!(row_status[7], 0); /* 8 is above the range */

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_fffrow_within_range_no_parens() {
        /* the docs recommend parentheses to avoid confusing other parsers,
        but the bare form is legal grammar as well */
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));
            let mut n_good_rows = 0;
            let mut row_status = [0 as c_char; 10];

            fits_find_rows(
                &mut f,
                &cc("INTCOL = 3:7"),
                1,
                10,
                &mut n_good_rows,
                &mut row_status,
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(n_good_rows, 5); /* rows 3,4,5,6,7 */

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_fffrow_within_range_equivalent_expression() {
        /* "(x=a:b)" is shorthand for "((a.le.x).and.(x.le.b))" - verify the
        two forms select exactly the same rows */
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));
            let mut n_short = 0;
            let mut short_status = [0 as c_char; 10];
            let mut n_long = 0;
            let mut long_status = [0 as c_char; 10];

            fits_find_rows(
                &mut f,
                &cc("(FLOATCOL = 2:8)"),
                1,
                10,
                &mut n_short,
                &mut short_status,
                &mut status,
            );
            assert_eq!(status, 0);
            fits_find_rows(
                &mut f,
                &cc("((2 .le. FLOATCOL) .and. (FLOATCOL .le. 8))"),
                1,
                10,
                &mut n_long,
                &mut long_status,
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(n_short, n_long);
            assert_eq!(short_status, long_status);
            assert_eq!(n_short, 6); /* 2.5,3.5,4.5,5.5,6.5,7.5 */

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_fffrow_within_range_reversed() {
        /* an inverted range can never be satisfied: a.le.x .and. x.le.b */
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));
            let mut n_good_rows = 0;
            let mut row_status = [0 as c_char; 10];

            fits_find_rows(
                &mut f,
                &cc("(INTCOL = 7:3)"),
                1,
                10,
                &mut n_good_rows,
                &mut row_status,
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(n_good_rows, 0);

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_fffrow_within_range_promotes_types() {
        /* integer column tested against a floating point range */
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));
            let mut n_good_rows = 0;
            let mut row_status = [0 as c_char; 10];

            fits_find_rows(
                &mut f,
                &cc("(INTCOL = 2.5:5.5)"),
                1,
                10,
                &mut n_good_rows,
                &mut row_status,
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(n_good_rows, 3); /* rows 3,4,5 */
            assert_eq!(row_status[1], 0);
            assert_eq!(row_status[2], 1);
            assert_eq!(row_status[4], 1);
            assert_eq!(row_status[5], 0);

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_fffrow_within_range_column_bounds() {
        /* x, a and b may all be columns or expressions:
        FLOATCOL is always INTCOL+0.5, so it lies in [INTCOL, INTCOL+1] */
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));
            let mut n_good_rows = 0;
            let mut row_status = [0 as c_char; 10];

            fits_find_rows(
                &mut f,
                &cc("(FLOATCOL = INTCOL:INTCOL+1)"),
                1,
                10,
                &mut n_good_rows,
                &mut row_status,
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(n_good_rows, 10);

            /* ... while INTCOL is always just below FLOATCOL */
            fits_find_rows(
                &mut f,
                &cc("(INTCOL = FLOATCOL:100)"),
                1,
                10,
                &mut n_good_rows,
                &mut row_status,
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(n_good_rows, 0);

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_fffrow_within_range_compound() {
        /* the result is an ordinary boolean, usable in a larger expression */
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));
            let mut n_good_rows = 0;
            let mut row_status = [0 as c_char; 10];

            fits_find_rows(
                &mut f,
                &cc("(INTCOL = 3:7) && BOOLCOL"),
                1,
                10,
                &mut n_good_rows,
                &mut row_status,
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(n_good_rows, 3); /* rows 3,5,7 - BOOLCOL is T,F,T,F,... */
            assert_eq!(row_status[2], 1);
            assert_eq!(row_status[3], 0);
            assert_eq!(row_status[4], 1);

            /* and it can be negated */
            fits_find_rows(
                &mut f,
                &cc("!(INTCOL = 3:7)"),
                1,
                10,
                &mut n_good_rows,
                &mut row_status,
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(n_good_rows, 5); /* rows 1,2,8,9,10 */

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_fffrow_within_range_row_number() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));
            let mut n_good_rows = 0;
            let mut row_status = [0 as c_char; 10];

            fits_find_rows(
                &mut f,
                &cc("(#row = 4:6)"),
                1,
                10,
                &mut n_good_rows,
                &mut row_status,
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(n_good_rows, 3);
            assert_eq!(row_status[3], 1);
            assert_eq!(row_status[5], 1);
            assert_eq!(row_status[6], 0);

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_fftexp_within_range() {
        /* "The result value is a boolean true or false." */
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));
            let mut datatype = 0;
            let mut naxis = 0;
            let mut nelem = 0;
            let mut naxes = [0 as c_long; 10];

            fits_test_expr(
                &mut f,
                &cc("(INTCOL = 3:7)"),
                10,
                &mut datatype,
                &mut nelem,
                &mut naxis,
                &mut naxes,
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(datatype, TLOGICAL);
            assert_eq!(nelem, 1);

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_ffcrow_within_range() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));
            let mut results = [0 as c_char; 10];
            let mut anynul = 0;

            fits_calc_rows(
                &mut f,
                TLOGICAL,
                &cc("(INTCOL = 4:6)"),
                1,
                10,
                core::ptr::null(),
                as_bytes_mut(&mut results),
                &mut anynul,
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(results, [0, 0, 0, 1, 1, 1, 0, 0, 0, 0]);

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_ffcrow_within_range_vector() {
        /* a vector x produces a vector boolean result */
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let vecdata: [f64; 5] = [1.0, 2.0, 3.0, 4.0, 5.0];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            create_tbl(
                f.as_deref_mut().unwrap(),
                BINARY_TBL,
                1,
                &["VECCOL"],
                &["5D"],
                None,
                &mut status,
            );
            fits_write_col_dbl(f.as_deref_mut().unwrap(), 1, 1, 1, 5, &vecdata, &mut status);
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);

            let mut datatype = 0;
            let mut naxis = 0;
            let mut nelem = 0;
            let mut naxes = [0 as c_long; 10];
            fits_test_expr(
                f.as_deref_mut().unwrap(),
                &cc("(VECCOL = 2:4)"),
                10,
                &mut datatype,
                &mut nelem,
                &mut naxis,
                &mut naxes,
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(datatype, TLOGICAL);
            assert_eq!(nelem, 5);

            let mut results = [0 as c_char; 5];
            let mut anynul = 0;
            fits_calc_rows(
                f.as_deref_mut().unwrap(),
                TLOGICAL,
                &cc("(VECCOL = 2:4)"),
                1,
                5,
                core::ptr::null(),
                as_bytes_mut(&mut results),
                &mut anynul,
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(results, [0, 1, 1, 1, 0]);

            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_ffcalc_within_range() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));
            let mut results = [0 as c_char; 10];
            let mut nullarr = [0 as c_char; 10];
            let mut anynul = 0;

            let fp_self: *mut fitsfile = &raw mut *f;
            fits_calculator(
                unsafe { &mut *fp_self },
                &cc("(INTCOL = 3:7)"),
                unsafe { &mut *fp_self },
                &cc("INRANGE"),
                &cc("1L"),
                &mut status,
            );
            assert_eq!(status, 0);

            fits_read_colnull_log(
                &mut f,
                5,
                1,
                1,
                10,
                &mut results,
                &mut nullarr,
                Some(&mut anynul),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(results, [0, 0, 1, 1, 1, 1, 1, 0, 0, 0]);

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_ffsrow_within_range() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            {
                let f = create_test_table(&name);
                fits_close_file(f, &mut status);
            }

            let mut f: Option<Box<fitsfile>> = None;
            fits_open_file(&mut f, &name, READWRITE, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);

            let fp: *mut fitsfile = f.as_deref_mut().unwrap();
            fits_select_rows(
                unsafe { &mut *fp },
                unsafe { &mut *fp },
                &cc("(INTCOL = 3:7)"),
                &mut status,
            );
            assert_eq!(status, 0);

            let mut nrows = 0;
            fits_get_num_rows(f.as_deref_mut().unwrap(), &mut nrows, &mut status);
            assert_eq!(nrows, 5);

            let mut intdata = [0 as c_long; 5];
            let mut anynull = 0;
            fits_read_col_lng(
                f.as_deref_mut().unwrap(),
                1,
                1,
                1,
                5,
                0,
                &mut intdata,
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(intdata, [3, 4, 5, 6, 7]);

            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    // ===================== operator / function coverage =====================
    //
    // The expression syntax accepted by fffrow/ffcrow/ffcalc/ffsrow is
    // documented in the "Row Filtering Specification" section of the CFITSIO
    // manual.  The tests below walk through that operator and function table.

    /// Evaluate a boolean expression over all 10 rows of `create_test_table`,
    /// returning (number of true rows, per-row flags).
    fn count_rows(f: &mut fitsfile, expr: &str) -> (c_long, [c_char; 10]) {
        let mut status = 0;
        let mut n_good_rows: c_long = 0;
        let mut row_status = [0 as c_char; 10];
        fits_find_rows(
            f,
            &cc(expr),
            1,
            10,
            &mut n_good_rows,
            &mut row_status,
            &mut status,
        );
        assert_eq!(status, 0, "fits_find_rows failed for expression: {expr}");
        (n_good_rows, row_status)
    }

    /// Evaluate an expression, returning N elements as doubles.
    fn eval_dbl<const N: usize>(f: &mut fitsfile, expr: &str) -> [f64; N] {
        let mut status = 0;
        let mut results = [0.0f64; N];
        let mut anynul = 0;
        fits_calc_rows(
            f,
            TDOUBLE,
            &cc(expr),
            1,
            N as c_long,
            core::ptr::null(),
            as_bytes_mut(&mut results),
            &mut anynul,
            &mut status,
        );
        assert_eq!(status, 0, "fits_calc_rows failed for expression: {expr}");
        results
    }

    /// Evaluate an expression, returning N elements as longs.
    fn eval_lng<const N: usize>(f: &mut fitsfile, expr: &str) -> [c_long; N] {
        let mut status = 0;
        let mut results = [0 as c_long; N];
        let mut anynul = 0;
        fits_calc_rows(
            f,
            TLONG,
            &cc(expr),
            1,
            N as c_long,
            core::ptr::null(),
            as_bytes_mut(&mut results),
            &mut anynul,
            &mut status,
        );
        assert_eq!(status, 0, "fits_calc_rows failed for expression: {expr}");
        results
    }

    /// Evaluate an expression, returning N elements as logicals.
    fn eval_log<const N: usize>(f: &mut fitsfile, expr: &str) -> [c_char; N] {
        let mut status = 0;
        let mut results = [0 as c_char; N];
        let mut anynul = 0;
        fits_calc_rows(
            f,
            TLOGICAL,
            &cc(expr),
            1,
            N as c_long,
            core::ptr::null(),
            as_bytes_mut(&mut results),
            &mut anynul,
            &mut status,
        );
        assert_eq!(status, 0, "fits_calc_rows failed for expression: {expr}");
        results
    }

    #[test]
    fn test_fffrow_fortran_boolean_spellings() {
        /* "Boolean operators can be used in the expression in either their
        Fortran or C forms." */
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));

            assert_eq!(count_rows(&mut f, "INTCOL .gt. 5").0, 5);
            assert_eq!(count_rows(&mut f, "INTCOL .GT. 5").0, 5);
            assert_eq!(count_rows(&mut f, "INTCOL gt. 5").0, 5);
            assert_eq!(count_rows(&mut f, "INTCOL .lt. 5").0, 4);
            assert_eq!(count_rows(&mut f, "INTCOL .LT. 5").0, 4);
            assert_eq!(count_rows(&mut f, "INTCOL .ge. 5").0, 6);
            assert_eq!(count_rows(&mut f, "INTCOL .le. 5").0, 5);
            assert_eq!(count_rows(&mut f, "INTCOL .eq. 5").0, 1);
            assert_eq!(count_rows(&mut f, "INTCOL .EQ. 5").0, 1);
            assert_eq!(count_rows(&mut f, "INTCOL .ne. 5").0, 9);

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_fffrow_fortran_logical_spellings() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));

            assert_eq!(
                count_rows(&mut f, "(INTCOL .gt. 3) .and. (INTCOL .lt. 6)").0,
                2
            );
            assert_eq!(
                count_rows(&mut f, "(INTCOL .lt. 3) .or. (INTCOL .gt. 8)").0,
                4
            );
            assert_eq!(count_rows(&mut f, ".not. (INTCOL .gt. 5)").0, 5);
            assert_eq!(count_rows(&mut f, "NOT. (INTCOL .GT. 5)").0, 5);
            assert_eq!(count_rows(&mut f, "(INTCOL .GT. 3) .AND. BOOLCOL").0, 3);

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_fffrow_alternate_comparison_spellings() {
        /* "=>" and "=<" are accepted as synonyms for ">=" and "<=" */
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));

            assert_eq!(count_rows(&mut f, "INTCOL => 5").0, 6);
            assert_eq!(count_rows(&mut f, "INTCOL =< 5").0, 5);

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_ffcrow_power_caret() {
        /* "exponentiation"  ** ^ */
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));

            let r = eval_dbl::<4>(&mut f, "INTCOL ^ 2");
            assert!((r[0] - 1.0).abs() < 1e-9);
            assert!((r[3] - 16.0).abs() < 1e-9);

            let r = eval_dbl::<1>(&mut f, "2 ** 3");
            assert!((r[0] - 8.0).abs() < 1e-9);

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_ffcrow_bitwise_xor() {
        /* "bitwise XOR"  x ^^ y  (32-bit int only) */
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));

            let r = eval_lng::<4>(&mut f, "INTCOL ^^ 3");
            assert_eq!(r, [2, 1, 0, 7]);

            let r = eval_lng::<4>(&mut f, "INTCOL .xor. 3");
            assert_eq!(r, [2, 1, 0, 7]);

            let r = eval_lng::<4>(&mut f, "INTCOL .XOR. 3");
            assert_eq!(r, [2, 1, 0, 7]);

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_ffcrow_type_casts() {
        /* "real to integer"  (int) x     "integer to real"  (float) i */
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));

            /* FLOATCOL is 1.5, 2.5, ... - the cast truncates toward zero */
            let r = eval_lng::<4>(&mut f, "(int) FLOATCOL");
            assert_eq!(r, [1, 2, 3, 4]);
            let r = eval_lng::<4>(&mut f, "(INT) FLOATCOL");
            assert_eq!(r, [1, 2, 3, 4]);
            let r = eval_lng::<1>(&mut f, "(int) -2.7");
            assert_eq!(r, [-2]);

            /* the integer to real cast promotes to double precision */
            let r = eval_dbl::<4>(&mut f, "(float) INTCOL + 0.25");
            assert!((r[0] - 1.25).abs() < 1e-9);
            let r = eval_dbl::<1>(&mut f, "(FLOAT) 7");
            assert!((r[0] - 7.0).abs() < 1e-9);
            let r = eval_dbl::<1>(&mut f, "(double) 7");
            assert!((r[0] - 7.0).abs() < 1e-9);

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_ffcrow_builtin_constants() {
        /* #pi, #e and #deg (#row is covered by test_ffcrow_row_number) */
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));

            let r = eval_dbl::<1>(&mut f, "#pi");
            assert!((r[0] - core::f64::consts::PI).abs() < 1e-12);
            let r = eval_dbl::<1>(&mut f, "#PI");
            assert!((r[0] - core::f64::consts::PI).abs() < 1e-12);
            let r = eval_dbl::<1>(&mut f, "#e");
            assert!((r[0] - core::f64::consts::E).abs() < 1e-12);
            let r = eval_dbl::<1>(&mut f, "#deg");
            assert!((r[0] - core::f64::consts::PI / 180.0).abs() < 1e-12);

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_ffcrow_integer_constant_bases() {
        /* 13245 decimal, 0x12f3 hex, 0o1373 octal, 0b01001 binary */
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));

            assert_eq!(eval_lng::<1>(&mut f, "13245"), [13245]);
            assert_eq!(eval_lng::<1>(&mut f, "0x12f3"), [0x12f3]);
            assert_eq!(eval_lng::<1>(&mut f, "0o1373"), [0o1373]);
            assert_eq!(eval_lng::<1>(&mut f, "0b01001"), [0b01001]);
            /* usable in ordinary arithmetic */
            assert_eq!(eval_lng::<3>(&mut f, "INTCOL + 0x10"), [17, 18, 19]);

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_fftexp_rejects_mismatched_operand_types() {
        /* An operator applied across incompatible operand sorts is a syntax
        error: the grammar has no production for it, so bison reports one and
        recovers through `line: error '\n'`.

        The recovery is what made this fragile. Resolving a column name during
        recovery runs find_column, which calls fits_parser_allocateCol; that
        used to be handed a fresh `tstatus = 0` whose value was written back
        over lParse.status, discarding the syntax error and leaving the parse
        looking successful. fits_test_expr then returned 0 and handed back the
        *left* operand -- BOOLCOL + INTCOL evaluated to BOOLCOL. */
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));

            for expr in [
                "BOOLCOL + INTCOL",
                "BOOLCOL - FLOATCOL",
                "BOOLCOL && INTCOL",
                "STRCOL - INTCOL",
                "STRCOL > INTCOL",
            ] {
                let mut datatype = 0;
                let mut nelem: c_long = 0;
                let mut naxis = 0;
                let mut naxes = [0 as c_long; 5];
                let mut st = 0;
                fits_test_expr(
                    &mut f,
                    &cc(expr),
                    5,
                    &mut datatype,
                    &mut nelem,
                    &mut naxis,
                    &mut naxes,
                    &mut st,
                );
                assert_eq!(st, PARSE_SYNTAX_ERR, "expr should be rejected: {expr}");
                fits_clear_errmsg();
            }

            /* the compatible forms still parse */
            for expr in [
                "INTCOL + FLOATCOL",
                "BOOLCOL && BOOLCOL",
                "INTCOL > FLOATCOL",
            ] {
                let mut datatype = 0;
                let mut nelem: c_long = 0;
                let mut naxis = 0;
                let mut naxes = [0 as c_long; 5];
                let mut st = 0;
                fits_test_expr(
                    &mut f,
                    &cc(expr),
                    5,
                    &mut datatype,
                    &mut nelem,
                    &mut naxis,
                    &mut naxes,
                    &mut st,
                );
                assert_eq!(st, 0, "expr should parse: {expr}");
            }

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_fffrow_string_keyword_reference() {
        /* '#KEY' where the keyword's value is a quoted string resolves to a
        string constant.  ffdtyp tested cval[0] against a backslash instead of
        a single quote, so quoted values were classified as integers and
        find_keywd reported BAD_C2I. */
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));
            fits_write_key_str(&mut f, &cc("STRKEY"), &cc("gamma"), None, &mut status);
            assert_eq!(status, 0);

            let mut datatype = 0;
            let mut nelem: c_long = 0;
            let mut naxis = 0;
            let mut naxes = [0 as c_long; 5];
            fits_test_expr(
                &mut f,
                &cc("#STRKEY"),
                5,
                &mut datatype,
                &mut nelem,
                &mut naxis,
                &mut naxes,
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(datatype, TSTRING);

            /* and it is usable as a string operand */
            assert_eq!(count_rows(&mut f, "STRCOL == #STRKEY").0, 1);

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_fffrow_empty_variable_name() {
        /* '$$' lexes as a variable with an empty name.  The lookup walks the
        header with namelen == 0, which used to underflow `namelen - 1` in
        ffgcrd_safe and panic; it must simply fail to find the column. */
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));
            let mut n_good_rows = 0;
            let mut row_status = [0 as c_char; 10];

            for expr in ["$$", "$$ + 1", "$$ > 1"] {
                status = 0;
                fits_find_rows(
                    &mut f,
                    &cc(expr),
                    1,
                    10,
                    &mut n_good_rows,
                    &mut row_status,
                    &mut status,
                );
                assert_ne!(status, 0, "expr should be rejected: {expr}");
            }

            status = 0;
            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_ffcrow_hex_constant_case() {
        /* 0x literals must be case-insensitive.  The lexer used to convert
        non-digits with (*p - 'a' + 10) without folding case, so uppercase
        A-F produced negative digit values and 0x1F evaluated to -1. */
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));

            for (expr, want) in [
                ("0x1f", 31),
                ("0x1F", 31),
                ("0xff", 255),
                ("0xFF", 255),
                ("0xFf", 255),
                ("0xa", 10),
                ("0xA", 10),
                ("0xabcdef", 11259375),
                ("0xABCDEF", 11259375),
                ("0x10", 16),
                ("0x0", 0),
            ] {
                assert_eq!(eval_lng::<1>(&mut f, expr), [want], "expr: {expr}");
            }

            /* the other radix prefixes must keep working */
            assert_eq!(eval_lng::<1>(&mut f, "0b1011"), [11]);
            assert_eq!(eval_lng::<1>(&mut f, "0o17"), [15]);

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_fftexp_bare_dot_rejected() {
        /* A lone '.' must not lex as the double 0.0.  The {real} pattern's
        third alternative used to allow zero digits before the point. */
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));

            for expr in [".", "..", "INTCOL + .", ". + 1", ".e5"] {
                let mut datatype = 0;
                let mut nelem: c_long = 0;
                let mut naxis = 0;
                let mut naxes = [0 as c_long; 5];
                let mut st = 0;
                fits_test_expr(
                    &mut f,
                    &cc(expr),
                    5,
                    &mut datatype,
                    &mut nelem,
                    &mut naxis,
                    &mut naxes,
                    &mut st,
                );
                assert_eq!(st, PARSE_SYNTAX_ERR, "expr should be rejected: {expr}");
            }

            /* valid reals are unaffected */
            for (expr, want) in [
                ("1.", 1.0),
                ("12.", 12.0),
                (".5", 0.5),
                ("1.5", 1.5),
                ("1.5e3", 1500.0),
                ("1.5E3", 1500.0),
                ("0.", 0.0),
            ] {
                assert_eq!(eval_dbl::<1>(&mut f, expr), [want], "expr: {expr}");
            }

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_fffrow_isnull_and_setnull() {
        /* "declare certain value null" SETNULL(x,y) - if x==y a NULL is
        returned, otherwise y;  "a null value?" ISNULL(x) */
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));

            assert_eq!(count_rows(&mut f, "ISNULL(INTCOL)").0, 0);

            let (n, flags) = count_rows(&mut f, "ISNULL(SETNULL(3,INTCOL))");
            assert_eq!(n, 1);
            assert_eq!(flags[2], 1);

            /* DEFNULL substitutes for the value made NULL by SETNULL */
            let r = eval_lng::<5>(&mut f, "DEFNULL(SETNULL(3,INTCOL), -99)");
            assert_eq!(r, [1, 2, -99, 4, 5]);

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_fffrow_null_constant() {
        /* "#null  undefined value" - useful for conditionally setting values
        to a NULL, eg. "col1==-99 ? #NULL : col1" */
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));

            let (n, flags) = count_rows(&mut f, "ISNULL(INTCOL==3 ? #NULL : INTCOL)");
            assert_eq!(n, 1);
            assert_eq!(flags[2], 1);

            let r = eval_lng::<4>(&mut f, "DEFNULL(INTCOL==3 ? #null : INTCOL, -1)");
            assert_eq!(r, [1, 2, -1, 4]);

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_fffrow_snull_constant() {
        /* "#snull  undefined string" */
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));

            let (n, flags) = count_rows(&mut f, "ISNULL(INTCOL==2 ? #snull : STRCOL)");
            assert_eq!(n, 1);
            assert_eq!(flags[1], 1);

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_fffrow_near_function() {
        /* near(value_1, value_2, tolerance) */
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));

            let (n, flags) = count_rows(&mut f, "NEAR(FLOATCOL, 3.5, 0.01)");
            assert_eq!(n, 1);
            assert_eq!(flags[2], 1);

            /* the tolerance test is strict: |x-y| < tolerance, so a tolerance
            of exactly 1.0 still excludes 2.5 and 4.5 */
            assert_eq!(count_rows(&mut f, "NEAR(FLOATCOL, 3.5, 1.0)").0, 1);
            assert_eq!(count_rows(&mut f, "NEAR(FLOATCOL, 3.5, 1.001)").0, 3);
            assert_eq!(count_rows(&mut f, "near(INTCOL, 5, 0.5)").0, 1);
            assert_eq!(count_rows(&mut f, "near(INTCOL, 5, 0)").0, 0);

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_ffcrow_hyperbolic_functions() {
        /* SINH(x), COSH(x), TANH(x) */
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));

            let r = eval_dbl::<1>(&mut f, "SINH(0)");
            assert!(r[0].abs() < 1e-12);
            let r = eval_dbl::<1>(&mut f, "COSH(0)");
            assert!((r[0] - 1.0).abs() < 1e-12);
            let r = eval_dbl::<1>(&mut f, "TANH(0)");
            assert!(r[0].abs() < 1e-12);
            let r = eval_dbl::<1>(&mut f, "SINH(1)");
            assert!((r[0] - 1.0_f64.sinh()).abs() < 1e-12);
            let r = eval_dbl::<1>(&mut f, "COSH(1)");
            assert!((r[0] - 1.0_f64.cosh()).abs() < 1e-12);
            let r = eval_dbl::<1>(&mut f, "TANH(1)");
            assert!((r[0] - 1.0_f64.tanh()).abs() < 1e-12);

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_ffcrow_arc_function_spellings() {
        /* ARCSIN/ARCCOS/ARCTAN are synonyms of ASIN/ACOS/ATAN */
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));

            let r = eval_dbl::<1>(&mut f, "ARCSIN(1)");
            assert!((r[0] - core::f64::consts::FRAC_PI_2).abs() < 1e-12);
            let r = eval_dbl::<1>(&mut f, "ARCCOS(0)");
            assert!((r[0] - core::f64::consts::FRAC_PI_2).abs() < 1e-12);
            let r = eval_dbl::<1>(&mut f, "ARCTAN(1)");
            assert!((r[0] - core::f64::consts::FRAC_PI_4).abs() < 1e-12);

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_ffcrow_angsep_function() {
        /* "angular separation"  angsep(ra1,dec1,ra2,dec2)  (all in degrees) */
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));

            /* 0*INTCOL keeps the arguments non-constant, so these go
            through Do_Func rather than being folded at parse time - see
            test_ffcrow_angsep_constant_arguments below */
            let r = eval_dbl::<3>(&mut f, "ANGSEP(0.0*INTCOL, 0.0, 0.0, 90.0)");
            for v in r {
                assert!((v - 90.0).abs() < 1e-9, "got {v}");
            }
            let r = eval_dbl::<1>(&mut f, "ANGSEP(0.0*INTCOL, 0.0, 90.0, 0.0)");
            assert!((r[0] - 90.0).abs() < 1e-9);
            let r = eval_dbl::<1>(&mut f, "ANGSEP(10.0 + 0.0*INTCOL, 20.0, 10.0, 20.0)");
            assert!(r[0].abs() < 1e-9);
            /* 1 degree of RA at declination 60 is 0.5 degrees on the sky */
            let r = eval_dbl::<1>(&mut f, "ANGSEP(0.0*INTCOL, 60.0, 1.0, 60.0)");
            assert!((r[0] - 0.5).abs() < 1e-3);
            /* the declination separation grows with the row number */
            let r = eval_dbl::<3>(&mut f, "ANGSEP(0.0, 0.0, 0.0, 9.0*INTCOL)");
            assert!((r[0] - 9.0).abs() < 1e-9);
            assert!((r[1] - 18.0).abs() < 1e-9);
            assert!((r[2] - 27.0).abs() < 1e-9);

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_ffcrow_angsep_constant_arguments() {
        /* When all four arguments are constants ANGSEP is folded at parse
        time, and the folded result must agree with the per-row one computed
        by test_ffcrow_angsep_function above.

        DEVIATION from CFITSIO 4.7.0: there "case angsep_fct:" is missing its
        "break;" and falls through into "case min1_fct:", so the folded form
        returns its first argument instead of the separation.  Fix submitted
        upstream; see the comment in New_Func in eval_y.rs. */
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));

            let r = eval_dbl::<1>(&mut f, "ANGSEP(0,0,0,90)");
            assert!((r[0] - 90.0).abs() < 1e-9, "got {}", r[0]);
            let r = eval_dbl::<1>(&mut f, "ANGSEP(7.5,0,0,90)");
            assert!((r[0] - 90.0).abs() < 1e-9, "got {}", r[0]);
            let r = eval_dbl::<1>(&mut f, "ANGSEP(0,0,90,0)");
            assert!((r[0] - 90.0).abs() < 1e-9, "got {}", r[0]);
            let r = eval_dbl::<1>(&mut f, "ANGSEP(10,20,10,20)");
            assert!(r[0].abs() < 1e-9, "got {}", r[0]);
            /* 1 degree of RA at declination 60 is 0.5 degrees on the sky */
            let r = eval_dbl::<1>(&mut f, "ANGSEP(0,60,1,60)");
            assert!((r[0] - 0.5).abs() < 1e-3, "got {}", r[0]);

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_ffcrow_random_functions() {
        /* random(), randomn() and randomp(x) */
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));

            let r = eval_dbl::<10>(&mut f, "RANDOM()");
            for v in r {
                assert!((0.0..1.0).contains(&v), "random() out of range: {v}");
            }
            assert!(r.iter().any(|v| *v != r[0]), "random() returned constants");

            /* a normal deviate - just check it is finite and not always zero */
            let r = eval_dbl::<10>(&mut f, "RANDOMN()");
            for v in r {
                assert!(v.is_finite());
            }
            assert!(r.iter().any(|v| *v != 0.0));

            /* Poisson deviates are non-negative integers */
            let r = eval_lng::<10>(&mut f, "RANDOMP(5)");
            for v in r {
                assert!(v >= 0);
            }

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_fffrow_string_functions() {
        /* "substring" strmid(s,p,n), "string search" strstr(s,r) */
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));

            /* STRCOL row 1 is 'alpha' */
            let (n, flags) = count_rows(&mut f, "STRMID(STRCOL,1,3) == 'alp'");
            assert_eq!(n, 1);
            assert_eq!(flags[0], 1);
            assert_eq!(count_rows(&mut f, "STRMID(STRCOL,2,3) == 'lph'").0, 1);

            /* position of the first 'a' in each string; 'epsilon' has none,
            so STRSTR returns NULL there */
            let r = eval_lng::<10>(&mut f, "DEFNULL(STRSTR(STRCOL,'a'), -1)");
            assert_eq!(r, [1, 4, 2, 5, -1, 4, 3, 5, 4, 2]);
            assert_eq!(count_rows(&mut f, "ISNULL(STRSTR(STRCOL,'a'))").0, 1);
            /* 'ta' starts at position 3 of beta, zeta and iota */
            assert_eq!(count_rows(&mut f, "STRSTR(STRCOL,'ta') == 3").0, 3);

            /* strings can be concatenated with '+' */
            let (n, flags) = count_rows(&mut f, "(STRCOL + '!') == 'alpha!'");
            assert_eq!(n, 1);
            assert_eq!(flags[0], 1);

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_fffrow_row_offset() {
        /* "PHA{-3} will evaluate to the value of column PHA, 3 rows above the
        row currently being processed.  Rows that fall outside the table will
        be treated as undefined, or NULLs." */
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));

            let r = eval_lng::<5>(&mut f, "DEFNULL(INTCOL{-1}, -99)");
            assert_eq!(r, [-99, 1, 2, 3, 4]);

            let r = eval_lng::<5>(&mut f, "DEFNULL(INTCOL{2}, -99)");
            assert_eq!(r, [3, 4, 5, 6, 7]);

            let (n, flags) = count_rows(&mut f, "ISNULL(INTCOL{-2})");
            assert_eq!(n, 2);
            assert_eq!(flags[0], 1);
            assert_eq!(flags[1], 1);
            assert_eq!(flags[2], 0);

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_ffcrow_accum_and_seqdiff() {
        /* "cumulative sum" accum(x), "sequential difference" seqdiff(x);
        the two are functional inverses */
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));

            /* over rows, accum runs through the whole column */
            let r = eval_lng::<5>(&mut f, "ACCUM(INTCOL)");
            assert_eq!(r, [1, 3, 6, 10, 15]);

            let r = eval_lng::<5>(&mut f, "SEQDIFF(ACCUM(INTCOL))");
            assert_eq!(r, [1, 2, 3, 4, 5]);

            let r = eval_dbl::<3>(&mut f, "ACCUM(FLOATCOL)");
            assert!((r[0] - 1.5).abs() < 1e-6);
            assert!((r[1] - 4.0).abs() < 1e-6);
            assert!((r[2] - 7.5).abs() < 1e-6);

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_fffrow_approx_equal() {
        /* "approx. equal(1e-7)"  ~  -- true when |x-y| < 1.0e-7 */
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));

            let (n, flags) = count_rows(&mut f, "FLOATCOL ~ 1.5");
            assert_eq!(n, 1);
            assert_eq!(flags[0], 1);

            let (n, flags) = count_rows(&mut f, "INTCOL ~ 5");
            assert_eq!(n, 1);
            assert_eq!(flags[4], 1);

            /* a difference below the 1e-7 threshold still compares equal ... */
            assert_eq!(count_rows(&mut f, "FLOATCOL ~ (1.5 + 1.0e-9)").0, 1);
            /* ... but one above it does not */
            assert_eq!(count_rows(&mut f, "FLOATCOL ~ (1.5 + 1.0e-5)").0, 0);

            fits_close_file(f, &mut status);
        });
    }

    /// Single-row table exercising the vector, bit and dereference syntax.
    ///
    /// VECCOL is the 5-vector [1,2,3,4,5], IVEC the 5-vector [10,20,30,40,50],
    /// MATRIX a 2x3 array holding 1..6, BITS an 8-bit column set to 11110000,
    /// and "MY COL" a scalar whose name needs the $...$ quoting syntax.
    fn create_vector_table(name: &[c_char]) -> Box<fitsfile> {
        let mut status = 0;
        let vecdata: [f64; 5] = [1.0, 2.0, 3.0, 4.0, 5.0];
        let ivecdata: [c_long; 5] = [10, 20, 30, 40, 50];
        let matdata: [c_long; 6] = [1, 2, 3, 4, 5, 6];
        let tdim: [c_long; 2] = [2, 3];
        let bits: [c_char; 8] = [1, 1, 1, 1, 0, 0, 0, 0];
        let spaced: [c_long; 1] = [42];

        let mut f: Option<Box<fitsfile>> = None;
        fits_create_file(&mut f, name, &mut status);
        fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
        create_tbl(
            f.as_deref_mut().unwrap(),
            BINARY_TBL,
            1,
            &["VECCOL", "IVEC", "MATRIX", "BITS", "MY COL"],
            &["5D", "5J", "6J", "8X", "1J"],
            None,
            &mut status,
        );
        fits_write_tdim(f.as_deref_mut().unwrap(), 3, 2, &tdim, &mut status);
        fits_write_col_dbl(f.as_deref_mut().unwrap(), 1, 1, 1, 5, &vecdata, &mut status);
        fits_write_col_lng(
            f.as_deref_mut().unwrap(),
            2,
            1,
            1,
            5,
            &ivecdata,
            &mut status,
        );
        fits_write_col_lng(f.as_deref_mut().unwrap(), 3, 1, 1, 6, &matdata, &mut status);
        fits_write_col_bit(f.as_deref_mut().unwrap(), 4, 1, 1, 8, &bits, &mut status);
        fits_write_col_lng(f.as_deref_mut().unwrap(), 5, 1, 1, 1, &spaced, &mut status);
        assert_eq!(status, 0, "create_vector_table setup failed");
        f.unwrap()
    }

    #[test]
    fn test_ffcrow_vector_functions() {
        /* MIN/MAX/AVERAGE/MEDIAN/SUM/STDDEV/NELEM/NVALID/NAXIS/NAXES */
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_vector_table(&to_buf(filename));

            let r = eval_dbl::<1>(&mut f, "MIN(VECCOL)");
            assert!((r[0] - 1.0).abs() < 1e-9);
            let r = eval_dbl::<1>(&mut f, "MAX(VECCOL)");
            assert!((r[0] - 5.0).abs() < 1e-9);
            let r = eval_dbl::<1>(&mut f, "SUM(VECCOL)");
            assert!((r[0] - 15.0).abs() < 1e-9);
            let r = eval_dbl::<1>(&mut f, "AVERAGE(VECCOL)");
            assert!((r[0] - 3.0).abs() < 1e-9);
            let r = eval_dbl::<1>(&mut f, "MEDIAN(VECCOL)");
            assert!((r[0] - 3.0).abs() < 1e-9);
            /* the sample standard deviation, i.e. 1/sqrt(N-1) */
            let r = eval_dbl::<1>(&mut f, "STDDEV(VECCOL)");
            assert!((r[0] - 2.5_f64.sqrt()).abs() < 1e-9, "got {}", r[0]);

            assert_eq!(eval_lng::<1>(&mut f, "NELEM(VECCOL)"), [5]);
            assert_eq!(eval_lng::<1>(&mut f, "NVALID(VECCOL)"), [5]);
            assert_eq!(eval_lng::<1>(&mut f, "NAXIS(VECCOL)"), [1]);
            assert_eq!(eval_lng::<1>(&mut f, "NAXES(VECCOL,1)"), [5]);

            /* the same functions over an integer vector */
            assert_eq!(eval_lng::<1>(&mut f, "MIN(IVEC)"), [10]);
            assert_eq!(eval_lng::<1>(&mut f, "MAX(IVEC)"), [50]);
            assert_eq!(eval_lng::<1>(&mut f, "SUM(IVEC)"), [150]);

            /* a 2x3 array */
            assert_eq!(eval_lng::<1>(&mut f, "NAXIS(MATRIX)"), [2]);
            assert_eq!(eval_lng::<1>(&mut f, "NAXES(MATRIX,1)"), [2]);
            assert_eq!(eval_lng::<1>(&mut f, "NAXES(MATRIX,2)"), [3]);
            assert_eq!(eval_lng::<1>(&mut f, "NELEM(MATRIX)"), [6]);

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_ffcrow_vector_boolean_sum() {
        /* "If V is a boolean vector, SUM returns the number of TRUE elements";
        the documented all-elements idiom is SUM(A>B) == NELEM(A) */
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_vector_table(&to_buf(filename));

            assert_eq!(eval_lng::<1>(&mut f, "SUM(VECCOL > 2)"), [3]);
            assert_eq!(
                eval_log::<1>(&mut f, "SUM(VECCOL > 0) == NELEM(VECCOL)"),
                [1]
            );
            assert_eq!(
                eval_log::<1>(&mut f, "SUM(VECCOL > 2) == NELEM(VECCOL)"),
                [0]
            );
            /* elementwise comparison of two vectors */
            assert_eq!(eval_lng::<1>(&mut f, "SUM(IVEC > VECCOL)"), [5]);

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_ffcrow_vector_element_position_functions() {
        /* ELEMENTNUM(V) and AXISELEM(V,n) return vectors of the same size */
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_vector_table(&to_buf(filename));

            assert_eq!(eval_lng::<5>(&mut f, "ELEMENTNUM(VECCOL)"), [1, 2, 3, 4, 5]);
            assert_eq!(eval_lng::<5>(&mut f, "AXISELEM(VECCOL,1)"), [1, 2, 3, 4, 5]);
            /* for the 2x3 array the axis positions cycle per axis */
            assert_eq!(
                eval_lng::<6>(&mut f, "ELEMENTNUM(MATRIX)"),
                [1, 2, 3, 4, 5, 6]
            );
            assert_eq!(
                eval_lng::<6>(&mut f, "AXISELEM(MATRIX,1)"),
                [1, 2, 1, 2, 1, 2]
            );
            assert_eq!(
                eval_lng::<6>(&mut f, "AXISELEM(MATRIX,2)"),
                [1, 1, 2, 2, 3, 3]
            );

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_ffcrow_array_function() {
        /* "promote to array"  ARRAY(X,d) */
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_vector_table(&to_buf(filename));

            /* a scalar promoted to a 4-vector */
            assert_eq!(eval_lng::<4>(&mut f, "ARRAY(3,4)"), [3, 3, 3, 3]);
            assert_eq!(eval_lng::<1>(&mut f, "NELEM(ARRAY(3,4))"), [4]);
            /* a multi-dimensional shape given as a vector of dimensions */
            assert_eq!(
                eval_lng::<6>(&mut f, "ARRAY(0,{2,3,1})"),
                [0, 0, 0, 0, 0, 0]
            );
            assert_eq!(eval_lng::<1>(&mut f, "NAXIS(ARRAY(0,{2,3,1}))"), [3]);
            assert_eq!(eval_lng::<1>(&mut f, "NAXES(ARRAY(0,{2,3,1}),2)"), [3]);
            /* a column promoted to an array */
            assert_eq!(eval_lng::<3>(&mut f, "ARRAY($MY COL$,3)"), [42, 42, 42]);
            /* re-dimensioning an existing vector */
            assert_eq!(eval_lng::<1>(&mut f, "NAXIS(ARRAY(MATRIX,6))"), [1]);
            assert_eq!(eval_lng::<6>(&mut f, "ARRAY(MATRIX,6)"), [1, 2, 3, 4, 5, 6]);

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_ffcrow_vector_literal() {
        /* "Vectors can be manually constructed ... using a comma-separated
        list of elements surrounded by curly braces" */
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_vector_table(&to_buf(filename));

            assert_eq!(eval_lng::<4>(&mut f, "{1,3,6,1}"), [1, 3, 6, 1]);
            assert_eq!(eval_lng::<1>(&mut f, "SUM({1,3,6,1})"), [11]);
            assert_eq!(eval_lng::<1>(&mut f, "NELEM({1,3,6,1})"), [4]);
            assert_eq!(eval_lng::<1>(&mut f, "MAX({1,3,6,1})"), [6]);
            /* "Any elements which are themselves vectors will be expanded
            out with each of its elements becoming an element" */
            assert_eq!(eval_lng::<1>(&mut f, "NELEM({IVEC,0})"), [6]);
            assert_eq!(eval_lng::<6>(&mut f, "{IVEC,0}"), [10, 20, 30, 40, 50, 0]);
            /* "elements will be promoted to the highest data type present" */
            let r = eval_dbl::<3>(&mut f, "{1,2.5,3}");
            assert!((r[1] - 2.5).abs() < 1e-9);

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_ffcrow_vector_dereference() {
        /* single elements are selected with a bracketed index list; a C-like
        reversed-index form is also accepted */
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_vector_table(&to_buf(filename));

            assert_eq!(eval_lng::<1>(&mut f, "IVEC[2]"), [20]);
            /* the index may itself be an expression */
            assert_eq!(eval_lng::<1>(&mut f, "IVEC[NELEM(IVEC)]"), [50]);
            assert_eq!(eval_lng::<1>(&mut f, "IVEC[1+1]"), [20]);

            /* MATRIX is 2x3 holding 1..6 in column-major (Fortran) order */
            assert_eq!(eval_lng::<1>(&mut f, "MATRIX[1,2]"), [3]);
            assert_eq!(eval_lng::<1>(&mut f, "MATRIX[2,3]"), [6]);
            /* the C-like syntax reverses the index order */
            assert_eq!(eval_lng::<1>(&mut f, "MATRIX[2][1]"), [3]);
            assert_eq!(eval_lng::<1>(&mut f, "MATRIX[3][2]"), [6]);
            /* with fewer indices a slice of the array is extracted */
            assert_eq!(eval_lng::<2>(&mut f, "MATRIX[3]"), [5, 6]);

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_ffcrow_dollar_quoted_column_name() {
        /* "Vector columns which contain spaces or arithmetic operators must
        have their names enclosed in '$' characters" */
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_vector_table(&to_buf(filename));

            assert_eq!(eval_lng::<1>(&mut f, "$MY COL$"), [42]);
            assert_eq!(eval_lng::<1>(&mut f, "$MY COL$ * 2"), [84]);
            assert_eq!(eval_log::<1>(&mut f, "$MY COL$ > 40"), [1]);

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_ffcrow_vector_functions_with_nulls() {
        /* "The first 6 of these functions ignore any null values"; NELEM
        counts all elements while NVALID counts the non-null ones */
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let mut f = create_vector_table(&name);

            /* make the 3rd element of VECCOL undefined */
            fits_write_col_null(&mut f, 1, 1, 3, 1, &mut status);
            assert_eq!(status, 0);

            assert_eq!(eval_lng::<1>(&mut f, "NELEM(VECCOL)"), [5]);
            assert_eq!(eval_lng::<1>(&mut f, "NVALID(VECCOL)"), [4]);
            let r = eval_dbl::<1>(&mut f, "SUM(VECCOL)");
            assert!((r[0] - 12.0).abs() < 1e-9, "got {}", r[0]);
            let r = eval_dbl::<1>(&mut f, "AVERAGE(VECCOL)");
            assert!((r[0] - 3.0).abs() < 1e-9, "got {}", r[0]);
            let r = eval_dbl::<1>(&mut f, "MIN(VECCOL)");
            assert!((r[0] - 1.0).abs() < 1e-9);
            let r = eval_dbl::<1>(&mut f, "MAX(VECCOL)");
            assert!((r[0] - 5.0).abs() < 1e-9);

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_ffcrow_bit_masks() {
        /* Bit masks in binary, octal and hex form, with 'x' wildcards */
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_vector_table(&to_buf(filename));

            assert_eq!(eval_log::<1>(&mut f, "BITS == b11110000"), [1]);
            assert_eq!(eval_log::<1>(&mut f, "BITS .eq. b11110000"), [1]);
            assert_eq!(eval_log::<1>(&mut f, "BITS == b00001111"), [0]);
            assert_eq!(eval_log::<1>(&mut f, "BITS != b00001111"), [1]);
            /* 'x' is a wildcard bit */
            assert_eq!(eval_log::<1>(&mut f, "BITS == bxxxx0000"), [1]);
            assert_eq!(eval_log::<1>(&mut f, "BITS == b1111xxxx"), [1]);
            assert_eq!(eval_log::<1>(&mut f, "BITS == bxxxx1xxx"), [0]);
            /* octal and hex masks describe the same pattern */
            assert_eq!(eval_log::<1>(&mut f, "BITS == o360"), [1]);
            assert_eq!(eval_log::<1>(&mut f, "BITS == hF0"), [1]);
            assert_eq!(eval_log::<1>(&mut f, "BITS == hf0"), [1]);
            /* NELEM returns the column width for bit columns */
            assert_eq!(eval_lng::<1>(&mut f, "NELEM(BITS)"), [8]);

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_ffcrow_bit_mask_ordering() {
        /* "It is also possible to test if a range of bits is less than, less
        than equal, greater than and greater than equal to a ... value" */
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_vector_table(&to_buf(filename));

            assert_eq!(eval_log::<1>(&mut f, "BITS <= b11111111"), [1]);
            assert_eq!(eval_log::<1>(&mut f, "BITS >= b11110000"), [1]);
            assert_eq!(eval_log::<1>(&mut f, "BITS > b00000000"), [1]);
            assert_eq!(eval_log::<1>(&mut f, "BITS < b11111111"), [1]);
            assert_eq!(eval_log::<1>(&mut f, "BITS > b11111111"), [0]);
            /* a wildcard limits which bits take part in the comparison */
            assert_eq!(eval_log::<1>(&mut f, "BITS .gt. bxxx100xx"), [0]);
            assert_eq!(eval_log::<1>(&mut f, "BITS .le. b1xxxxxxx"), [1]);

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_ffcrow_bit_operators() {
        /* "Bit wise AND, OR and NOT operations are also possible ... All of
        these operators result in a bit field"; bit fields can also be
        appended with '+' */
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_vector_table(&to_buf(filename));

            assert_eq!(eval_log::<1>(&mut f, "(!BITS) == b00001111"), [1]);
            assert_eq!(
                eval_log::<1>(&mut f, "(BITS & b10000001) == b10000000"),
                [1]
            );
            assert_eq!(
                eval_log::<1>(&mut f, "(BITS | b00000001) == b11110001"),
                [1]
            );
            assert_eq!(
                eval_log::<1>(&mut f, "(BITS + BITS) == b1111000011110000"),
                [1]
            );
            assert_eq!(eval_lng::<1>(&mut f, "NELEM(BITS + BITS)"), [16]);

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_fffrow_region_shape_functions() {
        /* CIRCLE(xcen,ycen,rad,x,y), BOX(xcen,ycen,xwid,ywid,rot,x,y) and
        ELLIPSE(xcen,ycen,xrad,yrad,rot,x,y) */
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));

            /* points (1,1.5) ... (10,10.5) against a circle of radius 5
            centred on the origin: rows 1..3 are inside */
            let (n, flags) = count_rows(&mut f, "CIRCLE(0, 0, 5, INTCOL, FLOATCOL)");
            assert_eq!(n, 3, "flags {flags:?}");
            assert_eq!(flags[2], 1);
            assert_eq!(flags[3], 0);

            /* an unrotated 8x8 box centred on the origin holds rows 1..3 */
            let (n, _) = count_rows(&mut f, "BOX(0, 0, 8, 8, 0, INTCOL, FLOATCOL)");
            assert_eq!(n, 3);
            /* rotating the box by 90 degrees leaves it unchanged */
            assert_eq!(
                count_rows(&mut f, "BOX(0, 0, 8, 8, 90, INTCOL, FLOATCOL)").0,
                3
            );
            /* a wide, flat box only admits the first row */
            assert_eq!(
                count_rows(&mut f, "BOX(0, 0, 20, 4, 0, INTCOL, FLOATCOL)").0,
                1
            );

            /* an ellipse with semi-axes 5 and 5 behaves like the circle */
            assert_eq!(
                count_rows(&mut f, "ELLIPSE(0, 0, 5, 5, 0, INTCOL, FLOATCOL)").0,
                3
            );
            assert_eq!(
                count_rows(&mut f, "ELLIPSE(0, 0, 20, 2, 0, INTCOL, FLOATCOL)").0,
                1
            );

            fits_close_file(f, &mut status);
        });
    }

    /// Build a file holding a 10-row event table (TIME = 0.5, 1.5, ... 9.5)
    /// followed by a GTI extension covering 0-2, 4-6 and 8-10, and leave the
    /// file positioned on the event table.
    fn create_gti_table(name: &[c_char]) -> Box<fitsfile> {
        let mut status = 0;
        let times: [f64; 10] = [0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5];
        let start: [f64; 3] = [0.0, 4.0, 8.0];
        let stop: [f64; 3] = [2.0, 6.0, 10.0];

        let mut f: Option<Box<fitsfile>> = None;
        fits_create_file(&mut f, name, &mut status);
        fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);

        create_tbl(
            f.as_deref_mut().unwrap(),
            BINARY_TBL,
            10,
            &["TIME"],
            &["1D"],
            None,
            &mut status,
        );
        fits_write_col_dbl(f.as_deref_mut().unwrap(), 1, 1, 1, 10, &times, &mut status);

        create_tbl(
            f.as_deref_mut().unwrap(),
            BINARY_TBL,
            3,
            &["START", "STOP"],
            &["1D", "1D"],
            None,
            &mut status,
        );
        let extname = cc("STDGTI");
        fits_update_key_str(
            f.as_deref_mut().unwrap(),
            &cc("EXTNAME"),
            &extname,
            None,
            &mut status,
        );
        fits_write_col_dbl(f.as_deref_mut().unwrap(), 1, 1, 1, 3, &start, &mut status);
        fits_write_col_dbl(f.as_deref_mut().unwrap(), 2, 1, 1, 3, &stop, &mut status);

        /* leave the file on the event table */
        fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
        assert_eq!(status, 0, "create_gti_table setup failed");
        f.unwrap()
    }

    #[test]
    fn test_fffrow_gtifilter() {
        /* gtifilter( [ "gtifile" [, expr [, "STARTCOL", "STOPCOL" ] ] ] ) -
        "gtifilter()" is equivalent to
        gtifilter( "", TIME, "*START*", "*STOP*" ) */
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_gti_table(&to_buf(filename));

            let (n, flags) = count_rows(&mut f, "gtifilter()");
            assert_eq!(n, 6, "flags {flags:?}");
            assert_eq!(flags, [1, 1, 0, 0, 1, 1, 0, 0, 1, 1]);

            /* the same, spelled out */
            let (n, flags) = count_rows(&mut f, "gtifilter('', TIME, '*START*', '*STOP*')");
            assert_eq!(n, 6);
            assert_eq!(flags, [1, 1, 0, 0, 1, 1, 0, 0, 1, 1]);

            /* an arbitrary time expression may be filtered: TIME+2 is in a
            GTI for rows 3,4 (4.5,5.5) and 7,8 (8.5,9.5) only */
            let (n, flags) = count_rows(&mut f, "gtifilter('', TIME + 2)");
            assert_eq!(n, 4, "flags {flags:?}");
            assert_eq!(flags, [0, 0, 1, 1, 0, 0, 1, 1, 0, 0]);

            /* the GTI extension may be named explicitly */
            assert_eq!(count_rows(&mut f, "gtifilter('[STDGTI]')").0, 6);
            assert_eq!(
                count_rows(&mut f, "gtifilter('', TIME, 'START', 'STOP')").0,
                6
            );

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_ffcrow_gtifind() {
        /* "gtifind() returns the row number in the GTI table that matches the
        time sample, or -1 if the time sample is not within any GTI" */
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_gti_table(&to_buf(filename));

            /* rows outside every GTI are returned as undefined, which
            without a null value becomes 0 */
            let r = eval_lng::<10>(&mut f, "gtifind('', TIME)");
            assert_eq!(r, [1, 1, 0, 0, 2, 2, 0, 0, 3, 3]);

            /* supplying a null value shows which rows are undefined; note
            that gtifind() is a bexpr in the grammar, so it cannot be nested
            inside DEFNULL()/ISNULL() even though it yields an integer */
            let mut results = [0 as c_long; 10];
            let nulval: c_long = -1;
            let mut anynul = 0;
            fits_calc_rows(
                &mut f,
                TLONG,
                &cc("gtifind('', TIME)"),
                1,
                10,
                (&raw const nulval).cast::<c_void>(),
                as_bytes_mut(&mut results),
                &mut anynul,
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(anynul, 1);
            assert_eq!(results, [1, 1, -1, -1, 2, 2, -1, -1, 3, 3]);

            let r = eval_lng::<10>(&mut f, "gtifind('', TIME, 'START', 'STOP')");
            assert_eq!(r, [1, 1, 0, 0, 2, 2, 0, 0, 3, 3]);

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_ffcrow_gtioverlap() {
        /* gtioverlap( "gtifile", startExpr, stopExpr [, "STARTCOL", "STOPCOL"] )
        computes the overlap between the requested range and the GTIs */
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_gti_table(&to_buf(filename));

            /* one second bins starting at each TIME value */
            let r = eval_dbl::<10>(&mut f, "gtioverlap('', TIME, TIME + 1)");
            let expect = [1.0, 0.5, 0.0, 0.5, 1.0, 0.5, 0.0, 0.5, 1.0, 0.5];
            for (i, e) in expect.iter().enumerate() {
                assert!(
                    (r[i] - e).abs() < 1e-9,
                    "row {}: got {} want {e}",
                    i + 1,
                    r[i]
                );
            }

            /* a range spanning every interval sees the full 6 seconds of GTI */
            let r = eval_dbl::<1>(&mut f, "gtioverlap('', 0*TIME, 0*TIME + 20)");
            assert!((r[0] - 6.0).abs() < 1e-9, "got {}", r[0]);

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_fffrow_regfilter() {
        /* regfilter( "regfilename" [, Xexpr, Yexpr [, "wcs cols" ] ] ) */
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));

            let dir = std::path::Path::new(filename).parent().unwrap();
            let regpath = dir.join("test_region.reg");
            std::fs::write(&regpath, "circle(5,5,3)\n").unwrap();
            let regpath = regpath.to_str().unwrap();

            /* the points are (INTCOL, FLOATCOL) = (1,1.5) ... (10,10.5);
            rows 3 to 6 lie within 3 units of (5,5) */
            let (n, flags) =
                count_rows(&mut f, &format!("regfilter('{regpath}', INTCOL, FLOATCOL)"));
            assert_eq!(n, 4, "flags {flags:?}");
            assert_eq!(flags, [0, 0, 1, 1, 1, 1, 0, 0, 0, 0]);

            /* an excluded region is subtracted from the accepted area */
            std::fs::write(regpath, "circle(5,5,3)\n-circle(5,5,1)\n").unwrap();
            let (n, flags) =
                count_rows(&mut f, &format!("regfilter('{regpath}', INTCOL, FLOATCOL)"));
            assert_eq!(n, 3, "flags {flags:?}");
            assert_eq!(flags[4], 0); /* (5,5.5) is inside the excluded circle */

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_fffrow_boolean_literals() {
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));

            /* note that for a constant expression fffrow reports every row
            as "good" regardless of the value, so check the flags */
            assert_eq!(count_rows(&mut f, "T").1, [1; 10]);
            assert_eq!(count_rows(&mut f, "F").1, [0; 10]);
            assert_eq!(count_rows(&mut f, "t").1, [1; 10]);
            assert_eq!(count_rows(&mut f, "f").1, [0; 10]);
            /* BOOLCOL alternates true/false starting with true */
            assert_eq!(count_rows(&mut f, "BOOLCOL == T").0, 5);
            assert_eq!(count_rows(&mut f, "BOOLCOL != T").0, 5);
            assert_eq!(count_rows(&mut f, "BOOLCOL .and. T").0, 5);
            assert_eq!(count_rows(&mut f, "BOOLCOL .or. T").0, 10);

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_fffrow_string_comparisons() {
        /* the comparison operators also apply to string expressions */
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));

            let (n, flags) = count_rows(&mut f, "STRCOL == 'gamma'");
            assert_eq!(n, 1);
            assert_eq!(flags[2], 1);
            assert_eq!(count_rows(&mut f, "STRCOL != 'gamma'").0, 9);
            /* only 'alpha' sorts before 'beta' */
            assert_eq!(count_rows(&mut f, "STRCOL < 'beta'").0, 1);
            assert_eq!(count_rows(&mut f, "STRCOL <= 'beta'").0, 2);
            /* 'theta' and 'zeta' sort after 'kappa' */
            assert_eq!(count_rows(&mut f, "STRCOL > 'kappa'").0, 2);
            assert_eq!(count_rows(&mut f, "STRCOL >= 'kappa'").0, 3);
            /* double quotes may be used instead of single quotes */
            assert_eq!(count_rows(&mut f, "STRCOL == \"gamma\"").0, 1);

            /* NELEM of a string column is its declared width */
            assert_eq!(eval_lng::<1>(&mut f, "NELEM(STRCOL)"), [10]);

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_fffrow_string_conditional() {
        /* "x and y can be any scalar data type (including string)" */
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));

            let (n, flags) = count_rows(&mut f, "(INTCOL > 5 ? 'big' : 'small') == 'big'");
            assert_eq!(n, 5);
            assert_eq!(flags[4], 0);
            assert_eq!(flags[5], 1);
            /* string concatenation inside the branches */
            assert_eq!(
                count_rows(&mut f, "(BOOLCOL ? STRCOL + '!' : STRCOL) == 'alpha!'").0,
                1
            );

            fits_close_file(f, &mut status);
        });
    }

    /// View a typed const slice as a byte slice for fits_write_col(TLONG, ...).
    fn as_bytes_const<T>(s: &[T]) -> &[u8] {
        unsafe { core::slice::from_raw_parts(s.as_ptr().cast::<u8>(), core::mem::size_of_val(s)) }
    }
}
