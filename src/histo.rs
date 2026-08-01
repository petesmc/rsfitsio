/*   Globally defined histogram parameters */

use core::ffi::{CStr, c_void};
use core::slice;
use core::{array, ptr};

use crate::bb;
use crate::c_types::{c_char, c_int, c_long, c_short};
use crate::cfileio::ffclos_safe;
use crate::cfileio::ffimport_file_safe;
use crate::cfileio::ffinit_safe;
use crate::cfileio::fits_find_match_delim;
use crate::cfileio::fits_get_token2_safe;
use crate::eval_defs::MAXDIMS;
use crate::eval_defs::Node;
use crate::eval_defs::ParseStatusVariables;
use crate::eval_f::ffcprs;
use crate::eval_f::ffiprs;
use crate::eval_f::fits_parser_set_temporary_col;
use crate::eval_f::fits_parser_workfn_safe;
use crate::fitscore::ffgcno_safe;
use crate::fitscore::ffgnrw_safe;
use crate::fitscore::ffmahd_safe;
use crate::fitscore::ffpmsg_slice;
use crate::fitscore::ffpmsg_str;
use crate::fitscore::fits_copy_pixlist2image_safe;
use crate::fitscore::fits_recalloc;
use crate::putcol::ffiter_safe;
use crate::putcol::fits_iter_set_by_num_safe;

use crate::aliases::rust_api::*;
use crate::putcol::fits_iter_get_array_safe;
use crate::putcol::fits_iter_set_datatype_safe;
use crate::putcol::fits_iter_set_file_safe;
use crate::putcol::fits_iter_set_iotype_safe;

// Alias for fits_iterate_data
use crate::putcol::ffiter_safe as fits_iterate_data;
use crate::putkey::ffcrim_safe;
use crate::wrappers::strcat_safe;
use crate::wrappers::strcpy_safe;
use crate::wrappers::strcspn_safe;
use crate::wrappers::strlen_safe;
use crate::wrappers::strncat_safe;
use crate::wrappers::strtod_safe;
use crate::wrappers::{isdigit_safe, strlen};

use bytemuck::{cast_slice, cast_slice_mut};

use crate::cs;
use crate::eval_defs::{ParseData, parseInfo};
use crate::fitscore::ffkeyn_safe;
use crate::fitsio::*;
use crate::getcol::ffgcv_safe;
use crate::getkey::ffgky_safe;
use crate::putkey::ffpky_safe;
use crate::{KeywordDatatype, KeywordDatatypeMut, NullValue, raw_to_slice};

/*  Structure holding all the histogramming information   */
#[derive(Default)]
pub(crate) struct HistType<'a> {
    hist: HistUnion,
    tblptr: *mut fitsfile,
    haxis: c_int,
    hcolnum: [c_int; 4],
    himagetype: c_int,
    haxis1: c_long,
    haxis2: c_long,
    haxis3: c_long,
    haxis4: c_long,
    amin1: f64,
    amin2: f64,
    amin3: f64,
    amin4: f64,
    maxbin1: f64,
    maxbin2: f64,
    maxbin3: f64,
    maxbin4: f64,
    binsize1: f64,
    binsize2: f64,
    binsize3: f64,
    binsize4: f64,
    incr: [c_long; 5],
    wtrecip: c_int,
    wtcolnum: c_int,
    wtexpr: *mut c_char,
    weight: f64,
    rowselector: *mut c_char,
    rowselector_cur: *mut c_char,
    repeat: c_long,
    startCols: [c_int; 5],
    numIterCols: c_int,
    iterCols: &'a mut [iteratorCol],
    parsers: &'a mut [ParseData],
    infos: &'a mut [parseInfo],
}

/*
The iterator work functions (ffwritehist, ffcalchist)
need to do their job... passed via *userPointer.
*/
pub(crate) union HistUnion {
    b: *mut c_char,
    i: *mut c_short,
    j: *mut c_int,
    r: *mut f32,
    d: *mut f64,
}

impl Default for HistUnion {
    fn default() -> Self {
        HistUnion { b: ptr::null_mut() }
    }
}

/*--------------------------------------------------------------------------*/
/// Parse the extended input binning specification string, returning
/// the binning parameters.  Supports up to 4 dimensions.  The binspec
/// string has one of these forms:
///
/// bin binsize                  - 2D histogram with binsize on each axis
/// bin xcol                     - 1D histogram on column xcol
/// bin (xcol, ycol) = binsize   - 2D histogram with binsize on each axis
/// bin x=min:max:size, y=min:max:size, z..., t...
/// bin x=:max, y=::size
/// bin x=size, y=min::size
/// bin x(expr), y(expr)=min:max:size, ...
///
/// most other reasonable combinations are supported. The (expr) is an
/// optional expression that will be calculated on the fly instead of
/// a table column name.  The name is still used for the output pixel
/// array metadata.
///
/// If expr == 0, then expressions are forbidden.  The caller does not
/// expect expressions.  
///
/// If exprs is non-zero, then upon return an array of expressions is
/// passed back to the caller.  Storage may be allocated by this routine,
/// If *exprs is non-zero upon return, the caller is responsible to
/// free(*exprs).  Upon return, the contains of exprs is,
///     (*exprs)[0] = expression for column 1 (or 0 if none)
///     (*exprs)[1] = expression for column 2 (or 0 if none)
///     (*exprs)[2] = expression for column 3 (or 0 if none)
///     (*exprs)[3] = expression for column 4 (or 0 if none)
///     (*exprs)[4] = expression for weighting (or 0 if none)
///
/// If the user specifies a column name and not an expression for bin
/// axis i, then the corresponding (*exprs)[i] will be a null pointer.
///
/// To be recognized as an expression, the weighting expression must be
/// enclosed in parentheses.
///
/// Expressions are never allowed using the bin (xcol,ycol) notation.
pub(crate) fn ffbinse(
    binspec: &[c_char],                      /* I - binning specification */
    imagetype: &mut c_int,                   /* O - image type, TINT or TSHORT */
    histaxis: &mut c_int,                    /* O - no. of axes in the histogram */
    colname: &mut [[c_char; FLEN_VALUE]; 4], /* column name for axis */
    minin: &mut [f64; 4],                    /* minimum value for each axis */
    maxin: &mut [f64; 4],                    /* maximum value for each axis */
    binsizein: &mut [f64; 4],                /* size of bins on each axis */
    minname: &mut [[c_char; FLEN_VALUE]; 4], /* keyword name for min */
    maxname: &mut [[c_char; FLEN_VALUE]; 4], /* keyword name for max */
    binname: &mut [[c_char; FLEN_VALUE]; 4], /* keyword name for binsize */
    weight: &mut f64,                        /* weighting factor          */
    wtname: &mut [c_char; 71],               /* keyword or column name for weight */
    recip: &mut c_int,                       /* the reciprocal of the weight? */
    exprs: Option<&mut [Box<[c_char]>; 5]>,  /* returned with expressions (or 0) */
    status: &mut c_int,
) -> c_int {
    let ii: c_int = 0;
    let mut slen: c_int;
    let mut defaulttype: c_int;
    let mut tmpname: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let file_expr: *mut c_char = core::ptr::null_mut();
    let mut dummy: f64 = 0.0;
    let mut exprbeg: [*mut c_char; 5] = [core::ptr::null_mut(); 5];
    let mut exprend: [*mut c_char; 5] = [core::ptr::null_mut(); 5];
    let mut has_exprs = false;

    if *status > 0 {
        return *status;
    }

    /* set the default values */
    *histaxis = 2;
    *imagetype = TINT;
    defaulttype = 1;
    *weight = 1.0;
    *recip = 0;
    wtname[0] = 0;

    /* set default values */
    for ii in 0..4 {
        colname[ii][0] = 0;
        minname[ii][0] = 0;
        maxname[ii][0] = 0;
        binname[ii][0] = 0;
        minin[ii] = DOUBLENULLVALUE; /* undefined values */
        maxin[ii] = DOUBLENULLVALUE;
        binsizein[ii] = DOUBLENULLVALUE;
    }

    let mut ptr = &binspec[3..]; /* skip over 'bin' */

    if ptr[0] == bb(b'i')
    /* bini */
    {
        *imagetype = TSHORT;
        defaulttype = 0;
        ptr = &ptr[1..];
    } else if ptr[0] == bb(b'j')
    /* binj; same as default */
    {
        defaulttype = 0;
        ptr = &ptr[1..];
    } else if ptr[0] == bb(b'r')
    /* binr */
    {
        *imagetype = TFLOAT;
        defaulttype = 0;
        ptr = &ptr[1..];
    } else if ptr[0] == bb(b'd')
    /* bind */
    {
        *imagetype = TDOUBLE;
        defaulttype = 0;
        ptr = &ptr[1..];
    } else if ptr[0] == bb(b'b')
    /* binb */
    {
        *imagetype = TBYTE;
        defaulttype = 0;
        ptr = &ptr[1..];
    }

    if ptr[0] == 0 {
        /* use all defaults for other parameters */
        return *status;
    } else if ptr[0] != bb(b' ')
    /* must be at least one blank */
    {
        ffpmsg_str("binning specification syntax error:");
        ffpmsg_slice(binspec);
        *status = URL_PARSE_ERROR;
        return *status;
    }

    while ptr[0] == bb(b' ') {
        /* skip over blanks */
        ptr = &ptr[1..];
    }

    if ptr[0] == 0 {
        /* no other parameters; use defaults */
        return *status;
    }

    /* Check if need to import expression from a file */
    let mut file_expr = None;
    let c: Vec<c_char>;
    if ptr[0] == bb(b'@') {
        if ffimport_file_safe(&ptr[1..], &mut file_expr, status) != 0 {
            return *status;
        }

        c = file_expr.unwrap();

        ptr = &c;
        while ptr[0] == bb(b' ') {
            ptr = &ptr[1..]; /* skip leading white space... again */
        }
    }

    // This look to used to emulate the goto getweight
    'getweight: loop {
        if ptr[0] == bb(b'(') {
            /* this must be the opening parenthesis around a list of column */
            /* names, optionally followed by a '=' and the binning spec. */

            for ii in 0..4 {
                ptr = &ptr[1..]; /* skip over the '(', ',', or ' ') */
                while ptr[0] == bb(b' ') {
                    /* skip over blanks */
                    ptr = &ptr[1..];
                }

                let slen = strcspn_safe(ptr, cs!(c" ,)"));
                if slen >= FLEN_VALUE {
                    ffpmsg_str("column name too long in binning specification");
                    ffpmsg_slice(binspec);
                    *status = URL_PARSE_ERROR;
                    return *status;
                }

                strncat_safe(&mut colname[ii], ptr, slen); /* copy 1st column name */

                ptr = &ptr[slen..];
                while ptr[0] == bb(b' ') {
                    /* skip over blanks */
                    ptr = &ptr[1..];
                }

                if ptr[0] == bb(b')')
                /* end of the list of names */
                {
                    *histaxis = (ii + 1) as c_int;
                    break;
                }
            }

            if ii == 4
            /* too many names in the list , or missing ')'  */
            {
                ffpmsg_str(
                    "binning specification has too many column names or is missing closing ')':",
                );
                ffpmsg_slice(binspec);
                *status = URL_PARSE_ERROR;
                return *status;
            }

            ptr = &ptr[1..]; /* skip over the closing parenthesis */
            while ptr[0] == bb(b' ') {
                /* skip over blanks */
                ptr = &ptr[1..];
            }

            if ptr[0] == 0 {
                return *status; /* parsed the entire string */
            } else if ptr[0] != bb(b'=')
            /* must be an equals sign now*/
            {
                ffpmsg_str("illegal binning specification in URL:");
                ffpmsg_str(" an equals sign '=' must follow the column names");
                ffpmsg_slice(binspec);
                *status = URL_PARSE_ERROR;
                return *status;
            }

            ptr = &ptr[1..]; /* skip over the equals sign */
            while ptr[0] == bb(b' ') {
                /* skip over blanks */
                ptr = &ptr[1..];
            }

            /* get the single range specification for all the columns */
            /* Note that the extended syntax is not allowed here */
            ffbinr_safe(
                &mut ptr,
                &mut tmpname,
                &mut minin[0],
                &mut maxin[0],
                &mut binsizein[0],
                &mut minname[0],
                &mut maxname[0],
                &mut binname[0],
                status,
            );
            if *status > 0 {
                ffpmsg_str("illegal binning specification in URL:");
                ffpmsg_slice(binspec);
                return *status;
            }

            for ii in 1..(*histaxis as usize) {
                minin[ii] = minin[0];
                maxin[ii] = maxin[0];
                binsizein[ii] = binsizein[0];
                // Read first to avoid borrow checker issues. This copies the data, so not ideal.
                let minname0 = minname[0];
                let maxname0 = maxname[0];
                let binname0 = binname[0];
                strcpy_safe(&mut minname[ii], &minname0);
                strcpy_safe(&mut maxname[ii], &maxname0);
                strcpy_safe(&mut binname[ii], &binname0);
            }

            while ptr[0] == bb(b' ') {
                /* skip over blanks */
                ptr = &ptr[1..];
            }

            if ptr[0] == bb(b';') {
                break 'getweight; /* a weighting factor is specified */
            }

            if ptr[0] != 0
            /* must have reached end of string */
            {
                ffpmsg_str("illegal syntax after binning range specification in URL:");
                ffpmsg_slice(binspec);
                *status = URL_PARSE_ERROR;
                return *status;
            }

            return *status;
        } /* end of case with list of column names in ( )  */

        /* if we've reached this point, then the binning specification */
        /* must be of the form: XCOL = min:max:binsize, YCOL = ...     */
        /* where the column name followed by '=' are optional.         */
        /* If the column name is not specified, then use the default name */

        let mut ii: usize = 0;
        while ii < 4 {
            /* allow up to 4 histogram dimensions */
            exprend[ii] = core::ptr::null_mut();
            exprbeg[ii] = core::ptr::null_mut();
            let mut exprbeg_idx: usize = 0;
            let mut exprend_idx: usize = 0;

            // Save ptr position BEFORE calling ffbinre
            let ptr_before_ffbinre = ptr;

            ffbinre(
                &mut ptr,
                &mut colname[ii],
                Some(&mut exprbeg_idx),
                Some(&mut exprend_idx),
                &mut minin[ii],
                &mut maxin[ii],
                &mut binsizein[ii],
                &mut minname[ii],
                &mut maxname[ii],
                &mut binname[ii],
                status,
            );

            /* Check for expressions */
            if exprbeg_idx != 0 {
                has_exprs = true;
                // Convert indices to pointers (relative to binspec, not ptr)
                // Calculate offset from binspec to ptr position BEFORE ffbinre was called
                let offset_from_binspec =
                    unsafe { ptr_before_ffbinre.as_ptr().offset_from(binspec.as_ptr()) as usize };

                exprbeg[ii] = binspec
                    .as_ptr()
                    .wrapping_add(offset_from_binspec + exprbeg_idx)
                    as *mut c_char;
                exprend[ii] = binspec
                    .as_ptr()
                    .wrapping_add(offset_from_binspec + exprend_idx)
                    as *mut c_char;
            }

            if *status > 0 {
                ffpmsg_str("illegal syntax in binning range specification in URL:");
                ffpmsg_slice(binspec);
                return *status;
            }

            if ptr[0] == 0 || ptr[0] == bb(b';') {
                break; /* reached the end of the string */
            }

            if ptr[0] == bb(b' ') {
                while ptr[0] == bb(b' ') {
                    /* skip over blanks */
                    ptr = &ptr[1..];
                }

                if ptr[0] == 0 || ptr[0] == bb(b';') {
                    break; /* reached the end of the string */
                }

                if ptr[0] == bb(b',') {
                    ptr = &ptr[1..]; /* comma separates the next column specification */
                }
            } else if ptr[0] == bb(b',') {
                ptr = &ptr[1..]; /* comma separates the next column specification */
            } else {
                ffpmsg_str("illegal characters following binning specification in URL:");
                ffpmsg_slice(binspec);
                *status = URL_PARSE_ERROR;
                return *status;
            }

            ii += 1;
        }

        if ii == 4 {
            /* there are yet more characters in the string */
            ffpmsg_str("illegal binning specification in URL:");
            ffpmsg_str("apparently greater than 4 histogram dimensions");
            ffpmsg_slice(binspec);
            *status = URL_PARSE_ERROR;
            return *status;
        } else {
            *histaxis = (ii + 1) as c_int;
        }

        /* special case: if a single number was entered it should be      */
        /* interpreted as the binning factor for the default X and Y axes */

        if *histaxis == 1
            && colname[0][0] == 0
            && minin[0] == DOUBLENULLVALUE
            && maxin[0] == DOUBLENULLVALUE
        {
            *histaxis = 2;
            binsizein[1] = binsizein[0];
        }

        break;
    } // Loop specifically to break to `getweight`

    // getweight:
    if ptr[0] == bb(b';') {
        /* looks like a weighting factor is given */
        ptr = &ptr[1..];

        while ptr[0] == bb(b' ') {
            /* skip over blanks */
            ptr = &ptr[1..];
        }

        *recip = 0;
        if ptr[0] == bb(b'/') {
            *recip = 1; /* the reciprocal of the weight is entered */
            ptr = &ptr[1..];

            while ptr[0] == bb(b' ') {
                /* skip over blanks */
                ptr = &ptr[1..];
            }
        }

        /* parse the weight as though it were a binrange. */
        /* either a column name or a numerical value will be returned */
        exprend[4] = core::ptr::null_mut();
        exprbeg[4] = core::ptr::null_mut();
        let mut exprbeg_idx: usize = 0;
        let mut exprend_idx: usize = 0;
        let mut dummy2: f64 = 0.0;
        let mut tmpname2: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
        let mut tmpname3: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
        ffbinre(
            &mut ptr,
            wtname,
            Some(&mut exprbeg_idx),
            Some(&mut exprend_idx),
            &mut dummy,
            &mut dummy2,
            weight,
            &mut tmpname,
            &mut tmpname2,
            &mut tmpname3,
            status,
        );
        if exprbeg_idx != 0 {
            has_exprs = true;
            // Convert indices to pointers
            exprbeg[4] = ptr.as_ptr().wrapping_add(exprbeg_idx) as *mut c_char;
            exprend[4] = ptr.as_ptr().wrapping_add(exprend_idx) as *mut c_char;
        }

        if *status > 0 {
            ffpmsg_str("illegal binning weight specification in URL:");
            ffpmsg_slice(binspec);
            return *status;
        }

        /* creat a float datatype histogram by default, if weight */
        /* factor is not = 1.0  */

        if (defaulttype != 0 && *weight != 1.0)
            || (defaulttype != 0 && wtname[0] != 0)
            || (defaulttype != 0 && !exprbeg[4].is_null())
        {
            *imagetype = TFLOAT;
        }
    }

    while ptr[0] == bb(b' ') {
        /* skip over blanks */
        ptr = &ptr[1..];
    }

    if ptr[0] != 0
    /* should have reached the end of string */
    {
        ffpmsg_str("illegal syntax after binning weight specification in URL:");
        ffpmsg_slice(binspec);
        *status = URL_PARSE_ERROR;
    }

    /* If we found expressions, this is where we accumulate them into
    something to be returned to the caller.  The start and end of
    each expression will be found in exprbeg[] and exprend[], with
    the 5th entry being the weight expression if any */
    if has_exprs && let Some(exprs_arr) = exprs {
        // For each dimension (0-3) and the weight (4), extract the expression string
        for ii in 0..5 {
            let expr_ptr: *const c_char;
            let expr_len: usize;

            if exprbeg[ii].is_null() || exprend[ii].is_null() {
                // No expression for this dimension
                expr_len = 0;
                expr_ptr = core::ptr::null();
            } else {
                // exprbeg and exprend are pointers into the original binspec string
                unsafe {
                    expr_ptr = exprbeg[ii];
                    // Calculate length by pointer subtraction (exprend - exprbeg)
                    expr_len = exprend[ii].offset_from(exprbeg[ii]) as usize;
                }
            }

            // Allocate and copy the expression string
            if expr_len > 0 {
                let mut expr_vec = vec![0 as c_char; expr_len + 1]; // +1 for null terminator
                unsafe {
                    core::ptr::copy_nonoverlapping(expr_ptr, expr_vec.as_mut_ptr(), expr_len);
                }
                expr_vec[expr_len] = 0; // Null terminator

                exprs_arr[ii] = expr_vec.into_boxed_slice();
            } else {
                // Empty expression
                exprs_arr[ii] = vec![0 as c_char].into_boxed_slice();
            }
        }
    }

    *status
}

/*--------------------------------------------------------------------------*/
/// Parse non-extended expression, but otherwise the same as ffbinse()
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffbins(
    binspec: *const c_char,                  /* I - binning specification */
    imagetype: *mut c_int,                   /* O - image type, TINT or TSHORT */
    histaxis: *mut c_int,                    /* O - no. of axes in the histogram */
    colname: *mut [[c_char; FLEN_VALUE]; 4], /* column name for axis */
    minin: *mut f64,                         /* minimum value for each axis */
    maxin: *mut f64,                         /* maximum value for each axis */
    binsizein: *mut f64,                     /* size of bins on each axis */
    minname: *mut [[c_char; FLEN_VALUE]; 4], /* keyword name for min */
    maxname: *mut [[c_char; FLEN_VALUE]; 4], /* keyword name for max */
    binname: *mut [[c_char; FLEN_VALUE]; 4], /* keyword name for binsize */
    weight: *mut f64,                        /* weighting factor          */
    wtname: *mut [c_char; FLEN_VALUE],       /* keyword or column name for weight */
    recip: *mut c_int,                       /* the reciprocal of the weight? */
    status: *mut c_int,
) -> c_int {
    unsafe {
        raw_to_slice!(binspec);

        let imagetype = imagetype.as_mut().expect(NULL_MSG);
        let histaxis = histaxis.as_mut().expect(NULL_MSG);
        let colname = colname.as_mut().expect(NULL_MSG);
        let minin: &mut [f64; 4] = slice::from_raw_parts_mut(minin, 4).try_into().unwrap();
        let maxin: &mut [f64; 4] = slice::from_raw_parts_mut(maxin, 4).try_into().unwrap();
        let binsizein: &mut [f64; 4] = slice::from_raw_parts_mut(binsizein, 4).try_into().unwrap();
        let minname = minname.as_mut().expect(NULL_MSG);
        let maxname = maxname.as_mut().expect(NULL_MSG);
        let binname = binname.as_mut().expect(NULL_MSG);
        let weight = weight.as_mut().expect(NULL_MSG);
        let wtname = wtname.as_mut().expect(NULL_MSG);

        let recip = recip.as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);

        ffbins_safe(
            binspec, imagetype, histaxis, colname, minin, maxin, binsizein, minname, maxname,
            binname, weight, wtname, recip, status,
        )
    }
}

/*--------------------------------------------------------------------------*/
/// Parse non-extended expression, but otherwise the same as ffbinse()
pub fn ffbins_safe(
    binspec: &[c_char],                      /* I - binning specification */
    imagetype: &mut c_int,                   /* O - image type, TINT or TSHORT */
    histaxis: &mut c_int,                    /* O - no. of axes in the histogram */
    colname: &mut [[c_char; FLEN_VALUE]; 4], /* column name for axis */
    minin: &mut [f64; 4],                    /* minimum value for each axis */
    maxin: &mut [f64; 4],                    /* maximum value for each axis */
    binsizein: &mut [f64; 4],                /* size of bins on each axis */
    minname: &mut [[c_char; FLEN_VALUE]; 4], /* keyword name for min */
    maxname: &mut [[c_char; FLEN_VALUE]; 4], /* keyword name for max */
    binname: &mut [[c_char; FLEN_VALUE]; 4], /* keyword name for binsize */
    weight: &mut f64,                        /* weighting factor          */
    wtcol: &mut [c_char; FLEN_VALUE],        /* keyword or column name for weight */
    recip: &mut c_int,                       /* the reciprocal of the weight? */
    status: &mut c_int,
) -> c_int {
    ffbinse(
        binspec, imagetype, histaxis, colname, minin, maxin, binsizein, minname, maxname, binname,
        weight, wtcol, recip, None, /* No exprs pointer */
        status,
    )
}

/*--------------------------------------------------------------------------*/
/// Parse the input binning range specification string, returning
/// the column name, histogram min and max values, and bin size.
///
/// This is the "extended" binning syntax that allows for an expression
/// of the form XCOL(expr).  The expression must be enclosed in parentheses.
///
/// The start and end of the expression are returned in *exprbeg and *exprend.
/// If exprbeg and exprend are null pointers then the expression is forbidden.
pub(crate) fn ffbinre(
    ptr: &mut &[c_char],
    colname: &mut [c_char],
    mut exprbeg: Option<&mut usize>,
    mut exprend: Option<&mut usize>,
    minin: &mut f64,
    maxin: &mut f64,
    binsizein: &mut f64,
    minname: &mut [c_char],
    maxname: &mut [c_char],
    binname: &mut [c_char],
    status: &mut c_int,
) -> c_int {
    let mut slen: c_int;
    let mut isanumber: c_int = 0;
    let mut token: Option<Vec<c_char>> = None;

    if *status > 0 {
        return *status;
    }

    // Initialize expression indices to 0 (no expression)
    if let Some(exprbeg_ref) = exprbeg.as_deref_mut() {
        *exprbeg_ref = 0;
    }
    if let Some(exprend_ref) = exprend.as_deref_mut() {
        *exprend_ref = 0;
    }

    let mut ptr_index = 0;
    slen = fits_get_token2_safe(
        ptr,
        &mut ptr_index,
        cs!(c" ,=:;("),
        &mut token,
        Some(&mut isanumber),
        status,
    ); /* get 1st token */

    if (*status != 0)
        || (slen == 0
            && ((*ptr)[ptr_index] == 0
                || (*ptr)[ptr_index] == bb(b',')
                || (*ptr)[ptr_index] == bb(b';')))
    {
        *ptr = &(*ptr)[ptr_index..];
        return *status; /* a null range string */
    }

    if isanumber == 0 && (*ptr)[ptr_index] != bb(b':') {
        /* this looks like the column name */

        /* Check for case where col name string is empty but '='
        is still there (indicating a following specification string).
        Musn't enter this block as token would not have been allocated. */
        if let Some(ref token_vec) = token {
            if strlen_safe(token_vec) > FLEN_VALUE - 1 {
                ffpmsg_str("column name too long (ffbinr)");

                *status = PARSE_SYNTAX_ERR;
                return *status;
            }
            if token_vec[0] == bb(b'#') && isdigit_safe(token_vec[1]) {
                /* omit the leading '#' in the column number */
                strcpy_safe(colname, &token_vec[1..]);
            } else {
                strcpy_safe(colname, token_vec);
            }

            token = None;
        }
        while (*ptr)[ptr_index] == bb(b' ') {
            /* skip over blanks */
            ptr_index += 1;
        }

        /* An optional expression of the form XCOL(expr) is allowed here, but only if exprbeg and exprend are non-null */
        if (*ptr)[ptr_index] == bb(b'(')
            && let Some(exprbeg) = exprbeg
            && let Some(exprend) = exprend
        {
            *exprbeg = ptr_index;
            let tmp = fits_find_match_delim(&ptr[ptr_index + 1..], bb(b')'));
            match tmp {
                None => {
                    /* find ')' */
                    ffpmsg_str("bin expression syntax error (ffbinr)");
                    *status = PARSE_SYNTAX_ERR;
                    return *status;
                }
                Some(tmp) => {
                    *exprend = ptr_index + tmp + 1; /* +1 to include the closing paren */
                    ptr_index += tmp + 1; /* Advance pointer past delimiter */
                }
            };
        }

        while (*ptr)[ptr_index] == bb(b' ') {
            ptr_index += 1; /* skip over more possible blanks */
        }

        if (*ptr)[ptr_index] != bb(b'=') {
            *ptr = &(*ptr)[ptr_index..];
            return *status; /* reached the end */
        }

        ptr_index += 1; /* skip over the = sign */

        while (*ptr)[ptr_index] == bb(b' ') {
            /* skip over blanks */
            ptr_index += 1;
        }

        /* get specification info */
        let mut ptr_index2 = 0;
        let ptr_sub = &(*ptr)[ptr_index..];
        slen = fits_get_token2_safe(
            ptr_sub,
            &mut ptr_index2,
            cs!(c" ,:;"),
            &mut token,
            Some(&mut isanumber),
            status,
        );
        ptr_index += ptr_index2;

        if *status != 0 {
            *ptr = &(*ptr)[ptr_index..];
            return *status;
        }
    }

    if (*ptr)[ptr_index] != bb(b':') {
        /* This is the first token, and since it is not followed by
        a ':' this must be the binsize token. Or it could be empty. */
        if let Some(ref token_vec) = token {
            if isanumber == 0 {
                if strlen_safe(token_vec) > FLEN_VALUE - 1 {
                    ffpmsg_str("binname too long (ffbinr)");

                    *status = PARSE_SYNTAX_ERR;
                    return *status;
                }
                strcpy_safe(binname, token_vec);
            } else {
                let mut endp = 0;
                *binsizein = strtod_safe(token_vec, &mut endp);
            }
        }

        *ptr = &(*ptr)[ptr_index..];
        return *status; /* reached the end */
    } else {
        /* the token contains the min value */
        if slen > 0
            && let Some(ref token_vec) = token
        {
            if isanumber == 0 {
                if strlen_safe(token_vec) > FLEN_VALUE - 1 {
                    ffpmsg_str("minname too long (ffbinr)");

                    *status = PARSE_SYNTAX_ERR;
                    return *status;
                }
                strcpy_safe(minname, token_vec);
            } else {
                let mut endp = 0;
                *minin = strtod_safe(token_vec, &mut endp);
            }

            token = None;
        }
    }

    ptr_index += 1; /* skip the colon between the min and max values */
    let mut ptr_index2 = 0;
    let ptr_sub = &(*ptr)[ptr_index..];
    slen = fits_get_token2_safe(
        ptr_sub,
        &mut ptr_index2,
        cs!(c" ,:;"),
        &mut token,
        Some(&mut isanumber),
        status,
    ); /* get token */
    ptr_index += ptr_index2;
    if *status != 0 {
        *ptr = &(*ptr)[ptr_index..];
        return *status;
    }

    /* the token contains the max value */
    if slen > 0
        && let Some(ref token_vec) = token
    {
        if isanumber == 0 {
            if strlen_safe(token_vec) > FLEN_VALUE - 1 {
                ffpmsg_str("maxname too long (ffbinr)");

                *status = PARSE_SYNTAX_ERR;
                return *status;
            }
            strcpy_safe(maxname, token_vec);
        } else {
            let mut endp = 0;
            *maxin = strtod_safe(token_vec, &mut endp);
        }

        token = None;
    }

    if (*ptr)[ptr_index] != bb(b':') {
        *ptr = &(*ptr)[ptr_index..];
        return *status; /* reached the end; no binsize token */
    }

    ptr_index += 1; /* skip the colon between the max and binsize values */
    let mut ptr_index2 = 0;
    let ptr_sub = &(*ptr)[ptr_index..];
    slen = fits_get_token2_safe(
        ptr_sub,
        &mut ptr_index2,
        cs!(c" ,:;"),
        &mut token,
        Some(&mut isanumber),
        status,
    ); /* get token */

    ptr_index += ptr_index2;

    if *status != 0 {
        *ptr = &(*ptr)[ptr_index..];
        return *status;
    }

    /* the token contains the binsize value */
    if slen > 0
        && let Some(ref token_vec) = token
    {
        if isanumber == 0 {
            if strlen_safe(token_vec) > FLEN_VALUE - 1 {
                ffpmsg_str("binname too long (ffbinr)");

                *status = PARSE_SYNTAX_ERR;
                return *status;
            }
            strcpy_safe(binname, token_vec);
        } else {
            let mut endp = 0;
            *binsizein = strtod_safe(token_vec, &mut endp);
        }
    }

    /* Update the pointer to point past what we've parsed */
    *ptr = &(*ptr)[ptr_index..];

    *status
}
/*--------------------------------------------------------------------------*/
/// Parse the input binning range specification string, returning
/// the column name, histogram min and max values, and bin size.
///
/// This is the non-extended version of the parser which disallows
/// binning expressions.  Only column names are allowed.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffbinr(
    ptr: *mut *const c_char,
    colname: *mut c_char,
    minin: *mut f64,
    maxin: *mut f64,
    binsizein: *mut f64,
    minname: *mut c_char,
    maxname: *mut c_char,
    binname: *mut c_char,
    status: *mut c_int,
) -> c_int {
    unsafe {
        let minin = minin.as_mut().expect(NULL_MSG);
        let maxin = maxin.as_mut().expect(NULL_MSG);
        let binsizein = binsizein.as_mut().expect(NULL_MSG);

        let colname = slice::from_raw_parts_mut(colname, FLEN_VALUE)
            .try_into()
            .unwrap();
        let minname = slice::from_raw_parts_mut(minname, FLEN_VALUE)
            .try_into()
            .unwrap();
        let maxname = slice::from_raw_parts_mut(maxname, FLEN_VALUE)
            .try_into()
            .unwrap();
        let binname = slice::from_raw_parts_mut(binname, FLEN_VALUE)
            .try_into()
            .unwrap();

        let status = status.as_mut().expect(NULL_MSG);

        let ptr = ptr.as_mut().expect(NULL_MSG);

        let ptr_slice = slice::from_raw_parts(*ptr, strlen(*ptr));

        ffbinr_safe(
            &mut &ptr_slice[..],
            colname,
            minin,
            maxin,
            binsizein,
            minname,
            maxname,
            binname,
            status,
        )
    }
}

/*--------------------------------------------------------------------------*/
/// Parse the input binning range specification string, returning
/// the column name, histogram min and max values, and bin size.
///
/// This is the non-extended version of the parser which disallows
/// binning expressions.  Only column names are allowed.
pub fn ffbinr_safe(
    ptr: &mut &[c_char],
    colname: &mut [c_char; FLEN_VALUE],
    minin: &mut f64,
    maxin: &mut f64,
    binsizein: &mut f64,
    minname: &mut [c_char; FLEN_VALUE],
    maxname: &mut [c_char; FLEN_VALUE],
    binname: &mut [c_char; FLEN_VALUE],
    status: &mut c_int,
) -> c_int {
    ffbinre(
        ptr, colname, None, None, minin, maxin, binsizein, minname, maxname, binname, status,
    )
}

/*--------------------------------------------------------------------------*/
pub(crate) fn ffhist2e(
    fptr: &mut Option<Box<fitsfile>>, /* IO - pointer to table with X and Y cols, on output, points to histogram image    */
    outfile: &[c_char],               /* I - name for the output histogram file      */
    imagetype: c_int,                 /* I - datatype for image: TINT, TSHORT, etc   */
    naxis: c_int,                     /* I - number of axes in the histogram image   */
    colname: &[[c_char; FLEN_VALUE]; 4], /* I - column names               */
    colexpr: Option<&[Option<&[c_char]>; 4]>, /* I - optionally, expression intead of colum  */
    minin: &[f64],                    /* I - minimum histogram value, for each axis */
    maxin: &[f64],                    /* I - maximum histogram value, for each axis */
    binsizein: &[f64],                /* I - bin size along each axis               */
    minname: &[[c_char; FLEN_VALUE]; 4], /* I - optional keywords for min    */
    maxname: &[[c_char; FLEN_VALUE]; 4], /* I - optional keywords for max    */
    binname: &[[c_char; FLEN_VALUE]; 4], /* I - optional keywords for binsize */
    weightin: f64,                    /* I - binning weighting factor          */
    wtcol: &[c_char; FLEN_VALUE],     /* I - optional keyword or col for weight*/
    wtexpr: Option<&[c_char]>,        /* I - optionally, weight expression     */
    recip: c_int,                     /* I - use reciprocal of the weight?     */
    selectrow: Option<&[c_char]>,     /* I - optional array (length = no. of   */
    /* rows in the table).  If the element is true */
    /* then the corresponding row of the table will*/
    /* be included in the histogram, otherwise the */
    /* row will be skipped.  Ingnored if *selectrow*/
    /* is equal to NULL.                           */
    status: &mut c_int,
) -> c_int {
    let mut histptr = None;
    let mut bitpix: c_int = 0;
    let mut colnum: [c_int; 4] = [0; 4];
    let mut wtcolnum: c_int = 0;
    let mut haxes: [c_long; 4] = [0; 4];
    let mut amin: [f64; 4] = [0.; 4];
    let mut amax: [f64; 4] = [0.; 4];
    let mut binsize: [f64; 4] = [0.; 4];
    let mut weight: f64 = 0.;
    let numIterCols: c_int = 0;
    let mut datatypes: [c_int; 4] = [0; 4];
    let mut wtdatatype: c_int = 0;
    let repeat: *mut c_long = core::ptr::null_mut();
    let mut wtrepeat: c_long = 0;
    let errmsg: [c_char; FLEN_ERRMSG] = [0; FLEN_ERRMSG];
    let mut vectorRepeat: c_long = 0;

    if *status > 0 {
        return *status;
    }

    if naxis > 4 {
        ffpmsg_str("histogram has more than 4 dimensions");
        *status = BAD_DIMEN;
        return *status;
    }

    if let Some(fptr_inner) = fptr.as_deref_mut() {
        /* reset position to the correct HDU if necessary */
        if fptr_inner.HDUposition != (fptr_inner.Fptr).curhdu {
            let hdu_pos = (fptr_inner.HDUposition) + 1;
            ffmahd_safe(fptr_inner, hdu_pos, None, status);
        }

        if imagetype == TBYTE {
            bitpix = BYTE_IMG;
        } else if imagetype == TSHORT {
            bitpix = SHORT_IMG;
        } else if imagetype == TINT {
            bitpix = LONG_IMG;
        } else if imagetype == TFLOAT {
            bitpix = FLOAT_IMG;
        } else if imagetype == TDOUBLE {
            bitpix = DOUBLE_IMG;
        } else {
            *status = BAD_DATATYPE;
            return *status;
        }

        /*    Calculate the binning parameters:    */
        /*   columm numbers, axes length, min values,  max values, and binsizes.  */

        // Need to make local mutable copies since fits_calc_binningde may modify them
        let mut minin_mut: [f64; 4] = [minin[0], minin[1], minin[2], minin[3]];
        let mut maxin_mut: [f64; 4] = [maxin[0], maxin[1], maxin[2], maxin[3]];
        let mut binsizein_mut: [f64; 4] = [binsizein[0], binsizein[1], binsizein[2], binsizein[3]];

        if fits_calc_binningde(
            fptr_inner,
            naxis,
            &mut colname.clone(),
            colexpr,
            &mut minin_mut,
            &mut maxin_mut,
            &mut binsizein_mut,
            minname,
            maxname,
            binname,
            &mut colnum,
            Some(&mut datatypes),
            &mut haxes,
            &mut amin,
            &mut amax,
            &mut binsize,
            Some(&mut vectorRepeat),
            status,
        ) > 0
        {
            ffpmsg_str("failed to determine binning parameters");
            return *status;
        }

        /* get the histogramming weighting factor, if any */
        if wtcol[0] != 0 {
            /* first, look for a keyword with the weight value */
            if ffgky_safe(
                fptr_inner,
                KeywordDatatypeMut::TDOUBLE(&mut weight),
                wtcol,
                None,
                status,
            ) == 0
            {
                /* Data type if keyword was found */
                wtdatatype = TDOUBLE;
                wtrepeat = 1;
            } else {
                /* not a keyword, so look for column with this name */
                *status = 0;

                /* get the column number in the table */
                if ffgcno_safe(fptr_inner, CASEINSEN as c_int, wtcol, &mut wtcolnum, status) > 0 {
                    ffpmsg_str("keyword or column for histogram weights doesn't exist: ");
                    ffpmsg_slice(wtcol);
                    return *status;
                }

                /* get the datatype of the column */
                fits_get_eqcoltype(
                    fptr_inner,
                    wtcolnum,
                    Some(&mut wtdatatype),
                    Some(&mut wtrepeat),
                    None,
                    status,
                );

                weight = DOUBLENULLVALUE;
            }
        } else if let Some(wtexpr) = wtexpr
            && wtexpr[0] != 0
        {
            /* A weighting expression - always TDOUBLE */

            /* Initialize the parser so that we can determine the datatype
            of the returned type as well as the vector dimensions.  The
            parsers is kept allocated so we can assemble an iterator that
            uses it below.
             */

            let mut naxis1: c_int = 0;
            let mut nelem: c_long = 0;
            let mut naxes: [c_long; MAXDIMS as usize] = [0; MAXDIMS as usize];
            let mut lParse: ParseData = ParseData::default();

            ffiprs(
                fptr_inner,
                0,
                wtexpr,
                MAXDIMS,
                &mut wtdatatype,
                &mut nelem,
                &mut naxis1,
                &mut naxes,
                &mut lParse,
                status,
            );
            ffcprs(&mut lParse);
            if nelem < 0 {
                nelem = 1; /* If it's a constant expression */
            }

            weight = DOUBLENULLVALUE;
            wtrepeat = nelem;
            wtdatatype = wtdatatype;
        } else {
            weight = weightin;
            wtrepeat = vectorRepeat;
            wtdatatype = TDOUBLE;
        }

        /* Make sure weighting column is not an un-binnable data type */
        if wtdatatype < 0 || wtdatatype == TSTRING || wtdatatype == TBIT || wtdatatype == TLOGICAL {
            ffpmsg_str("Invalid datatype for bin weighting factor");
            *status = BAD_DATATYPE;
            return *status;
        }

        /* And dimensions of weighting must agree with input column data */
        if wtrepeat != vectorRepeat {
            ffpmsg_str("Vector dimensions of weighting do not agree with binning columns");
            *status = BAD_DIMEN;
            return *status;
        }

        if weight <= 0. && weight != DOUBLENULLVALUE {
            ffpmsg_str("Illegal histogramming weighting factor <= 0.");
            *status = URL_PARSE_ERROR;
            return *status;
        }

        if recip != 0 && weight != DOUBLENULLVALUE {
            /* take reciprocal of weight */
            weight = 1.0 / weight;
        }

        /* size of histogram is now known, so create temp output file */
        if fits_create_file(&mut histptr, outfile, status) > 0 {
            ffpmsg_str("failed to create temp output file for histogram");
            return *status;
        }

        let histptr = histptr.as_mut().unwrap();

        /* create output FITS image HDU */
        if ffcrim_safe(histptr, bitpix, naxis, &haxes, status) > 0 {
            ffpmsg_str("failed to create output histogram FITS image");
            return *status;
        }

        /* copy header keywords, converting pixel list WCS keywords to image WCS form */
        if fits_copy_pixlist2image_safe(fptr_inner, histptr, 9, naxis, &colnum, status) > 0 {
            ffpmsg_str("failed to copy pixel list keywords to new histogram header");
            return *status;
        }

        /* if the table columns have no WCS keywords, then write default keywords */
        fits_write_keys_histoe(
            fptr_inner,
            histptr,
            naxis,
            &colnum,
            Some(colname),
            colexpr,
            status,
        );

        /* update the WCS keywords for the ref. pixel location, and pixel size */
        fits_rebin_wcsd_safe(histptr, naxis, &mut amin, &mut binsize, status);

        /* now compute the output image by binning the column values */
        if fits_make_histde(
            fptr_inner,
            histptr,
            Some(&mut datatypes),
            bitpix,
            naxis,
            &haxes,
            &colnum,
            colexpr,
            &amin,
            &amax,
            &binsize,
            weight,
            wtcolnum,
            wtexpr,
            recip,
            selectrow,
            status,
        ) > 0
        {
            ffpmsg_str("failed to calculate new histogram values");
            return *status;
        }
    }
    /* finally, close the original file and return ptr to the new image */
    let tmp = fptr.take().unwrap();
    ffclos_safe(tmp, status);

    *fptr = histptr;

    // cleanup:
    *status
}

/*--------------------------------------------------------------------------*/
/// Non-extended-syntax version of ffhist2e()
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffhist2(
    fptr: *mut Option<Box<fitsfile>>, /* IO - pointer to table with X and Y cols;    */
    /*     on output, points to histogram image    */
    outfile: *const c_char, /* I - name for the output histogram file      */
    imagetype: c_int,       /* I - datatype for image: TINT, TSHORT, etc   */
    naxis: c_int,           /* I - number of axes in the histogram image   */
    colname: *const [[c_char; FLEN_VALUE]; 4], /* I - column names               */
    minin: *const f64,      /* I - minimum histogram value, for each axis */
    maxin: *const f64,      /* I - maximum histogram value, for each axis */
    binsizein: *const f64,  /* I - bin size along each axis               */
    minname: *const [[c_char; FLEN_VALUE]; 4], /* I - optional keywords for min    */
    maxname: *const [[c_char; FLEN_VALUE]; 4], /* I - optional keywords for max    */
    binname: *const [[c_char; FLEN_VALUE]; 4], /* I - optional keywords for binsize */
    weightin: f64,          /* I - binning weighting factor          */
    wtcol: *const [c_char; FLEN_VALUE], /* I - optional keyword or col for weight*/
    recip: c_int,           /* I - use reciprocal of the weight?     */
    selectrow: *const c_char, /* I - optional array (length = no. of   */
    /* rows in the table).  If the element is true */
    /* then the corresponding row of the table will*/
    /* be included in the histogram, otherwise the */
    /* row will be skipped.  Ingnored if *selectrow*/
    /* is equal to NULL.                           */
    status: *mut c_int,
) -> c_int {
    unsafe {
        let fptr = (fptr).as_mut().expect(NULL_MSG);
        raw_to_slice!(outfile);
        let colname = colname.as_ref().expect(NULL_MSG);
        let minin = slice::from_raw_parts(minin, naxis as usize);
        let maxin = slice::from_raw_parts(maxin, naxis as usize);
        let binsizein = slice::from_raw_parts(binsizein, naxis as usize);
        let minname = minname.as_ref().expect(NULL_MSG);
        let maxname = maxname.as_ref().expect(NULL_MSG);
        let binname = binname.as_ref().expect(NULL_MSG);
        let wtcol = wtcol.as_ref().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);

        let mut nrows = 0;

        let f = fptr.as_mut().unwrap();
        ffgnrw_safe(f, &mut nrows, status); /* no. of rows */

        let selectrow = if selectrow.is_null() {
            None
        } else {
            Some(slice::from_raw_parts(selectrow, nrows as usize))
        };

        ffhist2_safe(
            fptr, outfile, imagetype, naxis, colname, minin, maxin, binsizein, minname, maxname,
            binname, weightin, wtcol, recip, selectrow, status,
        )
    }
}

/*--------------------------------------------------------------------------*/
/// Non-extended-syntax version of ffhist2e()
pub fn ffhist2_safe(
    fptr: &mut Option<Box<fitsfile>>, /* IO - pointer to table with X and Y cols;    */
    /*     on output, points to histogram image    */
    outfile: &[c_char], /* I - name for the output histogram file      */
    imagetype: c_int,   /* I - datatype for image: TINT, TSHORT, etc   */
    naxis: c_int,       /* I - number of axes in the histogram image   */
    colname: &[[c_char; FLEN_VALUE]; 4], /* I - column names               */
    minin: &[f64],      /* I - minimum histogram value, for each axis */
    maxin: &[f64],      /* I - maximum histogram value, for each axis */
    binsizein: &[f64],  /* I - bin size along each axis               */
    minname: &[[c_char; FLEN_VALUE]; 4], /* I - optional keywords for min    */
    maxname: &[[c_char; FLEN_VALUE]; 4], /* I - optional keywords for max    */
    binname: &[[c_char; FLEN_VALUE]; 4], /* I - optional keywords for binsize */
    weightin: f64,      /* I - binning weighting factor          */
    wtcol: &[c_char; FLEN_VALUE], /* I - optional keyword or col for weight*/
    recip: c_int,       /* I - use reciprocal of the weight?     */
    selectrow: Option<&[c_char]>, /* I - optional array (length = no. of   */
    /* rows in the table).  If the element is true */
    /* then the corresponding row of the table will*/
    /* be included in the histogram, otherwise the */
    /* row will be skipped.  Ingnored if *selectrow*/
    /* is equal to NULL.                           */
    status: &mut c_int,
) -> c_int {
    ffhist2e(
        fptr, outfile, imagetype, naxis, colname, None, minin, maxin, binsizein, minname, maxname,
        binname, weightin, wtcol, None, recip, selectrow, status,
    )
}

/*--------------------------------------------------------------------------*/
/// ffhist3: same as ffhist2, but does not close the original file
///  and/or replace the original file pointer
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffhist3(
    fptr: *const *mut fitsfile, /* I - pointer to table with X and Y cols;    */
    outfile: *const c_char,     /* I - name for the output histogram file      */
    imagetype: c_int,           /* I - datatype for image: TINT, TSHORT, etc   */
    naxis: c_int,               /* I - number of axes in the histogram image   */
    colname: *const [[c_char; FLEN_VALUE]; 4], /* I - column names               */
    minin: *const f64,          /* I - minimum histogram value, for each axis */
    maxin: *const f64,          /* I - maximum histogram value, for each axis */
    binsizein: *const f64,      /* I - bin size along each axis               */
    minname: *const [[c_char; FLEN_VALUE]; 4], /* I - optional keywords for min    */
    maxname: *const [[c_char; FLEN_VALUE]; 4], /* I - optional keywords for max    */
    binname: *const [[c_char; FLEN_VALUE]; 4], /* I - optional keywords for binsize */
    weightin: f64,              /* I - binning weighting factor          */
    wtcol: *const [c_char; FLEN_VALUE], /* I - optional keyword or col for weight*/
    recip: c_int,               /* I - use reciprocal of the weight?     */
    selectrow: *const c_char,   /* I - optional array (length = no. of   */
    /* rows in the table).  If the element is true */
    /* then the corresponding row of the table will*/
    /* be included in the histogram, otherwise the */
    /* row will be skipped.  Ingnored if *selectrow*/
    /* is equal to NULL.                           */
    status: *mut c_int,
) -> *mut fitsfile {
    unsafe {
        let fptr = &mut ((*fptr).as_mut()).expect(NULL_MSG);
        raw_to_slice!(outfile);
        let colname = colname.as_ref().unwrap();
        let minin = slice::from_raw_parts(minin, naxis as usize);
        let maxin = slice::from_raw_parts(maxin, naxis as usize);
        let binsizein = slice::from_raw_parts(binsizein, naxis as usize);
        let minname = minname.as_ref().unwrap();
        let maxname = maxname.as_ref().unwrap();
        let binname = binname.as_ref().unwrap();
        let wtcol = wtcol.as_ref().unwrap();
        let status = status.as_mut().expect(NULL_MSG);

        let mut nrows = 0;
        ffgnrw_safe(fptr, &mut nrows, status); /* no. of rows */

        let selectrow = if selectrow.is_null() {
            None
        } else {
            Some(slice::from_raw_parts(selectrow, nrows as usize))
        };

        let res = ffhist3_safe(
            fptr, outfile, imagetype, naxis, colname, minin, maxin, binsizein, minname, maxname,
            binname, weightin, wtcol, recip, selectrow, status,
        );

        match res {
            Some(boxed) => Box::into_raw(boxed),
            None => core::ptr::null_mut(),
        }
    }
}

/*--------------------------------------------------------------------------*/
/// ffhist3: same as ffhist2, but does not close the original file
///  and/or replace the original file pointer
pub fn ffhist3_safe(
    fptr: &mut fitsfile, /* IO - pointer to table with X and Y cols;    */
    outfile: &[c_char],  /* I - name for the output histogram file      */
    imagetype: c_int,    /* I - datatype for image: TINT, TSHORT, etc   */
    naxis: c_int,        /* I - number of axes in the histogram image   */
    colname: &[[c_char; FLEN_VALUE]; 4], /* I - column names               */
    minin: &[f64],       /* I - minimum histogram value, for each axis */
    maxin: &[f64],       /* I - maximum histogram value, for each axis */
    binsizein: &[f64],   /* I - bin size along each axis               */
    minname: &[[c_char; FLEN_VALUE]; 4], /* I - optional keywords for min    */
    maxname: &[[c_char; FLEN_VALUE]; 4], /* I - optional keywords for max    */
    binname: &[[c_char; FLEN_VALUE]; 4], /* I - optional keywords for binsize */
    weightin: f64,       /* I - binning weighting factor          */
    wtcol: &[c_char; FLEN_VALUE], /* I - optional keyword or col for weight*/
    recip: c_int,        /* I - use reciprocal of the weight?     */
    selectrow: Option<&[c_char]>, /* I - optional array (length = no. of   */
    /* rows in the table).  If the element is true */
    /* then the corresponding row of the table will*/
    /* be included in the histogram, otherwise the */
    /* row will be skipped.  Ingnored if *selectrow*/
    /* is equal to NULL.                           */
    status: &mut c_int,
) -> Option<Box<fitsfile>> {
    let mut histptr: Option<Box<fitsfile>> = None;
    let bitpix: c_int;
    let mut colnum: [c_int; 4] = [0; 4];
    let mut wtcolnum: c_int = 0;
    let mut haxes: [c_long; 4] = [0; 4];
    let mut amin: [f64; 4] = [0.0; 4];
    let mut amax: [f64; 4] = [0.0; 4];
    let mut binsize: [f64; 4] = [0.0; 4];
    let mut weight: f64 = 0.0;

    if *status > 0 {
        return None;
    }

    if naxis > 4 {
        ffpmsg_str("histogram has more than 4 dimensions");
        *status = BAD_DIMEN;
        return None;
    }

    /* reset position to the correct HDU if necessary */
    if (fptr).HDUposition != ((fptr).Fptr).curhdu {
        ffmahd_safe(fptr, ((fptr).HDUposition) + 1, None, status);
    }

    if imagetype == TBYTE {
        bitpix = BYTE_IMG;
    } else if imagetype == TSHORT {
        bitpix = SHORT_IMG;
    } else if imagetype == TINT {
        bitpix = LONG_IMG;
    } else if imagetype == TFLOAT {
        bitpix = FLOAT_IMG;
    } else if imagetype == TDOUBLE {
        bitpix = DOUBLE_IMG;
    } else {
        *status = BAD_DATATYPE;
        return None;
    }

    /*    Calculate the binning parameters:    */
    /*   columm numbers, axes length, min values,  max values, and binsizes.  */

    // Make mutable copies for fits_calc_binningd_safe
    let mut colname_mut = *colname;
    let mut minin_mut: [f64; 4] = [minin[0], minin[1], minin[2], minin[3]];
    let mut maxin_mut: [f64; 4] = [maxin[0], maxin[1], maxin[2], maxin[3]];
    let mut binsizein_mut: [f64; 4] = [binsizein[0], binsizein[1], binsizein[2], binsizein[3]];

    if fits_calc_binningd_safe(
        fptr,
        naxis,
        &mut colname_mut,
        &mut minin_mut,
        &mut maxin_mut,
        &mut binsizein_mut,
        minname,
        maxname,
        binname,
        &mut colnum,
        &mut haxes,
        &mut amin,
        &mut amax,
        &mut binsize,
        status,
    ) > 0
    {
        ffpmsg_str("failed to determine binning parameters");
        return None;
    }

    /* get the histogramming weighting factor, if any */
    if wtcol[0] != 0 {
        /* first, look for a keyword with the weight value */
        if fits_read_key(
            fptr,
            KeywordDatatypeMut::TDOUBLE(&mut weight),
            wtcol,
            None,
            status,
        ) != 0
        {
            /* not a keyword, so look for column with this name */
            *status = 0;

            /* get the column number in the table */
            if ffgcno_safe(
                fptr,
                CASEINSEN.try_into().unwrap(),
                wtcol,
                &mut wtcolnum,
                status,
            ) > 0
            {
                ffpmsg_str("keyword or column for histogram weights doesn't exist: ");
                ffpmsg_slice(wtcol);
                return None;
            }

            weight = DOUBLENULLVALUE;
        }
    } else {
        weight = weightin;
    }

    if weight <= 0. && weight != DOUBLENULLVALUE {
        ffpmsg_str("Illegal histogramming weighting factor <= 0.");
        *status = URL_PARSE_ERROR;
        return None;
    }

    if recip != 0 && weight != DOUBLENULLVALUE {
        /* take reciprocal of weight */
        weight = 1.0 / weight;
    }

    /* size of histogram is now known, so create temp output file */
    if fits_create_file(&mut histptr, outfile, status) > 0 {
        ffpmsg_str("failed to create temp output file for histogram");
        return None;
    }

    let mut histptr = histptr.unwrap();

    /* create output FITS image HDU */
    if ffcrim_safe(&mut histptr, bitpix, naxis, &haxes, status) > 0 {
        ffpmsg_str("failed to create output histogram FITS image");
        return None;
    }

    /* copy header keywords, converting pixel list WCS keywords to image WCS */
    if fits_copy_pixlist2image_safe(fptr, &mut histptr, 9, naxis, &colnum, status) > 0 {
        ffpmsg_str("failed to copy pixel list keywords to new histogram header");
        return None;
    }

    /* if the table columns have no WCS keywords, then write default keywords */
    fits_write_keys_histo_safe(fptr, &mut histptr, naxis, &colnum, status);

    /* update the WCS keywords for the ref. pixel location, and pixel size */
    fits_rebin_wcsd_safe(&mut histptr, naxis, &mut amin, &mut binsize, status);

    /* now compute the output image by binning the column values */
    if fits_make_histd_safe(
        fptr,
        &mut histptr,
        bitpix,
        naxis,
        &haxes,
        &colnum,
        &amin,
        &amax,
        &binsize,
        weight,
        wtcolnum,
        recip,
        selectrow,
        status,
    ) > 0
    {
        ffpmsg_str("failed to calculate new histogram values");
        return None;
    }

    Some(histptr)
}

/*--------------------------------------------------------------------------*/
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffhist(
    fptr: *mut Option<Box<fitsfile>>, /* I - pointer to table with X and Y cols; on output, points to histogram image    */
    outfile: *const c_char,           /* I - name for the output histogram file      */
    imagetype: c_int,                 /* I - datatype for image: TINT, TSHORT, etc   */
    naxis: c_int,                     /* I - number of axes in the histogram image   */
    colname: *mut [[c_char; FLEN_VALUE]; 4], /* I - column names               */
    minin: *mut f64,                  /* I - minimum histogram value, for each axis */
    maxin: *mut f64,                  /* I - maximum histogram value, for each axis */
    binsizein: *mut f64,              /* I - bin size along each axis               */
    minname: *const [[c_char; FLEN_VALUE]; 4], /* I - optional keywords for min    */
    maxname: *const [[c_char; FLEN_VALUE]; 4], /* I - optional keywords for max    */
    binname: *const [[c_char; FLEN_VALUE]; 4], /* I - optional keywords for binsize */
    weightin: f64,                    /* I - binning weighting factor          */
    wtcol: *const [c_char; FLEN_VALUE], /* I - optional keyword or col for weight*/
    recip: c_int,                     /* I - use reciprocal of the weight?     */
    selectrow: *const c_char,         /* I - optional array (length = no. of   */
    /* rows in the table).  If the element is true */
    /* then the corresponding row of the table will*/
    /* be included in the histogram, otherwise the */
    /* row will be skipped.  Ingnored if *selectrow*/
    /* is equal to NULL.                           */
    status: *mut c_int,
) -> c_int {
    unsafe {
        let fptr = (fptr).as_mut().expect(NULL_MSG);
        raw_to_slice!(outfile);
        let colname = colname.as_mut().expect(NULL_MSG);
        let minin = slice::from_raw_parts_mut(minin, naxis as usize);
        let maxin = slice::from_raw_parts_mut(maxin, naxis as usize);
        let binsizein = slice::from_raw_parts_mut(binsizein, naxis as usize);
        let minname = minname.as_ref().expect(NULL_MSG);
        let maxname = maxname.as_ref().expect(NULL_MSG);
        let binname = binname.as_ref().expect(NULL_MSG);
        let wtcol = wtcol.as_ref().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);

        let mut nrows = 0;

        ffgnrw_safe(fptr.as_deref_mut().unwrap(), &mut nrows, status); /* no. of rows */

        let selectrow = if selectrow.is_null() {
            None
        } else {
            Some(slice::from_raw_parts(selectrow, nrows as usize))
        };

        ffhist_safe(
            fptr, outfile, imagetype, naxis, colname, minin, maxin, binsizein, minname, maxname,
            binname, weightin, wtcol, recip, selectrow, status,
        )
    }
}

/*--------------------------------------------------------------------------*/
pub fn ffhist_safe(
    fptr: &mut Option<Box<fitsfile>>, /* I - pointer to table with X and Y cols; on output, points to histogram image    */
    outfile: &[c_char],               /* I - name for the output histogram file      */
    imagetype: c_int,                 /* I - datatype for image: TINT, TSHORT, etc   */
    naxis: c_int,                     /* I - number of axes in the histogram image   */
    colname: &mut [[c_char; FLEN_VALUE]; 4], /* I - column names               */
    minin: &mut [f64],                /* I - minimum histogram value, for each axis */
    maxin: &mut [f64],                /* I - maximum histogram value, for each axis */
    binsizein: &mut [f64],            /* I - bin size along each axis               */
    minname: &[[c_char; FLEN_VALUE]; 4], /* I - optional keywords for min    */
    maxname: &[[c_char; FLEN_VALUE]; 4], /* I - optional keywords for max    */
    binname: &[[c_char; FLEN_VALUE]; 4], /* I - optional keywords for binsize */
    weightin: f64,                    /* I - binning weighting factor          */
    wtcol: &[c_char; FLEN_VALUE],     /* I - optional keyword or col for weight*/
    recip: c_int,                     /* I - use reciprocal of the weight?     */
    selectrow: Option<&[c_char]>,     /* I - optional array (length = no. of   */
    /* rows in the table).  If the element is true */
    /* then the corresponding row of the table will*/
    /* be included in the histogram, otherwise the */
    /* row will be skipped.  Ingnored if *selectrow*/
    /* is equal to NULL.                           */
    status: &mut c_int,
) -> c_int {
    let mut ii: c_int;
    let mut datatype: c_int = 0;
    let mut repeat: c_int;
    let mut imin: c_int;
    let mut imax: c_int;
    let mut ibin: c_int;
    let bitpix: c_int;
    let mut tstatus: c_int;
    let mut use_datamax: c_int = 0;
    let mut haxes: [c_long; 4] = [0; 4];
    let mut histptr: Option<Box<fitsfile>> = None;
    let mut errmsg: [c_char; FLEN_ERRMSG] = [0; FLEN_ERRMSG];
    let mut keyname: [c_char; FLEN_KEYWORD] = [0; FLEN_KEYWORD];
    let mut card: [c_char; FLEN_CARD] = [0; FLEN_CARD];
    let mut comm: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];
    let mut colptr: &tcolumn;
    let mut imagepars: [iteratorCol; 1] = [iteratorCol::default(); 1];
    let n_cols: c_int = 1;
    let mut nkeys: c_int = 0;
    let offset: c_long = 0;
    let n_per_loop: c_long = -1; /* force whole array to be passed at one time */
    let mut histData: HistType = HistType::default(); /* Structure holding histogram info for iterator */

    let mut amin: [f64; 4] = [0.0; 4];
    let mut amax: [f64; 4] = [0.0; 4];
    let mut binsize: [f64; 4] = [0.0; 4];
    let mut maxbin: [f64; 4] = [0.0; 4];
    let mut datamin: f64 = DOUBLENULLVALUE;
    let mut datamax: f64 = DOUBLENULLVALUE;
    let mut svalue: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut dvalue: f64 = 0.0;
    let mut cpref: [[c_char; FLEN_VALUE]; 4] = [[0; FLEN_VALUE]; 4];
    let mut cptr: *mut c_char;

    if *status > 0 {
        return *status;
    }

    if naxis > 4 {
        ffpmsg_str("histogram has more than 4 dimensions");
        *status = BAD_DIMEN;
        return *status;
    }

    if let Some(fptr_inner) = fptr.as_deref_mut() {
        /* reset position to the correct HDU if necessary */
        if fptr_inner.HDUposition != fptr_inner.Fptr.curhdu {
            let hdu_pos = fptr_inner.HDUposition + 1;
            ffmahd_safe(fptr_inner, hdu_pos, None, status);
        }

        histData.tblptr = fptr_inner;
        histData.himagetype = imagetype;
        histData.haxis = naxis;
        histData.rowselector = selectrow
            .map(|s| s.as_ptr() as *mut c_char)
            .unwrap_or(ptr::null_mut());

        if imagetype == TBYTE {
            bitpix = BYTE_IMG;
        } else if imagetype == TSHORT {
            bitpix = SHORT_IMG;
        } else if imagetype == TINT {
            bitpix = LONG_IMG;
        } else if imagetype == TFLOAT {
            bitpix = FLOAT_IMG;
        } else if imagetype == TDOUBLE {
            bitpix = DOUBLE_IMG;
        } else {
            *status = BAD_DATATYPE;
            return *status;
        }

        /* The CPREF keyword, if it exists, gives the preferred columns. */
        /* Otherwise, assume "X", "Y", "Z", and "T"  */

        let (cpref0, other) = cpref.split_at_mut(1);
        let cpref0 = &mut cpref0[0];
        let (cpref1, other) = other.split_at_mut(1);
        let cpref1 = &mut cpref1[0];
        let (cpref2, other) = other.split_at_mut(1);
        let cpref2 = &mut cpref2[0];
        let cpref3 = &mut other[0];

        tstatus = 0;
        ffgky_safe(
            fptr_inner,
            KeywordDatatypeMut::TSTRING(cpref0),
            cs!(c"CPREF"),
            None,
            &mut tstatus,
        );

        let mut cptr: &mut [c_char];
        if tstatus == 0 {
            /* Preferred column names are given;  separate them */
            cptr = cpref0;

            /* the first preferred axis... */
            while cptr[0] != bb(b',') && cptr[0] != 0 {
                cptr = &mut cptr[1..];
            }

            if cptr[0] != 0 {
                cptr[0] = 0;
                cptr = &mut cptr[1..];
                while cptr[0] == bb(b' ') {
                    cptr = &mut cptr[1..];
                }

                strcpy_safe(cpref1, cptr);
                cptr = cpref1;

                /* the second preferred axis... */
                while cptr[0] != bb(b',') && cptr[0] != 0 {
                    cptr = &mut cptr[1..];
                }

                if cptr[0] != 0 {
                    cptr[0] = 0;
                    cptr = &mut cptr[1..];
                    while cptr[0] == bb(b' ') {
                        cptr = &mut cptr[1..];
                    }

                    strcpy_safe(cpref2, cptr);
                    cptr = cpref2;

                    /* the third preferred axis... */
                    while cptr[0] != bb(b',') && cptr[0] != 0 {
                        cptr = &mut cptr[1..];
                    }

                    if cptr[0] != 0 {
                        cptr[0] = 0;
                        cptr = &mut cptr[1..];
                        while cptr[0] == bb(b' ') {
                            cptr = &mut cptr[1..];
                        }

                        strcpy_safe(cpref3, cptr);
                    }
                }
            }
        }

        for ii in 0..naxis as usize {
            /* get the min, max, and binsize values from keywords, if specified */

            if minname[ii][0] != 0
                && ffgky_safe(
                    fptr_inner,
                    KeywordDatatypeMut::TDOUBLE(&mut minin[ii]),
                    &minname[ii],
                    None,
                    status,
                ) != 0
            {
                ffpmsg_str("error reading histogramming minimum keyword");
                ffpmsg_slice(&minname[ii]);
                return *status;
            }

            if maxname[ii][0] != 0
                && ffgky_safe(
                    fptr_inner,
                    KeywordDatatypeMut::TDOUBLE(&mut maxin[ii]),
                    &maxname[ii],
                    None,
                    status,
                ) != 0
            {
                ffpmsg_str("error reading histogramming maximum keyword");
                ffpmsg_slice(&maxname[ii]);
                return *status;
            }

            if binname[ii][0] != 0
                && ffgky_safe(
                    fptr_inner,
                    KeywordDatatypeMut::TDOUBLE(&mut binsizein[ii]),
                    &binname[ii],
                    None,
                    status,
                ) != 0
            {
                ffpmsg_str("error reading histogramming binsize keyword");
                ffpmsg_slice(&binname[ii]);
                return *status;
            }

            if binsizein[ii] == 0. {
                ffpmsg_str("error: histogram binsize = 0");
                *status = ZERO_SCALE;
                return *status;
            }

            if colname[ii][0] == 0 {
                strcpy_safe(&mut colname[ii], &cpref[ii]); /* try using the preferred column */
                if colname[ii][0] == 0 {
                    if ii == 0 {
                        strcpy_safe(&mut colname[ii], cs!(c"X"));
                    } else if ii == 1 {
                        strcpy_safe(&mut colname[ii], cs!(c"Y"));
                    } else if ii == 2 {
                        strcpy_safe(&mut colname[ii], cs!(c"Z"));
                    } else if ii == 3 {
                        strcpy_safe(&mut colname[ii], cs!(c"T"));
                    }
                }
            }

            /* get the column number in the table */
            if ffgcno_safe(
                fptr_inner,
                CASEINSEN as c_int,
                &colname[ii],
                &mut histData.hcolnum[ii],
                status,
            ) > 0
            {
                strcpy_safe(
                    &mut errmsg,
                    cs!(c"column for histogram axis doesn't exist: "),
                );
                let errmsg_len = strlen_safe(&errmsg);
                strncat_safe(&mut errmsg, &colname[ii], FLEN_ERRMSG - errmsg_len - 1);
                ffpmsg_slice(&errmsg);
                return *status;
            }

            colptr = &fptr_inner.Fptr.get_tableptr_as_slice()[(histData.hcolnum[ii] - 1) as usize];

            repeat = colptr.trepeat as c_int; /* vector repeat factor of the column */
            if repeat > 1 {
                strcpy_safe(&mut errmsg, cs!(c"Can't bin a vector column: "));
                let errmsg_len = strlen_safe(&errmsg);
                strncat_safe(&mut errmsg, &colname[ii], FLEN_ERRMSG - errmsg_len - 1);
                ffpmsg_slice(&errmsg);
                *status = BAD_DATATYPE;
                return *status;
            }

            /* get the datatype of the column */
            fits_get_eqcoltype(
                fptr_inner,
                histData.hcolnum[ii],
                Some(&mut datatype),
                None,
                None,
                status,
            );

            if datatype < 0 || datatype == TSTRING {
                strcpy_safe(
                    &mut errmsg,
                    cs!(c"Inappropriate datatype; can't bin this column: "),
                );
                let errmsg_len = strlen_safe(&errmsg);
                strncat_safe(&mut errmsg, &colname[ii], FLEN_ERRMSG - errmsg_len - 1);
                ffpmsg_slice(&errmsg);
                *status = BAD_DATATYPE;
                return *status;
            }

            /* use TLMINn and TLMAXn keyword values if min and max were not given */
            /* else use actual data min and max if TLMINn and TLMAXn don't exist */

            if minin[ii] == DOUBLENULLVALUE {
                ffkeyn_safe(cs!(c"TLMIN"), histData.hcolnum[ii], &mut keyname, status);
                if ffgky_safe(
                    fptr_inner,
                    KeywordDatatypeMut::TDOUBLE(&mut amin[ii]),
                    &keyname,
                    None,
                    status,
                ) > 0
                {
                    /* use actual data minimum value for the histogram minimum */
                    *status = 0;
                    if fits_get_col_minmax(
                        fptr_inner,
                        histData.hcolnum[ii],
                        &mut amin[ii],
                        &mut datamax,
                        status,
                    ) > 0
                    {
                        strcpy_safe(
                            &mut errmsg,
                            cs!(c"Error calculating datamin and datamax for column: "),
                        );
                        let errmsg_len = strlen_safe(&errmsg);
                        strncat_safe(&mut errmsg, &colname[ii], FLEN_ERRMSG - errmsg_len - 1);
                        ffpmsg_slice(&errmsg);
                        return *status;
                    }
                }
            } else {
                amin[ii] = minin[ii];
            }

            if maxin[ii] == DOUBLENULLVALUE {
                ffkeyn_safe(cs!(c"TLMAX"), histData.hcolnum[ii], &mut keyname, status);
                if ffgky_safe(
                    fptr_inner,
                    KeywordDatatypeMut::TDOUBLE(&mut amax[ii]),
                    &keyname,
                    None,
                    status,
                ) > 0
                {
                    *status = 0;
                    if datamax != DOUBLENULLVALUE
                    /* already computed max value */
                    {
                        amax[ii] = datamax;
                    } else {
                        /* use actual data maximum value for the histogram maximum */
                        if fits_get_col_minmax(
                            fptr_inner,
                            histData.hcolnum[ii],
                            &mut datamin,
                            &mut amax[ii],
                            status,
                        ) > 0
                        {
                            strcpy_safe(
                                &mut errmsg,
                                cs!(c"Error calculating datamin and datamax for column: "),
                            );
                            let errmsg_len = strlen_safe(&errmsg);
                            strncat_safe(&mut errmsg, &colname[ii], FLEN_ERRMSG - errmsg_len - 1);
                            ffpmsg_slice(&errmsg);
                            return *status;
                        }
                    }
                }
                use_datamax = 1; /* flag that the max was determined by the data values */
            /* and not specifically set by the calling program */
            } else {
                amax[ii] = maxin[ii];
            }

            /* use TDBINn keyword or else 1 if bin size is not given */
            if binsizein[ii] == DOUBLENULLVALUE {
                tstatus = 0;
                ffkeyn_safe(
                    cs!(c"TDBIN"),
                    histData.hcolnum[ii],
                    &mut keyname,
                    &mut tstatus,
                );

                if ffgky_safe(
                    fptr_inner,
                    KeywordDatatypeMut::TDOUBLE(&mut binsizein[ii]),
                    &keyname,
                    None,
                    &mut tstatus,
                ) > 0
                {
                    /* make at least 10 bins */
                    binsizein[ii] = (amax[ii] - amin[ii]) / 10.;
                    if binsizein[ii] > 1. {
                        binsizein[ii] = 1.; /* use default bin size */
                    }
                }
            }

            if (amin[ii] > amax[ii] && binsizein[ii] > 0.)
                || (amin[ii] < amax[ii] && binsizein[ii] < 0.)
            {
                binsize[ii] = -binsizein[ii]; /* reverse the sign of binsize */
            } else {
                binsize[ii] = binsizein[ii]; /* binsize has the correct sign */
            }

            ibin = binsize[ii] as c_int;
            imin = amin[ii] as c_int;
            imax = amax[ii] as c_int;

            /* Determine the range and number of bins in the histogram. This  */
            /* depends on whether the input columns are integer or floats, so */
            /* treat each case separately.                                    */

            if datatype <= TLONG
                && f64::from(imin) == amin[ii]
                && f64::from(imax) == amax[ii]
                && f64::from(ibin) == binsize[ii]
            {
                /* This is an integer column and integer limits were entered. */
                /* Shift the lower and upper histogramming limits by 0.5, so that */
                /* the values fall in the center of the bin, not on the edge. */

                haxes[ii] = c_long::from((imax - imin) / ibin + 1); /* last bin may only */
                /* be partially full */
                maxbin[ii] = haxes[ii] as f64 + 1.0; /* add 1. instead of 0.5 to avoid roundoff */

                if amin[ii] < amax[ii] {
                    amin[ii] -= 0.5;
                    amax[ii] += 0.5;
                } else {
                    amin[ii] += 0.5;
                    amax[ii] -= 0.5;
                }
            } else if use_datamax != 0 {
                /* Either the column datatype and/or the limits are floating point, */
                /* and the histogram limits are being defined by the min and max */
                /* values of the array.  Add 1 to the number of histogram bins to */
                /* make sure that pixels that are equal to the maximum or are */
                /* in the last partial bin are included.  */

                maxbin[ii] = (amax[ii] - amin[ii]) / binsize[ii];
                haxes[ii] = (maxbin[ii] as c_long + 1) as c_long;
            } else {
                /*  float datatype column and/or limits, and the maximum value to */
                /*  include in the histogram is specified by the calling program. */
                /*  The lower limit is inclusive, but upper limit is exclusive    */
                maxbin[ii] = (amax[ii] - amin[ii]) / binsize[ii];
                haxes[ii] = maxbin[ii] as c_long;

                if amin[ii] < amax[ii] {
                    if amin[ii] + (haxes[ii] as f64 * binsize[ii]) < amax[ii] {
                        haxes[ii] += 1; /* need to include another partial bin */
                    }
                } else if amin[ii] + (haxes[ii] as f64 * binsize[ii]) > amax[ii] {
                    haxes[ii] += 1; /* need to include another partial bin */
                }
            }
        }

        /* get the histogramming weighting factor */
        if wtcol[0] != 0 {
            /* first, look for a keyword with the weight value */
            if ffgky_safe(
                fptr_inner,
                KeywordDatatypeMut::TDOUBLE(&mut histData.weight),
                wtcol,
                None,
                status,
            ) != 0
            {
                /* not a keyword, so look for column with this name */
                *status = 0;

                /* get the column number in the table */
                if ffgcno_safe(
                    fptr_inner,
                    CASEINSEN as c_int,
                    wtcol,
                    &mut histData.wtcolnum,
                    status,
                ) > 0
                {
                    ffpmsg_str("keyword or column for histogram weights doesn't exist: ");
                    ffpmsg_slice(wtcol);
                    return *status;
                }

                histData.weight = DOUBLENULLVALUE;
            }
        } else {
            histData.weight = weightin;
        }

        if histData.weight <= 0. && histData.weight != DOUBLENULLVALUE {
            ffpmsg_str("Illegal histogramming weighting factor <= 0.");
            *status = URL_PARSE_ERROR;
            return *status;
        }

        if recip != 0 && histData.weight != DOUBLENULLVALUE {
            /* take reciprocal of weight */
            histData.weight = 1.0 / histData.weight;
        }

        histData.wtrecip = recip;

        /* size of histogram is now known, so create temp output file */
        if ffinit_safe(&mut histptr, outfile, status) > 0 {
            ffpmsg_str("failed to create temp output file for histogram");
            return *status;
        }

        let histptr_inner = histptr.as_deref_mut().unwrap();

        if ffcrim_safe(histptr_inner, bitpix, histData.haxis, &haxes, status) > 0 {
            ffpmsg_str("failed to create primary array histogram in temp file");
            ffclos_safe(histptr.unwrap(), status);
            return *status;
        }

        /* copy all non-structural keywords from the table to the image */
        fits_get_hdrspace(fptr_inner, Some(&mut nkeys), None, status);
        for ii in 1..=nkeys {
            fits_read_record(fptr_inner, ii, Some(&mut card), status);
            if fits_get_keyclass(&card) >= 120 {
                fits_write_record(histptr_inner, &card, status);
            }
        }

        /* Set global variables with histogram parameter values.    */
        /* Use separate scalar variables rather than arrays because */
        /* it is more efficient when computing the histogram.       */

        histData.amin1 = amin[0];
        histData.maxbin1 = maxbin[0];
        histData.binsize1 = binsize[0];
        histData.haxis1 = haxes[0];

        if histData.haxis > 1 {
            histData.amin2 = amin[1];
            histData.maxbin2 = maxbin[1];
            histData.binsize2 = binsize[1];
            histData.haxis2 = haxes[1];

            if histData.haxis > 2 {
                histData.amin3 = amin[2];
                histData.maxbin3 = maxbin[2];
                histData.binsize3 = binsize[2];
                histData.haxis3 = haxes[2];

                if histData.haxis > 3 {
                    histData.amin4 = amin[3];
                    histData.maxbin4 = maxbin[3];
                    histData.binsize4 = binsize[3];
                    histData.haxis4 = haxes[3];
                }
            }
        }

        /* define parameters of image for the iterator function */
        fits_iter_set_file_safe(&mut imagepars[0], histptr_inner); /* pointer to image */
        fits_iter_set_datatype_safe(&mut imagepars[0], imagetype); /* image datatype   */
        fits_iter_set_iotype_safe(&mut imagepars[0], OUTPUT_COL); /* image is output  */

        /* call the iterator function to write out the histogram image */
        if fits_iterate_data(
            n_cols,
            &mut imagepars,
            offset,
            n_per_loop,
            ffwritehisto,
            ((&mut histData) as *mut HistType).cast::<c_void>(),
            status,
        ) != 0
        {
            return *status;
        }

        /* write the World Coordinate System (WCS) keywords */
        /* create default values if WCS keywords are not present in the table */
        for ii in 0..histData.haxis as usize {
            /*  CTYPEn  */
            tstatus = 0;
            ffkeyn_safe(
                cs!(c"TCTYP"),
                histData.hcolnum[ii],
                &mut keyname,
                &mut tstatus,
            );
            ffgky_safe(
                fptr_inner,
                KeywordDatatypeMut::TSTRING(&mut svalue),
                &keyname,
                None,
                &mut tstatus,
            );
            if tstatus != 0 {
                /* just use column name as the type */
                tstatus = 0;
                ffkeyn_safe(
                    cs!(c"TTYPE"),
                    histData.hcolnum[ii],
                    &mut keyname,
                    &mut tstatus,
                );
                ffgky_safe(
                    fptr_inner,
                    KeywordDatatypeMut::TSTRING(&mut svalue),
                    &keyname,
                    None,
                    &mut tstatus,
                );
            }

            if tstatus == 0 {
                ffkeyn_safe(cs!(c"CTYPE"), (ii + 1) as c_int, &mut keyname, &mut tstatus);
                ffpky_safe(
                    histptr_inner,
                    KeywordDatatype::TSTRING(&svalue),
                    &keyname,
                    Some(cs!(c"Coordinate Type")),
                    &mut tstatus,
                );
            } else {
                tstatus = 0;
            }

            /*  CUNITn  */
            ffkeyn_safe(
                cs!(c"TCUNI"),
                histData.hcolnum[ii],
                &mut keyname,
                &mut tstatus,
            );
            ffgky_safe(
                fptr_inner,
                KeywordDatatypeMut::TSTRING(&mut svalue),
                &keyname,
                None,
                &mut tstatus,
            );
            if tstatus != 0 {
                /* use the column units */
                tstatus = 0;
                ffkeyn_safe(
                    cs!(c"TUNIT"),
                    histData.hcolnum[ii],
                    &mut keyname,
                    &mut tstatus,
                );
                ffgky_safe(
                    fptr_inner,
                    KeywordDatatypeMut::TSTRING(&mut svalue),
                    &keyname,
                    None,
                    &mut tstatus,
                );
            }

            if tstatus == 0 {
                ffkeyn_safe(cs!(c"CUNIT"), (ii + 1) as c_int, &mut keyname, &mut tstatus);
                ffpky_safe(
                    histptr_inner,
                    KeywordDatatype::TSTRING(&svalue),
                    &keyname,
                    Some(cs!(c"Coordinate Units")),
                    &mut tstatus,
                );
            } else {
                tstatus = 0;
            }

            /*  CRPIXn  - Reference Pixel  */
            ffkeyn_safe(
                cs!(c"TCRPX"),
                histData.hcolnum[ii],
                &mut keyname,
                &mut tstatus,
            );
            ffgky_safe(
                fptr_inner,
                KeywordDatatypeMut::TDOUBLE(&mut dvalue),
                &keyname,
                None,
                &mut tstatus,
            );
            if tstatus != 0 {
                dvalue = 1.0; /* choose first pixel in new image as ref. pix. */
                tstatus = 0;
            } else {
                /* calculate locate of the ref. pix. in the new image */
                dvalue = (dvalue - amin[ii]) / binsize[ii] + 0.5;
            }

            ffkeyn_safe(cs!(c"CRPIX"), (ii + 1) as c_int, &mut keyname, &mut tstatus);
            ffpky_safe(
                histptr_inner,
                KeywordDatatype::TDOUBLE(&dvalue),
                &keyname,
                Some(cs!(c"Reference Pixel")),
                &mut tstatus,
            );

            /*  CRVALn - Value at the location of the reference pixel */
            ffkeyn_safe(
                cs!(c"TCRVL"),
                histData.hcolnum[ii],
                &mut keyname,
                &mut tstatus,
            );
            ffgky_safe(
                fptr_inner,
                KeywordDatatypeMut::TDOUBLE(&mut dvalue),
                &keyname,
                None,
                &mut tstatus,
            );
            if tstatus != 0 {
                /* calculate value at ref. pix. location (at center of 1st pixel) */
                dvalue = amin[ii] + binsize[ii] / 2.;
                tstatus = 0;
            }

            ffkeyn_safe(cs!(c"CRVAL"), (ii + 1) as c_int, &mut keyname, &mut tstatus);
            ffpky_safe(
                histptr_inner,
                KeywordDatatype::TDOUBLE(&dvalue),
                &keyname,
                Some(cs!(c"Reference Value")),
                &mut tstatus,
            );

            /*  CDELTn - unit size of pixels  */
            ffkeyn_safe(
                cs!(c"TCDLT"),
                histData.hcolnum[ii],
                &mut keyname,
                &mut tstatus,
            );
            ffgky_safe(
                fptr_inner,
                KeywordDatatypeMut::TDOUBLE(&mut dvalue),
                &keyname,
                None,
                &mut tstatus,
            );
            if tstatus != 0 {
                dvalue = 1.0; /* use default pixel size */
                tstatus = 0;
            }

            dvalue *= binsize[ii];
            ffkeyn_safe(cs!(c"CDELT"), (ii + 1) as c_int, &mut keyname, &mut tstatus);
            ffpky_safe(
                histptr_inner,
                KeywordDatatype::TDOUBLE(&dvalue),
                &keyname,
                Some(cs!(c"Pixel size")),
                &mut tstatus,
            );

            /*  CROTAn - Rotation angle (degrees CCW)  */
            /*  There should only be a CROTA2 keyword, and only for 2+ D images */
            if ii == 1 {
                ffkeyn_safe(
                    cs!(c"TCROT"),
                    histData.hcolnum[ii],
                    &mut keyname,
                    &mut tstatus,
                );
                ffgky_safe(
                    fptr_inner,
                    KeywordDatatypeMut::TDOUBLE(&mut dvalue),
                    &keyname,
                    None,
                    &mut tstatus,
                );
                if tstatus == 0 && dvalue != 0.
                /* only write keyword if angle != 0 */
                {
                    ffkeyn_safe(cs!(c"CROTA"), (ii + 1) as c_int, &mut keyname, &mut tstatus);
                    ffpky_safe(
                        histptr_inner,
                        KeywordDatatype::TDOUBLE(&dvalue),
                        &keyname,
                        Some(cs!(c"Rotation angle")),
                        &mut tstatus,
                    );
                } else {
                    /* didn't find CROTA for the 2nd axis, so look for one */
                    /* on the first axis */
                    tstatus = 0;
                    ffkeyn_safe(
                        cs!(c"TCROT"),
                        histData.hcolnum[0],
                        &mut keyname,
                        &mut tstatus,
                    );
                    ffgky_safe(
                        fptr_inner,
                        KeywordDatatypeMut::TDOUBLE(&mut dvalue),
                        &keyname,
                        None,
                        &mut tstatus,
                    );
                    if tstatus == 0 && dvalue != 0.
                    /* only write keyword if angle != 0 */
                    {
                        dvalue *= -1.; /* negate the value, because mirror image */
                        ffkeyn_safe(cs!(c"CROTA"), (ii + 1) as c_int, &mut keyname, &mut tstatus);
                        ffpky_safe(
                            histptr_inner,
                            KeywordDatatype::TDOUBLE(&dvalue),
                            &keyname,
                            Some(cs!(c"Rotation angle")),
                            &mut tstatus,
                        );
                    }
                }
            }
        }

        /* convert any TPn_k keywords to PCi_j; the value remains unchanged */
        /* also convert any TCn_k to CDi_j; the value is modified by n binning size */
        /* This is a bit of a kludge, and only works for 2D WCS */

        if histData.haxis == 2 {
            /* PC1_1 */
            tstatus = 0;
            ffkeyn_safe(cs!(c"TP"), histData.hcolnum[0], &mut card, &mut tstatus);
            strcat_safe(&mut card, cs!(c"_"));
            ffkeyn_safe(&card, histData.hcolnum[0], &mut keyname, &mut tstatus);
            ffgky_safe(
                fptr_inner,
                KeywordDatatypeMut::TDOUBLE(&mut dvalue),
                &keyname,
                Some(&mut comm),
                &mut tstatus,
            );
            if tstatus == 0 {
                ffpky_safe(
                    histptr_inner,
                    KeywordDatatype::TDOUBLE(&dvalue),
                    cs!(c"PC1_1"),
                    Some(&comm),
                    &mut tstatus,
                );
            }

            tstatus = 0;
            keyname[1] = bb(b'C');
            ffgky_safe(
                fptr_inner,
                KeywordDatatypeMut::TDOUBLE(&mut dvalue),
                &keyname,
                Some(&mut comm),
                &mut tstatus,
            );
            if tstatus == 0 {
                dvalue *= binsize[0];
                ffpky_safe(
                    histptr_inner,
                    KeywordDatatype::TDOUBLE(&dvalue),
                    cs!(c"CD1_1"),
                    Some(&comm),
                    &mut tstatus,
                );
            }

            /* PC1_2 */
            tstatus = 0;
            ffkeyn_safe(cs!(c"TP"), histData.hcolnum[0], &mut card, &mut tstatus);
            strcat_safe(&mut card, cs!(c"_"));
            ffkeyn_safe(&card, histData.hcolnum[1], &mut keyname, &mut tstatus);
            ffgky_safe(
                fptr_inner,
                KeywordDatatypeMut::TDOUBLE(&mut dvalue),
                &keyname,
                Some(&mut comm),
                &mut tstatus,
            );
            if tstatus == 0 {
                ffpky_safe(
                    histptr_inner,
                    KeywordDatatype::TDOUBLE(&dvalue),
                    cs!(c"PC1_2"),
                    Some(&comm),
                    &mut tstatus,
                );
            }

            tstatus = 0;
            keyname[1] = bb(b'C');
            ffgky_safe(
                fptr_inner,
                KeywordDatatypeMut::TDOUBLE(&mut dvalue),
                &keyname,
                Some(&mut comm),
                &mut tstatus,
            );
            if tstatus == 0 {
                dvalue *= binsize[0];
                ffpky_safe(
                    histptr_inner,
                    KeywordDatatype::TDOUBLE(&dvalue),
                    cs!(c"CD1_2"),
                    Some(&comm),
                    &mut tstatus,
                );
            }

            /* PC2_1 */
            tstatus = 0;
            ffkeyn_safe(cs!(c"TP"), histData.hcolnum[1], &mut card, &mut tstatus);
            strcat_safe(&mut card, cs!(c"_"));
            ffkeyn_safe(&card, histData.hcolnum[0], &mut keyname, &mut tstatus);
            ffgky_safe(
                fptr_inner,
                KeywordDatatypeMut::TDOUBLE(&mut dvalue),
                &keyname,
                Some(&mut comm),
                &mut tstatus,
            );
            if tstatus == 0 {
                ffpky_safe(
                    histptr_inner,
                    KeywordDatatype::TDOUBLE(&dvalue),
                    cs!(c"PC2_1"),
                    Some(&comm),
                    &mut tstatus,
                );
            }

            tstatus = 0;
            keyname[1] = bb(b'C');
            ffgky_safe(
                fptr_inner,
                KeywordDatatypeMut::TDOUBLE(&mut dvalue),
                &keyname,
                Some(&mut comm),
                &mut tstatus,
            );
            if tstatus == 0 {
                dvalue *= binsize[1];
                ffpky_safe(
                    histptr_inner,
                    KeywordDatatype::TDOUBLE(&dvalue),
                    cs!(c"CD2_1"),
                    Some(&comm),
                    &mut tstatus,
                );
            }

            /* PC2_2 */
            tstatus = 0;
            ffkeyn_safe(cs!(c"TP"), histData.hcolnum[1], &mut card, &mut tstatus);
            strcat_safe(&mut card, cs!(c"_"));
            ffkeyn_safe(&card, histData.hcolnum[1], &mut keyname, &mut tstatus);
            ffgky_safe(
                fptr_inner,
                KeywordDatatypeMut::TDOUBLE(&mut dvalue),
                &keyname,
                Some(&mut comm),
                &mut tstatus,
            );
            if tstatus == 0 {
                ffpky_safe(
                    histptr_inner,
                    KeywordDatatype::TDOUBLE(&dvalue),
                    cs!(c"PC2_2"),
                    Some(&comm),
                    &mut tstatus,
                );
            }

            tstatus = 0;
            keyname[1] = bb(b'C');
            ffgky_safe(
                fptr_inner,
                KeywordDatatypeMut::TDOUBLE(&mut dvalue),
                &keyname,
                Some(&mut comm),
                &mut tstatus,
            );
            if tstatus == 0 {
                dvalue *= binsize[1];
                ffpky_safe(
                    histptr_inner,
                    KeywordDatatype::TDOUBLE(&dvalue),
                    cs!(c"CD2_2"),
                    Some(&comm),
                    &mut tstatus,
                );
            }
        }
    }
    /* finally, close the original file and return ptr to the new image */
    let tmp = fptr.take().unwrap();
    ffclos_safe(tmp, status);
    *fptr = histptr;

    *status
}

/*--------------------------------------------------------------------------*/
/// Single-precision version
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn fits_calc_binning(
    fptr: *mut fitsfile, /* IO - pointer to table to be binned      ;       */
    naxis: c_int,        /* I - number of axes/columns in the binned image  */
    colname: *mut [[c_char; FLEN_VALUE]; 4], /* I - optional column names         */
    minin: *mut f64,     /* I - optional lower bound value for each axis  */
    maxin: *mut f64,     /* I - optional upper bound value, for each axis */
    binsizein: *mut f64, /* I - optional bin size along each axis         */
    minname: *mut [[c_char; FLEN_VALUE]; 4], /* I - optional keywords for min       */
    maxname: *mut [[c_char; FLEN_VALUE]; 4], /* I - optional keywords for max       */
    binname: *mut [[c_char; FLEN_VALUE]; 4], /* I - optional keywords for binsize   */

    /* The returned parameters for each axis of the n-dimensional histogram are */
    colnum: *mut c_int, /* O - column numbers, to be binned */
    haxes: *mut c_long, /* O - number of bins in each histogram axis */
    amin: *mut f32,     /* O - lower bound of the histogram axes */
    amax: *mut f32,     /* O - upper bound of the histogram axes */
    binsize: *mut f32,  /* O - width of histogram bins/pixels on each axis */
    status: *mut c_int,
) -> c_int {
    unsafe {
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let colname = colname.as_mut().expect(NULL_MSG);
        let minin = slice::from_raw_parts_mut(minin.as_mut().expect(NULL_MSG), naxis as usize);
        let maxin = slice::from_raw_parts_mut(maxin.as_mut().expect(NULL_MSG), naxis as usize);
        let binsizein =
            slice::from_raw_parts_mut(binsizein.as_mut().expect(NULL_MSG), naxis as usize);
        let minname = minname.as_mut().expect(NULL_MSG);
        let maxname = maxname.as_mut().expect(NULL_MSG);
        let binname = binname.as_mut().expect(NULL_MSG);
        let colnum = slice::from_raw_parts_mut(colnum, naxis as usize);
        let haxes = slice::from_raw_parts_mut(haxes, naxis as usize);
        let amin = slice::from_raw_parts_mut(amin, naxis as usize);
        let amax = slice::from_raw_parts_mut(amax, naxis as usize);
        let binsize = slice::from_raw_parts_mut(binsize, naxis as usize);
        let status = status.as_mut().expect(NULL_MSG);

        fits_calc_binning_safe(
            fptr, naxis, colname, minin, maxin, binsizein, minname, maxname, binname, colnum,
            haxes, amin, amax, binsize, status,
        )
    }
}

/*--------------------------------------------------------------------------*/
/// Single-precision version
pub fn fits_calc_binning_safe(
    fptr: &mut fitsfile, /* IO - pointer to table to be binned      ;       */
    naxis: c_int,        /* I - number of axes/columns in the binned image  */
    colname: &mut [[c_char; FLEN_VALUE]; 4], /* I - optional column names         */
    minin: &mut [f64],   /* I - optional lower bound value for each axis  */
    maxin: &mut [f64],   /* I - optional upper bound value, for each axis */
    binsizein: &mut [f64], /* I - optional bin size along each axis         */
    minname: &mut [[c_char; FLEN_VALUE]; 4], /* I - optional keywords for min       */
    maxname: &mut [[c_char; FLEN_VALUE]; 4], /* I - optional keywords for max       */
    binname: &mut [[c_char; FLEN_VALUE]; 4], /* I - optional keywords for binsize   */

    /* The returned parameters for each axis of the n-dimensional histogram are */
    colnum: &mut [c_int], /* O - column numbers, to be binned */
    haxes: &mut [c_long], /* O - number of bins in each histogram axis */
    amin: &mut [f32],     /* O - lower bound of the histogram axes */
    amax: &mut [f32],     /* O - upper bound of the histogram axes */
    binsize: &mut [f32],  /* O - width of histogram bins/pixels on each axis */
    status: &mut c_int,
) -> c_int {
    let mut amind: [f64; 4] = [0.0; 4];
    let mut amaxd: [f64; 4] = [0.0; 4];
    let mut binsized: [f64; 4] = [0.0; 4];

    fits_calc_binningd_safe(
        fptr,
        naxis,
        colname,
        minin,
        maxin,
        binsizein,
        minname,
        maxname,
        binname,
        colnum,
        haxes,
        &mut amind,
        &mut amaxd,
        &mut binsized,
        status,
    );

    /* Copy double precision values into single precision */
    if *status == 0 {
        let mut naxis1 = 4;
        if naxis < naxis1 {
            naxis1 = naxis;
        }
        for i in 0..naxis1 as usize {
            amin[i] = amind[i] as f32;
            amax[i] = amaxd[i] as f32;
            binsize[i] = binsized[i] as f32;
        }
    }

    *status
}

/*--------------------------------------------------------------------------*/
/// Double-precision version with extended syntax
/// Calculate the actual binning parameters, based on various user inputoptions.
///
/// Note: caller is responsible to free parsers[*] upon return using ffcprs()
pub(crate) fn fits_calc_binningde(
    fptr: &mut fitsfile, /* IO - pointer to table to be binned      ;       */
    naxis: c_int,        /* I - number of axes/columns in the binned image  */
    colname: &mut [[c_char; FLEN_VALUE]; 4], /* I - optional column names         */
    colexpr: Option<&[Option<&[c_char]>; 4]>, /* I - optional column expressions instead of name    */
    minin: &mut [f64],                        /* IO - optional lower bound value for each axis  */
    maxin: &mut [f64],                        /* IO - optional upper bound value, for each axis */
    binsizein: &mut [f64],                    /* IO - optional bin size along each axis         */
    minname: &[[c_char; FLEN_VALUE]; 4],      /* I - optional keywords for min       */
    maxname: &[[c_char; FLEN_VALUE]; 4],      /* I - optional keywords for max       */
    binname: &[[c_char; FLEN_VALUE]; 4],      /* I - optional keywords for binsize   */

    /* The returned parameters for each axis of the n-dimensional histogram are */
    colnum: &mut [c_int],                /* O - column numbers, to be binned */
    mut datatypes: Option<&mut [c_int]>, /* O - datatype for each column */
    haxes: &mut [c_long],                /* O - number of bins in each histogram axis */
    amin: &mut [f64],                    /* O - lower bound of the histogram axes */
    amax: &mut [f64],                    /* O - upper bound of the histogram axes */
    binsize: &mut [f64],                 /* O - width of histogram bins/pixels on each axis */
    mut repeat: Option<&mut c_long>,     /* O - vector repeat of input columns */
    status: &mut c_int,
) -> c_int {
    let mut colptr: &mut tcolumn;
    let mut cpref: [[c_char; FLEN_VALUE]; 4] = [[0; FLEN_VALUE]; 4];
    let mut errmsg: [c_char; FLEN_ERRMSG] = [0; FLEN_ERRMSG];
    let mut keyname: [c_char; FLEN_KEYWORD] = [0; FLEN_KEYWORD];
    let mut tstatus: c_int;
    let mut ii: c_int;
    let mut datatype: c_int = 0;
    let mut imin: c_int;
    let mut imax: c_int;
    let mut ibin: c_int;
    let mut use_datamax: c_int = 0;
    let mut repeat1: c_long = 0;
    let mut datamin: f64;
    let mut datamax: f64;
    let mut ncols: c_int;

    /* check inputs */

    if *status > 0 {
        return *status;
    }

    /* Initialize the number of iterator columns required */
    if let Some(repeat) = repeat.as_deref_mut() {
        (*repeat) = 0;
    }

    if naxis > 4 {
        ffpmsg_str("histograms with more than 4 dimensions are not supported");
        *status = BAD_DIMEN;
        return *status;
    }

    /* reset position to the correct HDU if necessary */
    if (fptr).HDUposition != ((fptr).Fptr).curhdu {
        ffmahd_safe(fptr, ((fptr).HDUposition) + 1, None, status);
    }

    /* ============================================================= */
    /* The CPREF keyword, if it exists, gives the preferred columns. */
    /* Otherwise, assume "X", "Y", "Z", and "T"  */

    cpref[0][0] = 0;
    cpref[1][0] = 0;
    cpref[2][0] = 0;
    cpref[3][0] = 0;

    tstatus = 0;
    ffgky_safe(
        fptr,
        KeywordDatatypeMut::TSTRING(&mut cpref[0]),
        cs!(c"CPREF"),
        None,
        &mut tstatus,
    );

    // Work with indices instead of references to avoid borrow checker issues
    if tstatus == 0 {
        /* Preferred column names are given;  separate them */
        let mut idx = 0;

        /* the first preferred axis... */
        while cpref[0][idx] != bb(b',') && cpref[0][idx] != 0 {
            idx += 1;
        }

        if cpref[0][idx] != 0 {
            cpref[0][idx] = 0;
            idx += 1;
            while cpref[0][idx] == bb(b' ') {
                idx += 1;
            }

            // Copy from cpref[0][idx..] to cpref[1]
            let mut j = 0;
            while cpref[0][idx + j] != 0 && j < FLEN_VALUE - 1 {
                cpref[1][j] = cpref[0][idx + j];
                j += 1;
            }
            cpref[1][j] = 0;
            idx = 0;

            /* the second preferred axis... */
            while cpref[1][idx] != bb(b',') && cpref[1][idx] != 0 {
                idx += 1;
            }

            if cpref[1][idx] != 0 {
                cpref[1][idx] = 0;
                idx += 1;
                while cpref[1][idx] == bb(b' ') {
                    idx += 1;
                }

                // Copy from cpref[1][idx..] to cpref[2]
                let mut j = 0;
                while cpref[1][idx + j] != 0 && j < FLEN_VALUE - 1 {
                    cpref[2][j] = cpref[1][idx + j];
                    j += 1;
                }
                cpref[2][j] = 0;
                idx = 0;

                /* the third preferred axis... */
                while cpref[2][idx] != bb(b',') && cpref[2][idx] != 0 {
                    idx += 1;
                }

                if cpref[2][idx] != 0 {
                    cpref[2][idx] = 0;
                    idx += 1;
                    while cpref[2][idx] == bb(b' ') {
                        idx += 1;
                    }

                    // Copy from cpref[2][idx..] to cpref[3]
                    let mut j = 0;
                    while cpref[2][idx + j] != 0 && j < FLEN_VALUE - 1 {
                        cpref[3][j] = cpref[2][idx + j];
                        j += 1;
                    }
                    cpref[3][j] = 0;
                }
            }
        }
    }

    /* ============================================================= */
    /* Main Loop for calculating parameters for each column          */

    for ii in 0..naxis as usize {
        /* =========================================================== */
        /* Determine column Number, based on, in order of priority,
            1  input column name, or
        2  name given by CPREF keyword, or
        3  assume X, Y, Z and T for the name
         */

        let cond = colexpr.is_none()
            || colexpr.as_ref().unwrap()[ii].is_none()
            || colexpr.as_ref().unwrap()[ii].as_ref().unwrap()[0] == 0;
        if colname[ii][0] == 0 && cond {
            strcpy_safe(&mut colname[ii], &cpref[ii]); /* try using the preferred column */
            if colname[ii][0] == 0 {
                if ii == 0 {
                    strcpy_safe(&mut colname[ii], cs!(c"X"));
                } else if ii == 1 {
                    strcpy_safe(&mut colname[ii], cs!(c"Y"));
                } else if ii == 2 {
                    strcpy_safe(&mut colname[ii], cs!(c"Z"));
                } else if ii == 3 {
                    strcpy_safe(&mut colname[ii], cs!(c"T"));
                }
            }
        }

        /* get the column number in the table */
        colnum[ii] = 0;
        if cond {
            if ffgcno_safe(
                fptr,
                CASEINSEN as c_int,
                &colname[ii],
                &mut colnum[ii],
                status,
            ) > 0
            {
                strcpy_safe(
                    &mut errmsg,
                    cs!(c"column for histogram axis doesn't exist: "),
                );
                let errmsg_len = strlen_safe(&errmsg);
                strncat_safe(&mut errmsg, &colname[ii], FLEN_ERRMSG - errmsg_len - 1);
                ffpmsg_slice(&errmsg);
                return *status;
            }

            /* ================================================================ */
            /* check tha column is not a vector or a string                     */

            /* get the datatype of the column */
            fits_get_eqcoltype(
                fptr,
                colnum[ii],
                Some(&mut datatype),
                Some(&mut repeat1),
                None,
                status,
            );

            ncols = 1; /* Require only one iterator column, the actual column */
        } else {
            /* column expression: use parse to determine datatype and dimensions */

            let mut nelem: c_long = 0;
            let mut naxes: [c_long; MAXDIMS as usize] = [0; MAXDIMS as usize];
            let mut naxis: c_int = 0;
            let mut lParse: ParseData = ParseData::default();

            /* Initialize the parser so that we can determine the datatype
            of the returned type as well as the vector dimensions */
            let ce = colexpr.unwrap();
            if ffiprs(
                fptr,
                0,
                ce[ii].unwrap(),
                MAXDIMS,
                &mut datatype,
                &mut nelem,
                &mut naxis,
                &mut naxes,
                &mut lParse,
                status,
            ) != 0
            {
                strcpy_safe(&mut errmsg, cs!(c"Parser error of binning expression: "));
                let errmsg_len = strlen_safe(&errmsg);
                strncat_safe(&mut errmsg, ce[ii].unwrap(), FLEN_ERRMSG - errmsg_len - 1);
                ffpmsg_slice(&errmsg);
                ffcprs(&mut lParse);
                return *status;
            }
            if nelem < 0 {
                nelem = 1; /* If it's a constant expression */
            }

            repeat1 = nelem;

            /* We require lParse.nCols columns to be read from input,
            plus one for the Temporary calculator result */
            ncols = lParse.nCols + 1;
            ffcprs(&mut lParse);
        }

        /* Not sure why this repeat limitation is here -- CM
             The iterator system can handle vector columns just fine
        `       */
        if datatype < 0 || datatype == TSTRING {
            strcpy_safe(
                &mut errmsg,
                cs!(c"Inappropriate datatype; can't bin this column: "),
            );
            let errmsg_len = strlen_safe(&errmsg);
            strncat_safe(&mut errmsg, &colname[ii], FLEN_ERRMSG - errmsg_len - 1);
            ffpmsg_slice(&errmsg);
            *status = BAD_DATATYPE;
            return *status;
        }

        /* Store repeat value for future use */
        if let Some(repeat) = repeat.as_deref_mut() {
            if ii == 0 {
                *repeat = repeat1; /* First time around save the repeat value */
            } else if *repeat != repeat1 {
                /* later dimensions, keep same dims */

                strcpy_safe(
                    &mut errmsg,
                    cs!(c"Vector repeat of input columns do not agree"),
                );
                ffpmsg_slice(&errmsg);
                *status = BAD_DIMEN;
                return *status;
            }
        }

        if let Some(datatypes) = datatypes.as_deref_mut() {
            datatypes[ii] = datatype;
        }

        /* ================================================================ */
        /* get the minimum value */

        datamin = DOUBLENULLVALUE;
        datamax = DOUBLENULLVALUE;

        if minname[ii][0] != 0
            && ffgky_safe(
                fptr,
                KeywordDatatypeMut::TDOUBLE(&mut minin[ii]),
                &minname[ii],
                None,
                status,
            ) != 0
        {
            ffpmsg_str("error reading histogramming minimum keyword");
            ffpmsg_slice(&minname[ii]);
            return *status;
        }

        let cond = colexpr.is_none()
            || colexpr.as_ref().unwrap()[ii].is_none()
            || colexpr.as_ref().unwrap()[ii].as_ref().unwrap()[0] == 0;

        if minin[ii] != DOUBLENULLVALUE {
            amin[ii] = minin[ii];
        } else if cond {
            ffkeyn_safe(cs!(c"TLMIN"), colnum[ii], &mut keyname, status);
            if ffgky_safe(
                fptr,
                KeywordDatatypeMut::TDOUBLE(&mut amin[ii]),
                &keyname,
                None,
                status,
            ) > 0
            {
                /* use actual data minimum value for the histogram minimum */
                *status = 0;
                if fits_get_col_minmax(fptr, colnum[ii], &mut amin[ii], &mut datamax, status) > 0 {
                    strcpy_safe(
                        &mut errmsg,
                        cs!(c"Error calculating datamin and datamax for column: "),
                    );
                    let errmsg_len = strlen_safe(&errmsg);
                    strncat_safe(&mut errmsg, &colname[ii], FLEN_ERRMSG - errmsg_len - 1);
                    ffpmsg_slice(&errmsg);
                    return *status;
                }
            }
        } else {
            /* it's an expression */
            let ce = colexpr.unwrap();
            if fits_get_expr_minmax(
                fptr,
                ce[ii].unwrap(),
                Some(&mut amin[ii]),
                Some(&mut datamax),
                None,
                status,
            ) > 0
            {
                strcpy_safe(
                    &mut errmsg,
                    cs!(c"Error calculating datamin and datamax for expression: "),
                );
                ffpmsg_slice(&errmsg);
                ffpmsg_slice(ce[ii].unwrap());
                return *status;
            }
            if amin[ii] == DOUBLENULLVALUE {
                amin[ii] = 0.0;
            }
        }

        /* ================================================================ */
        /* get the maximum value */

        if maxname[ii][0] != 0
            && ffgky_safe(
                fptr,
                KeywordDatatypeMut::TDOUBLE(&mut maxin[ii]),
                &maxname[ii],
                None,
                status,
            ) != 0
        {
            ffpmsg_str("error reading histogramming maximum keyword");
            ffpmsg_slice(&maxname[ii]);
            return *status;
        }

        let cond = colexpr.is_none()
            || colexpr.as_ref().unwrap()[ii].is_none()
            || colexpr.as_ref().unwrap()[ii].as_ref().unwrap()[0] == 0;

        if maxin[ii] != DOUBLENULLVALUE {
            amax[ii] = maxin[ii];
        } else if cond {
            ffkeyn_safe(cs!(c"TLMAX"), colnum[ii], &mut keyname, status);
            if ffgky_safe(
                fptr,
                KeywordDatatypeMut::TDOUBLE(&mut amax[ii]),
                &keyname,
                None,
                status,
            ) > 0
            {
                *status = 0;
                if datamax != DOUBLENULLVALUE
                /* already computed max value */
                {
                    amax[ii] = datamax;
                } else {
                    /* use actual data maximum value for the histogram maximum */
                    if fits_get_col_minmax(fptr, colnum[ii], &mut datamin, &mut amax[ii], status)
                        > 0
                    {
                        strcpy_safe(
                            &mut errmsg,
                            cs!(c"Error calculating datamin and datamax for column: "),
                        );
                        let errmsg_len = strlen_safe(&errmsg);
                        strncat_safe(&mut errmsg, &colname[ii], FLEN_ERRMSG - errmsg_len - 1);
                        ffpmsg_slice(&errmsg);
                        return *status;
                    }
                }
            }
            use_datamax = 1; /* flag that the max was determined by the data values */
        /* and not specifically set by the calling program */
        } else {
            /* it's an expression */
            let ce = colexpr.unwrap();
            if fits_get_expr_minmax(
                fptr,
                ce[ii].unwrap(),
                Some(&mut datamin),
                Some(&mut amax[ii]),
                None,
                status,
            ) > 0
            {
                strcpy_safe(
                    &mut errmsg,
                    cs!(c"Error calculating datamin and datamax for expression: "),
                );
                ffpmsg_slice(&errmsg);
                ffpmsg_slice(ce[ii].unwrap());
                return *status;
            }
            if amax[ii] == DOUBLENULLVALUE {
                amin[ii] = 1.0;
            }
            use_datamax = 1;
        }

        /* ================================================================ */
        /* determine binning size and range                                 */

        if binname[ii][0] != 0
            && ffgky_safe(
                fptr,
                KeywordDatatypeMut::TDOUBLE(&mut binsizein[ii]),
                &binname[ii],
                None,
                status,
            ) != 0
        {
            ffpmsg_str("error reading histogramming binsize keyword");
            ffpmsg_slice(&binname[ii]);
            return *status;
        }

        if binsizein[ii] == 0. {
            ffpmsg_str("error: histogram binsize = 0");
            *status = ZERO_SCALE;
            return *status;
        }

        /* use TDBINn keyword or else 1 if bin size is not given */
        if binsizein[ii] != DOUBLENULLVALUE {
            binsize[ii] = binsizein[ii];
        } else {
            tstatus = 0;

            let cond = colexpr.is_none()
                || colexpr.as_ref().unwrap()[ii].is_none()
                || colexpr.as_ref().unwrap()[ii].as_ref().unwrap()[0] == 0;
            if cond {
                ffkeyn_safe(cs!(c"TDBIN"), colnum[ii], &mut keyname, &mut tstatus);
                ffgky_safe(
                    fptr,
                    KeywordDatatypeMut::TDOUBLE(&mut binsizein[ii]),
                    &keyname,
                    None,
                    &mut tstatus,
                );
            }

            let cond = colexpr.is_some()
                && colexpr.unwrap()[ii].is_some()
                && colexpr.unwrap()[ii].as_ref().unwrap()[0] != 0;
            if tstatus != 0 || cond {
                /* make at least 10 bins */
                binsize[ii] = (amax[ii] - amin[ii]) / 10.;
                if binsize[ii] > 1. {
                    binsize[ii] = 1.; /* use default bin size */
                }
            }
        }

        /* ================================================================ */
        /* if the min is greater than the max, make the binsize negative */
        if (amin[ii] > amax[ii] && binsize[ii] > 0.) || (amin[ii] < amax[ii] && binsize[ii] < 0.) {
            binsize[ii] = -binsize[ii]; /* reverse the sign of binsize */
        }

        if binsize[ii] == 0.0 {
            ffpmsg_str("error: computed histogram binsize = 0");

            let cond = colexpr.is_some()
                && colexpr.unwrap()[ii].is_some()
                && colexpr.unwrap()[ii].as_ref().unwrap()[0] != 0;
            if cond {
                ffpmsg_str("binning expression:");
                ffpmsg_slice(colexpr.unwrap()[ii].as_ref().unwrap());
            } else if colname[ii][0] != 0 {
                ffpmsg_str("binning column:");
                ffpmsg_slice(&colname[ii]);
            }
            *status = ZERO_SCALE;
            return *status;
        }

        ibin = binsize[ii] as c_int;
        imin = amin[ii] as c_int;
        imax = amax[ii] as c_int;

        /* Determine the range and number of bins in the histogram. This  */
        /* depends on whether the input columns are integer or floats, so */
        /* treat each case separately.                                    */

        if datatype <= TLONG
            && f64::from(imin) == amin[ii]
            && f64::from(imax) == amax[ii]
            && f64::from(ibin) == binsize[ii]
        {
            /* This is an integer column and integer limits were entered. */
            /* Shift the lower and upper histogramming limits by 0.5, so that */
            /* the values fall in the center of the bin, not on the edge. */

            haxes[ii] = c_long::from((imax - imin) / ibin + 1); /* last bin may only */
            /* be partially full */
            if amin[ii] < amax[ii] {
                amin[ii] -= 0.5;
                amax[ii] += 0.5;
            } else {
                amin[ii] += 0.5;
                amax[ii] -= 0.5;
            }
        } else if use_datamax != 0 {
            /* Either the column datatype and/or the limits are floating point, */
            /* and the histogram limits are being defined by the min and max */
            /* values of the array.  Add 1 to the number of histogram bins to */
            /* make sure that pixels that are equal to the maximum or are */
            /* in the last partial bin are included.  */

            haxes[ii] = (((amax[ii] - amin[ii]) / binsize[ii]) + 1.) as c_long;
        } else {
            /*  float datatype column and/or limits, and the maximum value to */
            /*  include in the histogram is specified by the calling program. */
            /*  The lower limit is inclusive, but upper limit is exclusive    */
            haxes[ii] = ((amax[ii] - amin[ii]) / binsize[ii]) as c_long;

            if amin[ii] < amax[ii] {
                if amin[ii] + (haxes[ii] as f64 * binsize[ii]) < amax[ii] {
                    haxes[ii] += 1; /* need to include another partial bin */
                }
            } else if amin[ii] + (haxes[ii] as f64 * binsize[ii]) > amax[ii] {
                haxes[ii] += 1; /* need to include another partial bin */
            }
        }
    }

    *status
}

/*--------------------------------------------------------------------------*/
/// Double-precision version, with non-extended syntax
/// Calculate the actual binning parameters, based on various user input options.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn fits_calc_binningd(
    fptr: *mut fitsfile, /* IO - pointer to table to be binned      ;       */
    naxis: c_int,        /* I - number of axes/columns in the binned image  */
    colname: *mut [[c_char; FLEN_VALUE]; 4], /* I - optional column names         */
    minin: *mut f64,     /* I - optional lower bound value for each axis  */
    maxin: *mut f64,     /* I - optional upper bound value, for each axis */
    binsizein: *mut f64, /* I - optional bin size along each axis         */
    minname: *const [[c_char; FLEN_VALUE]; 4], /* I - optional keywords for min       */
    maxname: *const [[c_char; FLEN_VALUE]; 4], /* I - optional keywords for max       */
    binname: *const [[c_char; FLEN_VALUE]; 4], /* I - optional keywords for binsize   */

    /* The returned parameters for each axis of the n-dimensional histogram are */
    colnum: *mut c_int, /* O - column numbers, to be binned */
    haxes: *mut c_long, /* O - number of bins in each histogram axis */
    amin: *mut f64,     /* O - lower bound of the histogram axes */
    amax: *mut f64,     /* O - upper bound of the histogram axes */
    binsize: *mut f64,  /* O - width of histogram bins/pixels on each axis */
    status: *mut c_int,
) -> c_int {
    unsafe {
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let colname = colname
            .cast::<[[c_char; FLEN_VALUE]; 4]>()
            .as_mut()
            .unwrap();
        let minin = slice::from_raw_parts_mut(minin.cast::<f64>(), naxis as usize);
        let maxin = slice::from_raw_parts_mut(maxin.cast::<f64>(), naxis as usize);
        let binsizein = slice::from_raw_parts_mut(binsizein.cast::<f64>(), naxis as usize);

        let minname = minname.as_ref().unwrap();
        let maxname = maxname.as_ref().unwrap();
        let binname = binname.as_ref().unwrap();
        let colnum = slice::from_raw_parts_mut(colnum, naxis as usize);
        let haxes = slice::from_raw_parts_mut(haxes, naxis as usize);
        let amin = slice::from_raw_parts_mut(amin, naxis as usize);
        let amax = slice::from_raw_parts_mut(amax, naxis as usize);
        let binsize = slice::from_raw_parts_mut(binsize, naxis as usize);
        let status = status.as_mut().expect(NULL_MSG);

        fits_calc_binningd_safe(
            fptr, naxis, colname, minin, maxin, binsizein, minname, maxname, binname, colnum,
            haxes, amin, amax, binsize, status,
        )
    }
}

/*--------------------------------------------------------------------------*/
/// Double-precision version, with non-extended syntax
/// Calculate the actual binning parameters, based on various user inputoptions.
pub fn fits_calc_binningd_safe(
    fptr: &mut fitsfile, /* IO - pointer to table to be binned      ;       */
    naxis: c_int,        /* I - number of axes/columns in the binned image  */
    colname: &mut [[c_char; FLEN_VALUE]; 4], /* I - optional column names         */
    minin: &mut [f64],   /* IO - optional lower bound value for each axis  */
    maxin: &mut [f64],   /* IO - optional upper bound value, for each axis */
    binsizein: &mut [f64], /* IO - optional bin size along each axis         */
    minname: &[[c_char; FLEN_VALUE]; 4], /* I - optional keywords for min       */
    maxname: &[[c_char; FLEN_VALUE]; 4], /* I - optional keywords for max       */
    binname: &[[c_char; FLEN_VALUE]; 4], /* I - optional keywords for binsize   */

    /* The returned parameters for each axis of the n-dimensional histogram are */
    colnum: &mut [c_int], /* O - column numbers, to be binned */
    haxes: &mut [c_long], /* O - number of bins in each histogram axis */
    amin: &mut [f64],     /* O - lower bound of the histogram axes */
    amax: &mut [f64],     /* O - upper bound of the histogram axes */
    binsize: &mut [f64],  /* O - width of histogram bins/pixels on each axis */
    status: &mut c_int,
) -> c_int {
    fits_calc_binningde(
        fptr, naxis, colname, None, minin, maxin, binsizein, minname, maxname, binname, colnum,
        None, haxes, amin, amax, binsize, None, status,
    )
}

/*--------------------------------------------------------------------------*/
/// Write default WCS keywords in the output histogram image header
/// if the keywords do not already exist.
pub(crate) fn fits_write_keys_histoe(
    fptr: &mut fitsfile,    /* I - pointer to table to be binned              */
    histptr: &mut fitsfile, /* I - pointer to output histogram image HDU      */
    naxis: c_int,           /* I - number of axes in the histogram image      */
    colnum: &[c_int],       /* I - column numbers (array length = naxis)      */
    colname: Option<&[[c_char; FLEN_VALUE]; 4]>, /* I - if expression, then column name to use */
    colexpr: Option<&[Option<&[c_char]>; 4]>, /* I - if expression, then column name to use */
    status: &mut c_int,
) -> c_int {
    let mut tstatus: c_int = 0;
    let mut keyname: [c_char; FLEN_KEYWORD] = [0; FLEN_KEYWORD];
    let mut svalue: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut dvalue: f64 = 0.0;

    if *status > 0 {
        return *status;
    }
    for ii in 0..(naxis as usize) {
        /*  CTYPEn  */
        tstatus = 0;

        if let Some(colexpr) = colexpr
            && let Some(colexpr_item) = colexpr[ii].as_ref()
            && colexpr_item[0] != 0
        {
            let colname = colname.as_ref().expect(NULL_MSG);

            // && colexpr[ii] && colexpr[ii][0] && colname[ii]) {
            /* Column expression: we need to put the column name from the binning expression */
            ffkeyn_safe(cs!(c"CTYPE"), (ii + 1) as c_int, &mut keyname, &mut tstatus);
            fits_write_key_str(
                histptr,
                &keyname,
                &colname[ii],
                Some(cs!(c"Coordinate Type")),
                &mut tstatus,
            );
        } else {
            /* Column name */
            tstatus = 0;
            ffkeyn_safe(cs!(c"CTYPE"), (ii + 1) as c_int, &mut keyname, &mut tstatus);
            fits_read_key_str(histptr, &keyname, &mut svalue, None, &mut tstatus);

            if tstatus == 0 {
                continue; /* keyword already exists, so skip to next axis */
            }

            /* use column name as the axis name */
            tstatus = 0;
            ffkeyn_safe(cs!(c"TTYPE"), colnum[ii], &mut keyname, &mut tstatus);
            fits_read_key_str(fptr, &keyname, &mut svalue, None, &mut tstatus);

            if tstatus == 0 {
                ffkeyn_safe(cs!(c"CTYPE"), (ii + 1) as c_int, &mut keyname, &mut tstatus);
                fits_write_key_str(
                    histptr,
                    &keyname,
                    &svalue,
                    Some(cs!(c"Coordinate Type")),
                    &mut tstatus,
                );
            }

            /*  CUNITn,  use the column units */
            tstatus = 0;
            ffkeyn_safe(cs!(c"TUNIT"), colnum[ii], &mut keyname, &mut tstatus);
            fits_read_key_str(fptr, &keyname, &mut svalue, None, &mut tstatus);

            if tstatus == 0 {
                ffkeyn_safe(cs!(c"CUNIT"), (ii + 1) as c_int, &mut keyname, &mut tstatus);
                fits_write_key_str(
                    histptr,
                    &keyname,
                    &svalue,
                    Some(cs!(c"Coordinate Units")),
                    &mut tstatus,
                );
            }
        }

        /*  CRPIXn  - Reference Pixel choose first pixel in new image as ref. pix. */
        dvalue = 1.0;
        tstatus = 0;
        ffkeyn_safe(cs!(c"CRPIX"), (ii + 1) as c_int, &mut keyname, &mut tstatus);
        ffpky_safe(
            histptr,
            KeywordDatatype::TDOUBLE(&dvalue),
            &keyname,
            Some(cs!(c"Reference Pixel")),
            &mut tstatus,
        );

        /*  CRVALn - Value at the location of the reference pixel */
        dvalue = 1.0;
        tstatus = 0;
        ffkeyn_safe(cs!(c"CRVAL"), (ii + 1) as c_int, &mut keyname, &mut tstatus);
        ffpky_safe(
            histptr,
            KeywordDatatype::TDOUBLE(&dvalue),
            &keyname,
            Some(cs!(c"Reference Value")),
            &mut tstatus,
        );

        /*  CDELTn - unit size of pixels  */
        dvalue = 1.0;
        tstatus = 0;
        dvalue = 1.0;
        ffkeyn_safe(cs!(c"CDELT"), (ii + 1) as c_int, &mut keyname, &mut tstatus);
        ffpky_safe(
            histptr,
            KeywordDatatype::TDOUBLE(&dvalue),
            &keyname,
            Some(cs!(c"Pixel size")),
            &mut tstatus,
        );
    }
    *status
}

/*--------------------------------------------------------------------------*/
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn fits_write_keys_histo(
    fptr: *mut fitsfile,    /* I - pointer to table to be binned              */
    histptr: *mut fitsfile, /* I - pointer to output histogram image HDU      */
    naxis: c_int,           /* I - number of axes in the histogram image      */
    colnum: *const c_int,   /* I - column numbers (array length = naxis)      */
    status: *mut c_int,
) -> c_int {
    unsafe {
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let histptr = histptr.as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);

        let colnum = slice::from_raw_parts(colnum, naxis as usize);

        fits_write_keys_histo_safe(fptr, histptr, naxis, colnum, status)
    }
}

/*--------------------------------------------------------------------------*/
pub fn fits_write_keys_histo_safe(
    fptr: &mut fitsfile,    /* I - pointer to table to be binned              */
    histptr: &mut fitsfile, /* I - pointer to output histogram image HDU      */
    naxis: c_int,           /* I - number of axes in the histogram image      */
    colnum: &[c_int],       /* I - column numbers (array length = naxis)      */
    status: &mut c_int,
) -> c_int {
    fits_write_keys_histoe(fptr, histptr, naxis, colnum, None, None, status)
}

/*--------------------------------------------------------------------------*/
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn fits_rebin_wcs(
    fptr: *mut fitsfile, /* I - pointer to table to be binned           */
    naxis: c_int,        /* I - number of axes in the histogram image   */
    amin: *mut f32,      /* I - first pixel include in each axis        */
    binsize: *mut f32,   /* I - binning factor for each axis            */
    status: *mut c_int,
) -> c_int {
    unsafe {
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let amin = slice::from_raw_parts_mut(amin, naxis as usize);
        let binsize = slice::from_raw_parts_mut(binsize, naxis as usize);
        let status = status.as_mut().expect(NULL_MSG);

        fits_rebin_wcs_safe(fptr, naxis, amin, binsize, status)
    }
}

/*--------------------------------------------------------------------------*/
pub fn fits_rebin_wcs_safe(
    fptr: &mut fitsfile, /* I - pointer to table to be binned           */
    naxis: c_int,        /* I - number of axes in the histogram image   */
    amin: &mut [f32],    /* I - first pixel include in each axis        */
    binsize: &mut [f32], /* I - binning factor for each axis            */
    status: &mut c_int,
) -> c_int {
    let mut amind: [f64; 4] = [0.0; 4];
    let mut binsized: [f64; 4] = [0.0; 4];

    /* Copy single precision values into double precision */
    if *status == 0 {
        let naxis1: usize = if naxis < 4 { naxis as usize } else { 4 };
        for i in 0..naxis1 {
            amind[i] = f64::from(amin[i]);
            binsized[i] = f64::from(binsize[i]);
        }

        fits_rebin_wcsd_safe(fptr, naxis, &mut amind, &mut binsized, status);
    }
    *status
}

/*--------------------------------------------------------------------------*/
/// Double precision version
/// Update the  WCS keywords that define the location of the reference
/// pixel, and the pixel size, along each axis.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn fits_rebin_wcsd(
    fptr: *mut fitsfile, /* I - pointer to table to be binned           */
    naxis: c_int,        /* I - number of axes in the histogram image   */
    amin: *mut f64,      /* I - first pixel include in each axis        */
    binsize: *mut f64,   /* I - binning factor for each axis            */
    status: *mut c_int,
) -> c_int {
    unsafe {
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let amin = slice::from_raw_parts_mut(amin, naxis as usize);
        let binsize = slice::from_raw_parts_mut(binsize, naxis as usize);
        let status = status.as_mut().expect(NULL_MSG);

        fits_rebin_wcsd_safe(fptr, naxis, amin, binsize, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Double precision version
/// Update the  WCS keywords that define the location of the reference
/// pixel, and the pixel size, along each axis.
pub fn fits_rebin_wcsd_safe(
    fptr: &mut fitsfile, /* I - pointer to table to be binned           */
    naxis: c_int,        /* I - number of axes in the histogram image   */
    amin: &mut [f64],    /* I - first pixel include in each axis        */
    binsize: &mut [f64], /* I - binning factor for each axis            */
    status: &mut c_int,
) -> c_int {
    let mut ii: c_int;
    let mut jj: c_int;
    let mut tstatus: c_int;
    let mut reset: c_int;
    let mut keyname: [c_char; FLEN_KEYWORD] = [0; FLEN_KEYWORD];
    let mut svalue: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut dvalue: f64 = 0.0;

    if *status > 0 {
        return *status;
    }

    for ii in 0..naxis as usize {
        reset = 0; /* flag to reset the reference pixel */
        tstatus = 0;
        ffkeyn_safe(cs!(c"CRVAL"), (ii + 1) as c_int, &mut keyname, &mut tstatus);
        /* get previous (pre-binning) value */
        ffgky_safe(
            fptr,
            KeywordDatatypeMut::TDOUBLE(&mut dvalue),
            &keyname,
            None,
            &mut tstatus,
        );
        if tstatus == 0 && dvalue == 1.0 {
            reset = 1;
        }

        tstatus = 0;
        /*  CRPIXn - update location of the ref. pix. in the binned image */
        ffkeyn_safe(cs!(c"CRPIX"), (ii + 1) as c_int, &mut keyname, &mut tstatus);

        /* get previous (pre-binning) value */
        ffgky_safe(
            fptr,
            KeywordDatatypeMut::TDOUBLE(&mut dvalue),
            &keyname,
            None,
            &mut tstatus,
        );

        if tstatus == 0 {
            if dvalue != 1.0 {
                reset = 0;
            }

            /* updated value to give pixel location after binning */
            dvalue = (dvalue - amin[ii]) / binsize[ii] + 0.5;

            fits_modify_key_dbl(fptr, &keyname, dvalue, -14, None, &mut tstatus);
        } else {
            reset = 0;
        }

        /*  CDELTn - update unit size of pixels  */
        tstatus = 0;
        ffkeyn_safe(cs!(c"CDELT"), (ii + 1) as c_int, &mut keyname, &mut tstatus);

        /* get previous (pre-binning) value */
        ffgky_safe(
            fptr,
            KeywordDatatypeMut::TDOUBLE(&mut dvalue),
            &keyname,
            None,
            &mut tstatus,
        );

        if tstatus == 0 {
            if dvalue != 1.0 {
                reset = 0;
            }

            /* updated to give post-binning value */
            dvalue *= binsize[ii];

            fits_modify_key_dbl(fptr, &keyname, dvalue, -14, None, &mut tstatus);
        } else {
            /* no CDELTn keyword, so look for a CDij keywords */
            reset = 0;

            for jj in 0..naxis {
                tstatus = 0;
                ffkeyn_safe(cs!(c"CD"), jj + 1, &mut svalue, &mut tstatus);
                strcat_safe(&mut svalue, cs!(c"_"));
                ffkeyn_safe(&svalue, (ii + 1) as c_int, &mut keyname, &mut tstatus);

                /* get previous (pre-binning) value */
                ffgky_safe(
                    fptr,
                    KeywordDatatypeMut::TDOUBLE(&mut dvalue),
                    &keyname,
                    None,
                    &mut tstatus,
                );

                if tstatus == 0 {
                    /* updated to give post-binning value */
                    dvalue *= binsize[ii];

                    fits_modify_key_dbl(fptr, &keyname, dvalue, -14, None, &mut tstatus);
                }
            }
        }

        if reset != 0 {
            /* the original CRPIX, CRVAL, and CDELT keywords were all = 1.0 */
            /* In this special case, reset the reference pixel to be the */
            /* first pixel in the array (instead of possibly far off the array) */

            dvalue = 1.0;
            ffkeyn_safe(cs!(c"CRPIX"), (ii + 1) as c_int, &mut keyname, &mut tstatus);
            fits_modify_key_dbl(fptr, &keyname, dvalue, -14, None, &mut tstatus);

            ffkeyn_safe(cs!(c"CRVAL"), (ii + 1) as c_int, &mut keyname, &mut tstatus);
            dvalue = amin[ii] + (binsize[ii] / 2.0);
            fits_modify_key_dbl(fptr, &keyname, dvalue, -14, None, &mut tstatus);
        }
    }
    *status
}

/*--------------------------------------------------------------------------*/
/// Single-precision version
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn fits_make_hist(
    fptr: *mut fitsfile,      /* IO - pointer to table with X and Y cols; */
    histptr: *mut fitsfile,   /* I - pointer to output FITS image      */
    bitpix: c_int,            /* I - datatype for image: 16, 32, -32, etc    */
    naxis: c_int,             /* I - number of axes in the histogram image   */
    naxes: *const c_long,     /* I - size of axes in the histogram image   */
    colnum: *const c_int,     /* I - column numbers (array length = naxis)   */
    amin: *const f32,         /* I - minimum histogram value, for each axis */
    amax: *const f32,         /* I - maximum histogram value, for each axis */
    binsize: *const f32,      /* I - bin size along each axis               */
    weight: f32,              /* I - binning weighting factor          */
    wtcolnum: c_int,          /* I - optional keyword or col for weight*/
    recip: c_int,             /* I - use reciprocal of the weight?     */
    selectrow: *const c_char, /* I - optional array (length = no. of   */
    /* rows in the table).  If the element is true */
    /* then the corresponding row of the table will*/
    /* be included in the histogram, otherwise the */
    /* row will be skipped.  Ingnored if *selectrow*/
    /* is equal to NULL.                           */
    status: *mut c_int,
) -> c_int {
    unsafe {
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let histptr = histptr.as_mut().expect(NULL_MSG);
        let naxes = slice::from_raw_parts(naxes, naxis as usize);
        let colnum = slice::from_raw_parts(colnum, naxis as usize);
        let amin = slice::from_raw_parts(amin, naxis as usize);
        let amax = slice::from_raw_parts(amax, naxis as usize);
        let binsize = slice::from_raw_parts(binsize, naxis as usize);
        let status = status.as_mut().expect(NULL_MSG);

        let mut nrows = 0;
        ffgnrw_safe(fptr, &mut nrows, status); /* no. of rows */

        let selectrow = if selectrow.is_null() {
            None
        } else {
            Some(slice::from_raw_parts(selectrow, nrows as usize))
        };

        fits_make_hist_safe(
            fptr, histptr, bitpix, naxis, naxes, colnum, amin, amax, binsize, weight, wtcolnum,
            recip, selectrow, status,
        )
    }
}

/*--------------------------------------------------------------------------*/
/// Single-precision version
pub fn fits_make_hist_safe(
    fptr: &mut fitsfile,          /* IO - pointer to table with X and Y cols; */
    histptr: &mut fitsfile,       /* I - pointer to output FITS image      */
    bitpix: c_int,                /* I - datatype for image: 16, 32, -32, etc    */
    naxis: c_int,                 /* I - number of axes in the histogram image   */
    naxes: &[c_long],             /* I - size of axes in the histogram image   */
    colnum: &[c_int],             /* I - column numbers (array length = naxis)   */
    amin: &[f32],                 /* I - minimum histogram value, for each axis */
    amax: &[f32],                 /* I - maximum histogram value, for each axis */
    binsize: &[f32],              /* I - bin size along each axis               */
    weight: f32,                  /* I - binning weighting factor          */
    wtcolnum: c_int,              /* I - optional keyword or col for weight*/
    recip: c_int,                 /* I - use reciprocal of the weight?     */
    selectrow: Option<&[c_char]>, /* I - optional array (length = no. of   */
    /* rows in the table).  If the element is true */
    /* then the corresponding row of the table will*/
    /* be included in the histogram, otherwise the */
    /* row will be skipped.  Ingnored if *selectrow*/
    /* is equal to NULL.                           */
    status: &mut c_int,
) -> c_int {
    let mut amind: [f64; 4] = [0.0; 4];
    let mut amaxd: [f64; 4] = [0.0; 4];
    let mut binsized: [f64; 4] = [0.0; 4];
    let weightd: f64 = f64::from(weight);

    /* Copy single precision values into double precision */
    if *status == 0 {
        let naxis1: usize = if naxis < 4 { naxis as usize } else { 4 };
        for i in 0..naxis1 {
            amind[i] = f64::from(amin[i]);
            amaxd[i] = f64::from(amax[i]);
            binsized[i] = f64::from(binsize[i]);
        }

        fits_make_histd_safe(
            fptr, histptr, bitpix, naxis, naxes, colnum, &amind, &amaxd, &binsized, weightd,
            wtcolnum, recip, selectrow, status,
        );
    }

    *status
}

/*--------------------------------------------------------------------------*/
/// Double-precision version
pub(crate) fn fits_make_histde(
    fptr: &mut fitsfile,    /* IO - pointer to table with X and Y cols; */
    histptr: &mut fitsfile, /* I - pointer to output FITS image      */
    mut datatypes: Option<&mut [c_int]>, /*  I - datatype of input (or 0 for auto) */
    bitpix: c_int,          /* I - datatype for image: 16, 32, -32, etc    */
    naxis: c_int,           /* I - number of axes in the histogram image   */
    naxes: &[c_long],       /* I - size of axes in the histogram image   */
    colnum: &[c_int],       /* I - column numbers (array length = naxis)   */
    colexpr: Option<&[Option<&[c_char]>; 4]>, /* I - optional expression instead of column */
    amin: &[f64],           /* I - minimum histogram value, for each axis */
    amax: &[f64],           /* I - maximum histogram value, for each axis */
    binsize: &[f64],        /* I - bin size along each axis               */
    mut weight: f64,        /* I - binning weighting factor (0 or DOUBLENULLVALUE means null) */
    wtcolnum: c_int,        /* I - optional keyword or col for weight*/
    wtexpr: Option<&[c_char]>, /* I - optional weighting expression */
    /*  disambiguation of weight values */
    /*    non-null weight: use that value */
    /*    null weight: use wtexpr if non-null, else wtcolnum */
    recip: c_int,                 /* I - use reciprocal of the weight?     */
    selectrow: Option<&[c_char]>, /* I - optional array (length = no. of   */
    /* rows in the table).  If the element is true */
    /* then the corresponding row of the table will*/
    /* be included in the histogram, otherwise the */
    /* row will be skipped.  Ingnored if *selectrow*/
    /* is equal to NULL.                           */
    status: &mut c_int,
) -> c_int {
    let mut ii: usize;
    let imagetype: c_int;
    let n_cols: c_int = 1;
    let mut imin: c_long;
    let mut imax: c_long;
    let mut ibin: c_long;
    let offset: c_long = 0;
    let n_per_loop: c_long = -1; /* force whole array to be passed at one time */
    let mut taxes: [f64; 4] = [0.0; 4];
    let mut tmin: [f64; 4] = [0.0; 4];
    let mut tmax: [f64; 4] = [0.0; 4];
    let mut tbin: [f64; 4] = [0.0; 4];
    let mut maxbin: [f64; 4] = [0.0; 4];
    let mut histData: HistType = HistType::default();
    let mut imagepars: [iteratorCol; 1] = [iteratorCol::default(); 1];
    let mut nrows: c_long = -1;
    let mut parsers: [ParseData; 5] = array::from_fn(|_| ParseData::default());
    let mut infos: [parseInfo; 5] = array::from_fn(|_| parseInfo::default());
    let mut numAllocCols: c_int = 0;
    let mut startCol: c_int = -1;
    let mut numIterCols: c_int = 0;
    let mut iterCols: *mut iteratorCol = core::ptr::null_mut::<iteratorCol>();
    let mut double_nulval: f64 = DOUBLENULLVALUE;
    let mut repeat: c_long = 0;
    let mut wtrepeat: c_long = 0;

    /* check inputs */

    if *status > 0 {
        return *status;
    }

    if naxis > 4 {
        ffpmsg_str("histogram has more than 4 dimensions");
        *status = BAD_DIMEN;
        return *status;
    }

    if bitpix == BYTE_IMG {
        imagetype = TBYTE;
    } else if bitpix == SHORT_IMG {
        imagetype = TSHORT;
    } else if bitpix == LONG_IMG {
        imagetype = TINT;
    } else if bitpix == FLOAT_IMG {
        imagetype = TFLOAT;
    } else if bitpix == DOUBLE_IMG {
        imagetype = TDOUBLE;
    } else {
        *status = BAD_DATATYPE;
        return *status;
    }

    /* reset position to the correct HDU if necessary */
    if (fptr).HDUposition != ((fptr).Fptr).curhdu {
        ffmahd_safe(fptr, ((fptr).HDUposition) + 1, None, status);
    }

    /* Resolve the conflict between wtexpr, wtcolnum, and weight */
    if ((wtcolnum > 0) || (wtexpr.is_some() && wtexpr.unwrap()[0] != 0)) && weight == 0.0 {
        weight = DOUBLENULLVALUE;
    }
    histData.weight = weight;
    histData.wtcolnum = wtcolnum;
    histData.wtexpr = wtexpr
        .map(|s| s.as_ptr() as *mut c_char)
        .unwrap_or(ptr::null_mut());
    histData.wtrecip = recip;
    histData.tblptr = fptr;
    histData.himagetype = imagetype;
    histData.haxis = naxis;
    histData.rowselector = selectrow
        .map(|s| s.as_ptr() as *mut c_char)
        .unwrap_or(ptr::null_mut());

    'cleanup: loop {
        /* Now make iterator columns for input, as well as any calculated values */
        numAllocCols = 5;
        iterCols = unsafe {
            fits_recalloc(
                ptr::null_mut(),
                0,
                numAllocCols as usize,
                core::mem::size_of::<iteratorCol>(),
            )
            .cast::<iteratorCol>()
        };
        if iterCols.is_null() {
            ffpmsg_str("memory allocation failure (fits_make_histde)");
            *status = MEMORY_ALLOCATION;
            break 'cleanup;
        }

        /* We fill the iterCols in order, starting from column 1 through 4, and
        then moving on to the weighting column */
        for ii in 0..5 {
            histData.startCols[ii] = -1;
        }
        startCol = 0;

        /* Loop through each axis and recheck the binning parameters */
        ii = 0;
        while ii < naxis as usize {
            let mut colrepeat: c_long = 0;
            let mut datatype: c_int = 0;
            histData.startCols[ii] = startCol;

            taxes[ii] = naxes[ii] as f64;
            tmin[ii] = amin[ii];
            tmax[ii] = amax[ii];
            if (amin[ii] > amax[ii] && binsize[ii] > 0.)
                || (amin[ii] < amax[ii] && binsize[ii] < 0.)
            {
                tbin[ii] = -binsize[ii]; /* reverse the sign of binsize */
            } else {
                tbin[ii] = binsize[ii]; /* binsize has the correct sign */
            }

            imin = tmin[ii] as c_long;
            imax = tmax[ii] as c_long;
            ibin = tbin[ii] as c_long;

            /* get the datatype of the column and repeat */
            let cond = colexpr.is_some()
                && colexpr.unwrap()[ii].is_some()
                && colexpr.unwrap()[ii].as_ref().unwrap()[0] != 0;
            if !cond {
                fits_get_eqcoltype(
                    fptr,
                    colnum[ii],
                    Some(&mut datatype),
                    Some(&mut colrepeat),
                    None,
                    status,
                );
            }

            /* If caller specified datatype, use that */
            if let Some(dt) = datatypes.as_deref_mut()
                && dt[ii] != 0
            {
                datatype = dt[ii];
            }

            if datatype <= TLONG
                && (imin as f64) == tmin[ii]
                && (imax as f64) == tmax[ii]
                && (ibin as f64) == tbin[ii]
            {
                /* This is an integer column and integer limits were entered. */
                /* Shift the lower and upper histogramming limits by 0.5, so that */
                /* the values fall in the center of the bin, not on the edge. */

                maxbin[ii] = taxes[ii] + 1.0; /* add 1. instead of 0.5 to avoid roundoff */

                if tmin[ii] < tmax[ii] {
                    tmin[ii] -= 0.5;
                    tmax[ii] += 0.5;
                } else {
                    tmin[ii] += 0.5;
                    tmax[ii] -= 0.5;
                }
            } else {
                /* not an integer column with integer limits */
                maxbin[ii] = (tmax[ii] - tmin[ii]) / tbin[ii];
            }

            /* This is a column expression.  Here is where we allocate the
            parser for it during the actual evaluation. */
            let cond = colexpr.is_some()
                && colexpr.unwrap()[ii].is_some()
                && colexpr.unwrap()[ii].as_ref().unwrap()[0] != 0;
            if cond {
                let mut datatype: c_int = 0;
                let mut naxis1: c_int = 0;
                let mut nelem: c_long = 0;
                let mut naxes: [c_long; MAXDIMS as usize] = [0; MAXDIMS as usize];
                let jj: c_int = 0;

                /* Initialize the parser for this binning expression */
                let ce = colexpr.unwrap();
                ffiprs(
                    fptr,
                    0,
                    ce[ii].unwrap(),
                    MAXDIMS,
                    &mut datatype,
                    &mut nelem,
                    &mut naxis1,
                    &mut naxes,
                    &mut (parsers[ii]),
                    status,
                );
                if *status != 0 {
                    break 'cleanup;
                }
                if nelem < 0 {
                    nelem = 1; /* If it's a constant expression */
                }

                colrepeat = nelem;

                /* Set up the parser data for evaluation to a TemporaryCol */
                fits_get_num_rows(fptr, &mut nrows, status);
                if fits_parser_set_temporary_col(
                    &mut (parsers[ii]),
                    &mut (infos[ii]),
                    nrows,
                    (&mut (double_nulval) as *mut f64).cast::<c_void>(),
                    status,
                ) != 0
                {
                    break 'cleanup;
                }

                /* Copy iterator columns from the parser to the master iterator columns */
                iterCols = unsafe {
                    fits_recalloc(
                        iterCols.cast::<c_void>(),
                        numAllocCols as usize,
                        (numAllocCols + parsers[ii].nCols) as usize,
                        core::mem::size_of::<iteratorCol>(),
                    ) as *mut iteratorCol
                };
                if iterCols.is_null() {
                    *status = MEMORY_ALLOCATION;
                    break 'cleanup;
                }
                numAllocCols += parsers[ii].nCols;
                for jj in 0..parsers[ii].nCols as usize {
                    unsafe { *iterCols.add(startCol as usize) = parsers[ii].colData[jj] };
                    startCol += 1;
                }
            } else {
                /* Just a "regular" column name, we already have enough allocated for these */
                fits_iter_set_by_num_safe(
                    unsafe { &mut *iterCols.add(startCol as usize) },
                    fptr,
                    colnum[ii],
                    TDOUBLE,
                    INPUT_COL,
                );
                startCol += 1;
            }

            /* Check that all the vector dimensions agree */
            if repeat == 0 {
                repeat = colrepeat;
            } else if repeat != colrepeat {
                ffpmsg_str("vector dimensions of binning values do not agree");
                *status = BAD_DIMEN;
                break 'cleanup;
            }
            ii += 1;
        } /* End of loop over columns */

        /* Now initialize the iterator column data for the weighting */
        if wtexpr.is_some() && wtexpr.unwrap()[0] != 0 && weight == DOUBLENULLVALUE {
            let mut wtdatatype: c_int = 0;
            let mut wtnaxis: c_int = 0;
            let mut wtnaxes: [c_long; MAXDIMS as usize] = [0; MAXDIMS as usize];
            let jj: c_int = 0;

            histData.startCols[4] = startCol;
            ffiprs(
                fptr,
                0,
                wtexpr.unwrap(),
                MAXDIMS,
                &mut wtdatatype,
                &mut wtrepeat,
                &mut wtnaxis,
                &mut wtnaxes,
                &mut (parsers[4]),
                status,
            );
            if *status != 0 {
                break 'cleanup;
            }
            if wtrepeat < 0 {
                wtrepeat = 1; /* If it's a constant expression */
            }

            /* Set up the parser data for evaluation to a TemporaryCol */
            /* It's a weighting expression, set that up and ... */
            fits_get_num_rows(fptr, &mut nrows, status);
            if fits_parser_set_temporary_col(
                &mut (parsers[4]),
                &mut (infos[4]),
                nrows,
                (&mut (double_nulval) as *mut f64).cast::<c_void>(),
                status,
            ) != 0
            {
                break 'cleanup;
            }

            /* Copy iterator columns from the parser to the master iterator columns */
            iterCols = unsafe {
                fits_recalloc(
                    iterCols.cast::<c_void>(),
                    numAllocCols as usize,
                    (numAllocCols + parsers[4].nCols) as usize,
                    core::mem::size_of::<iteratorCol>(),
                ) as *mut iteratorCol
            };
            if iterCols.is_null() {
                *status = MEMORY_ALLOCATION;
                break 'cleanup;
            }
            numAllocCols += parsers[ii].nCols;
            for jj in 0..parsers[4].nCols as usize {
                unsafe { *iterCols.add(startCol as usize) = parsers[4].colData[jj] };
                startCol += 1;
            }
        } else if weight == DOUBLENULLVALUE {
            let mut wtdatatype: c_int = 0;

            /* It's a "regular" weighting column */
            fits_get_eqcoltype(
                fptr,
                wtcolnum,
                Some(&mut wtdatatype),
                Some(&mut wtrepeat),
                None,
                status,
            );

            histData.startCols[4] = startCol;
            fits_iter_set_by_num_safe(
                unsafe { &mut *iterCols.add(startCol as usize) },
                fptr,
                wtcolnum,
                TDOUBLE,
                INPUT_COL,
            );
            startCol += 1;
        } else {
            /* In case of explicit numerical value, we can just use that number
            in the vector expression, so the vector repeat of the weighting can
            be set to that of the input */
            wtrepeat = repeat;
        }

        /* Vector dimension of weighting must agree with binning */
        if wtrepeat != 0 && repeat != 0 && wtrepeat != repeat {
            ffpmsg_str("vector dimensions of weights do not agree with bins");
            *status = BAD_DIMEN;
            break 'cleanup;
        }

        /* We now know he number of iterator columns */
        numIterCols = startCol;

        /* Fill in iterator information for the parser*/
        histData.numIterCols = numIterCols;
        histData.iterCols = unsafe { slice::from_raw_parts_mut(iterCols, numIterCols as usize) };
        histData.parsers = &mut parsers;
        histData.infos = &mut infos;
        histData.repeat = repeat;

        /* Set global variables with histogram parameter values.    */
        /* Use separate scalar variables rather than arrays because */
        /* it is more efficient when computing the histogram.       */

        histData.hcolnum[0] = colnum[0];
        histData.amin1 = tmin[0];
        histData.maxbin1 = maxbin[0];
        histData.binsize1 = tbin[0];
        histData.haxis1 = taxes[0] as c_long;
        histData.incr[0] = 1;

        if histData.haxis > 1 {
            histData.hcolnum[1] = colnum[1];
            histData.amin2 = tmin[1];
            histData.maxbin2 = maxbin[1];
            histData.binsize2 = tbin[1];
            histData.haxis2 = taxes[1] as c_long;
            histData.incr[1] = histData.incr[0] * histData.haxis1;

            if histData.haxis > 2 {
                histData.hcolnum[2] = colnum[2];
                histData.amin3 = tmin[2];
                histData.maxbin3 = maxbin[2];
                histData.binsize3 = tbin[2];
                histData.haxis3 = taxes[2] as c_long;
                histData.incr[2] = histData.incr[1] * histData.haxis2;

                if histData.haxis > 3 {
                    histData.hcolnum[3] = colnum[3];
                    histData.amin4 = tmin[3];
                    histData.maxbin4 = maxbin[3];
                    histData.binsize4 = tbin[3];
                    histData.haxis4 = taxes[3] as c_long;
                    histData.incr[3] = histData.incr[2] * histData.haxis3;
                }
            }
        }

        /* define parameters of image for the iterator function */
        let histptr_mut = histptr;
        fits_iter_set_file_safe(&mut imagepars[0], histptr_mut); /* pointer to image */
        fits_iter_set_datatype_safe(&mut imagepars[0], imagetype); /* image datatype   */
        fits_iter_set_iotype_safe(&mut imagepars[0], OUTPUT_COL); /* image is output  */

        /* call the iterator function to write out the histogram image */
        fits_iterate_data(
            n_cols,
            &mut imagepars,
            offset,
            n_per_loop,
            ffwritehisto,
            (&mut histData as *mut HistType).cast::<c_void>(),
            status,
        );

        break 'cleanup;
    } // End of cleanup loop

    // cleanup:
    /* Free any allocated memory ... */
    if !iterCols.is_null() {
        // SAFETY / TODO
        // unsafe { free(iterCols as *mut c_void) };
    }
    /* ... and parsers */
    for ii in 0..=4 {
        if parsers[ii].nCols > 0 {
            ffcprs(&mut (parsers[ii]));
        }
    }
    *status
}

/*--------------------------------------------------------------------------*/
/// Double-precision version, non-extended syntax
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn fits_make_histd(
    fptr: *mut fitsfile,      /* IO - pointer to table with X and Y cols; */
    histptr: *mut fitsfile,   /* I - pointer to output FITS image      */
    bitpix: c_int,            /* I - datatype for image: 16, 32, -32, etc    */
    naxis: c_int,             /* I - number of axes in the histogram image   */
    naxes: *const c_long,     /* I - size of axes in the histogram image   */
    colnum: *const c_int,     /* I - column numbers (array length = naxis)   */
    amin: *const f64,         /* I - minimum histogram value, for each axis */
    amax: *const f64,         /* I - maximum histogram value, for each axis */
    binsize: *const f64,      /* I - bin size along each axis               */
    weight: f64,              /* I - binning weighting factor          */
    wtcolnum: c_int,          /* I - optional keyword or col for weight*/
    recip: c_int,             /* I - use reciprocal of the weight?     */
    selectrow: *const c_char, /* I - optional array (length = no. of   */
    /* rows in the table).  If the element is true */
    /* then the corresponding row of the table will*/
    /* be included in the histogram, otherwise the */
    /* row will be skipped.  Ingnored if *selectrow*/
    /* is equal to NULL.                           */
    status: *mut c_int,
) -> c_int {
    unsafe {
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let histptr = histptr.as_mut().expect(NULL_MSG);
        let naxes = slice::from_raw_parts(naxes, naxis as usize);
        let colnum = slice::from_raw_parts(colnum, naxis as usize);
        let amin = slice::from_raw_parts(amin, naxis as usize);
        let amax = slice::from_raw_parts(amax, naxis as usize);
        let binsize = slice::from_raw_parts(binsize, naxis as usize);
        let status = status.as_mut().expect(NULL_MSG);

        let mut nrows = 0;
        ffgnrw_safe(fptr, &mut nrows, status); /* no. of rows */

        let selectrow = if selectrow.is_null() {
            None
        } else {
            Some(slice::from_raw_parts(selectrow, nrows as usize))
        };

        fits_make_histd_safe(
            fptr, histptr, bitpix, naxis, naxes, colnum, amin, amax, binsize, weight, wtcolnum,
            recip, selectrow, status,
        )
    }
}

/*--------------------------------------------------------------------------*/
/// Double-precision version, non-extended syntax
pub fn fits_make_histd_safe(
    fptr: &mut fitsfile,          /* IO - pointer to table with X and Y cols; */
    histptr: &mut fitsfile,       /* I - pointer to output FITS image      */
    bitpix: c_int,                /* I - datatype for image: 16, 32, -32, etc    */
    naxis: c_int,                 /* I - number of axes in the histogram image   */
    naxes: &[c_long],             /* I - size of axes in the histogram image   */
    colnum: &[c_int],             /* I - column numbers (array length = naxis)   */
    amin: &[f64],                 /* I - minimum histogram value, for each axis */
    amax: &[f64],                 /* I - maximum histogram value, for each axis */
    binsize: &[f64],              /* I - bin size along each axis               */
    weight: f64,     /* I - binning weighting factor (0 or DOUBLENULLVALUE means null) */
    wtcolnum: c_int, /* I - optional keyword or col for weight*/
    recip: c_int,    /* I - use reciprocal of the weight?     */
    selectrow: Option<&[c_char]>, /* I - optional array (length = no. of   */
    /* rows in the table).  If the element is true */
    /* then the corresponding row of the table will*/
    /* be included in the histogram, otherwise the */
    /* row will be skipped.  Ingnored if *selectrow*/
    /* is equal to NULL.                           */
    status: &mut c_int,
) -> c_int {
    fits_make_histde(
        fptr, histptr, None, bitpix, naxis, naxes, colnum, None, amin, amax, binsize, weight,
        wtcolnum, None, recip, selectrow, status,
    )
}

/*--------------------------------------------------------------------------*/
/// Simple utility routine to compute the min and max value in a column
pub(crate) fn fits_get_col_minmax(
    fptr: &mut fitsfile,
    colnum: c_int,
    datamin: &mut f64,
    datamax: &mut f64,
    status: &mut c_int,
) -> c_int {
    let mut anynul: c_int = 0;
    let mut nrows: c_long = 0;
    let mut ntodo: c_long = 0;
    let mut firstrow: c_long = 1;
    let mut array: [f64; 1000] = [0.0; 1000];
    let nulval = NullValue::Double(DOUBLENULLVALUE);

    ffgky_safe(
        fptr,
        crate::KeywordDatatypeMut::TLONG(&mut nrows),
        cs!(c"NAXIS2"),
        None,
        status,
    ); /* no. of rows */

    *datamin = 9.0E36;
    *datamax = -9.0E36;

    while nrows > 0 {
        ntodo = c_long::min(nrows, 100);
        ffgcv_safe(
            fptr,
            TDOUBLE,
            colnum,
            firstrow as LONGLONG,
            1,
            ntodo as LONGLONG,
            Some(nulval.clone()),
            cast_slice_mut(&mut array),
            Some(&mut anynul),
            status,
        );

        for ii in 0..ntodo {
            if array[ii as usize] != nulval.get_value_as_f64() {
                *datamin = f64::min(*datamin, array[ii as usize]);
                *datamax = f64::max(*datamax, array[ii as usize]);
            }
        }

        nrows -= ntodo;
        firstrow += ntodo;
    }

    *status
}

#[derive(Default)]
struct histo_minmax_workfn_struct {
    Info: *mut parseInfo,
    datamin: f64,
    datamax: f64,
    ntotal: c_long,
    ngood: c_long,
}

/*---------------------------------------------------------------------------*/
/// Iterator work function which evaluates a parser result and computes
/// min max value
extern "C" fn histo_minmax_expr_workfn(
    totalrows: c_long,         /* I - Total rows to be processed     */
    offset: c_long,            /* I - Number of rows skipped at start*/
    firstrow: c_long,          /* I - First row of this iteration    */
    nrows: c_long,             /* I - Number of rows in this iter    */
    nCols: c_int,              /* I - Number of columns in use       */
    colData: *mut iteratorCol, /* IO- Column information/data        */
    userPtr: *mut c_void,      /* I - Data handling instructions     */
) -> c_int {
    let colData = unsafe { slice::from_raw_parts_mut(colData, nCols as usize) };
    let userPtr = unsafe { &mut *userPtr.cast::<histo_minmax_workfn_struct>() };
    let mut status: c_int = 0;
    let mut i: c_long;
    let mut data: *mut f64;

    let wf = userPtr;
    let pv: &mut ParseStatusVariables = unsafe { &mut ((*(wf.Info)).parseVariables) };

    /* Call calculator work function.  Result is put in final column of colData as a TemporaryCol */
    status = fits_parser_workfn_safe(
        totalrows,
        offset,
        firstrow,
        nrows,
        nCols,
        colData,
        wf.Info.cast::<c_void>(),
    );

    let outcol: &mut iteratorCol = &mut (colData[nCols as usize - 1]);

    /* The result of the calculation is in pv->Data, and null value in pv->Null */
    let data = unsafe {
        slice::from_raw_parts(outcol.array.cast::<f64>(), (nrows * pv.repeat + 1) as usize)
    };
    let nulval: f64 = unsafe { *(*(wf.Info)).nullPtr.cast::<f64>() };

    for i in 1..=(nrows * pv.repeat) as usize {
        /* Note that data[0] == 0 indicates no null values at all!!! */
        if data[0] == 0.0 || data[i] != nulval {
            if data[i] < wf.datamin || wf.datamin == DOUBLENULLVALUE {
                wf.datamin = data[i];
            }
            if data[i] > wf.datamax || wf.datamax == DOUBLENULLVALUE {
                wf.datamax = data[i];
            }
            wf.ngood += 1;
        }
        wf.ntotal += 1;
    }

    status
}

/*--------------------------------------------------------------------------*/
/// Simple utility routine to compute the min and max value in an expression
fn fits_get_expr_minmax(
    fptr: &mut fitsfile,
    expr: &[c_char],
    datamin: Option<&mut f64>,
    datamax: Option<&mut f64>,
    mut datatype: Option<&mut c_int>,
    status: &mut c_int,
) -> c_int {
    let mut Info: parseInfo = parseInfo::default();
    let mut lParse: ParseData = ParseData::default();
    let mut minmaxWorkFn: histo_minmax_workfn_struct = histo_minmax_workfn_struct::default();
    let mut naxis: c_int = 0;
    let constant: c_int = 0;
    let typecode: c_int = 0;
    let newNullKwd: c_int = 0;
    let mut nelem: c_long = 0;
    let mut naxes: [c_long; MAXDIMS as usize] = [0; MAXDIMS as usize];
    let repeat: c_long = 0;
    let width: c_long = 0;
    let mut nrows: c_long = 0;
    let col_cnt: c_int = 0;
    let colNo: c_int = 0;
    let result: &mut Node;
    let card: [c_char; 81] = [0; 81];
    let tform: [c_char; 16] = [0; 16];
    let nullKwd: [c_char; 9] = [0; 9];
    let tdimKwd: [c_char; 9] = [0; 9];

    let mut double_nulval: f64 = DOUBLENULLVALUE;

    if *status != 0 {
        return *status;
    }

    if let Some(ref mut dt) = datatype {
        **dt = 0;
    }

    ffgky_safe(
        fptr,
        KeywordDatatypeMut::TLONG(&mut nrows),
        cs!(c"NAXIS2"),
        None,
        status,
    ); /* no. of rows */

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

    if let Some(ref mut dt) = datatype {
        **dt = Info.datatype;
    }

    if nelem < 0 {
        /* Constant already computed */
        result = &mut lParse.Nodes[lParse.resultNode as usize];
        match Info.datatype {
            TDOUBLE => {
                let tmp = unsafe { result.value.data.dbl };
                if let Some(dm) = datamax {
                    *dm = tmp;
                }
                if let Some(dm) = datamin {
                    *dm = tmp;
                }
            }
            TLONG => {
                let tmp = unsafe { result.value.data.lng as f64 };
                if let Some(dm) = datamax {
                    *dm = tmp;
                }
                if let Some(dm) = datamin {
                    *dm = tmp;
                }
            }
            TLOGICAL => {
                let tmp = unsafe { if result.value.data.log == 1 { 1.0 } else { 0.0 } };
                if let Some(dm) = datamax {
                    *dm = tmp;
                }
                if let Some(dm) = datamin {
                    *dm = tmp;
                }
            }
            TBIT => {
                let tmp = unsafe {
                    if result.value.data.astr[0] != 0 {
                        1.0
                    } else {
                        0.0
                    }
                };
                if let Some(dm) = datamax {
                    *dm = tmp;
                }
                if let Some(dm) = datamin {
                    *dm = tmp;
                }
            }
            _ => {
                // For any other data types, we don't handle constants
                ffcprs(&mut lParse);
                return *status;
            }
        }
        ffcprs(&mut lParse);
        return *status;
    }

    Info.parseData = &mut lParse;

    /* Add a temporary column which contains the expression value */
    if fits_parser_set_temporary_col(
        &mut lParse,
        &mut Info,
        nrows,
        (&mut double_nulval as *mut f64).cast::<c_void>(),
        status,
    ) != 0
    {
        ffcprs(&mut lParse);
        return *status;
    }

    /* Initialize the work function computing min/max */
    minmaxWorkFn.Info = &mut Info;
    minmaxWorkFn.datamax = DOUBLENULLVALUE;
    minmaxWorkFn.datamin = DOUBLENULLVALUE;
    minmaxWorkFn.ngood = 0;
    minmaxWorkFn.ntotal = 0;

    if ffiter_safe(
        lParse.nCols,
        &mut lParse.colData,
        0,
        0,
        histo_minmax_expr_workfn
            as extern "C" fn(
                c_long,
                c_long,
                c_long,
                c_long,
                c_int,
                *mut iteratorCol,
                *mut c_void,
            ) -> c_int,
        (&mut minmaxWorkFn as *mut histo_minmax_workfn_struct).cast::<c_void>(),
        status,
    ) == -1
    {
        *status = 0; /* -1 indicates exitted without error before end... OK */
    }

    if let Some(datamin) = datamin {
        *datamin = minmaxWorkFn.datamin;
    }
    if let Some(datamax) = datamax {
        *datamax = minmaxWorkFn.datamax;
    }

    ffcprs(&mut lParse);
    *status
}

/*--------------------------------------------------------------------------*/
/// Interator work function that writes out the histogram.
/// The histogram values are calculated by another work function, ffcalchisto.
/// This work function only gets called once, and totaln = nvalues.
extern "C" fn ffwritehisto(
    totaln: c_long,
    pixoffset: c_long,
    firstn: c_long,
    nvalues: c_long,
    narrays: c_int,
    imagepars: *mut iteratorCol,
    userPointer: *mut c_void,
) -> c_int {
    unsafe {
        let imagepars = imagepars.as_mut().expect(NULL_MSG);

        let userPointer = userPointer.cast::<HistType>();
        let userPointer = userPointer.as_mut().expect(NULL_MSG);

        ffwritehisto_safe(
            totaln,
            pixoffset,
            firstn,
            nvalues,
            narrays,
            imagepars,
            userPointer,
        )
    }
}

/*--------------------------------------------------------------------------*/
/// Interator work function that writes out the histogram.
/// The histogram values are calculated by another work function, ffcalchisto.
/// This work function only gets called once, and totaln = nvalues.
fn ffwritehisto_safe(
    totaln: c_long,
    pixoffset: c_long,
    firstn: c_long,
    nvalues: c_long,
    narrays: c_int,
    imagepars: &mut iteratorCol,
    userPointer: &mut HistType,
) -> c_int {
    let mut colpars: &mut [iteratorCol];
    let mut ii: c_int;
    let mut status: c_int = 0;
    let mut ncols: c_int;
    let rows_per_loop: c_long = 0;
    let offset: c_long = 0;

    let histData = userPointer;

    /* store pointer to the histogram array, and initialize to zero */

    match histData.himagetype {
        TBYTE => {
            histData.hist.b = fits_iter_get_array_safe(imagepars) as *mut c_char;
        }
        TSHORT => {
            histData.hist.i = fits_iter_get_array_safe(imagepars) as *mut c_short;
        }
        TINT => {
            histData.hist.j = fits_iter_get_array_safe(imagepars) as *mut c_int;
        }
        TFLOAT => {
            histData.hist.r = fits_iter_get_array_safe(imagepars) as *mut f32;
        }
        TDOUBLE => {
            histData.hist.d = fits_iter_get_array_safe(imagepars) as *mut f64;
        }
        _ => {
            // Unsupported type
        }
    }

    /* call iterator function to calc the histogram pixel values */

    /* must lock this call in multithreaded environoments because */
    /* the ffcalchist work routine uses static vaiables that would */
    /* get clobbered if multiple threads were running at the same time */
    let n_cols = histData.numIterCols;

    // TODO / WARNING this is super unsafe
    unsafe {
        let iterCols = histData.iterCols.as_mut_ptr();
        let iterCols = slice::from_raw_parts_mut(iterCols, n_cols as usize);

        ffiter_safe(
            n_cols,
            iterCols,
            offset,
            rows_per_loop,
            ffcalchist
                as extern "C" fn(
                    c_long,
                    c_long,
                    c_long,
                    c_long,
                    c_int,
                    *mut iteratorCol,
                    *mut c_void,
                ) -> c_int,
            (histData as *mut HistType).cast::<c_void>(),
            &mut status,
        );
    }

    status
}

/*--------------------------------------------------------------------------*/
/// Interator work function that calculates values for the 2D histogram.
extern "C" fn ffcalchist(
    totalrows: c_long,
    offset: c_long,
    firstrow: c_long,
    nrows: c_long,
    ncols: c_int,
    colpars: *mut iteratorCol,
    userPointer: *mut c_void,
) -> c_int {
    let colpars = unsafe { &mut *colpars };
    let userPointer = unsafe { &mut *userPointer.cast::<HistType>() };
    let mut ii: c_long;
    let mut ipix: c_long;
    let mut iaxisbin: c_long;
    let mut pix: f64;
    let mut axisbin: f64;
    let mut rowselect: *const c_char;
    let mut colptr: [*mut f64; 5] = [ptr::null_mut(); 5];
    let mut status: c_int = 0;
    let irow: c_long = 0;
    let mut adjustedRepeat: c_long = 0;

    let histData = userPointer;

    if firstrow == 1 {
        histData.rowselector_cur = histData.rowselector;
    }

    rowselect = histData.rowselector_cur;

    for ii in 0..=4 {
        let startCol: c_int = histData.startCols[ii];
        let mut outcol: Option<&mut iteratorCol> = None;
        /* Call calculator work function.  Result is put in final column of colData as a TemporaryCol */
        colptr[ii] = core::ptr::null_mut();

        /* Do not process unspecified axes (but do process weight column) */
        if (ii >= histData.haxis as usize && ii != 4) || histData.startCols[ii] < 0 {
            continue;
        }

        /* We have a parser for this, evaluate it */
        if histData.parsers[ii].nCols > 0 {
            let nCols: c_int = histData.parsers[ii].nCols;
            let colData_slice: &mut [iteratorCol] =
                &mut histData.iterCols[startCol as usize..][..nCols as usize];

            status = fits_parser_workfn_safe(
                totalrows,
                offset,
                firstrow,
                nrows,
                nCols,
                colData_slice,
                ((&mut (histData.infos[ii])) as *mut parseInfo).cast::<c_void>(),
            );
            if status != 0 {
                return status;
            }
            /* Output column is last iterator column, which better be a TemporaryCol */
            outcol = Some(&mut histData.iterCols[(startCol + nCols - 1) as usize]);
        } else {
            outcol = Some(&mut histData.iterCols[startCol as usize]);
        }

        if outcol.is_some() {
            /* Note that the 0th array element returned by the iterator is
            actually the null value!  This is actually rather a big
            undocumented "feature" of the iterator. However, "ii" below
            starts at a value of 1 which skips over the null value */
            colptr[ii] = fits_iter_get_array_safe(outcol.unwrap()) as *mut f64;
        }
    }

    /* Main loop over rows
        For tables:
         irow = row counter (1 .. nrows)
         elem = counter of element (1 .. histData->repeat) for each row
         ii = counts up from 1 (see note below) used to index colptr[]'s
        For images:
         irow = pixel counter (1 .. totalnpix)
         elem = 1  (not applicable)
    */
    if !histData.tblptr.is_null() && (unsafe { &*histData.tblptr }.Fptr.hdutype != IMAGE_HDU) {
        adjustedRepeat = histData.repeat;
    } else {
        adjustedRepeat = 1;
    }

    /* Note that ii starts at 1 because position [0] in the
    column data arrays is for the "null" value! */
    let mut ii: usize = 1;
    for irow in 1..=nrows {
        // for (ii = 1, irow = 1; irow <= nrows; irow++) {

        if !rowselect.is_null() {
            /* if a row selector array is supplied... */

            if unsafe { *rowselect } != 0 {
                rowselect = unsafe { rowselect.add(1) }; /* this row is included in the histogram */
            } else {
                rowselect = unsafe { rowselect.add(1) }; /* this row is excluded from the histogram */

                ii += adjustedRepeat as usize; /* skip this portion of data */
                continue;
            }
        }

        /* Loop over elements in each row, increment ii after each element */
        for elem in 1..=adjustedRepeat {
            // for (elem = 1; elem <= histData.repeat; elem++, ii++) {
            if unsafe { *colptr[0].add(ii) } == DOUBLENULLVALUE {
                /* test for null value */
                ii += 1;
                continue;
            }
            if !colptr[4].is_null() && unsafe { *colptr[4].add(ii) } == DOUBLENULLVALUE {
                /* and null weight */
                ii += 1;
                continue;
            }

            pix = (unsafe { *colptr[0].add(ii) } - histData.amin1) / histData.binsize1;
            ipix = (pix + 1.) as c_long; /* add 1 because the 1st pixel is the null value */

            /* test if bin is within range */
            if ipix < 1 || ipix > histData.haxis1 || pix > histData.maxbin1 {
                ii += 1;
                continue;
            }

            if histData.haxis > 1 {
                if unsafe { *colptr[1].add(ii) } == DOUBLENULLVALUE {
                    ii += 1;
                    continue;
                }

                axisbin = (unsafe { *colptr[1].add(ii) } - histData.amin2) / histData.binsize2;
                iaxisbin = axisbin as c_long;

                if axisbin < 0. || iaxisbin >= histData.haxis2 || axisbin > histData.maxbin2 {
                    ii += 1;
                    continue;
                }

                ipix += iaxisbin * histData.incr[1];

                if histData.haxis > 2 {
                    if unsafe { *colptr[2].add(ii) } == DOUBLENULLVALUE {
                        ii += 1;
                        continue;
                    }

                    axisbin = (unsafe { *colptr[2].add(ii) } - histData.amin3) / histData.binsize3;
                    iaxisbin = axisbin as c_long;
                    if axisbin < 0. || iaxisbin >= histData.haxis3 || axisbin > histData.maxbin3 {
                        ii += 1;
                        continue;
                    }

                    ipix += iaxisbin * histData.incr[2];

                    if histData.haxis > 3 {
                        if unsafe { *colptr[3].add(ii) } == DOUBLENULLVALUE {
                            ii += 1;
                            continue;
                        }

                        axisbin =
                            (unsafe { *colptr[3].add(ii) } - histData.amin4) / histData.binsize4;
                        iaxisbin = axisbin as c_long;
                        if axisbin < 0. || iaxisbin >= histData.haxis4 || axisbin > histData.maxbin4
                        {
                            ii += 1;
                            continue;
                        }

                        ipix += iaxisbin * histData.incr[3];
                    } /* end of haxis > 3 case */
                } /* end of haxis > 2 case */
            } /* end of haxis > 1 case */

            /* increment the histogram pixel */
            unsafe {
                if histData.weight != DOUBLENULLVALUE {
                    /* constant weight factor */

                    /* Note that if wtrecip == 1, the reciprocal was precomputed above */
                    if histData.himagetype == TINT {
                        *histData.hist.j.add(ipix as usize) += histData.weight as c_int;
                    } else if histData.himagetype == TSHORT {
                        *histData.hist.i.add(ipix as usize) += histData.weight as c_short;
                    } else if histData.himagetype == TFLOAT {
                        *histData.hist.r.add(ipix as usize) += histData.weight as f32;
                    } else if histData.himagetype == TDOUBLE {
                        *histData.hist.d.add(ipix as usize) += histData.weight;
                    } else if histData.himagetype == TBYTE {
                        *histData.hist.b.add(ipix as usize) += histData.weight as c_char;
                    }
                } else if histData.wtrecip != 0 {
                    /* use reciprocal of the weight */
                    if histData.himagetype == TINT {
                        *histData.hist.j.add(ipix as usize) += (1. / *colptr[4].add(ii)) as c_int;
                    } else if histData.himagetype == TSHORT {
                        *histData.hist.i.add(ipix as usize) += (1. / *colptr[4].add(ii)) as c_short;
                    } else if histData.himagetype == TFLOAT {
                        *histData.hist.r.add(ipix as usize) += (1. / *colptr[4].add(ii)) as f32;
                    } else if histData.himagetype == TDOUBLE {
                        *histData.hist.d.add(ipix as usize) += 1. / *colptr[4].add(ii);
                    } else if histData.himagetype == TBYTE {
                        *histData.hist.b.add(ipix as usize) += (1. / *colptr[4].add(ii)) as c_char;
                    }
                } else {
                    /* no weights */
                    if histData.himagetype == TINT {
                        *histData.hist.j.add(ipix as usize) += *colptr[4].add(ii) as c_int;
                    } else if histData.himagetype == TSHORT {
                        *histData.hist.i.add(ipix as usize) += *colptr[4].add(ii) as c_short;
                    } else if histData.himagetype == TFLOAT {
                        *histData.hist.r.add(ipix as usize) += *colptr[4].add(ii) as f32;
                    } else if histData.himagetype == TDOUBLE {
                        *histData.hist.d.add(ipix as usize) += *colptr[4].add(ii);
                    } else if histData.himagetype == TBYTE {
                        *histData.hist.b.add(ipix as usize) += *colptr[4].add(ii) as c_char;
                    }
                }
            }
        } /* end of loop over elements per row */
    } /* end of main loop over all rows */

    histData.rowselector_cur = rowselect as *mut c_char; /* Save row pointer for next go-round */
    status
}

#[cfg(test)]
mod tests {
    use super::*;

    // Helper function to compare C strings
    fn c_str_eq(c_arr: &[c_char], expected: &str) -> bool {
        if c_arr.is_empty() || c_arr[0] == 0 {
            return expected.is_empty();
        }

        // Find null terminator
        let len = c_arr.iter().position(|&c| c == 0).unwrap_or(c_arr.len());
        let slice = &c_arr[..len];

        // Convert to bytes for comparison
        let bytes: &[u8] = unsafe { core::slice::from_raw_parts(slice.as_ptr().cast::<u8>(), len) };

        expected.as_bytes() == bytes
    }

    /*--------------------------------------------------------------------------*/
    /* Test ffbins: Parse basic binning specification strings */
    /*--------------------------------------------------------------------------*/

    #[test]
    fn test_ffbins_simple_binsize() {
        let mut status = 0;
        let mut imagetype = 0;
        let mut histaxis = 0;
        let mut recip = 0;
        let mut colname = [[0; FLEN_VALUE]; 4];
        let mut minin = [0.0; 4];
        let mut maxin = [0.0; 4];
        let mut binsizein = [0.0; 4];
        let mut minname = [[0; FLEN_VALUE]; 4];
        let mut maxname = [[0; FLEN_VALUE]; 4];
        let mut binname = [[0; FLEN_VALUE]; 4];
        let mut wt = 0.0;
        let mut wtname = [0; FLEN_VALUE];

        let binspec = cs!(c"bin 10");

        ffbins_safe(
            binspec,
            &mut imagetype,
            &mut histaxis,
            &mut colname,
            &mut minin,
            &mut maxin,
            &mut binsizein,
            &mut minname,
            &mut maxname,
            &mut binname,
            &mut wt,
            &mut wtname,
            &mut recip,
            &mut status,
        );

        assert_eq!(status, 0);
        assert_eq!(histaxis, 2);
        assert_eq!(imagetype, TINT);
        assert_eq!(binsizein[0], 10.0);
        assert_eq!(binsizein[1], 10.0);
        assert_eq!(wt, 1.0);
        assert_eq!(recip, 0);
    }

    #[test]
    fn test_ffbins_single_column() {
        let mut status = 0;
        let mut imagetype = 0;
        let mut histaxis = 0;
        let mut recip = 0;
        let mut colname = [[0; FLEN_VALUE]; 4];
        let mut minin = [0.0; 4];
        let mut maxin = [0.0; 4];
        let mut binsizein = [0.0; 4];
        let mut minname = [[0; FLEN_VALUE]; 4];
        let mut maxname = [[0; FLEN_VALUE]; 4];
        let mut binname = [[0; FLEN_VALUE]; 4];
        let mut wt = 0.0;
        let mut wtname = [0; FLEN_VALUE];

        let binspec = cs!(c"bin X");

        ffbins_safe(
            binspec,
            &mut imagetype,
            &mut histaxis,
            &mut colname,
            &mut minin,
            &mut maxin,
            &mut binsizein,
            &mut minname,
            &mut maxname,
            &mut binname,
            &mut wt,
            &mut wtname,
            &mut recip,
            &mut status,
        );

        assert_eq!(status, 0);
        assert_eq!(histaxis, 1);
        assert!(c_str_eq(&colname[0], "X"));
    }

    #[test]
    fn test_ffbins_two_columns_with_binsize() {
        let mut status = 0;
        let mut imagetype = 0;
        let mut histaxis = 0;
        let mut recip = 0;
        let mut colname = [[0; FLEN_VALUE]; 4];
        let mut minin = [0.0; 4];
        let mut maxin = [0.0; 4];
        let mut binsizein = [0.0; 4];
        let mut minname = [[0; FLEN_VALUE]; 4];
        let mut maxname = [[0; FLEN_VALUE]; 4];
        let mut binname = [[0; FLEN_VALUE]; 4];
        let mut wt = 0.0;
        let mut wtname = [0; FLEN_VALUE];

        let binspec = cs!(c"bin (X,Y) = 5");

        ffbins_safe(
            binspec,
            &mut imagetype,
            &mut histaxis,
            &mut colname,
            &mut minin,
            &mut maxin,
            &mut binsizein,
            &mut minname,
            &mut maxname,
            &mut binname,
            &mut wt,
            &mut wtname,
            &mut recip,
            &mut status,
        );

        assert_eq!(status, 0);
        assert_eq!(histaxis, 2);
        assert!(c_str_eq(&colname[0], "X"));
        assert!(c_str_eq(&colname[1], "Y"));
        assert_eq!(binsizein[0], 5.0);
        assert_eq!(binsizein[1], 5.0);
    }

    #[test]
    fn test_ffbins_full_specification() {
        let mut status = 0;
        let mut imagetype = 0;
        let mut histaxis = 0;
        let mut recip = 0;
        let mut colname = [[0; FLEN_VALUE]; 4];
        let mut minin = [0.0; 4];
        let mut maxin = [0.0; 4];
        let mut binsizein = [0.0; 4];
        let mut minname = [[0; FLEN_VALUE]; 4];
        let mut maxname = [[0; FLEN_VALUE]; 4];
        let mut binname = [[0; FLEN_VALUE]; 4];
        let mut wt = 0.0;
        let mut wtname = [0; FLEN_VALUE];

        let binspec = cs!(c"bin X=0:100:2, Y=10:200:5");

        ffbins_safe(
            binspec,
            &mut imagetype,
            &mut histaxis,
            &mut colname,
            &mut minin,
            &mut maxin,
            &mut binsizein,
            &mut minname,
            &mut maxname,
            &mut binname,
            &mut wt,
            &mut wtname,
            &mut recip,
            &mut status,
        );

        assert_eq!(status, 0);
        assert_eq!(histaxis, 2);
        assert!(c_str_eq(&colname[0], "X"));
        assert!(c_str_eq(&colname[1], "Y"));
        assert_eq!(minin[0], 0.0);
        assert_eq!(maxin[0], 100.0);
        assert_eq!(binsizein[0], 2.0);
        assert_eq!(minin[1], 10.0);
        assert_eq!(maxin[1], 200.0);
        assert_eq!(binsizein[1], 5.0);
    }

    #[test]
    fn test_ffbins_partial_specification() {
        let mut status = 0;
        let mut imagetype = 0;
        let mut histaxis = 0;
        let mut recip = 0;
        let mut colname = [[0; FLEN_VALUE]; 4];
        let mut minin = [0.0; 4];
        let mut maxin = [0.0; 4];
        let mut binsizein = [0.0; 4];
        let mut minname = [[0; FLEN_VALUE]; 4];
        let mut maxname = [[0; FLEN_VALUE]; 4];
        let mut binname = [[0; FLEN_VALUE]; 4];
        let mut wt = 0.0;
        let mut wtname = [0; FLEN_VALUE];

        let binspec = cs!(c"bin X=:100, Y=::5");

        ffbins_safe(
            binspec,
            &mut imagetype,
            &mut histaxis,
            &mut colname,
            &mut minin,
            &mut maxin,
            &mut binsizein,
            &mut minname,
            &mut maxname,
            &mut binname,
            &mut wt,
            &mut wtname,
            &mut recip,
            &mut status,
        );

        assert_eq!(status, 0);
        assert_eq!(histaxis, 2);
        assert!(c_str_eq(&colname[0], "X"));
        assert!(c_str_eq(&colname[1], "Y"));
        // minin[0] should be DOUBLENULLVALUE (undefined)
        assert_eq!(maxin[0], 100.0);
        // binsizein[0] should be DOUBLENULLVALUE (undefined)
        assert_eq!(binsizein[1], 5.0);
    }

    #[test]
    fn test_ffbins_image_type_short() {
        let mut status = 0;
        let mut imagetype = 0;
        let mut histaxis = 0;
        let mut recip = 0;
        let mut colname = [[0; FLEN_VALUE]; 4];
        let mut minin = [0.0; 4];
        let mut maxin = [0.0; 4];
        let mut binsizein = [0.0; 4];
        let mut minname = [[0; FLEN_VALUE]; 4];
        let mut maxname = [[0; FLEN_VALUE]; 4];
        let mut binname = [[0; FLEN_VALUE]; 4];
        let mut wt = 0.0;
        let mut wtname = [0; FLEN_VALUE];

        let binspec = cs!(c"bini X");

        ffbins_safe(
            binspec,
            &mut imagetype,
            &mut histaxis,
            &mut colname,
            &mut minin,
            &mut maxin,
            &mut binsizein,
            &mut minname,
            &mut maxname,
            &mut binname,
            &mut wt,
            &mut wtname,
            &mut recip,
            &mut status,
        );

        assert_eq!(status, 0);
        assert_eq!(imagetype, TSHORT);
        assert!(c_str_eq(&colname[0], "X"));
    }

    #[test]
    fn test_ffbins_image_type_int() {
        let mut status = 0;
        let mut imagetype = 0;
        let mut histaxis = 0;
        let mut recip = 0;
        let mut colname = [[0; FLEN_VALUE]; 4];
        let mut minin = [0.0; 4];
        let mut maxin = [0.0; 4];
        let mut binsizein = [0.0; 4];
        let mut minname = [[0; FLEN_VALUE]; 4];
        let mut maxname = [[0; FLEN_VALUE]; 4];
        let mut binname = [[0; FLEN_VALUE]; 4];
        let mut wt = 0.0;
        let mut wtname = [0; FLEN_VALUE];

        let binspec = cs!(c"binj X");

        ffbins_safe(
            binspec,
            &mut imagetype,
            &mut histaxis,
            &mut colname,
            &mut minin,
            &mut maxin,
            &mut binsizein,
            &mut minname,
            &mut maxname,
            &mut binname,
            &mut wt,
            &mut wtname,
            &mut recip,
            &mut status,
        );

        assert_eq!(status, 0);
        assert_eq!(imagetype, TINT);
        assert!(c_str_eq(&colname[0], "X"));
    }

    #[test]
    fn test_ffbins_four_dimensions() {
        let mut status = 0;
        let mut imagetype = 0;
        let mut histaxis = 0;
        let mut recip = 0;
        let mut colname = [[0; FLEN_VALUE]; 4];
        let mut minin = [0.0; 4];
        let mut maxin = [0.0; 4];
        let mut binsizein = [0.0; 4];
        let mut minname = [[0; FLEN_VALUE]; 4];
        let mut maxname = [[0; FLEN_VALUE]; 4];
        let mut binname = [[0; FLEN_VALUE]; 4];
        let mut wt = 0.0;
        let mut wtname = [0; FLEN_VALUE];

        let binspec = cs!(c"bin X=0:10:1, Y=0:20:2, Z=0:30:3, T=0:40:4");

        ffbins_safe(
            binspec,
            &mut imagetype,
            &mut histaxis,
            &mut colname,
            &mut minin,
            &mut maxin,
            &mut binsizein,
            &mut minname,
            &mut maxname,
            &mut binname,
            &mut wt,
            &mut wtname,
            &mut recip,
            &mut status,
        );

        assert_eq!(status, 0);
        assert_eq!(histaxis, 4);
        assert!(c_str_eq(&colname[0], "X"));
        assert!(c_str_eq(&colname[1], "Y"));
        assert!(c_str_eq(&colname[2], "Z"));
        assert!(c_str_eq(&colname[3], "T"));
        assert_eq!(binsizein[0], 1.0);
        assert_eq!(binsizein[1], 2.0);
        assert_eq!(binsizein[2], 3.0);
        assert_eq!(binsizein[3], 4.0);
    }

    #[test]
    fn test_ffbins_keyword_names() {
        let mut status = 0;
        let mut imagetype = 0;
        let mut histaxis = 0;
        let mut recip = 0;
        let mut colname = [[0; FLEN_VALUE]; 4];
        let mut minin = [0.0; 4];
        let mut maxin = [0.0; 4];
        let mut binsizein = [0.0; 4];
        let mut minname = [[0; FLEN_VALUE]; 4];
        let mut maxname = [[0; FLEN_VALUE]; 4];
        let mut binname = [[0; FLEN_VALUE]; 4];
        let mut wt = 0.0;
        let mut wtname = [0; FLEN_VALUE];

        let binspec = cs!(c"bin X=XMIN:XMAX:XBIN");

        ffbins_safe(
            binspec,
            &mut imagetype,
            &mut histaxis,
            &mut colname,
            &mut minin,
            &mut maxin,
            &mut binsizein,
            &mut minname,
            &mut maxname,
            &mut binname,
            &mut wt,
            &mut wtname,
            &mut recip,
            &mut status,
        );

        assert_eq!(status, 0);
        assert!(c_str_eq(&colname[0], "X"));
        assert!(c_str_eq(&minname[0], "XMIN"));
        assert!(c_str_eq(&maxname[0], "XMAX"));
        assert!(c_str_eq(&binname[0], "XBIN"));
    }

    #[test]
    fn test_ffbins_three_dimensions() {
        let mut status = 0;
        let mut imagetype = 0;
        let mut histaxis = 0;
        let mut recip = 0;
        let mut colname = [[0; FLEN_VALUE]; 4];
        let mut minin = [0.0; 4];
        let mut maxin = [0.0; 4];
        let mut binsizein = [0.0; 4];
        let mut minname = [[0; FLEN_VALUE]; 4];
        let mut maxname = [[0; FLEN_VALUE]; 4];
        let mut binname = [[0; FLEN_VALUE]; 4];
        let mut wt = 0.0;
        let mut wtname = [0; FLEN_VALUE];

        let binspec = cs!(c"bin X=0:10:1, Y=0:20:2, Z=0:30:3");

        ffbins_safe(
            binspec,
            &mut imagetype,
            &mut histaxis,
            &mut colname,
            &mut minin,
            &mut maxin,
            &mut binsizein,
            &mut minname,
            &mut maxname,
            &mut binname,
            &mut wt,
            &mut wtname,
            &mut recip,
            &mut status,
        );

        assert_eq!(status, 0);
        assert_eq!(histaxis, 3);
        assert!(c_str_eq(&colname[0], "X"));
        assert!(c_str_eq(&colname[1], "Y"));
        assert!(c_str_eq(&colname[2], "Z"));
        assert_eq!(binsizein[0], 1.0);
        assert_eq!(binsizein[1], 2.0);
        assert_eq!(binsizein[2], 3.0);
    }

    #[test]
    fn test_ffbins_with_spaces() {
        let mut status = 0;
        let mut imagetype = 0;
        let mut histaxis = 0;
        let mut recip = 0;
        let mut colname = [[0; FLEN_VALUE]; 4];
        let mut minin = [0.0; 4];
        let mut maxin = [0.0; 4];
        let mut binsizein = [0.0; 4];
        let mut minname = [[0; FLEN_VALUE]; 4];
        let mut maxname = [[0; FLEN_VALUE]; 4];
        let mut binname = [[0; FLEN_VALUE]; 4];
        let mut wt = 0.0;
        let mut wtname = [0; FLEN_VALUE];

        let binspec = cs!(c"bin X = 0:100:2, Y = 10:200:5");

        ffbins_safe(
            binspec,
            &mut imagetype,
            &mut histaxis,
            &mut colname,
            &mut minin,
            &mut maxin,
            &mut binsizein,
            &mut minname,
            &mut maxname,
            &mut binname,
            &mut wt,
            &mut wtname,
            &mut recip,
            &mut status,
        );

        assert_eq!(status, 0);
        assert_eq!(histaxis, 2);
        assert!(c_str_eq(&colname[0], "X"));
        assert!(c_str_eq(&colname[1], "Y"));
        assert_eq!(minin[0], 0.0);
        assert_eq!(maxin[0], 100.0);
        assert_eq!(binsizein[0], 2.0);
    }

    /*--------------------------------------------------------------------------*/
    /* Test ffbinse: Extended binning specification with expressions */
    /*--------------------------------------------------------------------------*/

    #[test]
    fn test_ffbinse_with_expression() {
        let mut status = 0;
        let mut imagetype = 0;
        let mut histaxis = 0;
        let mut recip = 0;
        let mut colname = [[0; FLEN_VALUE]; 4];
        let mut minin = [0.0; 4];
        let mut maxin = [0.0; 4];
        let mut binsizein = [0.0; 4];
        let mut minname = [[0; FLEN_VALUE]; 4];
        let mut maxname = [[0; FLEN_VALUE]; 4];
        let mut binname = [[0; FLEN_VALUE]; 4];
        let mut wt = 0.0;
        let mut wtname = [0; FLEN_VALUE];
        let mut exprs: [Box<[c_char]>; 5] = Default::default();

        let binspec = cs!(c"bin X(X*2)=0:100:2");

        ffbinse(
            binspec,
            &mut imagetype,
            &mut histaxis,
            &mut colname,
            &mut minin,
            &mut maxin,
            &mut binsizein,
            &mut minname,
            &mut maxname,
            &mut binname,
            &mut wt,
            &mut wtname,
            &mut recip,
            Some(&mut exprs),
            &mut status,
        );

        assert_eq!(status, 0);
        assert_eq!(histaxis, 1);
        assert!(c_str_eq(&colname[0], "X"));
        // Expression should be "(X*2)"
        assert!(c_str_eq(&exprs[0], "(X*2)"));
    }

    #[test]
    fn test_ffbinse_two_expressions() {
        let mut status = 0;
        let mut imagetype = 0;
        let mut histaxis = 0;
        let mut recip = 0;
        let mut colname = [[0; FLEN_VALUE]; 4];
        let mut minin = [0.0; 4];
        let mut maxin = [0.0; 4];
        let mut binsizein = [0.0; 4];
        let mut minname = [[0; FLEN_VALUE]; 4];
        let mut maxname = [[0; FLEN_VALUE]; 4];
        let mut binname = [[0; FLEN_VALUE]; 4];
        let mut wt = 0.0;
        let mut wtname = [0; FLEN_VALUE];
        let mut exprs: [Box<[c_char]>; 5] = Default::default();

        let binspec = cs!(c"bin X(X+Y)=0:10:1, Y(Y*2)=0:20:2");

        ffbinse(
            binspec,
            &mut imagetype,
            &mut histaxis,
            &mut colname,
            &mut minin,
            &mut maxin,
            &mut binsizein,
            &mut minname,
            &mut maxname,
            &mut binname,
            &mut wt,
            &mut wtname,
            &mut recip,
            Some(&mut exprs),
            &mut status,
        );

        assert_eq!(status, 0);
        assert_eq!(histaxis, 2);
        // Expressions should be "(X+Y)" and "(Y*2)"
        assert!(c_str_eq(&exprs[0], "(X+Y)"));
        assert!(c_str_eq(&exprs[1], "(Y*2)"));
    }

    #[test]
    fn test_ffbinse_mixed_column_and_expression() {
        let mut status = 0;
        let mut imagetype = 0;
        let mut histaxis = 0;
        let mut recip = 0;
        let mut colname = [[0; FLEN_VALUE]; 4];
        let mut minin = [0.0; 4];
        let mut maxin = [0.0; 4];
        let mut binsizein = [0.0; 4];
        let mut minname = [[0; FLEN_VALUE]; 4];
        let mut maxname = [[0; FLEN_VALUE]; 4];
        let mut binname = [[0; FLEN_VALUE]; 4];
        let mut wt = 0.0;
        let mut wtname = [0; FLEN_VALUE];
        let mut exprs: [Box<[c_char]>; 5] = Default::default();

        let binspec = cs!(c"bin X=0:10:1, Y(sqrt(Y))=0:20:2");

        ffbinse(
            binspec,
            &mut imagetype,
            &mut histaxis,
            &mut colname,
            &mut minin,
            &mut maxin,
            &mut binsizein,
            &mut minname,
            &mut maxname,
            &mut binname,
            &mut wt,
            &mut wtname,
            &mut recip,
            Some(&mut exprs),
            &mut status,
        );

        assert_eq!(status, 0);
        assert_eq!(histaxis, 2);
        assert!(c_str_eq(&colname[0], "X"));
        assert!(c_str_eq(&colname[1], "Y"));
        // exprs[0] should be empty (no expression for X)
        assert!(c_str_eq(&exprs[0], ""));
        // exprs[1] should contain "(sqrt(Y))"
        assert!(c_str_eq(&exprs[1], "(sqrt(Y))"));
    }

    /*--------------------------------------------------------------------------*/
    /* Test ffbinr: Parse binning range specification */
    /*--------------------------------------------------------------------------*/

    #[test]
    fn test_ffbinr_column_only() {
        let mut status = 0;
        let mut colname = [0; FLEN_VALUE];
        let mut minin = 0.0;
        let mut maxin = 0.0;
        let mut binsizein = 0.0;
        let mut minname = [0; FLEN_VALUE];
        let mut maxname = [0; FLEN_VALUE];
        let mut binname = [0; FLEN_VALUE];

        let spec = cs!(c"XCOL");
        let mut ptr = spec;

        ffbinr_safe(
            &mut ptr,
            &mut colname,
            &mut minin,
            &mut maxin,
            &mut binsizein,
            &mut minname,
            &mut maxname,
            &mut binname,
            &mut status,
        );

        assert_eq!(status, 0);
        assert!(c_str_eq(&colname, "XCOL"));
    }

    #[test]
    fn test_ffbinr_column_with_binsize() {
        let mut status = 0;
        let mut colname = [0; FLEN_VALUE];
        let mut minin = 0.0;
        let mut maxin = 0.0;
        let mut binsizein = 0.0;
        let mut minname = [0; FLEN_VALUE];
        let mut maxname = [0; FLEN_VALUE];
        let mut binname = [0; FLEN_VALUE];

        let spec = cs!(c"XCOL=5");
        let mut ptr = spec;

        ffbinr_safe(
            &mut ptr,
            &mut colname,
            &mut minin,
            &mut maxin,
            &mut binsizein,
            &mut minname,
            &mut maxname,
            &mut binname,
            &mut status,
        );

        assert_eq!(status, 0);
        assert!(c_str_eq(&colname, "XCOL"));
        assert_eq!(binsizein, 5.0);
    }

    #[test]
    fn test_ffbinr_full_range() {
        let mut status = 0;
        let mut colname = [0; FLEN_VALUE];
        let mut minin = 0.0;
        let mut maxin = 0.0;
        let mut binsizein = 0.0;
        let mut minname = [0; FLEN_VALUE];
        let mut maxname = [0; FLEN_VALUE];
        let mut binname = [0; FLEN_VALUE];

        let spec = cs!(c"XCOL=10:100:2");
        let mut ptr = spec;

        ffbinr_safe(
            &mut ptr,
            &mut colname,
            &mut minin,
            &mut maxin,
            &mut binsizein,
            &mut minname,
            &mut maxname,
            &mut binname,
            &mut status,
        );

        assert_eq!(status, 0);
        assert!(c_str_eq(&colname, "XCOL"));
        assert_eq!(minin, 10.0);
        assert_eq!(maxin, 100.0);
        assert_eq!(binsizein, 2.0);
    }

    #[test]
    fn test_ffbinr_partial_range_max_only() {
        let mut status = 0;
        let mut colname = [0; FLEN_VALUE];
        let mut minin = 0.0;
        let mut maxin = 0.0;
        let mut binsizein = 0.0;
        let mut minname = [0; FLEN_VALUE];
        let mut maxname = [0; FLEN_VALUE];
        let mut binname = [0; FLEN_VALUE];

        let spec = cs!(c"XCOL=:100");
        let mut ptr = spec;

        ffbinr_safe(
            &mut ptr,
            &mut colname,
            &mut minin,
            &mut maxin,
            &mut binsizein,
            &mut minname,
            &mut maxname,
            &mut binname,
            &mut status,
        );

        assert_eq!(status, 0);
        assert!(c_str_eq(&colname, "XCOL"));
        assert_eq!(maxin, 100.0);
    }

    #[test]
    fn test_ffbinr_partial_range_binsize_only() {
        let mut status = 0;
        let mut colname = [0; FLEN_VALUE];
        let mut minin = 0.0;
        let mut maxin = 0.0;
        let mut binsizein = 0.0;
        let mut minname = [0; FLEN_VALUE];
        let mut maxname = [0; FLEN_VALUE];
        let mut binname = [0; FLEN_VALUE];

        let spec = cs!(c"XCOL=::5");
        let mut ptr = spec;

        ffbinr_safe(
            &mut ptr,
            &mut colname,
            &mut minin,
            &mut maxin,
            &mut binsizein,
            &mut minname,
            &mut maxname,
            &mut binname,
            &mut status,
        );

        assert_eq!(status, 0);
        assert!(c_str_eq(&colname, "XCOL"));
        assert_eq!(binsizein, 5.0);
    }

    #[test]
    fn test_ffbinr_column_number() {
        let mut status = 0;
        let mut colname = [0; FLEN_VALUE];
        let mut minin = 0.0;
        let mut maxin = 0.0;
        let mut binsizein = 0.0;
        let mut minname = [0; FLEN_VALUE];
        let mut maxname = [0; FLEN_VALUE];
        let mut binname = [0; FLEN_VALUE];

        let spec = cs!(c"#5=0:100:2");
        let mut ptr = spec;

        ffbinr_safe(
            &mut ptr,
            &mut colname,
            &mut minin,
            &mut maxin,
            &mut binsizein,
            &mut minname,
            &mut maxname,
            &mut binname,
            &mut status,
        );

        assert_eq!(status, 0);
        // The '#' should be stripped, leaving just "5"
        assert!(c_str_eq(&colname, "5"));
        assert_eq!(minin, 0.0);
        assert_eq!(maxin, 100.0);
        assert_eq!(binsizein, 2.0);
    }

    #[test]
    fn test_ffbinr_with_keyword_names() {
        let mut status = 0;
        let mut colname = [0; FLEN_VALUE];
        let mut minin = 0.0;
        let mut maxin = 0.0;
        let mut binsizein = 0.0;
        let mut minname = [0; FLEN_VALUE];
        let mut maxname = [0; FLEN_VALUE];
        let mut binname = [0; FLEN_VALUE];

        let spec = cs!(c"XCOL=XMIN:XMAX:XBIN");
        let mut ptr = spec;

        ffbinr_safe(
            &mut ptr,
            &mut colname,
            &mut minin,
            &mut maxin,
            &mut binsizein,
            &mut minname,
            &mut maxname,
            &mut binname,
            &mut status,
        );

        assert_eq!(status, 0);
        assert!(c_str_eq(&colname, "XCOL"));
        assert!(c_str_eq(&minname, "XMIN"));
        assert!(c_str_eq(&maxname, "XMAX"));
        assert!(c_str_eq(&binname, "XBIN"));
    }

    #[test]
    fn test_ffbinr_min_and_max_only() {
        let mut status = 0;
        let mut colname = [0; FLEN_VALUE];
        let mut minin = 0.0;
        let mut maxin = 0.0;
        let mut binsizein = 0.0;
        let mut minname = [0; FLEN_VALUE];
        let mut maxname = [0; FLEN_VALUE];
        let mut binname = [0; FLEN_VALUE];

        let spec = cs!(c"XCOL=10:100");
        let mut ptr = spec;

        ffbinr_safe(
            &mut ptr,
            &mut colname,
            &mut minin,
            &mut maxin,
            &mut binsizein,
            &mut minname,
            &mut maxname,
            &mut binname,
            &mut status,
        );

        assert_eq!(status, 0);
        assert!(c_str_eq(&colname, "XCOL"));
        assert_eq!(minin, 10.0);
        assert_eq!(maxin, 100.0);
    }

    /*--------------------------------------------------------------------------*/
    /* Test ffbinre: Extended range parser with expressions */
    /*--------------------------------------------------------------------------*/

    #[test]
    fn test_ffbinre_with_expression() {
        let mut status = 0;
        let mut colname = [0; FLEN_VALUE];
        let mut minin = 0.0;
        let mut maxin = 0.0;
        let mut binsizein = 0.0;
        let mut minname = [0; FLEN_VALUE];
        let mut maxname = [0; FLEN_VALUE];
        let mut binname = [0; FLEN_VALUE];
        let mut exprbeg: usize = 0;
        let mut exprend: usize = 0;

        let spec = cs!(c"XCOL(X*2)=0:100:2");
        let mut ptr = spec;

        ffbinre(
            &mut ptr,
            &mut colname,
            Some(&mut exprbeg),
            Some(&mut exprend),
            &mut minin,
            &mut maxin,
            &mut binsizein,
            &mut minname,
            &mut maxname,
            &mut binname,
            &mut status,
        );

        assert_eq!(status, 0);
        assert!(c_str_eq(&colname, "XCOL"));

        // Expression should be "(X*2)" which is 5 chars
        let expr_len = exprend - exprbeg;
        assert_eq!(expr_len, 5);
        assert_eq!(minin, 0.0);
        assert_eq!(maxin, 100.0);
        assert_eq!(binsizein, 2.0);
    }

    #[test]
    fn test_ffbinre_without_expression() {
        let mut status = 0;
        let mut colname = [0; FLEN_VALUE];
        let mut minin = 0.0;
        let mut maxin = 0.0;
        let mut binsizein = 0.0;
        let mut minname = [0; FLEN_VALUE];
        let mut maxname = [0; FLEN_VALUE];
        let mut binname = [0; FLEN_VALUE];
        let mut exprbeg: usize = 999; // Set to non-zero to test it gets cleared
        let mut exprend: usize = 999;

        let spec = cs!(c"XCOL=0:100:2");
        let mut ptr = spec;

        ffbinre(
            &mut ptr,
            &mut colname,
            Some(&mut exprbeg),
            Some(&mut exprend),
            &mut minin,
            &mut maxin,
            &mut binsizein,
            &mut minname,
            &mut maxname,
            &mut binname,
            &mut status,
        );

        assert_eq!(status, 0);
        assert!(c_str_eq(&colname, "XCOL"));
        // No expression, so both should be 0
        assert_eq!(exprbeg, 0);
        assert_eq!(exprend, 0);
        assert_eq!(minin, 0.0);
        assert_eq!(maxin, 100.0);
        assert_eq!(binsizein, 2.0);
    }
}
