/*   Globally defined histogram parameters */

use core::slice;
use std::ffi::CStr;
use std::ffi::c_void;
use std::ptr;

use crate::c_types::{c_char, c_int, c_long, c_short};

use bytemuck::{cast_slice, cast_slice_mut};
use cstr::cstr;

use crate::aliases::ffgnrw_safe;
use crate::aliases::safer::{fits_read_key_str, fits_write_key_str};

use crate::cs;
use crate::eval_defs::{ParseData, parseInfo};
use crate::fitscore::ffkeyn_safe;
use crate::fitsio::*;
use crate::getcol::ffgcv_safe;
use crate::getkey::ffgky_safe;
use crate::putkey::ffpky_safe;
use crate::{KeywordDatatype, NullValue, raw_to_slice};

/*  Structure holding all the histogramming information   */
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
    iterCols: &'a mut iteratorCol,
    parsers: &'a mut ParseData<'a>,
    infos: &'a mut parseInfo<'a>,
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
    binspec: &[c_char],                             /* I - binning specification */
    imagetype: &mut c_int,                          /* O - image type, TINT or TSHORT */
    histaxis: &mut c_int,                           /* O - no. of axes in the histogram */
    colname: &[[c_char; FLEN_VALUE]; 4],            /* column name for axis */
    minin: &[f64; 4],                               /* minimum value for each axis */
    maxin: &[f64; 4],                               /* maximum value for each axis */
    binsizein: &[f64; 4],                           /* size of bins on each axis */
    minname: &[[c_char; FLEN_VALUE]; 4],            /* keyword name for min */
    maxname: &[[c_char; FLEN_VALUE]; 4],            /* keyword name for max */
    binname: &[[c_char; FLEN_VALUE]; 4],            /* keyword name for binsize */
    weight: &mut f64,                               /* weighting factor          */
    wtcol: &[c_char; 71],                           /* keyword or column name for weight */
    recip: &mut c_int,                              /* the reciprocal of the weight? */
    exprs: Option<&mut Option<[Box<[c_char]>; 5]>>, /* returned with expressions (or 0) */
    status: &mut c_int,
) -> c_int {
    todo!()
}

/*--------------------------------------------------------------------------*/
/// Parse non-extended expression, but otherwise the same as ffbinse()
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffbins(
    binspec: *const c_char,                    /* I - binning specification */
    imagetype: *mut c_int,                     /* O - image type, TINT or TSHORT */
    histaxis: *mut c_int,                      /* O - no. of axes in the histogram */
    colname: *const [[c_char; FLEN_VALUE]; 4], /* column name for axis */
    minin: *const f64,                         /* minimum value for each axis */
    maxin: *const f64,                         /* maximum value for each axis */
    binsizein: *const f64,                     /* size of bins on each axis */
    minname: *mut [[c_char; FLEN_VALUE]; 4],   /* keyword name for min */
    maxname: *mut [[c_char; FLEN_VALUE]; 4],   /* keyword name for max */
    binname: *mut [[c_char; FLEN_VALUE]; 4],   /* keyword name for binsize */
    wt: *mut f64,                              /* weighting factor          */
    wtname: *mut [c_char; FLEN_VALUE],         /* keyword or column name for weight */
    recip: *mut c_int,                         /* the reciprocal of the weight? */
    status: *mut c_int,
) -> c_int {
    unsafe {
        raw_to_slice!(binspec);

        let imagetype = imagetype.as_mut().expect(NULL_MSG);
        let histaxis = histaxis.as_mut().expect(NULL_MSG);
        let colname = colname.as_ref().expect(NULL_MSG);
        let minin: &[f64; 4] = slice::from_raw_parts(minin, 4).try_into().unwrap();
        let maxin: &[f64; 4] = slice::from_raw_parts(maxin, 4).try_into().unwrap();
        let binsizein: &[f64; 4] = slice::from_raw_parts(binsizein, 4).try_into().unwrap();
        let minname = minname.as_mut().expect(NULL_MSG);
        let maxname = maxname.as_mut().expect(NULL_MSG);
        let binname = binname.as_mut().expect(NULL_MSG);
        let weight = wt.as_mut().expect(NULL_MSG);
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
pub(crate) fn ffbins_safe(
    binspec: &[c_char],                  /* I - binning specification */
    imagetype: &mut c_int,               /* O - image type, TINT or TSHORT */
    histaxis: &mut c_int,                /* O - no. of axes in the histogram */
    colname: &[[c_char; FLEN_VALUE]; 4], /* column name for axis */
    minin: &[f64; 4],                    /* minimum value for each axis */
    maxin: &[f64; 4],                    /* maximum value for each axis */
    binsizein: &[f64; 4],                /* size of bins on each axis */
    minname: &[[c_char; FLEN_VALUE]; 4], /* keyword name for min */
    maxname: &[[c_char; FLEN_VALUE]; 4], /* keyword name for max */
    binname: &[[c_char; FLEN_VALUE]; 4], /* keyword name for binsize */
    weight: &mut f64,                    /* weighting factor          */
    wtcol: &[c_char; FLEN_VALUE],        /* keyword or column name for weight */
    recip: &mut c_int,                   /* the reciprocal of the weight? */
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
    ptr: *mut *mut c_char,
    colname: *mut c_char,
    exprbeg: *mut *mut c_char,
    exprend: *mut *mut c_char,
    minin: *mut f64,
    maxin: *mut f64,
    binsizein: *mut f64,
    minname: *mut c_char,
    maxname: *mut c_char,
    binname: *mut c_char,
    status: *mut c_int,
) -> c_int {
    todo!();
}

/*--------------------------------------------------------------------------*/
/// Parse the input binning range specification string, returning
/// the column name, histogram min and max values, and bin size.
///
/// This is the non-extended version of the parser which disallows
/// binning expressions.  Only column names are allowed.
pub(crate) fn ffbinr(
    ptr: *mut *mut c_char,
    colname: *mut c_char,
    minin: *mut f64,
    maxin: *mut f64,
    binsizein: *mut f64,
    minname: *mut c_char,
    maxname: *mut c_char,
    binname: *mut c_char,
    status: *mut c_int,
) -> c_int {
    ffbinre(
        ptr,
        colname,
        ptr::null_mut(),
        ptr::null_mut(),
        minin,
        maxin,
        binsizein,
        minname,
        maxname,
        binname,
        status,
    )
}

/*--------------------------------------------------------------------------*/
pub(crate) fn ffhist2e(
    fptr: &mut Option<Box<fitsfile>>, /* IO - pointer to table with X and Y cols, on output, points to histogram image    */
    outfile: &[c_char],               /* I - name for the output histogram file      */
    imagetype: c_int,                 /* I - datatype for image: TINT, TSHORT, etc   */
    naxis: c_int,                     /* I - number of axes in the histogram image   */
    colname: &[[c_char; FLEN_VALUE]; 4], /* I - column names               */
    colexpr: Option<&[&[c_char]; 4]>, /* I - optionally, expression intead of colum  */
    minin: &[f64],                    /* I - minimum histogram value, for each axis */
    maxin: &[f64],                    /* I - maximum histogram value, for each axis */
    binsizein: &[f64],                /* I - bin size along each axis               */
    minname: Option<&[[c_char; FLEN_VALUE]; 4]>, /* I - optional keywords for min    */
    maxname: Option<&[[c_char; FLEN_VALUE]; 4]>, /* I - optional keywords for max    */
    binname: Option<&[[c_char; FLEN_VALUE]; 4]>, /* I - optional keywords for binsize */
    weightin: f64,                    /* I - binning weighting factor          */
    wtcol: Option<&[c_char; FLEN_VALUE]>, /* I - optional keyword or col for weight*/
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
    todo!()
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
        let minname = minname.as_ref();
        let maxname = maxname.as_ref();
        let binname = binname.as_ref();
        let wtcol = wtcol.as_ref();
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
pub(crate) fn ffhist2_safe(
    fptr: &mut Option<Box<fitsfile>>, /* IO - pointer to table with X and Y cols;    */
    /*     on output, points to histogram image    */
    outfile: &[c_char], /* I - name for the output histogram file      */
    imagetype: c_int,   /* I - datatype for image: TINT, TSHORT, etc   */
    naxis: c_int,       /* I - number of axes in the histogram image   */
    colname: &[[c_char; FLEN_VALUE]; 4], /* I - column names               */
    minin: &[f64],      /* I - minimum histogram value, for each axis */
    maxin: &[f64],      /* I - maximum histogram value, for each axis */
    binsizein: &[f64],  /* I - bin size along each axis               */
    minname: Option<&[[c_char; FLEN_VALUE]; 4]>, /* I - optional keywords for min    */
    maxname: Option<&[[c_char; FLEN_VALUE]; 4]>, /* I - optional keywords for max    */
    binname: Option<&[[c_char; FLEN_VALUE]; 4]>, /* I - optional keywords for binsize */
    weightin: f64,      /* I - binning weighting factor          */
    wtcol: Option<&[c_char; FLEN_VALUE]>, /* I - optional keyword or col for weight*/
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
) -> c_int {
    unsafe {
        let fptr = &mut ((*fptr).as_mut()).expect(NULL_MSG);
        raw_to_slice!(outfile);
        let colname = colname.as_ref().unwrap();
        let minin = slice::from_raw_parts(minin, naxis as usize);
        let maxin = slice::from_raw_parts(maxin, naxis as usize);
        let binsizein = slice::from_raw_parts(binsizein, naxis as usize);
        let minname = minname.as_ref();
        let maxname = maxname.as_ref();
        let binname = binname.as_ref();
        let wtcol = wtcol.as_ref();
        let status = status.as_mut().expect(NULL_MSG);

        let mut nrows = 0;
        ffgnrw_safe(fptr, &mut nrows, status); /* no. of rows */

        let selectrow = if selectrow.is_null() {
            None
        } else {
            Some(slice::from_raw_parts(selectrow, nrows as usize))
        };

        ffhist3_safe(
            fptr, outfile, imagetype, naxis, colname, minin, maxin, binsizein, minname, maxname,
            binname, weightin, wtcol, recip, selectrow, status,
        )
    }
}

/*--------------------------------------------------------------------------*/
/// ffhist3: same as ffhist2, but does not close the original file
///  and/or replace the original file pointer
pub(crate) fn ffhist3_safe(
    fptr: &&mut fitsfile, /* IO - pointer to table with X and Y cols;    */
    outfile: &[c_char],   /* I - name for the output histogram file      */
    imagetype: c_int,     /* I - datatype for image: TINT, TSHORT, etc   */
    naxis: c_int,         /* I - number of axes in the histogram image   */
    colname: &[[c_char; FLEN_VALUE]; 4], /* I - column names               */
    minin: &[f64],        /* I - minimum histogram value, for each axis */
    maxin: &[f64],        /* I - maximum histogram value, for each axis */
    binsizein: &[f64],    /* I - bin size along each axis               */
    minname: Option<&[[c_char; FLEN_VALUE]; 4]>, /* I - optional keywords for min    */
    maxname: Option<&[[c_char; FLEN_VALUE]; 4]>, /* I - optional keywords for max    */
    binname: Option<&[[c_char; FLEN_VALUE]; 4]>, /* I - optional keywords for binsize */
    weightin: f64,        /* I - binning weighting factor          */
    wtcol: Option<&[c_char; FLEN_VALUE]>, /* I - optional keyword or col for weight*/
    recip: c_int,         /* I - use reciprocal of the weight?     */
    selectrow: Option<&[c_char]>, /* I - optional array (length = no. of   */
    /* rows in the table).  If the element is true */
    /* then the corresponding row of the table will*/
    /* be included in the histogram, otherwise the */
    /* row will be skipped.  Ingnored if *selectrow*/
    /* is equal to NULL.                           */
    status: &mut c_int,
) -> c_int {
    todo!();
}

/*--------------------------------------------------------------------------*/
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffhist(
    fptr: *mut *mut fitsfile, /* I - pointer to table with X and Y cols; on output, points to histogram image    */
    outfile: *const c_char,   /* I - name for the output histogram file      */
    imagetype: c_int,         /* I - datatype for image: TINT, TSHORT, etc   */
    naxis: c_int,             /* I - number of axes in the histogram image   */
    colname: *const [[c_char; FLEN_VALUE]; 4], /* I - column names               */
    minin: *const f64,        /* I - minimum histogram value, for each axis */
    maxin: *const f64,        /* I - maximum histogram value, for each axis */
    binsizein: *const f64,    /* I - bin size along each axis               */
    minname: *const [[c_char; FLEN_VALUE]; 4], /* I - optional keywords for min    */
    maxname: *const [[c_char; FLEN_VALUE]; 4], /* I - optional keywords for max    */
    binname: *const [[c_char; FLEN_VALUE]; 4], /* I - optional keywords for binsize */
    weightin: f64,            /* I - binning weighting factor          */
    wtcol: *const [c_char; FLEN_VALUE], /* I - optional keyword or col for weight*/
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
        let fptr = &mut (*fptr).as_mut().expect(NULL_MSG);
        raw_to_slice!(outfile);
        let colname = colname.as_ref().expect(NULL_MSG);
        let minin = slice::from_raw_parts(minin, naxis as usize);
        let maxin = slice::from_raw_parts(maxin, naxis as usize);
        let binsizein = slice::from_raw_parts(binsizein, naxis as usize);
        let minname = minname.as_ref();
        let maxname = maxname.as_ref();
        let binname = binname.as_ref();
        let wtcol = wtcol.as_ref();
        let status = status.as_mut().expect(NULL_MSG);

        let mut nrows = 0;
        ffgnrw_safe(fptr, &mut nrows, status); /* no. of rows */

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
pub(crate) fn ffhist_safe(
    fptr: &mut &mut fitsfile, /* I - pointer to table with X and Y cols; on output, points to histogram image    */
    outfile: &[c_char],       /* I - name for the output histogram file      */
    imagetype: c_int,         /* I - datatype for image: TINT, TSHORT, etc   */
    naxis: c_int,             /* I - number of axes in the histogram image   */
    colname: &[[c_char; FLEN_VALUE]; 4], /* I - column names               */
    minin: &[f64],            /* I - minimum histogram value, for each axis */
    maxin: &[f64],            /* I - maximum histogram value, for each axis */
    binsizein: &[f64],        /* I - bin size along each axis               */
    minname: Option<&[[c_char; FLEN_VALUE]; 4]>, /* I - optional keywords for min    */
    maxname: Option<&[[c_char; FLEN_VALUE]; 4]>, /* I - optional keywords for max    */
    binname: Option<&[[c_char; FLEN_VALUE]; 4]>, /* I - optional keywords for binsize */
    weightin: f64,            /* I - binning weighting factor          */
    wtcol: Option<&[c_char; FLEN_VALUE]>, /* I - optional keyword or col for weight*/
    recip: c_int,             /* I - use reciprocal of the weight?     */
    selectrow: Option<&[c_char]>, /* I - optional array (length = no. of   */
    /* rows in the table).  If the element is true */
    /* then the corresponding row of the table will*/
    /* be included in the histogram, otherwise the */
    /* row will be skipped.  Ingnored if *selectrow*/
    /* is equal to NULL.                           */
    status: &mut c_int,
) -> c_int {
    todo!();
}

/*--------------------------------------------------------------------------*/
/// Single-precision version
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn fits_calc_binning(
    fptr: *mut fitsfile, /* IO - pointer to table to be binned      ;       */
    naxis: c_int,        /* I - number of axes/columns in the binned image  */
    colname: *const [[c_char; FLEN_VALUE]; 4], /* I - optional column names         */
    minin: *const f64,   /* I - optional lower bound value for each axis  */
    maxin: *const f64,   /* I - optional upper bound value, for each axis */
    binsizein: *const f64, /* I - optional bin size along each axis         */
    minname: *const [[c_char; FLEN_VALUE]; 4], /* I - optional keywords for min       */
    maxname: *const [[c_char; FLEN_VALUE]; 4], /* I - optional keywords for max       */
    binname: *const [[c_char; FLEN_VALUE]; 4], /* I - optional keywords for binsize   */

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
        let colname = colname.as_ref();
        let minin = if minin.is_null() {
            None
        } else {
            Some(slice::from_raw_parts(minin, naxis as usize))
        };
        let maxin = if maxin.is_null() {
            None
        } else {
            Some(slice::from_raw_parts(maxin, naxis as usize))
        };
        let binsizein = if binsizein.is_null() {
            None
        } else {
            Some(slice::from_raw_parts(binsizein, naxis as usize))
        };
        let minname = minname.as_ref();
        let maxname = maxname.as_ref();
        let binname = binname.as_ref();
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
pub(crate) fn fits_calc_binning_safe(
    fptr: &mut fitsfile, /* IO - pointer to table to be binned      ;       */
    naxis: c_int,        /* I - number of axes/columns in the binned image  */
    colname: Option<&[[c_char; FLEN_VALUE]; 4]>, /* I - optional column names         */
    minin: Option<&[f64]>, /* I - optional lower bound value for each axis  */
    maxin: Option<&[f64]>, /* I - optional upper bound value, for each axis */
    binsizein: Option<&[f64]>, /* I - optional bin size along each axis         */
    minname: Option<&[[c_char; FLEN_VALUE]; 4]>, /* I - optional keywords for min       */
    maxname: Option<&[[c_char; FLEN_VALUE]; 4]>, /* I - optional keywords for max       */
    binname: Option<&[[c_char; FLEN_VALUE]; 4]>, /* I - optional keywords for binsize   */

    /* The returned parameters for each axis of the n-dimensional histogram are */
    colnum: &mut [c_int], /* O - column numbers, to be binned */
    haxes: &mut [c_long], /* O - number of bins in each histogram axis */
    amin: &mut [f32],     /* O - lower bound of the histogram axes */
    amax: &mut [f32],     /* O - upper bound of the histogram axes */
    binsize: &mut [f32],  /* O - width of histogram bins/pixels on each axis */
    status: &mut c_int,
) -> c_int {
    todo!();
}

/*--------------------------------------------------------------------------*/
/// Double-precision version with extended syntax
/// Calculate the actual binning parameters, based on various user inputoptions.
///
/// Note: caller is responsible to free parsers[*] upon return using ffcprs()
pub(crate) fn fits_calc_binningde(
    fptr: &mut fitsfile, /* IO - pointer to table to be binned      ;       */
    naxis: c_int,        /* I - number of axes/columns in the binned image  */
    colname: Option<&[[c_char; FLEN_VALUE]; 4]>, /* I - optional column names         */
    colexpr: Option<&[&[c_char]; 4]>, /* I - optional column expressions instead of name    */
    minin: Option<&[f64]>, /* I - optional lower bound value for each axis  */
    maxin: Option<&[f64]>, /* I - optional upper bound value, for each axis */
    binsizein: Option<&[f64]>, /* I - optional bin size along each axis         */
    minname: Option<&[[c_char; FLEN_VALUE]; 4]>, /* I - optional keywords for min       */
    maxname: Option<&[[c_char; FLEN_VALUE]; 4]>, /* I - optional keywords for max       */
    binname: Option<&[[c_char; FLEN_VALUE]; 4]>, /* I - optional keywords for binsize   */

    /* The returned parameters for each axis of the n-dimensional histogram are */
    colnum: &mut [c_int],            /* O - column numbers, to be binned */
    datatypes: Option<&mut [c_int]>, /* O - datatype for each column */
    haxes: &mut [c_long],            /* O - number of bins in each histogram axis */
    amin: &mut [f64],                /* O - lower bound of the histogram axes */
    amax: &mut [f64],                /* O - upper bound of the histogram axes */
    binsize: &mut [f64],             /* O - width of histogram bins/pixels on each axis */
    repeat: Option<&mut [c_long]>,   /* O - vector repeat of input columns */
    status: &mut c_int,
) -> c_int {
    todo!();
}

/*--------------------------------------------------------------------------*/
/// Double-precision version, with non-extended syntax
/// Calculate the actual binning parameters, based on various user input options.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn fits_calc_binningd(
    fptr: *mut fitsfile, /* IO - pointer to table to be binned      ;       */
    naxis: c_int,        /* I - number of axes/columns in the binned image  */
    colname: *const [[c_char; FLEN_VALUE]; 4], /* I - optional column names         */
    minin: *const f64,   /* I - optional lower bound value for each axis  */
    maxin: *const f64,   /* I - optional upper bound value, for each axis */
    binsizein: *const f64, /* I - optional bin size along each axis         */
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
        let colname = colname.as_ref();
        let minin = if minin.is_null() {
            None
        } else {
            Some(slice::from_raw_parts(minin, naxis as usize))
        };
        let maxin = if maxin.is_null() {
            None
        } else {
            Some(slice::from_raw_parts(maxin, naxis as usize))
        };
        let binsizein = if binsizein.is_null() {
            None
        } else {
            Some(slice::from_raw_parts(binsizein, naxis as usize))
        };
        let minname = minname.as_ref();
        let maxname = maxname.as_ref();
        let binname = binname.as_ref();
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
pub(crate) fn fits_calc_binningd_safe(
    fptr: &mut fitsfile, /* IO - pointer to table to be binned      ;       */
    naxis: c_int,        /* I - number of axes/columns in the binned image  */
    colname: Option<&[[c_char; FLEN_VALUE]; 4]>, /* I - optional column names         */
    minin: Option<&[f64]>, /* I - optional lower bound value for each axis  */
    maxin: Option<&[f64]>, /* I - optional upper bound value, for each axis */
    binsizein: Option<&[f64]>, /* I - optional bin size along each axis         */
    minname: Option<&[[c_char; FLEN_VALUE]; 4]>, /* I - optional keywords for min       */
    maxname: Option<&[[c_char; FLEN_VALUE]; 4]>, /* I - optional keywords for max       */
    binname: Option<&[[c_char; FLEN_VALUE]; 4]>, /* I - optional keywords for binsize   */

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
    colexpr: Option<&[&[c_char]; 4]>, /* I - if expression, then column name to use */
    status: &mut c_int,
) -> c_int {
    let ii: c_int = 0;
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
            && colexpr[ii][0] != 0
        {
            let colname = colname.as_ref().expect(NULL_MSG);

            // && colexpr[ii] && colexpr[ii][0] && colname[ii]) {
            /* Column expression: we need to put the column name from the binning expression */
            ffkeyn_safe(cs!("CTYPE"), (ii + 1) as c_int, &mut keyname, &mut tstatus);
            fits_write_key_str(
                histptr,
                &keyname,
                &colname[ii],
                Some(cs!("Coordinate Type")),
                &mut tstatus,
            );
        } else {
            /* Column name */
            tstatus = 0;
            ffkeyn_safe(cs!("CTYPE"), (ii + 1) as c_int, &mut keyname, &mut tstatus);
            fits_read_key_str(histptr, &keyname, &mut svalue, None, &mut tstatus);

            if tstatus == 0 {
                continue; /* keyword already exists, so skip to next axis */
            }

            /* use column name as the axis name */
            tstatus = 0;
            ffkeyn_safe(cs!("TTYPE"), colnum[ii], &mut keyname, &mut tstatus);
            fits_read_key_str(fptr, &keyname, &mut svalue, None, &mut tstatus);

            if tstatus == 0 {
                ffkeyn_safe(cs!("CTYPE"), (ii + 1) as c_int, &mut keyname, &mut tstatus);
                fits_write_key_str(
                    histptr,
                    &keyname,
                    &svalue,
                    Some(cs!("Coordinate Type")),
                    &mut tstatus,
                );
            }

            /*  CUNITn,  use the column units */
            tstatus = 0;
            ffkeyn_safe(cs!("TUNIT"), colnum[ii], &mut keyname, &mut tstatus);
            fits_read_key_str(fptr, &keyname, &mut svalue, None, &mut tstatus);

            if tstatus == 0 {
                ffkeyn_safe(cs!("CUNIT"), (ii + 1) as c_int, &mut keyname, &mut tstatus);
                fits_write_key_str(
                    histptr,
                    &keyname,
                    &svalue,
                    Some(cs!("Coordinate Units")),
                    &mut tstatus,
                );
            }
        }

        /*  CRPIXn  - Reference Pixel choose first pixel in new image as ref. pix. */
        dvalue = 1.0;
        tstatus = 0;
        ffkeyn_safe(cs!("CRPIX"), (ii + 1) as c_int, &mut keyname, &mut tstatus);
        ffpky_safe(
            histptr,
            KeywordDatatype::TDOUBLE(&dvalue),
            &keyname,
            Some(cs!("Reference Pixel")),
            &mut tstatus,
        );

        /*  CRVALn - Value at the location of the reference pixel */
        dvalue = 1.0;
        tstatus = 0;
        ffkeyn_safe(cs!("CRVAL"), (ii + 1) as c_int, &mut keyname, &mut tstatus);
        ffpky_safe(
            histptr,
            KeywordDatatype::TDOUBLE(&dvalue),
            &keyname,
            Some(cs!("Reference Value")),
            &mut tstatus,
        );

        /*  CDELTn - unit size of pixels  */
        dvalue = 1.0;
        tstatus = 0;
        dvalue = 1.0;
        ffkeyn_safe(cs!("CDELT"), (ii + 1) as c_int, &mut keyname, &mut tstatus);
        ffpky_safe(
            histptr,
            KeywordDatatype::TDOUBLE(&dvalue),
            &keyname,
            Some(cs!("Pixel size")),
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
pub(crate) fn fits_write_keys_histo_safe(
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
pub(crate) fn fits_rebin_wcs_safe(
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
pub(crate) fn fits_rebin_wcsd_safe(
    fptr: &mut fitsfile, /* I - pointer to table to be binned           */
    naxis: c_int,        /* I - number of axes in the histogram image   */
    amin: &mut [f64],    /* I - first pixel include in each axis        */
    binsize: &mut [f64], /* I - binning factor for each axis            */
    status: &mut c_int,
) -> c_int {
    todo!();
}

/*--------------------------------------------------------------------------*/
/// Single-precision version
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn fits_make_hist(
    fptr: *mut fitsfile,      /* IO - pointer to table with X and Y cols; */
    histptr: *const fitsfile, /* I - pointer to output FITS image      */
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
        let histptr = histptr.as_ref().expect(NULL_MSG);
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
pub(crate) fn fits_make_hist_safe(
    fptr: &mut fitsfile,          /* IO - pointer to table with X and Y cols; */
    histptr: &fitsfile,           /* I - pointer to output FITS image      */
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
    fptr: &mut fitsfile,              /* IO - pointer to table with X and Y cols; */
    histptr: &fitsfile,               /* I - pointer to output FITS image      */
    datatypes: Option<&mut [c_int]>,  /*  I - datatype of input (or 0 for auto) */
    bitpix: c_int,                    /* I - datatype for image: 16, 32, -32, etc    */
    naxis: c_int,                     /* I - number of axes in the histogram image   */
    naxes: &[c_long],                 /* I - size of axes in the histogram image   */
    colnum: &[c_int],                 /* I - column numbers (array length = naxis)   */
    colexpr: Option<&[&[c_char]; 4]>, /* I - optional expression instead of column */
    amin: &[f64],                     /* I - minimum histogram value, for each axis */
    amax: &[f64],                     /* I - maximum histogram value, for each axis */
    binsize: &[f64],                  /* I - bin size along each axis               */
    weight: f64,     /* I - binning weighting factor (0 or DOUBLENULLVALUE means null) */
    wtcolnum: c_int, /* I - optional keyword or col for weight*/
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
    todo!();
}

/*--------------------------------------------------------------------------*/
/// Double-precision version, non-extended syntax
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn fits_make_histd(
    fptr: *mut fitsfile,      /* IO - pointer to table with X and Y cols; */
    histptr: *const fitsfile, /* I - pointer to output FITS image      */
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
        let histptr = histptr.as_ref().expect(NULL_MSG);
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
pub(crate) fn fits_make_histd_safe(
    fptr: &mut fitsfile,          /* IO - pointer to table with X and Y cols; */
    histptr: &fitsfile,           /* I - pointer to output FITS image      */
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
        cs!("NAXIS2"),
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

/*---------------------------------------------------------------------------*/
/// Iterator work function which evaluates a parser result and computes
/// min max value
fn histo_minmax_expr_workfn(
    totalrows: c_long,         /* I - Total rows to be processed     */
    offset: c_long,            /* I - Number of rows skipped at start*/
    firstrow: c_long,          /* I - First row of this iteration    */
    nrows: c_long,             /* I - Number of rows in this iter    */
    nCols: c_int,              /* I - Number of columns in use       */
    colData: &mut iteratorCol, /* IO- Column information/data        */
    userPtr: *mut c_void,      /* I - Data handling instructions     */
) -> c_int {
    todo!()
}

/*--------------------------------------------------------------------------*/
/// Simple utility routine to compute the min and max value in an expression
fn fits_get_expr_minmax(
    fptr: &mut fitsfile,
    expr: &mut [c_char],
    datamin: &mut f64,
    datamax: &mut f64,
    datatype: &mut c_int,
    status: &mut c_int,
) -> c_int {
    todo!()
}

/*--------------------------------------------------------------------------*/
/// Interator work function that writes out the histogram.
/// The histogram values are calculated by another work function, ffcalchisto.
/// This work function only gets called once, and totaln = nvalues.
fn ffwritehisto(
    totaln: c_long,
    pixoffset: c_long,
    firstn: c_long,
    nvalues: c_long,
    narrays: c_int,
    imagepars: &mut iteratorCol,
    userPointer: *mut c_void,
) -> c_int {
    todo!();
}

/*--------------------------------------------------------------------------*/
/// Interator work function that calculates values for the 2D histogram.
fn ffcalchist(
    totalrows: c_long,
    offset: c_long,
    firstrow: c_long,
    nrows: c_long,
    ncols: c_int,
    colpars: &mut iteratorCol,
    userPointer: *mut c_void,
) -> c_int {
    todo!()
}
