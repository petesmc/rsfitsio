use std::ffi::c_void;

use crate::c_types::{c_char, c_int, c_long};

use crate::fitsio::{LONGLONG, PixelFilter, fitsfile, iteratorCol};

pub const MAXDIMS: c_int = 5;
pub const MAXSUBS: c_int = 10;
pub const MAXVARNAME: c_int = 80;
pub const CONST_OP: c_int = -1000;
pub const P_ERROR: c_int = -1;
pub const MAX_STRLEN: c_int = 256;
pub const MAX_STRLEN_S: &str = "255";

pub struct DataInfo {
    pub name: [c_char; MAXVARNAME as usize + 1],
    pub dtype: c_int,
    pub nelem: c_long,
    pub naxis: c_int,
    pub naxes: [c_long; MAXDIMS as usize],
    pub undef: *mut c_char,
    pub data: *mut c_void,
}

#[derive(Copy, Clone)]
pub union data_union {
    pub dbl: f64,
    pub lng: c_long,
    pub log: c_char,
    pub astr: [c_char; MAX_STRLEN as usize],
    pub dblptr: *mut f64,
    pub lngptr: *mut c_long,
    pub logptr: *mut c_char,
    pub strptr: *mut *mut c_char,
    pub ptr: *mut c_void,
}

impl Default for data_union {
    fn default() -> Self {
        data_union {
            ptr: std::ptr::null_mut(),
        }
    }
}

#[derive(Default, Copy, Clone)]
pub struct lval {
    pub nelem: c_long,
    pub naxis: c_int,
    pub naxes: [c_long; MAXDIMS as usize],
    pub undef: *mut c_char,
    pub data: data_union,
}

pub struct Node {
    pub operation: c_int,
    pub DoOp: Option<unsafe fn(p: *mut ParseData, this: *mut Node)>,
    pub nSubNodes: c_int,
    pub SubNodes: [c_int; MAXSUBS as usize],
    pub ntype: c_int,
    pub value: lval,
}

pub struct ParseData {
    pub def_fptr: *mut fitsfile,
    pub getData:
        Option<fn(p: *mut ParseData, dataName: *mut c_char, dataValue: *mut c_void) -> c_int>,
    pub loadData: Option<
        fn(
            p: *mut ParseData,
            varNum: c_int,
            fRow: c_long,
            nRows: c_long,
            data: *mut c_void,
            undef: *mut c_char,
        ) -> c_int,
    >,
    pub compressed: c_int,
    pub timeCol: c_int,
    pub parCol: c_int,
    pub valCol: c_int,
    pub expr: *mut c_char,
    pub index: c_int,
    pub is_eobuf: c_int,
    pub Nodes: *mut Node,
    pub nNodes: c_int,
    pub nNodesAlloc: c_int,
    pub resultNode: c_int,
    pub firstRow: c_long,
    pub nRows: c_long,
    pub nCols: c_int,
    pub nElements: c_long,
    pub nAxis: c_int,
    pub nAxes: [c_long; MAXDIMS as usize],
    pub colData: *mut iteratorCol, // This is a list
    pub varData: *mut DataInfo,
    pub pixFilter: *mut PixelFilter,
    pub firstDataRow: c_long,
    pub nDataRows: c_long,
    pub totalRows: c_long,
    pub nPrevDataRows: c_long,
    pub datatype: c_int,
    pub hdutype: c_int,
    pub status: c_int,
}

enum funcOp {
    rnd_fct = 1001,
    sum_fct,
    nelem_fct,
    sin_fct,
    cos_fct,
    tan_fct,
    asin_fct,
    acos_fct,
    atan_fct,
    sinh_fct,
    cosh_fct,
    tanh_fct,
    exp_fct,
    log_fct,
    log10_fct,
    sqrt_fct,
    abs_fct,
    atan2_fct,
    ceil_fct,
    floor_fct,
    round_fct,
    min1_fct,
    min2_fct,
    max1_fct,
    max2_fct,
    near_fct,
    circle_fct,
    box_fct,
    elps_fct,
    isnull_fct,
    defnull_fct,
    gtifilt_fct,
    regfilt_fct,
    ifthenelse_fct,
    row_fct,
    null_fct,
    median_fct,
    average_fct,
    stddev_fct,
    nonnull_fct,
    angsep_fct,
    gasrnd_fct,
    poirnd_fct,
    strmid_fct,
    strpos_fct,
    setnull_fct,
    gtiover_fct,
    gtifind_fct,
    elemnum_fct,
    axiselem_fct,
    array_fct,
}

pub struct ParseStatusVariables<'a> {
    /* These variables were 'static' in fits_parse_workfn() */
    Data: &'a c_void,
    Null: &'a c_void,
    datasize: c_int,
    lastRow: c_long,
    repeat: c_long,
    resDataSize: c_long,
    jnull: LONGLONG,
    userInfo: &'a parseInfo<'a>,
    zeros: [c_long; 4],
}

pub struct parseInfo<'a> {
    datatype: c_int,           /* Data type to cast parse results into for user       */
    dataPtr: &'a c_void,       /* Pointer to array of results, NULL if to use iterCol */
    nullPtr: &'a c_void,       /* Pointer to nulval, use zero if NULL                 */
    maxRows: c_long,           /* Max No. of rows to process, -1=all, 0=1 iteration   */
    anyNull: c_int,            /* Flag indicating at least 1 undef value encountered  */
    parseData: *mut ParseData, /* Pointer to parser configuration */
    parseVariables: ParseStatusVariables<'a>,
}

/* Not sure why this is needed but it is */
// pub type YYSTYPE  = FITS_PARSER_YYSTYPE;
/* How ParseData is accessed from the lexer, i.e. by yyextra */
//pub type YY_EXTRA_TYPE =  ParseData;
