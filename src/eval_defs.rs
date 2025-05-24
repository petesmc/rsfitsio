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

pub struct DataInfo<'a> {
    name: [c_char; MAXVARNAME as usize + 1],
    r#type: c_int,
    nelem: c_long,
    naxis: c_int,
    naxes: [c_long; MAXDIMS as usize],
    undef: &'a c_char,
    data: &'a c_void,
}

union data_union {
    dbl: f64,
    lng: c_long,
    log: c_char,
    str: [c_char; MAX_STRLEN as usize],
    dblptr: *mut f64,
    lngptr: *mut c_long,
    logptr: *mut c_char,
    strptr: *mut c_char,
    ptr: *mut c_void,
}

pub struct lval<'a> {
    nelem: c_long,
    naxis: c_int,
    naxes: [c_long; MAXDIMS as usize],
    undef: &'a c_char,
    data: data_union,
}

pub struct Node<'a> {
    operation: c_int,
    DoOp: fn(p: ParseData, this: Node) -> c_void,
    nSubNodes: c_int,
    SubNodes: [c_int; MAXSUBS as usize],
    r#type: c_int,
    value: &'a lval<'a>,
}

pub struct ParseData<'a> {
    def_fptr: &'a fitsfile,
    getData: fn(p: ParseData, dataName: c_char, dataValue: c_void) -> c_int,
    loadData: fn(
        p: ParseData,
        varNum: c_int,
        fRow: c_long,
        nRows: c_long,
        data: c_void,
        undef: c_char,
    ) -> c_int,
    compressed: c_int,
    timeCol: c_int,
    parCol: c_int,
    valCol: c_int,
    expr: c_char,
    index: c_int,
    is_eobuf: c_int,
    Nodes: &'a Node<'a>,
    nNodes: c_int,
    nNodesAlloc: c_int,
    resultNode: c_int,
    firstRow: c_long,
    nRows: c_long,
    nCols: c_int,
    nElements: c_long,
    nAxis: c_int,
    nAxes: [c_long; MAXDIMS as usize],
    colData: &'a iteratorCol,
    varData: &'a DataInfo<'a>,
    pixFilter: &'a PixelFilter,
    firstDataRow: c_long,
    nDataRows: c_long,
    totalRows: c_long,
    nPrevDataRows: c_long,
    datatype: c_int,
    hdutype: c_int,
    status: c_int,
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
    datatype: c_int,     /* Data type to cast parse results into for user       */
    dataPtr: &'a c_void, /* Pointer to array of results, NULL if to use iterCol */
    nullPtr: &'a c_void, /* Pointer to nulval, use zero if NULL                 */
    maxRows: c_long,     /* Max No. of rows to process, -1=all, 0=1 iteration   */
    anyNull: c_int,      /* Flag indicating at least 1 undef value encountered  */
    parseData: &'a ParseData<'a>, /* Pointer to parser configuration */
    parseVariables: ParseStatusVariables<'a>,
}

/* Not sure why this is needed but it is */
// pub type YYSTYPE  = FITS_PARSER_YYSTYPE;
/* How ParseData is accessed from the lexer, i.e. by yyextra */
//pub type YY_EXTRA_TYPE =  ParseData;
