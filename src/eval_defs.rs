use std::ffi::c_void;

use crate::c_types::{c_char, c_int, c_long};

use crate::eval_l::yyguts_t;
use crate::eval_tab::FITS_PARSER_YYSTYPE;
use crate::fitsio::{LONGLONG, PixelFilter, fitsfile, iteratorCol};

pub const MAXDIMS: c_int = 5;
pub const MAXSUBS: c_int = 10;
pub const MAXVARNAME: usize = 80;
pub const CONST_OP: c_int = -1000;
pub const P_ERROR: c_int = -1;
pub const MAX_STRLEN: c_int = 256;
pub const MAX_STRLEN_S: &str = "255";

pub(crate) type yyscan_t<'a> = &'a mut yyguts_t;

#[derive(Debug, Copy, Clone)]
pub struct DataInfo {
    pub name: [c_char; MAXVARNAME + 1],
    pub dtype: c_int,
    pub nelem: c_long,
    pub naxis: c_int,
    pub naxes: [c_long; MAXDIMS as usize],
    pub undef: *mut c_char,
    pub data: *mut c_void,
}

impl Default for DataInfo {
    fn default() -> Self {
        DataInfo {
            name: [0; MAXVARNAME + 1],
            dtype: 0,
            nelem: 0,
            naxis: 0,
            naxes: [0; MAXDIMS as usize],
            undef: std::ptr::null_mut(),
            data: std::ptr::null_mut(),
        }
    }
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

impl std::fmt::Debug for data_union {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "data_union {{ long: {:?} }}", unsafe { self.lng })
    }
}

impl Default for data_union {
    fn default() -> Self {
        data_union {
            ptr: std::ptr::null_mut(),
        }
    }
}

#[derive(Default, Debug, Copy, Clone)]
pub struct lval {
    pub nelem: c_long,
    pub naxis: c_int,
    pub naxes: [c_long; MAXDIMS as usize],
    pub undef: *mut c_char,
    pub data: data_union,
}

#[derive(Default, Debug, Copy, Clone)]
pub struct Node {
    pub operation: c_int,
    pub DoOp: Option<fn(p: &mut ParseData, this_node_idx: usize)>,
    pub nSubNodes: c_int,
    pub SubNodes: [usize; MAXSUBS as usize],
    pub ntype: c_int,
    pub value: lval,
}

#[derive(Default)]
pub struct ParseData {
    pub def_fptr: *mut fitsfile,
    pub getData: Option<
        fn(p: &mut ParseData, dataName: &[c_char], dataValue: &mut FITS_PARSER_YYSTYPE) -> c_int,
    >,
    pub loadData: Option<
        fn(
            p: &mut ParseData,
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
    pub expr: Option<Box<[u8]>>,
    pub index: c_int,
    pub is_eobuf: c_int,
    pub Nodes: Vec<Node>,
    pub nNodes: c_int,
    pub nNodesAlloc: c_int,
    pub resultNode: c_int,
    pub firstRow: c_long,
    pub nRows: c_long,
    pub nCols: c_int,
    pub nElements: c_long,
    pub nAxis: c_int,
    pub nAxes: [c_long; MAXDIMS as usize],
    pub colData: Vec<iteratorCol>, // This is a list
    pub varData: Vec<DataInfo>,
    pub pixFilter: *mut PixelFilter,
    pub firstDataRow: c_long,
    pub nDataRows: c_long,
    pub totalRows: c_long,
    pub nPrevDataRows: c_long,
    pub datatype: c_int,
    pub hdutype: c_int,
    pub status: c_int,
}

impl ParseData {
    pub(crate) fn reset(&mut self) {
        // Reset all fields
        self.def_fptr = Default::default();
        self.getData = Default::default();
        self.loadData = Default::default();
        self.compressed = Default::default();
        self.timeCol = Default::default();
        self.parCol = Default::default();
        self.valCol = Default::default();
        self.expr = Default::default();
        self.index = Default::default();
        self.is_eobuf = Default::default();
        self.Nodes = Default::default();
        self.nNodes = Default::default();
        self.nNodesAlloc = Default::default();
        self.resultNode = Default::default();
        self.firstRow = Default::default();
        self.nRows = Default::default();
        self.nCols = Default::default();
        self.nElements = Default::default();
        self.nAxis = Default::default();
        self.nAxes = Default::default();
        self.colData = Default::default();
        self.varData = Default::default();
        self.pixFilter = Default::default();
        self.firstDataRow = Default::default();
        self.nDataRows = Default::default();
        self.totalRows = Default::default();
        self.nPrevDataRows = Default::default();
        self.datatype = Default::default();
        self.hdutype = Default::default();
        self.status = Default::default();
    }
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

#[derive(Default)]
pub(crate) struct ParseStatusVariables {
    /* These variables were 'static' in fits_parse_workfn() */
    pub(crate) Data: *mut c_void,
    pub(crate) Null: *mut c_void,
    pub(crate) datasize: c_int,
    pub(crate) lastRow: c_long,
    pub(crate) repeat: c_long,
    pub(crate) resDataSize: c_long,
    pub(crate) jnull: LONGLONG,
    pub(crate) userInfo: *mut parseInfo,
    pub(crate) zeros: [c_long; 4],
}

#[derive(Default)]
pub(crate) struct parseInfo {
    pub(crate) datatype: c_int, /* Data type to cast parse results into for user       */
    pub(crate) dataPtr: *mut c_void, /* Pointer to array of results, NULL if to use iterCol */
    pub(crate) nullPtr: *mut c_void, /* Pointer to nulval, use zero if NULL                 */
    pub(crate) maxRows: c_long, /* Max No. of rows to process, -1=all, 0=1 iteration   */
    pub(crate) anyNull: c_int,  /* Flag indicating at least 1 undef value encountered  */
    pub(crate) parseData: *mut ParseData, /* Pointer to parser configuration */
    pub(crate) parseVariables: ParseStatusVariables,
}
