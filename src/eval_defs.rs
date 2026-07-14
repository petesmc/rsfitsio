use core::ffi::c_void;

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

#[derive(Debug)]
pub struct DataInfo {
    pub name: [c_char; MAXVARNAME + 1],
    pub dtype: c_int,
    pub nelem: c_long,
    pub naxis: c_int,
    pub naxes: [c_long; MAXDIMS as usize],
    pub undef: Option<Box<[c_char]>>,
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
            undef: None,
            data: core::ptr::null_mut(),
        }
    }
}

/// A single zeroed, 8-byte-aligned raw heap allocation. Node result buffers are
/// kept as raw allocations (rather than `Box`/`Vec`, which are `noalias`) so the
/// union's raw *views* into them stay valid under Stacked Borrows. Freed on drop.
struct RawBuf {
    ptr: *mut u8,
    /// Number of `u64` words backing the allocation.
    words: usize,
}

impl RawBuf {
    /// Allocate a zeroed buffer of at least `len` bytes, 8-byte aligned (buffers
    /// store f64/c_long and arrays of pointers, needing the alignment the former
    /// `calloc` provided), or `None` on allocation failure.
    fn new(len: usize) -> Option<Self> {
        let words = len.div_ceil(8).max(1);
        let mut v: Vec<u64> = Vec::new();
        v.try_reserve_exact(words).ok()?;
        v.resize(words, 0);
        let raw = Box::into_raw(v.into_boxed_slice());
        Some(RawBuf {
            ptr: raw as *mut u8,
            words,
        })
    }

    fn as_ptr(&self) -> *mut u8 {
        self.ptr
    }
}

impl Drop for RawBuf {
    fn drop(&mut self) {
        // Reconstitute the boxed [u64] slice and drop it, freeing via the same
        // (Rust global) allocator it was created with.
        unsafe {
            drop(Box::from_raw(core::ptr::slice_from_raw_parts_mut(
                self.ptr.cast::<u64>(),
                self.words,
            )));
        }
    }
}

/// Owns the heap storage for one node's result buffer. Numeric nodes use only
/// `primary` (the data-plus-undef buffer). STRING/BITSTR nodes use `primary` as
/// the array of row pointers (`strptr`) and `secondary` as the single backing
/// character buffer those pointers index into. Everything is freed on drop.
pub struct NodeBuf {
    primary: RawBuf,
    secondary: Option<RawBuf>,
}

impl NodeBuf {
    /// Numeric node buffer: `len` bytes of data followed by the undef flags.
    pub fn new_numeric(len: usize) -> Option<Self> {
        Some(NodeBuf {
            primary: RawBuf::new(len)?,
            secondary: None,
        })
    }

    /// String node buffer: `ptrs_len` bytes for the row-pointer array plus a
    /// `backing_len`-byte character buffer they point into.
    pub fn new_string(ptrs_len: usize, backing_len: usize) -> Option<Self> {
        let primary = RawBuf::new(ptrs_len)?;
        let secondary = RawBuf::new(backing_len)?;
        Some(NodeBuf {
            primary,
            secondary: Some(secondary),
        })
    }

    /// Raw view of the primary buffer (numeric data, or the `strptr` array).
    pub fn as_ptr(&self) -> *mut u8 {
        self.primary.as_ptr()
    }

    /// Raw view of the string backing buffer (only for string nodes).
    pub fn backing_ptr(&self) -> *mut u8 {
        self.secondary
            .as_ref()
            .expect("backing_ptr on a numeric NodeBuf")
            .as_ptr()
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

impl core::fmt::Debug for data_union {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "data_union {{ long: {:?} }}", unsafe { self.lng })
    }
}

impl Default for data_union {
    fn default() -> Self {
        data_union {
            ptr: core::ptr::null_mut(),
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
    /// Rust-owned backing store for the per-node result buffers of numeric
    /// (DOUBLE/LONG/BOOLEAN) nodes, indexed by node number. When a slot is
    /// `Some`, the node's `value.data.ptr`/`value.undef` are raw *views* into
    /// this box (data followed by the undef flags), and the memory is owned and
    /// freed here rather than via malloc/free. A `None` slot means the node's
    /// pointer (if any) is owned elsewhere (legacy malloc, e.g. GTI data).
    /// See `Allocate_Ptrs` / `free_node_data` in eval_y.
    pub node_buffers: Vec<Option<NodeBuf>>,
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
        self.node_buffers = Default::default();
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
