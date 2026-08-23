//! Shared definitions for the expression parser.
//!
//! The types the three parser files agree on: the evaluation tree's `Node` and
//! `NodeValue`, the `ParseData` state the parse builds up, and the column and
//! keyword descriptions it collects. All are crate-private. Corresponds to
//! CFITSIO's `eval_defs.h`.
//!
//! See [`crate::eval_f`] for how the three files fit together.
#![warn(missing_docs)]

use alloc::alloc::{alloc_zeroed, dealloc};
use alloc::rc::Rc;
use core::alloc::Layout;
use core::ffi::c_void;

use crate::c_types::{c_char, c_int, c_long};

use crate::eval_l::yyguts_t;
/* Not sure why this is needed but it is */
use crate::eval_tab::FITS_PARSER_YYSTYPE;
use crate::fitsio::{LONGLONG, PixelFilter, fitsfile, iteratorCol};
use crate::region::SAORegion;

/// Maximum number of dimensions a parsed vector column may have.
pub const MAXDIMS: c_int = 5;
/// Maximum number of subscripts in a column reference.
pub const MAXSUBS: c_int = 10;
/// Maximum length of a variable name in an expression.
pub const MAXVARNAME: usize = 80;
/// Operation code marking a node as a constant rather than a computation.
pub const CONST_OP: c_int = -1000;
/// Parser error indicator.
pub const P_ERROR: c_int = -1;
/// Maximum length of a string value in an expression, including the
/// terminator.
pub const MAX_STRLEN: c_int = 256;
/// [`MAX_STRLEN`] less the terminator, as a string, for building format
/// specifiers.
pub const MAX_STRLEN_S: &str = "255";

/// The lexer's scanner state. An opaque pointer in the C.
pub(crate) type yyscan_t<'a> = &'a mut yyguts_t;

/// The sort of a parser value: which of `eval.y`'s four nonterminals a node
/// belongs to.
///
/// The transpilation stored this as a bare `c_int` holding a lexer token id, so
/// every test read `node.ntype == fits_parser_yytokentype::DOUBLE as c_int`.
///
/// The discriminants are the original token numbers, and the derived `Ord` is
/// the numeric promotion order that `eval.y` depended on -- its `%token`
/// block carries the comment *"First 3 must be in order of increasing
/// promotion"*, which is why `PROMOTE` could be a `>` comparison. Keeping the
/// order here makes that explicit rather than accidental.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, PartialOrd, Ord)]
#[repr(i32)]
pub(crate) enum ValueSort {
    #[default]
    Boolean = 258,
    Long = 259,
    Double = 260,
    String = 261,
    Bits = 262,
}

impl ValueSort {
    /// The token number, for the places that still round-trip through `c_int`.
    pub(crate) fn code(self) -> c_int {
        self as c_int
    }

    /// The sort for a token number, or `None` if it names no sort.
    pub(crate) fn from_code(v: c_int) -> Option<ValueSort> {
        Some(match v {
            258 => ValueSort::Boolean,
            259 => ValueSort::Long,
            260 => ValueSort::Double,
            261 => ValueSort::String,
            262 => ValueSort::Bits,
            _ => return None,
        })
    }
}

#[derive(Debug)]
pub(crate) struct DataInfo {
    pub name: [c_char; MAXVARNAME + 1],
    pub dtype: ValueSort,
    pub nelem: c_long,
    pub naxis: c_int,
    pub naxes: [c_long; MAXDIMS as usize],
    pub undef: Option<Vec<c_char>>,
    pub data: *mut c_void,
    /// `Bits` columns only: the row-pointer array `data` points at, and the
    /// single character block those pointers index into.
    ///
    /// `Setup_DataArrays` malloc'd both and `ffcprs` freed them. They are owned
    /// here instead, so they go when the `ParseData` does; `data` is still the
    /// raw `*mut c_void` the nodes read, taken from `bit_rows`.
    pub bit_rows: Vec<*mut c_char>,
    /// The character block `bit_rows` indexes into. See [`DataInfo::bit_rows`].
    pub bit_block: Vec<c_char>,
}

impl Default for DataInfo {
    fn default() -> Self {
        DataInfo {
            name: [0; MAXVARNAME + 1],
            dtype: ValueSort::default(),
            nelem: 0,
            naxis: 0,
            naxes: [0; MAXDIMS as usize],
            undef: None,
            data: core::ptr::null_mut(),
            bit_rows: Vec::new(),
            bit_block: Vec::new(),
        }
    }
}

/// Which element type a [`NodeValue::Buffer`] holds.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum BufferKind {
    Long,
    Double,
    Logical,
    /// An array of row pointers into one contiguous character block.
    Text,
}

impl BufferKind {
    /// The alignment a buffer of this kind must have.
    ///
    /// `Allocate_Ptrs` sizes its blocks in bytes -- the numeric one packs the
    /// data and the trailing `undef` flags into a single allocation -- so the
    /// alignment cannot come from a `Layout::array::<T>` and has to be stated.
    const fn align(self) -> usize {
        match self {
            BufferKind::Long => align_of::<c_long>(),
            BufferKind::Double => align_of::<f64>(),
            BufferKind::Logical => align_of::<c_char>(),
            BufferKind::Text => align_of::<*mut c_char>(),
        }
    }
}

/// Allocate `bytes` zeroed bytes aligned for `kind`, returning the block and
/// the [`Layout`] it must be released with.
///
/// This is the C's `malloc`/`calloc` for the expression engine's row buffers.
/// It is `alloc_zeroed` rather than `libc::calloc` so the blocks come from
/// whatever global allocator the program installed, and the layout travels with
/// the pointer because Rust's `dealloc` needs it back -- C's `free` did not.
///
/// A zero-length request is rounded up to one byte, matching the C, where
/// `malloc(0)` returns a non-null pointer that the callers treat as success.
/// Returns `None` on allocation failure, which is the null the C checked for.
pub(crate) fn alloc_row_bytes(kind: BufferKind, bytes: usize) -> Option<(*mut c_void, Layout)> {
    let layout = Layout::from_size_align(bytes.max(1), kind.align()).ok()?;
    // SAFETY: the layout has a non-zero size.
    let p = unsafe { alloc_zeroed(layout) };
    if p.is_null() {
        None
    } else {
        Some((p.cast::<c_void>(), layout))
    }
}

/// Release a block from [`alloc_row_bytes`].
///
/// # Safety
/// `ptr` must have come from [`alloc_row_bytes`] with this exact `layout`, and
/// must not be used again.
pub(crate) unsafe fn free_row_bytes(ptr: *mut c_void, layout: Layout) {
    if !ptr.is_null() {
        unsafe { dealloc(ptr.cast::<u8>(), layout) };
    }
}

/// The value carried by a [`Node`].
///
/// A node holds either a constant folded at parse time, or a buffer with one
/// element per (row x element) that the `Do_*` routines fill in.
///
/// This was an untagged union, which is how `ACCUM(BOOLCOL)` came to read a
/// `char` buffer through a `*long`: nothing checked that the arm being read
/// matched the node's sort. The accessors below assert it, so that class of
/// mistake is a panic in a debug build instead of undefined behaviour.
///
/// The buffer is still a raw allocation rather than a `Vec`. Two things stand
/// in the way, neither of them `Copy`:
///
/// * The `Do_*` routines read two operand buffers while writing a third, all
///   indexed out of `ParseData::Nodes`. That needs one `&mut` and two `&` at
///   distinct indices -- `split_at_mut` territory, since the two operands may
///   be the same node.
/// * Only *computed* nodes own their buffer. `Evaluate_Parser` re-points every
///   *column* node at `varData[column].data` before each evaluation, and for a
///   String column that is the iterator's own array, laid out by `ffiter`. An
///   owning arm would need a borrowed one beside it, and the ownership question
///   moves to `ffiter`.
///
/// What the allocation *is* no longer depends on the C runtime: the block comes
/// from [`alloc_row_bytes`], i.e. the global Rust allocator, and the variant
/// carries the [`Layout`] to release it with. `Owned::None` is how a column
/// node's borrowed view into `varData` is distinguished from a buffer the node
/// must free, which the callers used to have to know out of band.
///
/// The `Text` arm makes this 264 bytes where the others need 16. That is the
/// union's own footprint, so it is not a regression.
#[allow(clippy::large_enum_variant)]
#[derive(Clone, Default)]
pub(crate) enum NodeValue {
    /// A freshly allocated node, or one whose buffer has been freed.
    #[default]
    Empty,
    Long(c_long),
    Double(f64),
    Logical(c_char),
    /// A string or bit-string scalar, NUL-terminated within the array. Also
    /// used as the scratch buffer for a scalar string result.
    Text([c_char; MAX_STRLEN as usize]),
    /// A row buffer, allocated by `Allocate_Ptrs` and released by
    /// `free_node_buffer`.
    Buffer {
        kind: BufferKind,
        ptr: *mut c_void,
        /// How `ptr` (and, for `Text`, the character block behind it) is
        /// released.
        own: Owned,
    },
    /// A region, by index into [`ParseData::regions`].
    ///
    /// The region itself is an `Rc` held there rather than in the node,
    /// because `NodeValue` is `Copy` and an `Rc` arm would end that. The node
    /// carrying an index instead is what lets several nodes share one region
    /// without any of them owning it.
    Region(usize),
}

impl core::fmt::Debug for NodeValue {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            NodeValue::Empty => write!(f, "Empty"),
            NodeValue::Long(v) => write!(f, "Long({v})"),
            NodeValue::Double(v) => write!(f, "Double({v})"),
            NodeValue::Logical(v) => write!(f, "Logical({v})"),
            NodeValue::Text(_) => write!(f, "Text(..)"),
            NodeValue::Buffer { kind, ptr, own } => {
                write!(f, "Buffer({kind:?}, {ptr:?}, {own:?})")
            }
            NodeValue::Region(i) => write!(f, "Region({i})"),
        }
    }
}

impl NodeValue {
    /// The scalar integer value. Panics in debug if the node holds something
    /// else, which means an operator was applied to the wrong sort.
    pub(crate) fn lng(&self) -> c_long {
        match self {
            NodeValue::Long(v) => *v,
            other => wrong_arm("Long", other),
        }
    }

    pub(crate) fn dbl(&self) -> f64 {
        match self {
            NodeValue::Double(v) => *v,
            other => wrong_arm("Double", other),
        }
    }

    pub(crate) fn log(&self) -> c_char {
        match self {
            NodeValue::Logical(v) => *v,
            other => wrong_arm("Logical", other),
        }
    }

    /// The scalar string, as the fixed-size array the C code indexes.
    pub(crate) fn text(&self) -> &[c_char; MAX_STRLEN as usize] {
        match self {
            NodeValue::Text(v) => v,
            other => wrong_arm("Text", other),
        }
    }

    /// A writable view of the scalar string buffer, for the `*_safe` string
    /// helpers.
    ///
    /// The C wrote into `data.str` of a node that had no value yet, so this
    /// installs an empty `Text` first when the node holds something else.
    pub(crate) fn text_mut(&mut self) -> &mut [c_char; MAX_STRLEN as usize] {
        if !matches!(self, NodeValue::Text(_)) {
            *self = NodeValue::Text([0; MAX_STRLEN as usize]);
        }
        match self {
            NodeValue::Text(v) => v,
            _ => unreachable!(),
        }
    }

    /// A writable pointer to the scalar string buffer.
    ///
    /// The C wrote into `data.str` of a node that had no value yet, so this
    /// installs an empty `Text` first when the node holds something else.
    pub(crate) fn text_mut_ptr(&mut self) -> *mut c_char {
        self.text_mut().as_mut_ptr()
    }

    fn buffer(&self, want: BufferKind) -> *mut c_void {
        match self {
            NodeValue::Buffer { kind, ptr, .. } if *kind == want => *ptr,
            /* an unallocated node reads as a null buffer, as the union did */
            NodeValue::Empty => core::ptr::null_mut(),
            other => wrong_arm(&format!("Buffer({want:?})"), other),
        }
    }

    pub(crate) fn lng_buf(&self) -> *mut c_long {
        self.buffer(BufferKind::Long).cast()
    }
    pub(crate) fn dbl_buf(&self) -> *mut f64 {
        self.buffer(BufferKind::Double).cast()
    }
    pub(crate) fn log_buf(&self) -> *mut c_char {
        self.buffer(BufferKind::Logical).cast()
    }
    pub(crate) fn str_buf(&self) -> *mut *mut c_char {
        self.buffer(BufferKind::Text).cast()
    }

    /// The buffer pointer whatever its element type, for the paths that only
    /// allocate or free it.
    pub(crate) fn raw(&self) -> *mut c_void {
        match self {
            NodeValue::Buffer { ptr, .. } => *ptr,
            _ => core::ptr::null_mut(),
        }
    }

    /// The index into [`ParseData::regions`] this node refers to.
    pub(crate) fn region(&self) -> usize {
        match self {
            NodeValue::Region(i) => *i,
            other => wrong_arm("Region", other),
        }
    }

    /// Point the node at a buffer somebody else owns.
    ///
    /// `Evaluate_Parser` uses this for column nodes, whose data lives in
    /// `varData` (and, for a String column, in the iterator's own array).
    /// Nothing is freed through such a node.
    pub(crate) fn set_borrowed_buffer(&mut self, kind: BufferKind, ptr: *mut c_void) {
        *self = NodeValue::Buffer {
            kind,
            ptr,
            own: Owned::None,
        };
    }

    /// Give the node a buffer from [`alloc_row_bytes`] to own.
    pub(crate) fn set_owned_buffer(&mut self, kind: BufferKind, ptr: *mut c_void, layout: Layout) {
        *self = NodeValue::Buffer {
            kind,
            ptr,
            own: Owned::Block(layout),
        };
    }

    /// Record the single character block a `Text` buffer's row pointers index
    /// into, which `Allocate_Ptrs` stores at `ptr[0]`.
    ///
    /// # Panics
    /// If the node does not hold an owned buffer.
    pub(crate) fn set_text_block(&mut self, block: *mut c_char, block_layout: Layout) {
        match self {
            NodeValue::Buffer {
                own: own @ Owned::Block(_),
                ..
            } => {
                let Owned::Block(rows) = *own else {
                    unreachable!()
                };
                *own = Owned::BlockAndText {
                    rows,
                    block,
                    block_layout,
                };
            }
            other => wrong_arm("an owned Buffer", other),
        }
    }

    /// A pointer to the scalar payload, for the conversion helpers that take a
    /// `void *` and a datatype.
    ///
    /// The union put the payload at offset 0, so the C could simply cast the
    /// whole thing; an enum has a discriminant in front, so the address has to
    /// come from the variant. Borrow-wise this is a pointer *into `self`*, so
    /// it lives exactly as long as the borrow.
    pub(crate) fn scalar_ptr(&self) -> *const c_void {
        match self {
            NodeValue::Long(v) => core::ptr::from_ref::<c_long>(v).cast(),
            NodeValue::Double(v) => core::ptr::from_ref::<f64>(v).cast(),
            NodeValue::Logical(v) => core::ptr::from_ref::<c_char>(v).cast(),
            NodeValue::Text(v) => v.as_ptr().cast(),
            NodeValue::Empty | NodeValue::Buffer { .. } | NodeValue::Region(_) => core::ptr::null(),
        }
    }

    /// Release whatever this node owns and mark it empty.
    ///
    /// A borrowed buffer (`Owned::None`) is left alone; that is what makes it
    /// safe to call on a column node, which points into `varData`. For a
    /// `Text` buffer both allocations go -- the character block first, then the
    /// row-pointer array -- so callers no longer have to free the block by hand
    /// before calling this.
    ///
    /// # Safety
    /// Nothing else may still refer to the buffer.
    pub(crate) unsafe fn free_buffer(&mut self) {
        if let NodeValue::Buffer { ptr, own, .. } = *self {
            match own {
                Owned::None => {}
                Owned::Block(layout) => unsafe { free_row_bytes(ptr, layout) },
                Owned::BlockAndText {
                    rows,
                    block,
                    block_layout,
                } => unsafe {
                    /* the row pointers all index one block, allocated at [0] */
                    free_row_bytes(block.cast::<c_void>(), block_layout);
                    free_row_bytes(ptr, rows);
                },
            }
        }
        *self = NodeValue::Empty;
    }
}

/// How a [`NodeValue::Buffer`]'s storage is released.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) enum Owned {
    /// The node only borrows the block; something else frees it. This is a
    /// column node pointing into `varData`.
    #[default]
    None,
    /// The node owns one block, allocated with this layout.
    Block(Layout),
    /// A `Text` buffer: a row-pointer array with the `rows` layout, plus the
    /// single character block at `ptr[0]` that every row pointer indexes into.
    BlockAndText {
        rows: Layout,
        block: *mut c_char,
        block_layout: Layout,
    },
}

#[cold]
#[track_caller]
fn wrong_arm(want: &str, got: &NodeValue) -> ! {
    panic!("node value is {got:?}, expected {want}");
}

#[derive(Default, Debug, Clone)]
pub(crate) struct lval {
    pub nelem: c_long,
    pub naxis: c_int,
    pub naxes: [c_long; MAXDIMS as usize],
    pub undef: *mut c_char,
    pub data: NodeValue,
}

impl lval {
    /// The raw start of row `row` of a String/Bits row buffer.
    ///
    /// # Safety
    /// See [`lval::str_row`].
    pub(crate) unsafe fn str_row_ptr(&self, row: c_long) -> *mut c_char {
        unsafe { *self.data.str_buf().offset(row as isize) }
    }

    /// Row `row` of a String/Bits row buffer.
    ///
    /// `Allocate_Ptrs` reserves `nelem + 1` characters per row -- the string
    /// and its terminator -- so that is the length of the slice. This is what
    /// lets the `Do_*` routines reach a row with the `*_safe` string helpers
    /// instead of an unbounded `strcpy` on a bare `*mut c_char`.
    ///
    /// # Safety
    /// - The node must hold a live `Text` row buffer from `Allocate_Ptrs`, and
    ///   `row` must be below the `nRows` it was allocated for.
    /// - The returned lifetime is unbounded: the slice borrows the heap block,
    ///   not `self`, which is what lets a caller hold a destination row and two
    ///   operand rows at the same time. The caller must not create two
    ///   overlapping `&mut` rows. In the engine the destination is always a
    ///   freshly allocated node, so it cannot alias an operand.
    pub(crate) unsafe fn str_row<'a>(&self, row: c_long) -> &'a [c_char] {
        unsafe { core::slice::from_raw_parts(self.str_row_ptr(row), self.nelem as usize + 1) }
    }

    /// Row `row` of a String/Bits row buffer, writable.
    ///
    /// # Safety
    /// See [`lval::str_row`].
    pub(crate) unsafe fn str_row_mut<'a>(&self, row: c_long) -> &'a mut [c_char] {
        unsafe { core::slice::from_raw_parts_mut(self.str_row_ptr(row), self.nelem as usize + 1) }
    }
}

#[derive(Default, Debug, Clone)]
pub(crate) struct Node {
    pub operation: c_int,
    pub DoOp: Option<fn(p: &mut ParseData, this_node_idx: usize)>,
    pub nSubNodes: c_int,
    pub SubNodes: [usize; MAXSUBS as usize],
    pub ntype: ValueSort,
    /// Whether a GTI node's START/STOP columns are fully time-ordered.
    ///
    /// `New_GTI` records this and `Do_GTI` reads it back to pick the binary or
    /// the linear search. The C kept it in `this->type` of the START/STOP data
    /// node -- an `int` field otherwise holding a value sort -- which only
    /// worked because nothing typed that field. It gets its own home here.
    pub gtiOrdered: bool,
    pub value: lval,
}

impl Node {
    /// A constant folded at parse time.
    pub(crate) fn is_const(&self) -> bool {
        self.operation == CONST_OP
    }

    /// A reference to a table column, rather than something computed.
    ///
    /// `New_Column` stores the column as `-ColNum`, and `ColNum` indexes
    /// `varData` from zero -- so column 0 is `operation == 0`, and the test
    /// has to be `<= 0` rather than `< 0`. Getting that wrong classifies the
    /// first column as something computed, and the engine then reads a buffer
    /// the node never had.
    pub(crate) fn is_column(&self) -> bool {
        self.operation <= 0 && self.operation != CONST_OP
    }

    /// The index into [`ParseData::varData`] this column node refers to.
    pub(crate) fn column(&self) -> c_int {
        debug_assert!(self.is_column(), "not a column node: {}", self.operation);
        -self.operation
    }

    /// Whether this node computes something -- an operator, a cast or a
    /// function -- as opposed to being a constant or a column.
    ///
    /// `Node::operation` packs four encodings into one `c_int`: a function
    /// code (1001 and up), an operator (an ASCII character or a bison token,
    /// both positive), `CONST_OP`, or a negated column number. The engine
    /// tests `operation > 0` to mean "has operands to evaluate and a buffer of
    /// its own to free"; this says so.
    pub(crate) fn is_computed(&self) -> bool {
        self.operation > 0
    }
}

/// Fetches a single named value for the parser (ffcalc's keyword lookup).
pub(crate) type GetDataFn =
    fn(p: &mut ParseData, dataName: &[c_char], dataValue: &mut FITS_PARSER_YYSTYPE) -> c_int;

/// Loads a row range of one parser variable into the caller's buffers.
pub(crate) type LoadDataFn = fn(
    p: &mut ParseData,
    varNum: c_int,
    fRow: c_long,
    nRows: c_long,
    data: *mut c_void,
    undef: *mut c_char,
) -> c_int;

#[derive(Default)]
pub(crate) struct ParseData {
    pub def_fptr: *mut fitsfile,
    pub getData: Option<GetDataFn>,
    pub loadData: Option<LoadDataFn>,
    pub compressed: c_int,
    pub timeCol: c_int,
    pub parCol: c_int,
    pub valCol: c_int,
    pub expr: Option<Box<[u8]>>,
    pub index: c_int,
    pub is_eobuf: c_int,
    pub Nodes: Vec<Node>,
    /// Regions read by `New_REG`, shared with the nodes that filter on them.
    pub regions: Vec<Rc<SAORegion>>,
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
        self.regions = Default::default();
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
