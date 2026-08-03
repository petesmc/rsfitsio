use core::ffi::c_void;

use crate::c_types::{c_char, c_double, c_int, c_long};

use crate::fitsio::{LONGLONG, PixelFilter, fitsfile, iteratorCol};

pub const MAXDIMS: c_int = 5;
pub const MAXSUBS: c_int = 10;
pub const MAXVARNAME: usize = 80;
pub const P_ERROR: c_int = -1;
pub const MAX_STRLEN: c_int = 256;
pub const MAX_STRLEN_S: &str = "255";

#[derive(Debug)]
pub(crate) struct DataInfo {
    pub name: [c_char; MAXVARNAME + 1],
    pub dtype: ValueSort,
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
            dtype: ValueSort::default(),
            nelem: 0,
            naxis: 0,
            naxes: [0; MAXDIMS as usize],
            undef: None,
            data: core::ptr::null_mut(),
        }
    }
}

/// The function a `Node` computes, stored in `Node::operation` for nodes whose
/// `DoOp` is `Do_Func`.
///
/// The C had these as a `typedef enum FuncOp` in `eval_defs.h`; the
/// transpilation flattened them to bare integers, so `Do_Func` was a match over
/// fifty numeric literals. The values are unchanged from the C.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[repr(i32)]
pub(crate) enum FuncOp {
    /// `rnd_fct`
    Rnd = 1001,
    /// `sum_fct`
    Sum = 1002,
    /// `nelem_fct`
    Nelem = 1003,
    /// `sin_fct`
    Sin = 1004,
    /// `cos_fct`
    Cos = 1005,
    /// `tan_fct`
    Tan = 1006,
    /// `asin_fct`
    Asin = 1007,
    /// `acos_fct`
    Acos = 1008,
    /// `atan_fct`
    Atan = 1009,
    /// `sinh_fct`
    Sinh = 1010,
    /// `cosh_fct`
    Cosh = 1011,
    /// `tanh_fct`
    Tanh = 1012,
    /// `exp_fct`
    Exp = 1013,
    /// `log_fct`
    Log = 1014,
    /// `log10_fct`
    Log10 = 1015,
    /// `sqrt_fct`
    Sqrt = 1016,
    /// `abs_fct`
    Abs = 1017,
    /// `atan2_fct`
    Atan2 = 1018,
    /// `ceil_fct`
    Ceil = 1019,
    /// `floor_fct`
    Floor = 1020,
    /// `round_fct`
    Round = 1021,
    /// `min1_fct`
    Min1 = 1022,
    /// `min2_fct`
    Min2 = 1023,
    /// `max1_fct`
    Max1 = 1024,
    /// `max2_fct`
    Max2 = 1025,
    /// `near_fct`
    Near = 1026,
    /// `circle_fct`
    Circle = 1027,
    /// `box_fct`
    Box = 1028,
    /// `elps_fct`
    Ellipse = 1029,
    /// `isnull_fct`
    IsNull = 1030,
    /// `defnull_fct`
    DefNull = 1031,
    /// `gtifilt_fct`
    GtiFilt = 1032,
    /// `regfilt_fct`
    RegFilt = 1033,
    /// `ifthenelse_fct`
    IfThenElse = 1034,
    /// `row_fct`
    Row = 1035,
    /// `null_fct`
    Null = 1036,
    /// `median_fct`
    Median = 1037,
    /// `average_fct`
    Average = 1038,
    /// `stddev_fct`
    Stddev = 1039,
    /// `nonnull_fct`
    NonNull = 1040,
    /// `angsep_fct`
    AngSep = 1041,
    /// `gasrnd_fct`
    GasRnd = 1042,
    /// `poirnd_fct`
    PoiRnd = 1043,
    /// `strmid_fct`
    StrMid = 1044,
    /// `strpos_fct`
    StrPos = 1045,
    /// `setnull_fct`
    SetNull = 1046,
    /// `gtiover_fct`
    GtiOver = 1047,
    /// `gtifind_fct`
    GtiFind = 1048,
    /// `elemnum_fct`
    ElemNum = 1049,
    /// `axiselem_fct`
    AxisElem = 1050,
    /// `array_fct`
    Array = 1051,
}

impl FuncOp {
    /// The `Node::operation` encoding.
    pub(crate) fn code(self) -> c_int {
        self as c_int
    }

    /// Recover the function from a `Node::operation` value.
    ///
    /// `operation` also encodes column indices, `CONST_OP` and operator
    /// characters, so this returns `None` for anything that is not a function.
    pub(crate) fn from_code(v: c_int) -> Option<FuncOp> {
        Some(match v {
            1001 => FuncOp::Rnd,
            1002 => FuncOp::Sum,
            1003 => FuncOp::Nelem,
            1004 => FuncOp::Sin,
            1005 => FuncOp::Cos,
            1006 => FuncOp::Tan,
            1007 => FuncOp::Asin,
            1008 => FuncOp::Acos,
            1009 => FuncOp::Atan,
            1010 => FuncOp::Sinh,
            1011 => FuncOp::Cosh,
            1012 => FuncOp::Tanh,
            1013 => FuncOp::Exp,
            1014 => FuncOp::Log,
            1015 => FuncOp::Log10,
            1016 => FuncOp::Sqrt,
            1017 => FuncOp::Abs,
            1018 => FuncOp::Atan2,
            1019 => FuncOp::Ceil,
            1020 => FuncOp::Floor,
            1021 => FuncOp::Round,
            1022 => FuncOp::Min1,
            1023 => FuncOp::Min2,
            1024 => FuncOp::Max1,
            1025 => FuncOp::Max2,
            1026 => FuncOp::Near,
            1027 => FuncOp::Circle,
            1028 => FuncOp::Box,
            1029 => FuncOp::Ellipse,
            1030 => FuncOp::IsNull,
            1031 => FuncOp::DefNull,
            1032 => FuncOp::GtiFilt,
            1033 => FuncOp::RegFilt,
            1034 => FuncOp::IfThenElse,
            1035 => FuncOp::Row,
            1036 => FuncOp::Null,
            1037 => FuncOp::Median,
            1038 => FuncOp::Average,
            1039 => FuncOp::Stddev,
            1040 => FuncOp::NonNull,
            1041 => FuncOp::AngSep,
            1042 => FuncOp::GasRnd,
            1043 => FuncOp::PoiRnd,
            1044 => FuncOp::StrMid,
            1045 => FuncOp::StrPos,
            1046 => FuncOp::SetNull,
            1047 => FuncOp::GtiOver,
            1048 => FuncOp::GtiFind,
            1049 => FuncOp::ElemNum,
            1050 => FuncOp::AxisElem,
            1051 => FuncOp::Array,
            _ => return None,
        })
    }
}
/// What a computed node does, stored in [`Operation::Op`].
///
/// `eval.y` used the operator character for the arithmetic and bitwise cases
/// and a token number for the rest, and the transpilation left both as bare
/// integers -- so `Do_BinOp_lng` matched on `43` for `+` and `279` for `==`.
/// The discriminants are unchanged.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[repr(i32)]
pub(crate) enum OpCode {
    /* spelled as the character eval.y wrote */
    Mod = 37,
    BitAnd = 38,
    Mul = 42,
    Add = 43,
    Sub = 45,
    Div = 47,
    BitXor = 94,
    BitOr = 124,
    /// `~`, approximate equality
    Approx = 126,
    /// `[`, a subscript node built by `New_Deref`
    Deref = 91,
    /// `{`, a vector literal built by `New_Vector`
    Vector = 123,

    /* a conversion's opcode is its target sort */
    ToBoolean = 258,
    ToLong = 259,
    ToDouble = 260,

    Or = 277,
    And = 278,
    Eq = 279,
    Ne = 280,
    Gt = 281,
    Lt = 282,
    Lte = 283,
    Gte = 284,
    Power = 286,
    Not = 287,
    IntCast = 288,
    FltCast = 289,
    UMinus = 290,
    /// `ACCUM()`, a running total
    Accum = 291,
    /// `SEQDIFF()`, a running difference
    Diff = 292,
}

impl OpCode {
    pub(crate) fn code(self) -> c_int {
        self as c_int
    }

    pub(crate) fn from_code(v: c_int) -> Option<OpCode> {
        use OpCode::*;
        Some(match v {
            37 => Mod,
            38 => BitAnd,
            42 => Mul,
            43 => Add,
            45 => Sub,
            47 => Div,
            94 => BitXor,
            124 => BitOr,
            126 => Approx,
            258 => ToBoolean,
            259 => ToLong,
            260 => ToDouble,
            277 => Or,
            278 => And,
            279 => Eq,
            280 => Ne,
            281 => Gt,
            282 => Lt,
            283 => Lte,
            284 => Gte,
            286 => Power,
            287 => Not,
            288 => IntCast,
            289 => FltCast,
            290 => UMinus,
            291 => Accum,
            292 => Diff,
            _ => return None,
        })
    }
}

/// What a [`Node`] is.
///
/// The C packed all four cases into one `int`: a constant was the sentinel
/// -1000, a column was the *negated* index of its entry in `colData`, and
/// anything positive was either an operator character, a token number or a
/// `funcOp`. Tests like `operation > 0` meaning "has sub-nodes to evaluate"
/// only made sense with that layout in mind.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) enum Operation {
    /// A constant folded at parse time; was `CONST_OP`.
    #[default]
    Const,
    /// A leaf reading `colData[index]`; was a negated index.
    Column(c_int),
    /// A unary or binary operator.
    Op(OpCode),
    /// A function call.
    Func(FuncOp),
}

impl Operation {
    /// True for nodes with sub-nodes to evaluate first; was `operation > 0`.
    pub(crate) fn is_computed(self) -> bool {
        matches!(self, Operation::Op(_) | Operation::Func(_))
    }

    /// The column index, for a column leaf.
    pub(crate) fn column(self) -> Option<c_int> {
        match self {
            Operation::Column(i) => Some(i),
            _ => None,
        }
    }

    pub(crate) fn op(self) -> Option<OpCode> {
        match self {
            Operation::Op(o) => Some(o),
            _ => None,
        }
    }
}

/// The sort of a parser value: which of `eval.y`'s four nonterminals a node
/// belongs to.
///
/// The transpilation stored this as a bare `c_int` holding a lexer token id, so
/// every test read `node.ntype == <token> as c_int`.
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

    /// True for the sorts `eval.y` calls `expr`; `Boolean` is a `bexpr`.
    pub(crate) fn is_expr(self) -> bool {
        matches!(self, ValueSort::Long | ValueSort::Double)
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
    /// Not an element array at all: a single owned object the engine stashes
    /// in the node, such as the `SAORegion` built by `New_REG`. Only ever read
    /// back through [`NodeValue::raw`].
    Opaque,
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
/// The buffer is still a raw allocation rather than a `Vec`. The `Do_*`
/// routines read two operand buffers while writing a third, all indexed out of
/// `ParseData::Nodes`; owning the storage would make that a borrow conflict,
/// and untangling it is a separate change from tagging the union.
/// The `Text` arm makes this 264 bytes where the others need 16. That is the
/// union's own footprint, so it is not a regression, and boxing it is not an
/// option while `lval` and `Node` are `Copy` and get copied by value all
/// through the engine.
#[allow(clippy::large_enum_variant)]
#[derive(Copy, Clone, Default)]
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
    },
}

impl core::fmt::Debug for NodeValue {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            NodeValue::Empty => write!(f, "Empty"),
            NodeValue::Long(v) => write!(f, "Long({v})"),
            NodeValue::Double(v) => write!(f, "Double({v})"),
            NodeValue::Logical(v) => write!(f, "Logical({v})"),
            NodeValue::Text(_) => write!(f, "Text(..)"),
            NodeValue::Buffer { kind, ptr } => write!(f, "Buffer({kind:?}, {ptr:?})"),
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

    /// A writable pointer to the scalar string buffer.
    ///
    /// The C wrote into `data.str` of a node that had no value yet, so this
    /// installs an empty `Text` first when the node holds something else.
    pub(crate) fn text_mut_ptr(&mut self) -> *mut c_char {
        if !matches!(self, NodeValue::Text(_)) {
            *self = NodeValue::Text([0; MAX_STRLEN as usize]);
        }
        match self {
            NodeValue::Text(v) => v.as_mut_ptr(),
            _ => unreachable!(),
        }
    }

    fn buffer(&self, want: BufferKind) -> *mut c_void {
        match self {
            NodeValue::Buffer { kind, ptr } if *kind == want => *ptr,
            /* an unallocated node reads as a null buffer, as the union did */
            NodeValue::Empty => core::ptr::null_mut(),
            other => wrong_arm("Buffer", other),
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

    pub(crate) fn set_buffer(&mut self, kind: BufferKind, ptr: *mut c_void) {
        *self = NodeValue::Buffer { kind, ptr };
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
            NodeValue::Long(v) => (v as *const c_long).cast(),
            NodeValue::Double(v) => (v as *const f64).cast(),
            NodeValue::Logical(v) => (v as *const c_char).cast(),
            NodeValue::Text(v) => v.as_ptr().cast(),
            NodeValue::Empty | NodeValue::Buffer { .. } => core::ptr::null(),
        }
    }

    /// Release a row buffer and mark the node empty.
    ///
    /// # Safety
    /// The buffer must have come from `malloc`/`calloc` and must not be
    /// referenced elsewhere.
    pub(crate) unsafe fn free_buffer(&mut self) {
        if let NodeValue::Buffer { ptr, .. } = *self
            && !ptr.is_null()
        {
            unsafe { libc::free(ptr) };
        }
        *self = NodeValue::Empty;
    }
}

#[cold]
#[track_caller]
fn wrong_arm(want: &str, got: &NodeValue) -> ! {
    panic!("node value is {got:?}, expected {want}");
}

#[derive(Default, Debug, Copy, Clone)]
pub(crate) struct lval {
    pub nelem: c_long,
    pub naxis: c_int,
    pub naxes: [c_long; MAXDIMS as usize],
    pub undef: *mut c_char,
    pub data: NodeValue,
}

#[derive(Default, Debug, Copy, Clone)]
pub(crate) struct Node {
    pub operation: Operation,
    pub DoOp: Option<fn(p: &mut ParseData, this_node_idx: usize)>,
    pub nSubNodes: c_int,
    pub SubNodes: [usize; MAXSUBS as usize],
    pub ntype: ValueSort,
    /// Only meaningful on the START/STOP node `New_GTI` builds: whether the
    /// GTI rows are fully time-ordered, which `GTIOVERLAP` requires and
    /// `Search_GTI` uses to pick a binary search over a linear one.
    ///
    /// The C stashed this in `ntype`, which every other node uses for its
    /// value sort.
    pub gti_ordered: bool,
    pub value: lval,
}

/// The sort of a table column, as the parser sees it.
///
/// These were the `COLUMN` / `BCOLUMN` / `SCOLUMN` / `BITCOL` lexer tokens.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum ColumnSort {
    /// Integer or floating point.
    Numeric,
    Boolean,
    String,
    Bits,
}

/// What a name in an expression resolves to.
///
/// This replaces bison's `%union YYSTYPE` together with the parallel `c_int`
/// token kind that [`GetDataFn`] used to return beside it. Two things follow
/// from making it an enum: the variant *is* the kind, so the value and the tag
/// can no longer disagree, and reading it needs no `unsafe`.
///
/// The union also had a `Node` arm for the parser's own stack, which no longer
/// exists, and could in principle hold a bit-string constant, which no resolver
/// ever produced. Neither is representable here.
#[derive(Clone, Debug, PartialEq)]
pub(crate) enum ParserValue {
    /// A table column, by index into [`ParseData::varData`].
    Column { index: c_int, sort: ColumnSort },
    /// A header keyword that resolved to an integer.
    Long(c_long),
    /// A header keyword that resolved to a float.
    Double(c_double),
    /// A header keyword that resolved to `T` or `F`.
    Boolean(bool),
    /// A header keyword that resolved to a quoted string, NUL-terminated.
    Str(Vec<c_char>),
}

/// Fetches a single named value for the parser (ffcalc's keyword lookup).
///
/// Returns `None` when the name cannot be resolved, having already recorded a
/// status in [`ParseData::status`] and a message on the error stack.
pub(crate) type GetDataFn = fn(p: &mut ParseData, dataName: &[c_char]) -> Option<ParserValue>;

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
    /// The expression lowered for the new columnar evaluator, when it could
    /// handle every construct in it. `None` means the `Node` arena evaluates
    /// this expression, which is still the case for most of them.
    #[cfg(feature = "new-eval")]
    pub expr_tree: Option<crate::eval::expr::Expr>,
    /// The running values of the expression's `ACCUM`/`SEQDIFF` nodes, which
    /// have to survive from one batch to the next.
    pub accum_state: Vec<crate::eval::expr::AccumState>,
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
        #[cfg(feature = "new-eval")]
        {
            self.expr_tree = None;
            self.accum_state = Vec::new();
        }
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
