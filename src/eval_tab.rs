use crate::{
    c_types::{c_char, c_double, c_int, c_long},
    eval_defs::MAX_STRLEN,
};

#[repr(i32)]
pub(crate) enum fits_parser_yytokentype {
    FITS_PARSER_YYEMPTY = -2,
    FITS_PARSER_YYEOF = 0,     /* "end of file"  */
    FITS_PARSER_YYerror = 256, /* error  */
    FITS_PARSER_YYUNDEF = 257, /* "invalid token"  */
    BOOLEAN = 258,             /* First 3 must be in order of        */
    LONG = 259,                /* increasing promotion for later use */
    DOUBLE = 260,              /* DOUBLE  */
    STRING = 261,              /* STRING  */
    BITSTR = 262,              /* BITSTR  */
    FUNCTION = 263,            /* FUNCTION  */
    BFUNCTION = 264,           /* Bit function */
    IFUNCTION = 265,           /* Integer function */
    GTIFILTER = 266,           /* GTIFILTER  */
    GTIOVERLAP = 267,          /* GTIOVERLAP  */
    GTIFIND = 268,             /* GTIFIND  */
    REGFILTER = 269,           /* REGFILTER  */
    COLUMN = 270,              /* COLUMN  */
    BCOLUMN = 271,             /* BCOLUMN  */
    SCOLUMN = 272,             /* SCOLUMN  */
    BITCOL = 273,              /* BITCOL  */
    ROWREF = 274,              /* ROWREF  */
    NULLREF = 275,             /* NULLREF  */
    SNULLREF = 276,            /* SNULLREF  */
    OR = 277,                  /* OR  */
    AND = 278,                 /* AND  */
    EQ = 279,                  /* EQ  */
    NE = 280,                  /* NE  */
    GT = 281,                  /* GT  */
    LT = 282,                  /* LT  */
    LTE = 283,                 /* LTE  */
    GTE = 284,                 /* GTE  */
    XOR = 285,                 /* XOR  */
    POWER = 286,               /* POWER  */
    NOT = 287,                 /* NOT  */
    INTCAST = 288,             /* INTCAST  */
    FLTCAST = 289,             /* FLTCAST  */
    UMINUS = 290,              /* UMINUS  */
    ACCUM = 291,               /* ACCUM  */
    DIFF = 292,                /* DIFF  */
}

impl From<c_int> for fits_parser_yytokentype {
    fn from(value: c_int) -> Self {
        match value {
            -2 => Self::FITS_PARSER_YYEMPTY,
            0 => Self::FITS_PARSER_YYEOF,
            256 => Self::FITS_PARSER_YYerror,
            257 => Self::FITS_PARSER_YYUNDEF,
            258 => Self::BOOLEAN,
            259 => Self::LONG,
            260 => Self::DOUBLE,
            261 => Self::STRING,
            262 => Self::BITSTR,
            263 => Self::FUNCTION,
            264 => Self::BFUNCTION,
            265 => Self::IFUNCTION,
            266 => Self::GTIFILTER,
            267 => Self::GTIOVERLAP,
            268 => Self::GTIFIND,
            269 => Self::REGFILTER,
            270 => Self::COLUMN,
            271 => Self::BCOLUMN,
            272 => Self::SCOLUMN,
            273 => Self::BITCOL,
            274 => Self::ROWREF,
            275 => Self::NULLREF,
            276 => Self::SNULLREF,
            277 => Self::OR,
            278 => Self::AND,
            279 => Self::EQ,
            280 => Self::NE,
            281 => Self::GT,
            282 => Self::LT,
            283 => Self::LTE,
            284 => Self::GTE,
            285 => Self::XOR,
            286 => Self::POWER,
            287 => Self::NOT,
            288 => Self::INTCAST,
            289 => Self::FLTCAST,
            290 => Self::UMINUS,
            291 => Self::ACCUM,
            292 => Self::DIFF,
            _ => panic!(),
        }
    }
}

/// The semantic value bison carries on its value stack.
///
/// This was `%union YYSTYPE`. Making it an enum means the variant *is* the
/// kind, so a slot written as one arm cannot be read back as another without
/// saying so -- the same guarantee [`NodeValue`](crate::eval_defs::NodeValue)
/// gives node values.
///
/// It stays `Copy`: the stack is a `[FITS_PARSER_YYSTYPE; YYINITDEPTH]` array
/// that bison assigns through by value on every shift and reduce.
#[allow(clippy::large_enum_variant)]
#[derive(Copy, Clone, Default)]
pub(crate) enum FITS_PARSER_YYSTYPE {
    /// A stack slot bison has not written yet.
    #[default]
    Empty,
    Node(c_int),                         /* Index of Node */
    Double(c_double),                    /* real value    */
    Long(c_long),                        /* integer value */
    Logical(c_char),                     /* logical value */
    Text([c_char; MAX_STRLEN as usize]), /* string value  */
}

impl core::fmt::Debug for FITS_PARSER_YYSTYPE {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            FITS_PARSER_YYSTYPE::Empty => write!(f, "Empty"),
            FITS_PARSER_YYSTYPE::Node(v) => write!(f, "Node({v})"),
            FITS_PARSER_YYSTYPE::Double(v) => write!(f, "Double({v})"),
            FITS_PARSER_YYSTYPE::Long(v) => write!(f, "Long({v})"),
            FITS_PARSER_YYSTYPE::Logical(v) => write!(f, "Logical({v})"),
            FITS_PARSER_YYSTYPE::Text(_) => write!(f, "Text(..)"),
        }
    }
}

impl FITS_PARSER_YYSTYPE {
    /// The node index this slot carries.
    ///
    /// An unwritten slot reads as node 0, which is what the union's zeroed
    /// stack gave; bison's `$$ = $1` default copies slots that no action has
    /// filled in, and the grammar never uses such a value.
    pub(crate) fn node(&self) -> c_int {
        match self {
            FITS_PARSER_YYSTYPE::Node(v) => *v,
            FITS_PARSER_YYSTYPE::Empty => 0,
            other => wrong_arm("Node", other),
        }
    }

    pub(crate) fn dbl(&self) -> c_double {
        match self {
            FITS_PARSER_YYSTYPE::Double(v) => *v,
            other => wrong_arm("Double", other),
        }
    }

    pub(crate) fn lng(&self) -> c_long {
        match self {
            FITS_PARSER_YYSTYPE::Long(v) => *v,
            other => wrong_arm("Long", other),
        }
    }

    pub(crate) fn log(&self) -> c_char {
        match self {
            FITS_PARSER_YYSTYPE::Logical(v) => *v,
            other => wrong_arm("Logical", other),
        }
    }

    /// The string value, as the fixed-size array the lexer and grammar index.
    pub(crate) fn text(&self) -> &[c_char; MAX_STRLEN as usize] {
        match self {
            FITS_PARSER_YYSTYPE::Text(v) => v,
            other => wrong_arm("Text", other),
        }
    }

    /// A writable view of the string buffer for the `*_safe` string helpers,
    /// installing an empty one when the slot holds something else -- the lexer
    /// writes into `yylval->str` of a slot whose previous contents are
    /// irrelevant.
    pub(crate) fn text_mut(&mut self) -> &mut [c_char; MAX_STRLEN as usize] {
        if !matches!(self, FITS_PARSER_YYSTYPE::Text(_)) {
            *self = FITS_PARSER_YYSTYPE::Text([0; MAX_STRLEN as usize]);
        }
        match self {
            FITS_PARSER_YYSTYPE::Text(v) => v,
            _ => unreachable!(),
        }
    }

    /// A writable pointer to the string buffer, installing an empty one when
    /// the slot holds something else -- the lexer writes into `yylval->str` of
    /// a slot whose previous contents are irrelevant.
    pub(crate) fn text_mut_ptr(&mut self) -> *mut c_char {
        self.text_mut().as_mut_ptr()
    }
}

#[cold]
#[track_caller]
fn wrong_arm(want: &str, got: &FITS_PARSER_YYSTYPE) -> ! {
    panic!("parser value is {got:?}, expected {want}");
}
