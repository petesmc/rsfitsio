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
    BOOLEAN = 258,             /* BOOLEAN  */
    LONG = 259,                /* LONG  */
    DOUBLE = 260,              /* DOUBLE  */
    STRING = 261,              /* STRING  */
    BITSTR = 262,              /* BITSTR  */
    FUNCTION = 263,            /* FUNCTION  */
    BFUNCTION = 264,           /* BFUNCTION  */
    IFUNCTION = 265,           /* IFUNCTION  */
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

#[derive(Copy, Clone)]
pub(crate) union FITS_PARSER_YYSTYPE {
    pub(crate) Node: c_int,                         /* Index of Node */
    pub(crate) dbl: c_double,                       /* real value    */
    pub(crate) lng: c_long,                         /* integer value */
    pub(crate) log: c_char,                         /* logical value */
    pub(crate) astr: [c_char; MAX_STRLEN as usize], /* string value  */
}
