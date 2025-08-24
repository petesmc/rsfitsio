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

#[derive(Copy, Clone)]
pub(crate) union FITS_PARSER_YYSTYPE {
    pub(crate) Node: c_int,                         /* Index of Node */
    pub(crate) dbl: c_double,                       /* real value    */
    pub(crate) lng: c_long,                         /* integer value */
    pub(crate) log: c_char,                         /* logical value */
    pub(crate) astr: [c_char; MAX_STRLEN as usize], /* string value  */
}
