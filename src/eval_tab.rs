//! Numeric tags shared between the parser and the evaluation engine.
//!
//! This was `eval_tab.h`, bison's token-number header. The parser no longer
//! produces tokens — `crate::parser::token::Tok` does, with no numbering — so
//! what survives here are the two roles those numbers still play inside
//! `Node`:
//!
//! * **sort tags** (`BOOLEAN` … `BITSTR`), stored in `Node::ntype` and
//!   `DataInfo::dtype`. The first three are ordered so that numeric promotion
//!   is a `>` comparison, which is why the values matter and must not be
//!   renumbered.
//! * **operator codes** (`OR` … `DIFF`), stored in `Node::operation` and
//!   dispatched on by the `Do_*` routines.
//!
//! The lexer-token variants that used to sit between them — `FUNCTION`,
//! `GTIFILTER`, `COLUMN`, `ROWREF` and the rest — are gone: function
//! classification is now `parser::token::CallKind`, and a resolved column is
//! `eval_defs::ColumnSort`. Their numbers are left as gaps rather than reused,
//! so the surviving values still line up with upstream's `eval_tab.h`.

use crate::c_types::c_int;

#[repr(i32)]
pub(crate) enum fits_parser_yytokentype {
    /* value sorts -- ordered for promotion */
    BOOLEAN = 258,
    LONG = 259,
    DOUBLE = 260,
    STRING = 261,
    BITSTR = 262,

    /* operator codes stored in Node::operation */
    OR = 277,
    AND = 278,
    EQ = 279,
    NE = 280,
    GT = 281,
    LT = 282,
    LTE = 283,
    GTE = 284,
    XOR = 285,
    POWER = 286,
    NOT = 287,
    INTCAST = 288,
    FLTCAST = 289,
    UMINUS = 290,
    ACCUM = 291,
    DIFF = 292,
}

impl From<c_int> for fits_parser_yytokentype {
    /// Callers convert both `Node::ntype` / `DataInfo::dtype` (a sort) and
    /// `Node::operation` (an operator code), so every surviving variant has to
    /// round-trip. `Node::operation` also holds raw operator characters and
    /// `funcOp` values, which are matched before this is reached.
    fn from(value: c_int) -> Self {
        match value {
            258 => Self::BOOLEAN,
            259 => Self::LONG,
            260 => Self::DOUBLE,
            261 => Self::STRING,
            262 => Self::BITSTR,
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
            _ => panic!("not a parser sort or operator code: {value}"),
        }
    }
}
