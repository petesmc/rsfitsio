//! Name resolution pre-pass.
//!
//! The flex scanner resolved every `{variable}` and `{constant}` token *as it
//! lexed it*, by calling back into `ParseData::getData`. That has two
//! observable consequences the new front end has to reproduce:
//!
//! 1. **Registration order.** Columns are appended to `varData` / `colData` in
//!    the order their names appear in the text, and `ffiprs` reports them in
//!    that order.
//!
//! 2. **Truncation on failure.** A failed lookup returned `pERROR` (-1). Bison
//!    treats any token `<= YYEOF` as end-of-input, so the parse silently
//!    stopped at that point. `NOSUCHCOL` therefore parsed as the *empty*
//!    expression and `ffiprs` reported the resolver's status (`COL_NOT_FOUND`),
//!    while `1 + NOSUCHCOL` truncated to `1 +`, which is a syntax error, and
//!    reported `PARSE_SYNTAX_ERR`.
//!
//! Doing this as a separate pass over the token stream keeps the tokenizer a
//! pure function of its bytes, while resolving exactly the tokens flex would
//! have resolved, in the same order.

use std::collections::HashMap;

use super::token::{Spanned, Tok};
use crate::c_types::{c_char, c_int};
use crate::eval_tab::{FITS_PARSER_YYSTYPE, fits_parser_yytokentype as T};
use crate::eval_defs::ParseData;
use crate::eval_y::fits_parser_yyGetVariable;
use crate::fitsio::PARSE_SYNTAX_ERR;

/// What a name resolved to: the token kind the lexer would have returned, plus
/// its semantic value. Keyed by the name token's byte offset, which is unique.
pub(crate) type Resolutions = HashMap<usize, (c_int, FITS_PARSER_YYSTYPE)>;

fn cstr(b: &[u8]) -> Vec<c_char> {
    let mut v: Vec<c_char> = b.iter().map(|&c| c as c_char).collect();
    v.push(0);
    v
}

/// Is `dtype` one of the token kinds a name may resolve to?
pub(crate) fn token_of(v: c_int) -> Option<T> {
    match v {
        x if x == T::BOOLEAN as c_int => Some(T::BOOLEAN),
        x if x == T::LONG as c_int => Some(T::LONG),
        x if x == T::DOUBLE as c_int => Some(T::DOUBLE),
        x if x == T::STRING as c_int => Some(T::STRING),
        x if x == T::BITSTR as c_int => Some(T::BITSTR),
        x if x == T::COLUMN as c_int => Some(T::COLUMN),
        x if x == T::BCOLUMN as c_int => Some(T::BCOLUMN),
        x if x == T::SCOLUMN as c_int => Some(T::SCOLUMN),
        x if x == T::BITCOL as c_int => Some(T::BITCOL),
        _ => None,
    }
}

/// Resolve every name token, left to right.
///
/// On the first failure the token stream is truncated there, mirroring bison's
/// treatment of `pERROR` as end-of-input. The resolver has already recorded its
/// own status in `lParse.status`.
pub(crate) fn resolve_names(lParse: &mut ParseData, toks: &mut Vec<Spanned>) -> Resolutions {
    let mut out = Resolutions::new();

    for i in 0..toks.len() {
        let (name, bare) = match &toks[i].tok {
            Tok::Ident(n) => (n.clone(), true),
            Tok::Keyword(n) => {
                let mut full = vec![b'#'];
                full.extend_from_slice(n);
                (full, false)
            }
            _ => continue,
        };

        let cname = cstr(&name);
        let mut lval = FITS_PARSER_YYSTYPE { Node: 0 };
        let dtype = if bare {
            /* `{variable}`: registered columns first, then getData */
            fits_parser_yyGetVariable(lParse, &cname, &mut lval)
        } else if let Some(get) = lParse.getData {
            /* `{constant}`: straight to the keyword lookup */
            get(lParse, &cname, &mut lval)
        } else {
            if lParse.status == 0 {
                lParse.status = PARSE_SYNTAX_ERR;
            }
            -1
        };

        if token_of(dtype).is_none() {
            /* the resolver has set lParse.status; the parse stops here */
            toks.truncate(i);
            break;
        }
        out.insert(toks[i].at, (dtype, lval));
    }

    out
}
