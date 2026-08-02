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
use crate::c_types::c_char;
use crate::eval_defs::{ParseData, ParserValue};
use crate::eval_y::fits_parser_yyGetVariable;
use crate::fitsio::PARSE_SYNTAX_ERR;

/// What each name resolved to, keyed by the name token's byte offset (unique,
/// since two tokens cannot start at the same place).
pub(crate) type Resolutions = HashMap<usize, ParserValue>;

fn cstr(b: &[u8]) -> Vec<c_char> {
    let mut v: Vec<c_char> = b.iter().map(|&c| c as c_char).collect();
    v.push(0);
    v
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
        let resolved = if bare {
            /* `{variable}`: registered columns first, then getData */
            fits_parser_yyGetVariable(lParse, &cname)
        } else if let Some(get) = lParse.getData {
            /* `{constant}`: straight to the keyword lookup */
            get(lParse, &cname)
        } else {
            if lParse.status == 0 {
                lParse.status = PARSE_SYNTAX_ERR;
            }
            None
        };

        let Some(value) = resolved else {
            /* the resolver has set lParse.status; the parse stops here */
            toks.truncate(i);
            break;
        };
        out.insert(toks[i].at, value);
    }

    out
}
