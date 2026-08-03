//! Hand-written front end for the CFITSIO expression language.
//!
//! Replaces the flex scanner and the bison LALR driver that used to live in
//! `src/eval_l.rs` and the first half of `src/eval_y.rs`. The pipeline is:
//!
//! ```text
//! bytes -> lexer::tokenize -> [Spanned]  (nom patterns + flex's longest-match rule)
//!       -> grammar::parse  -> Ast        (Pratt / precedence climbing, untyped)
//!       -> lower::lower    -> Node arena (identifier resolution, sorts, promotions)
//! ```
//!
//! The split matters: `eval.y` fused syntax, sort-checking and node building
//! into one LALR machine with four mutually-recursive nonterminals
//! (`expr`/`bexpr`/`sexpr`/`bits`). Here the grammar is untyped and every
//! type-directed decision happens in `lower`, which is why 135 productions
//! collapse to about 30 AST variants. See `PARSER_SPEC.md` and
//! `PARSER_MIGRATION.md`.

pub(crate) mod ast;
pub(crate) mod error;
pub(crate) mod grammar;
pub(crate) mod lexer;
pub(crate) mod lower;
pub(crate) mod resolve;
pub(crate) mod token;

use crate::c_types::c_int;
#[cfg(feature = "new-eval")]
use crate::c_types::c_long;
use crate::eval_defs::ParseData;
#[cfg(feature = "new-eval")]
use crate::eval_defs::ValueSort;
use crate::fitsio::PARSE_SYNTAX_ERR;

/// Parse `lParse.expr`, filling `lParse.Nodes` and setting `lParse.resultNode`.
///
/// Errors are reported through `lParse.status` and the CFITSIO error stack; the
/// return value is always 0, and `ffiprs` surfaces `lParse.status`.
///
/// That indirection matters for compatibility. When the flex scanner failed to
/// resolve a name, `find_column` set `lParse.status` (to `COL_NOT_FOUND`, say)
/// and returned `pERROR`; bison treats any token `<= 0` as end-of-input, so
/// `yyparse` reduced the empty `lines` rule and returned *success*, and
/// `ffiprs` then reported the resolver's status rather than a syntax error.
/// So `NOSUCHCOL` is a 202, not a 431, and the GTI/region builders likewise
/// surface their own file-open failures.
///
/// A blank expression is not an error here — it produces no nodes, and
/// `ffiprs` reports "Blank expression" from `nNodes == 0`, exactly as before.
pub(crate) fn parse_expression(lParse: &mut ParseData) -> c_int {
    let src = match lParse.expr.as_deref() {
        Some(s) => s.to_vec(),
        None => return 0,
    };

    /// Report a failure. A status already set by `getData` or by one of the
    /// node builders wins; only a pure syntax error becomes `PARSE_SYNTAX_ERR`.
    fn fail(lParse: &mut ParseData, e: error::ParseError) -> c_int {
        if lParse.status == 0 {
            e.report();
            lParse.status = PARSE_SYNTAX_ERR;
        }
        0
    }

    let mut toks = match lexer::tokenize(&src) {
        Ok(t) => t,
        Err(e) => {
            fail(lParse, e);
            return PARSE_SYNTAX_ERR;
        }
    };

    /* resolve names in source order, truncating at the first failure */
    let names = resolve::resolve_names(lParse, &mut toks);

    let tree = match grammar::parse(&toks) {
        Ok(t) => t,
        Err(e) => {
            /* a syntax error outranks whatever the resolver recorded, which is
            how `1 + NOSUCHCOL` reports PARSE_SYNTAX_ERR and a bare
            `NOSUCHCOL` reports COL_NOT_FOUND */
            fail(lParse, e);
            return PARSE_SYNTAX_ERR;
        }
    };
    let Some(tree) = tree else {
        /* nothing but blank lines; ffiprs turns nNodes == 0 into an error */
        return 0;
    };

    let mut lowerer = lower::Lowerer {
        p: lParse,
        names: &names,
    };
    let status = match lowerer.lower(&tree) {
        Ok(node) => {
            lParse.resultNode = node;
            0
        }
        Err(e) => fail(lParse, e),
    };

    /* The arena is still built either way: `ffiprs` reads the result node for
    the expression's datatype and shape, and the new evaluator only replaces
    the per-row computation. */
    #[cfg(feature = "new-eval")]
    if status == 0 && lParse.status == 0 {
        /* the per-column shapes let the lowering decide whether a subscript
        names a single element or a slice */
        let shapes: Vec<(c_long, Vec<c_long>)> = lParse
            .varData
            .iter()
            .map(|v| {
                (
                    v.nelem,
                    v.naxes[..(v.naxis.max(0) as usize).min(v.naxes.len())]
                        .iter()
                        .map(|&n| n.max(0))
                        .collect(),
                )
            })
            .collect();
        let cols = crate::eval::lower::Columns {
            shapes,
            sorts: lParse.varData.iter().map(|v| v.dtype).collect(),
        };
        lParse.expr_tree = crate::eval::lower::lower(&tree, &names, &cols)
            .ok()
            /* A bit-string *result* is never retrievable -- the engine reports
            432 for a row-varying one and 433 for a constant -- so leave those
            with the arena rather than reproduce the error. Bit-valued
            subexpressions are fine; only the top-level sort matters. */
            .filter(|t| {
                t.sort(&|i| cols.sorts.get(i).copied().unwrap_or(ValueSort::Long))
                    != ValueSort::Bits
            });
    }

    status
}
