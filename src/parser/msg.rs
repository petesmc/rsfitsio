//! What CFITSIO says when it turns an expression down.
//!
//! The wording of a rejection is part of the library's behaviour, not an
//! implementation detail: a caller can read the error stack with
//! `fits_read_errmsg` and compare it against the C library's. So these are
//! transcribed from `eval.y` rather than written, and checked against the real
//! thing with `ORACLE_MSGS=1 tests/oracle/oracle` -- see PARSER_MIGRATION.md
//! §10.6.
//!
//! They live here, beside the grammar they describe, rather than in
//! `eval::lower` which raises them.

use super::error::ParseError;
use crate::eval_defs::ValueSort;

/// True for the sorts `eval.y` calls `expr` -- the numbers. A boolean is a
/// `bexpr`, a different nonterminal, and most productions take only one.
fn is_expr(t: ValueSort) -> bool {
    t.is_expr()
}

/// A function called with argument sorts it has no form for.
///
/// CFITSIO's grammar is type-stratified -- `expr`, `bexpr`, `sexpr` and `bits`
/// are separate nonterminals -- so a call is a production over a *fixed* list
/// of argument sorts, and the message is a literal attached to whichever
/// production reduced. Two consequences, both of which the wording here has to
/// follow:
///
///   * The message names the production, not the arguments as written. That is
///     why `MIN(BOOLCOL,BOOLCOL)` reports `Boolean Function(expr,expr)` and not
///     anything mentioning `bool` -- it reduced by `bexpr: FUNCTION bexpr ','
///     bexpr ')'`, whose literal says `expr,expr`.
///
///   * A shape with no production at all never reaches an action, so bison
///     fails first and the message is a bare `syntax error`. `MIN(1,2,3)` is
///     one of these: there is no three-argument numeric production, though
///     there is a four-argument one.
///
/// Which layer a call belongs to is the lexer's doing, not the grammar's --
/// `BOX`/`CIRCLE`/`ELLIPSE`/`NEAR`/`ISNULL` come back as `BFUNCTION` and get
/// the `Boolean Function...` wordings, everything else as `FUNCTION`. That is
/// why `NEAR(1)` and `NOSUCHFUNC(1)` differ despite the identical shape.
///
/// The name itself only decides whether the production's action *accepts* the
/// call; when it does not, the production's literal is what reaches the stack.
/// Derived from `eval.y` and checked against the library with
/// `ORACLE_MSGS=1 tests/oracle/oracle`.
pub(crate) fn unsupported(name: &[u8], args: &[ValueSort]) -> ParseError {
    /* no production matched, so bison reported the error itself. The offset is
    not part of a rendered message -- see `ParseError::report` -- so nothing is
    lost by not threading one down to here. */
    let syntax = || ParseError::syntax("syntax error", 0, b"");
    let semantic = |m: &str| ParseError::semantic(m);

    let num = |t: &ValueSort| is_expr(*t);
    let all_num = args.iter().all(num);

    if matches!(name, b"BOX" | b"CIRCLE" | b"ELLIPSE" | b"NEAR" | b"ISNULL") {
        return match args.len() {
            /* `BFUNCTION expr ')'`, `BFUNCTION bexpr ')'` and `BFUNCTION
            sexpr ')'` all exist and all say the same thing; there is no
            `BFUNCTION bits ')'`, so ISNULL(BITS) is a syntax error */
            1 if args[0] != ValueSort::Bits => semantic("Boolean Function(expr) not supported"),
            3 | 5 if all_num => semantic("Boolean Function not supported"),
            7 if all_num => semantic("SAO Image Function not supported"),
            _ => syntax(),
        };
    }
    /* STRSTR is the only IFUNCTION and its one production is always valid, so
    any shape reaching here matched nothing. */
    if name == b"STRSTR" {
        return syntax();
    }

    match args {
        [] => semantic("Function() not supported"),
        [ValueSort::Boolean] => semantic("Function(bool) not supported"),
        [ValueSort::String] => semantic("Function(str) not supported"),
        [ValueSort::Bits] => semantic("Function(bits) not supported"),
        [t] if num(t) => semantic("Function(expr) not supported"),
        [ValueSort::Boolean, b] if num(b) => semantic("Function(bool,expr) not supported"),
        [ValueSort::Boolean, ValueSort::Boolean] => {
            semantic("Boolean Function(expr,expr) not supported")
        }
        [ValueSort::String, ValueSort::String] => semantic("Function(string,string) not supported"),
        [a, b] if num(a) && num(b) => semantic("Function(expr,expr) not supported"),
        [ValueSort::String, b, c] if num(b) && num(c) => {
            semantic("Function(string,expr,expr) not supported")
        }
        _ if args.len() == 4 && all_num => semantic("Function(expr,expr,expr,expr) not supported"),
        _ => syntax(),
    }
}
