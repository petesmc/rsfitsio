//! Token type produced by [`super::lexer`].
//!
//! One variant per `eval.l` rule (see `PARSER_SPEC.md` §2.4). Unlike the flex
//! lexer, identifiers are *not* resolved here: `Ident` and `Keyword` carry the
//! raw name and are looked up during lowering, so the tokenizer is a pure
//! function of the input bytes.

use crate::c_types::c_long;

/// Which of the six function-token classes the lexer assigned to a name.
///
/// The classification is by name alone (`PARSER_SPEC.md` §4); which concrete
/// function is meant is decided during lowering, from the argument sorts.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum CallKind {
    /// Everything not named below.
    Function,
    /// `BOX(`, `CIRCLE(`, `ELLIPSE(`, `NEAR(`, `ISNULL(` — always boolean.
    BFunction,
    /// `STRSTR(` — returns an integer.
    IFunction,
    GtiFilter,
    GtiOverlap,
    GtiFind,
    RegFilter,
}

#[derive(Clone, Debug, PartialEq)]
pub(crate) enum Tok {
    /* ---- literals ---- */
    Long(c_long),
    Double(f64),
    Boolean(bool),
    /// String literal with the surrounding quotes removed.
    Str(Vec<u8>),
    /// Bit string, already expanded to ASCII `'0'` / `'1'` / `'x'`.
    BitStr(Vec<u8>),

    /* ---- names, resolved during lowering ---- */
    /// A column or variable name (`$…$` delimiters already stripped).
    Ident(Vec<u8>),
    /// A `#KEYWORD` header reference, *without* the leading `#`.
    Keyword(Vec<u8>),
    /// `#ROW`
    RowRef,
    /// `#NULL`
    NullRef,
    /// `#SNULL`
    SNullRef,
    /// A function name plus its opening parenthesis. The name is upper-cased
    /// and does not include the `(`.
    Call(CallKind, Vec<u8>),

    /* ---- multi-character operators ---- */
    /// `^` or `**`
    Power,
    /// `!`, `.not.`, `.NOT.`, `not.`, `NOT.`
    Not,
    /// `||`, `.or.`, …
    Or,
    /// `&&`, `.and.`, …
    And,
    /// `==`, `.eq.`, …
    Eq,
    /// `!=`, `.ne.`, …
    Ne,
    /// `>`, `.gt.`, …
    Gt,
    /// `<`, `.lt.`, …
    Lt,
    /// `>=`, `=>`, `.ge.`, …
    Gte,
    /// `<=`, `=<`, `.le.`, …
    Lte,
    /// `^^`, `.xor.`, `.XOR.`
    Xor,
    /// `(int)` / `(INT)`
    IntCast,
    /// `(float)` / `(FLOAT)` / `(double)` / `(DOUBLE)`
    FltCast,

    /* ---- single characters and the terminator ---- */
    /// Any other single byte: `( ) [ ] { } , : ? + - * / % & | ~ =` and junk.
    Char(u8),
    /// The `\n` that terminates every expression.
    Newline,
}

impl Tok {
    /// The single byte this token stands for, if it is a bare character.
    pub(crate) fn as_char(&self) -> Option<u8> {
        match self {
            Tok::Char(c) => Some(*c),
            _ => None,
        }
    }

    /// True when this token is the bare character `c`.
    pub(crate) fn is_char(&self, c: u8) -> bool {
        self.as_char() == Some(c)
    }
}

/// A token plus the byte offset in the source expression at which it started.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct Spanned {
    pub(crate) tok: Tok,
    pub(crate) at: usize,
}
