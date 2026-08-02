//! Parse diagnostics.
//!
//! The flex/bison front end reported nearly everything as a bare
//! `"syntax error"` plus whatever `ffpmsg` calls the semantic actions happened
//! to make. This carries the same information but keeps the byte offset, so the
//! message can point at the offending token.

use crate::fitscore::ffpmsg_str;

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) enum ParseErrorKind {
    /// A literal exceeded `MAX_STRLEN`. The payload names the kind of literal,
    /// matching the wording of the C messages ("Bit string", "Variable", …).
    TooLong(String),
    /// The token stream did not match the grammar.
    Syntax(String),
    /// The expression parsed but is not well typed, or violates one of the
    /// structural constraints in `PARSER_SPEC.md` §3.4.
    Semantic(String),
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct ParseError {
    pub(crate) kind: ParseErrorKind,
    /// Byte offset into the expression at which the problem was detected.
    pub(crate) at: usize,
    /// A short excerpt of the offending text, for the message.
    pub(crate) text: Vec<u8>,
}

impl ParseError {
    pub(crate) fn new(kind: ParseErrorKind, at: usize, text: &[u8]) -> Self {
        ParseError {
            kind,
            at,
            text: text.iter().copied().take(20).collect(),
        }
    }

    /// A syntax error at `at`, quoting `text`.
    pub(crate) fn syntax(msg: impl Into<String>, at: usize, text: &[u8]) -> Self {
        Self::new(ParseErrorKind::Syntax(msg.into()), at, text)
    }

    /// A type or structural error. These have no useful excerpt.
    pub(crate) fn semantic(msg: impl Into<String>) -> Self {
        Self::new(ParseErrorKind::Semantic(msg.into()), 0, b"")
    }

    /// Render to the CFITSIO error stack, the way the C did.
    pub(crate) fn report(&self) {
        let excerpt = String::from_utf8_lossy(&self.text);
        match &self.kind {
            ParseErrorKind::TooLong(what) => {
                ffpmsg_str(&format!("{what} exceeds maximum length: '{excerpt}...'"));
            }
            ParseErrorKind::Syntax(msg) => {
                if self.text.is_empty() {
                    ffpmsg_str(msg);
                } else {
                    ffpmsg_str(&format!("{msg} at '{excerpt}'"));
                }
            }
            ParseErrorKind::Semantic(msg) => ffpmsg_str(msg),
        }
    }
}
