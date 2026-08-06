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
    ///
    /// The payload is a *reason*, for whoever is reading this code or a
    /// backtrace -- it is deliberately not what gets reported. Bison had no
    /// such detail to give: a token sequence with no production simply failed,
    /// and `yyerror` pushed the bare string "syntax error". Since the nom front
    /// end knows exactly what it wanted and where, it would be easy to say so,
    /// but a caller diffing rsfitsio's error stack against CFITSIO's has to see
    /// the same text. See `report`.
    Syntax(String),
    /// The expression parsed but is not well typed, or violates one of the
    /// structural constraints in `PARSER_SPEC.md` §3.4.
    Semantic(String),
    /// A node builder failed and has already said why.
    ///
    /// `#define TEST(a) if( (a)<0 ) YYERROR` -- the check the C makes after
    /// almost every `New_*` call is silent. The builder pushed its own message
    /// and set the status; `YYERROR` only unwinds. So `INTCOL[1]` reports
    /// "Cannot index a scalar value" and nothing after it.
    ///
    /// "Couldn't build node structure: out of memory?" is *not* this case: it
    /// belongs to the top-level `line:` rule, for a node that came back
    /// negative without any production having tested it.
    Aborted,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct ParseError {
    pub(crate) kind: ParseErrorKind,
    /// Byte offset into the expression at which the problem was detected.
    pub(crate) at: usize,
    /// A short excerpt of the offending text, for the message.
    pub(crate) text: Vec<u8>,
}

/// Push a parser message the way the bison `yyerror` did.
///
/// ```c
/// static void yyerror(yyscan_t scanner, ParseData *lParse, char *s)
/// {
///     char msg[80];
///     strncpy(msg, s, 80);
///     msg[79] = '\0';
///     ffpmsg(msg);
/// }
/// ```
///
/// The fixed buffer means every message the parser raises is cut to 79
/// characters *before* it reaches the stack, so the longer ones -- the bitwise
/// and combined-string-length complaints -- lose their tails. `ffpmsg` itself
/// splits at 80 and would otherwise carry the remainder onto a second entry,
/// which is why the truncation has to happen here rather than there.
fn yyerror(msg: &str) {
    /* byte-wise, as `strncpy` is; the messages are all ASCII literals, and a
    char boundary is only in question if one ever stops being */
    let clamped = match msg.char_indices().nth(79) {
        Some((i, _)) => &msg[..i],
        None => msg,
    };
    ffpmsg_str(clamped);
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

    /// A builder failed having already reported. See [`ParseErrorKind::Aborted`].
    pub(crate) fn aborted() -> Self {
        Self::new(ParseErrorKind::Aborted, 0, b"")
    }

    /// Render to the CFITSIO error stack, the way the C did.
    pub(crate) fn report(&self) {
        let excerpt = String::from_utf8_lossy(&self.text);
        match &self.kind {
            ParseErrorKind::TooLong(what) => {
                yyerror(&format!("{what} exceeds maximum length: '{excerpt}...'"));
            }
            /* Everything the grammar rejects reaches the stack as this one
            string, whatever the `Syntax` payload says and wherever it was
            raised -- see the variant's own note. No position or excerpt is
            appended either, though `at` and `text` carry one, because bison
            pushed the bare message and nothing else. Both are still kept on the
            error for debugging and for callers that want to point at a token. */
            ParseErrorKind::Syntax(_) => yyerror("syntax error"),
            ParseErrorKind::Semantic(msg) => yyerror(msg),
            /* the builder has already reported; TEST() adds nothing */
            ParseErrorKind::Aborted => {}
        }
    }
}
