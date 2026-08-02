//! Tokenizer for the CFITSIO expression language.
//!
//! Replaces the flex DFA that used to live in `src/eval_l.rs`. The pattern for
//! each rule is written as a `nom` parser; flex's *disambiguation* — longest
//! match wins, ties broken by declaration order — is supplied by [`next_token`],
//! because that is precisely what a PEG-style combinator library does not do.
//!
//! The rule order in [`RULES`] is the order of the rules in `eval.l` and is
//! semantically significant: `b101` is a bit string but `b1x0y` is a variable,
//! `T` is a boolean but `TT` is a variable. See `PARSER_SPEC.md` §2.6.

use nom::{
    Parser,
    bytes::complete::{tag, take_while, take_while1},
    character::complete::{char as nchar, digit0, digit1, one_of},
    combinator::recognize,
};

use super::error::{ParseError, ParseErrorKind};
use super::token::{CallKind, Spanned, Tok};
use crate::c_types::c_long;
use crate::eval_defs::MAX_STRLEN;

/// Digit expansions used by the octal and hexadecimal bit-string forms.
const OCT_BITS: [&[u8]; 8] = [
    b"000", b"001", b"010", b"011", b"100", b"101", b"110", b"111",
];
const HEX_BITS: [&[u8]; 16] = [
    b"0000", b"0001", b"0010", b"0011", b"0100", b"0101", b"0110", b"0111", b"1000", b"1001",
    b"1010", b"1011", b"1100", b"1101", b"1110", b"1111",
];

/// Largest bit string the evaluation engine can hold, from `eval.l`'s
/// `bitstring[256]` with a NUL terminator.
const MAX_BITS: usize = 255;

// ---------------------------------------------------------------------------
// Pattern matchers
//
// Each returns the length of the longest match at the start of `i`, or None.
// The bodies transcribe the `%{ ... %}` definitions block of eval.l verbatim;
// the definition is quoted above each one.
// ---------------------------------------------------------------------------

/// Length of the match produced by a `nom` recognizer, if any.
fn len_of<'a, P>(mut p: P, i: &'a [u8]) -> Option<usize>
where
    P: Parser<&'a [u8], Output = &'a [u8], Error = nom::error::Error<&'a [u8]>>,
{
    match p.parse(i) {
        Ok((rest, _)) => Some(i.len() - rest.len()),
        Err(_) => None,
    }
}

/// A pattern matcher: the length of the longest match at the start of the
/// input, or `None`.
type Matcher = fn(&[u8]) -> Option<usize>;

/// Longest match among several recognizers.
///
/// `alt` is an *ordered* choice, which is not what a flex alternation means:
/// `{real}`'s three branches overlap and the longest has to win (`1.5e3` is one
/// token, not `1.5` followed by `e3`).
fn longest(alts: &[Matcher], i: &[u8]) -> Option<usize> {
    alts.iter().filter_map(|f| f(i)).max()
}

fn is_bit_digit(c: u8) -> bool {
    matches!(c, b'0' | b'1' | b'x' | b'X')
}
fn is_oct_digit(c: u8) -> bool {
    c.is_ascii_digit() && c < b'8' || c == b'x' || c == b'X'
}
fn is_hex_digit(c: u8) -> bool {
    c.is_ascii_hexdigit() || c == b'x' || c == b'X'
}
fn is_ident_start(c: u8) -> bool {
    c.is_ascii_alphabetic() || c == b'_'
}
fn is_ident_cont(c: u8) -> bool {
    c.is_ascii_alphanumeric() || c == b'_'
}

/// `bit ([bB][01xX]+)`
fn m_bit(i: &[u8]) -> Option<usize> {
    len_of(recognize((one_of("bB"), take_while1(is_bit_digit))), i)
}

/// `oct ([oO][01234567xX]+)`
fn m_oct(i: &[u8]) -> Option<usize> {
    len_of(recognize((one_of("oO"), take_while1(is_oct_digit))), i)
}

/// `hex ([hH][0123456789aAbBcCdDeEfFxX]+)`
fn m_hex(i: &[u8]) -> Option<usize> {
    len_of(recognize((one_of("hH"), take_while1(is_hex_digit))), i)
}

/// `bitconst (0b[01]+)`
fn m_bitconst(i: &[u8]) -> Option<usize> {
    len_of(
        recognize((tag("0b"), take_while1(|c| matches!(c, b'0' | b'1')))),
        i,
    )
}

/// `octconst (0o[01234567]+)`
fn m_octconst(i: &[u8]) -> Option<usize> {
    len_of(
        recognize((
            tag("0o"),
            take_while1(|c: u8| c.is_ascii_digit() && c < b'8'),
        )),
        i,
    )
}

/// `hexconst (0x[0123456789aAbBcCdDeEfF]+)`
fn m_hexconst(i: &[u8]) -> Option<usize> {
    len_of(
        recognize((tag("0x"), take_while1(|c: u8| c.is_ascii_hexdigit()))),
        i,
    )
}

/// `integer [0-9]+`
fn m_integer(i: &[u8]) -> Option<usize> {
    len_of(recognize(digit1), i)
}

/// `boolean (t|f|T|F)`
fn m_boolean(i: &[u8]) -> Option<usize> {
    len_of(recognize(one_of("tfTF")), i)
}

/// `real ([0-9]*"."[0-9]+)|([0-9]*"."*[0-9]+[eEdD][+-]?[0-9]+)|([0-9]+".")`
///
/// The third branch requires at least one digit; allowing zero digits is what
/// made a lone `.` lex as `0.0` (`PARSER_SPEC.md` §6.2).
fn m_real(i: &[u8]) -> Option<usize> {
    fn frac(i: &[u8]) -> Option<usize> {
        len_of(recognize((digit0, nchar('.'), digit1)), i)
    }
    /// `[0-9]*"."*[0-9]+[eEdD][+-]?[0-9]+`
    ///
    /// Written by hand rather than with `digit0` + `digit1`, because those are
    /// greedy and nom does not backtrack: `digit0` would eat the `1` of
    /// `1e308`, leaving `digit1` with nothing and failing the whole branch.
    /// A regex engine backtracks and matches. The mantissa is therefore split
    /// explicitly into digits / dots / digits.
    fn expo(i: &[u8]) -> Option<usize> {
        let digits = |from: usize| {
            let mut k = from;
            while k < i.len() && i[k].is_ascii_digit() {
                k += 1;
            }
            k
        };
        let lead = digits(0);
        let mut dots = lead;
        while dots < i.len() && i[dots] == b'.' {
            dots += 1;
        }
        let mant = digits(dots);
        /* the trailing `[0-9]+` is mandatory: with dots it must come after
        them, without dots the leading run supplies it */
        if dots > lead {
            if mant == dots {
                return None;
            }
        } else if lead == 0 {
            return None;
        }

        let mut k = mant;
        if k >= i.len() || !matches!(i[k], b'e' | b'E' | b'd' | b'D') {
            return None;
        }
        k += 1;
        if k < i.len() && (i[k] == b'+' || i[k] == b'-') {
            k += 1;
        }
        let end = digits(k);
        (end > k).then_some(end)
    }
    fn trailing_dot(i: &[u8]) -> Option<usize> {
        len_of(recognize((digit1, nchar('.'))), i)
    }
    longest(&[frac, expo, trailing_dot], i)
}

/// `constant ("#"[a-zA-Z0-9_]+)|("#""$"[^\n]*"$")`
fn m_constant(i: &[u8]) -> Option<usize> {
    fn plain(i: &[u8]) -> Option<usize> {
        len_of(recognize((nchar('#'), take_while1(is_ident_cont))), i)
    }
    /// `#$…$` — `[^\n]*` is greedy, so this runs to the *last* `$` on the line.
    fn dollar(i: &[u8]) -> Option<usize> {
        if !i.starts_with(b"#$") {
            return None;
        }
        let line_end = i.iter().position(|&c| c == b'\n').unwrap_or(i.len());
        i[2..line_end]
            .iter()
            .rposition(|&c| c == b'$')
            .map(|p| p + 3)
    }
    longest(&[plain, dollar], i)
}

/// `string ([\"][^\"\n]*[\"])|([\'][^\'\n]*[\'])`
fn m_string(i: &[u8]) -> Option<usize> {
    fn quoted(q: u8, i: &[u8]) -> Option<usize> {
        if i.first() != Some(&q) {
            return None;
        }
        i[1..]
            .iter()
            .position(|&c| c == q || c == b'\n')
            .filter(|&p| i[1 + p] == q)
            .map(|p| p + 2)
    }
    fn dq(i: &[u8]) -> Option<usize> {
        quoted(b'"', i)
    }
    fn sq(i: &[u8]) -> Option<usize> {
        quoted(b'\'', i)
    }
    longest(&[dq, sq], i)
}

/// `variable ([a-zA-Z_][a-zA-Z0-9_]*)|("$"[^$\n]*"$")`
fn m_variable(i: &[u8]) -> Option<usize> {
    fn plain(i: &[u8]) -> Option<usize> {
        len_of(
            recognize((
                nom::character::complete::satisfy(|c| is_ident_start(c as u8)),
                take_while(is_ident_cont),
            )),
            i,
        )
    }
    /// `$…$` — the body excludes `$`, so this stops at the *first* closing `$`.
    fn dollar(i: &[u8]) -> Option<usize> {
        if i.first() != Some(&b'$') {
            return None;
        }
        i[1..]
            .iter()
            .position(|&c| c == b'$' || c == b'\n')
            .filter(|&p| i[1 + p] == b'$')
            .map(|p| p + 2)
    }
    longest(&[plain, dollar], i)
}

/// `function [a-zA-Z][a-zA-Z0-9]+"("`
///
/// Note: at least *two* name characters, no underscore, and no space before
/// the parenthesis — `X(` is a variable followed by `(`.
fn m_function(i: &[u8]) -> Option<usize> {
    len_of(
        recognize((
            one_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"),
            take_while1(|c: u8| c.is_ascii_alphanumeric()),
            nchar('('),
        )),
        i,
    )
}

/// Matcher for a fixed set of literal spellings, longest first.
fn m_literals(alts: &[&str], i: &[u8]) -> Option<usize> {
    alts.iter()
        .filter(|a| i.starts_with(a.as_bytes()))
        .map(|a| a.len())
        .max()
}

fn m_intcast(i: &[u8]) -> Option<usize> {
    m_literals(&["(int)", "(INT)"], i)
}
fn m_fltcast(i: &[u8]) -> Option<usize> {
    m_literals(&["(float)", "(FLOAT)", "(double)", "(DOUBLE)"], i)
}
fn m_power(i: &[u8]) -> Option<usize> {
    m_literals(&["^", "**"], i)
}
fn m_not(i: &[u8]) -> Option<usize> {
    m_literals(&["!", ".not.", ".NOT.", "not.", "NOT."], i)
}
fn m_or(i: &[u8]) -> Option<usize> {
    m_literals(&["||", ".or.", ".OR.", "or.", "OR."], i)
}
fn m_and(i: &[u8]) -> Option<usize> {
    m_literals(&["&&", ".and.", ".AND.", "and.", "AND."], i)
}
fn m_eq(i: &[u8]) -> Option<usize> {
    m_literals(&["==", ".eq.", ".EQ.", "eq.", "EQ."], i)
}
fn m_ne(i: &[u8]) -> Option<usize> {
    m_literals(&["!=", ".ne.", ".NE.", "ne.", "NE."], i)
}
fn m_gt(i: &[u8]) -> Option<usize> {
    m_literals(&[">", ".gt.", ".GT.", "gt.", "GT."], i)
}
fn m_lt(i: &[u8]) -> Option<usize> {
    m_literals(&["<", ".lt.", ".LT.", "lt.", "LT."], i)
}
fn m_gte(i: &[u8]) -> Option<usize> {
    m_literals(&[">=", "=>", ".ge.", ".GE.", "ge.", "GE."], i)
}
fn m_lte(i: &[u8]) -> Option<usize> {
    m_literals(&["<=", "=<", ".le.", ".LE.", "le.", "LE."], i)
}
fn m_xor(i: &[u8]) -> Option<usize> {
    m_literals(&["^^", ".xor.", ".XOR."], i)
}
fn m_newline(i: &[u8]) -> Option<usize> {
    (i.first() == Some(&b'\n')).then_some(1)
}

/// Identifies which `eval.l` rule matched, so that the token can be built and
/// so that ties in match length are broken by declaration order.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Rule {
    Bit,
    Oct,
    Hex,
    BitConst,
    OctConst,
    HexConst,
    Integer,
    Boolean,
    Real,
    Constant,
    Str,
    Variable,
    Function,
    IntCast,
    FltCast,
    Power,
    Not,
    Or,
    And,
    Eq,
    Ne,
    Gt,
    Lt,
    Gte,
    Lte,
    Xor,
    Newline,
}

/// The rules of `eval.l`, in declaration order. Earlier entries win ties.
const RULES: &[(Rule, Matcher)] = &[
    (Rule::Bit, m_bit),
    (Rule::Oct, m_oct),
    (Rule::Hex, m_hex),
    (Rule::BitConst, m_bitconst),
    (Rule::OctConst, m_octconst),
    (Rule::HexConst, m_hexconst),
    (Rule::Integer, m_integer),
    (Rule::Boolean, m_boolean),
    (Rule::Real, m_real),
    (Rule::Constant, m_constant),
    (Rule::Str, m_string),
    (Rule::Variable, m_variable),
    (Rule::Function, m_function),
    (Rule::IntCast, m_intcast),
    (Rule::FltCast, m_fltcast),
    (Rule::Power, m_power),
    (Rule::Not, m_not),
    (Rule::Or, m_or),
    (Rule::And, m_and),
    (Rule::Eq, m_eq),
    (Rule::Ne, m_ne),
    (Rule::Gt, m_gt),
    (Rule::Lt, m_lt),
    (Rule::Gte, m_gte),
    (Rule::Lte, m_lte),
    (Rule::Xor, m_xor),
    (Rule::Newline, m_newline),
];

/// Pick the rule flex would pick: longest match, ties to the earliest rule.
fn next_rule(i: &[u8]) -> Option<(Rule, usize)> {
    let mut best: Option<(Rule, usize)> = None;
    for &(rule, m) in RULES {
        if let Some(n) = m(i)
            && n > 0
            && best.is_none_or(|(_, bn)| n > bn)
        {
            best = Some((rule, n));
        }
    }
    best
}

// ---------------------------------------------------------------------------
// Token construction
// ---------------------------------------------------------------------------

/// `atol`: leading digits, saturating at the `c_long` range like glibc.
fn atol(text: &[u8]) -> c_long {
    let mut v: i128 = 0;
    for &c in text {
        if !c.is_ascii_digit() {
            break;
        }
        v = v * 10 + i128::from(c - b'0');
        if v > i128::from(c_long::MAX) {
            return c_long::MAX;
        }
    }
    v as c_long
}

/// `atof`: parses a leading C double. Note that the Fortran `d`/`D` exponent
/// spellings accepted by `{real}` are *not* understood here, exactly as in C —
/// `1.5d3` evaluates to 1.5 (`PARSER_SPEC.md` §6.4).
fn atof(text: &[u8]) -> f64 {
    let mut end = 0;
    let b = text;
    if end < b.len() && (b[end] == b'+' || b[end] == b'-') {
        end += 1;
    }
    while end < b.len() && b[end].is_ascii_digit() {
        end += 1;
    }
    if end < b.len() && b[end] == b'.' {
        end += 1;
        while end < b.len() && b[end].is_ascii_digit() {
            end += 1;
        }
    }
    if end < b.len() && (b[end] == b'e' || b[end] == b'E') {
        let save = end;
        let mut j = end + 1;
        if j < b.len() && (b[j] == b'+' || b[j] == b'-') {
            j += 1;
        }
        if j < b.len() && b[j].is_ascii_digit() {
            while j < b.len() && b[j].is_ascii_digit() {
                j += 1;
            }
            end = j;
        } else {
            end = save;
        }
    }
    core::str::from_utf8(&b[..end])
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(0.0)
}

/// Expand an octal or hexadecimal bit-string literal.
///
/// `x`/`X` become a run of `x` (a "don't care" bit). Returns `None` if the
/// result would exceed the engine's 255-bit buffer.
fn expand_bits(body: &[u8], hex: bool) -> Option<Vec<u8>> {
    let width = if hex { 4 } else { 3 };
    let mut out = Vec::with_capacity(body.len() * width);
    for &c in body {
        let chunk: &[u8] = if c == b'x' || c == b'X' {
            if hex { b"xxxx" } else { b"xxx" }
        } else if hex {
            HEX_BITS[(c as char).to_digit(16)? as usize]
        } else {
            OCT_BITS[usize::from(c - b'0')]
        };
        if out.len() + chunk.len() > MAX_BITS {
            return None;
        }
        out.extend_from_slice(chunk);
    }
    Some(out)
}

/// Value of a `0x` / `0o` / `0b` integer literal. Digits are shifted in and
/// allowed to wrap, matching the C, which is how `0xffffffffffffffff` is -1.
fn radix_const(digits: &[u8], bits: u32) -> c_long {
    let mut v: u64 = 0;
    for &c in digits {
        let d = u64::from((c as char).to_digit(16).unwrap_or(0));
        v = (v << bits) | d;
    }
    v as c_long
}

/// Classify a function name into one of the six call tokens.
fn call_kind(name: &[u8]) -> CallKind {
    match name {
        b"BOX" | b"CIRCLE" | b"ELLIPSE" | b"NEAR" | b"ISNULL" => CallKind::BFunction,
        b"GTIFILTER" => CallKind::GtiFilter,
        b"GTIOVERLAP" => CallKind::GtiOverlap,
        b"GTIFIND" => CallKind::GtiFind,
        b"REGFILTER" => CallKind::RegFilter,
        b"STRSTR" => CallKind::IFunction,
        _ => CallKind::Function,
    }
}

fn too_long(what: &str, text: &[u8], at: usize) -> ParseError {
    ParseError::new(ParseErrorKind::TooLong(what.into()), at, text)
}

/// Build the token for `rule` from the matched text.
fn build(rule: Rule, text: &[u8], at: usize) -> Result<Tok, ParseError> {
    let max = MAX_STRLEN as usize;
    Ok(match rule {
        Rule::Bit => {
            let body = &text[1..];
            if body.len() >= max {
                return Err(too_long("Bit string", text, at));
            }
            Tok::BitStr(body.to_vec())
        }
        Rule::Oct | Rule::Hex => {
            let hex = rule == Rule::Hex;
            let what = if hex { "Hex string" } else { "Bit string" };
            if text.len() >= 256 {
                return Err(too_long(what, text, at));
            }
            match expand_bits(&text[1..], hex) {
                Some(bits) => Tok::BitStr(bits),
                None => return Err(too_long(what, text, at)),
            }
        }
        Rule::BitConst => Tok::Long(radix_const(&text[2..], 1)),
        Rule::OctConst => Tok::Long(radix_const(&text[2..], 3)),
        Rule::HexConst => Tok::Long(radix_const(&text[2..], 4)),
        Rule::Integer => Tok::Long(atol(text)),
        Rule::Boolean => Tok::Boolean(text[0] == b't' || text[0] == b'T'),
        Rule::Real => Tok::Double(atof(text)),
        Rule::Constant => {
            /* Builtin constants first, case-insensitively. */
            let upper: Vec<u8> = text.iter().map(|c| c.to_ascii_uppercase()).collect();
            match upper.as_slice() {
                b"#PI" => Tok::Double(4.0 * 1.0_f64.atan()),
                b"#E" => Tok::Double(1.0_f64.exp()),
                b"#DEG" => Tok::Double(4.0 * 1.0_f64.atan() / 180.0),
                b"#ROW" => Tok::RowRef,
                b"#NULL" => Tok::NullRef,
                b"#SNULL" => Tok::SNullRef,
                _ => {
                    if text.len() > 1 && text[1] == b'$' {
                        /* '#$name$' -> '#name'; C checks len against
                        MAX_STRLEN-1 because it prepends the '#'. */
                        let body = &text[2..text.len() - 1];
                        if body.len() >= max - 1 {
                            return Err(too_long("Keyword string", text, at));
                        }
                        Tok::Keyword(body.to_vec())
                    } else {
                        Tok::Keyword(text[1..].to_vec())
                    }
                }
            }
        }
        Rule::Str => {
            let body = &text[1..text.len() - 1];
            if body.len() >= max {
                return Err(too_long("String", text, at));
            }
            Tok::Str(body.to_vec())
        }
        Rule::Variable => {
            if text.len() >= max {
                return Err(too_long("Variable", text, at));
            }
            let name = if text[0] == b'$' {
                &text[1..text.len() - 1]
            } else {
                text
            };
            Tok::Ident(name.to_vec())
        }
        Rule::Function => {
            if text.len() >= max {
                return Err(too_long("Function", text, at));
            }
            /* drop the '(' and upper-case, as the C does before FSTRCMP */
            let name: Vec<u8> = text[..text.len() - 1]
                .iter()
                .map(|c| c.to_ascii_uppercase())
                .collect();
            Tok::Call(call_kind(&name), name)
        }
        Rule::IntCast => Tok::IntCast,
        Rule::FltCast => Tok::FltCast,
        Rule::Power => Tok::Power,
        Rule::Not => Tok::Not,
        Rule::Or => Tok::Or,
        Rule::And => Tok::And,
        Rule::Eq => Tok::Eq,
        Rule::Ne => Tok::Ne,
        Rule::Gt => Tok::Gt,
        Rule::Lt => Tok::Lt,
        Rule::Gte => Tok::Gte,
        Rule::Lte => Tok::Lte,
        Rule::Xor => Tok::Xor,
        Rule::Newline => Tok::Newline,
    })
}

/// Tokenize a whole expression.
///
/// The input is the raw expression bytes; a trailing `\n` is supplied by the
/// caller (`ffiprs` appends one) and becomes a [`Tok::Newline`].
pub(crate) fn tokenize(input: &[u8]) -> Result<Vec<Spanned>, ParseError> {
    let mut out = Vec::new();
    let mut pos = 0usize;
    while pos < input.len() {
        let c = input[pos];
        if c == b' ' || c == b'\t' {
            pos += 1;
            continue;
        }
        if c == 0 {
            break;
        }
        match next_rule(&input[pos..]) {
            Some((rule, n)) => {
                let text = &input[pos..pos + n];
                out.push(Spanned {
                    tok: build(rule, text, pos)?,
                    at: pos,
                });
                pos += n;
            }
            None => {
                /* flex's final `.` rule: emit the byte itself */
                out.push(Spanned {
                    tok: Tok::Char(c),
                    at: pos,
                });
                pos += 1;
            }
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn toks(s: &str) -> Vec<Tok> {
        tokenize(s.as_bytes())
            .unwrap()
            .into_iter()
            .map(|t| t.tok)
            .collect()
    }

    fn ident(s: &str) -> Tok {
        Tok::Ident(s.as_bytes().to_vec())
    }

    #[test]
    fn longest_match_beats_rule_order() {
        /* `bit` is declared before `variable` and wins the 4-char tie ... */
        assert_eq!(toks("b101"), [Tok::BitStr(b"101".to_vec())]);
        /* ... but `variable` matches one char more here, so it wins */
        assert_eq!(toks("b1x0y"), [ident("b1x0y")]);
        assert_eq!(toks("B01x"), [Tok::BitStr(b"01x".to_vec())]);
    }

    #[test]
    fn boolean_only_when_a_single_letter() {
        assert_eq!(toks("T"), [Tok::Boolean(true)]);
        assert_eq!(toks("f"), [Tok::Boolean(false)]);
        assert_eq!(toks("TT"), [ident("TT")]);
        assert_eq!(toks("false"), [ident("false")]);
    }

    #[test]
    fn bare_dot_is_not_a_number() {
        assert_eq!(toks("."), [Tok::Char(b'.')]);
        assert_eq!(toks(".."), [Tok::Char(b'.'), Tok::Char(b'.')]);
        assert_eq!(toks("1."), [Tok::Double(1.0)]);
        assert_eq!(toks(".5"), [Tok::Double(0.5)]);
        assert_eq!(toks("1.2.3"), [Tok::Double(1.2), Tok::Double(0.3)]);
    }

    #[test]
    fn real_exponent_forms() {
        assert_eq!(toks("1.5e3"), [Tok::Double(1500.0)]);
        assert_eq!(toks("1.5e-3"), [Tok::Double(0.0015)]);
        /* the Fortran exponent is matched by {real} but not understood by
        atof, so the exponent is dropped -- as in C */
        assert_eq!(toks("1.5d3"), [Tok::Double(1.5)]);
        /* the exponent branch needs a digit before [eEdD] */
        assert_eq!(toks("1.e5"), [Tok::Double(1.0), ident("e5")]);
        assert_eq!(toks("1.5e"), [Tok::Double(1.5), ident("e")]);
    }

    #[test]
    fn hex_literals_are_case_insensitive() {
        assert_eq!(toks("0x1f"), [Tok::Long(31)]);
        assert_eq!(toks("0x1F"), [Tok::Long(31)]);
        assert_eq!(toks("0xFF"), [Tok::Long(255)]);
        assert_eq!(toks("0xffffffffffffffff"), [Tok::Long(-1)]);
        /* the prefixes need at least one digit */
        assert_eq!(toks("0x"), [Tok::Long(0), ident("x")]);
    }

    #[test]
    fn function_needs_two_name_chars_and_no_space() {
        assert_eq!(toks("ab("), [Tok::Call(CallKind::Function, b"AB".to_vec())]);
        assert_eq!(
            toks("X(1)"),
            [ident("X"), Tok::Char(b'('), Tok::Long(1), Tok::Char(b')')]
        );
        assert_eq!(toks("a ("), [ident("a"), Tok::Char(b'(')]);
        assert_eq!(
            toks("ISNULL("),
            [Tok::Call(CallKind::BFunction, b"ISNULL".to_vec())]
        );
        assert_eq!(
            toks("strstr("),
            [Tok::Call(CallKind::IFunction, b"STRSTR".to_vec())]
        );
    }

    #[test]
    fn dollar_forms_differ_in_greediness() {
        /* '#$…$' runs to the last '$'; '$…$' stops at the first */
        assert_eq!(toks("#$a$b$"), [Tok::Keyword(b"a$b".to_vec())]);
        assert_eq!(toks("$a$b$"), [ident("a"), ident("b"), Tok::Char(b'$')]);
        assert_eq!(toks("$MY COL$"), [ident("MY COL")]);
    }

    #[test]
    fn operator_spellings() {
        assert_eq!(toks("a=>b"), [ident("a"), Tok::Gte, ident("b")]);
        assert_eq!(toks("a=<b"), [ident("a"), Tok::Lte, ident("b")]);
        assert_eq!(toks("a**b"), [ident("a"), Tok::Power, ident("b")]);
        assert_eq!(
            toks("a***b"),
            [ident("a"), Tok::Power, Tok::Char(b'*'), ident("b")]
        );
        assert_eq!(toks("a.NE.b"), [ident("a"), Tok::Ne, ident("b")]);
        /* mixed case is not one of the listed spellings */
        assert_eq!(
            toks("a.Ne.b"),
            [
                ident("a"),
                Tok::Char(b'.'),
                ident("Ne"),
                Tok::Char(b'.'),
                ident("b")
            ]
        );
        assert_eq!(toks("(int)x"), [Tok::IntCast, ident("x")]);
        assert_eq!(
            toks("(int )x"),
            [Tok::Char(b'('), ident("int"), Tok::Char(b')'), ident("x")]
        );
    }

    #[test]
    fn bit_string_expansions() {
        assert_eq!(toks("h1f"), [Tok::BitStr(b"00011111".to_vec())]);
        assert_eq!(toks("o17"), [Tok::BitStr(b"001111".to_vec())]);
        assert_eq!(toks("hxx"), [Tok::BitStr(b"xxxxxxxx".to_vec())]);
        assert_eq!(toks("oxx"), [Tok::BitStr(b"xxxxxx".to_vec())]);
    }

    #[test]
    fn unterminated_string_degrades_to_characters() {
        assert_eq!(toks("\"abc"), [Tok::Char(b'"'), ident("abc")]);
    }

    #[test]
    fn builtin_constants() {
        assert!(matches!(toks("#PI")[0], Tok::Double(_)));
        assert!(matches!(toks("#pi")[0], Tok::Double(_)));
        assert_eq!(toks("#ROW"), [Tok::RowRef]);
        assert_eq!(toks("#null"), [Tok::NullRef]);
        assert_eq!(toks("#SNULL"), [Tok::SNullRef]);
        assert_eq!(toks("#MYKEY"), [Tok::Keyword(b"MYKEY".to_vec())]);
    }

    #[test]
    fn overlong_literals_are_rejected() {
        let long = "x".repeat(300);
        assert!(tokenize(format!("'{long}'").as_bytes()).is_err());
        assert!(tokenize(long.as_bytes()).is_err());
        assert!(tokenize(format!("b{}", "1".repeat(300)).as_bytes()).is_err());
    }
}
