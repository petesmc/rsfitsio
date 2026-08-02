//! Precedence-climbing (Pratt) parser over the token stream.
//!
//! This replaces bison's LALR driver. The binding-power table below is a direct
//! transcription of the `%left` / `%right` block of `eval.y`
//! (`PARSER_SPEC.md` §3.1): level *n* counted from the bottom becomes
//! `(2n, 2n+1)` for a left-associative operator and `(2n+1, 2n)` for a
//! right-associative one.
//!
//! Two levels are deliberately unlike C and must not be "fixed":
//!
//! * bitwise `| & XOR` bind **tighter** than `* /`, so `1|2*3` is `(1|2)*3` = 9;
//! * unary minus binds **tighter** than `**`, so `-2**2` is `(-2)**2` = +4.
//!
//! The grammar built here is untyped; sorts are resolved in `super::lower`.

use super::ast::{Ast, AstKind, BinOp, UnOp};
use super::error::ParseError;
use super::token::{Spanned, Tok};

/// Bison's `%left ',' '=' ':' '{' '}'` is level 1 and `%right '?'` is level 2,
/// so the ternary's right operand is parsed at 4 and its left binding power
/// is 5.
const BP_TERNARY: (u8, u8) = (5, 4);

/// `NOT` sits at bison level 11, above everything but the casts, unary minus
/// and `[`. A single binding power cannot reproduce that, because `NOT`'s only
/// productions are `NOT bexpr` and `NOT bits`: when the operand is numeric the
/// LALR parser cannot reduce and keeps shifting until a boolean has been built.
///
/// Binding the operand at 10 — above `AND` (8) and `OR` (6), at or below the
/// comparison levels (10, 12) — reproduces every discriminating case:
///
/// | input | parse | why |
/// |---|---|---|
/// | `!B \|\| B` | `(!B) \|\| B` | `\|\|` is 6 < 10, stop |
/// | `!B && 1>2` | `(!B) && (1>2)` | `&&` is 8 < 10, stop |
/// | `!X > 5` | `!(X > 5)` | `>` is 12 >= 10, absorb |
/// | `!X + 1 > 5` | `!((X+1) > 5)` | `+` is 14, `>` is 12, both absorb |
///
/// The single known divergence is `!B == B`, which the LALR parser groups as
/// `(!B) == B` and this groups as `!(B == B)`. For booleans the two are
/// equivalent (`!(a==b)` is `(!a)==b`), so only the node tree differs.
const BP_NOT_OPERAND: u8 = 10;

/// Casts, bison level 12.
const BP_CAST_OPERAND: u8 = 24;
/// Unary `+` / `-` (`%prec UMINUS`), bison level 13.
const BP_UNARY_OPERAND: u8 = 26;
/// Subscripting, bison level 14 — binds tighter than everything.
const BP_SUBSCRIPT: u8 = 28;

/// Maximum nesting depth.
///
/// The bison parser grew a heap stack, so nesting cost it almost nothing;
/// recursive descent uses the machine stack, and an unbounded parser simply
/// crashes on adversarial input like 10,000 open parentheses.
///
/// Only *right* nesting counts: `1+1+1+…` is folded iteratively by the loop in
/// [`Parser::expr`] and stays at depth 2 however long it is. Parentheses,
/// prefix operators, call arguments and the right operand of `**` do nest.
///
/// The limit is measured, not guessed: an unbounded build overflows between 175
/// and 200 levels in a debug `cargo test` thread (2 MiB stack). 100 leaves a
/// ~75-level margin there and far more in release, while being an order of
/// magnitude deeper than any hand-written expression. It also happens to match
/// bison's `YYINITDEPTH`.
const MAX_DEPTH: u32 = 100;

/// Left and right binding powers of an infix operator, or `None` if the token
/// is not infix (`,` `:` `)` `]` `}` are delimiters, not operators).
fn infix_bp(t: &Tok) -> Option<(u8, u8, BinOp)> {
    let r = match t {
        Tok::Or => (6, 7, BinOp::Or),
        Tok::And => (8, 9, BinOp::And),
        Tok::Eq => (10, 11, BinOp::Eq),
        Tok::Ne => (10, 11, BinOp::Ne),
        Tok::Gt => (12, 13, BinOp::Gt),
        Tok::Lt => (12, 13, BinOp::Lt),
        Tok::Gte => (12, 13, BinOp::Gte),
        Tok::Lte => (12, 13, BinOp::Lte),
        Tok::Xor => (18, 19, BinOp::BitXor),
        /* right-associative */
        Tok::Power => (21, 20, BinOp::Pow),
        Tok::Char(c) => match c {
            b'~' => (10, 11, BinOp::Approx),
            b'+' => (14, 15, BinOp::Add),
            b'-' => (14, 15, BinOp::Sub),
            b'%' => (14, 15, BinOp::Mod),
            b'*' => (16, 17, BinOp::Mul),
            b'/' => (16, 17, BinOp::Div),
            b'|' => (18, 19, BinOp::BitOr),
            b'&' => (18, 19, BinOp::BitAnd),
            _ => return None,
        },
        _ => return None,
    };
    Some(r)
}

pub(crate) struct Parser<'a> {
    toks: &'a [Spanned],
    pos: usize,
    depth: u32,
}

impl<'a> Parser<'a> {
    fn peek(&self) -> Option<&'a Tok> {
        self.toks.get(self.pos).map(|s| &s.tok)
    }

    fn at(&self) -> usize {
        self.toks
            .get(self.pos)
            .or_else(|| self.toks.last())
            .map_or(0, |s| s.at)
    }

    fn bump(&mut self) {
        self.pos += 1;
    }

    /// Describe the current token for an error message.
    fn here(&self) -> Vec<u8> {
        match self.peek() {
            None => b"end of expression".to_vec(),
            Some(Tok::Newline) => b"end of expression".to_vec(),
            Some(Tok::Char(c)) => vec![*c],
            Some(Tok::Ident(n)) | Some(Tok::Call(_, n)) | Some(Tok::Keyword(n)) => n.clone(),
            Some(_) => b"operator".to_vec(),
        }
    }

    fn err(&self, msg: impl Into<String>) -> ParseError {
        let here = self.here();
        ParseError::syntax(msg, self.at(), &here)
    }

    /// Consume the bare character `c`, or fail.
    fn expect_char(&mut self, c: u8) -> Result<(), ParseError> {
        match self.peek() {
            Some(t) if t.is_char(c) => {
                self.bump();
                Ok(())
            }
            _ => Err(self.err(format!("expected '{}'", c as char))),
        }
    }

    fn eat_char(&mut self, c: u8) -> bool {
        match self.peek() {
            Some(t) if t.is_char(c) => {
                self.bump();
                true
            }
            _ => false,
        }
    }

    fn enter(&mut self) -> Result<(), ParseError> {
        self.depth += 1;
        if self.depth > MAX_DEPTH {
            return Err(self.err("expression nested too deeply"));
        }
        Ok(())
    }

    fn leave(&mut self) {
        self.depth -= 1;
    }

    /// Parse a full expression with the given minimum binding power.
    fn expr(&mut self, min_bp: u8) -> Result<Ast, ParseError> {
        self.enter()?;
        let mut lhs = self.prefix()?;

        loop {
            /* postfix: subscripting binds tighter than every infix operator */
            if BP_SUBSCRIPT >= min_bp
                && self.peek().is_some_and(|t| t.is_char(b'['))
            {
                lhs = self.subscript(lhs)?;
                continue;
            }

            let Some(tok) = self.peek() else { break };

            /* the two ternaries share level-1/2 handling */
            if tok.is_char(b'?') {
                let (lbp, rbp) = BP_TERNARY;
                if lbp < min_bp {
                    break;
                }
                let at = self.at();
                self.bump();
                let then = self.expr(0)?;
                self.expect_char(b':')?;
                let els = self.expr(rbp)?;
                lhs = Ast::new(
                    AstKind::Ternary {
                        cond: Box::new(lhs),
                        then: Box::new(then),
                        els: Box::new(els),
                    },
                    at,
                );
                continue;
            }
            if tok.is_char(b'=') {
                /* `val = lo : hi`; '=' and ':' are bison level 1 */
                if 2 < min_bp {
                    break;
                }
                let at = self.at();
                self.bump();
                let lo = self.expr(3)?;
                self.expect_char(b':')?;
                let hi = self.expr(3)?;
                lhs = Ast::new(
                    AstKind::Range {
                        val: Box::new(lhs),
                        lo: Box::new(lo),
                        hi: Box::new(hi),
                    },
                    at,
                );
                continue;
            }

            let Some((lbp, rbp, op)) = infix_bp(tok) else {
                break;
            };
            if lbp < min_bp {
                break;
            }
            let at = self.at();
            self.bump();
            let rhs = self.expr(rbp)?;
            lhs = Ast::new(
                AstKind::Binary {
                    op,
                    lhs: Box::new(lhs),
                    rhs: Box::new(rhs),
                },
                at,
            );
        }

        self.leave();
        Ok(lhs)
    }

    /// `base [ i1, … i5 ]`
    fn subscript(&mut self, base: Ast) -> Result<Ast, ParseError> {
        let at = self.at();
        self.expect_char(b'[')?;
        let mut idx = Vec::new();
        loop {
            idx.push(self.expr(0)?);
            if !self.eat_char(b',') {
                break;
            }
        }
        self.expect_char(b']')?;
        if idx.len() > 5 {
            return Err(ParseError::syntax(
                "at most 5 subscripts are allowed",
                at,
                b"[",
            ));
        }
        Ok(Ast::new(
            AstKind::Deref {
                base: Box::new(base),
                idx,
            },
            at,
        ))
    }

    fn prefix(&mut self) -> Result<Ast, ParseError> {
        let at = self.at();
        let Some(tok) = self.peek().cloned() else {
            return Err(self.err("unexpected end of expression"));
        };

        let kind = match tok {
            Tok::Long(v) => {
                self.bump();
                AstKind::Long(v)
            }
            Tok::Double(v) => {
                self.bump();
                AstKind::Double(v)
            }
            Tok::Boolean(v) => {
                self.bump();
                AstKind::Boolean(v)
            }
            Tok::Str(s) => {
                self.bump();
                AstKind::Str(s)
            }
            Tok::BitStr(s) => {
                self.bump();
                AstKind::BitStr(s)
            }
            Tok::RowRef => {
                self.bump();
                AstKind::RowRef
            }
            Tok::NullRef => {
                self.bump();
                AstKind::NullRef
            }
            Tok::SNullRef => {
                self.bump();
                AstKind::SNullRef
            }
            Tok::Keyword(name) => {
                self.bump();
                AstKind::Keyword(name)
            }

            Tok::Ident(name) => {
                self.bump();
                /* `NAME{offset}` is only spelled directly on a name, never on a
                parenthesised or computed expression. */
                if self.peek().is_some_and(|t| t.is_char(b'{')) {
                    self.bump();
                    let off = self.expr(0)?;
                    self.expect_char(b'}')?;
                    AstKind::Offset {
                        name,
                        off: Box::new(off),
                    }
                } else {
                    AstKind::Ident(name)
                }
            }

            Tok::Not => {
                self.bump();
                AstKind::Unary {
                    op: UnOp::Not,
                    arg: Box::new(self.expr(BP_NOT_OPERAND)?),
                }
            }
            Tok::IntCast => {
                self.bump();
                AstKind::Unary {
                    op: UnOp::IntCast,
                    arg: Box::new(self.expr(BP_CAST_OPERAND)?),
                }
            }
            Tok::FltCast => {
                self.bump();
                AstKind::Unary {
                    op: UnOp::FltCast,
                    arg: Box::new(self.expr(BP_CAST_OPERAND)?),
                }
            }

            Tok::Call(kind, name) => {
                self.bump();
                let mut args = Vec::new();
                if !self.eat_char(b')') {
                    loop {
                        args.push(self.expr(0)?);
                        if !self.eat_char(b',') {
                            break;
                        }
                    }
                    self.expect_char(b')')?;
                }
                AstKind::Call { kind, name, args }
            }

            Tok::Char(b'(') => {
                self.bump();
                let inner = self.expr(0)?;
                self.expect_char(b')')?;
                /* parentheses are not represented in the tree; `eval.y`'s
                `'(' expr ')'` rule is `$$ = $2` */
                return Ok(inner);
            }
            Tok::Char(b'{') => {
                self.bump();
                let mut items = Vec::new();
                loop {
                    items.push(self.expr(0)?);
                    if !self.eat_char(b',') {
                        break;
                    }
                }
                self.expect_char(b'}')?;
                AstKind::Vector(items)
            }
            /* `%prec UMINUS`: unary plus is the identity, unary minus negates */
            Tok::Char(b'+') => {
                self.bump();
                return self.expr(BP_UNARY_OPERAND);
            }
            Tok::Char(b'-') => {
                self.bump();
                AstKind::Unary {
                    op: UnOp::Neg,
                    arg: Box::new(self.expr(BP_UNARY_OPERAND)?),
                }
            }

            _ => return Err(self.err("expected an expression")),
        };

        Ok(Ast::new(kind, at))
    }
}

/// Parse a token stream into a syntax tree.
///
/// `eval.y`'s start symbol is `lines: ε | lines line`, where each `line` is an
/// expression terminated by `\n` and `lParse->resultNode` is overwritten each
/// time — so with several lines the last non-empty one wins, and an input of
/// nothing but newlines yields no result at all.
pub(crate) fn parse(toks: &[Spanned]) -> Result<Option<Ast>, ParseError> {
    let mut p = Parser {
        toks,
        pos: 0,
        depth: 0,
    };
    let mut result = None;

    loop {
        /* skip blank lines */
        while p.peek() == Some(&Tok::Newline) {
            p.bump();
        }
        if p.peek().is_none() {
            break;
        }
        result = Some(p.expr(0)?);
        match p.peek() {
            Some(Tok::Newline) => p.bump(),
            /* `line: expr '\n'` requires the terminator. Running out of tokens
            instead is a syntax error, which is what makes `1.e5` a 431: the
            resolution pre-pass truncated the stream at the unresolvable `e5`,
            taking the newline with it. */
            _ => return Err(p.err("unexpected trailing input")),
        }
    }
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::super::lexer::tokenize;
    use super::*;

    /// Render a tree in fully-parenthesised prefix form, so that grouping is
    /// visible in the assertions.
    fn show(a: &Ast) -> String {
        match &a.kind {
            AstKind::Long(v) => format!("{v}"),
            AstKind::Double(v) => format!("{v}"),
            AstKind::Boolean(v) => format!("{v}"),
            AstKind::Str(s) => format!("'{}'", String::from_utf8_lossy(s)),
            AstKind::BitStr(s) => format!("b{}", String::from_utf8_lossy(s)),
            AstKind::RowRef => "#ROW".into(),
            AstKind::NullRef => "#NULL".into(),
            AstKind::SNullRef => "#SNULL".into(),
            AstKind::Keyword(n) => format!("#{}", String::from_utf8_lossy(n)),
            AstKind::Ident(n) => String::from_utf8_lossy(n).into_owned(),
            AstKind::Offset { name, off } => {
                format!("{}{{{}}}", String::from_utf8_lossy(name), show(off))
            }
            AstKind::Unary { op, arg } => format!("({op:?} {})", show(arg)),
            AstKind::Binary { op, lhs, rhs } => {
                format!("({op:?} {} {})", show(lhs), show(rhs))
            }
            AstKind::Ternary { cond, then, els } => {
                format!("(?: {} {} {})", show(cond), show(then), show(els))
            }
            AstKind::Range { val, lo, hi } => {
                format!("(range {} {} {})", show(val), show(lo), show(hi))
            }
            AstKind::Deref { base, idx } => {
                let parts: Vec<String> = idx.iter().map(show).collect();
                format!("(idx {} [{}])", show(base), parts.join(" "))
            }
            AstKind::Vector(items) => {
                let parts: Vec<String> = items.iter().map(show).collect();
                format!("(vec {})", parts.join(" "))
            }
            AstKind::Call { name, args, .. } => {
                let parts: Vec<String> = args.iter().map(show).collect();
                format!("({} {})", String::from_utf8_lossy(name), parts.join(" "))
            }
        }
    }

    /// `ffiprs` appends the terminating newline before parsing, and the
    /// grammar requires it, so the helpers do the same.
    fn p(src: &str) -> String {
        let toks = tokenize(format!("{src}\n").as_bytes()).expect("lex");
        show(&parse(&toks).expect("parse").expect("non-empty"))
    }

    fn fails(src: &str) -> bool {
        match tokenize(format!("{src}\n").as_bytes()) {
            Err(_) => true,
            Ok(toks) => parse(&toks).is_err(),
        }
    }

    #[test]
    fn arithmetic_precedence() {
        assert_eq!(p("2+3*4"), "(Add 2 (Mul 3 4))");
        assert_eq!(p("2*3+4"), "(Add (Mul 2 3) 4)");
        assert_eq!(p("10-3-2"), "(Sub (Sub 10 3) 2)");
        assert_eq!(p("10/2/5"), "(Div (Div 10 2) 5)");
    }

    #[test]
    fn power_is_right_associative() {
        assert_eq!(p("2**3**2"), "(Pow 2 (Pow 3 2))");
        assert_eq!(p("2^3^2"), "(Pow 2 (Pow 3 2))");
    }

    #[test]
    fn unary_minus_binds_tighter_than_power() {
        /* -2**2 is (-2)**2 = +4, not -(2**2) */
        assert_eq!(p("-2**2"), "(Pow (Neg 2) 2)");
        assert_eq!(p("-2*3"), "(Mul (Neg 2) 3)");
    }

    #[test]
    fn bitwise_binds_tighter_than_multiplication() {
        assert_eq!(p("1|2*3"), "(Mul (BitOr 1 2) 3)");
        assert_eq!(p("1+2|4"), "(Add 1 (BitOr 2 4))");
        assert_eq!(p("2*3|1"), "(Mul 2 (BitOr 3 1))");
    }

    #[test]
    fn casts_bind_tighter_than_multiplication() {
        assert_eq!(p("(int)2.7*2"), "(Mul (IntCast 2.7) 2)");
        assert_eq!(p("(int)(2.7*2)"), "(IntCast (Mul 2.7 2))");
    }

    #[test]
    fn not_reaches_through_numeric_operands() {
        assert_eq!(p("!B||B"), "(Or (Not B) B)");
        assert_eq!(p("!B&&B"), "(And (Not B) B)");
        assert_eq!(p("!X>5"), "(Not (Gt X 5))");
        assert_eq!(p("!X>10"), "(Not (Gt X 10))");
        assert_eq!(p("!X+1>5"), "(Not (Gt (Add X 1) 5))");
        assert_eq!(p("!X>5||B"), "(Or (Not (Gt X 5)) B)");
        assert_eq!(p("!X>5&&B"), "(And (Not (Gt X 5)) B)");
    }

    #[test]
    fn ternary_and_range() {
        assert_eq!(p("B?1:2"), "(?: B 1 2)");
        assert_eq!(p("B?B?1:2:3"), "(?: B (?: B 1 2) 3)");
        assert_eq!(p("X=1:2"), "(range X 1 2)");
        assert_eq!(p("B?1=2:3:4"), "(?: B (range 1 2 3) 4)");
        assert_eq!(p("X=Y:Y+1"), "(range X Y (Add Y 1))");
        /* ranges do not chain */
        assert!(fails("X=1:2:3"));
        assert!(fails("1:2"));
    }

    #[test]
    fn offsets_only_on_bare_names() {
        assert_eq!(p("COL{1}"), "COL{1}");
        assert!(fails("(COL){1}"));
        assert!(fails("(COL+1){1}"));
    }

    #[test]
    fn vectors_and_subscripts() {
        assert_eq!(p("{1,2,3}"), "(vec 1 2 3)");
        assert_eq!(p("{1,2}[2]"), "(idx (vec 1 2) [2])");
        assert_eq!(p("V[1,2]"), "(idx V [1 2])");
        assert_eq!(p("V[1+1]"), "(idx V [(Add 1 1)])");
        assert!(fails("V[1,2,3,4,5,6]"));
        assert!(fails("{1,2"));
        assert!(fails("V["));
    }

    #[test]
    fn calls() {
        assert_eq!(p("SIN(X)"), "(SIN X)");
        assert_eq!(p("RANDOM()"), "(RANDOM )");
        assert_eq!(p("MIN(X,3)"), "(MIN X 3)");
        assert_eq!(p("ANGSEP(1,2,3,4)"), "(ANGSEP 1 2 3 4)");
        assert!(fails("SIN(X"));
        assert!(fails("SIN(,)"));
    }

    #[test]
    fn rejects_trailing_and_missing_operands() {
        assert!(fails("X+"));
        assert!(fails("+"));
        assert!(fails("(X"));
        assert!(fails("X)"));
        assert!(fails("&&"));
        assert!(fails("?:"));
    }

    #[test]
    fn depth_is_bounded() {
        /* deep nesting is refused rather than overflowing the stack */
        let deep = format!("{}1{}", "(".repeat(500), ")".repeat(500));
        assert!(fails(&deep));
        assert!(fails(&format!("{}1", "!".repeat(500))));
        assert!(fails(&format!("{}1{}", "SIN(".repeat(500), ")".repeat(500))));

        /* ordinary nesting still works ... */
        let ok = format!("{}1{}", "(".repeat(50), ")".repeat(50));
        assert_eq!(p(&ok), "1");

        /* ... and left-associative chains cost no depth at all, so a very long
        sum parses fine */
        let chain = format!("{}\n", vec!["1"; 5000].join("+"));
        assert!(tokenize(chain.as_bytes()).is_ok_and(|t| parse(&t).is_ok()));
    }

    #[test]
    fn blank_input_yields_no_tree() {
        let toks = tokenize(b"\n").unwrap();
        assert_eq!(parse(&toks).unwrap(), None);
        let toks = tokenize(b"").unwrap();
        assert_eq!(parse(&toks).unwrap(), None);
        /* an expression with no terminator is a syntax error */
        let toks = tokenize(b"1").unwrap();
        assert!(parse(&toks).is_err());
    }

    #[test]
    fn last_line_wins() {
        let toks = tokenize(b"1\n2\n").unwrap();
        assert_eq!(show(&parse(&toks).unwrap().unwrap()), "2");
    }
}
