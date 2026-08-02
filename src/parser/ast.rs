//! Untyped syntax tree.
//!
//! `eval.y` stratifies its grammar into four mutually-recursive nonterminals —
//! `expr`, `bexpr`, `sexpr`, `bits` — so that the LALR machine can pick the
//! right semantic action from the operand *sorts*. That is what makes the
//! grammar conflict-free, and it is also why there are 135 productions for a
//! language with about 30 distinct constructs.
//!
//! This tree is deliberately sort-agnostic: `a + b` is one [`AstKind::Binary`]
//! whatever `a` and `b` turn out to be. Every type-directed decision — which of
//! the three `+` rules applies, whether a `{…}` literal is boolean or numeric,
//! which overload of `NELEM(` is meant — is made in `super::lower`.

use super::token::CallKind;
use crate::c_types::c_long;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum UnOp {
    /// Unary `-`
    Neg,
    /// `!`, `.not.`
    Not,
    /// `(int)`
    IntCast,
    /// `(float)`, `(double)`
    FltCast,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum BinOp {
    Add,
    Sub,
    Mul,
    Div,
    Mod,
    Pow,
    /// `&` — bitwise on integers, AND on bit strings
    BitAnd,
    /// `|` — bitwise on integers, OR on bit strings
    BitOr,
    /// `^^`, `.xor.`
    BitXor,
    And,
    Or,
    Eq,
    Ne,
    Gt,
    Lt,
    Gte,
    Lte,
    /// `~` — approximate equality
    Approx,
}

#[derive(Clone, Debug, PartialEq)]
pub(crate) enum AstKind {
    Long(c_long),
    Double(f64),
    Boolean(bool),
    Str(Vec<u8>),
    BitStr(Vec<u8>),
    RowRef,
    NullRef,
    SNullRef,
    /// `#KEYWORD`, without the `#`.
    Keyword(Vec<u8>),
    /// A column or variable name.
    Ident(Vec<u8>),
    /// `NAME{offset}`. The base is a bare name by construction: the parser only
    /// builds this immediately after an identifier, which is how
    /// `(INTCOL){1}` stays a syntax error (`PARSER_SPEC.md` §3.4 item 1).
    Offset { name: Vec<u8>, off: Box<Ast> },
    Unary {
        op: UnOp,
        arg: Box<Ast>,
    },
    Binary {
        op: BinOp,
        lhs: Box<Ast>,
        rhs: Box<Ast>,
    },
    /// `cond ? then : els`
    Ternary {
        cond: Box<Ast>,
        then: Box<Ast>,
        els: Box<Ast>,
    },
    /// `val = lo : hi`, the range test.
    Range {
        val: Box<Ast>,
        lo: Box<Ast>,
        hi: Box<Ast>,
    },
    /// `base[i1, … i5]`, 1 to 5 subscripts.
    Deref { base: Box<Ast>, idx: Vec<Ast> },
    /// `{ e1, e2, … }`
    Vector(Vec<Ast>),
    Call {
        kind: CallKind,
        /// Upper-cased, without the `(`.
        name: Vec<u8>,
        args: Vec<Ast>,
    },
}

#[derive(Clone, Debug, PartialEq)]
pub(crate) struct Ast {
    pub(crate) kind: AstKind,
    /// Byte offset of the construct in the source expression.
    pub(crate) at: usize,
}

impl Ast {
    pub(crate) fn new(kind: AstKind, at: usize) -> Self {
        Ast { kind, at }
    }

    pub(crate) fn boxed(kind: AstKind, at: usize) -> Box<Self> {
        Box::new(Ast::new(kind, at))
    }

    /// True for the literal forms that `New_Const` turns into constant nodes.
    /// Used by the lowering pass where `eval.y` tests `OPER(x) == CONST_OP`.
    pub(crate) fn is_literal(&self) -> bool {
        matches!(
            self.kind,
            AstKind::Long(_) | AstKind::Double(_) | AstKind::Boolean(_) | AstKind::Str(_)
        )
    }
}
