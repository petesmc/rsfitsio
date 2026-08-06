//! Lowering: syntax tree to `Node` arena.
//!
//! This is where the sort information that `eval.y` encoded in its four
//! nonterminals reappears. Each [`Ast`] node is lowered to an index into
//! `ParseData::Nodes`, and the *sort* of that index (`Nodes[i].ntype`) selects
//! which of the original productions applies — so the four `NELEM(`
//! implementations become one `match`, and the twenty-five subscript rules
//! become one.
//!
//! Everything here transcribes the semantic actions of `eval.y`; the `New_*`,
//! `Test_Dims` and `Copy_Dims` helpers it calls are the originals, unchanged,
//! in `crate::eval_y`.
//!
//! The nine constraints in `PARSER_SPEC.md` §3.4 that the stratified grammar
//! enforced by *not having a production* must be enforced explicitly here.
//! Each is marked with a `spec 3.4 (n)` comment.

use super::ast::{Ast, AstKind, BinOp, UnOp};
use super::error::ParseError;
use super::msg::unsupported;
use std::collections::HashMap;

use crate::region::SAORegion;

use super::resolve::Resolutions;
use super::token::CallKind;
use crate::c_types::{c_char, c_double, c_int, c_long};
use crate::eval_defs::{
    FuncOp, MAX_STRLEN, MAXSUBS, NodeValue, OpCode, Operation, ParseData, ParserValue, ValueSort,
};
use crate::eval_y::{
    Close_Vec, Copy_Dims, New_Array, New_BinOp, New_Column, New_Const, New_Deref, New_Func,
    New_FuncSize, New_GTI, New_Offset, New_REG, New_Unary, New_Vector, Test_Dims,
};

/// A node index in `ParseData::Nodes`.
type NodeId = c_int;
type LRes = Result<NodeId, ParseError>;

/// Sentinel the GTI/region builders use for an omitted node argument.
const NO_NODE: NodeId = -99;

/// True for the sorts that `eval.y` calls `expr`. `Boolean` is deliberately
/// excluded: a boolean-valued node is a `bexpr`, a *different* nonterminal, and
/// most arithmetic and comparison productions accept only `expr`. That is why
/// `INTCOL + BOOLCOL` and `1<2==1` are syntax errors.
fn is_expr(t: ValueSort) -> bool {
    t.is_expr()
}

/// `expr` or `bexpr` — the sorts the casts and vector literals accept.
fn is_numeric(t: ValueSort) -> bool {
    t == ValueSort::Boolean || is_expr(t)
}

/// Convert a byte slice to a NUL-terminated `c_char` buffer.
fn cstr(b: &[u8]) -> Vec<c_char> {
    let mut v: Vec<c_char> = b.iter().map(|&c| c as c_char).collect();
    v.push(0);
    v
}

/// The good-time intervals a `GTIFILTER`/`GTIFIND` call loaded, and whether
/// they came back in ascending order.
///
/// `New_GTI` reads these from the file while building the arena. Recording
/// them here lets the columnar evaluator use the same load rather than
/// repeating four hundred lines of HDU navigation, in the same way
/// [`Resolutions`] hands it the already-resolved column names.
#[derive(Clone, Debug, Default)]
pub(crate) struct GtiData {
    pub(crate) start: Vec<c_double>,
    pub(crate) stop: Vec<c_double>,
    pub(crate) ordered: bool,
    /// The column the times come from when the call did not name an
    /// expression -- `gtifilter()` means the TIME column, which `New_GTI`
    /// resolves and the source text does not mention.
    pub(crate) time_col: Option<c_int>,
}

/// The GTI loads of one expression, keyed by the call's byte offset.
pub(crate) type GtiLoads = HashMap<usize, GtiData>;

/// The regions of one expression, keyed by the call's byte offset.
pub(crate) type RegionLoads = HashMap<usize, RegionData>;

/// The region a `REGFILTER` call parsed, and the coordinate columns it used.
#[derive(Clone, Debug)]
pub(crate) struct RegionData {
    pub(crate) region: SAORegion,
    /// The X and Y columns `New_REG` resolved when the call did not name
    /// expressions, as `regfilter("f.reg")` does not.
    pub(crate) x_col: Option<c_int>,
    pub(crate) y_col: Option<c_int>,
}

pub(crate) struct Lowerer<'a> {
    pub(crate) p: &'a mut ParseData,
    /// Names already resolved by [`super::resolve`], keyed by byte offset.
    pub(crate) names: &'a Resolutions,
    /// What each GTI call read, for the columnar lowering that follows.
    pub(crate) gti: GtiLoads,
    /// And likewise the region each `REGFILTER` parsed.
    pub(crate) regions: RegionLoads,
}

impl Lowerer<'_> {
    // -----------------------------------------------------------------
    // Node accessors, mirroring eval.y's TYPE / OPER / SIZE macros
    // -----------------------------------------------------------------

    fn ntype(&self, n: NodeId) -> ValueSort {
        self.p.Nodes[n as usize].ntype
    }
    fn set_ntype(&mut self, n: NodeId, t: ValueSort) {
        self.p.Nodes[n as usize].ntype = t;
    }
    fn size(&self, n: NodeId) -> c_long {
        self.p.Nodes[n as usize].value.nelem
    }
    fn set_size(&mut self, n: NodeId, v: c_long) {
        self.p.Nodes[n as usize].value.nelem = v;
    }
    fn is_const(&self, n: NodeId) -> bool {
        self.p.Nodes[n as usize].operation == Operation::Const
    }
    /// The value of a node already known to be a `ValueSort::Long` constant.
    fn const_lng(&self, n: NodeId) -> c_long {
        self.p.Nodes[n as usize].value.data.lng()
    }

    /// `TEST(a)`: a negative index means the builder failed.
    fn test(&self, n: NodeId) -> LRes {
        if n < 0 {
            /* `#define TEST(a) if( (a)<0 ) YYERROR` -- silent. Whichever
            builder returned the negative index has already pushed its own
            message and set the status. */
            Err(ParseError::aborted())
        } else {
            Ok(n)
        }
    }

    /// `PROMOTE(a,b)`: raise whichever side is lower on ValueSort::Boolean < ValueSort::Long < ValueSort::Double.
    fn promote(&mut self, a: &mut NodeId, b: &mut NodeId) {
        let (ta, tb) = (self.ntype(*a), self.ntype(*b));
        if ta > tb {
            *b = New_Unary(self.p, ta, None, *b);
        } else if ta < tb {
            *a = New_Unary(self.p, tb, None, *a);
        }
    }

    fn dims_ok(&mut self, a: NodeId, b: NodeId) -> bool {
        Test_Dims(self.p, a, b) != 0
    }

    /// `if( SIZE(a) < SIZE(b) ) Copy_Dims(dst, b);`
    fn copy_dims_if(&mut self, dst: NodeId, a: NodeId, b: NodeId) {
        if self.size(a) < self.size(b) {
            Copy_Dims(self.p, dst, b);
        }
    }

    fn func(&mut self, ret: Option<ValueSort>, op: FuncOp, args: &[NodeId]) -> NodeId {
        let mut a = [0 as NodeId; 7];
        a[..args.len()].copy_from_slice(args);
        New_Func(
            self.p,
            ret,
            op,
            args.len() as c_int,
            a[0],
            a[1],
            a[2],
            a[3],
            a[4],
            a[5],
            a[6],
        )
    }

    fn func_size(
        &mut self,
        ret: Option<ValueSort>,
        op: FuncOp,
        args: &[NodeId],
        size: c_int,
    ) -> NodeId {
        let mut a = [0 as NodeId; 7];
        a[..args.len()].copy_from_slice(args);
        New_FuncSize(
            self.p,
            ret,
            op,
            args.len() as c_int,
            a[0],
            a[1],
            a[2],
            a[3],
            a[4],
            a[5],
            a[6],
            size,
        )
    }

    fn const_long(&mut self, v: c_long) -> NodeId {
        New_Const(self.p, ValueSort::Long, NodeValue::Long(v))
    }

    fn const_double(&mut self, v: c_double) -> NodeId {
        New_Const(self.p, ValueSort::Double, NodeValue::Double(v))
    }

    fn const_bool(&mut self, v: bool) -> NodeId {
        New_Const(
            self.p,
            ValueSort::Boolean,
            NodeValue::Logical(c_char::from(v)),
        )
    }

    /// A string or bit-string constant. `text` is the content without a NUL.
    fn const_text(&mut self, tag: ValueSort, text: &[c_char]) -> LRes {
        let mut buf = [0 as c_char; MAX_STRLEN as usize];
        let len = text.len().min(buf.len() - 1);
        buf[..len].copy_from_slice(&text[..len]);
        let n = New_Const(self.p, tag, NodeValue::Text(buf));
        let n = self.test(n)?;
        self.set_size(n, len as c_long);
        Ok(n)
    }

    /// `if( TYPE(n) != ValueSort::Double ) n = New_Unary(ValueSort::Double, 0, n);`
    fn as_double(&mut self, n: NodeId) -> NodeId {
        if self.ntype(n) != ValueSort::Double {
            New_Unary(self.p, ValueSort::Double, None, n)
        } else {
            n
        }
    }

    /// `if( TYPE(n) != ValueSort::Long ) n = New_Unary(ValueSort::Long, 0, n);`
    fn as_long(&mut self, n: NodeId) -> NodeId {
        if self.ntype(n) != ValueSort::Long {
            New_Unary(self.p, ValueSort::Long, None, n)
        } else {
            n
        }
    }

    // -----------------------------------------------------------------
    // Identifier resolution
    // -----------------------------------------------------------------

    /// The resolution recorded for the name token at `at`.
    ///
    /// The pre-pass resolved every name before parsing began and truncated the
    /// stream at the first failure, so a name that reaches lowering always has
    /// an entry.
    fn resolved(&self, name: &[u8], at: usize) -> Result<ParserValue, ParseError> {
        self.names
            .get(&at)
            .cloned()
            .ok_or_else(|| ParseError::syntax("unknown column or keyword", at, name))
    }

    /// Build the node the grammar would have built for a resolved name.
    fn resolve(&mut self, name: &[u8], at: usize) -> LRes {
        let n = match self.resolved(name, at)? {
            ParserValue::Column { index, .. } => New_Column(self.p, index),
            ParserValue::Long(v) => self.const_long(v),
            ParserValue::Double(v) => self.const_double(v),
            ParserValue::Boolean(v) => self.const_bool(v),
            ParserValue::Str(text) => {
                /* `Str` is NUL-terminated; the node size excludes it */
                return self.const_text(ValueSort::String, &text[..text.len() - 1]);
            }
        };
        self.test(n)
    }

    /// Resolve a name that must be a *column*, for the `NAME{offset}` form.
    fn resolve_column(&mut self, name: &[u8], at: usize) -> Result<c_int, ParseError> {
        /* spec 3.4 (1): the offset form exists only for columns, not for a
        keyword that happens to resolve to a constant */
        match self.resolved(name, at)? {
            ParserValue::Column { index, .. } => Ok(index),
            _ => Err(ParseError::syntax(
                "row offsets are only allowed on a column",
                at,
                name,
            )),
        }
    }

    // -----------------------------------------------------------------
    // Main dispatch
    // -----------------------------------------------------------------

    pub(crate) fn lower(&mut self, a: &Ast) -> LRes {
        match &a.kind {
            AstKind::Long(v) => {
                let v = *v;
                let n = self.const_long(v);
                self.test(n)
            }
            AstKind::Double(v) => {
                let v = *v;
                let n = self.const_double(v);
                self.test(n)
            }
            AstKind::Boolean(v) => {
                let n = self.const_bool(*v);
                self.test(n)
            }
            AstKind::Str(s) => self.new_text_const(s, ValueSort::String),
            AstKind::BitStr(s) => self.new_text_const(s, ValueSort::Bits),
            AstKind::RowRef => {
                let n = self.func(Some(ValueSort::Long), FuncOp::Row, &[]);
                self.test(n)
            }
            AstKind::NullRef => {
                let n = self.func(Some(ValueSort::Long), FuncOp::Null, &[]);
                self.test(n)
            }
            AstKind::SNullRef => {
                let n = self.func(Some(ValueSort::String), FuncOp::Null, &[]);
                self.test(n)
            }
            AstKind::Ident(name) => self.resolve(name, a.at),
            AstKind::Keyword(name) => self.resolve(name, a.at),
            AstKind::Offset { name, off } => {
                let col = self.resolve_column(name, a.at)?;
                let o = self.lower(off)?;
                /* `COLUMN '{' expr '}'` takes an `expr`, so a boolean or string
                offset matches no production; the action's own complaint below
                is only ever about a number that is the wrong kind of number */
                if !is_expr(self.ntype(o)) {
                    return Err(ParseError::syntax("offset must be numeric", a.at, b"{"));
                }
                if self.ntype(o) != ValueSort::Long || !self.is_const(o) {
                    return Err(ParseError::semantic(
                        "Offset argument must be a constant integer",
                    ));
                }
                let n = New_Offset(self.p, col, o);
                self.test(n)
            }
            AstKind::Unary { op, arg } => self.lower_unary(*op, arg, a.at),
            AstKind::Binary { op, lhs, rhs } => self.lower_binary(*op, lhs, rhs, a.at),
            AstKind::Ternary { cond, then, els } => self.lower_ternary(cond, then, els, a.at),
            AstKind::Range { val, lo, hi } => self.lower_range(val, lo, hi),
            AstKind::Deref { base, idx } => self.lower_deref(base, idx, a.at),
            AstKind::Vector(items) => self.lower_vector(items, a.at),
            AstKind::Call { kind, name, args } => self.lower_call(*kind, name, args, a.at),
        }
    }

    /// `New_Const` for a string or bit-string literal.
    fn new_text_const(&mut self, s: &[u8], tag: ValueSort) -> LRes {
        let text: Vec<c_char> = s.iter().map(|&c| c as c_char).collect();
        self.const_text(tag, &text)
    }

    fn lower_unary(&mut self, op: UnOp, arg: &Ast, at: usize) -> LRes {
        let n = self.lower(arg)?;
        let t = self.ntype(n);
        let r = match op {
            UnOp::Neg => {
                if !is_expr(t) {
                    return Err(ParseError::syntax("unary '-' needs a number", at, b"-"));
                }
                New_Unary(self.p, t, Some(OpCode::UMinus), n)
            }
            /* `'+' expr { $$ = $2; }` -- no node, but the sort still matters */
            UnOp::Plus => {
                if !is_expr(t) {
                    return Err(ParseError::syntax("unary '+' needs a number", at, b"+"));
                }
                return Ok(n);
            }
            /* spec 3.4 (3): `NOT bexpr` and `NOT bits` only */
            UnOp::Not => match t {
                ValueSort::Boolean => New_Unary(self.p, ValueSort::Boolean, Some(OpCode::Not), n),
                ValueSort::Bits => New_Unary(self.p, ValueSort::Bits, Some(OpCode::Not), n),
                _ => {
                    return Err(ParseError::syntax("syntax error", at, b"!"));
                }
            },
            UnOp::IntCast => {
                if !is_numeric(t) {
                    return Err(ParseError::syntax("(int) needs a number", at, b"(int)"));
                }
                New_Unary(self.p, ValueSort::Long, Some(OpCode::IntCast), n)
            }
            UnOp::FltCast => {
                if !is_numeric(t) {
                    return Err(ParseError::syntax("(float) needs a number", at, b"(float)"));
                }
                New_Unary(self.p, ValueSort::Double, Some(OpCode::FltCast), n)
            }
        };
        self.test(r)
    }

    /// `expr '=' expr ':' expr` — desugars to `(lo <= val) && (val <= hi)`.
    fn lower_range(&mut self, val: &Ast, lo: &Ast, hi: &Ast) -> LRes {
        let mut v = self.lower(val)?;
        let mut l = self.lower(lo)?;
        let mut h = self.lower(hi)?;
        for n in [v, l, h] {
            if !is_expr(self.ntype(n)) {
                /* `bexpr: expr '=' expr ':' expr` is the only range production,
                so a non-numeric operand matches nothing and bison reports it */
                return Err(ParseError::syntax(
                    "range test needs numeric operands",
                    0,
                    b"",
                ));
            }
        }
        self.promote(&mut v, &mut l);
        self.promote(&mut v, &mut h);
        self.promote(&mut l, &mut h);
        let lo_le = New_BinOp(self.p, ValueSort::Boolean, l, OpCode::Lte, v);
        let le_hi = New_BinOp(self.p, ValueSort::Boolean, v, OpCode::Lte, h);
        let r = New_BinOp(self.p, ValueSort::Boolean, lo_le, OpCode::And, le_hi);
        self.test(r)
    }

    fn lower_deref(&mut self, base: &Ast, idx: &[Ast], at: usize) -> LRes {
        let b = self.lower(base)?;
        /* spec 3.4 (2): expr, bexpr and bits subscript; sexpr does not */
        if self.ntype(b) == ValueSort::String {
            return Err(ParseError::syntax(
                "strings cannot be subscripted",
                at,
                b"[",
            ));
        }
        let mut dims = [0 as NodeId; 5];
        for (slot, e) in dims.iter_mut().zip(idx) {
            *slot = self.lower(e)?;
            /* every subscript production takes `expr` in the index positions,
            so a boolean or string index matches none of them and never reaches
            New_Deref -- which is what reports "Index value must be an integer
            type", and only ever sees the numbers that did parse */
            if !is_expr(self.ntype(*slot)) {
                return Err(ParseError::syntax(
                    "subscript index must be numeric",
                    at,
                    b"[",
                ));
            }
        }
        let r = New_Deref(
            self.p,
            b,
            idx.len() as c_int,
            dims[0],
            dims[1],
            dims[2],
            dims[3],
            dims[4],
        );
        self.test(r)
    }

    /// `{ e1, e2, … }`
    ///
    /// A literal is a *boolean* vector while every element so far has been
    /// boolean; the first numeric element promotes the whole thing, exactly as
    /// the `bvector` / `vector` productions do. Elements are chunked at
    /// `MAXSUBS`, closing and re-wrapping the node when it fills up.
    fn lower_vector(&mut self, items: &[Ast], at: usize) -> LRes {
        let mut all_bool = true;
        let mut vec_node: Option<NodeId> = None;

        for item in items {
            let e = self.lower(item)?;
            let te = self.ntype(e);
            if !is_numeric(te) {
                return Err(ParseError::syntax(
                    "vector elements must be numeric or boolean",
                    at,
                    b"{",
                ));
            }

            let Some(mut v) = vec_node else {
                let v = New_Vector(self.p, e);
                let v = self.test(v)?;
                all_bool = te == ValueSort::Boolean;
                vec_node = Some(v);
                continue;
            };

            if te != ValueSort::Boolean {
                if all_bool {
                    /* `bvector ',' expr`: the literal turns numeric */
                    self.set_ntype(v, te);
                    all_bool = false;
                } else if self.ntype(v) < te {
                    /* `vector ',' expr`: widen to the larger numeric type */
                    self.set_ntype(v, te);
                }
            }

            if self.p.Nodes[v as usize].nSubNodes >= MAXSUBS {
                let closed = Close_Vec(self.p, v);
                let closed = self.test(closed)?;
                let fresh = New_Vector(self.p, closed);
                v = self.test(fresh)?;
            }
            let k = self.p.Nodes[v as usize].nSubNodes as usize;
            self.p.Nodes[v as usize].SubNodes[k] = e as usize;
            self.p.Nodes[v as usize].nSubNodes += 1;
            vec_node = Some(v);
        }

        let Some(v) = vec_node else {
            return Err(ParseError::syntax("empty vector literal", at, b"{"));
        };
        let r = Close_Vec(self.p, v);
        self.test(r)
    }

    // -----------------------------------------------------------------
    // Binary operators
    // -----------------------------------------------------------------

    fn lower_binary(&mut self, op: BinOp, lhs: &Ast, rhs: &Ast, at: usize) -> LRes {
        let mut a = self.lower(lhs)?;
        let mut b = self.lower(rhs)?;
        let (ta, tb) = (self.ntype(a), self.ntype(b));

        let r = match op {
            BinOp::Add => match (ta, tb) {
                (ValueSort::String, ValueSort::String) | (ValueSort::Bits, ValueSort::Bits) => {
                    let total = self.size(a) + self.size(b);
                    if total >= c_long::from(MAX_STRLEN) {
                        return Err(ParseError::semantic(if ta == ValueSort::String {
                            "Combined string size exceeds 255 characters"
                        } else {
                            "Combined bit string size exceeds 255 bits"
                        }));
                    }
                    let n = New_BinOp(self.p, ta, a, OpCode::Add, b);
                    let n = self.test(n)?;
                    self.set_size(n, total);
                    return Ok(n);
                }
                _ if is_expr(ta) && is_expr(tb) => self.arith(a, b, OpCode::Add)?,
                _ => return Err(self.bad_operands("+", at)),
            },
            BinOp::Sub => self.numeric_only(a, b, OpCode::Sub, at)?,
            BinOp::Div => self.numeric_only(a, b, OpCode::Div, at)?,
            BinOp::Mod => self.numeric_only(a, b, OpCode::Mod, at)?,
            BinOp::Mul => match (ta, tb) {
                /* `expr '*' bexpr` and `bexpr '*' expr` coerce the boolean */
                (x, ValueSort::Boolean) if is_expr(x) => {
                    b = New_Unary(self.p, ta, None, b);
                    New_BinOp(self.p, ta, a, OpCode::Mul, b)
                }
                (ValueSort::Boolean, y) if is_expr(y) => {
                    a = New_Unary(self.p, tb, None, a);
                    New_BinOp(self.p, tb, a, OpCode::Mul, b)
                }
                _ => self.numeric_only(a, b, OpCode::Mul, at)?,
            },
            BinOp::Pow => {
                if !is_expr(ta) || !is_expr(tb) {
                    return Err(self.bad_operands("**", at));
                }
                self.promote(&mut a, &mut b);
                let t = self.ntype(a);
                New_BinOp(self.p, t, a, OpCode::Power, b)
            }

            /* `bits & bits` / `bits | bits`, or integer bitwise */
            BinOp::BitAnd | BinOp::BitOr => {
                let ch = if op == BinOp::BitAnd {
                    OpCode::BitAnd
                } else {
                    OpCode::BitOr
                };
                if ta == ValueSort::Bits && tb == ValueSort::Bits {
                    let n = New_BinOp(self.p, ValueSort::Bits, a, ch, b);
                    let n = self.test(n)?;
                    let sz = self.size(a).max(self.size(b));
                    self.set_size(n, sz);
                    return Ok(n);
                }
                /* The C grammar has only `bits OP bits` and `expr OP expr`, so
                a boolean, string or mixed pair matches no production at all and
                bison reports a bare syntax error before any action runs. The
                message below belongs to the numeric rule, which accepts the
                pair and then rejects it for not being integral -- so it is
                reached only when both sides are already numbers. */
                if !is_expr(ta) || !is_expr(tb) {
                    return Err(self.bad_operands(if op == BinOp::BitAnd { "&" } else { "|" }, at));
                }
                if ta != ValueSort::Long || tb != ValueSort::Long {
                    return Err(ParseError::semantic(
                        "Bitwise operations with incompatible types; only (bit OP bit) and (int OP int) are allowed",
                    ));
                }
                New_BinOp(self.p, ta, a, ch, b)
            }
            /* `bits` has `&`, `|` and `+` but no XOR, so unlike the two above
            this one has a numeric production and nothing else. */
            BinOp::BitXor => {
                if !is_expr(ta) || !is_expr(tb) {
                    return Err(self.bad_operands("^^", at));
                }
                if ta != ValueSort::Long || tb != ValueSort::Long {
                    return Err(ParseError::semantic(
                        "Bitwise operations with incompatible types; only (bit OP bit) and (int OP int) are allowed",
                    ));
                }
                New_BinOp(self.p, ta, a, OpCode::BitXor, b)
            }

            BinOp::And | BinOp::Or => {
                if ta != ValueSort::Boolean || tb != ValueSort::Boolean {
                    return Err(self.bad_operands(if op == BinOp::And { "&&" } else { "||" }, at));
                }
                let code = if op == BinOp::And {
                    OpCode::And
                } else {
                    OpCode::Or
                };
                New_BinOp(self.p, ValueSort::Boolean, a, code, b)
            }

            BinOp::Eq | BinOp::Ne | BinOp::Gt | BinOp::Lt | BinOp::Gte | BinOp::Lte => {
                return self.lower_compare(op, a, b, at);
            }
            BinOp::Approx => {
                if !is_expr(ta) || !is_expr(tb) {
                    return Err(self.bad_operands("~", at));
                }
                self.promote(&mut a, &mut b);
                New_BinOp(self.p, ValueSort::Boolean, a, OpCode::Approx, b)
            }
        };
        self.test(r)
    }

    /// The sorts of a call's arguments, in order.
    ///
    /// `unsupported` needs the whole list, not just the argument that failed:
    /// the message names the production that would have reduced, and that is
    /// chosen by the full signature.
    fn sorts(&self, ns: &[NodeId]) -> Vec<ValueSort> {
        ns.iter().map(|&n| self.ntype(n)).collect()
    }

    /// `FUNCTION expr …` productions accept only `expr` operands. Without this
    /// a string or bit-string argument reaches `Do_Func`, which reads the
    /// buffer through a `*long` and segfaults.
    fn require_expr(&self, n: NodeId, name: &[u8], all: &[ValueSort]) -> Result<(), ParseError> {
        if is_expr(self.ntype(n)) {
            Ok(())
        } else {
            Err(unsupported(name, all))
        }
    }

    /// `FUNCTION (expr|bexpr) …`: the first argument of AXISELEM / NAXES /
    /// ARRAY may also be boolean.
    fn require_numeric(&self, n: NodeId, name: &[u8], all: &[ValueSort]) -> Result<(), ParseError> {
        if is_numeric(self.ntype(n)) {
            Ok(())
        } else {
            Err(unsupported(name, all))
        }
    }

    /// The `expr` positions of the GTI and region calls.
    ///
    /// `GTIFILTER STRING ',' expr ')'` and its neighbours are productions in
    /// their own right, over their own tokens, and none of them carries a
    /// "not supported" wording -- so an argument of the wrong sort matches no
    /// production and bison reports a bare syntax error.
    fn require_expr_arg(&self, n: NodeId) -> Result<(), ParseError> {
        if is_expr(self.ntype(n)) {
            Ok(())
        } else {
            Err(ParseError::syntax("syntax error", 0, b""))
        }
    }

    /// Operands an operator is not defined for.
    ///
    /// CFITSIO's grammar is type-stratified -- `expr`, `bexpr`, `sexpr` and
    /// `bits` are separate nonterminals -- so it does not diagnose these as
    /// type errors at all: no production matches and bison reports `syntax
    /// error`. Checked against the real library with
    /// `ORACLE_MSGS=1 tests/oracle/oracle`, which answers `syntax error` to
    /// `1 && 2`, `'ab' > 1` and every other combination this rejects.
    ///
    /// The message a caller reads back with `fits_read_errmsg` is part of the
    /// library's behaviour, so it says what CFITSIO says. `op` and `at` still
    /// place the error.
    fn bad_operands(&self, op: &str, at: usize) -> ParseError {
        ParseError::syntax("syntax error", at, op.as_bytes())
    }

    /// `expr OP expr` with numeric promotion, the shape of most `eval.y`
    /// arithmetic actions.
    fn arith(&mut self, mut a: NodeId, mut b: NodeId, ch: OpCode) -> LRes {
        self.promote(&mut a, &mut b);
        let t = self.ntype(a);
        Ok(New_BinOp(self.p, t, a, ch, b))
    }

    fn numeric_only(&mut self, a: NodeId, b: NodeId, ch: OpCode, at: usize) -> LRes {
        if !is_expr(self.ntype(a)) || !is_expr(self.ntype(b)) {
            return Err(self.bad_operands(&format!("{ch:?}"), at));
        }
        self.arith(a, b, ch)
    }

    fn lower_compare(&mut self, op: BinOp, mut a: NodeId, mut b: NodeId, at: usize) -> LRes {
        let (ta, tb) = (self.ntype(a), self.ntype(b));
        let code = match op {
            BinOp::Eq => OpCode::Eq,
            BinOp::Ne => OpCode::Ne,
            BinOp::Gt => OpCode::Gt,
            BinOp::Lt => OpCode::Lt,
            BinOp::Gte => OpCode::Gte,
            BinOp::Lte => OpCode::Lte,
            _ => unreachable!(),
        };

        /* `bits CMP bits` and `sexpr CMP sexpr` force a scalar result */
        if (ta == ValueSort::Bits && tb == ValueSort::Bits)
            || (ta == ValueSort::String && tb == ValueSort::String)
        {
            let n = New_BinOp(self.p, ValueSort::Boolean, a, code, b);
            let n = self.test(n)?;
            self.set_size(n, 1);
            return Ok(n);
        }
        /* `bexpr EQ bexpr` and `bexpr NE bexpr` — no ordering on booleans */
        if ta == ValueSort::Boolean
            && tb == ValueSort::Boolean
            && matches!(op, BinOp::Eq | BinOp::Ne)
        {
            let n = New_BinOp(self.p, ValueSort::Boolean, a, code, b);
            return self.test(n);
        }
        if !is_expr(ta) || !is_expr(tb) {
            return Err(self.bad_operands("comparison", at));
        }
        self.promote(&mut a, &mut b);
        let n = New_BinOp(self.p, ValueSort::Boolean, a, code, b);
        self.test(n)
    }

    // -----------------------------------------------------------------
    // Conditional
    // -----------------------------------------------------------------

    fn lower_ternary(&mut self, cond: &Ast, then: &Ast, els: &Ast, at: usize) -> LRes {
        let c = self.lower(cond)?;
        if self.ntype(c) != ValueSort::Boolean {
            return Err(ParseError::syntax(
                "the condition of '?:' must be boolean",
                at,
                b"?",
            ));
        }
        let mut x = self.lower(then)?;
        let mut y = self.lower(els)?;
        let (tx, ty) = (self.ntype(x), self.ntype(y));

        /* spec 3.4 (5): there is no `bexpr '?' bits ':' bits` */
        if tx == ValueSort::Bits || ty == ValueSort::Bits {
            return Err(ParseError::syntax("syntax error", at, b"?"));
        }

        if tx == ValueSort::String || ty == ValueSort::String {
            if tx != ty {
                return Err(ParseError::syntax(
                    "'?:' branches have incompatible types",
                    at,
                    b"?",
                ));
            }
            if self.size(c) != 1 {
                return Err(ParseError::semantic("Cannot have a vector string column"));
            }
            /* the output size must be known up front to avoid an overflow */
            let out = self.size(x).max(self.size(y)) as c_int;
            let n = self.func_size(None, FuncOp::IfThenElse, &[x, y, c], out);
            let n = self.test(n)?;
            self.copy_dims_if(n, x, y);
            return Ok(n);
        }

        let both_bool = tx == ValueSort::Boolean && ty == ValueSort::Boolean;
        if !both_bool {
            self.promote(&mut x, &mut y);
        }
        if !self.dims_ok(x, y) {
            return Err(ParseError::semantic(
                "Incompatible dimensions in '?:' arguments",
            ));
        }
        let n = self.func(None, FuncOp::IfThenElse, &[x, y, c]);
        let n = self.test(n)?;
        self.copy_dims_if(n, x, y);

        /* the condition is temporarily retyped so that Test_Dims compares the
        shapes rather than the sorts, exactly as `eval.y` does */
        if both_bool {
            if !self.dims_ok(c, n) {
                return Err(ParseError::semantic(
                    "Incompatible dimensions in '?:' condition",
                ));
            }
        } else {
            let tx2 = self.ntype(x);
            self.set_ntype(c, tx2);
            let ok = self.dims_ok(c, n);
            self.set_ntype(c, ValueSort::Boolean);
            if !ok {
                return Err(ParseError::semantic(
                    "Incompatible dimensions in '?:' condition",
                ));
            }
        }
        if self.size(n) < self.size(c) {
            Copy_Dims(self.p, n, c);
        }
        Ok(n)
    }

    // -----------------------------------------------------------------
    // Calls
    // -----------------------------------------------------------------

    fn lower_call(&mut self, kind: CallKind, name: &[u8], args: &[Ast], at: usize) -> LRes {
        match kind {
            CallKind::GtiFilter => self.lower_gti(FuncOp::GtiFilt, name, args, at),
            CallKind::GtiFind => self.lower_gti(FuncOp::GtiFind, name, args, at),
            CallKind::GtiOverlap => self.lower_gti_overlap(name, args, at),
            CallKind::RegFilter => self.lower_regfilter(name, args, at),
            CallKind::IFunction => {
                /* STRSTR(sexpr, sexpr) */
                let [s, sub] = self.lower_n::<2>(args, name, at)?;
                if self.ntype(s) != ValueSort::String || self.ntype(sub) != ValueSort::String {
                    return Err(unsupported(name, &[self.ntype(s), self.ntype(sub)]));
                }
                let n = self.func(Some(ValueSort::Long), FuncOp::StrPos, &[s, sub]);
                self.test(n)
            }
            CallKind::BFunction => self.lower_bfunction(name, args, at),
            CallKind::Function => self.lower_function(name, args, at),
        }
    }

    /// Lower exactly `N` arguments, or report the arity mismatch.
    fn lower_n<const N: usize>(
        &mut self,
        args: &[Ast],
        name: &[u8],
        at: usize,
    ) -> Result<[NodeId; N], ParseError> {
        if args.len() != N {
            return Err(ParseError::syntax(
                format!(
                    "{}() takes {N} argument(s), got {}",
                    String::from_utf8_lossy(name),
                    args.len()
                ),
                at,
                name,
            ));
        }
        let mut out = [0 as NodeId; N];
        for (slot, a) in out.iter_mut().zip(args) {
            *slot = self.lower(a)?;
        }
        Ok(out)
    }

    fn lower_all(&mut self, args: &[Ast]) -> Result<Vec<NodeId>, ParseError> {
        args.iter().map(|a| self.lower(a)).collect()
    }

    fn lower_bfunction(&mut self, name: &[u8], args: &[Ast], at: usize) -> LRes {
        let ns = self.lower_all(args)?;
        match (name, ns.len()) {
            (b"ISNULL", 1) => {
                let a = ns[0];
                match self.ntype(a) {
                    ValueSort::String => {
                        let n = self.func(Some(ValueSort::Boolean), FuncOp::IsNull, &[a]);
                        self.test(n)
                    }
                    /* spec 3.4 (7): ISNULL is not defined for bit strings */
                    ValueSort::Bits => Err(unsupported(name, &[ValueSort::Bits])),
                    _ => {
                        let n = self.func(None, FuncOp::IsNull, &[a]);
                        let n = self.test(n)?;
                        /* keep the argument's size but return ValueSort::Boolean */
                        self.set_ntype(n, ValueSort::Boolean);
                        Ok(n)
                    }
                }
            }
            (b"NEAR", 3) => self.region_fct(FuncOp::Near, &ns, name, "NEAR"),
            (b"CIRCLE", 5) => self.region_fct(FuncOp::Circle, &ns, name, "CIRCLE"),
            (b"BOX", 7) => self.region_fct(FuncOp::Box, &ns, name, "BOX or ELLIPSE"),
            (b"ELLIPSE", 7) => self.region_fct(FuncOp::Ellipse, &ns, name, "BOX or ELLIPSE"),
            _ => {
                let ts = self.sorts(&ns);
                Err(unsupported(name, &ts))
            }
        }
    }

    /// `NEAR`, `CIRCLE`, `BOX` and `ELLIPSE`: all arguments to ValueSort::Double, all
    /// dimensions pairwise compatible, then propagate the largest shape.
    fn region_fct(&mut self, op: FuncOp, ns: &[NodeId], name: &[u8], what: &str) -> LRes {
        let ts = self.sorts(ns);
        for &n in ns {
            self.require_expr(n, name, &ts)?;
        }
        let ds: Vec<NodeId> = ns.iter().map(|&n| self.as_double(n)).collect();
        for w in ds.windows(2) {
            if !self.dims_ok(w[0], w[1]) {
                return Err(ParseError::semantic(format!(
                    "Dimensions of {what} arguments are not compatible"
                )));
            }
        }
        let n = self.func(Some(ValueSort::Boolean), op, &ds);
        let n = self.test(n)?;
        self.copy_dims_if(n, n, ds[0]);
        for w in ds.windows(2) {
            self.copy_dims_if(n, w[0], w[1]);
        }
        Ok(n)
    }

    fn lower_function(&mut self, name: &[u8], args: &[Ast], at: usize) -> LRes {
        let ns = self.lower_all(args)?;
        match ns.len() {
            0 => match name {
                b"RANDOM" => {
                    let n = self.func(Some(ValueSort::Double), FuncOp::Rnd, &[]);
                    self.test(n)
                }
                b"RANDOMN" => {
                    let n = self.func(Some(ValueSort::Double), FuncOp::GasRnd, &[]);
                    self.test(n)
                }
                _ => Err(unsupported(name, &[])),
            },
            1 => self.func1(name, ns[0]),
            2 => self.func2(name, ns[0], ns[1]),
            3 => self.func3(name, ns[0], ns[1], ns[2]),
            4 => self.func4(name, &ns),
            _ => {
                let ts = self.sorts(&ns);
                Err(unsupported(name, &ts))
            }
        }
    }

    fn func1(&mut self, name: &[u8], a: NodeId) -> LRes {
        let t = self.ntype(a);
        let sz = self.size(a);

        /* forms that exist for every argument sort */
        if name == b"NELEM" {
            let n = self.const_long(sz);
            return self.test(n);
        }

        match t {
            ValueSort::Boolean => match name {
                b"SUM" => {
                    let n = self.func(Some(ValueSort::Long), FuncOp::Sum, &[a]);
                    self.test(n)
                }
                b"ACCUM" => self.accum(a, ValueSort::Long, OpCode::Accum),
                /* `FUNCTION bexpr ')'` supports only SUM, NELEM and ACCUM.
                Falling through to the numeric list would accept `ABS(BOOLCOL)`
                and then hand Do_Func a char buffer to read through a *long. */
                _ => Err(unsupported(name, &[t])),
            },
            ValueSort::String => match name {
                b"NVALID" => {
                    let n = self.func(Some(ValueSort::Long), FuncOp::NonNull, &[a]);
                    self.test(n)
                }
                _ => Err(unsupported(name, &[t])),
            },
            ValueSort::Bits => match name {
                /* bit arrays have no NULLs, so NVALID is the element count */
                b"NVALID" => {
                    let n = self.const_long(sz);
                    self.test(n)
                }
                b"SUM" => {
                    let n = self.func(Some(ValueSort::Long), FuncOp::Sum, &[a]);
                    self.test(n)
                }
                b"MIN" | b"MAX" => {
                    let op = if name == b"MIN" {
                        FuncOp::Min1
                    } else {
                        FuncOp::Max1
                    };
                    let n = self.func(Some(t), op, &[a]);
                    let n = self.test(n)?;
                    /* a is a vector, so the result is never constant and it is
                    safe to force the size */
                    self.set_size(n, 1);
                    Ok(n)
                }
                b"ACCUM" => self.accum(a, ValueSort::Long, OpCode::Accum),
                _ => Err(unsupported(name, &[t])),
            },
            _ => self.func1_numeric(name, a),
        }
    }

    fn accum(&mut self, a: NodeId, ret: ValueSort, code: OpCode) -> LRes {
        let zero = if ret == ValueSort::Double {
            self.const_double(0.0)
        } else {
            self.const_long(0)
        };
        let n = New_BinOp(self.p, ret, a, code, zero);
        self.test(n)
    }

    fn func1_numeric(&mut self, name: &[u8], a: NodeId) -> LRes {
        let t = self.ntype(a);
        let sz = self.size(a);
        let n = match name {
            b"SUM" => self.func(Some(t), FuncOp::Sum, &[a]),
            b"AVERAGE" => self.func(Some(ValueSort::Double), FuncOp::Average, &[a]),
            b"STDDEV" => self.func(Some(ValueSort::Double), FuncOp::Stddev, &[a]),
            b"MEDIAN" => self.func(Some(t), FuncOp::Median, &[a]),
            b"NVALID" => self.func(Some(ValueSort::Long), FuncOp::NonNull, &[a]),
            b"ACCUM" if t == ValueSort::Long => {
                return self.accum(a, ValueSort::Long, OpCode::Accum);
            }
            b"ACCUM" if t == ValueSort::Double => {
                return self.accum(a, ValueSort::Double, OpCode::Accum);
            }
            b"SEQDIFF" if t == ValueSort::Long => {
                return self.accum(a, ValueSort::Long, OpCode::Diff);
            }
            b"SEQDIFF" if t == ValueSort::Double => {
                return self.accum(a, ValueSort::Double, OpCode::Diff);
            }
            b"ABS" => self.func(None, FuncOp::Abs, &[a]),
            b"MIN" => self.func(Some(t), FuncOp::Min1, &[a]),
            b"MAX" => self.func(Some(t), FuncOp::Max1, &[a]),
            b"RANDOM" | b"RANDOMN" => {
                let op = if name == b"RANDOM" {
                    FuncOp::Rnd
                } else {
                    FuncOp::GasRnd
                };
                let n = self.func(None, op, &[a]);
                let n = self.test(n)?;
                self.set_ntype(n, ValueSort::Double);
                return Ok(n);
            }
            b"ELEMENTNUM" => {
                if self.is_const(a) {
                    let n = self.const_long(1);
                    return self.test(n);
                }
                let n = self.func(None, FuncOp::ElemNum, &[a]);
                let n = self.test(n)?;
                self.set_ntype(n, ValueSort::Long);
                return Ok(n);
            }
            b"NAXIS" => {
                /* constant-folded: a constant is always one-dimensional */
                let v = if self.is_const(a) {
                    1
                } else {
                    c_long::from(self.p.Nodes[a as usize].value.naxis)
                };
                let n = self.const_long(v);
                return self.test(n);
            }
            _ => {
                /* everything else takes a ValueSort::Double */
                let d = self.as_double(a);
                let op = match name {
                    b"SIN" => FuncOp::Sin,
                    b"COS" => FuncOp::Cos,
                    b"TAN" => FuncOp::Tan,
                    b"ARCSIN" | b"ASIN" => FuncOp::Asin,
                    b"ARCCOS" | b"ACOS" => FuncOp::Acos,
                    b"ARCTAN" | b"ATAN" => FuncOp::Atan,
                    b"SINH" => FuncOp::Sinh,
                    b"COSH" => FuncOp::Cosh,
                    b"TANH" => FuncOp::Tanh,
                    b"EXP" => FuncOp::Exp,
                    b"LOG" => FuncOp::Log,
                    b"LOG10" => FuncOp::Log10,
                    b"SQRT" => FuncOp::Sqrt,
                    b"ROUND" => FuncOp::Round,
                    b"FLOOR" => FuncOp::Floor,
                    b"CEIL" => FuncOp::Ceil,
                    b"RANDOMP" => {
                        let n = self.func(None, FuncOp::PoiRnd, &[d]);
                        let n = self.test(n)?;
                        self.set_ntype(n, ValueSort::Long);
                        return Ok(n);
                    }
                    _ => return Err(unsupported(name, &[t])),
                };
                self.func(None, op, &[d])
            }
        };
        let _ = sz;
        self.test(n)
    }

    fn func2(&mut self, name: &[u8], a: NodeId, b: NodeId) -> LRes {
        let (ta, tb) = (self.ntype(a), self.ntype(b));

        match name {
            b"DEFNULL" => {
                /* spec 3.4 (8): no DEFNULL over bit strings */
                if ta == ValueSort::Bits || tb == ValueSort::Bits {
                    /* CFITSIO's grammar has no `FUNCTION bits ',' bits`
                    production at all, so it reports a syntax error */
                    return Err(ParseError::semantic("syntax error"));
                }
                if ta == ValueSort::String || tb == ValueSort::String {
                    if ta != tb {
                        return Err(unsupported(name, &[ta, tb]));
                    }
                    let out = self.size(a).max(self.size(b)) as c_int;
                    let n = self.func_size(None, FuncOp::DefNull, &[a, b], out);
                    let n = self.test(n)?;
                    if self.size(b) > self.size(a) {
                        let sb = self.size(b);
                        self.set_size(n, sb);
                    }
                    return Ok(n);
                }
                self.require_numeric(a, name, &[ta, tb])?;
                self.require_numeric(b, name, &[ta, tb])?;
                /* `FUNCTION expr ',' expr` and `FUNCTION bexpr ',' bexpr` are
                separate productions; there is no mixed form */
                if (ta == ValueSort::Boolean) != (tb == ValueSort::Boolean) {
                    return Err(unsupported(name, &[ta, tb]));
                }
                if !(self.size(a) >= self.size(b) && self.dims_ok(a, b)) {
                    return Err(ParseError::semantic(
                        "Dimensions of DEFNULL arguments are not compatible",
                    ));
                }
                let (mut a, mut b) = (a, b);
                if !(ta == ValueSort::Boolean && tb == ValueSort::Boolean) {
                    self.promote(&mut a, &mut b);
                }
                let n = self.func(None, FuncOp::DefNull, &[a, b]);
                self.test(n)
            }
            b"ARCTAN2" => {
                self.require_expr(a, name, &[ta, tb])?;
                self.require_expr(b, name, &[ta, tb])?;
                let a = self.as_double(a);
                let b = self.as_double(b);
                if !self.dims_ok(a, b) {
                    return Err(ParseError::semantic(
                        "Dimensions of arctan2 arguments are not compatible",
                    ));
                }
                let n = self.func(None, FuncOp::Atan2, &[a, b]);
                let n = self.test(n)?;
                self.copy_dims_if(n, a, b);
                Ok(n)
            }
            b"MIN" | b"MAX" => {
                self.require_expr(a, name, &[ta, tb])?;
                self.require_expr(b, name, &[ta, tb])?;
                let (mut a, mut b) = (a, b);
                self.promote(&mut a, &mut b);
                if !self.dims_ok(a, b) {
                    return Err(ParseError::semantic(format!(
                        "Dimensions of {}(a,b) arguments are not compatible",
                        String::from_utf8_lossy(name).to_lowercase()
                    )));
                }
                let op = if name == b"MIN" {
                    FuncOp::Min2
                } else {
                    FuncOp::Max2
                };
                let n = self.func(None, op, &[a, b]);
                let n = self.test(n)?;
                self.copy_dims_if(n, a, b);
                Ok(n)
            }
            b"SETNULL" => {
                self.require_expr(a, name, &[ta, tb])?;
                self.require_expr(b, name, &[ta, tb])?;
                if !self.is_const(a) || self.size(a) != 1 {
                    return Err(ParseError::semantic(
                        "SETNULL first argument must be a scalar constant",
                    ));
                }
                /* make the sentinel the same sort as the value */
                let a = if ta != tb {
                    New_Unary(self.p, tb, None, a)
                } else {
                    a
                };
                let n = self.func(None, FuncOp::SetNull, &[b, a]);
                self.test(n)
            }
            b"AXISELEM" => {
                self.require_numeric(a, name, &[ta, tb])?;
                self.require_expr(b, name, &[ta, tb])?;
                if !self.is_const(b) || self.size(b) != 1 {
                    return Err(ParseError::semantic(
                        "AXISELEM second argument must be a scalar constant",
                    ));
                }
                if self.is_const(a) {
                    let n = self.const_long(1);
                    return self.test(n);
                }
                let b = self.as_long(b);
                let n = self.func(None, FuncOp::AxisElem, &[a, b]);
                let n = self.test(n)?;
                self.set_ntype(n, ValueSort::Long);
                Ok(n)
            }
            b"NAXES" => {
                self.require_numeric(a, name, &[ta, tb])?;
                self.require_expr(b, name, &[ta, tb])?;
                if !self.is_const(b) || self.size(b) != 1 {
                    return Err(ParseError::semantic(
                        "NAXES second argument must be a scalar constant",
                    ));
                }
                if self.is_const(a) {
                    let n = self.const_long(1);
                    return self.test(n);
                }
                let b = self.as_long(b);
                /* both are constant now, so fold the axis length directly */
                let iaxis = self.const_lng(b);
                let naxis = c_long::from(self.p.Nodes[a as usize].value.naxis);
                let v = if iaxis == 0 {
                    naxis
                } else if iaxis <= naxis {
                    self.p.Nodes[a as usize].value.naxes[(iaxis - 1) as usize]
                } else {
                    1
                };
                let n = self.const_long(v);
                self.test(n)
            }
            b"ARRAY" => {
                self.require_numeric(a, name, &[ta, tb])?;
                self.require_expr(b, name, &[ta, tb])?;
                let n = New_Array(self.p, a, b);
                self.test(n)
            }
            _ => Err(unsupported(name, &[ta, tb])),
        }
    }

    fn func3(&mut self, name: &[u8], s: NodeId, pos: NodeId, len: NodeId) -> LRes {
        if name != b"STRMID" {
            return Err(unsupported(name, &self.sorts(&[s, pos, len])));
        }
        if self.ntype(s) != ValueSort::String {
            return Err(unsupported(name, &self.sorts(&[s, pos, len])));
        }
        if self.ntype(pos) != ValueSort::Long
            || self.size(pos) != 1
            || self.ntype(len) != ValueSort::Long
            || self.size(len) != 1
        {
            return Err(ParseError::semantic(
                "When using STRMID(S,P,N), P and N must be integers (and not vector columns)",
            ));
        }
        let n_chars = if self.is_const(len) {
            self.const_lng(len)
        } else {
            /* variable length: reserve the largest the source could give */
            self.size(s)
        };
        if n_chars <= 0 || n_chars >= c_long::from(MAX_STRLEN) {
            return Err(ParseError::semantic("STRMID(S,P,N), N must be 1-255"));
        }
        let n = self.func_size(None, FuncOp::StrMid, &[s, pos, len], n_chars as c_int);
        self.test(n)
    }

    fn func4(&mut self, name: &[u8], ns: &[NodeId]) -> LRes {
        if name != b"ANGSEP" {
            return Err(unsupported(name, &self.sorts(ns)));
        }
        let ts = self.sorts(ns);
        for &n in ns {
            self.require_expr(n, name, &ts)?;
        }
        let ds: Vec<NodeId> = ns.iter().map(|&n| self.as_double(n)).collect();
        for w in ds.windows(2) {
            if !self.dims_ok(w[0], w[1]) {
                return Err(ParseError::semantic(
                    "Dimensions of ANGSEP arguments are not compatible",
                ));
            }
        }
        let n = self.func(None, FuncOp::AngSep, &ds);
        let n = self.test(n)?;
        for w in ds.windows(2) {
            self.copy_dims_if(n, w[0], w[1]);
        }
        Ok(n)
    }

    // -----------------------------------------------------------------
    // GTI and region filters
    // -----------------------------------------------------------------

    /// spec 3.4 (4): the file name must be a literal string, not a string
    /// expression — `GTIFILTER(STRCOL)` is a syntax error.
    fn literal_str(&self, a: &Ast, name: &[u8]) -> Result<Vec<c_char>, ParseError> {
        match &a.kind {
            AstKind::Str(s) => Ok(cstr(s)),
            _ => Err(ParseError::syntax(
                format!(
                    "{}() needs a literal string here",
                    String::from_utf8_lossy(name)
                ),
                a.at,
                name,
            )),
        }
    }

    fn lower_gti(&mut self, op: FuncOp, name: &[u8], args: &[Ast], at: usize) -> LRes {
        let mut fname = cstr(b"");
        let mut start = cstr(b"*START*");
        let mut stop = cstr(b"*STOP*");
        let mut node1 = NO_NODE;

        match args.len() {
            0 => {}
            1 => fname = self.literal_str(&args[0], name)?,
            2 => {
                fname = self.literal_str(&args[0], name)?;
                node1 = self.lower(&args[1])?;
                self.require_expr_arg(node1)?;
            }
            4 => {
                fname = self.literal_str(&args[0], name)?;
                node1 = self.lower(&args[1])?;
                self.require_expr_arg(node1)?;
                start = self.literal_str(&args[2], name)?;
                stop = self.literal_str(&args[3], name)?;
            }
            n => {
                return Err(ParseError::syntax(
                    format!(
                        "{}() takes 0, 1, 2 or 4 arguments, got {n}",
                        String::from_utf8_lossy(name)
                    ),
                    at,
                    name,
                ));
            }
        }
        let n = New_GTI(
            self.p,
            op,
            fname.as_mut_ptr(),
            node1,
            NO_NODE,
            start.as_mut_ptr(),
            stop.as_mut_ptr(),
        );
        let n = self.test(n)?;
        self.record_gti(n, at);
        Ok(n)
    }

    /// Copy the region `New_REG` just parsed out of the node holding it.
    ///
    /// The node owns a raw `SAORegion` for the arena's own use, so this takes
    /// a clone rather than the pointer: two owners of one leaked box is not
    /// worth the saving.
    fn record_region(&mut self, n: NodeId, at: usize) {
        let held = &self.p.Nodes[self.p.Nodes[n as usize].SubNodes[0]];
        let raw = held.value.data.raw().cast::<SAORegion>();
        if raw.is_null() {
            return;
        }
        let region = unsafe { (*raw).clone() };
        let node = &self.p.Nodes[n as usize];
        let x_col = self.p.Nodes[node.SubNodes[1]].operation.column();
        let y_col = self.p.Nodes[node.SubNodes[2]].operation.column();
        self.regions.insert(
            at,
            RegionData {
                region,
                x_col,
                y_col,
            },
        );
    }

    /// Copy the intervals `New_GTI` just loaded out of the node holding them.
    ///
    /// They live in the first subnode as one buffer of `2 * nGTI` doubles:
    /// the starts, then the stops.
    fn record_gti(&mut self, n: NodeId, at: usize) {
        let node = &self.p.Nodes[n as usize];
        let times = node.SubNodes[0];
        let held = &self.p.Nodes[times];
        let count = held.value.nelem.max(0) as usize;
        let time_col0 = self.p.Nodes[node.SubNodes[1]].operation.column();
        if count == 0 {
            self.gti.insert(
                at,
                GtiData {
                    time_col: time_col0,
                    ..Default::default()
                },
            );
            return;
        }
        let ordered = held.gti_ordered;
        /* the time subnode, for the forms that did not spell one out */
        let time_col = self.p.Nodes[node.SubNodes[1]].operation.column();
        let buf = unsafe {
            let p = held.value.data.dbl_buf();
            if p.is_null() {
                &[][..]
            } else {
                core::slice::from_raw_parts(p, count * 2)
            }
        };
        self.gti.insert(
            at,
            GtiData {
                start: buf[..count].to_vec(),
                stop: buf[count..].to_vec(),
                ordered,
                time_col,
            },
        );
    }

    fn lower_gti_overlap(&mut self, name: &[u8], args: &[Ast], at: usize) -> LRes {
        let (mut start, mut stop) = (cstr(b"*START*"), cstr(b"*STOP*"));
        if args.len() != 3 && args.len() != 5 {
            return Err(ParseError::syntax(
                format!(
                    "{}() takes 3 or 5 arguments, got {}",
                    String::from_utf8_lossy(name),
                    args.len()
                ),
                at,
                name,
            ));
        }
        let mut fname = self.literal_str(&args[0], name)?;
        let n1 = self.lower(&args[1])?;
        let n2 = self.lower(&args[2])?;
        self.require_expr_arg(n1)?;
        self.require_expr_arg(n2)?;
        if args.len() == 5 {
            start = self.literal_str(&args[3], name)?;
            stop = self.literal_str(&args[4], name)?;
        }
        let n = New_GTI(
            self.p,
            FuncOp::GtiOver,
            fname.as_mut_ptr(),
            n1,
            n2,
            start.as_mut_ptr(),
            stop.as_mut_ptr(),
        );
        let n = self.test(n)?;
        self.record_gti(n, at);
        Ok(n)
    }

    fn lower_regfilter(&mut self, name: &[u8], args: &[Ast], at: usize) -> LRes {
        let mut cols = cstr(b"");
        let (mut nx, mut ny) = (NO_NODE, NO_NODE);
        if args.is_empty() || args.len() == 2 || args.len() > 4 {
            return Err(ParseError::syntax(
                format!(
                    "{}() takes 1, 3 or 4 arguments, got {}",
                    String::from_utf8_lossy(name),
                    args.len()
                ),
                at,
                name,
            ));
        }
        let mut fname = self.literal_str(&args[0], name)?;
        if args.len() >= 3 {
            nx = self.lower(&args[1])?;
            ny = self.lower(&args[2])?;
            self.require_expr_arg(nx)?;
            self.require_expr_arg(ny)?;
        }
        if args.len() == 4 {
            cols = self.literal_str(&args[3], name)?;
        }
        let n = New_REG(self.p, fname.as_mut_ptr(), nx, ny, cols.as_mut_ptr());
        let n = self.test(n)?;
        self.record_region(n, at);
        Ok(n)
    }
}
