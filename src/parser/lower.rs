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
use super::resolve::{Resolutions, token_of};
use super::token::CallKind;
use crate::c_types::{c_char, c_double, c_int, c_long, c_void};
use crate::eval_defs::{CONST_OP, MAX_STRLEN, MAXSUBS, ParseData};
use crate::eval_tab::{FITS_PARSER_YYSTYPE, fits_parser_yytokentype as T};
use crate::eval_y::{
    ABS_FCT, ACOS_FCT, ANGSEP_FCT, ASIN_FCT, ATAN2_FCT, ATAN_FCT, AVERAGE_FCT,
    AXISELEM_FCT, BOX_FCT, CEIL_FCT, CIRCLE_FCT, COS_FCT, COSH_FCT, Close_Vec, Copy_Dims,
    DEFNULL_FCT, ELEMNUM_FCT, ELPS_FCT, EXP_FCT, FLOOR_FCT, GASRND_FCT, GTIFILT_FCT, GTIFIND_FCT,
    GTIOVER_FCT, IFTHENELSE_FCT, ISNULL_FCT, LOG10_FCT, LOG_FCT, MAX1_FCT, MAX2_FCT, MEDIAN_FCT,
    MIN1_FCT, MIN2_FCT, NEAR_FCT, NONNULL_FCT, NULL_FCT, New_Array, New_BinOp, New_Column,
    New_Const, New_Deref, New_Func, New_FuncSize, New_GTI, New_Offset, New_REG, New_Unary,
    New_Vector, POIRND_FCT, RND_FCT, ROUND_FCT, ROW_FCT, SETNULL_FCT, SIN_FCT, SINH_FCT, SQRT_FCT,
    STDDEV_FCT, STRMID_FCT, STRPOS_FCT, SUM_FCT, TAN_FCT, TANH_FCT, Test_Dims,
    funcOp,
};

/// A node index in `ParseData::Nodes`.
type NodeId = c_int;
type LRes = Result<NodeId, ParseError>;

/* Sort tags. These double as token ids and as `Node::ntype`, and the first
three are ordered so that numeric promotion is a `>` comparison. */
const BOOLEAN: c_int = T::BOOLEAN as c_int;
const LONG: c_int = T::LONG as c_int;
const DOUBLE: c_int = T::DOUBLE as c_int;
const STRING: c_int = T::STRING as c_int;
const BITSTR: c_int = T::BITSTR as c_int;

/// Sentinel the GTI/region builders use for an omitted node argument.
const NO_NODE: NodeId = -99;

/// True for the sorts that `eval.y` calls `expr`. `BOOLEAN` is deliberately
/// excluded: a boolean-valued node is a `bexpr`, a *different* nonterminal, and
/// most arithmetic and comparison productions accept only `expr`. That is why
/// `INTCOL + BOOLCOL` and `1<2==1` are syntax errors.
fn is_expr(t: c_int) -> bool {
    t == LONG || t == DOUBLE
}

/// `expr` or `bexpr` — the sorts the casts and vector literals accept.
fn is_numeric(t: c_int) -> bool {
    t == BOOLEAN || is_expr(t)
}

/// Convert a byte slice to a NUL-terminated `c_char` buffer.
fn cstr(b: &[u8]) -> Vec<c_char> {
    let mut v: Vec<c_char> = b.iter().map(|&c| c as c_char).collect();
    v.push(0);
    v
}

fn unsupported(what: &str, name: &[u8]) -> ParseError {
    ParseError::semantic(format!(
        "Function({what}) not supported: {}(",
        String::from_utf8_lossy(name)
    ))
}

pub(crate) struct Lowerer<'a> {
    pub(crate) p: &'a mut ParseData,
    /// Names already resolved by [`super::resolve`], keyed by byte offset.
    pub(crate) names: &'a Resolutions,
}

impl Lowerer<'_> {
    // -----------------------------------------------------------------
    // Node accessors, mirroring eval.y's TYPE / OPER / SIZE macros
    // -----------------------------------------------------------------

    fn ntype(&self, n: NodeId) -> c_int {
        self.p.Nodes[n as usize].ntype
    }
    fn set_ntype(&mut self, n: NodeId, t: c_int) {
        self.p.Nodes[n as usize].ntype = t;
    }
    fn size(&self, n: NodeId) -> c_long {
        self.p.Nodes[n as usize].value.nelem
    }
    fn set_size(&mut self, n: NodeId, v: c_long) {
        self.p.Nodes[n as usize].value.nelem = v;
    }
    fn is_const(&self, n: NodeId) -> bool {
        self.p.Nodes[n as usize].operation == CONST_OP
    }
    /// The value of a node already known to be a `LONG` constant.
    fn const_lng(&self, n: NodeId) -> c_long {
        unsafe { self.p.Nodes[n as usize].value.data.lng }
    }

    /// `TEST(a)`: a negative index means the builder failed.
    fn test(&self, n: NodeId) -> LRes {
        if n < 0 {
            Err(ParseError::semantic(
                "Couldn't build node structure: out of memory?",
            ))
        } else {
            Ok(n)
        }
    }

    /// `PROMOTE(a,b)`: raise whichever side is lower on BOOLEAN < LONG < DOUBLE.
    fn promote(&mut self, a: &mut NodeId, b: &mut NodeId) {
        let (ta, tb) = (self.ntype(*a), self.ntype(*b));
        if ta > tb {
            *b = New_Unary(self.p, ta, 0, *b);
        } else if ta < tb {
            *a = New_Unary(self.p, tb, 0, *a);
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

    fn func(&mut self, ret: c_int, op: funcOp, args: &[NodeId]) -> NodeId {
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

    fn func_size(&mut self, ret: c_int, op: funcOp, args: &[NodeId], size: c_int) -> NodeId {
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
        New_Const(
            self.p,
            LONG,
            (&raw const v).cast::<c_void>(),
            size_of::<c_long>() as c_long,
        )
    }

    fn const_double(&mut self, v: c_double) -> NodeId {
        New_Const(
            self.p,
            DOUBLE,
            (&raw const v).cast::<c_void>(),
            size_of::<c_double>() as c_long,
        )
    }

    /// `if( TYPE(n) != DOUBLE ) n = New_Unary(DOUBLE, 0, n);`
    fn as_double(&mut self, n: NodeId) -> NodeId {
        if self.ntype(n) != DOUBLE {
            New_Unary(self.p, DOUBLE, 0, n)
        } else {
            n
        }
    }

    /// `if( TYPE(n) != LONG ) n = New_Unary(LONG, 0, n);`
    fn as_long(&mut self, n: NodeId) -> NodeId {
        if self.ntype(n) != LONG {
            New_Unary(self.p, LONG, 0, n)
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
    fn resolve_token(
        &mut self,
        name: &[u8],
        at: usize,
    ) -> Result<(c_int, FITS_PARSER_YYSTYPE), ParseError> {
        match self.names.get(&at) {
            Some(&(dtype, lval)) => Ok((dtype, lval)),
            None => Err(ParseError::syntax("unknown column or keyword", at, name)),
        }
    }

    /// Build the node the grammar would have built for a resolved name.
    fn resolve(&mut self, name: &[u8], at: usize) -> LRes {
        let (dtype, lval) = self.resolve_token(name, at)?;
        unsafe {
            let n = match token_of(dtype) {
                Some(T::COLUMN | T::BCOLUMN | T::SCOLUMN | T::BITCOL) => {
                    New_Column(self.p, lval.lng as c_int)
                }
                Some(T::LONG) => self.const_long(lval.lng),
                Some(T::DOUBLE) => self.const_double(lval.dbl),
                Some(T::BOOLEAN) => New_Const(
                    self.p,
                    BOOLEAN,
                    (&raw const lval.log).cast::<c_void>(),
                    size_of::<c_char>() as c_long,
                ),
                Some(t @ (T::STRING | T::BITSTR)) => {
                    let tag = if matches!(t, T::STRING) { STRING } else { BITSTR };
                    let len = lval.astr.iter().position(|&c| c == 0).unwrap_or(0);
                    let n = New_Const(
                        self.p,
                        tag,
                        lval.astr.as_ptr().cast::<c_void>(),
                        len as c_long + 1,
                    );
                    let n = self.test(n)?;
                    self.set_size(n, len as c_long);
                    n
                }
                /* getData has already set lParse.status and pushed a message */
                _ => return Err(ParseError::syntax("unknown column or keyword", at, name)),
            };
            self.test(n)
        }
    }

    /// Resolve a name that must be a *column*, for the `NAME{offset}` form.
    fn resolve_column(&mut self, name: &[u8], at: usize) -> Result<c_int, ParseError> {
        let (dtype, lval) = self.resolve_token(name, at)?;
        match token_of(dtype) {
            /* spec 3.4 (1): the offset form exists only for column tokens */
            Some(T::COLUMN | T::BCOLUMN | T::SCOLUMN | T::BITCOL) => unsafe { Ok(lval.lng as c_int) },
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
                let b: c_char = c_char::from(*v);
                let n = New_Const(
                    self.p,
                    BOOLEAN,
                    (&raw const b).cast::<c_void>(),
                    size_of::<c_char>() as c_long,
                );
                self.test(n)
            }
            AstKind::Str(s) => self.new_text_const(s, STRING),
            AstKind::BitStr(s) => self.new_text_const(s, BITSTR),
            AstKind::RowRef => {
                let n = self.func(LONG, ROW_FCT, &[]);
                self.test(n)
            }
            AstKind::NullRef => {
                let n = self.func(LONG, NULL_FCT, &[]);
                self.test(n)
            }
            AstKind::SNullRef => {
                let n = self.func(STRING, NULL_FCT, &[]);
                self.test(n)
            }
            AstKind::Ident(name) => self.resolve(name, a.at),
            AstKind::Keyword(name) => self.resolve(name, a.at),
            AstKind::Offset { name, off } => {
                let col = self.resolve_column(name, a.at)?;
                let o = self.lower(off)?;
                if self.ntype(o) != LONG || !self.is_const(o) {
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

    /// `New_Const` for a string or bit-string literal, with `SIZE` set to the
    /// length rather than the length-with-NUL that `New_Const` records.
    fn new_text_const(&mut self, s: &[u8], tag: c_int) -> LRes {
        let buf = cstr(s);
        let n = New_Const(
            self.p,
            tag,
            buf.as_ptr().cast::<c_void>(),
            s.len() as c_long + 1,
        );
        let n = self.test(n)?;
        self.set_size(n, s.len() as c_long);
        Ok(n)
    }

    fn lower_unary(&mut self, op: UnOp, arg: &Ast, at: usize) -> LRes {
        let n = self.lower(arg)?;
        let t = self.ntype(n);
        let r = match op {
            UnOp::Neg => {
                if !is_expr(t) {
                    return Err(ParseError::syntax("unary '-' needs a number", at, b"-"));
                }
                New_Unary(self.p, t, T::UMINUS as c_int, n)
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
                BOOLEAN => New_Unary(self.p, BOOLEAN, T::NOT as c_int, n),
                BITSTR => New_Unary(self.p, BITSTR, T::NOT as c_int, n),
                _ => {
                    return Err(ParseError::syntax(
                        "'!' needs a boolean or bit-string operand",
                        at,
                        b"!",
                    ));
                }
            },
            UnOp::IntCast => {
                if !is_numeric(t) {
                    return Err(ParseError::syntax("(int) needs a number", at, b"(int)"));
                }
                New_Unary(self.p, LONG, T::INTCAST as c_int, n)
            }
            UnOp::FltCast => {
                if !is_numeric(t) {
                    return Err(ParseError::syntax("(float) needs a number", at, b"(float)"));
                }
                New_Unary(self.p, DOUBLE, T::FLTCAST as c_int, n)
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
                return Err(ParseError::semantic("range test needs numeric operands"));
            }
        }
        self.promote(&mut v, &mut l);
        self.promote(&mut v, &mut h);
        self.promote(&mut l, &mut h);
        let lo_le = New_BinOp(self.p, BOOLEAN, l, T::LTE as c_int, v);
        let le_hi = New_BinOp(self.p, BOOLEAN, v, T::LTE as c_int, h);
        let r = New_BinOp(self.p, BOOLEAN, lo_le, T::AND as c_int, le_hi);
        self.test(r)
    }

    fn lower_deref(&mut self, base: &Ast, idx: &[Ast], at: usize) -> LRes {
        let b = self.lower(base)?;
        /* spec 3.4 (2): expr, bexpr and bits subscript; sexpr does not */
        if self.ntype(b) == STRING {
            return Err(ParseError::syntax("strings cannot be subscripted", at, b"["));
        }
        let mut dims = [0 as NodeId; 5];
        for (slot, e) in dims.iter_mut().zip(idx) {
            *slot = self.lower(e)?;
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
                all_bool = te == BOOLEAN;
                vec_node = Some(v);
                continue;
            };

            if te != BOOLEAN {
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
                (STRING, STRING) | (BITSTR, BITSTR) => {
                    let total = self.size(a) + self.size(b);
                    if total >= c_long::from(MAX_STRLEN) {
                        return Err(ParseError::semantic(if ta == STRING {
                            "Combined string size exceeds 255 characters"
                        } else {
                            "Combined bit string size exceeds 255 bits"
                        }));
                    }
                    let n = New_BinOp(self.p, ta, a, i32::from(b'+'), b);
                    let n = self.test(n)?;
                    self.set_size(n, total);
                    return Ok(n);
                }
                _ if is_expr(ta) && is_expr(tb) => self.arith(a, b, b'+')?,
                _ => return Err(self.bad_operands("+", at)),
            },
            BinOp::Sub => self.numeric_only(a, b, b'-', at)?,
            BinOp::Div => self.numeric_only(a, b, b'/', at)?,
            BinOp::Mod => self.numeric_only(a, b, b'%', at)?,
            BinOp::Mul => match (ta, tb) {
                /* `expr '*' bexpr` and `bexpr '*' expr` coerce the boolean */
                (x, BOOLEAN) if is_expr(x) => {
                    b = New_Unary(self.p, ta, 0, b);
                    New_BinOp(self.p, ta, a, i32::from(b'*'), b)
                }
                (BOOLEAN, y) if is_expr(y) => {
                    a = New_Unary(self.p, tb, 0, a);
                    New_BinOp(self.p, tb, a, i32::from(b'*'), b)
                }
                _ => self.numeric_only(a, b, b'*', at)?,
            },
            BinOp::Pow => {
                if !is_expr(ta) || !is_expr(tb) {
                    return Err(self.bad_operands("**", at));
                }
                self.promote(&mut a, &mut b);
                let t = self.ntype(a);
                New_BinOp(self.p, t, a, T::POWER as c_int, b)
            }

            /* `bits & bits` / `bits | bits`, or integer bitwise */
            BinOp::BitAnd | BinOp::BitOr => {
                let ch = if op == BinOp::BitAnd { b'&' } else { b'|' };
                if ta == BITSTR && tb == BITSTR {
                    let n = New_BinOp(self.p, BITSTR, a, i32::from(ch), b);
                    let n = self.test(n)?;
                    let sz = self.size(a).max(self.size(b));
                    self.set_size(n, sz);
                    return Ok(n);
                }
                if ta != LONG || tb != LONG {
                    return Err(ParseError::semantic(
                        "Bitwise operations with incompatible types; only (bit OP bit) and (int OP int) are allowed",
                    ));
                }
                New_BinOp(self.p, ta, a, i32::from(ch), b)
            }
            BinOp::BitXor => {
                if ta != LONG || tb != LONG {
                    return Err(ParseError::semantic(
                        "Bitwise operations with incompatible types; only (bit OP bit) and (int OP int) are allowed",
                    ));
                }
                New_BinOp(self.p, ta, a, i32::from(b'^'), b)
            }

            BinOp::And | BinOp::Or => {
                if ta != BOOLEAN || tb != BOOLEAN {
                    return Err(self.bad_operands(
                        if op == BinOp::And { "&&" } else { "||" },
                        at,
                    ));
                }
                let code = if op == BinOp::And {
                    T::AND as c_int
                } else {
                    T::OR as c_int
                };
                New_BinOp(self.p, BOOLEAN, a, code, b)
            }

            BinOp::Eq | BinOp::Ne | BinOp::Gt | BinOp::Lt | BinOp::Gte | BinOp::Lte => {
                return self.lower_compare(op, a, b, at);
            }
            BinOp::Approx => {
                if !is_expr(ta) || !is_expr(tb) {
                    return Err(self.bad_operands("~", at));
                }
                self.promote(&mut a, &mut b);
                New_BinOp(self.p, BOOLEAN, a, i32::from(b'~'), b)
            }
        };
        self.test(r)
    }

    /// `FUNCTION expr …` productions accept only `expr` operands. Without this
    /// a string or bit-string argument reaches `Do_Func`, which reads the
    /// buffer through a `*long` and segfaults.
    fn require_expr(&self, n: NodeId, name: &[u8]) -> Result<(), ParseError> {
        if is_expr(self.ntype(n)) {
            Ok(())
        } else {
            Err(unsupported("expr", name))
        }
    }

    /// `FUNCTION (expr|bexpr) …`: the first argument of AXISELEM / NAXES /
    /// ARRAY may also be boolean.
    fn require_numeric(&self, n: NodeId, name: &[u8]) -> Result<(), ParseError> {
        if is_numeric(self.ntype(n)) {
            Ok(())
        } else {
            Err(unsupported("expr", name))
        }
    }

    fn bad_operands(&self, op: &str, at: usize) -> ParseError {
        ParseError::syntax(
            format!("operands of '{op}' have incompatible types"),
            at,
            op.as_bytes(),
        )
    }

    /// `expr OP expr` with numeric promotion, the shape of most `eval.y`
    /// arithmetic actions.
    fn arith(&mut self, mut a: NodeId, mut b: NodeId, ch: u8) -> LRes {
        self.promote(&mut a, &mut b);
        let t = self.ntype(a);
        Ok(New_BinOp(self.p, t, a, i32::from(ch), b))
    }

    fn numeric_only(&mut self, a: NodeId, b: NodeId, ch: u8, at: usize) -> LRes {
        if !is_expr(self.ntype(a)) || !is_expr(self.ntype(b)) {
            return Err(self.bad_operands(
                core::str::from_utf8(core::slice::from_ref(&ch)).unwrap_or("?"),
                at,
            ));
        }
        self.arith(a, b, ch)
    }

    fn lower_compare(&mut self, op: BinOp, mut a: NodeId, mut b: NodeId, at: usize) -> LRes {
        let (ta, tb) = (self.ntype(a), self.ntype(b));
        let code = match op {
            BinOp::Eq => T::EQ as c_int,
            BinOp::Ne => T::NE as c_int,
            BinOp::Gt => T::GT as c_int,
            BinOp::Lt => T::LT as c_int,
            BinOp::Gte => T::GTE as c_int,
            BinOp::Lte => T::LTE as c_int,
            _ => unreachable!(),
        };

        /* `bits CMP bits` and `sexpr CMP sexpr` force a scalar result */
        if (ta == BITSTR && tb == BITSTR) || (ta == STRING && tb == STRING) {
            let n = New_BinOp(self.p, BOOLEAN, a, code, b);
            let n = self.test(n)?;
            self.set_size(n, 1);
            return Ok(n);
        }
        /* `bexpr EQ bexpr` and `bexpr NE bexpr` — no ordering on booleans */
        if ta == BOOLEAN && tb == BOOLEAN && matches!(op, BinOp::Eq | BinOp::Ne) {
            let n = New_BinOp(self.p, BOOLEAN, a, code, b);
            return self.test(n);
        }
        if !is_expr(ta) || !is_expr(tb) {
            return Err(self.bad_operands("comparison", at));
        }
        self.promote(&mut a, &mut b);
        let n = New_BinOp(self.p, BOOLEAN, a, code, b);
        self.test(n)
    }

    // -----------------------------------------------------------------
    // Conditional
    // -----------------------------------------------------------------

    fn lower_ternary(&mut self, cond: &Ast, then: &Ast, els: &Ast, at: usize) -> LRes {
        let c = self.lower(cond)?;
        if self.ntype(c) != BOOLEAN {
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
        if tx == BITSTR || ty == BITSTR {
            return Err(ParseError::syntax(
                "'?:' is not defined for bit strings",
                at,
                b"?",
            ));
        }

        if tx == STRING || ty == STRING {
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
            let n = self.func_size(0, IFTHENELSE_FCT, &[x, y, c], out);
            let n = self.test(n)?;
            self.copy_dims_if(n, x, y);
            return Ok(n);
        }

        let both_bool = tx == BOOLEAN && ty == BOOLEAN;
        if !both_bool {
            self.promote(&mut x, &mut y);
        }
        if !self.dims_ok(x, y) {
            return Err(ParseError::semantic(
                "Incompatible dimensions in '?:' arguments",
            ));
        }
        let n = self.func(0, IFTHENELSE_FCT, &[x, y, c]);
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
            self.set_ntype(c, BOOLEAN);
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
            CallKind::GtiFilter => self.lower_gti(GTIFILT_FCT, name, args, at),
            CallKind::GtiFind => self.lower_gti(GTIFIND_FCT, name, args, at),
            CallKind::GtiOverlap => self.lower_gti_overlap(name, args, at),
            CallKind::RegFilter => self.lower_regfilter(name, args, at),
            CallKind::IFunction => {
                /* STRSTR(sexpr, sexpr) */
                let [s, sub] = self.lower_n::<2>(args, name, at)?;
                if self.ntype(s) != STRING || self.ntype(sub) != STRING {
                    return Err(unsupported("string,string", name));
                }
                let n = self.func(LONG, STRPOS_FCT, &[s, sub]);
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
                    STRING => {
                        let n = self.func(BOOLEAN, ISNULL_FCT, &[a]);
                        self.test(n)
                    }
                    /* spec 3.4 (7): ISNULL is not defined for bit strings */
                    BITSTR => Err(unsupported("bits", name)),
                    _ => {
                        let n = self.func(0, ISNULL_FCT, &[a]);
                        let n = self.test(n)?;
                        /* keep the argument's size but return BOOLEAN */
                        self.set_ntype(n, BOOLEAN);
                        Ok(n)
                    }
                }
            }
            (b"NEAR", 3) => self.region_fct(NEAR_FCT, &ns, "NEAR"),
            (b"CIRCLE", 5) => self.region_fct(CIRCLE_FCT, &ns, "CIRCLE"),
            (b"BOX", 7) => self.region_fct(BOX_FCT, &ns, "BOX or ELLIPSE"),
            (b"ELLIPSE", 7) => self.region_fct(ELPS_FCT, &ns, "BOX or ELLIPSE"),
            _ => Err(ParseError::syntax(
                format!(
                    "Boolean Function {}() with {} argument(s) not supported",
                    String::from_utf8_lossy(name),
                    ns.len()
                ),
                at,
                name,
            )),
        }
    }

    /// `NEAR`, `CIRCLE`, `BOX` and `ELLIPSE`: all arguments to DOUBLE, all
    /// dimensions pairwise compatible, then propagate the largest shape.
    fn region_fct(&mut self, op: funcOp, ns: &[NodeId], what: &str) -> LRes {
        for &n in ns {
            self.require_expr(n, what.as_bytes())?;
        }
        let ds: Vec<NodeId> = ns.iter().map(|&n| self.as_double(n)).collect();
        for w in ds.windows(2) {
            if !self.dims_ok(w[0], w[1]) {
                return Err(ParseError::semantic(format!(
                    "Dimensions of {what} arguments are not compatible"
                )));
            }
        }
        let n = self.func(BOOLEAN, op, &ds);
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
                    let n = self.func(DOUBLE, RND_FCT, &[]);
                    self.test(n)
                }
                b"RANDOMN" => {
                    let n = self.func(DOUBLE, GASRND_FCT, &[]);
                    self.test(n)
                }
                _ => Err(unsupported("", name)),
            },
            1 => self.func1(name, ns[0]),
            2 => self.func2(name, ns[0], ns[1]),
            3 => self.func3(name, ns[0], ns[1], ns[2]),
            4 => self.func4(name, &ns),
            _ => Err(ParseError::syntax(
                format!(
                    "{}() with {} arguments is not supported",
                    String::from_utf8_lossy(name),
                    ns.len()
                ),
                at,
                name,
            )),
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
            BOOLEAN => match name {
                b"SUM" => {
                    let n = self.func(LONG, SUM_FCT, &[a]);
                    self.test(n)
                }
                b"ACCUM" => self.accum(a, LONG, T::ACCUM as c_int),
                /* `FUNCTION bexpr ')'` supports only SUM, NELEM and ACCUM.
                Falling through to the numeric list would accept `ABS(BOOLCOL)`
                and then hand Do_Func a char buffer to read through a *long. */
                _ => Err(unsupported("bool", name)),
            },
            STRING => match name {
                b"NVALID" => {
                    let n = self.func(LONG, NONNULL_FCT, &[a]);
                    self.test(n)
                }
                _ => Err(unsupported("str", name)),
            },
            BITSTR => match name {
                /* bit arrays have no NULLs, so NVALID is the element count */
                b"NVALID" => {
                    let n = self.const_long(sz);
                    self.test(n)
                }
                b"SUM" => {
                    let n = self.func(LONG, SUM_FCT, &[a]);
                    self.test(n)
                }
                b"MIN" | b"MAX" => {
                    let op = if name == b"MIN" { MIN1_FCT } else { MAX1_FCT };
                    let n = self.func(t, op, &[a]);
                    let n = self.test(n)?;
                    /* a is a vector, so the result is never constant and it is
                    safe to force the size */
                    self.set_size(n, 1);
                    Ok(n)
                }
                b"ACCUM" => self.accum(a, LONG, T::ACCUM as c_int),
                _ => Err(unsupported("bits", name)),
            },
            _ => self.func1_numeric(name, a),
        }
    }

    fn accum(&mut self, a: NodeId, ret: c_int, code: c_int) -> LRes {
        let zero = if ret == DOUBLE {
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
            b"SUM" => self.func(t, SUM_FCT, &[a]),
            b"AVERAGE" => self.func(DOUBLE, AVERAGE_FCT, &[a]),
            b"STDDEV" => self.func(DOUBLE, STDDEV_FCT, &[a]),
            b"MEDIAN" => self.func(t, MEDIAN_FCT, &[a]),
            b"NVALID" => self.func(LONG, NONNULL_FCT, &[a]),
            b"ACCUM" if t == LONG => return self.accum(a, LONG, T::ACCUM as c_int),
            b"ACCUM" if t == DOUBLE => return self.accum(a, DOUBLE, T::ACCUM as c_int),
            b"SEQDIFF" if t == LONG => return self.accum(a, LONG, T::DIFF as c_int),
            b"SEQDIFF" if t == DOUBLE => return self.accum(a, DOUBLE, T::DIFF as c_int),
            b"ABS" => self.func(0, ABS_FCT, &[a]),
            b"MIN" => self.func(t, MIN1_FCT, &[a]),
            b"MAX" => self.func(t, MAX1_FCT, &[a]),
            b"RANDOM" | b"RANDOMN" => {
                let op = if name == b"RANDOM" { RND_FCT } else { GASRND_FCT };
                let n = self.func(0, op, &[a]);
                let n = self.test(n)?;
                self.set_ntype(n, DOUBLE);
                return Ok(n);
            }
            b"ELEMENTNUM" => {
                if self.is_const(a) {
                    let n = self.const_long(1);
                    return self.test(n);
                }
                let n = self.func(0, ELEMNUM_FCT, &[a]);
                let n = self.test(n)?;
                self.set_ntype(n, LONG);
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
                /* everything else takes a DOUBLE */
                let d = self.as_double(a);
                let op = match name {
                    b"SIN" => SIN_FCT,
                    b"COS" => COS_FCT,
                    b"TAN" => TAN_FCT,
                    b"ARCSIN" | b"ASIN" => ASIN_FCT,
                    b"ARCCOS" | b"ACOS" => ACOS_FCT,
                    b"ARCTAN" | b"ATAN" => ATAN_FCT,
                    b"SINH" => SINH_FCT,
                    b"COSH" => COSH_FCT,
                    b"TANH" => TANH_FCT,
                    b"EXP" => EXP_FCT,
                    b"LOG" => LOG_FCT,
                    b"LOG10" => LOG10_FCT,
                    b"SQRT" => SQRT_FCT,
                    b"ROUND" => ROUND_FCT,
                    b"FLOOR" => FLOOR_FCT,
                    b"CEIL" => CEIL_FCT,
                    b"RANDOMP" => {
                        let n = self.func(0, POIRND_FCT, &[d]);
                        let n = self.test(n)?;
                        self.set_ntype(n, LONG);
                        return Ok(n);
                    }
                    _ => return Err(unsupported("expr", name)),
                };
                self.func(0, op, &[d])
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
                if ta == BITSTR || tb == BITSTR {
                    return Err(unsupported("bits,bits", name));
                }
                if ta == STRING || tb == STRING {
                    if ta != tb {
                        return Err(unsupported("string,expr", name));
                    }
                    let out = self.size(a).max(self.size(b)) as c_int;
                    let n = self.func_size(0, DEFNULL_FCT, &[a, b], out);
                    let n = self.test(n)?;
                    if self.size(b) > self.size(a) {
                        let sb = self.size(b);
                        self.set_size(n, sb);
                    }
                    return Ok(n);
                }
                self.require_numeric(a, name)?;
                self.require_numeric(b, name)?;
                /* `FUNCTION expr ',' expr` and `FUNCTION bexpr ',' bexpr` are
                separate productions; there is no mixed form */
                if (ta == BOOLEAN) != (tb == BOOLEAN) {
                    return Err(unsupported("expr,expr", name));
                }
                if !(self.size(a) >= self.size(b) && self.dims_ok(a, b)) {
                    return Err(ParseError::semantic(
                        "Dimensions of DEFNULL arguments are not compatible",
                    ));
                }
                let (mut a, mut b) = (a, b);
                if !(ta == BOOLEAN && tb == BOOLEAN) {
                    self.promote(&mut a, &mut b);
                }
                let n = self.func(0, DEFNULL_FCT, &[a, b]);
                self.test(n)
            }
            b"ARCTAN2" => {
                self.require_expr(a, name)?;
                self.require_expr(b, name)?;
                let a = self.as_double(a);
                let b = self.as_double(b);
                if !self.dims_ok(a, b) {
                    return Err(ParseError::semantic(
                        "Dimensions of arctan2 arguments are not compatible",
                    ));
                }
                let n = self.func(0, ATAN2_FCT, &[a, b]);
                let n = self.test(n)?;
                self.copy_dims_if(n, a, b);
                Ok(n)
            }
            b"MIN" | b"MAX" => {
                self.require_expr(a, name)?;
                self.require_expr(b, name)?;
                let (mut a, mut b) = (a, b);
                self.promote(&mut a, &mut b);
                if !self.dims_ok(a, b) {
                    return Err(ParseError::semantic(format!(
                        "Dimensions of {}(a,b) arguments are not compatible",
                        String::from_utf8_lossy(name).to_lowercase()
                    )));
                }
                let op = if name == b"MIN" { MIN2_FCT } else { MAX2_FCT };
                let n = self.func(0, op, &[a, b]);
                let n = self.test(n)?;
                self.copy_dims_if(n, a, b);
                Ok(n)
            }
            b"SETNULL" => {
                self.require_expr(a, name)?;
                self.require_expr(b, name)?;
                if !self.is_const(a) || self.size(a) != 1 {
                    return Err(ParseError::semantic(
                        "SETNULL first argument must be a scalar constant",
                    ));
                }
                /* make the sentinel the same sort as the value */
                let a = if ta != tb {
                    New_Unary(self.p, tb, 0, a)
                } else {
                    a
                };
                let n = self.func(0, SETNULL_FCT, &[b, a]);
                self.test(n)
            }
            b"AXISELEM" => {
                self.require_numeric(a, name)?;
                self.require_expr(b, name)?;
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
                let n = self.func(0, AXISELEM_FCT, &[a, b]);
                let n = self.test(n)?;
                self.set_ntype(n, LONG);
                Ok(n)
            }
            b"NAXES" => {
                self.require_numeric(a, name)?;
                self.require_expr(b, name)?;
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
                self.require_numeric(a, name)?;
                self.require_expr(b, name)?;
                let n = New_Array(self.p, a, b);
                self.test(n)
            }
            _ => Err(unsupported("expr,expr", name)),
        }
    }

    fn func3(&mut self, name: &[u8], s: NodeId, pos: NodeId, len: NodeId) -> LRes {
        if name != b"STRMID" {
            return Err(unsupported("expr,expr,expr", name));
        }
        if self.ntype(s) != STRING {
            return Err(unsupported("string,expr,expr", name));
        }
        if self.ntype(pos) != LONG
            || self.size(pos) != 1
            || self.ntype(len) != LONG
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
        let n = self.func_size(0, STRMID_FCT, &[s, pos, len], n_chars as c_int);
        self.test(n)
    }

    fn func4(&mut self, name: &[u8], ns: &[NodeId]) -> LRes {
        if name != b"ANGSEP" {
            return Err(unsupported("expr,expr,expr,expr", name));
        }
        for &n in ns {
            self.require_expr(n, name)?;
        }
        let ds: Vec<NodeId> = ns.iter().map(|&n| self.as_double(n)).collect();
        for w in ds.windows(2) {
            if !self.dims_ok(w[0], w[1]) {
                return Err(ParseError::semantic(
                    "Dimensions of ANGSEP arguments are not compatible",
                ));
            }
        }
        let n = self.func(0, ANGSEP_FCT, &ds);
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

    fn lower_gti(&mut self, op: funcOp, name: &[u8], args: &[Ast], at: usize) -> LRes {
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
                self.require_expr(node1, name)?;
            }
            4 => {
                fname = self.literal_str(&args[0], name)?;
                node1 = self.lower(&args[1])?;
                self.require_expr(node1, name)?;
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
        self.test(n)
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
        self.require_expr(n1, name)?;
        self.require_expr(n2, name)?;
        if args.len() == 5 {
            start = self.literal_str(&args[3], name)?;
            stop = self.literal_str(&args[4], name)?;
        }
        let n = New_GTI(
            self.p,
            GTIOVER_FCT,
            fname.as_mut_ptr(),
            n1,
            n2,
            start.as_mut_ptr(),
            stop.as_mut_ptr(),
        );
        self.test(n)
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
            self.require_expr(nx, name)?;
            self.require_expr(ny, name)?;
        }
        if args.len() == 4 {
            cols = self.literal_str(&args[3], name)?;
        }
        let n = New_REG(
            self.p,
            fname.as_mut_ptr(),
            nx,
            ny,
            cols.as_mut_ptr(),
        );
        self.test(n)
    }
}

