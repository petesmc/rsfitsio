//! Reading what a `GTIFILTER` or `REGFILTER` call refers to.
//!
//! These are the two places where parsing an expression touches the filesystem:
//! a good-time-interval file and an SAO region file, both named by a literal
//! string inside the expression. `eval.y` did the reading inside `New_GTI` and
//! `New_REG`, as the arena was built, and `parser::lower` still gets its
//! [`GtiLoads`] and [`RegionLoads`] by calling those and reading the result
//! back off the node.
//!
//! That is the last thing tying the columnar lowering to the arena. It walks
//! the tree with a shared borrow and the file work needs `&mut ParseData`, so
//! the read has to happen before either lowering rather than during one. Doing
//! it here, as a pass over the `Ast`, gives both the same table -- keyed by the
//! call's byte offset, which is what `Columns` already looks them up by.
//!
//! The file read is not all `New_GTI` and `New_REG` do, and the rest has to
//! come along or it is silently lost: resolving `TIME`, `X` and `Y` when the
//! call leaves them out, and splitting `REGFILTER`'s column-name pair. See
//! PARSER_MIGRATION.md §10.7.

use super::ast::{Ast, AstKind};
use super::error::ParseError;
use super::resolve::Resolutions;
use super::token::CallKind;
use crate::c_types::{c_char, c_double, c_int};
use crate::eval_defs::{ColumnSort, FuncOp, ParseData, ParserValue};
use crate::eval_y::{load_gti, load_region};
use crate::fitscore::ffgcno_safe;
use crate::fitscore::ffpmsg_str;
use crate::region::SAORegion;
use bytemuck::cast_slice;
use std::collections::HashMap;

/// The good-time intervals a `GTIFILTER`/`GTIFIND` call loaded, and whether
/// they came back in ascending order.
///
/// [`load_all`] reads these once, before lowering, and the evaluator uses
/// that load rather than repeating four hundred lines of HDU navigation -- in
/// the same way [`Resolutions`] hands it the already-resolved column names.
#[derive(Clone, Debug, Default, PartialEq)]
pub(crate) struct GtiData {
    pub(crate) start: Vec<c_double>,
    pub(crate) stop: Vec<c_double>,
    pub(crate) ordered: bool,
    /// The column the times come from when the call did not name an
    /// expression -- `gtifilter()` means the TIME column, which is resolved
    /// here and which the source text does not mention.
    pub(crate) time_col: Option<c_int>,
}

/// The GTI loads of one expression, keyed by the call's byte offset.
pub(crate) type GtiLoads = HashMap<usize, GtiData>;

/// The regions of one expression, keyed by the call's byte offset.
pub(crate) type RegionLoads = HashMap<usize, RegionData>;

/// The region a `REGFILTER` call parsed, and the coordinate columns it used.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct RegionData {
    pub(crate) region: SAORegion,
    /// The X and Y columns resolved when the call did not name expressions,
    /// as `regfilter("f.reg")` does not.
    pub(crate) x_col: Option<c_int>,
    pub(crate) y_col: Option<c_int>,
}

/// A NUL-terminated buffer, as the loaders take.
fn cstr(b: &[u8]) -> Vec<c_char> {
    let mut v: Vec<c_char> = b.iter().map(|&c| c as c_char).collect();
    v.push(0);
    v
}

/// The literal string an argument has to be.
///
/// `GTIFILTER STRING ')'` and its neighbours take a `STRING` token, not an
/// expression, so anything else there matches no production.
fn literal(a: &Ast) -> Result<Vec<c_char>, ParseError> {
    match &a.kind {
        AstKind::Str(s) => Ok(cstr(s)),
        _ => Err(ParseError::syntax(
            "a literal string is required here",
            a.at,
            b"",
        )),
    }
}

/// The column an argument names, when it plainly names one.
///
/// `record_gti` and `record_region` take this off the node the argument
/// lowered to, which reports a column only when the node *is* a column. An
/// expression over one is not a column, and the evaluator has to evaluate it.
fn column_of(a: &Ast, names: &Resolutions) -> Option<c_int> {
    match &a.kind {
        AstKind::Ident(_) | AstKind::Keyword(_) => match names.get(&a.at) {
            Some(ParserValue::Column { index, .. }) => Some(*index),
            _ => None,
        },
        _ => None,
    }
}

/// The numeric column a bare name refers to, the way the C looks one up while
/// building a node for it.
macro_rules! named {
    ($p:expr, $n:literal) => {
        match crate::eval_y::fits_parser_yyGetVariable($p, crate::cs!($n)) {
            Some(ParserValue::Column {
                index,
                sort: ColumnSort::Numeric,
            }) => Some(index),
            _ => None,
        }
    };
}

/// Read every GTI and region file the expression names.
///
/// Errors are the ones `New_GTI` and `New_REG` raise before they build
/// anything: a missing default column, a `GTIOVERLAP` without both ends, an
/// unusable column-name pair. A file that will not open is *not* one of them --
/// the loaders push their own message and set the status, and return `None`,
/// which is reported the same way here as `TEST` reported it there.
pub(crate) fn load_all(
    p: &mut ParseData,
    ast: &Ast,
    names: &Resolutions,
    gti: &mut GtiLoads,
    regions: &mut RegionLoads,
) -> Result<(), ParseError> {
    if let AstKind::Call { kind, name, args } = &ast.kind {
        match kind {
            CallKind::GtiFilter | CallKind::GtiFind | CallKind::GtiOverlap => {
                load_one_gti(p, ast.at, *kind, name, args, names, gti)?;
            }
            CallKind::RegFilter => {
                load_one_region(p, ast.at, args, names, regions)?;
            }
            _ => {}
        }
    }
    for child in children(ast) {
        load_all(p, child, names, gti, regions)?;
    }
    Ok(())
}

/// Every sub-expression of `ast`, so the walk reaches a nested call.
fn children(ast: &Ast) -> Vec<&Ast> {
    match &ast.kind {
        AstKind::Unary { arg, .. } => vec![arg],
        AstKind::Binary { lhs, rhs, .. } => vec![lhs, rhs],
        AstKind::Ternary { cond, then, els } => vec![cond, then, els],
        AstKind::Range { val, lo, hi } => vec![val, lo, hi],
        AstKind::Call { args, .. } | AstKind::Vector(args) => args.iter().collect(),
        AstKind::Deref { base, idx } => {
            let mut v = vec![&**base];
            v.extend(idx.iter());
            v
        }
        AstKind::Offset { off, .. } => vec![off],
        _ => Vec::new(),
    }
}

fn load_one_gti(
    p: &mut ParseData,
    at: usize,
    kind: CallKind,
    name: &[u8],
    args: &[Ast],
    names: &Resolutions,
    out: &mut GtiLoads,
) -> Result<(), ParseError> {
    let op = match kind {
        CallKind::GtiFind => FuncOp::GtiFind,
        CallKind::GtiOverlap => FuncOp::GtiOver,
        _ => FuncOp::GtiFilt,
    };
    let mut fname = cstr(b"");
    let mut start = cstr(b"*START*");
    let mut stop = cstr(b"*STOP*");
    /* which argument, if any, gave the times: GTIOVERLAP takes a pair and the
    other two take one, and the `*START*`/`*STOP*` names come last */
    let time_arg: Option<&Ast> = if op == FuncOp::GtiOver {
        match args.len() {
            3 | 5 => {
                fname = literal(&args[0])?;
                if args.len() == 5 {
                    start = literal(&args[3])?;
                    stop = literal(&args[4])?;
                }
                Some(&args[1])
            }
            _ => {
                /* `startExpr and stopExpr values must be defined for GTIOVERLAP`
                is New_GTI's, but it only reaches it for a shape the grammar
                admits; anything else never reduced */
                return Err(ParseError::syntax("GTIOVERLAP arity", at, name));
            }
        }
    } else {
        match args.len() {
            0 => None,
            1 => {
                fname = literal(&args[0])?;
                None
            }
            2 | 4 => {
                fname = literal(&args[0])?;
                if args.len() == 4 {
                    start = literal(&args[2])?;
                    stop = literal(&args[3])?;
                }
                Some(&args[1])
            }
            _ => return Err(ParseError::syntax("GTIFILTER arity", at, name)),
        }
    };

    /* with no time expression the C resolves TIME itself, and says so when
    there is not one */
    let time_col = match time_arg {
        Some(a) => column_of(a, names),
        None => match named!(p, c"TIME") {
            Some(i) => Some(i),
            None => {
                ffpmsg_str("Could not build TIME column for GTIFILTER/GTIFIND");
                return Err(ParseError::aborted());
            }
        },
    };

    let Some((nrows, loaded)) = load_gti(
        p,
        op,
        fname.as_mut_ptr(),
        start.as_mut_ptr(),
        stop.as_mut_ptr(),
    ) else {
        /* the loader has already said why -- a file that would not open, a
        missing column -- so this adds nothing, as TEST did not */
        return Err(ParseError::aborted());
    };

    let count = nrows.max(0) as usize;
    let data = match loaded {
        Some((times, ordered)) => GtiData {
            start: times[..count].to_vec(),
            stop: times[count..count * 2].to_vec(),
            ordered,
            time_col,
        },
        None => GtiData {
            time_col,
            ..Default::default()
        },
    };
    out.insert(at, data);
    Ok(())
}

fn load_one_region(
    p: &mut ParseData,
    at: usize,
    args: &[Ast],
    names: &Resolutions,
    out: &mut RegionLoads,
) -> Result<(), ParseError> {
    if args.is_empty() || args.len() == 2 || args.len() > 4 {
        return Err(ParseError::syntax("REGFILTER arity", at, b""));
    }
    let mut fname = literal(&args[0])?;

    /* the coordinates, either written out or the X and Y the C looks up */
    let (x_col, y_col) = if args.len() >= 3 {
        (column_of(&args[1], names), column_of(&args[2], names))
    } else {
        let x = match named!(p, c"X") {
            Some(i) => i,
            None => {
                ffpmsg_str("Could not build X column for REGFILTER");
                return Err(ParseError::aborted());
            }
        };
        let y = match named!(p, c"Y") {
            Some(i) => i,
            None => {
                ffpmsg_str("Could not build Y column for REGFILTER");
                return Err(ParseError::aborted());
            }
        };
        (Some(x), Some(y))
    };

    /* A fourth argument names the pair of *table* columns the region's
    coordinates are in, which is a different thing from the expressions above:
    `load_region` wants their numbers, for the WCS. Zero means "not given". */
    let (mut xnum, mut ynum) = (0 as c_int, 0 as c_int);
    if args.len() == 4 {
        let spec = literal(&args[3])?;
        let text: Vec<u8> = spec
            .iter()
            .take_while(|&&c| c != 0)
            .map(|&c| c as u8)
            .collect();
        let mut parts = text
            .split(|&c| c == b' ' || c == b',')
            .filter(|s| !s.is_empty());
        let (Some(cx), Some(cy)) = (parts.next(), parts.next()) else {
            ffpmsg_str("Could not extract valid pair of column names from REGFILTER");
            return Err(ParseError::aborted());
        };
        /* `def_fptr` is the raw handle the parser was given; the C reaches
        through it the same way here */
        let fptr = unsafe { p.def_fptr.as_mut() }.unwrap();
        let mut status = 0;
        ffgcno_safe(fptr, 0, &cstr(cx), &mut xnum, &mut status);
        ffgcno_safe(fptr, 0, &cstr(cy), &mut ynum, &mut status);
        p.status = status;
    }

    let Some(region) = load_region(p, fname.as_mut_ptr(), xnum, ynum) else {
        return Err(ParseError::aborted());
    };
    out.insert(
        at,
        RegionData {
            region: *region,
            x_col,
            y_col,
        },
    );
    Ok(())
}
