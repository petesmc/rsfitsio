//! Per-function kernels for `Do_Func`.
//!
//! `Do_Func` used to be one 3,300-line function holding two `match`es over
//! fifty unrelated operations -- a constant-folding branch and a row-loop
//! branch. Each arm is now a named function here, which makes the operations
//! greppable, gives each one somewhere to hang a comment, and is the shape a
//! kernel-per-operation evaluator would want.
//!
//! Every kernel takes the node it is writing and the arguments `Do_Func`
//! gathered; the `_const` forms compute a single folded value, the `_rows`
//! forms fill the node's row buffer.

use super::*;

/// A kernel: it fills in `Nodes[this]` from the gathered arguments.
pub(crate) type Kernel = fn(&mut ParseData, usize, &mut FuncArgs);

/// The arguments of a function node, gathered once by [`Do_Func`].
#[derive(Clone, Copy)]
pub(crate) struct FuncArgs {
    /// Node index of each argument.
    pub params: [usize; 10],
    /// Element count of each argument, or 0 when it is a constant.
    pub vector: [c_int; 10],
    /// Value of each constant argument.
    pub vals: [lval; 10],
    /// Null flag of each constant argument.
    pub nulls: [c_char; 10],
}

impl Default for FuncArgs {
    fn default() -> Self {
        FuncArgs {
            params: [0; 10],
            vector: [0; 10],
            vals: [lval {
                nelem: 0,
                naxis: 0,
                naxes: [0; 5],
                undef: core::ptr::null_mut::<c_char>(),
                data: NodeValue::Double(0.0),
            }; 10],
            nulls: [0; 10],
        }
    }
}

/// `FuncOp::Sum`
pub(crate) fn sum_const(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    unsafe {
        if (lParse.Nodes[a.params[0]]).ntype == ValueSort::Boolean {
            (lParse.Nodes[this_node_idx]).value.data =
                NodeValue::Long(c_long::from(if c_int::from(a.vals[0].data.log()) != 0 {
                    1
                } else {
                    0
                }));
        } else if (lParse.Nodes[a.params[0]]).ntype == ValueSort::Long {
            (lParse.Nodes[this_node_idx]).value.data = NodeValue::Long(a.vals[0].data.lng());
        } else if (lParse.Nodes[a.params[0]]).ntype == ValueSort::Double {
            (lParse.Nodes[this_node_idx]).value.data = NodeValue::Double(a.vals[0].data.dbl());
        } else if (lParse.Nodes[a.params[0]]).ntype == ValueSort::Bits {
            strcpy(
                (lParse.Nodes[this_node_idx]).value.data.text_mut_ptr(),
                a.vals[0].data.text_mut_ptr(),
            );
        }
    }
}

/// `FuncOp::Average`
pub(crate) fn average_const(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    if (lParse.Nodes[a.params[0]]).ntype == ValueSort::Long {
        (lParse.Nodes[this_node_idx]).value.data =
            NodeValue::Double(a.vals[0].data.lng() as c_double);
    } else if (lParse.Nodes[a.params[0]]).ntype == ValueSort::Double {
        (lParse.Nodes[this_node_idx]).value.data = NodeValue::Double(a.vals[0].data.dbl());
    }
}

/// `FuncOp::Stddev`
pub(crate) fn stddev_const(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    (lParse.Nodes[this_node_idx]).value.data = NodeValue::Double(0.0);
}

/// `FuncOp::Median`
pub(crate) fn median_const(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    if (lParse.Nodes[a.params[0]]).ntype == ValueSort::Boolean {
        (lParse.Nodes[this_node_idx]).value.data =
            NodeValue::Long(c_long::from(if c_int::from(a.vals[0].data.log()) != 0 {
                1
            } else {
                0
            }));
    } else if (lParse.Nodes[a.params[0]]).ntype == ValueSort::Long {
        (lParse.Nodes[this_node_idx]).value.data = NodeValue::Long(a.vals[0].data.lng());
    } else {
        (lParse.Nodes[this_node_idx]).value.data = NodeValue::Double(a.vals[0].data.dbl());
    }
}

/// `FuncOp::PoiRnd`
pub(crate) fn poi_rnd_const(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    if (lParse.Nodes[a.params[0]]).ntype == ValueSort::Double {
        (lParse.Nodes[this_node_idx]).value.data =
            NodeValue::Long(c_long::from(simplerng_getpoisson(a.vals[0].data.dbl())));
    } else {
        (lParse.Nodes[this_node_idx]).value.data = NodeValue::Long(c_long::from(
            simplerng_getpoisson(a.vals[0].data.lng() as c_double),
        ));
    }
}

/// `FuncOp::Abs`
pub(crate) fn abs_const(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    let mut ival: c_long = 0;
    let mut dval: c_double = 0.0;

    if (lParse.Nodes[a.params[0]]).ntype == ValueSort::Double {
        dval = a.vals[0].data.dbl();
        (lParse.Nodes[this_node_idx]).value.data =
            NodeValue::Double(if dval > 0.0 { dval } else { -dval });
    } else {
        ival = a.vals[0].data.lng();
        (lParse.Nodes[this_node_idx]).value.data =
            NodeValue::Long(if ival > 0 { ival } else { -ival });
    }
}

/// `FuncOp::NonNull`
pub(crate) fn non_null_const(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    (lParse.Nodes[this_node_idx]).value.data = NodeValue::Long(1);
}

/// `FuncOp::IsNull`
pub(crate) fn is_null_const(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    (lParse.Nodes[this_node_idx]).value.data = NodeValue::Logical(0);
}

/// `FuncOp::DefNull`
pub(crate) fn def_null_const(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    unsafe {
        if (lParse.Nodes[this_node_idx]).ntype == ValueSort::Boolean {
            (lParse.Nodes[this_node_idx]).value.data = NodeValue::Logical(a.vals[0].data.log());
        } else if (lParse.Nodes[this_node_idx]).ntype == ValueSort::Long {
            (lParse.Nodes[this_node_idx]).value.data = NodeValue::Long(a.vals[0].data.lng());
        } else if (lParse.Nodes[this_node_idx]).ntype == ValueSort::Double {
            (lParse.Nodes[this_node_idx]).value.data = NodeValue::Double(a.vals[0].data.dbl());
        } else if (lParse.Nodes[this_node_idx]).ntype == ValueSort::String {
            strcpy(
                (lParse.Nodes[this_node_idx]).value.data.text_mut_ptr(),
                a.vals[0].data.text_mut_ptr(),
            );
        }
    }
}

/// `FuncOp::SetNull`
pub(crate) fn set_null_const(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    if (lParse.Nodes[this_node_idx]).ntype == ValueSort::Long {
        (lParse.Nodes[this_node_idx]).value.data = NodeValue::Long(a.vals[0].data.lng());
    } else if (lParse.Nodes[this_node_idx]).ntype == ValueSort::Double {
        (lParse.Nodes[this_node_idx]).value.data = NodeValue::Double(a.vals[0].data.dbl());
    }
}

/// `FuncOp::Sin`
pub(crate) fn sin_const(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    (lParse.Nodes[this_node_idx]).value.data = NodeValue::Double(sin(a.vals[0].data.dbl()));
}

/// `FuncOp::Cos`
pub(crate) fn cos_const(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    (lParse.Nodes[this_node_idx]).value.data = NodeValue::Double(cos(a.vals[0].data.dbl()));
}

/// `FuncOp::Tan`
pub(crate) fn tan_const(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    (lParse.Nodes[this_node_idx]).value.data = NodeValue::Double(tan(a.vals[0].data.dbl()));
}

/// `FuncOp::Asin`
pub(crate) fn asin_const(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    let mut dval: c_double = 0.0;

    dval = a.vals[0].data.dbl();
    if dval < -1.0 || dval > 1.0 {
        fits_parser_yyerror(lParse, cs!(c"Out of range argument to arcsin"));
    } else {
        (lParse.Nodes[this_node_idx]).value.data = NodeValue::Double(asin(dval));
    }
}

/// `FuncOp::Acos`
pub(crate) fn acos_const(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    let mut dval: c_double = 0.0;

    dval = a.vals[0].data.dbl();
    if dval < -1.0 || dval > 1.0 {
        fits_parser_yyerror(lParse, cs!(c"Out of range argument to arccos"));
    } else {
        (lParse.Nodes[this_node_idx]).value.data = NodeValue::Double(acos(dval));
    }
}

/// `FuncOp::Atan`
pub(crate) fn atan_const(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    (lParse.Nodes[this_node_idx]).value.data = NodeValue::Double(atan(a.vals[0].data.dbl()));
}

/// `FuncOp::Sinh`
pub(crate) fn sinh_const(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    (lParse.Nodes[this_node_idx]).value.data = NodeValue::Double(sinh(a.vals[0].data.dbl()));
}

/// `FuncOp::Cosh`
pub(crate) fn cosh_const(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    (lParse.Nodes[this_node_idx]).value.data = NodeValue::Double(cosh(a.vals[0].data.dbl()));
}

/// `FuncOp::Tanh`
pub(crate) fn tanh_const(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    (lParse.Nodes[this_node_idx]).value.data = NodeValue::Double(tanh(a.vals[0].data.dbl()));
}

/// `FuncOp::Exp`
pub(crate) fn exp_const(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    (lParse.Nodes[this_node_idx]).value.data = NodeValue::Double(exp(a.vals[0].data.dbl()));
}

/// `FuncOp::Log`
pub(crate) fn log_const(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    let mut dval: c_double = 0.0;

    dval = a.vals[0].data.dbl();
    if dval <= 0.0 {
        fits_parser_yyerror(lParse, cs!(c"Out of range argument to log"));
    } else {
        (lParse.Nodes[this_node_idx]).value.data = NodeValue::Double(log(dval));
    }
}

/// `FuncOp::Log10`
pub(crate) fn log10_const(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    let mut dval: c_double = 0.0;

    dval = a.vals[0].data.dbl();
    if dval <= 0.0 {
        fits_parser_yyerror(lParse, cs!(c"Out of range argument to log10"));
    } else {
        (lParse.Nodes[this_node_idx]).value.data = NodeValue::Double(log10(dval));
    }
}

/// `FuncOp::Sqrt`
pub(crate) fn sqrt_const(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    let mut dval: c_double = 0.0;

    dval = a.vals[0].data.dbl();
    if dval < 0.0 {
        fits_parser_yyerror(lParse, cs!(c"Out of range argument to sqrt"));
    } else {
        (lParse.Nodes[this_node_idx]).value.data = NodeValue::Double(sqrt(dval));
    }
}

/// `FuncOp::Ceil`
pub(crate) fn ceil_const(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    (lParse.Nodes[this_node_idx]).value.data = NodeValue::Double(ceil(a.vals[0].data.dbl()));
}

/// `FuncOp::Floor`
pub(crate) fn floor_const(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    (lParse.Nodes[this_node_idx]).value.data = NodeValue::Double(floor(a.vals[0].data.dbl()));
}

/// `FuncOp::Round`
pub(crate) fn round_const(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    (lParse.Nodes[this_node_idx]).value.data = NodeValue::Double(floor(a.vals[0].data.dbl() + 0.5));
}

/// `FuncOp::Atan2`
pub(crate) fn atan2_const(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    (lParse.Nodes[this_node_idx]).value.data =
        NodeValue::Double(atan2(a.vals[0].data.dbl(), a.vals[1].data.dbl()));
}

/// `FuncOp::AngSep`
pub(crate) fn ang_sep_const(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    (lParse.Nodes[this_node_idx]).value.data = NodeValue::Double(angsep_calc(
        a.vals[0].data.dbl(),
        a.vals[1].data.dbl(),
        a.vals[2].data.dbl(),
        a.vals[3].data.dbl(),
    ));
    /* DEVIATION from CFITSIO 4.7.0: "case angsep_fct:" there
    is missing its "break;" and falls through into
    "case min1_fct:", which overwrites the separation just
    computed with a.vals[0] - so a constant-folded ANGSEP
    returns its first argument.  Fix submitted upstream; the
    non-constant path in Do_Func has always been correct. */
}

/// `FuncOp::Min1`
pub(crate) fn min1_const(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    unsafe {
        if (lParse.Nodes[this_node_idx]).ntype == ValueSort::Double {
            (lParse.Nodes[this_node_idx]).value.data = NodeValue::Double(a.vals[0].data.dbl());
        } else if (lParse.Nodes[this_node_idx]).ntype == ValueSort::Long {
            (lParse.Nodes[this_node_idx]).value.data = NodeValue::Long(a.vals[0].data.lng());
        } else if (lParse.Nodes[this_node_idx]).ntype == ValueSort::Bits {
            strcpy(
                (lParse.Nodes[this_node_idx]).value.data.text_mut_ptr(),
                a.vals[0].data.text_mut_ptr(),
            );
        }
    }
}

/// `FuncOp::Min2`
pub(crate) fn min2_const(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    if (lParse.Nodes[this_node_idx]).ntype == ValueSort::Double {
        (lParse.Nodes[this_node_idx]).value.data =
            NodeValue::Double(if a.vals[0].data.dbl() < a.vals[1].data.dbl() {
                a.vals[0].data.dbl()
            } else {
                a.vals[1].data.dbl()
            });
    } else if (lParse.Nodes[this_node_idx]).ntype == ValueSort::Long {
        (lParse.Nodes[this_node_idx]).value.data =
            NodeValue::Long(if a.vals[0].data.lng() < a.vals[1].data.lng() {
                a.vals[0].data.lng()
            } else {
                a.vals[1].data.lng()
            });
    }
}

/// `FuncOp::Max1`
pub(crate) fn max1_const(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    unsafe {
        if (lParse.Nodes[this_node_idx]).ntype == ValueSort::Double {
            (lParse.Nodes[this_node_idx]).value.data = NodeValue::Double(a.vals[0].data.dbl());
        } else if (lParse.Nodes[this_node_idx]).ntype == ValueSort::Long {
            (lParse.Nodes[this_node_idx]).value.data = NodeValue::Long(a.vals[0].data.lng());
        } else if (lParse.Nodes[this_node_idx]).ntype == ValueSort::Bits {
            strcpy(
                (lParse.Nodes[this_node_idx]).value.data.text_mut_ptr(),
                a.vals[0].data.text_mut_ptr(),
            );
        }
    }
}

/// `FuncOp::Max2`
pub(crate) fn max2_const(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    if (lParse.Nodes[this_node_idx]).ntype == ValueSort::Double {
        (lParse.Nodes[this_node_idx]).value.data =
            NodeValue::Double(if a.vals[0].data.dbl() > a.vals[1].data.dbl() {
                a.vals[0].data.dbl()
            } else {
                a.vals[1].data.dbl()
            });
    } else if (lParse.Nodes[this_node_idx]).ntype == ValueSort::Long {
        (lParse.Nodes[this_node_idx]).value.data =
            NodeValue::Long(if a.vals[0].data.lng() > a.vals[1].data.lng() {
                a.vals[0].data.lng()
            } else {
                a.vals[1].data.lng()
            });
    }
}

/// `FuncOp::Near`
pub(crate) fn near_const(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    (lParse.Nodes[this_node_idx]).value.data = NodeValue::Logical(bnear(
        a.vals[0].data.dbl(),
        a.vals[1].data.dbl(),
        a.vals[2].data.dbl(),
    ));
}

/// `FuncOp::Circle`
pub(crate) fn circle_const(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    (lParse.Nodes[this_node_idx]).value.data = NodeValue::Logical(circle(
        a.vals[0].data.dbl(),
        a.vals[1].data.dbl(),
        a.vals[2].data.dbl(),
        a.vals[3].data.dbl(),
        a.vals[4].data.dbl(),
    ));
}

/// `FuncOp::Box`
pub(crate) fn box_const(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    (lParse.Nodes[this_node_idx]).value.data = NodeValue::Logical(saobox(
        a.vals[0].data.dbl(),
        a.vals[1].data.dbl(),
        a.vals[2].data.dbl(),
        a.vals[3].data.dbl(),
        a.vals[4].data.dbl(),
        a.vals[5].data.dbl(),
        a.vals[6].data.dbl(),
    ));
}

/// `FuncOp::Ellipse`
pub(crate) fn ellipse_const(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    (lParse.Nodes[this_node_idx]).value.data = NodeValue::Logical(ellipse(
        a.vals[0].data.dbl(),
        a.vals[1].data.dbl(),
        a.vals[2].data.dbl(),
        a.vals[3].data.dbl(),
        a.vals[4].data.dbl(),
        a.vals[5].data.dbl(),
        a.vals[6].data.dbl(),
    ));
}

/// `FuncOp::IfThenElse`
pub(crate) fn if_then_else_const(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    unsafe {
        match (lParse.Nodes[this_node_idx]).ntype {
            ValueSort::Boolean => {
                (lParse.Nodes[this_node_idx]).value.data = NodeValue::Logical(
                    (if c_int::from(a.vals[2].data.log()) != 0 {
                        c_int::from(a.vals[0].data.log())
                    } else {
                        c_int::from(a.vals[1].data.log())
                    }) as c_char,
                );
            }
            ValueSort::Long => {
                (lParse.Nodes[this_node_idx]).value.data =
                    NodeValue::Long(if c_int::from(a.vals[2].data.log()) != 0 {
                        a.vals[0].data.lng()
                    } else {
                        a.vals[1].data.lng()
                    });
            }
            ValueSort::Double => {
                (lParse.Nodes[this_node_idx]).value.data =
                    NodeValue::Double(if c_int::from(a.vals[2].data.log()) != 0 {
                        a.vals[0].data.dbl()
                    } else {
                        a.vals[1].data.dbl()
                    });
            }
            ValueSort::String => {
                strcpy(
                    (lParse.Nodes[this_node_idx]).value.data.text_mut_ptr(),
                    if c_int::from(a.vals[2].data.log()) != 0 {
                        a.vals[0].data.text_mut_ptr()
                    } else {
                        a.vals[1].data.text_mut_ptr()
                    },
                );
            }
            _ => {}
        }
    }
}

/// `FuncOp::StrMid`
pub(crate) fn str_mid_const(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    let nelem: c_long = 0;

    let dest_str = (lParse.Nodes[this_node_idx]).value.data.text_mut_ptr();
    let dest_len = (lParse.Nodes[this_node_idx]).value.nelem as c_int;
    cstrmid(
        lParse,
        dest_str,
        dest_len,
        a.vals[0].data.text_mut_ptr(),
        a.vals[0].nelem as c_int,
        a.vals[1].data.lng() as c_int,
    );
}

/// `FuncOp::StrPos`
pub(crate) fn str_pos_const(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    unsafe {
        let res: *mut c_char = strstr(a.vals[0].data.text_mut_ptr(), a.vals[1].data.text_mut_ptr());
        if res.is_null() {
            (lParse.Nodes[this_node_idx]).value.data = NodeValue::Long(0);
        } else {
            (lParse.Nodes[this_node_idx]).value.data =
                NodeValue::Long(res.offset_from(a.vals[0].data.text_mut_ptr()) as c_long + 1);
        }
    }
}

/// `FuncOp::Row`
pub(crate) fn row_rows(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    unsafe {
        let mut row: c_long = lParse.nRows;

        loop {
            let fresh53 = row;
            row -= 1;
            if fresh53 == 0 {
                break;
            }
            *((lParse.Nodes[this_node_idx]).value.data.lng_buf()).offset(row as isize) =
                lParse.firstRow + row;
            *((lParse.Nodes[this_node_idx]).value.undef).offset(row as isize) = 0;
        }
    }
}

/// `FuncOp::Null`
pub(crate) fn null_rows(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    unsafe {
        let mut row: c_long = lParse.nRows;

        if (lParse.Nodes[this_node_idx]).ntype == ValueSort::Long {
            loop {
                let fresh54 = row;
                row -= 1;
                if fresh54 == 0 {
                    break;
                }
                *((lParse.Nodes[this_node_idx]).value.data.lng_buf()).offset(row as isize) = 0;
                *((lParse.Nodes[this_node_idx]).value.undef).offset(row as isize) = 1;
            }
        } else if (lParse.Nodes[this_node_idx]).ntype == ValueSort::String {
            loop {
                let fresh55 = row;
                row -= 1;
                if fresh55 == 0 {
                    break;
                }
                *(*((lParse.Nodes[this_node_idx]).value.data.str_buf()).offset(row as isize))
                    .offset(0) = 0;
                *((lParse.Nodes[this_node_idx]).value.undef).offset(row as isize) = 1;
            }
        }
    }
}

/// `FuncOp::AxisElem`
pub(crate) fn axis_elem_rows(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    unsafe {
        let elem: c_long = lParse.nRows * (lParse.Nodes[this_node_idx]).value.nelem;

        let mut ielem: c_long = 0;
        let mut iaxis: [c_long; 5] = [1, 1, 1, 1, 1];
        let ipos: c_long = a.vals[1].data.lng() - 1;
        let naxis: c_int = (lParse.Nodes[this_node_idx]).value.naxis;
        let mut j: c_int = 0;
        if ipos < 0 || ipos >= 5 as c_long {
            fits_parser_yyerror(
                lParse,
                cs!(c"AXISELEM(V,n) n value exceeded maximum dimension"),
            );
            free((lParse.Nodes[this_node_idx]).value.data.raw());
        } else {
            ielem = 0;
            while ielem < elem {
                *((lParse.Nodes[this_node_idx]).value.data.lng_buf()).offset(ielem as isize) =
                    iaxis[ipos as usize];
                *((lParse.Nodes[this_node_idx]).value.undef).offset(ielem as isize) = 0;
                iaxis[0] += 1;
                j = 0;
                while j < naxis {
                    if iaxis[j as usize] <= (lParse.Nodes[this_node_idx]).value.naxes[j as usize] {
                        break;
                    }
                    iaxis[j as usize] = 1;
                    if j < naxis - 1 {
                        iaxis[(j + 1) as usize] += 1;
                    }
                    j += 1;
                }
                ielem += 1;
            }
        }
    }
}

/// `FuncOp::ElemNum`
pub(crate) fn elem_num_rows(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    unsafe {
        let nelem: c_long = 0;
        let elem: c_long = lParse.nRows * (lParse.Nodes[this_node_idx]).value.nelem;

        let mut ielem_0: c_long = 0;
        let mut elemnum: c_long = 1;
        let j_0: c_int = 0;
        ielem_0 = 0;
        while ielem_0 < elem {
            *((lParse.Nodes[this_node_idx]).value.data.lng_buf()).offset(ielem_0 as isize) =
                elemnum;
            *((lParse.Nodes[this_node_idx]).value.undef).offset(ielem_0 as isize) = 0;
            elemnum += 1;
            if elemnum > (lParse.Nodes[this_node_idx]).value.nelem {
                elemnum = 1;
            }
            ielem_0 += 1;
        }
    }
}

/// `FuncOp::Rnd`
pub(crate) fn rnd_rows(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    unsafe {
        let mut elem: c_long = lParse.nRows * (lParse.Nodes[this_node_idx]).value.nelem;

        loop {
            let fresh56 = elem;
            elem -= 1;
            if fresh56 == 0 {
                break;
            }
            *((lParse.Nodes[this_node_idx]).value.data.dbl_buf()).offset(elem as isize) =
                simplerng_getuniform();
            *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize) = 0;
        }
    }
}

/// `FuncOp::GasRnd`
pub(crate) fn gas_rnd_rows(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    unsafe {
        let mut elem: c_long = lParse.nRows * (lParse.Nodes[this_node_idx]).value.nelem;

        loop {
            let fresh57 = elem;
            elem -= 1;
            if fresh57 == 0 {
                break;
            }
            *((lParse.Nodes[this_node_idx]).value.data.dbl_buf()).offset(elem as isize) =
                simplerng_getnorm();
            *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize) = 0;
        }
    }
}

/// `FuncOp::PoiRnd`
pub(crate) fn poi_rnd_rows(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    unsafe {
        let mut elem: c_long = lParse.nRows * (lParse.Nodes[this_node_idx]).value.nelem;

        if (lParse.Nodes[a.params[0]]).ntype == ValueSort::Double {
            if (lParse.Nodes[a.params[0]]).operation == Operation::Const {
                loop {
                    let fresh58 = elem;
                    elem -= 1;
                    if fresh58 == 0 {
                        break;
                    }
                    *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize) =
                        if a.vals[0].data.dbl() < 0.0 { 1 } else { 0 };
                    if *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize) == 0 {
                        *((lParse.Nodes[this_node_idx]).value.data.lng_buf())
                            .offset(elem as isize) =
                            c_long::from(simplerng_getpoisson(a.vals[0].data.dbl()));
                    }
                }
            } else {
                loop {
                    let fresh59 = elem;
                    elem -= 1;
                    if fresh59 == 0 {
                        break;
                    }
                    *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize) =
                        *((lParse.Nodes[a.params[0]]).value.undef).offset(elem as isize);
                    if *((lParse.Nodes[a.params[0]]).value.data.dbl_buf()).offset(elem as isize)
                        < 0.0
                    {
                        *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize) = 1;
                    }
                    if *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize) == 0 {
                        *((lParse.Nodes[this_node_idx]).value.data.lng_buf())
                            .offset(elem as isize) = c_long::from(simplerng_getpoisson(
                            *((lParse.Nodes[a.params[0]]).value.data.dbl_buf())
                                .offset(elem as isize),
                        ));
                    }
                }
            }
        } else if (lParse.Nodes[a.params[0]]).operation == Operation::Const {
            loop {
                let fresh60 = elem;
                elem -= 1;
                if fresh60 == 0 {
                    break;
                }
                *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize) =
                    if a.vals[0].data.lng() < 0 { 1 } else { 0 };
                if *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize) == 0 {
                    *((lParse.Nodes[this_node_idx]).value.data.lng_buf()).offset(elem as isize) =
                        c_long::from(simplerng_getpoisson(a.vals[0].data.lng() as c_double));
                }
            }
        } else {
            loop {
                let fresh61 = elem;
                elem -= 1;
                if fresh61 == 0 {
                    break;
                }
                *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize) =
                    *((lParse.Nodes[a.params[0]]).value.undef).offset(elem as isize);
                if *((lParse.Nodes[a.params[0]]).value.data.lng_buf()).offset(elem as isize) < 0 {
                    *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize) = 1;
                }
                if *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize) == 0 {
                    *((lParse.Nodes[this_node_idx]).value.data.lng_buf()).offset(elem as isize) =
                        c_long::from(simplerng_getpoisson(
                            *((lParse.Nodes[a.params[0]]).value.data.lng_buf())
                                .offset(elem as isize) as c_double,
                        ));
                }
            }
        }
    }
}

/// `FuncOp::Sum`
pub(crate) fn sum_rows(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    unsafe {
        let mut nelem: c_long = 0;
        let mut row: c_long = lParse.nRows;
        let mut elem: c_long = lParse.nRows * (lParse.Nodes[this_node_idx]).value.nelem;

        elem = row * (lParse.Nodes[a.params[0]]).value.nelem;
        if (lParse.Nodes[a.params[0]]).ntype == ValueSort::Boolean {
            loop {
                let fresh62 = row;
                row -= 1;
                if fresh62 == 0 {
                    break;
                }
                *((lParse.Nodes[this_node_idx]).value.data.lng_buf()).offset(row as isize) = 0;
                *((lParse.Nodes[this_node_idx]).value.undef).offset(row as isize) = 1;
                nelem = (lParse.Nodes[a.params[0]]).value.nelem;
                loop {
                    let fresh63 = nelem;
                    nelem -= 1;
                    if fresh63 == 0 {
                        break;
                    }
                    elem -= 1;

                    if *((lParse.Nodes[a.params[0]]).value.undef).offset(elem as isize) == 0 {
                        *((lParse.Nodes[this_node_idx]).value.data.lng_buf())
                            .offset(row as isize) += c_long::from(
                            if c_int::from(
                                *((lParse.Nodes[a.params[0]]).value.data.log_buf())
                                    .offset(elem as isize),
                            ) != 0
                            {
                                1
                            } else {
                                0
                            },
                        );
                        *((lParse.Nodes[this_node_idx]).value.undef).offset(row as isize) = 0;
                    }
                }
            }
        } else if (lParse.Nodes[a.params[0]]).ntype == ValueSort::Long {
            loop {
                let fresh64 = row;
                row -= 1;
                if fresh64 == 0 {
                    break;
                }
                *((lParse.Nodes[this_node_idx]).value.data.lng_buf()).offset(row as isize) = 0;
                *((lParse.Nodes[this_node_idx]).value.undef).offset(row as isize) = 1;
                nelem = (lParse.Nodes[a.params[0]]).value.nelem;
                loop {
                    let fresh65 = nelem;
                    nelem -= 1;
                    if fresh65 == 0 {
                        break;
                    }
                    elem -= 1;

                    if *((lParse.Nodes[a.params[0]]).value.undef).offset(elem as isize) == 0 {
                        *((lParse.Nodes[this_node_idx]).value.data.lng_buf())
                            .offset(row as isize) +=
                            *((lParse.Nodes[a.params[0]]).value.data.lng_buf())
                                .offset(elem as isize);
                        *((lParse.Nodes[this_node_idx]).value.undef).offset(row as isize) = 0;
                    }
                }
            }
        } else if (lParse.Nodes[a.params[0]]).ntype == ValueSort::Double {
            loop {
                let fresh66 = row;
                row -= 1;
                if fresh66 == 0 {
                    break;
                }
                *((lParse.Nodes[this_node_idx]).value.data.dbl_buf()).offset(row as isize) = 0.0;
                *((lParse.Nodes[this_node_idx]).value.undef).offset(row as isize) = 1;
                nelem = (lParse.Nodes[a.params[0]]).value.nelem;
                loop {
                    let fresh67 = nelem;
                    nelem -= 1;
                    if fresh67 == 0 {
                        break;
                    }
                    elem -= 1;

                    if *((lParse.Nodes[a.params[0]]).value.undef).offset(elem as isize) == 0 {
                        *((lParse.Nodes[this_node_idx]).value.data.dbl_buf())
                            .offset(row as isize) +=
                            *((lParse.Nodes[a.params[0]]).value.data.dbl_buf())
                                .offset(elem as isize);
                        *((lParse.Nodes[this_node_idx]).value.undef).offset(row as isize) = 0;
                    }
                }
            }
        } else {
            nelem = (lParse.Nodes[a.params[0]]).value.nelem;
            loop {
                let fresh68 = row;
                row -= 1;
                if fresh68 == 0 {
                    break;
                }
                let mut sptr1: *mut c_char =
                    *((lParse.Nodes[a.params[0]]).value.data.str_buf()).offset(row as isize);
                *((lParse.Nodes[this_node_idx]).value.data.lng_buf()).offset(row as isize) = 0;
                *((lParse.Nodes[this_node_idx]).value.undef).offset(row as isize) = 0;
                while *sptr1 != 0 {
                    if c_int::from(*sptr1) == '1' as i32 {
                        let fresh69 = &mut *((lParse.Nodes[this_node_idx]).value.data.lng_buf())
                            .offset(row as isize);
                        *fresh69 += 1;
                    }
                    sptr1 = sptr1.offset(1);
                }
            }
        }
    }
}

/// `FuncOp::Average`
pub(crate) fn average_rows(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    unsafe {
        let mut nelem: c_long = 0;
        let mut row: c_long = lParse.nRows;
        let mut elem: c_long = lParse.nRows * (lParse.Nodes[this_node_idx]).value.nelem;

        elem = row * (lParse.Nodes[a.params[0]]).value.nelem;
        if (lParse.Nodes[a.params[0]]).ntype == ValueSort::Long {
            loop {
                let fresh70 = row;
                row -= 1;
                if fresh70 == 0 {
                    break;
                }
                let mut count: c_int = 0;
                *((lParse.Nodes[this_node_idx]).value.data.dbl_buf()).offset(row as isize) = 0.0;
                nelem = (lParse.Nodes[a.params[0]]).value.nelem;
                loop {
                    let fresh71 = nelem;
                    nelem -= 1;
                    if fresh71 == 0 {
                        break;
                    }
                    elem -= 1;

                    if c_int::from(*((lParse.Nodes[a.params[0]]).value.undef).offset(elem as isize))
                        == 0
                    {
                        *((lParse.Nodes[this_node_idx]).value.data.dbl_buf())
                            .offset(row as isize) +=
                            *((lParse.Nodes[a.params[0]]).value.data.lng_buf())
                                .offset(elem as isize) as c_double;
                        count += 1;
                    }
                }
                if count == 0 {
                    *((lParse.Nodes[this_node_idx]).value.undef).offset(row as isize) = 1;
                } else {
                    *((lParse.Nodes[this_node_idx]).value.undef).offset(row as isize) = 0;
                    *((lParse.Nodes[this_node_idx]).value.data.dbl_buf()).offset(row as isize) /=
                        c_double::from(count);
                }
            }
        } else if (lParse.Nodes[a.params[0]]).ntype == ValueSort::Double {
            loop {
                let fresh72 = row;
                row -= 1;
                if fresh72 == 0 {
                    break;
                }
                let mut count_0: c_int = 0;
                *((lParse.Nodes[this_node_idx]).value.data.dbl_buf()).offset(row as isize) = 0.0;
                nelem = (lParse.Nodes[a.params[0]]).value.nelem;
                loop {
                    let fresh73 = nelem;
                    nelem -= 1;
                    if fresh73 == 0 {
                        break;
                    }
                    elem -= 1;

                    if c_int::from(*((lParse.Nodes[a.params[0]]).value.undef).offset(elem as isize))
                        == 0
                    {
                        *((lParse.Nodes[this_node_idx]).value.data.dbl_buf())
                            .offset(row as isize) +=
                            *((lParse.Nodes[a.params[0]]).value.data.dbl_buf())
                                .offset(elem as isize);
                        count_0 += 1;
                    }
                }
                if count_0 == 0 {
                    *((lParse.Nodes[this_node_idx]).value.undef).offset(row as isize) = 1;
                } else {
                    *((lParse.Nodes[this_node_idx]).value.undef).offset(row as isize) = 0;
                    *((lParse.Nodes[this_node_idx]).value.data.dbl_buf()).offset(row as isize) /=
                        c_double::from(count_0);
                }
            }
        }
    }
}

/// `FuncOp::Stddev`
pub(crate) fn stddev_rows(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    unsafe {
        let mut nelem: c_long = 0;
        let mut row: c_long = lParse.nRows;
        let mut elem: c_long = lParse.nRows * (lParse.Nodes[this_node_idx]).value.nelem;

        elem = row * (lParse.Nodes[a.params[0]]).value.nelem;
        if (lParse.Nodes[a.params[0]]).ntype == ValueSort::Long {
            loop {
                let fresh74 = row;
                row -= 1;
                if fresh74 == 0 {
                    break;
                }
                let mut count_1: c_int = 0;
                let mut sum: c_double = 0.0;
                let mut sum2: c_double = 0.0;
                nelem = (lParse.Nodes[a.params[0]]).value.nelem;
                loop {
                    let fresh75 = nelem;
                    nelem -= 1;
                    if fresh75 == 0 {
                        break;
                    }
                    elem -= 1;

                    if c_int::from(*((lParse.Nodes[a.params[0]]).value.undef).offset(elem as isize))
                        == 0
                    {
                        sum += *((lParse.Nodes[a.params[0]]).value.data.lng_buf())
                            .offset(elem as isize) as c_double;
                        count_1 += 1;
                    }
                }
                if count_1 > 1 {
                    sum /= c_double::from(count_1);
                    nelem = (lParse.Nodes[a.params[0]]).value.nelem;
                    elem += nelem;
                    loop {
                        let fresh76 = nelem;
                        nelem -= 1;
                        if fresh76 == 0 {
                            break;
                        }
                        elem -= 1;

                        if c_int::from(
                            *((lParse.Nodes[a.params[0]]).value.undef).offset(elem as isize),
                        ) == 0
                        {
                            let dx: c_double = *((lParse.Nodes[a.params[0]]).value.data.lng_buf())
                                .offset(elem as isize)
                                as c_double
                                - sum;
                            sum2 += dx * dx;
                        }
                    }
                    sum2 /= c_double::from(count_1) - c_double::from(1);
                    *((lParse.Nodes[this_node_idx]).value.undef).offset(row as isize) = 0;
                    *((lParse.Nodes[this_node_idx]).value.data.dbl_buf()).offset(row as isize) =
                        sqrt(sum2);
                } else {
                    *((lParse.Nodes[this_node_idx]).value.undef).offset(row as isize) = 0;
                    *((lParse.Nodes[this_node_idx]).value.data.dbl_buf()).offset(row as isize) =
                        0.0;
                }
            }
        } else if (lParse.Nodes[a.params[0]]).ntype == ValueSort::Double {
            loop {
                let fresh77 = row;
                row -= 1;
                if fresh77 == 0 {
                    break;
                }
                let mut count_2: c_int = 0;
                let mut sum_0: c_double = 0.0;
                let mut sum2_0: c_double = 0.0;
                nelem = (lParse.Nodes[a.params[0]]).value.nelem;
                loop {
                    let fresh78 = nelem;
                    nelem -= 1;
                    if fresh78 == 0 {
                        break;
                    }
                    elem -= 1;

                    if c_int::from(*((lParse.Nodes[a.params[0]]).value.undef).offset(elem as isize))
                        == 0
                    {
                        sum_0 += *((lParse.Nodes[a.params[0]]).value.data.dbl_buf())
                            .offset(elem as isize);
                        count_2 += 1;
                    }
                }
                if count_2 > 1 {
                    sum_0 /= c_double::from(count_2);
                    nelem = (lParse.Nodes[a.params[0]]).value.nelem;
                    elem += nelem;
                    loop {
                        let fresh79 = nelem;
                        nelem -= 1;
                        if fresh79 == 0 {
                            break;
                        }
                        elem -= 1;

                        if c_int::from(
                            *((lParse.Nodes[a.params[0]]).value.undef).offset(elem as isize),
                        ) == 0
                        {
                            let dx_0: c_double =
                                *((lParse.Nodes[a.params[0]]).value.data.dbl_buf())
                                    .offset(elem as isize)
                                    - sum_0;
                            sum2_0 += dx_0 * dx_0;
                        }
                    }
                    sum2_0 /= c_double::from(count_2) - c_double::from(1);
                    *((lParse.Nodes[this_node_idx]).value.undef).offset(row as isize) = 0;
                    *((lParse.Nodes[this_node_idx]).value.data.dbl_buf()).offset(row as isize) =
                        sqrt(sum2_0);
                } else {
                    *((lParse.Nodes[this_node_idx]).value.undef).offset(row as isize) = 0;
                    *((lParse.Nodes[this_node_idx]).value.data.dbl_buf()).offset(row as isize) =
                        0.0;
                }
            }
        }
    }
}

/// `FuncOp::Median`
pub(crate) fn median_rows(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    unsafe {
        let mut nelem: c_long = 0;
        let row: c_long = lParse.nRows;
        let mut elem: c_long = lParse.nRows * (lParse.Nodes[this_node_idx]).value.nelem;

        elem = row * (lParse.Nodes[a.params[0]]).value.nelem;
        nelem = (lParse.Nodes[a.params[0]]).value.nelem;
        if (lParse.Nodes[a.params[0]]).ntype == ValueSort::Long {
            let mut dptr: *mut c_long = (lParse.Nodes[a.params[0]]).value.data.lng_buf();
            let mut uptr: *mut c_char = (lParse.Nodes[a.params[0]]).value.undef;
            let mptr: *mut c_long = malloc(
                (::core::mem::size_of::<c_long>() as c_ulong)
                    .wrapping_mul(nelem as c_ulong)
                    .try_into()
                    .unwrap(),
            )
            .cast::<c_long>();
            let mut irow: c_int = 0;
            if mptr.is_null() {
                fits_parser_yyerror(
                    lParse,
                    cs!(c"Could not allocate temporary memory in median function"),
                );
                free((lParse.Nodes[this_node_idx]).value.data.raw());
            } else {
                irow = 0;
                while c_long::from(irow) < row {
                    let mut p: *mut c_long = mptr;
                    let mut nelem1: c_int = nelem as c_int;
                    loop {
                        let fresh80 = nelem1;
                        nelem1 -= 1;
                        if fresh80 == 0 {
                            break;
                        }
                        if c_int::from(*uptr) == 0 {
                            let fresh81 = p;
                            p = p.offset(1);
                            *fresh81 = *dptr;
                        }
                        dptr = dptr.offset(1);
                        uptr = uptr.offset(1);
                    }
                    nelem1 = p.offset_from(mptr) as c_long as c_int;
                    if nelem1 > 0 {
                        *((lParse.Nodes[this_node_idx]).value.undef).offset(irow as isize) = 0;
                        *((lParse.Nodes[this_node_idx]).value.data.lng_buf())
                            .offset(irow as isize) =
                            qselect_median(slice::from_raw_parts_mut(mptr, nelem1 as usize));
                    } else {
                        *((lParse.Nodes[this_node_idx]).value.undef).offset(irow as isize) = 1;
                        *((lParse.Nodes[this_node_idx]).value.data.lng_buf())
                            .offset(irow as isize) = 0;
                    }
                    irow += 1;
                }
                free(mptr.cast::<c_void>());
            }
        } else {
            let mut dptr_0: *mut c_double = (lParse.Nodes[a.params[0]]).value.data.dbl_buf();
            let mut uptr_0: *mut c_char = (lParse.Nodes[a.params[0]]).value.undef;
            let mptr_0: *mut c_double = malloc(
                (::core::mem::size_of::<c_double>() as c_ulong)
                    .wrapping_mul(nelem as c_ulong)
                    .try_into()
                    .unwrap(),
            )
            .cast::<c_double>();
            let mut irow_0: c_int = 0;
            if mptr_0.is_null() {
                fits_parser_yyerror(
                    lParse,
                    cs!(c"Could not allocate temporary memory in median function"),
                );
                free((lParse.Nodes[this_node_idx]).value.data.raw());
            } else {
                irow_0 = 0;
                while c_long::from(irow_0) < row {
                    let mut p_0: *mut c_double = mptr_0;
                    let mut nelem1_0: c_int = nelem as c_int;
                    loop {
                        let fresh82 = nelem1_0;
                        nelem1_0 -= 1;
                        if fresh82 == 0 {
                            break;
                        }
                        if c_int::from(*uptr_0) == 0 {
                            let fresh83 = p_0;
                            p_0 = p_0.offset(1);
                            *fresh83 = *dptr_0;
                        }
                        dptr_0 = dptr_0.offset(1);
                        uptr_0 = uptr_0.offset(1);
                    }
                    nelem1_0 = p_0.offset_from(mptr_0) as c_long as c_int;
                    if nelem1_0 > 0 {
                        *((lParse.Nodes[this_node_idx]).value.undef).offset(irow_0 as isize) = 0;
                        *((lParse.Nodes[this_node_idx]).value.data.dbl_buf())
                            .offset(irow_0 as isize) =
                            qselect_median(slice::from_raw_parts_mut(mptr_0, nelem1_0 as usize));
                    } else {
                        *((lParse.Nodes[this_node_idx]).value.undef).offset(irow_0 as isize) = 1;
                        *((lParse.Nodes[this_node_idx]).value.data.dbl_buf())
                            .offset(irow_0 as isize) = 0.0;
                    }
                    irow_0 += 1;
                }
                free(mptr_0.cast::<c_void>());
            }
        }
    }
}

/// `FuncOp::Abs`
pub(crate) fn abs_rows(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    unsafe {
        let mut ival: c_long = 0;
        let mut dval: c_double = 0.0;
        let mut elem: c_long = lParse.nRows * (lParse.Nodes[this_node_idx]).value.nelem;

        if (lParse.Nodes[a.params[0]]).ntype == ValueSort::Double {
            loop {
                let fresh84 = elem;
                elem -= 1;
                if fresh84 == 0 {
                    break;
                }
                dval = *((lParse.Nodes[a.params[0]]).value.data.dbl_buf()).offset(elem as isize);
                *((lParse.Nodes[this_node_idx]).value.data.dbl_buf()).offset(elem as isize) =
                    if dval > 0.0 { dval } else { -dval };
                *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize) =
                    *((lParse.Nodes[a.params[0]]).value.undef).offset(elem as isize);
            }
        } else {
            loop {
                let fresh85 = elem;
                elem -= 1;
                if fresh85 == 0 {
                    break;
                }
                ival = *((lParse.Nodes[a.params[0]]).value.data.lng_buf()).offset(elem as isize);
                *((lParse.Nodes[this_node_idx]).value.data.lng_buf()).offset(elem as isize) =
                    if ival > 0 { ival } else { -ival };
                *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize) =
                    *((lParse.Nodes[a.params[0]]).value.undef).offset(elem as isize);
            }
        }
    }
}

/// `FuncOp::NonNull`
pub(crate) fn non_null_rows(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    unsafe {
        let mut nelem: c_long = 0;
        let mut row: c_long = lParse.nRows;
        let mut elem: c_long = lParse.nRows * (lParse.Nodes[this_node_idx]).value.nelem;

        nelem = (lParse.Nodes[a.params[0]]).value.nelem;
        if (lParse.Nodes[a.params[0]]).ntype == ValueSort::String {
            nelem = 1;
        }
        elem = row * nelem;
        loop {
            let fresh86 = row;
            row -= 1;
            if fresh86 == 0 {
                break;
            }
            let mut nelem1_1: c_int = nelem as c_int;
            *((lParse.Nodes[this_node_idx]).value.undef).offset(row as isize) = 0;
            *((lParse.Nodes[this_node_idx]).value.data.lng_buf()).offset(row as isize) = 0;
            loop {
                let fresh87 = nelem1_1;
                nelem1_1 -= 1;
                if fresh87 == 0 {
                    break;
                }
                elem -= 1;

                if c_int::from(*((lParse.Nodes[a.params[0]]).value.undef).offset(elem as isize))
                    == 0
                {
                    let fresh88 = &mut *((lParse.Nodes[this_node_idx]).value.data.lng_buf())
                        .offset(row as isize);
                    *fresh88 += 1;
                }
            }
        }
    }
}

/// `FuncOp::IsNull`
pub(crate) fn is_null_rows(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    unsafe {
        let row: c_long = lParse.nRows;
        let mut elem: c_long = lParse.nRows * (lParse.Nodes[this_node_idx]).value.nelem;

        if (lParse.Nodes[a.params[0]]).ntype == ValueSort::String {
            elem = row;
        }
        loop {
            let fresh89 = elem;
            elem -= 1;
            if fresh89 == 0 {
                break;
            }
            *((lParse.Nodes[this_node_idx]).value.data.log_buf()).offset(elem as isize) =
                *((lParse.Nodes[a.params[0]]).value.undef).offset(elem as isize);
            *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize) = 0;
        }
    }
}

/// `FuncOp::DefNull`
pub(crate) fn def_null_rows(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    unsafe {
        let mut i: c_int = 0;
        let mut nelem: c_long = 0;
        let mut row: c_long = lParse.nRows;
        let mut elem: c_long = lParse.nRows * (lParse.Nodes[this_node_idx]).value.nelem;

        match (lParse.Nodes[this_node_idx]).ntype {
            ValueSort::Boolean => loop {
                let fresh90 = row;
                row -= 1;
                if fresh90 == 0 {
                    break;
                }
                nelem = (lParse.Nodes[this_node_idx]).value.nelem;
                loop {
                    let fresh91 = nelem;
                    nelem -= 1;
                    if fresh91 == 0 {
                        break;
                    }
                    elem -= 1;

                    i = 2;
                    loop {
                        let fresh92 = i;
                        i -= 1;
                        if fresh92 == 0 {
                            break;
                        }
                        if a.vector[i as usize] > 1 {
                            a.nulls[i as usize] =
                                *((lParse.Nodes[a.params[i as usize]]).value.undef)
                                    .offset(elem as isize);
                            a.vals[i as usize].data = NodeValue::Logical(
                                *((lParse.Nodes[a.params[i as usize]]).value.data.log_buf())
                                    .offset(elem as isize),
                            );
                        } else if a.vector[i as usize] != 0 {
                            a.nulls[i as usize] =
                                *((lParse.Nodes[a.params[i as usize]]).value.undef)
                                    .offset(row as isize);
                            a.vals[i as usize].data = NodeValue::Logical(
                                *((lParse.Nodes[a.params[i as usize]]).value.data.log_buf())
                                    .offset(row as isize),
                            );
                        }
                    }
                    if a.nulls[0] != 0 {
                        *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize) =
                            a.nulls[1];
                        *((lParse.Nodes[this_node_idx]).value.data.log_buf())
                            .offset(elem as isize) = a.vals[1].data.log();
                    } else {
                        *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize) = 0;
                        *((lParse.Nodes[this_node_idx]).value.data.log_buf())
                            .offset(elem as isize) = a.vals[0].data.log();
                    }
                }
            },
            ValueSort::Long => loop {
                let fresh93 = row;
                row -= 1;
                if fresh93 == 0 {
                    break;
                }
                nelem = (lParse.Nodes[this_node_idx]).value.nelem;
                loop {
                    let fresh94 = nelem;
                    nelem -= 1;
                    if fresh94 == 0 {
                        break;
                    }
                    elem -= 1;

                    i = 2;
                    loop {
                        let fresh95 = i;
                        i -= 1;
                        if fresh95 == 0 {
                            break;
                        }
                        if a.vector[i as usize] > 1 {
                            a.nulls[i as usize] =
                                *((lParse.Nodes[a.params[i as usize]]).value.undef)
                                    .offset(elem as isize);
                            a.vals[i as usize].data = NodeValue::Long(
                                *((lParse.Nodes[a.params[i as usize]]).value.data.lng_buf())
                                    .offset(elem as isize),
                            );
                        } else if a.vector[i as usize] != 0 {
                            a.nulls[i as usize] =
                                *((lParse.Nodes[a.params[i as usize]]).value.undef)
                                    .offset(row as isize);
                            a.vals[i as usize].data = NodeValue::Long(
                                *((lParse.Nodes[a.params[i as usize]]).value.data.lng_buf())
                                    .offset(row as isize),
                            );
                        }
                    }
                    if a.nulls[0] != 0 {
                        *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize) =
                            a.nulls[1];
                        *((lParse.Nodes[this_node_idx]).value.data.lng_buf())
                            .offset(elem as isize) = a.vals[1].data.lng();
                    } else {
                        *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize) = 0;
                        *((lParse.Nodes[this_node_idx]).value.data.lng_buf())
                            .offset(elem as isize) = a.vals[0].data.lng();
                    }
                }
            },
            ValueSort::Double => loop {
                let fresh96 = row;
                row -= 1;
                if fresh96 == 0 {
                    break;
                }
                nelem = (lParse.Nodes[this_node_idx]).value.nelem;
                loop {
                    let fresh97 = nelem;
                    nelem -= 1;
                    if fresh97 == 0 {
                        break;
                    }
                    elem -= 1;

                    i = 2;
                    loop {
                        let fresh98 = i;
                        i -= 1;
                        if fresh98 == 0 {
                            break;
                        }
                        if a.vector[i as usize] > 1 {
                            a.nulls[i as usize] =
                                *((lParse.Nodes[a.params[i as usize]]).value.undef)
                                    .offset(elem as isize);
                            a.vals[i as usize].data = NodeValue::Double(
                                *((lParse.Nodes[a.params[i as usize]]).value.data.dbl_buf())
                                    .offset(elem as isize),
                            );
                        } else if a.vector[i as usize] != 0 {
                            a.nulls[i as usize] =
                                *((lParse.Nodes[a.params[i as usize]]).value.undef)
                                    .offset(row as isize);
                            a.vals[i as usize].data = NodeValue::Double(
                                *((lParse.Nodes[a.params[i as usize]]).value.data.dbl_buf())
                                    .offset(row as isize),
                            );
                        }
                    }
                    if a.nulls[0] != 0 {
                        *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize) =
                            a.nulls[1];
                        *((lParse.Nodes[this_node_idx]).value.data.dbl_buf())
                            .offset(elem as isize) = a.vals[1].data.dbl();
                    } else {
                        *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize) = 0;
                        *((lParse.Nodes[this_node_idx]).value.data.dbl_buf())
                            .offset(elem as isize) = a.vals[0].data.dbl();
                    }
                }
            },
            ValueSort::String => loop {
                let fresh99 = row;
                row -= 1;
                if fresh99 == 0 {
                    break;
                }
                i = 2;
                loop {
                    let fresh100 = i;
                    i -= 1;
                    if fresh100 == 0 {
                        break;
                    }
                    if a.vector[i as usize] != 0 {
                        a.nulls[i as usize] = *((lParse.Nodes[a.params[i as usize]]).value.undef)
                            .offset(row as isize);
                        strcpy(
                            a.vals[i as usize].data.text_mut_ptr(),
                            *((lParse.Nodes[a.params[i as usize]]).value.data.str_buf())
                                .offset(row as isize),
                        );
                    }
                }
                if a.nulls[0] != 0 {
                    *((lParse.Nodes[this_node_idx]).value.undef).offset(row as isize) = a.nulls[1];
                    strcpy(
                        *((lParse.Nodes[this_node_idx]).value.data.str_buf()).offset(row as isize),
                        a.vals[1].data.text_mut_ptr(),
                    );
                } else {
                    *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize) = 0;
                    strcpy(
                        *((lParse.Nodes[this_node_idx]).value.data.str_buf()).offset(row as isize),
                        a.vals[0].data.text_mut_ptr(),
                    );
                }
            },
            _ => {}
        }
    }
}

/// `FuncOp::SetNull`
pub(crate) fn set_null_rows(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    unsafe {
        let mut elem: c_long = lParse.nRows * (lParse.Nodes[this_node_idx]).value.nelem;

        match (lParse.Nodes[this_node_idx]).ntype {
            ValueSort::Long => loop {
                let fresh101 = elem;
                elem -= 1;
                if fresh101 == 0 {
                    break;
                }
                if (lParse.Nodes[a.params[1]]).value.data.lng()
                    == *((lParse.Nodes[a.params[0]]).value.data.lng_buf()).offset(elem as isize)
                {
                    *((lParse.Nodes[this_node_idx]).value.data.lng_buf()).offset(elem as isize) = 0;
                    *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize) = 1;
                } else {
                    *((lParse.Nodes[this_node_idx]).value.data.lng_buf()).offset(elem as isize) =
                        *((lParse.Nodes[a.params[0]]).value.data.lng_buf()).offset(elem as isize);
                    *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize) =
                        *((lParse.Nodes[a.params[0]]).value.undef).offset(elem as isize);
                }
            },
            ValueSort::Double => loop {
                let fresh102 = elem;
                elem -= 1;
                if fresh102 == 0 {
                    break;
                }
                if (lParse.Nodes[a.params[1]]).value.data.dbl()
                    == *((lParse.Nodes[a.params[0]]).value.data.dbl_buf()).offset(elem as isize)
                {
                    *((lParse.Nodes[this_node_idx]).value.data.dbl_buf()).offset(elem as isize) =
                        0.0;
                    *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize) = 1;
                } else {
                    *((lParse.Nodes[this_node_idx]).value.data.dbl_buf()).offset(elem as isize) =
                        *((lParse.Nodes[a.params[0]]).value.data.dbl_buf()).offset(elem as isize);
                    *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize) =
                        *((lParse.Nodes[a.params[0]]).value.undef).offset(elem as isize);
                }
            },
            _ => {}
        }
    }
}

/// `FuncOp::Sin`
pub(crate) fn sin_rows(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    unsafe {
        let mut elem: c_long = lParse.nRows * (lParse.Nodes[this_node_idx]).value.nelem;

        loop {
            let fresh103 = elem;
            elem -= 1;
            if fresh103 == 0 {
                break;
            }
            let fresh104 = &mut *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize);
            *fresh104 = *((lParse.Nodes[a.params[0]]).value.undef).offset(elem as isize);
            if *fresh104 == 0 {
                *((lParse.Nodes[this_node_idx]).value.data.dbl_buf()).offset(elem as isize) =
                    sin(*((lParse.Nodes[a.params[0]]).value.data.dbl_buf()).offset(elem as isize));
            }
        }
    }
}

/// `FuncOp::Cos`
pub(crate) fn cos_rows(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    unsafe {
        let mut elem: c_long = lParse.nRows * (lParse.Nodes[this_node_idx]).value.nelem;

        loop {
            let fresh105 = elem;
            elem -= 1;
            if fresh105 == 0 {
                break;
            }
            let fresh106 = &mut *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize);
            *fresh106 = *((lParse.Nodes[a.params[0]]).value.undef).offset(elem as isize);
            if *fresh106 == 0 {
                *((lParse.Nodes[this_node_idx]).value.data.dbl_buf()).offset(elem as isize) =
                    cos(*((lParse.Nodes[a.params[0]]).value.data.dbl_buf()).offset(elem as isize));
            }
        }
    }
}

/// `FuncOp::Tan`
pub(crate) fn tan_rows(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    unsafe {
        let mut elem: c_long = lParse.nRows * (lParse.Nodes[this_node_idx]).value.nelem;

        loop {
            let fresh107 = elem;
            elem -= 1;
            if fresh107 == 0 {
                break;
            }
            let fresh108 = &mut *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize);
            *fresh108 = *((lParse.Nodes[a.params[0]]).value.undef).offset(elem as isize);
            if *fresh108 == 0 {
                *((lParse.Nodes[this_node_idx]).value.data.dbl_buf()).offset(elem as isize) =
                    tan(*((lParse.Nodes[a.params[0]]).value.data.dbl_buf()).offset(elem as isize));
            }
        }
    }
}

/// `FuncOp::Asin`
pub(crate) fn asin_rows(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    unsafe {
        let mut dval: c_double = 0.0;
        let mut elem: c_long = lParse.nRows * (lParse.Nodes[this_node_idx]).value.nelem;

        loop {
            let fresh109 = elem;
            elem -= 1;
            if fresh109 == 0 {
                break;
            }
            let fresh110 = &mut *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize);
            *fresh110 = *((lParse.Nodes[a.params[0]]).value.undef).offset(elem as isize);
            if *fresh110 == 0 {
                dval = *((lParse.Nodes[a.params[0]]).value.data.dbl_buf()).offset(elem as isize);
                if dval < -1.0 || dval > 1.0 {
                    *((lParse.Nodes[this_node_idx]).value.data.dbl_buf()).offset(elem as isize) =
                        0.0;
                    *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize) = 1;
                } else {
                    *((lParse.Nodes[this_node_idx]).value.data.dbl_buf()).offset(elem as isize) =
                        asin(dval);
                }
            }
        }
    }
}

/// `FuncOp::Acos`
pub(crate) fn acos_rows(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    unsafe {
        let mut dval: c_double = 0.0;
        let mut elem: c_long = lParse.nRows * (lParse.Nodes[this_node_idx]).value.nelem;

        loop {
            let fresh111 = elem;
            elem -= 1;
            if fresh111 == 0 {
                break;
            }
            let fresh112 = &mut *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize);
            *fresh112 = *((lParse.Nodes[a.params[0]]).value.undef).offset(elem as isize);
            if *fresh112 == 0 {
                dval = *((lParse.Nodes[a.params[0]]).value.data.dbl_buf()).offset(elem as isize);
                if dval < -1.0 || dval > 1.0 {
                    *((lParse.Nodes[this_node_idx]).value.data.dbl_buf()).offset(elem as isize) =
                        0.0;
                    *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize) = 1;
                } else {
                    *((lParse.Nodes[this_node_idx]).value.data.dbl_buf()).offset(elem as isize) =
                        acos(dval);
                }
            }
        }
    }
}

/// `FuncOp::Atan`
pub(crate) fn atan_rows(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    unsafe {
        let mut dval: c_double = 0.0;
        let mut elem: c_long = lParse.nRows * (lParse.Nodes[this_node_idx]).value.nelem;

        loop {
            let fresh113 = elem;
            elem -= 1;
            if fresh113 == 0 {
                break;
            }
            let fresh114 = &mut *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize);
            *fresh114 = *((lParse.Nodes[a.params[0]]).value.undef).offset(elem as isize);
            if *fresh114 == 0 {
                dval = *((lParse.Nodes[a.params[0]]).value.data.dbl_buf()).offset(elem as isize);
                *((lParse.Nodes[this_node_idx]).value.data.dbl_buf()).offset(elem as isize) =
                    atan(dval);
            }
        }
    }
}

/// `FuncOp::Sinh`
pub(crate) fn sinh_rows(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    unsafe {
        let mut elem: c_long = lParse.nRows * (lParse.Nodes[this_node_idx]).value.nelem;

        loop {
            let fresh115 = elem;
            elem -= 1;
            if fresh115 == 0 {
                break;
            }
            let fresh116 = &mut *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize);
            *fresh116 = *((lParse.Nodes[a.params[0]]).value.undef).offset(elem as isize);
            if *fresh116 == 0 {
                *((lParse.Nodes[this_node_idx]).value.data.dbl_buf()).offset(elem as isize) =
                    sinh(*((lParse.Nodes[a.params[0]]).value.data.dbl_buf()).offset(elem as isize));
            }
        }
    }
}

/// `FuncOp::Cosh`
pub(crate) fn cosh_rows(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    unsafe {
        let mut elem: c_long = lParse.nRows * (lParse.Nodes[this_node_idx]).value.nelem;

        loop {
            let fresh117 = elem;
            elem -= 1;
            if fresh117 == 0 {
                break;
            }
            let fresh118 = &mut *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize);
            *fresh118 = *((lParse.Nodes[a.params[0]]).value.undef).offset(elem as isize);
            if *fresh118 == 0 {
                *((lParse.Nodes[this_node_idx]).value.data.dbl_buf()).offset(elem as isize) =
                    cosh(*((lParse.Nodes[a.params[0]]).value.data.dbl_buf()).offset(elem as isize));
            }
        }
    }
}

/// `FuncOp::Tanh`
pub(crate) fn tanh_rows(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    unsafe {
        let mut elem: c_long = lParse.nRows * (lParse.Nodes[this_node_idx]).value.nelem;

        loop {
            let fresh119 = elem;
            elem -= 1;
            if fresh119 == 0 {
                break;
            }
            let fresh120 = &mut *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize);
            *fresh120 = *((lParse.Nodes[a.params[0]]).value.undef).offset(elem as isize);
            if *fresh120 == 0 {
                *((lParse.Nodes[this_node_idx]).value.data.dbl_buf()).offset(elem as isize) =
                    tanh(*((lParse.Nodes[a.params[0]]).value.data.dbl_buf()).offset(elem as isize));
            }
        }
    }
}

/// `FuncOp::Exp`
pub(crate) fn exp_rows(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    unsafe {
        let mut dval: c_double = 0.0;
        let mut elem: c_long = lParse.nRows * (lParse.Nodes[this_node_idx]).value.nelem;

        loop {
            let fresh121 = elem;
            elem -= 1;
            if fresh121 == 0 {
                break;
            }
            let fresh122 = &mut *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize);
            *fresh122 = *((lParse.Nodes[a.params[0]]).value.undef).offset(elem as isize);
            if *fresh122 == 0 {
                dval = *((lParse.Nodes[a.params[0]]).value.data.dbl_buf()).offset(elem as isize);
                *((lParse.Nodes[this_node_idx]).value.data.dbl_buf()).offset(elem as isize) =
                    exp(dval);
            }
        }
    }
}

/// `FuncOp::Log`
pub(crate) fn log_rows(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    unsafe {
        let mut dval: c_double = 0.0;
        let mut elem: c_long = lParse.nRows * (lParse.Nodes[this_node_idx]).value.nelem;

        loop {
            let fresh123 = elem;
            elem -= 1;
            if fresh123 == 0 {
                break;
            }
            let fresh124 = &mut *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize);
            *fresh124 = *((lParse.Nodes[a.params[0]]).value.undef).offset(elem as isize);
            if *fresh124 == 0 {
                dval = *((lParse.Nodes[a.params[0]]).value.data.dbl_buf()).offset(elem as isize);
                if dval <= 0.0 {
                    *((lParse.Nodes[this_node_idx]).value.data.dbl_buf()).offset(elem as isize) =
                        0.0;
                    *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize) = 1;
                } else {
                    *((lParse.Nodes[this_node_idx]).value.data.dbl_buf()).offset(elem as isize) =
                        log(dval);
                }
            }
        }
    }
}

/// `FuncOp::Log10`
pub(crate) fn log10_rows(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    unsafe {
        let mut dval: c_double = 0.0;
        let mut elem: c_long = lParse.nRows * (lParse.Nodes[this_node_idx]).value.nelem;

        loop {
            let fresh125 = elem;
            elem -= 1;
            if fresh125 == 0 {
                break;
            }
            let fresh126 = &mut *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize);
            *fresh126 = *((lParse.Nodes[a.params[0]]).value.undef).offset(elem as isize);
            if *fresh126 == 0 {
                dval = *((lParse.Nodes[a.params[0]]).value.data.dbl_buf()).offset(elem as isize);
                if dval <= 0.0 {
                    *((lParse.Nodes[this_node_idx]).value.data.dbl_buf()).offset(elem as isize) =
                        0.0;
                    *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize) = 1;
                } else {
                    *((lParse.Nodes[this_node_idx]).value.data.dbl_buf()).offset(elem as isize) =
                        log10(dval);
                }
            }
        }
    }
}

/// `FuncOp::Sqrt`
pub(crate) fn sqrt_rows(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    unsafe {
        let mut dval: c_double = 0.0;
        let mut elem: c_long = lParse.nRows * (lParse.Nodes[this_node_idx]).value.nelem;

        loop {
            let fresh127 = elem;
            elem -= 1;
            if fresh127 == 0 {
                break;
            }
            let fresh128 = &mut *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize);
            *fresh128 = *((lParse.Nodes[a.params[0]]).value.undef).offset(elem as isize);
            if *fresh128 == 0 {
                dval = *((lParse.Nodes[a.params[0]]).value.data.dbl_buf()).offset(elem as isize);
                if dval < 0.0 {
                    *((lParse.Nodes[this_node_idx]).value.data.dbl_buf()).offset(elem as isize) =
                        0.0;
                    *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize) = 1;
                } else {
                    *((lParse.Nodes[this_node_idx]).value.data.dbl_buf()).offset(elem as isize) =
                        sqrt(dval);
                }
            }
        }
    }
}

/// `FuncOp::Ceil`
pub(crate) fn ceil_rows(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    unsafe {
        let mut elem: c_long = lParse.nRows * (lParse.Nodes[this_node_idx]).value.nelem;

        loop {
            let fresh129 = elem;
            elem -= 1;
            if fresh129 == 0 {
                break;
            }
            let fresh130 = &mut *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize);
            *fresh130 = *((lParse.Nodes[a.params[0]]).value.undef).offset(elem as isize);
            if *fresh130 == 0 {
                *((lParse.Nodes[this_node_idx]).value.data.dbl_buf()).offset(elem as isize) =
                    ceil(*((lParse.Nodes[a.params[0]]).value.data.dbl_buf()).offset(elem as isize));
            }
        }
    }
}

/// `FuncOp::Floor`
pub(crate) fn floor_rows(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    unsafe {
        let mut elem: c_long = lParse.nRows * (lParse.Nodes[this_node_idx]).value.nelem;

        loop {
            let fresh131 = elem;
            elem -= 1;
            if fresh131 == 0 {
                break;
            }
            let fresh132 = &mut *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize);
            *fresh132 = *((lParse.Nodes[a.params[0]]).value.undef).offset(elem as isize);
            if *fresh132 == 0 {
                *((lParse.Nodes[this_node_idx]).value.data.dbl_buf()).offset(elem as isize) = floor(
                    *((lParse.Nodes[a.params[0]]).value.data.dbl_buf()).offset(elem as isize),
                );
            }
        }
    }
}

/// `FuncOp::Round`
pub(crate) fn round_rows(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    unsafe {
        let mut elem: c_long = lParse.nRows * (lParse.Nodes[this_node_idx]).value.nelem;

        loop {
            let fresh133 = elem;
            elem -= 1;
            if fresh133 == 0 {
                break;
            }
            let fresh134 = &mut *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize);
            *fresh134 = *((lParse.Nodes[a.params[0]]).value.undef).offset(elem as isize);
            if *fresh134 == 0 {
                *((lParse.Nodes[this_node_idx]).value.data.dbl_buf()).offset(elem as isize) = floor(
                    *((lParse.Nodes[a.params[0]]).value.data.dbl_buf()).offset(elem as isize) + 0.5,
                );
            }
        }
    }
}

/// `FuncOp::Atan2`
pub(crate) fn atan2_rows(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    unsafe {
        let mut i: c_int = 0;
        let mut nelem: c_long = 0;
        let mut row: c_long = lParse.nRows;
        let mut elem: c_long = lParse.nRows * (lParse.Nodes[this_node_idx]).value.nelem;

        loop {
            let fresh135 = row;
            row -= 1;
            if fresh135 == 0 {
                break;
            }
            nelem = (lParse.Nodes[this_node_idx]).value.nelem;
            loop {
                let fresh136 = nelem;
                nelem -= 1;
                if fresh136 == 0 {
                    break;
                }
                elem -= 1;

                i = 2;
                loop {
                    let fresh137 = i;
                    i -= 1;
                    if fresh137 == 0 {
                        break;
                    }
                    if a.vector[i as usize] > 1 {
                        a.vals[i as usize].data = NodeValue::Double(
                            *((lParse.Nodes[a.params[i as usize]]).value.data.dbl_buf())
                                .offset(elem as isize),
                        );
                        a.nulls[i as usize] = *((lParse.Nodes[a.params[i as usize]]).value.undef)
                            .offset(elem as isize);
                    } else if a.vector[i as usize] != 0 {
                        a.vals[i as usize].data = NodeValue::Double(
                            *((lParse.Nodes[a.params[i as usize]]).value.data.dbl_buf())
                                .offset(row as isize),
                        );
                        a.nulls[i as usize] = *((lParse.Nodes[a.params[i as usize]]).value.undef)
                            .offset(row as isize);
                    }
                }
                let fresh138 =
                    &mut *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize);
                *fresh138 = if c_int::from(a.nulls[0]) != 0 || c_int::from(a.nulls[1]) != 0 {
                    1
                } else {
                    0
                };
                if *fresh138 == 0 {
                    *((lParse.Nodes[this_node_idx]).value.data.dbl_buf()).offset(elem as isize) =
                        atan2(a.vals[0].data.dbl(), a.vals[1].data.dbl());
                }
            }
        }
    }
}

/// `FuncOp::AngSep`
pub(crate) fn ang_sep_rows(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    unsafe {
        let mut i: c_int = 0;
        let mut nelem: c_long = 0;
        let mut row: c_long = lParse.nRows;
        let mut elem: c_long = lParse.nRows * (lParse.Nodes[this_node_idx]).value.nelem;

        loop {
            let fresh139 = row;
            row -= 1;
            if fresh139 == 0 {
                break;
            }
            nelem = (lParse.Nodes[this_node_idx]).value.nelem;
            loop {
                let fresh140 = nelem;
                nelem -= 1;
                if fresh140 == 0 {
                    break;
                }
                elem -= 1;
                i = 4 as c_int;
                loop {
                    let fresh141 = i;
                    i -= 1;
                    if fresh141 == 0 {
                        break;
                    }
                    if a.vector[i as usize] > 1 {
                        a.vals[i as usize].data = NodeValue::Double(
                            *((lParse.Nodes[a.params[i as usize]]).value.data.dbl_buf())
                                .offset(elem as isize),
                        );
                        a.nulls[i as usize] = *((lParse.Nodes[a.params[i as usize]]).value.undef)
                            .offset(elem as isize);
                    } else if a.vector[i as usize] != 0 {
                        a.vals[i as usize].data = NodeValue::Double(
                            *((lParse.Nodes[a.params[i as usize]]).value.data.dbl_buf())
                                .offset(row as isize),
                        );
                        a.nulls[i as usize] = *((lParse.Nodes[a.params[i as usize]]).value.undef)
                            .offset(row as isize);
                    }
                }
                let fresh142 =
                    &mut *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize);
                *fresh142 = c_int::from(
                    c_int::from(a.nulls[0]) != 0
                        || c_int::from(a.nulls[1]) != 0
                        || c_int::from(a.nulls[2]) != 0
                        || c_int::from(a.nulls[3]) != 0,
                ) as c_char;
                if *fresh142 == 0 {
                    *((lParse.Nodes[this_node_idx]).value.data.dbl_buf()).offset(elem as isize) =
                        angsep_calc(
                            a.vals[0].data.dbl(),
                            a.vals[1].data.dbl(),
                            a.vals[2].data.dbl(),
                            a.vals[3].data.dbl(),
                        );
                }
            }
        }
    }
}

/// `FuncOp::Min1`
pub(crate) fn min1_rows(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    unsafe {
        let mut valInit: c_int = 0;
        let mut nelem: c_long = 0;
        let mut row: c_long = lParse.nRows;
        let mut elem: c_long = lParse.nRows * (lParse.Nodes[this_node_idx]).value.nelem;

        elem = row * (lParse.Nodes[a.params[0]]).value.nelem;
        if (lParse.Nodes[this_node_idx]).ntype == ValueSort::Long {
            let mut minVal: c_long = 0;
            loop {
                let fresh143 = row;
                row -= 1;
                if fresh143 == 0 {
                    break;
                }
                valInit = 1;
                *((lParse.Nodes[this_node_idx]).value.undef).offset(row as isize) = 1;
                nelem = (lParse.Nodes[a.params[0]]).value.nelem;
                loop {
                    let fresh144 = nelem;
                    nelem -= 1;
                    if fresh144 == 0 {
                        break;
                    }
                    elem -= 1;
                    if *((lParse.Nodes[a.params[0]]).value.undef).offset(elem as isize) == 0 {
                        if valInit != 0 {
                            valInit = 0;
                            minVal = *((lParse.Nodes[a.params[0]]).value.data.lng_buf())
                                .offset(elem as isize);
                        } else {
                            minVal = if minVal
                                < *((lParse.Nodes[a.params[0]]).value.data.lng_buf())
                                    .offset(elem as isize)
                            {
                                minVal
                            } else {
                                *((lParse.Nodes[a.params[0]]).value.data.lng_buf())
                                    .offset(elem as isize)
                            };
                        }
                        *((lParse.Nodes[this_node_idx]).value.undef).offset(row as isize) = 0;
                    }
                }
                *((lParse.Nodes[this_node_idx]).value.data.lng_buf()).offset(row as isize) = minVal;
            }
        } else if (lParse.Nodes[this_node_idx]).ntype == ValueSort::Double {
            let mut minVal_0: c_double = 0.0;
            loop {
                let fresh145 = row;
                row -= 1;
                if fresh145 == 0 {
                    break;
                }
                valInit = 1;
                *((lParse.Nodes[this_node_idx]).value.undef).offset(row as isize) = 1;
                nelem = (lParse.Nodes[a.params[0]]).value.nelem;
                loop {
                    let fresh146 = nelem;
                    nelem -= 1;
                    if fresh146 == 0 {
                        break;
                    }
                    elem -= 1;
                    if *((lParse.Nodes[a.params[0]]).value.undef).offset(elem as isize) == 0 {
                        if valInit != 0 {
                            valInit = 0;
                            minVal_0 = *((lParse.Nodes[a.params[0]]).value.data.dbl_buf())
                                .offset(elem as isize);
                        } else {
                            minVal_0 = if minVal_0
                                < *((lParse.Nodes[a.params[0]]).value.data.dbl_buf())
                                    .offset(elem as isize)
                            {
                                minVal_0
                            } else {
                                *((lParse.Nodes[a.params[0]]).value.data.dbl_buf())
                                    .offset(elem as isize)
                            };
                        }
                        *((lParse.Nodes[this_node_idx]).value.undef).offset(row as isize) = 0;
                    }
                }
                *((lParse.Nodes[this_node_idx]).value.data.dbl_buf()).offset(row as isize) =
                    minVal_0;
            }
        } else if (lParse.Nodes[this_node_idx]).ntype == ValueSort::Bits {
            let mut minVal_1: c_char = 0;
            loop {
                let fresh147 = row;
                row -= 1;
                if fresh147 == 0 {
                    break;
                }
                let mut sptr1_0: *mut c_char =
                    *((lParse.Nodes[a.params[0]]).value.data.str_buf()).offset(row as isize);
                minVal_1 = b'1' as c_char;
                while *sptr1_0 != 0 {
                    if c_int::from(*sptr1_0) == '0' as i32 {
                        minVal_1 = b'0' as c_char;
                    }
                    sptr1_0 = sptr1_0.offset(1);
                }
                *(*((lParse.Nodes[this_node_idx]).value.data.str_buf()).offset(row as isize))
                    .offset(0) = minVal_1;
                *(*((lParse.Nodes[this_node_idx]).value.data.str_buf()).offset(row as isize))
                    .offset(1) = 0;
            }
        }
    }
}

/// `FuncOp::Min2`
pub(crate) fn min2_rows(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    unsafe {
        let mut i: c_int = 0;
        let mut nelem: c_long = 0;
        let mut row: c_long = lParse.nRows;
        let mut elem: c_long = lParse.nRows * (lParse.Nodes[this_node_idx]).value.nelem;

        if (lParse.Nodes[this_node_idx]).ntype == ValueSort::Long {
            loop {
                let fresh148 = row;
                row -= 1;
                if fresh148 == 0 {
                    break;
                }
                nelem = (lParse.Nodes[this_node_idx]).value.nelem;
                loop {
                    let fresh149 = nelem;
                    nelem -= 1;
                    if fresh149 == 0 {
                        break;
                    }
                    elem -= 1;
                    i = 2;
                    loop {
                        let fresh150 = i;
                        i -= 1;
                        if fresh150 == 0 {
                            break;
                        }
                        if a.vector[i as usize] > 1 {
                            a.vals[i as usize].data = NodeValue::Long(
                                *((lParse.Nodes[a.params[i as usize]]).value.data.lng_buf())
                                    .offset(elem as isize),
                            );
                            a.nulls[i as usize] =
                                *((lParse.Nodes[a.params[i as usize]]).value.undef)
                                    .offset(elem as isize);
                        } else if a.vector[i as usize] != 0 {
                            a.vals[i as usize].data = NodeValue::Long(
                                *((lParse.Nodes[a.params[i as usize]]).value.data.lng_buf())
                                    .offset(row as isize),
                            );
                            a.nulls[i as usize] =
                                *((lParse.Nodes[a.params[i as usize]]).value.undef)
                                    .offset(row as isize);
                        }
                    }
                    if c_int::from(a.nulls[0]) != 0 && c_int::from(a.nulls[1]) != 0 {
                        *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize) = 1;
                        *((lParse.Nodes[this_node_idx]).value.data.lng_buf())
                            .offset(elem as isize) = 0;
                    } else if a.nulls[0] != 0 {
                        *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize) = 0;
                        *((lParse.Nodes[this_node_idx]).value.data.lng_buf())
                            .offset(elem as isize) = a.vals[1].data.lng();
                    } else if a.nulls[1] != 0 {
                        *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize) = 0;
                        *((lParse.Nodes[this_node_idx]).value.data.lng_buf())
                            .offset(elem as isize) = a.vals[0].data.lng();
                    } else {
                        *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize) = 0;
                        *((lParse.Nodes[this_node_idx]).value.data.lng_buf())
                            .offset(elem as isize) = if a.vals[0].data.lng() < a.vals[1].data.lng()
                        {
                            a.vals[0].data.lng()
                        } else {
                            a.vals[1].data.lng()
                        };
                    }
                }
            }
        } else if (lParse.Nodes[this_node_idx]).ntype == ValueSort::Double {
            loop {
                let fresh151 = row;
                row -= 1;
                if fresh151 == 0 {
                    break;
                }
                nelem = (lParse.Nodes[this_node_idx]).value.nelem;
                loop {
                    let fresh152 = nelem;
                    nelem -= 1;
                    if fresh152 == 0 {
                        break;
                    }
                    elem -= 1;
                    i = 2;
                    loop {
                        let fresh153 = i;
                        i -= 1;
                        if fresh153 == 0 {
                            break;
                        }
                        if a.vector[i as usize] > 1 {
                            a.vals[i as usize].data = NodeValue::Double(
                                *((lParse.Nodes[a.params[i as usize]]).value.data.dbl_buf())
                                    .offset(elem as isize),
                            );
                            a.nulls[i as usize] =
                                *((lParse.Nodes[a.params[i as usize]]).value.undef)
                                    .offset(elem as isize);
                        } else if a.vector[i as usize] != 0 {
                            a.vals[i as usize].data = NodeValue::Double(
                                *((lParse.Nodes[a.params[i as usize]]).value.data.dbl_buf())
                                    .offset(row as isize),
                            );
                            a.nulls[i as usize] =
                                *((lParse.Nodes[a.params[i as usize]]).value.undef)
                                    .offset(row as isize);
                        }
                    }
                    if c_int::from(a.nulls[0]) != 0 && c_int::from(a.nulls[1]) != 0 {
                        *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize) = 1;
                        *((lParse.Nodes[this_node_idx]).value.data.dbl_buf())
                            .offset(elem as isize) = 0.0;
                    } else if a.nulls[0] != 0 {
                        *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize) = 0;
                        *((lParse.Nodes[this_node_idx]).value.data.dbl_buf())
                            .offset(elem as isize) = a.vals[1].data.dbl();
                    } else if a.nulls[1] != 0 {
                        *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize) = 0;
                        *((lParse.Nodes[this_node_idx]).value.data.dbl_buf())
                            .offset(elem as isize) = a.vals[0].data.dbl();
                    } else {
                        *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize) = 0;
                        *((lParse.Nodes[this_node_idx]).value.data.dbl_buf())
                            .offset(elem as isize) = if a.vals[0].data.dbl() < a.vals[1].data.dbl()
                        {
                            a.vals[0].data.dbl()
                        } else {
                            a.vals[1].data.dbl()
                        };
                    }
                }
            }
        }
    }
}

/// `FuncOp::Max1`
pub(crate) fn max1_rows(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    unsafe {
        let mut valInit: c_int = 0;
        let mut nelem: c_long = 0;
        let mut row: c_long = lParse.nRows;
        let mut elem: c_long = lParse.nRows * (lParse.Nodes[this_node_idx]).value.nelem;

        elem = row * (lParse.Nodes[a.params[0]]).value.nelem;
        if (lParse.Nodes[this_node_idx]).ntype == ValueSort::Long {
            let mut maxVal: c_long = 0;
            loop {
                let fresh154 = row;
                row -= 1;
                if fresh154 == 0 {
                    break;
                }
                valInit = 1;
                *((lParse.Nodes[this_node_idx]).value.undef).offset(row as isize) = 1;
                nelem = (lParse.Nodes[a.params[0]]).value.nelem;
                loop {
                    let fresh155 = nelem;
                    nelem -= 1;
                    if fresh155 == 0 {
                        break;
                    }
                    elem -= 1;
                    if *((lParse.Nodes[a.params[0]]).value.undef).offset(elem as isize) == 0 {
                        if valInit != 0 {
                            valInit = 0;
                            maxVal = *((lParse.Nodes[a.params[0]]).value.data.lng_buf())
                                .offset(elem as isize);
                        } else {
                            maxVal = if maxVal
                                > *((lParse.Nodes[a.params[0]]).value.data.lng_buf())
                                    .offset(elem as isize)
                            {
                                maxVal
                            } else {
                                *((lParse.Nodes[a.params[0]]).value.data.lng_buf())
                                    .offset(elem as isize)
                            };
                        }
                        *((lParse.Nodes[this_node_idx]).value.undef).offset(row as isize) = 0;
                    }
                }
                *((lParse.Nodes[this_node_idx]).value.data.lng_buf()).offset(row as isize) = maxVal;
            }
        } else if (lParse.Nodes[this_node_idx]).ntype == ValueSort::Double {
            let mut maxVal_0: c_double = 0.0;
            loop {
                let fresh156 = row;
                row -= 1;
                if fresh156 == 0 {
                    break;
                }
                valInit = 1;
                *((lParse.Nodes[this_node_idx]).value.undef).offset(row as isize) = 1;
                nelem = (lParse.Nodes[a.params[0]]).value.nelem;
                loop {
                    let fresh157 = nelem;
                    nelem -= 1;
                    if fresh157 == 0 {
                        break;
                    }
                    elem -= 1;
                    if *((lParse.Nodes[a.params[0]]).value.undef).offset(elem as isize) == 0 {
                        if valInit != 0 {
                            valInit = 0;
                            maxVal_0 = *((lParse.Nodes[a.params[0]]).value.data.dbl_buf())
                                .offset(elem as isize);
                        } else {
                            maxVal_0 = if maxVal_0
                                > *((lParse.Nodes[a.params[0]]).value.data.dbl_buf())
                                    .offset(elem as isize)
                            {
                                maxVal_0
                            } else {
                                *((lParse.Nodes[a.params[0]]).value.data.dbl_buf())
                                    .offset(elem as isize)
                            };
                        }
                        *((lParse.Nodes[this_node_idx]).value.undef).offset(row as isize) = 0;
                    }
                }
                *((lParse.Nodes[this_node_idx]).value.data.dbl_buf()).offset(row as isize) =
                    maxVal_0;
            }
        } else if (lParse.Nodes[this_node_idx]).ntype == ValueSort::Bits {
            let mut maxVal_1: c_char = 0;
            loop {
                let fresh158 = row;
                row -= 1;
                if fresh158 == 0 {
                    break;
                }
                let mut sptr1_1: *mut c_char =
                    *((lParse.Nodes[a.params[0]]).value.data.str_buf()).offset(row as isize);
                maxVal_1 = b'0' as c_char;
                while *sptr1_1 != 0 {
                    if c_int::from(*sptr1_1) == '1' as i32 {
                        maxVal_1 = b'1' as c_char;
                    }
                    sptr1_1 = sptr1_1.offset(1);
                }
                *(*((lParse.Nodes[this_node_idx]).value.data.str_buf()).offset(row as isize))
                    .offset(0) = maxVal_1;
                *(*((lParse.Nodes[this_node_idx]).value.data.str_buf()).offset(row as isize))
                    .offset(1) = 0;
            }
        }
    }
}

/// `FuncOp::Max2`
pub(crate) fn max2_rows(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    unsafe {
        let mut i: c_int = 0;
        let mut nelem: c_long = 0;
        let mut row: c_long = lParse.nRows;
        let mut elem: c_long = lParse.nRows * (lParse.Nodes[this_node_idx]).value.nelem;

        if (lParse.Nodes[this_node_idx]).ntype == ValueSort::Long {
            loop {
                let fresh159 = row;
                row -= 1;
                if fresh159 == 0 {
                    break;
                }
                nelem = (lParse.Nodes[this_node_idx]).value.nelem;
                loop {
                    let fresh160 = nelem;
                    nelem -= 1;
                    if fresh160 == 0 {
                        break;
                    }
                    elem -= 1;
                    i = 2;
                    loop {
                        let fresh161 = i;
                        i -= 1;
                        if fresh161 == 0 {
                            break;
                        }
                        if a.vector[i as usize] > 1 {
                            a.vals[i as usize].data = NodeValue::Long(
                                *((lParse.Nodes[a.params[i as usize]]).value.data.lng_buf())
                                    .offset(elem as isize),
                            );
                            a.nulls[i as usize] =
                                *((lParse.Nodes[a.params[i as usize]]).value.undef)
                                    .offset(elem as isize);
                        } else if a.vector[i as usize] != 0 {
                            a.vals[i as usize].data = NodeValue::Long(
                                *((lParse.Nodes[a.params[i as usize]]).value.data.lng_buf())
                                    .offset(row as isize),
                            );
                            a.nulls[i as usize] =
                                *((lParse.Nodes[a.params[i as usize]]).value.undef)
                                    .offset(row as isize);
                        }
                    }
                    if c_int::from(a.nulls[0]) != 0 && c_int::from(a.nulls[1]) != 0 {
                        *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize) = 1;
                        *((lParse.Nodes[this_node_idx]).value.data.lng_buf())
                            .offset(elem as isize) = 0;
                    } else if a.nulls[0] != 0 {
                        *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize) = 0;
                        *((lParse.Nodes[this_node_idx]).value.data.lng_buf())
                            .offset(elem as isize) = a.vals[1].data.lng();
                    } else if a.nulls[1] != 0 {
                        *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize) = 0;
                        *((lParse.Nodes[this_node_idx]).value.data.lng_buf())
                            .offset(elem as isize) = a.vals[0].data.lng();
                    } else {
                        *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize) = 0;
                        *((lParse.Nodes[this_node_idx]).value.data.lng_buf())
                            .offset(elem as isize) = if a.vals[0].data.lng() > a.vals[1].data.lng()
                        {
                            a.vals[0].data.lng()
                        } else {
                            a.vals[1].data.lng()
                        };
                    }
                }
            }
        } else if (lParse.Nodes[this_node_idx]).ntype == ValueSort::Double {
            loop {
                let fresh162 = row;
                row -= 1;
                if fresh162 == 0 {
                    break;
                }
                nelem = (lParse.Nodes[this_node_idx]).value.nelem;
                loop {
                    let fresh163 = nelem;
                    nelem -= 1;
                    if fresh163 == 0 {
                        break;
                    }
                    elem -= 1;
                    i = 2;
                    loop {
                        let fresh164 = i;
                        i -= 1;
                        if fresh164 == 0 {
                            break;
                        }
                        if a.vector[i as usize] > 1 {
                            a.vals[i as usize].data = NodeValue::Double(
                                *((lParse.Nodes[a.params[i as usize]]).value.data.dbl_buf())
                                    .offset(elem as isize),
                            );
                            a.nulls[i as usize] =
                                *((lParse.Nodes[a.params[i as usize]]).value.undef)
                                    .offset(elem as isize);
                        } else if a.vector[i as usize] != 0 {
                            a.vals[i as usize].data = NodeValue::Double(
                                *((lParse.Nodes[a.params[i as usize]]).value.data.dbl_buf())
                                    .offset(row as isize),
                            );
                            a.nulls[i as usize] =
                                *((lParse.Nodes[a.params[i as usize]]).value.undef)
                                    .offset(row as isize);
                        }
                    }
                    if c_int::from(a.nulls[0]) != 0 && c_int::from(a.nulls[1]) != 0 {
                        *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize) = 1;
                        *((lParse.Nodes[this_node_idx]).value.data.dbl_buf())
                            .offset(elem as isize) = 0.0;
                    } else if a.nulls[0] != 0 {
                        *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize) = 0;
                        *((lParse.Nodes[this_node_idx]).value.data.dbl_buf())
                            .offset(elem as isize) = a.vals[1].data.dbl();
                    } else if a.nulls[1] != 0 {
                        *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize) = 0;
                        *((lParse.Nodes[this_node_idx]).value.data.dbl_buf())
                            .offset(elem as isize) = a.vals[0].data.dbl();
                    } else {
                        *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize) = 0;
                        *((lParse.Nodes[this_node_idx]).value.data.dbl_buf())
                            .offset(elem as isize) = if a.vals[0].data.dbl() > a.vals[1].data.dbl()
                        {
                            a.vals[0].data.dbl()
                        } else {
                            a.vals[1].data.dbl()
                        };
                    }
                }
            }
        }
    }
}

/// `FuncOp::Near`
pub(crate) fn near_rows(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    unsafe {
        let mut i: c_int = 0;
        let mut nelem: c_long = 0;
        let mut row: c_long = lParse.nRows;
        let mut elem: c_long = lParse.nRows * (lParse.Nodes[this_node_idx]).value.nelem;

        loop {
            let fresh165 = row;
            row -= 1;
            if fresh165 == 0 {
                break;
            }
            nelem = (lParse.Nodes[this_node_idx]).value.nelem;
            loop {
                let fresh166 = nelem;
                nelem -= 1;
                if fresh166 == 0 {
                    break;
                }
                elem -= 1;
                i = 3 as c_int;
                loop {
                    let fresh167 = i;
                    i -= 1;
                    if fresh167 == 0 {
                        break;
                    }
                    if a.vector[i as usize] > 1 {
                        a.vals[i as usize].data = NodeValue::Double(
                            *((lParse.Nodes[a.params[i as usize]]).value.data.dbl_buf())
                                .offset(elem as isize),
                        );
                        a.nulls[i as usize] = *((lParse.Nodes[a.params[i as usize]]).value.undef)
                            .offset(elem as isize);
                    } else if a.vector[i as usize] != 0 {
                        a.vals[i as usize].data = NodeValue::Double(
                            *((lParse.Nodes[a.params[i as usize]]).value.data.dbl_buf())
                                .offset(row as isize),
                        );
                        a.nulls[i as usize] = *((lParse.Nodes[a.params[i as usize]]).value.undef)
                            .offset(row as isize);
                    }
                }
                let fresh168 =
                    &mut *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize);
                *fresh168 = c_int::from(
                    c_int::from(a.nulls[0]) != 0
                        || c_int::from(a.nulls[1]) != 0
                        || c_int::from(a.nulls[2]) != 0,
                ) as c_char;
                if *fresh168 == 0 {
                    *((lParse.Nodes[this_node_idx]).value.data.log_buf()).offset(elem as isize) =
                        bnear(
                            a.vals[0].data.dbl(),
                            a.vals[1].data.dbl(),
                            a.vals[2].data.dbl(),
                        );
                }
            }
        }
    }
}

/// `FuncOp::Circle`
pub(crate) fn circle_rows(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    unsafe {
        let mut i: c_int = 0;
        let mut nelem: c_long = 0;
        let mut row: c_long = lParse.nRows;
        let mut elem: c_long = lParse.nRows * (lParse.Nodes[this_node_idx]).value.nelem;

        loop {
            let fresh169 = row;
            row -= 1;
            if fresh169 == 0 {
                break;
            }
            nelem = (lParse.Nodes[this_node_idx]).value.nelem;
            loop {
                let fresh170 = nelem;
                nelem -= 1;
                if fresh170 == 0 {
                    break;
                }
                elem -= 1;
                i = 5 as c_int;
                loop {
                    let fresh171 = i;
                    i -= 1;
                    if fresh171 == 0 {
                        break;
                    }
                    if a.vector[i as usize] > 1 {
                        a.vals[i as usize].data = NodeValue::Double(
                            *((lParse.Nodes[a.params[i as usize]]).value.data.dbl_buf())
                                .offset(elem as isize),
                        );
                        a.nulls[i as usize] = *((lParse.Nodes[a.params[i as usize]]).value.undef)
                            .offset(elem as isize);
                    } else if a.vector[i as usize] != 0 {
                        a.vals[i as usize].data = NodeValue::Double(
                            *((lParse.Nodes[a.params[i as usize]]).value.data.dbl_buf())
                                .offset(row as isize),
                        );
                        a.nulls[i as usize] = *((lParse.Nodes[a.params[i as usize]]).value.undef)
                            .offset(row as isize);
                    }
                }
                let fresh172 =
                    &mut *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize);
                *fresh172 = c_int::from(
                    c_int::from(a.nulls[0]) != 0
                        || c_int::from(a.nulls[1]) != 0
                        || c_int::from(a.nulls[2]) != 0
                        || c_int::from(a.nulls[3]) != 0
                        || c_int::from(a.nulls[4]) != 0,
                ) as c_char;
                if *fresh172 == 0 {
                    *((lParse.Nodes[this_node_idx]).value.data.log_buf()).offset(elem as isize) =
                        circle(
                            a.vals[0].data.dbl(),
                            a.vals[1].data.dbl(),
                            a.vals[2].data.dbl(),
                            a.vals[3].data.dbl(),
                            a.vals[4].data.dbl(),
                        );
                }
            }
        }
    }
}

/// `FuncOp::Box`
pub(crate) fn box_rows(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    unsafe {
        let mut i: c_int = 0;
        let mut nelem: c_long = 0;
        let mut row: c_long = lParse.nRows;
        let mut elem: c_long = lParse.nRows * (lParse.Nodes[this_node_idx]).value.nelem;

        loop {
            let fresh173 = row;
            row -= 1;
            if fresh173 == 0 {
                break;
            }
            nelem = (lParse.Nodes[this_node_idx]).value.nelem;
            loop {
                let fresh174 = nelem;
                nelem -= 1;
                if fresh174 == 0 {
                    break;
                }
                elem -= 1;
                i = 7 as c_int;
                loop {
                    let fresh175 = i;
                    i -= 1;
                    if fresh175 == 0 {
                        break;
                    }
                    if a.vector[i as usize] > 1 {
                        a.vals[i as usize].data = NodeValue::Double(
                            *((lParse.Nodes[a.params[i as usize]]).value.data.dbl_buf())
                                .offset(elem as isize),
                        );
                        a.nulls[i as usize] = *((lParse.Nodes[a.params[i as usize]]).value.undef)
                            .offset(elem as isize);
                    } else if a.vector[i as usize] != 0 {
                        a.vals[i as usize].data = NodeValue::Double(
                            *((lParse.Nodes[a.params[i as usize]]).value.data.dbl_buf())
                                .offset(row as isize),
                        );
                        a.nulls[i as usize] = *((lParse.Nodes[a.params[i as usize]]).value.undef)
                            .offset(row as isize);
                    }
                }
                let fresh176 =
                    &mut *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize);
                *fresh176 = c_int::from(
                    c_int::from(a.nulls[0]) != 0
                        || c_int::from(a.nulls[1]) != 0
                        || c_int::from(a.nulls[2]) != 0
                        || c_int::from(a.nulls[3]) != 0
                        || c_int::from(a.nulls[4]) != 0
                        || c_int::from(a.nulls[5]) != 0
                        || c_int::from(a.nulls[6]) != 0,
                ) as c_char;
                if *fresh176 == 0 {
                    *((lParse.Nodes[this_node_idx]).value.data.log_buf()).offset(elem as isize) =
                        saobox(
                            a.vals[0].data.dbl(),
                            a.vals[1].data.dbl(),
                            a.vals[2].data.dbl(),
                            a.vals[3].data.dbl(),
                            a.vals[4].data.dbl(),
                            a.vals[5].data.dbl(),
                            a.vals[6].data.dbl(),
                        );
                }
            }
        }
    }
}

/// `FuncOp::Ellipse`
pub(crate) fn ellipse_rows(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    unsafe {
        let mut i: c_int = 0;
        let mut nelem: c_long = 0;
        let mut row: c_long = lParse.nRows;
        let mut elem: c_long = lParse.nRows * (lParse.Nodes[this_node_idx]).value.nelem;

        loop {
            let fresh177 = row;
            row -= 1;
            if fresh177 == 0 {
                break;
            }
            nelem = (lParse.Nodes[this_node_idx]).value.nelem;
            loop {
                let fresh178 = nelem;
                nelem -= 1;
                if fresh178 == 0 {
                    break;
                }
                elem -= 1;
                i = 7 as c_int;
                loop {
                    let fresh179 = i;
                    i -= 1;
                    if fresh179 == 0 {
                        break;
                    }
                    if a.vector[i as usize] > 1 {
                        a.vals[i as usize].data = NodeValue::Double(
                            *((lParse.Nodes[a.params[i as usize]]).value.data.dbl_buf())
                                .offset(elem as isize),
                        );
                        a.nulls[i as usize] = *((lParse.Nodes[a.params[i as usize]]).value.undef)
                            .offset(elem as isize);
                    } else if a.vector[i as usize] != 0 {
                        a.vals[i as usize].data = NodeValue::Double(
                            *((lParse.Nodes[a.params[i as usize]]).value.data.dbl_buf())
                                .offset(row as isize),
                        );
                        a.nulls[i as usize] = *((lParse.Nodes[a.params[i as usize]]).value.undef)
                            .offset(row as isize);
                    }
                }
                let fresh180 =
                    &mut *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize);
                *fresh180 = c_int::from(
                    c_int::from(a.nulls[0]) != 0
                        || c_int::from(a.nulls[1]) != 0
                        || c_int::from(a.nulls[2]) != 0
                        || c_int::from(a.nulls[3]) != 0
                        || c_int::from(a.nulls[4]) != 0
                        || c_int::from(a.nulls[5]) != 0
                        || c_int::from(a.nulls[6]) != 0,
                ) as c_char;
                if *fresh180 == 0 {
                    *((lParse.Nodes[this_node_idx]).value.data.log_buf()).offset(elem as isize) =
                        ellipse(
                            a.vals[0].data.dbl(),
                            a.vals[1].data.dbl(),
                            a.vals[2].data.dbl(),
                            a.vals[3].data.dbl(),
                            a.vals[4].data.dbl(),
                            a.vals[5].data.dbl(),
                            a.vals[6].data.dbl(),
                        );
                }
            }
        }
    }
}

/// `FuncOp::IfThenElse`
pub(crate) fn if_then_else_rows(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    unsafe {
        let mut i: c_int = 0;
        let mut nelem: c_long = 0;
        let mut row: c_long = lParse.nRows;
        let mut elem: c_long = lParse.nRows * (lParse.Nodes[this_node_idx]).value.nelem;

        match (lParse.Nodes[this_node_idx]).ntype {
            ValueSort::Boolean => loop {
                let fresh181 = row;
                row -= 1;
                if fresh181 == 0 {
                    break;
                }
                nelem = (lParse.Nodes[this_node_idx]).value.nelem;
                loop {
                    let fresh182 = nelem;
                    nelem -= 1;
                    if fresh182 == 0 {
                        break;
                    }
                    elem -= 1;
                    if a.vector[2] > 1 {
                        a.vals[2].data = NodeValue::Logical(
                            *((lParse.Nodes[a.params[2]]).value.data.log_buf())
                                .offset(elem as isize),
                        );
                        a.nulls[2] =
                            *((lParse.Nodes[a.params[2]]).value.undef).offset(elem as isize);
                    } else if a.vector[2] != 0 {
                        a.vals[2].data = NodeValue::Logical(
                            *((lParse.Nodes[a.params[2]]).value.data.log_buf())
                                .offset(row as isize),
                        );
                        a.nulls[2] =
                            *((lParse.Nodes[a.params[2]]).value.undef).offset(row as isize);
                    }
                    i = 2;
                    loop {
                        let fresh183 = i;
                        i -= 1;
                        if fresh183 == 0 {
                            break;
                        }
                        if a.vector[i as usize] > 1 {
                            a.vals[i as usize].data = NodeValue::Logical(
                                *((lParse.Nodes[a.params[i as usize]]).value.data.log_buf())
                                    .offset(elem as isize),
                            );
                            a.nulls[i as usize] =
                                *((lParse.Nodes[a.params[i as usize]]).value.undef)
                                    .offset(elem as isize);
                        } else if a.vector[i as usize] != 0 {
                            a.vals[i as usize].data = NodeValue::Logical(
                                *((lParse.Nodes[a.params[i as usize]]).value.data.log_buf())
                                    .offset(row as isize),
                            );
                            a.nulls[i as usize] =
                                *((lParse.Nodes[a.params[i as usize]]).value.undef)
                                    .offset(row as isize);
                        }
                    }
                    let fresh184 =
                        &mut *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize);
                    *fresh184 = a.nulls[2];
                    if *fresh184 == 0 {
                        if a.vals[2].data.log() != 0 {
                            *((lParse.Nodes[this_node_idx]).value.data.log_buf())
                                .offset(elem as isize) = a.vals[0].data.log();
                            *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize) =
                                a.nulls[0];
                        } else {
                            *((lParse.Nodes[this_node_idx]).value.data.log_buf())
                                .offset(elem as isize) = a.vals[1].data.log();
                            *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize) =
                                a.nulls[1];
                        }
                    }
                }
            },
            ValueSort::Long => loop {
                let fresh185 = row;
                row -= 1;
                if fresh185 == 0 {
                    break;
                }
                nelem = (lParse.Nodes[this_node_idx]).value.nelem;
                loop {
                    let fresh186 = nelem;
                    nelem -= 1;
                    if fresh186 == 0 {
                        break;
                    }
                    elem -= 1;
                    if a.vector[2] > 1 {
                        a.vals[2].data = NodeValue::Logical(
                            *((lParse.Nodes[a.params[2]]).value.data.log_buf())
                                .offset(elem as isize),
                        );
                        a.nulls[2] =
                            *((lParse.Nodes[a.params[2]]).value.undef).offset(elem as isize);
                    } else if a.vector[2] != 0 {
                        a.vals[2].data = NodeValue::Logical(
                            *((lParse.Nodes[a.params[2]]).value.data.log_buf())
                                .offset(row as isize),
                        );
                        a.nulls[2] =
                            *((lParse.Nodes[a.params[2]]).value.undef).offset(row as isize);
                    }
                    i = 2;
                    loop {
                        let fresh187 = i;
                        i -= 1;
                        if fresh187 == 0 {
                            break;
                        }
                        if a.vector[i as usize] > 1 {
                            a.vals[i as usize].data = NodeValue::Long(
                                *((lParse.Nodes[a.params[i as usize]]).value.data.lng_buf())
                                    .offset(elem as isize),
                            );
                            a.nulls[i as usize] =
                                *((lParse.Nodes[a.params[i as usize]]).value.undef)
                                    .offset(elem as isize);
                        } else if a.vector[i as usize] != 0 {
                            a.vals[i as usize].data = NodeValue::Long(
                                *((lParse.Nodes[a.params[i as usize]]).value.data.lng_buf())
                                    .offset(row as isize),
                            );
                            a.nulls[i as usize] =
                                *((lParse.Nodes[a.params[i as usize]]).value.undef)
                                    .offset(row as isize);
                        }
                    }
                    let fresh188 =
                        &mut *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize);
                    *fresh188 = a.nulls[2];
                    if *fresh188 == 0 {
                        if a.vals[2].data.log() != 0 {
                            *((lParse.Nodes[this_node_idx]).value.data.lng_buf())
                                .offset(elem as isize) = a.vals[0].data.lng();
                            *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize) =
                                a.nulls[0];
                        } else {
                            *((lParse.Nodes[this_node_idx]).value.data.lng_buf())
                                .offset(elem as isize) = a.vals[1].data.lng();
                            *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize) =
                                a.nulls[1];
                        }
                    }
                }
            },
            ValueSort::Double => loop {
                let fresh189 = row;
                row -= 1;
                if fresh189 == 0 {
                    break;
                }
                nelem = (lParse.Nodes[this_node_idx]).value.nelem;
                loop {
                    let fresh190 = nelem;
                    nelem -= 1;
                    if fresh190 == 0 {
                        break;
                    }
                    elem -= 1;
                    if a.vector[2] > 1 {
                        a.vals[2].data = NodeValue::Logical(
                            *((lParse.Nodes[a.params[2]]).value.data.log_buf())
                                .offset(elem as isize),
                        );
                        a.nulls[2] =
                            *((lParse.Nodes[a.params[2]]).value.undef).offset(elem as isize);
                    } else if a.vector[2] != 0 {
                        a.vals[2].data = NodeValue::Logical(
                            *((lParse.Nodes[a.params[2]]).value.data.log_buf())
                                .offset(row as isize),
                        );
                        a.nulls[2] =
                            *((lParse.Nodes[a.params[2]]).value.undef).offset(row as isize);
                    }
                    i = 2;
                    loop {
                        let fresh191 = i;
                        i -= 1;
                        if fresh191 == 0 {
                            break;
                        }
                        if a.vector[i as usize] > 1 {
                            a.vals[i as usize].data = NodeValue::Double(
                                *((lParse.Nodes[a.params[i as usize]]).value.data.dbl_buf())
                                    .offset(elem as isize),
                            );
                            a.nulls[i as usize] =
                                *((lParse.Nodes[a.params[i as usize]]).value.undef)
                                    .offset(elem as isize);
                        } else if a.vector[i as usize] != 0 {
                            a.vals[i as usize].data = NodeValue::Double(
                                *((lParse.Nodes[a.params[i as usize]]).value.data.dbl_buf())
                                    .offset(row as isize),
                            );
                            a.nulls[i as usize] =
                                *((lParse.Nodes[a.params[i as usize]]).value.undef)
                                    .offset(row as isize);
                        }
                    }
                    let fresh192 =
                        &mut *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize);
                    *fresh192 = a.nulls[2];
                    if *fresh192 == 0 {
                        if a.vals[2].data.log() != 0 {
                            *((lParse.Nodes[this_node_idx]).value.data.dbl_buf())
                                .offset(elem as isize) = a.vals[0].data.dbl();
                            *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize) =
                                a.nulls[0];
                        } else {
                            *((lParse.Nodes[this_node_idx]).value.data.dbl_buf())
                                .offset(elem as isize) = a.vals[1].data.dbl();
                            *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize) =
                                a.nulls[1];
                        }
                    }
                }
            },
            ValueSort::String => loop {
                let fresh193 = row;
                row -= 1;
                if fresh193 == 0 {
                    break;
                }
                if a.vector[2] != 0 {
                    a.vals[2].data = NodeValue::Logical(
                        *((lParse.Nodes[a.params[2]]).value.data.log_buf()).offset(row as isize),
                    );
                    a.nulls[2] = *((lParse.Nodes[a.params[2]]).value.undef).offset(row as isize);
                }
                i = 2;
                loop {
                    let fresh194 = i;
                    i -= 1;
                    if fresh194 == 0 {
                        break;
                    }
                    if a.vector[i as usize] != 0 {
                        strcpy(
                            a.vals[i as usize].data.text_mut_ptr(),
                            *((lParse.Nodes[a.params[i as usize]]).value.data.str_buf())
                                .offset(row as isize),
                        );
                        a.nulls[i as usize] = *((lParse.Nodes[a.params[i as usize]]).value.undef)
                            .offset(row as isize);
                    }
                }
                let fresh195 =
                    &mut *((lParse.Nodes[this_node_idx]).value.undef).offset(row as isize);
                *fresh195 = a.nulls[2];
                if *fresh195 == 0 {
                    if a.vals[2].data.log() != 0 {
                        strcpy(
                            *((lParse.Nodes[this_node_idx]).value.data.str_buf())
                                .offset(row as isize),
                            a.vals[0].data.text_mut_ptr(),
                        );
                        *((lParse.Nodes[this_node_idx]).value.undef).offset(row as isize) =
                            a.nulls[0];
                    } else {
                        strcpy(
                            *((lParse.Nodes[this_node_idx]).value.data.str_buf())
                                .offset(row as isize),
                            a.vals[1].data.text_mut_ptr(),
                        );
                        *((lParse.Nodes[this_node_idx]).value.undef).offset(row as isize) =
                            a.nulls[1];
                    }
                } else {
                    *(*((lParse.Nodes[this_node_idx]).value.data.str_buf()).offset(row as isize))
                        .offset(0) = 0;
                }
            },
            _ => {}
        }
    }
}

/// `FuncOp::StrMid`
pub(crate) fn str_mid_rows(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    unsafe {
        let nelem: c_long = 0;
        let mut row: c_long = lParse.nRows;

        let strconst: c_int =
            c_int::from((lParse.Nodes[a.params[0]]).operation == Operation::Const);
        let posconst: c_int =
            c_int::from((lParse.Nodes[a.params[1]]).operation == Operation::Const);
        let lenconst: c_int =
            c_int::from((lParse.Nodes[a.params[2]]).operation == Operation::Const);
        let dest_len: c_int = (lParse.Nodes[this_node_idx]).value.nelem as c_int;
        let mut src_len: c_int = (lParse.Nodes[a.params[0]]).value.nelem as c_int;
        loop {
            let fresh196 = row;
            row -= 1;
            if fresh196 == 0 {
                break;
            }
            let mut pos: c_int = 0;
            let mut len: c_int = 0;
            let mut str: *mut c_char = ptr::null_mut();
            let mut undef: c_int = 0;
            if posconst != 0 {
                pos = (lParse.Nodes[a.params[1]]).value.data.lng() as c_int;
            } else {
                pos = *((lParse.Nodes[a.params[1]]).value.data.lng_buf()).offset(row as isize)
                    as c_int;
                if *((lParse.Nodes[a.params[1]]).value.undef).offset(row as isize) != 0 {
                    undef = 1;
                }
            }
            if strconst != 0 {
                str = (lParse.Nodes[a.params[0]]).value.data.text_mut_ptr();
                if src_len == 0 {
                    src_len = strlen(str) as c_int;
                }
            } else {
                str = *((lParse.Nodes[a.params[0]]).value.data.str_buf()).offset(row as isize);
                if *((lParse.Nodes[a.params[0]]).value.undef).offset(row as isize) != 0 {
                    undef = 1;
                }
            }
            if lenconst != 0 {
                len = dest_len;
            } else {
                len = *((lParse.Nodes[a.params[2]]).value.data.lng_buf()).offset(row as isize)
                    as c_int;
                if *((lParse.Nodes[a.params[2]]).value.undef).offset(row as isize) != 0 {
                    undef = 1;
                }
            }
            *(*((lParse.Nodes[this_node_idx]).value.data.str_buf()).offset(row as isize))
                .offset(0) = 0;
            if pos == 0 {
                undef = 1;
            }
            if undef == 0
                && cstrmid(
                    lParse,
                    *((lParse.Nodes[this_node_idx]).value.data.str_buf()).offset(row as isize),
                    len,
                    str,
                    src_len,
                    pos,
                ) < 0
            {
                break;
            }
            *((lParse.Nodes[this_node_idx]).value.undef).offset(row as isize) = undef as c_char;
        }
    }
}

/// `FuncOp::StrPos`
pub(crate) fn str_pos_rows(lParse: &mut ParseData, this_node_idx: usize, a: &mut FuncArgs) {
    unsafe {
        let mut row: c_long = lParse.nRows;

        let const1: c_int = c_int::from((lParse.Nodes[a.params[0]]).operation == Operation::Const);
        let const2: c_int = c_int::from((lParse.Nodes[a.params[1]]).operation == Operation::Const);
        loop {
            let fresh197 = row;
            row -= 1;
            if fresh197 == 0 {
                break;
            }
            let mut str1: *mut c_char = ptr::null_mut();
            let mut str2: *mut c_char = ptr::null_mut();
            let mut undef_0: c_int = 0;
            if const1 != 0 {
                str1 = (lParse.Nodes[a.params[0]]).value.data.text_mut_ptr();
            } else {
                str1 = *((lParse.Nodes[a.params[0]]).value.data.str_buf()).offset(row as isize);
                if *((lParse.Nodes[a.params[0]]).value.undef).offset(row as isize) != 0 {
                    undef_0 = 1;
                }
            }
            if const2 != 0 {
                str2 = (lParse.Nodes[a.params[1]]).value.data.text_mut_ptr();
            } else {
                str2 = *((lParse.Nodes[a.params[1]]).value.data.str_buf()).offset(row as isize);
                if *((lParse.Nodes[a.params[1]]).value.undef).offset(row as isize) != 0 {
                    undef_0 = 1;
                }
            }
            *((lParse.Nodes[this_node_idx]).value.data.lng_buf()).offset(row as isize) = 0;
            if undef_0 == 0 {
                let res_0: *mut c_char = strstr(str1, str2);
                if res_0.is_null() {
                    undef_0 = 1;
                    *((lParse.Nodes[this_node_idx]).value.data.lng_buf()).offset(row as isize) = 0;
                } else {
                    *((lParse.Nodes[this_node_idx]).value.data.lng_buf()).offset(row as isize) =
                        res_0.offset_from(str1) as c_long + 1;
                }
            }
            *((lParse.Nodes[this_node_idx]).value.undef).offset(row as isize) = undef_0 as c_char;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Build a node holding a constant, and return its index.
    fn push_const(p: &mut ParseData, sort: ValueSort, v: NodeValue) -> usize {
        let n = New_Const(p, sort, v);
        assert!(n >= 0);
        n as usize
    }

    /// Run a `_const` kernel over the given constant arguments and hand back
    /// the value it folded.
    ///
    /// This is the shape the split bought: a kernel is now a plain function
    /// over a node and its arguments, so it can be exercised without going
    /// anywhere near a FITS file, a parser or an iterator.
    fn fold(kernel: Kernel, ret: ValueSort, args: &[(ValueSort, NodeValue)]) -> NodeValue {
        let mut p = ParseData::default();
        let mut a = FuncArgs::default();
        for (i, &(sort, v)) in args.iter().enumerate() {
            a.params[i] = push_const(&mut p, sort, v);
            a.vals[i].data = v;
        }
        let this = push_const(&mut p, ret, NodeValue::Empty);
        kernel(&mut p, this, &mut a);
        p.Nodes[this].value.data
    }

    fn dbl(v: f64) -> (ValueSort, NodeValue) {
        (ValueSort::Double, NodeValue::Double(v))
    }
    fn lng(v: c_long) -> (ValueSort, NodeValue) {
        (ValueSort::Long, NodeValue::Long(v))
    }

    fn as_dbl(v: NodeValue) -> f64 {
        v.dbl()
    }

    #[test]
    fn constant_folding_of_the_transcendental_kernels() {
        let cases: &[(Kernel, f64, f64)] = &[
            (sin_const, 0.0, 0.0),
            (cos_const, 0.0, 1.0),
            (tan_const, 0.0, 0.0),
            (asin_const, 1.0, core::f64::consts::FRAC_PI_2),
            (acos_const, 1.0, 0.0),
            (atan_const, 0.0, 0.0),
            (sinh_const, 0.0, 0.0),
            (cosh_const, 0.0, 1.0),
            (tanh_const, 0.0, 0.0),
            (exp_const, 0.0, 1.0),
            (log_const, 1.0, 0.0),
            (log10_const, 100.0, 2.0),
            (sqrt_const, 9.0, 3.0),
            (ceil_const, 2.1, 3.0),
            (floor_const, 2.9, 2.0),
            (round_const, 2.5, 3.0),
        ];
        for &(kernel, arg, want) in cases {
            let got = as_dbl(fold(kernel, ValueSort::Double, &[dbl(arg)]));
            assert!(
                (got - want).abs() < 1e-12,
                "kernel({arg}) = {got}, want {want}"
            );
        }
    }

    #[test]
    fn abs_kernel_covers_both_numeric_sorts() {
        assert_eq!(
            fold(abs_const, ValueSort::Long, &[lng(-7)]).lng(),
            7,
            "abs of a long"
        );
        assert!(
            (as_dbl(fold(abs_const, ValueSort::Double, &[dbl(-2.5)])) - 2.5).abs() < 1e-12,
            "abs of a double"
        );
    }

    #[test]
    fn min_and_max_of_two_arguments() {
        assert_eq!(
            fold(min2_const, ValueSort::Long, &[lng(3), lng(7)]).lng(),
            3
        );
        assert_eq!(
            fold(max2_const, ValueSort::Long, &[lng(3), lng(7)]).lng(),
            7
        );
        assert_eq!(
            fold(min2_const, ValueSort::Long, &[lng(7), lng(3)]).lng(),
            3
        );
        assert!(
            (as_dbl(fold(min2_const, ValueSort::Double, &[dbl(1.5), dbl(-0.5)])) + 0.5).abs()
                < 1e-12
        );
    }

    #[test]
    fn atan2_and_angsep() {
        let q = as_dbl(fold(atan2_const, ValueSort::Double, &[dbl(1.0), dbl(1.0)]));
        assert!((q - core::f64::consts::FRAC_PI_4).abs() < 1e-12, "got {q}");

        /* a quarter turn along a meridian is 90 degrees */
        let sep = as_dbl(fold(
            ang_sep_const,
            ValueSort::Double,
            &[dbl(0.0), dbl(0.0), dbl(0.0), dbl(90.0)],
        ));
        assert!((sep - 90.0).abs() < 1e-9, "got {sep}");

        /* and a point is zero degrees from itself */
        let zero = as_dbl(fold(
            ang_sep_const,
            ValueSort::Double,
            &[dbl(10.0), dbl(20.0), dbl(10.0), dbl(20.0)],
        ));
        assert!(zero.abs() < 1e-9, "got {zero}");
    }

    #[test]
    fn null_predicates_on_constants() {
        /* a constant is never null */
        assert_eq!(fold(is_null_const, ValueSort::Boolean, &[lng(1)]).log(), 0);
        /* and NVALID counts its one defined element */
        assert_eq!(fold(non_null_const, ValueSort::Long, &[lng(1)]).lng(), 1);
    }
}
