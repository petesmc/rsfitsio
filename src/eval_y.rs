#![allow(
    dead_code,
    mutable_transmutes,
    non_camel_case_types,
    non_snake_case,
    non_upper_case_globals,
    unused_assignments,
    unused_mut,
    deprecated
)]

use std::ffi::CStr;
use std::ptr;

use bytemuck::cast_slice;
use libc::{calloc, free, malloc, memcpy, memset, realloc, sprintf, time, time_t};

// Math function wrappers
fn sqrt(x: f64) -> f64 {
    x.sqrt()
}
fn sin(x: f64) -> f64 {
    x.sin()
}
fn cos(x: f64) -> f64 {
    x.cos()
}
fn ceil(x: f64) -> f64 {
    x.ceil()
}
fn floor(x: f64) -> f64 {
    x.floor()
}
fn fabs(x: f64) -> f64 {
    x.abs()
}
fn acos(x: f64) -> f64 {
    x.acos()
}
fn asin(x: f64) -> f64 {
    x.asin()
}
fn atan(x: f64) -> f64 {
    x.atan()
}
fn atan2(y: f64, x: f64) -> f64 {
    y.atan2(x)
}
fn cosh(x: f64) -> f64 {
    x.cosh()
}
fn exp(x: f64) -> f64 {
    x.exp()
}
fn log(x: f64) -> f64 {
    x.ln()
}
fn log10(x: f64) -> f64 {
    x.log10()
}
fn pow(x: f64, y: f64) -> f64 {
    x.powf(y)
}
fn sinh(x: f64) -> f64 {
    x.sinh()
}
fn tan(x: f64) -> f64 {
    x.tan()
}
fn tanh(x: f64) -> f64 {
    x.tanh()
}

use crate::atoi;
use crate::c_types::{
    c_char, c_double, c_int, c_long, c_schar, c_short, c_uchar, c_uint, c_ulong, c_void,
};
use crate::cfileio::{ffclos, ffexts, ffopen};
use crate::eval_defs::{Node, ParseData, data_union, lval};
use crate::eval_l::{fits_parser_yyGetVariable, fits_parser_yylex};
use crate::eval_tab::{FITS_PARSER_YYSTYPE, fits_parser_yytokentype};
use crate::fitscore::ffpmsg;
use crate::fitscore::{ffgcno, ffghdn, ffmahd, ffmnhd, ffupch};
use crate::fitsio::{LONGLONG, fitsfile};
use crate::getcold::ffgcvd;
use crate::getkey::{ffgkyd, ffgkyj, ffgkys};
use crate::region::{SAORegion, WCSdata, fits_in_region, fits_read_rgnfile};
use crate::simplerng::{
    simplerng_getnorm, simplerng_getpoisson, simplerng_getuniform, simplerng_srand,
};
use crate::wcssub::ffgtcs;
use crate::wrappers::strncpy;
use crate::wrappers::{strcat, strcmp, strcpy, strlen, strstr};

pub type yyscan_t = *mut c_void;

pub type yy_state_t = yytype_int16;
pub type yytype_int16 = c_short;
pub type yysymbol_kind_t = c_int;
pub const YYSYMBOL_sexpr: yysymbol_kind_t = 65;
pub const YYSYMBOL_bits: yysymbol_kind_t = 64;
pub const YYSYMBOL_bexpr: yysymbol_kind_t = 63;
pub const YYSYMBOL_expr: yysymbol_kind_t = 62;
pub const YYSYMBOL_vector: yysymbol_kind_t = 61;
pub const YYSYMBOL_bvector: yysymbol_kind_t = 60;
pub const YYSYMBOL_line: yysymbol_kind_t = 59;
pub const YYSYMBOL_lines: yysymbol_kind_t = 58;
pub const YYSYMBOL_YYACCEPT: yysymbol_kind_t = 57;
pub const YYSYMBOL_56_: yysymbol_kind_t = 56;
pub const YYSYMBOL_55_: yysymbol_kind_t = 55;
pub const YYSYMBOL_54_: yysymbol_kind_t = 54;
pub const YYSYMBOL_53_n_: yysymbol_kind_t = 53;
pub const YYSYMBOL_DIFF: yysymbol_kind_t = 52;
pub const YYSYMBOL_ACCUM: yysymbol_kind_t = 51;
pub const YYSYMBOL_50_: yysymbol_kind_t = 50;
pub const YYSYMBOL_UMINUS: yysymbol_kind_t = 49;
pub const YYSYMBOL_FLTCAST: yysymbol_kind_t = 48;
pub const YYSYMBOL_INTCAST: yysymbol_kind_t = 47;
pub const YYSYMBOL_NOT: yysymbol_kind_t = 46;
pub const YYSYMBOL_POWER: yysymbol_kind_t = 45;
pub const YYSYMBOL_XOR: yysymbol_kind_t = 44;
pub const YYSYMBOL_43_: yysymbol_kind_t = 43;
pub const YYSYMBOL_42_: yysymbol_kind_t = 42;
pub const YYSYMBOL_41_: yysymbol_kind_t = 41;
pub const YYSYMBOL_40_: yysymbol_kind_t = 40;
pub const YYSYMBOL_39_: yysymbol_kind_t = 39;
pub const YYSYMBOL_38_: yysymbol_kind_t = 38;
pub const YYSYMBOL_37_: yysymbol_kind_t = 37;
pub const YYSYMBOL_GTE: yysymbol_kind_t = 36;
pub const YYSYMBOL_LTE: yysymbol_kind_t = 35;
pub const YYSYMBOL_LT: yysymbol_kind_t = 34;
pub const YYSYMBOL_GT: yysymbol_kind_t = 33;
pub const YYSYMBOL_32_: yysymbol_kind_t = 32;
pub const YYSYMBOL_NE: yysymbol_kind_t = 31;
pub const YYSYMBOL_EQ: yysymbol_kind_t = 30;
pub const YYSYMBOL_AND: yysymbol_kind_t = 29;
pub const YYSYMBOL_OR: yysymbol_kind_t = 28;
pub const YYSYMBOL_27_: yysymbol_kind_t = 27;
pub const YYSYMBOL_26_: yysymbol_kind_t = 26;
pub const YYSYMBOL_25_: yysymbol_kind_t = 25;
pub const YYSYMBOL_24_: yysymbol_kind_t = 24;
pub const YYSYMBOL_23_: yysymbol_kind_t = 23;
pub const YYSYMBOL_22_: yysymbol_kind_t = 22;
pub const YYSYMBOL_SNULLREF: yysymbol_kind_t = 21;
pub const YYSYMBOL_NULLREF: yysymbol_kind_t = 20;
pub const YYSYMBOL_ROWREF: yysymbol_kind_t = 19;
pub const YYSYMBOL_BITCOL: yysymbol_kind_t = 18;
pub const YYSYMBOL_SCOLUMN: yysymbol_kind_t = 17;
pub const YYSYMBOL_BCOLUMN: yysymbol_kind_t = 16;
pub const YYSYMBOL_COLUMN: yysymbol_kind_t = 15;
pub const YYSYMBOL_REGFILTER: yysymbol_kind_t = 14;
pub const YYSYMBOL_GTIFIND: yysymbol_kind_t = 13;
pub const YYSYMBOL_GTIOVERLAP: yysymbol_kind_t = 12;
pub const YYSYMBOL_GTIFILTER: yysymbol_kind_t = 11;
pub const YYSYMBOL_IFUNCTION: yysymbol_kind_t = 10;
pub const YYSYMBOL_BFUNCTION: yysymbol_kind_t = 9;
pub const YYSYMBOL_FUNCTION: yysymbol_kind_t = 8;
pub const YYSYMBOL_BITSTR: yysymbol_kind_t = 7;
pub const YYSYMBOL_STRING: yysymbol_kind_t = 6;
pub const YYSYMBOL_DOUBLE: yysymbol_kind_t = 5;
pub const YYSYMBOL_LONG: yysymbol_kind_t = 4;
pub const YYSYMBOL_BOOLEAN: yysymbol_kind_t = 3;
pub const YYSYMBOL_YYUNDEF: yysymbol_kind_t = 2;
pub const YYSYMBOL_YYerror: yysymbol_kind_t = 1;
pub const YYSYMBOL_YYEOF: yysymbol_kind_t = 0;
pub const YYSYMBOL_YYEMPTY: yysymbol_kind_t = -2;
pub type yytype_int8 = c_schar;
pub type yy_state_fast_t = c_int;
pub type funcOp = c_uint;
pub const array_fct: funcOp = 1051;
pub const axiselem_fct: funcOp = 1050;
pub const elemnum_fct: funcOp = 1049;
pub const gtifind_fct: funcOp = 1048;
pub const gtiover_fct: funcOp = 1047;
pub const setnull_fct: funcOp = 1046;
pub const strpos_fct: funcOp = 1045;
pub const strmid_fct: funcOp = 1044;
pub const poirnd_fct: funcOp = 1043;
pub const gasrnd_fct: funcOp = 1042;
pub const angsep_fct: funcOp = 1041;
pub const nonnull_fct: funcOp = 1040;
pub const stddev_fct: funcOp = 1039;
pub const average_fct: funcOp = 1038;
pub const median_fct: funcOp = 1037;
pub const null_fct: funcOp = 1036;
pub const row_fct: funcOp = 1035;
pub const ifthenelse_fct: funcOp = 1034;
pub const regfilt_fct: funcOp = 1033;
pub const gtifilt_fct: funcOp = 1032;
pub const defnull_fct: funcOp = 1031;
pub const isnull_fct: funcOp = 1030;
pub const elps_fct: funcOp = 1029;
pub const box_fct: funcOp = 1028;
pub const circle_fct: funcOp = 1027;
pub const near_fct: funcOp = 1026;
pub const max2_fct: funcOp = 1025;
pub const max1_fct: funcOp = 1024;
pub const min2_fct: funcOp = 1023;
pub const min1_fct: funcOp = 1022;
pub const round_fct: funcOp = 1021;
pub const floor_fct: funcOp = 1020;
pub const ceil_fct: funcOp = 1019;
pub const atan2_fct: funcOp = 1018;
pub const abs_fct: funcOp = 1017;
pub const sqrt_fct: funcOp = 1016;
pub const log10_fct: funcOp = 1015;
pub const log_fct: funcOp = 1014;
pub const exp_fct: funcOp = 1013;
pub const tanh_fct: funcOp = 1012;
pub const cosh_fct: funcOp = 1011;
pub const sinh_fct: funcOp = 1010;
pub const atan_fct: funcOp = 1009;
pub const acos_fct: funcOp = 1008;
pub const asin_fct: funcOp = 1007;
pub const tan_fct: funcOp = 1006;
pub const cos_fct: funcOp = 1005;
pub const sin_fct: funcOp = 1004;
pub const nelem_fct: funcOp = 1003;
pub const sum_fct: funcOp = 1002;
pub const rnd_fct: funcOp = 1001;

pub type shapeType = c_uint;
pub const bpanda_rgn: shapeType = 14;
pub const epanda_rgn: shapeType = 13;
pub const panda_rgn: shapeType = 12;
pub const poly_rgn: shapeType = 11;
pub const sector_rgn: shapeType = 10;
pub const diamond_rgn: shapeType = 9;
pub const rectangle_rgn: shapeType = 8;
pub const boxannulus_rgn: shapeType = 7;
pub const box_rgn: shapeType = 6;
pub const elliptannulus_rgn: shapeType = 5;
pub const ellipse_rgn: shapeType = 4;
pub const annulus_rgn: shapeType = 3;
pub const circle_rgn: shapeType = 2;
pub const line_rgn: shapeType = 1;
pub const point_rgn: shapeType = 0;
pub type yytype_uint8 = c_uchar;
#[derive(Copy, Clone)]
#[repr(C)]
pub union yyalloc {
    pub yyss_alloc: yy_state_t,
    pub yyvs_alloc: FITS_PARSER_YYSTYPE,
}

/* YYTRANSLATE[TOKEN-NUM] -- Symbol number corresponding to TOKEN-NUM
as returned by yylex.  */
static yytranslate: [yytype_int8; 293] = [
    0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 53, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
    2, 2, 2, 2, 2, 2, 39, 43, 2, 55, 56, 40, 37, 22, 38, 2, 41, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 24,
    2, 2, 23, 2, 27, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
    2, 50, 2, 54, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
    2, 2, 25, 42, 26, 32, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
    2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
    20, 21, 28, 29, 30, 31, 33, 34, 35, 36, 44, 45, 46, 47, 48, 49, 51, 52,
];

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
STATE-NUM.  */
static yypact: [yytype_int16; 322] = [
    -41, 316, -41, -40, -41, -41, -41, -41, -41, 369, 423, 423, -5, 15, -4, 27, 36, 38, 40, 41,
    -41, -41, -41, 423, 423, 423, 423, 423, 423, -41, 423, -41, -7, 10, 1226, 81, 1646, 83, -41,
    -41, 450, 116, 309, 12, 479, 185, 152, 222, 1593, 1673, 1675, -19, -41, 13, -18, -41, 6, 423,
    423, 423, 423, 1593, 1673, 1684, 17, 17, 19, 24, 17, 19, 17, 19, 710, 1253, 1611, 365, 423,
    -41, 423, -41, 423, 423, 423, 423, 423, 423, 423, 423, 423, 423, 423, 423, 423, 423, 423, 423,
    423, 423, -41, 423, 423, 423, 423, 423, 423, 423, -41, -2, -2, -2, -2, -2, -2, -2, -2, -2, 423,
    -41, 423, 423, 423, 423, 423, 423, 423, -41, 423, -41, 423, -41, -41, 423, -41, 423, -41, -41,
    -41, 423, 423, -41, 423, 423, -41, 423, -41, 1455, 1478, 1501, 1524, -41, -41, -41, -41, 1593,
    1673, 1593, 1673, 1547, 1712, 1712, 1712, 1726, 1726, 1726, 1726, 368, 368, 368, 28, 19, 28, 5,
    5, 5, 5, 851, 1570, 425, 260, 128, -20, 14, 14, 28, 876, -2, -2, -25, -25, -25, -25, -25, -25,
    -36, 24, 24, 901, 140, 140, 39, 39, 39, 39, -41, 508, 738, 1258, 1288, 1629, 1312, 1638, 537,
    1336, 566, 1360, -41, -41, -41, -41, 423, 423, -41, 423, 423, 423, 423, -41, 24, 189, 423, -41,
    423, -41, -41, -41, 423, -41, 423, -41, 93, -41, 423, 94, -41, 423, 1694, 926, 1694, 1673,
    1694, 1673, 1684, 951, 976, 1384, 766, 595, 79, 624, 80, 653, 423, -41, 423, -41, 423, -41,
    423, -41, 423, -41, 100, 101, -41, 117, 118, -41, 1001, 1026, 1051, 794, 1408, 72, 111, 85, 99,
    423, -41, 423, -41, 423, -41, -41, 423, -41, 129, -41, -41, 1076, 1101, 1126, 682, 104, 423,
    -41, 423, -41, 423, -41, 423, -41, -41, 1151, 1176, 1201, 1432, -41, -41, -41, 423, 822, -41,
];

/* YYDEFACT[STATE-NUM] -- Default reduction number in state STATE-NUM.
Performed when YYTABLE does not specify something else to do.  Zero
means the default is an error.  */
static yydefact: [yytype_uint8; 322] = [
    2, 0, 1, 0, 74, 31, 32, 127, 18, 0, 0, 0, 0, 0, 0, 0, 33, 75, 128, 19, 35, 36, 130, 0, 0, 0, 0,
    0, 0, 4, 0, 3, 0, 0, 0, 0, 0, 0, 9, 54, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 109, 0, 0, 113, 0,
    0, 0, 0, 0, 12, 10, 0, 46, 47, 125, 29, 70, 71, 72, 73, 0, 0, 0, 0, 0, 17, 0, 16, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 8, 0, 0, 0, 0, 0, 0, 0, 7, 0, 59, 0, 55, 58, 0, 57, 0, 102, 103, 104, 0, 0, 110, 0, 0, 114,
    0, 117, 0, 0, 0, 0, 48, 126, 30, 131, 15, 11, 13, 14, 0, 88, 89, 87, 83, 84, 86, 85, 38, 39,
    37, 40, 49, 41, 43, 42, 44, 45, 0, 0, 0, 0, 97, 96, 98, 99, 50, 0, 0, 0, 77, 78, 81, 79, 80,
    82, 23, 22, 21, 0, 90, 91, 92, 94, 95, 93, 132, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 34, 76, 129,
    20, 0, 0, 65, 0, 0, 0, 0, 120, 29, 0, 0, 24, 0, 61, 56, 105, 0, 134, 0, 60, 0, 111, 0, 0, 115,
    0, 100, 0, 51, 53, 52, 101, 133, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 66, 0, 121, 0, 25, 0, 135, 0,
    106, 0, 0, 63, 0, 0, 118, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 67, 0, 122, 0, 26, 62, 0, 112, 0, 116,
    119, 0, 0, 0, 0, 0, 0, 68, 0, 123, 0, 27, 0, 107, 64, 0, 0, 0, 0, 69, 124, 28, 0, 0, 108,
];

/* YYPGOTO[NTERM-NUM].  */
static yypgoto: [yytype_int16; 9] = [-41, -41, -41, -41, -41, -1, 170, 96, 30];

/* YYDEFGOTO[NTERM-NUM].  */
static yydefgoto: [yytype_int8; 9] = [0, 1, 31, 32, 33, 48, 49, 46, 63];

/* YYTABLE[YYPACT[STATE-NUM]] -- What to do in state STATE-NUM.  If
positive, shift that token.  If negative, reduce the rule whose
number is the opposite.  If YYTABLE_NINF, syntax error.  */
static yytable: [yytype_int16; 1777] = [
    34, 51, 54, 138, 141, 8, 114, 115, 40, 44, 102, 103, 113, 38, 116, 76, 19, 114, 115, 77, 104,
    53, 61, 64, 65, 116, 68, 70, 143, 72, 105, 37, 78, 56, 131, 140, 79, 139, 142, 43, 47, 50, 118,
    119, 185, 120, 121, 122, 123, 124, 96, 52, 55, 186, 104, 97, 145, 146, 147, 148, 75, 57, 144,
    58, 105, 59, 60, 97, 132, 105, 93, 94, 95, 96, 116, 153, 124, 155, 97, 157, 158, 159, 160, 161,
    162, 163, 164, 165, 166, 167, 168, 170, 171, 172, 173, 174, 175, 36, 176, 257, 259, 271, 274,
    183, 184, 42, 282, 283, 99, 100, 101, 102, 103, 118, 119, 196, 120, 121, 122, 123, 124, 104,
    67, 284, 285, 204, 74, 205, 294, 178, 207, 105, 209, 295, 106, 302, 125, 211, 128, 212, 213,
    296, 214, 99, 100, 101, 102, 103, 197, 198, 199, 200, 201, 202, 203, 297, 104, 101, 102, 103,
    311, 208, 0, 0, 0, 0, 105, 210, 104, 0, 0, 35, 129, 120, 121, 122, 123, 124, 105, 41, 45, 0,
    107, 108, 0, 109, 110, 111, 112, 113, 0, 0, 0, 62, 114, 115, 66, 69, 71, 0, 73, 0, 116, 187,
    188, 189, 190, 191, 192, 193, 194, 195, 99, 100, 101, 102, 103, 0, 245, 246, 0, 247, 249, 0,
    252, 104, 113, 0, 253, 0, 254, 114, 115, 0, 255, 105, 256, 0, 0, 116, 258, 135, 0, 260, 0, 151,
    154, 0, 156, 0, 0, 0, 118, 119, 251, 120, 121, 122, 123, 124, 277, 169, 278, 0, 279, 0, 280, 0,
    281, 177, 179, 180, 181, 182, 0, 0, 0, 0, 136, 0, 0, 227, 228, 0, 224, 298, 0, 299, 0, 300,
    118, 119, 301, 120, 121, 122, 123, 124, 206, 0, 0, 0, 312, 0, 313, 0, 314, 0, 315, 0, 0, 0, 0,
    0, 0, 0, 2, 3, 320, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 0,
    107, 108, 23, 109, 110, 111, 112, 113, 0, 0, 0, 0, 114, 115, 24, 25, 0, 0, 0, 0, 116, 0, 0, 26,
    27, 28, 130, 0, 0, 0, 29, 0, 30, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
    21, 22, 0, 248, 250, 23, 118, 119, 0, 120, 121, 122, 123, 124, 0, 0, 0, 24, 25, 91, 92, 93, 94,
    95, 96, 0, 26, 27, 28, 97, 0, 0, 152, 0, 0, 30, 39, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
    16, 17, 18, 19, 20, 21, 22, 0, 0, 0, 23, 223, 0, 0, 99, 100, 101, 102, 103, 0, 0, 0, 24, 25, 0,
    0, 0, 104, 0, 0, 0, 26, 27, 28, 126, 80, 0, 105, 0, 0, 30, 0, 81, 82, 83, 84, 85, 86, 87, 88,
    89, 90, 91, 92, 93, 94, 95, 96, 0, 0, 0, 0, 97, 133, 80, 0, 0, 0, 127, 0, 0, 81, 82, 83, 84,
    85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 0, 0, 0, 0, 97, 231, 80, 0, 0, 0, 134, 0, 0,
    81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 0, 0, 0, 0, 97, 239, 80, 0, 0,
    0, 232, 0, 0, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 0, 0, 0, 0, 97,
    242, 80, 0, 0, 0, 240, 0, 0, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 0,
    0, 0, 0, 97, 269, 80, 0, 0, 0, 243, 0, 0, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93,
    94, 95, 96, 0, 0, 0, 0, 97, 272, 80, 0, 0, 0, 270, 0, 0, 81, 82, 83, 84, 85, 86, 87, 88, 89,
    90, 91, 92, 93, 94, 95, 96, 0, 0, 0, 0, 97, 275, 80, 0, 0, 0, 273, 0, 0, 81, 82, 83, 84, 85,
    86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 0, 0, 0, 0, 97, 309, 80, 0, 0, 0, 276, 0, 0, 81,
    82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 0, 0, 0, 0, 97, 80, 0, 0, 0, 0,
    310, 0, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 0, 0, 0, 0, 97, 80, 0,
    0, 0, 0, 149, 0, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 0, 0, 0, 0,
    97, 80, 0, 0, 0, 0, 233, 0, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 0,
    0, 0, 0, 97, 80, 0, 0, 0, 0, 268, 0, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94,
    95, 96, 0, 0, 0, 0, 97, 80, 0, 0, 0, 0, 292, 0, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92,
    93, 94, 95, 96, 0, 0, 0, 0, 97, 220, 80, 0, 0, 0, 321, 0, 0, 81, 82, 83, 84, 85, 86, 87, 88,
    89, 90, 91, 92, 93, 94, 95, 96, 0, 225, 80, 0, 97, 0, 0, 0, 221, 81, 82, 83, 84, 85, 86, 87,
    88, 89, 90, 91, 92, 93, 94, 95, 96, 0, 229, 80, 0, 97, 0, 0, 0, 226, 81, 82, 83, 84, 85, 86,
    87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 0, 261, 80, 0, 97, 0, 0, 0, 230, 81, 82, 83, 84, 85,
    86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 0, 263, 80, 0, 97, 0, 0, 0, 262, 81, 82, 83, 84,
    85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 0, 265, 80, 0, 97, 0, 0, 0, 264, 81, 82, 83,
    84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 0, 286, 80, 0, 97, 0, 0, 0, 266, 81, 82,
    83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 0, 288, 80, 0, 97, 0, 0, 0, 287, 81,
    82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 0, 290, 80, 0, 97, 0, 0, 0, 289,
    81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 0, 303, 80, 0, 97, 0, 0, 0,
    291, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 0, 305, 80, 0, 97, 0, 0,
    0, 304, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 0, 307, 80, 0, 97, 0,
    0, 0, 306, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 0, 0, 80, 0, 97, 0,
    0, 0, 308, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 0, 0, 80, 0, 97, 0,
    0, 0, 316, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 0, 0, 80, 0, 97, 0,
    0, 0, 317, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 0, 0, 80, 0, 97, 0,
    0, 0, 318, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 0, 0, 0, 0, 97, 0,
    0, 98, 99, 100, 101, 102, 103, 99, 100, 101, 102, 103, 0, 0, 0, 104, 0, 0, 0, 0, 104, 0, 0, 0,
    0, 105, 0, 0, 0, 0, 105, 150, 235, 80, 0, 0, 234, 0, 0, 0, 81, 82, 83, 84, 85, 86, 87, 88, 89,
    90, 91, 92, 93, 94, 95, 96, 237, 80, 0, 0, 97, 0, 0, 0, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90,
    91, 92, 93, 94, 95, 96, 241, 80, 0, 0, 97, 0, 0, 0, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91,
    92, 93, 94, 95, 96, 244, 80, 0, 0, 97, 0, 0, 0, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92,
    93, 94, 95, 96, 267, 80, 0, 0, 97, 0, 0, 0, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93,
    94, 95, 96, 293, 80, 0, 0, 97, 0, 0, 0, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94,
    95, 96, 319, 80, 0, 0, 97, 0, 0, 0, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95,
    96, 80, 0, 0, 215, 97, 0, 0, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96,
    80, 0, 0, 216, 97, 0, 0, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 80, 0,
    0, 217, 97, 0, 0, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 80, 0, 0,
    218, 97, 0, 0, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 80, 219, 0, 0,
    97, 0, 0, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 80, 222, 0, 0, 97, 0,
    0, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 80, 0, 0, 0, 97, 0, 0, 81,
    82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 0, 0, 107, 108, 97, 109, 110, 111,
    112, 113, 0, 0, 0, 0, 114, 115, 0, 0, 0, 0, 118, 119, 116, 120, 121, 122, 123, 124, 151, 118,
    119, 0, 120, 121, 122, 123, 124, 107, 108, 0, 109, 110, 111, 112, 113, 0, 236, 0, 0, 114, 115,
    0, 0, 0, 0, 238, 0, 116, 137, 0, 117, 99, 100, 101, 102, 103, 118, 119, 0, 120, 121, 122, 123,
    124, 104, 118, 119, 0, 120, 121, 122, 123, 124, 0, 105, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90,
    91, 92, 93, 94, 95, 96, 0, 0, 0, 0, 97, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 0,
    0, 0, 0, 97, 88, 89, 90, 91, 92, 93, 94, 95, 96, 0, 0, 0, 0, 97,
];

static yycheck: [yytype_int16; 1777] = [
    1, 6, 6, 22, 22, 7, 42, 43, 9, 10, 30, 31, 37, 53, 50, 22, 18, 42, 43, 26, 40, 6, 23, 24, 25,
    50, 27, 28, 22, 30, 50, 1, 22, 6, 22, 22, 26, 56, 56, 9, 10, 11, 30, 31, 46, 33, 34, 35, 36,
    37, 45, 56, 56, 55, 40, 50, 57, 58, 59, 60, 30, 25, 56, 25, 50, 25, 25, 50, 56, 50, 42, 43, 44,
    45, 50, 76, 37, 78, 50, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97,
    1, 99, 6, 6, 22, 22, 104, 105, 9, 6, 6, 27, 28, 29, 30, 31, 30, 31, 116, 33, 34, 35, 36, 37,
    40, 26, 6, 6, 126, 30, 128, 56, 99, 131, 50, 133, 22, 53, 6, 53, 138, 22, 140, 141, 56, 143,
    27, 28, 29, 30, 31, 118, 119, 120, 121, 122, 123, 124, 56, 40, 29, 30, 31, 56, 131, -1, -1, -1,
    -1, 50, 137, 40, -1, -1, 1, 56, 33, 34, 35, 36, 37, 50, 9, 10, -1, 30, 31, -1, 33, 34, 35, 36,
    37, -1, -1, -1, 23, 42, 43, 26, 27, 28, -1, 30, -1, 50, 107, 108, 109, 110, 111, 112, 113, 114,
    115, 27, 28, 29, 30, 31, -1, 219, 220, -1, 222, 223, -1, 225, 40, 37, -1, 229, -1, 231, 42, 43,
    -1, 235, 50, 237, -1, -1, 50, 241, 56, -1, 244, -1, 56, 76, -1, 78, -1, -1, -1, 30, 31, 224,
    33, 34, 35, 36, 37, 261, 91, 263, -1, 265, -1, 267, -1, 269, 99, 100, 101, 102, 103, -1, -1,
    -1, -1, 56, -1, -1, 185, 186, -1, 24, 286, -1, 288, -1, 290, 30, 31, 293, 33, 34, 35, 36, 37,
    128, -1, -1, -1, 303, -1, 305, -1, 307, -1, 309, -1, -1, -1, -1, -1, -1, -1, 0, 1, 319, 3, 4,
    5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, -1, 30, 31, 25, 33, 34, 35, 36,
    37, -1, -1, -1, -1, 42, 43, 37, 38, -1, -1, -1, -1, 50, -1, -1, 46, 47, 48, 56, -1, -1, -1, 53,
    -1, 55, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, -1, 222, 223, 25,
    30, 31, -1, 33, 34, 35, 36, 37, -1, -1, -1, 37, 38, 40, 41, 42, 43, 44, 45, -1, 46, 47, 48, 50,
    -1, -1, 56, -1, -1, 55, 56, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
    21, -1, -1, -1, 25, 24, -1, -1, 27, 28, 29, 30, 31, -1, -1, -1, 37, 38, -1, -1, -1, 40, -1, -1,
    -1, 46, 47, 48, 22, 23, -1, 50, -1, -1, 55, -1, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41,
    42, 43, 44, 45, -1, -1, -1, -1, 50, 22, 23, -1, -1, -1, 56, -1, -1, 30, 31, 32, 33, 34, 35, 36,
    37, 38, 39, 40, 41, 42, 43, 44, 45, -1, -1, -1, -1, 50, 22, 23, -1, -1, -1, 56, -1, -1, 30, 31,
    32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, -1, -1, -1, -1, 50, 22, 23, -1, -1, -1,
    56, -1, -1, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, -1, -1, -1, -1, 50,
    22, 23, -1, -1, -1, 56, -1, -1, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45,
    -1, -1, -1, -1, 50, 22, 23, -1, -1, -1, 56, -1, -1, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40,
    41, 42, 43, 44, 45, -1, -1, -1, -1, 50, 22, 23, -1, -1, -1, 56, -1, -1, 30, 31, 32, 33, 34, 35,
    36, 37, 38, 39, 40, 41, 42, 43, 44, 45, -1, -1, -1, -1, 50, 22, 23, -1, -1, -1, 56, -1, -1, 30,
    31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, -1, -1, -1, -1, 50, 22, 23, -1, -1,
    -1, 56, -1, -1, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, -1, -1, -1, -1,
    50, 23, -1, -1, -1, -1, 56, -1, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45,
    -1, -1, -1, -1, 50, 23, -1, -1, -1, -1, 56, -1, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41,
    42, 43, 44, 45, -1, -1, -1, -1, 50, 23, -1, -1, -1, -1, 56, -1, 30, 31, 32, 33, 34, 35, 36, 37,
    38, 39, 40, 41, 42, 43, 44, 45, -1, -1, -1, -1, 50, 23, -1, -1, -1, -1, 56, -1, 30, 31, 32, 33,
    34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, -1, -1, -1, -1, 50, 23, -1, -1, -1, -1, 56, -1,
    30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, -1, -1, -1, -1, 50, 22, 23, -1,
    -1, -1, 56, -1, -1, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, -1, 22, 23,
    -1, 50, -1, -1, -1, 54, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, -1, 22,
    23, -1, 50, -1, -1, -1, 54, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, -1,
    22, 23, -1, 50, -1, -1, -1, 54, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45,
    -1, 22, 23, -1, 50, -1, -1, -1, 54, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44,
    45, -1, 22, 23, -1, 50, -1, -1, -1, 54, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43,
    44, 45, -1, 22, 23, -1, 50, -1, -1, -1, 54, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42,
    43, 44, 45, -1, 22, 23, -1, 50, -1, -1, -1, 54, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41,
    42, 43, 44, 45, -1, 22, 23, -1, 50, -1, -1, -1, 54, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40,
    41, 42, 43, 44, 45, -1, 22, 23, -1, 50, -1, -1, -1, 54, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39,
    40, 41, 42, 43, 44, 45, -1, 22, 23, -1, 50, -1, -1, -1, 54, 30, 31, 32, 33, 34, 35, 36, 37, 38,
    39, 40, 41, 42, 43, 44, 45, -1, 22, 23, -1, 50, -1, -1, -1, 54, 30, 31, 32, 33, 34, 35, 36, 37,
    38, 39, 40, 41, 42, 43, 44, 45, -1, -1, 23, -1, 50, -1, -1, -1, 54, 30, 31, 32, 33, 34, 35, 36,
    37, 38, 39, 40, 41, 42, 43, 44, 45, -1, -1, 23, -1, 50, -1, -1, -1, 54, 30, 31, 32, 33, 34, 35,
    36, 37, 38, 39, 40, 41, 42, 43, 44, 45, -1, -1, 23, -1, 50, -1, -1, -1, 54, 30, 31, 32, 33, 34,
    35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, -1, -1, 23, -1, 50, -1, -1, -1, 54, 30, 31, 32, 33,
    34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, -1, -1, -1, -1, 50, -1, -1, 53, 27, 28, 29, 30,
    31, 27, 28, 29, 30, 31, -1, -1, -1, 40, -1, -1, -1, -1, 40, -1, -1, -1, -1, 50, -1, -1, -1, -1,
    50, 56, 22, 23, -1, -1, 56, -1, -1, -1, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43,
    44, 45, 22, 23, -1, -1, 50, -1, -1, -1, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43,
    44, 45, 22, 23, -1, -1, 50, -1, -1, -1, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43,
    44, 45, 22, 23, -1, -1, 50, -1, -1, -1, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43,
    44, 45, 22, 23, -1, -1, 50, -1, -1, -1, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43,
    44, 45, 22, 23, -1, -1, 50, -1, -1, -1, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43,
    44, 45, 22, 23, -1, -1, 50, -1, -1, -1, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43,
    44, 45, 23, -1, -1, 26, 50, -1, -1, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44,
    45, 23, -1, -1, 26, 50, -1, -1, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45,
    23, -1, -1, 26, 50, -1, -1, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 23,
    -1, -1, 26, 50, -1, -1, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 23, 24,
    -1, -1, 50, -1, -1, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 23, 24, -1,
    -1, 50, -1, -1, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 23, -1, -1, -1,
    50, -1, -1, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, -1, -1, 30, 31, 50,
    33, 34, 35, 36, 37, -1, -1, -1, -1, 42, 43, -1, -1, -1, -1, 30, 31, 50, 33, 34, 35, 36, 37, 56,
    30, 31, -1, 33, 34, 35, 36, 37, 30, 31, -1, 33, 34, 35, 36, 37, -1, 56, -1, -1, 42, 43, -1, -1,
    -1, -1, 56, -1, 50, 22, -1, 53, 27, 28, 29, 30, 31, 30, 31, -1, 33, 34, 35, 36, 37, 40, 30, 31,
    -1, 33, 34, 35, 36, 37, -1, 50, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45,
    -1, -1, -1, -1, 50, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, -1, -1, -1, -1, 50, 37,
    38, 39, 40, 41, 42, 43, 44, 45, -1, -1, -1, -1, 50,
];

/* YYSTOS[STATE-NUM] -- The symbol kind of the accessing symbol of
state STATE-NUM.  */
static yystos: [yytype_int8; 322] = [
    0, 58, 0, 1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 25, 37, 38,
    46, 47, 48, 53, 55, 59, 60, 61, 62, 63, 64, 65, 53, 56, 62, 63, 64, 65, 62, 63, 64, 65, 62, 63,
    65, 6, 56, 6, 6, 56, 6, 25, 25, 25, 25, 62, 63, 65, 62, 62, 63, 64, 62, 63, 62, 63, 62, 63, 64,
    65, 22, 26, 22, 26, 23, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 50, 53,
    27, 28, 29, 30, 31, 40, 50, 53, 30, 31, 33, 34, 35, 36, 37, 42, 43, 50, 53, 30, 31, 33, 34, 35,
    36, 37, 53, 22, 56, 22, 56, 56, 22, 56, 22, 56, 56, 56, 22, 22, 56, 22, 22, 56, 22, 56, 62, 62,
    62, 62, 56, 56, 56, 56, 62, 63, 62, 63, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 63, 62,
    62, 62, 62, 62, 62, 62, 63, 65, 63, 63, 63, 63, 62, 62, 46, 55, 64, 64, 64, 64, 64, 64, 64, 64,
    64, 62, 65, 65, 65, 65, 65, 65, 65, 62, 62, 63, 62, 65, 62, 65, 62, 62, 62, 62, 26, 26, 26, 26,
    24, 22, 54, 24, 24, 24, 22, 54, 64, 64, 22, 54, 22, 56, 56, 56, 22, 56, 22, 56, 22, 56, 22, 22,
    56, 22, 62, 62, 62, 63, 62, 63, 65, 62, 62, 62, 62, 62, 6, 62, 6, 62, 22, 54, 22, 54, 22, 54,
    22, 56, 22, 56, 22, 22, 56, 22, 22, 56, 62, 62, 62, 62, 62, 6, 6, 6, 6, 22, 54, 22, 54, 22, 54,
    56, 22, 56, 22, 56, 56, 62, 62, 62, 62, 6, 22, 54, 22, 54, 22, 54, 22, 56, 56, 62, 62, 62, 62,
    54, 54, 54, 22, 62, 56,
];

/* YYR1[RULE-NUM] -- Symbol kind of the left-hand side of rule RULE-NUM.  */
static yyr1: [yytype_int8; 136] = [
    0, 57, 58, 58, 59, 59, 59, 59, 59, 59, 60, 60, 61, 61, 61, 61, 62, 63, 64, 64, 64, 64, 64, 64,
    64, 64, 64, 64, 64, 64, 64, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62,
    62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62,
    62, 62, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63,
    63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63,
    63, 63, 63, 63, 63, 63, 63, 65, 65, 65, 65, 65, 65, 65, 65, 65,
];

/* YYR2[RULE-NUM] -- Number of symbols on the right-hand side of rule RULE-NUM.  */
static yyr2: [yytype_int8; 136] = [
    0, 2, 0, 2, 1, 2, 2, 2, 2, 2, 2, 3, 2, 3, 3, 3, 2, 2, 1, 1, 4, 3, 3, 3, 4, 6, 8, 10, 12, 2, 3,
    1, 1, 1, 4, 1, 1, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 3, 3, 3, 5, 5, 5, 2, 3, 5, 3, 3, 3, 5, 5, 9,
    7, 11, 4, 6, 8, 10, 12, 2, 2, 2, 2, 1, 1, 4, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 3, 3, 5, 5, 3, 3, 3, 5, 7, 11, 15, 2, 3, 5, 9, 2, 3, 5, 9, 3, 7, 9, 4, 6, 8, 10,
    12, 2, 3, 1, 1, 4, 1, 3, 3, 5, 5, 7,
];

unsafe fn Alloc_Node(mut lParse: *mut ParseData) -> c_int {
    unsafe {
        let mut newNodePtr: *mut Node = std::ptr::null_mut::<Node>();
        if (*lParse).nNodes == (*lParse).nNodesAlloc {
            if !((*lParse).Nodes).is_null() {
                (*lParse).nNodesAlloc += (*lParse).nNodesAlloc;
                newNodePtr = realloc(
                    (*lParse).Nodes as *mut c_void,
                    (::core::mem::size_of::<Node>() as c_ulong)
                        .wrapping_mul((*lParse).nNodesAlloc as c_ulong)
                        .try_into()
                        .unwrap(),
                ) as *mut Node;
            } else {
                (*lParse).nNodesAlloc = 100 as c_int;
                newNodePtr = malloc(
                    (::core::mem::size_of::<Node>() as c_ulong)
                        .wrapping_mul((*lParse).nNodesAlloc as c_ulong)
                        .try_into()
                        .unwrap(),
                ) as *mut Node;
            }
            if !newNodePtr.is_null() {
                (*lParse).Nodes = newNodePtr;
            } else {
                (*lParse).status = 113 as c_int;
                return -(1);
            }
        }
        let fresh0 = (*lParse).nNodes;
        (*lParse).nNodes += 1;
        fresh0
    }
}
unsafe fn Free_Last_Node(mut lParse: *mut ParseData) {
    unsafe {
        if (*lParse).nNodes != 0 {
            (*lParse).nNodes -= 1;
            (*lParse).nNodes;
        }
    }
}
unsafe fn New_Const(
    mut lParse: *mut ParseData,
    mut returnType: c_int,
    mut value: *mut c_void,
    mut len: c_long,
) -> c_int {
    unsafe {
        let mut this: *mut Node = std::ptr::null_mut::<Node>();
        let mut n: c_int = 0;
        n = Alloc_Node(lParse);
        if n >= 0 {
            this = ((*lParse).Nodes).offset(n as isize);
            (*this).operation = -(1000 as c_int);
            (*this).DoOp = None;
            (*this).nSubNodes = 0;
            (*this).ntype = returnType;
            memcpy(
                &mut (*this).value.data as *const _ as *mut c_void,
                value,
                (len as c_ulong).try_into().unwrap(),
            );
            (*this).value.undef = ptr::null_mut();
            (*this).value.nelem = 1;
            (*this).value.naxis = 1;
            (*this).value.naxes[0] = 1;
        }
        n
    }
}
unsafe fn New_Column(mut lParse: *mut ParseData, mut ColNum: c_int) -> c_int {
    unsafe {
        let mut this: *mut Node = std::ptr::null_mut::<Node>();
        let mut n: c_int = 0;
        let mut i: c_int = 0;
        n = Alloc_Node(lParse);
        if n >= 0 {
            this = ((*lParse).Nodes).offset(n as isize);
            (*this).operation = -ColNum;
            (*this).DoOp = None;
            (*this).nSubNodes = 0;
            (*this).ntype = (*((*lParse).varData).offset(ColNum as isize)).dtype;
            (*this).value.nelem = (*((*lParse).varData).offset(ColNum as isize)).nelem;
            (*this).value.naxis = (*((*lParse).varData).offset(ColNum as isize)).naxis;
            i = 0;
            while i < (*((*lParse).varData).offset(ColNum as isize)).naxis {
                (*this).value.naxes[i as usize] =
                    (*((*lParse).varData).offset(ColNum as isize)).naxes[i as usize];
                i += 1;
                i;
            }
        }
        n
    }
}
unsafe fn New_Offset(
    mut lParse: *mut ParseData,
    mut ColNum: c_int,
    mut offsetNode: c_int,
) -> c_int {
    unsafe {
        let mut this: *mut Node = std::ptr::null_mut::<Node>();
        let mut n: c_int = 0;
        let mut i: c_int = 0;
        let mut colNode: c_int = 0;
        colNode = New_Column(lParse, ColNum);
        if colNode < 0 {
            return -(1);
        }
        n = Alloc_Node(lParse);
        if n >= 0 {
            this = ((*lParse).Nodes).offset(n as isize);
            (*this).operation = b'{' as i32;
            (*this).DoOp = Some(Do_Offset);
            (*this).nSubNodes = 2 as c_int;
            (*this).SubNodes[0] = colNode;
            (*this).SubNodes[1] = offsetNode;
            (*this).ntype = (*((*lParse).varData).offset(ColNum as isize)).dtype;
            (*this).value.nelem = (*((*lParse).varData).offset(ColNum as isize)).nelem;
            (*this).value.naxis = (*((*lParse).varData).offset(ColNum as isize)).naxis;
            i = 0;
            while i < (*((*lParse).varData).offset(ColNum as isize)).naxis {
                (*this).value.naxes[i as usize] =
                    (*((*lParse).varData).offset(ColNum as isize)).naxes[i as usize];
                i += 1;
                i;
            }
        }
        n
    }
}
unsafe fn New_Unary(
    mut lParse: *mut ParseData,
    mut returnType: c_int,
    mut Op: c_int,
    mut Node1: c_int,
) -> c_int {
    unsafe {
        let mut this: *mut Node = std::ptr::null_mut::<Node>();
        let mut that: *mut Node = std::ptr::null_mut::<Node>();
        let mut i: c_int = 0;
        let mut n: c_int = 0;
        if Node1 < 0 {
            return -(1);
        }
        that = ((*lParse).Nodes).offset(Node1 as isize);
        if Op == 0 {
            Op = returnType;
        }
        if (Op == fits_parser_yytokentype::DOUBLE as c_int
            || Op == fits_parser_yytokentype::FLTCAST as c_int)
            && (*that).ntype == fits_parser_yytokentype::DOUBLE as c_int
        {
            return Node1;
        }
        if (Op == fits_parser_yytokentype::LONG as c_int
            || Op == fits_parser_yytokentype::INTCAST as c_int)
            && (*that).ntype == fits_parser_yytokentype::LONG as c_int
        {
            return Node1;
        }
        if Op == fits_parser_yytokentype::BOOLEAN as c_int
            && (*that).ntype == fits_parser_yytokentype::BOOLEAN as c_int
        {
            return Node1;
        }
        n = Alloc_Node(lParse);
        if n >= 0 {
            this = ((*lParse).Nodes).offset(n as isize);
            (*this).operation = Op;
            (*this).DoOp = Some(Do_Unary);
            (*this).nSubNodes = 1;
            (*this).SubNodes[0] = Node1;
            (*this).ntype = returnType;
            that = ((*lParse).Nodes).offset(Node1 as isize);
            (*this).value.nelem = (*that).value.nelem;
            (*this).value.naxis = (*that).value.naxis;
            i = 0;
            while i < (*that).value.naxis {
                (*this).value.naxes[i as usize] = (*that).value.naxes[i as usize];
                i += 1;
                i;
            }
            if (*that).operation == -(1000 as c_int) {
                ((*this).DoOp).expect("non-null function pointer")(lParse, this);
            }
        }
        n
    }
}
unsafe fn New_BinOp(
    mut lParse: *mut ParseData,
    mut returnType: c_int,
    mut Node1: c_int,
    mut Op: c_int,
    mut Node2: c_int,
) -> c_int {
    unsafe {
        let mut this: *mut Node = std::ptr::null_mut::<Node>();
        let mut that1: *mut Node = std::ptr::null_mut::<Node>();
        let mut that2: *mut Node = std::ptr::null_mut::<Node>();
        let mut n: c_int = 0;
        let mut i: c_int = 0;
        let mut constant: c_int = 0;
        if Node1 < 0 || Node2 < 0 {
            return -(1);
        }
        n = Alloc_Node(lParse);
        if n >= 0 {
            this = ((*lParse).Nodes).offset(n as isize);
            (*this).operation = Op;
            (*this).nSubNodes = 2 as c_int;
            (*this).SubNodes[0] = Node1;
            (*this).SubNodes[1] = Node2;
            (*this).ntype = returnType;
            that1 = ((*lParse).Nodes).offset(Node1 as isize);
            that2 = ((*lParse).Nodes).offset(Node2 as isize);
            constant = ((*that1).operation == -(1000 as c_int)
                && (*that2).operation == -(1000 as c_int)) as c_int;
            if (*that1).ntype != fits_parser_yytokentype::STRING as c_int
                && (*that1).ntype != fits_parser_yytokentype::BITSTR as c_int
                && Test_Dims(lParse, Node1, Node2) == 0
            {
                Free_Last_Node(lParse);
                fits_parser_yyerror(
                    0 as yyscan_t,
                    lParse,
                    b"Array sizes/dims do not match for binary operator\0" as *const u8
                        as *const c_char as *mut c_char,
                );
                return -(1);
            }
            if (*that1).value.nelem == 1 {
                that1 = that2;
            }
            (*this).value.nelem = (*that1).value.nelem;
            (*this).value.naxis = (*that1).value.naxis;
            i = 0;
            while i < (*that1).value.naxis {
                (*this).value.naxes[i as usize] = (*that1).value.naxes[i as usize];
                i += 1;
                i;
            }
            if Op == fits_parser_yytokentype::ACCUM as c_int
                && (*that1).ntype == fits_parser_yytokentype::BITSTR as c_int
            {
                (*this).value.nelem = 1;
                (*this).value.naxis = 1;
                (*this).value.naxes[0] = 1;
            }
            match (*that1).ntype {
                262 => {
                    (*this).DoOp = Some(Do_BinOp_bit);
                }
                261 => {
                    (*this).DoOp = Some(Do_BinOp_str);
                }
                258 => {
                    (*this).DoOp = Some(Do_BinOp_log);
                }
                259 => {
                    (*this).DoOp = Some(Do_BinOp_lng);
                }
                260 => {
                    (*this).DoOp = Some(Do_BinOp_dbl);
                }
                _ => {}
            }
            if constant != 0 {
                ((*this).DoOp).expect("non-null function pointer")(lParse, this);
            }
        }
        n
    }
}
unsafe fn New_Func(
    mut lParse: *mut ParseData,
    mut returnType: c_int,
    mut Op: funcOp,
    mut nNodes: c_int,
    mut Node1: c_int,
    mut Node2: c_int,
    mut Node3: c_int,
    mut Node4: c_int,
    mut Node5: c_int,
    mut Node6: c_int,
    mut Node7: c_int,
) -> c_int {
    unsafe {
        New_FuncSize(
            lParse, returnType, Op, nNodes, Node1, Node2, Node3, Node4, Node5, Node6, Node7, 0,
        )
    }
}
unsafe fn yydestruct(
    mut yymsg: *const c_char,
    mut yykind: yysymbol_kind_t,
    mut yyvaluep: *mut FITS_PARSER_YYSTYPE,
    mut scanner: yyscan_t,
    mut lParse: *mut ParseData,
) {
    if yymsg.is_null() {
        yymsg = b"Deleting\0" as *const u8 as *const c_char;
    }
}
unsafe fn New_FuncSize(
    mut lParse: *mut ParseData,
    mut returnType: c_int,
    mut Op: funcOp,
    mut nNodes: c_int,
    mut Node1: c_int,
    mut Node2: c_int,
    mut Node3: c_int,
    mut Node4: c_int,
    mut Node5: c_int,
    mut Node6: c_int,
    mut Node7: c_int,
    mut Size: c_int,
) -> c_int {
    unsafe {
        let mut this: *mut Node = std::ptr::null_mut::<Node>();
        let mut that: *mut Node = std::ptr::null_mut::<Node>();
        let mut i: c_int = 0;
        let mut n: c_int = 0;
        let mut constant: c_int = 0;
        if Node1 < 0 || Node2 < 0 || Node3 < 0 || Node4 < 0 || Node5 < 0 || Node6 < 0 || Node7 < 0 {
            return -(1);
        }
        n = Alloc_Node(lParse);
        if n >= 0 {
            this = ((*lParse).Nodes).offset(n as isize);
            (*this).operation = Op as c_int;
            (*this).DoOp = Some(Do_Func);
            (*this).nSubNodes = nNodes;
            (*this).SubNodes[0] = Node1;
            (*this).SubNodes[1] = Node2;
            (*this).SubNodes[2] = Node3;
            (*this).SubNodes[3] = Node4;
            (*this).SubNodes[4] = Node5;
            (*this).SubNodes[5] = Node6;
            (*this).SubNodes[6] = Node7;
            constant = nNodes;
            i = constant;
            if Op as c_uint == poirnd_fct as c_int as c_uint {
                constant = 0;
            }
            loop {
                let fresh1 = i;
                i -= 1;
                if fresh1 == 0 {
                    break;
                }
                constant = (constant != 0
                    && (*((*lParse).Nodes).offset((*this).SubNodes[i as usize] as isize)).operation
                        == -(1000 as c_int)) as c_int;
            }
            if returnType != 0 {
                (*this).ntype = returnType;
                (*this).value.nelem = 1;
                (*this).value.naxis = 1;
                (*this).value.naxes[0] = 1;
            } else {
                that = ((*lParse).Nodes).offset(Node1 as isize);
                (*this).ntype = (*that).ntype;
                (*this).value.nelem = (*that).value.nelem;
                (*this).value.naxis = (*that).value.naxis;
                i = 0;
                while i < (*that).value.naxis {
                    (*this).value.naxes[i as usize] = (*that).value.naxes[i as usize];
                    i += 1;
                    i;
                }
            }
            if Size > 0 {
                (*this).value.nelem = Size as c_long;
            }
            if constant != 0 {
                ((*this).DoOp).expect("non-null function pointer")(lParse, this);
            }
        }
        n
    }
}

pub unsafe fn fits_parser_yyparse(mut scanner: yyscan_t, mut lParse: *mut ParseData) -> c_int {
    unsafe {
        let mut current_block: u64;
        let mut yychar: c_int = 0;
        static mut yyval_default: FITS_PARSER_YYSTYPE = FITS_PARSER_YYSTYPE { Node: 0 };
        let mut yylval: FITS_PARSER_YYSTYPE = yyval_default;
        let mut fits_parser_yynerrs: c_int = 0;
        let mut yystate: yy_state_fast_t = 0;
        let mut yyerrstatus: c_int = 0;
        let mut yystacksize: c_long = 100 as c_long;
        let mut yyssa: [yy_state_t; 100] = [0; 100];
        let mut yyss: *mut yy_state_t = yyssa.as_mut_ptr();
        let mut yyssp: *mut yy_state_t = yyss;
        let mut yyvsa: [FITS_PARSER_YYSTYPE; 100] = [FITS_PARSER_YYSTYPE { Node: 0 }; 100];
        let mut yyvs: *mut FITS_PARSER_YYSTYPE = yyvsa.as_mut_ptr();
        let mut yyvsp: *mut FITS_PARSER_YYSTYPE = yyvs;
        let mut yyn: c_int = 0;
        let mut yyresult: c_int = 0;
        let mut yytoken: yysymbol_kind_t = YYSYMBOL_YYEMPTY;
        let mut yyval: FITS_PARSER_YYSTYPE = FITS_PARSER_YYSTYPE { Node: 0 };
        let mut yylen: c_int = 0;
        yychar = fits_parser_yytokentype::FITS_PARSER_YYEMPTY as c_int;
        's_54: loop {
            (0 as c_int != 0 && (0 as c_int <= yystate && yystate < 322 as c_int)) as c_int;
            *yyssp = yystate as yy_state_t;
            if yyss
                .offset(yystacksize as isize)
                .offset(-(1 as c_int as isize))
                <= yyssp
            {
                let mut yysize: c_long = yyssp.offset_from(yyss) as c_long + 1;
                if 10000 as c_long <= yystacksize {
                    current_block = 11794367917084412820;
                    break;
                }
                yystacksize *= 2 as c_long;
                if (10000 as c_long) < yystacksize {
                    yystacksize = 10000 as c_long;
                }
                let mut yyss1: *mut yy_state_t = yyss;
                let mut yyptr: *mut yyalloc = malloc(
                    ((yystacksize
                        * (::core::mem::size_of::<yy_state_t>() as c_ulong as c_long
                            + ::core::mem::size_of::<FITS_PARSER_YYSTYPE>() as c_ulong as c_long)
                        + (::core::mem::size_of::<yyalloc>() as c_ulong as c_long - 1))
                        as c_ulong)
                        .try_into()
                        .unwrap(),
                ) as *mut yyalloc;
                if yyptr.is_null() {
                    current_block = 11794367917084412820;
                    break;
                }
                let mut yynewbytes: c_long = 0;
                libc::memcpy(
                    &mut (*yyptr).yyss_alloc as *mut yy_state_t as *mut c_void,
                    yyss as *const c_void,
                    (yysize as c_ulong)
                        .wrapping_mul(::core::mem::size_of::<yy_state_t>() as c_ulong)
                        as libc::size_t,
                );
                yyss = &mut (*yyptr).yyss_alloc;
                yynewbytes = yystacksize
                    * ::core::mem::size_of::<yy_state_t>() as c_ulong as c_long
                    + (::core::mem::size_of::<yyalloc>() as c_ulong as c_long - 1);
                yyptr = yyptr.offset(
                    (yynewbytes / ::core::mem::size_of::<yyalloc>() as c_ulong as c_long) as isize,
                );
                let mut yynewbytes_0: c_long = 0;
                libc::memcpy(
                    &mut (*yyptr).yyvs_alloc as *mut FITS_PARSER_YYSTYPE as *mut c_void,
                    yyvs as *const c_void,
                    (yysize as c_ulong)
                        .wrapping_mul(::core::mem::size_of::<FITS_PARSER_YYSTYPE>() as c_ulong)
                        as libc::size_t,
                );
                yyvs = &mut (*yyptr).yyvs_alloc;
                yynewbytes_0 = yystacksize
                    * ::core::mem::size_of::<FITS_PARSER_YYSTYPE>() as c_ulong as c_long
                    + (::core::mem::size_of::<yyalloc>() as c_ulong as c_long - 1);
                yyptr = yyptr.offset(
                    (yynewbytes_0 / ::core::mem::size_of::<yyalloc>() as c_ulong as c_long)
                        as isize,
                );
                if yyss1 != yyssa.as_mut_ptr() {
                    free(yyss1 as *mut c_void);
                }
                yyssp = yyss.offset(yysize as isize).offset(-(1 as c_int as isize));
                yyvsp = yyvs.offset(yysize as isize).offset(-(1 as c_int as isize));
                if yyss
                    .offset(yystacksize as isize)
                    .offset(-(1 as c_int as isize))
                    <= yyssp
                {
                    current_block = 3964311021479492664;
                    break;
                }
            }
            if yystate == 2 as c_int {
                yyresult = 0;
                current_block = 15864720325503947191;
                break;
            } else {
                yyn = yypact[yystate as usize] as c_int;
                if yyn == -(41 as c_int) {
                    current_block = 5937473999264333383;
                } else {
                    if yychar == fits_parser_yytokentype::FITS_PARSER_YYEMPTY as c_int {
                        yychar =
                            fits_parser_yylex(&mut yylval as *mut FITS_PARSER_YYSTYPE, scanner);
                    }
                    if yychar <= fits_parser_yytokentype::FITS_PARSER_YYEOF as c_int {
                        yychar = fits_parser_yytokentype::FITS_PARSER_YYEOF as c_int;
                        yytoken = YYSYMBOL_YYEOF;
                        current_block = 1924505913685386279;
                    } else if yychar == fits_parser_yytokentype::FITS_PARSER_YYerror as c_int {
                        yychar = fits_parser_yytokentype::FITS_PARSER_YYUNDEF as c_int;
                        yytoken = YYSYMBOL_YYerror;
                        current_block = 1774893048582444437;
                    } else {
                        yytoken = (if 0 <= yychar && yychar <= 292 as c_int {
                            yytranslate[yychar as usize] as yysymbol_kind_t as c_int
                        } else {
                            YYSYMBOL_YYUNDEF as c_int
                        }) as yysymbol_kind_t;
                        current_block = 1924505913685386279;
                    }
                    match current_block {
                        1774893048582444437 => {}
                        _ => {
                            yyn += yytoken as c_int;
                            if yyn < 0
                                || (1776 as c_int) < yyn
                                || yycheck[yyn as usize] as c_int != yytoken as c_int
                            {
                                current_block = 5937473999264333383;
                            } else {
                                yyn = yytable[yyn as usize] as c_int;
                                if yyn <= 0 {
                                    yyn = -yyn;
                                    current_block = 670225253387957849;
                                } else {
                                    if yyerrstatus != 0 {
                                        yyerrstatus -= 1;
                                        yyerrstatus;
                                    }
                                    yystate = yyn;
                                    yyvsp = yyvsp.offset(1);
                                    *yyvsp = yylval;
                                    yychar = fits_parser_yytokentype::FITS_PARSER_YYEMPTY as c_int;
                                    current_block = 7872030484262409139;
                                }
                            }
                        }
                    }
                }
                if current_block == 5937473999264333383 {
                    yyn = yydefact[yystate as usize] as c_int;
                    if yyn == 0 {
                        yytoken =
                            (if yychar == fits_parser_yytokentype::FITS_PARSER_YYEMPTY as c_int {
                                YYSYMBOL_YYEMPTY as c_int
                            } else if 0 <= yychar && yychar <= 292 as c_int {
                                yytranslate[yychar as usize] as yysymbol_kind_t as c_int
                            } else {
                                YYSYMBOL_YYUNDEF as c_int
                            }) as yysymbol_kind_t;
                        if yyerrstatus == 0 {
                            fits_parser_yynerrs += 1;
                            fits_parser_yynerrs;
                            fits_parser_yyerror(
                                scanner,
                                lParse,
                                b"syntax error\0" as *const u8 as *const c_char as *mut c_char,
                            );
                        }
                        if yyerrstatus == 3 as c_int {
                            if yychar <= fits_parser_yytokentype::FITS_PARSER_YYEOF as c_int {
                                if yychar == fits_parser_yytokentype::FITS_PARSER_YYEOF as c_int {
                                    current_block = 3964311021479492664;
                                    break;
                                }
                            } else {
                                yydestruct(
                                    b"Error: discarding\0" as *const u8 as *const c_char,
                                    yytoken,
                                    &mut yylval,
                                    scanner,
                                    lParse,
                                );
                                yychar = fits_parser_yytokentype::FITS_PARSER_YYEMPTY as c_int;
                            }
                        }
                        current_block = 1774893048582444437;
                    } else {
                        current_block = 670225253387957849;
                    }
                }
                if current_block == 670225253387957849 {
                    yylen = yyr2[yyn as usize] as c_int;
                    yyval = *yyvsp.offset((1 as c_int - yylen) as isize);
                    match yyn {
                        4 => {
                            current_block = 17353983478346836848;
                        }
                        5 => {
                            if (*yyvsp.offset(-1_isize)).Node < 0 {
                                fits_parser_yyerror(
                                    scanner,
                                    lParse,
                                    b"Couldn't build node structure: out of memory?\0" as *const u8
                                        as *const c_char
                                        as *mut c_char,
                                );
                                current_block = 4830776507462815627;
                            } else {
                                (*lParse).resultNode = (*yyvsp.offset(-1_isize)).Node;
                                current_block = 17353983478346836848;
                            }
                        }
                        6 => {
                            if (*yyvsp.offset(-1_isize)).Node < 0 {
                                fits_parser_yyerror(
                                    scanner,
                                    lParse,
                                    b"Couldn't build node structure: out of memory?\0" as *const u8
                                        as *const c_char
                                        as *mut c_char,
                                );
                                current_block = 4830776507462815627;
                            } else {
                                (*lParse).resultNode = (*yyvsp.offset(-1_isize)).Node;
                                current_block = 17353983478346836848;
                            }
                        }
                        7 => {
                            if (*yyvsp.offset(-1_isize)).Node < 0 {
                                fits_parser_yyerror(
                                    scanner,
                                    lParse,
                                    b"Couldn't build node structure: out of memory?\0" as *const u8
                                        as *const c_char
                                        as *mut c_char,
                                );
                                current_block = 4830776507462815627;
                            } else {
                                (*lParse).resultNode = (*yyvsp.offset(-1_isize)).Node;
                                current_block = 17353983478346836848;
                            }
                        }
                        8 => {
                            if (*yyvsp.offset(-1_isize)).Node < 0 {
                                fits_parser_yyerror(
                                    scanner,
                                    lParse,
                                    b"Couldn't build node structure: out of memory?\0" as *const u8
                                        as *const c_char
                                        as *mut c_char,
                                );
                                current_block = 4830776507462815627;
                            } else {
                                (*lParse).resultNode = (*yyvsp.offset(-1_isize)).Node;
                                current_block = 17353983478346836848;
                            }
                        }
                        9 => {
                            yyerrstatus = 0;
                            current_block = 17353983478346836848;
                        }
                        10 => {
                            yyval.Node = New_Vector(lParse, (*yyvsp.offset(0)).Node);
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        11 => {
                            if (*((*lParse).Nodes).offset((*yyvsp.offset(-2_isize)).Node as isize))
                                .nSubNodes
                                >= 10 as c_int
                            {
                                (*yyvsp.offset(-2_isize)).Node =
                                    Close_Vec(lParse, (*yyvsp.offset(-2_isize)).Node);
                                if (*yyvsp.offset(-2_isize)).Node < 0 {
                                    current_block = 4830776507462815627;
                                } else {
                                    yyval.Node = New_Vector(lParse, (*yyvsp.offset(-2_isize)).Node);
                                    if yyval.Node < 0 {
                                        current_block = 4830776507462815627;
                                    } else {
                                        current_block = 2616667235040759262;
                                    }
                                }
                            } else {
                                yyval.Node = (*yyvsp.offset(-2_isize)).Node;
                                current_block = 2616667235040759262;
                            }
                            match current_block {
                                4830776507462815627 => {}
                                _ => {
                                    let fresh2 = &mut (*((*lParse).Nodes)
                                        .offset(yyval.Node as isize))
                                    .nSubNodes;
                                    let fresh3 = *fresh2;
                                    *fresh2 += 1;
                                    (*((*lParse).Nodes).offset(yyval.Node as isize)).SubNodes
                                        [fresh3 as usize] = (*yyvsp.offset(0)).Node;
                                    current_block = 17353983478346836848;
                                }
                            }
                        }
                        12 => {
                            yyval.Node = New_Vector(lParse, (*yyvsp.offset(0)).Node);
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        13 => {
                            if (*((*lParse).Nodes).offset((*yyvsp.offset(-2_isize)).Node as isize))
                                .ntype
                                < (*((*lParse).Nodes).offset((*yyvsp.offset(0)).Node as isize))
                                    .ntype
                            {
                                (*((*lParse).Nodes)
                                    .offset((*yyvsp.offset(-2_isize)).Node as isize))
                                .ntype = (*((*lParse).Nodes)
                                    .offset((*yyvsp.offset(0)).Node as isize))
                                .ntype;
                            }
                            if (*((*lParse).Nodes).offset((*yyvsp.offset(-2_isize)).Node as isize))
                                .nSubNodes
                                >= 10 as c_int
                            {
                                (*yyvsp.offset(-2_isize)).Node =
                                    Close_Vec(lParse, (*yyvsp.offset(-2_isize)).Node);
                                if (*yyvsp.offset(-2_isize)).Node < 0 {
                                    current_block = 4830776507462815627;
                                } else {
                                    yyval.Node = New_Vector(lParse, (*yyvsp.offset(-2_isize)).Node);
                                    if yyval.Node < 0 {
                                        current_block = 4830776507462815627;
                                    } else {
                                        current_block = 2904036176499606090;
                                    }
                                }
                            } else {
                                yyval.Node = (*yyvsp.offset(-2_isize)).Node;
                                current_block = 2904036176499606090;
                            }
                            match current_block {
                                4830776507462815627 => {}
                                _ => {
                                    let fresh4 = &mut (*((*lParse).Nodes)
                                        .offset(yyval.Node as isize))
                                    .nSubNodes;
                                    let fresh5 = *fresh4;
                                    *fresh4 += 1;
                                    (*((*lParse).Nodes).offset(yyval.Node as isize)).SubNodes
                                        [fresh5 as usize] = (*yyvsp.offset(0)).Node;
                                    current_block = 17353983478346836848;
                                }
                            }
                        }
                        14 => {
                            if (*((*lParse).Nodes).offset((*yyvsp.offset(-2_isize)).Node as isize))
                                .nSubNodes
                                >= 10 as c_int
                            {
                                (*yyvsp.offset(-2_isize)).Node =
                                    Close_Vec(lParse, (*yyvsp.offset(-2_isize)).Node);
                                if (*yyvsp.offset(-2_isize)).Node < 0 {
                                    current_block = 4830776507462815627;
                                } else {
                                    yyval.Node = New_Vector(lParse, (*yyvsp.offset(-2_isize)).Node);
                                    if yyval.Node < 0 {
                                        current_block = 4830776507462815627;
                                    } else {
                                        current_block = 17702298541784679949;
                                    }
                                }
                            } else {
                                yyval.Node = (*yyvsp.offset(-2_isize)).Node;
                                current_block = 17702298541784679949;
                            }
                            match current_block {
                                4830776507462815627 => {}
                                _ => {
                                    let fresh6 = &mut (*((*lParse).Nodes)
                                        .offset(yyval.Node as isize))
                                    .nSubNodes;
                                    let fresh7 = *fresh6;
                                    *fresh6 += 1;
                                    (*((*lParse).Nodes).offset(yyval.Node as isize)).SubNodes
                                        [fresh7 as usize] = (*yyvsp.offset(0)).Node;
                                    current_block = 17353983478346836848;
                                }
                            }
                        }
                        15 => {
                            (*((*lParse).Nodes).offset((*yyvsp.offset(-2_isize)).Node as isize))
                                .ntype =
                                (*((*lParse).Nodes).offset((*yyvsp.offset(0)).Node as isize)).ntype;
                            if (*((*lParse).Nodes).offset((*yyvsp.offset(-2_isize)).Node as isize))
                                .nSubNodes
                                >= 10 as c_int
                            {
                                (*yyvsp.offset(-2_isize)).Node =
                                    Close_Vec(lParse, (*yyvsp.offset(-2_isize)).Node);
                                if (*yyvsp.offset(-2_isize)).Node < 0 {
                                    current_block = 4830776507462815627;
                                } else {
                                    yyval.Node = New_Vector(lParse, (*yyvsp.offset(-2_isize)).Node);
                                    if yyval.Node < 0 {
                                        current_block = 4830776507462815627;
                                    } else {
                                        current_block = 1069630499025798221;
                                    }
                                }
                            } else {
                                yyval.Node = (*yyvsp.offset(-2_isize)).Node;
                                current_block = 1069630499025798221;
                            }
                            match current_block {
                                4830776507462815627 => {}
                                _ => {
                                    let fresh8 = &mut (*((*lParse).Nodes)
                                        .offset(yyval.Node as isize))
                                    .nSubNodes;
                                    let fresh9 = *fresh8;
                                    *fresh8 += 1;
                                    (*((*lParse).Nodes).offset(yyval.Node as isize)).SubNodes
                                        [fresh9 as usize] = (*yyvsp.offset(0)).Node;
                                    current_block = 17353983478346836848;
                                }
                            }
                        }
                        16 => {
                            yyval.Node = Close_Vec(lParse, (*yyvsp.offset(-1_isize)).Node);
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        17 => {
                            yyval.Node = Close_Vec(lParse, (*yyvsp.offset(-1_isize)).Node);
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        18 => {
                            yyval.Node = New_Const(
                                lParse,
                                fits_parser_yytokentype::BITSTR as c_int,
                                ((*yyvsp.offset(0)).astr).as_mut_ptr() as *mut c_void,
                                (strlen(((*yyvsp.offset(0)).astr).as_mut_ptr()))
                                    .wrapping_add((1 as c_int as c_ulong).try_into().unwrap())
                                    as c_long,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                ((*((*lParse).Nodes).offset(yyval.Node as isize)).value).nelem =
                                    strlen(((*yyvsp.offset(0)).astr).as_mut_ptr()) as c_long;
                                current_block = 17353983478346836848;
                            }
                        }
                        19 => {
                            yyval.Node = New_Column(lParse, (*yyvsp.offset(0)).lng as c_int);
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        20 => {
                            if (*((*lParse).Nodes).offset((*yyvsp.offset(-1_isize)).Node as isize))
                                .ntype
                                != fits_parser_yytokentype::LONG as c_int
                                || (*((*lParse).Nodes)
                                    .offset((*yyvsp.offset(-1_isize)).Node as isize))
                                .operation
                                    != -(1000 as c_int)
                            {
                                fits_parser_yyerror(
                                    scanner,
                                    lParse,
                                    b"Offset argument must be a constant integer\0" as *const u8
                                        as *const c_char
                                        as *mut c_char,
                                );
                                current_block = 4830776507462815627;
                            } else {
                                yyval.Node = New_Offset(
                                    lParse,
                                    (*yyvsp.offset(-3_isize)).lng as c_int,
                                    (*yyvsp.offset(-1_isize)).Node,
                                );
                                if yyval.Node < 0 {
                                    current_block = 4830776507462815627;
                                } else {
                                    current_block = 17353983478346836848;
                                }
                            }
                        }
                        21 => {
                            yyval.Node = New_BinOp(
                                lParse,
                                fits_parser_yytokentype::BITSTR as c_int,
                                (*yyvsp.offset(-2_isize)).Node,
                                '&' as i32,
                                (*yyvsp.offset(0)).Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                ((*((*lParse).Nodes).offset(yyval.Node as isize)).value).nelem =
                                    if ((*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-2_isize)).Node as isize))
                                    .value)
                                        .nelem
                                        > (*((*lParse).Nodes)
                                            .offset((*yyvsp.offset(0)).Node as isize))
                                        .value
                                        .nelem
                                    {
                                        (*((*lParse).Nodes)
                                            .offset((*yyvsp.offset(-2_isize)).Node as isize))
                                        .value
                                        .nelem
                                    } else {
                                        (*((*lParse).Nodes)
                                            .offset((*yyvsp.offset(0)).Node as isize))
                                        .value
                                        .nelem
                                    };
                                current_block = 17353983478346836848;
                            }
                        }
                        22 => {
                            yyval.Node = New_BinOp(
                                lParse,
                                fits_parser_yytokentype::BITSTR as c_int,
                                (*yyvsp.offset(-2_isize)).Node,
                                '|' as i32,
                                (*yyvsp.offset(0)).Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                ((*((*lParse).Nodes).offset(yyval.Node as isize)).value).nelem =
                                    if ((*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-2_isize)).Node as isize))
                                    .value)
                                        .nelem
                                        > (*((*lParse).Nodes)
                                            .offset((*yyvsp.offset(0)).Node as isize))
                                        .value
                                        .nelem
                                    {
                                        (*((*lParse).Nodes)
                                            .offset((*yyvsp.offset(-2_isize)).Node as isize))
                                        .value
                                        .nelem
                                    } else {
                                        (*((*lParse).Nodes)
                                            .offset((*yyvsp.offset(0)).Node as isize))
                                        .value
                                        .nelem
                                    };
                                current_block = 17353983478346836848;
                            }
                        }
                        23 => {
                            if (*((*lParse).Nodes).offset((*yyvsp.offset(-2_isize)).Node as isize))
                                .value
                                .nelem
                                + (*((*lParse).Nodes).offset((*yyvsp.offset(0)).Node as isize))
                                    .value
                                    .nelem
                                >= 256 as c_long
                            {
                                fits_parser_yyerror(
                                    scanner,
                                    lParse,
                                    b"Combined bit string size exceeds 255 bits\0" as *const u8
                                        as *const c_char
                                        as *mut c_char,
                                );
                                current_block = 4830776507462815627;
                            } else {
                                yyval.Node = New_BinOp(
                                    lParse,
                                    fits_parser_yytokentype::BITSTR as c_int,
                                    (*yyvsp.offset(-2_isize)).Node,
                                    '+' as i32,
                                    (*yyvsp.offset(0)).Node,
                                );
                                if yyval.Node < 0 {
                                    current_block = 4830776507462815627;
                                } else {
                                    ((*((*lParse).Nodes).offset(yyval.Node as isize)).value)
                                        .nelem = (*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-2_isize)).Node as isize))
                                    .value
                                    .nelem
                                        + (*((*lParse).Nodes)
                                            .offset((*yyvsp.offset(0)).Node as isize))
                                        .value
                                        .nelem;
                                    current_block = 17353983478346836848;
                                }
                            }
                        }
                        24 => {
                            yyval.Node = New_Deref(
                                lParse,
                                (*yyvsp.offset(-3_isize)).Node,
                                1,
                                (*yyvsp.offset(-1_isize)).Node,
                                0,
                                0,
                                0,
                                0,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        25 => {
                            yyval.Node = New_Deref(
                                lParse,
                                (*yyvsp.offset(-(5))).Node,
                                2 as c_int,
                                (*yyvsp.offset(-3_isize)).Node,
                                (*yyvsp.offset(-1_isize)).Node,
                                0,
                                0,
                                0,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        26 => {
                            yyval.Node = New_Deref(
                                lParse,
                                (*yyvsp.offset(-(7 as c_int) as isize)).Node,
                                3 as c_int,
                                (*yyvsp.offset(-(5))).Node,
                                (*yyvsp.offset(-3_isize)).Node,
                                (*yyvsp.offset(-1_isize)).Node,
                                0,
                                0,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        27 => {
                            yyval.Node = New_Deref(
                                lParse,
                                (*yyvsp.offset(-(9 as c_int) as isize)).Node,
                                4 as c_int,
                                (*yyvsp.offset(-(7 as c_int) as isize)).Node,
                                (*yyvsp.offset(-(5))).Node,
                                (*yyvsp.offset(-3_isize)).Node,
                                (*yyvsp.offset(-1_isize)).Node,
                                0,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        28 => {
                            yyval.Node = New_Deref(
                                lParse,
                                (*yyvsp.offset(-(11 as c_int) as isize)).Node,
                                5 as c_int,
                                (*yyvsp.offset(-(9 as c_int) as isize)).Node,
                                (*yyvsp.offset(-(7 as c_int) as isize)).Node,
                                (*yyvsp.offset(-(5))).Node,
                                (*yyvsp.offset(-3_isize)).Node,
                                (*yyvsp.offset(-1_isize)).Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        29 => {
                            yyval.Node = New_Unary(
                                lParse,
                                fits_parser_yytokentype::BITSTR as c_int,
                                fits_parser_yytokentype::NOT as c_int,
                                (*yyvsp.offset(0)).Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        30 => {
                            yyval.Node = (*yyvsp.offset(-1_isize)).Node;
                            current_block = 17353983478346836848;
                        }
                        31 => {
                            yyval.Node = New_Const(
                                lParse,
                                fits_parser_yytokentype::LONG as c_int,
                                &mut (*yyvsp.offset(0)).lng as *mut c_long as *mut c_void,
                                ::core::mem::size_of::<c_long>() as c_ulong as c_long,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        32 => {
                            yyval.Node = New_Const(
                                lParse,
                                fits_parser_yytokentype::DOUBLE as c_int,
                                &mut (*yyvsp.offset(0)).dbl as *mut c_double as *mut c_void,
                                ::core::mem::size_of::<c_double>() as c_ulong as c_long,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        33 => {
                            yyval.Node = New_Column(lParse, (*yyvsp.offset(0)).lng as c_int);
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        34 => {
                            if (*((*lParse).Nodes).offset((*yyvsp.offset(-1_isize)).Node as isize))
                                .ntype
                                != fits_parser_yytokentype::LONG as c_int
                                || (*((*lParse).Nodes)
                                    .offset((*yyvsp.offset(-1_isize)).Node as isize))
                                .operation
                                    != -(1000 as c_int)
                            {
                                fits_parser_yyerror(
                                    scanner,
                                    lParse,
                                    b"Offset argument must be a constant integer\0" as *const u8
                                        as *const c_char
                                        as *mut c_char,
                                );
                                current_block = 4830776507462815627;
                            } else {
                                yyval.Node = New_Offset(
                                    lParse,
                                    (*yyvsp.offset(-3_isize)).lng as c_int,
                                    (*yyvsp.offset(-1_isize)).Node,
                                );
                                if yyval.Node < 0 {
                                    current_block = 4830776507462815627;
                                } else {
                                    current_block = 17353983478346836848;
                                }
                            }
                        }
                        35 => {
                            yyval.Node = New_Func(
                                lParse,
                                fits_parser_yytokentype::LONG as c_int,
                                row_fct,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                            );
                            current_block = 17353983478346836848;
                        }
                        36 => {
                            yyval.Node = New_Func(
                                lParse,
                                fits_parser_yytokentype::LONG as c_int,
                                null_fct,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                            );
                            current_block = 17353983478346836848;
                        }
                        37 => {
                            if (*((*lParse).Nodes).offset((*yyvsp.offset(-2_isize)).Node as isize))
                                .ntype
                                > (*((*lParse).Nodes).offset((*yyvsp.offset(0)).Node as isize))
                                    .ntype
                            {
                                (*yyvsp.offset(0)).Node = New_Unary(
                                    lParse,
                                    (*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-2_isize)).Node as isize))
                                    .ntype,
                                    0,
                                    (*yyvsp.offset(0)).Node,
                                );
                            } else if (*((*lParse).Nodes)
                                .offset((*yyvsp.offset(-2_isize)).Node as isize))
                            .ntype
                                < (*((*lParse).Nodes).offset((*yyvsp.offset(0)).Node as isize))
                                    .ntype
                            {
                                (*yyvsp.offset(-2_isize)).Node = New_Unary(
                                    lParse,
                                    (*((*lParse).Nodes).offset((*yyvsp.offset(0)).Node as isize))
                                        .ntype,
                                    0,
                                    (*yyvsp.offset(-2_isize)).Node,
                                );
                            }
                            yyval.Node = New_BinOp(
                                lParse,
                                (*((*lParse).Nodes)
                                    .offset((*yyvsp.offset(-2_isize)).Node as isize))
                                .ntype,
                                (*yyvsp.offset(-2_isize)).Node,
                                '%' as i32,
                                (*yyvsp.offset(0)).Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        38 => {
                            if (*((*lParse).Nodes).offset((*yyvsp.offset(-2_isize)).Node as isize))
                                .ntype
                                > (*((*lParse).Nodes).offset((*yyvsp.offset(0)).Node as isize))
                                    .ntype
                            {
                                (*yyvsp.offset(0)).Node = New_Unary(
                                    lParse,
                                    (*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-2_isize)).Node as isize))
                                    .ntype,
                                    0,
                                    (*yyvsp.offset(0)).Node,
                                );
                            } else if (*((*lParse).Nodes)
                                .offset((*yyvsp.offset(-2_isize)).Node as isize))
                            .ntype
                                < (*((*lParse).Nodes).offset((*yyvsp.offset(0)).Node as isize))
                                    .ntype
                            {
                                (*yyvsp.offset(-2_isize)).Node = New_Unary(
                                    lParse,
                                    (*((*lParse).Nodes).offset((*yyvsp.offset(0)).Node as isize))
                                        .ntype,
                                    0,
                                    (*yyvsp.offset(-2_isize)).Node,
                                );
                            }
                            yyval.Node = New_BinOp(
                                lParse,
                                (*((*lParse).Nodes)
                                    .offset((*yyvsp.offset(-2_isize)).Node as isize))
                                .ntype,
                                (*yyvsp.offset(-2_isize)).Node,
                                '+' as i32,
                                (*yyvsp.offset(0)).Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        39 => {
                            if (*((*lParse).Nodes).offset((*yyvsp.offset(-2_isize)).Node as isize))
                                .ntype
                                > (*((*lParse).Nodes).offset((*yyvsp.offset(0)).Node as isize))
                                    .ntype
                            {
                                (*yyvsp.offset(0)).Node = New_Unary(
                                    lParse,
                                    (*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-2_isize)).Node as isize))
                                    .ntype,
                                    0,
                                    (*yyvsp.offset(0)).Node,
                                );
                            } else if (*((*lParse).Nodes)
                                .offset((*yyvsp.offset(-2_isize)).Node as isize))
                            .ntype
                                < (*((*lParse).Nodes).offset((*yyvsp.offset(0)).Node as isize))
                                    .ntype
                            {
                                (*yyvsp.offset(-2_isize)).Node = New_Unary(
                                    lParse,
                                    (*((*lParse).Nodes).offset((*yyvsp.offset(0)).Node as isize))
                                        .ntype,
                                    0,
                                    (*yyvsp.offset(-2_isize)).Node,
                                );
                            }
                            yyval.Node = New_BinOp(
                                lParse,
                                (*((*lParse).Nodes)
                                    .offset((*yyvsp.offset(-2_isize)).Node as isize))
                                .ntype,
                                (*yyvsp.offset(-2_isize)).Node,
                                '-' as i32,
                                (*yyvsp.offset(0)).Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        40 => {
                            if (*((*lParse).Nodes).offset((*yyvsp.offset(-2_isize)).Node as isize))
                                .ntype
                                > (*((*lParse).Nodes).offset((*yyvsp.offset(0)).Node as isize))
                                    .ntype
                            {
                                (*yyvsp.offset(0)).Node = New_Unary(
                                    lParse,
                                    (*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-2_isize)).Node as isize))
                                    .ntype,
                                    0,
                                    (*yyvsp.offset(0)).Node,
                                );
                            } else if (*((*lParse).Nodes)
                                .offset((*yyvsp.offset(-2_isize)).Node as isize))
                            .ntype
                                < (*((*lParse).Nodes).offset((*yyvsp.offset(0)).Node as isize))
                                    .ntype
                            {
                                (*yyvsp.offset(-2_isize)).Node = New_Unary(
                                    lParse,
                                    (*((*lParse).Nodes).offset((*yyvsp.offset(0)).Node as isize))
                                        .ntype,
                                    0,
                                    (*yyvsp.offset(-2_isize)).Node,
                                );
                            }
                            yyval.Node = New_BinOp(
                                lParse,
                                (*((*lParse).Nodes)
                                    .offset((*yyvsp.offset(-2_isize)).Node as isize))
                                .ntype,
                                (*yyvsp.offset(-2_isize)).Node,
                                '*' as i32,
                                (*yyvsp.offset(0)).Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        41 => {
                            if (*((*lParse).Nodes).offset((*yyvsp.offset(-2_isize)).Node as isize))
                                .ntype
                                > (*((*lParse).Nodes).offset((*yyvsp.offset(0)).Node as isize))
                                    .ntype
                            {
                                (*yyvsp.offset(0)).Node = New_Unary(
                                    lParse,
                                    (*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-2_isize)).Node as isize))
                                    .ntype,
                                    0,
                                    (*yyvsp.offset(0)).Node,
                                );
                            } else if (*((*lParse).Nodes)
                                .offset((*yyvsp.offset(-2_isize)).Node as isize))
                            .ntype
                                < (*((*lParse).Nodes).offset((*yyvsp.offset(0)).Node as isize))
                                    .ntype
                            {
                                (*yyvsp.offset(-2_isize)).Node = New_Unary(
                                    lParse,
                                    (*((*lParse).Nodes).offset((*yyvsp.offset(0)).Node as isize))
                                        .ntype,
                                    0,
                                    (*yyvsp.offset(-2_isize)).Node,
                                );
                            }
                            yyval.Node = New_BinOp(
                                lParse,
                                (*((*lParse).Nodes)
                                    .offset((*yyvsp.offset(-2_isize)).Node as isize))
                                .ntype,
                                (*yyvsp.offset(-2_isize)).Node,
                                '/' as i32,
                                (*yyvsp.offset(0)).Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        42 => {
                            if (*((*lParse).Nodes).offset((*yyvsp.offset(-2_isize)).Node as isize))
                                .ntype
                                != fits_parser_yytokentype::LONG as c_int
                                || (*((*lParse).Nodes).offset((*yyvsp.offset(0)).Node as isize))
                                    .ntype
                                    != fits_parser_yytokentype::LONG as c_int
                            {
                                fits_parser_yyerror(
                                scanner,
                                lParse,
                                b"Bitwise operations with incompatible types; only (bit OP bit) and (int OP int) are allowed\0"
                                    as *const u8 as *const c_char
                                    as *mut c_char,
                            );
                                current_block = 4830776507462815627;
                            } else {
                                yyval.Node = New_BinOp(
                                    lParse,
                                    (*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-2_isize)).Node as isize))
                                    .ntype,
                                    (*yyvsp.offset(-2_isize)).Node,
                                    '&' as i32,
                                    (*yyvsp.offset(0)).Node,
                                );
                                current_block = 17353983478346836848;
                            }
                        }
                        43 => {
                            if (*((*lParse).Nodes).offset((*yyvsp.offset(-2_isize)).Node as isize))
                                .ntype
                                != fits_parser_yytokentype::LONG as c_int
                                || (*((*lParse).Nodes).offset((*yyvsp.offset(0)).Node as isize))
                                    .ntype
                                    != fits_parser_yytokentype::LONG as c_int
                            {
                                fits_parser_yyerror(
                                scanner,
                                lParse,
                                b"Bitwise operations with incompatible types; only (bit OP bit) and (int OP int) are allowed\0"
                                    as *const u8 as *const c_char
                                    as *mut c_char,
                            );
                                current_block = 4830776507462815627;
                            } else {
                                yyval.Node = New_BinOp(
                                    lParse,
                                    (*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-2_isize)).Node as isize))
                                    .ntype,
                                    (*yyvsp.offset(-2_isize)).Node,
                                    '|' as i32,
                                    (*yyvsp.offset(0)).Node,
                                );
                                current_block = 17353983478346836848;
                            }
                        }
                        44 => {
                            if (*((*lParse).Nodes).offset((*yyvsp.offset(-2_isize)).Node as isize))
                                .ntype
                                != fits_parser_yytokentype::LONG as c_int
                                || (*((*lParse).Nodes).offset((*yyvsp.offset(0)).Node as isize))
                                    .ntype
                                    != fits_parser_yytokentype::LONG as c_int
                            {
                                fits_parser_yyerror(
                                scanner,
                                lParse,
                                b"Bitwise operations with incompatible types; only (bit OP bit) and (int OP int) are allowed\0"
                                    as *const u8 as *const c_char
                                    as *mut c_char,
                            );
                                current_block = 4830776507462815627;
                            } else {
                                yyval.Node = New_BinOp(
                                    lParse,
                                    (*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-2_isize)).Node as isize))
                                    .ntype,
                                    (*yyvsp.offset(-2_isize)).Node,
                                    '^' as i32,
                                    (*yyvsp.offset(0)).Node,
                                );
                                current_block = 17353983478346836848;
                            }
                        }
                        45 => {
                            if (*((*lParse).Nodes).offset((*yyvsp.offset(-2_isize)).Node as isize))
                                .ntype
                                > (*((*lParse).Nodes).offset((*yyvsp.offset(0)).Node as isize))
                                    .ntype
                            {
                                (*yyvsp.offset(0)).Node = New_Unary(
                                    lParse,
                                    (*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-2_isize)).Node as isize))
                                    .ntype,
                                    0,
                                    (*yyvsp.offset(0)).Node,
                                );
                            } else if (*((*lParse).Nodes)
                                .offset((*yyvsp.offset(-2_isize)).Node as isize))
                            .ntype
                                < (*((*lParse).Nodes).offset((*yyvsp.offset(0)).Node as isize))
                                    .ntype
                            {
                                (*yyvsp.offset(-2_isize)).Node = New_Unary(
                                    lParse,
                                    (*((*lParse).Nodes).offset((*yyvsp.offset(0)).Node as isize))
                                        .ntype,
                                    0,
                                    (*yyvsp.offset(-2_isize)).Node,
                                );
                            }
                            yyval.Node = New_BinOp(
                                lParse,
                                (*((*lParse).Nodes)
                                    .offset((*yyvsp.offset(-2_isize)).Node as isize))
                                .ntype,
                                (*yyvsp.offset(-2_isize)).Node,
                                fits_parser_yytokentype::POWER as c_int,
                                (*yyvsp.offset(0)).Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        46 => {
                            yyval.Node = (*yyvsp.offset(0)).Node;
                            current_block = 17353983478346836848;
                        }
                        47 => {
                            yyval.Node = New_Unary(
                                lParse,
                                (*((*lParse).Nodes).offset((*yyvsp.offset(0)).Node as isize)).ntype,
                                fits_parser_yytokentype::UMINUS as c_int,
                                (*yyvsp.offset(0)).Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        48 => {
                            yyval.Node = (*yyvsp.offset(-1_isize)).Node;
                            current_block = 17353983478346836848;
                        }
                        49 => {
                            (*yyvsp.offset(0)).Node = New_Unary(
                                lParse,
                                (*((*lParse).Nodes)
                                    .offset((*yyvsp.offset(-2_isize)).Node as isize))
                                .ntype,
                                0,
                                (*yyvsp.offset(0)).Node,
                            );
                            yyval.Node = New_BinOp(
                                lParse,
                                (*((*lParse).Nodes)
                                    .offset((*yyvsp.offset(-2_isize)).Node as isize))
                                .ntype,
                                (*yyvsp.offset(-2_isize)).Node,
                                '*' as i32,
                                (*yyvsp.offset(0)).Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        50 => {
                            (*yyvsp.offset(-2_isize)).Node = New_Unary(
                                lParse,
                                (*((*lParse).Nodes).offset((*yyvsp.offset(0)).Node as isize)).ntype,
                                0,
                                (*yyvsp.offset(-2_isize)).Node,
                            );
                            yyval.Node = New_BinOp(
                                lParse,
                                (*((*lParse).Nodes).offset((*yyvsp.offset(0)).Node as isize)).ntype,
                                (*yyvsp.offset(-2_isize)).Node,
                                '*' as i32,
                                (*yyvsp.offset(0)).Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        51 => {
                            if (*((*lParse).Nodes).offset((*yyvsp.offset(-2_isize)).Node as isize))
                                .ntype
                                > (*((*lParse).Nodes).offset((*yyvsp.offset(0)).Node as isize))
                                    .ntype
                            {
                                (*yyvsp.offset(0)).Node = New_Unary(
                                    lParse,
                                    (*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-2_isize)).Node as isize))
                                    .ntype,
                                    0,
                                    (*yyvsp.offset(0)).Node,
                                );
                            } else if (*((*lParse).Nodes)
                                .offset((*yyvsp.offset(-2_isize)).Node as isize))
                            .ntype
                                < (*((*lParse).Nodes).offset((*yyvsp.offset(0)).Node as isize))
                                    .ntype
                            {
                                (*yyvsp.offset(-2_isize)).Node = New_Unary(
                                    lParse,
                                    (*((*lParse).Nodes).offset((*yyvsp.offset(0)).Node as isize))
                                        .ntype,
                                    0,
                                    (*yyvsp.offset(-2_isize)).Node,
                                );
                            }
                            if Test_Dims(
                                lParse,
                                (*yyvsp.offset(-2_isize)).Node,
                                (*yyvsp.offset(0)).Node,
                            ) == 0
                            {
                                fits_parser_yyerror(
                                    scanner,
                                    lParse,
                                    b"Incompatible dimensions in '?:' arguments\0" as *const u8
                                        as *const c_char
                                        as *mut c_char,
                                );
                                current_block = 4830776507462815627;
                            } else {
                                yyval.Node = New_Func(
                                    lParse,
                                    0,
                                    ifthenelse_fct,
                                    3 as c_int,
                                    (*yyvsp.offset(-2_isize)).Node,
                                    (*yyvsp.offset(0)).Node,
                                    (*yyvsp.offset(-4_isize)).Node,
                                    0,
                                    0,
                                    0,
                                    0,
                                );
                                if yyval.Node < 0 {
                                    current_block = 4830776507462815627;
                                } else {
                                    if (*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-2_isize)).Node as isize))
                                    .value
                                    .nelem
                                        < (*((*lParse).Nodes)
                                            .offset((*yyvsp.offset(0)).Node as isize))
                                        .value
                                        .nelem
                                    {
                                        Copy_Dims(lParse, yyval.Node, (*yyvsp.offset(0)).Node);
                                    }
                                    (*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-4_isize)).Node as isize))
                                    .ntype = (*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-2_isize)).Node as isize))
                                    .ntype;
                                    if Test_Dims(lParse, (*yyvsp.offset(-4_isize)).Node, yyval.Node)
                                        == 0
                                    {
                                        fits_parser_yyerror(
                                            scanner,
                                            lParse,
                                            b"Incompatible dimensions in '?:' condition\0"
                                                as *const u8
                                                as *const c_char
                                                as *mut c_char,
                                        );
                                        current_block = 4830776507462815627;
                                    } else {
                                        (*((*lParse).Nodes)
                                            .offset((*yyvsp.offset(-4_isize)).Node as isize))
                                        .ntype = fits_parser_yytokentype::BOOLEAN as c_int;
                                        if ((*((*lParse).Nodes).offset(yyval.Node as isize)).value)
                                            .nelem
                                            < (*((*lParse).Nodes)
                                                .offset((*yyvsp.offset(-4_isize)).Node as isize))
                                            .value
                                            .nelem
                                        {
                                            Copy_Dims(
                                                lParse,
                                                yyval.Node,
                                                (*yyvsp.offset(-4_isize)).Node,
                                            );
                                        }
                                        current_block = 17353983478346836848;
                                    }
                                }
                            }
                        }
                        52 => {
                            if (*((*lParse).Nodes).offset((*yyvsp.offset(-2_isize)).Node as isize))
                                .ntype
                                > (*((*lParse).Nodes).offset((*yyvsp.offset(0)).Node as isize))
                                    .ntype
                            {
                                (*yyvsp.offset(0)).Node = New_Unary(
                                    lParse,
                                    (*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-2_isize)).Node as isize))
                                    .ntype,
                                    0,
                                    (*yyvsp.offset(0)).Node,
                                );
                            } else if (*((*lParse).Nodes)
                                .offset((*yyvsp.offset(-2_isize)).Node as isize))
                            .ntype
                                < (*((*lParse).Nodes).offset((*yyvsp.offset(0)).Node as isize))
                                    .ntype
                            {
                                (*yyvsp.offset(-2_isize)).Node = New_Unary(
                                    lParse,
                                    (*((*lParse).Nodes).offset((*yyvsp.offset(0)).Node as isize))
                                        .ntype,
                                    0,
                                    (*yyvsp.offset(-2_isize)).Node,
                                );
                            }
                            if Test_Dims(
                                lParse,
                                (*yyvsp.offset(-2_isize)).Node,
                                (*yyvsp.offset(0)).Node,
                            ) == 0
                            {
                                fits_parser_yyerror(
                                    scanner,
                                    lParse,
                                    b"Incompatible dimensions in '?:' arguments\0" as *const u8
                                        as *const c_char
                                        as *mut c_char,
                                );
                                current_block = 4830776507462815627;
                            } else {
                                yyval.Node = New_Func(
                                    lParse,
                                    0,
                                    ifthenelse_fct,
                                    3 as c_int,
                                    (*yyvsp.offset(-2_isize)).Node,
                                    (*yyvsp.offset(0)).Node,
                                    (*yyvsp.offset(-4_isize)).Node,
                                    0,
                                    0,
                                    0,
                                    0,
                                );
                                if yyval.Node < 0 {
                                    current_block = 4830776507462815627;
                                } else {
                                    if (*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-2_isize)).Node as isize))
                                    .value
                                    .nelem
                                        < (*((*lParse).Nodes)
                                            .offset((*yyvsp.offset(0)).Node as isize))
                                        .value
                                        .nelem
                                    {
                                        Copy_Dims(lParse, yyval.Node, (*yyvsp.offset(0)).Node);
                                    }
                                    (*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-4_isize)).Node as isize))
                                    .ntype = (*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-2_isize)).Node as isize))
                                    .ntype;
                                    if Test_Dims(lParse, (*yyvsp.offset(-4_isize)).Node, yyval.Node)
                                        == 0
                                    {
                                        fits_parser_yyerror(
                                            scanner,
                                            lParse,
                                            b"Incompatible dimensions in '?:' condition\0"
                                                as *const u8
                                                as *const c_char
                                                as *mut c_char,
                                        );
                                        current_block = 4830776507462815627;
                                    } else {
                                        (*((*lParse).Nodes)
                                            .offset((*yyvsp.offset(-4_isize)).Node as isize))
                                        .ntype = fits_parser_yytokentype::BOOLEAN as c_int;
                                        if ((*((*lParse).Nodes).offset(yyval.Node as isize)).value)
                                            .nelem
                                            < (*((*lParse).Nodes)
                                                .offset((*yyvsp.offset(-4_isize)).Node as isize))
                                            .value
                                            .nelem
                                        {
                                            Copy_Dims(
                                                lParse,
                                                yyval.Node,
                                                (*yyvsp.offset(-4_isize)).Node,
                                            );
                                        }
                                        current_block = 17353983478346836848;
                                    }
                                }
                            }
                        }
                        53 => {
                            if (*((*lParse).Nodes).offset((*yyvsp.offset(-2_isize)).Node as isize))
                                .ntype
                                > (*((*lParse).Nodes).offset((*yyvsp.offset(0)).Node as isize))
                                    .ntype
                            {
                                (*yyvsp.offset(0)).Node = New_Unary(
                                    lParse,
                                    (*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-2_isize)).Node as isize))
                                    .ntype,
                                    0,
                                    (*yyvsp.offset(0)).Node,
                                );
                            } else if (*((*lParse).Nodes)
                                .offset((*yyvsp.offset(-2_isize)).Node as isize))
                            .ntype
                                < (*((*lParse).Nodes).offset((*yyvsp.offset(0)).Node as isize))
                                    .ntype
                            {
                                (*yyvsp.offset(-2_isize)).Node = New_Unary(
                                    lParse,
                                    (*((*lParse).Nodes).offset((*yyvsp.offset(0)).Node as isize))
                                        .ntype,
                                    0,
                                    (*yyvsp.offset(-2_isize)).Node,
                                );
                            }
                            if Test_Dims(
                                lParse,
                                (*yyvsp.offset(-2_isize)).Node,
                                (*yyvsp.offset(0)).Node,
                            ) == 0
                            {
                                fits_parser_yyerror(
                                    scanner,
                                    lParse,
                                    b"Incompatible dimensions in '?:' arguments\0" as *const u8
                                        as *const c_char
                                        as *mut c_char,
                                );
                                current_block = 4830776507462815627;
                            } else {
                                yyval.Node = New_Func(
                                    lParse,
                                    0,
                                    ifthenelse_fct,
                                    3 as c_int,
                                    (*yyvsp.offset(-2_isize)).Node,
                                    (*yyvsp.offset(0)).Node,
                                    (*yyvsp.offset(-4_isize)).Node,
                                    0,
                                    0,
                                    0,
                                    0,
                                );
                                if yyval.Node < 0 {
                                    current_block = 4830776507462815627;
                                } else {
                                    if (*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-2_isize)).Node as isize))
                                    .value
                                    .nelem
                                        < (*((*lParse).Nodes)
                                            .offset((*yyvsp.offset(0)).Node as isize))
                                        .value
                                        .nelem
                                    {
                                        Copy_Dims(lParse, yyval.Node, (*yyvsp.offset(0)).Node);
                                    }
                                    (*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-4_isize)).Node as isize))
                                    .ntype = (*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-2_isize)).Node as isize))
                                    .ntype;
                                    if Test_Dims(lParse, (*yyvsp.offset(-4_isize)).Node, yyval.Node)
                                        == 0
                                    {
                                        fits_parser_yyerror(
                                            scanner,
                                            lParse,
                                            b"Incompatible dimensions in '?:' condition\0"
                                                as *const u8
                                                as *const c_char
                                                as *mut c_char,
                                        );
                                        current_block = 4830776507462815627;
                                    } else {
                                        (*((*lParse).Nodes)
                                            .offset((*yyvsp.offset(-4_isize)).Node as isize))
                                        .ntype = fits_parser_yytokentype::BOOLEAN as c_int;
                                        if ((*((*lParse).Nodes).offset(yyval.Node as isize)).value)
                                            .nelem
                                            < (*((*lParse).Nodes)
                                                .offset((*yyvsp.offset(-4_isize)).Node as isize))
                                            .value
                                            .nelem
                                        {
                                            Copy_Dims(
                                                lParse,
                                                yyval.Node,
                                                (*yyvsp.offset(-4_isize)).Node,
                                            );
                                        }
                                        current_block = 17353983478346836848;
                                    }
                                }
                            }
                        }
                        54 => {
                            if (if ((*yyvsp.offset(-1_isize)).astr[0] as c_int)
                                < (*::core::mem::transmute::<&[u8; 8], &[c_char; 8]>(b"RANDOM(\0"))
                                    [0] as c_int
                            {
                                -(1)
                            } else if (*yyvsp.offset(-1_isize)).astr[0] as c_int
                                > (*::core::mem::transmute::<&[u8; 8], &[c_char; 8]>(b"RANDOM(\0"))
                                    [0] as c_int
                            {
                                1
                            } else {
                                strcmp(
                                    ((*yyvsp.offset(-1_isize)).astr).as_mut_ptr(),
                                    b"RANDOM(\0" as *const u8 as *const c_char,
                                )
                            }) == 0
                            {
                                yyval.Node = New_Func(
                                    lParse,
                                    fits_parser_yytokentype::DOUBLE as c_int,
                                    rnd_fct,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                );
                                current_block = 1297461190301222800;
                            } else if (if ((*yyvsp.offset(-1_isize)).astr[0] as c_int)
                                < (*::core::mem::transmute::<&[u8; 9], &[c_char; 9]>(b"RANDOMN(\0"))
                                    [0] as c_int
                            {
                                -(1)
                            } else if (*yyvsp.offset(-1_isize)).astr[0] as c_int
                                > (*::core::mem::transmute::<&[u8; 9], &[c_char; 9]>(b"RANDOMN(\0"))
                                    [0] as c_int
                            {
                                1
                            } else {
                                strcmp(
                                    ((*yyvsp.offset(-1_isize)).astr).as_mut_ptr(),
                                    b"RANDOMN(\0" as *const u8 as *const c_char,
                                )
                            }) == 0
                            {
                                yyval.Node = New_Func(
                                    lParse,
                                    fits_parser_yytokentype::DOUBLE as c_int,
                                    gasrnd_fct,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                );
                                current_block = 1297461190301222800;
                            } else {
                                fits_parser_yyerror(
                                    scanner,
                                    lParse,
                                    b"Function() not supported\0" as *const u8 as *const c_char
                                        as *mut c_char,
                                );
                                current_block = 4830776507462815627;
                            }
                            match current_block {
                                4830776507462815627 => {}
                                _ => {
                                    if yyval.Node < 0 {
                                        current_block = 4830776507462815627;
                                    } else {
                                        current_block = 17353983478346836848;
                                    }
                                }
                            }
                        }
                        55 => {
                            if (if ((*yyvsp.offset(-2_isize)).astr[0] as c_int)
                                < (*::core::mem::transmute::<&[u8; 5], &[c_char; 5]>(b"SUM(\0"))[0]
                                    as c_int
                            {
                                -(1)
                            } else if (*yyvsp.offset(-2_isize)).astr[0] as c_int
                                > (*::core::mem::transmute::<&[u8; 5], &[c_char; 5]>(b"SUM(\0"))[0]
                                    as c_int
                            {
                                1
                            } else {
                                strcmp(
                                    ((*yyvsp.offset(-2_isize)).astr).as_mut_ptr(),
                                    b"SUM(\0" as *const u8 as *const c_char,
                                )
                            }) == 0
                            {
                                yyval.Node = New_Func(
                                    lParse,
                                    fits_parser_yytokentype::LONG as c_int,
                                    sum_fct,
                                    1,
                                    (*yyvsp.offset(-1_isize)).Node,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                );
                                current_block = 10848699504537784535;
                            } else if (if ((*yyvsp.offset(-2_isize)).astr[0] as c_int)
                                < (*::core::mem::transmute::<&[u8; 7], &[c_char; 7]>(b"NELEM(\0"))
                                    [0] as c_int
                            {
                                -(1)
                            } else if (*yyvsp.offset(-2_isize)).astr[0] as c_int
                                > (*::core::mem::transmute::<&[u8; 7], &[c_char; 7]>(b"NELEM(\0"))
                                    [0] as c_int
                            {
                                1
                            } else {
                                strcmp(
                                    ((*yyvsp.offset(-2_isize)).astr).as_mut_ptr(),
                                    b"NELEM(\0" as *const u8 as *const c_char,
                                )
                            }) == 0
                            {
                                yyval.Node = New_Const(
                                    lParse,
                                    fits_parser_yytokentype::LONG as c_int,
                                    &mut (*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-1_isize)).Node as isize))
                                    .value
                                    .nelem as *mut c_long
                                        as *mut c_void,
                                    ::core::mem::size_of::<c_long>() as c_ulong as c_long,
                                );
                                current_block = 10848699504537784535;
                            } else if (if ((*yyvsp.offset(-2_isize)).astr[0] as c_int)
                                < (*::core::mem::transmute::<&[u8; 7], &[c_char; 7]>(b"ACCUM(\0"))
                                    [0] as c_int
                            {
                                -(1)
                            } else if (*yyvsp.offset(-2_isize)).astr[0] as c_int
                                > (*::core::mem::transmute::<&[u8; 7], &[c_char; 7]>(b"ACCUM(\0"))
                                    [0] as c_int
                            {
                                1
                            } else {
                                strcmp(
                                    ((*yyvsp.offset(-2_isize)).astr).as_mut_ptr(),
                                    b"ACCUM(\0" as *const u8 as *const c_char,
                                )
                            }) == 0
                            {
                                let mut zero: c_long = 0;
                                yyval.Node = New_BinOp(
                                    lParse,
                                    fits_parser_yytokentype::LONG as c_int,
                                    (*yyvsp.offset(-1_isize)).Node,
                                    fits_parser_yytokentype::ACCUM as c_int,
                                    New_Const(
                                        lParse,
                                        fits_parser_yytokentype::LONG as c_int,
                                        &mut zero as *mut c_long as *mut c_void,
                                        ::core::mem::size_of::<c_long>() as c_ulong as c_long,
                                    ),
                                );
                                current_block = 10848699504537784535;
                            } else {
                                fits_parser_yyerror(
                                    scanner,
                                    lParse,
                                    b"Function(bool) not supported\0" as *const u8 as *const c_char
                                        as *mut c_char,
                                );
                                current_block = 4830776507462815627;
                            }
                            match current_block {
                                4830776507462815627 => {}
                                _ => {
                                    if yyval.Node < 0 {
                                        current_block = 4830776507462815627;
                                    } else {
                                        current_block = 17353983478346836848;
                                    }
                                }
                            }
                        }
                        56 => {
                            if (if ((*yyvsp.offset(-4_isize)).astr[0] as c_int)
                                < (*::core::mem::transmute::<&[u8; 10], &[c_char; 10]>(
                                    b"AXISELEM(\0",
                                ))[0] as c_int
                            {
                                -(1)
                            } else if (*yyvsp.offset(-4_isize)).astr[0] as c_int
                                > (*::core::mem::transmute::<&[u8; 10], &[c_char; 10]>(
                                    b"AXISELEM(\0",
                                ))[0] as c_int
                            {
                                1
                            } else {
                                strcmp(
                                    ((*yyvsp.offset(-4_isize)).astr).as_mut_ptr(),
                                    b"AXISELEM(\0" as *const u8 as *const c_char,
                                )
                            }) == 0
                            {
                                if (*((*lParse).Nodes)
                                    .offset((*yyvsp.offset(-1_isize)).Node as isize))
                                .operation
                                    != -(1000 as c_int)
                                    || (*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-1_isize)).Node as isize))
                                    .value
                                    .nelem
                                        != 1
                                {
                                    fits_parser_yyerror(
                                        scanner,
                                        lParse,
                                        b"AXISELEM second argument must be a scalar constant\0"
                                            as *const u8
                                            as *const c_char
                                            as *mut c_char,
                                    );
                                    current_block = 4830776507462815627;
                                } else if (*((*lParse).Nodes)
                                    .offset((*yyvsp.offset(-3_isize)).Node as isize))
                                .operation
                                    == -(1000 as c_int)
                                {
                                    let mut one: c_long = 1;
                                    yyval.Node = New_Const(
                                        lParse,
                                        fits_parser_yytokentype::LONG as c_int,
                                        &mut one as *mut c_long as *mut c_void,
                                        ::core::mem::size_of::<c_long>() as c_ulong as c_long,
                                    );
                                    current_block = 13755523488868872559;
                                } else {
                                    if (*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-1_isize)).Node as isize))
                                    .ntype
                                        != fits_parser_yytokentype::LONG as c_int
                                    {
                                        (*yyvsp.offset(-1_isize)).Node = New_Unary(
                                            lParse,
                                            fits_parser_yytokentype::LONG as c_int,
                                            0,
                                            (*yyvsp.offset(-1_isize)).Node,
                                        );
                                    }
                                    yyval.Node = New_Func(
                                        lParse,
                                        0,
                                        axiselem_fct,
                                        2 as c_int,
                                        (*yyvsp.offset(-3_isize)).Node,
                                        (*yyvsp.offset(-1_isize)).Node,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                    );
                                    if yyval.Node < 0 {
                                        current_block = 4830776507462815627;
                                    } else {
                                        (*((*lParse).Nodes).offset(yyval.Node as isize)).ntype =
                                            fits_parser_yytokentype::LONG as c_int;
                                        current_block = 13755523488868872559;
                                    }
                                }
                            } else if (if ((*yyvsp.offset(-4_isize)).astr[0] as c_int)
                                < (*::core::mem::transmute::<&[u8; 7], &[c_char; 7]>(b"NAXES(\0"))
                                    [0] as c_int
                            {
                                -(1)
                            } else if (*yyvsp.offset(-4_isize)).astr[0] as c_int
                                > (*::core::mem::transmute::<&[u8; 7], &[c_char; 7]>(b"NAXES(\0"))
                                    [0] as c_int
                            {
                                1
                            } else {
                                strcmp(
                                    ((*yyvsp.offset(-4_isize)).astr).as_mut_ptr(),
                                    b"NAXES(\0" as *const u8 as *const c_char,
                                )
                            }) == 0
                            {
                                if (*((*lParse).Nodes)
                                    .offset((*yyvsp.offset(-1_isize)).Node as isize))
                                .operation
                                    != -(1000 as c_int)
                                    || (*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-1_isize)).Node as isize))
                                    .value
                                    .nelem
                                        != 1
                                {
                                    fits_parser_yyerror(
                                        scanner,
                                        lParse,
                                        b"NAXES second argument must be a scalar constant\0"
                                            as *const u8
                                            as *const c_char
                                            as *mut c_char,
                                    );
                                    current_block = 4830776507462815627;
                                } else if (*((*lParse).Nodes)
                                    .offset((*yyvsp.offset(-3_isize)).Node as isize))
                                .operation
                                    == -(1000 as c_int)
                                {
                                    let mut one_0: c_long = 1;
                                    yyval.Node = New_Const(
                                        lParse,
                                        fits_parser_yytokentype::LONG as c_int,
                                        &mut one_0 as *mut c_long as *mut c_void,
                                        ::core::mem::size_of::<c_long>() as c_ulong as c_long,
                                    );
                                    current_block = 13755523488868872559;
                                } else {
                                    let mut iaxis: c_long = 0;
                                    let mut naxis: c_int = 0;
                                    if (*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-1_isize)).Node as isize))
                                    .ntype
                                        != fits_parser_yytokentype::LONG as c_int
                                    {
                                        (*yyvsp.offset(-1_isize)).Node = New_Unary(
                                            lParse,
                                            fits_parser_yytokentype::LONG as c_int,
                                            0,
                                            (*yyvsp.offset(-1_isize)).Node,
                                        );
                                    }
                                    iaxis = (*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-1_isize)).Node as isize))
                                    .value
                                    .data
                                    .lng;
                                    naxis = (*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-3_isize)).Node as isize))
                                    .value
                                    .naxis;
                                    if iaxis == 0 {
                                        iaxis = naxis as c_long;
                                    } else if iaxis <= naxis as c_long {
                                        iaxis = (*((*lParse).Nodes)
                                            .offset((*yyvsp.offset(-3_isize)).Node as isize))
                                        .value
                                        .naxes[(iaxis - 1) as usize];
                                    } else {
                                        iaxis = 1;
                                    }
                                    yyval.Node = New_Const(
                                        lParse,
                                        fits_parser_yytokentype::LONG as c_int,
                                        &mut iaxis as *mut c_long as *mut c_void,
                                        ::core::mem::size_of::<c_long>() as c_ulong as c_long,
                                    );
                                    if yyval.Node < 0 {
                                        current_block = 4830776507462815627;
                                    } else {
                                        current_block = 13755523488868872559;
                                    }
                                }
                            } else if (if ((*yyvsp.offset(-4_isize)).astr[0] as c_int)
                                < (*::core::mem::transmute::<&[u8; 7], &[c_char; 7]>(b"ARRAY(\0"))
                                    [0] as c_int
                            {
                                -(1)
                            } else if (*yyvsp.offset(-4_isize)).astr[0] as c_int
                                > (*::core::mem::transmute::<&[u8; 7], &[c_char; 7]>(b"ARRAY(\0"))
                                    [0] as c_int
                            {
                                1
                            } else {
                                strcmp(
                                    ((*yyvsp.offset(-4_isize)).astr).as_mut_ptr(),
                                    b"ARRAY(\0" as *const u8 as *const c_char,
                                )
                            }) == 0
                            {
                                yyval.Node = New_Array(
                                    lParse,
                                    (*yyvsp.offset(-3_isize)).Node,
                                    (*yyvsp.offset(-1_isize)).Node,
                                );
                                if yyval.Node < 0 {
                                    current_block = 4830776507462815627;
                                } else {
                                    current_block = 13755523488868872559;
                                }
                            } else {
                                fits_parser_yyerror(
                                    scanner,
                                    lParse,
                                    b"Function(bool,expr) not supported\0" as *const u8
                                        as *const c_char
                                        as *mut c_char,
                                );
                                current_block = 4830776507462815627;
                            }
                            match current_block {
                                4830776507462815627 => {}
                                _ => {
                                    if yyval.Node < 0 {
                                        current_block = 4830776507462815627;
                                    } else {
                                        current_block = 17353983478346836848;
                                    }
                                }
                            }
                        }
                        57 => {
                            if (if ((*yyvsp.offset(-2_isize)).astr[0] as c_int)
                                < (*::core::mem::transmute::<&[u8; 7], &[c_char; 7]>(b"NELEM(\0"))
                                    [0] as c_int
                            {
                                -(1)
                            } else if (*yyvsp.offset(-2_isize)).astr[0] as c_int
                                > (*::core::mem::transmute::<&[u8; 7], &[c_char; 7]>(b"NELEM(\0"))
                                    [0] as c_int
                            {
                                1
                            } else {
                                strcmp(
                                    ((*yyvsp.offset(-2_isize)).astr).as_mut_ptr(),
                                    b"NELEM(\0" as *const u8 as *const c_char,
                                )
                            }) == 0
                            {
                                yyval.Node = New_Const(
                                    lParse,
                                    fits_parser_yytokentype::LONG as c_int,
                                    &mut (*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-1_isize)).Node as isize))
                                    .value
                                    .nelem as *mut c_long
                                        as *mut c_void,
                                    ::core::mem::size_of::<c_long>() as c_ulong as c_long,
                                );
                                current_block = 15752106442776732052;
                            } else if (if ((*yyvsp.offset(-2_isize)).astr[0] as c_int)
                                < (*::core::mem::transmute::<&[u8; 8], &[c_char; 8]>(b"NVALID(\0"))
                                    [0] as c_int
                            {
                                -(1)
                            } else if (*yyvsp.offset(-2_isize)).astr[0] as c_int
                                > (*::core::mem::transmute::<&[u8; 8], &[c_char; 8]>(b"NVALID(\0"))
                                    [0] as c_int
                            {
                                1
                            } else {
                                strcmp(
                                    ((*yyvsp.offset(-2_isize)).astr).as_mut_ptr(),
                                    b"NVALID(\0" as *const u8 as *const c_char,
                                )
                            }) == 0
                            {
                                yyval.Node = New_Func(
                                    lParse,
                                    fits_parser_yytokentype::LONG as c_int,
                                    nonnull_fct,
                                    1,
                                    (*yyvsp.offset(-1_isize)).Node,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                );
                                current_block = 15752106442776732052;
                            } else {
                                fits_parser_yyerror(
                                    scanner,
                                    lParse,
                                    b"Function(str) not supported\0" as *const u8 as *const c_char
                                        as *mut c_char,
                                );
                                current_block = 4830776507462815627;
                            }
                            match current_block {
                                4830776507462815627 => {}
                                _ => {
                                    if yyval.Node < 0 {
                                        current_block = 4830776507462815627;
                                    } else {
                                        current_block = 17353983478346836848;
                                    }
                                }
                            }
                        }
                        58 => {
                            if (if ((*yyvsp.offset(-2_isize)).astr[0] as c_int)
                                < (*::core::mem::transmute::<&[u8; 7], &[c_char; 7]>(b"NELEM(\0"))
                                    [0] as c_int
                            {
                                -(1)
                            } else if (*yyvsp.offset(-2_isize)).astr[0] as c_int
                                > (*::core::mem::transmute::<&[u8; 7], &[c_char; 7]>(b"NELEM(\0"))
                                    [0] as c_int
                            {
                                1
                            } else {
                                strcmp(
                                    ((*yyvsp.offset(-2_isize)).astr).as_mut_ptr(),
                                    b"NELEM(\0" as *const u8 as *const c_char,
                                )
                            }) == 0
                            {
                                yyval.Node = New_Const(
                                    lParse,
                                    fits_parser_yytokentype::LONG as c_int,
                                    &mut (*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-1_isize)).Node as isize))
                                    .value
                                    .nelem as *mut c_long
                                        as *mut c_void,
                                    ::core::mem::size_of::<c_long>() as c_ulong as c_long,
                                );
                                current_block = 494012601817399562;
                            } else if (if ((*yyvsp.offset(-2_isize)).astr[0] as c_int)
                                < (*::core::mem::transmute::<&[u8; 8], &[c_char; 8]>(b"NVALID(\0"))
                                    [0] as c_int
                            {
                                -(1)
                            } else if (*yyvsp.offset(-2_isize)).astr[0] as c_int
                                > (*::core::mem::transmute::<&[u8; 8], &[c_char; 8]>(b"NVALID(\0"))
                                    [0] as c_int
                            {
                                1
                            } else {
                                strcmp(
                                    ((*yyvsp.offset(-2_isize)).astr).as_mut_ptr(),
                                    b"NVALID(\0" as *const u8 as *const c_char,
                                )
                            }) == 0
                            {
                                yyval.Node = New_Const(
                                    lParse,
                                    fits_parser_yytokentype::LONG as c_int,
                                    &mut (*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-1_isize)).Node as isize))
                                    .value
                                    .nelem as *mut c_long
                                        as *mut c_void,
                                    ::core::mem::size_of::<c_long>() as c_ulong as c_long,
                                );
                                current_block = 494012601817399562;
                            } else if (if ((*yyvsp.offset(-2_isize)).astr[0] as c_int)
                                < (*::core::mem::transmute::<&[u8; 5], &[c_char; 5]>(b"SUM(\0"))[0]
                                    as c_int
                            {
                                -(1)
                            } else if (*yyvsp.offset(-2_isize)).astr[0] as c_int
                                > (*::core::mem::transmute::<&[u8; 5], &[c_char; 5]>(b"SUM(\0"))[0]
                                    as c_int
                            {
                                1
                            } else {
                                strcmp(
                                    ((*yyvsp.offset(-2_isize)).astr).as_mut_ptr(),
                                    b"SUM(\0" as *const u8 as *const c_char,
                                )
                            }) == 0
                            {
                                yyval.Node = New_Func(
                                    lParse,
                                    fits_parser_yytokentype::LONG as c_int,
                                    sum_fct,
                                    1,
                                    (*yyvsp.offset(-1_isize)).Node,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                );
                                current_block = 494012601817399562;
                            } else if (if ((*yyvsp.offset(-2_isize)).astr[0] as c_int)
                                < (*::core::mem::transmute::<&[u8; 5], &[c_char; 5]>(b"MIN(\0"))[0]
                                    as c_int
                            {
                                -(1)
                            } else if (*yyvsp.offset(-2_isize)).astr[0] as c_int
                                > (*::core::mem::transmute::<&[u8; 5], &[c_char; 5]>(b"MIN(\0"))[0]
                                    as c_int
                            {
                                1
                            } else {
                                strcmp(
                                    ((*yyvsp.offset(-2_isize)).astr).as_mut_ptr(),
                                    b"MIN(\0" as *const u8 as *const c_char,
                                )
                            }) == 0
                            {
                                yyval.Node = New_Func(
                                    lParse,
                                    (*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-1_isize)).Node as isize))
                                    .ntype,
                                    min1_fct,
                                    1,
                                    (*yyvsp.offset(-1_isize)).Node,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                );
                                ((*((*lParse).Nodes).offset(yyval.Node as isize)).value).nelem = 1;
                                current_block = 494012601817399562;
                            } else if (if ((*yyvsp.offset(-2_isize)).astr[0] as c_int)
                                < (*::core::mem::transmute::<&[u8; 7], &[c_char; 7]>(b"ACCUM(\0"))
                                    [0] as c_int
                            {
                                -(1)
                            } else if (*yyvsp.offset(-2_isize)).astr[0] as c_int
                                > (*::core::mem::transmute::<&[u8; 7], &[c_char; 7]>(b"ACCUM(\0"))
                                    [0] as c_int
                            {
                                1
                            } else {
                                strcmp(
                                    ((*yyvsp.offset(-2_isize)).astr).as_mut_ptr(),
                                    b"ACCUM(\0" as *const u8 as *const c_char,
                                )
                            }) == 0
                            {
                                let mut zero_0: c_long = 0;
                                yyval.Node = New_BinOp(
                                    lParse,
                                    fits_parser_yytokentype::LONG as c_int,
                                    (*yyvsp.offset(-1_isize)).Node,
                                    fits_parser_yytokentype::ACCUM as c_int,
                                    New_Const(
                                        lParse,
                                        fits_parser_yytokentype::LONG as c_int,
                                        &mut zero_0 as *mut c_long as *mut c_void,
                                        ::core::mem::size_of::<c_long>() as c_ulong as c_long,
                                    ),
                                );
                                current_block = 494012601817399562;
                            } else if (if ((*yyvsp.offset(-2_isize)).astr[0] as c_int)
                                < (*::core::mem::transmute::<&[u8; 5], &[c_char; 5]>(b"MAX(\0"))[0]
                                    as c_int
                            {
                                -(1)
                            } else if (*yyvsp.offset(-2_isize)).astr[0] as c_int
                                > (*::core::mem::transmute::<&[u8; 5], &[c_char; 5]>(b"MAX(\0"))[0]
                                    as c_int
                            {
                                1
                            } else {
                                strcmp(
                                    ((*yyvsp.offset(-2_isize)).astr).as_mut_ptr(),
                                    b"MAX(\0" as *const u8 as *const c_char,
                                )
                            }) == 0
                            {
                                yyval.Node = New_Func(
                                    lParse,
                                    (*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-1_isize)).Node as isize))
                                    .ntype,
                                    max1_fct,
                                    1,
                                    (*yyvsp.offset(-1_isize)).Node,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                );
                                ((*((*lParse).Nodes).offset(yyval.Node as isize)).value).nelem = 1;
                                current_block = 494012601817399562;
                            } else {
                                fits_parser_yyerror(
                                    scanner,
                                    lParse,
                                    b"Function(bits) not supported\0" as *const u8 as *const c_char
                                        as *mut c_char,
                                );
                                current_block = 4830776507462815627;
                            }
                            match current_block {
                                4830776507462815627 => {}
                                _ => {
                                    if yyval.Node < 0 {
                                        current_block = 4830776507462815627;
                                    } else {
                                        current_block = 17353983478346836848;
                                    }
                                }
                            }
                        }
                        59 => {
                            if (if ((*yyvsp.offset(-2_isize)).astr[0] as c_int)
                                < (*::core::mem::transmute::<&[u8; 5], &[c_char; 5]>(b"SUM(\0"))[0]
                                    as c_int
                            {
                                -(1)
                            } else if (*yyvsp.offset(-2_isize)).astr[0] as c_int
                                > (*::core::mem::transmute::<&[u8; 5], &[c_char; 5]>(b"SUM(\0"))[0]
                                    as c_int
                            {
                                1
                            } else {
                                strcmp(
                                    ((*yyvsp.offset(-2_isize)).astr).as_mut_ptr(),
                                    b"SUM(\0" as *const u8 as *const c_char,
                                )
                            }) == 0
                            {
                                yyval.Node = New_Func(
                                    lParse,
                                    (*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-1_isize)).Node as isize))
                                    .ntype,
                                    sum_fct,
                                    1,
                                    (*yyvsp.offset(-1_isize)).Node,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                );
                                current_block = 7600445499126923600;
                            } else if (if ((*yyvsp.offset(-2_isize)).astr[0] as c_int)
                                < (*::core::mem::transmute::<&[u8; 9], &[c_char; 9]>(b"AVERAGE(\0"))
                                    [0] as c_int
                            {
                                -(1)
                            } else if (*yyvsp.offset(-2_isize)).astr[0] as c_int
                                > (*::core::mem::transmute::<&[u8; 9], &[c_char; 9]>(b"AVERAGE(\0"))
                                    [0] as c_int
                            {
                                1
                            } else {
                                strcmp(
                                    ((*yyvsp.offset(-2_isize)).astr).as_mut_ptr(),
                                    b"AVERAGE(\0" as *const u8 as *const c_char,
                                )
                            }) == 0
                            {
                                yyval.Node = New_Func(
                                    lParse,
                                    fits_parser_yytokentype::DOUBLE as c_int,
                                    average_fct,
                                    1,
                                    (*yyvsp.offset(-1_isize)).Node,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                );
                                current_block = 7600445499126923600;
                            } else if (if ((*yyvsp.offset(-2_isize)).astr[0] as c_int)
                                < (*::core::mem::transmute::<&[u8; 8], &[c_char; 8]>(b"STDDEV(\0"))
                                    [0] as c_int
                            {
                                -(1)
                            } else if (*yyvsp.offset(-2_isize)).astr[0] as c_int
                                > (*::core::mem::transmute::<&[u8; 8], &[c_char; 8]>(b"STDDEV(\0"))
                                    [0] as c_int
                            {
                                1
                            } else {
                                strcmp(
                                    ((*yyvsp.offset(-2_isize)).astr).as_mut_ptr(),
                                    b"STDDEV(\0" as *const u8 as *const c_char,
                                )
                            }) == 0
                            {
                                yyval.Node = New_Func(
                                    lParse,
                                    fits_parser_yytokentype::DOUBLE as c_int,
                                    stddev_fct,
                                    1,
                                    (*yyvsp.offset(-1_isize)).Node,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                );
                                current_block = 7600445499126923600;
                            } else if (if ((*yyvsp.offset(-2_isize)).astr[0] as c_int)
                                < (*::core::mem::transmute::<&[u8; 8], &[c_char; 8]>(b"MEDIAN(\0"))
                                    [0] as c_int
                            {
                                -(1)
                            } else if (*yyvsp.offset(-2_isize)).astr[0] as c_int
                                > (*::core::mem::transmute::<&[u8; 8], &[c_char; 8]>(b"MEDIAN(\0"))
                                    [0] as c_int
                            {
                                1
                            } else {
                                strcmp(
                                    ((*yyvsp.offset(-2_isize)).astr).as_mut_ptr(),
                                    b"MEDIAN(\0" as *const u8 as *const c_char,
                                )
                            }) == 0
                            {
                                yyval.Node = New_Func(
                                    lParse,
                                    (*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-1_isize)).Node as isize))
                                    .ntype,
                                    median_fct,
                                    1,
                                    (*yyvsp.offset(-1_isize)).Node,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                );
                                current_block = 7600445499126923600;
                            } else if (if ((*yyvsp.offset(-2_isize)).astr[0] as c_int)
                                < (*::core::mem::transmute::<&[u8; 7], &[c_char; 7]>(b"NELEM(\0"))
                                    [0] as c_int
                            {
                                -(1)
                            } else if (*yyvsp.offset(-2_isize)).astr[0] as c_int
                                > (*::core::mem::transmute::<&[u8; 7], &[c_char; 7]>(b"NELEM(\0"))
                                    [0] as c_int
                            {
                                1
                            } else {
                                strcmp(
                                    ((*yyvsp.offset(-2_isize)).astr).as_mut_ptr(),
                                    b"NELEM(\0" as *const u8 as *const c_char,
                                )
                            }) == 0
                            {
                                yyval.Node = New_Const(
                                    lParse,
                                    fits_parser_yytokentype::LONG as c_int,
                                    &mut (*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-1_isize)).Node as isize))
                                    .value
                                    .nelem as *mut c_long
                                        as *mut c_void,
                                    ::core::mem::size_of::<c_long>() as c_ulong as c_long,
                                );
                                current_block = 7600445499126923600;
                            } else if (if ((*yyvsp.offset(-2_isize)).astr[0] as c_int)
                                < (*::core::mem::transmute::<&[u8; 8], &[c_char; 8]>(b"NVALID(\0"))
                                    [0] as c_int
                            {
                                -(1)
                            } else if (*yyvsp.offset(-2_isize)).astr[0] as c_int
                                > (*::core::mem::transmute::<&[u8; 8], &[c_char; 8]>(b"NVALID(\0"))
                                    [0] as c_int
                            {
                                1
                            } else {
                                strcmp(
                                    ((*yyvsp.offset(-2_isize)).astr).as_mut_ptr(),
                                    b"NVALID(\0" as *const u8 as *const c_char,
                                )
                            }) == 0
                            {
                                yyval.Node = New_Func(
                                    lParse,
                                    fits_parser_yytokentype::LONG as c_int,
                                    nonnull_fct,
                                    1,
                                    (*yyvsp.offset(-1_isize)).Node,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                );
                                current_block = 7600445499126923600;
                            } else if (if ((*yyvsp.offset(-2_isize)).astr[0] as c_int)
                                < (*::core::mem::transmute::<&[u8; 7], &[c_char; 7]>(b"ACCUM(\0"))
                                    [0] as c_int
                            {
                                -(1)
                            } else if (*yyvsp.offset(-2_isize)).astr[0] as c_int
                                > (*::core::mem::transmute::<&[u8; 7], &[c_char; 7]>(b"ACCUM(\0"))
                                    [0] as c_int
                            {
                                1
                            } else {
                                strcmp(
                                    ((*yyvsp.offset(-2_isize)).astr).as_mut_ptr(),
                                    b"ACCUM(\0" as *const u8 as *const c_char,
                                )
                            }) == 0
                                && (*((*lParse).Nodes)
                                    .offset((*yyvsp.offset(-1_isize)).Node as isize))
                                .ntype
                                    == fits_parser_yytokentype::LONG as c_int
                            {
                                let mut zero_1: c_long = 0;
                                yyval.Node = New_BinOp(
                                    lParse,
                                    fits_parser_yytokentype::LONG as c_int,
                                    (*yyvsp.offset(-1_isize)).Node,
                                    fits_parser_yytokentype::ACCUM as c_int,
                                    New_Const(
                                        lParse,
                                        fits_parser_yytokentype::LONG as c_int,
                                        &mut zero_1 as *mut c_long as *mut c_void,
                                        ::core::mem::size_of::<c_long>() as c_ulong as c_long,
                                    ),
                                );
                                current_block = 7600445499126923600;
                            } else if (if ((*yyvsp.offset(-2_isize)).astr[0] as c_int)
                                < (*::core::mem::transmute::<&[u8; 7], &[c_char; 7]>(b"ACCUM(\0"))
                                    [0] as c_int
                            {
                                -(1)
                            } else if (*yyvsp.offset(-2_isize)).astr[0] as c_int
                                > (*::core::mem::transmute::<&[u8; 7], &[c_char; 7]>(b"ACCUM(\0"))
                                    [0] as c_int
                            {
                                1
                            } else {
                                strcmp(
                                    ((*yyvsp.offset(-2_isize)).astr).as_mut_ptr(),
                                    b"ACCUM(\0" as *const u8 as *const c_char,
                                )
                            }) == 0
                                && (*((*lParse).Nodes)
                                    .offset((*yyvsp.offset(-1_isize)).Node as isize))
                                .ntype
                                    == fits_parser_yytokentype::DOUBLE as c_int
                            {
                                let mut zero_2: c_double = 0.0;
                                yyval.Node = New_BinOp(
                                    lParse,
                                    fits_parser_yytokentype::DOUBLE as c_int,
                                    (*yyvsp.offset(-1_isize)).Node,
                                    fits_parser_yytokentype::ACCUM as c_int,
                                    New_Const(
                                        lParse,
                                        fits_parser_yytokentype::DOUBLE as c_int,
                                        &mut zero_2 as *mut c_double as *mut c_void,
                                        ::core::mem::size_of::<c_double>() as c_ulong as c_long,
                                    ),
                                );
                                current_block = 7600445499126923600;
                            } else if (if ((*yyvsp.offset(-2_isize)).astr[0] as c_int)
                                < (*::core::mem::transmute::<&[u8; 9], &[c_char; 9]>(b"SEQDIFF(\0"))
                                    [0] as c_int
                            {
                                -(1)
                            } else if (*yyvsp.offset(-2_isize)).astr[0] as c_int
                                > (*::core::mem::transmute::<&[u8; 9], &[c_char; 9]>(b"SEQDIFF(\0"))
                                    [0] as c_int
                            {
                                1
                            } else {
                                strcmp(
                                    ((*yyvsp.offset(-2_isize)).astr).as_mut_ptr(),
                                    b"SEQDIFF(\0" as *const u8 as *const c_char,
                                )
                            }) == 0
                                && (*((*lParse).Nodes)
                                    .offset((*yyvsp.offset(-1_isize)).Node as isize))
                                .ntype
                                    == fits_parser_yytokentype::LONG as c_int
                            {
                                let mut zero_3: c_long = 0;
                                yyval.Node = New_BinOp(
                                    lParse,
                                    fits_parser_yytokentype::LONG as c_int,
                                    (*yyvsp.offset(-1_isize)).Node,
                                    fits_parser_yytokentype::DIFF as c_int,
                                    New_Const(
                                        lParse,
                                        fits_parser_yytokentype::LONG as c_int,
                                        &mut zero_3 as *mut c_long as *mut c_void,
                                        ::core::mem::size_of::<c_long>() as c_ulong as c_long,
                                    ),
                                );
                                current_block = 7600445499126923600;
                            } else if (if ((*yyvsp.offset(-2_isize)).astr[0] as c_int)
                                < (*::core::mem::transmute::<&[u8; 9], &[c_char; 9]>(b"SEQDIFF(\0"))
                                    [0] as c_int
                            {
                                -(1)
                            } else if (*yyvsp.offset(-2_isize)).astr[0] as c_int
                                > (*::core::mem::transmute::<&[u8; 9], &[c_char; 9]>(b"SEQDIFF(\0"))
                                    [0] as c_int
                            {
                                1
                            } else {
                                strcmp(
                                    ((*yyvsp.offset(-2_isize)).astr).as_mut_ptr(),
                                    b"SEQDIFF(\0" as *const u8 as *const c_char,
                                )
                            }) == 0
                                && (*((*lParse).Nodes)
                                    .offset((*yyvsp.offset(-1_isize)).Node as isize))
                                .ntype
                                    == fits_parser_yytokentype::DOUBLE as c_int
                            {
                                let mut zero_4: c_double = 0.0;
                                yyval.Node = New_BinOp(
                                    lParse,
                                    fits_parser_yytokentype::DOUBLE as c_int,
                                    (*yyvsp.offset(-1_isize)).Node,
                                    fits_parser_yytokentype::DIFF as c_int,
                                    New_Const(
                                        lParse,
                                        fits_parser_yytokentype::DOUBLE as c_int,
                                        &mut zero_4 as *mut c_double as *mut c_void,
                                        ::core::mem::size_of::<c_double>() as c_ulong as c_long,
                                    ),
                                );
                                current_block = 7600445499126923600;
                            } else if (if ((*yyvsp.offset(-2_isize)).astr[0] as c_int)
                                < (*::core::mem::transmute::<&[u8; 5], &[c_char; 5]>(b"ABS(\0"))[0]
                                    as c_int
                            {
                                -(1)
                            } else if (*yyvsp.offset(-2_isize)).astr[0] as c_int
                                > (*::core::mem::transmute::<&[u8; 5], &[c_char; 5]>(b"ABS(\0"))[0]
                                    as c_int
                            {
                                1
                            } else {
                                strcmp(
                                    ((*yyvsp.offset(-2_isize)).astr).as_mut_ptr(),
                                    b"ABS(\0" as *const u8 as *const c_char,
                                )
                            }) == 0
                            {
                                yyval.Node = New_Func(
                                    lParse,
                                    0,
                                    abs_fct,
                                    1,
                                    (*yyvsp.offset(-1_isize)).Node,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                );
                                current_block = 7600445499126923600;
                            } else if (if ((*yyvsp.offset(-2_isize)).astr[0] as c_int)
                                < (*::core::mem::transmute::<&[u8; 5], &[c_char; 5]>(b"MIN(\0"))[0]
                                    as c_int
                            {
                                -(1)
                            } else if (*yyvsp.offset(-2_isize)).astr[0] as c_int
                                > (*::core::mem::transmute::<&[u8; 5], &[c_char; 5]>(b"MIN(\0"))[0]
                                    as c_int
                            {
                                1
                            } else {
                                strcmp(
                                    ((*yyvsp.offset(-2_isize)).astr).as_mut_ptr(),
                                    b"MIN(\0" as *const u8 as *const c_char,
                                )
                            }) == 0
                            {
                                yyval.Node = New_Func(
                                    lParse,
                                    (*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-1_isize)).Node as isize))
                                    .ntype,
                                    min1_fct,
                                    1,
                                    (*yyvsp.offset(-1_isize)).Node,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                );
                                current_block = 7600445499126923600;
                            } else if (if ((*yyvsp.offset(-2_isize)).astr[0] as c_int)
                                < (*::core::mem::transmute::<&[u8; 5], &[c_char; 5]>(b"MAX(\0"))[0]
                                    as c_int
                            {
                                -(1)
                            } else if (*yyvsp.offset(-2_isize)).astr[0] as c_int
                                > (*::core::mem::transmute::<&[u8; 5], &[c_char; 5]>(b"MAX(\0"))[0]
                                    as c_int
                            {
                                1
                            } else {
                                strcmp(
                                    ((*yyvsp.offset(-2_isize)).astr).as_mut_ptr(),
                                    b"MAX(\0" as *const u8 as *const c_char,
                                )
                            }) == 0
                            {
                                yyval.Node = New_Func(
                                    lParse,
                                    (*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-1_isize)).Node as isize))
                                    .ntype,
                                    max1_fct,
                                    1,
                                    (*yyvsp.offset(-1_isize)).Node,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                );
                                current_block = 7600445499126923600;
                            } else if (if ((*yyvsp.offset(-2_isize)).astr[0] as c_int)
                                < (*::core::mem::transmute::<&[u8; 8], &[c_char; 8]>(b"RANDOM(\0"))
                                    [0] as c_int
                            {
                                -(1)
                            } else if (*yyvsp.offset(-2_isize)).astr[0] as c_int
                                > (*::core::mem::transmute::<&[u8; 8], &[c_char; 8]>(b"RANDOM(\0"))
                                    [0] as c_int
                            {
                                1
                            } else {
                                strcmp(
                                    ((*yyvsp.offset(-2_isize)).astr).as_mut_ptr(),
                                    b"RANDOM(\0" as *const u8 as *const c_char,
                                )
                            }) == 0
                            {
                                yyval.Node = New_Func(
                                    lParse,
                                    0,
                                    rnd_fct,
                                    1,
                                    (*yyvsp.offset(-1_isize)).Node,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                );
                                if yyval.Node < 0 {
                                    current_block = 4830776507462815627;
                                } else {
                                    (*((*lParse).Nodes).offset(yyval.Node as isize)).ntype =
                                        fits_parser_yytokentype::DOUBLE as c_int;
                                    current_block = 7600445499126923600;
                                }
                            } else if (if ((*yyvsp.offset(-2_isize)).astr[0] as c_int)
                                < (*::core::mem::transmute::<&[u8; 9], &[c_char; 9]>(b"RANDOMN(\0"))
                                    [0] as c_int
                            {
                                -(1)
                            } else if (*yyvsp.offset(-2_isize)).astr[0] as c_int
                                > (*::core::mem::transmute::<&[u8; 9], &[c_char; 9]>(b"RANDOMN(\0"))
                                    [0] as c_int
                            {
                                1
                            } else {
                                strcmp(
                                    ((*yyvsp.offset(-2_isize)).astr).as_mut_ptr(),
                                    b"RANDOMN(\0" as *const u8 as *const c_char,
                                )
                            }) == 0
                            {
                                yyval.Node = New_Func(
                                    lParse,
                                    0,
                                    gasrnd_fct,
                                    1,
                                    (*yyvsp.offset(-1_isize)).Node,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                );
                                if yyval.Node < 0 {
                                    current_block = 4830776507462815627;
                                } else {
                                    (*((*lParse).Nodes).offset(yyval.Node as isize)).ntype =
                                        fits_parser_yytokentype::DOUBLE as c_int;
                                    current_block = 7600445499126923600;
                                }
                            } else if (if ((*yyvsp.offset(-2_isize)).astr[0] as c_int)
                                < (*::core::mem::transmute::<&[u8; 12], &[c_char; 12]>(
                                    b"ELEMENTNUM(\0",
                                ))[0] as c_int
                            {
                                -(1)
                            } else if (*yyvsp.offset(-2_isize)).astr[0] as c_int
                                > (*::core::mem::transmute::<&[u8; 12], &[c_char; 12]>(
                                    b"ELEMENTNUM(\0",
                                ))[0] as c_int
                            {
                                1
                            } else {
                                strcmp(
                                    ((*yyvsp.offset(-2_isize)).astr).as_mut_ptr(),
                                    b"ELEMENTNUM(\0" as *const u8 as *const c_char,
                                )
                            }) == 0
                            {
                                if (*((*lParse).Nodes)
                                    .offset((*yyvsp.offset(-1_isize)).Node as isize))
                                .operation
                                    == -(1000 as c_int)
                                {
                                    let mut one_1: c_long = 1;
                                    yyval.Node = New_Const(
                                        lParse,
                                        fits_parser_yytokentype::LONG as c_int,
                                        &mut one_1 as *mut c_long as *mut c_void,
                                        ::core::mem::size_of::<c_long>() as c_ulong as c_long,
                                    );
                                    current_block = 7600445499126923600;
                                } else {
                                    yyval.Node = New_Func(
                                        lParse,
                                        0,
                                        elemnum_fct,
                                        1,
                                        (*yyvsp.offset(-1_isize)).Node,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                    );
                                    if yyval.Node < 0 {
                                        current_block = 4830776507462815627;
                                    } else {
                                        (*((*lParse).Nodes).offset(yyval.Node as isize)).ntype =
                                            fits_parser_yytokentype::LONG as c_int;
                                        current_block = 7600445499126923600;
                                    }
                                }
                            } else if (if ((*yyvsp.offset(-2_isize)).astr[0] as c_int)
                                < (*::core::mem::transmute::<&[u8; 7], &[c_char; 7]>(b"NAXIS(\0"))
                                    [0] as c_int
                            {
                                -(1)
                            } else if (*yyvsp.offset(-2_isize)).astr[0] as c_int
                                > (*::core::mem::transmute::<&[u8; 7], &[c_char; 7]>(b"NAXIS(\0"))
                                    [0] as c_int
                            {
                                1
                            } else {
                                strcmp(
                                    ((*yyvsp.offset(-2_isize)).astr).as_mut_ptr(),
                                    b"NAXIS(\0" as *const u8 as *const c_char,
                                )
                            }) == 0
                            {
                                if (*((*lParse).Nodes)
                                    .offset((*yyvsp.offset(-1_isize)).Node as isize))
                                .operation
                                    == -(1000 as c_int)
                                {
                                    let mut one_2: c_long = 1;
                                    yyval.Node = New_Const(
                                        lParse,
                                        fits_parser_yytokentype::LONG as c_int,
                                        &mut one_2 as *mut c_long as *mut c_void,
                                        ::core::mem::size_of::<c_long>() as c_ulong as c_long,
                                    );
                                    current_block = 7600445499126923600;
                                } else {
                                    let mut naxis_0: c_long = (*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-1_isize)).Node as isize))
                                    .value
                                    .naxis
                                        as c_long;
                                    yyval.Node = New_Const(
                                        lParse,
                                        fits_parser_yytokentype::LONG as c_int,
                                        &mut naxis_0 as *mut c_long as *mut c_void,
                                        ::core::mem::size_of::<c_long>() as c_ulong as c_long,
                                    );
                                    if yyval.Node < 0 {
                                        current_block = 4830776507462815627;
                                    } else {
                                        current_block = 7600445499126923600;
                                    }
                                }
                            } else {
                                if (*((*lParse).Nodes)
                                    .offset((*yyvsp.offset(-1_isize)).Node as isize))
                                .ntype
                                    != fits_parser_yytokentype::DOUBLE as c_int
                                {
                                    (*yyvsp.offset(-1_isize)).Node = New_Unary(
                                        lParse,
                                        fits_parser_yytokentype::DOUBLE as c_int,
                                        0,
                                        (*yyvsp.offset(-1_isize)).Node,
                                    );
                                }
                                if (if ((*yyvsp.offset(-2_isize)).astr[0] as c_int)
                                    < (*::core::mem::transmute::<&[u8; 5], &[c_char; 5]>(b"SIN(\0"))
                                        [0] as c_int
                                {
                                    -(1)
                                } else if (*yyvsp.offset(-2_isize)).astr[0] as c_int
                                    > (*::core::mem::transmute::<&[u8; 5], &[c_char; 5]>(b"SIN(\0"))
                                        [0] as c_int
                                {
                                    1
                                } else {
                                    strcmp(
                                        ((*yyvsp.offset(-2_isize)).astr).as_mut_ptr(),
                                        b"SIN(\0" as *const u8 as *const c_char,
                                    )
                                }) == 0
                                {
                                    yyval.Node = New_Func(
                                        lParse,
                                        0,
                                        sin_fct,
                                        1,
                                        (*yyvsp.offset(-1_isize)).Node,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                    );
                                    current_block = 7600445499126923600;
                                } else if (if ((*yyvsp.offset(-2_isize)).astr[0] as c_int)
                                    < (*::core::mem::transmute::<&[u8; 5], &[c_char; 5]>(b"COS(\0"))
                                        [0] as c_int
                                {
                                    -(1)
                                } else if (*yyvsp.offset(-2_isize)).astr[0] as c_int
                                    > (*::core::mem::transmute::<&[u8; 5], &[c_char; 5]>(b"COS(\0"))
                                        [0] as c_int
                                {
                                    1
                                } else {
                                    strcmp(
                                        ((*yyvsp.offset(-2_isize)).astr).as_mut_ptr(),
                                        b"COS(\0" as *const u8 as *const c_char,
                                    )
                                }) == 0
                                {
                                    yyval.Node = New_Func(
                                        lParse,
                                        0,
                                        cos_fct,
                                        1,
                                        (*yyvsp.offset(-1_isize)).Node,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                    );
                                    current_block = 7600445499126923600;
                                } else if (if ((*yyvsp.offset(-2_isize)).astr[0] as c_int)
                                    < (*::core::mem::transmute::<&[u8; 5], &[c_char; 5]>(b"TAN(\0"))
                                        [0] as c_int
                                {
                                    -(1)
                                } else if (*yyvsp.offset(-2_isize)).astr[0] as c_int
                                    > (*::core::mem::transmute::<&[u8; 5], &[c_char; 5]>(b"TAN(\0"))
                                        [0] as c_int
                                {
                                    1
                                } else {
                                    strcmp(
                                        ((*yyvsp.offset(-2_isize)).astr).as_mut_ptr(),
                                        b"TAN(\0" as *const u8 as *const c_char,
                                    )
                                }) == 0
                                {
                                    yyval.Node = New_Func(
                                        lParse,
                                        0,
                                        tan_fct,
                                        1,
                                        (*yyvsp.offset(-1_isize)).Node,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                    );
                                    current_block = 7600445499126923600;
                                } else if (if ((*yyvsp.offset(-2_isize)).astr[0] as c_int)
                                    < (*::core::mem::transmute::<&[u8; 8], &[c_char; 8]>(
                                        b"ARCSIN(\0",
                                    ))[0] as c_int
                                {
                                    -(1)
                                } else if (*yyvsp.offset(-2_isize)).astr[0] as c_int
                                    > (*::core::mem::transmute::<&[u8; 8], &[c_char; 8]>(
                                        b"ARCSIN(\0",
                                    ))[0] as c_int
                                {
                                    1
                                } else {
                                    strcmp(
                                        ((*yyvsp.offset(-2_isize)).astr).as_mut_ptr(),
                                        b"ARCSIN(\0" as *const u8 as *const c_char,
                                    )
                                }) == 0
                                    || (if ((*yyvsp.offset(-2_isize)).astr[0] as c_int)
                                        < (*::core::mem::transmute::<&[u8; 6], &[c_char; 6]>(
                                            b"ASIN(\0",
                                        ))[0] as c_int
                                    {
                                        -(1)
                                    } else if (*yyvsp.offset(-2_isize)).astr[0] as c_int
                                        > (*::core::mem::transmute::<&[u8; 6], &[c_char; 6]>(
                                            b"ASIN(\0",
                                        ))[0] as c_int
                                    {
                                        1
                                    } else {
                                        strcmp(
                                            ((*yyvsp.offset(-2_isize)).astr).as_mut_ptr(),
                                            b"ASIN(\0" as *const u8 as *const c_char,
                                        )
                                    }) == 0
                                {
                                    yyval.Node = New_Func(
                                        lParse,
                                        0,
                                        asin_fct,
                                        1,
                                        (*yyvsp.offset(-1_isize)).Node,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                    );
                                    current_block = 7600445499126923600;
                                } else if (if ((*yyvsp.offset(-2_isize)).astr[0] as c_int)
                                    < (*::core::mem::transmute::<&[u8; 8], &[c_char; 8]>(
                                        b"ARCCOS(\0",
                                    ))[0] as c_int
                                {
                                    -(1)
                                } else if (*yyvsp.offset(-2_isize)).astr[0] as c_int
                                    > (*::core::mem::transmute::<&[u8; 8], &[c_char; 8]>(
                                        b"ARCCOS(\0",
                                    ))[0] as c_int
                                {
                                    1
                                } else {
                                    strcmp(
                                        ((*yyvsp.offset(-2_isize)).astr).as_mut_ptr(),
                                        b"ARCCOS(\0" as *const u8 as *const c_char,
                                    )
                                }) == 0
                                    || (if ((*yyvsp.offset(-2_isize)).astr[0] as c_int)
                                        < (*::core::mem::transmute::<&[u8; 6], &[c_char; 6]>(
                                            b"ACOS(\0",
                                        ))[0] as c_int
                                    {
                                        -(1)
                                    } else if (*yyvsp.offset(-2_isize)).astr[0] as c_int
                                        > (*::core::mem::transmute::<&[u8; 6], &[c_char; 6]>(
                                            b"ACOS(\0",
                                        ))[0] as c_int
                                    {
                                        1
                                    } else {
                                        strcmp(
                                            ((*yyvsp.offset(-2_isize)).astr).as_mut_ptr(),
                                            b"ACOS(\0" as *const u8 as *const c_char,
                                        )
                                    }) == 0
                                {
                                    yyval.Node = New_Func(
                                        lParse,
                                        0,
                                        acos_fct,
                                        1,
                                        (*yyvsp.offset(-1_isize)).Node,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                    );
                                    current_block = 7600445499126923600;
                                } else if (if ((*yyvsp.offset(-2_isize)).astr[0] as c_int)
                                    < (*::core::mem::transmute::<&[u8; 8], &[c_char; 8]>(
                                        b"ARCTAN(\0",
                                    ))[0] as c_int
                                {
                                    -(1)
                                } else if (*yyvsp.offset(-2_isize)).astr[0] as c_int
                                    > (*::core::mem::transmute::<&[u8; 8], &[c_char; 8]>(
                                        b"ARCTAN(\0",
                                    ))[0] as c_int
                                {
                                    1
                                } else {
                                    strcmp(
                                        ((*yyvsp.offset(-2_isize)).astr).as_mut_ptr(),
                                        b"ARCTAN(\0" as *const u8 as *const c_char,
                                    )
                                }) == 0
                                    || (if ((*yyvsp.offset(-2_isize)).astr[0] as c_int)
                                        < (*::core::mem::transmute::<&[u8; 6], &[c_char; 6]>(
                                            b"ATAN(\0",
                                        ))[0] as c_int
                                    {
                                        -(1)
                                    } else if (*yyvsp.offset(-2_isize)).astr[0] as c_int
                                        > (*::core::mem::transmute::<&[u8; 6], &[c_char; 6]>(
                                            b"ATAN(\0",
                                        ))[0] as c_int
                                    {
                                        1
                                    } else {
                                        strcmp(
                                            ((*yyvsp.offset(-2_isize)).astr).as_mut_ptr(),
                                            b"ATAN(\0" as *const u8 as *const c_char,
                                        )
                                    }) == 0
                                {
                                    yyval.Node = New_Func(
                                        lParse,
                                        0,
                                        atan_fct,
                                        1,
                                        (*yyvsp.offset(-1_isize)).Node,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                    );
                                    current_block = 7600445499126923600;
                                } else if (if ((*yyvsp.offset(-2_isize)).astr[0] as c_int)
                                    < (*::core::mem::transmute::<&[u8; 6], &[c_char; 6]>(
                                        b"SINH(\0",
                                    ))[0] as c_int
                                {
                                    -(1)
                                } else if (*yyvsp.offset(-2_isize)).astr[0] as c_int
                                    > (*::core::mem::transmute::<&[u8; 6], &[c_char; 6]>(
                                        b"SINH(\0",
                                    ))[0] as c_int
                                {
                                    1
                                } else {
                                    strcmp(
                                        ((*yyvsp.offset(-2_isize)).astr).as_mut_ptr(),
                                        b"SINH(\0" as *const u8 as *const c_char,
                                    )
                                }) == 0
                                {
                                    yyval.Node = New_Func(
                                        lParse,
                                        0,
                                        sinh_fct,
                                        1,
                                        (*yyvsp.offset(-1_isize)).Node,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                    );
                                    current_block = 7600445499126923600;
                                } else if (if ((*yyvsp.offset(-2_isize)).astr[0] as c_int)
                                    < (*::core::mem::transmute::<&[u8; 6], &[c_char; 6]>(
                                        b"COSH(\0",
                                    ))[0] as c_int
                                {
                                    -(1)
                                } else if (*yyvsp.offset(-2_isize)).astr[0] as c_int
                                    > (*::core::mem::transmute::<&[u8; 6], &[c_char; 6]>(
                                        b"COSH(\0",
                                    ))[0] as c_int
                                {
                                    1
                                } else {
                                    strcmp(
                                        ((*yyvsp.offset(-2_isize)).astr).as_mut_ptr(),
                                        b"COSH(\0" as *const u8 as *const c_char,
                                    )
                                }) == 0
                                {
                                    yyval.Node = New_Func(
                                        lParse,
                                        0,
                                        cosh_fct,
                                        1,
                                        (*yyvsp.offset(-1_isize)).Node,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                    );
                                    current_block = 7600445499126923600;
                                } else if (if ((*yyvsp.offset(-2_isize)).astr[0] as c_int)
                                    < (*::core::mem::transmute::<&[u8; 6], &[c_char; 6]>(
                                        b"TANH(\0",
                                    ))[0] as c_int
                                {
                                    -(1)
                                } else if (*yyvsp.offset(-2_isize)).astr[0] as c_int
                                    > (*::core::mem::transmute::<&[u8; 6], &[c_char; 6]>(
                                        b"TANH(\0",
                                    ))[0] as c_int
                                {
                                    1
                                } else {
                                    strcmp(
                                        ((*yyvsp.offset(-2_isize)).astr).as_mut_ptr(),
                                        b"TANH(\0" as *const u8 as *const c_char,
                                    )
                                }) == 0
                                {
                                    yyval.Node = New_Func(
                                        lParse,
                                        0,
                                        tanh_fct,
                                        1,
                                        (*yyvsp.offset(-1_isize)).Node,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                    );
                                    current_block = 7600445499126923600;
                                } else if (if ((*yyvsp.offset(-2_isize)).astr[0] as c_int)
                                    < (*::core::mem::transmute::<&[u8; 5], &[c_char; 5]>(b"EXP(\0"))
                                        [0] as c_int
                                {
                                    -(1)
                                } else if (*yyvsp.offset(-2_isize)).astr[0] as c_int
                                    > (*::core::mem::transmute::<&[u8; 5], &[c_char; 5]>(b"EXP(\0"))
                                        [0] as c_int
                                {
                                    1
                                } else {
                                    strcmp(
                                        ((*yyvsp.offset(-2_isize)).astr).as_mut_ptr(),
                                        b"EXP(\0" as *const u8 as *const c_char,
                                    )
                                }) == 0
                                {
                                    yyval.Node = New_Func(
                                        lParse,
                                        0,
                                        exp_fct,
                                        1,
                                        (*yyvsp.offset(-1_isize)).Node,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                    );
                                    current_block = 7600445499126923600;
                                } else if (if ((*yyvsp.offset(-2_isize)).astr[0] as c_int)
                                    < (*::core::mem::transmute::<&[u8; 5], &[c_char; 5]>(b"LOG(\0"))
                                        [0] as c_int
                                {
                                    -(1)
                                } else if (*yyvsp.offset(-2_isize)).astr[0] as c_int
                                    > (*::core::mem::transmute::<&[u8; 5], &[c_char; 5]>(b"LOG(\0"))
                                        [0] as c_int
                                {
                                    1
                                } else {
                                    strcmp(
                                        ((*yyvsp.offset(-2_isize)).astr).as_mut_ptr(),
                                        b"LOG(\0" as *const u8 as *const c_char,
                                    )
                                }) == 0
                                {
                                    yyval.Node = New_Func(
                                        lParse,
                                        0,
                                        log_fct,
                                        1,
                                        (*yyvsp.offset(-1_isize)).Node,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                    );
                                    current_block = 7600445499126923600;
                                } else if (if ((*yyvsp.offset(-2_isize)).astr[0] as c_int)
                                    < (*::core::mem::transmute::<&[u8; 7], &[c_char; 7]>(
                                        b"LOG10(\0",
                                    ))[0] as c_int
                                {
                                    -(1)
                                } else if (*yyvsp.offset(-2_isize)).astr[0] as c_int
                                    > (*::core::mem::transmute::<&[u8; 7], &[c_char; 7]>(
                                        b"LOG10(\0",
                                    ))[0] as c_int
                                {
                                    1
                                } else {
                                    strcmp(
                                        ((*yyvsp.offset(-2_isize)).astr).as_mut_ptr(),
                                        b"LOG10(\0" as *const u8 as *const c_char,
                                    )
                                }) == 0
                                {
                                    yyval.Node = New_Func(
                                        lParse,
                                        0,
                                        log10_fct,
                                        1,
                                        (*yyvsp.offset(-1_isize)).Node,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                    );
                                    current_block = 7600445499126923600;
                                } else if (if ((*yyvsp.offset(-2_isize)).astr[0] as c_int)
                                    < (*::core::mem::transmute::<&[u8; 6], &[c_char; 6]>(
                                        b"SQRT(\0",
                                    ))[0] as c_int
                                {
                                    -(1)
                                } else if (*yyvsp.offset(-2_isize)).astr[0] as c_int
                                    > (*::core::mem::transmute::<&[u8; 6], &[c_char; 6]>(
                                        b"SQRT(\0",
                                    ))[0] as c_int
                                {
                                    1
                                } else {
                                    strcmp(
                                        ((*yyvsp.offset(-2_isize)).astr).as_mut_ptr(),
                                        b"SQRT(\0" as *const u8 as *const c_char,
                                    )
                                }) == 0
                                {
                                    yyval.Node = New_Func(
                                        lParse,
                                        0,
                                        sqrt_fct,
                                        1,
                                        (*yyvsp.offset(-1_isize)).Node,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                    );
                                    current_block = 7600445499126923600;
                                } else if (if ((*yyvsp.offset(-2_isize)).astr[0] as c_int)
                                    < (*::core::mem::transmute::<&[u8; 7], &[c_char; 7]>(
                                        b"ROUND(\0",
                                    ))[0] as c_int
                                {
                                    -(1)
                                } else if (*yyvsp.offset(-2_isize)).astr[0] as c_int
                                    > (*::core::mem::transmute::<&[u8; 7], &[c_char; 7]>(
                                        b"ROUND(\0",
                                    ))[0] as c_int
                                {
                                    1
                                } else {
                                    strcmp(
                                        ((*yyvsp.offset(-2_isize)).astr).as_mut_ptr(),
                                        b"ROUND(\0" as *const u8 as *const c_char,
                                    )
                                }) == 0
                                {
                                    yyval.Node = New_Func(
                                        lParse,
                                        0,
                                        round_fct,
                                        1,
                                        (*yyvsp.offset(-1_isize)).Node,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                    );
                                    current_block = 7600445499126923600;
                                } else if (if ((*yyvsp.offset(-2_isize)).astr[0] as c_int)
                                    < (*::core::mem::transmute::<&[u8; 7], &[c_char; 7]>(
                                        b"FLOOR(\0",
                                    ))[0] as c_int
                                {
                                    -(1)
                                } else if (*yyvsp.offset(-2_isize)).astr[0] as c_int
                                    > (*::core::mem::transmute::<&[u8; 7], &[c_char; 7]>(
                                        b"FLOOR(\0",
                                    ))[0] as c_int
                                {
                                    1
                                } else {
                                    strcmp(
                                        ((*yyvsp.offset(-2_isize)).astr).as_mut_ptr(),
                                        b"FLOOR(\0" as *const u8 as *const c_char,
                                    )
                                }) == 0
                                {
                                    yyval.Node = New_Func(
                                        lParse,
                                        0,
                                        floor_fct,
                                        1,
                                        (*yyvsp.offset(-1_isize)).Node,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                    );
                                    current_block = 7600445499126923600;
                                } else if (if ((*yyvsp.offset(-2_isize)).astr[0] as c_int)
                                    < (*::core::mem::transmute::<&[u8; 6], &[c_char; 6]>(
                                        b"CEIL(\0",
                                    ))[0] as c_int
                                {
                                    -(1)
                                } else if (*yyvsp.offset(-2_isize)).astr[0] as c_int
                                    > (*::core::mem::transmute::<&[u8; 6], &[c_char; 6]>(
                                        b"CEIL(\0",
                                    ))[0] as c_int
                                {
                                    1
                                } else {
                                    strcmp(
                                        ((*yyvsp.offset(-2_isize)).astr).as_mut_ptr(),
                                        b"CEIL(\0" as *const u8 as *const c_char,
                                    )
                                }) == 0
                                {
                                    yyval.Node = New_Func(
                                        lParse,
                                        0,
                                        ceil_fct,
                                        1,
                                        (*yyvsp.offset(-1_isize)).Node,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                    );
                                    current_block = 7600445499126923600;
                                } else if (if ((*yyvsp.offset(-2_isize)).astr[0] as c_int)
                                    < (*::core::mem::transmute::<&[u8; 9], &[c_char; 9]>(
                                        b"RANDOMP(\0",
                                    ))[0] as c_int
                                {
                                    -(1)
                                } else if (*yyvsp.offset(-2_isize)).astr[0] as c_int
                                    > (*::core::mem::transmute::<&[u8; 9], &[c_char; 9]>(
                                        b"RANDOMP(\0",
                                    ))[0] as c_int
                                {
                                    1
                                } else {
                                    strcmp(
                                        ((*yyvsp.offset(-2_isize)).astr).as_mut_ptr(),
                                        b"RANDOMP(\0" as *const u8 as *const c_char,
                                    )
                                }) == 0
                                {
                                    yyval.Node = New_Func(
                                        lParse,
                                        0,
                                        poirnd_fct,
                                        1,
                                        (*yyvsp.offset(-1_isize)).Node,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                    );
                                    (*((*lParse).Nodes).offset(yyval.Node as isize)).ntype =
                                        fits_parser_yytokentype::LONG as c_int;
                                    current_block = 7600445499126923600;
                                } else {
                                    fits_parser_yyerror(
                                        scanner,
                                        lParse,
                                        b"Function(expr) not supported\0" as *const u8
                                            as *const c_char
                                            as *mut c_char,
                                    );
                                    current_block = 4830776507462815627;
                                }
                            }
                            match current_block {
                                4830776507462815627 => {}
                                _ => {
                                    if yyval.Node < 0 {
                                        current_block = 4830776507462815627;
                                    } else {
                                        current_block = 17353983478346836848;
                                    }
                                }
                            }
                        }
                        60 => {
                            if (if ((*yyvsp.offset(-4_isize)).astr[0] as c_int)
                                < (*::core::mem::transmute::<&[u8; 8], &[c_char; 8]>(b"STRSTR(\0"))
                                    [0] as c_int
                            {
                                -(1)
                            } else if (*yyvsp.offset(-4_isize)).astr[0] as c_int
                                > (*::core::mem::transmute::<&[u8; 8], &[c_char; 8]>(b"STRSTR(\0"))
                                    [0] as c_int
                            {
                                1
                            } else {
                                strcmp(
                                    ((*yyvsp.offset(-4_isize)).astr).as_mut_ptr(),
                                    b"STRSTR(\0" as *const u8 as *const c_char,
                                )
                            }) == 0
                            {
                                yyval.Node = New_Func(
                                    lParse,
                                    fits_parser_yytokentype::LONG as c_int,
                                    strpos_fct,
                                    2 as c_int,
                                    (*yyvsp.offset(-3_isize)).Node,
                                    (*yyvsp.offset(-1_isize)).Node,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                );
                                if yyval.Node < 0 {
                                    current_block = 4830776507462815627;
                                } else {
                                    current_block = 12473999534294722905;
                                }
                            } else {
                                current_block = 12473999534294722905;
                            }
                            match current_block {
                                4830776507462815627 => {}
                                _ => {
                                    current_block = 17353983478346836848;
                                }
                            }
                        }
                        61 => {
                            if (if ((*yyvsp.offset(-4_isize)).astr[0] as c_int)
                                < (*::core::mem::transmute::<&[u8; 9], &[c_char; 9]>(b"DEFNULL(\0"))
                                    [0] as c_int
                            {
                                -(1)
                            } else if (*yyvsp.offset(-4_isize)).astr[0] as c_int
                                > (*::core::mem::transmute::<&[u8; 9], &[c_char; 9]>(b"DEFNULL(\0"))
                                    [0] as c_int
                            {
                                1
                            } else {
                                strcmp(
                                    ((*yyvsp.offset(-4_isize)).astr).as_mut_ptr(),
                                    b"DEFNULL(\0" as *const u8 as *const c_char,
                                )
                            }) == 0
                            {
                                if (*((*lParse).Nodes)
                                    .offset((*yyvsp.offset(-3_isize)).Node as isize))
                                .value
                                .nelem
                                    >= (*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-1_isize)).Node as isize))
                                    .value
                                    .nelem
                                    && Test_Dims(
                                        lParse,
                                        (*yyvsp.offset(-3_isize)).Node,
                                        (*yyvsp.offset(-1_isize)).Node,
                                    ) != 0
                                {
                                    if (*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-3_isize)).Node as isize))
                                    .ntype
                                        > (*((*lParse).Nodes)
                                            .offset((*yyvsp.offset(-1_isize)).Node as isize))
                                        .ntype
                                    {
                                        (*yyvsp.offset(-1_isize)).Node = New_Unary(
                                            lParse,
                                            (*((*lParse).Nodes)
                                                .offset((*yyvsp.offset(-3_isize)).Node as isize))
                                            .ntype,
                                            0,
                                            (*yyvsp.offset(-1_isize)).Node,
                                        );
                                    } else if (*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-3_isize)).Node as isize))
                                    .ntype
                                        < (*((*lParse).Nodes)
                                            .offset((*yyvsp.offset(-1_isize)).Node as isize))
                                        .ntype
                                    {
                                        (*yyvsp.offset(-3_isize)).Node = New_Unary(
                                            lParse,
                                            (*((*lParse).Nodes)
                                                .offset((*yyvsp.offset(-1_isize)).Node as isize))
                                            .ntype,
                                            0,
                                            (*yyvsp.offset(-3_isize)).Node,
                                        );
                                    }
                                    yyval.Node = New_Func(
                                        lParse,
                                        0,
                                        defnull_fct,
                                        2 as c_int,
                                        (*yyvsp.offset(-3_isize)).Node,
                                        (*yyvsp.offset(-1_isize)).Node,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                    );
                                    if yyval.Node < 0 {
                                        current_block = 4830776507462815627;
                                    } else {
                                        current_block = 9966817879908499150;
                                    }
                                } else {
                                    fits_parser_yyerror(
                                        scanner,
                                        lParse,
                                        b"Dimensions of DEFNULL arguments are not compatible\0"
                                            as *const u8
                                            as *const c_char
                                            as *mut c_char,
                                    );
                                    current_block = 4830776507462815627;
                                }
                            } else if (if ((*yyvsp.offset(-4_isize)).astr[0] as c_int)
                                < (*::core::mem::transmute::<&[u8; 9], &[c_char; 9]>(b"ARCTAN2(\0"))
                                    [0] as c_int
                            {
                                -(1)
                            } else if (*yyvsp.offset(-4_isize)).astr[0] as c_int
                                > (*::core::mem::transmute::<&[u8; 9], &[c_char; 9]>(b"ARCTAN2(\0"))
                                    [0] as c_int
                            {
                                1
                            } else {
                                strcmp(
                                    ((*yyvsp.offset(-4_isize)).astr).as_mut_ptr(),
                                    b"ARCTAN2(\0" as *const u8 as *const c_char,
                                )
                            }) == 0
                            {
                                if (*((*lParse).Nodes)
                                    .offset((*yyvsp.offset(-3_isize)).Node as isize))
                                .ntype
                                    != fits_parser_yytokentype::DOUBLE as c_int
                                {
                                    (*yyvsp.offset(-3_isize)).Node = New_Unary(
                                        lParse,
                                        fits_parser_yytokentype::DOUBLE as c_int,
                                        0,
                                        (*yyvsp.offset(-3_isize)).Node,
                                    );
                                }
                                if (*((*lParse).Nodes)
                                    .offset((*yyvsp.offset(-1_isize)).Node as isize))
                                .ntype
                                    != fits_parser_yytokentype::DOUBLE as c_int
                                {
                                    (*yyvsp.offset(-1_isize)).Node = New_Unary(
                                        lParse,
                                        fits_parser_yytokentype::DOUBLE as c_int,
                                        0,
                                        (*yyvsp.offset(-1_isize)).Node,
                                    );
                                }
                                if Test_Dims(
                                    lParse,
                                    (*yyvsp.offset(-3_isize)).Node,
                                    (*yyvsp.offset(-1_isize)).Node,
                                ) != 0
                                {
                                    yyval.Node = New_Func(
                                        lParse,
                                        0,
                                        atan2_fct,
                                        2 as c_int,
                                        (*yyvsp.offset(-3_isize)).Node,
                                        (*yyvsp.offset(-1_isize)).Node,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                    );
                                    if yyval.Node < 0 {
                                        current_block = 4830776507462815627;
                                    } else {
                                        if (*((*lParse).Nodes)
                                            .offset((*yyvsp.offset(-3_isize)).Node as isize))
                                        .value
                                        .nelem
                                            < (*((*lParse).Nodes)
                                                .offset((*yyvsp.offset(-1_isize)).Node as isize))
                                            .value
                                            .nelem
                                        {
                                            Copy_Dims(
                                                lParse,
                                                yyval.Node,
                                                (*yyvsp.offset(-1_isize)).Node,
                                            );
                                        }
                                        current_block = 9966817879908499150;
                                    }
                                } else {
                                    fits_parser_yyerror(
                                        scanner,
                                        lParse,
                                        b"Dimensions of arctan2 arguments are not compatible\0"
                                            as *const u8
                                            as *const c_char
                                            as *mut c_char,
                                    );
                                    current_block = 4830776507462815627;
                                }
                            } else if (if ((*yyvsp.offset(-4_isize)).astr[0] as c_int)
                                < (*::core::mem::transmute::<&[u8; 5], &[c_char; 5]>(b"MIN(\0"))[0]
                                    as c_int
                            {
                                -(1)
                            } else if (*yyvsp.offset(-4_isize)).astr[0] as c_int
                                > (*::core::mem::transmute::<&[u8; 5], &[c_char; 5]>(b"MIN(\0"))[0]
                                    as c_int
                            {
                                1
                            } else {
                                strcmp(
                                    ((*yyvsp.offset(-4_isize)).astr).as_mut_ptr(),
                                    b"MIN(\0" as *const u8 as *const c_char,
                                )
                            }) == 0
                            {
                                if (*((*lParse).Nodes)
                                    .offset((*yyvsp.offset(-3_isize)).Node as isize))
                                .ntype
                                    > (*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-1_isize)).Node as isize))
                                    .ntype
                                {
                                    (*yyvsp.offset(-1_isize)).Node = New_Unary(
                                        lParse,
                                        (*((*lParse).Nodes)
                                            .offset((*yyvsp.offset(-3_isize)).Node as isize))
                                        .ntype,
                                        0,
                                        (*yyvsp.offset(-1_isize)).Node,
                                    );
                                } else if (*((*lParse).Nodes)
                                    .offset((*yyvsp.offset(-3_isize)).Node as isize))
                                .ntype
                                    < (*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-1_isize)).Node as isize))
                                    .ntype
                                {
                                    (*yyvsp.offset(-3_isize)).Node = New_Unary(
                                        lParse,
                                        (*((*lParse).Nodes)
                                            .offset((*yyvsp.offset(-1_isize)).Node as isize))
                                        .ntype,
                                        0,
                                        (*yyvsp.offset(-3_isize)).Node,
                                    );
                                }
                                if Test_Dims(
                                    lParse,
                                    (*yyvsp.offset(-3_isize)).Node,
                                    (*yyvsp.offset(-1_isize)).Node,
                                ) != 0
                                {
                                    yyval.Node = New_Func(
                                        lParse,
                                        0,
                                        min2_fct,
                                        2 as c_int,
                                        (*yyvsp.offset(-3_isize)).Node,
                                        (*yyvsp.offset(-1_isize)).Node,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                    );
                                    if yyval.Node < 0 {
                                        current_block = 4830776507462815627;
                                    } else {
                                        if (*((*lParse).Nodes)
                                            .offset((*yyvsp.offset(-3_isize)).Node as isize))
                                        .value
                                        .nelem
                                            < (*((*lParse).Nodes)
                                                .offset((*yyvsp.offset(-1_isize)).Node as isize))
                                            .value
                                            .nelem
                                        {
                                            Copy_Dims(
                                                lParse,
                                                yyval.Node,
                                                (*yyvsp.offset(-1_isize)).Node,
                                            );
                                        }
                                        current_block = 9966817879908499150;
                                    }
                                } else {
                                    fits_parser_yyerror(
                                        scanner,
                                        lParse,
                                        b"Dimensions of min(a,b) arguments are not compatible\0"
                                            as *const u8
                                            as *const c_char
                                            as *mut c_char,
                                    );
                                    current_block = 4830776507462815627;
                                }
                            } else if (if ((*yyvsp.offset(-4_isize)).astr[0] as c_int)
                                < (*::core::mem::transmute::<&[u8; 5], &[c_char; 5]>(b"MAX(\0"))[0]
                                    as c_int
                            {
                                -(1)
                            } else if (*yyvsp.offset(-4_isize)).astr[0] as c_int
                                > (*::core::mem::transmute::<&[u8; 5], &[c_char; 5]>(b"MAX(\0"))[0]
                                    as c_int
                            {
                                1
                            } else {
                                strcmp(
                                    ((*yyvsp.offset(-4_isize)).astr).as_mut_ptr(),
                                    b"MAX(\0" as *const u8 as *const c_char,
                                )
                            }) == 0
                            {
                                if (*((*lParse).Nodes)
                                    .offset((*yyvsp.offset(-3_isize)).Node as isize))
                                .ntype
                                    > (*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-1_isize)).Node as isize))
                                    .ntype
                                {
                                    (*yyvsp.offset(-1_isize)).Node = New_Unary(
                                        lParse,
                                        (*((*lParse).Nodes)
                                            .offset((*yyvsp.offset(-3_isize)).Node as isize))
                                        .ntype,
                                        0,
                                        (*yyvsp.offset(-1_isize)).Node,
                                    );
                                } else if (*((*lParse).Nodes)
                                    .offset((*yyvsp.offset(-3_isize)).Node as isize))
                                .ntype
                                    < (*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-1_isize)).Node as isize))
                                    .ntype
                                {
                                    (*yyvsp.offset(-3_isize)).Node = New_Unary(
                                        lParse,
                                        (*((*lParse).Nodes)
                                            .offset((*yyvsp.offset(-1_isize)).Node as isize))
                                        .ntype,
                                        0,
                                        (*yyvsp.offset(-3_isize)).Node,
                                    );
                                }
                                if Test_Dims(
                                    lParse,
                                    (*yyvsp.offset(-3_isize)).Node,
                                    (*yyvsp.offset(-1_isize)).Node,
                                ) != 0
                                {
                                    yyval.Node = New_Func(
                                        lParse,
                                        0,
                                        max2_fct,
                                        2 as c_int,
                                        (*yyvsp.offset(-3_isize)).Node,
                                        (*yyvsp.offset(-1_isize)).Node,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                    );
                                    if yyval.Node < 0 {
                                        current_block = 4830776507462815627;
                                    } else {
                                        if (*((*lParse).Nodes)
                                            .offset((*yyvsp.offset(-3_isize)).Node as isize))
                                        .value
                                        .nelem
                                            < (*((*lParse).Nodes)
                                                .offset((*yyvsp.offset(-1_isize)).Node as isize))
                                            .value
                                            .nelem
                                        {
                                            Copy_Dims(
                                                lParse,
                                                yyval.Node,
                                                (*yyvsp.offset(-1_isize)).Node,
                                            );
                                        }
                                        current_block = 9966817879908499150;
                                    }
                                } else {
                                    fits_parser_yyerror(
                                        scanner,
                                        lParse,
                                        b"Dimensions of max(a,b) arguments are not compatible\0"
                                            as *const u8
                                            as *const c_char
                                            as *mut c_char,
                                    );
                                    current_block = 4830776507462815627;
                                }
                            } else if (if ((*yyvsp.offset(-4_isize)).astr[0] as c_int)
                                < (*::core::mem::transmute::<&[u8; 9], &[c_char; 9]>(b"SETNULL(\0"))
                                    [0] as c_int
                            {
                                -(1)
                            } else if (*yyvsp.offset(-4_isize)).astr[0] as c_int
                                > (*::core::mem::transmute::<&[u8; 9], &[c_char; 9]>(b"SETNULL(\0"))
                                    [0] as c_int
                            {
                                1
                            } else {
                                strcmp(
                                    ((*yyvsp.offset(-4_isize)).astr).as_mut_ptr(),
                                    b"SETNULL(\0" as *const u8 as *const c_char,
                                )
                            }) == 0
                            {
                                if (*((*lParse).Nodes)
                                    .offset((*yyvsp.offset(-3_isize)).Node as isize))
                                .operation
                                    != -(1000 as c_int)
                                    || (*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-3_isize)).Node as isize))
                                    .value
                                    .nelem
                                        != 1
                                {
                                    fits_parser_yyerror(
                                        scanner,
                                        lParse,
                                        b"SETNULL first argument must be a scalar constant\0"
                                            as *const u8
                                            as *const c_char
                                            as *mut c_char,
                                    );
                                    current_block = 4830776507462815627;
                                } else {
                                    if (*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-3_isize)).Node as isize))
                                    .ntype
                                        != (*((*lParse).Nodes)
                                            .offset((*yyvsp.offset(-1_isize)).Node as isize))
                                        .ntype
                                    {
                                        (*yyvsp.offset(-3_isize)).Node = New_Unary(
                                            lParse,
                                            (*((*lParse).Nodes)
                                                .offset((*yyvsp.offset(-1_isize)).Node as isize))
                                            .ntype,
                                            0,
                                            (*yyvsp.offset(-3_isize)).Node,
                                        );
                                    }
                                    yyval.Node = New_Func(
                                        lParse,
                                        0,
                                        setnull_fct,
                                        2 as c_int,
                                        (*yyvsp.offset(-1_isize)).Node,
                                        (*yyvsp.offset(-3_isize)).Node,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                    );
                                    current_block = 9966817879908499150;
                                }
                            } else if (if ((*yyvsp.offset(-4_isize)).astr[0] as c_int)
                                < (*::core::mem::transmute::<&[u8; 10], &[c_char; 10]>(
                                    b"AXISELEM(\0",
                                ))[0] as c_int
                            {
                                -(1)
                            } else if (*yyvsp.offset(-4_isize)).astr[0] as c_int
                                > (*::core::mem::transmute::<&[u8; 10], &[c_char; 10]>(
                                    b"AXISELEM(\0",
                                ))[0] as c_int
                            {
                                1
                            } else {
                                strcmp(
                                    ((*yyvsp.offset(-4_isize)).astr).as_mut_ptr(),
                                    b"AXISELEM(\0" as *const u8 as *const c_char,
                                )
                            }) == 0
                            {
                                if (*((*lParse).Nodes)
                                    .offset((*yyvsp.offset(-1_isize)).Node as isize))
                                .operation
                                    != -(1000 as c_int)
                                    || (*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-1_isize)).Node as isize))
                                    .value
                                    .nelem
                                        != 1
                                {
                                    fits_parser_yyerror(
                                        scanner,
                                        lParse,
                                        b"AXISELEM second argument must be a scalar constant\0"
                                            as *const u8
                                            as *const c_char
                                            as *mut c_char,
                                    );
                                    current_block = 4830776507462815627;
                                } else if (*((*lParse).Nodes)
                                    .offset((*yyvsp.offset(-3_isize)).Node as isize))
                                .operation
                                    == -(1000 as c_int)
                                {
                                    let mut one_3: c_long = 1;
                                    yyval.Node = New_Const(
                                        lParse,
                                        fits_parser_yytokentype::LONG as c_int,
                                        &mut one_3 as *mut c_long as *mut c_void,
                                        ::core::mem::size_of::<c_long>() as c_ulong as c_long,
                                    );
                                    current_block = 9966817879908499150;
                                } else {
                                    if (*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-1_isize)).Node as isize))
                                    .ntype
                                        != fits_parser_yytokentype::LONG as c_int
                                    {
                                        (*yyvsp.offset(-1_isize)).Node = New_Unary(
                                            lParse,
                                            fits_parser_yytokentype::LONG as c_int,
                                            0,
                                            (*yyvsp.offset(-1_isize)).Node,
                                        );
                                    }
                                    yyval.Node = New_Func(
                                        lParse,
                                        0,
                                        axiselem_fct,
                                        2 as c_int,
                                        (*yyvsp.offset(-3_isize)).Node,
                                        (*yyvsp.offset(-1_isize)).Node,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                    );
                                    if yyval.Node < 0 {
                                        current_block = 4830776507462815627;
                                    } else {
                                        (*((*lParse).Nodes).offset(yyval.Node as isize)).ntype =
                                            fits_parser_yytokentype::LONG as c_int;
                                        current_block = 9966817879908499150;
                                    }
                                }
                            } else if (if ((*yyvsp.offset(-4_isize)).astr[0] as c_int)
                                < (*::core::mem::transmute::<&[u8; 7], &[c_char; 7]>(b"NAXES(\0"))
                                    [0] as c_int
                            {
                                -(1)
                            } else if (*yyvsp.offset(-4_isize)).astr[0] as c_int
                                > (*::core::mem::transmute::<&[u8; 7], &[c_char; 7]>(b"NAXES(\0"))
                                    [0] as c_int
                            {
                                1
                            } else {
                                strcmp(
                                    ((*yyvsp.offset(-4_isize)).astr).as_mut_ptr(),
                                    b"NAXES(\0" as *const u8 as *const c_char,
                                )
                            }) == 0
                            {
                                if (*((*lParse).Nodes)
                                    .offset((*yyvsp.offset(-1_isize)).Node as isize))
                                .operation
                                    != -(1000 as c_int)
                                    || (*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-1_isize)).Node as isize))
                                    .value
                                    .nelem
                                        != 1
                                {
                                    fits_parser_yyerror(
                                        scanner,
                                        lParse,
                                        b"NAXES second argument must be a scalar constant\0"
                                            as *const u8
                                            as *const c_char
                                            as *mut c_char,
                                    );
                                    current_block = 4830776507462815627;
                                } else if (*((*lParse).Nodes)
                                    .offset((*yyvsp.offset(-3_isize)).Node as isize))
                                .operation
                                    == -(1000 as c_int)
                                {
                                    let mut one_4: c_long = 1;
                                    yyval.Node = New_Const(
                                        lParse,
                                        fits_parser_yytokentype::LONG as c_int,
                                        &mut one_4 as *mut c_long as *mut c_void,
                                        ::core::mem::size_of::<c_long>() as c_ulong as c_long,
                                    );
                                    current_block = 9966817879908499150;
                                } else {
                                    let mut iaxis_0: c_long = 0;
                                    let mut naxis_1: c_int = 0;
                                    if (*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-1_isize)).Node as isize))
                                    .ntype
                                        != fits_parser_yytokentype::LONG as c_int
                                    {
                                        (*yyvsp.offset(-1_isize)).Node = New_Unary(
                                            lParse,
                                            fits_parser_yytokentype::LONG as c_int,
                                            0,
                                            (*yyvsp.offset(-1_isize)).Node,
                                        );
                                    }
                                    iaxis_0 = (*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-1_isize)).Node as isize))
                                    .value
                                    .data
                                    .lng;
                                    naxis_1 = (*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-3_isize)).Node as isize))
                                    .value
                                    .naxis;
                                    if iaxis_0 == 0 {
                                        iaxis_0 = naxis_1 as c_long;
                                    } else if iaxis_0 <= naxis_1 as c_long {
                                        iaxis_0 = (*((*lParse).Nodes)
                                            .offset((*yyvsp.offset(-3_isize)).Node as isize))
                                        .value
                                        .naxes[(iaxis_0 - 1) as usize];
                                    } else {
                                        iaxis_0 = 1;
                                    }
                                    yyval.Node = New_Const(
                                        lParse,
                                        fits_parser_yytokentype::LONG as c_int,
                                        &mut iaxis_0 as *mut c_long as *mut c_void,
                                        ::core::mem::size_of::<c_long>() as c_ulong as c_long,
                                    );
                                    if yyval.Node < 0 {
                                        current_block = 4830776507462815627;
                                    } else {
                                        current_block = 9966817879908499150;
                                    }
                                }
                            } else if (if ((*yyvsp.offset(-4_isize)).astr[0] as c_int)
                                < (*::core::mem::transmute::<&[u8; 7], &[c_char; 7]>(b"ARRAY(\0"))
                                    [0] as c_int
                            {
                                -(1)
                            } else if (*yyvsp.offset(-4_isize)).astr[0] as c_int
                                > (*::core::mem::transmute::<&[u8; 7], &[c_char; 7]>(b"ARRAY(\0"))
                                    [0] as c_int
                            {
                                1
                            } else {
                                strcmp(
                                    ((*yyvsp.offset(-4_isize)).astr).as_mut_ptr(),
                                    b"ARRAY(\0" as *const u8 as *const c_char,
                                )
                            }) == 0
                            {
                                yyval.Node = New_Array(
                                    lParse,
                                    (*yyvsp.offset(-3_isize)).Node,
                                    (*yyvsp.offset(-1_isize)).Node,
                                );
                                if yyval.Node < 0 {
                                    current_block = 4830776507462815627;
                                } else {
                                    current_block = 9966817879908499150;
                                }
                            } else {
                                fits_parser_yyerror(
                                    scanner,
                                    lParse,
                                    b"Function(expr,expr) not supported\0" as *const u8
                                        as *const c_char
                                        as *mut c_char,
                                );
                                current_block = 4830776507462815627;
                            }
                            match current_block {
                                4830776507462815627 => {}
                                _ => {
                                    current_block = 17353983478346836848;
                                }
                            }
                        }
                        62 => {
                            if (if ((*yyvsp.offset(-(8 as c_int) as isize)).astr[0] as c_int)
                                < (*::core::mem::transmute::<&[u8; 8], &[c_char; 8]>(b"ANGSEP(\0"))
                                    [0] as c_int
                            {
                                -(1)
                            } else if (*yyvsp.offset(-(8 as c_int) as isize)).astr[0] as c_int
                                > (*::core::mem::transmute::<&[u8; 8], &[c_char; 8]>(b"ANGSEP(\0"))
                                    [0] as c_int
                            {
                                1
                            } else {
                                strcmp(
                                    ((*yyvsp.offset(-(8 as c_int) as isize)).astr).as_mut_ptr(),
                                    b"ANGSEP(\0" as *const u8 as *const c_char,
                                )
                            }) == 0
                            {
                                if (*((*lParse).Nodes)
                                    .offset((*yyvsp.offset(-(7 as c_int) as isize)).Node as isize))
                                .ntype
                                    != fits_parser_yytokentype::DOUBLE as c_int
                                {
                                    (*yyvsp.offset(-(7 as c_int) as isize)).Node = New_Unary(
                                        lParse,
                                        fits_parser_yytokentype::DOUBLE as c_int,
                                        0,
                                        (*yyvsp.offset(-(7 as c_int) as isize)).Node,
                                    );
                                }
                                if (*((*lParse).Nodes).offset((*yyvsp.offset(-(5))).Node as isize))
                                    .ntype
                                    != fits_parser_yytokentype::DOUBLE as c_int
                                {
                                    (*yyvsp.offset(-(5))).Node = New_Unary(
                                        lParse,
                                        fits_parser_yytokentype::DOUBLE as c_int,
                                        0,
                                        (*yyvsp.offset(-(5))).Node,
                                    );
                                }
                                if (*((*lParse).Nodes)
                                    .offset((*yyvsp.offset(-3_isize)).Node as isize))
                                .ntype
                                    != fits_parser_yytokentype::DOUBLE as c_int
                                {
                                    (*yyvsp.offset(-3_isize)).Node = New_Unary(
                                        lParse,
                                        fits_parser_yytokentype::DOUBLE as c_int,
                                        0,
                                        (*yyvsp.offset(-3_isize)).Node,
                                    );
                                }
                                if (*((*lParse).Nodes)
                                    .offset((*yyvsp.offset(-1_isize)).Node as isize))
                                .ntype
                                    != fits_parser_yytokentype::DOUBLE as c_int
                                {
                                    (*yyvsp.offset(-1_isize)).Node = New_Unary(
                                        lParse,
                                        fits_parser_yytokentype::DOUBLE as c_int,
                                        0,
                                        (*yyvsp.offset(-1_isize)).Node,
                                    );
                                }
                                if Test_Dims(
                                    lParse,
                                    (*yyvsp.offset(-(7 as c_int) as isize)).Node,
                                    (*yyvsp.offset(-(5))).Node,
                                ) != 0
                                    && Test_Dims(
                                        lParse,
                                        (*yyvsp.offset(-(5))).Node,
                                        (*yyvsp.offset(-3_isize)).Node,
                                    ) != 0
                                    && Test_Dims(
                                        lParse,
                                        (*yyvsp.offset(-3_isize)).Node,
                                        (*yyvsp.offset(-1_isize)).Node,
                                    ) != 0
                                {
                                    yyval.Node = New_Func(
                                        lParse,
                                        0,
                                        angsep_fct,
                                        4 as c_int,
                                        (*yyvsp.offset(-(7 as c_int) as isize)).Node,
                                        (*yyvsp.offset(-(5))).Node,
                                        (*yyvsp.offset(-3_isize)).Node,
                                        (*yyvsp.offset(-1_isize)).Node,
                                        0,
                                        0,
                                        0,
                                    );
                                    if yyval.Node < 0 {
                                        current_block = 4830776507462815627;
                                    } else {
                                        if (*((*lParse).Nodes).offset(
                                            (*yyvsp.offset(-(7 as c_int) as isize)).Node as isize,
                                        ))
                                        .value
                                        .nelem
                                            < (*((*lParse).Nodes)
                                                .offset((*yyvsp.offset(-(5))).Node as isize))
                                            .value
                                            .nelem
                                        {
                                            Copy_Dims(
                                                lParse,
                                                yyval.Node,
                                                (*yyvsp.offset(-(5))).Node,
                                            );
                                        }
                                        if (*((*lParse).Nodes)
                                            .offset((*yyvsp.offset(-(5))).Node as isize))
                                        .value
                                        .nelem
                                            < (*((*lParse).Nodes)
                                                .offset((*yyvsp.offset(-3_isize)).Node as isize))
                                            .value
                                            .nelem
                                        {
                                            Copy_Dims(
                                                lParse,
                                                yyval.Node,
                                                (*yyvsp.offset(-3_isize)).Node,
                                            );
                                        }
                                        if (*((*lParse).Nodes)
                                            .offset((*yyvsp.offset(-3_isize)).Node as isize))
                                        .value
                                        .nelem
                                            < (*((*lParse).Nodes)
                                                .offset((*yyvsp.offset(-1_isize)).Node as isize))
                                            .value
                                            .nelem
                                        {
                                            Copy_Dims(
                                                lParse,
                                                yyval.Node,
                                                (*yyvsp.offset(-1_isize)).Node,
                                            );
                                        }
                                        current_block = 17353983478346836848;
                                    }
                                } else {
                                    fits_parser_yyerror(
                                        scanner,
                                        lParse,
                                        b"Dimensions of ANGSEP arguments are not compatible\0"
                                            as *const u8
                                            as *const c_char
                                            as *mut c_char,
                                    );
                                    current_block = 4830776507462815627;
                                }
                            } else {
                                fits_parser_yyerror(
                                    scanner,
                                    lParse,
                                    b"Function(expr,expr,expr,expr) not supported\0" as *const u8
                                        as *const c_char
                                        as *mut c_char,
                                );
                                current_block = 4830776507462815627;
                            }
                        }
                        63 => {
                            yyval.Node = New_GTI(
                                lParse,
                                gtiover_fct,
                                ((*yyvsp.offset(-(5))).astr).as_mut_ptr(),
                                (*yyvsp.offset(-3_isize)).Node,
                                (*yyvsp.offset(-1_isize)).Node,
                                b"*START*\0" as *const u8 as *const c_char as *mut c_char,
                                b"*STOP*\0" as *const u8 as *const c_char as *mut c_char,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        64 => {
                            yyval.Node = New_GTI(
                                lParse,
                                gtiover_fct,
                                ((*yyvsp.offset(-(9 as c_int) as isize)).astr).as_mut_ptr(),
                                (*yyvsp.offset(-(7 as c_int) as isize)).Node,
                                (*yyvsp.offset(-(5))).Node,
                                ((*yyvsp.offset(-3_isize)).astr).as_mut_ptr(),
                                ((*yyvsp.offset(-1_isize)).astr).as_mut_ptr(),
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        65 => {
                            yyval.Node = New_Deref(
                                lParse,
                                (*yyvsp.offset(-3_isize)).Node,
                                1,
                                (*yyvsp.offset(-1_isize)).Node,
                                0,
                                0,
                                0,
                                0,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        66 => {
                            yyval.Node = New_Deref(
                                lParse,
                                (*yyvsp.offset(-(5))).Node,
                                2 as c_int,
                                (*yyvsp.offset(-3_isize)).Node,
                                (*yyvsp.offset(-1_isize)).Node,
                                0,
                                0,
                                0,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        67 => {
                            yyval.Node = New_Deref(
                                lParse,
                                (*yyvsp.offset(-(7 as c_int) as isize)).Node,
                                3 as c_int,
                                (*yyvsp.offset(-(5))).Node,
                                (*yyvsp.offset(-3_isize)).Node,
                                (*yyvsp.offset(-1_isize)).Node,
                                0,
                                0,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        68 => {
                            yyval.Node = New_Deref(
                                lParse,
                                (*yyvsp.offset(-(9 as c_int) as isize)).Node,
                                4 as c_int,
                                (*yyvsp.offset(-(7 as c_int) as isize)).Node,
                                (*yyvsp.offset(-(5))).Node,
                                (*yyvsp.offset(-3_isize)).Node,
                                (*yyvsp.offset(-1_isize)).Node,
                                0,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        69 => {
                            yyval.Node = New_Deref(
                                lParse,
                                (*yyvsp.offset(-(11 as c_int) as isize)).Node,
                                5 as c_int,
                                (*yyvsp.offset(-(9 as c_int) as isize)).Node,
                                (*yyvsp.offset(-(7 as c_int) as isize)).Node,
                                (*yyvsp.offset(-(5))).Node,
                                (*yyvsp.offset(-3_isize)).Node,
                                (*yyvsp.offset(-1_isize)).Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        70 => {
                            yyval.Node = New_Unary(
                                lParse,
                                fits_parser_yytokentype::LONG as c_int,
                                fits_parser_yytokentype::INTCAST as c_int,
                                (*yyvsp.offset(0)).Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        71 => {
                            yyval.Node = New_Unary(
                                lParse,
                                fits_parser_yytokentype::LONG as c_int,
                                fits_parser_yytokentype::INTCAST as c_int,
                                (*yyvsp.offset(0)).Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        72 => {
                            yyval.Node = New_Unary(
                                lParse,
                                fits_parser_yytokentype::DOUBLE as c_int,
                                fits_parser_yytokentype::FLTCAST as c_int,
                                (*yyvsp.offset(0)).Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        73 => {
                            yyval.Node = New_Unary(
                                lParse,
                                fits_parser_yytokentype::DOUBLE as c_int,
                                fits_parser_yytokentype::FLTCAST as c_int,
                                (*yyvsp.offset(0)).Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        74 => {
                            yyval.Node = New_Const(
                                lParse,
                                fits_parser_yytokentype::BOOLEAN as c_int,
                                &mut (*yyvsp.offset(0)).log as *mut c_char as *mut c_void,
                                ::core::mem::size_of::<c_char>() as c_ulong as c_long,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        75 => {
                            yyval.Node = New_Column(lParse, (*yyvsp.offset(0)).lng as c_int);
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        76 => {
                            if (*((*lParse).Nodes).offset((*yyvsp.offset(-1_isize)).Node as isize))
                                .ntype
                                != fits_parser_yytokentype::LONG as c_int
                                || (*((*lParse).Nodes)
                                    .offset((*yyvsp.offset(-1_isize)).Node as isize))
                                .operation
                                    != -(1000 as c_int)
                            {
                                fits_parser_yyerror(
                                    scanner,
                                    lParse,
                                    b"Offset argument must be a constant integer\0" as *const u8
                                        as *const c_char
                                        as *mut c_char,
                                );
                                current_block = 4830776507462815627;
                            } else {
                                yyval.Node = New_Offset(
                                    lParse,
                                    (*yyvsp.offset(-3_isize)).lng as c_int,
                                    (*yyvsp.offset(-1_isize)).Node,
                                );
                                if yyval.Node < 0 {
                                    current_block = 4830776507462815627;
                                } else {
                                    current_block = 17353983478346836848;
                                }
                            }
                        }
                        77 => {
                            yyval.Node = New_BinOp(
                                lParse,
                                fits_parser_yytokentype::BOOLEAN as c_int,
                                (*yyvsp.offset(-2_isize)).Node,
                                fits_parser_yytokentype::EQ as c_int,
                                (*yyvsp.offset(0)).Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                ((*((*lParse).Nodes).offset(yyval.Node as isize)).value).nelem = 1;
                                current_block = 17353983478346836848;
                            }
                        }
                        78 => {
                            yyval.Node = New_BinOp(
                                lParse,
                                fits_parser_yytokentype::BOOLEAN as c_int,
                                (*yyvsp.offset(-2_isize)).Node,
                                fits_parser_yytokentype::NE as c_int,
                                (*yyvsp.offset(0)).Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                ((*((*lParse).Nodes).offset(yyval.Node as isize)).value).nelem = 1;
                                current_block = 17353983478346836848;
                            }
                        }
                        79 => {
                            yyval.Node = New_BinOp(
                                lParse,
                                fits_parser_yytokentype::BOOLEAN as c_int,
                                (*yyvsp.offset(-2_isize)).Node,
                                fits_parser_yytokentype::LT as c_int,
                                (*yyvsp.offset(0)).Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                ((*((*lParse).Nodes).offset(yyval.Node as isize)).value).nelem = 1;
                                current_block = 17353983478346836848;
                            }
                        }
                        80 => {
                            yyval.Node = New_BinOp(
                                lParse,
                                fits_parser_yytokentype::BOOLEAN as c_int,
                                (*yyvsp.offset(-2_isize)).Node,
                                fits_parser_yytokentype::LTE as c_int,
                                (*yyvsp.offset(0)).Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                ((*((*lParse).Nodes).offset(yyval.Node as isize)).value).nelem = 1;
                                current_block = 17353983478346836848;
                            }
                        }
                        81 => {
                            yyval.Node = New_BinOp(
                                lParse,
                                fits_parser_yytokentype::BOOLEAN as c_int,
                                (*yyvsp.offset(-2_isize)).Node,
                                fits_parser_yytokentype::GT as c_int,
                                (*yyvsp.offset(0)).Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                ((*((*lParse).Nodes).offset(yyval.Node as isize)).value).nelem = 1;
                                current_block = 17353983478346836848;
                            }
                        }
                        82 => {
                            yyval.Node = New_BinOp(
                                lParse,
                                fits_parser_yytokentype::BOOLEAN as c_int,
                                (*yyvsp.offset(-2_isize)).Node,
                                fits_parser_yytokentype::GTE as c_int,
                                (*yyvsp.offset(0)).Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                ((*((*lParse).Nodes).offset(yyval.Node as isize)).value).nelem = 1;
                                current_block = 17353983478346836848;
                            }
                        }
                        83 => {
                            if (*((*lParse).Nodes).offset((*yyvsp.offset(-2_isize)).Node as isize))
                                .ntype
                                > (*((*lParse).Nodes).offset((*yyvsp.offset(0)).Node as isize))
                                    .ntype
                            {
                                (*yyvsp.offset(0)).Node = New_Unary(
                                    lParse,
                                    (*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-2_isize)).Node as isize))
                                    .ntype,
                                    0,
                                    (*yyvsp.offset(0)).Node,
                                );
                            } else if (*((*lParse).Nodes)
                                .offset((*yyvsp.offset(-2_isize)).Node as isize))
                            .ntype
                                < (*((*lParse).Nodes).offset((*yyvsp.offset(0)).Node as isize))
                                    .ntype
                            {
                                (*yyvsp.offset(-2_isize)).Node = New_Unary(
                                    lParse,
                                    (*((*lParse).Nodes).offset((*yyvsp.offset(0)).Node as isize))
                                        .ntype,
                                    0,
                                    (*yyvsp.offset(-2_isize)).Node,
                                );
                            }
                            yyval.Node = New_BinOp(
                                lParse,
                                fits_parser_yytokentype::BOOLEAN as c_int,
                                (*yyvsp.offset(-2_isize)).Node,
                                fits_parser_yytokentype::GT as c_int,
                                (*yyvsp.offset(0)).Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        84 => {
                            if (*((*lParse).Nodes).offset((*yyvsp.offset(-2_isize)).Node as isize))
                                .ntype
                                > (*((*lParse).Nodes).offset((*yyvsp.offset(0)).Node as isize))
                                    .ntype
                            {
                                (*yyvsp.offset(0)).Node = New_Unary(
                                    lParse,
                                    (*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-2_isize)).Node as isize))
                                    .ntype,
                                    0,
                                    (*yyvsp.offset(0)).Node,
                                );
                            } else if (*((*lParse).Nodes)
                                .offset((*yyvsp.offset(-2_isize)).Node as isize))
                            .ntype
                                < (*((*lParse).Nodes).offset((*yyvsp.offset(0)).Node as isize))
                                    .ntype
                            {
                                (*yyvsp.offset(-2_isize)).Node = New_Unary(
                                    lParse,
                                    (*((*lParse).Nodes).offset((*yyvsp.offset(0)).Node as isize))
                                        .ntype,
                                    0,
                                    (*yyvsp.offset(-2_isize)).Node,
                                );
                            }
                            yyval.Node = New_BinOp(
                                lParse,
                                fits_parser_yytokentype::BOOLEAN as c_int,
                                (*yyvsp.offset(-2_isize)).Node,
                                fits_parser_yytokentype::LT as c_int,
                                (*yyvsp.offset(0)).Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        85 => {
                            if (*((*lParse).Nodes).offset((*yyvsp.offset(-2_isize)).Node as isize))
                                .ntype
                                > (*((*lParse).Nodes).offset((*yyvsp.offset(0)).Node as isize))
                                    .ntype
                            {
                                (*yyvsp.offset(0)).Node = New_Unary(
                                    lParse,
                                    (*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-2_isize)).Node as isize))
                                    .ntype,
                                    0,
                                    (*yyvsp.offset(0)).Node,
                                );
                            } else if (*((*lParse).Nodes)
                                .offset((*yyvsp.offset(-2_isize)).Node as isize))
                            .ntype
                                < (*((*lParse).Nodes).offset((*yyvsp.offset(0)).Node as isize))
                                    .ntype
                            {
                                (*yyvsp.offset(-2_isize)).Node = New_Unary(
                                    lParse,
                                    (*((*lParse).Nodes).offset((*yyvsp.offset(0)).Node as isize))
                                        .ntype,
                                    0,
                                    (*yyvsp.offset(-2_isize)).Node,
                                );
                            }
                            yyval.Node = New_BinOp(
                                lParse,
                                fits_parser_yytokentype::BOOLEAN as c_int,
                                (*yyvsp.offset(-2_isize)).Node,
                                fits_parser_yytokentype::GTE as c_int,
                                (*yyvsp.offset(0)).Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        86 => {
                            if (*((*lParse).Nodes).offset((*yyvsp.offset(-2_isize)).Node as isize))
                                .ntype
                                > (*((*lParse).Nodes).offset((*yyvsp.offset(0)).Node as isize))
                                    .ntype
                            {
                                (*yyvsp.offset(0)).Node = New_Unary(
                                    lParse,
                                    (*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-2_isize)).Node as isize))
                                    .ntype,
                                    0,
                                    (*yyvsp.offset(0)).Node,
                                );
                            } else if (*((*lParse).Nodes)
                                .offset((*yyvsp.offset(-2_isize)).Node as isize))
                            .ntype
                                < (*((*lParse).Nodes).offset((*yyvsp.offset(0)).Node as isize))
                                    .ntype
                            {
                                (*yyvsp.offset(-2_isize)).Node = New_Unary(
                                    lParse,
                                    (*((*lParse).Nodes).offset((*yyvsp.offset(0)).Node as isize))
                                        .ntype,
                                    0,
                                    (*yyvsp.offset(-2_isize)).Node,
                                );
                            }
                            yyval.Node = New_BinOp(
                                lParse,
                                fits_parser_yytokentype::BOOLEAN as c_int,
                                (*yyvsp.offset(-2_isize)).Node,
                                fits_parser_yytokentype::LTE as c_int,
                                (*yyvsp.offset(0)).Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        87 => {
                            if (*((*lParse).Nodes).offset((*yyvsp.offset(-2_isize)).Node as isize))
                                .ntype
                                > (*((*lParse).Nodes).offset((*yyvsp.offset(0)).Node as isize))
                                    .ntype
                            {
                                (*yyvsp.offset(0)).Node = New_Unary(
                                    lParse,
                                    (*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-2_isize)).Node as isize))
                                    .ntype,
                                    0,
                                    (*yyvsp.offset(0)).Node,
                                );
                            } else if (*((*lParse).Nodes)
                                .offset((*yyvsp.offset(-2_isize)).Node as isize))
                            .ntype
                                < (*((*lParse).Nodes).offset((*yyvsp.offset(0)).Node as isize))
                                    .ntype
                            {
                                (*yyvsp.offset(-2_isize)).Node = New_Unary(
                                    lParse,
                                    (*((*lParse).Nodes).offset((*yyvsp.offset(0)).Node as isize))
                                        .ntype,
                                    0,
                                    (*yyvsp.offset(-2_isize)).Node,
                                );
                            }
                            yyval.Node = New_BinOp(
                                lParse,
                                fits_parser_yytokentype::BOOLEAN as c_int,
                                (*yyvsp.offset(-2_isize)).Node,
                                '~' as i32,
                                (*yyvsp.offset(0)).Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        88 => {
                            if (*((*lParse).Nodes).offset((*yyvsp.offset(-2_isize)).Node as isize))
                                .ntype
                                > (*((*lParse).Nodes).offset((*yyvsp.offset(0)).Node as isize))
                                    .ntype
                            {
                                (*yyvsp.offset(0)).Node = New_Unary(
                                    lParse,
                                    (*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-2_isize)).Node as isize))
                                    .ntype,
                                    0,
                                    (*yyvsp.offset(0)).Node,
                                );
                            } else if (*((*lParse).Nodes)
                                .offset((*yyvsp.offset(-2_isize)).Node as isize))
                            .ntype
                                < (*((*lParse).Nodes).offset((*yyvsp.offset(0)).Node as isize))
                                    .ntype
                            {
                                (*yyvsp.offset(-2_isize)).Node = New_Unary(
                                    lParse,
                                    (*((*lParse).Nodes).offset((*yyvsp.offset(0)).Node as isize))
                                        .ntype,
                                    0,
                                    (*yyvsp.offset(-2_isize)).Node,
                                );
                            }
                            yyval.Node = New_BinOp(
                                lParse,
                                fits_parser_yytokentype::BOOLEAN as c_int,
                                (*yyvsp.offset(-2_isize)).Node,
                                fits_parser_yytokentype::EQ as c_int,
                                (*yyvsp.offset(0)).Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        89 => {
                            if (*((*lParse).Nodes).offset((*yyvsp.offset(-2_isize)).Node as isize))
                                .ntype
                                > (*((*lParse).Nodes).offset((*yyvsp.offset(0)).Node as isize))
                                    .ntype
                            {
                                (*yyvsp.offset(0)).Node = New_Unary(
                                    lParse,
                                    (*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-2_isize)).Node as isize))
                                    .ntype,
                                    0,
                                    (*yyvsp.offset(0)).Node,
                                );
                            } else if (*((*lParse).Nodes)
                                .offset((*yyvsp.offset(-2_isize)).Node as isize))
                            .ntype
                                < (*((*lParse).Nodes).offset((*yyvsp.offset(0)).Node as isize))
                                    .ntype
                            {
                                (*yyvsp.offset(-2_isize)).Node = New_Unary(
                                    lParse,
                                    (*((*lParse).Nodes).offset((*yyvsp.offset(0)).Node as isize))
                                        .ntype,
                                    0,
                                    (*yyvsp.offset(-2_isize)).Node,
                                );
                            }
                            yyval.Node = New_BinOp(
                                lParse,
                                fits_parser_yytokentype::BOOLEAN as c_int,
                                (*yyvsp.offset(-2_isize)).Node,
                                fits_parser_yytokentype::NE as c_int,
                                (*yyvsp.offset(0)).Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        90 => {
                            yyval.Node = New_BinOp(
                                lParse,
                                fits_parser_yytokentype::BOOLEAN as c_int,
                                (*yyvsp.offset(-2_isize)).Node,
                                fits_parser_yytokentype::EQ as c_int,
                                (*yyvsp.offset(0)).Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                ((*((*lParse).Nodes).offset(yyval.Node as isize)).value).nelem = 1;
                                current_block = 17353983478346836848;
                            }
                        }
                        91 => {
                            yyval.Node = New_BinOp(
                                lParse,
                                fits_parser_yytokentype::BOOLEAN as c_int,
                                (*yyvsp.offset(-2_isize)).Node,
                                fits_parser_yytokentype::NE as c_int,
                                (*yyvsp.offset(0)).Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                ((*((*lParse).Nodes).offset(yyval.Node as isize)).value).nelem = 1;
                                current_block = 17353983478346836848;
                            }
                        }
                        92 => {
                            yyval.Node = New_BinOp(
                                lParse,
                                fits_parser_yytokentype::BOOLEAN as c_int,
                                (*yyvsp.offset(-2_isize)).Node,
                                fits_parser_yytokentype::GT as c_int,
                                (*yyvsp.offset(0)).Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                ((*((*lParse).Nodes).offset(yyval.Node as isize)).value).nelem = 1;
                                current_block = 17353983478346836848;
                            }
                        }
                        93 => {
                            yyval.Node = New_BinOp(
                                lParse,
                                fits_parser_yytokentype::BOOLEAN as c_int,
                                (*yyvsp.offset(-2_isize)).Node,
                                fits_parser_yytokentype::GTE as c_int,
                                (*yyvsp.offset(0)).Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                ((*((*lParse).Nodes).offset(yyval.Node as isize)).value).nelem = 1;
                                current_block = 17353983478346836848;
                            }
                        }
                        94 => {
                            yyval.Node = New_BinOp(
                                lParse,
                                fits_parser_yytokentype::BOOLEAN as c_int,
                                (*yyvsp.offset(-2_isize)).Node,
                                fits_parser_yytokentype::LT as c_int,
                                (*yyvsp.offset(0)).Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                ((*((*lParse).Nodes).offset(yyval.Node as isize)).value).nelem = 1;
                                current_block = 17353983478346836848;
                            }
                        }
                        95 => {
                            yyval.Node = New_BinOp(
                                lParse,
                                fits_parser_yytokentype::BOOLEAN as c_int,
                                (*yyvsp.offset(-2_isize)).Node,
                                fits_parser_yytokentype::LTE as c_int,
                                (*yyvsp.offset(0)).Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                ((*((*lParse).Nodes).offset(yyval.Node as isize)).value).nelem = 1;
                                current_block = 17353983478346836848;
                            }
                        }
                        96 => {
                            yyval.Node = New_BinOp(
                                lParse,
                                fits_parser_yytokentype::BOOLEAN as c_int,
                                (*yyvsp.offset(-2_isize)).Node,
                                fits_parser_yytokentype::AND as c_int,
                                (*yyvsp.offset(0)).Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        97 => {
                            yyval.Node = New_BinOp(
                                lParse,
                                fits_parser_yytokentype::BOOLEAN as c_int,
                                (*yyvsp.offset(-2_isize)).Node,
                                fits_parser_yytokentype::OR as c_int,
                                (*yyvsp.offset(0)).Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        98 => {
                            yyval.Node = New_BinOp(
                                lParse,
                                fits_parser_yytokentype::BOOLEAN as c_int,
                                (*yyvsp.offset(-2_isize)).Node,
                                fits_parser_yytokentype::EQ as c_int,
                                (*yyvsp.offset(0)).Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        99 => {
                            yyval.Node = New_BinOp(
                                lParse,
                                fits_parser_yytokentype::BOOLEAN as c_int,
                                (*yyvsp.offset(-2_isize)).Node,
                                fits_parser_yytokentype::NE as c_int,
                                (*yyvsp.offset(0)).Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        100 => {
                            if (*((*lParse).Nodes).offset((*yyvsp.offset(-4_isize)).Node as isize))
                                .ntype
                                > (*((*lParse).Nodes)
                                    .offset((*yyvsp.offset(-2_isize)).Node as isize))
                                .ntype
                            {
                                (*yyvsp.offset(-2_isize)).Node = New_Unary(
                                    lParse,
                                    (*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-4_isize)).Node as isize))
                                    .ntype,
                                    0,
                                    (*yyvsp.offset(-2_isize)).Node,
                                );
                            } else if (*((*lParse).Nodes)
                                .offset((*yyvsp.offset(-4_isize)).Node as isize))
                            .ntype
                                < (*((*lParse).Nodes)
                                    .offset((*yyvsp.offset(-2_isize)).Node as isize))
                                .ntype
                            {
                                (*yyvsp.offset(-4_isize)).Node = New_Unary(
                                    lParse,
                                    (*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-2_isize)).Node as isize))
                                    .ntype,
                                    0,
                                    (*yyvsp.offset(-4_isize)).Node,
                                );
                            }
                            if (*((*lParse).Nodes).offset((*yyvsp.offset(-4_isize)).Node as isize))
                                .ntype
                                > (*((*lParse).Nodes).offset((*yyvsp.offset(0)).Node as isize))
                                    .ntype
                            {
                                (*yyvsp.offset(0)).Node = New_Unary(
                                    lParse,
                                    (*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-4_isize)).Node as isize))
                                    .ntype,
                                    0,
                                    (*yyvsp.offset(0)).Node,
                                );
                            } else if (*((*lParse).Nodes)
                                .offset((*yyvsp.offset(-4_isize)).Node as isize))
                            .ntype
                                < (*((*lParse).Nodes).offset((*yyvsp.offset(0)).Node as isize))
                                    .ntype
                            {
                                (*yyvsp.offset(-4_isize)).Node = New_Unary(
                                    lParse,
                                    (*((*lParse).Nodes).offset((*yyvsp.offset(0)).Node as isize))
                                        .ntype,
                                    0,
                                    (*yyvsp.offset(-4_isize)).Node,
                                );
                            }
                            if (*((*lParse).Nodes).offset((*yyvsp.offset(-2_isize)).Node as isize))
                                .ntype
                                > (*((*lParse).Nodes).offset((*yyvsp.offset(0)).Node as isize))
                                    .ntype
                            {
                                (*yyvsp.offset(0)).Node = New_Unary(
                                    lParse,
                                    (*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-2_isize)).Node as isize))
                                    .ntype,
                                    0,
                                    (*yyvsp.offset(0)).Node,
                                );
                            } else if (*((*lParse).Nodes)
                                .offset((*yyvsp.offset(-2_isize)).Node as isize))
                            .ntype
                                < (*((*lParse).Nodes).offset((*yyvsp.offset(0)).Node as isize))
                                    .ntype
                            {
                                (*yyvsp.offset(-2_isize)).Node = New_Unary(
                                    lParse,
                                    (*((*lParse).Nodes).offset((*yyvsp.offset(0)).Node as isize))
                                        .ntype,
                                    0,
                                    (*yyvsp.offset(-2_isize)).Node,
                                );
                            }
                            (*yyvsp.offset(-2_isize)).Node = New_BinOp(
                                lParse,
                                fits_parser_yytokentype::BOOLEAN as c_int,
                                (*yyvsp.offset(-2_isize)).Node,
                                fits_parser_yytokentype::LTE as c_int,
                                (*yyvsp.offset(-4_isize)).Node,
                            );
                            (*yyvsp.offset(0)).Node = New_BinOp(
                                lParse,
                                fits_parser_yytokentype::BOOLEAN as c_int,
                                (*yyvsp.offset(-4_isize)).Node,
                                fits_parser_yytokentype::LTE as c_int,
                                (*yyvsp.offset(0)).Node,
                            );
                            yyval.Node = New_BinOp(
                                lParse,
                                fits_parser_yytokentype::BOOLEAN as c_int,
                                (*yyvsp.offset(-2_isize)).Node,
                                fits_parser_yytokentype::AND as c_int,
                                (*yyvsp.offset(0)).Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        101 => {
                            if Test_Dims(
                                lParse,
                                (*yyvsp.offset(-2_isize)).Node,
                                (*yyvsp.offset(0)).Node,
                            ) == 0
                            {
                                fits_parser_yyerror(
                                    scanner,
                                    lParse,
                                    b"Incompatible dimensions in '?:' arguments\0" as *const u8
                                        as *const c_char
                                        as *mut c_char,
                                );
                                current_block = 4830776507462815627;
                            } else {
                                yyval.Node = New_Func(
                                    lParse,
                                    0,
                                    ifthenelse_fct,
                                    3 as c_int,
                                    (*yyvsp.offset(-2_isize)).Node,
                                    (*yyvsp.offset(0)).Node,
                                    (*yyvsp.offset(-4_isize)).Node,
                                    0,
                                    0,
                                    0,
                                    0,
                                );
                                if yyval.Node < 0 {
                                    current_block = 4830776507462815627;
                                } else {
                                    if (*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-2_isize)).Node as isize))
                                    .value
                                    .nelem
                                        < (*((*lParse).Nodes)
                                            .offset((*yyvsp.offset(0)).Node as isize))
                                        .value
                                        .nelem
                                    {
                                        Copy_Dims(lParse, yyval.Node, (*yyvsp.offset(0)).Node);
                                    }
                                    if Test_Dims(lParse, (*yyvsp.offset(-4_isize)).Node, yyval.Node)
                                        == 0
                                    {
                                        fits_parser_yyerror(
                                            scanner,
                                            lParse,
                                            b"Incompatible dimensions in '?:' condition\0"
                                                as *const u8
                                                as *const c_char
                                                as *mut c_char,
                                        );
                                        current_block = 4830776507462815627;
                                    } else {
                                        if ((*((*lParse).Nodes).offset(yyval.Node as isize)).value)
                                            .nelem
                                            < (*((*lParse).Nodes)
                                                .offset((*yyvsp.offset(-4_isize)).Node as isize))
                                            .value
                                            .nelem
                                        {
                                            Copy_Dims(
                                                lParse,
                                                yyval.Node,
                                                (*yyvsp.offset(-4_isize)).Node,
                                            );
                                        }
                                        current_block = 17353983478346836848;
                                    }
                                }
                            }
                        }
                        102 => {
                            if (if ((*yyvsp.offset(-2_isize)).astr[0] as c_int)
                                < (*::core::mem::transmute::<&[u8; 8], &[c_char; 8]>(b"ISNULL(\0"))
                                    [0] as c_int
                            {
                                -(1)
                            } else if (*yyvsp.offset(-2_isize)).astr[0] as c_int
                                > (*::core::mem::transmute::<&[u8; 8], &[c_char; 8]>(b"ISNULL(\0"))
                                    [0] as c_int
                            {
                                1
                            } else {
                                strcmp(
                                    ((*yyvsp.offset(-2_isize)).astr).as_mut_ptr(),
                                    b"ISNULL(\0" as *const u8 as *const c_char,
                                )
                            }) == 0
                            {
                                yyval.Node = New_Func(
                                    lParse,
                                    0,
                                    isnull_fct,
                                    1,
                                    (*yyvsp.offset(-1_isize)).Node,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                );
                                if yyval.Node < 0 {
                                    current_block = 4830776507462815627;
                                } else {
                                    (*((*lParse).Nodes).offset(yyval.Node as isize)).ntype =
                                        fits_parser_yytokentype::BOOLEAN as c_int;
                                    current_block = 17353983478346836848;
                                }
                            } else {
                                fits_parser_yyerror(
                                    scanner,
                                    lParse,
                                    b"Boolean Function(expr) not supported\0" as *const u8
                                        as *const c_char
                                        as *mut c_char,
                                );
                                current_block = 4830776507462815627;
                            }
                        }
                        103 => {
                            if (if ((*yyvsp.offset(-2_isize)).astr[0] as c_int)
                                < (*::core::mem::transmute::<&[u8; 8], &[c_char; 8]>(b"ISNULL(\0"))
                                    [0] as c_int
                            {
                                -(1)
                            } else if (*yyvsp.offset(-2_isize)).astr[0] as c_int
                                > (*::core::mem::transmute::<&[u8; 8], &[c_char; 8]>(b"ISNULL(\0"))
                                    [0] as c_int
                            {
                                1
                            } else {
                                strcmp(
                                    ((*yyvsp.offset(-2_isize)).astr).as_mut_ptr(),
                                    b"ISNULL(\0" as *const u8 as *const c_char,
                                )
                            }) == 0
                            {
                                yyval.Node = New_Func(
                                    lParse,
                                    0,
                                    isnull_fct,
                                    1,
                                    (*yyvsp.offset(-1_isize)).Node,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                );
                                if yyval.Node < 0 {
                                    current_block = 4830776507462815627;
                                } else {
                                    (*((*lParse).Nodes).offset(yyval.Node as isize)).ntype =
                                        fits_parser_yytokentype::BOOLEAN as c_int;
                                    current_block = 17353983478346836848;
                                }
                            } else {
                                fits_parser_yyerror(
                                    scanner,
                                    lParse,
                                    b"Boolean Function(expr) not supported\0" as *const u8
                                        as *const c_char
                                        as *mut c_char,
                                );
                                current_block = 4830776507462815627;
                            }
                        }
                        104 => {
                            if (if ((*yyvsp.offset(-2_isize)).astr[0] as c_int)
                                < (*::core::mem::transmute::<&[u8; 8], &[c_char; 8]>(b"ISNULL(\0"))
                                    [0] as c_int
                            {
                                -(1)
                            } else if (*yyvsp.offset(-2_isize)).astr[0] as c_int
                                > (*::core::mem::transmute::<&[u8; 8], &[c_char; 8]>(b"ISNULL(\0"))
                                    [0] as c_int
                            {
                                1
                            } else {
                                strcmp(
                                    ((*yyvsp.offset(-2_isize)).astr).as_mut_ptr(),
                                    b"ISNULL(\0" as *const u8 as *const c_char,
                                )
                            }) == 0
                            {
                                yyval.Node = New_Func(
                                    lParse,
                                    fits_parser_yytokentype::BOOLEAN as c_int,
                                    isnull_fct,
                                    1,
                                    (*yyvsp.offset(-1_isize)).Node,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                );
                                if yyval.Node < 0 {
                                    current_block = 4830776507462815627;
                                } else {
                                    current_block = 17353983478346836848;
                                }
                            } else {
                                fits_parser_yyerror(
                                    scanner,
                                    lParse,
                                    b"Boolean Function(expr) not supported\0" as *const u8
                                        as *const c_char
                                        as *mut c_char,
                                );
                                current_block = 4830776507462815627;
                            }
                        }
                        105 => {
                            if (if ((*yyvsp.offset(-4_isize)).astr[0] as c_int)
                                < (*::core::mem::transmute::<&[u8; 9], &[c_char; 9]>(b"DEFNULL(\0"))
                                    [0] as c_int
                            {
                                -(1)
                            } else if (*yyvsp.offset(-4_isize)).astr[0] as c_int
                                > (*::core::mem::transmute::<&[u8; 9], &[c_char; 9]>(b"DEFNULL(\0"))
                                    [0] as c_int
                            {
                                1
                            } else {
                                strcmp(
                                    ((*yyvsp.offset(-4_isize)).astr).as_mut_ptr(),
                                    b"DEFNULL(\0" as *const u8 as *const c_char,
                                )
                            }) == 0
                            {
                                if (*((*lParse).Nodes)
                                    .offset((*yyvsp.offset(-3_isize)).Node as isize))
                                .value
                                .nelem
                                    >= (*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-1_isize)).Node as isize))
                                    .value
                                    .nelem
                                    && Test_Dims(
                                        lParse,
                                        (*yyvsp.offset(-3_isize)).Node,
                                        (*yyvsp.offset(-1_isize)).Node,
                                    ) != 0
                                {
                                    yyval.Node = New_Func(
                                        lParse,
                                        0,
                                        defnull_fct,
                                        2 as c_int,
                                        (*yyvsp.offset(-3_isize)).Node,
                                        (*yyvsp.offset(-1_isize)).Node,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                    );
                                    if yyval.Node < 0 {
                                        current_block = 4830776507462815627;
                                    } else {
                                        current_block = 17353983478346836848;
                                    }
                                } else {
                                    fits_parser_yyerror(
                                        scanner,
                                        lParse,
                                        b"Dimensions of DEFNULL arguments are not compatible\0"
                                            as *const u8
                                            as *const c_char
                                            as *mut c_char,
                                    );
                                    current_block = 4830776507462815627;
                                }
                            } else {
                                fits_parser_yyerror(
                                    scanner,
                                    lParse,
                                    b"Boolean Function(expr,expr) not supported\0" as *const u8
                                        as *const c_char
                                        as *mut c_char,
                                );
                                current_block = 4830776507462815627;
                            }
                        }
                        106 => {
                            if (*((*lParse).Nodes).offset((*yyvsp.offset(-(5))).Node as isize))
                                .ntype
                                != fits_parser_yytokentype::DOUBLE as c_int
                            {
                                (*yyvsp.offset(-(5))).Node = New_Unary(
                                    lParse,
                                    fits_parser_yytokentype::DOUBLE as c_int,
                                    0,
                                    (*yyvsp.offset(-(5))).Node,
                                );
                            }
                            if (*((*lParse).Nodes).offset((*yyvsp.offset(-3_isize)).Node as isize))
                                .ntype
                                != fits_parser_yytokentype::DOUBLE as c_int
                            {
                                (*yyvsp.offset(-3_isize)).Node = New_Unary(
                                    lParse,
                                    fits_parser_yytokentype::DOUBLE as c_int,
                                    0,
                                    (*yyvsp.offset(-3_isize)).Node,
                                );
                            }
                            if (*((*lParse).Nodes).offset((*yyvsp.offset(-1_isize)).Node as isize))
                                .ntype
                                != fits_parser_yytokentype::DOUBLE as c_int
                            {
                                (*yyvsp.offset(-1_isize)).Node = New_Unary(
                                    lParse,
                                    fits_parser_yytokentype::DOUBLE as c_int,
                                    0,
                                    (*yyvsp.offset(-1_isize)).Node,
                                );
                            }
                            if !(Test_Dims(
                                lParse,
                                (*yyvsp.offset(-(5))).Node,
                                (*yyvsp.offset(-3_isize)).Node,
                            ) != 0
                                && Test_Dims(
                                    lParse,
                                    (*yyvsp.offset(-3_isize)).Node,
                                    (*yyvsp.offset(-1_isize)).Node,
                                ) != 0)
                            {
                                fits_parser_yyerror(
                                    scanner,
                                    lParse,
                                    b"Dimensions of NEAR arguments are not compatible\0"
                                        as *const u8
                                        as *const c_char
                                        as *mut c_char,
                                );
                                current_block = 4830776507462815627;
                            } else if (if ((*yyvsp.offset(-(6 as c_int) as isize)).astr[0] as c_int)
                                < (*::core::mem::transmute::<&[u8; 6], &[c_char; 6]>(b"NEAR(\0"))[0]
                                    as c_int
                            {
                                -(1)
                            } else if (*yyvsp.offset(-(6 as c_int) as isize)).astr[0] as c_int
                                > (*::core::mem::transmute::<&[u8; 6], &[c_char; 6]>(b"NEAR(\0"))[0]
                                    as c_int
                            {
                                1
                            } else {
                                strcmp(
                                    ((*yyvsp.offset(-(6 as c_int) as isize)).astr).as_mut_ptr(),
                                    b"NEAR(\0" as *const u8 as *const c_char,
                                )
                            }) == 0
                            {
                                yyval.Node = New_Func(
                                    lParse,
                                    fits_parser_yytokentype::BOOLEAN as c_int,
                                    near_fct,
                                    3 as c_int,
                                    (*yyvsp.offset(-(5))).Node,
                                    (*yyvsp.offset(-3_isize)).Node,
                                    (*yyvsp.offset(-1_isize)).Node,
                                    0,
                                    0,
                                    0,
                                    0,
                                );
                                if yyval.Node < 0 {
                                    current_block = 4830776507462815627;
                                } else {
                                    if ((*((*lParse).Nodes).offset(yyval.Node as isize)).value)
                                        .nelem
                                        < (*((*lParse).Nodes)
                                            .offset((*yyvsp.offset(-(5))).Node as isize))
                                        .value
                                        .nelem
                                    {
                                        Copy_Dims(lParse, yyval.Node, (*yyvsp.offset(-(5))).Node);
                                    }
                                    if (*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-(5))).Node as isize))
                                    .value
                                    .nelem
                                        < (*((*lParse).Nodes)
                                            .offset((*yyvsp.offset(-3_isize)).Node as isize))
                                        .value
                                        .nelem
                                    {
                                        Copy_Dims(
                                            lParse,
                                            yyval.Node,
                                            (*yyvsp.offset(-3_isize)).Node,
                                        );
                                    }
                                    if (*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-3_isize)).Node as isize))
                                    .value
                                    .nelem
                                        < (*((*lParse).Nodes)
                                            .offset((*yyvsp.offset(-1_isize)).Node as isize))
                                        .value
                                        .nelem
                                    {
                                        Copy_Dims(
                                            lParse,
                                            yyval.Node,
                                            (*yyvsp.offset(-1_isize)).Node,
                                        );
                                    }
                                    current_block = 17353983478346836848;
                                }
                            } else {
                                fits_parser_yyerror(
                                    scanner,
                                    lParse,
                                    b"Boolean Function not supported\0" as *const u8
                                        as *const c_char
                                        as *mut c_char,
                                );
                                current_block = 4830776507462815627;
                            }
                        }
                        107 => {
                            if (*((*lParse).Nodes)
                                .offset((*yyvsp.offset(-(9 as c_int) as isize)).Node as isize))
                            .ntype
                                != fits_parser_yytokentype::DOUBLE as c_int
                            {
                                (*yyvsp.offset(-(9 as c_int) as isize)).Node = New_Unary(
                                    lParse,
                                    fits_parser_yytokentype::DOUBLE as c_int,
                                    0,
                                    (*yyvsp.offset(-(9 as c_int) as isize)).Node,
                                );
                            }
                            if (*((*lParse).Nodes)
                                .offset((*yyvsp.offset(-(7 as c_int) as isize)).Node as isize))
                            .ntype
                                != fits_parser_yytokentype::DOUBLE as c_int
                            {
                                (*yyvsp.offset(-(7 as c_int) as isize)).Node = New_Unary(
                                    lParse,
                                    fits_parser_yytokentype::DOUBLE as c_int,
                                    0,
                                    (*yyvsp.offset(-(7 as c_int) as isize)).Node,
                                );
                            }
                            if (*((*lParse).Nodes).offset((*yyvsp.offset(-(5))).Node as isize))
                                .ntype
                                != fits_parser_yytokentype::DOUBLE as c_int
                            {
                                (*yyvsp.offset(-(5))).Node = New_Unary(
                                    lParse,
                                    fits_parser_yytokentype::DOUBLE as c_int,
                                    0,
                                    (*yyvsp.offset(-(5))).Node,
                                );
                            }
                            if (*((*lParse).Nodes).offset((*yyvsp.offset(-3_isize)).Node as isize))
                                .ntype
                                != fits_parser_yytokentype::DOUBLE as c_int
                            {
                                (*yyvsp.offset(-3_isize)).Node = New_Unary(
                                    lParse,
                                    fits_parser_yytokentype::DOUBLE as c_int,
                                    0,
                                    (*yyvsp.offset(-3_isize)).Node,
                                );
                            }
                            if (*((*lParse).Nodes).offset((*yyvsp.offset(-1_isize)).Node as isize))
                                .ntype
                                != fits_parser_yytokentype::DOUBLE as c_int
                            {
                                (*yyvsp.offset(-1_isize)).Node = New_Unary(
                                    lParse,
                                    fits_parser_yytokentype::DOUBLE as c_int,
                                    0,
                                    (*yyvsp.offset(-1_isize)).Node,
                                );
                            }
                            if !(Test_Dims(
                                lParse,
                                (*yyvsp.offset(-(9 as c_int) as isize)).Node,
                                (*yyvsp.offset(-(7 as c_int) as isize)).Node,
                            ) != 0
                                && Test_Dims(
                                    lParse,
                                    (*yyvsp.offset(-(7 as c_int) as isize)).Node,
                                    (*yyvsp.offset(-(5))).Node,
                                ) != 0
                                && Test_Dims(
                                    lParse,
                                    (*yyvsp.offset(-(5))).Node,
                                    (*yyvsp.offset(-3_isize)).Node,
                                ) != 0
                                && Test_Dims(
                                    lParse,
                                    (*yyvsp.offset(-3_isize)).Node,
                                    (*yyvsp.offset(-1_isize)).Node,
                                ) != 0)
                            {
                                fits_parser_yyerror(
                                    scanner,
                                    lParse,
                                    b"Dimensions of CIRCLE arguments are not compatible\0"
                                        as *const u8
                                        as *const c_char
                                        as *mut c_char,
                                );
                                current_block = 4830776507462815627;
                            } else if (if ((*yyvsp.offset(-(10 as c_int) as isize)).astr[0]
                                as c_int)
                                < (*::core::mem::transmute::<&[u8; 8], &[c_char; 8]>(b"CIRCLE(\0"))
                                    [0] as c_int
                            {
                                -(1)
                            } else if (*yyvsp.offset(-(10 as c_int) as isize)).astr[0] as c_int
                                > (*::core::mem::transmute::<&[u8; 8], &[c_char; 8]>(b"CIRCLE(\0"))
                                    [0] as c_int
                            {
                                1
                            } else {
                                strcmp(
                                    ((*yyvsp.offset(-(10 as c_int) as isize)).astr).as_mut_ptr(),
                                    b"CIRCLE(\0" as *const u8 as *const c_char,
                                )
                            }) == 0
                            {
                                yyval.Node = New_Func(
                                    lParse,
                                    fits_parser_yytokentype::BOOLEAN as c_int,
                                    circle_fct,
                                    5 as c_int,
                                    (*yyvsp.offset(-(9 as c_int) as isize)).Node,
                                    (*yyvsp.offset(-(7 as c_int) as isize)).Node,
                                    (*yyvsp.offset(-(5))).Node,
                                    (*yyvsp.offset(-3_isize)).Node,
                                    (*yyvsp.offset(-1_isize)).Node,
                                    0,
                                    0,
                                );
                                if yyval.Node < 0 {
                                    current_block = 4830776507462815627;
                                } else {
                                    if ((*((*lParse).Nodes).offset(yyval.Node as isize)).value)
                                        .nelem
                                        < (*((*lParse).Nodes).offset(
                                            (*yyvsp.offset(-(9 as c_int) as isize)).Node as isize,
                                        ))
                                        .value
                                        .nelem
                                    {
                                        Copy_Dims(
                                            lParse,
                                            yyval.Node,
                                            (*yyvsp.offset(-(9 as c_int) as isize)).Node,
                                        );
                                    }
                                    if (*((*lParse).Nodes).offset(
                                        (*yyvsp.offset(-(9 as c_int) as isize)).Node as isize,
                                    ))
                                    .value
                                    .nelem
                                        < (*((*lParse).Nodes).offset(
                                            (*yyvsp.offset(-(7 as c_int) as isize)).Node as isize,
                                        ))
                                        .value
                                        .nelem
                                    {
                                        Copy_Dims(
                                            lParse,
                                            yyval.Node,
                                            (*yyvsp.offset(-(7 as c_int) as isize)).Node,
                                        );
                                    }
                                    if (*((*lParse).Nodes).offset(
                                        (*yyvsp.offset(-(7 as c_int) as isize)).Node as isize,
                                    ))
                                    .value
                                    .nelem
                                        < (*((*lParse).Nodes)
                                            .offset((*yyvsp.offset(-(5))).Node as isize))
                                        .value
                                        .nelem
                                    {
                                        Copy_Dims(lParse, yyval.Node, (*yyvsp.offset(-(5))).Node);
                                    }
                                    if (*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-(5))).Node as isize))
                                    .value
                                    .nelem
                                        < (*((*lParse).Nodes)
                                            .offset((*yyvsp.offset(-3_isize)).Node as isize))
                                        .value
                                        .nelem
                                    {
                                        Copy_Dims(
                                            lParse,
                                            yyval.Node,
                                            (*yyvsp.offset(-3_isize)).Node,
                                        );
                                    }
                                    if (*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-3_isize)).Node as isize))
                                    .value
                                    .nelem
                                        < (*((*lParse).Nodes)
                                            .offset((*yyvsp.offset(-1_isize)).Node as isize))
                                        .value
                                        .nelem
                                    {
                                        Copy_Dims(
                                            lParse,
                                            yyval.Node,
                                            (*yyvsp.offset(-1_isize)).Node,
                                        );
                                    }
                                    current_block = 17353983478346836848;
                                }
                            } else {
                                fits_parser_yyerror(
                                    scanner,
                                    lParse,
                                    b"Boolean Function not supported\0" as *const u8
                                        as *const c_char
                                        as *mut c_char,
                                );
                                current_block = 4830776507462815627;
                            }
                        }
                        108 => {
                            if (*((*lParse).Nodes)
                                .offset((*yyvsp.offset(-(13 as c_int) as isize)).Node as isize))
                            .ntype
                                != fits_parser_yytokentype::DOUBLE as c_int
                            {
                                (*yyvsp.offset(-(13 as c_int) as isize)).Node = New_Unary(
                                    lParse,
                                    fits_parser_yytokentype::DOUBLE as c_int,
                                    0,
                                    (*yyvsp.offset(-(13 as c_int) as isize)).Node,
                                );
                            }
                            if (*((*lParse).Nodes)
                                .offset((*yyvsp.offset(-(11 as c_int) as isize)).Node as isize))
                            .ntype
                                != fits_parser_yytokentype::DOUBLE as c_int
                            {
                                (*yyvsp.offset(-(11 as c_int) as isize)).Node = New_Unary(
                                    lParse,
                                    fits_parser_yytokentype::DOUBLE as c_int,
                                    0,
                                    (*yyvsp.offset(-(11 as c_int) as isize)).Node,
                                );
                            }
                            if (*((*lParse).Nodes)
                                .offset((*yyvsp.offset(-(9 as c_int) as isize)).Node as isize))
                            .ntype
                                != fits_parser_yytokentype::DOUBLE as c_int
                            {
                                (*yyvsp.offset(-(9 as c_int) as isize)).Node = New_Unary(
                                    lParse,
                                    fits_parser_yytokentype::DOUBLE as c_int,
                                    0,
                                    (*yyvsp.offset(-(9 as c_int) as isize)).Node,
                                );
                            }
                            if (*((*lParse).Nodes)
                                .offset((*yyvsp.offset(-(7 as c_int) as isize)).Node as isize))
                            .ntype
                                != fits_parser_yytokentype::DOUBLE as c_int
                            {
                                (*yyvsp.offset(-(7 as c_int) as isize)).Node = New_Unary(
                                    lParse,
                                    fits_parser_yytokentype::DOUBLE as c_int,
                                    0,
                                    (*yyvsp.offset(-(7 as c_int) as isize)).Node,
                                );
                            }
                            if (*((*lParse).Nodes).offset((*yyvsp.offset(-(5))).Node as isize))
                                .ntype
                                != fits_parser_yytokentype::DOUBLE as c_int
                            {
                                (*yyvsp.offset(-(5))).Node = New_Unary(
                                    lParse,
                                    fits_parser_yytokentype::DOUBLE as c_int,
                                    0,
                                    (*yyvsp.offset(-(5))).Node,
                                );
                            }
                            if (*((*lParse).Nodes).offset((*yyvsp.offset(-3_isize)).Node as isize))
                                .ntype
                                != fits_parser_yytokentype::DOUBLE as c_int
                            {
                                (*yyvsp.offset(-3_isize)).Node = New_Unary(
                                    lParse,
                                    fits_parser_yytokentype::DOUBLE as c_int,
                                    0,
                                    (*yyvsp.offset(-3_isize)).Node,
                                );
                            }
                            if (*((*lParse).Nodes).offset((*yyvsp.offset(-1_isize)).Node as isize))
                                .ntype
                                != fits_parser_yytokentype::DOUBLE as c_int
                            {
                                (*yyvsp.offset(-1_isize)).Node = New_Unary(
                                    lParse,
                                    fits_parser_yytokentype::DOUBLE as c_int,
                                    0,
                                    (*yyvsp.offset(-1_isize)).Node,
                                );
                            }
                            if !(Test_Dims(
                                lParse,
                                (*yyvsp.offset(-(13 as c_int) as isize)).Node,
                                (*yyvsp.offset(-(11 as c_int) as isize)).Node,
                            ) != 0
                                && Test_Dims(
                                    lParse,
                                    (*yyvsp.offset(-(11 as c_int) as isize)).Node,
                                    (*yyvsp.offset(-(9 as c_int) as isize)).Node,
                                ) != 0
                                && Test_Dims(
                                    lParse,
                                    (*yyvsp.offset(-(9 as c_int) as isize)).Node,
                                    (*yyvsp.offset(-(7 as c_int) as isize)).Node,
                                ) != 0
                                && Test_Dims(
                                    lParse,
                                    (*yyvsp.offset(-(7 as c_int) as isize)).Node,
                                    (*yyvsp.offset(-(5))).Node,
                                ) != 0
                                && Test_Dims(
                                    lParse,
                                    (*yyvsp.offset(-(5))).Node,
                                    (*yyvsp.offset(-3_isize)).Node,
                                ) != 0
                                && Test_Dims(
                                    lParse,
                                    (*yyvsp.offset(-3_isize)).Node,
                                    (*yyvsp.offset(-1_isize)).Node,
                                ) != 0)
                            {
                                fits_parser_yyerror(
                                    scanner,
                                    lParse,
                                    b"Dimensions of BOX or ELLIPSE arguments are not compatible\0"
                                        as *const u8
                                        as *const c_char
                                        as *mut c_char,
                                );
                                current_block = 4830776507462815627;
                            } else {
                                if (if ((*yyvsp.offset(-(14 as c_int) as isize)).astr[0] as c_int)
                                    < (*::core::mem::transmute::<&[u8; 5], &[c_char; 5]>(b"BOX(\0"))
                                        [0] as c_int
                                {
                                    -(1)
                                } else if (*yyvsp.offset(-(14 as c_int) as isize)).astr[0] as c_int
                                    > (*::core::mem::transmute::<&[u8; 5], &[c_char; 5]>(b"BOX(\0"))
                                        [0] as c_int
                                {
                                    1
                                } else {
                                    strcmp(
                                        ((*yyvsp.offset(-(14 as c_int) as isize)).astr)
                                            .as_mut_ptr(),
                                        b"BOX(\0" as *const u8 as *const c_char,
                                    )
                                }) == 0
                                {
                                    yyval.Node = New_Func(
                                        lParse,
                                        fits_parser_yytokentype::BOOLEAN as c_int,
                                        box_fct,
                                        7 as c_int,
                                        (*yyvsp.offset(-(13 as c_int) as isize)).Node,
                                        (*yyvsp.offset(-(11 as c_int) as isize)).Node,
                                        (*yyvsp.offset(-(9 as c_int) as isize)).Node,
                                        (*yyvsp.offset(-(7 as c_int) as isize)).Node,
                                        (*yyvsp.offset(-(5))).Node,
                                        (*yyvsp.offset(-3_isize)).Node,
                                        (*yyvsp.offset(-1_isize)).Node,
                                    );
                                    current_block = 3023179740610631044;
                                } else if (if ((*yyvsp.offset(-(14 as c_int) as isize)).astr[0]
                                    as c_int)
                                    < (*::core::mem::transmute::<&[u8; 9], &[c_char; 9]>(
                                        b"ELLIPSE(\0",
                                    ))[0] as c_int
                                {
                                    -(1)
                                } else if (*yyvsp.offset(-(14 as c_int) as isize)).astr[0] as c_int
                                    > (*::core::mem::transmute::<&[u8; 9], &[c_char; 9]>(
                                        b"ELLIPSE(\0",
                                    ))[0] as c_int
                                {
                                    1
                                } else {
                                    strcmp(
                                        ((*yyvsp.offset(-(14 as c_int) as isize)).astr)
                                            .as_mut_ptr(),
                                        b"ELLIPSE(\0" as *const u8 as *const c_char,
                                    )
                                }) == 0
                                {
                                    yyval.Node = New_Func(
                                        lParse,
                                        fits_parser_yytokentype::BOOLEAN as c_int,
                                        elps_fct,
                                        7 as c_int,
                                        (*yyvsp.offset(-(13 as c_int) as isize)).Node,
                                        (*yyvsp.offset(-(11 as c_int) as isize)).Node,
                                        (*yyvsp.offset(-(9 as c_int) as isize)).Node,
                                        (*yyvsp.offset(-(7 as c_int) as isize)).Node,
                                        (*yyvsp.offset(-(5))).Node,
                                        (*yyvsp.offset(-3_isize)).Node,
                                        (*yyvsp.offset(-1_isize)).Node,
                                    );
                                    current_block = 3023179740610631044;
                                } else {
                                    fits_parser_yyerror(
                                        scanner,
                                        lParse,
                                        b"SAO Image Function not supported\0" as *const u8
                                            as *const c_char
                                            as *mut c_char,
                                    );
                                    current_block = 4830776507462815627;
                                }
                                match current_block {
                                    4830776507462815627 => {}
                                    _ => {
                                        if yyval.Node < 0 {
                                            current_block = 4830776507462815627;
                                        } else {
                                            if ((*((*lParse).Nodes).offset(yyval.Node as isize))
                                                .value)
                                                .nelem
                                                < (*((*lParse).Nodes).offset(
                                                    (*yyvsp.offset(-(13 as c_int) as isize)).Node
                                                        as isize,
                                                ))
                                                .value
                                                .nelem
                                            {
                                                Copy_Dims(
                                                    lParse,
                                                    yyval.Node,
                                                    (*yyvsp.offset(-(13 as c_int) as isize)).Node,
                                                );
                                            }
                                            if (*((*lParse).Nodes).offset(
                                                (*yyvsp.offset(-(13 as c_int) as isize)).Node
                                                    as isize,
                                            ))
                                            .value
                                            .nelem
                                                < (*((*lParse).Nodes).offset(
                                                    (*yyvsp.offset(-(11 as c_int) as isize)).Node
                                                        as isize,
                                                ))
                                                .value
                                                .nelem
                                            {
                                                Copy_Dims(
                                                    lParse,
                                                    yyval.Node,
                                                    (*yyvsp.offset(-(11 as c_int) as isize)).Node,
                                                );
                                            }
                                            if (*((*lParse).Nodes).offset(
                                                (*yyvsp.offset(-(11 as c_int) as isize)).Node
                                                    as isize,
                                            ))
                                            .value
                                            .nelem
                                                < (*((*lParse).Nodes).offset(
                                                    (*yyvsp.offset(-(9 as c_int) as isize)).Node
                                                        as isize,
                                                ))
                                                .value
                                                .nelem
                                            {
                                                Copy_Dims(
                                                    lParse,
                                                    yyval.Node,
                                                    (*yyvsp.offset(-(9 as c_int) as isize)).Node,
                                                );
                                            }
                                            if (*((*lParse).Nodes).offset(
                                                (*yyvsp.offset(-(9 as c_int) as isize)).Node
                                                    as isize,
                                            ))
                                            .value
                                            .nelem
                                                < (*((*lParse).Nodes).offset(
                                                    (*yyvsp.offset(-(7 as c_int) as isize)).Node
                                                        as isize,
                                                ))
                                                .value
                                                .nelem
                                            {
                                                Copy_Dims(
                                                    lParse,
                                                    yyval.Node,
                                                    (*yyvsp.offset(-(7 as c_int) as isize)).Node,
                                                );
                                            }
                                            if (*((*lParse).Nodes).offset(
                                                (*yyvsp.offset(-(7 as c_int) as isize)).Node
                                                    as isize,
                                            ))
                                            .value
                                            .nelem
                                                < (*((*lParse).Nodes)
                                                    .offset((*yyvsp.offset(-(5))).Node as isize))
                                                .value
                                                .nelem
                                            {
                                                Copy_Dims(
                                                    lParse,
                                                    yyval.Node,
                                                    (*yyvsp.offset(-(5))).Node,
                                                );
                                            }
                                            if (*((*lParse).Nodes)
                                                .offset((*yyvsp.offset(-(5))).Node as isize))
                                            .value
                                            .nelem
                                                < (*((*lParse).Nodes).offset(
                                                    (*yyvsp.offset(-3_isize)).Node as isize,
                                                ))
                                                .value
                                                .nelem
                                            {
                                                Copy_Dims(
                                                    lParse,
                                                    yyval.Node,
                                                    (*yyvsp.offset(-3_isize)).Node,
                                                );
                                            }
                                            if (*((*lParse).Nodes)
                                                .offset((*yyvsp.offset(-3_isize)).Node as isize))
                                            .value
                                            .nelem
                                                < (*((*lParse).Nodes).offset(
                                                    (*yyvsp.offset(-1_isize)).Node as isize,
                                                ))
                                                .value
                                                .nelem
                                            {
                                                Copy_Dims(
                                                    lParse,
                                                    yyval.Node,
                                                    (*yyvsp.offset(-1_isize)).Node,
                                                );
                                            }
                                            current_block = 17353983478346836848;
                                        }
                                    }
                                }
                            }
                        }
                        109 => {
                            yyval.Node = New_GTI(
                                lParse,
                                gtifilt_fct,
                                b"\0" as *const u8 as *const c_char as *mut c_char,
                                -(99 as c_int),
                                -(99 as c_int),
                                b"*START*\0" as *const u8 as *const c_char as *mut c_char,
                                b"*STOP*\0" as *const u8 as *const c_char as *mut c_char,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        110 => {
                            yyval.Node = New_GTI(
                                lParse,
                                gtifilt_fct,
                                ((*yyvsp.offset(-1_isize)).astr).as_mut_ptr(),
                                -(99 as c_int),
                                -(99 as c_int),
                                b"*START*\0" as *const u8 as *const c_char as *mut c_char,
                                b"*STOP*\0" as *const u8 as *const c_char as *mut c_char,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        111 => {
                            yyval.Node = New_GTI(
                                lParse,
                                gtifilt_fct,
                                ((*yyvsp.offset(-3_isize)).astr).as_mut_ptr(),
                                (*yyvsp.offset(-1_isize)).Node,
                                -(99 as c_int),
                                b"*START*\0" as *const u8 as *const c_char as *mut c_char,
                                b"*STOP*\0" as *const u8 as *const c_char as *mut c_char,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        112 => {
                            yyval.Node = New_GTI(
                                lParse,
                                gtifilt_fct,
                                ((*yyvsp.offset(-(7 as c_int) as isize)).astr).as_mut_ptr(),
                                (*yyvsp.offset(-(5))).Node,
                                -(99 as c_int),
                                ((*yyvsp.offset(-3_isize)).astr).as_mut_ptr(),
                                ((*yyvsp.offset(-1_isize)).astr).as_mut_ptr(),
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        113 => {
                            yyval.Node = New_GTI(
                                lParse,
                                gtifind_fct,
                                b"\0" as *const u8 as *const c_char as *mut c_char,
                                -(99 as c_int),
                                -(99 as c_int),
                                b"*START*\0" as *const u8 as *const c_char as *mut c_char,
                                b"*STOP*\0" as *const u8 as *const c_char as *mut c_char,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        114 => {
                            yyval.Node = New_GTI(
                                lParse,
                                gtifind_fct,
                                ((*yyvsp.offset(-1_isize)).astr).as_mut_ptr(),
                                -(99 as c_int),
                                -(99 as c_int),
                                b"*START*\0" as *const u8 as *const c_char as *mut c_char,
                                b"*STOP*\0" as *const u8 as *const c_char as *mut c_char,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        115 => {
                            yyval.Node = New_GTI(
                                lParse,
                                gtifind_fct,
                                ((*yyvsp.offset(-3_isize)).astr).as_mut_ptr(),
                                (*yyvsp.offset(-1_isize)).Node,
                                -(99 as c_int),
                                b"*START*\0" as *const u8 as *const c_char as *mut c_char,
                                b"*STOP*\0" as *const u8 as *const c_char as *mut c_char,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        116 => {
                            yyval.Node = New_GTI(
                                lParse,
                                gtifind_fct,
                                ((*yyvsp.offset(-(7 as c_int) as isize)).astr).as_mut_ptr(),
                                (*yyvsp.offset(-(5))).Node,
                                -(99 as c_int),
                                ((*yyvsp.offset(-3_isize)).astr).as_mut_ptr(),
                                ((*yyvsp.offset(-1_isize)).astr).as_mut_ptr(),
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        117 => {
                            yyval.Node = New_REG(
                                lParse,
                                ((*yyvsp.offset(-1_isize)).astr).as_mut_ptr(),
                                -(99 as c_int),
                                -(99 as c_int),
                                b"\0" as *const u8 as *const c_char as *mut c_char,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        118 => {
                            yyval.Node = New_REG(
                                lParse,
                                ((*yyvsp.offset(-(5))).astr).as_mut_ptr(),
                                (*yyvsp.offset(-3_isize)).Node,
                                (*yyvsp.offset(-1_isize)).Node,
                                b"\0" as *const u8 as *const c_char as *mut c_char,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        119 => {
                            yyval.Node = New_REG(
                                lParse,
                                ((*yyvsp.offset(-(7 as c_int) as isize)).astr).as_mut_ptr(),
                                (*yyvsp.offset(-(5))).Node,
                                (*yyvsp.offset(-3_isize)).Node,
                                ((*yyvsp.offset(-1_isize)).astr).as_mut_ptr(),
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        120 => {
                            yyval.Node = New_Deref(
                                lParse,
                                (*yyvsp.offset(-3_isize)).Node,
                                1,
                                (*yyvsp.offset(-1_isize)).Node,
                                0,
                                0,
                                0,
                                0,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        121 => {
                            yyval.Node = New_Deref(
                                lParse,
                                (*yyvsp.offset(-(5))).Node,
                                2 as c_int,
                                (*yyvsp.offset(-3_isize)).Node,
                                (*yyvsp.offset(-1_isize)).Node,
                                0,
                                0,
                                0,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        122 => {
                            yyval.Node = New_Deref(
                                lParse,
                                (*yyvsp.offset(-(7 as c_int) as isize)).Node,
                                3 as c_int,
                                (*yyvsp.offset(-(5))).Node,
                                (*yyvsp.offset(-3_isize)).Node,
                                (*yyvsp.offset(-1_isize)).Node,
                                0,
                                0,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        123 => {
                            yyval.Node = New_Deref(
                                lParse,
                                (*yyvsp.offset(-(9 as c_int) as isize)).Node,
                                4 as c_int,
                                (*yyvsp.offset(-(7 as c_int) as isize)).Node,
                                (*yyvsp.offset(-(5))).Node,
                                (*yyvsp.offset(-3_isize)).Node,
                                (*yyvsp.offset(-1_isize)).Node,
                                0,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        124 => {
                            yyval.Node = New_Deref(
                                lParse,
                                (*yyvsp.offset(-(11 as c_int) as isize)).Node,
                                5 as c_int,
                                (*yyvsp.offset(-(9 as c_int) as isize)).Node,
                                (*yyvsp.offset(-(7 as c_int) as isize)).Node,
                                (*yyvsp.offset(-(5))).Node,
                                (*yyvsp.offset(-3_isize)).Node,
                                (*yyvsp.offset(-1_isize)).Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        125 => {
                            yyval.Node = New_Unary(
                                lParse,
                                fits_parser_yytokentype::BOOLEAN as c_int,
                                fits_parser_yytokentype::NOT as c_int,
                                (*yyvsp.offset(0)).Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        126 => {
                            yyval.Node = (*yyvsp.offset(-1_isize)).Node;
                            current_block = 17353983478346836848;
                        }
                        127 => {
                            yyval.Node = New_Const(
                                lParse,
                                fits_parser_yytokentype::STRING as c_int,
                                ((*yyvsp.offset(0)).astr).as_mut_ptr() as *mut c_void,
                                (strlen(((*yyvsp.offset(0)).astr).as_mut_ptr()))
                                    .wrapping_add((1 as c_int as c_ulong).try_into().unwrap())
                                    as c_long,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                ((*((*lParse).Nodes).offset(yyval.Node as isize)).value).nelem =
                                    strlen(((*yyvsp.offset(0)).astr).as_mut_ptr()) as c_long;
                                current_block = 17353983478346836848;
                            }
                        }
                        128 => {
                            yyval.Node = New_Column(lParse, (*yyvsp.offset(0)).lng as c_int);
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        129 => {
                            if (*((*lParse).Nodes).offset((*yyvsp.offset(-1_isize)).Node as isize))
                                .ntype
                                != fits_parser_yytokentype::LONG as c_int
                                || (*((*lParse).Nodes)
                                    .offset((*yyvsp.offset(-1_isize)).Node as isize))
                                .operation
                                    != -(1000 as c_int)
                            {
                                fits_parser_yyerror(
                                    scanner,
                                    lParse,
                                    b"Offset argument must be a constant integer\0" as *const u8
                                        as *const c_char
                                        as *mut c_char,
                                );
                                current_block = 4830776507462815627;
                            } else {
                                yyval.Node = New_Offset(
                                    lParse,
                                    (*yyvsp.offset(-3_isize)).lng as c_int,
                                    (*yyvsp.offset(-1_isize)).Node,
                                );
                                if yyval.Node < 0 {
                                    current_block = 4830776507462815627;
                                } else {
                                    current_block = 17353983478346836848;
                                }
                            }
                        }
                        130 => {
                            yyval.Node = New_Func(
                                lParse,
                                fits_parser_yytokentype::STRING as c_int,
                                null_fct,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                            );
                            current_block = 17353983478346836848;
                        }
                        131 => {
                            yyval.Node = (*yyvsp.offset(-1_isize)).Node;
                            current_block = 17353983478346836848;
                        }
                        132 => {
                            if (*((*lParse).Nodes).offset((*yyvsp.offset(-2_isize)).Node as isize))
                                .value
                                .nelem
                                + (*((*lParse).Nodes).offset((*yyvsp.offset(0)).Node as isize))
                                    .value
                                    .nelem
                                >= 256 as c_long
                            {
                                fits_parser_yyerror(
                                    scanner,
                                    lParse,
                                    b"Combined string size exceeds 255 characters\0" as *const u8
                                        as *const c_char
                                        as *mut c_char,
                                );
                                current_block = 4830776507462815627;
                            } else {
                                yyval.Node = New_BinOp(
                                    lParse,
                                    fits_parser_yytokentype::STRING as c_int,
                                    (*yyvsp.offset(-2_isize)).Node,
                                    '+' as i32,
                                    (*yyvsp.offset(0)).Node,
                                );
                                if yyval.Node < 0 {
                                    current_block = 4830776507462815627;
                                } else {
                                    ((*((*lParse).Nodes).offset(yyval.Node as isize)).value)
                                        .nelem = (*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-2_isize)).Node as isize))
                                    .value
                                    .nelem
                                        + (*((*lParse).Nodes)
                                            .offset((*yyvsp.offset(0)).Node as isize))
                                        .value
                                        .nelem;
                                    current_block = 17353983478346836848;
                                }
                            }
                        }
                        133 => {
                            let mut outSize: c_int = 0;
                            if (*((*lParse).Nodes).offset((*yyvsp.offset(-4_isize)).Node as isize))
                                .value
                                .nelem
                                != 1
                            {
                                fits_parser_yyerror(
                                    scanner,
                                    lParse,
                                    b"Cannot have a vector string column\0" as *const u8
                                        as *const c_char
                                        as *mut c_char,
                                );
                                current_block = 4830776507462815627;
                            } else {
                                outSize = ((*((*lParse).Nodes)
                                    .offset((*yyvsp.offset(-2_isize)).Node as isize))
                                .value)
                                    .nelem as c_int;
                                if (*((*lParse).Nodes).offset((*yyvsp.offset(0)).Node as isize))
                                    .value
                                    .nelem
                                    > outSize as c_long
                                {
                                    outSize = (*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(0)).Node as isize))
                                    .value
                                    .nelem as c_int;
                                }
                                yyval.Node = New_FuncSize(
                                    lParse,
                                    0,
                                    ifthenelse_fct,
                                    3 as c_int,
                                    (*yyvsp.offset(-2_isize)).Node,
                                    (*yyvsp.offset(0)).Node,
                                    (*yyvsp.offset(-4_isize)).Node,
                                    0,
                                    0,
                                    0,
                                    0,
                                    outSize,
                                );
                                if yyval.Node < 0 {
                                    current_block = 4830776507462815627;
                                } else {
                                    if (*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-2_isize)).Node as isize))
                                    .value
                                    .nelem
                                        < (*((*lParse).Nodes)
                                            .offset((*yyvsp.offset(0)).Node as isize))
                                        .value
                                        .nelem
                                    {
                                        Copy_Dims(lParse, yyval.Node, (*yyvsp.offset(0)).Node);
                                    }
                                    current_block = 17353983478346836848;
                                }
                            }
                        }
                        134 => {
                            if (if ((*yyvsp.offset(-4_isize)).astr[0] as c_int)
                                < (*::core::mem::transmute::<&[u8; 9], &[c_char; 9]>(b"DEFNULL(\0"))
                                    [0] as c_int
                            {
                                -(1)
                            } else if (*yyvsp.offset(-4_isize)).astr[0] as c_int
                                > (*::core::mem::transmute::<&[u8; 9], &[c_char; 9]>(b"DEFNULL(\0"))
                                    [0] as c_int
                            {
                                1
                            } else {
                                strcmp(
                                    ((*yyvsp.offset(-4_isize)).astr).as_mut_ptr(),
                                    b"DEFNULL(\0" as *const u8 as *const c_char,
                                )
                            }) == 0
                            {
                                let mut outSize_0: c_int = 0;
                                outSize_0 = (*((*lParse).Nodes)
                                    .offset((*yyvsp.offset(-3_isize)).Node as isize))
                                .value
                                .nelem as c_int;
                                if (*((*lParse).Nodes)
                                    .offset((*yyvsp.offset(-1_isize)).Node as isize))
                                .value
                                .nelem
                                    > outSize_0 as c_long
                                {
                                    outSize_0 = (*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-1_isize)).Node as isize))
                                    .value
                                    .nelem as c_int;
                                }
                                yyval.Node = New_FuncSize(
                                    lParse,
                                    0,
                                    defnull_fct,
                                    2 as c_int,
                                    (*yyvsp.offset(-3_isize)).Node,
                                    (*yyvsp.offset(-1_isize)).Node,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                    outSize_0,
                                );
                                if yyval.Node < 0 {
                                    current_block = 4830776507462815627;
                                } else {
                                    if (*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-1_isize)).Node as isize))
                                    .value
                                    .nelem
                                        > (*((*lParse).Nodes)
                                            .offset((*yyvsp.offset(-3_isize)).Node as isize))
                                        .value
                                        .nelem
                                    {
                                        ((*((*lParse).Nodes).offset(yyval.Node as isize)).value)
                                            .nelem = (*((*lParse).Nodes)
                                            .offset((*yyvsp.offset(-1_isize)).Node as isize))
                                        .value
                                        .nelem;
                                    }
                                    current_block = 17353983478346836848;
                                }
                            } else {
                                fits_parser_yyerror(
                                    scanner,
                                    lParse,
                                    b"Function(string,string) not supported\0" as *const u8
                                        as *const c_char
                                        as *mut c_char,
                                );
                                current_block = 4830776507462815627;
                            }
                        }
                        135 => {
                            if (if ((*yyvsp.offset(-(6 as c_int) as isize)).astr[0] as c_int)
                                < (*::core::mem::transmute::<&[u8; 8], &[c_char; 8]>(b"STRMID(\0"))
                                    [0] as c_int
                            {
                                -(1)
                            } else if (*yyvsp.offset(-(6 as c_int) as isize)).astr[0] as c_int
                                > (*::core::mem::transmute::<&[u8; 8], &[c_char; 8]>(b"STRMID(\0"))
                                    [0] as c_int
                            {
                                1
                            } else {
                                strcmp(
                                    ((*yyvsp.offset(-(6 as c_int) as isize)).astr).as_mut_ptr(),
                                    b"STRMID(\0" as *const u8 as *const c_char,
                                )
                            }) == 0
                            {
                                let mut len: c_int = 0;
                                if (*((*lParse).Nodes)
                                    .offset((*yyvsp.offset(-3_isize)).Node as isize))
                                .ntype
                                    != fits_parser_yytokentype::LONG as c_int
                                    || (*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-3_isize)).Node as isize))
                                    .value
                                    .nelem
                                        != 1
                                    || (*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-1_isize)).Node as isize))
                                    .ntype
                                        != fits_parser_yytokentype::LONG as c_int
                                    || (*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-1_isize)).Node as isize))
                                    .value
                                    .nelem
                                        != 1
                                {
                                    fits_parser_yyerror(
                                    scanner,
                                    lParse,
                                    b"When using STRMID(S,P,N), P and N must be integers (and not vector columns)\0"
                                        as *const u8 as *const c_char
                                        as *mut c_char,
                                );
                                    current_block = 4830776507462815627;
                                } else {
                                    if (*((*lParse).Nodes)
                                        .offset((*yyvsp.offset(-1_isize)).Node as isize))
                                    .operation
                                        == -(1000 as c_int)
                                    {
                                        len = (*((*lParse).Nodes)
                                            .offset((*yyvsp.offset(-1_isize)).Node as isize))
                                        .value
                                        .data
                                        .lng as c_int;
                                    } else {
                                        len = (*((*lParse).Nodes)
                                            .offset((*yyvsp.offset(-(5))).Node as isize))
                                        .value
                                        .nelem
                                            as c_int;
                                    }
                                    if len <= 0 || len >= 256 as c_int {
                                        fits_parser_yyerror(
                                            scanner,
                                            lParse,
                                            b"STRMID(S,P,N), N must be 1-255\0" as *const u8
                                                as *const c_char
                                                as *mut c_char,
                                        );
                                        current_block = 4830776507462815627;
                                    } else {
                                        yyval.Node = New_FuncSize(
                                            lParse,
                                            0,
                                            strmid_fct,
                                            3 as c_int,
                                            (*yyvsp.offset(-(5))).Node,
                                            (*yyvsp.offset(-3_isize)).Node,
                                            (*yyvsp.offset(-1_isize)).Node,
                                            0,
                                            0,
                                            0,
                                            0,
                                            len,
                                        );
                                        if yyval.Node < 0 {
                                            current_block = 4830776507462815627;
                                        } else {
                                            current_block = 17353983478346836848;
                                        }
                                    }
                                }
                            } else {
                                fits_parser_yyerror(
                                    scanner,
                                    lParse,
                                    b"Function(string,expr,expr) not supported\0" as *const u8
                                        as *const c_char
                                        as *mut c_char,
                                );
                                current_block = 4830776507462815627;
                            }
                        }
                        _ => {
                            current_block = 17353983478346836848;
                        }
                    }
                    match current_block {
                        4830776507462815627 => {
                            fits_parser_yynerrs += 1;
                            fits_parser_yynerrs;
                            yyvsp = yyvsp.offset(-(yylen as isize));
                            yyssp = yyssp.offset(-(yylen as isize));
                            yylen = 0;
                            yystate = *yyssp as yy_state_fast_t;
                            current_block = 1774893048582444437;
                        }
                        _ => {
                            yyvsp = yyvsp.offset(-(yylen as isize));
                            yyssp = yyssp.offset(-(yylen as isize));
                            yylen = 0;
                            yyvsp = yyvsp.offset(1);
                            *yyvsp = yyval;
                            let yylhs: c_int = yyr1[yyn as usize] as c_int - 57 as c_int;
                            let yyi: c_int = yypgoto[yylhs as usize] as c_int + *yyssp as c_int;
                            yystate = if 0 <= yyi
                                && yyi <= 1776 as c_int
                                && yycheck[yyi as usize] as c_int == *yyssp as c_int
                            {
                                yytable[yyi as usize] as c_int
                            } else {
                                yydefgoto[yylhs as usize] as c_int
                            };
                            current_block = 7872030484262409139;
                        }
                    }
                }
                if current_block == 1774893048582444437 {
                    yyerrstatus = 3 as c_int;
                    loop {
                        yyn = yypact[yystate as usize] as c_int;
                        if yyn != -(41 as c_int) {
                            yyn += YYSYMBOL_YYerror as c_int;
                            if 0 <= yyn
                                && yyn <= 1776 as c_int
                                && yycheck[yyn as usize] as c_int == YYSYMBOL_YYerror as c_int
                            {
                                yyn = yytable[yyn as usize] as c_int;
                                if (0 as c_int) < yyn {
                                    break;
                                }
                            }
                        }
                        if yyssp == yyss {
                            current_block = 3964311021479492664;
                            break 's_54;
                        }
                        yydestruct(
                            b"Error: popping\0" as *const u8 as *const c_char,
                            yystos[yystate as usize] as yysymbol_kind_t,
                            yyvsp,
                            scanner,
                            lParse,
                        );
                        yyvsp = yyvsp.offset(-(1 as c_int as isize));
                        yyssp = yyssp.offset(-(1 as c_int as isize));
                        yystate = *yyssp as yy_state_fast_t;
                    }
                    yyvsp = yyvsp.offset(1);
                    *yyvsp = yylval;
                    yystate = yyn;
                }
                yyssp = yyssp.offset(1);
                yyssp;
            }
        }
        match current_block {
            11794367917084412820 => {
                fits_parser_yyerror(
                    scanner,
                    lParse,
                    b"memory exhausted\0" as *const u8 as *const c_char as *mut c_char,
                );
                yyresult = 2 as c_int;
            }
            3964311021479492664 => {
                yyresult = 1;
            }
            _ => {}
        }
        if yychar != fits_parser_yytokentype::FITS_PARSER_YYEMPTY as c_int {
            yytoken = (if 0 <= yychar && yychar <= 292 as c_int {
                yytranslate[yychar as usize] as yysymbol_kind_t as c_int
            } else {
                YYSYMBOL_YYUNDEF as c_int
            }) as yysymbol_kind_t;
            yydestruct(
                b"Cleanup: discarding lookahead\0" as *const u8 as *const c_char,
                yytoken,
                &mut yylval,
                scanner,
                lParse,
            );
        }
        yyvsp = yyvsp.offset(-(yylen as isize));
        yyssp = yyssp.offset(-(yylen as isize));
        while yyssp != yyss {
            yydestruct(
                b"Cleanup: popping\0" as *const u8 as *const c_char,
                yystos[*yyssp as usize] as yysymbol_kind_t,
                yyvsp,
                scanner,
                lParse,
            );
            yyvsp = yyvsp.offset(-(1 as c_int as isize));
            yyssp = yyssp.offset(-(1 as c_int as isize));
        }
        if yyss != yyssa.as_mut_ptr() {
            free(yyss as *mut c_void);
        }
        yyresult
    }
}
unsafe fn New_Deref(
    mut lParse: *mut ParseData,
    mut Var: c_int,
    mut nDim: c_int,
    mut Dim1: c_int,
    mut Dim2: c_int,
    mut Dim3: c_int,
    mut Dim4: c_int,
    mut Dim5: c_int,
) -> c_int {
    unsafe {
        let mut n: c_int = 0;
        let mut idx: c_int = 0;
        let mut constant: c_int = 0;
        let mut elem: c_long = 0;
        let mut this: *mut Node = std::ptr::null_mut::<Node>();
        let mut theVar: *mut Node = std::ptr::null_mut::<Node>();
        let mut theDim: [*mut Node; 5] = [std::ptr::null_mut::<Node>(); 5];
        if Var < 0 || Dim1 < 0 || Dim2 < 0 || Dim3 < 0 || Dim4 < 0 || Dim5 < 0 {
            return -(1);
        }
        theVar = ((*lParse).Nodes).offset(Var as isize);
        if (*theVar).operation == -(1000 as c_int) || (*theVar).value.nelem == 1 {
            fits_parser_yyerror(
                0 as yyscan_t,
                lParse,
                b"Cannot index a scalar value\0" as *const u8 as *const c_char as *mut c_char,
            );
            return -(1);
        }
        n = Alloc_Node(lParse);
        if n >= 0 {
            this = ((*lParse).Nodes).offset(n as isize);
            (*this).nSubNodes = nDim + 1;
            (*this).SubNodes[0] = Var;
            theVar = ((*lParse).Nodes).offset((*this).SubNodes[0] as isize);
            (*this).SubNodes[1] = Dim1;
            theDim[0] = ((*lParse).Nodes).offset((*this).SubNodes[1] as isize);
            (*this).SubNodes[2] = Dim2;
            theDim[1] = ((*lParse).Nodes).offset((*this).SubNodes[2] as isize);
            (*this).SubNodes[3] = Dim3;
            theDim[2] = ((*lParse).Nodes).offset((*this).SubNodes[3] as isize);
            (*this).SubNodes[4] = Dim4;
            theDim[3] = ((*lParse).Nodes).offset((*this).SubNodes[4] as isize);
            (*this).SubNodes[5] = Dim5;
            theDim[4] = ((*lParse).Nodes).offset((*this).SubNodes[5] as isize);
            constant = ((*theVar).operation == -(1000 as c_int)) as c_int;
            idx = 0;
            while idx < nDim {
                constant = (constant != 0 && (*theDim[idx as usize]).operation == -(1000 as c_int))
                    as c_int;
                idx += 1;
                idx;
            }
            idx = 0;
            while idx < nDim {
                if (*theDim[idx as usize]).value.nelem > 1 {
                    Free_Last_Node(lParse);
                    fits_parser_yyerror(
                        0 as yyscan_t,
                        lParse,
                        b"Cannot use an array as an index value\0" as *const u8 as *const c_char
                            as *mut c_char,
                    );
                    return -(1);
                } else if (*theDim[idx as usize]).ntype != fits_parser_yytokentype::LONG as c_int {
                    Free_Last_Node(lParse);
                    fits_parser_yyerror(
                        0 as yyscan_t,
                        lParse,
                        b"Index value must be an integer type\0" as *const u8 as *const c_char
                            as *mut c_char,
                    );
                    return -(1);
                }
                idx += 1;
                idx;
            }
            (*this).operation = '[' as i32;
            (*this).DoOp = Some(Do_Deref);
            (*this).ntype = (*theVar).ntype;
            if (*theVar).value.naxis == nDim {
                (*this).value.nelem = 1;
                (*this).value.naxis = 1;
                (*this).value.naxes[0] = 1;
            } else if nDim == 1 {
                elem = 1;
                (*this).value.naxis = (*theVar).value.naxis - 1;
                idx = 0;
                while idx < (*this).value.naxis {
                    (*this).value.naxes[idx as usize] = (*theVar).value.naxes[idx as usize];
                    elem *= (*this).value.naxes[idx as usize];
                    idx += 1;
                    idx;
                }
                (*this).value.nelem = elem;
            } else {
                Free_Last_Node(lParse);
                fits_parser_yyerror(
                    0 as yyscan_t,
                    lParse,
                    b"Must specify just one or all indices for vector\0" as *const u8
                        as *const c_char as *mut c_char,
                );
                return -(1);
            }
            if constant != 0 {
                ((*this).DoOp).expect("non-null function pointer")(lParse, this);
            }
        }
        n
    }
}
unsafe fn New_GTI(
    mut lParse: *mut ParseData,
    mut Op: funcOp,
    mut fname: *mut c_char,
    mut Node1: c_int,
    mut Node2: c_int,
    mut start: *mut c_char,
    mut stop: *mut c_char,
) -> c_int {
    unsafe {
        let mut fptr: *mut fitsfile = std::ptr::null_mut();
        let mut this: *mut Node = std::ptr::null_mut::<Node>();
        let mut that0: *mut Node = std::ptr::null_mut::<Node>();
        let mut that1: *mut Node = std::ptr::null_mut::<Node>();
        let mut that2: *mut Node = std::ptr::null_mut::<Node>();
        let mut type_0: c_int = 0;
        let mut i: c_int = 0;
        let mut n: c_int = 0;
        let mut startCol: c_int = 0;
        let mut stopCol: c_int = 0;
        let mut Node0: c_int = 0;
        let mut hdutype: c_int = 0;
        let mut hdunum: c_int = 0;
        let mut evthdu: c_int = 0;
        let mut samefile: c_int = 0;
        let mut extvers: c_int = 0;
        let mut movetotype: c_int = 0;
        let mut tstat: c_int = 0;
        let mut extname: [c_char; 100] = [0; 100];
        let mut nrows: c_long = 0;
        let mut timeZeroI: [c_double; 2] = [0.; 2];
        let mut timeZeroF: [c_double; 2] = [0.; 2];
        let mut dt: c_double = 0.;
        let mut timeSpan: c_double = 0.;
        let mut xcol: [c_char; 20] = [0; 20];
        let mut xexpr: [c_char; 20] = [0; 20];
        let mut colVal: FITS_PARSER_YYSTYPE = FITS_PARSER_YYSTYPE { Node: 0 };
        if (Op as c_uint == gtifilt_fct as c_int as c_uint
            || Op as c_uint == gtifind_fct as c_int as c_uint)
            && Node1 == -(99 as c_int)
        {
            type_0 = fits_parser_yyGetVariable(
                lParse,
                b"TIME\0" as *const u8 as *const c_char as *mut c_char,
                &mut colVal,
            );
            if type_0 == fits_parser_yytokentype::COLUMN as c_int {
                Node1 = New_Column(lParse, colVal.lng as c_int);
            } else {
                fits_parser_yyerror(
                    0 as yyscan_t,
                    lParse,
                    b"Could not build TIME column for GTIFILTER/GTIFIND\0" as *const u8
                        as *const c_char as *mut c_char,
                );
                return -(1);
            }
        }
        if Op as c_uint == gtiover_fct as c_int as c_uint {
            if Node1 == -(99 as c_int) || Node2 == -(99 as c_int) {
                fits_parser_yyerror(
                    0 as yyscan_t,
                    lParse,
                    b"startExpr and stopExpr values must be defined for GTIOVERLAP\0" as *const u8
                        as *const c_char as *mut c_char,
                );
                return -(1);
            }
            Node2 = New_Unary(lParse, fits_parser_yytokentype::DOUBLE as c_int, 0, Node2);
            if Node2 < 0 {
                return -(1);
            }
        }
        Node1 = New_Unary(lParse, fits_parser_yytokentype::DOUBLE as c_int, 0, Node1);
        Node0 = Alloc_Node(lParse);
        if Node1 < 0 || Node0 < 0 {
            return -(1);
        }
        fptr = (*lParse).def_fptr;
        ffghdn(fptr, &mut evthdu);
        tstat = 0;
        if ffgkyd(
            fptr,
            b"TIMEZERO\0" as *const u8 as *const c_char,
            timeZeroI.as_mut_ptr(),
            std::ptr::null_mut::<c_char>(),
            &mut tstat,
        ) != 0
        {
            tstat = 0;
            if ffgkyd(
                fptr,
                b"TIMEZERI\0" as *const u8 as *const c_char,
                timeZeroI.as_mut_ptr(),
                std::ptr::null_mut::<c_char>(),
                &mut tstat,
            ) != 0
            {
                timeZeroF[0] = 0.0f64;
                timeZeroI[0] = timeZeroF[0];
            } else if ffgkyd(
                fptr,
                b"TIMEZERF\0" as *const u8 as *const c_char,
                timeZeroF.as_mut_ptr(),
                std::ptr::null_mut::<c_char>(),
                &mut tstat,
            ) != 0
            {
                timeZeroF[0] = 0.0f64;
            }
        } else {
            timeZeroF[0] = 0.0f64;
        }
        match *fname.offset(0) as c_int {
            0 => {
                samefile = 1;
                hdunum = 1;
            }
            91 => {
                samefile = 1;
                i = 1;
                while *fname.offset(i as isize) as c_int != 0_i32
                    && *fname.offset(i as isize) as c_int != ']' as i32
                {
                    i += 1;
                    i;
                }
                if *fname.offset(i as isize) != 0 {
                    *fname.offset(i as isize) = 0;
                    fname = fname.offset(1);
                    fname;
                    ffexts(
                        fname,
                        &mut hdunum,
                        extname.as_mut_ptr(),
                        &mut extvers,
                        &mut movetotype,
                        xcol.as_mut_ptr(),
                        xexpr.as_mut_ptr(),
                        &mut (*lParse).status,
                    );
                    if *extname.as_mut_ptr() != 0 {
                        ffmnhd(
                            fptr,
                            movetotype,
                            extname.as_mut_ptr(),
                            extvers,
                            &mut (*lParse).status,
                        );
                        ffghdn(fptr, &mut hdunum);
                    } else if hdunum != 0 {
                        hdunum += 1;
                        ffmahd(fptr, hdunum, &mut hdutype, &mut (*lParse).status);
                    } else if (*lParse).status == 0 {
                        fits_parser_yyerror(
                            0 as yyscan_t,
                            lParse,
                            b"Cannot use primary array for GTI filter\0" as *const u8
                                as *const c_char as *mut c_char,
                        );
                        return -(1);
                    }
                } else {
                    fits_parser_yyerror(
                        0 as yyscan_t,
                        lParse,
                        b"File extension specifier lacks closing ']'\0" as *const u8
                            as *const c_char as *mut c_char,
                    );
                    return -(1);
                }
            }
            43 => {
                samefile = 1;
                hdunum =
                    atoi(std::ffi::CStr::from_ptr(fname).to_str().unwrap_or("0")).unwrap_or(0) + 1;
                if hdunum > 1 {
                    ffmahd(fptr, hdunum, &mut hdutype, &mut (*lParse).status);
                } else {
                    fits_parser_yyerror(
                        0 as yyscan_t,
                        lParse,
                        b"Cannot use primary array for GTI filter / GTIFIND\0" as *const u8
                            as *const c_char as *mut c_char,
                    );
                    return -(1);
                }
            }
            _ => {
                samefile = 0;
                let mut fptr_tmp: *mut Option<Box<fitsfile>> =
                    &mut None as *mut Option<Box<fitsfile>>;
                if ffopen(fptr_tmp, fname, 0, &mut (*lParse).status) == 0 {
                    ffghdn(fptr, &mut hdunum);
                }
            }
        }
        if (*lParse).status != 0 {
            return -(1);
        }
        if hdunum == 1 {
            loop {
                hdunum += 1;
                hdunum;
                if ffmahd(fptr, hdunum, &mut hdutype, &mut (*lParse).status) != 0 {
                    break;
                }
                if hdutype == 0 {
                    continue;
                }
                tstat = 0;
                if ffgkys(
                    fptr,
                    b"EXTNAME\0" as *const u8 as *const c_char,
                    extname.as_mut_ptr(),
                    std::ptr::null_mut::<c_char>(),
                    &mut tstat,
                ) != 0
                {
                    continue;
                }
                ffupch(extname.as_mut_ptr());
                if !(strstr(extname.as_mut_ptr(), b"GTI\0" as *const u8 as *const c_char)).is_null()
                {
                    break;
                }
            }
            if (*lParse).status != 0 {
                if (*lParse).status == 107 as c_int {
                    fits_parser_yyerror(
                        0 as yyscan_t,
                        lParse,
                        b"GTI extension not found in this file\0" as *const u8 as *const c_char
                            as *mut c_char,
                    );
                }
                return -(1);
            }
        }
        ffgcno(fptr, 0, start, &mut startCol, &mut (*lParse).status);
        ffgcno(fptr, 0, stop, &mut stopCol, &mut (*lParse).status);
        if (*lParse).status != 0 {
            return -(1);
        }
        tstat = 0;
        if ffgkyd(
            fptr,
            b"TIMEZERO\0" as *const u8 as *const c_char,
            timeZeroI.as_mut_ptr().offset(1 as c_int as isize),
            std::ptr::null_mut::<c_char>(),
            &mut tstat,
        ) != 0
        {
            tstat = 0;
            if ffgkyd(
                fptr,
                b"TIMEZERI\0" as *const u8 as *const c_char,
                timeZeroI.as_mut_ptr().offset(1 as c_int as isize),
                std::ptr::null_mut::<c_char>(),
                &mut tstat,
            ) != 0
            {
                timeZeroF[1] = 0.0f64;
                timeZeroI[1] = timeZeroF[1];
            } else if ffgkyd(
                fptr,
                b"TIMEZERF\0" as *const u8 as *const c_char,
                timeZeroF.as_mut_ptr().offset(1 as c_int as isize),
                std::ptr::null_mut::<c_char>(),
                &mut tstat,
            ) != 0
            {
                timeZeroF[1] = 0.0f64;
            }
        } else {
            timeZeroF[1] = 0.0f64;
        }
        n = Alloc_Node(lParse);
        if n >= 0 {
            this = ((*lParse).Nodes).offset(n as isize);
            (*this).SubNodes[1] = Node1;
            (*this).operation = Op as c_int;
            if Op as c_uint == gtifilt_fct as c_int as c_uint {
                (*this).nSubNodes = 2 as c_int;
                (*this).DoOp = Some(Do_GTI);
                (*this).ntype = fits_parser_yytokentype::BOOLEAN as c_int;
            } else if Op as c_uint == gtifind_fct as c_int as c_uint {
                (*this).nSubNodes = 2 as c_int;
                (*this).DoOp = Some(Do_GTI);
                (*this).ntype = fits_parser_yytokentype::LONG as c_int;
            } else {
                (*this).nSubNodes = 3 as c_int;
                (*this).DoOp = Some(Do_GTI_Over);
                (*this).ntype = fits_parser_yytokentype::DOUBLE as c_int;
            }
            that1 = ((*lParse).Nodes).offset(Node1 as isize);
            (*this).value.nelem = (*that1).value.nelem;
            (*this).value.naxis = (*that1).value.naxis;
            i = 0;
            while i < (*that1).value.naxis {
                (*this).value.naxes[i as usize] = (*that1).value.naxes[i as usize];
                i += 1;
                i;
            }
            if Op as c_uint == gtiover_fct as c_int as c_uint {
                (*this).SubNodes[2] = Node2;
                that2 = ((*lParse).Nodes).offset(Node2 as isize);
                if (*that1).value.nelem != (*that2).value.nelem {
                    fits_parser_yyerror(
                        0 as yyscan_t,
                        lParse,
                        b"Dimensions of TIME and TIME_STOP must match for GTIOVERLAP\0" as *const u8
                            as *const c_char as *mut c_char,
                    );
                    return -(1);
                }
            }
            (*this).SubNodes[0] = Node0;
            that0 = ((*lParse).Nodes).offset(Node0 as isize);
            (*that0).operation = -(1000 as c_int);
            (*that0).DoOp = None;
            (*that0).value.data.ptr = std::ptr::null_mut::<c_void>();
            if ffgkyj(
                fptr,
                b"NAXIS2\0" as *const u8 as *const c_char,
                &mut nrows,
                std::ptr::null_mut::<c_char>(),
                &mut (*lParse).status,
            ) != 0
            {
                return -(1);
            }
            (*that0).value.nelem = nrows;
            if nrows != 0 {
                let mut startptr: *mut c_double = std::ptr::null_mut::<c_double>();
                let mut stopptr: *mut c_double = std::ptr::null_mut::<c_double>();
                (*that0).value.data.dblptr = malloc(
                    ((2 as c_long * nrows) as c_ulong)
                        .wrapping_mul(::core::mem::size_of::<c_double>() as c_ulong)
                        .try_into()
                        .unwrap(),
                ) as *mut c_double;
                if ((*that0).value.data.dblptr).is_null() {
                    (*lParse).status = 113 as c_int;
                    return -(1);
                }
                startptr = (*that0).value.data.dblptr;
                stopptr = ((*that0).value.data.dblptr).offset(nrows as isize);
                ffgcvd(
                    fptr,
                    startCol,
                    1 as LONGLONG,
                    1 as LONGLONG,
                    nrows as LONGLONG,
                    0.0f64,
                    startptr,
                    &mut i,
                    &mut (*lParse).status,
                );
                ffgcvd(
                    fptr,
                    stopCol,
                    1 as LONGLONG,
                    1 as LONGLONG,
                    nrows as LONGLONG,
                    0.0f64,
                    stopptr,
                    &mut i,
                    &mut (*lParse).status,
                );
                if (*lParse).status != 0 {
                    free((*that0).value.data.dblptr as *mut c_void);
                    return -(1);
                }
                (*that0).ntype = 1;
                i = nrows as c_int;
                loop {
                    i -= 1;
                    if i == 0 {
                        break;
                    }
                    if !(*startptr.offset(i as isize) > *stopptr.offset(i as isize)
                        || *startptr.offset(i as isize) < *stopptr.offset((i - 1) as isize))
                    {
                        continue;
                    }
                    (*that0).ntype = 0;
                    break;
                }
                if (*that0).ntype != 1 && Op as c_uint == gtiover_fct as c_int as c_uint {
                    let mut errmsg: [c_char; 120] = [0; 120];
                    sprintf(
                        errmsg.as_mut_ptr(),
                        b"Input GTI must be time-ordered for GTIOVERLAP (row %ld)\0" as *const u8
                            as *const c_char,
                        i + 1,
                    );
                    fits_parser_yyerror(0 as yyscan_t, lParse, errmsg.as_mut_ptr());
                    return -(1);
                }
                dt = timeZeroI[1] - timeZeroI[0] + (timeZeroF[1] - timeZeroF[0]);
                timeSpan = *stopptr.offset((nrows - 1) as isize) - *startptr.offset(0);
                if timeSpan == 0.0 {
                    timeSpan = 1.0f64;
                }
                if fabs(dt / timeSpan) > 1e-12f64 {
                    i = 0;
                    while (i as c_long) < nrows {
                        *startptr.offset(i as isize) += dt;
                        *stopptr.offset(i as isize) += dt;
                        i += 1;
                        i;
                    }
                }
            }
            if (*((*lParse).Nodes).offset(Node1 as isize)).operation == -(1000 as c_int)
                && (Op as c_uint == gtifilt_fct as c_int as c_uint
                    || (*((*lParse).Nodes).offset(Node2 as isize)).operation == -(1000 as c_int))
            {
                ((*this).DoOp).expect("non-null function pointer")(lParse, this);
            }
        }
        if samefile != 0 {
            ffmahd(fptr, evthdu, &mut hdutype, &mut (*lParse).status);
        } else {
            ffclos(Some(Box::from_raw(fptr)), &mut (*lParse).status);
        }
        n
    }
}

unsafe fn New_REG(
    mut lParse: *mut ParseData,
    mut fname: *mut c_char,
    mut NodeX: c_int,
    mut NodeY: c_int,
    mut colNames: *mut c_char,
) -> c_int {
    unsafe {
        let mut this: *mut Node = std::ptr::null_mut::<Node>();
        let mut that0: *mut Node = std::ptr::null_mut::<Node>();
        let mut type_0: c_int = 0;
        let mut n: c_int = 0;
        let mut Node0: c_int = 0;
        let mut Xcol: c_int = 0;
        let mut Ycol: c_int = 0;
        let mut tstat: c_int = 0;
        let mut wcs: WCSdata = WCSdata {
            exists: false,
            xrefval: 0.,
            yrefval: 0.,
            xrefpix: 0.,
            yrefpix: 0.,
            xinc: 0.,
            yinc: 0.,
            rot: 0.,
            dtype: [0; 5],
        };
        let mut Rgn: *mut SAORegion = ptr::null_mut();
        let mut cX: *mut c_char = ptr::null_mut();
        let mut cY: *mut c_char = ptr::null_mut();
        let mut colVal: FITS_PARSER_YYSTYPE = FITS_PARSER_YYSTYPE { Node: 0 };
        if NodeX == -(99 as c_int) {
            type_0 = fits_parser_yyGetVariable(
                lParse,
                b"X\0" as *const u8 as *const c_char as *mut c_char,
                &mut colVal,
            );
            if type_0 == fits_parser_yytokentype::COLUMN as c_int {
                NodeX = New_Column(lParse, colVal.lng as c_int);
            } else {
                fits_parser_yyerror(
                    0 as yyscan_t,
                    lParse,
                    b"Could not build X column for REGFILTER\0" as *const u8 as *const c_char
                        as *mut c_char,
                );
                return -(1);
            }
        }
        if NodeY == -(99 as c_int) {
            type_0 = fits_parser_yyGetVariable(
                lParse,
                b"Y\0" as *const u8 as *const c_char as *mut c_char,
                &mut colVal,
            );
            if type_0 == fits_parser_yytokentype::COLUMN as c_int {
                NodeY = New_Column(lParse, colVal.lng as c_int);
            } else {
                fits_parser_yyerror(
                    0 as yyscan_t,
                    lParse,
                    b"Could not build Y column for REGFILTER\0" as *const u8 as *const c_char
                        as *mut c_char,
                );
                return -(1);
            }
        }
        NodeX = New_Unary(lParse, fits_parser_yytokentype::DOUBLE as c_int, 0, NodeX);
        NodeY = New_Unary(lParse, fits_parser_yytokentype::DOUBLE as c_int, 0, NodeY);
        Node0 = Alloc_Node(lParse);
        if NodeX < 0 || NodeY < 0 || Node0 < 0 {
            return -(1);
        }
        if Test_Dims(lParse, NodeX, NodeY) == 0 {
            fits_parser_yyerror(
                0 as yyscan_t,
                lParse,
                b"Dimensions of REGFILTER arguments are not compatible\0" as *const u8
                    as *const c_char as *mut c_char,
            );
            return -(1);
        }
        n = Alloc_Node(lParse);
        if n >= 0 {
            this = ((*lParse).Nodes).offset(n as isize);
            (*this).nSubNodes = 3 as c_int;
            (*this).SubNodes[0] = Node0;
            (*this).SubNodes[1] = NodeX;
            (*this).SubNodes[2] = NodeY;
            (*this).operation = regfilt_fct as c_int;
            (*this).DoOp = Some(Do_REG);
            (*this).ntype = fits_parser_yytokentype::BOOLEAN as c_int;
            (*this).value.nelem = 1;
            (*this).value.naxis = 1;
            (*this).value.naxes[0] = 1;
            Copy_Dims(lParse, n, NodeX);
            if ((*((*lParse).Nodes).offset(NodeX as isize)).value).nelem
                < ((*((*lParse).Nodes).offset(NodeY as isize)).value).nelem
            {
                Copy_Dims(lParse, n, NodeY);
            }
            that0 = ((*lParse).Nodes).offset(Node0 as isize);
            (*that0).operation = -(1000 as c_int);
            (*that0).DoOp = None;
            Ycol = 0;
            Xcol = Ycol;
            if *colNames != 0 {
                while *colNames as c_int == ' ' as i32 {
                    colNames = colNames.offset(1);
                    colNames;
                }
                cY = colNames;
                cX = cY;
                while *cY as c_int != 0 && *cY as c_int != ' ' as i32 && *cY as c_int != ',' as i32
                {
                    cY = cY.offset(1);
                    cY;
                }
                if *cY != 0 {
                    let fresh10 = cY;
                    cY = cY.offset(1);
                    *fresh10 = 0;
                }
                while *cY as c_int == ' ' as i32 {
                    cY = cY.offset(1);
                    cY;
                }
                if *cY == 0 {
                    fits_parser_yyerror(
                        0 as yyscan_t,
                        lParse,
                        b"Could not extract valid pair of column names from REGFILTER\0"
                            as *const u8 as *const c_char as *mut c_char,
                    );
                    Free_Last_Node(lParse);
                    return -(1);
                }
                ffgcno((*lParse).def_fptr, 0, cX, &mut Xcol, &mut (*lParse).status);
                ffgcno((*lParse).def_fptr, 0, cY, &mut Ycol, &mut (*lParse).status);
                if (*lParse).status != 0 {
                    fits_parser_yyerror(
                        0 as yyscan_t,
                        lParse,
                        b"Could not locate columns indicated for WCS info\0" as *const u8
                            as *const c_char as *mut c_char,
                    );
                    Free_Last_Node(lParse);
                    return -(1);
                }
            } else {
                Xcol = Locate_Col(lParse, ((*lParse).Nodes).offset(NodeX as isize));
                Ycol = Locate_Col(lParse, ((*lParse).Nodes).offset(NodeY as isize));
                if Xcol < 0 || Ycol < 0 {
                    fits_parser_yyerror(
                        0 as yyscan_t,
                        lParse,
                        b"Found multiple X/Y column references in REGFILTER\0" as *const u8
                            as *const c_char as *mut c_char,
                    );
                    Free_Last_Node(lParse);
                    return -(1);
                }
            }
            wcs.exists = false;
            if Xcol > 0 && Ycol > 0 {
                tstat = 0;
                ffgtcs(
                    (*lParse).def_fptr,
                    Xcol,
                    Ycol,
                    &mut wcs.xrefval,
                    &mut wcs.yrefval,
                    &mut wcs.xrefpix,
                    &mut wcs.yrefpix,
                    &mut wcs.xinc,
                    &mut wcs.yinc,
                    &mut wcs.rot,
                    &mut (wcs.dtype),
                    &mut tstat,
                );
                if tstat == 505 as c_int {
                    wcs.exists = false;
                } else if tstat != 0 {
                    (*lParse).status = tstat;
                    Free_Last_Node(lParse);
                    return -(1);
                } else {
                    wcs.exists = true;
                }
            }

            let fname_slice = CStr::from_ptr(fname);
            let rgn_box = Box::from_raw(Rgn);
            fits_read_rgnfile(
                cast_slice(fname_slice.to_bytes_with_nul()),
                &mut wcs,
                &mut Some(rgn_box),
                &mut (*lParse).status,
            );
            if (*lParse).status != 0 {
                Free_Last_Node(lParse);
                return -(1);
            }
            (*that0).value.data.ptr = Rgn as *mut c_void;
            if (*((*lParse).Nodes).offset(NodeX as isize)).operation == -(1000 as c_int)
                && (*((*lParse).Nodes).offset(NodeY as isize)).operation == -(1000 as c_int)
            {
                ((*this).DoOp).expect("non-null function pointer")(lParse, this);
            }
        }
        n
    }
}
unsafe fn New_Vector(mut lParse: *mut ParseData, mut subNode: c_int) -> c_int {
    unsafe {
        let mut this: *mut Node = std::ptr::null_mut::<Node>();
        let mut that: *mut Node = std::ptr::null_mut::<Node>();
        let mut n: c_int = 0;
        n = Alloc_Node(lParse);
        if n >= 0 {
            this = ((*lParse).Nodes).offset(n as isize);
            that = ((*lParse).Nodes).offset(subNode as isize);
            (*this).ntype = (*that).ntype;
            (*this).nSubNodes = 1;
            (*this).SubNodes[0] = subNode;
            (*this).operation = '{' as i32;
            (*this).DoOp = Some(Do_Vector);
        }
        n
    }
}
unsafe fn Close_Vec(mut lParse: *mut ParseData, mut vecNode: c_int) -> c_int {
    unsafe {
        let mut this: *mut Node = std::ptr::null_mut::<Node>();
        let mut n: c_int = 0;
        let mut nelem: c_int = 0;
        this = ((*lParse).Nodes).offset(vecNode as isize);
        n = 0;
        while n < (*this).nSubNodes {
            if (*((*lParse).Nodes).offset((*this).SubNodes[n as usize] as isize)).ntype
                != (*this).ntype
            {
                (*this).SubNodes[n as usize] =
                    New_Unary(lParse, (*this).ntype, 0, (*this).SubNodes[n as usize]);
                if (*this).SubNodes[n as usize] < 0 {
                    return -(1);
                }
            }
            nelem = (nelem as c_long
                + (*((*lParse).Nodes).offset((*this).SubNodes[n as usize] as isize))
                    .value
                    .nelem) as c_int;
            n += 1;
            n;
        }
        (*this).value.naxis = 1;
        (*this).value.nelem = nelem as c_long;
        (*this).value.naxes[0] = nelem as c_long;
        vecNode
    }
}
unsafe fn New_Array(mut lParse: *mut ParseData, mut valueNode: c_int, mut dimNode: c_int) -> c_int {
    unsafe {
        let mut dims: *mut Node = std::ptr::null_mut::<Node>();
        let mut naxis: c_long = 0;
        let mut nelem: c_long = 0;
        let mut naxes: [c_long; 5] = [0; 5];
        let mut this: *mut Node = std::ptr::null_mut::<Node>();
        let mut n: c_int = 0;
        let mut i: c_int = 0;
        if valueNode < 0 || dimNode < 0 {
            return -(1);
        }
        dims = &mut *((*lParse).Nodes).offset(dimNode as isize) as *mut Node;
        i = 0;
        while i < 5 as c_int {
            naxes[i as usize] = 1;
            i += 1;
            i;
        }
        if (*((*lParse).Nodes).offset(dimNode as isize)).operation == -(1000 as c_int) {
            if (*((*lParse).Nodes).offset(dimNode as isize)).ntype
                != fits_parser_yytokentype::LONG as c_int
            {
                dimNode = New_Unary(lParse, fits_parser_yytokentype::LONG as c_int, 0, dimNode);
            }
            if dimNode < 0 {
                return -(1);
            }
            naxis = 1;
            naxes[0] = ((*((*lParse).Nodes).offset(dimNode as isize)).value)
                .data
                .lng;
        } else if (*((*lParse).Nodes).offset(dimNode as isize)).operation == '{' as i32 {
            if (*dims).nSubNodes > 5 as c_int {
                fits_parser_yyerror(
                    0 as yyscan_t,
                    lParse,
                    b"ARRAY(V,{...}) number of dimensions must not exceed 5\0" as *const u8
                        as *const c_char as *mut c_char,
                );
                return -(1);
            }
            naxis = (*dims).nSubNodes as c_long;
            i = 0;
            while i < (*dims).nSubNodes {
                if (*((*lParse).Nodes).offset((*dims).SubNodes[i as usize] as isize)).ntype
                    != fits_parser_yytokentype::LONG as c_int
                {
                    (*dims).SubNodes[i as usize] = New_Unary(
                        lParse,
                        fits_parser_yytokentype::LONG as c_int,
                        0,
                        (*dims).SubNodes[i as usize],
                    );
                    if (*dims).SubNodes[i as usize] < 0 {
                        return -(1);
                    }
                }
                naxes[i as usize] = (*((*lParse).Nodes)
                    .offset((*dims).SubNodes[i as usize] as isize))
                .value
                .data
                .lng;
                i += 1;
                i;
            }
        } else {
            fits_parser_yyerror(
                0 as yyscan_t,
                lParse,
                b"ARRAY(V,dims) dims must be either integer or const vector\0" as *const u8
                    as *const c_char as *mut c_char,
            );
            return -(1);
        }
        nelem = 1;
        i = 0;
        while (i as c_long) < naxis {
            if naxes[i as usize] <= 0 {
                fits_parser_yyerror(
                    0 as yyscan_t,
                    lParse,
                    b"ARRAY(V,dims) must have positive dimensions\0" as *const u8 as *const c_char
                        as *mut c_char,
                );
                return -(1);
            }
            nelem *= naxes[i as usize];
            i += 1;
            i;
        }
        if !(((*((*lParse).Nodes).offset(valueNode as isize)).value).nelem == nelem && nelem > 1) {
            if ((*((*lParse).Nodes).offset(valueNode as isize)).value).nelem > 1 && nelem > 1 {
                fits_parser_yyerror(
                    0 as yyscan_t,
                    lParse,
                    b"ARRAY(V,d) mismatch between number of elements in V and d\0" as *const u8
                        as *const c_char as *mut c_char,
                );
                return -(1);
            } else if ((*((*lParse).Nodes).offset(valueNode as isize)).value).nelem > 1 {
                fits_parser_yyerror(
                    0 as yyscan_t,
                    lParse,
                    b"ARRAY(V,n) value V must have vector dimension of 1\0" as *const u8
                        as *const c_char as *mut c_char,
                );
                return -(1);
            }
        }
        n = Alloc_Node(lParse);
        if n >= 0 {
            this = ((*lParse).Nodes).offset(n as isize);
            (*this).operation = array_fct as c_int;
            (*this).nSubNodes = 1;
            (*this).SubNodes[0] = valueNode;
            (*this).ntype = (*((*lParse).Nodes).offset(valueNode as isize)).ntype;
            (*this).value.nelem = nelem;
            (*this).value.naxis = naxis as c_int;
            i = 0;
            while (i as c_long) < naxis {
                (*this).value.naxes[i as usize] = naxes[i as usize];
                i += 1;
                i;
            }
            (*this).DoOp = Some(Do_Array);
        }
        n
    }
}
unsafe fn Locate_Col(mut lParse: *mut ParseData, mut this: *mut Node) -> c_int {
    unsafe {
        let mut that: *mut Node = std::ptr::null_mut::<Node>();
        let mut i: c_int = 0;
        let mut col: c_int = 0;
        let mut newCol: c_int = 0;
        let mut nfound: c_int = 0;
        if (*this).nSubNodes == 0 && (*this).operation <= 0 && (*this).operation != -(1000 as c_int)
        {
            return (*((*lParse).colData).offset(-(*this).operation as isize)).colnum;
        }
        i = 0;
        while i < (*this).nSubNodes {
            that = ((*lParse).Nodes).offset((*this).SubNodes[i as usize] as isize);
            if (*that).operation > 0 {
                newCol = Locate_Col(lParse, that);
                if newCol <= 0 {
                    nfound += -newCol;
                } else if nfound == 0 {
                    col = newCol;
                    nfound += 1;
                    nfound;
                } else if col != newCol {
                    nfound += 1;
                    nfound;
                }
            } else if (*that).operation != -(1000 as c_int) {
                newCol = (*((*lParse).colData).offset(-(*that).operation as isize)).colnum;
                if nfound == 0 {
                    col = newCol;
                    nfound += 1;
                    nfound;
                } else if col != newCol {
                    nfound += 1;
                    nfound;
                }
            }
            i += 1;
            i;
        }
        if nfound != 1 { -nfound } else { col }
    }
}
unsafe fn Test_Dims(mut lParse: *mut ParseData, mut Node1: c_int, mut Node2: c_int) -> c_int {
    unsafe {
        let mut that1: *mut Node = std::ptr::null_mut::<Node>();
        let mut that2: *mut Node = std::ptr::null_mut::<Node>();
        let mut valid: c_int = 0;
        let mut i: c_int = 0;
        if Node1 < 0 || Node2 < 0 {
            return 0;
        }
        that1 = ((*lParse).Nodes).offset(Node1 as isize);
        that2 = ((*lParse).Nodes).offset(Node2 as isize);
        if (*that1).value.nelem == 1 || (*that2).value.nelem == 1 {
            valid = 1;
        } else if (*that1).ntype == (*that2).ntype
            && (*that1).value.nelem == (*that2).value.nelem
            && (*that1).value.naxis == (*that2).value.naxis
        {
            valid = 1;
            i = 0;
            while i < (*that1).value.naxis {
                if (*that1).value.naxes[i as usize] != (*that2).value.naxes[i as usize] {
                    valid = 0;
                }
                i += 1;
                i;
            }
        } else {
            valid = 0;
        }
        valid
    }
}
unsafe fn Copy_Dims(mut lParse: *mut ParseData, mut Node1: c_int, mut Node2: c_int) {
    unsafe {
        let mut that1: *mut Node = std::ptr::null_mut::<Node>();
        let mut that2: *mut Node = std::ptr::null_mut::<Node>();
        let mut i: c_int = 0;
        if Node1 < 0 || Node2 < 0 {
            return;
        }
        that1 = ((*lParse).Nodes).offset(Node1 as isize);
        that2 = ((*lParse).Nodes).offset(Node2 as isize);
        (*that1).value.nelem = (*that2).value.nelem;
        (*that1).value.naxis = (*that2).value.naxis;
        i = 0;
        while i < (*that2).value.naxis {
            (*that1).value.naxes[i as usize] = (*that2).value.naxes[i as usize];
            i += 1;
            i;
        }
    }
}

pub unsafe fn Evaluate_Parser(mut lParse: *mut ParseData, mut firstRow: c_long, mut nRows: c_long) {
    unsafe {
        let mut i: c_int = 0;
        let mut column: c_int = 0;
        let mut offset: c_long = 0;
        let mut rowOffset: c_long = 0;
        static mut rand_initialized: c_int = 0;
        if rand_initialized == 0 {
            simplerng_srand(time(std::ptr::null_mut::<time_t>()) as c_uint);
            rand_initialized = 1;
        }
        (*lParse).firstRow = firstRow;
        (*lParse).nRows = nRows;
        rowOffset = firstRow - (*lParse).firstDataRow;
        i = 0;
        while i < (*lParse).nNodes {
            if !((*((*lParse).Nodes).offset(i as isize)).operation > 0
                || (*((*lParse).Nodes).offset(i as isize)).operation == -(1000 as c_int))
            {
                column = -(*((*lParse).Nodes).offset(i as isize)).operation;
                offset = (*((*lParse).varData).offset(column as isize)).nelem * rowOffset;
                let fresh11 = &mut ((*((*lParse).Nodes).offset(i as isize)).value).undef;
                *fresh11 =
                    ((*((*lParse).varData).offset(column as isize)).undef).offset(offset as isize);
                match (*((*lParse).Nodes).offset(i as isize)).ntype {
                    262 => {
                        let fresh12 =
                            &mut ((*((*lParse).Nodes).offset(i as isize)).value).data.strptr;
                        *fresh12 = ((*((*lParse).varData).offset(column as isize)).data
                            as *mut *mut c_char)
                            .offset(rowOffset as isize);
                        let fresh13 = &mut ((*((*lParse).Nodes).offset(i as isize)).value).undef;
                        *fresh13 = ptr::null_mut();
                    }
                    261 => {
                        let fresh14 =
                            &mut ((*((*lParse).Nodes).offset(i as isize)).value).data.strptr;
                        *fresh14 = ((*((*lParse).varData).offset(column as isize)).data
                            as *mut *mut c_char)
                            .offset(rowOffset as isize);
                        let fresh15 = &mut ((*((*lParse).Nodes).offset(i as isize)).value).undef;
                        *fresh15 = ((*((*lParse).varData).offset(column as isize)).undef)
                            .offset(rowOffset as isize);
                    }
                    258 => {
                        let fresh16 =
                            &mut ((*((*lParse).Nodes).offset(i as isize)).value).data.logptr;
                        *fresh16 = ((*((*lParse).varData).offset(column as isize)).data as *const _
                            as *mut c_char)
                            .offset(offset as isize);
                    }
                    259 => {
                        let fresh17 =
                            &mut ((*((*lParse).Nodes).offset(i as isize)).value).data.lngptr;
                        *fresh17 = ((*((*lParse).varData).offset(column as isize)).data as *const _
                            as *mut c_long)
                            .offset(offset as isize);
                    }
                    260 => {
                        let fresh18 =
                            &mut ((*((*lParse).Nodes).offset(i as isize)).value).data.dblptr;
                        *fresh18 = ((*((*lParse).varData).offset(column as isize)).data as *const _
                            as *mut c_double)
                            .offset(offset as isize);
                    }
                    _ => {}
                }
            }
            i += 1;
            i;
        }
        Evaluate_Node(lParse, (*lParse).resultNode);
    }
}
unsafe fn Evaluate_Node(mut lParse: *mut ParseData, mut thisNode: c_int) {
    unsafe {
        let mut this: *mut Node = std::ptr::null_mut::<Node>();
        let mut i: c_int = 0;
        if (*lParse).status != 0 {
            return;
        }
        this = ((*lParse).Nodes).offset(thisNode as isize);
        if (*this).operation > 0 {
            i = (*this).nSubNodes;
            loop {
                let fresh19 = i;
                i -= 1;
                if fresh19 == 0 {
                    break;
                }
                Evaluate_Node(lParse, (*this).SubNodes[i as usize]);
                if (*lParse).status != 0 {
                    return;
                }
            }
            ((*this).DoOp).expect("non-null function pointer")(lParse, this);
        }
    }
}
unsafe fn Allocate_Ptrs(mut lParse: *mut ParseData, mut this: *mut Node) {
    unsafe {
        let mut elem: c_long = 0;
        let mut row: c_long = 0;
        let mut size: c_long = 0;
        if (*this).ntype == fits_parser_yytokentype::BITSTR as c_int
            || (*this).ntype == fits_parser_yytokentype::STRING as c_int
        {
            (*this).value.data.strptr = malloc(
                ((*lParse).nRows as c_ulong)
                    .wrapping_mul(::core::mem::size_of::<*mut c_char>() as c_ulong)
                    .try_into()
                    .unwrap(),
            ) as *mut *mut c_char;
            if !((*this).value.data.strptr).is_null() {
                let fresh20 = &mut *((*this).value.data.strptr).offset(0);
                *fresh20 = malloc(
                    (((*lParse).nRows * ((*this).value.nelem + 2 as c_long)) as c_ulong)
                        .wrapping_mul(::core::mem::size_of::<c_char>() as c_ulong)
                        .try_into()
                        .unwrap(),
                ) as *mut c_char;
                if !(*((*this).value.data.strptr).offset(0)).is_null() {
                    row = 0;
                    loop {
                        row += 1;
                        if row >= (*lParse).nRows {
                            break;
                        }
                        let fresh21 = &mut *((*this).value.data.strptr).offset(row as isize);
                        *fresh21 = (*((*this).value.data.strptr).offset((row - 1) as isize))
                            .offset((*this).value.nelem as isize)
                            .offset(1 as c_int as isize);
                    }
                    if (*this).ntype == fits_parser_yytokentype::STRING as c_int {
                        (*this).value.undef = (*((*this).value.data.strptr)
                            .offset((row - 1) as isize))
                        .offset((*this).value.nelem as isize)
                        .offset(1 as c_int as isize);
                    } else {
                        (*this).value.undef = ptr::null_mut();
                    }
                } else {
                    (*lParse).status = 113 as c_int;
                    free((*this).value.data.strptr as *mut c_void);
                }
            } else {
                (*lParse).status = 113 as c_int;
            }
        } else {
            elem = (*this).value.nelem * (*lParse).nRows;
            match (*this).ntype {
                260 => {
                    size = ::core::mem::size_of::<c_double>() as c_ulong as c_long;
                }
                259 => {
                    size = ::core::mem::size_of::<c_long>() as c_ulong as c_long;
                }
                258 => {
                    size = ::core::mem::size_of::<c_char>() as c_ulong as c_long;
                }
                _ => {
                    size = 1;
                }
            }
            (*this).value.data.ptr = calloc(
                ((size + 1) as c_ulong).try_into().unwrap(),
                (elem as c_ulong).try_into().unwrap(),
            );
            if ((*this).value.data.ptr).is_null() {
                (*lParse).status = 113 as c_int;
            } else {
                (*this).value.undef =
                    ((*this).value.data.ptr as *mut c_char).offset((elem * size) as isize);
            }
        };
    }
}
unsafe fn Do_Unary(mut lParse: *mut ParseData, mut this: *mut Node) {
    unsafe {
        let mut that: *mut Node = std::ptr::null_mut::<Node>();
        let mut elem: c_long = 0;
        that = ((*lParse).Nodes).offset((*this).SubNodes[0] as isize);
        if (*that).operation == -(1000 as c_int) {
            match (*this).operation {
                x if x == fits_parser_yytokentype::DOUBLE as c_int
                    || x == fits_parser_yytokentype::FLTCAST as c_int =>
                {
                    if (*that).ntype == fits_parser_yytokentype::LONG as c_int {
                        (*this).value.data.dbl = (*that).value.data.lng as c_double;
                    } else if (*that).ntype == fits_parser_yytokentype::BOOLEAN as c_int {
                        (*this).value.data.dbl = if (*that).value.data.log as c_int != 0 {
                            1.0f64
                        } else {
                            0.0f64
                        };
                    }
                }
                x if x == fits_parser_yytokentype::LONG as c_int
                    || x == fits_parser_yytokentype::INTCAST as c_int =>
                {
                    if (*that).ntype == fits_parser_yytokentype::DOUBLE as c_int {
                        (*this).value.data.lng = (*that).value.data.dbl as c_long;
                    } else if (*that).ntype == fits_parser_yytokentype::BOOLEAN as c_int {
                        (*this).value.data.lng = if (*that).value.data.log as c_int != 0 {
                            1
                        } else {
                            0
                        };
                    }
                }
                x if x == fits_parser_yytokentype::BOOLEAN as c_int => {
                    if (*that).ntype == fits_parser_yytokentype::DOUBLE as c_int {
                        (*this).value.data.log = if (*that).value.data.dbl != 0.0f64 {
                            1
                        } else {
                            0
                        };
                    } else if (*that).ntype == fits_parser_yytokentype::LONG as c_int {
                        (*this).value.data.log = if (*that).value.data.lng != 0 { 1 } else { 0 };
                    }
                }
                x if x == fits_parser_yytokentype::UMINUS as c_int => {
                    if (*that).ntype == fits_parser_yytokentype::DOUBLE as c_int {
                        (*this).value.data.dbl = -(*that).value.data.dbl;
                    } else if (*that).ntype == fits_parser_yytokentype::LONG as c_int {
                        (*this).value.data.lng = -(*that).value.data.lng;
                    }
                }
                x if x == fits_parser_yytokentype::NOT as c_int => {
                    if (*that).ntype == fits_parser_yytokentype::BOOLEAN as c_int {
                        (*this).value.data.log = if (*that).value.data.log == 0 { 1 } else { 0 };
                    } else if (*that).ntype == fits_parser_yytokentype::BITSTR as c_int {
                        bitnot(
                            ((*this).value.data.astr).as_mut_ptr(),
                            ((*that).value.data.astr).as_mut_ptr(),
                        );
                    }
                }
                _ => {}
            }
            (*this).operation = -(1000 as c_int);
        } else {
            Allocate_Ptrs(lParse, this);
            if (*lParse).status == 0 {
                if (*this).ntype != fits_parser_yytokentype::BITSTR as c_int {
                    elem = (*lParse).nRows;
                    if (*this).ntype != fits_parser_yytokentype::STRING as c_int {
                        elem *= (*this).value.nelem;
                    }
                    loop {
                        let fresh22 = elem;
                        elem -= 1;
                        if fresh22 == 0 {
                            break;
                        }
                        *((*this).value.undef).offset(elem as isize) =
                            *((*that).value.undef).offset(elem as isize);
                    }
                }
                elem = (*lParse).nRows * (*this).value.nelem;
                match (*this).operation {
                    258 => {
                        if (*that).ntype == fits_parser_yytokentype::DOUBLE as c_int {
                            loop {
                                let fresh23 = elem;
                                elem -= 1;
                                if fresh23 == 0 {
                                    break;
                                }
                                *((*this).value.data.logptr).offset(elem as isize) =
                                    if *((*that).value.data.dblptr).offset(elem as isize) != 0.0f64
                                    {
                                        1
                                    } else {
                                        0
                                    };
                            }
                        } else if (*that).ntype == fits_parser_yytokentype::LONG as c_int {
                            loop {
                                let fresh24 = elem;
                                elem -= 1;
                                if fresh24 == 0 {
                                    break;
                                }
                                *((*this).value.data.logptr).offset(elem as isize) =
                                    if *((*that).value.data.lngptr).offset(elem as isize) != 0 {
                                        1
                                    } else {
                                        0
                                    };
                            }
                        }
                    }
                    260 | 289 => {
                        if (*that).ntype == fits_parser_yytokentype::LONG as c_int {
                            loop {
                                let fresh25 = elem;
                                elem -= 1;
                                if fresh25 == 0 {
                                    break;
                                }
                                *((*this).value.data.dblptr).offset(elem as isize) =
                                    *((*that).value.data.lngptr).offset(elem as isize) as c_double;
                            }
                        } else if (*that).ntype == fits_parser_yytokentype::BOOLEAN as c_int {
                            loop {
                                let fresh26 = elem;
                                elem -= 1;
                                if fresh26 == 0 {
                                    break;
                                }
                                *((*this).value.data.dblptr).offset(elem as isize) =
                                    if *((*that).value.data.logptr).offset(elem as isize) as c_int
                                        != 0
                                    {
                                        1.0f64
                                    } else {
                                        0.0f64
                                    };
                            }
                        }
                    }
                    259 | 288 => {
                        if (*that).ntype == fits_parser_yytokentype::DOUBLE as c_int {
                            loop {
                                let fresh27 = elem;
                                elem -= 1;
                                if fresh27 == 0 {
                                    break;
                                }
                                *((*this).value.data.lngptr).offset(elem as isize) =
                                    *((*that).value.data.dblptr).offset(elem as isize) as c_long;
                            }
                        } else if (*that).ntype == fits_parser_yytokentype::BOOLEAN as c_int {
                            loop {
                                let fresh28 = elem;
                                elem -= 1;
                                if fresh28 == 0 {
                                    break;
                                }
                                *((*this).value.data.lngptr).offset(elem as isize) =
                                    if *((*that).value.data.logptr).offset(elem as isize) as c_int
                                        != 0
                                    {
                                        1
                                    } else {
                                        0
                                    };
                            }
                        }
                    }
                    290 => {
                        if (*that).ntype == fits_parser_yytokentype::DOUBLE as c_int {
                            loop {
                                let fresh29 = elem;
                                elem -= 1;
                                if fresh29 == 0 {
                                    break;
                                }
                                *((*this).value.data.dblptr).offset(elem as isize) =
                                    -*((*that).value.data.dblptr).offset(elem as isize);
                            }
                        } else if (*that).ntype == fits_parser_yytokentype::LONG as c_int {
                            loop {
                                let fresh30 = elem;
                                elem -= 1;
                                if fresh30 == 0 {
                                    break;
                                }
                                *((*this).value.data.lngptr).offset(elem as isize) =
                                    -*((*that).value.data.lngptr).offset(elem as isize);
                            }
                        }
                    }
                    287 => {
                        if (*that).ntype == fits_parser_yytokentype::BOOLEAN as c_int {
                            loop {
                                let fresh31 = elem;
                                elem -= 1;
                                if fresh31 == 0 {
                                    break;
                                }
                                *((*this).value.data.logptr).offset(elem as isize) =
                                    (*((*that).value.data.logptr).offset(elem as isize) == 0)
                                        as c_int as c_char;
                            }
                        } else if (*that).ntype == fits_parser_yytokentype::BITSTR as c_int {
                            elem = (*lParse).nRows;
                            loop {
                                let fresh32 = elem;
                                elem -= 1;
                                if fresh32 == 0 {
                                    break;
                                }
                                bitnot(
                                    *((*this).value.data.strptr).offset(elem as isize),
                                    *((*that).value.data.strptr).offset(elem as isize),
                                );
                            }
                        }
                    }
                    _ => {}
                }
            }
        }
        if (*that).operation > 0 {
            free((*that).value.data.ptr);
        }
    }
}
unsafe fn Do_Offset(mut lParse: *mut ParseData, mut this: *mut Node) {
    unsafe {
        let mut col: *mut Node = std::ptr::null_mut::<Node>();
        let mut fRow: c_long = 0;
        let mut nRowOverlap: c_long = 0;
        let mut nRowReload: c_long = 0;
        let mut rowOffset: c_long = 0;
        let mut nelem: c_long = 0;
        let mut elem: c_long = 0;
        let mut offset: c_long = 0;
        let mut nRealElem: c_long = 0;
        let mut status: c_int = 0;
        col = ((*lParse).Nodes).offset((*this).SubNodes[0] as isize);
        rowOffset = (*((*lParse).Nodes).offset((*this).SubNodes[1] as isize))
            .value
            .data
            .lng;
        Allocate_Ptrs(lParse, this);
        fRow = (*lParse).firstRow + rowOffset;
        if (*this).ntype == fits_parser_yytokentype::STRING as c_int
            || (*this).ntype == fits_parser_yytokentype::BITSTR as c_int
        {
            nRealElem = 1;
        } else {
            nRealElem = (*this).value.nelem;
        }
        nelem = nRealElem;
        if fRow < (*lParse).firstDataRow {
            nRowReload = (*lParse).firstDataRow - fRow;
            if nRowReload > (*lParse).nRows {
                nRowReload = (*lParse).nRows;
            }
            nRowOverlap = (*lParse).nRows - nRowReload;
            offset = 0;
            while fRow < 1 && nRowReload > 0 {
                if (*this).ntype == fits_parser_yytokentype::BITSTR as c_int {
                    nelem = (*this).value.nelem;
                    *(*((*this).value.data.strptr).offset(offset as isize))
                        .offset(nelem as isize) = 0;
                    loop {
                        let fresh33 = nelem;
                        nelem -= 1;
                        if fresh33 == 0 {
                            break;
                        }
                        *(*((*this).value.data.strptr).offset(offset as isize))
                            .offset(nelem as isize) = b'0' as c_char;
                    }
                    offset += 1;
                    offset;
                } else {
                    loop {
                        let fresh34 = nelem;
                        nelem -= 1;
                        if fresh34 == 0 {
                            break;
                        }
                        let fresh35 = offset;
                        offset += 1;
                        *((*this).value.undef).offset(fresh35 as isize) = 1;
                    }
                }
                nelem = nRealElem;
                fRow += 1;
                fRow;
                nRowReload -= 1;
                nRowReload;
            }
        } else if fRow + (*lParse).nRows > (*lParse).firstDataRow + (*lParse).nDataRows {
            nRowReload = fRow + (*lParse).nRows - ((*lParse).firstDataRow + (*lParse).nDataRows);
            if nRowReload > (*lParse).nRows {
                nRowReload = (*lParse).nRows;
            } else {
                fRow = (*lParse).firstDataRow + (*lParse).nDataRows;
            }
            nRowOverlap = (*lParse).nRows - nRowReload;
            offset = nRowOverlap * nelem;
            elem = (*lParse).nRows * nelem;
            while fRow + nRowReload > (*lParse).totalRows && nRowReload > 0 {
                if (*this).ntype == fits_parser_yytokentype::BITSTR as c_int {
                    nelem = (*this).value.nelem;
                    elem -= 1;
                    elem;
                    *(*((*this).value.data.strptr).offset(elem as isize)).offset(nelem as isize) =
                        0;
                    loop {
                        let fresh36 = nelem;
                        nelem -= 1;
                        if fresh36 == 0 {
                            break;
                        }
                        *(*((*this).value.data.strptr).offset(elem as isize))
                            .offset(nelem as isize) = b'0' as c_char;
                    }
                } else {
                    loop {
                        let fresh37 = nelem;
                        nelem -= 1;
                        if fresh37 == 0 {
                            break;
                        }
                        elem -= 1;
                        *((*this).value.undef).offset(elem as isize) = 1;
                    }
                }
                nelem = nRealElem;
                nRowReload -= 1;
                nRowReload;
            }
        } else {
            nRowReload = 0;
            nRowOverlap = (*lParse).nRows;
            offset = 0;
        }
        if nRowReload > 0 {
            match (*this).ntype {
                262 | 261 => {
                    status = ((*lParse).loadData).expect("non-null function pointer")(
                        lParse,
                        -(*col).operation,
                        fRow,
                        nRowReload,
                        ((*this).value.data.strptr).offset(offset as isize) as *mut c_void,
                        ((*this).value.undef).offset(offset as isize),
                    );
                }
                258 => {
                    status = ((*lParse).loadData).expect("non-null function pointer")(
                        lParse,
                        -(*col).operation,
                        fRow,
                        nRowReload,
                        ((*this).value.data.logptr).offset(offset as isize) as *mut c_void,
                        ((*this).value.undef).offset(offset as isize),
                    );
                }
                259 => {
                    status = ((*lParse).loadData).expect("non-null function pointer")(
                        lParse,
                        -(*col).operation,
                        fRow,
                        nRowReload,
                        ((*this).value.data.lngptr).offset(offset as isize) as *mut c_void,
                        ((*this).value.undef).offset(offset as isize),
                    );
                }
                260 => {
                    status = ((*lParse).loadData).expect("non-null function pointer")(
                        lParse,
                        -(*col).operation,
                        fRow,
                        nRowReload,
                        ((*this).value.data.dblptr).offset(offset as isize) as *mut c_void,
                        ((*this).value.undef).offset(offset as isize),
                    );
                }
                _ => {}
            }
        }
        if nRowOverlap <= 0 {
            return;
        }
        if rowOffset > 0 {
            elem = nRowOverlap * nelem;
        } else {
            elem = (*lParse).nRows * nelem;
        }
        offset = nelem * rowOffset;
        loop {
            let fresh38 = nRowOverlap;
            nRowOverlap -= 1;
            if !(fresh38 != 0 && (*lParse).status == 0) {
                break;
            }
            loop {
                let fresh39 = nelem;
                nelem -= 1;
                if !(fresh39 != 0 && (*lParse).status == 0) {
                    break;
                }
                elem -= 1;
                elem;
                if (*this).ntype != fits_parser_yytokentype::BITSTR as c_int {
                    *((*this).value.undef).offset(elem as isize) =
                        *((*col).value.undef).offset((elem + offset) as isize);
                }
                match (*this).ntype {
                    262 => {
                        strcpy(
                            *((*this).value.data.strptr).offset(elem as isize),
                            *((*col).value.data.strptr).offset((elem + offset) as isize),
                        );
                    }
                    261 => {
                        strcpy(
                            *((*this).value.data.strptr).offset(elem as isize),
                            *((*col).value.data.strptr).offset((elem + offset) as isize),
                        );
                    }
                    258 => {
                        *((*this).value.data.logptr).offset(elem as isize) =
                            *((*col).value.data.logptr).offset((elem + offset) as isize);
                    }
                    259 => {
                        *((*this).value.data.lngptr).offset(elem as isize) =
                            *((*col).value.data.lngptr).offset((elem + offset) as isize);
                    }
                    260 => {
                        *((*this).value.data.dblptr).offset(elem as isize) =
                            *((*col).value.data.dblptr).offset((elem + offset) as isize);
                    }
                    _ => {}
                }
            }
            nelem = nRealElem;
        }
    }
}
unsafe fn Do_BinOp_bit(mut lParse: *mut ParseData, mut this: *mut Node) {
    unsafe {
        let mut that1: *mut Node = std::ptr::null_mut::<Node>();
        let mut that2: *mut Node = std::ptr::null_mut::<Node>();
        let mut sptr1: *mut c_char = ptr::null_mut();
        let mut sptr2: *mut c_char = ptr::null_mut();
        let mut const1: c_int = 0;
        let mut const2: c_int = 0;
        let mut rows: c_long = 0;
        that1 = ((*lParse).Nodes).offset((*this).SubNodes[0] as isize);
        that2 = ((*lParse).Nodes).offset((*this).SubNodes[1] as isize);
        const1 = ((*that1).operation == -(1000 as c_int)) as c_int;
        const2 = ((*that2).operation == -(1000 as c_int)) as c_int;
        sptr1 = if const1 != 0 {
            ((*that1).value.data.astr).as_mut_ptr()
        } else {
            std::ptr::null_mut::<c_char>()
        };
        sptr2 = if const2 != 0 {
            ((*that2).value.data.astr).as_mut_ptr()
        } else {
            std::ptr::null_mut::<c_char>()
        };
        if const1 != 0 && const2 != 0 {
            match (*this).operation {
                280 => {
                    (*this).value.data.log = if bitcmp(sptr1, sptr2) == 0 { 1 } else { 0 };
                }
                279 => {
                    (*this).value.data.log = bitcmp(sptr1, sptr2);
                }
                281..=284 => {
                    (*this).value.data.log = bitlgte(sptr1, (*this).operation, sptr2);
                }
                124 => {
                    bitor(((*this).value.data.astr).as_mut_ptr(), sptr1, sptr2);
                }
                38 => {
                    bitand(((*this).value.data.astr).as_mut_ptr(), sptr1, sptr2);
                }
                43 => {
                    strcpy(((*this).value.data.astr).as_mut_ptr(), sptr1);
                    strcat(((*this).value.data.astr).as_mut_ptr(), sptr2);
                }
                291 => {
                    (*this).value.data.lng = 0;
                    while *sptr1 != 0 {
                        if *sptr1 as c_int == '1' as i32 {
                            (*this).value.data.lng += 1;
                            (*this).value.data.lng;
                        }
                        sptr1 = sptr1.offset(1);
                        sptr1;
                    }
                }
                _ => {}
            }
            (*this).operation = -(1000 as c_int);
        } else {
            Allocate_Ptrs(lParse, this);
            if (*lParse).status == 0 {
                rows = (*lParse).nRows;
                match (*this).operation {
                    279..=284 => loop {
                        let fresh40 = rows;
                        rows -= 1;
                        if fresh40 == 0 {
                            break;
                        }
                        if const1 == 0 {
                            sptr1 = *((*that1).value.data.strptr).offset(rows as isize);
                        }
                        if const2 == 0 {
                            sptr2 = *((*that2).value.data.strptr).offset(rows as isize);
                        }
                        match (*this).operation {
                            280 => {
                                *((*this).value.data.logptr).offset(rows as isize) =
                                    if bitcmp(sptr1, sptr2) == 0 { 1 } else { 0 };
                            }
                            279 => {
                                *((*this).value.data.logptr).offset(rows as isize) =
                                    bitcmp(sptr1, sptr2);
                            }
                            281..=284 => {
                                *((*this).value.data.logptr).offset(rows as isize) =
                                    bitlgte(sptr1, (*this).operation, sptr2);
                            }
                            _ => {}
                        }
                        *((*this).value.undef).offset(rows as isize) = 0;
                    },
                    124 | 38 | 43 => loop {
                        let fresh41 = rows;
                        rows -= 1;
                        if fresh41 == 0 {
                            break;
                        }
                        if const1 == 0 {
                            sptr1 = *((*that1).value.data.strptr).offset(rows as isize);
                        }
                        if const2 == 0 {
                            sptr2 = *((*that2).value.data.strptr).offset(rows as isize);
                        }
                        if (*this).operation == '|' as i32 {
                            bitor(
                                *((*this).value.data.strptr).offset(rows as isize),
                                sptr1,
                                sptr2,
                            );
                        } else if (*this).operation == '&' as i32 {
                            bitand(
                                *((*this).value.data.strptr).offset(rows as isize),
                                sptr1,
                                sptr2,
                            );
                        } else {
                            strcpy(*((*this).value.data.strptr).offset(rows as isize), sptr1);
                            strcat(*((*this).value.data.strptr).offset(rows as isize), sptr2);
                        }
                    },
                    291 => {
                        let mut i: c_long = 0;
                        let mut previous: c_long = 0;
                        let mut curr: c_long = 0;
                        previous = (*that2).value.data.lng;
                        i = 0;
                        while i < rows {
                            sptr1 = *((*that1).value.data.strptr).offset(i as isize);
                            curr = 0;
                            while *sptr1 != 0 {
                                if *sptr1 as c_int == '1' as i32 {
                                    curr += 1;
                                    curr;
                                }
                                sptr1 = sptr1.offset(1);
                                sptr1;
                            }
                            previous += curr;
                            *((*this).value.data.lngptr).offset(i as isize) = previous;
                            *((*this).value.undef).offset(i as isize) = 0;
                            i += 1;
                            i;
                        }
                        (*that2).value.data.lng = previous;
                    }
                    _ => {}
                }
            }
        }
        if (*that1).operation > 0 {
            free(*((*that1).value.data.strptr).offset(0) as *mut c_void);
            free((*that1).value.data.strptr as *mut c_void);
        }
        if (*that2).operation > 0 {
            free(*((*that2).value.data.strptr).offset(0) as *mut c_void);
            free((*that2).value.data.strptr as *mut c_void);
        }
    }
}
unsafe fn Do_BinOp_str(mut lParse: *mut ParseData, mut this: *mut Node) {
    unsafe {
        let mut that1: *mut Node = std::ptr::null_mut::<Node>();
        let mut that2: *mut Node = std::ptr::null_mut::<Node>();
        let mut sptr1: *mut c_char = ptr::null_mut();
        let mut sptr2: *mut c_char = ptr::null_mut();
        let mut null1: c_char = 0;
        let mut null2: c_char = 0;
        let mut const1: c_int = 0;
        let mut const2: c_int = 0;
        let mut val: c_int = 0;
        let mut rows: c_long = 0;
        that1 = ((*lParse).Nodes).offset((*this).SubNodes[0] as isize);
        that2 = ((*lParse).Nodes).offset((*this).SubNodes[1] as isize);
        const1 = ((*that1).operation == -(1000 as c_int)) as c_int;
        const2 = ((*that2).operation == -(1000 as c_int)) as c_int;
        sptr1 = if const1 != 0 {
            ((*that1).value.data.astr).as_mut_ptr()
        } else {
            std::ptr::null_mut::<c_char>()
        };
        sptr2 = if const2 != 0 {
            ((*that2).value.data.astr).as_mut_ptr()
        } else {
            std::ptr::null_mut::<c_char>()
        };
        if const1 != 0 && const2 != 0 {
            match (*this).operation {
                280 | 279 => {
                    val = ((if (*sptr1.offset(0) as c_int) < *sptr2.offset(0) as c_int {
                        -(1)
                    } else if *sptr1.offset(0) as c_int > *sptr2.offset(0) as c_int {
                        1
                    } else {
                        strcmp(sptr1, sptr2)
                    }) == 0) as c_int;
                    (*this).value.data.log =
                        (if (*this).operation == fits_parser_yytokentype::EQ as c_int {
                            val
                        } else {
                            (val == 0) as c_int
                        }) as c_char;
                }
                281 => {
                    (*this).value.data.log =
                        if (if (*sptr1.offset(0) as c_int) < *sptr2.offset(0) as c_int {
                            -(1)
                        } else if *sptr1.offset(0) as c_int > *sptr2.offset(0) as c_int {
                            1
                        } else {
                            strcmp(sptr1, sptr2)
                        }) > 0
                        {
                            1
                        } else {
                            0
                        };
                }
                282 => {
                    (*this).value.data.log =
                        if (if (*sptr1.offset(0) as c_int) < *sptr2.offset(0) as c_int {
                            -(1)
                        } else if *sptr1.offset(0) as c_int > *sptr2.offset(0) as c_int {
                            1
                        } else {
                            strcmp(sptr1, sptr2)
                        }) < 0
                        {
                            1
                        } else {
                            0
                        };
                }
                284 => {
                    (*this).value.data.log =
                        if (if (*sptr1.offset(0) as c_int) < *sptr2.offset(0) as c_int {
                            -(1)
                        } else if *sptr1.offset(0) as c_int > *sptr2.offset(0) as c_int {
                            1
                        } else {
                            strcmp(sptr1, sptr2)
                        }) >= 0
                        {
                            1
                        } else {
                            0
                        };
                }
                283 => {
                    (*this).value.data.log =
                        if (if (*sptr1.offset(0) as c_int) < *sptr2.offset(0) as c_int {
                            -(1)
                        } else if *sptr1.offset(0) as c_int > *sptr2.offset(0) as c_int {
                            1
                        } else {
                            strcmp(sptr1, sptr2)
                        }) <= 0
                        {
                            1
                        } else {
                            0
                        };
                }
                43 => {
                    strcpy(((*this).value.data.astr).as_mut_ptr(), sptr1);
                    strcat(((*this).value.data.astr).as_mut_ptr(), sptr2);
                }
                _ => {}
            }
            (*this).operation = -(1000 as c_int);
        } else {
            Allocate_Ptrs(lParse, this);
            if (*lParse).status == 0 {
                rows = (*lParse).nRows;
                match (*this).operation {
                    280 | 279 => loop {
                        let fresh42 = rows;
                        rows -= 1;
                        if fresh42 == 0 {
                            break;
                        }
                        if const1 == 0 {
                            null1 = *((*that1).value.undef).offset(rows as isize);
                        }
                        if const2 == 0 {
                            null2 = *((*that2).value.undef).offset(rows as isize);
                        }
                        *((*this).value.undef).offset(rows as isize) =
                            if null1 as c_int != 0 || null2 as c_int != 0 {
                                1
                            } else {
                                0
                            };
                        if *((*this).value.undef).offset(rows as isize) == 0 {
                            if const1 == 0 {
                                sptr1 = *((*that1).value.data.strptr).offset(rows as isize);
                            }
                            if const2 == 0 {
                                sptr2 = *((*that2).value.data.strptr).offset(rows as isize);
                            }
                            val = ((if (*sptr1.offset(0) as c_int) < *sptr2.offset(0) as c_int {
                                -(1)
                            } else if *sptr1.offset(0) as c_int > *sptr2.offset(0) as c_int {
                                1
                            } else {
                                strcmp(sptr1, sptr2)
                            }) == 0) as c_int;
                            *((*this).value.data.logptr).offset(rows as isize) =
                                (if (*this).operation == fits_parser_yytokentype::EQ as c_int {
                                    val
                                } else {
                                    (val == 0) as c_int
                                }) as c_char;
                        }
                    },
                    281 | 282 => loop {
                        let fresh43 = rows;
                        rows -= 1;
                        if fresh43 == 0 {
                            break;
                        }
                        if const1 == 0 {
                            null1 = *((*that1).value.undef).offset(rows as isize);
                        }
                        if const2 == 0 {
                            null2 = *((*that2).value.undef).offset(rows as isize);
                        }
                        *((*this).value.undef).offset(rows as isize) =
                            if null1 as c_int != 0 || null2 as c_int != 0 {
                                1
                            } else {
                                0
                            };
                        if *((*this).value.undef).offset(rows as isize) == 0 {
                            if const1 == 0 {
                                sptr1 = *((*that1).value.data.strptr).offset(rows as isize);
                            }
                            if const2 == 0 {
                                sptr2 = *((*that2).value.data.strptr).offset(rows as isize);
                            }
                            val = if (*sptr1.offset(0) as c_int) < *sptr2.offset(0) as c_int {
                                -(1)
                            } else if *sptr1.offset(0) as c_int > *sptr2.offset(0) as c_int {
                                1
                            } else {
                                strcmp(sptr1, sptr2)
                            };
                            *((*this).value.data.logptr).offset(rows as isize) =
                                (if (*this).operation == fits_parser_yytokentype::GT as c_int {
                                    (val > 0) as c_int
                                } else {
                                    (val < 0) as c_int
                                }) as c_char;
                        }
                    },
                    284 | 283 => loop {
                        let fresh44 = rows;
                        rows -= 1;
                        if fresh44 == 0 {
                            break;
                        }
                        if const1 == 0 {
                            null1 = *((*that1).value.undef).offset(rows as isize);
                        }
                        if const2 == 0 {
                            null2 = *((*that2).value.undef).offset(rows as isize);
                        }
                        *((*this).value.undef).offset(rows as isize) =
                            if null1 as c_int != 0 || null2 as c_int != 0 {
                                1
                            } else {
                                0
                            };
                        if *((*this).value.undef).offset(rows as isize) == 0 {
                            if const1 == 0 {
                                sptr1 = *((*that1).value.data.strptr).offset(rows as isize);
                            }
                            if const2 == 0 {
                                sptr2 = *((*that2).value.data.strptr).offset(rows as isize);
                            }
                            val = if (*sptr1.offset(0) as c_int) < *sptr2.offset(0) as c_int {
                                -(1)
                            } else if *sptr1.offset(0) as c_int > *sptr2.offset(0) as c_int {
                                1
                            } else {
                                strcmp(sptr1, sptr2)
                            };
                            *((*this).value.data.logptr).offset(rows as isize) =
                                (if (*this).operation == fits_parser_yytokentype::GTE as c_int {
                                    (val >= 0) as c_int
                                } else {
                                    (val <= 0) as c_int
                                }) as c_char;
                        }
                    },
                    43 => loop {
                        let fresh45 = rows;
                        rows -= 1;
                        if fresh45 == 0 {
                            break;
                        }
                        if const1 == 0 {
                            null1 = *((*that1).value.undef).offset(rows as isize);
                        }
                        if const2 == 0 {
                            null2 = *((*that2).value.undef).offset(rows as isize);
                        }
                        *((*this).value.undef).offset(rows as isize) =
                            if null1 as c_int != 0 || null2 as c_int != 0 {
                                1
                            } else {
                                0
                            };
                        if *((*this).value.undef).offset(rows as isize) == 0 {
                            if const1 == 0 {
                                sptr1 = *((*that1).value.data.strptr).offset(rows as isize);
                            }
                            if const2 == 0 {
                                sptr2 = *((*that2).value.data.strptr).offset(rows as isize);
                            }
                            strcpy(*((*this).value.data.strptr).offset(rows as isize), sptr1);
                            strcat(*((*this).value.data.strptr).offset(rows as isize), sptr2);
                        }
                    },
                    _ => {}
                }
            }
        }
        if (*that1).operation > 0 {
            free(*((*that1).value.data.strptr).offset(0) as *mut c_void);
            free((*that1).value.data.strptr as *mut c_void);
        }
        if (*that2).operation > 0 {
            free(*((*that2).value.data.strptr).offset(0) as *mut c_void);
            free((*that2).value.data.strptr as *mut c_void);
        }
    }
}
unsafe fn Do_BinOp_log(mut lParse: *mut ParseData, mut this: *mut Node) {
    unsafe {
        let mut that1: *mut Node = std::ptr::null_mut::<Node>();
        let mut that2: *mut Node = std::ptr::null_mut::<Node>();
        let mut vector1: c_int = 0;
        let mut vector2: c_int = 0;
        let mut val1: c_char = 0;
        let mut val2: c_char = 0;
        let mut null1: c_char = 0;
        let mut null2: c_char = 0;
        let mut rows: c_long = 0;
        let mut nelem: c_long = 0;
        let mut elem: c_long = 0;
        that1 = ((*lParse).Nodes).offset((*this).SubNodes[0] as isize);
        that2 = ((*lParse).Nodes).offset((*this).SubNodes[1] as isize);
        vector1 = ((*that1).operation != -(1000 as c_int)) as c_int;
        if vector1 != 0 {
            vector1 = (*that1).value.nelem as c_int;
        } else {
            val1 = (*that1).value.data.log;
        }
        vector2 = ((*that2).operation != -(1000 as c_int)) as c_int;
        if vector2 != 0 {
            vector2 = (*that2).value.nelem as c_int;
        } else {
            val2 = (*that2).value.data.log;
        }
        if vector1 == 0 && vector2 == 0 {
            match (*this).operation {
                277 => {
                    (*this).value.data.log = if val1 as c_int != 0 || val2 as c_int != 0 {
                        1
                    } else {
                        0
                    };
                }
                278 => {
                    (*this).value.data.log = if val1 as c_int != 0 && val2 as c_int != 0 {
                        1
                    } else {
                        0
                    };
                }
                279 => {
                    (*this).value.data.log = (val1 as c_int != 0 && val2 as c_int != 0
                        || val1 == 0 && val2 == 0)
                        as c_int as c_char;
                }
                280 => {
                    (*this).value.data.log =
                        if val1 as c_int != 0 && val2 == 0 || val1 == 0 && val2 as c_int != 0 {
                            1
                        } else {
                            0
                        };
                }
                291 => {
                    (*this).value.data.lng = val1 as c_long;
                }
                _ => {}
            }
            (*this).operation = -(1000 as c_int);
        } else if (*this).operation == fits_parser_yytokentype::ACCUM as c_int {
            let mut i: c_long = 0;
            let mut previous: c_long = 0;
            let mut curr: c_long = 0;
            rows = (*lParse).nRows;
            nelem = (*this).value.nelem;
            elem = (*this).value.nelem * rows;
            Allocate_Ptrs(lParse, this);
            if (*lParse).status == 0 {
                previous = (*that2).value.data.lng;
                i = 0;
                while i < elem {
                    if *((*that1).value.undef).offset(i as isize) == 0 {
                        curr = *((*that1).value.data.logptr).offset(i as isize) as c_long;
                        previous += curr;
                    }
                    *((*this).value.data.lngptr).offset(i as isize) = previous;
                    *((*this).value.undef).offset(i as isize) = 0;
                    i += 1;
                    i;
                }
                (*that2).value.data.lng = previous;
            }
        } else {
            rows = (*lParse).nRows;
            nelem = (*this).value.nelem;
            elem = (*this).value.nelem * rows;
            Allocate_Ptrs(lParse, this);
            if (*lParse).status == 0 {
                if (*this).operation == fits_parser_yytokentype::ACCUM as c_int {
                    let mut i_0: c_long = 0;
                    let mut previous_0: c_long = 0;
                    let mut curr_0: c_long = 0;
                    previous_0 = (*that2).value.data.lng;
                    i_0 = 0;
                    while i_0 < elem {
                        if *((*that1).value.undef).offset(i_0 as isize) == 0 {
                            curr_0 = *((*that1).value.data.logptr).offset(i_0 as isize) as c_long;
                            previous_0 += curr_0;
                        }
                        *((*this).value.data.lngptr).offset(i_0 as isize) = previous_0;
                        *((*this).value.undef).offset(i_0 as isize) = 0;
                        i_0 += 1;
                        i_0;
                    }
                    (*that2).value.data.lng = previous_0;
                }
                loop {
                    let fresh46 = rows;
                    rows -= 1;
                    if fresh46 == 0 {
                        break;
                    }
                    loop {
                        let fresh47 = nelem;
                        nelem -= 1;
                        if fresh47 == 0 {
                            break;
                        }
                        elem -= 1;
                        elem;
                        if vector1 > 1 {
                            val1 = *((*that1).value.data.logptr).offset(elem as isize);
                            null1 = *((*that1).value.undef).offset(elem as isize);
                        } else if vector1 != 0 {
                            val1 = *((*that1).value.data.logptr).offset(rows as isize);
                            null1 = *((*that1).value.undef).offset(rows as isize);
                        }
                        if vector2 > 1 {
                            val2 = *((*that2).value.data.logptr).offset(elem as isize);
                            null2 = *((*that2).value.undef).offset(elem as isize);
                        } else if vector2 != 0 {
                            val2 = *((*that2).value.data.logptr).offset(rows as isize);
                            null2 = *((*that2).value.undef).offset(rows as isize);
                        }
                        *((*this).value.undef).offset(elem as isize) =
                            if null1 as c_int != 0 || null2 as c_int != 0 {
                                1
                            } else {
                                0
                            };
                        match (*this).operation {
                            277 => {
                                if null1 == 0 && null2 == 0 {
                                    *((*this).value.data.logptr).offset(elem as isize) =
                                        if val1 as c_int != 0 || val2 as c_int != 0 {
                                            1
                                        } else {
                                            0
                                        };
                                } else if null1 as c_int != 0 && null2 == 0 && val2 as c_int != 0
                                    || null1 == 0 && null2 as c_int != 0 && val1 as c_int != 0
                                {
                                    *((*this).value.data.logptr).offset(elem as isize) = 1;
                                    *((*this).value.undef).offset(elem as isize) = 0;
                                }
                            }
                            278 => {
                                if null1 == 0 && null2 == 0 {
                                    *((*this).value.data.logptr).offset(elem as isize) =
                                        if val1 as c_int != 0 && val2 as c_int != 0 {
                                            1
                                        } else {
                                            0
                                        };
                                } else if null1 as c_int != 0 && null2 == 0 && val2 == 0
                                    || null1 == 0 && null2 as c_int != 0 && val1 == 0
                                {
                                    *((*this).value.data.logptr).offset(elem as isize) = 0;
                                    *((*this).value.undef).offset(elem as isize) = 0;
                                }
                            }
                            279 => {
                                *((*this).value.data.logptr).offset(elem as isize) =
                                    (val1 as c_int != 0 && val2 as c_int != 0
                                        || val1 == 0 && val2 == 0)
                                        as c_int as c_char;
                            }
                            280 => {
                                *((*this).value.data.logptr).offset(elem as isize) =
                                    (val1 as c_int != 0 && val2 == 0
                                        || val1 == 0 && val2 as c_int != 0)
                                        as c_int as c_char;
                            }
                            _ => {}
                        }
                    }
                    nelem = (*this).value.nelem;
                }
            }
        }
        if (*that1).operation > 0 {
            free((*that1).value.data.ptr);
        }
        if (*that2).operation > 0 {
            free((*that2).value.data.ptr);
        }
    }
}
unsafe fn Do_BinOp_lng(mut lParse: *mut ParseData, mut this: *mut Node) {
    unsafe {
        let mut that1: *mut Node = std::ptr::null_mut::<Node>();
        let mut that2: *mut Node = std::ptr::null_mut::<Node>();
        let mut vector1: c_int = 0;
        let mut vector2: c_int = 0;
        let mut val1: c_long = 0;
        let mut val2: c_long = 0;
        let mut null1: c_char = 0;
        let mut null2: c_char = 0;
        let mut rows: c_long = 0;
        let mut nelem: c_long = 0;
        let mut elem: c_long = 0;
        that1 = ((*lParse).Nodes).offset((*this).SubNodes[0] as isize);
        that2 = ((*lParse).Nodes).offset((*this).SubNodes[1] as isize);
        vector1 = ((*that1).operation != -(1000 as c_int)) as c_int;
        if vector1 != 0 {
            vector1 = (*that1).value.nelem as c_int;
        } else {
            val1 = (*that1).value.data.lng;
        }
        vector2 = ((*that2).operation != -(1000 as c_int)) as c_int;
        if vector2 != 0 {
            vector2 = (*that2).value.nelem as c_int;
        } else {
            val2 = (*that2).value.data.lng;
        }
        if vector1 == 0 && vector2 == 0 {
            match (*this).operation {
                126 | 279 => {
                    (*this).value.data.log = if val1 == val2 { 1 } else { 0 };
                }
                280 => {
                    (*this).value.data.log = if val1 != val2 { 1 } else { 0 };
                }
                281 => {
                    (*this).value.data.log = if val1 > val2 { 1 } else { 0 };
                }
                282 => {
                    (*this).value.data.log = if val1 < val2 { 1 } else { 0 };
                }
                283 => {
                    (*this).value.data.log = if val1 <= val2 { 1 } else { 0 };
                }
                284 => {
                    (*this).value.data.log = if val1 >= val2 { 1 } else { 0 };
                }
                43 => {
                    (*this).value.data.lng = val1 + val2;
                }
                45 => {
                    (*this).value.data.lng = val1 - val2;
                }
                42 => {
                    (*this).value.data.lng = val1 * val2;
                }
                38 => {
                    (*this).value.data.lng = val1 & val2;
                }
                124 => {
                    (*this).value.data.lng = val1 | val2;
                }
                94 => {
                    (*this).value.data.lng = val1 ^ val2;
                }
                37 => {
                    if val2 != 0 {
                        (*this).value.data.lng = val1 % val2;
                    } else {
                        fits_parser_yyerror(
                            0 as yyscan_t,
                            lParse,
                            b"Divide by Zero\0" as *const u8 as *const c_char as *mut c_char,
                        );
                    }
                }
                47 => {
                    if val2 != 0 {
                        (*this).value.data.lng = val1 / val2;
                    } else {
                        fits_parser_yyerror(
                            0 as yyscan_t,
                            lParse,
                            b"Divide by Zero\0" as *const u8 as *const c_char as *mut c_char,
                        );
                    }
                }
                286 => {
                    (*this).value.data.lng = pow(val1 as c_double, val2 as c_double) as c_long;
                }
                291 => {
                    (*this).value.data.lng = val1;
                }
                292 => {
                    (*this).value.data.lng = 0;
                }
                _ => {}
            }
            (*this).operation = -(1000 as c_int);
        } else if (*this).operation == fits_parser_yytokentype::ACCUM as c_int
            || (*this).operation == fits_parser_yytokentype::DIFF as c_int
        {
            let mut i: c_long = 0;
            let mut previous: c_long = 0;
            let mut curr: c_long = 0;
            let mut undef: c_long = 0;
            rows = (*lParse).nRows;
            nelem = (*this).value.nelem;
            elem = (*this).value.nelem * rows;
            Allocate_Ptrs(lParse, this);
            if (*lParse).status == 0 {
                previous = (*that2).value.data.lng;
                undef = *(*that2).value.undef as c_long;
                if (*this).operation == fits_parser_yytokentype::ACCUM as c_int {
                    i = 0;
                    while i < elem {
                        if *((*that1).value.undef).offset(i as isize) == 0 {
                            curr = *((*that1).value.data.lngptr).offset(i as isize);
                            previous += curr;
                        }
                        *((*this).value.data.lngptr).offset(i as isize) = previous;
                        *((*this).value.undef).offset(i as isize) = 0;
                        i += 1;
                        i;
                    }
                } else {
                    i = 0;
                    while i < elem {
                        curr = *((*that1).value.data.lngptr).offset(i as isize);
                        if *((*that1).value.undef).offset(i as isize) as c_int != 0 || undef != 0 {
                            *((*this).value.data.lngptr).offset(i as isize) = 0;
                            *((*this).value.undef).offset(i as isize) = 1;
                        } else {
                            *((*this).value.data.lngptr).offset(i as isize) = curr - previous;
                            *((*this).value.undef).offset(i as isize) = 0;
                        }
                        previous = curr;
                        undef = *((*that1).value.undef).offset(i as isize) as c_long;
                        i += 1;
                        i;
                    }
                }
                (*that2).value.data.lng = previous;
                (*that2).value.undef = undef as *mut c_char;
            }
        } else {
            rows = (*lParse).nRows;
            nelem = (*this).value.nelem;
            elem = (*this).value.nelem * rows;
            Allocate_Ptrs(lParse, this);
            loop {
                let fresh48 = rows;
                rows -= 1;
                if !(fresh48 != 0 && (*lParse).status == 0) {
                    break;
                }
                loop {
                    let fresh49 = nelem;
                    nelem -= 1;
                    if !(fresh49 != 0 && (*lParse).status == 0) {
                        break;
                    }
                    elem -= 1;
                    elem;
                    if vector1 > 1 {
                        val1 = *((*that1).value.data.lngptr).offset(elem as isize);
                        null1 = *((*that1).value.undef).offset(elem as isize);
                    } else if vector1 != 0 {
                        val1 = *((*that1).value.data.lngptr).offset(rows as isize);
                        null1 = *((*that1).value.undef).offset(rows as isize);
                    }
                    if vector2 > 1 {
                        val2 = *((*that2).value.data.lngptr).offset(elem as isize);
                        null2 = *((*that2).value.undef).offset(elem as isize);
                    } else if vector2 != 0 {
                        val2 = *((*that2).value.data.lngptr).offset(rows as isize);
                        null2 = *((*that2).value.undef).offset(rows as isize);
                    }
                    *((*this).value.undef).offset(elem as isize) =
                        if null1 as c_int != 0 || null2 as c_int != 0 {
                            1
                        } else {
                            0
                        };
                    match (*this).operation {
                        126 | 279 => {
                            *((*this).value.data.logptr).offset(elem as isize) =
                                if val1 == val2 { 1 } else { 0 };
                        }
                        280 => {
                            *((*this).value.data.logptr).offset(elem as isize) =
                                if val1 != val2 { 1 } else { 0 };
                        }
                        281 => {
                            *((*this).value.data.logptr).offset(elem as isize) =
                                if val1 > val2 { 1 } else { 0 };
                        }
                        282 => {
                            *((*this).value.data.logptr).offset(elem as isize) =
                                if val1 < val2 { 1 } else { 0 };
                        }
                        283 => {
                            *((*this).value.data.logptr).offset(elem as isize) =
                                if val1 <= val2 { 1 } else { 0 };
                        }
                        284 => {
                            *((*this).value.data.logptr).offset(elem as isize) =
                                if val1 >= val2 { 1 } else { 0 };
                        }
                        43 => {
                            *((*this).value.data.lngptr).offset(elem as isize) = val1 + val2;
                        }
                        45 => {
                            *((*this).value.data.lngptr).offset(elem as isize) = val1 - val2;
                        }
                        42 => {
                            *((*this).value.data.lngptr).offset(elem as isize) = val1 * val2;
                        }
                        38 => {
                            *((*this).value.data.lngptr).offset(elem as isize) = val1 & val2;
                        }
                        124 => {
                            *((*this).value.data.lngptr).offset(elem as isize) = val1 | val2;
                        }
                        94 => {
                            *((*this).value.data.lngptr).offset(elem as isize) = val1 ^ val2;
                        }
                        37 => {
                            if val2 != 0 {
                                *((*this).value.data.lngptr).offset(elem as isize) = val1 % val2;
                            } else {
                                *((*this).value.data.lngptr).offset(elem as isize) = 0;
                                *((*this).value.undef).offset(elem as isize) = 1;
                            }
                        }
                        47 => {
                            if val2 != 0 {
                                *((*this).value.data.lngptr).offset(elem as isize) = val1 / val2;
                            } else {
                                *((*this).value.data.lngptr).offset(elem as isize) = 0;
                                *((*this).value.undef).offset(elem as isize) = 1;
                            }
                        }
                        286 => {
                            *((*this).value.data.lngptr).offset(elem as isize) =
                                pow(val1 as c_double, val2 as c_double) as c_long;
                        }
                        _ => {}
                    }
                }
                nelem = (*this).value.nelem;
            }
        }
        if (*that1).operation > 0 {
            free((*that1).value.data.ptr);
        }
        if (*that2).operation > 0 {
            free((*that2).value.data.ptr);
        }
    }
}
unsafe fn Do_BinOp_dbl(mut lParse: *mut ParseData, mut this: *mut Node) {
    unsafe {
        let mut that1: *mut Node = std::ptr::null_mut::<Node>();
        let mut that2: *mut Node = std::ptr::null_mut::<Node>();
        let mut vector1: c_int = 0;
        let mut vector2: c_int = 0;
        let mut val1: c_double = 0.0f64;
        let mut val2: c_double = 0.0f64;
        let mut null1: c_char = 0;
        let mut null2: c_char = 0;
        let mut rows: c_long = 0;
        let mut nelem: c_long = 0;
        let mut elem: c_long = 0;
        that1 = ((*lParse).Nodes).offset((*this).SubNodes[0] as isize);
        that2 = ((*lParse).Nodes).offset((*this).SubNodes[1] as isize);
        vector1 = ((*that1).operation != -(1000 as c_int)) as c_int;
        if vector1 != 0 {
            vector1 = (*that1).value.nelem as c_int;
        } else {
            val1 = (*that1).value.data.dbl;
        }
        vector2 = ((*that2).operation != -(1000 as c_int)) as c_int;
        if vector2 != 0 {
            vector2 = (*that2).value.nelem as c_int;
        } else {
            val2 = (*that2).value.data.dbl;
        }
        if vector1 == 0 && vector2 == 0 {
            match (*this).operation {
                126 => {
                    (*this).value.data.log = if fabs(val1 - val2) < 1.0e-7f64 { 1 } else { 0 };
                }
                279 => {
                    (*this).value.data.log = if val1 == val2 { 1 } else { 0 };
                }
                280 => {
                    (*this).value.data.log = if val1 != val2 { 1 } else { 0 };
                }
                281 => {
                    (*this).value.data.log = if val1 > val2 { 1 } else { 0 };
                }
                282 => {
                    (*this).value.data.log = if val1 < val2 { 1 } else { 0 };
                }
                283 => {
                    (*this).value.data.log = if val1 <= val2 { 1 } else { 0 };
                }
                284 => {
                    (*this).value.data.log = if val1 >= val2 { 1 } else { 0 };
                }
                43 => {
                    (*this).value.data.dbl = val1 + val2;
                }
                45 => {
                    (*this).value.data.dbl = val1 - val2;
                }
                42 => {
                    (*this).value.data.dbl = val1 * val2;
                }
                37 => {
                    if val2 != 0. {
                        (*this).value.data.dbl = val1 - val2 * (val1 / val2) as c_int as c_double;
                    } else {
                        fits_parser_yyerror(
                            0 as yyscan_t,
                            lParse,
                            b"Divide by Zero\0" as *const u8 as *const c_char as *mut c_char,
                        );
                    }
                }
                47 => {
                    if val2 != 0. {
                        (*this).value.data.dbl = val1 / val2;
                    } else {
                        fits_parser_yyerror(
                            0 as yyscan_t,
                            lParse,
                            b"Divide by Zero\0" as *const u8 as *const c_char as *mut c_char,
                        );
                    }
                }
                286 => {
                    (*this).value.data.dbl = pow(val1, val2);
                }
                291 => {
                    (*this).value.data.dbl = val1;
                }
                292 => {
                    (*this).value.data.dbl = 0.0;
                }
                _ => {}
            }
            (*this).operation = -(1000 as c_int);
        } else if (*this).operation == fits_parser_yytokentype::ACCUM as c_int
            || (*this).operation == fits_parser_yytokentype::DIFF as c_int
        {
            let mut i: c_long = 0;
            let mut undef: c_long = 0;
            let mut previous: c_double = 0.;
            let mut curr: c_double = 0.;
            rows = (*lParse).nRows;
            nelem = (*this).value.nelem;
            elem = (*this).value.nelem * rows;
            Allocate_Ptrs(lParse, this);
            if (*lParse).status == 0 {
                previous = (*that2).value.data.dbl;
                undef = *(*that2).value.undef as c_long;
                if (*this).operation == fits_parser_yytokentype::ACCUM as c_int {
                    i = 0;
                    while i < elem {
                        if *((*that1).value.undef).offset(i as isize) == 0 {
                            curr = *((*that1).value.data.dblptr).offset(i as isize);
                            previous += curr;
                        }
                        *((*this).value.data.dblptr).offset(i as isize) = previous;
                        *((*this).value.undef).offset(i as isize) = 0;
                        i += 1;
                        i;
                    }
                } else {
                    i = 0;
                    while i < elem {
                        curr = *((*that1).value.data.dblptr).offset(i as isize);
                        if *((*that1).value.undef).offset(i as isize) as c_int != 0 || undef != 0 {
                            *((*this).value.data.dblptr).offset(i as isize) = 0.0;
                            *((*this).value.undef).offset(i as isize) = 1;
                        } else {
                            *((*this).value.data.dblptr).offset(i as isize) = curr - previous;
                            *((*this).value.undef).offset(i as isize) = 0;
                        }
                        previous = curr;
                        undef = *((*that1).value.undef).offset(i as isize) as c_long;
                        i += 1;
                        i;
                    }
                }
                (*that2).value.data.dbl = previous;
                (*that2).value.undef = undef as *mut c_char;
            }
        } else {
            rows = (*lParse).nRows;
            nelem = (*this).value.nelem;
            elem = (*this).value.nelem * rows;
            Allocate_Ptrs(lParse, this);
            loop {
                let fresh50 = rows;
                rows -= 1;
                if !(fresh50 != 0 && (*lParse).status == 0) {
                    break;
                }
                loop {
                    let fresh51 = nelem;
                    nelem -= 1;
                    if !(fresh51 != 0 && (*lParse).status == 0) {
                        break;
                    }
                    elem -= 1;
                    elem;
                    if vector1 > 1 {
                        val1 = *((*that1).value.data.dblptr).offset(elem as isize);
                        null1 = *((*that1).value.undef).offset(elem as isize);
                    } else if vector1 != 0 {
                        val1 = *((*that1).value.data.dblptr).offset(rows as isize);
                        null1 = *((*that1).value.undef).offset(rows as isize);
                    }
                    if vector2 > 1 {
                        val2 = *((*that2).value.data.dblptr).offset(elem as isize);
                        null2 = *((*that2).value.undef).offset(elem as isize);
                    } else if vector2 != 0 {
                        val2 = *((*that2).value.data.dblptr).offset(rows as isize);
                        null2 = *((*that2).value.undef).offset(rows as isize);
                    }
                    *((*this).value.undef).offset(elem as isize) =
                        if null1 as c_int != 0 || null2 as c_int != 0 {
                            1
                        } else {
                            0
                        };
                    match (*this).operation {
                        126 => {
                            *((*this).value.data.logptr).offset(elem as isize) =
                                if fabs(val1 - val2) < 1.0e-7f64 { 1 } else { 0 };
                        }
                        279 => {
                            *((*this).value.data.logptr).offset(elem as isize) =
                                if val1 == val2 { 1 } else { 0 };
                        }
                        280 => {
                            *((*this).value.data.logptr).offset(elem as isize) =
                                if val1 != val2 { 1 } else { 0 };
                        }
                        281 => {
                            *((*this).value.data.logptr).offset(elem as isize) =
                                if val1 > val2 { 1 } else { 0 };
                        }
                        282 => {
                            *((*this).value.data.logptr).offset(elem as isize) =
                                if val1 < val2 { 1 } else { 0 };
                        }
                        283 => {
                            *((*this).value.data.logptr).offset(elem as isize) =
                                if val1 <= val2 { 1 } else { 0 };
                        }
                        284 => {
                            *((*this).value.data.logptr).offset(elem as isize) =
                                if val1 >= val2 { 1 } else { 0 };
                        }
                        43 => {
                            *((*this).value.data.dblptr).offset(elem as isize) = val1 + val2;
                        }
                        45 => {
                            *((*this).value.data.dblptr).offset(elem as isize) = val1 - val2;
                        }
                        42 => {
                            *((*this).value.data.dblptr).offset(elem as isize) = val1 * val2;
                        }
                        37 => {
                            if val2 != 0. {
                                *((*this).value.data.dblptr).offset(elem as isize) =
                                    val1 - val2 * (val1 / val2) as c_int as c_double;
                            } else {
                                *((*this).value.data.dblptr).offset(elem as isize) = 0.0f64;
                                *((*this).value.undef).offset(elem as isize) = 1;
                            }
                        }
                        47 => {
                            if val2 != 0. {
                                *((*this).value.data.dblptr).offset(elem as isize) = val1 / val2;
                            } else {
                                *((*this).value.data.dblptr).offset(elem as isize) = 0.0f64;
                                *((*this).value.undef).offset(elem as isize) = 1;
                            }
                        }
                        286 => {
                            *((*this).value.data.dblptr).offset(elem as isize) = pow(val1, val2);
                        }
                        _ => {}
                    }
                }
                nelem = (*this).value.nelem;
            }
        }
        if (*that1).operation > 0 {
            free((*that1).value.data.ptr);
        }
        if (*that2).operation > 0 {
            free((*that2).value.data.ptr);
        }
    }
}

pub unsafe fn qselect_median_lng(mut arr: *mut c_long, mut n: c_int) -> c_long {
    unsafe {
        let mut low: c_int = 0;
        let mut high: c_int = 0;
        let mut median: c_int = 0;
        let mut middle: c_int = 0;
        let mut ll: c_int = 0;
        let mut hh: c_int = 0;
        low = 0;
        high = n - 1;
        median = (low + high) / 2 as c_int;
        loop {
            if high <= low {
                return *arr.offset(median as isize);
            }
            if high == low + 1 {
                if *arr.offset(low as isize) > *arr.offset(high as isize) {
                    let mut t: c_long = *arr.offset(low as isize);
                    *arr.offset(low as isize) = *arr.offset(high as isize);
                    *arr.offset(high as isize) = t;
                }
                return *arr.offset(median as isize);
            }
            middle = (low + high) / 2 as c_int;
            if *arr.offset(middle as isize) > *arr.offset(high as isize) {
                let mut t_0: c_long = *arr.offset(middle as isize);
                *arr.offset(middle as isize) = *arr.offset(high as isize);
                *arr.offset(high as isize) = t_0;
            }
            if *arr.offset(low as isize) > *arr.offset(high as isize) {
                let mut t_1: c_long = *arr.offset(low as isize);
                *arr.offset(low as isize) = *arr.offset(high as isize);
                *arr.offset(high as isize) = t_1;
            }
            if *arr.offset(middle as isize) > *arr.offset(low as isize) {
                let mut t_2: c_long = *arr.offset(middle as isize);
                *arr.offset(middle as isize) = *arr.offset(low as isize);
                *arr.offset(low as isize) = t_2;
            }
            let mut t_3: c_long = *arr.offset(middle as isize);
            *arr.offset(middle as isize) = *arr.offset((low + 1) as isize);
            *arr.offset((low + 1) as isize) = t_3;
            ll = low + 1;
            hh = high;
            loop {
                loop {
                    ll += 1;
                    ll;
                    if *arr.offset(low as isize) <= *arr.offset(ll as isize) {
                        break;
                    }
                }
                loop {
                    hh -= 1;
                    hh;
                    if *arr.offset(hh as isize) <= *arr.offset(low as isize) {
                        break;
                    }
                }
                if hh < ll {
                    break;
                }
                let mut t_4: c_long = *arr.offset(ll as isize);
                *arr.offset(ll as isize) = *arr.offset(hh as isize);
                *arr.offset(hh as isize) = t_4;
            }
            let mut t_5: c_long = *arr.offset(low as isize);
            *arr.offset(low as isize) = *arr.offset(hh as isize);
            *arr.offset(hh as isize) = t_5;
            if hh <= median {
                low = ll;
            }
            if hh >= median {
                high = hh - 1;
            }
        }
    }
}

pub unsafe fn qselect_median_dbl(mut arr: *mut c_double, mut n: c_int) -> c_double {
    unsafe {
        let mut low: c_int = 0;
        let mut high: c_int = 0;
        let mut median: c_int = 0;
        let mut middle: c_int = 0;
        let mut ll: c_int = 0;
        let mut hh: c_int = 0;
        low = 0;
        high = n - 1;
        median = (low + high) / 2 as c_int;
        loop {
            if high <= low {
                return *arr.offset(median as isize);
            }
            if high == low + 1 {
                if *arr.offset(low as isize) > *arr.offset(high as isize) {
                    let mut t: c_double = *arr.offset(low as isize);
                    *arr.offset(low as isize) = *arr.offset(high as isize);
                    *arr.offset(high as isize) = t;
                }
                return *arr.offset(median as isize);
            }
            middle = (low + high) / 2 as c_int;
            if *arr.offset(middle as isize) > *arr.offset(high as isize) {
                let mut t_0: c_double = *arr.offset(middle as isize);
                *arr.offset(middle as isize) = *arr.offset(high as isize);
                *arr.offset(high as isize) = t_0;
            }
            if *arr.offset(low as isize) > *arr.offset(high as isize) {
                let mut t_1: c_double = *arr.offset(low as isize);
                *arr.offset(low as isize) = *arr.offset(high as isize);
                *arr.offset(high as isize) = t_1;
            }
            if *arr.offset(middle as isize) > *arr.offset(low as isize) {
                let mut t_2: c_double = *arr.offset(middle as isize);
                *arr.offset(middle as isize) = *arr.offset(low as isize);
                *arr.offset(low as isize) = t_2;
            }
            let mut t_3: c_double = *arr.offset(middle as isize);
            *arr.offset(middle as isize) = *arr.offset((low + 1) as isize);
            *arr.offset((low + 1) as isize) = t_3;
            ll = low + 1;
            hh = high;
            loop {
                loop {
                    ll += 1;
                    ll;
                    if !(*arr.offset(low as isize) > *arr.offset(ll as isize)) {
                        break;
                    }
                }
                loop {
                    hh -= 1;
                    hh;
                    if !(*arr.offset(hh as isize) > *arr.offset(low as isize)) {
                        break;
                    }
                }
                if hh < ll {
                    break;
                }
                let mut t_4: c_double = *arr.offset(ll as isize);
                *arr.offset(ll as isize) = *arr.offset(hh as isize);
                *arr.offset(hh as isize) = t_4;
            }
            let mut t_5: c_double = *arr.offset(low as isize);
            *arr.offset(low as isize) = *arr.offset(hh as isize);
            *arr.offset(hh as isize) = t_5;
            if hh <= median {
                low = ll;
            }
            if hh >= median {
                high = hh - 1;
            }
        }
    }
}

pub unsafe fn angsep_calc(
    mut ra1: c_double,
    mut dec1: c_double,
    mut ra2: c_double,
    mut dec2: c_double,
) -> c_double {
    unsafe {
        static mut deg: c_double = 0.0;
        let mut a: c_double = 0.;
        let mut sdec: c_double = 0.;
        let mut sra: c_double = 0.;
        if deg == 0.0 {
            deg = 4.0 * atan(1.0) / 180.0;
        }
        sra = ((ra2 - ra1) * deg / 2.0).sin();
        sdec = ((dec2 - dec1) * deg / 2.0).sin();
        a = sdec * sdec + (dec1 * deg).cos() * cos(dec2 * deg) * sra * sra;
        if a < 0.0 {
            a = 0.0;
        }
        if a > 1 as c_double {
            a = 1 as c_double;
        }
        2.0f64 * atan2((a).sqrt(), (1.0f64 - a).sqrt()) / deg
    }
}
unsafe fn Do_Func(mut lParse: *mut ParseData, mut this: *mut Node) {
    unsafe {
        let mut theParams: [*mut Node; 10] = [std::ptr::null_mut::<Node>(); 10];
        let mut vector: [c_int; 10] = [0; 10];
        let mut allConst: c_int = 0;
        let mut pVals: [lval; 10] = [lval {
            nelem: 0,
            naxis: 0,
            naxes: [0; 5],
            undef: std::ptr::null_mut::<c_char>(),
            data: data_union { dbl: 0. },
        }; 10];
        let mut pNull: [c_char; 10] = [0; 10];
        let mut ival: c_long = 0;
        let mut dval: c_double = 0.;
        let mut i: c_int = 0;
        let mut valInit: c_int = 0;
        let mut row: c_long = 0;
        let mut elem: c_long = 0;
        let mut nelem: c_long = 0;
        i = (*this).nSubNodes;
        allConst = 1;
        loop {
            let fresh52 = i;
            i -= 1;
            if fresh52 == 0 {
                break;
            }
            theParams[i as usize] = ((*lParse).Nodes).offset((*this).SubNodes[i as usize] as isize);
            vector[i as usize] = ((*theParams[i as usize]).operation != -(1000 as c_int)) as c_int;
            if vector[i as usize] != 0 {
                allConst = 0;
                vector[i as usize] = (*theParams[i as usize]).value.nelem as c_int;
            } else {
                if (*theParams[i as usize]).ntype == fits_parser_yytokentype::DOUBLE as c_int {
                    pVals[i as usize].data.dbl = (*theParams[i as usize]).value.data.dbl;
                } else if (*theParams[i as usize]).ntype == fits_parser_yytokentype::LONG as c_int {
                    pVals[i as usize].data.lng = (*theParams[i as usize]).value.data.lng;
                } else if (*theParams[i as usize]).ntype
                    == fits_parser_yytokentype::BOOLEAN as c_int
                {
                    pVals[i as usize].data.log = (*theParams[i as usize]).value.data.log;
                } else {
                    strcpy(
                        (pVals[i as usize].data.astr).as_mut_ptr(),
                        ((*theParams[i as usize]).value.data.astr).as_mut_ptr(),
                    );
                }
                pNull[i as usize] = 0;
            }
        }
        if (*this).nSubNodes == 0 {
            allConst = 0;
        }
        if (*this).operation == poirnd_fct as c_int {
            allConst = 0;
        }
        if (*this).operation == gasrnd_fct as c_int {
            allConst = 0;
        }
        if (*this).operation == rnd_fct as c_int {
            allConst = 0;
        }
        if allConst != 0 {
            let mut current_block_139: u64;
            match (*this).operation {
                1002 => {
                    if (*theParams[0]).ntype == fits_parser_yytokentype::BOOLEAN as c_int {
                        (*this).value.data.lng = (if pVals[0].data.log as c_int != 0 {
                            1
                        } else {
                            0
                        }) as c_long;
                    } else if (*theParams[0]).ntype == fits_parser_yytokentype::LONG as c_int {
                        (*this).value.data.lng = pVals[0].data.lng;
                    } else if (*theParams[0]).ntype == fits_parser_yytokentype::DOUBLE as c_int {
                        (*this).value.data.dbl = pVals[0].data.dbl;
                    } else if (*theParams[0]).ntype == fits_parser_yytokentype::BITSTR as c_int {
                        strcpy(
                            ((*this).value.data.astr).as_mut_ptr(),
                            (pVals[0].data.astr).as_mut_ptr(),
                        );
                    }
                    current_block_139 = 7627602990488000394;
                }
                1038 => {
                    if (*theParams[0]).ntype == fits_parser_yytokentype::LONG as c_int {
                        (*this).value.data.dbl = pVals[0].data.lng as c_double;
                    } else if (*theParams[0]).ntype == fits_parser_yytokentype::DOUBLE as c_int {
                        (*this).value.data.dbl = pVals[0].data.dbl;
                    }
                    current_block_139 = 7627602990488000394;
                }
                1039 => {
                    (*this).value.data.dbl = 0.0;
                    current_block_139 = 7627602990488000394;
                }
                1037 => {
                    if (*theParams[0]).ntype == fits_parser_yytokentype::BOOLEAN as c_int {
                        (*this).value.data.lng = (if pVals[0].data.log as c_int != 0 {
                            1
                        } else {
                            0
                        }) as c_long;
                    } else if (*theParams[0]).ntype == fits_parser_yytokentype::LONG as c_int {
                        (*this).value.data.lng = pVals[0].data.lng;
                    } else {
                        (*this).value.data.dbl = pVals[0].data.dbl;
                    }
                    current_block_139 = 7627602990488000394;
                }
                1043 => {
                    if (*theParams[0]).ntype == fits_parser_yytokentype::DOUBLE as c_int {
                        (*this).value.data.lng = simplerng_getpoisson(pVals[0].data.dbl) as c_long;
                    } else {
                        (*this).value.data.lng =
                            simplerng_getpoisson(pVals[0].data.lng as c_double) as c_long;
                    }
                    current_block_139 = 7627602990488000394;
                }
                1017 => {
                    if (*theParams[0]).ntype == fits_parser_yytokentype::DOUBLE as c_int {
                        dval = pVals[0].data.dbl;
                        (*this).value.data.dbl = if dval > 0.0f64 { dval } else { -dval };
                    } else {
                        ival = pVals[0].data.lng;
                        (*this).value.data.lng = if ival > 0 { ival } else { -ival };
                    }
                    current_block_139 = 7627602990488000394;
                }
                1040 => {
                    (*this).value.data.lng = 1;
                    current_block_139 = 7627602990488000394;
                }
                1030 => {
                    (*this).value.data.log = 0;
                    current_block_139 = 7627602990488000394;
                }
                1031 => {
                    if (*this).ntype == fits_parser_yytokentype::BOOLEAN as c_int {
                        (*this).value.data.log = pVals[0].data.log;
                    } else if (*this).ntype == fits_parser_yytokentype::LONG as c_int {
                        (*this).value.data.lng = pVals[0].data.lng;
                    } else if (*this).ntype == fits_parser_yytokentype::DOUBLE as c_int {
                        (*this).value.data.dbl = pVals[0].data.dbl;
                    } else if (*this).ntype == fits_parser_yytokentype::STRING as c_int {
                        strcpy(
                            ((*this).value.data.astr).as_mut_ptr(),
                            (pVals[0].data.astr).as_mut_ptr(),
                        );
                    }
                    current_block_139 = 7627602990488000394;
                }
                1046 => {
                    if (*this).ntype == fits_parser_yytokentype::LONG as c_int {
                        (*this).value.data.lng = pVals[0].data.lng;
                    } else if (*this).ntype == fits_parser_yytokentype::DOUBLE as c_int {
                        (*this).value.data.dbl = pVals[0].data.dbl;
                    }
                    current_block_139 = 7627602990488000394;
                }
                1004 => {
                    (*this).value.data.dbl = sin(pVals[0].data.dbl);
                    current_block_139 = 7627602990488000394;
                }
                1005 => {
                    (*this).value.data.dbl = cos(pVals[0].data.dbl);
                    current_block_139 = 7627602990488000394;
                }
                1006 => {
                    (*this).value.data.dbl = tan(pVals[0].data.dbl);
                    current_block_139 = 7627602990488000394;
                }
                1007 => {
                    dval = pVals[0].data.dbl;
                    if dval < -1.0f64 || dval > 1.0f64 {
                        fits_parser_yyerror(
                            0 as yyscan_t,
                            lParse,
                            b"Out of range argument to arcsin\0" as *const u8 as *const c_char
                                as *mut c_char,
                        );
                    } else {
                        (*this).value.data.dbl = asin(dval);
                    }
                    current_block_139 = 7627602990488000394;
                }
                1008 => {
                    dval = pVals[0].data.dbl;
                    if dval < -1.0f64 || dval > 1.0f64 {
                        fits_parser_yyerror(
                            0 as yyscan_t,
                            lParse,
                            b"Out of range argument to arccos\0" as *const u8 as *const c_char
                                as *mut c_char,
                        );
                    } else {
                        (*this).value.data.dbl = acos(dval);
                    }
                    current_block_139 = 7627602990488000394;
                }
                1009 => {
                    (*this).value.data.dbl = atan(pVals[0].data.dbl);
                    current_block_139 = 7627602990488000394;
                }
                1010 => {
                    (*this).value.data.dbl = sinh(pVals[0].data.dbl);
                    current_block_139 = 7627602990488000394;
                }
                1011 => {
                    (*this).value.data.dbl = cosh(pVals[0].data.dbl);
                    current_block_139 = 7627602990488000394;
                }
                1012 => {
                    (*this).value.data.dbl = tanh(pVals[0].data.dbl);
                    current_block_139 = 7627602990488000394;
                }
                1013 => {
                    (*this).value.data.dbl = exp(pVals[0].data.dbl);
                    current_block_139 = 7627602990488000394;
                }
                1014 => {
                    dval = pVals[0].data.dbl;
                    if dval <= 0.0f64 {
                        fits_parser_yyerror(
                            0 as yyscan_t,
                            lParse,
                            b"Out of range argument to log\0" as *const u8 as *const c_char
                                as *mut c_char,
                        );
                    } else {
                        (*this).value.data.dbl = log(dval);
                    }
                    current_block_139 = 7627602990488000394;
                }
                1015 => {
                    dval = pVals[0].data.dbl;
                    if dval <= 0.0f64 {
                        fits_parser_yyerror(
                            0 as yyscan_t,
                            lParse,
                            b"Out of range argument to log10\0" as *const u8 as *const c_char
                                as *mut c_char,
                        );
                    } else {
                        (*this).value.data.dbl = log10(dval);
                    }
                    current_block_139 = 7627602990488000394;
                }
                1016 => {
                    dval = pVals[0].data.dbl;
                    if dval < 0.0f64 {
                        fits_parser_yyerror(
                            0 as yyscan_t,
                            lParse,
                            b"Out of range argument to sqrt\0" as *const u8 as *const c_char
                                as *mut c_char,
                        );
                    } else {
                        (*this).value.data.dbl = sqrt(dval);
                    }
                    current_block_139 = 7627602990488000394;
                }
                1019 => {
                    (*this).value.data.dbl = ceil(pVals[0].data.dbl);
                    current_block_139 = 7627602990488000394;
                }
                1020 => {
                    (*this).value.data.dbl = floor(pVals[0].data.dbl);
                    current_block_139 = 7627602990488000394;
                }
                1021 => {
                    (*this).value.data.dbl = floor(pVals[0].data.dbl + 0.5f64);
                    current_block_139 = 7627602990488000394;
                }
                1018 => {
                    (*this).value.data.dbl = atan2(pVals[0].data.dbl, pVals[1].data.dbl);
                    current_block_139 = 7627602990488000394;
                }
                1041 => {
                    (*this).value.data.dbl = angsep_calc(
                        pVals[0].data.dbl,
                        pVals[1].data.dbl,
                        pVals[2].data.dbl,
                        pVals[3].data.dbl,
                    );
                    current_block_139 = 15934000668868306918;
                }
                1022 => {
                    current_block_139 = 15934000668868306918;
                }
                1023 => {
                    if (*this).ntype == fits_parser_yytokentype::DOUBLE as c_int {
                        (*this).value.data.dbl = if pVals[0].data.dbl < pVals[1].data.dbl {
                            pVals[0].data.dbl
                        } else {
                            pVals[1].data.dbl
                        };
                    } else if (*this).ntype == fits_parser_yytokentype::LONG as c_int {
                        (*this).value.data.lng = if pVals[0].data.lng < pVals[1].data.lng {
                            pVals[0].data.lng
                        } else {
                            pVals[1].data.lng
                        };
                    }
                    current_block_139 = 7627602990488000394;
                }
                1024 => {
                    if (*this).ntype == fits_parser_yytokentype::DOUBLE as c_int {
                        (*this).value.data.dbl = pVals[0].data.dbl;
                    } else if (*this).ntype == fits_parser_yytokentype::LONG as c_int {
                        (*this).value.data.lng = pVals[0].data.lng;
                    } else if (*this).ntype == fits_parser_yytokentype::BITSTR as c_int {
                        strcpy(
                            ((*this).value.data.astr).as_mut_ptr(),
                            (pVals[0].data.astr).as_mut_ptr(),
                        );
                    }
                    current_block_139 = 7627602990488000394;
                }
                1025 => {
                    if (*this).ntype == fits_parser_yytokentype::DOUBLE as c_int {
                        (*this).value.data.dbl = if pVals[0].data.dbl > pVals[1].data.dbl {
                            pVals[0].data.dbl
                        } else {
                            pVals[1].data.dbl
                        };
                    } else if (*this).ntype == fits_parser_yytokentype::LONG as c_int {
                        (*this).value.data.lng = if pVals[0].data.lng > pVals[1].data.lng {
                            pVals[0].data.lng
                        } else {
                            pVals[1].data.lng
                        };
                    }
                    current_block_139 = 7627602990488000394;
                }
                1026 => {
                    (*this).value.data.log =
                        bnear(pVals[0].data.dbl, pVals[1].data.dbl, pVals[2].data.dbl);
                    current_block_139 = 7627602990488000394;
                }
                1027 => {
                    (*this).value.data.log = circle(
                        pVals[0].data.dbl,
                        pVals[1].data.dbl,
                        pVals[2].data.dbl,
                        pVals[3].data.dbl,
                        pVals[4].data.dbl,
                    );
                    current_block_139 = 7627602990488000394;
                }
                1028 => {
                    (*this).value.data.log = saobox(
                        pVals[0].data.dbl,
                        pVals[1].data.dbl,
                        pVals[2].data.dbl,
                        pVals[3].data.dbl,
                        pVals[4].data.dbl,
                        pVals[5].data.dbl,
                        pVals[6].data.dbl,
                    );
                    current_block_139 = 7627602990488000394;
                }
                1029 => {
                    (*this).value.data.log = ellipse(
                        pVals[0].data.dbl,
                        pVals[1].data.dbl,
                        pVals[2].data.dbl,
                        pVals[3].data.dbl,
                        pVals[4].data.dbl,
                        pVals[5].data.dbl,
                        pVals[6].data.dbl,
                    );
                    current_block_139 = 7627602990488000394;
                }
                1034 => {
                    match (*this).ntype {
                        258 => {
                            (*this).value.data.log = (if pVals[2].data.log as c_int != 0 {
                                pVals[0].data.log as c_int
                            } else {
                                pVals[1].data.log as c_int
                            }) as c_char;
                        }
                        259 => {
                            (*this).value.data.lng = if pVals[2].data.log as c_int != 0 {
                                pVals[0].data.lng
                            } else {
                                pVals[1].data.lng
                            };
                        }
                        260 => {
                            (*this).value.data.dbl = if pVals[2].data.log as c_int != 0 {
                                pVals[0].data.dbl
                            } else {
                                pVals[1].data.dbl
                            };
                        }
                        261 => {
                            strcpy(
                                ((*this).value.data.astr).as_mut_ptr(),
                                if pVals[2].data.log as c_int != 0 {
                                    (pVals[0].data.astr).as_mut_ptr()
                                } else {
                                    (pVals[1].data.astr).as_mut_ptr()
                                },
                            );
                        }
                        _ => {}
                    }
                    current_block_139 = 7627602990488000394;
                }
                1044 => {
                    cstrmid(
                        lParse,
                        ((*this).value.data.astr).as_mut_ptr(),
                        (*this).value.nelem as c_int,
                        (pVals[0].data.astr).as_mut_ptr(),
                        pVals[0].nelem as c_int,
                        pVals[1].data.lng as c_int,
                    );
                    current_block_139 = 7627602990488000394;
                }
                1045 => {
                    let mut res: *mut c_char = strstr(
                        (pVals[0].data.astr).as_mut_ptr(),
                        (pVals[1].data.astr).as_mut_ptr(),
                    );
                    if res.is_null() {
                        (*this).value.data.lng = 0;
                    } else {
                        (*this).value.data.lng =
                            res.offset_from((pVals[0].data.astr).as_mut_ptr()) as c_long + 1;
                    }
                    current_block_139 = 7627602990488000394;
                }
                _ => {
                    current_block_139 = 7627602990488000394;
                }
            }
            if current_block_139 == 15934000668868306918 {
                if (*this).ntype == fits_parser_yytokentype::DOUBLE as c_int {
                    (*this).value.data.dbl = pVals[0].data.dbl;
                } else if (*this).ntype == fits_parser_yytokentype::LONG as c_int {
                    (*this).value.data.lng = pVals[0].data.lng;
                } else if (*this).ntype == fits_parser_yytokentype::BITSTR as c_int {
                    strcpy(
                        ((*this).value.data.astr).as_mut_ptr(),
                        (pVals[0].data.astr).as_mut_ptr(),
                    );
                }
            }
            (*this).operation = -(1000 as c_int);
        } else {
            Allocate_Ptrs(lParse, this);
            row = (*lParse).nRows;
            elem = row * (*this).value.nelem;
            if (*lParse).status == 0 {
                match (*this).operation {
                    1035 => loop {
                        let fresh53 = row;
                        row -= 1;
                        if fresh53 == 0 {
                            break;
                        }
                        *((*this).value.data.lngptr).offset(row as isize) =
                            (*lParse).firstRow + row;
                        *((*this).value.undef).offset(row as isize) = 0;
                    },
                    1036 => {
                        if (*this).ntype == fits_parser_yytokentype::LONG as c_int {
                            loop {
                                let fresh54 = row;
                                row -= 1;
                                if fresh54 == 0 {
                                    break;
                                }
                                *((*this).value.data.lngptr).offset(row as isize) = 0;
                                *((*this).value.undef).offset(row as isize) = 1;
                            }
                        } else if (*this).ntype == fits_parser_yytokentype::STRING as c_int {
                            loop {
                                let fresh55 = row;
                                row -= 1;
                                if fresh55 == 0 {
                                    break;
                                }
                                *(*((*this).value.data.strptr).offset(row as isize)).offset(0) = 0;
                                *((*this).value.undef).offset(row as isize) = 1;
                            }
                        }
                    }
                    1050 => {
                        let mut ielem: c_long = 0;
                        let mut iaxis: [c_long; 5] = [1, 1, 1, 1, 1];
                        let mut ipos: c_long = pVals[1].data.lng - 1;
                        let mut naxis: c_int = (*this).value.naxis;
                        let mut j: c_int = 0;
                        if ipos < 0 || ipos >= 5 as c_long {
                            fits_parser_yyerror(
                                0 as yyscan_t,
                                lParse,
                                b"AXISELEM(V,n) n value exceeded maximum dimension\0" as *const u8
                                    as *const c_char as *mut c_char,
                            );
                            free((*this).value.data.ptr);
                        } else {
                            ielem = 0;
                            while ielem < elem {
                                *((*this).value.data.lngptr).offset(ielem as isize) =
                                    iaxis[ipos as usize];
                                *((*this).value.undef).offset(ielem as isize) = 0;
                                iaxis[0] += 1;
                                iaxis[0];
                                j = 0;
                                while j < naxis {
                                    if iaxis[j as usize] <= (*this).value.naxes[j as usize] {
                                        break;
                                    }
                                    iaxis[j as usize] = 1;
                                    if j < naxis - 1 {
                                        iaxis[(j + 1) as usize] += 1;
                                        iaxis[(j + 1) as usize];
                                    }
                                    j += 1;
                                    j;
                                }
                                ielem += 1;
                                ielem;
                            }
                        }
                    }
                    1049 => {
                        let mut ielem_0: c_long = 0;
                        let mut elemnum: c_long = 1;
                        let mut j_0: c_int = 0;
                        ielem_0 = 0;
                        while ielem_0 < elem {
                            *((*this).value.data.lngptr).offset(ielem_0 as isize) = elemnum;
                            *((*this).value.undef).offset(ielem_0 as isize) = 0;
                            elemnum += 1;
                            elemnum;
                            if elemnum > (*this).value.nelem {
                                elemnum = 1;
                            }
                            ielem_0 += 1;
                            ielem_0;
                        }
                    }
                    1001 => loop {
                        let fresh56 = elem;
                        elem -= 1;
                        if fresh56 == 0 {
                            break;
                        }
                        *((*this).value.data.dblptr).offset(elem as isize) = simplerng_getuniform();
                        *((*this).value.undef).offset(elem as isize) = 0;
                    },
                    1042 => loop {
                        let fresh57 = elem;
                        elem -= 1;
                        if fresh57 == 0 {
                            break;
                        }
                        *((*this).value.data.dblptr).offset(elem as isize) = simplerng_getnorm();
                        *((*this).value.undef).offset(elem as isize) = 0;
                    },
                    1043 => {
                        if (*theParams[0]).ntype == fits_parser_yytokentype::DOUBLE as c_int {
                            if (*theParams[0]).operation == -(1000 as c_int) {
                                loop {
                                    let fresh58 = elem;
                                    elem -= 1;
                                    if fresh58 == 0 {
                                        break;
                                    }
                                    *((*this).value.undef).offset(elem as isize) =
                                        if pVals[0].data.dbl < 0.0 { 1 } else { 0 };
                                    if *((*this).value.undef).offset(elem as isize) == 0 {
                                        *((*this).value.data.lngptr).offset(elem as isize) =
                                            simplerng_getpoisson(pVals[0].data.dbl) as c_long;
                                    }
                                }
                            } else {
                                loop {
                                    let fresh59 = elem;
                                    elem -= 1;
                                    if fresh59 == 0 {
                                        break;
                                    }
                                    *((*this).value.undef).offset(elem as isize) =
                                        *((*theParams[0]).value.undef).offset(elem as isize);
                                    if *((*theParams[0]).value.data.dblptr).offset(elem as isize)
                                        < 0.0
                                    {
                                        *((*this).value.undef).offset(elem as isize) = 1;
                                    }
                                    if *((*this).value.undef).offset(elem as isize) == 0 {
                                        *((*this).value.data.lngptr).offset(elem as isize) =
                                            simplerng_getpoisson(
                                                *((*theParams[0]).value.data.dblptr)
                                                    .offset(elem as isize),
                                            ) as c_long;
                                    }
                                }
                            }
                        } else if (*theParams[0]).operation == -(1000 as c_int) {
                            loop {
                                let fresh60 = elem;
                                elem -= 1;
                                if fresh60 == 0 {
                                    break;
                                }
                                *((*this).value.undef).offset(elem as isize) =
                                    if pVals[0].data.lng < 0 { 1 } else { 0 };
                                if *((*this).value.undef).offset(elem as isize) == 0 {
                                    *((*this).value.data.lngptr).offset(elem as isize) =
                                        simplerng_getpoisson(pVals[0].data.lng as c_double)
                                            as c_long;
                                }
                            }
                        } else {
                            loop {
                                let fresh61 = elem;
                                elem -= 1;
                                if fresh61 == 0 {
                                    break;
                                }
                                *((*this).value.undef).offset(elem as isize) =
                                    *((*theParams[0]).value.undef).offset(elem as isize);
                                if *((*theParams[0]).value.data.lngptr).offset(elem as isize) < 0 {
                                    *((*this).value.undef).offset(elem as isize) = 1;
                                }
                                if *((*this).value.undef).offset(elem as isize) == 0 {
                                    *((*this).value.data.lngptr).offset(elem as isize) =
                                        simplerng_getpoisson(
                                            *((*theParams[0]).value.data.lngptr)
                                                .offset(elem as isize)
                                                as c_double,
                                        ) as c_long;
                                }
                            }
                        }
                    }
                    1002 => {
                        elem = row * (*theParams[0]).value.nelem;
                        if (*theParams[0]).ntype == fits_parser_yytokentype::BOOLEAN as c_int {
                            loop {
                                let fresh62 = row;
                                row -= 1;
                                if fresh62 == 0 {
                                    break;
                                }
                                *((*this).value.data.lngptr).offset(row as isize) = 0;
                                *((*this).value.undef).offset(row as isize) = 1;
                                nelem = (*theParams[0]).value.nelem;
                                loop {
                                    let fresh63 = nelem;
                                    nelem -= 1;
                                    if fresh63 == 0 {
                                        break;
                                    }
                                    elem -= 1;
                                    elem;
                                    if *((*theParams[0]).value.undef).offset(elem as isize) == 0 {
                                        *((*this).value.data.lngptr).offset(row as isize) +=
                                            (if *((*theParams[0]).value.data.logptr)
                                                .offset(elem as isize)
                                                as c_int
                                                != 0
                                            {
                                                1
                                            } else {
                                                0
                                            })
                                                as c_long;
                                        *((*this).value.undef).offset(row as isize) = 0;
                                    }
                                }
                            }
                        } else if (*theParams[0]).ntype == fits_parser_yytokentype::LONG as c_int {
                            loop {
                                let fresh64 = row;
                                row -= 1;
                                if fresh64 == 0 {
                                    break;
                                }
                                *((*this).value.data.lngptr).offset(row as isize) = 0;
                                *((*this).value.undef).offset(row as isize) = 1;
                                nelem = (*theParams[0]).value.nelem;
                                loop {
                                    let fresh65 = nelem;
                                    nelem -= 1;
                                    if fresh65 == 0 {
                                        break;
                                    }
                                    elem -= 1;
                                    elem;
                                    if *((*theParams[0]).value.undef).offset(elem as isize) == 0 {
                                        *((*this).value.data.lngptr).offset(row as isize) +=
                                            *((*theParams[0]).value.data.lngptr)
                                                .offset(elem as isize);
                                        *((*this).value.undef).offset(row as isize) = 0;
                                    }
                                }
                            }
                        } else if (*theParams[0]).ntype == fits_parser_yytokentype::DOUBLE as c_int
                        {
                            loop {
                                let fresh66 = row;
                                row -= 1;
                                if fresh66 == 0 {
                                    break;
                                }
                                *((*this).value.data.dblptr).offset(row as isize) = 0.0f64;
                                *((*this).value.undef).offset(row as isize) = 1;
                                nelem = (*theParams[0]).value.nelem;
                                loop {
                                    let fresh67 = nelem;
                                    nelem -= 1;
                                    if fresh67 == 0 {
                                        break;
                                    }
                                    elem -= 1;
                                    elem;
                                    if *((*theParams[0]).value.undef).offset(elem as isize) == 0 {
                                        *((*this).value.data.dblptr).offset(row as isize) +=
                                            *((*theParams[0]).value.data.dblptr)
                                                .offset(elem as isize);
                                        *((*this).value.undef).offset(row as isize) = 0;
                                    }
                                }
                            }
                        } else {
                            nelem = (*theParams[0]).value.nelem;
                            loop {
                                let fresh68 = row;
                                row -= 1;
                                if fresh68 == 0 {
                                    break;
                                }
                                let mut sptr1: *mut c_char =
                                    *((*theParams[0]).value.data.strptr).offset(row as isize);
                                *((*this).value.data.lngptr).offset(row as isize) = 0;
                                *((*this).value.undef).offset(row as isize) = 0;
                                while *sptr1 != 0 {
                                    if *sptr1 as c_int == '1' as i32 {
                                        let fresh69 =
                                            &mut *((*this).value.data.lngptr).offset(row as isize);
                                        *fresh69 += 1;
                                        *fresh69;
                                    }
                                    sptr1 = sptr1.offset(1);
                                    sptr1;
                                }
                            }
                        }
                    }
                    1038 => {
                        elem = row * (*theParams[0]).value.nelem;
                        if (*theParams[0]).ntype == fits_parser_yytokentype::LONG as c_int {
                            loop {
                                let fresh70 = row;
                                row -= 1;
                                if fresh70 == 0 {
                                    break;
                                }
                                let mut count: c_int = 0;
                                *((*this).value.data.dblptr).offset(row as isize) = 0.0;
                                nelem = (*theParams[0]).value.nelem;
                                loop {
                                    let fresh71 = nelem;
                                    nelem -= 1;
                                    if fresh71 == 0 {
                                        break;
                                    }
                                    elem -= 1;
                                    elem;
                                    if *((*theParams[0]).value.undef).offset(elem as isize) as c_int
                                        == 0
                                    {
                                        *((*this).value.data.dblptr).offset(row as isize) +=
                                            *((*theParams[0]).value.data.lngptr)
                                                .offset(elem as isize)
                                                as c_double;
                                        count += 1;
                                        count;
                                    }
                                }
                                if count == 0 {
                                    *((*this).value.undef).offset(row as isize) = 1;
                                } else {
                                    *((*this).value.undef).offset(row as isize) = 0;
                                    *((*this).value.data.dblptr).offset(row as isize) /=
                                        count as c_double;
                                }
                            }
                        } else if (*theParams[0]).ntype == fits_parser_yytokentype::DOUBLE as c_int
                        {
                            loop {
                                let fresh72 = row;
                                row -= 1;
                                if fresh72 == 0 {
                                    break;
                                }
                                let mut count_0: c_int = 0;
                                *((*this).value.data.dblptr).offset(row as isize) = 0.0;
                                nelem = (*theParams[0]).value.nelem;
                                loop {
                                    let fresh73 = nelem;
                                    nelem -= 1;
                                    if fresh73 == 0 {
                                        break;
                                    }
                                    elem -= 1;
                                    elem;
                                    if *((*theParams[0]).value.undef).offset(elem as isize) as c_int
                                        == 0
                                    {
                                        *((*this).value.data.dblptr).offset(row as isize) +=
                                            *((*theParams[0]).value.data.dblptr)
                                                .offset(elem as isize);
                                        count_0 += 1;
                                        count_0;
                                    }
                                }
                                if count_0 == 0 {
                                    *((*this).value.undef).offset(row as isize) = 1;
                                } else {
                                    *((*this).value.undef).offset(row as isize) = 0;
                                    *((*this).value.data.dblptr).offset(row as isize) /=
                                        count_0 as c_double;
                                }
                            }
                        }
                    }
                    1039 => {
                        elem = row * (*theParams[0]).value.nelem;
                        if (*theParams[0]).ntype == fits_parser_yytokentype::LONG as c_int {
                            loop {
                                let fresh74 = row;
                                row -= 1;
                                if fresh74 == 0 {
                                    break;
                                }
                                let mut count_1: c_int = 0;
                                let mut sum: c_double = 0.0;
                                let mut sum2: c_double = 0.0;
                                nelem = (*theParams[0]).value.nelem;
                                loop {
                                    let fresh75 = nelem;
                                    nelem -= 1;
                                    if fresh75 == 0 {
                                        break;
                                    }
                                    elem -= 1;
                                    elem;
                                    if *((*theParams[0]).value.undef).offset(elem as isize) as c_int
                                        == 0
                                    {
                                        sum += *((*theParams[0]).value.data.lngptr)
                                            .offset(elem as isize)
                                            as c_double;
                                        count_1 += 1;
                                        count_1;
                                    }
                                }
                                if count_1 > 1 {
                                    sum /= count_1 as c_double;
                                    nelem = (*theParams[0]).value.nelem;
                                    elem += nelem;
                                    loop {
                                        let fresh76 = nelem;
                                        nelem -= 1;
                                        if fresh76 == 0 {
                                            break;
                                        }
                                        elem -= 1;
                                        elem;
                                        if *((*theParams[0]).value.undef).offset(elem as isize)
                                            as c_int
                                            == 0
                                        {
                                            let mut dx: c_double =
                                                *((*theParams[0]).value.data.lngptr)
                                                    .offset(elem as isize)
                                                    as c_double
                                                    - sum;
                                            sum2 += dx * dx;
                                        }
                                    }
                                    sum2 /= count_1 as c_double - 1 as c_double;
                                    *((*this).value.undef).offset(row as isize) = 0;
                                    *((*this).value.data.dblptr).offset(row as isize) = sqrt(sum2);
                                } else {
                                    *((*this).value.undef).offset(row as isize) = 0;
                                    *((*this).value.data.dblptr).offset(row as isize) = 0.0;
                                }
                            }
                        } else if (*theParams[0]).ntype == fits_parser_yytokentype::DOUBLE as c_int
                        {
                            loop {
                                let fresh77 = row;
                                row -= 1;
                                if fresh77 == 0 {
                                    break;
                                }
                                let mut count_2: c_int = 0;
                                let mut sum_0: c_double = 0.0;
                                let mut sum2_0: c_double = 0.0;
                                nelem = (*theParams[0]).value.nelem;
                                loop {
                                    let fresh78 = nelem;
                                    nelem -= 1;
                                    if fresh78 == 0 {
                                        break;
                                    }
                                    elem -= 1;
                                    elem;
                                    if *((*theParams[0]).value.undef).offset(elem as isize) as c_int
                                        == 0
                                    {
                                        sum_0 += *((*theParams[0]).value.data.dblptr)
                                            .offset(elem as isize);
                                        count_2 += 1;
                                        count_2;
                                    }
                                }
                                if count_2 > 1 {
                                    sum_0 /= count_2 as c_double;
                                    nelem = (*theParams[0]).value.nelem;
                                    elem += nelem;
                                    loop {
                                        let fresh79 = nelem;
                                        nelem -= 1;
                                        if fresh79 == 0 {
                                            break;
                                        }
                                        elem -= 1;
                                        elem;
                                        if *((*theParams[0]).value.undef).offset(elem as isize)
                                            as c_int
                                            == 0
                                        {
                                            let mut dx_0: c_double =
                                                *((*theParams[0]).value.data.dblptr)
                                                    .offset(elem as isize)
                                                    - sum_0;
                                            sum2_0 += dx_0 * dx_0;
                                        }
                                    }
                                    sum2_0 /= count_2 as c_double - 1 as c_double;
                                    *((*this).value.undef).offset(row as isize) = 0;
                                    *((*this).value.data.dblptr).offset(row as isize) =
                                        sqrt(sum2_0);
                                } else {
                                    *((*this).value.undef).offset(row as isize) = 0;
                                    *((*this).value.data.dblptr).offset(row as isize) = 0.0;
                                }
                            }
                        }
                    }
                    1037 => {
                        elem = row * (*theParams[0]).value.nelem;
                        nelem = (*theParams[0]).value.nelem;
                        if (*theParams[0]).ntype == fits_parser_yytokentype::LONG as c_int {
                            let mut dptr: *mut c_long = (*theParams[0]).value.data.lngptr;
                            let mut uptr: *mut c_char = (*theParams[0]).value.undef;
                            let mut mptr: *mut c_long = malloc(
                                (::core::mem::size_of::<c_long>() as c_ulong)
                                    .wrapping_mul(nelem as c_ulong)
                                    .try_into()
                                    .unwrap(),
                            )
                                as *mut c_long;
                            let mut irow: c_int = 0;
                            if mptr.is_null() {
                                fits_parser_yyerror(
                                    0 as yyscan_t,
                                    lParse,
                                    b"Could not allocate temporary memory in median function\0"
                                        as *const u8
                                        as *const c_char
                                        as *mut c_char,
                                );
                                free((*this).value.data.ptr);
                            } else {
                                irow = 0;
                                while (irow as c_long) < row {
                                    let mut p: *mut c_long = mptr;
                                    let mut nelem1: c_int = nelem as c_int;
                                    loop {
                                        let fresh80 = nelem1;
                                        nelem1 -= 1;
                                        if fresh80 == 0 {
                                            break;
                                        }
                                        if *uptr as c_int == 0 {
                                            let fresh81 = p;
                                            p = p.offset(1);
                                            *fresh81 = *dptr;
                                        }
                                        dptr = dptr.offset(1);
                                        dptr;
                                        uptr = uptr.offset(1);
                                        uptr;
                                    }
                                    nelem1 = p.offset_from(mptr) as c_long as c_int;
                                    if nelem1 > 0 {
                                        *((*this).value.undef).offset(irow as isize) = 0;
                                        *((*this).value.data.lngptr).offset(irow as isize) =
                                            qselect_median_lng(mptr, nelem1);
                                    } else {
                                        *((*this).value.undef).offset(irow as isize) = 1;
                                        *((*this).value.data.lngptr).offset(irow as isize) = 0;
                                    }
                                    irow += 1;
                                    irow;
                                }
                                free(mptr as *mut c_void);
                            }
                        } else {
                            let mut dptr_0: *mut c_double = (*theParams[0]).value.data.dblptr;
                            let mut uptr_0: *mut c_char = (*theParams[0]).value.undef;
                            let mut mptr_0: *mut c_double = malloc(
                                (::core::mem::size_of::<c_double>() as c_ulong)
                                    .wrapping_mul(nelem as c_ulong)
                                    .try_into()
                                    .unwrap(),
                            )
                                as *mut c_double;
                            let mut irow_0: c_int = 0;
                            if mptr_0.is_null() {
                                fits_parser_yyerror(
                                    0 as yyscan_t,
                                    lParse,
                                    b"Could not allocate temporary memory in median function\0"
                                        as *const u8
                                        as *const c_char
                                        as *mut c_char,
                                );
                                free((*this).value.data.ptr);
                            } else {
                                irow_0 = 0;
                                while (irow_0 as c_long) < row {
                                    let mut p_0: *mut c_double = mptr_0;
                                    let mut nelem1_0: c_int = nelem as c_int;
                                    loop {
                                        let fresh82 = nelem1_0;
                                        nelem1_0 -= 1;
                                        if fresh82 == 0 {
                                            break;
                                        }
                                        if *uptr_0 as c_int == 0 {
                                            let fresh83 = p_0;
                                            p_0 = p_0.offset(1);
                                            *fresh83 = *dptr_0;
                                        }
                                        dptr_0 = dptr_0.offset(1);
                                        dptr_0;
                                        uptr_0 = uptr_0.offset(1);
                                        uptr_0;
                                    }
                                    nelem1_0 = p_0.offset_from(mptr_0) as c_long as c_int;
                                    if nelem1_0 > 0 {
                                        *((*this).value.undef).offset(irow_0 as isize) = 0;
                                        *((*this).value.data.dblptr).offset(irow_0 as isize) =
                                            qselect_median_dbl(mptr_0, nelem1_0);
                                    } else {
                                        *((*this).value.undef).offset(irow_0 as isize) = 1;
                                        *((*this).value.data.dblptr).offset(irow_0 as isize) = 0.0;
                                    }
                                    irow_0 += 1;
                                    irow_0;
                                }
                                free(mptr_0 as *mut c_void);
                            }
                        }
                    }
                    1017 => {
                        if (*theParams[0]).ntype == fits_parser_yytokentype::DOUBLE as c_int {
                            loop {
                                let fresh84 = elem;
                                elem -= 1;
                                if fresh84 == 0 {
                                    break;
                                }
                                dval = *((*theParams[0]).value.data.dblptr).offset(elem as isize);
                                *((*this).value.data.dblptr).offset(elem as isize) =
                                    if dval > 0.0f64 { dval } else { -dval };
                                *((*this).value.undef).offset(elem as isize) =
                                    *((*theParams[0]).value.undef).offset(elem as isize);
                            }
                        } else {
                            loop {
                                let fresh85 = elem;
                                elem -= 1;
                                if fresh85 == 0 {
                                    break;
                                }
                                ival = *((*theParams[0]).value.data.lngptr).offset(elem as isize);
                                *((*this).value.data.lngptr).offset(elem as isize) =
                                    if ival > 0 { ival } else { -ival };
                                *((*this).value.undef).offset(elem as isize) =
                                    *((*theParams[0]).value.undef).offset(elem as isize);
                            }
                        }
                    }
                    1040 => {
                        nelem = (*theParams[0]).value.nelem;
                        if (*theParams[0]).ntype == fits_parser_yytokentype::STRING as c_int {
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
                            *((*this).value.undef).offset(row as isize) = 0;
                            *((*this).value.data.lngptr).offset(row as isize) = 0;
                            loop {
                                let fresh87 = nelem1_1;
                                nelem1_1 -= 1;
                                if fresh87 == 0 {
                                    break;
                                }
                                elem -= 1;
                                elem;
                                if *((*theParams[0]).value.undef).offset(elem as isize) as c_int
                                    == 0
                                {
                                    let fresh88 =
                                        &mut *((*this).value.data.lngptr).offset(row as isize);
                                    *fresh88 += 1;
                                    *fresh88;
                                }
                            }
                        }
                    }
                    1030 => {
                        if (*theParams[0]).ntype == fits_parser_yytokentype::STRING as c_int {
                            elem = row;
                        }
                        loop {
                            let fresh89 = elem;
                            elem -= 1;
                            if fresh89 == 0 {
                                break;
                            }
                            *((*this).value.data.logptr).offset(elem as isize) =
                                *((*theParams[0]).value.undef).offset(elem as isize);
                            *((*this).value.undef).offset(elem as isize) = 0;
                        }
                    }
                    1031 => match (*this).ntype {
                        258 => loop {
                            let fresh90 = row;
                            row -= 1;
                            if fresh90 == 0 {
                                break;
                            }
                            nelem = (*this).value.nelem;
                            loop {
                                let fresh91 = nelem;
                                nelem -= 1;
                                if fresh91 == 0 {
                                    break;
                                }
                                elem -= 1;
                                elem;
                                i = 2 as c_int;
                                loop {
                                    let fresh92 = i;
                                    i -= 1;
                                    if fresh92 == 0 {
                                        break;
                                    }
                                    if vector[i as usize] > 1 {
                                        pNull[i as usize] = *((*theParams[i as usize]).value.undef)
                                            .offset(elem as isize);
                                        pVals[i as usize].data.log =
                                            *((*theParams[i as usize]).value.data.logptr)
                                                .offset(elem as isize);
                                    } else if vector[i as usize] != 0 {
                                        pNull[i as usize] = *((*theParams[i as usize]).value.undef)
                                            .offset(row as isize);
                                        pVals[i as usize].data.log =
                                            *((*theParams[i as usize]).value.data.logptr)
                                                .offset(row as isize);
                                    }
                                }
                                if pNull[0] != 0 {
                                    *((*this).value.undef).offset(elem as isize) = pNull[1];
                                    *((*this).value.data.logptr).offset(elem as isize) =
                                        pVals[1].data.log;
                                } else {
                                    *((*this).value.undef).offset(elem as isize) = 0;
                                    *((*this).value.data.logptr).offset(elem as isize) =
                                        pVals[0].data.log;
                                }
                            }
                        },
                        259 => loop {
                            let fresh93 = row;
                            row -= 1;
                            if fresh93 == 0 {
                                break;
                            }
                            nelem = (*this).value.nelem;
                            loop {
                                let fresh94 = nelem;
                                nelem -= 1;
                                if fresh94 == 0 {
                                    break;
                                }
                                elem -= 1;
                                elem;
                                i = 2 as c_int;
                                loop {
                                    let fresh95 = i;
                                    i -= 1;
                                    if fresh95 == 0 {
                                        break;
                                    }
                                    if vector[i as usize] > 1 {
                                        pNull[i as usize] = *((*theParams[i as usize]).value.undef)
                                            .offset(elem as isize);
                                        pVals[i as usize].data.lng =
                                            *((*theParams[i as usize]).value.data.lngptr)
                                                .offset(elem as isize);
                                    } else if vector[i as usize] != 0 {
                                        pNull[i as usize] = *((*theParams[i as usize]).value.undef)
                                            .offset(row as isize);
                                        pVals[i as usize].data.lng =
                                            *((*theParams[i as usize]).value.data.lngptr)
                                                .offset(row as isize);
                                    }
                                }
                                if pNull[0] != 0 {
                                    *((*this).value.undef).offset(elem as isize) = pNull[1];
                                    *((*this).value.data.lngptr).offset(elem as isize) =
                                        pVals[1].data.lng;
                                } else {
                                    *((*this).value.undef).offset(elem as isize) = 0;
                                    *((*this).value.data.lngptr).offset(elem as isize) =
                                        pVals[0].data.lng;
                                }
                            }
                        },
                        260 => loop {
                            let fresh96 = row;
                            row -= 1;
                            if fresh96 == 0 {
                                break;
                            }
                            nelem = (*this).value.nelem;
                            loop {
                                let fresh97 = nelem;
                                nelem -= 1;
                                if fresh97 == 0 {
                                    break;
                                }
                                elem -= 1;
                                elem;
                                i = 2 as c_int;
                                loop {
                                    let fresh98 = i;
                                    i -= 1;
                                    if fresh98 == 0 {
                                        break;
                                    }
                                    if vector[i as usize] > 1 {
                                        pNull[i as usize] = *((*theParams[i as usize]).value.undef)
                                            .offset(elem as isize);
                                        pVals[i as usize].data.dbl =
                                            *((*theParams[i as usize]).value.data.dblptr)
                                                .offset(elem as isize);
                                    } else if vector[i as usize] != 0 {
                                        pNull[i as usize] = *((*theParams[i as usize]).value.undef)
                                            .offset(row as isize);
                                        pVals[i as usize].data.dbl =
                                            *((*theParams[i as usize]).value.data.dblptr)
                                                .offset(row as isize);
                                    }
                                }
                                if pNull[0] != 0 {
                                    *((*this).value.undef).offset(elem as isize) = pNull[1];
                                    *((*this).value.data.dblptr).offset(elem as isize) =
                                        pVals[1].data.dbl;
                                } else {
                                    *((*this).value.undef).offset(elem as isize) = 0;
                                    *((*this).value.data.dblptr).offset(elem as isize) =
                                        pVals[0].data.dbl;
                                }
                            }
                        },
                        261 => loop {
                            let fresh99 = row;
                            row -= 1;
                            if fresh99 == 0 {
                                break;
                            }
                            i = 2 as c_int;
                            loop {
                                let fresh100 = i;
                                i -= 1;
                                if fresh100 == 0 {
                                    break;
                                }
                                if vector[i as usize] != 0 {
                                    pNull[i as usize] = *((*theParams[i as usize]).value.undef)
                                        .offset(row as isize);
                                    strcpy(
                                        (pVals[i as usize].data.astr).as_mut_ptr(),
                                        *((*theParams[i as usize]).value.data.strptr)
                                            .offset(row as isize),
                                    );
                                }
                            }
                            if pNull[0] != 0 {
                                *((*this).value.undef).offset(row as isize) = pNull[1];
                                strcpy(
                                    *((*this).value.data.strptr).offset(row as isize),
                                    (pVals[1].data.astr).as_mut_ptr(),
                                );
                            } else {
                                *((*this).value.undef).offset(elem as isize) = 0;
                                strcpy(
                                    *((*this).value.data.strptr).offset(row as isize),
                                    (pVals[0].data.astr).as_mut_ptr(),
                                );
                            }
                        },
                        _ => {}
                    },
                    1046 => match (*this).ntype {
                        259 => loop {
                            let fresh101 = elem;
                            elem -= 1;
                            if fresh101 == 0 {
                                break;
                            }
                            if (*theParams[1]).value.data.lng
                                == *((*theParams[0]).value.data.lngptr).offset(elem as isize)
                            {
                                *((*this).value.data.lngptr).offset(elem as isize) = 0;
                                *((*this).value.undef).offset(elem as isize) = 1;
                            } else {
                                *((*this).value.data.lngptr).offset(elem as isize) =
                                    *((*theParams[0]).value.data.lngptr).offset(elem as isize);
                                *((*this).value.undef).offset(elem as isize) =
                                    *((*theParams[0]).value.undef).offset(elem as isize);
                            }
                        },
                        260 => loop {
                            let fresh102 = elem;
                            elem -= 1;
                            if fresh102 == 0 {
                                break;
                            }
                            if (*theParams[1]).value.data.dbl
                                == *((*theParams[0]).value.data.dblptr).offset(elem as isize)
                            {
                                *((*this).value.data.dblptr).offset(elem as isize) = 0.0;
                                *((*this).value.undef).offset(elem as isize) = 1;
                            } else {
                                *((*this).value.data.dblptr).offset(elem as isize) =
                                    *((*theParams[0]).value.data.dblptr).offset(elem as isize);
                                *((*this).value.undef).offset(elem as isize) =
                                    *((*theParams[0]).value.undef).offset(elem as isize);
                            }
                        },
                        _ => {}
                    },
                    1004 => loop {
                        let fresh103 = elem;
                        elem -= 1;
                        if fresh103 == 0 {
                            break;
                        }
                        let fresh104 = &mut *((*this).value.undef).offset(elem as isize);
                        *fresh104 = *((*theParams[0]).value.undef).offset(elem as isize);
                        if *fresh104 == 0 {
                            *((*this).value.data.dblptr).offset(elem as isize) =
                                sin(*((*theParams[0]).value.data.dblptr).offset(elem as isize));
                        }
                    },
                    1005 => loop {
                        let fresh105 = elem;
                        elem -= 1;
                        if fresh105 == 0 {
                            break;
                        }
                        let fresh106 = &mut *((*this).value.undef).offset(elem as isize);
                        *fresh106 = *((*theParams[0]).value.undef).offset(elem as isize);
                        if *fresh106 == 0 {
                            *((*this).value.data.dblptr).offset(elem as isize) =
                                cos(*((*theParams[0]).value.data.dblptr).offset(elem as isize));
                        }
                    },
                    1006 => loop {
                        let fresh107 = elem;
                        elem -= 1;
                        if fresh107 == 0 {
                            break;
                        }
                        let fresh108 = &mut *((*this).value.undef).offset(elem as isize);
                        *fresh108 = *((*theParams[0]).value.undef).offset(elem as isize);
                        if *fresh108 == 0 {
                            *((*this).value.data.dblptr).offset(elem as isize) =
                                tan(*((*theParams[0]).value.data.dblptr).offset(elem as isize));
                        }
                    },
                    1007 => loop {
                        let fresh109 = elem;
                        elem -= 1;
                        if fresh109 == 0 {
                            break;
                        }
                        let fresh110 = &mut *((*this).value.undef).offset(elem as isize);
                        *fresh110 = *((*theParams[0]).value.undef).offset(elem as isize);
                        if *fresh110 == 0 {
                            dval = *((*theParams[0]).value.data.dblptr).offset(elem as isize);
                            if dval < -1.0f64 || dval > 1.0f64 {
                                *((*this).value.data.dblptr).offset(elem as isize) = 0.0f64;
                                *((*this).value.undef).offset(elem as isize) = 1;
                            } else {
                                *((*this).value.data.dblptr).offset(elem as isize) = asin(dval);
                            }
                        }
                    },
                    1008 => loop {
                        let fresh111 = elem;
                        elem -= 1;
                        if fresh111 == 0 {
                            break;
                        }
                        let fresh112 = &mut *((*this).value.undef).offset(elem as isize);
                        *fresh112 = *((*theParams[0]).value.undef).offset(elem as isize);
                        if *fresh112 == 0 {
                            dval = *((*theParams[0]).value.data.dblptr).offset(elem as isize);
                            if dval < -1.0f64 || dval > 1.0f64 {
                                *((*this).value.data.dblptr).offset(elem as isize) = 0.0f64;
                                *((*this).value.undef).offset(elem as isize) = 1;
                            } else {
                                *((*this).value.data.dblptr).offset(elem as isize) = acos(dval);
                            }
                        }
                    },
                    1009 => loop {
                        let fresh113 = elem;
                        elem -= 1;
                        if fresh113 == 0 {
                            break;
                        }
                        let fresh114 = &mut *((*this).value.undef).offset(elem as isize);
                        *fresh114 = *((*theParams[0]).value.undef).offset(elem as isize);
                        if *fresh114 == 0 {
                            dval = *((*theParams[0]).value.data.dblptr).offset(elem as isize);
                            *((*this).value.data.dblptr).offset(elem as isize) = atan(dval);
                        }
                    },
                    1010 => loop {
                        let fresh115 = elem;
                        elem -= 1;
                        if fresh115 == 0 {
                            break;
                        }
                        let fresh116 = &mut *((*this).value.undef).offset(elem as isize);
                        *fresh116 = *((*theParams[0]).value.undef).offset(elem as isize);
                        if *fresh116 == 0 {
                            *((*this).value.data.dblptr).offset(elem as isize) =
                                sinh(*((*theParams[0]).value.data.dblptr).offset(elem as isize));
                        }
                    },
                    1011 => loop {
                        let fresh117 = elem;
                        elem -= 1;
                        if fresh117 == 0 {
                            break;
                        }
                        let fresh118 = &mut *((*this).value.undef).offset(elem as isize);
                        *fresh118 = *((*theParams[0]).value.undef).offset(elem as isize);
                        if *fresh118 == 0 {
                            *((*this).value.data.dblptr).offset(elem as isize) =
                                cosh(*((*theParams[0]).value.data.dblptr).offset(elem as isize));
                        }
                    },
                    1012 => loop {
                        let fresh119 = elem;
                        elem -= 1;
                        if fresh119 == 0 {
                            break;
                        }
                        let fresh120 = &mut *((*this).value.undef).offset(elem as isize);
                        *fresh120 = *((*theParams[0]).value.undef).offset(elem as isize);
                        if *fresh120 == 0 {
                            *((*this).value.data.dblptr).offset(elem as isize) =
                                tanh(*((*theParams[0]).value.data.dblptr).offset(elem as isize));
                        }
                    },
                    1013 => loop {
                        let fresh121 = elem;
                        elem -= 1;
                        if fresh121 == 0 {
                            break;
                        }
                        let fresh122 = &mut *((*this).value.undef).offset(elem as isize);
                        *fresh122 = *((*theParams[0]).value.undef).offset(elem as isize);
                        if *fresh122 == 0 {
                            dval = *((*theParams[0]).value.data.dblptr).offset(elem as isize);
                            *((*this).value.data.dblptr).offset(elem as isize) = exp(dval);
                        }
                    },
                    1014 => loop {
                        let fresh123 = elem;
                        elem -= 1;
                        if fresh123 == 0 {
                            break;
                        }
                        let fresh124 = &mut *((*this).value.undef).offset(elem as isize);
                        *fresh124 = *((*theParams[0]).value.undef).offset(elem as isize);
                        if *fresh124 == 0 {
                            dval = *((*theParams[0]).value.data.dblptr).offset(elem as isize);
                            if dval <= 0.0f64 {
                                *((*this).value.data.dblptr).offset(elem as isize) = 0.0f64;
                                *((*this).value.undef).offset(elem as isize) = 1;
                            } else {
                                *((*this).value.data.dblptr).offset(elem as isize) = log(dval);
                            }
                        }
                    },
                    1015 => loop {
                        let fresh125 = elem;
                        elem -= 1;
                        if fresh125 == 0 {
                            break;
                        }
                        let fresh126 = &mut *((*this).value.undef).offset(elem as isize);
                        *fresh126 = *((*theParams[0]).value.undef).offset(elem as isize);
                        if *fresh126 == 0 {
                            dval = *((*theParams[0]).value.data.dblptr).offset(elem as isize);
                            if dval <= 0.0f64 {
                                *((*this).value.data.dblptr).offset(elem as isize) = 0.0f64;
                                *((*this).value.undef).offset(elem as isize) = 1;
                            } else {
                                *((*this).value.data.dblptr).offset(elem as isize) = log10(dval);
                            }
                        }
                    },
                    1016 => loop {
                        let fresh127 = elem;
                        elem -= 1;
                        if fresh127 == 0 {
                            break;
                        }
                        let fresh128 = &mut *((*this).value.undef).offset(elem as isize);
                        *fresh128 = *((*theParams[0]).value.undef).offset(elem as isize);
                        if *fresh128 == 0 {
                            dval = *((*theParams[0]).value.data.dblptr).offset(elem as isize);
                            if dval < 0.0f64 {
                                *((*this).value.data.dblptr).offset(elem as isize) = 0.0f64;
                                *((*this).value.undef).offset(elem as isize) = 1;
                            } else {
                                *((*this).value.data.dblptr).offset(elem as isize) = sqrt(dval);
                            }
                        }
                    },
                    1019 => loop {
                        let fresh129 = elem;
                        elem -= 1;
                        if fresh129 == 0 {
                            break;
                        }
                        let fresh130 = &mut *((*this).value.undef).offset(elem as isize);
                        *fresh130 = *((*theParams[0]).value.undef).offset(elem as isize);
                        if *fresh130 == 0 {
                            *((*this).value.data.dblptr).offset(elem as isize) =
                                ceil(*((*theParams[0]).value.data.dblptr).offset(elem as isize));
                        }
                    },
                    1020 => loop {
                        let fresh131 = elem;
                        elem -= 1;
                        if fresh131 == 0 {
                            break;
                        }
                        let fresh132 = &mut *((*this).value.undef).offset(elem as isize);
                        *fresh132 = *((*theParams[0]).value.undef).offset(elem as isize);
                        if *fresh132 == 0 {
                            *((*this).value.data.dblptr).offset(elem as isize) =
                                floor(*((*theParams[0]).value.data.dblptr).offset(elem as isize));
                        }
                    },
                    1021 => loop {
                        let fresh133 = elem;
                        elem -= 1;
                        if fresh133 == 0 {
                            break;
                        }
                        let fresh134 = &mut *((*this).value.undef).offset(elem as isize);
                        *fresh134 = *((*theParams[0]).value.undef).offset(elem as isize);
                        if *fresh134 == 0 {
                            *((*this).value.data.dblptr).offset(elem as isize) = floor(
                                *((*theParams[0]).value.data.dblptr).offset(elem as isize) + 0.5f64,
                            );
                        }
                    },
                    1018 => loop {
                        let fresh135 = row;
                        row -= 1;
                        if fresh135 == 0 {
                            break;
                        }
                        nelem = (*this).value.nelem;
                        loop {
                            let fresh136 = nelem;
                            nelem -= 1;
                            if fresh136 == 0 {
                                break;
                            }
                            elem -= 1;
                            elem;
                            i = 2 as c_int;
                            loop {
                                let fresh137 = i;
                                i -= 1;
                                if fresh137 == 0 {
                                    break;
                                }
                                if vector[i as usize] > 1 {
                                    pVals[i as usize].data.dbl =
                                        *((*theParams[i as usize]).value.data.dblptr)
                                            .offset(elem as isize);
                                    pNull[i as usize] = *((*theParams[i as usize]).value.undef)
                                        .offset(elem as isize);
                                } else if vector[i as usize] != 0 {
                                    pVals[i as usize].data.dbl =
                                        *((*theParams[i as usize]).value.data.dblptr)
                                            .offset(row as isize);
                                    pNull[i as usize] = *((*theParams[i as usize]).value.undef)
                                        .offset(row as isize);
                                }
                            }
                            let fresh138 = &mut *((*this).value.undef).offset(elem as isize);
                            *fresh138 = if pNull[0] as c_int != 0 || pNull[1] as c_int != 0 {
                                1
                            } else {
                                0
                            };
                            if *fresh138 == 0 {
                                *((*this).value.data.dblptr).offset(elem as isize) =
                                    atan2(pVals[0].data.dbl, pVals[1].data.dbl);
                            }
                        }
                    },
                    1041 => loop {
                        let fresh139 = row;
                        row -= 1;
                        if fresh139 == 0 {
                            break;
                        }
                        nelem = (*this).value.nelem;
                        loop {
                            let fresh140 = nelem;
                            nelem -= 1;
                            if fresh140 == 0 {
                                break;
                            }
                            elem -= 1;
                            elem;
                            i = 4 as c_int;
                            loop {
                                let fresh141 = i;
                                i -= 1;
                                if fresh141 == 0 {
                                    break;
                                }
                                if vector[i as usize] > 1 {
                                    pVals[i as usize].data.dbl =
                                        *((*theParams[i as usize]).value.data.dblptr)
                                            .offset(elem as isize);
                                    pNull[i as usize] = *((*theParams[i as usize]).value.undef)
                                        .offset(elem as isize);
                                } else if vector[i as usize] != 0 {
                                    pVals[i as usize].data.dbl =
                                        *((*theParams[i as usize]).value.data.dblptr)
                                            .offset(row as isize);
                                    pNull[i as usize] = *((*theParams[i as usize]).value.undef)
                                        .offset(row as isize);
                                }
                            }
                            let fresh142 = &mut *((*this).value.undef).offset(elem as isize);
                            *fresh142 = (pNull[0] as c_int != 0
                                || pNull[1] as c_int != 0
                                || pNull[2] as c_int != 0
                                || pNull[3] as c_int != 0)
                                as c_int as c_char;
                            if *fresh142 == 0 {
                                *((*this).value.data.dblptr).offset(elem as isize) = angsep_calc(
                                    pVals[0].data.dbl,
                                    pVals[1].data.dbl,
                                    pVals[2].data.dbl,
                                    pVals[3].data.dbl,
                                );
                            }
                        }
                    },
                    1022 => {
                        elem = row * (*theParams[0]).value.nelem;
                        if (*this).ntype == fits_parser_yytokentype::LONG as c_int {
                            let mut minVal: c_long = 0;
                            loop {
                                let fresh143 = row;
                                row -= 1;
                                if fresh143 == 0 {
                                    break;
                                }
                                valInit = 1;
                                *((*this).value.undef).offset(row as isize) = 1;
                                nelem = (*theParams[0]).value.nelem;
                                loop {
                                    let fresh144 = nelem;
                                    nelem -= 1;
                                    if fresh144 == 0 {
                                        break;
                                    }
                                    elem -= 1;
                                    elem;
                                    if *((*theParams[0]).value.undef).offset(elem as isize) == 0 {
                                        if valInit != 0 {
                                            valInit = 0;
                                            minVal = *((*theParams[0]).value.data.lngptr)
                                                .offset(elem as isize);
                                        } else {
                                            minVal = if minVal
                                                < *((*theParams[0]).value.data.lngptr)
                                                    .offset(elem as isize)
                                            {
                                                minVal
                                            } else {
                                                *((*theParams[0]).value.data.lngptr)
                                                    .offset(elem as isize)
                                            };
                                        }
                                        *((*this).value.undef).offset(row as isize) = 0;
                                    }
                                }
                                *((*this).value.data.lngptr).offset(row as isize) = minVal;
                            }
                        } else if (*this).ntype == fits_parser_yytokentype::DOUBLE as c_int {
                            let mut minVal_0: c_double = 0.0f64;
                            loop {
                                let fresh145 = row;
                                row -= 1;
                                if fresh145 == 0 {
                                    break;
                                }
                                valInit = 1;
                                *((*this).value.undef).offset(row as isize) = 1;
                                nelem = (*theParams[0]).value.nelem;
                                loop {
                                    let fresh146 = nelem;
                                    nelem -= 1;
                                    if fresh146 == 0 {
                                        break;
                                    }
                                    elem -= 1;
                                    elem;
                                    if *((*theParams[0]).value.undef).offset(elem as isize) == 0 {
                                        if valInit != 0 {
                                            valInit = 0;
                                            minVal_0 = *((*theParams[0]).value.data.dblptr)
                                                .offset(elem as isize);
                                        } else {
                                            minVal_0 = if minVal_0
                                                < *((*theParams[0]).value.data.dblptr)
                                                    .offset(elem as isize)
                                            {
                                                minVal_0
                                            } else {
                                                *((*theParams[0]).value.data.dblptr)
                                                    .offset(elem as isize)
                                            };
                                        }
                                        *((*this).value.undef).offset(row as isize) = 0;
                                    }
                                }
                                *((*this).value.data.dblptr).offset(row as isize) = minVal_0;
                            }
                        } else if (*this).ntype == fits_parser_yytokentype::BITSTR as c_int {
                            let mut minVal_1: c_char = 0;
                            loop {
                                let fresh147 = row;
                                row -= 1;
                                if fresh147 == 0 {
                                    break;
                                }
                                let mut sptr1_0: *mut c_char =
                                    *((*theParams[0]).value.data.strptr).offset(row as isize);
                                minVal_1 = b'1' as c_char;
                                while *sptr1_0 != 0 {
                                    if *sptr1_0 as c_int == '0' as i32 {
                                        minVal_1 = b'0' as c_char;
                                    }
                                    sptr1_0 = sptr1_0.offset(1);
                                    sptr1_0;
                                }
                                *(*((*this).value.data.strptr).offset(row as isize)).offset(0) =
                                    minVal_1;
                                *(*((*this).value.data.strptr).offset(row as isize))
                                    .offset(1 as c_int as isize) = 0;
                            }
                        }
                    }
                    1023 => {
                        if (*this).ntype == fits_parser_yytokentype::LONG as c_int {
                            loop {
                                let fresh148 = row;
                                row -= 1;
                                if fresh148 == 0 {
                                    break;
                                }
                                nelem = (*this).value.nelem;
                                loop {
                                    let fresh149 = nelem;
                                    nelem -= 1;
                                    if fresh149 == 0 {
                                        break;
                                    }
                                    elem -= 1;
                                    elem;
                                    i = 2 as c_int;
                                    loop {
                                        let fresh150 = i;
                                        i -= 1;
                                        if fresh150 == 0 {
                                            break;
                                        }
                                        if vector[i as usize] > 1 {
                                            pVals[i as usize].data.lng =
                                                *((*theParams[i as usize]).value.data.lngptr)
                                                    .offset(elem as isize);
                                            pNull[i as usize] =
                                                *((*theParams[i as usize]).value.undef)
                                                    .offset(elem as isize);
                                        } else if vector[i as usize] != 0 {
                                            pVals[i as usize].data.lng =
                                                *((*theParams[i as usize]).value.data.lngptr)
                                                    .offset(row as isize);
                                            pNull[i as usize] =
                                                *((*theParams[i as usize]).value.undef)
                                                    .offset(row as isize);
                                        }
                                    }
                                    if pNull[0] as c_int != 0 && pNull[1] as c_int != 0 {
                                        *((*this).value.undef).offset(elem as isize) = 1;
                                        *((*this).value.data.lngptr).offset(elem as isize) = 0;
                                    } else if pNull[0] != 0 {
                                        *((*this).value.undef).offset(elem as isize) = 0;
                                        *((*this).value.data.lngptr).offset(elem as isize) =
                                            pVals[1].data.lng;
                                    } else if pNull[1] != 0 {
                                        *((*this).value.undef).offset(elem as isize) = 0;
                                        *((*this).value.data.lngptr).offset(elem as isize) =
                                            pVals[0].data.lng;
                                    } else {
                                        *((*this).value.undef).offset(elem as isize) = 0;
                                        *((*this).value.data.lngptr).offset(elem as isize) =
                                            if pVals[0].data.lng < pVals[1].data.lng {
                                                pVals[0].data.lng
                                            } else {
                                                pVals[1].data.lng
                                            };
                                    }
                                }
                            }
                        } else if (*this).ntype == fits_parser_yytokentype::DOUBLE as c_int {
                            loop {
                                let fresh151 = row;
                                row -= 1;
                                if fresh151 == 0 {
                                    break;
                                }
                                nelem = (*this).value.nelem;
                                loop {
                                    let fresh152 = nelem;
                                    nelem -= 1;
                                    if fresh152 == 0 {
                                        break;
                                    }
                                    elem -= 1;
                                    elem;
                                    i = 2 as c_int;
                                    loop {
                                        let fresh153 = i;
                                        i -= 1;
                                        if fresh153 == 0 {
                                            break;
                                        }
                                        if vector[i as usize] > 1 {
                                            pVals[i as usize].data.dbl =
                                                *((*theParams[i as usize]).value.data.dblptr)
                                                    .offset(elem as isize);
                                            pNull[i as usize] =
                                                *((*theParams[i as usize]).value.undef)
                                                    .offset(elem as isize);
                                        } else if vector[i as usize] != 0 {
                                            pVals[i as usize].data.dbl =
                                                *((*theParams[i as usize]).value.data.dblptr)
                                                    .offset(row as isize);
                                            pNull[i as usize] =
                                                *((*theParams[i as usize]).value.undef)
                                                    .offset(row as isize);
                                        }
                                    }
                                    if pNull[0] as c_int != 0 && pNull[1] as c_int != 0 {
                                        *((*this).value.undef).offset(elem as isize) = 1;
                                        *((*this).value.data.dblptr).offset(elem as isize) = 0.0;
                                    } else if pNull[0] != 0 {
                                        *((*this).value.undef).offset(elem as isize) = 0;
                                        *((*this).value.data.dblptr).offset(elem as isize) =
                                            pVals[1].data.dbl;
                                    } else if pNull[1] != 0 {
                                        *((*this).value.undef).offset(elem as isize) = 0;
                                        *((*this).value.data.dblptr).offset(elem as isize) =
                                            pVals[0].data.dbl;
                                    } else {
                                        *((*this).value.undef).offset(elem as isize) = 0;
                                        *((*this).value.data.dblptr).offset(elem as isize) =
                                            if pVals[0].data.dbl < pVals[1].data.dbl {
                                                pVals[0].data.dbl
                                            } else {
                                                pVals[1].data.dbl
                                            };
                                    }
                                }
                            }
                        }
                    }
                    1024 => {
                        elem = row * (*theParams[0]).value.nelem;
                        if (*this).ntype == fits_parser_yytokentype::LONG as c_int {
                            let mut maxVal: c_long = 0;
                            loop {
                                let fresh154 = row;
                                row -= 1;
                                if fresh154 == 0 {
                                    break;
                                }
                                valInit = 1;
                                *((*this).value.undef).offset(row as isize) = 1;
                                nelem = (*theParams[0]).value.nelem;
                                loop {
                                    let fresh155 = nelem;
                                    nelem -= 1;
                                    if fresh155 == 0 {
                                        break;
                                    }
                                    elem -= 1;
                                    elem;
                                    if *((*theParams[0]).value.undef).offset(elem as isize) == 0 {
                                        if valInit != 0 {
                                            valInit = 0;
                                            maxVal = *((*theParams[0]).value.data.lngptr)
                                                .offset(elem as isize);
                                        } else {
                                            maxVal = if maxVal
                                                > *((*theParams[0]).value.data.lngptr)
                                                    .offset(elem as isize)
                                            {
                                                maxVal
                                            } else {
                                                *((*theParams[0]).value.data.lngptr)
                                                    .offset(elem as isize)
                                            };
                                        }
                                        *((*this).value.undef).offset(row as isize) = 0;
                                    }
                                }
                                *((*this).value.data.lngptr).offset(row as isize) = maxVal;
                            }
                        } else if (*this).ntype == fits_parser_yytokentype::DOUBLE as c_int {
                            let mut maxVal_0: c_double = 0.0f64;
                            loop {
                                let fresh156 = row;
                                row -= 1;
                                if fresh156 == 0 {
                                    break;
                                }
                                valInit = 1;
                                *((*this).value.undef).offset(row as isize) = 1;
                                nelem = (*theParams[0]).value.nelem;
                                loop {
                                    let fresh157 = nelem;
                                    nelem -= 1;
                                    if fresh157 == 0 {
                                        break;
                                    }
                                    elem -= 1;
                                    elem;
                                    if *((*theParams[0]).value.undef).offset(elem as isize) == 0 {
                                        if valInit != 0 {
                                            valInit = 0;
                                            maxVal_0 = *((*theParams[0]).value.data.dblptr)
                                                .offset(elem as isize);
                                        } else {
                                            maxVal_0 = if maxVal_0
                                                > *((*theParams[0]).value.data.dblptr)
                                                    .offset(elem as isize)
                                            {
                                                maxVal_0
                                            } else {
                                                *((*theParams[0]).value.data.dblptr)
                                                    .offset(elem as isize)
                                            };
                                        }
                                        *((*this).value.undef).offset(row as isize) = 0;
                                    }
                                }
                                *((*this).value.data.dblptr).offset(row as isize) = maxVal_0;
                            }
                        } else if (*this).ntype == fits_parser_yytokentype::BITSTR as c_int {
                            let mut maxVal_1: c_char = 0;
                            loop {
                                let fresh158 = row;
                                row -= 1;
                                if fresh158 == 0 {
                                    break;
                                }
                                let mut sptr1_1: *mut c_char =
                                    *((*theParams[0]).value.data.strptr).offset(row as isize);
                                maxVal_1 = b'0' as c_char;
                                while *sptr1_1 != 0 {
                                    if *sptr1_1 as c_int == '1' as i32 {
                                        maxVal_1 = b'1' as c_char;
                                    }
                                    sptr1_1 = sptr1_1.offset(1);
                                    sptr1_1;
                                }
                                *(*((*this).value.data.strptr).offset(row as isize)).offset(0) =
                                    maxVal_1;
                                *(*((*this).value.data.strptr).offset(row as isize))
                                    .offset(1 as c_int as isize) = 0;
                            }
                        }
                    }
                    1025 => {
                        if (*this).ntype == fits_parser_yytokentype::LONG as c_int {
                            loop {
                                let fresh159 = row;
                                row -= 1;
                                if fresh159 == 0 {
                                    break;
                                }
                                nelem = (*this).value.nelem;
                                loop {
                                    let fresh160 = nelem;
                                    nelem -= 1;
                                    if fresh160 == 0 {
                                        break;
                                    }
                                    elem -= 1;
                                    elem;
                                    i = 2 as c_int;
                                    loop {
                                        let fresh161 = i;
                                        i -= 1;
                                        if fresh161 == 0 {
                                            break;
                                        }
                                        if vector[i as usize] > 1 {
                                            pVals[i as usize].data.lng =
                                                *((*theParams[i as usize]).value.data.lngptr)
                                                    .offset(elem as isize);
                                            pNull[i as usize] =
                                                *((*theParams[i as usize]).value.undef)
                                                    .offset(elem as isize);
                                        } else if vector[i as usize] != 0 {
                                            pVals[i as usize].data.lng =
                                                *((*theParams[i as usize]).value.data.lngptr)
                                                    .offset(row as isize);
                                            pNull[i as usize] =
                                                *((*theParams[i as usize]).value.undef)
                                                    .offset(row as isize);
                                        }
                                    }
                                    if pNull[0] as c_int != 0 && pNull[1] as c_int != 0 {
                                        *((*this).value.undef).offset(elem as isize) = 1;
                                        *((*this).value.data.lngptr).offset(elem as isize) = 0;
                                    } else if pNull[0] != 0 {
                                        *((*this).value.undef).offset(elem as isize) = 0;
                                        *((*this).value.data.lngptr).offset(elem as isize) =
                                            pVals[1].data.lng;
                                    } else if pNull[1] != 0 {
                                        *((*this).value.undef).offset(elem as isize) = 0;
                                        *((*this).value.data.lngptr).offset(elem as isize) =
                                            pVals[0].data.lng;
                                    } else {
                                        *((*this).value.undef).offset(elem as isize) = 0;
                                        *((*this).value.data.lngptr).offset(elem as isize) =
                                            if pVals[0].data.lng > pVals[1].data.lng {
                                                pVals[0].data.lng
                                            } else {
                                                pVals[1].data.lng
                                            };
                                    }
                                }
                            }
                        } else if (*this).ntype == fits_parser_yytokentype::DOUBLE as c_int {
                            loop {
                                let fresh162 = row;
                                row -= 1;
                                if fresh162 == 0 {
                                    break;
                                }
                                nelem = (*this).value.nelem;
                                loop {
                                    let fresh163 = nelem;
                                    nelem -= 1;
                                    if fresh163 == 0 {
                                        break;
                                    }
                                    elem -= 1;
                                    elem;
                                    i = 2 as c_int;
                                    loop {
                                        let fresh164 = i;
                                        i -= 1;
                                        if fresh164 == 0 {
                                            break;
                                        }
                                        if vector[i as usize] > 1 {
                                            pVals[i as usize].data.dbl =
                                                *((*theParams[i as usize]).value.data.dblptr)
                                                    .offset(elem as isize);
                                            pNull[i as usize] =
                                                *((*theParams[i as usize]).value.undef)
                                                    .offset(elem as isize);
                                        } else if vector[i as usize] != 0 {
                                            pVals[i as usize].data.dbl =
                                                *((*theParams[i as usize]).value.data.dblptr)
                                                    .offset(row as isize);
                                            pNull[i as usize] =
                                                *((*theParams[i as usize]).value.undef)
                                                    .offset(row as isize);
                                        }
                                    }
                                    if pNull[0] as c_int != 0 && pNull[1] as c_int != 0 {
                                        *((*this).value.undef).offset(elem as isize) = 1;
                                        *((*this).value.data.dblptr).offset(elem as isize) = 0.0;
                                    } else if pNull[0] != 0 {
                                        *((*this).value.undef).offset(elem as isize) = 0;
                                        *((*this).value.data.dblptr).offset(elem as isize) =
                                            pVals[1].data.dbl;
                                    } else if pNull[1] != 0 {
                                        *((*this).value.undef).offset(elem as isize) = 0;
                                        *((*this).value.data.dblptr).offset(elem as isize) =
                                            pVals[0].data.dbl;
                                    } else {
                                        *((*this).value.undef).offset(elem as isize) = 0;
                                        *((*this).value.data.dblptr).offset(elem as isize) =
                                            if pVals[0].data.dbl > pVals[1].data.dbl {
                                                pVals[0].data.dbl
                                            } else {
                                                pVals[1].data.dbl
                                            };
                                    }
                                }
                            }
                        }
                    }
                    1026 => loop {
                        let fresh165 = row;
                        row -= 1;
                        if fresh165 == 0 {
                            break;
                        }
                        nelem = (*this).value.nelem;
                        loop {
                            let fresh166 = nelem;
                            nelem -= 1;
                            if fresh166 == 0 {
                                break;
                            }
                            elem -= 1;
                            elem;
                            i = 3 as c_int;
                            loop {
                                let fresh167 = i;
                                i -= 1;
                                if fresh167 == 0 {
                                    break;
                                }
                                if vector[i as usize] > 1 {
                                    pVals[i as usize].data.dbl =
                                        *((*theParams[i as usize]).value.data.dblptr)
                                            .offset(elem as isize);
                                    pNull[i as usize] = *((*theParams[i as usize]).value.undef)
                                        .offset(elem as isize);
                                } else if vector[i as usize] != 0 {
                                    pVals[i as usize].data.dbl =
                                        *((*theParams[i as usize]).value.data.dblptr)
                                            .offset(row as isize);
                                    pNull[i as usize] = *((*theParams[i as usize]).value.undef)
                                        .offset(row as isize);
                                }
                            }
                            let fresh168 = &mut *((*this).value.undef).offset(elem as isize);
                            *fresh168 = (pNull[0] as c_int != 0
                                || pNull[1] as c_int != 0
                                || pNull[2] as c_int != 0)
                                as c_int as c_char;
                            if *fresh168 == 0 {
                                *((*this).value.data.logptr).offset(elem as isize) =
                                    bnear(pVals[0].data.dbl, pVals[1].data.dbl, pVals[2].data.dbl);
                            }
                        }
                    },
                    1027 => loop {
                        let fresh169 = row;
                        row -= 1;
                        if fresh169 == 0 {
                            break;
                        }
                        nelem = (*this).value.nelem;
                        loop {
                            let fresh170 = nelem;
                            nelem -= 1;
                            if fresh170 == 0 {
                                break;
                            }
                            elem -= 1;
                            elem;
                            i = 5 as c_int;
                            loop {
                                let fresh171 = i;
                                i -= 1;
                                if fresh171 == 0 {
                                    break;
                                }
                                if vector[i as usize] > 1 {
                                    pVals[i as usize].data.dbl =
                                        *((*theParams[i as usize]).value.data.dblptr)
                                            .offset(elem as isize);
                                    pNull[i as usize] = *((*theParams[i as usize]).value.undef)
                                        .offset(elem as isize);
                                } else if vector[i as usize] != 0 {
                                    pVals[i as usize].data.dbl =
                                        *((*theParams[i as usize]).value.data.dblptr)
                                            .offset(row as isize);
                                    pNull[i as usize] = *((*theParams[i as usize]).value.undef)
                                        .offset(row as isize);
                                }
                            }
                            let fresh172 = &mut *((*this).value.undef).offset(elem as isize);
                            *fresh172 = (pNull[0] as c_int != 0
                                || pNull[1] as c_int != 0
                                || pNull[2] as c_int != 0
                                || pNull[3] as c_int != 0
                                || pNull[4] as c_int != 0)
                                as c_int as c_char;
                            if *fresh172 == 0 {
                                *((*this).value.data.logptr).offset(elem as isize) = circle(
                                    pVals[0].data.dbl,
                                    pVals[1].data.dbl,
                                    pVals[2].data.dbl,
                                    pVals[3].data.dbl,
                                    pVals[4].data.dbl,
                                );
                            }
                        }
                    },
                    1028 => loop {
                        let fresh173 = row;
                        row -= 1;
                        if fresh173 == 0 {
                            break;
                        }
                        nelem = (*this).value.nelem;
                        loop {
                            let fresh174 = nelem;
                            nelem -= 1;
                            if fresh174 == 0 {
                                break;
                            }
                            elem -= 1;
                            elem;
                            i = 7 as c_int;
                            loop {
                                let fresh175 = i;
                                i -= 1;
                                if fresh175 == 0 {
                                    break;
                                }
                                if vector[i as usize] > 1 {
                                    pVals[i as usize].data.dbl =
                                        *((*theParams[i as usize]).value.data.dblptr)
                                            .offset(elem as isize);
                                    pNull[i as usize] = *((*theParams[i as usize]).value.undef)
                                        .offset(elem as isize);
                                } else if vector[i as usize] != 0 {
                                    pVals[i as usize].data.dbl =
                                        *((*theParams[i as usize]).value.data.dblptr)
                                            .offset(row as isize);
                                    pNull[i as usize] = *((*theParams[i as usize]).value.undef)
                                        .offset(row as isize);
                                }
                            }
                            let fresh176 = &mut *((*this).value.undef).offset(elem as isize);
                            *fresh176 = (pNull[0] as c_int != 0
                                || pNull[1] as c_int != 0
                                || pNull[2] as c_int != 0
                                || pNull[3] as c_int != 0
                                || pNull[4] as c_int != 0
                                || pNull[5] as c_int != 0
                                || pNull[6] as c_int != 0)
                                as c_int as c_char;
                            if *fresh176 == 0 {
                                *((*this).value.data.logptr).offset(elem as isize) = saobox(
                                    pVals[0].data.dbl,
                                    pVals[1].data.dbl,
                                    pVals[2].data.dbl,
                                    pVals[3].data.dbl,
                                    pVals[4].data.dbl,
                                    pVals[5].data.dbl,
                                    pVals[6].data.dbl,
                                );
                            }
                        }
                    },
                    1029 => loop {
                        let fresh177 = row;
                        row -= 1;
                        if fresh177 == 0 {
                            break;
                        }
                        nelem = (*this).value.nelem;
                        loop {
                            let fresh178 = nelem;
                            nelem -= 1;
                            if fresh178 == 0 {
                                break;
                            }
                            elem -= 1;
                            elem;
                            i = 7 as c_int;
                            loop {
                                let fresh179 = i;
                                i -= 1;
                                if fresh179 == 0 {
                                    break;
                                }
                                if vector[i as usize] > 1 {
                                    pVals[i as usize].data.dbl =
                                        *((*theParams[i as usize]).value.data.dblptr)
                                            .offset(elem as isize);
                                    pNull[i as usize] = *((*theParams[i as usize]).value.undef)
                                        .offset(elem as isize);
                                } else if vector[i as usize] != 0 {
                                    pVals[i as usize].data.dbl =
                                        *((*theParams[i as usize]).value.data.dblptr)
                                            .offset(row as isize);
                                    pNull[i as usize] = *((*theParams[i as usize]).value.undef)
                                        .offset(row as isize);
                                }
                            }
                            let fresh180 = &mut *((*this).value.undef).offset(elem as isize);
                            *fresh180 = (pNull[0] as c_int != 0
                                || pNull[1] as c_int != 0
                                || pNull[2] as c_int != 0
                                || pNull[3] as c_int != 0
                                || pNull[4] as c_int != 0
                                || pNull[5] as c_int != 0
                                || pNull[6] as c_int != 0)
                                as c_int as c_char;
                            if *fresh180 == 0 {
                                *((*this).value.data.logptr).offset(elem as isize) = ellipse(
                                    pVals[0].data.dbl,
                                    pVals[1].data.dbl,
                                    pVals[2].data.dbl,
                                    pVals[3].data.dbl,
                                    pVals[4].data.dbl,
                                    pVals[5].data.dbl,
                                    pVals[6].data.dbl,
                                );
                            }
                        }
                    },
                    1034 => match (*this).ntype {
                        258 => loop {
                            let fresh181 = row;
                            row -= 1;
                            if fresh181 == 0 {
                                break;
                            }
                            nelem = (*this).value.nelem;
                            loop {
                                let fresh182 = nelem;
                                nelem -= 1;
                                if fresh182 == 0 {
                                    break;
                                }
                                elem -= 1;
                                elem;
                                if vector[2] > 1 {
                                    pVals[2].data.log =
                                        *((*theParams[2]).value.data.logptr).offset(elem as isize);
                                    pNull[2] = *((*theParams[2]).value.undef).offset(elem as isize);
                                } else if vector[2] != 0 {
                                    pVals[2].data.log =
                                        *((*theParams[2]).value.data.logptr).offset(row as isize);
                                    pNull[2] = *((*theParams[2]).value.undef).offset(row as isize);
                                }
                                i = 2 as c_int;
                                loop {
                                    let fresh183 = i;
                                    i -= 1;
                                    if fresh183 == 0 {
                                        break;
                                    }
                                    if vector[i as usize] > 1 {
                                        pVals[i as usize].data.log =
                                            *((*theParams[i as usize]).value.data.logptr)
                                                .offset(elem as isize);
                                        pNull[i as usize] = *((*theParams[i as usize]).value.undef)
                                            .offset(elem as isize);
                                    } else if vector[i as usize] != 0 {
                                        pVals[i as usize].data.log =
                                            *((*theParams[i as usize]).value.data.logptr)
                                                .offset(row as isize);
                                        pNull[i as usize] = *((*theParams[i as usize]).value.undef)
                                            .offset(row as isize);
                                    }
                                }
                                let fresh184 = &mut *((*this).value.undef).offset(elem as isize);
                                *fresh184 = pNull[2];
                                if *fresh184 == 0 {
                                    if pVals[2].data.log != 0 {
                                        *((*this).value.data.logptr).offset(elem as isize) =
                                            pVals[0].data.log;
                                        *((*this).value.undef).offset(elem as isize) = pNull[0];
                                    } else {
                                        *((*this).value.data.logptr).offset(elem as isize) =
                                            pVals[1].data.log;
                                        *((*this).value.undef).offset(elem as isize) = pNull[1];
                                    }
                                }
                            }
                        },
                        259 => loop {
                            let fresh185 = row;
                            row -= 1;
                            if fresh185 == 0 {
                                break;
                            }
                            nelem = (*this).value.nelem;
                            loop {
                                let fresh186 = nelem;
                                nelem -= 1;
                                if fresh186 == 0 {
                                    break;
                                }
                                elem -= 1;
                                elem;
                                if vector[2] > 1 {
                                    pVals[2].data.log =
                                        *((*theParams[2]).value.data.logptr).offset(elem as isize);
                                    pNull[2] = *((*theParams[2]).value.undef).offset(elem as isize);
                                } else if vector[2] != 0 {
                                    pVals[2].data.log =
                                        *((*theParams[2]).value.data.logptr).offset(row as isize);
                                    pNull[2] = *((*theParams[2]).value.undef).offset(row as isize);
                                }
                                i = 2 as c_int;
                                loop {
                                    let fresh187 = i;
                                    i -= 1;
                                    if fresh187 == 0 {
                                        break;
                                    }
                                    if vector[i as usize] > 1 {
                                        pVals[i as usize].data.lng =
                                            *((*theParams[i as usize]).value.data.lngptr)
                                                .offset(elem as isize);
                                        pNull[i as usize] = *((*theParams[i as usize]).value.undef)
                                            .offset(elem as isize);
                                    } else if vector[i as usize] != 0 {
                                        pVals[i as usize].data.lng =
                                            *((*theParams[i as usize]).value.data.lngptr)
                                                .offset(row as isize);
                                        pNull[i as usize] = *((*theParams[i as usize]).value.undef)
                                            .offset(row as isize);
                                    }
                                }
                                let fresh188 = &mut *((*this).value.undef).offset(elem as isize);
                                *fresh188 = pNull[2];
                                if *fresh188 == 0 {
                                    if pVals[2].data.log != 0 {
                                        *((*this).value.data.lngptr).offset(elem as isize) =
                                            pVals[0].data.lng;
                                        *((*this).value.undef).offset(elem as isize) = pNull[0];
                                    } else {
                                        *((*this).value.data.lngptr).offset(elem as isize) =
                                            pVals[1].data.lng;
                                        *((*this).value.undef).offset(elem as isize) = pNull[1];
                                    }
                                }
                            }
                        },
                        260 => loop {
                            let fresh189 = row;
                            row -= 1;
                            if fresh189 == 0 {
                                break;
                            }
                            nelem = (*this).value.nelem;
                            loop {
                                let fresh190 = nelem;
                                nelem -= 1;
                                if fresh190 == 0 {
                                    break;
                                }
                                elem -= 1;
                                elem;
                                if vector[2] > 1 {
                                    pVals[2].data.log =
                                        *((*theParams[2]).value.data.logptr).offset(elem as isize);
                                    pNull[2] = *((*theParams[2]).value.undef).offset(elem as isize);
                                } else if vector[2] != 0 {
                                    pVals[2].data.log =
                                        *((*theParams[2]).value.data.logptr).offset(row as isize);
                                    pNull[2] = *((*theParams[2]).value.undef).offset(row as isize);
                                }
                                i = 2 as c_int;
                                loop {
                                    let fresh191 = i;
                                    i -= 1;
                                    if fresh191 == 0 {
                                        break;
                                    }
                                    if vector[i as usize] > 1 {
                                        pVals[i as usize].data.dbl =
                                            *((*theParams[i as usize]).value.data.dblptr)
                                                .offset(elem as isize);
                                        pNull[i as usize] = *((*theParams[i as usize]).value.undef)
                                            .offset(elem as isize);
                                    } else if vector[i as usize] != 0 {
                                        pVals[i as usize].data.dbl =
                                            *((*theParams[i as usize]).value.data.dblptr)
                                                .offset(row as isize);
                                        pNull[i as usize] = *((*theParams[i as usize]).value.undef)
                                            .offset(row as isize);
                                    }
                                }
                                let fresh192 = &mut *((*this).value.undef).offset(elem as isize);
                                *fresh192 = pNull[2];
                                if *fresh192 == 0 {
                                    if pVals[2].data.log != 0 {
                                        *((*this).value.data.dblptr).offset(elem as isize) =
                                            pVals[0].data.dbl;
                                        *((*this).value.undef).offset(elem as isize) = pNull[0];
                                    } else {
                                        *((*this).value.data.dblptr).offset(elem as isize) =
                                            pVals[1].data.dbl;
                                        *((*this).value.undef).offset(elem as isize) = pNull[1];
                                    }
                                }
                            }
                        },
                        261 => loop {
                            let fresh193 = row;
                            row -= 1;
                            if fresh193 == 0 {
                                break;
                            }
                            if vector[2] != 0 {
                                pVals[2].data.log =
                                    *((*theParams[2]).value.data.logptr).offset(row as isize);
                                pNull[2] = *((*theParams[2]).value.undef).offset(row as isize);
                            }
                            i = 2 as c_int;
                            loop {
                                let fresh194 = i;
                                i -= 1;
                                if fresh194 == 0 {
                                    break;
                                }
                                if vector[i as usize] != 0 {
                                    strcpy(
                                        (pVals[i as usize].data.astr).as_mut_ptr(),
                                        *((*theParams[i as usize]).value.data.strptr)
                                            .offset(row as isize),
                                    );
                                    pNull[i as usize] = *((*theParams[i as usize]).value.undef)
                                        .offset(row as isize);
                                }
                            }
                            let fresh195 = &mut *((*this).value.undef).offset(row as isize);
                            *fresh195 = pNull[2];
                            if *fresh195 == 0 {
                                if pVals[2].data.log != 0 {
                                    strcpy(
                                        *((*this).value.data.strptr).offset(row as isize),
                                        (pVals[0].data.astr).as_mut_ptr(),
                                    );
                                    *((*this).value.undef).offset(row as isize) = pNull[0];
                                } else {
                                    strcpy(
                                        *((*this).value.data.strptr).offset(row as isize),
                                        (pVals[1].data.astr).as_mut_ptr(),
                                    );
                                    *((*this).value.undef).offset(row as isize) = pNull[1];
                                }
                            } else {
                                *(*((*this).value.data.strptr).offset(row as isize)).offset(0) = 0;
                            }
                        },
                        _ => {}
                    },
                    1044 => {
                        let mut strconst: c_int =
                            ((*theParams[0]).operation == -(1000 as c_int)) as c_int;
                        let mut posconst: c_int =
                            ((*theParams[1]).operation == -(1000 as c_int)) as c_int;
                        let mut lenconst: c_int =
                            ((*theParams[2]).operation == -(1000 as c_int)) as c_int;
                        let mut dest_len: c_int = (*this).value.nelem as c_int;
                        let mut src_len: c_int = (*theParams[0]).value.nelem as c_int;
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
                                pos = (*theParams[1]).value.data.lng as c_int;
                            } else {
                                pos = *((*theParams[1]).value.data.lngptr).offset(row as isize)
                                    as c_int;
                                if *((*theParams[1]).value.undef).offset(row as isize) != 0 {
                                    undef = 1;
                                }
                            }
                            if strconst != 0 {
                                str = ((*theParams[0]).value.data.astr).as_mut_ptr();
                                if src_len == 0 {
                                    src_len = strlen(str) as c_int;
                                }
                            } else {
                                str = *((*theParams[0]).value.data.strptr).offset(row as isize);
                                if *((*theParams[0]).value.undef).offset(row as isize) != 0 {
                                    undef = 1;
                                }
                            }
                            if lenconst != 0 {
                                len = dest_len;
                            } else {
                                len = *((*theParams[2]).value.data.lngptr).offset(row as isize)
                                    as c_int;
                                if *((*theParams[2]).value.undef).offset(row as isize) != 0 {
                                    undef = 1;
                                }
                            }
                            *(*((*this).value.data.strptr).offset(row as isize)).offset(0) = 0;
                            if pos == 0 {
                                undef = 1;
                            }
                            if undef == 0
                                && cstrmid(
                                    lParse,
                                    *((*this).value.data.strptr).offset(row as isize),
                                    len,
                                    str,
                                    src_len,
                                    pos,
                                ) < 0
                            {
                                break;
                            }
                            *((*this).value.undef).offset(row as isize) = undef as c_char;
                        }
                    }
                    1045 => {
                        let mut const1: c_int =
                            ((*theParams[0]).operation == -(1000 as c_int)) as c_int;
                        let mut const2: c_int =
                            ((*theParams[1]).operation == -(1000 as c_int)) as c_int;
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
                                str1 = ((*theParams[0]).value.data.astr).as_mut_ptr();
                            } else {
                                str1 = *((*theParams[0]).value.data.strptr).offset(row as isize);
                                if *((*theParams[0]).value.undef).offset(row as isize) != 0 {
                                    undef_0 = 1;
                                }
                            }
                            if const2 != 0 {
                                str2 = ((*theParams[1]).value.data.astr).as_mut_ptr();
                            } else {
                                str2 = *((*theParams[1]).value.data.strptr).offset(row as isize);
                                if *((*theParams[1]).value.undef).offset(row as isize) != 0 {
                                    undef_0 = 1;
                                }
                            }
                            *((*this).value.data.lngptr).offset(row as isize) = 0;
                            if undef_0 == 0 {
                                let mut res_0: *mut c_char = strstr(str1, str2);
                                if res_0.is_null() {
                                    undef_0 = 1;
                                    *((*this).value.data.lngptr).offset(row as isize) = 0;
                                } else {
                                    *((*this).value.data.lngptr).offset(row as isize) =
                                        res_0.offset_from(str1) as c_long + 1;
                                }
                            }
                            *((*this).value.undef).offset(row as isize) = undef_0 as c_char;
                        }
                    }
                    _ => {}
                }
            }
        }
        i = (*this).nSubNodes;
        loop {
            let fresh198 = i;
            i -= 1;
            if fresh198 == 0 {
                break;
            }
            if (*theParams[i as usize]).operation > 0 {
                free((*theParams[i as usize]).value.data.ptr);
            }
        }
    }
}
unsafe fn Do_Deref(mut lParse: *mut ParseData, mut this: *mut Node) {
    unsafe {
        let mut theVar: *mut Node = std::ptr::null_mut::<Node>();
        let mut theDims: [*mut Node; 5] = [std::ptr::null_mut::<Node>(); 5];
        let mut isConst: [c_int; 5] = [0; 5];
        let mut allConst: c_int = 0;
        let mut dimVals: [c_long; 5] = [0; 5];
        let mut i: c_int = 0;
        let mut nDims: c_int = 0;
        let mut row: c_long = 0;
        let mut elem: c_long = 0;
        let mut dsize: c_long = 0;
        theVar = ((*lParse).Nodes).offset((*this).SubNodes[0] as isize);
        nDims = (*this).nSubNodes - 1;
        i = nDims;
        allConst = 1;
        loop {
            let fresh199 = i;
            i -= 1;
            if fresh199 == 0 {
                break;
            }
            theDims[i as usize] =
                ((*lParse).Nodes).offset((*this).SubNodes[(i + 1) as usize] as isize);
            isConst[i as usize] = ((*theDims[i as usize]).operation == -(1000 as c_int)) as c_int;
            if isConst[i as usize] != 0 {
                dimVals[i as usize] = (*theDims[i as usize]).value.data.lng;
            } else {
                allConst = 0;
            }
        }
        if (*this).ntype == fits_parser_yytokentype::DOUBLE as c_int {
            dsize = ::core::mem::size_of::<c_double>() as c_ulong as c_long;
        } else if (*this).ntype == fits_parser_yytokentype::LONG as c_int {
            dsize = ::core::mem::size_of::<c_long>() as c_ulong as c_long;
        } else if (*this).ntype == fits_parser_yytokentype::BOOLEAN as c_int {
            dsize = ::core::mem::size_of::<c_char>() as c_ulong as c_long;
        } else {
            dsize = 0;
        }
        Allocate_Ptrs(lParse, this);
        if (*lParse).status == 0 {
            if allConst != 0 && (*theVar).value.naxis == nDims {
                elem = 0;
                i = nDims;
                loop {
                    let fresh200 = i;
                    i -= 1;
                    if fresh200 == 0 {
                        break;
                    }
                    if dimVals[i as usize] < 1
                        || dimVals[i as usize] > (*theVar).value.naxes[i as usize]
                    {
                        break;
                    }
                    elem = (*theVar).value.naxes[i as usize] * elem + dimVals[i as usize] - 1;
                }
                if i < 0 {
                    row = 0;
                    while row < (*lParse).nRows {
                        if (*this).ntype == fits_parser_yytokentype::STRING as c_int {
                            *((*this).value.undef).offset(row as isize) =
                                *((*theVar).value.undef).offset(row as isize);
                        } else if (*this).ntype != fits_parser_yytokentype::BITSTR as c_int {
                            *((*this).value.undef).offset(row as isize) =
                                *((*theVar).value.undef).offset(elem as isize);
                        }
                        if (*this).ntype == fits_parser_yytokentype::DOUBLE as c_int {
                            *((*this).value.data.dblptr).offset(row as isize) =
                                *((*theVar).value.data.dblptr).offset(elem as isize);
                        } else if (*this).ntype == fits_parser_yytokentype::LONG as c_int {
                            *((*this).value.data.lngptr).offset(row as isize) =
                                *((*theVar).value.data.lngptr).offset(elem as isize);
                        } else if (*this).ntype == fits_parser_yytokentype::BOOLEAN as c_int {
                            *((*this).value.data.logptr).offset(row as isize) =
                                *((*theVar).value.data.logptr).offset(elem as isize);
                        } else {
                            *(*((*this).value.data.strptr).offset(row as isize)).offset(0) =
                                *(*((*theVar).value.data.strptr).offset(0))
                                    .offset((elem + row) as isize);
                            *(*((*this).value.data.strptr).offset(row as isize))
                                .offset(1 as c_int as isize) = 0;
                        }
                        elem += (*theVar).value.nelem;
                        row += 1;
                        row;
                    }
                } else {
                    fits_parser_yyerror(
                        0 as yyscan_t,
                        lParse,
                        b"Index out of range\0" as *const u8 as *const c_char as *mut c_char,
                    );
                    free((*this).value.data.ptr);
                }
            } else if allConst != 0 && nDims == 1 {
                if dimVals[0] < 1
                    || dimVals[0] > (*theVar).value.naxes[((*theVar).value.naxis - 1) as usize]
                {
                    fits_parser_yyerror(
                        0 as yyscan_t,
                        lParse,
                        b"Index out of range\0" as *const u8 as *const c_char as *mut c_char,
                    );
                    free((*this).value.data.ptr);
                } else if (*this).ntype == fits_parser_yytokentype::BITSTR as c_int
                    || (*this).ntype == fits_parser_yytokentype::STRING as c_int
                {
                    elem = (*this).value.nelem * (dimVals[0] - 1);
                    row = 0;
                    while row < (*lParse).nRows {
                        if !((*this).value.undef).is_null() {
                            *((*this).value.undef).offset(row as isize) =
                                *((*theVar).value.undef).offset(row as isize);
                        }
                        memcpy(
                            (*((*this).value.data.strptr).offset(0)).offset(
                                (row as c_ulong)
                                    .wrapping_mul(::core::mem::size_of::<c_char>() as c_ulong)
                                    .wrapping_mul(((*this).value.nelem + 1) as c_ulong)
                                    as isize,
                            ) as *mut c_void,
                            (*((*theVar).value.data.strptr).offset(0)).offset(
                                (elem as c_ulong)
                                    .wrapping_mul(::core::mem::size_of::<c_char>() as c_ulong)
                                    as isize,
                            ) as *const c_void,
                            ((*this).value.nelem as c_ulong)
                                .wrapping_mul(::core::mem::size_of::<c_char>() as c_ulong)
                                .try_into()
                                .unwrap(),
                        );
                        *(*((*this).value.data.strptr).offset(row as isize))
                            .offset((*this).value.nelem as isize) = 0;
                        elem += (*theVar).value.nelem + 1;
                        row += 1;
                        row;
                    }
                } else {
                    elem = (*this).value.nelem * (dimVals[0] - 1);
                    row = 0;
                    while row < (*lParse).nRows {
                        memcpy(
                            ((*this).value.undef).offset((row * (*this).value.nelem) as isize)
                                as *mut c_void,
                            ((*theVar).value.undef).offset(elem as isize) as *const c_void,
                            ((*this).value.nelem as c_ulong)
                                .wrapping_mul(::core::mem::size_of::<c_char>() as c_ulong)
                                .try_into()
                                .unwrap(),
                        );
                        memcpy(
                            ((*this).value.data.ptr as *mut c_char)
                                .offset((row * dsize * (*this).value.nelem) as isize)
                                as *mut c_void,
                            ((*theVar).value.data.ptr as *mut c_char)
                                .offset((elem * dsize) as isize)
                                as *const c_void,
                            (((*this).value.nelem * dsize) as c_ulong)
                                .try_into()
                                .unwrap(),
                        );
                        elem += (*theVar).value.nelem;
                        row += 1;
                        row;
                    }
                }
            } else if (*theVar).value.naxis == nDims {
                row = 0;
                while row < (*lParse).nRows {
                    i = 0;
                    while i < nDims {
                        if isConst[i as usize] == 0 {
                            if *((*theDims[i as usize]).value.undef).offset(row as isize) != 0 {
                                fits_parser_yyerror(
                                    0 as yyscan_t,
                                    lParse,
                                    b"Null encountered as vector index\0" as *const u8
                                        as *const c_char
                                        as *mut c_char,
                                );
                                free((*this).value.data.ptr);
                                break;
                            } else {
                                dimVals[i as usize] = *((*theDims[i as usize]).value.data.lngptr)
                                    .offset(row as isize);
                            }
                        }
                        i += 1;
                        i;
                    }
                    if (*lParse).status != 0 {
                        break;
                    }
                    elem = 0;
                    i = nDims;
                    loop {
                        let fresh201 = i;
                        i -= 1;
                        if fresh201 == 0 {
                            break;
                        }
                        if dimVals[i as usize] < 1
                            || dimVals[i as usize] > (*theVar).value.naxes[i as usize]
                        {
                            break;
                        }
                        elem = (*theVar).value.naxes[i as usize] * elem + dimVals[i as usize] - 1;
                    }
                    if i < 0 {
                        elem += row * (*theVar).value.nelem;
                        if (*this).ntype == fits_parser_yytokentype::STRING as c_int {
                            *((*this).value.undef).offset(row as isize) =
                                *((*theVar).value.undef).offset(row as isize);
                        } else if (*this).ntype != fits_parser_yytokentype::BITSTR as c_int {
                            *((*this).value.undef).offset(row as isize) =
                                *((*theVar).value.undef).offset(elem as isize);
                        }
                        if (*this).ntype == fits_parser_yytokentype::DOUBLE as c_int {
                            *((*this).value.data.dblptr).offset(row as isize) =
                                *((*theVar).value.data.dblptr).offset(elem as isize);
                        } else if (*this).ntype == fits_parser_yytokentype::LONG as c_int {
                            *((*this).value.data.lngptr).offset(row as isize) =
                                *((*theVar).value.data.lngptr).offset(elem as isize);
                        } else if (*this).ntype == fits_parser_yytokentype::BOOLEAN as c_int {
                            *((*this).value.data.logptr).offset(row as isize) =
                                *((*theVar).value.data.logptr).offset(elem as isize);
                        } else {
                            *(*((*this).value.data.strptr).offset(row as isize)).offset(0) =
                                *(*((*theVar).value.data.strptr).offset(0))
                                    .offset((elem + row) as isize);
                            *(*((*this).value.data.strptr).offset(row as isize))
                                .offset(1 as c_int as isize) = 0;
                        }
                    } else {
                        fits_parser_yyerror(
                            0 as yyscan_t,
                            lParse,
                            b"Index out of range\0" as *const u8 as *const c_char as *mut c_char,
                        );
                        free((*this).value.data.ptr);
                    }
                    row += 1;
                    row;
                }
            } else {
                row = 0;
                while row < (*lParse).nRows {
                    if *((*theDims[0]).value.undef).offset(row as isize) != 0 {
                        fits_parser_yyerror(
                            0 as yyscan_t,
                            lParse,
                            b"Null encountered as vector index\0" as *const u8 as *const c_char
                                as *mut c_char,
                        );
                        free((*this).value.data.ptr);
                        break;
                    } else {
                        dimVals[0] = *((*theDims[0]).value.data.lngptr).offset(row as isize);
                        if dimVals[0] < 1
                            || dimVals[0]
                                > (*theVar).value.naxes[((*theVar).value.naxis - 1) as usize]
                        {
                            fits_parser_yyerror(
                                0 as yyscan_t,
                                lParse,
                                b"Index out of range\0" as *const u8 as *const c_char
                                    as *mut c_char,
                            );
                            free((*this).value.data.ptr);
                        } else if (*this).ntype == fits_parser_yytokentype::BITSTR as c_int
                            || (*this).ntype == fits_parser_yytokentype::STRING as c_int
                        {
                            elem = (*this).value.nelem * (dimVals[0] - 1);
                            elem += row * ((*theVar).value.nelem + 1);
                            if !((*this).value.undef).is_null() {
                                *((*this).value.undef).offset(row as isize) =
                                    *((*theVar).value.undef).offset(row as isize);
                            }
                            memcpy(
                                (*((*this).value.data.strptr).offset(0)).offset(
                                    (row as c_ulong)
                                        .wrapping_mul(::core::mem::size_of::<c_char>() as c_ulong)
                                        .wrapping_mul(((*this).value.nelem + 1) as c_ulong)
                                        as isize,
                                ) as *mut c_void,
                                (*((*theVar).value.data.strptr).offset(0)).offset(
                                    (elem as c_ulong)
                                        .wrapping_mul(::core::mem::size_of::<c_char>() as c_ulong)
                                        as isize,
                                ) as *const c_void,
                                ((*this).value.nelem as c_ulong)
                                    .wrapping_mul(::core::mem::size_of::<c_char>() as c_ulong)
                                    .try_into()
                                    .unwrap(),
                            );
                            *(*((*this).value.data.strptr).offset(row as isize))
                                .offset((*this).value.nelem as isize) = 0;
                        } else {
                            elem = (*this).value.nelem * (dimVals[0] - 1);
                            elem += row * (*theVar).value.nelem;
                            memcpy(
                                ((*this).value.undef).offset((row * (*this).value.nelem) as isize)
                                    as *mut c_void,
                                ((*theVar).value.undef).offset(elem as isize) as *const c_void,
                                ((*this).value.nelem as c_ulong)
                                    .wrapping_mul(::core::mem::size_of::<c_char>() as c_ulong)
                                    .try_into()
                                    .unwrap(),
                            );
                            memcpy(
                                ((*this).value.data.ptr as *mut c_char)
                                    .offset((row * dsize * (*this).value.nelem) as isize)
                                    as *mut c_void,
                                ((*theVar).value.data.ptr as *mut c_char)
                                    .offset((elem * dsize) as isize)
                                    as *const c_void,
                                (((*this).value.nelem * dsize) as c_ulong)
                                    .try_into()
                                    .unwrap(),
                            );
                        }
                        row += 1;
                        row;
                    }
                }
            }
        }
        if (*theVar).operation > 0 {
            if (*theVar).ntype == fits_parser_yytokentype::STRING as c_int
                || (*theVar).ntype == fits_parser_yytokentype::BITSTR as c_int
            {
                free(*((*theVar).value.data.strptr).offset(0) as *mut c_void);
            } else {
                free((*theVar).value.data.ptr);
            }
        }
        i = 0;
        while i < nDims {
            if (*theDims[i as usize]).operation > 0 {
                free((*theDims[i as usize]).value.data.ptr);
            }
            i += 1;
            i;
        }
    }
}
unsafe fn Do_GTI(mut lParse: *mut ParseData, mut this: *mut Node) {
    unsafe {
        let mut theExpr: *mut Node = std::ptr::null_mut::<Node>();
        let mut theTimes: *mut Node = std::ptr::null_mut::<Node>();
        let mut start: *mut c_double = std::ptr::null_mut::<c_double>();
        let mut stop: *mut c_double = std::ptr::null_mut::<c_double>();
        let mut times: *mut c_double = std::ptr::null_mut::<c_double>();
        let mut elem: c_long = 0;
        let mut nGTI: c_long = 0;
        let mut gti: c_long = 0;
        let mut ordered: c_int = 0;
        let mut dorow: c_int = ((*this).operation == gtifind_fct as c_int) as c_int;
        theTimes = ((*lParse).Nodes).offset((*this).SubNodes[0] as isize);
        theExpr = ((*lParse).Nodes).offset((*this).SubNodes[1] as isize);
        nGTI = (*theTimes).value.nelem;
        start = (*theTimes).value.data.dblptr;
        stop = ((*theTimes).value.data.dblptr).offset(nGTI as isize);
        ordered = (*theTimes).ntype;
        if (*theExpr).operation == -(1000 as c_int) {
            gti = Search_GTI(
                (*theExpr).value.data.dbl,
                nGTI,
                start,
                stop,
                ordered,
                std::ptr::null_mut::<c_long>(),
            );
            if dorow != 0 {
                (*this).value.data.lng = if gti >= 0 { gti + 1 } else { -(1) as c_long };
            } else {
                (*this).value.data.log = if gti >= 0 { 1 } else { 0 };
            }
            (*this).operation = -(1000 as c_int);
        } else {
            Allocate_Ptrs(lParse, this);
            times = (*theExpr).value.data.dblptr;
            if (*lParse).status == 0 {
                elem = (*lParse).nRows * (*this).value.nelem;
                if nGTI != 0 {
                    gti = -(1) as c_long;
                    loop {
                        let fresh202 = elem;
                        elem -= 1;
                        if fresh202 == 0 {
                            break;
                        }
                        let fresh203 = &mut *((*this).value.undef).offset(elem as isize);
                        *fresh203 = *((*theExpr).value.undef).offset(elem as isize);
                        if *fresh203 != 0 {
                            continue;
                        }
                        if gti < 0
                            || *times.offset(elem as isize) < *start.offset(gti as isize)
                            || *times.offset(elem as isize) > *stop.offset(gti as isize)
                        {
                            gti = Search_GTI(
                                *times.offset(elem as isize),
                                nGTI,
                                start,
                                stop,
                                ordered,
                                std::ptr::null_mut::<c_long>(),
                            );
                        }
                        if dorow != 0 {
                            *((*this).value.data.lngptr).offset(elem as isize) =
                                if gti >= 0 { gti + 1 } else { -(1) as c_long };
                            *((*this).value.undef).offset(elem as isize) =
                                (if gti >= 0 { 0 } else { 1 }) as c_char;
                        } else {
                            *((*this).value.data.logptr).offset(elem as isize) =
                                if gti >= 0 { 1 } else { 0 };
                        }
                    }
                } else if dorow != 0 {
                    loop {
                        let fresh204 = elem;
                        elem -= 1;
                        if fresh204 == 0 {
                            break;
                        }
                        *((*this).value.undef).offset(elem as isize) = 1;
                    }
                } else {
                    loop {
                        let fresh205 = elem;
                        elem -= 1;
                        if fresh205 == 0 {
                            break;
                        }
                        *((*this).value.data.logptr).offset(elem as isize) = 0;
                        *((*this).value.undef).offset(elem as isize) = 0;
                    }
                }
            }
        }
        if (*theExpr).operation > 0 {
            free((*theExpr).value.data.ptr);
        }
    }
}
unsafe fn Do_GTI_Over(mut lParse: *mut ParseData, mut this: *mut Node) {
    unsafe {
        let mut theTimes: *mut Node = std::ptr::null_mut::<Node>();
        let mut theStart: *mut Node = std::ptr::null_mut::<Node>();
        let mut theStop: *mut Node = std::ptr::null_mut::<Node>();
        let mut gtiStart: *mut c_double = std::ptr::null_mut::<c_double>();
        let mut gtiStop: *mut c_double = std::ptr::null_mut::<c_double>();
        let mut evtStart: *mut c_double = std::ptr::null_mut::<c_double>();
        let mut evtStop: *mut c_double = std::ptr::null_mut::<c_double>();
        let mut elem: c_long = 0;
        let mut nGTI: c_long = 0;
        let mut gti: c_long = 0;
        let mut nextGTI: c_long = 0;
        let mut ordered: c_int = 0;
        theTimes = ((*lParse).Nodes).offset((*this).SubNodes[0] as isize);
        theStop = ((*lParse).Nodes).offset((*this).SubNodes[2] as isize);
        theStart = ((*lParse).Nodes).offset((*this).SubNodes[1] as isize);
        nGTI = (*theTimes).value.nelem;
        gtiStart = (*theTimes).value.data.dblptr;
        gtiStop = ((*theTimes).value.data.dblptr).offset(nGTI as isize);
        if (*theStart).operation == -(1000 as c_int) && (*theStop).operation == -(1000 as c_int) {
            (*this).value.data.dbl = GTI_Over(
                (*theStart).value.data.dbl,
                (*theStop).value.data.dbl,
                nGTI,
                gtiStart,
                gtiStop,
                &mut gti,
            );
            (*this).operation = -(1000 as c_int);
        } else {
            let mut undefStart: c_char = 0;
            let mut undefStop: c_char = 0;
            let mut uStart: c_double = 0.;
            let mut uStop: c_double = 0.;
            if (*theStart).operation == -(1000 as c_int) {
                uStart = (*theStart).value.data.dbl;
            }
            if (*theStop).operation == -(1000 as c_int) {
                uStop = (*theStop).value.data.dbl;
            }
            Allocate_Ptrs(lParse, this);
            evtStart = (*theStart).value.data.dblptr;
            evtStop = (*theStop).value.data.dblptr;
            if (*lParse).status == 0 {
                elem = (*lParse).nRows * (*this).value.nelem;
                if nGTI != 0 {
                    let mut toverlap: c_double = 0.0f64;
                    gti = -(1) as c_long;
                    loop {
                        let fresh206 = elem;
                        elem -= 1;
                        if fresh206 == 0 {
                            break;
                        }
                        if (*theStart).operation != -(1000 as c_int) {
                            undefStart = *((*theStart).value.undef).offset(elem as isize);
                            uStart = *evtStart.offset(elem as isize);
                        }
                        if (*theStop).operation != -(1000 as c_int) {
                            undefStop = *((*theStop).value.undef).offset(elem as isize);
                            uStop = *evtStop.offset(elem as isize);
                        }
                        let fresh207 = &mut *((*this).value.undef).offset(elem as isize);
                        *fresh207 = if undefStart as c_int != 0 || undefStop as c_int != 0 {
                            1
                        } else {
                            0
                        };
                        if *fresh207 != 0 {
                            continue;
                        }
                        if gti < 0
                            || uStart < *gtiStart.offset(gti as isize)
                            || uStart > *gtiStop.offset(gti as isize)
                            || uStop < *gtiStart.offset(gti as isize)
                            || uStop > *gtiStop.offset(gti as isize)
                        {
                            toverlap = GTI_Over(uStart, uStop, nGTI, gtiStart, gtiStop, &mut gti);
                        } else {
                            toverlap = uStop - uStart;
                        }
                        *((*this).value.data.dblptr).offset(elem as isize) = toverlap;
                    }
                } else {
                    loop {
                        let fresh208 = elem;
                        elem -= 1;
                        if fresh208 == 0 {
                            break;
                        }
                        *((*this).value.data.dblptr).offset(elem as isize) = 0.0f64;
                        *((*this).value.undef).offset(elem as isize) = 0;
                    }
                }
            }
        }
        if (*theStart).operation > 0 {
            free((*theStart).value.data.ptr);
        }
        if (*theStop).operation > 0 {
            free((*theStop).value.data.ptr);
        }
    }
}
unsafe fn GTI_Over(
    mut evtStart: c_double,
    mut evtStop: c_double,
    mut nGTI: c_long,
    mut start: *mut c_double,
    mut stop: *mut c_double,
    mut gtiout: *mut c_long,
) -> c_double {
    unsafe {
        let mut gti1: c_long = 0;
        let mut gti2: c_long = 0;
        let mut nextGTI1: c_long = 0;
        let mut nextGTI2: c_long = 0;
        let mut gti: c_long = 0;
        let mut nMax: c_long = 0;
        let mut overlap: c_double = 0.0f64;
        *gtiout = -(1 as c_long);
        if evtStop <= evtStart {
            return 0.0f64;
        }
        gti1 = Search_GTI(evtStart, nGTI, start, stop, 1, &mut nextGTI1);
        gti2 = Search_GTI(evtStop, nGTI, start, stop, 1, &mut nextGTI2);
        if gti1 >= 0 {
            *gtiout = gti1;
        }
        if nextGTI1 < 0 && nextGTI2 < 0 {
            return 0.0f64;
        }
        if gti1 < 0 && gti2 < 0 && nextGTI1 == nextGTI2 {
            return 0.0f64;
        }
        if gti1 >= 0 && gti1 == gti2 {
            return evtStop - evtStart;
        }
        if nextGTI2 < 0 {
            nMax = nGTI - 1;
        } else if gti2 >= 0 {
            nMax = nextGTI2;
        } else {
            nMax = nextGTI2 - 1;
        }
        gti = nextGTI1;
        while gti <= nMax {
            let mut starti: c_double = *start.offset(gti as isize);
            let mut stopi: c_double = *stop.offset(gti as isize);
            if evtStart > starti {
                starti = evtStart;
            }
            if evtStop < stopi {
                stopi = evtStop;
            }
            overlap += stopi - starti;
            gti += 1;
            gti;
        }
        overlap
    }
}
unsafe fn Search_GTI(
    mut evtTime: c_double,
    mut nGTI: c_long,
    mut start: *mut c_double,
    mut stop: *mut c_double,
    mut ordered: c_int,
    mut nextGTI0: *mut c_long,
) -> c_long {
    unsafe {
        let mut gti: c_long = 0;
        let mut nextGTI: c_long = -(1 as c_long);
        let mut step: c_long = 0;
        if ordered != 0 && nGTI > 15 as c_long {
            if evtTime >= *start.offset(0) && evtTime <= *stop.offset((nGTI - 1) as isize) {
                step = nGTI >> 1;
                gti = step;
                loop {
                    if step > 1 {
                        step >>= 1;
                    }
                    if evtTime > *stop.offset(gti as isize) {
                        if evtTime >= *start.offset((gti + 1) as isize) {
                            gti += step;
                        } else {
                            nextGTI = gti + 1;
                            gti = -(1 as c_long);
                            break;
                        }
                    } else if evtTime < *start.offset(gti as isize) {
                        if evtTime <= *stop.offset((gti - 1) as isize) {
                            gti -= step;
                        } else {
                            nextGTI = gti;
                            gti = -(1 as c_long);
                            break;
                        }
                    } else {
                        nextGTI = gti;
                        break;
                    }
                }
            } else {
                if *start.offset(0) > evtTime {
                    nextGTI = 0;
                }
                gti = -(1 as c_long);
            }
        } else {
            gti = nGTI;
            loop {
                let fresh209 = gti;
                gti -= 1;
                if fresh209 == 0 {
                    break;
                }
                if *stop.offset(gti as isize) >= evtTime {
                    nextGTI = gti;
                }
                if evtTime >= *start.offset(gti as isize) && evtTime <= *stop.offset(gti as isize) {
                    break;
                }
            }
        }
        if nextGTI >= nGTI {
            nextGTI = -(1) as c_long;
        }
        if !nextGTI0.is_null() {
            *nextGTI0 = nextGTI;
        }
        gti
    }
}
unsafe fn Do_REG(mut lParse: *mut ParseData, mut this: *mut Node) {
    unsafe {
        let mut theRegion: *mut Node = std::ptr::null_mut::<Node>();
        let mut theX: *mut Node = std::ptr::null_mut::<Node>();
        let mut theY: *mut Node = std::ptr::null_mut::<Node>();
        let mut Xval: c_double = 0.0f64;
        let mut Yval: c_double = 0.0f64;
        let mut Xnull: c_char = 0;
        let mut Ynull: c_char = 0;
        let mut Xvector: c_int = 0;
        let mut Yvector: c_int = 0;
        let mut nelem: c_long = 0;
        let mut elem: c_long = 0;
        let mut rows: c_long = 0;
        theRegion = ((*lParse).Nodes).offset((*this).SubNodes[0] as isize);
        theX = ((*lParse).Nodes).offset((*this).SubNodes[1] as isize);
        theY = ((*lParse).Nodes).offset((*this).SubNodes[2] as isize);
        Xvector = ((*theX).operation != -(1000 as c_int)) as c_int;
        if Xvector != 0 {
            Xvector = (*theX).value.nelem as c_int;
        } else {
            Xval = (*theX).value.data.dbl;
        }
        Yvector = ((*theY).operation != -(1000 as c_int)) as c_int;
        if Yvector != 0 {
            Yvector = (*theY).value.nelem as c_int;
        } else {
            Yval = (*theY).value.data.dbl;
        }
        if Xvector == 0 && Yvector == 0 {
            (*this).value.data.log = if fits_in_region(
                Xval,
                Yval,
                &mut *((*theRegion).value.data.ptr as *mut SAORegion),
            ) != 0
            {
                1
            } else {
                0
            };
            (*this).operation = -(1000 as c_int);
        } else {
            Allocate_Ptrs(lParse, this);
            if (*lParse).status == 0 {
                rows = (*lParse).nRows;
                nelem = (*this).value.nelem;
                elem = rows * nelem;
                loop {
                    let fresh210 = rows;
                    rows -= 1;
                    if fresh210 == 0 {
                        break;
                    }
                    loop {
                        let fresh211 = nelem;
                        nelem -= 1;
                        if fresh211 == 0 {
                            break;
                        }
                        elem -= 1;
                        elem;
                        if Xvector > 1 {
                            Xval = *((*theX).value.data.dblptr).offset(elem as isize);
                            Xnull = *((*theX).value.undef).offset(elem as isize);
                        } else if Xvector != 0 {
                            Xval = *((*theX).value.data.dblptr).offset(rows as isize);
                            Xnull = *((*theX).value.undef).offset(rows as isize);
                        }
                        if Yvector > 1 {
                            Yval = *((*theY).value.data.dblptr).offset(elem as isize);
                            Ynull = *((*theY).value.undef).offset(elem as isize);
                        } else if Yvector != 0 {
                            Yval = *((*theY).value.data.dblptr).offset(rows as isize);
                            Ynull = *((*theY).value.undef).offset(rows as isize);
                        }
                        *((*this).value.undef).offset(elem as isize) =
                            if Xnull as c_int != 0 || Ynull as c_int != 0 {
                                1
                            } else {
                                0
                            };
                        if *((*this).value.undef).offset(elem as isize) != 0 {
                            continue;
                        }
                        *((*this).value.data.logptr).offset(elem as isize) = if fits_in_region(
                            Xval,
                            Yval,
                            &mut *((*theRegion).value.data.ptr as *mut SAORegion),
                        ) != 0
                        {
                            1
                        } else {
                            0
                        };
                    }
                    nelem = (*this).value.nelem;
                }
            }
        }
        if (*theX).operation > 0 {
            free((*theX).value.data.ptr);
        }
        if (*theY).operation > 0 {
            free((*theY).value.data.ptr);
        }
    }
}
unsafe fn Do_Vector(mut lParse: *mut ParseData, mut this: *mut Node) {
    unsafe {
        let mut that: *mut Node = std::ptr::null_mut::<Node>();
        let mut row: c_long = 0;
        let mut elem: c_long = 0;
        let mut idx: c_long = 0;
        let mut jdx: c_long = 0;
        let mut offset: c_long = 0;
        let mut node: c_int = 0;
        Allocate_Ptrs(lParse, this);
        if (*lParse).status == 0 {
            node = 0;
            while node < (*this).nSubNodes {
                that = ((*lParse).Nodes).offset((*this).SubNodes[node as usize] as isize);
                if (*that).operation == -(1000 as c_int) {
                    idx = (*lParse).nRows * (*this).value.nelem + offset;
                    loop {
                        idx -= (*this).value.nelem;
                        if idx < 0 {
                            break;
                        }
                        *((*this).value.undef).offset(idx as isize) = 0;
                        match (*this).ntype {
                            258 => {
                                *((*this).value.data.logptr).offset(idx as isize) =
                                    (*that).value.data.log;
                            }
                            259 => {
                                *((*this).value.data.lngptr).offset(idx as isize) =
                                    (*that).value.data.lng;
                            }
                            260 => {
                                *((*this).value.data.dblptr).offset(idx as isize) =
                                    (*that).value.data.dbl;
                            }
                            _ => {}
                        }
                    }
                } else {
                    row = (*lParse).nRows;
                    idx = row * (*that).value.nelem;
                    loop {
                        let fresh212 = row;
                        row -= 1;
                        if fresh212 == 0 {
                            break;
                        }
                        elem = (*that).value.nelem;
                        jdx = row * (*this).value.nelem + offset;
                        loop {
                            let fresh213 = elem;
                            elem -= 1;
                            if fresh213 == 0 {
                                break;
                            }
                            idx -= 1;
                            *((*this).value.undef).offset((jdx + elem) as isize) =
                                *((*that).value.undef).offset(idx as isize);
                            match (*this).ntype {
                                258 => {
                                    *((*this).value.data.logptr).offset((jdx + elem) as isize) =
                                        *((*that).value.data.logptr).offset(idx as isize);
                                }
                                259 => {
                                    *((*this).value.data.lngptr).offset((jdx + elem) as isize) =
                                        *((*that).value.data.lngptr).offset(idx as isize);
                                }
                                260 => {
                                    *((*this).value.data.dblptr).offset((jdx + elem) as isize) =
                                        *((*that).value.data.dblptr).offset(idx as isize);
                                }
                                _ => {}
                            }
                        }
                    }
                }
                offset += (*that).value.nelem;
                node += 1;
                node;
            }
        }
        node = 0;
        while node < (*this).nSubNodes {
            if (*((*lParse).Nodes).offset((*this).SubNodes[node as usize] as isize)).operation > 0 {
                free(
                    (*((*lParse).Nodes).offset((*this).SubNodes[node as usize] as isize))
                        .value
                        .data
                        .ptr,
                );
            }
            node += 1;
            node;
        }
    }
}

unsafe fn Do_Array(mut lParse: *mut ParseData, mut this: *mut Node) {
    unsafe {
        let mut that: *mut Node = std::ptr::null_mut::<Node>();
        let mut row: c_long = 0;
        let mut elem: c_long = 0;
        let mut idx: c_long = 0;
        let mut jdx: c_long = 0;
        let mut offset: c_long = 0;
        let mut node: c_int = 0;
        Allocate_Ptrs(lParse, this);
        if (*lParse).status == 0 {
            that = ((*lParse).Nodes).offset((*this).SubNodes[0] as isize);
            if (*that).operation == -(1000 as c_int) {
                idx = (*lParse).nRows * (*this).value.nelem + offset;
                loop {
                    let fresh214 = idx;
                    idx -= 1;
                    if fresh214 == 0 {
                        break;
                    }
                    *((*this).value.undef).offset(idx as isize) = 0;
                    match (*this).ntype {
                        258 => {
                            *((*this).value.data.logptr).offset(idx as isize) =
                                (*that).value.data.log;
                        }
                        259 => {
                            *((*this).value.data.lngptr).offset(idx as isize) =
                                (*that).value.data.lng;
                        }
                        260 => {
                            *((*this).value.data.dblptr).offset(idx as isize) =
                                (*that).value.data.dbl;
                        }
                        _ => {}
                    }
                }
            } else if (*that).value.nelem > 1 {
                idx = (*lParse).nRows * (*this).value.nelem;
                loop {
                    let fresh215 = idx;
                    idx -= 1;
                    if fresh215 == 0 {
                        break;
                    }
                    *((*this).value.undef).offset(idx as isize) =
                        *((*that).value.undef).offset(idx as isize);
                    match (*this).ntype {
                        258 => {
                            *((*this).value.data.logptr).offset(idx as isize) =
                                *((*that).value.data.logptr).offset(idx as isize);
                        }
                        259 => {
                            *((*this).value.data.lngptr).offset(idx as isize) =
                                *((*that).value.data.lngptr).offset(idx as isize);
                        }
                        260 => {
                            *((*this).value.data.dblptr).offset(idx as isize) =
                                *((*that).value.data.dblptr).offset(idx as isize);
                        }
                        _ => {}
                    }
                }
            } else {
                row = (*lParse).nRows;
                idx = row * (*this).value.nelem - 1;
                loop {
                    let fresh216 = row;
                    row -= 1;
                    if fresh216 == 0 {
                        break;
                    }
                    elem = (*this).value.nelem;
                    loop {
                        let fresh217 = elem;
                        elem -= 1;
                        if fresh217 == 0 {
                            break;
                        }
                        *((*this).value.undef).offset(idx as isize) =
                            *((*that).value.undef).offset(row as isize);
                        match (*this).ntype {
                            258 => {
                                *((*this).value.data.logptr).offset(idx as isize) =
                                    *((*that).value.data.logptr).offset(row as isize);
                            }
                            259 => {
                                *((*this).value.data.lngptr).offset(idx as isize) =
                                    *((*that).value.data.lngptr).offset(row as isize);
                            }
                            260 => {
                                *((*this).value.data.dblptr).offset(idx as isize) =
                                    *((*that).value.data.dblptr).offset(row as isize);
                            }
                            _ => {}
                        }
                        idx -= 1;
                        idx;
                    }
                }
            }
            if (*((*lParse).Nodes).offset((*this).SubNodes[0] as isize)).operation > 0 {
                free(
                    (*((*lParse).Nodes).offset((*this).SubNodes[0] as isize))
                        .value
                        .data
                        .ptr,
                );
            }
        }
    }
}

unsafe fn bitlgte(mut bits1: *mut c_char, mut oper: c_int, mut bits2: *mut c_char) -> c_char {
    unsafe {
        let mut val1: c_int = 0;
        let mut val2: c_int = 0;
        let mut nextbit: c_int = 0;
        let mut result: c_char = 0;
        let mut i: c_int = 0;
        let mut l1: c_int = 0;
        let mut l2: c_int = 0;
        let mut length: c_int = 0;
        let mut ldiff: c_int = 0;
        let mut stream: *mut c_char = ptr::null_mut();
        let mut chr1: c_char = 0;
        let mut chr2: c_char = 0;
        l1 = strlen(bits1) as c_int;
        l2 = strlen(bits2) as c_int;
        length = if l1 > l2 { l1 } else { l2 };
        stream = malloc(
            (::core::mem::size_of::<c_char>() as c_ulong)
                .wrapping_mul((length + 1) as c_ulong)
                .try_into()
                .unwrap(),
        ) as *mut c_char;
        if l1 < l2 {
            ldiff = l2 - l1;
            i = 0;
            loop {
                let fresh218 = ldiff;
                ldiff -= 1;
                if fresh218 == 0 {
                    break;
                }
                let fresh219 = i;
                i += 1;
                *stream.offset(fresh219 as isize) = b'0' as c_char;
            }
            loop {
                let fresh220 = l1;
                l1 -= 1;
                if fresh220 == 0 {
                    break;
                }
                let fresh221 = bits1;
                bits1 = bits1.offset(1);
                let fresh222 = i;
                i += 1;
                *stream.offset(fresh222 as isize) = *fresh221;
            }
            *stream.offset(i as isize) = 0;
            bits1 = stream;
        } else if l2 < l1 {
            ldiff = l1 - l2;
            i = 0;
            loop {
                let fresh223 = ldiff;
                ldiff -= 1;
                if fresh223 == 0 {
                    break;
                }
                let fresh224 = i;
                i += 1;
                *stream.offset(fresh224 as isize) = b'0' as c_char;
            }
            loop {
                let fresh225 = l2;
                l2 -= 1;
                if fresh225 == 0 {
                    break;
                }
                let fresh226 = bits2;
                bits2 = bits2.offset(1);
                let fresh227 = i;
                i += 1;
                *stream.offset(fresh227 as isize) = *fresh226;
            }
            *stream.offset(i as isize) = 0;
            bits2 = stream;
        }
        val2 = 0;
        val1 = val2;
        nextbit = 1;
        loop {
            let fresh228 = length;
            length -= 1;
            if fresh228 == 0 {
                break;
            }
            chr1 = *bits1.offset(length as isize);
            chr2 = *bits2.offset(length as isize);
            if chr1 as c_int != 'x' as i32
                && chr1 as c_int != 'X' as i32
                && chr2 as c_int != 'x' as i32
                && chr2 as c_int != 'X' as i32
            {
                if chr1 as c_int == '1' as i32 {
                    val1 += nextbit;
                }
                if chr2 as c_int == '1' as i32 {
                    val2 += nextbit;
                }
                nextbit *= 2 as c_int;
            }
        }
        result = 0;
        match oper {
            282 => {
                if val1 < val2 {
                    result = 1;
                }
            }
            283 => {
                if val1 <= val2 {
                    result = 1;
                }
            }
            281 => {
                if val1 > val2 {
                    result = 1;
                }
            }
            284 => {
                if val1 >= val2 {
                    result = 1;
                }
            }
            _ => {}
        }
        free(stream as *mut c_void);
        result
    }
}
unsafe fn bitand(mut result: *mut c_char, mut bitstrm1: *mut c_char, mut bitstrm2: *mut c_char) {
    unsafe {
        let mut i: c_int = 0;
        let mut l1: c_int = 0;
        let mut l2: c_int = 0;
        let mut ldiff: c_int = 0;
        let mut largestStream: c_int = 0;
        let mut stream: *mut c_char = ptr::null_mut();
        let mut chr1: c_char = 0;
        let mut chr2: c_char = 0;
        l1 = strlen(bitstrm1) as c_int;
        l2 = strlen(bitstrm2) as c_int;
        largestStream = if l1 > l2 { l1 } else { l2 };
        stream = malloc(
            (::core::mem::size_of::<c_char>() as c_ulong)
                .wrapping_mul((largestStream + 1) as c_ulong)
                .try_into()
                .unwrap(),
        ) as *mut c_char;
        if l1 < l2 {
            ldiff = l2 - l1;
            i = 0;
            loop {
                let fresh229 = ldiff;
                ldiff -= 1;
                if fresh229 == 0 {
                    break;
                }
                let fresh230 = i;
                i += 1;
                *stream.offset(fresh230 as isize) = b'0' as c_char;
            }
            loop {
                let fresh231 = l1;
                l1 -= 1;
                if fresh231 == 0 {
                    break;
                }
                let fresh232 = bitstrm1;
                bitstrm1 = bitstrm1.offset(1);
                let fresh233 = i;
                i += 1;
                *stream.offset(fresh233 as isize) = *fresh232;
            }
            *stream.offset(i as isize) = 0;
            bitstrm1 = stream;
        } else if l2 < l1 {
            ldiff = l1 - l2;
            i = 0;
            loop {
                let fresh234 = ldiff;
                ldiff -= 1;
                if fresh234 == 0 {
                    break;
                }
                let fresh235 = i;
                i += 1;
                *stream.offset(fresh235 as isize) = b'0' as c_char;
            }
            loop {
                let fresh236 = l2;
                l2 -= 1;
                if fresh236 == 0 {
                    break;
                }
                let fresh237 = bitstrm2;
                bitstrm2 = bitstrm2.offset(1);
                let fresh238 = i;
                i += 1;
                *stream.offset(fresh238 as isize) = *fresh237;
            }
            *stream.offset(i as isize) = 0;
            bitstrm2 = stream;
        }
        loop {
            let fresh239 = bitstrm1;
            bitstrm1 = bitstrm1.offset(1);
            chr1 = *fresh239;
            if chr1 == 0 {
                break;
            }
            let fresh240 = bitstrm2;
            bitstrm2 = bitstrm2.offset(1);
            chr2 = *fresh240;
            if chr1 as c_int == 'x' as i32 || chr2 as c_int == 'x' as i32 {
                *result = b'x' as c_char;
            } else if chr1 as c_int == '1' as i32 && chr2 as c_int == '1' as i32 {
                *result = b'1' as c_char;
            } else {
                *result = b'0' as c_char;
            }
            result = result.offset(1);
            result;
        }
        free(stream as *mut c_void);
        *result = 0;
    }
}
unsafe fn bitor(mut result: *mut c_char, mut bitstrm1: *mut c_char, mut bitstrm2: *mut c_char) {
    unsafe {
        let mut i: c_int = 0;
        let mut l1: c_int = 0;
        let mut l2: c_int = 0;
        let mut ldiff: c_int = 0;
        let mut largestStream: c_int = 0;
        let mut stream: *mut c_char = ptr::null_mut();
        let mut chr1: c_char = 0;
        let mut chr2: c_char = 0;
        l1 = strlen(bitstrm1) as c_int;
        l2 = strlen(bitstrm2) as c_int;
        largestStream = if l1 > l2 { l1 } else { l2 };
        stream = malloc(
            (::core::mem::size_of::<c_char>() as c_ulong)
                .wrapping_mul((largestStream + 1) as c_ulong)
                .try_into()
                .unwrap(),
        ) as *mut c_char;
        if l1 < l2 {
            ldiff = l2 - l1;
            i = 0;
            loop {
                let fresh241 = ldiff;
                ldiff -= 1;
                if fresh241 == 0 {
                    break;
                }
                let fresh242 = i;
                i += 1;
                *stream.offset(fresh242 as isize) = b'0' as c_char;
            }
            loop {
                let fresh243 = l1;
                l1 -= 1;
                if fresh243 == 0 {
                    break;
                }
                let fresh244 = bitstrm1;
                bitstrm1 = bitstrm1.offset(1);
                let fresh245 = i;
                i += 1;
                *stream.offset(fresh245 as isize) = *fresh244;
            }
            *stream.offset(i as isize) = 0;
            bitstrm1 = stream;
        } else if l2 < l1 {
            ldiff = l1 - l2;
            i = 0;
            loop {
                let fresh246 = ldiff;
                ldiff -= 1;
                if fresh246 == 0 {
                    break;
                }
                let fresh247 = i;
                i += 1;
                *stream.offset(fresh247 as isize) = b'0' as c_char;
            }
            loop {
                let fresh248 = l2;
                l2 -= 1;
                if fresh248 == 0 {
                    break;
                }
                let fresh249 = bitstrm2;
                bitstrm2 = bitstrm2.offset(1);
                let fresh250 = i;
                i += 1;
                *stream.offset(fresh250 as isize) = *fresh249;
            }
            *stream.offset(i as isize) = 0;
            bitstrm2 = stream;
        }
        loop {
            let fresh251 = bitstrm1;
            bitstrm1 = bitstrm1.offset(1);
            chr1 = *fresh251;
            if chr1 == 0 {
                break;
            }
            let fresh252 = bitstrm2;
            bitstrm2 = bitstrm2.offset(1);
            chr2 = *fresh252;
            if chr1 as c_int == '1' as i32 || chr2 as c_int == '1' as i32 {
                *result = b'1' as c_char;
            } else if chr1 as c_int == '0' as i32 || chr2 as c_int == '0' as i32 {
                *result = b'0' as c_char;
            } else {
                *result = b'x' as c_char;
            }
            result = result.offset(1);
            result;
        }
        free(stream as *mut c_void);
        *result = 0;
    }
}
unsafe fn bitnot(mut result: *mut c_char, mut bits: *mut c_char) {
    unsafe {
        let mut length: c_int = 0;
        let mut chr: c_char = 0;
        length = strlen(bits) as c_int;
        loop {
            let fresh253 = length;
            length -= 1;
            if fresh253 == 0 {
                break;
            }
            let fresh254 = bits;
            bits = bits.offset(1);
            chr = *fresh254;
            let fresh255 = result;
            result = result.offset(1);
            *fresh255 = (if chr as c_int == '1' as i32 {
                '0' as i32
            } else if chr as c_int == '0' as i32 {
                '1' as i32
            } else {
                chr as c_int
            }) as c_char;
        }
        *result = 0;
    }
}
unsafe fn bitcmp(mut bitstrm1: *mut c_char, mut bitstrm2: *mut c_char) -> c_char {
    unsafe {
        let mut i: c_int = 0;
        let mut l1: c_int = 0;
        let mut l2: c_int = 0;
        let mut ldiff: c_int = 0;
        let mut largestStream: c_int = 0;
        let mut stream: *mut c_char = ptr::null_mut();
        let mut chr1: c_char = 0;
        let mut chr2: c_char = 0;
        l1 = strlen(bitstrm1) as c_int;
        l2 = strlen(bitstrm2) as c_int;
        largestStream = if l1 > l2 { l1 } else { l2 };
        stream = malloc(
            (::core::mem::size_of::<c_char>() as c_ulong)
                .wrapping_mul((largestStream + 1) as c_ulong)
                .try_into()
                .unwrap(),
        ) as *mut c_char;
        if l1 < l2 {
            ldiff = l2 - l1;
            i = 0;
            loop {
                let fresh256 = ldiff;
                ldiff -= 1;
                if fresh256 == 0 {
                    break;
                }
                let fresh257 = i;
                i += 1;
                *stream.offset(fresh257 as isize) = b'0' as c_char;
            }
            loop {
                let fresh258 = l1;
                l1 -= 1;
                if fresh258 == 0 {
                    break;
                }
                let fresh259 = bitstrm1;
                bitstrm1 = bitstrm1.offset(1);
                let fresh260 = i;
                i += 1;
                *stream.offset(fresh260 as isize) = *fresh259;
            }
            *stream.offset(i as isize) = 0;
            bitstrm1 = stream;
        } else if l2 < l1 {
            ldiff = l1 - l2;
            i = 0;
            loop {
                let fresh261 = ldiff;
                ldiff -= 1;
                if fresh261 == 0 {
                    break;
                }
                let fresh262 = i;
                i += 1;
                *stream.offset(fresh262 as isize) = b'0' as c_char;
            }
            loop {
                let fresh263 = l2;
                l2 -= 1;
                if fresh263 == 0 {
                    break;
                }
                let fresh264 = bitstrm2;
                bitstrm2 = bitstrm2.offset(1);
                let fresh265 = i;
                i += 1;
                *stream.offset(fresh265 as isize) = *fresh264;
            }
            *stream.offset(i as isize) = 0;
            bitstrm2 = stream;
        }
        loop {
            let fresh266 = bitstrm1;
            bitstrm1 = bitstrm1.offset(1);
            chr1 = *fresh266;
            if chr1 == 0 {
                break;
            }
            let fresh267 = bitstrm2;
            bitstrm2 = bitstrm2.offset(1);
            chr2 = *fresh267;
            if chr1 as c_int == '0' as i32 && chr2 as c_int == '1' as i32
                || chr1 as c_int == '1' as i32 && chr2 as c_int == '0' as i32
            {
                free(stream as *mut c_void);
                return 0;
            }
        }
        free(stream as *mut c_void);
        1
    }
}
unsafe fn bnear(mut x: c_double, mut y: c_double, mut tolerance: c_double) -> c_char {
    if fabs(x - y) < tolerance { 1 } else { 0 }
}
unsafe fn saobox(
    mut xcen: c_double,
    mut ycen: c_double,
    mut xwid: c_double,
    mut ywid: c_double,
    mut rot: c_double,
    mut xcol: c_double,
    mut ycol: c_double,
) -> c_char {
    let mut x: c_double = 0.;
    let mut y: c_double = 0.;
    let mut xprime: c_double = 0.;
    let mut yprime: c_double = 0.;
    let mut xmin: c_double = 0.;
    let mut xmax: c_double = 0.;
    let mut ymin: c_double = 0.;
    let mut ymax: c_double = 0.;
    let mut theta: c_double = 0.;
    theta = rot / 180.0f64 * 3.141_592_653_589_793_f64;
    xprime = xcol - xcen;
    yprime = ycol - ycen;
    x = xprime * cos(theta) + yprime * sin(theta);
    y = -xprime * sin(theta) + yprime * cos(theta);
    xmin = -0.5f64 * xwid;
    xmax = 0.5f64 * xwid;
    ymin = -0.5f64 * ywid;
    ymax = 0.5f64 * ywid;
    if x >= xmin && x <= xmax && y >= ymin && y <= ymax {
        1
    } else {
        0
    }
}
unsafe fn circle(
    mut xcen: c_double,
    mut ycen: c_double,
    mut rad: c_double,
    mut xcol: c_double,
    mut ycol: c_double,
) -> c_char {
    let mut r2: c_double = 0.;
    let mut dx: c_double = 0.;
    let mut dy: c_double = 0.;
    let mut dlen: c_double = 0.;
    dx = xcol - xcen;
    dy = ycol - ycen;
    dx *= dx;
    dy *= dy;
    dlen = dx + dy;
    r2 = rad * rad;
    if dlen <= r2 { 1 } else { 0 }
}
unsafe fn ellipse(
    mut xcen: c_double,
    mut ycen: c_double,
    mut xrad: c_double,
    mut yrad: c_double,
    mut rot: c_double,
    mut xcol: c_double,
    mut ycol: c_double,
) -> c_char {
    let mut x: c_double = 0.;
    let mut y: c_double = 0.;
    let mut xprime: c_double = 0.;
    let mut yprime: c_double = 0.;
    let mut dx: c_double = 0.;
    let mut dy: c_double = 0.;
    let mut dlen: c_double = 0.;
    let mut theta: c_double = 0.;
    theta = rot / 180.0f64 * 3.141_592_653_589_793_f64;
    xprime = xcol - xcen;
    yprime = ycol - ycen;
    x = xprime * cos(theta) + yprime * sin(theta);
    y = -xprime * sin(theta) + yprime * cos(theta);
    dx = x / xrad;
    dy = y / yrad;
    dx *= dx;
    dy *= dy;
    dlen = dx + dy;
    if dlen <= 1.0f64 { 1 } else { 0 }
}
unsafe fn cstrmid(
    mut lParse: *mut ParseData,
    mut dest_str: *mut c_char,
    mut dest_len: c_int,
    mut src_str: *mut c_char,
    mut src_len: c_int,
    mut pos: c_int,
) -> c_int {
    unsafe {
        let mut fill_char: c_char = 0;
        if src_len == 0 {
            src_len = strlen(src_str) as c_int;
        }
        if pos < 0 {
            fits_parser_yyerror(
                0 as yyscan_t,
                lParse,
                b"STRMID(S,P,N) P must be 0 or greater\0" as *const u8 as *const c_char
                    as *mut c_char,
            );
            return -(1);
        }
        if pos > src_len || pos == 0 {
            memset(
                dest_str as *mut c_void,
                fill_char as c_int,
                (dest_len as c_ulong).try_into().unwrap(),
            );
        } else if pos + dest_len > src_len {
            let mut nsub: c_int = src_len - pos + 1;
            let mut npad: c_int = dest_len - nsub;
            memcpy(
                dest_str as *mut c_void,
                src_str.offset(pos as isize).offset(-(1 as c_int as isize)) as *const c_void,
                (nsub as c_ulong).try_into().unwrap(),
            );
            memset(
                dest_str.offset(nsub as isize) as *mut c_void,
                fill_char as c_int,
                (npad as c_ulong).try_into().unwrap(),
            );
        } else {
            memcpy(
                dest_str as *mut c_void,
                src_str.offset(pos as isize).offset(-(1 as c_int as isize)) as *const c_void,
                (dest_len as c_ulong).try_into().unwrap(),
            );
        }
        *dest_str.offset(dest_len as isize) = 0;
        0
    }
}
unsafe fn fits_parser_yyerror(
    mut scanner: yyscan_t,
    mut lParse: *mut ParseData,
    mut s: *mut c_char,
) {
    unsafe {
        let mut msg: [c_char; 80] = [0; 80];
        if (*lParse).status == 0 {
            (*lParse).status = 431 as c_int;
        }
        strncpy(
            msg.as_mut_ptr(),
            s,
            (80 as c_int as c_ulong).try_into().unwrap(),
        );
        msg[79] = 0;
        ffpmsg(msg.as_mut_ptr());
    }
}
