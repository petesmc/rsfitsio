use alloc::rc::Rc;
use core::ffi::CStr;
use core::slice;
use core::{cmp, ptr};

use bytemuck::{cast_slice, cast_slice_mut};
use libc::{calloc, free, malloc, memcpy, time, time_t};

const APPROX: f64 = 1.0e-7;

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

use crate::c_types::{
    c_char, c_double, c_int, c_long, c_schar, c_short, c_uchar, c_uint, c_ulong, c_void,
};
use crate::cfileio::{ffclos_safe, ffexts_safe, ffopen_safe};
use crate::eval_defs::{CONST_OP, MAXDIMS, MAXSUBS, Node, ParseData, data_union, lval, yyscan_t};
use crate::eval_l::{fits_parser_yyGetVariable, fits_parser_yylex, yyguts_t};
use crate::eval_tab::fits_parser_yytokentype::BITSTR;
use crate::eval_tab::{FITS_PARSER_YYSTYPE, fits_parser_yytokentype};
use crate::fitscore::ffpmsg_slice;
use crate::fitscore::{ffgcno_safe, ffghdn_safe, ffmahd_safe, ffmnhd_safe, ffupch_safe};
use crate::fitsio::{LONG_MAX, LONG_MIN, LONGLONG, MEMORY_ALLOCATION, PARSE_SYNTAX_ERR, fitsfile};
use crate::getcold::ffgcvd_safe;
use crate::getkey::{ffgkyd_safe, ffgkyj_safe, ffgkys_safe};
use crate::region::{MY_PI, SAORegion, WCSdata, fits_in_region, fits_read_rgnfile};
use crate::simplerng::{
    simplerng_getnorm, simplerng_getpoisson, simplerng_getuniform, simplerng_srand,
};
use crate::wcssub::ffgtcs_safe;
use crate::wrappers::strncpy_safe;
use crate::wrappers::{strcat, strcmp, strcpy, strlen, strstr};
use crate::{atoi, cs, int_snprintf};

pub type yy_state_t = yytype_int16;
pub type yytype_int16 = c_short;
pub type yysymbol_kind_t = c_int;
pub const YYSYMBOL_SEXPR: yysymbol_kind_t = 65;
pub const YYSYMBOL_BITS: yysymbol_kind_t = 64;
pub const YYSYMBOL_BEXPR: yysymbol_kind_t = 63;
pub const YYSYMBOL_EXPR: yysymbol_kind_t = 62;
pub const YYSYMBOL_VECTOR: yysymbol_kind_t = 61;
pub const YYSYMBOL_BVECTOR: yysymbol_kind_t = 60;
pub const YYSYMBOL_LINE: yysymbol_kind_t = 59;
pub const YYSYMBOL_LINES: yysymbol_kind_t = 58;
pub const YYSYMBOL_YYACCEPT: yysymbol_kind_t = 57;
pub const YYSYMBOL_56_: yysymbol_kind_t = 56;
pub const YYSYMBOL_55_: yysymbol_kind_t = 55;
pub const YYSYMBOL_54_: yysymbol_kind_t = 54;
pub const YYSYMBOL_53_N: yysymbol_kind_t = 53;
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
pub const YYSYMBOL_YYERROR: yysymbol_kind_t = 1;
pub const YYSYMBOL_YYEOF: yysymbol_kind_t = 0;
pub const YYSYMBOL_YYEMPTY: yysymbol_kind_t = -2;
pub type yytype_int8 = c_schar;
pub type yy_state_fast_t = c_int;
pub type funcOp = c_uint;
pub const ARRAY_FCT: funcOp = 1051;
pub const AXISELEM_FCT: funcOp = 1050;
pub const ELEMNUM_FCT: funcOp = 1049;
pub const GTIFIND_FCT: funcOp = 1048;
pub const GTIOVER_FCT: funcOp = 1047;
pub const SETNULL_FCT: funcOp = 1046;
pub const STRPOS_FCT: funcOp = 1045;
pub const STRMID_FCT: funcOp = 1044;
pub const POIRND_FCT: funcOp = 1043;
pub const GASRND_FCT: funcOp = 1042;
pub const ANGSEP_FCT: funcOp = 1041;
pub const NONNULL_FCT: funcOp = 1040;
pub const STDDEV_FCT: funcOp = 1039;
pub const AVERAGE_FCT: funcOp = 1038;
pub const MEDIAN_FCT: funcOp = 1037;
pub const NULL_FCT: funcOp = 1036;
pub const ROW_FCT: funcOp = 1035;
pub const IFTHENELSE_FCT: funcOp = 1034;
pub const REGFILT_FCT: funcOp = 1033;
pub const GTIFILT_FCT: funcOp = 1032;
pub const DEFNULL_FCT: funcOp = 1031;
pub const ISNULL_FCT: funcOp = 1030;
pub const ELPS_FCT: funcOp = 1029;
pub const BOX_FCT: funcOp = 1028;
pub const CIRCLE_FCT: funcOp = 1027;
pub const NEAR_FCT: funcOp = 1026;
pub const MAX2_FCT: funcOp = 1025;
pub const MAX1_FCT: funcOp = 1024;
pub const MIN2_FCT: funcOp = 1023;
pub const MIN1_FCT: funcOp = 1022;
pub const ROUND_FCT: funcOp = 1021;
pub const FLOOR_FCT: funcOp = 1020;
pub const CEIL_FCT: funcOp = 1019;
pub const ATAN2_FCT: funcOp = 1018;
pub const ABS_FCT: funcOp = 1017;
pub const SQRT_FCT: funcOp = 1016;
pub const LOG10_FCT: funcOp = 1015;
pub const LOG_FCT: funcOp = 1014;
pub const EXP_FCT: funcOp = 1013;
pub const TANH_FCT: funcOp = 1012;
pub const COSH_FCT: funcOp = 1011;
pub const SINH_FCT: funcOp = 1010;
pub const ATAN_FCT: funcOp = 1009;
pub const ACOS_FCT: funcOp = 1008;
pub const ASIN_FCT: funcOp = 1007;
pub const TAN_FCT: funcOp = 1006;
pub const COS_FCT: funcOp = 1005;
pub const SIN_FCT: funcOp = 1004;
pub const NELEM_FCT: funcOp = 1003;
pub const SUM_FCT: funcOp = 1002;
pub const RND_FCT: funcOp = 1001;

pub type shapeType = c_uint;
pub const BPANDA_RGN: shapeType = 14;
pub const EPANDA_RGN: shapeType = 13;
pub const PANDA_RGN: shapeType = 12;
pub const POLY_RGN: shapeType = 11;
pub const SECTOR_RGN: shapeType = 10;
pub const DIAMOND_RGN: shapeType = 9;
pub const RECTANGLE_RGN: shapeType = 8;
pub const BOXANNULUS_RGN: shapeType = 7;
pub const BOX_RGN: shapeType = 6;
pub const ELLIPTANNULUS_RGN: shapeType = 5;
pub const ELLIPSE_RGN: shapeType = 4;
pub const ANNULUS_RGN: shapeType = 3;
pub const CIRCLE_RGN: shapeType = 2;
pub const LINE_RGN: shapeType = 1;
pub const POINT_RGN: shapeType = 0;
pub type yytype_uint8 = c_uchar;

const PARSER_VECTOR_MIN_ADDR: usize = 0x1000;

#[derive(Copy, Clone)]
#[repr(C)]
pub(crate) union yyalloc {
    pub yyss_alloc: yy_state_t,
    pub yyvs_alloc: FITS_PARSER_YYSTYPE,
}

/* YYFINAL -- State number of the termination state.  */
const YYFINAL: usize = 2;
/* YYLAST -- Last index in YYTABLE.  */
const YYLAST: usize = 1776;

/* YYNTOKENS -- Number of terminals.  */
const YYNTOKENS: usize = 57;
/* YYNNTS -- Number of nonterminals.  */
const YYNNTS: usize = 9;
/* YYNRULES -- Number of rules.  */
const YYNRULES: usize = 135;
/* YYNSTATES -- Number of states.  */
const YYNSTATES: usize = 322;

/* YYMAXUTOK -- Last valid token kind.  */
const YYMAXUTOK: usize = 292;

/* YYTRANSLATE[TOKEN-NUM] -- Symbol number corresponding to TOKEN-NUM as returned by yylex.  */
static YYTRANSLATE: [yytype_int8; 293] = [
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

/* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
static YYRLINE: [yytype_int16; 136] = [
    0, 266, 266, 267, 270, 271, 277, 283, 289, 295, 298, 300, 313, 315, 328, 339, 353, 357, 361,
    365, 367, 376, 379, 382, 391, 393, 395, 397, 399, 401, 404, 408, 410, 412, 414, 423, 425, 427,
    430, 433, 436, 439, 442, 451, 460, 469, 472, 474, 476, 478, 482, 486, 505, 524, 543, 554, 568,
    617, 629, 660, 774, 782, 885, 911, 914, 918, 920, 922, 924, 926, 928, 930, 932, 934, 938, 940,
    942, 951, 954, 957, 960, 963, 966, 969, 972, 975, 978, 981, 984, 987, 990, 993, 996, 999, 1002,
    1005, 1008, 1010, 1012, 1014, 1017, 1024, 1041, 1054, 1067, 1078, 1094, 1118, 1146, 1183, 1187,
    1191, 1194, 1200, 1204, 1208, 1211, 1216, 1220, 1223, 1227, 1229, 1231, 1233, 1235, 1237, 1239,
    1243, 1246, 1248, 1257, 1259, 1261, 1270, 1289, 1308,
];

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
STATE-NUM.  */
static YYPACT: [yytype_int16; 322] = [
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
static YYDEFACT: [yytype_uint8; 322] = [
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
static YYPGOTO: [yytype_int16; 9] = [-41, -41, -41, -41, -41, -1, 170, 96, 30];

/* YYDEFGOTO[NTERM-NUM].  */
static YYDEFGOTO: [yytype_int8; 9] = [0, 1, 31, 32, 33, 48, 49, 46, 63];

/* YYTABLE[YYPACT[STATE-NUM]] -- What to do in state STATE-NUM.  If
positive, shift that token.  If negative, reduce the rule whose
number is the opposite.  If YYTABLE_NINF, syntax error.  */
static YYTABLE: [yytype_int16; 1777] = [
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

static YYCHECK: [yytype_int16; 1777] = [
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
static YYSTOS: [yytype_int8; 322] = [
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
static YYR1: [yytype_int8; 136] = [
    0, 57, 58, 58, 59, 59, 59, 59, 59, 59, 60, 60, 61, 61, 61, 61, 62, 63, 64, 64, 64, 64, 64, 64,
    64, 64, 64, 64, 64, 64, 64, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62,
    62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62,
    62, 62, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63,
    63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63,
    63, 63, 63, 63, 63, 63, 63, 65, 65, 65, 65, 65, 65, 65, 65, 65,
];

/* YYR2[RULE-NUM] -- Number of symbols on the right-hand side of rule RULE-NUM.  */
static YYR2: [yytype_int8; 136] = [
    0, 2, 0, 2, 1, 2, 2, 2, 2, 2, 2, 3, 2, 3, 3, 3, 2, 2, 1, 1, 4, 3, 3, 3, 4, 6, 8, 10, 12, 2, 3,
    1, 1, 1, 4, 1, 1, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 3, 3, 3, 5, 5, 5, 2, 3, 5, 3, 3, 3, 5, 5, 9,
    7, 11, 4, 6, 8, 10, 12, 2, 2, 2, 2, 1, 1, 4, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 3, 3, 5, 5, 3, 3, 3, 5, 7, 11, 15, 2, 3, 5, 9, 2, 3, 5, 9, 3, 7, 9, 4, 6, 8, 10,
    12, 2, 3, 1, 1, 4, 1, 3, 3, 5, 5, 7,
];

fn Alloc_Node(lParse: &mut ParseData) -> c_int {
    // If number of nodes == number allocated, then realloc
    if lParse.nNodes == lParse.nNodesAlloc {
        if lParse.nNodesAlloc == 0 {
            lParse.nNodesAlloc = 100;
            lParse.Nodes = Vec::with_capacity(lParse.nNodesAlloc as usize);
            lParse
                .Nodes
                .resize(lParse.nNodesAlloc as usize, Node::default());
        } else if lParse
            .Nodes
            .try_reserve_exact(lParse.nNodesAlloc as usize)
            .is_err()
        {
            lParse.status = MEMORY_ALLOCATION;
            return -(1);
        } else {
            lParse.nNodesAlloc *= 2;
            lParse
                .Nodes
                .resize(lParse.nNodesAlloc as usize, Node::default());
        }
    }

    lParse.nNodes += 1;
    lParse.nNodes - 1
}

fn Free_Last_Node(lParse: &mut ParseData) {
    if lParse.nNodes != 0 {
        lParse.nNodes -= 1;
        lParse.nNodes;
    }
}

fn New_Const(
    lParse: &mut ParseData,
    returnType: c_int,
    value: *const c_void,
    len: c_long,
) -> c_int {
    unsafe {
        let mut n: c_int = 0;
        n = Alloc_Node(lParse);
        if n >= 0 {
            let this_node = &mut lParse.Nodes[n as usize];
            this_node.operation = CONST_OP; /* Flag a constant */
            this_node.DoOp = None;
            this_node.nSubNodes = 0;
            this_node.ntype = returnType;
            memcpy(
                &mut this_node.value.data as *const _ as *mut c_void,
                value,
                (len as c_ulong).try_into().unwrap(),
            );
            this_node.value.undef = ptr::null_mut();
            this_node.value.nelem = 1;
            this_node.value.naxis = 1;
            this_node.value.naxes[0] = 1;
        }
        n
    }
}

fn New_Column(lParse: &mut ParseData, ColNum: c_int) -> c_int {
    let mut this_node_idx: usize;
    let mut n: c_int = 0;
    let mut i: c_int = 0;

    n = Alloc_Node(lParse);

    if n >= 0 {
        let this_node = &mut lParse.Nodes[n as usize];
        (this_node).operation = -ColNum;
        (this_node).DoOp = None;
        (this_node).nSubNodes = 0;
        (this_node).ntype = ((lParse.varData)[ColNum as usize]).dtype;
        (this_node).value.nelem = ((lParse.varData)[ColNum as usize]).nelem;
        (this_node).value.naxis = ((lParse.varData)[ColNum as usize]).naxis;

        while i < ((lParse.varData)[ColNum as usize]).naxis {
            (this_node).value.naxes[i as usize] =
                ((lParse.varData)[ColNum as usize]).naxes[i as usize];
            i += 1;
        }
    }

    n
}

fn New_Offset(lParse: &mut ParseData, ColNum: c_int, offsetNode: c_int) -> c_int {
    let mut this_node_idx: usize;
    let mut n: c_int = 0;
    let mut i: c_int = 0;
    let mut colNode: c_int = 0;
    colNode = New_Column(lParse, ColNum);
    if colNode < 0 {
        return -(1);
    }

    let colNode = colNode as usize;

    n = Alloc_Node(lParse);

    if n >= 0 {
        let this_node = &mut lParse.Nodes[n as usize];
        (this_node).operation = i32::from(b'{');
        (this_node).DoOp = Some(Do_Offset);
        (this_node).nSubNodes = 2;
        (this_node).SubNodes[0] = colNode;
        (this_node).SubNodes[1] = offsetNode.try_into().unwrap();
        (this_node).ntype = ((lParse.varData)[ColNum as usize]).dtype;
        (this_node).value.nelem = ((lParse.varData)[ColNum as usize]).nelem;
        (this_node).value.naxis = ((lParse.varData)[ColNum as usize]).naxis;
        i = 0;
        while i < ((lParse.varData)[ColNum as usize]).naxis {
            (this_node).value.naxes[i as usize] =
                ((lParse.varData)[ColNum as usize]).naxes[i as usize];
            i += 1;
        }
    }
    n
}

fn New_Unary(lParse: &mut ParseData, returnType: c_int, mut Op: c_int, Node1: c_int) -> c_int {
    let mut this_node_idx: usize;

    let mut i: c_int = 0;
    let mut n: c_int = 0;

    if Node1 < 0 {
        return -(1);
    }

    let that: &mut Node = &mut (lParse.Nodes)[Node1 as usize];

    if Op == 0 {
        Op = returnType;
    }

    if (Op == fits_parser_yytokentype::DOUBLE as c_int
        || Op == fits_parser_yytokentype::FLTCAST as c_int)
        && that.ntype == fits_parser_yytokentype::DOUBLE as c_int
    {
        return Node1;
    }

    if (Op == fits_parser_yytokentype::LONG as c_int
        || Op == fits_parser_yytokentype::INTCAST as c_int)
        && that.ntype == fits_parser_yytokentype::LONG as c_int
    {
        return Node1;
    }

    if Op == fits_parser_yytokentype::BOOLEAN as c_int
        && that.ntype == fits_parser_yytokentype::BOOLEAN as c_int
    {
        return Node1;
    }

    n = Alloc_Node(lParse);

    if n >= 0 {
        let (this_node, that_node) =
            get_this_that_nodes(&mut lParse.Nodes, n as usize, Node1 as usize);

        (this_node).operation = Op;
        (this_node).DoOp = Some(Do_Unary);
        (this_node).nSubNodes = 1;
        (this_node).SubNodes[0] = Node1.try_into().unwrap();
        (this_node).ntype = returnType;

        (this_node).value.nelem = (that_node).value.nelem;
        (this_node).value.naxis = (that_node).value.naxis;

        while i < (that_node).value.naxis {
            (this_node).value.naxes[i as usize] = (that_node).value.naxes[i as usize];
            i += 1;
        }

        if (that_node).operation == CONST_OP {
            ((this_node).DoOp).expect("non-null function pointer")(lParse, n as usize);
        }
    }
    n
}

fn New_BinOp(
    lParse: &mut ParseData,
    returnType: c_int,
    Node1: c_int,
    Op: c_int,
    Node2: c_int,
) -> c_int {
    let mut n: c_int = 0;
    let mut constant: c_int = 0;

    if Node1 < 0 || Node2 < 0 {
        return -(1);
    }

    n = Alloc_Node(lParse);

    if n >= 0 {
        let mut that1_idx = Node1 as usize;
        let that2_idx = Node2 as usize;
        let this_node_idx = n as usize;

        (lParse.Nodes[this_node_idx]).operation = Op;
        (lParse.Nodes[this_node_idx]).nSubNodes = 2;
        (lParse.Nodes[this_node_idx]).SubNodes[0] = Node1.try_into().unwrap();
        (lParse.Nodes[this_node_idx]).SubNodes[1] = Node2.try_into().unwrap();
        (lParse.Nodes[this_node_idx]).ntype = returnType;

        constant = c_int::from(
            (lParse.Nodes[that1_idx]).operation == CONST_OP
                && (lParse.Nodes[that2_idx]).operation == CONST_OP,
        );

        if (lParse.Nodes[that1_idx]).ntype != fits_parser_yytokentype::STRING as c_int
            && (lParse.Nodes[that1_idx]).ntype != fits_parser_yytokentype::BITSTR as c_int
            && Test_Dims(lParse, Node1, Node2) == 0
        {
            Free_Last_Node(lParse);
            fits_parser_yyerror(
                lParse,
                cs!(c"Array sizes/dims do not match for binary operator"),
            );
            return -(1);
        }

        if (lParse.Nodes[that1_idx]).value.nelem == 1 {
            that1_idx = that2_idx;
        }

        (lParse.Nodes[this_node_idx]).value.nelem = (lParse.Nodes[that1_idx]).value.nelem;
        (lParse.Nodes[this_node_idx]).value.naxis = (lParse.Nodes[that1_idx]).value.naxis;

        for i in 0..lParse.Nodes[that1_idx].value.naxis as usize {
            (lParse.Nodes[this_node_idx]).value.naxes[i] = (lParse.Nodes[that1_idx]).value.naxes[i];
        }

        if Op == fits_parser_yytokentype::ACCUM as c_int
            && (lParse.Nodes[that1_idx]).ntype == fits_parser_yytokentype::BITSTR as c_int
        {
            /* ACCUM is rank-reducing on bit strings */
            (lParse.Nodes[this_node_idx]).value.nelem = 1;
            (lParse.Nodes[this_node_idx]).value.naxis = 1;
            (lParse.Nodes[this_node_idx]).value.naxes[0] = 1;
        }

        /*  Both subnodes should be of same time  */
        match (lParse.Nodes[that1_idx]).ntype.into() {
            fits_parser_yytokentype::BITSTR => {
                (lParse.Nodes[this_node_idx]).DoOp = Some(Do_BinOp_bit);
            }
            fits_parser_yytokentype::STRING => {
                (lParse.Nodes[this_node_idx]).DoOp = Some(Do_BinOp_str);
            }
            fits_parser_yytokentype::BOOLEAN => {
                (lParse.Nodes[this_node_idx]).DoOp = Some(Do_BinOp_log);
            }
            fits_parser_yytokentype::LONG => {
                (lParse.Nodes[this_node_idx]).DoOp = Some(Do_BinOp_lng);
            }
            fits_parser_yytokentype::DOUBLE => {
                (lParse.Nodes[this_node_idx]).DoOp = Some(Do_BinOp_dbl);
            }
            _ => {}
        }

        if constant != 0 {
            ((lParse.Nodes[this_node_idx]).DoOp).expect("non-null function pointer")(
                lParse, n as usize,
            );
        }
    }
    n
}

fn New_Func(
    lParse: &mut ParseData,
    returnType: c_int,
    Op: funcOp,
    nNodes: c_int,
    Node1: c_int,
    Node2: c_int,
    Node3: c_int,
    Node4: c_int,
    Node5: c_int,
    Node6: c_int,
    Node7: c_int,
) -> c_int {
    New_FuncSize(
        lParse, returnType, Op, nNodes, Node1, Node2, Node3, Node4, Node5, Node6, Node7, 0,
    )
}

macro_rules! YY_SYMBOL_PRINT {
    ($title:expr, $kind:expr, $value:expr, $scanner:expr, $lParse:expr) => {{
        if YYDEBUG {
            eprintln!("{}", $title);
            yy_symbol_print($kind, $value, $scanner, $lParse);
            eprintln!();
        }
    }};
}

fn yysymbol_name(yykind: yysymbol_kind_t) -> &'static str {
    YYTNAME[yykind as usize]
}

/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static YYTNAME: [&str; 67] = [
    "\"end of file\"",
    "error",
    "\"invalid token\"",
    "BOOLEAN",
    "LONG",
    "DOUBLE",
    "STRING",
    "BITSTR",
    "FUNCTION",
    "BFUNCTION",
    "IFUNCTION",
    "GTIFILTER",
    "GTIOVERLAP",
    "GTIFIND",
    "REGFILTER",
    "COLUMN",
    "BCOLUMN",
    "SCOLUMN",
    "BITCOL",
    "ROWREF",
    "NULLREF",
    "SNULLREF",
    "','",
    "'='",
    "':'",
    "'{'",
    "'}'",
    "'?'",
    "OR",
    "AND",
    "EQ",
    "NE",
    "'~'",
    "GT",
    "LT",
    "LTE",
    "GTE",
    "'+'",
    "'-'",
    "'%'",
    "'*'",
    "'/'",
    "'|'",
    "'&'",
    "XOR",
    "POWER",
    "NOT",
    "INTCAST",
    "FLTCAST",
    "UMINUS",
    "'['",
    "ACCUM",
    "DIFF",
    "'\\n'",
    "']'",
    "'('",
    "')'",
    "$accept",
    "lines",
    "line",
    "bvector",
    "vector",
    "expr",
    "bexpr",
    "bits",
    "sexpr",
    "NULL PTR",
];

/*-----------------------------------.
| Print this symbol's value on YYO.  |
`-----------------------------------*/

fn yy_symbol_value_print(
    yykind: yysymbol_kind_t,
    yyvaluep: *const FITS_PARSER_YYSTYPE,
    _scanner: yyscan_t,
    _lParse: &mut ParseData,
) {
}

/*---------------------------.
| Print this symbol on YYO.  |
`---------------------------*/

fn yy_symbol_print(
    yykind: yysymbol_kind_t,
    yyvaluep: *const FITS_PARSER_YYSTYPE,
    scanner: yyscan_t,
    lParse: &mut ParseData,
) {
    eprint!(
        "{} {} (",
        if (yykind as usize) < YYNTOKENS {
            "token"
        } else {
            "nterm"
        },
        yysymbol_name(yykind)
    );
    yy_symbol_value_print(yykind, yyvaluep, scanner, lParse);
    eprintln!(")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

fn yy_stackprint(yybottom: &[yy_state_t], yytop: usize) {
    eprint!("Stack now");
    let mut ptr = 0;
    while ptr <= yytop {
        let yybot = yybottom[ptr];
        eprint!(" {}", yybot);
        ptr += 1;
    }
    eprintln!();
}

macro_rules! YY_STACK_PRINT {
    ($bottom:expr, $top:expr) => {{
        if YYDEBUG {
            yy_stackprint($bottom, $top);
        }
    }};
}

/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

fn yy_reduce_print(
    yyssp: *const yy_state_t,
    yyvsp: *const FITS_PARSER_YYSTYPE,
    yyrule: c_int,
    scanner: yyscan_t,
    lParse: &mut ParseData,
) {
    unsafe {
        let yylno = YYRLINE[yyrule as usize];
        let yynrhs = YYR2[yyrule as usize];
        eprintln!("Reducing stack by rule {} (line {}):", yyrule - 1, yylno);
        /* The symbols being reduced.  */
        for yyi in 0..yynrhs {
            eprint!("   ${} = ", yyi + 1);
            yy_symbol_print(
                (YYSTOS[(*yyssp.add((yyi + 1 - yynrhs) as usize)) as usize]).into(),
                yyvsp.add((yyi + 1 - yynrhs) as usize),
                scanner,
                lParse,
            );
            eprintln!();
        }
    }
}

const YYDEBUG: bool = false;

/* YYINITDEPTH -- initial size of the parser's stacks.  */
const YYINITDEPTH: usize = 200;

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
if the built-in stack extension method is used).

Do not make this value too large; the results are undefined if
YYSTACK_ALLOC_MAXIMUM < YYSTACK_BYTES (YYMAXDEPTH)
evaluated with infinite-precision integer arithmetic.  */

const YYMAXDEPTH: usize = 10000;

/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

fn yydestruct(
    mut yymsg: *const c_char,
    yykind: yysymbol_kind_t,
    yyvaluep: *mut FITS_PARSER_YYSTYPE,
    scanner: yyscan_t,
    lParse: &mut ParseData,
) {
    if yymsg.is_null() {
        yymsg = c"Deleting".as_ptr();
    }
}

/* If returnType==0 , use Node1's type and vector sizes as returnType, */
/* else return a single value of type returnType                       */
fn New_FuncSize(
    lParse: &mut ParseData,
    returnType: c_int,
    Op: funcOp,
    nNodes: c_int,
    Node1: c_int,
    Node2: c_int,
    Node3: c_int,
    Node4: c_int,
    Node5: c_int,
    Node6: c_int,
    Node7: c_int,
    Size: c_int,
) -> c_int {
    let mut this_node_idx: usize;
    let mut i: c_int = 0;
    let mut n: c_int = 0;
    let mut constant: c_int = 0;

    if Node1 < 0 || Node2 < 0 || Node3 < 0 || Node4 < 0 || Node5 < 0 || Node6 < 0 || Node7 < 0 {
        return -(1);
    }

    n = Alloc_Node(lParse);
    if n >= 0 {
        let this_node = &mut lParse.Nodes[n as usize];

        (this_node).operation = Op as c_int;
        (this_node).DoOp = Some(Do_Func);
        (this_node).nSubNodes = nNodes;
        (this_node).SubNodes[0] = Node1.try_into().unwrap();
        (this_node).SubNodes[1] = Node2.try_into().unwrap();
        (this_node).SubNodes[2] = Node3.try_into().unwrap();
        (this_node).SubNodes[3] = Node4.try_into().unwrap();
        (this_node).SubNodes[4] = Node5.try_into().unwrap();
        (this_node).SubNodes[5] = Node6.try_into().unwrap();
        (this_node).SubNodes[6] = Node7.try_into().unwrap();

        constant = nNodes;

        i = constant; /* Functions with zero params are not const */

        if Op as c_uint == POIRND_FCT as c_int as c_uint {
            constant = 0; /* Nor is Poisson deviate */
        }

        loop {
            let fresh1 = i;
            i -= 1;
            if fresh1 == 0 {
                break;
            }
            constant = c_int::from(
                constant != 0
                    && ((lParse.Nodes)[(lParse.Nodes[n as usize]).SubNodes[i as usize] as usize])
                        .operation
                        == CONST_OP,
            );
        }

        if returnType != 0 {
            (lParse.Nodes[n as usize]).ntype = returnType;
            (lParse.Nodes[n as usize]).value.nelem = 1;
            (lParse.Nodes[n as usize]).value.naxis = 1;
            (lParse.Nodes[n as usize]).value.naxes[0] = 1;
        } else {
            (lParse.Nodes[n as usize]).ntype = (lParse.Nodes[Node1 as usize]).ntype;
            (lParse.Nodes[n as usize]).value.nelem = (lParse.Nodes[Node1 as usize]).value.nelem;
            (lParse.Nodes[n as usize]).value.naxis = (lParse.Nodes[Node1 as usize]).value.naxis;
            i = 0;
            while i < (lParse.Nodes[Node1 as usize]).value.naxis {
                (lParse.Nodes[n as usize]).value.naxes[i as usize] =
                    (lParse.Nodes[Node1 as usize]).value.naxes[i as usize];
                i += 1;
            }
        }

        if Size > 0 {
            (lParse.Nodes[n as usize]).value.nelem = c_long::from(Size);
        }

        if constant != 0 {
            ((lParse.Nodes[n as usize]).DoOp).expect("non-null function pointer")(
                lParse, n as usize,
            );
        }
    }
    n
}

pub(crate) fn fits_parser_yyparse(scanner: &mut yyguts_t, lParse: &mut ParseData) -> c_int {
    unsafe {
        let mut current_block: u64;
        let mut yychar: c_int = 0;
        static mut YYVAL_DEFAULT: FITS_PARSER_YYSTYPE = FITS_PARSER_YYSTYPE { Node: 0 };
        let mut yylval: FITS_PARSER_YYSTYPE = YYVAL_DEFAULT;
        let mut fits_parser_yynerrs: c_int = 0;
        let mut yystate: yy_state_fast_t = 0;
        let mut yyerrstatus: c_int = 0;
        let mut yystacksize: c_long = 100 as c_long;
        let mut yyssa: [yy_state_t; 100] = [0; 100];
        let mut yyss: &mut [yy_state_t] = &mut yyssa;
        let mut yyssp: usize = 0;
        let mut yyvsa: [FITS_PARSER_YYSTYPE; 100] = [FITS_PARSER_YYSTYPE { Node: 0 }; 100];
        let mut yyvs: &mut [FITS_PARSER_YYSTYPE] = &mut yyvsa;
        let mut yyvsp: usize = 0;
        let mut yyn: c_int = 0;
        let mut yyresult: c_int = 0;
        let mut yytoken: yysymbol_kind_t = YYSYMBOL_YYEMPTY;
        let mut yyval: FITS_PARSER_YYSTYPE = FITS_PARSER_YYSTYPE { Node: 0 };
        let mut yylen: c_int = 0;
        yychar = fits_parser_yytokentype::FITS_PARSER_YYEMPTY as c_int;
        's_54: loop {
            /*--------------------------------------------------------------------.
            | yysetstate -- set current state (the top of the stack) to yystate.  |
            `--------------------------------------------------------------------*/
            if YYDEBUG {
                eprintln!("Entering state {}", yystate);
            }
            assert!(0 <= yystate && yystate < YYNSTATES as c_int);
            YY_STACK_PRINT!(yyss, yyssp);

            yyss[yyssp] = yystate as yy_state_t;
            if (yystacksize - 1) as usize <= yyssp {
                /* Get the current used size of the three stacks, in elements.  */
                let yysize: c_long = (yyssp + 1) as c_long;
                if YYMAXDEPTH as c_long <= yystacksize {
                    current_block = 11794367917084412820; // goto yyexhaustedlab;
                    break;
                }
                yystacksize *= 2 as c_long;
                if (YYMAXDEPTH as c_long) < yystacksize {
                    /* Extend the stack our own way.  */
                    yystacksize = YYMAXDEPTH as c_long;
                }
                let yyss1: *mut yy_state_t = yyss.as_mut_ptr();
                let mut yyptr: *mut yyalloc = malloc(
                    ((yystacksize
                        * (::core::mem::size_of::<yy_state_t>() as c_ulong as c_long
                            + ::core::mem::size_of::<FITS_PARSER_YYSTYPE>() as c_ulong as c_long)
                        + (::core::mem::size_of::<yyalloc>() as c_ulong as c_long - 1))
                        as c_ulong)
                        .try_into()
                        .unwrap(),
                )
                .cast::<yyalloc>();
                if yyptr.is_null() {
                    current_block = 11794367917084412820; // goto yyexhaustedlab;
                    break;
                }
                let mut yynewbytes: c_long = 0;
                libc::memcpy(
                    (&mut (*yyptr).yyss_alloc as *mut yy_state_t).cast::<c_void>(),
                    yyss.as_ptr().cast::<c_void>(),
                    (yysize as c_ulong)
                        .wrapping_mul(::core::mem::size_of::<yy_state_t>() as c_ulong)
                        as libc::size_t,
                );
                yyss = slice::from_raw_parts_mut(&mut (*yyptr).yyss_alloc, yystacksize as usize);
                yynewbytes = yystacksize
                    * ::core::mem::size_of::<yy_state_t>() as c_ulong as c_long
                    + (::core::mem::size_of::<yyalloc>() as c_ulong as c_long - 1);
                yyptr = yyptr.offset(
                    (yynewbytes / ::core::mem::size_of::<yyalloc>() as c_ulong as c_long) as isize,
                );
                let mut yynewbytes_0: c_long = 0;
                libc::memcpy(
                    (&mut (*yyptr).yyvs_alloc as *mut FITS_PARSER_YYSTYPE).cast::<c_void>(),
                    yyvs.as_ptr().cast::<c_void>(),
                    (yysize as c_ulong)
                        .wrapping_mul(::core::mem::size_of::<FITS_PARSER_YYSTYPE>() as c_ulong)
                        as libc::size_t,
                );
                yyvs = slice::from_raw_parts_mut(&mut (*yyptr).yyvs_alloc, yystacksize as usize);
                yynewbytes_0 = yystacksize
                    * ::core::mem::size_of::<FITS_PARSER_YYSTYPE>() as c_ulong as c_long
                    + (::core::mem::size_of::<yyalloc>() as c_ulong as c_long - 1);
                yyptr = yyptr.offset(
                    (yynewbytes_0 / ::core::mem::size_of::<yyalloc>() as c_ulong as c_long)
                        as isize,
                );
                if yyss1 != yyssa.as_mut_ptr() {
                    free(yyss1.cast::<c_void>());
                }
                yyssp = yysize as usize - 1;
                yyvsp = yysize as usize - 1;
                if (yystacksize - 1) as usize <= yyssp {
                    current_block = 3964311021479492664; // goto yyabortlab;
                    break;
                }
            }
            if yystate == YYFINAL as c_int {
                yyresult = 0;
                current_block = 15864720325503947191; // goto yyacceptlab;
                break;
            } else {
                /*-----------.
                | yybackup.  |
                `-----------*/
                /* Do appropriate processing given the current state.  Read a
                lookahead token if we need one and don't already have one.  */

                /* First try to decide what to do without reference to lookahead token.  */
                yyn = c_int::from(YYPACT[yystate as usize]);
                if yyn == -41 {
                    current_block = 5937473999264333383; // goto yydefault;
                } else {
                    /* Not known => get a lookahead token if don't already have one.  */

                    /* YYCHAR is either empty, or end-of-input, or a valid lookahead.  */

                    if yychar == fits_parser_yytokentype::FITS_PARSER_YYEMPTY as c_int {
                        if YYDEBUG {
                            eprintln!("Read a token");
                        }
                        yychar = fits_parser_yylex(&mut yylval, scanner);
                    }

                    if yychar <= fits_parser_yytokentype::FITS_PARSER_YYEOF as c_int {
                        yychar = fits_parser_yytokentype::FITS_PARSER_YYEOF as c_int;
                        yytoken = YYSYMBOL_YYEOF;
                        if YYDEBUG {
                            eprintln!("Now at end of input.");
                        }
                        current_block = 1924505913685386279;
                    } else if yychar == fits_parser_yytokentype::FITS_PARSER_YYerror as c_int {
                        /* The scanner already issued an error message, process directly
                        to error recovery.  But do not keep the error token as
                        lookahead, it is too special and may lead us to an endless
                        loop in error recovery. */
                        yychar = fits_parser_yytokentype::FITS_PARSER_YYUNDEF as c_int;
                        yytoken = YYSYMBOL_YYERROR;
                        current_block = 1774893048582444437; // goto yyerrlab1;
                    } else {
                        yytoken = (if 0 <= yychar && yychar <= 292 as c_int {
                            yysymbol_kind_t::from(YYTRANSLATE[yychar as usize]) as c_int
                        } else {
                            YYSYMBOL_YYUNDEF as c_int
                        }) as yysymbol_kind_t;
                        YY_SYMBOL_PRINT!("Next token is", yytoken, &yylval, scanner, lParse);
                        current_block = 1924505913685386279;
                    }
                    match current_block {
                        1774893048582444437 => {}
                        _ => {
                            yyn += yytoken as c_int;
                            if yyn < 0
                                || (YYLAST as c_int) < yyn
                                || c_int::from(YYCHECK[yyn as usize]) != yytoken as c_int
                            {
                                current_block = 5937473999264333383;
                            } else {
                                yyn = c_int::from(YYTABLE[yyn as usize]);
                                if yyn <= 0 {
                                    yyn = -yyn;
                                    current_block = 670225253387957849;
                                } else {
                                    /* Count tokens shifted since error; after three, turn off error status.  */
                                    if yyerrstatus != 0 {
                                        yyerrstatus -= 1;
                                    } /* Shift the lookahead token.  */
                                    YY_SYMBOL_PRINT!("Shifting", yytoken, &yylval, scanner, lParse);
                                    yystate = yyn;
                                    yyvsp += 1;
                                    yyvs[yyvsp] = yylval;
                                    yychar = fits_parser_yytokentype::FITS_PARSER_YYEMPTY as c_int; /* Discard the shifted token.  */
                                    current_block = 7872030484262409139;
                                }
                            }
                        }
                    }
                }

                if current_block == 5937473999264333383 {
                    /*-----------------------------------------------------------.
                    | yydefault -- do the default action for the current state.  |
                    `-----------------------------------------------------------*/

                    yyn = c_int::from(YYDEFACT[yystate as usize]);
                    if yyn == 0 {
                        // goto yyerrlab;
                        yytoken =
                            (if yychar == fits_parser_yytokentype::FITS_PARSER_YYEMPTY as c_int {
                                YYSYMBOL_YYEMPTY as c_int
                            } else if 0 <= yychar && yychar <= 292 as c_int {
                                yysymbol_kind_t::from(YYTRANSLATE[yychar as usize]) as c_int
                            } else {
                                YYSYMBOL_YYUNDEF as c_int
                            }) as yysymbol_kind_t;
                        if yyerrstatus == 0 {
                            fits_parser_yynerrs += 1;
                            fits_parser_yyerror(lParse, cs!(c"syntax error"));
                        }
                        if yyerrstatus == 3 as c_int {
                            if yychar <= fits_parser_yytokentype::FITS_PARSER_YYEOF as c_int {
                                if yychar == fits_parser_yytokentype::FITS_PARSER_YYEOF as c_int {
                                    current_block = 3964311021479492664;
                                    break;
                                }
                            } else {
                                yydestruct(
                                    c"Error: discarding".as_ptr(),
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
                        current_block = 670225253387957849; // goto yyreduce;
                    }
                }
                if current_block == 670225253387957849 {
                    /*-----------------------------.
                    | yyreduce -- do a reduction.  |
                    `-----------------------------*/

                    /* yyn is the number of a rule to reduce with.  */
                    yylen = c_int::from(YYR2[yyn as usize]);

                    /* If YYLEN is nonzero, implement the default value of the action:
                    '$$ = $1'.

                    Otherwise, the following line sets YYVAL to garbage.
                    This behavior is undocumented and Bison
                    users should not rely upon it.  Assigning to YYVAL
                    unconditionally makes the parser a bit smaller, and it avoids a
                    GCC warning that YYVAL may be used uninitialized.  */
                    yyval = yyvs[yyvsp + 1 - yylen as usize];

                    if YYDEBUG {
                        yy_reduce_print(&yyss[yyssp], &yyvs[yyvsp], yyn, scanner, lParse);
                    }

                    match yyn {
                        4 => {
                            /* line: '\n'  */
                            current_block = 17353983478346836848;
                        }
                        5 => {
                            /* line: expr '\n'  */
                            if yyvs[yyvsp - 1].Node < 0 {
                                fits_parser_yyerror(
                                    lParse,
                                    cs!(c"Couldn't build node structure: out of memory?"),
                                );
                                current_block = 4830776507462815627;
                            } else {
                                lParse.resultNode = yyvs[yyvsp - 1].Node;
                                current_block = 17353983478346836848;
                            }
                        }
                        6 => {
                            /* line: bexpr '\n'  */
                            if yyvs[yyvsp - 1].Node < 0 {
                                fits_parser_yyerror(
                                    lParse,
                                    cs!(c"Couldn't build node structure: out of memory?"),
                                );
                                current_block = 4830776507462815627;
                            } else {
                                lParse.resultNode = yyvs[yyvsp - 1].Node;
                                current_block = 17353983478346836848;
                            }
                        }
                        7 => {
                            /* line: sexpr '\n'  */
                            if yyvs[yyvsp - 1].Node < 0 {
                                fits_parser_yyerror(
                                    lParse,
                                    cs!(c"Couldn't build node structure: out of memory?"),
                                );
                                current_block = 4830776507462815627;
                            } else {
                                lParse.resultNode = yyvs[yyvsp - 1].Node;
                                current_block = 17353983478346836848;
                            }
                        }
                        8 => {
                            /* line: bits '\n'  */
                            if yyvs[yyvsp - 1].Node < 0 {
                                fits_parser_yyerror(
                                    lParse,
                                    cs!(c"Couldn't build node structure: out of memory?"),
                                );
                                current_block = 4830776507462815627;
                            } else {
                                lParse.resultNode = yyvs[yyvsp - 1].Node;
                                current_block = 17353983478346836848;
                            }
                        }
                        9 => {
                            /* line: error '\n'  */
                            yyerrstatus = 0;
                            current_block = 17353983478346836848;
                        }
                        10 => {
                            /* bvector: '{' bexpr  */
                            yyval.Node = New_Vector(lParse, yyvs[yyvsp].Node);
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        11 => {
                            /* bvector: bvector ',' bexpr  */
                            if (lParse.Nodes[yyvs[yyvsp - 2].Node as usize]).nSubNodes
                                >= 10 as c_int
                            {
                                yyvs[yyvsp - 2].Node = Close_Vec(lParse, yyvs[yyvsp - 2].Node);
                                if yyvs[yyvsp - 2].Node < 0 {
                                    current_block = 4830776507462815627;
                                } else {
                                    yyval.Node = New_Vector(lParse, yyvs[yyvsp - 2].Node);
                                    if yyval.Node < 0 {
                                        current_block = 4830776507462815627;
                                    } else {
                                        current_block = 2616667235040759262;
                                    }
                                }
                            } else {
                                yyval.Node = yyvs[yyvsp - 2].Node;
                                current_block = 2616667235040759262;
                            }
                            match current_block {
                                4830776507462815627 => {}
                                _ => {
                                    let fresh2 =
                                        &mut ((lParse.Nodes)[yyval.Node as usize]).nSubNodes;
                                    let fresh3 = *fresh2;
                                    *fresh2 += 1;
                                    ((lParse.Nodes)[yyval.Node as usize]).SubNodes
                                        [fresh3 as usize] = yyvs[yyvsp].Node.try_into().unwrap();
                                    current_block = 17353983478346836848;
                                }
                            }
                        }
                        12 => {
                            /* vector: '{' expr  */
                            yyval.Node = New_Vector(lParse, yyvs[yyvsp].Node);
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        13 => {
                            /* vector: vector ',' expr  */
                            if (lParse.Nodes[yyvs[yyvsp - 2].Node as usize]).ntype
                                < (lParse.Nodes[yyvs[yyvsp].Node as usize]).ntype
                            {
                                (lParse.Nodes[yyvs[yyvsp - 2].Node as usize]).ntype =
                                    ((lParse.Nodes)[yyvs[yyvsp].Node as usize]).ntype;
                            }
                            if (lParse.Nodes[yyvs[yyvsp - 2].Node as usize]).nSubNodes
                                >= 10 as c_int
                            {
                                yyvs[yyvsp - 2].Node = Close_Vec(lParse, yyvs[yyvsp - 2].Node);
                                if yyvs[yyvsp - 2].Node < 0 {
                                    current_block = 4830776507462815627;
                                } else {
                                    yyval.Node = New_Vector(lParse, yyvs[yyvsp - 2].Node);
                                    if yyval.Node < 0 {
                                        current_block = 4830776507462815627;
                                    } else {
                                        current_block = 2904036176499606090;
                                    }
                                }
                            } else {
                                yyval.Node = yyvs[yyvsp - 2].Node;
                                current_block = 2904036176499606090;
                            }
                            match current_block {
                                4830776507462815627 => {}
                                _ => {
                                    let fresh4 =
                                        &mut ((lParse.Nodes)[yyval.Node as usize]).nSubNodes;
                                    let fresh5 = *fresh4;
                                    *fresh4 += 1;
                                    ((lParse.Nodes)[yyval.Node as usize]).SubNodes
                                        [fresh5 as usize] = yyvs[yyvsp].Node.try_into().unwrap();
                                    current_block = 17353983478346836848;
                                }
                            }
                        }
                        14 => {
                            /* vector: vector ',' bexpr  */
                            if (lParse.Nodes[yyvs[yyvsp - 2].Node as usize]).nSubNodes
                                >= MAXSUBS as c_int
                            {
                                yyvs[yyvsp - 2].Node = Close_Vec(lParse, yyvs[yyvsp - 2].Node);
                                if yyvs[yyvsp - 2].Node < 0 {
                                    current_block = 4830776507462815627;
                                } else {
                                    yyval.Node = New_Vector(lParse, yyvs[yyvsp - 2].Node);
                                    if yyval.Node < 0 {
                                        current_block = 4830776507462815627;
                                    } else {
                                        current_block = 17702298541784679949;
                                    }
                                }
                            } else {
                                yyval.Node = yyvs[yyvsp - 2].Node;
                                current_block = 17702298541784679949;
                            }
                            match current_block {
                                4830776507462815627 => {}
                                _ => {
                                    let fresh6 =
                                        &mut ((lParse.Nodes)[yyval.Node as usize]).nSubNodes;
                                    let fresh7 = *fresh6;
                                    *fresh6 += 1;
                                    ((lParse.Nodes)[yyval.Node as usize]).SubNodes
                                        [fresh7 as usize] = yyvs[yyvsp].Node.try_into().unwrap();
                                    current_block = 17353983478346836848;
                                }
                            }
                        }
                        15 => {
                            /* vector: bvector ',' expr  */
                            (lParse.Nodes[yyvs[yyvsp - 2].Node as usize]).ntype =
                                (lParse.Nodes[yyvs[yyvsp].Node as usize]).ntype;
                            if (lParse.Nodes[yyvs[yyvsp - 2].Node as usize]).nSubNodes
                                >= 10 as c_int
                            {
                                yyvs[yyvsp - 2].Node = Close_Vec(lParse, yyvs[yyvsp - 2].Node);
                                if yyvs[yyvsp - 2].Node < 0 {
                                    current_block = 4830776507462815627;
                                } else {
                                    yyval.Node = New_Vector(lParse, yyvs[yyvsp - 2].Node);
                                    if yyval.Node < 0 {
                                        current_block = 4830776507462815627;
                                    } else {
                                        current_block = 1069630499025798221;
                                    }
                                }
                            } else {
                                yyval.Node = yyvs[yyvsp - 2].Node;
                                current_block = 1069630499025798221;
                            }
                            match current_block {
                                4830776507462815627 => {}
                                _ => {
                                    let fresh8 =
                                        &mut ((lParse.Nodes)[yyval.Node as usize]).nSubNodes;
                                    let fresh9 = *fresh8;
                                    *fresh8 += 1;
                                    ((lParse.Nodes)[yyval.Node as usize]).SubNodes
                                        [fresh9 as usize] = yyvs[yyvsp].Node.try_into().unwrap();
                                    current_block = 17353983478346836848;
                                }
                            }
                        }
                        16 => {
                            /* expr: vector '}'  */
                            yyval.Node = Close_Vec(lParse, yyvs[yyvsp - 1].Node);
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        17 => {
                            /* bexpr: bvector '}'  */
                            yyval.Node = Close_Vec(lParse, yyvs[yyvsp - 1].Node);
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        18 => {
                            /* bits: BITSTR  */
                            yyval.Node = New_Const(
                                lParse,
                                fits_parser_yytokentype::BITSTR as c_int,
                                (yyvs[yyvsp].astr).as_mut_ptr().cast::<c_void>(),
                                (strlen((yyvs[yyvsp].astr).as_mut_ptr()))
                                    .wrapping_add((1).try_into().unwrap())
                                    as c_long,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                (((lParse.Nodes)[yyval.Node as usize]).value).nelem =
                                    strlen((yyvs[yyvsp].astr).as_mut_ptr()) as c_long;
                                current_block = 17353983478346836848;
                            }
                        }
                        19 => {
                            /* bits: BITCOL  */
                            yyval.Node = New_Column(lParse, yyvs[yyvsp].lng as c_int);
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        20 => {
                            /* bits: BITCOL '{' expr '}'  */
                            if ((lParse.Nodes)[yyvs[yyvsp - 1].Node as usize]).ntype
                                != fits_parser_yytokentype::LONG as c_int
                                || ((lParse.Nodes)[yyvs[yyvsp - 1].Node as usize]).operation
                                    != CONST_OP
                            {
                                fits_parser_yyerror(
                                    lParse,
                                    cs!(c"Offset argument must be a constant integer"),
                                );
                                current_block = 4830776507462815627;
                            } else {
                                yyval.Node = New_Offset(
                                    lParse,
                                    yyvs[yyvsp - 3].lng as c_int,
                                    yyvs[yyvsp - 1].Node,
                                );
                                if yyval.Node < 0 {
                                    current_block = 4830776507462815627;
                                } else {
                                    current_block = 17353983478346836848;
                                }
                            }
                        }
                        21 => {
                            /* bits: bits '&' bits  */
                            yyval.Node = New_BinOp(
                                lParse,
                                fits_parser_yytokentype::BITSTR as c_int,
                                yyvs[yyvsp - 2].Node,
                                '&' as i32,
                                yyvs[yyvsp].Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                (((lParse.Nodes)[yyval.Node as usize]).value).nelem =
                                    if (((lParse.Nodes)[yyvs[yyvsp - 2].Node as usize]).value).nelem
                                        > (lParse.Nodes[yyvs[yyvsp].Node as usize]).value.nelem
                                    {
                                        ((lParse.Nodes)[yyvs[yyvsp - 2].Node as usize]).value.nelem
                                    } else {
                                        (lParse.Nodes[yyvs[yyvsp].Node as usize]).value.nelem
                                    };
                                current_block = 17353983478346836848;
                            }
                        }
                        22 => {
                            /* bits: bits '|' bits  */
                            yyval.Node = New_BinOp(
                                lParse,
                                fits_parser_yytokentype::BITSTR as c_int,
                                yyvs[yyvsp - 2].Node,
                                '|' as i32,
                                yyvs[yyvsp].Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                (((lParse.Nodes)[yyval.Node as usize]).value).nelem =
                                    if (((lParse.Nodes)[yyvs[yyvsp - 2].Node as usize]).value).nelem
                                        > (lParse.Nodes[yyvs[yyvsp].Node as usize]).value.nelem
                                    {
                                        ((lParse.Nodes)[yyvs[yyvsp - 2].Node as usize]).value.nelem
                                    } else {
                                        (lParse.Nodes[yyvs[yyvsp].Node as usize]).value.nelem
                                    };
                                current_block = 17353983478346836848;
                            }
                        }
                        23 => {
                            /* bits: bits '+' bits  */
                            if (lParse.Nodes[yyvs[yyvsp - 2].Node as usize]).value.nelem
                                + (lParse.Nodes[yyvs[yyvsp].Node as usize]).value.nelem
                                >= 256 as c_long
                            {
                                fits_parser_yyerror(
                                    lParse,
                                    cs!(c"Combined bit string size exceeds 255 bits"),
                                );
                                current_block = 4830776507462815627;
                            } else {
                                yyval.Node = New_BinOp(
                                    lParse,
                                    fits_parser_yytokentype::BITSTR as c_int,
                                    yyvs[yyvsp - 2].Node,
                                    '+' as i32,
                                    yyvs[yyvsp].Node,
                                );
                                if yyval.Node < 0 {
                                    current_block = 4830776507462815627;
                                } else {
                                    (((lParse.Nodes)[yyval.Node as usize]).value).nelem =
                                        ((lParse.Nodes)[yyvs[yyvsp - 2].Node as usize]).value.nelem
                                            + ((lParse.Nodes)[yyvs[yyvsp].Node as usize])
                                                .value
                                                .nelem;
                                    current_block = 17353983478346836848;
                                }
                            }
                        }
                        24 => {
                            /* bits: bits '[' expr ']'  */
                            yyval.Node = New_Deref(
                                lParse,
                                yyvs[yyvsp - 3].Node,
                                1,
                                yyvs[yyvsp - 1].Node,
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
                            /* bits: bits '[' expr ',' expr ']'  */
                            yyval.Node = New_Deref(
                                lParse,
                                yyvs[yyvsp - 5].Node,
                                2,
                                yyvs[yyvsp - 3].Node,
                                yyvs[yyvsp - 1].Node,
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
                            /* bits: bits '[' expr ',' expr ',' expr ']'  */
                            yyval.Node = New_Deref(
                                lParse,
                                yyvs[yyvsp - 7].Node,
                                3 as c_int,
                                yyvs[yyvsp - 5].Node,
                                yyvs[yyvsp - 3].Node,
                                yyvs[yyvsp - 1].Node,
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
                            /* bits: bits '[' expr ',' expr ',' expr ',' expr ']'  */
                            yyval.Node = New_Deref(
                                lParse,
                                yyvs[yyvsp - 9].Node,
                                4 as c_int,
                                yyvs[yyvsp - 7].Node,
                                yyvs[yyvsp - 5].Node,
                                yyvs[yyvsp - 3].Node,
                                yyvs[yyvsp - 1].Node,
                                0,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        28 => {
                            /* bits: bits '[' expr ',' expr ',' expr ',' expr ',' expr ']'  */
                            yyval.Node = New_Deref(
                                lParse,
                                yyvs[yyvsp - 11].Node,
                                5 as c_int,
                                yyvs[yyvsp - 9].Node,
                                yyvs[yyvsp - 7].Node,
                                yyvs[yyvsp - 5].Node,
                                yyvs[yyvsp - 3].Node,
                                yyvs[yyvsp - 1].Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        29 => {
                            /* bits: NOT bits  */
                            yyval.Node = New_Unary(
                                lParse,
                                fits_parser_yytokentype::BITSTR as c_int,
                                fits_parser_yytokentype::NOT as c_int,
                                yyvs[yyvsp].Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        30 => {
                            /* bits: '(' bits ')'  */
                            yyval.Node = yyvs[yyvsp - 1].Node;
                            current_block = 17353983478346836848;
                        }
                        31 => {
                            /* expr: LONG  */
                            yyval.Node = New_Const(
                                lParse,
                                fits_parser_yytokentype::LONG as c_int,
                                (&mut yyvs[yyvsp].lng as *mut c_long).cast::<c_void>(),
                                ::core::mem::size_of::<c_long>() as c_ulong as c_long,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        32 => {
                            /* expr: DOUBLE  */
                            yyval.Node = New_Const(
                                lParse,
                                fits_parser_yytokentype::DOUBLE as c_int,
                                (&mut yyvs[yyvsp].dbl as *mut c_double).cast::<c_void>(),
                                ::core::mem::size_of::<c_double>() as c_ulong as c_long,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        33 => {
                            /* expr: COLUMN  */
                            yyval.Node = New_Column(lParse, yyvs[yyvsp].lng as c_int);
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        34 => {
                            /* expr: COLUMN '{' expr '}'  */
                            if ((lParse.Nodes)[yyvs[yyvsp - 1].Node as usize]).ntype
                                != fits_parser_yytokentype::LONG as c_int
                                || ((lParse.Nodes)[yyvs[yyvsp - 1].Node as usize]).operation
                                    != CONST_OP
                            {
                                fits_parser_yyerror(
                                    lParse,
                                    cs!(c"Offset argument must be a constant integer"),
                                );
                                current_block = 4830776507462815627;
                            } else {
                                yyval.Node = New_Offset(
                                    lParse,
                                    yyvs[yyvsp - 3].lng as c_int,
                                    yyvs[yyvsp - 1].Node,
                                );
                                if yyval.Node < 0 {
                                    current_block = 4830776507462815627;
                                } else {
                                    current_block = 17353983478346836848;
                                }
                            }
                        }
                        35 => {
                            /* expr: ROWREF  */
                            yyval.Node = New_Func(
                                lParse,
                                fits_parser_yytokentype::LONG as c_int,
                                ROW_FCT,
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
                            /* expr: NULLREF  */
                            yyval.Node = New_Func(
                                lParse,
                                fits_parser_yytokentype::LONG as c_int,
                                NULL_FCT,
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
                            /* expr: expr '%' expr  */
                            if (lParse.Nodes[yyvs[yyvsp - 2].Node as usize]).ntype
                                > (lParse.Nodes[yyvs[yyvsp].Node as usize]).ntype
                            {
                                yyvs[yyvsp].Node = New_Unary(
                                    lParse,
                                    ((lParse.Nodes)[yyvs[yyvsp - 2].Node as usize]).ntype,
                                    0,
                                    yyvs[yyvsp].Node,
                                );
                            } else if ((lParse.Nodes)[yyvs[yyvsp - 2].Node as usize]).ntype
                                < (lParse.Nodes[yyvs[yyvsp].Node as usize]).ntype
                            {
                                yyvs[yyvsp - 2].Node = New_Unary(
                                    lParse,
                                    (lParse.Nodes[yyvs[yyvsp].Node as usize]).ntype,
                                    0,
                                    yyvs[yyvsp - 2].Node,
                                );
                            }
                            yyval.Node = New_BinOp(
                                lParse,
                                (lParse.Nodes[yyvs[yyvsp - 2].Node as usize]).ntype,
                                yyvs[yyvsp - 2].Node,
                                '%' as i32,
                                yyvs[yyvsp].Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        38 => {
                            /* expr: expr '+' expr  */
                            if (lParse.Nodes[yyvs[yyvsp - 2].Node as usize]).ntype
                                > (lParse.Nodes[yyvs[yyvsp].Node as usize]).ntype
                            {
                                yyvs[yyvsp].Node = New_Unary(
                                    lParse,
                                    ((lParse.Nodes)[yyvs[yyvsp - 2].Node as usize]).ntype,
                                    0,
                                    yyvs[yyvsp].Node,
                                );
                            } else if ((lParse.Nodes)[yyvs[yyvsp - 2].Node as usize]).ntype
                                < (lParse.Nodes[yyvs[yyvsp].Node as usize]).ntype
                            {
                                yyvs[yyvsp - 2].Node = New_Unary(
                                    lParse,
                                    (lParse.Nodes[yyvs[yyvsp].Node as usize]).ntype,
                                    0,
                                    yyvs[yyvsp - 2].Node,
                                );
                            }
                            yyval.Node = New_BinOp(
                                lParse,
                                (lParse.Nodes[yyvs[yyvsp - 2].Node as usize]).ntype,
                                yyvs[yyvsp - 2].Node,
                                '+' as i32,
                                yyvs[yyvsp].Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        39 => {
                            /* expr: expr '-' expr  */
                            if (lParse.Nodes[yyvs[yyvsp - 2].Node as usize]).ntype
                                > (lParse.Nodes[yyvs[yyvsp].Node as usize]).ntype
                            {
                                yyvs[yyvsp].Node = New_Unary(
                                    lParse,
                                    ((lParse.Nodes)[yyvs[yyvsp - 2].Node as usize]).ntype,
                                    0,
                                    yyvs[yyvsp].Node,
                                );
                            } else if ((lParse.Nodes)[yyvs[yyvsp - 2].Node as usize]).ntype
                                < (lParse.Nodes[yyvs[yyvsp].Node as usize]).ntype
                            {
                                yyvs[yyvsp - 2].Node = New_Unary(
                                    lParse,
                                    (lParse.Nodes[yyvs[yyvsp].Node as usize]).ntype,
                                    0,
                                    yyvs[yyvsp - 2].Node,
                                );
                            }
                            yyval.Node = New_BinOp(
                                lParse,
                                (lParse.Nodes[yyvs[yyvsp - 2].Node as usize]).ntype,
                                yyvs[yyvsp - 2].Node,
                                '-' as i32,
                                yyvs[yyvsp].Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        40 => {
                            /* expr: expr '*' expr  */
                            if (lParse.Nodes[yyvs[yyvsp - 2].Node as usize]).ntype
                                > (lParse.Nodes[yyvs[yyvsp].Node as usize]).ntype
                            {
                                yyvs[yyvsp].Node = New_Unary(
                                    lParse,
                                    ((lParse.Nodes)[yyvs[yyvsp - 2].Node as usize]).ntype,
                                    0,
                                    yyvs[yyvsp].Node,
                                );
                            } else if ((lParse.Nodes)[yyvs[yyvsp - 2].Node as usize]).ntype
                                < (lParse.Nodes[yyvs[yyvsp].Node as usize]).ntype
                            {
                                yyvs[yyvsp - 2].Node = New_Unary(
                                    lParse,
                                    (lParse.Nodes[yyvs[yyvsp].Node as usize]).ntype,
                                    0,
                                    yyvs[yyvsp - 2].Node,
                                );
                            }
                            yyval.Node = New_BinOp(
                                lParse,
                                (lParse.Nodes[yyvs[yyvsp - 2].Node as usize]).ntype,
                                yyvs[yyvsp - 2].Node,
                                '*' as i32,
                                yyvs[yyvsp].Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        41 => {
                            /* expr: expr '/' expr  */
                            if (lParse.Nodes[yyvs[yyvsp - 2].Node as usize]).ntype
                                > (lParse.Nodes[yyvs[yyvsp].Node as usize]).ntype
                            {
                                yyvs[yyvsp].Node = New_Unary(
                                    lParse,
                                    ((lParse.Nodes)[yyvs[yyvsp - 2].Node as usize]).ntype,
                                    0,
                                    yyvs[yyvsp].Node,
                                );
                            } else if ((lParse.Nodes)[yyvs[yyvsp - 2].Node as usize]).ntype
                                < (lParse.Nodes[yyvs[yyvsp].Node as usize]).ntype
                            {
                                yyvs[yyvsp - 2].Node = New_Unary(
                                    lParse,
                                    (lParse.Nodes[yyvs[yyvsp].Node as usize]).ntype,
                                    0,
                                    yyvs[yyvsp - 2].Node,
                                );
                            }
                            yyval.Node = New_BinOp(
                                lParse,
                                (lParse.Nodes[yyvs[yyvsp - 2].Node as usize]).ntype,
                                yyvs[yyvsp - 2].Node,
                                '/' as i32,
                                yyvs[yyvsp].Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        42 => {
                            /* expr: expr '&' expr  */
                            if (lParse.Nodes[yyvs[yyvsp - 2].Node as usize]).ntype
                                != fits_parser_yytokentype::LONG as c_int
                                || (lParse.Nodes[yyvs[yyvsp].Node as usize]).ntype
                                    != fits_parser_yytokentype::LONG as c_int
                            {
                                fits_parser_yyerror(lParse, cs!(c"Bitwise operations with incompatible types; only (bit OP bit) and (int OP int) are allowed"));
                                current_block = 4830776507462815627;
                            } else {
                                yyval.Node = New_BinOp(
                                    lParse,
                                    ((lParse.Nodes)[yyvs[yyvsp - 2].Node as usize]).ntype,
                                    yyvs[yyvsp - 2].Node,
                                    '&' as i32,
                                    yyvs[yyvsp].Node,
                                );
                                current_block = 17353983478346836848;
                            }
                        }
                        43 => {
                            /* expr: expr '|' expr  */
                            if (lParse.Nodes[yyvs[yyvsp - 2].Node as usize]).ntype
                                != fits_parser_yytokentype::LONG as c_int
                                || (lParse.Nodes[yyvs[yyvsp].Node as usize]).ntype
                                    != fits_parser_yytokentype::LONG as c_int
                            {
                                fits_parser_yyerror(lParse, cs!(c"Bitwise operations with incompatible types; only (bit OP bit) and (int OP int) are allowed"));
                                current_block = 4830776507462815627;
                            } else {
                                yyval.Node = New_BinOp(
                                    lParse,
                                    ((lParse.Nodes)[yyvs[yyvsp - 2].Node as usize]).ntype,
                                    yyvs[yyvsp - 2].Node,
                                    '|' as i32,
                                    yyvs[yyvsp].Node,
                                );
                                current_block = 17353983478346836848;
                            }
                        }
                        44 => {
                            /* expr: expr XOR expr  */
                            if (lParse.Nodes[yyvs[yyvsp - 2].Node as usize]).ntype
                                != fits_parser_yytokentype::LONG as c_int
                                || (lParse.Nodes[yyvs[yyvsp].Node as usize]).ntype
                                    != fits_parser_yytokentype::LONG as c_int
                            {
                                fits_parser_yyerror(lParse, cs!(c"Bitwise operations with incompatible types; only (bit OP bit) and (int OP int) are allowed"));
                                current_block = 4830776507462815627;
                            } else {
                                yyval.Node = New_BinOp(
                                    lParse,
                                    ((lParse.Nodes)[yyvs[yyvsp - 2].Node as usize]).ntype,
                                    yyvs[yyvsp - 2].Node,
                                    '^' as i32,
                                    yyvs[yyvsp].Node,
                                );
                                current_block = 17353983478346836848;
                            }
                        }
                        45 => {
                            /* expr: expr POWER expr  */
                            if (lParse.Nodes[yyvs[yyvsp - 2].Node as usize]).ntype
                                > (lParse.Nodes[yyvs[yyvsp].Node as usize]).ntype
                            {
                                yyvs[yyvsp].Node = New_Unary(
                                    lParse,
                                    ((lParse.Nodes)[yyvs[yyvsp - 2].Node as usize]).ntype,
                                    0,
                                    yyvs[yyvsp].Node,
                                );
                            } else if ((lParse.Nodes)[yyvs[yyvsp - 2].Node as usize]).ntype
                                < (lParse.Nodes[yyvs[yyvsp].Node as usize]).ntype
                            {
                                yyvs[yyvsp - 2].Node = New_Unary(
                                    lParse,
                                    (lParse.Nodes[yyvs[yyvsp].Node as usize]).ntype,
                                    0,
                                    yyvs[yyvsp - 2].Node,
                                );
                            }
                            yyval.Node = New_BinOp(
                                lParse,
                                (lParse.Nodes[yyvs[yyvsp - 2].Node as usize]).ntype,
                                yyvs[yyvsp - 2].Node,
                                fits_parser_yytokentype::POWER as c_int,
                                yyvs[yyvsp].Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        46 => {
                            /* expr: '+' expr  */
                            yyval.Node = yyvs[yyvsp].Node;
                            current_block = 17353983478346836848;
                        }
                        47 => {
                            /* expr: '-' expr  */
                            yyval.Node = New_Unary(
                                lParse,
                                (lParse.Nodes[yyvs[yyvsp].Node as usize]).ntype,
                                fits_parser_yytokentype::UMINUS as c_int,
                                yyvs[yyvsp].Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        48 => {
                            /* expr: '(' expr ')'  */
                            yyval.Node = yyvs[yyvsp - 1].Node;
                            current_block = 17353983478346836848;
                        }
                        49 => {
                            /* expr: expr '*' bexpr  */
                            yyvs[yyvsp].Node = New_Unary(
                                lParse,
                                (lParse.Nodes[yyvs[yyvsp - 2].Node as usize]).ntype,
                                0,
                                yyvs[yyvsp].Node,
                            );
                            yyval.Node = New_BinOp(
                                lParse,
                                (lParse.Nodes[yyvs[yyvsp - 2].Node as usize]).ntype,
                                yyvs[yyvsp - 2].Node,
                                '*' as i32,
                                yyvs[yyvsp].Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        50 => {
                            /* expr: bexpr '*' expr  */
                            yyvs[yyvsp - 2].Node = New_Unary(
                                lParse,
                                (lParse.Nodes[yyvs[yyvsp].Node as usize]).ntype,
                                0,
                                yyvs[yyvsp - 2].Node,
                            );
                            yyval.Node = New_BinOp(
                                lParse,
                                (lParse.Nodes[yyvs[yyvsp].Node as usize]).ntype,
                                yyvs[yyvsp - 2].Node,
                                '*' as i32,
                                yyvs[yyvsp].Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        51 => {
                            /* expr: bexpr '?' expr ':' expr  */
                            if (lParse.Nodes[yyvs[yyvsp - 2].Node as usize]).ntype
                                > (lParse.Nodes[yyvs[yyvsp].Node as usize]).ntype
                            {
                                yyvs[yyvsp].Node = New_Unary(
                                    lParse,
                                    ((lParse.Nodes)[yyvs[yyvsp - 2].Node as usize]).ntype,
                                    0,
                                    yyvs[yyvsp].Node,
                                );
                            } else if ((lParse.Nodes)[yyvs[yyvsp - 2].Node as usize]).ntype
                                < (lParse.Nodes[yyvs[yyvsp].Node as usize]).ntype
                            {
                                yyvs[yyvsp - 2].Node = New_Unary(
                                    lParse,
                                    (lParse.Nodes[yyvs[yyvsp].Node as usize]).ntype,
                                    0,
                                    yyvs[yyvsp - 2].Node,
                                );
                            }
                            if Test_Dims(lParse, yyvs[yyvsp - 2].Node, yyvs[yyvsp].Node) == 0 {
                                fits_parser_yyerror(
                                    lParse,
                                    cs!(c"Incompatible dimensions in '?:' arguments"),
                                );
                                current_block = 4830776507462815627;
                            } else {
                                yyval.Node = New_Func(
                                    lParse,
                                    0,
                                    IFTHENELSE_FCT,
                                    3 as c_int,
                                    yyvs[yyvsp - 2].Node,
                                    yyvs[yyvsp].Node,
                                    yyvs[yyvsp - 4].Node,
                                    0,
                                    0,
                                    0,
                                    0,
                                );
                                if yyval.Node < 0 {
                                    current_block = 4830776507462815627;
                                } else {
                                    if ((lParse.Nodes)[yyvs[yyvsp - 2].Node as usize]).value.nelem
                                        < (lParse.Nodes[yyvs[yyvsp].Node as usize]).value.nelem
                                    {
                                        Copy_Dims(lParse, yyval.Node, yyvs[yyvsp].Node);
                                    }
                                    ((lParse.Nodes)[yyvs[yyvsp - 4].Node as usize]).ntype =
                                        ((lParse.Nodes)[yyvs[yyvsp - 2].Node as usize]).ntype;
                                    if Test_Dims(lParse, yyvs[yyvsp - 4].Node, yyval.Node) == 0 {
                                        fits_parser_yyerror(
                                            lParse,
                                            cs!(c"Incompatible dimensions in '?:' condition"),
                                        );
                                        current_block = 4830776507462815627;
                                    } else {
                                        ((lParse.Nodes)[yyvs[yyvsp - 4].Node as usize]).ntype =
                                            fits_parser_yytokentype::BOOLEAN as c_int;
                                        if (((lParse.Nodes)[yyval.Node as usize]).value).nelem
                                            < ((lParse.Nodes)[yyvs[yyvsp - 4].Node as usize])
                                                .value
                                                .nelem
                                        {
                                            Copy_Dims(lParse, yyval.Node, yyvs[yyvsp - 4].Node);
                                        }
                                        current_block = 17353983478346836848;
                                    }
                                }
                            }
                        }
                        52 => {
                            /* expr: bexpr '?' bexpr ':' expr  */
                            if (lParse.Nodes[yyvs[yyvsp - 2].Node as usize]).ntype
                                > (lParse.Nodes[yyvs[yyvsp].Node as usize]).ntype
                            {
                                yyvs[yyvsp].Node = New_Unary(
                                    lParse,
                                    ((lParse.Nodes)[yyvs[yyvsp - 2].Node as usize]).ntype,
                                    0,
                                    yyvs[yyvsp].Node,
                                );
                            } else if ((lParse.Nodes)[yyvs[yyvsp - 2].Node as usize]).ntype
                                < (lParse.Nodes[yyvs[yyvsp].Node as usize]).ntype
                            {
                                yyvs[yyvsp - 2].Node = New_Unary(
                                    lParse,
                                    (lParse.Nodes[yyvs[yyvsp].Node as usize]).ntype,
                                    0,
                                    yyvs[yyvsp - 2].Node,
                                );
                            }
                            if Test_Dims(lParse, yyvs[yyvsp - 2].Node, yyvs[yyvsp].Node) == 0 {
                                fits_parser_yyerror(
                                    lParse,
                                    cs!(c"Incompatible dimensions in '?:' arguments"),
                                );
                                current_block = 4830776507462815627;
                            } else {
                                yyval.Node = New_Func(
                                    lParse,
                                    0,
                                    IFTHENELSE_FCT,
                                    3 as c_int,
                                    yyvs[yyvsp - 2].Node,
                                    yyvs[yyvsp].Node,
                                    yyvs[yyvsp - 4].Node,
                                    0,
                                    0,
                                    0,
                                    0,
                                );
                                if yyval.Node < 0 {
                                    current_block = 4830776507462815627;
                                } else {
                                    if ((lParse.Nodes)[yyvs[yyvsp - 2].Node as usize]).value.nelem
                                        < (lParse.Nodes[yyvs[yyvsp].Node as usize]).value.nelem
                                    {
                                        Copy_Dims(lParse, yyval.Node, yyvs[yyvsp].Node);
                                    }
                                    ((lParse.Nodes)[yyvs[yyvsp - 4].Node as usize]).ntype =
                                        ((lParse.Nodes)[yyvs[yyvsp - 2].Node as usize]).ntype;
                                    if Test_Dims(lParse, yyvs[yyvsp - 4].Node, yyval.Node) == 0 {
                                        fits_parser_yyerror(
                                            lParse,
                                            cs!(c"Incompatible dimensions in '?:' condition"),
                                        );
                                        current_block = 4830776507462815627;
                                    } else {
                                        ((lParse.Nodes)[yyvs[yyvsp - 4].Node as usize]).ntype =
                                            fits_parser_yytokentype::BOOLEAN as c_int;
                                        if (((lParse.Nodes)[yyval.Node as usize]).value).nelem
                                            < ((lParse.Nodes)[yyvs[yyvsp - 4].Node as usize])
                                                .value
                                                .nelem
                                        {
                                            Copy_Dims(lParse, yyval.Node, yyvs[yyvsp - 4].Node);
                                        }
                                        current_block = 17353983478346836848;
                                    }
                                }
                            }
                        }
                        53 => {
                            /* expr: bexpr '?' expr ':' bexpr  */
                            if (lParse.Nodes[yyvs[yyvsp - 2].Node as usize]).ntype
                                > (lParse.Nodes[yyvs[yyvsp].Node as usize]).ntype
                            {
                                yyvs[yyvsp].Node = New_Unary(
                                    lParse,
                                    ((lParse.Nodes)[yyvs[yyvsp - 2].Node as usize]).ntype,
                                    0,
                                    yyvs[yyvsp].Node,
                                );
                            } else if ((lParse.Nodes)[yyvs[yyvsp - 2].Node as usize]).ntype
                                < (lParse.Nodes[yyvs[yyvsp].Node as usize]).ntype
                            {
                                yyvs[yyvsp - 2].Node = New_Unary(
                                    lParse,
                                    (lParse.Nodes[yyvs[yyvsp].Node as usize]).ntype,
                                    0,
                                    yyvs[yyvsp - 2].Node,
                                );
                            }
                            if Test_Dims(lParse, yyvs[yyvsp - 2].Node, yyvs[yyvsp].Node) == 0 {
                                fits_parser_yyerror(
                                    lParse,
                                    cs!(c"Incompatible dimensions in '?:' arguments"),
                                );
                                current_block = 4830776507462815627;
                            } else {
                                yyval.Node = New_Func(
                                    lParse,
                                    0,
                                    IFTHENELSE_FCT,
                                    3 as c_int,
                                    yyvs[yyvsp - 2].Node,
                                    yyvs[yyvsp].Node,
                                    yyvs[yyvsp - 4].Node,
                                    0,
                                    0,
                                    0,
                                    0,
                                );
                                if yyval.Node < 0 {
                                    current_block = 4830776507462815627;
                                } else {
                                    if ((lParse.Nodes)[yyvs[yyvsp - 2].Node as usize]).value.nelem
                                        < (lParse.Nodes[yyvs[yyvsp].Node as usize]).value.nelem
                                    {
                                        Copy_Dims(lParse, yyval.Node, yyvs[yyvsp].Node);
                                    }
                                    ((lParse.Nodes)[yyvs[yyvsp - 4].Node as usize]).ntype =
                                        ((lParse.Nodes)[yyvs[yyvsp - 2].Node as usize]).ntype;
                                    if Test_Dims(lParse, yyvs[yyvsp - 4].Node, yyval.Node) == 0 {
                                        fits_parser_yyerror(
                                            lParse,
                                            cs!(c"Incompatible dimensions in '?:' condition"),
                                        );
                                        current_block = 4830776507462815627;
                                    } else {
                                        ((lParse.Nodes)[yyvs[yyvsp - 4].Node as usize]).ntype =
                                            fits_parser_yytokentype::BOOLEAN as c_int;
                                        if (((lParse.Nodes)[yyval.Node as usize]).value).nelem
                                            < ((lParse.Nodes)[yyvs[yyvsp - 4].Node as usize])
                                                .value
                                                .nelem
                                        {
                                            Copy_Dims(lParse, yyval.Node, yyvs[yyvsp - 4].Node);
                                        }
                                        current_block = 17353983478346836848;
                                    }
                                }
                            }
                        }
                        54 => {
                            /* expr: FUNCTION ')'  */
                            if (if c_int::from(yyvs[yyvsp - 1].astr[0])
                                < c_int::from((cs!(c"RANDOM("))[0])
                            {
                                -(1)
                            } else if c_int::from(yyvs[yyvsp - 1].astr[0])
                                > c_int::from((cs!(c"RANDOM("))[0])
                            {
                                1
                            } else {
                                strcmp((yyvs[yyvsp - 1].astr).as_mut_ptr(), c"RANDOM(".as_ptr())
                            }) == 0
                            {
                                yyval.Node = New_Func(
                                    lParse,
                                    fits_parser_yytokentype::DOUBLE as c_int,
                                    RND_FCT,
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
                            } else if (if c_int::from(yyvs[yyvsp - 1].astr[0])
                                < c_int::from((cs!(c"RANDOMN("))[0])
                            {
                                -(1)
                            } else if c_int::from(yyvs[yyvsp - 1].astr[0])
                                > c_int::from((cs!(c"RANDOMN("))[0])
                            {
                                1
                            } else {
                                strcmp((yyvs[yyvsp - 1].astr).as_mut_ptr(), c"RANDOMN(".as_ptr())
                            }) == 0
                            {
                                yyval.Node = New_Func(
                                    lParse,
                                    fits_parser_yytokentype::DOUBLE as c_int,
                                    GASRND_FCT,
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
                                fits_parser_yyerror(lParse, cs!(c"Function() not supported"));
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
                            /* expr: FUNCTION bexpr ')'  */
                            if (if c_int::from(yyvs[yyvsp - 2].astr[0])
                                < c_int::from((cs!(c"SUM("))[0])
                            {
                                -(1)
                            } else if c_int::from(yyvs[yyvsp - 2].astr[0])
                                > c_int::from((cs!(c"SUM("))[0])
                            {
                                1
                            } else {
                                strcmp((yyvs[yyvsp - 2].astr).as_mut_ptr(), c"SUM(".as_ptr())
                            }) == 0
                            {
                                yyval.Node = New_Func(
                                    lParse,
                                    fits_parser_yytokentype::LONG as c_int,
                                    SUM_FCT,
                                    1,
                                    yyvs[yyvsp - 1].Node,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                );
                                current_block = 10848699504537784535;
                            } else if (if c_int::from(yyvs[yyvsp - 2].astr[0])
                                < c_int::from((cs!(c"NELEM("))[0])
                            {
                                -(1)
                            } else if c_int::from(yyvs[yyvsp - 2].astr[0])
                                > c_int::from((cs!(c"NELEM("))[0])
                            {
                                1
                            } else {
                                strcmp((yyvs[yyvsp - 2].astr).as_mut_ptr(), c"NELEM(".as_ptr())
                            }) == 0
                            {
                                yyval.Node = New_Const(
                                    lParse,
                                    fits_parser_yytokentype::LONG as c_int,
                                    (&((lParse.Nodes)[yyvs[yyvsp - 1].Node as usize]).value.nelem
                                        as *const c_long)
                                        .cast::<c_void>(),
                                    ::core::mem::size_of::<c_long>() as c_ulong as c_long,
                                );
                                current_block = 10848699504537784535;
                            } else if (if c_int::from(yyvs[yyvsp - 2].astr[0])
                                < c_int::from((cs!(c"ACCUM("))[0])
                            {
                                -(1)
                            } else if c_int::from(yyvs[yyvsp - 2].astr[0])
                                > c_int::from((cs!(c"ACCUM("))[0])
                            {
                                1
                            } else {
                                strcmp((yyvs[yyvsp - 2].astr).as_mut_ptr(), c"ACCUM(".as_ptr())
                            }) == 0
                            {
                                let mut zero: c_long = 0;
                                let new_node = New_Const(
                                    lParse,
                                    fits_parser_yytokentype::LONG as c_int,
                                    (&mut zero as *mut c_long).cast::<c_void>(),
                                    ::core::mem::size_of::<c_long>() as c_ulong as c_long,
                                );

                                yyval.Node = New_BinOp(
                                    lParse,
                                    fits_parser_yytokentype::LONG as c_int,
                                    yyvs[yyvsp - 1].Node,
                                    fits_parser_yytokentype::ACCUM as c_int,
                                    new_node,
                                );
                                current_block = 10848699504537784535;
                            } else {
                                fits_parser_yyerror(lParse, cs!(c"Function(bool) not supported"));
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
                            /* expr: FUNCTION bexpr ',' expr ')'  */
                            if (if c_int::from(yyvs[yyvsp - 4].astr[0])
                                < c_int::from((cs!(c"AXISELEM"))[0])
                            {
                                -(1)
                            } else if c_int::from(yyvs[yyvsp - 4].astr[0])
                                > c_int::from((cs!(c"AXISELEM"))[0])
                            {
                                1
                            } else {
                                strcmp((yyvs[yyvsp - 4].astr).as_mut_ptr(), c"AXISELEM(".as_ptr())
                            }) == 0
                            {
                                if ((lParse.Nodes)[yyvs[yyvsp - 1].Node as usize]).operation
                                    != CONST_OP
                                    || ((lParse.Nodes)[yyvs[yyvsp - 1].Node as usize]).value.nelem
                                        != 1
                                {
                                    fits_parser_yyerror(
                                        lParse,
                                        cs!(c"AXISELEM second argument must be a scalar constant"),
                                    );
                                    current_block = 4830776507462815627;
                                } else if ((lParse.Nodes)[yyvs[yyvsp - 3].Node as usize]).operation
                                    == CONST_OP
                                {
                                    let mut one: c_long = 1;
                                    yyval.Node = New_Const(
                                        lParse,
                                        fits_parser_yytokentype::LONG as c_int,
                                        (&mut one as *mut c_long).cast::<c_void>(),
                                        ::core::mem::size_of::<c_long>() as c_ulong as c_long,
                                    );
                                    current_block = 13755523488868872559;
                                } else {
                                    if ((lParse.Nodes)[yyvs[yyvsp - 1].Node as usize]).ntype
                                        != fits_parser_yytokentype::LONG as c_int
                                    {
                                        yyvs[yyvsp - 1].Node = New_Unary(
                                            lParse,
                                            fits_parser_yytokentype::LONG as c_int,
                                            0,
                                            yyvs[yyvsp - 1].Node,
                                        );
                                    }
                                    yyval.Node = New_Func(
                                        lParse,
                                        0,
                                        AXISELEM_FCT,
                                        2,
                                        yyvs[yyvsp - 3].Node,
                                        yyvs[yyvsp - 1].Node,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                    );
                                    if yyval.Node < 0 {
                                        current_block = 4830776507462815627;
                                    } else {
                                        ((lParse.Nodes)[yyval.Node as usize]).ntype =
                                            fits_parser_yytokentype::LONG as c_int;
                                        current_block = 13755523488868872559;
                                    }
                                }
                            } else if (if c_int::from(yyvs[yyvsp - 4].astr[0])
                                < c_int::from((cs!(c"NAXES("))[0])
                            {
                                -(1)
                            } else if c_int::from(yyvs[yyvsp - 4].astr[0])
                                > c_int::from((cs!(c"NAXES("))[0])
                            {
                                1
                            } else {
                                strcmp((yyvs[yyvsp - 4].astr).as_mut_ptr(), c"NAXES(".as_ptr())
                            }) == 0
                            {
                                if ((lParse.Nodes)[yyvs[yyvsp - 1].Node as usize]).operation
                                    != CONST_OP
                                    || ((lParse.Nodes)[yyvs[yyvsp - 1].Node as usize]).value.nelem
                                        != 1
                                {
                                    fits_parser_yyerror(
                                        lParse,
                                        cs!(c"NAXES second argument must be a scalar constant"),
                                    );
                                    current_block = 4830776507462815627;
                                } else if ((lParse.Nodes)[yyvs[yyvsp - 3].Node as usize]).operation
                                    == CONST_OP
                                {
                                    let mut one_0: c_long = 1;
                                    yyval.Node = New_Const(
                                        lParse,
                                        fits_parser_yytokentype::LONG as c_int,
                                        (&mut one_0 as *mut c_long).cast::<c_void>(),
                                        ::core::mem::size_of::<c_long>() as c_ulong as c_long,
                                    );
                                    current_block = 13755523488868872559;
                                } else {
                                    let mut iaxis: c_long = 0;
                                    let mut naxis: c_int = 0;
                                    if ((lParse.Nodes)[yyvs[yyvsp - 1].Node as usize]).ntype
                                        != fits_parser_yytokentype::LONG as c_int
                                    {
                                        yyvs[yyvsp - 1].Node = New_Unary(
                                            lParse,
                                            fits_parser_yytokentype::LONG as c_int,
                                            0,
                                            yyvs[yyvsp - 1].Node,
                                        );
                                    }
                                    iaxis = ((lParse.Nodes)[yyvs[yyvsp - 1].Node as usize])
                                        .value
                                        .data
                                        .lng;
                                    naxis =
                                        ((lParse.Nodes)[yyvs[yyvsp - 3].Node as usize]).value.naxis;
                                    if iaxis == 0 {
                                        iaxis = c_long::from(naxis);
                                    } else if iaxis <= c_long::from(naxis) {
                                        iaxis = ((lParse.Nodes)[yyvs[yyvsp - 3].Node as usize])
                                            .value
                                            .naxes
                                            [(iaxis - 1) as usize];
                                    } else {
                                        iaxis = 1;
                                    }
                                    yyval.Node = New_Const(
                                        lParse,
                                        fits_parser_yytokentype::LONG as c_int,
                                        (&mut iaxis as *mut c_long).cast::<c_void>(),
                                        ::core::mem::size_of::<c_long>() as c_ulong as c_long,
                                    );
                                    if yyval.Node < 0 {
                                        current_block = 4830776507462815627;
                                    } else {
                                        current_block = 13755523488868872559;
                                    }
                                }
                            } else if (if c_int::from(yyvs[yyvsp - 4].astr[0])
                                < c_int::from((cs!(c"ARRAY("))[0])
                            {
                                -(1)
                            } else if c_int::from(yyvs[yyvsp - 4].astr[0])
                                > c_int::from((cs!(c"ARRAY("))[0])
                            {
                                1
                            } else {
                                strcmp((yyvs[yyvsp - 4].astr).as_mut_ptr(), c"ARRAY(".as_ptr())
                            }) == 0
                            {
                                yyval.Node =
                                    New_Array(lParse, yyvs[yyvsp - 3].Node, yyvs[yyvsp - 1].Node);
                                if yyval.Node < 0 {
                                    current_block = 4830776507462815627;
                                } else {
                                    current_block = 13755523488868872559;
                                }
                            } else {
                                fits_parser_yyerror(
                                    lParse,
                                    cs!(c"Function(bool,expr) not supported"),
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
                            /* expr: FUNCTION sexpr ')'  */
                            if (if c_int::from(yyvs[yyvsp - 2].astr[0])
                                < c_int::from((cs!(c"NELEM("))[0])
                            {
                                -(1)
                            } else if c_int::from(yyvs[yyvsp - 2].astr[0])
                                > c_int::from((cs!(c"NELEM("))[0])
                            {
                                1
                            } else {
                                strcmp((yyvs[yyvsp - 2].astr).as_mut_ptr(), c"NELEM(".as_ptr())
                            }) == 0
                            {
                                yyval.Node = New_Const(
                                    lParse,
                                    fits_parser_yytokentype::LONG as c_int,
                                    (&((lParse.Nodes)[yyvs[yyvsp - 1].Node as usize]).value.nelem
                                        as *const c_long)
                                        .cast::<c_void>(),
                                    ::core::mem::size_of::<c_long>() as c_ulong as c_long,
                                );
                                current_block = 15752106442776732052;
                            } else if (if c_int::from(yyvs[yyvsp - 2].astr[0])
                                < c_int::from((cs!(c"NVALID("))[0])
                            {
                                -(1)
                            } else if c_int::from(yyvs[yyvsp - 2].astr[0])
                                > c_int::from((cs!(c"NVALID("))[0])
                            {
                                1
                            } else {
                                strcmp((yyvs[yyvsp - 2].astr).as_mut_ptr(), c"NVALID(".as_ptr())
                            }) == 0
                            {
                                yyval.Node = New_Func(
                                    lParse,
                                    fits_parser_yytokentype::LONG as c_int,
                                    NONNULL_FCT,
                                    1,
                                    yyvs[yyvsp - 1].Node,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                );
                                current_block = 15752106442776732052;
                            } else {
                                fits_parser_yyerror(lParse, cs!(c"Function(str) not supported"));
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
                            /* expr: FUNCTION bits ')'  */
                            if (if c_int::from(yyvs[yyvsp - 2].astr[0])
                                < c_int::from((cs!(c"NELEM("))[0])
                            {
                                -(1)
                            } else if c_int::from(yyvs[yyvsp - 2].astr[0])
                                > c_int::from((cs!(c"NELEM("))[0])
                            {
                                1
                            } else {
                                strcmp((yyvs[yyvsp - 2].astr).as_mut_ptr(), c"NELEM(".as_ptr())
                            }) == 0
                            {
                                yyval.Node = New_Const(
                                    lParse,
                                    fits_parser_yytokentype::LONG as c_int,
                                    (&((lParse.Nodes)[yyvs[yyvsp - 1].Node as usize]).value.nelem
                                        as *const c_long)
                                        .cast::<c_void>(),
                                    ::core::mem::size_of::<c_long>() as c_ulong as c_long,
                                );
                                current_block = 494012601817399562;
                            } else if (if c_int::from(yyvs[yyvsp - 2].astr[0])
                                < c_int::from((cs!(c"NVALID("))[0])
                            {
                                -(1)
                            } else if c_int::from(yyvs[yyvsp - 2].astr[0])
                                > c_int::from((cs!(c"NVALID("))[0])
                            {
                                1
                            } else {
                                strcmp((yyvs[yyvsp - 2].astr).as_mut_ptr(), c"NVALID(".as_ptr())
                            }) == 0
                            {
                                yyval.Node = New_Const(
                                    lParse,
                                    fits_parser_yytokentype::LONG as c_int,
                                    (&((lParse.Nodes)[yyvs[yyvsp - 1].Node as usize]).value.nelem
                                        as *const c_long)
                                        .cast::<c_void>(),
                                    ::core::mem::size_of::<c_long>() as c_ulong as c_long,
                                );
                                current_block = 494012601817399562;
                            } else if (if c_int::from(yyvs[yyvsp - 2].astr[0])
                                < c_int::from((cs!(c"SUM("))[0])
                            {
                                -(1)
                            } else if c_int::from(yyvs[yyvsp - 2].astr[0])
                                > c_int::from((cs!(c"SUM("))[0])
                            {
                                1
                            } else {
                                strcmp((yyvs[yyvsp - 2].astr).as_mut_ptr(), c"SUM(".as_ptr())
                            }) == 0
                            {
                                yyval.Node = New_Func(
                                    lParse,
                                    fits_parser_yytokentype::LONG as c_int,
                                    SUM_FCT,
                                    1,
                                    yyvs[yyvsp - 1].Node,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                );
                                current_block = 494012601817399562;
                            } else if (if c_int::from(yyvs[yyvsp - 2].astr[0])
                                < c_int::from((cs!(c"MIN("))[0])
                            {
                                -(1)
                            } else if c_int::from(yyvs[yyvsp - 2].astr[0])
                                > c_int::from((cs!(c"MIN("))[0])
                            {
                                1
                            } else {
                                strcmp((yyvs[yyvsp - 2].astr).as_mut_ptr(), c"MIN(".as_ptr())
                            }) == 0
                            {
                                yyval.Node = New_Func(
                                    lParse,
                                    ((lParse.Nodes)[yyvs[yyvsp - 1].Node as usize]).ntype,
                                    MIN1_FCT,
                                    1,
                                    yyvs[yyvsp - 1].Node,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                );
                                (((lParse.Nodes)[yyval.Node as usize]).value).nelem = 1;
                                current_block = 494012601817399562;
                            } else if (if c_int::from(yyvs[yyvsp - 2].astr[0])
                                < c_int::from((cs!(c"ACCUM("))[0])
                            {
                                -(1)
                            } else if c_int::from(yyvs[yyvsp - 2].astr[0])
                                > c_int::from((cs!(c"ACCUM("))[0])
                            {
                                1
                            } else {
                                strcmp((yyvs[yyvsp - 2].astr).as_mut_ptr(), c"ACCUM(".as_ptr())
                            }) == 0
                            {
                                let mut zero_0: c_long = 0;
                                let new_node = New_Const(
                                    lParse,
                                    fits_parser_yytokentype::LONG as c_int,
                                    (&mut zero_0 as *mut c_long).cast::<c_void>(),
                                    ::core::mem::size_of::<c_long>() as c_ulong as c_long,
                                );

                                yyval.Node = New_BinOp(
                                    lParse,
                                    fits_parser_yytokentype::LONG as c_int,
                                    yyvs[yyvsp - 1].Node,
                                    fits_parser_yytokentype::ACCUM as c_int,
                                    new_node,
                                );
                                current_block = 494012601817399562;
                            } else if (if c_int::from(yyvs[yyvsp - 2].astr[0])
                                < c_int::from((cs!(c"MAX("))[0])
                            {
                                -(1)
                            } else if c_int::from(yyvs[yyvsp - 2].astr[0])
                                > c_int::from((cs!(c"MAX("))[0])
                            {
                                1
                            } else {
                                strcmp((yyvs[yyvsp - 2].astr).as_mut_ptr(), c"MAX(".as_ptr())
                            }) == 0
                            {
                                yyval.Node = New_Func(
                                    lParse,
                                    ((lParse.Nodes)[yyvs[yyvsp - 1].Node as usize]).ntype,
                                    MAX1_FCT,
                                    1,
                                    yyvs[yyvsp - 1].Node,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                );
                                (((lParse.Nodes)[yyval.Node as usize]).value).nelem = 1;
                                current_block = 494012601817399562;
                            } else {
                                fits_parser_yyerror(lParse, cs!(c"Function(bits) not supported"));
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
                            /* expr: FUNCTION expr ')'  */
                            if (if c_int::from(yyvs[yyvsp - 2].astr[0])
                                < c_int::from((cs!(c"SUM("))[0])
                            {
                                -(1)
                            } else if c_int::from(yyvs[yyvsp - 2].astr[0])
                                > c_int::from((cs!(c"SUM("))[0])
                            {
                                1
                            } else {
                                strcmp((yyvs[yyvsp - 2].astr).as_mut_ptr(), c"SUM(".as_ptr())
                            }) == 0
                            {
                                yyval.Node = New_Func(
                                    lParse,
                                    ((lParse.Nodes)[yyvs[yyvsp - 1].Node as usize]).ntype,
                                    SUM_FCT,
                                    1,
                                    yyvs[yyvsp - 1].Node,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                );
                                current_block = 7600445499126923600;
                            } else if (if c_int::from(yyvs[yyvsp - 2].astr[0])
                                < c_int::from((cs!(c"AVERAGE("))[0])
                            {
                                -(1)
                            } else if c_int::from(yyvs[yyvsp - 2].astr[0])
                                > c_int::from((cs!(c"AVERAGE("))[0])
                            {
                                1
                            } else {
                                strcmp((yyvs[yyvsp - 2].astr).as_mut_ptr(), c"AVERAGE(".as_ptr())
                            }) == 0
                            {
                                yyval.Node = New_Func(
                                    lParse,
                                    fits_parser_yytokentype::DOUBLE as c_int,
                                    AVERAGE_FCT,
                                    1,
                                    yyvs[yyvsp - 1].Node,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                );
                                current_block = 7600445499126923600;
                            } else if (if c_int::from(yyvs[yyvsp - 2].astr[0])
                                < c_int::from((cs!(c"STDDEV("))[0])
                            {
                                -(1)
                            } else if c_int::from(yyvs[yyvsp - 2].astr[0])
                                > c_int::from((cs!(c"STDDEV("))[0])
                            {
                                1
                            } else {
                                strcmp((yyvs[yyvsp - 2].astr).as_mut_ptr(), c"STDDEV(".as_ptr())
                            }) == 0
                            {
                                yyval.Node = New_Func(
                                    lParse,
                                    fits_parser_yytokentype::DOUBLE as c_int,
                                    STDDEV_FCT,
                                    1,
                                    yyvs[yyvsp - 1].Node,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                );
                                current_block = 7600445499126923600;
                            } else if (if c_int::from(yyvs[yyvsp - 2].astr[0])
                                < c_int::from((cs!(c"MEDIAN("))[0])
                            {
                                -(1)
                            } else if c_int::from(yyvs[yyvsp - 2].astr[0])
                                > c_int::from((cs!(c"MEDIAN("))[0])
                            {
                                1
                            } else {
                                strcmp((yyvs[yyvsp - 2].astr).as_mut_ptr(), c"MEDIAN(".as_ptr())
                            }) == 0
                            {
                                yyval.Node = New_Func(
                                    lParse,
                                    ((lParse.Nodes)[yyvs[yyvsp - 1].Node as usize]).ntype,
                                    MEDIAN_FCT,
                                    1,
                                    yyvs[yyvsp - 1].Node,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                );
                                current_block = 7600445499126923600;
                            } else if (if c_int::from(yyvs[yyvsp - 2].astr[0])
                                < c_int::from((cs!(c"NELEM("))[0])
                            {
                                -(1)
                            } else if c_int::from(yyvs[yyvsp - 2].astr[0])
                                > c_int::from((cs!(c"NELEM("))[0])
                            {
                                1
                            } else {
                                strcmp((yyvs[yyvsp - 2].astr).as_mut_ptr(), c"NELEM(".as_ptr())
                            }) == 0
                            {
                                yyval.Node = New_Const(
                                    lParse,
                                    fits_parser_yytokentype::LONG as c_int,
                                    (&((lParse.Nodes)[yyvs[yyvsp - 1].Node as usize]).value.nelem
                                        as *const c_long)
                                        .cast::<c_void>(),
                                    ::core::mem::size_of::<c_long>() as c_ulong as c_long,
                                );
                                current_block = 7600445499126923600;
                            } else if (if c_int::from(yyvs[yyvsp - 2].astr[0])
                                < c_int::from((cs!(c"NVALID("))[0])
                            {
                                -(1)
                            } else if c_int::from(yyvs[yyvsp - 2].astr[0])
                                > c_int::from((cs!(c"NVALID("))[0])
                            {
                                1
                            } else {
                                strcmp((yyvs[yyvsp - 2].astr).as_mut_ptr(), c"NVALID(".as_ptr())
                            }) == 0
                            {
                                yyval.Node = New_Func(
                                    lParse,
                                    fits_parser_yytokentype::LONG as c_int,
                                    NONNULL_FCT,
                                    1,
                                    yyvs[yyvsp - 1].Node,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                );
                                current_block = 7600445499126923600;
                            } else if (if c_int::from(yyvs[yyvsp - 2].astr[0])
                                < c_int::from((cs!(c"ACCUM("))[0])
                            {
                                -(1)
                            } else if c_int::from(yyvs[yyvsp - 2].astr[0])
                                > c_int::from((cs!(c"ACCUM("))[0])
                            {
                                1
                            } else {
                                strcmp((yyvs[yyvsp - 2].astr).as_mut_ptr(), c"ACCUM(".as_ptr())
                            }) == 0
                                && ((lParse.Nodes)[yyvs[yyvsp - 1].Node as usize]).ntype
                                    == fits_parser_yytokentype::LONG as c_int
                            {
                                let mut zero_1: c_long = 0;

                                let rc_parse: Rc<ParseData> = Rc::from_raw(lParse);

                                let new_node = New_Const(
                                    lParse,
                                    fits_parser_yytokentype::LONG as c_int,
                                    (&mut zero_1 as *mut c_long).cast::<c_void>(),
                                    ::core::mem::size_of::<c_long>() as c_ulong as c_long,
                                );

                                yyval.Node = New_BinOp(
                                    lParse,
                                    fits_parser_yytokentype::LONG as c_int,
                                    yyvs[yyvsp - 1].Node,
                                    fits_parser_yytokentype::ACCUM as c_int,
                                    new_node,
                                );
                                current_block = 7600445499126923600;
                            } else if (if c_int::from(yyvs[yyvsp - 2].astr[0])
                                < c_int::from((cs!(c"ACCUM("))[0])
                            {
                                -(1)
                            } else if c_int::from(yyvs[yyvsp - 2].astr[0])
                                > c_int::from((cs!(c"ACCUM("))[0])
                            {
                                1
                            } else {
                                strcmp((yyvs[yyvsp - 2].astr).as_mut_ptr(), c"ACCUM(".as_ptr())
                            }) == 0
                                && ((lParse.Nodes)[yyvs[yyvsp - 1].Node as usize]).ntype
                                    == fits_parser_yytokentype::DOUBLE as c_int
                            {
                                let mut zero_2: c_double = 0.0;
                                let new_node = New_Const(
                                    lParse,
                                    fits_parser_yytokentype::DOUBLE as c_int,
                                    (&mut zero_2 as *mut c_double).cast::<c_void>(),
                                    ::core::mem::size_of::<c_double>() as c_ulong as c_long,
                                );

                                yyval.Node = New_BinOp(
                                    lParse,
                                    fits_parser_yytokentype::DOUBLE as c_int,
                                    yyvs[yyvsp - 1].Node,
                                    fits_parser_yytokentype::ACCUM as c_int,
                                    new_node,
                                );
                                current_block = 7600445499126923600;
                            } else if (if c_int::from(yyvs[yyvsp - 2].astr[0])
                                < c_int::from((cs!(c"SEQDIFF("))[0])
                            {
                                -(1)
                            } else if c_int::from(yyvs[yyvsp - 2].astr[0])
                                > c_int::from((cs!(c"SEQDIFF("))[0])
                            {
                                1
                            } else {
                                strcmp((yyvs[yyvsp - 2].astr).as_mut_ptr(), c"SEQDIFF(".as_ptr())
                            }) == 0
                                && ((lParse.Nodes)[yyvs[yyvsp - 1].Node as usize]).ntype
                                    == fits_parser_yytokentype::LONG as c_int
                            {
                                let mut zero_3: c_long = 0;
                                let new_node = New_Const(
                                    lParse,
                                    fits_parser_yytokentype::LONG as c_int,
                                    (&mut zero_3 as *mut c_long).cast::<c_void>(),
                                    ::core::mem::size_of::<c_long>() as c_ulong as c_long,
                                );

                                yyval.Node = New_BinOp(
                                    lParse,
                                    fits_parser_yytokentype::LONG as c_int,
                                    yyvs[yyvsp - 1].Node,
                                    fits_parser_yytokentype::DIFF as c_int,
                                    new_node,
                                );
                                current_block = 7600445499126923600;
                            } else if (if c_int::from(yyvs[yyvsp - 2].astr[0])
                                < c_int::from((cs!(c"SEQDIFF("))[0])
                            {
                                -(1)
                            } else if c_int::from(yyvs[yyvsp - 2].astr[0])
                                > c_int::from((cs!(c"SEQDIFF("))[0])
                            {
                                1
                            } else {
                                strcmp((yyvs[yyvsp - 2].astr).as_mut_ptr(), c"SEQDIFF(".as_ptr())
                            }) == 0
                                && ((lParse.Nodes)[yyvs[yyvsp - 1].Node as usize]).ntype
                                    == fits_parser_yytokentype::DOUBLE as c_int
                            {
                                let mut zero_4: c_double = 0.0;
                                let new_node = New_Const(
                                    lParse,
                                    fits_parser_yytokentype::DOUBLE as c_int,
                                    (&mut zero_4 as *mut c_double).cast::<c_void>(),
                                    ::core::mem::size_of::<c_double>() as c_ulong as c_long,
                                );

                                yyval.Node = New_BinOp(
                                    lParse,
                                    fits_parser_yytokentype::DOUBLE as c_int,
                                    yyvs[yyvsp - 1].Node,
                                    fits_parser_yytokentype::DIFF as c_int,
                                    new_node,
                                );
                                current_block = 7600445499126923600;
                            } else if (if c_int::from(yyvs[yyvsp - 2].astr[0])
                                < c_int::from((cs!(c"ABS("))[0])
                            {
                                -(1)
                            } else if c_int::from(yyvs[yyvsp - 2].astr[0])
                                > c_int::from((cs!(c"ABS("))[0])
                            {
                                1
                            } else {
                                strcmp((yyvs[yyvsp - 2].astr).as_mut_ptr(), c"ABS(".as_ptr())
                            }) == 0
                            {
                                yyval.Node = New_Func(
                                    lParse,
                                    0,
                                    ABS_FCT,
                                    1,
                                    yyvs[yyvsp - 1].Node,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                );
                                current_block = 7600445499126923600;
                            } else if (if c_int::from(yyvs[yyvsp - 2].astr[0])
                                < c_int::from((cs!(c"MIN("))[0])
                            {
                                -(1)
                            } else if c_int::from(yyvs[yyvsp - 2].astr[0])
                                > c_int::from((cs!(c"MIN("))[0])
                            {
                                1
                            } else {
                                strcmp((yyvs[yyvsp - 2].astr).as_mut_ptr(), c"MIN(".as_ptr())
                            }) == 0
                            {
                                yyval.Node = New_Func(
                                    lParse,
                                    ((lParse.Nodes)[yyvs[yyvsp - 1].Node as usize]).ntype,
                                    MIN1_FCT,
                                    1,
                                    yyvs[yyvsp - 1].Node,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                );
                                current_block = 7600445499126923600;
                            } else if (if c_int::from(yyvs[yyvsp - 2].astr[0])
                                < c_int::from((cs!(c"MAX("))[0])
                            {
                                -(1)
                            } else if c_int::from(yyvs[yyvsp - 2].astr[0])
                                > c_int::from((cs!(c"MAX("))[0])
                            {
                                1
                            } else {
                                strcmp((yyvs[yyvsp - 2].astr).as_mut_ptr(), c"MAX(".as_ptr())
                            }) == 0
                            {
                                yyval.Node = New_Func(
                                    lParse,
                                    ((lParse.Nodes)[yyvs[yyvsp - 1].Node as usize]).ntype,
                                    MAX1_FCT,
                                    1,
                                    yyvs[yyvsp - 1].Node,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0,
                                );
                                current_block = 7600445499126923600;
                            } else if (if c_int::from(yyvs[yyvsp - 2].astr[0])
                                < c_int::from((cs!(c"RANDOM("))[0])
                            {
                                -(1)
                            } else if c_int::from(yyvs[yyvsp - 2].astr[0])
                                > c_int::from((cs!(c"RANDOM("))[0])
                            {
                                1
                            } else {
                                strcmp((yyvs[yyvsp - 2].astr).as_mut_ptr(), c"RANDOM(".as_ptr())
                            }) == 0
                            {
                                yyval.Node = New_Func(
                                    lParse,
                                    0,
                                    RND_FCT,
                                    1,
                                    yyvs[yyvsp - 1].Node,
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
                                    ((lParse.Nodes)[yyval.Node as usize]).ntype =
                                        fits_parser_yytokentype::DOUBLE as c_int;
                                    current_block = 7600445499126923600;
                                }
                            } else if (if c_int::from(yyvs[yyvsp - 2].astr[0])
                                < c_int::from((cs!(c"RANDOMN("))[0])
                            {
                                -(1)
                            } else if c_int::from(yyvs[yyvsp - 2].astr[0])
                                > c_int::from((cs!(c"RANDOMN("))[0])
                            {
                                1
                            } else {
                                strcmp((yyvs[yyvsp - 2].astr).as_mut_ptr(), c"RANDOMN(".as_ptr())
                            }) == 0
                            {
                                yyval.Node = New_Func(
                                    lParse,
                                    0,
                                    GASRND_FCT,
                                    1,
                                    yyvs[yyvsp - 1].Node,
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
                                    ((lParse.Nodes)[yyval.Node as usize]).ntype =
                                        fits_parser_yytokentype::DOUBLE as c_int;
                                    current_block = 7600445499126923600;
                                }
                            } else if (if c_int::from(yyvs[yyvsp - 2].astr[0])
                                < c_int::from((cs!(c"ELEMENTNUM"))[0])
                            {
                                -(1)
                            } else if c_int::from(yyvs[yyvsp - 2].astr[0])
                                > c_int::from((cs!(c"ELEMENTNUM"))[0])
                            {
                                1
                            } else {
                                strcmp((yyvs[yyvsp - 2].astr).as_mut_ptr(), c"ELEMENTNUM(".as_ptr())
                            }) == 0
                            {
                                if ((lParse.Nodes)[yyvs[yyvsp - 1].Node as usize]).operation
                                    == CONST_OP
                                {
                                    let mut one_1: c_long = 1;
                                    yyval.Node = New_Const(
                                        lParse,
                                        fits_parser_yytokentype::LONG as c_int,
                                        (&mut one_1 as *mut c_long).cast::<c_void>(),
                                        ::core::mem::size_of::<c_long>() as c_ulong as c_long,
                                    );
                                    current_block = 7600445499126923600;
                                } else {
                                    yyval.Node = New_Func(
                                        lParse,
                                        0,
                                        ELEMNUM_FCT,
                                        1,
                                        yyvs[yyvsp - 1].Node,
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
                                        ((lParse.Nodes)[yyval.Node as usize]).ntype =
                                            fits_parser_yytokentype::LONG as c_int;
                                        current_block = 7600445499126923600;
                                    }
                                }
                            } else if (if c_int::from(yyvs[yyvsp - 2].astr[0])
                                < c_int::from((cs!(c"NAXIS("))[0])
                            {
                                -(1)
                            } else if c_int::from(yyvs[yyvsp - 2].astr[0])
                                > c_int::from((cs!(c"NAXIS("))[0])
                            {
                                1
                            } else {
                                strcmp((yyvs[yyvsp - 2].astr).as_mut_ptr(), c"NAXIS(".as_ptr())
                            }) == 0
                            {
                                if ((lParse.Nodes)[yyvs[yyvsp - 1].Node as usize]).operation
                                    == CONST_OP
                                {
                                    let mut one_2: c_long = 1;
                                    yyval.Node = New_Const(
                                        lParse,
                                        fits_parser_yytokentype::LONG as c_int,
                                        (&mut one_2 as *mut c_long).cast::<c_void>(),
                                        ::core::mem::size_of::<c_long>() as c_ulong as c_long,
                                    );
                                    current_block = 7600445499126923600;
                                } else {
                                    let mut naxis_0: c_long = c_long::from(
                                        ((lParse.Nodes)[yyvs[yyvsp - 1].Node as usize]).value.naxis,
                                    );
                                    yyval.Node = New_Const(
                                        lParse,
                                        fits_parser_yytokentype::LONG as c_int,
                                        (&mut naxis_0 as *mut c_long).cast::<c_void>(),
                                        ::core::mem::size_of::<c_long>() as c_ulong as c_long,
                                    );
                                    if yyval.Node < 0 {
                                        current_block = 4830776507462815627;
                                    } else {
                                        current_block = 7600445499126923600;
                                    }
                                }
                            } else {
                                if ((lParse.Nodes)[yyvs[yyvsp - 1].Node as usize]).ntype
                                    != fits_parser_yytokentype::DOUBLE as c_int
                                {
                                    yyvs[yyvsp - 1].Node = New_Unary(
                                        lParse,
                                        fits_parser_yytokentype::DOUBLE as c_int,
                                        0,
                                        yyvs[yyvsp - 1].Node,
                                    );
                                }
                                if (if c_int::from(yyvs[yyvsp - 2].astr[0])
                                    < c_int::from((cs!(c"SIN("))[0])
                                {
                                    -(1)
                                } else if c_int::from(yyvs[yyvsp - 2].astr[0])
                                    > c_int::from((cs!(c"SIN("))[0])
                                {
                                    1
                                } else {
                                    strcmp((yyvs[yyvsp - 2].astr).as_mut_ptr(), c"SIN(".as_ptr())
                                }) == 0
                                {
                                    yyval.Node = New_Func(
                                        lParse,
                                        0,
                                        SIN_FCT,
                                        1,
                                        yyvs[yyvsp - 1].Node,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                    );
                                    current_block = 7600445499126923600;
                                } else if (if c_int::from(yyvs[yyvsp - 2].astr[0])
                                    < c_int::from((cs!(c"COS("))[0])
                                {
                                    -(1)
                                } else if c_int::from(yyvs[yyvsp - 2].astr[0])
                                    > c_int::from((cs!(c"COS("))[0])
                                {
                                    1
                                } else {
                                    strcmp((yyvs[yyvsp - 2].astr).as_mut_ptr(), c"COS(".as_ptr())
                                }) == 0
                                {
                                    yyval.Node = New_Func(
                                        lParse,
                                        0,
                                        COS_FCT,
                                        1,
                                        yyvs[yyvsp - 1].Node,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                    );
                                    current_block = 7600445499126923600;
                                } else if (if c_int::from(yyvs[yyvsp - 2].astr[0])
                                    < c_int::from((cs!(c"TAN("))[0])
                                {
                                    -(1)
                                } else if c_int::from(yyvs[yyvsp - 2].astr[0])
                                    > c_int::from((cs!(c"TAN("))[0])
                                {
                                    1
                                } else {
                                    strcmp((yyvs[yyvsp - 2].astr).as_mut_ptr(), c"TAN(".as_ptr())
                                }) == 0
                                {
                                    yyval.Node = New_Func(
                                        lParse,
                                        0,
                                        TAN_FCT,
                                        1,
                                        yyvs[yyvsp - 1].Node,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                    );
                                    current_block = 7600445499126923600;
                                } else if (if c_int::from(yyvs[yyvsp - 2].astr[0])
                                    < c_int::from((cs!(c"ARCSIN"))[0])
                                {
                                    -(1)
                                } else if c_int::from(yyvs[yyvsp - 2].astr[0])
                                    > c_int::from((cs!(c"ARCSIN"))[0])
                                {
                                    1
                                } else {
                                    strcmp((yyvs[yyvsp - 2].astr).as_mut_ptr(), c"ARCSIN(".as_ptr())
                                }) == 0
                                    || (if c_int::from(yyvs[yyvsp - 2].astr[0])
                                        < c_int::from((cs!(c"ASIN"))[0])
                                    {
                                        -(1)
                                    } else if c_int::from(yyvs[yyvsp - 2].astr[0])
                                        > c_int::from((cs!(c"ASIN"))[0])
                                    {
                                        1
                                    } else {
                                        strcmp(
                                            (yyvs[yyvsp - 2].astr).as_mut_ptr(),
                                            c"ASIN(".as_ptr(),
                                        )
                                    }) == 0
                                {
                                    yyval.Node = New_Func(
                                        lParse,
                                        0,
                                        ASIN_FCT,
                                        1,
                                        yyvs[yyvsp - 1].Node,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                    );
                                    current_block = 7600445499126923600;
                                } else if (if c_int::from(yyvs[yyvsp - 2].astr[0])
                                    < c_int::from((cs!(c"ARCCOS"))[0])
                                {
                                    -(1)
                                } else if c_int::from(yyvs[yyvsp - 2].astr[0])
                                    > c_int::from((cs!(c"ARCCOS"))[0])
                                {
                                    1
                                } else {
                                    strcmp((yyvs[yyvsp - 2].astr).as_mut_ptr(), c"ARCCOS(".as_ptr())
                                }) == 0
                                    || (if c_int::from(yyvs[yyvsp - 2].astr[0])
                                        < c_int::from((cs!(c"ACOS"))[0])
                                    {
                                        -(1)
                                    } else if c_int::from(yyvs[yyvsp - 2].astr[0])
                                        > c_int::from((cs!(c"ACOS"))[0])
                                    {
                                        1
                                    } else {
                                        strcmp(
                                            (yyvs[yyvsp - 2].astr).as_mut_ptr(),
                                            c"ACOS(".as_ptr(),
                                        )
                                    }) == 0
                                {
                                    yyval.Node = New_Func(
                                        lParse,
                                        0,
                                        ACOS_FCT,
                                        1,
                                        yyvs[yyvsp - 1].Node,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                    );
                                    current_block = 7600445499126923600;
                                } else if (if c_int::from(yyvs[yyvsp - 2].astr[0])
                                    < c_int::from((cs!(c"ARCTAN"))[0])
                                {
                                    -(1)
                                } else if c_int::from(yyvs[yyvsp - 2].astr[0])
                                    > c_int::from((cs!(c"ARCTAN"))[0])
                                {
                                    1
                                } else {
                                    strcmp((yyvs[yyvsp - 2].astr).as_mut_ptr(), c"ARCTAN(".as_ptr())
                                }) == 0
                                    || (if c_int::from(yyvs[yyvsp - 2].astr[0])
                                        < c_int::from((cs!(c"ATAN"))[0])
                                    {
                                        -(1)
                                    } else if c_int::from(yyvs[yyvsp - 2].astr[0])
                                        > c_int::from((cs!(c"ATAN"))[0])
                                    {
                                        1
                                    } else {
                                        strcmp(
                                            (yyvs[yyvsp - 2].astr).as_mut_ptr(),
                                            c"ATAN(".as_ptr(),
                                        )
                                    }) == 0
                                {
                                    yyval.Node = New_Func(
                                        lParse,
                                        0,
                                        ATAN_FCT,
                                        1,
                                        yyvs[yyvsp - 1].Node,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                    );
                                    current_block = 7600445499126923600;
                                } else if (if c_int::from(yyvs[yyvsp - 2].astr[0])
                                    < c_int::from((cs!(c"SINH"))[0])
                                {
                                    -(1)
                                } else if c_int::from(yyvs[yyvsp - 2].astr[0])
                                    > c_int::from((cs!(c"SINH"))[0])
                                {
                                    1
                                } else {
                                    strcmp((yyvs[yyvsp - 2].astr).as_mut_ptr(), c"SINH(".as_ptr())
                                }) == 0
                                {
                                    yyval.Node = New_Func(
                                        lParse,
                                        0,
                                        SINH_FCT,
                                        1,
                                        yyvs[yyvsp - 1].Node,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                    );
                                    current_block = 7600445499126923600;
                                } else if (if c_int::from(yyvs[yyvsp - 2].astr[0])
                                    < c_int::from((cs!(c"COSH"))[0])
                                {
                                    -(1)
                                } else if c_int::from(yyvs[yyvsp - 2].astr[0])
                                    > c_int::from((cs!(c"COSH"))[0])
                                {
                                    1
                                } else {
                                    strcmp((yyvs[yyvsp - 2].astr).as_mut_ptr(), c"COSH(".as_ptr())
                                }) == 0
                                {
                                    yyval.Node = New_Func(
                                        lParse,
                                        0,
                                        COSH_FCT,
                                        1,
                                        yyvs[yyvsp - 1].Node,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                    );
                                    current_block = 7600445499126923600;
                                } else if (if c_int::from(yyvs[yyvsp - 2].astr[0])
                                    < c_int::from((cs!(c"TANH"))[0])
                                {
                                    -(1)
                                } else if c_int::from(yyvs[yyvsp - 2].astr[0])
                                    > c_int::from((cs!(c"TANH"))[0])
                                {
                                    1
                                } else {
                                    strcmp((yyvs[yyvsp - 2].astr).as_mut_ptr(), c"TANH(".as_ptr())
                                }) == 0
                                {
                                    yyval.Node = New_Func(
                                        lParse,
                                        0,
                                        TANH_FCT,
                                        1,
                                        yyvs[yyvsp - 1].Node,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                    );
                                    current_block = 7600445499126923600;
                                } else if (if c_int::from(yyvs[yyvsp - 2].astr[0])
                                    < c_int::from((cs!(c"EXP("))[0])
                                {
                                    -(1)
                                } else if c_int::from(yyvs[yyvsp - 2].astr[0])
                                    > c_int::from((cs!(c"EXP("))[0])
                                {
                                    1
                                } else {
                                    strcmp((yyvs[yyvsp - 2].astr).as_mut_ptr(), c"EXP(".as_ptr())
                                }) == 0
                                {
                                    yyval.Node = New_Func(
                                        lParse,
                                        0,
                                        EXP_FCT,
                                        1,
                                        yyvs[yyvsp - 1].Node,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                    );
                                    current_block = 7600445499126923600;
                                } else if (if c_int::from(yyvs[yyvsp - 2].astr[0])
                                    < c_int::from((cs!(c"LOG("))[0])
                                {
                                    -(1)
                                } else if c_int::from(yyvs[yyvsp - 2].astr[0])
                                    > c_int::from((cs!(c"LOG("))[0])
                                {
                                    1
                                } else {
                                    strcmp((yyvs[yyvsp - 2].astr).as_mut_ptr(), c"LOG(".as_ptr())
                                }) == 0
                                {
                                    yyval.Node = New_Func(
                                        lParse,
                                        0,
                                        LOG_FCT,
                                        1,
                                        yyvs[yyvsp - 1].Node,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                    );
                                    current_block = 7600445499126923600;
                                } else if (if c_int::from(yyvs[yyvsp - 2].astr[0])
                                    < c_int::from((cs!(c"LOG10"))[0])
                                {
                                    -(1)
                                } else if c_int::from(yyvs[yyvsp - 2].astr[0])
                                    > c_int::from((cs!(c"LOG10"))[0])
                                {
                                    1
                                } else {
                                    strcmp((yyvs[yyvsp - 2].astr).as_mut_ptr(), c"LOG10(".as_ptr())
                                }) == 0
                                {
                                    yyval.Node = New_Func(
                                        lParse,
                                        0,
                                        LOG10_FCT,
                                        1,
                                        yyvs[yyvsp - 1].Node,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                    );
                                    current_block = 7600445499126923600;
                                } else if (if c_int::from(yyvs[yyvsp - 2].astr[0])
                                    < c_int::from((cs!(c"SQRT"))[0])
                                {
                                    -(1)
                                } else if c_int::from(yyvs[yyvsp - 2].astr[0])
                                    > c_int::from((cs!(c"SQRT"))[0])
                                {
                                    1
                                } else {
                                    strcmp((yyvs[yyvsp - 2].astr).as_mut_ptr(), c"SQRT(".as_ptr())
                                }) == 0
                                {
                                    yyval.Node = New_Func(
                                        lParse,
                                        0,
                                        SQRT_FCT,
                                        1,
                                        yyvs[yyvsp - 1].Node,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                    );
                                    current_block = 7600445499126923600;
                                } else if (if c_int::from(yyvs[yyvsp - 2].astr[0])
                                    < c_int::from((cs!(c"ROUND"))[0])
                                {
                                    -(1)
                                } else if c_int::from(yyvs[yyvsp - 2].astr[0])
                                    > c_int::from((cs!(c"ROUND"))[0])
                                {
                                    1
                                } else {
                                    strcmp((yyvs[yyvsp - 2].astr).as_mut_ptr(), c"ROUND(".as_ptr())
                                }) == 0
                                {
                                    yyval.Node = New_Func(
                                        lParse,
                                        0,
                                        ROUND_FCT,
                                        1,
                                        yyvs[yyvsp - 1].Node,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                    );
                                    current_block = 7600445499126923600;
                                } else if (if c_int::from(yyvs[yyvsp - 2].astr[0])
                                    < c_int::from((cs!(c"FLOOR"))[0])
                                {
                                    -(1)
                                } else if c_int::from(yyvs[yyvsp - 2].astr[0])
                                    > c_int::from((cs!(c"FLOOR"))[0])
                                {
                                    1
                                } else {
                                    strcmp((yyvs[yyvsp - 2].astr).as_mut_ptr(), c"FLOOR(".as_ptr())
                                }) == 0
                                {
                                    yyval.Node = New_Func(
                                        lParse,
                                        0,
                                        FLOOR_FCT,
                                        1,
                                        yyvs[yyvsp - 1].Node,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                    );
                                    current_block = 7600445499126923600;
                                } else if (if c_int::from(yyvs[yyvsp - 2].astr[0])
                                    < c_int::from((cs!(c"CEIL"))[0])
                                {
                                    -(1)
                                } else if c_int::from(yyvs[yyvsp - 2].astr[0])
                                    > c_int::from((cs!(c"CEIL"))[0])
                                {
                                    1
                                } else {
                                    strcmp((yyvs[yyvsp - 2].astr).as_mut_ptr(), c"CEIL(".as_ptr())
                                }) == 0
                                {
                                    yyval.Node = New_Func(
                                        lParse,
                                        0,
                                        CEIL_FCT,
                                        1,
                                        yyvs[yyvsp - 1].Node,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                    );
                                    current_block = 7600445499126923600;
                                } else if (if c_int::from(yyvs[yyvsp - 2].astr[0])
                                    < c_int::from((cs!(c"RANDOMP"))[0])
                                {
                                    -(1)
                                } else if c_int::from(yyvs[yyvsp - 2].astr[0])
                                    > c_int::from((cs!(c"RANDOMP"))[0])
                                {
                                    1
                                } else {
                                    strcmp(
                                        (yyvs[yyvsp - 2].astr).as_mut_ptr(),
                                        c"RANDOMP(".as_ptr(),
                                    )
                                }) == 0
                                {
                                    yyval.Node = New_Func(
                                        lParse,
                                        0,
                                        POIRND_FCT,
                                        1,
                                        yyvs[yyvsp - 1].Node,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                    );
                                    ((lParse.Nodes)[yyval.Node as usize]).ntype =
                                        fits_parser_yytokentype::LONG as c_int;
                                    current_block = 7600445499126923600;
                                } else {
                                    fits_parser_yyerror(
                                        lParse,
                                        cs!(c"Function(expr) not supported"),
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
                            /* expr: IFUNCTION sexpr ',' sexpr ')'  */
                            if (if c_int::from(yyvs[yyvsp - 4].astr[0])
                                < c_int::from((cs!(c"STRSTR("))[0])
                            {
                                -(1)
                            } else if c_int::from(yyvs[yyvsp - 4].astr[0])
                                > c_int::from((cs!(c"STRSTR("))[0])
                            {
                                1
                            } else {
                                strcmp((yyvs[yyvsp - 4].astr).as_mut_ptr(), c"STRSTR(".as_ptr())
                            }) == 0
                            {
                                yyval.Node = New_Func(
                                    lParse,
                                    fits_parser_yytokentype::LONG as c_int,
                                    STRPOS_FCT,
                                    2,
                                    yyvs[yyvsp - 3].Node,
                                    yyvs[yyvsp - 1].Node,
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
                            /* expr: FUNCTION expr ',' expr ')'  */
                            if (if c_int::from(yyvs[yyvsp - 4].astr[0])
                                < c_int::from((cs!(c"DEFNULL("))[0])
                            {
                                -(1)
                            } else if c_int::from(yyvs[yyvsp - 4].astr[0])
                                > c_int::from((cs!(c"DEFNULL("))[0])
                            {
                                1
                            } else {
                                strcmp((yyvs[yyvsp - 4].astr).as_mut_ptr(), c"DEFNULL(".as_ptr())
                            }) == 0
                            {
                                if ((lParse.Nodes)[yyvs[yyvsp - 3].Node as usize]).value.nelem
                                    >= ((lParse.Nodes)[yyvs[yyvsp - 1].Node as usize]).value.nelem
                                    && Test_Dims(lParse, yyvs[yyvsp - 3].Node, yyvs[yyvsp - 1].Node)
                                        != 0
                                {
                                    if ((lParse.Nodes)[yyvs[yyvsp - 3].Node as usize]).ntype
                                        > ((lParse.Nodes)[yyvs[yyvsp - 1].Node as usize]).ntype
                                    {
                                        yyvs[yyvsp - 1].Node = New_Unary(
                                            lParse,
                                            ((lParse.Nodes)[yyvs[yyvsp - 3].Node as usize]).ntype,
                                            0,
                                            yyvs[yyvsp - 1].Node,
                                        );
                                    } else if ((lParse.Nodes)[yyvs[yyvsp - 3].Node as usize]).ntype
                                        < ((lParse.Nodes)[yyvs[yyvsp - 1].Node as usize]).ntype
                                    {
                                        yyvs[yyvsp - 3].Node = New_Unary(
                                            lParse,
                                            ((lParse.Nodes)[yyvs[yyvsp - 1].Node as usize]).ntype,
                                            0,
                                            yyvs[yyvsp - 3].Node,
                                        );
                                    }
                                    yyval.Node = New_Func(
                                        lParse,
                                        0,
                                        DEFNULL_FCT,
                                        2,
                                        yyvs[yyvsp - 3].Node,
                                        yyvs[yyvsp - 1].Node,
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
                                        lParse,
                                        cs!(c"Dimensions of DEFNULL arguments are not compatible"),
                                    );
                                    current_block = 4830776507462815627;
                                }
                            } else if (if c_int::from(yyvs[yyvsp - 4].astr[0])
                                < c_int::from((cs!(c"ARCTAN2("))[0])
                            {
                                -(1)
                            } else if c_int::from(yyvs[yyvsp - 4].astr[0])
                                > c_int::from((cs!(c"ARCTAN2("))[0])
                            {
                                1
                            } else {
                                strcmp((yyvs[yyvsp - 4].astr).as_mut_ptr(), c"ARCTAN2(".as_ptr())
                            }) == 0
                            {
                                if ((lParse.Nodes)[yyvs[yyvsp - 3].Node as usize]).ntype
                                    != fits_parser_yytokentype::DOUBLE as c_int
                                {
                                    yyvs[yyvsp - 3].Node = New_Unary(
                                        lParse,
                                        fits_parser_yytokentype::DOUBLE as c_int,
                                        0,
                                        yyvs[yyvsp - 3].Node,
                                    );
                                }
                                if ((lParse.Nodes)[yyvs[yyvsp - 1].Node as usize]).ntype
                                    != fits_parser_yytokentype::DOUBLE as c_int
                                {
                                    yyvs[yyvsp - 1].Node = New_Unary(
                                        lParse,
                                        fits_parser_yytokentype::DOUBLE as c_int,
                                        0,
                                        yyvs[yyvsp - 1].Node,
                                    );
                                }
                                if Test_Dims(lParse, yyvs[yyvsp - 3].Node, yyvs[yyvsp - 1].Node)
                                    != 0
                                {
                                    yyval.Node = New_Func(
                                        lParse,
                                        0,
                                        ATAN2_FCT,
                                        2,
                                        yyvs[yyvsp - 3].Node,
                                        yyvs[yyvsp - 1].Node,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                    );
                                    if yyval.Node < 0 {
                                        current_block = 4830776507462815627;
                                    } else {
                                        if ((lParse.Nodes)[yyvs[yyvsp - 3].Node as usize])
                                            .value
                                            .nelem
                                            < ((lParse.Nodes)[yyvs[yyvsp - 1].Node as usize])
                                                .value
                                                .nelem
                                        {
                                            Copy_Dims(lParse, yyval.Node, yyvs[yyvsp - 1].Node);
                                        }
                                        current_block = 9966817879908499150;
                                    }
                                } else {
                                    fits_parser_yyerror(
                                        lParse,
                                        cs!(c"Dimensions of arctan2 arguments are not compatible"),
                                    );
                                    current_block = 4830776507462815627;
                                }
                            } else if (if c_int::from(yyvs[yyvsp - 4].astr[0])
                                < c_int::from((cs!(c"MIN("))[0])
                            {
                                -(1)
                            } else if c_int::from(yyvs[yyvsp - 4].astr[0])
                                > c_int::from((cs!(c"MIN("))[0])
                            {
                                1
                            } else {
                                strcmp((yyvs[yyvsp - 4].astr).as_mut_ptr(), c"MIN(".as_ptr())
                            }) == 0
                            {
                                if ((lParse.Nodes)[yyvs[yyvsp - 3].Node as usize]).ntype
                                    > ((lParse.Nodes)[yyvs[yyvsp - 1].Node as usize]).ntype
                                {
                                    yyvs[yyvsp - 1].Node = New_Unary(
                                        lParse,
                                        ((lParse.Nodes)[yyvs[yyvsp - 3].Node as usize]).ntype,
                                        0,
                                        yyvs[yyvsp - 1].Node,
                                    );
                                } else if ((lParse.Nodes)[yyvs[yyvsp - 3].Node as usize]).ntype
                                    < ((lParse.Nodes)[yyvs[yyvsp - 1].Node as usize]).ntype
                                {
                                    yyvs[yyvsp - 3].Node = New_Unary(
                                        lParse,
                                        ((lParse.Nodes)[yyvs[yyvsp - 1].Node as usize]).ntype,
                                        0,
                                        yyvs[yyvsp - 3].Node,
                                    );
                                }
                                if Test_Dims(lParse, yyvs[yyvsp - 3].Node, yyvs[yyvsp - 1].Node)
                                    != 0
                                {
                                    yyval.Node = New_Func(
                                        lParse,
                                        0,
                                        MIN2_FCT,
                                        2,
                                        yyvs[yyvsp - 3].Node,
                                        yyvs[yyvsp - 1].Node,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                    );
                                    if yyval.Node < 0 {
                                        current_block = 4830776507462815627;
                                    } else {
                                        if ((lParse.Nodes)[yyvs[yyvsp - 3].Node as usize])
                                            .value
                                            .nelem
                                            < ((lParse.Nodes)[yyvs[yyvsp - 1].Node as usize])
                                                .value
                                                .nelem
                                        {
                                            Copy_Dims(lParse, yyval.Node, yyvs[yyvsp - 1].Node);
                                        }
                                        current_block = 9966817879908499150;
                                    }
                                } else {
                                    fits_parser_yyerror(
                                        lParse,
                                        cs!(c"Dimensions of min(a,b) arguments are not compatible"),
                                    );
                                    current_block = 4830776507462815627;
                                }
                            } else if (if c_int::from(yyvs[yyvsp - 4].astr[0])
                                < c_int::from((cs!(c"MAX("))[0])
                            {
                                -(1)
                            } else if c_int::from(yyvs[yyvsp - 4].astr[0])
                                > c_int::from((cs!(c"MAX("))[0])
                            {
                                1
                            } else {
                                strcmp((yyvs[yyvsp - 4].astr).as_mut_ptr(), c"MAX(".as_ptr())
                            }) == 0
                            {
                                if ((lParse.Nodes)[yyvs[yyvsp - 3].Node as usize]).ntype
                                    > ((lParse.Nodes)[yyvs[yyvsp - 1].Node as usize]).ntype
                                {
                                    yyvs[yyvsp - 1].Node = New_Unary(
                                        lParse,
                                        ((lParse.Nodes)[yyvs[yyvsp - 3].Node as usize]).ntype,
                                        0,
                                        yyvs[yyvsp - 1].Node,
                                    );
                                } else if ((lParse.Nodes)[yyvs[yyvsp - 3].Node as usize]).ntype
                                    < ((lParse.Nodes)[yyvs[yyvsp - 1].Node as usize]).ntype
                                {
                                    yyvs[yyvsp - 3].Node = New_Unary(
                                        lParse,
                                        ((lParse.Nodes)[yyvs[yyvsp - 1].Node as usize]).ntype,
                                        0,
                                        yyvs[yyvsp - 3].Node,
                                    );
                                }
                                if Test_Dims(lParse, yyvs[yyvsp - 3].Node, yyvs[yyvsp - 1].Node)
                                    != 0
                                {
                                    yyval.Node = New_Func(
                                        lParse,
                                        0,
                                        MAX2_FCT,
                                        2,
                                        yyvs[yyvsp - 3].Node,
                                        yyvs[yyvsp - 1].Node,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                    );
                                    if yyval.Node < 0 {
                                        current_block = 4830776507462815627;
                                    } else {
                                        if ((lParse.Nodes)[yyvs[yyvsp - 3].Node as usize])
                                            .value
                                            .nelem
                                            < ((lParse.Nodes)[yyvs[yyvsp - 1].Node as usize])
                                                .value
                                                .nelem
                                        {
                                            Copy_Dims(lParse, yyval.Node, yyvs[yyvsp - 1].Node);
                                        }
                                        current_block = 9966817879908499150;
                                    }
                                } else {
                                    fits_parser_yyerror(
                                        lParse,
                                        cs!(c"Dimensions of max(a,b) arguments are not compatible"),
                                    );
                                    current_block = 4830776507462815627;
                                }
                            } else if (if c_int::from(yyvs[yyvsp - 4].astr[0])
                                < c_int::from((cs!(c"SETNULL("))[0])
                            {
                                -(1)
                            } else if c_int::from(yyvs[yyvsp - 4].astr[0])
                                > c_int::from((cs!(c"SETNULL("))[0])
                            {
                                1
                            } else {
                                strcmp((yyvs[yyvsp - 4].astr).as_mut_ptr(), c"SETNULL(".as_ptr())
                            }) == 0
                            {
                                if ((lParse.Nodes)[yyvs[yyvsp - 3].Node as usize]).operation
                                    != CONST_OP
                                    || ((lParse.Nodes)[yyvs[yyvsp - 3].Node as usize]).value.nelem
                                        != 1
                                {
                                    fits_parser_yyerror(
                                        lParse,
                                        cs!(c"SETNULL first argument must be a scalar constant"),
                                    );
                                    current_block = 4830776507462815627;
                                } else {
                                    if ((lParse.Nodes)[yyvs[yyvsp - 3].Node as usize]).ntype
                                        != ((lParse.Nodes)[yyvs[yyvsp - 1].Node as usize]).ntype
                                    {
                                        yyvs[yyvsp - 3].Node = New_Unary(
                                            lParse,
                                            ((lParse.Nodes)[yyvs[yyvsp - 1].Node as usize]).ntype,
                                            0,
                                            yyvs[yyvsp - 3].Node,
                                        );
                                    }
                                    yyval.Node = New_Func(
                                        lParse,
                                        0,
                                        SETNULL_FCT,
                                        2,
                                        yyvs[yyvsp - 1].Node,
                                        yyvs[yyvsp - 3].Node,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                    );
                                    current_block = 9966817879908499150;
                                }
                            } else if (if c_int::from(yyvs[yyvsp - 4].astr[0])
                                < c_int::from((cs!(c"AXISELEM"))[0])
                            {
                                -(1)
                            } else if c_int::from(yyvs[yyvsp - 4].astr[0])
                                > c_int::from((cs!(c"AXISELEM"))[0])
                            {
                                1
                            } else {
                                strcmp((yyvs[yyvsp - 4].astr).as_mut_ptr(), c"AXISELEM(".as_ptr())
                            }) == 0
                            {
                                if ((lParse.Nodes)[yyvs[yyvsp - 1].Node as usize]).operation
                                    != CONST_OP
                                    || ((lParse.Nodes)[yyvs[yyvsp - 1].Node as usize]).value.nelem
                                        != 1
                                {
                                    fits_parser_yyerror(
                                        lParse,
                                        cs!(c"AXISELEM second argument must be a scalar constant"),
                                    );
                                    current_block = 4830776507462815627;
                                } else if ((lParse.Nodes)[yyvs[yyvsp - 3].Node as usize]).operation
                                    == CONST_OP
                                {
                                    let mut one_3: c_long = 1;
                                    yyval.Node = New_Const(
                                        lParse,
                                        fits_parser_yytokentype::LONG as c_int,
                                        (&mut one_3 as *mut c_long).cast::<c_void>(),
                                        ::core::mem::size_of::<c_long>() as c_ulong as c_long,
                                    );
                                    current_block = 9966817879908499150;
                                } else {
                                    if ((lParse.Nodes)[yyvs[yyvsp - 1].Node as usize]).ntype
                                        != fits_parser_yytokentype::LONG as c_int
                                    {
                                        yyvs[yyvsp - 1].Node = New_Unary(
                                            lParse,
                                            fits_parser_yytokentype::LONG as c_int,
                                            0,
                                            yyvs[yyvsp - 1].Node,
                                        );
                                    }
                                    yyval.Node = New_Func(
                                        lParse,
                                        0,
                                        AXISELEM_FCT,
                                        2,
                                        yyvs[yyvsp - 3].Node,
                                        yyvs[yyvsp - 1].Node,
                                        0,
                                        0,
                                        0,
                                        0,
                                        0,
                                    );
                                    if yyval.Node < 0 {
                                        current_block = 4830776507462815627;
                                    } else {
                                        ((lParse.Nodes)[yyval.Node as usize]).ntype =
                                            fits_parser_yytokentype::LONG as c_int;
                                        current_block = 9966817879908499150;
                                    }
                                }
                            } else if (if c_int::from(yyvs[yyvsp - 4].astr[0])
                                < c_int::from((cs!(c"NAXES("))[0])
                            {
                                -(1)
                            } else if c_int::from(yyvs[yyvsp - 4].astr[0])
                                > c_int::from((cs!(c"NAXES("))[0])
                            {
                                1
                            } else {
                                strcmp((yyvs[yyvsp - 4].astr).as_mut_ptr(), c"NAXES(".as_ptr())
                            }) == 0
                            {
                                if ((lParse.Nodes)[yyvs[yyvsp - 1].Node as usize]).operation
                                    != CONST_OP
                                    || ((lParse.Nodes)[yyvs[yyvsp - 1].Node as usize]).value.nelem
                                        != 1
                                {
                                    fits_parser_yyerror(
                                        lParse,
                                        cs!(c"NAXES second argument must be a scalar constant"),
                                    );
                                    current_block = 4830776507462815627;
                                } else if ((lParse.Nodes)[yyvs[yyvsp - 3].Node as usize]).operation
                                    == CONST_OP
                                {
                                    let mut one_4: c_long = 1;
                                    yyval.Node = New_Const(
                                        lParse,
                                        fits_parser_yytokentype::LONG as c_int,
                                        (&mut one_4 as *mut c_long).cast::<c_void>(),
                                        ::core::mem::size_of::<c_long>() as c_ulong as c_long,
                                    );
                                    current_block = 9966817879908499150;
                                } else {
                                    let mut iaxis_0: c_long = 0;
                                    let mut naxis_1: c_int = 0;
                                    if ((lParse.Nodes)[yyvs[yyvsp - 1].Node as usize]).ntype
                                        != fits_parser_yytokentype::LONG as c_int
                                    {
                                        yyvs[yyvsp - 1].Node = New_Unary(
                                            lParse,
                                            fits_parser_yytokentype::LONG as c_int,
                                            0,
                                            yyvs[yyvsp - 1].Node,
                                        );
                                    }
                                    iaxis_0 = ((lParse.Nodes)[yyvs[yyvsp - 1].Node as usize])
                                        .value
                                        .data
                                        .lng;
                                    naxis_1 =
                                        ((lParse.Nodes)[yyvs[yyvsp - 3].Node as usize]).value.naxis;
                                    if iaxis_0 == 0 {
                                        iaxis_0 = c_long::from(naxis_1);
                                    } else if iaxis_0 <= c_long::from(naxis_1) {
                                        iaxis_0 = ((lParse.Nodes)[yyvs[yyvsp - 3].Node as usize])
                                            .value
                                            .naxes
                                            [(iaxis_0 - 1) as usize];
                                    } else {
                                        iaxis_0 = 1;
                                    }
                                    yyval.Node = New_Const(
                                        lParse,
                                        fits_parser_yytokentype::LONG as c_int,
                                        (&mut iaxis_0 as *mut c_long).cast::<c_void>(),
                                        ::core::mem::size_of::<c_long>() as c_ulong as c_long,
                                    );
                                    if yyval.Node < 0 {
                                        current_block = 4830776507462815627;
                                    } else {
                                        current_block = 9966817879908499150;
                                    }
                                }
                            } else if (if c_int::from(yyvs[yyvsp - 4].astr[0])
                                < c_int::from((cs!(c"ARRAY("))[0])
                            {
                                -(1)
                            } else if c_int::from(yyvs[yyvsp - 4].astr[0])
                                > c_int::from((cs!(c"ARRAY("))[0])
                            {
                                1
                            } else {
                                strcmp((yyvs[yyvsp - 4].astr).as_mut_ptr(), c"ARRAY(".as_ptr())
                            }) == 0
                            {
                                yyval.Node =
                                    New_Array(lParse, yyvs[yyvsp - 3].Node, yyvs[yyvsp - 1].Node);
                                if yyval.Node < 0 {
                                    current_block = 4830776507462815627;
                                } else {
                                    current_block = 9966817879908499150;
                                }
                            } else {
                                fits_parser_yyerror(
                                    lParse,
                                    cs!(c"Function(expr,expr) not supported"),
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
                            /* expr: FUNCTION expr ',' expr ',' expr ',' expr ')'  */
                            if (if c_int::from(yyvs[yyvsp - 8].astr[0])
                                < c_int::from((cs!(c"ANGSEP("))[0])
                            {
                                -(1)
                            } else if c_int::from(yyvs[yyvsp - 8].astr[0])
                                > c_int::from((cs!(c"ANGSEP("))[0])
                            {
                                1
                            } else {
                                strcmp((yyvs[yyvsp - 8].astr).as_mut_ptr(), c"ANGSEP(".as_ptr())
                            }) == 0
                            {
                                if ((lParse.Nodes)[yyvs[yyvsp - 7].Node as usize]).ntype
                                    != fits_parser_yytokentype::DOUBLE as c_int
                                {
                                    yyvs[yyvsp - 7].Node = New_Unary(
                                        lParse,
                                        fits_parser_yytokentype::DOUBLE as c_int,
                                        0,
                                        yyvs[yyvsp - 7].Node,
                                    );
                                }
                                if ((lParse.Nodes)[yyvs[yyvsp - 5].Node as usize]).ntype
                                    != fits_parser_yytokentype::DOUBLE as c_int
                                {
                                    yyvs[yyvsp - 5].Node = New_Unary(
                                        lParse,
                                        fits_parser_yytokentype::DOUBLE as c_int,
                                        0,
                                        yyvs[yyvsp - 5].Node,
                                    );
                                }
                                if ((lParse.Nodes)[yyvs[yyvsp - 3].Node as usize]).ntype
                                    != fits_parser_yytokentype::DOUBLE as c_int
                                {
                                    yyvs[yyvsp - 3].Node = New_Unary(
                                        lParse,
                                        fits_parser_yytokentype::DOUBLE as c_int,
                                        0,
                                        yyvs[yyvsp - 3].Node,
                                    );
                                }
                                if ((lParse.Nodes)[yyvs[yyvsp - 1].Node as usize]).ntype
                                    != fits_parser_yytokentype::DOUBLE as c_int
                                {
                                    yyvs[yyvsp - 1].Node = New_Unary(
                                        lParse,
                                        fits_parser_yytokentype::DOUBLE as c_int,
                                        0,
                                        yyvs[yyvsp - 1].Node,
                                    );
                                }
                                if Test_Dims(lParse, yyvs[yyvsp - 7].Node, yyvs[yyvsp - 5].Node)
                                    != 0
                                    && Test_Dims(lParse, yyvs[yyvsp - 5].Node, yyvs[yyvsp - 3].Node)
                                        != 0
                                    && Test_Dims(lParse, yyvs[yyvsp - 3].Node, yyvs[yyvsp - 1].Node)
                                        != 0
                                {
                                    yyval.Node = New_Func(
                                        lParse,
                                        0,
                                        ANGSEP_FCT,
                                        4 as c_int,
                                        yyvs[yyvsp - 7].Node,
                                        yyvs[yyvsp - 5].Node,
                                        yyvs[yyvsp - 3].Node,
                                        yyvs[yyvsp - 1].Node,
                                        0,
                                        0,
                                        0,
                                    );
                                    if yyval.Node < 0 {
                                        current_block = 4830776507462815627;
                                    } else {
                                        if ((lParse.Nodes)[yyvs[yyvsp - 7].Node as usize])
                                            .value
                                            .nelem
                                            < ((lParse.Nodes)[yyvs[yyvsp - 5].Node as usize])
                                                .value
                                                .nelem
                                        {
                                            Copy_Dims(lParse, yyval.Node, yyvs[yyvsp - 5].Node);
                                        }
                                        if ((lParse.Nodes)[yyvs[yyvsp - 5].Node as usize])
                                            .value
                                            .nelem
                                            < ((lParse.Nodes)[yyvs[yyvsp - 3].Node as usize])
                                                .value
                                                .nelem
                                        {
                                            Copy_Dims(lParse, yyval.Node, yyvs[yyvsp - 3].Node);
                                        }
                                        if ((lParse.Nodes)[yyvs[yyvsp - 3].Node as usize])
                                            .value
                                            .nelem
                                            < ((lParse.Nodes)[yyvs[yyvsp - 1].Node as usize])
                                                .value
                                                .nelem
                                        {
                                            Copy_Dims(lParse, yyval.Node, yyvs[yyvsp - 1].Node);
                                        }
                                        current_block = 17353983478346836848;
                                    }
                                } else {
                                    fits_parser_yyerror(
                                        lParse,
                                        cs!(c"Dimensions of ANGSEP arguments are not compatible"),
                                    );
                                    current_block = 4830776507462815627;
                                }
                            } else {
                                fits_parser_yyerror(
                                    lParse,
                                    cs!(c"Function(expr,expr,expr,expr) not supported"),
                                );
                                current_block = 4830776507462815627;
                            }
                        }
                        63 => {
                            /* expr: GTIOVERLAP STRING ',' expr ',' expr ')'  */
                            yyval.Node = New_GTI(
                                lParse,
                                GTIOVER_FCT,
                                (yyvs[yyvsp - 5].astr).as_mut_ptr(),
                                yyvs[yyvsp - 3].Node,
                                yyvs[yyvsp - 1].Node,
                                c"*START*".as_ptr() as *mut c_char,
                                c"*STOP*".as_ptr() as *mut c_char,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        64 => {
                            /* expr: GTIOVERLAP STRING ',' expr ',' expr ',' STRING ',' STRING ')'  */
                            yyval.Node = New_GTI(
                                lParse,
                                GTIOVER_FCT,
                                (yyvs[yyvsp - 9].astr).as_mut_ptr(),
                                yyvs[yyvsp - 7].Node,
                                yyvs[yyvsp - 5].Node,
                                (yyvs[yyvsp - 3].astr).as_mut_ptr(),
                                (yyvs[yyvsp - 1].astr).as_mut_ptr(),
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        65 => {
                            /* expr: expr '[' expr ']'  */
                            yyval.Node = New_Deref(
                                lParse,
                                yyvs[yyvsp - 3].Node,
                                1,
                                yyvs[yyvsp - 1].Node,
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
                            /* expr: expr '[' expr ',' expr ']'  */
                            yyval.Node = New_Deref(
                                lParse,
                                yyvs[yyvsp - 5].Node,
                                2,
                                yyvs[yyvsp - 3].Node,
                                yyvs[yyvsp - 1].Node,
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
                            /* expr: expr '[' expr ',' expr ',' expr ']'  */
                            yyval.Node = New_Deref(
                                lParse,
                                yyvs[yyvsp - 7].Node,
                                3 as c_int,
                                yyvs[yyvsp - 5].Node,
                                yyvs[yyvsp - 3].Node,
                                yyvs[yyvsp - 1].Node,
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
                            /* expr: expr '[' expr ',' expr ',' expr ',' expr ']'  */
                            yyval.Node = New_Deref(
                                lParse,
                                yyvs[yyvsp - 9].Node,
                                4 as c_int,
                                yyvs[yyvsp - 7].Node,
                                yyvs[yyvsp - 5].Node,
                                yyvs[yyvsp - 3].Node,
                                yyvs[yyvsp - 1].Node,
                                0,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        69 => {
                            /* expr: expr '[' expr ',' expr ',' expr ',' expr ',' expr ']'  */
                            yyval.Node = New_Deref(
                                lParse,
                                yyvs[yyvsp - 11].Node,
                                5 as c_int,
                                yyvs[yyvsp - 9].Node,
                                yyvs[yyvsp - 7].Node,
                                yyvs[yyvsp - 5].Node,
                                yyvs[yyvsp - 3].Node,
                                yyvs[yyvsp - 1].Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        70 => {
                            /* expr: INTCAST expr  */
                            yyval.Node = New_Unary(
                                lParse,
                                fits_parser_yytokentype::LONG as c_int,
                                fits_parser_yytokentype::INTCAST as c_int,
                                yyvs[yyvsp].Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        71 => {
                            /* expr: INTCAST bexpr  */
                            yyval.Node = New_Unary(
                                lParse,
                                fits_parser_yytokentype::LONG as c_int,
                                fits_parser_yytokentype::INTCAST as c_int,
                                yyvs[yyvsp].Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        72 => {
                            /* expr: FLTCAST expr  */
                            yyval.Node = New_Unary(
                                lParse,
                                fits_parser_yytokentype::DOUBLE as c_int,
                                fits_parser_yytokentype::FLTCAST as c_int,
                                yyvs[yyvsp].Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        73 => {
                            /* expr: FLTCAST bexpr  */
                            yyval.Node = New_Unary(
                                lParse,
                                fits_parser_yytokentype::DOUBLE as c_int,
                                fits_parser_yytokentype::FLTCAST as c_int,
                                yyvs[yyvsp].Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        74 => {
                            /* bexpr: BOOLEAN  */
                            yyval.Node = New_Const(
                                lParse,
                                fits_parser_yytokentype::BOOLEAN as c_int,
                                (&mut yyvs[yyvsp].log as *mut c_char).cast::<c_void>(),
                                ::core::mem::size_of::<c_char>() as c_ulong as c_long,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        75 => {
                            /* bexpr: BCOLUMN  */
                            yyval.Node = New_Column(lParse, yyvs[yyvsp].lng as c_int);
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        76 => {
                            /* bexpr: BCOLUMN '{' expr '}'  */
                            if ((lParse.Nodes)[yyvs[yyvsp - 1].Node as usize]).ntype
                                != fits_parser_yytokentype::LONG as c_int
                                || ((lParse.Nodes)[yyvs[yyvsp - 1].Node as usize]).operation
                                    != CONST_OP
                            {
                                fits_parser_yyerror(
                                    lParse,
                                    cs!(c"Offset argument must be a constant integer"),
                                );
                                current_block = 4830776507462815627;
                            } else {
                                yyval.Node = New_Offset(
                                    lParse,
                                    yyvs[yyvsp - 3].lng as c_int,
                                    yyvs[yyvsp - 1].Node,
                                );
                                if yyval.Node < 0 {
                                    current_block = 4830776507462815627;
                                } else {
                                    current_block = 17353983478346836848;
                                }
                            }
                        }
                        77 => {
                            /* bexpr: bits EQ bits  */
                            yyval.Node = New_BinOp(
                                lParse,
                                fits_parser_yytokentype::BOOLEAN as c_int,
                                yyvs[yyvsp - 2].Node,
                                fits_parser_yytokentype::EQ as c_int,
                                yyvs[yyvsp].Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                (((lParse.Nodes)[yyval.Node as usize]).value).nelem = 1;
                                current_block = 17353983478346836848;
                            }
                        }
                        78 => {
                            /* bexpr: bits NE bits  */
                            yyval.Node = New_BinOp(
                                lParse,
                                fits_parser_yytokentype::BOOLEAN as c_int,
                                yyvs[yyvsp - 2].Node,
                                fits_parser_yytokentype::NE as c_int,
                                yyvs[yyvsp].Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                (((lParse.Nodes)[yyval.Node as usize]).value).nelem = 1;
                                current_block = 17353983478346836848;
                            }
                        }
                        79 => {
                            /* bexpr: bits LT bits  */
                            yyval.Node = New_BinOp(
                                lParse,
                                fits_parser_yytokentype::BOOLEAN as c_int,
                                yyvs[yyvsp - 2].Node,
                                fits_parser_yytokentype::LT as c_int,
                                yyvs[yyvsp].Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                (((lParse.Nodes)[yyval.Node as usize]).value).nelem = 1;
                                current_block = 17353983478346836848;
                            }
                        }
                        80 => {
                            /* bexpr: bits LTE bits  */
                            yyval.Node = New_BinOp(
                                lParse,
                                fits_parser_yytokentype::BOOLEAN as c_int,
                                yyvs[yyvsp - 2].Node,
                                fits_parser_yytokentype::LTE as c_int,
                                yyvs[yyvsp].Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                (((lParse.Nodes)[yyval.Node as usize]).value).nelem = 1;
                                current_block = 17353983478346836848;
                            }
                        }
                        81 => {
                            /* bexpr: bits GT bits  */
                            yyval.Node = New_BinOp(
                                lParse,
                                fits_parser_yytokentype::BOOLEAN as c_int,
                                yyvs[yyvsp - 2].Node,
                                fits_parser_yytokentype::GT as c_int,
                                yyvs[yyvsp].Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                (((lParse.Nodes)[yyval.Node as usize]).value).nelem = 1;
                                current_block = 17353983478346836848;
                            }
                        }
                        82 => {
                            /* bexpr: bits GTE bits  */
                            yyval.Node = New_BinOp(
                                lParse,
                                fits_parser_yytokentype::BOOLEAN as c_int,
                                yyvs[yyvsp - 2].Node,
                                fits_parser_yytokentype::GTE as c_int,
                                yyvs[yyvsp].Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                (((lParse.Nodes)[yyval.Node as usize]).value).nelem = 1;
                                current_block = 17353983478346836848;
                            }
                        }
                        83 => {
                            /* bexpr: expr GT expr  */
                            if (lParse.Nodes[yyvs[yyvsp - 2].Node as usize]).ntype
                                > (lParse.Nodes[yyvs[yyvsp].Node as usize]).ntype
                            {
                                yyvs[yyvsp].Node = New_Unary(
                                    lParse,
                                    ((lParse.Nodes)[yyvs[yyvsp - 2].Node as usize]).ntype,
                                    0,
                                    yyvs[yyvsp].Node,
                                );
                            } else if ((lParse.Nodes)[yyvs[yyvsp - 2].Node as usize]).ntype
                                < (lParse.Nodes[yyvs[yyvsp].Node as usize]).ntype
                            {
                                yyvs[yyvsp - 2].Node = New_Unary(
                                    lParse,
                                    (lParse.Nodes[yyvs[yyvsp].Node as usize]).ntype,
                                    0,
                                    yyvs[yyvsp - 2].Node,
                                );
                            }
                            yyval.Node = New_BinOp(
                                lParse,
                                fits_parser_yytokentype::BOOLEAN as c_int,
                                yyvs[yyvsp - 2].Node,
                                fits_parser_yytokentype::GT as c_int,
                                yyvs[yyvsp].Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        84 => {
                            /* bexpr: expr LT expr  */
                            if (lParse.Nodes[yyvs[yyvsp - 2].Node as usize]).ntype
                                > (lParse.Nodes[yyvs[yyvsp].Node as usize]).ntype
                            {
                                yyvs[yyvsp].Node = New_Unary(
                                    lParse,
                                    ((lParse.Nodes)[yyvs[yyvsp - 2].Node as usize]).ntype,
                                    0,
                                    yyvs[yyvsp].Node,
                                );
                            } else if ((lParse.Nodes)[yyvs[yyvsp - 2].Node as usize]).ntype
                                < (lParse.Nodes[yyvs[yyvsp].Node as usize]).ntype
                            {
                                yyvs[yyvsp - 2].Node = New_Unary(
                                    lParse,
                                    (lParse.Nodes[yyvs[yyvsp].Node as usize]).ntype,
                                    0,
                                    yyvs[yyvsp - 2].Node,
                                );
                            }
                            yyval.Node = New_BinOp(
                                lParse,
                                fits_parser_yytokentype::BOOLEAN as c_int,
                                yyvs[yyvsp - 2].Node,
                                fits_parser_yytokentype::LT as c_int,
                                yyvs[yyvsp].Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        85 => {
                            /* bexpr: expr GTE expr  */
                            if (lParse.Nodes[yyvs[yyvsp - 2].Node as usize]).ntype
                                > (lParse.Nodes[yyvs[yyvsp].Node as usize]).ntype
                            {
                                yyvs[yyvsp].Node = New_Unary(
                                    lParse,
                                    ((lParse.Nodes)[yyvs[yyvsp - 2].Node as usize]).ntype,
                                    0,
                                    yyvs[yyvsp].Node,
                                );
                            } else if ((lParse.Nodes)[yyvs[yyvsp - 2].Node as usize]).ntype
                                < (lParse.Nodes[yyvs[yyvsp].Node as usize]).ntype
                            {
                                yyvs[yyvsp - 2].Node = New_Unary(
                                    lParse,
                                    (lParse.Nodes[yyvs[yyvsp].Node as usize]).ntype,
                                    0,
                                    yyvs[yyvsp - 2].Node,
                                );
                            }
                            yyval.Node = New_BinOp(
                                lParse,
                                fits_parser_yytokentype::BOOLEAN as c_int,
                                yyvs[yyvsp - 2].Node,
                                fits_parser_yytokentype::GTE as c_int,
                                yyvs[yyvsp].Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        86 => {
                            /* bexpr: expr LTE expr  */
                            if (lParse.Nodes[yyvs[yyvsp - 2].Node as usize]).ntype
                                > (lParse.Nodes[yyvs[yyvsp].Node as usize]).ntype
                            {
                                yyvs[yyvsp].Node = New_Unary(
                                    lParse,
                                    ((lParse.Nodes)[yyvs[yyvsp - 2].Node as usize]).ntype,
                                    0,
                                    yyvs[yyvsp].Node,
                                );
                            } else if ((lParse.Nodes)[yyvs[yyvsp - 2].Node as usize]).ntype
                                < (lParse.Nodes[yyvs[yyvsp].Node as usize]).ntype
                            {
                                yyvs[yyvsp - 2].Node = New_Unary(
                                    lParse,
                                    (lParse.Nodes[yyvs[yyvsp].Node as usize]).ntype,
                                    0,
                                    yyvs[yyvsp - 2].Node,
                                );
                            }
                            yyval.Node = New_BinOp(
                                lParse,
                                fits_parser_yytokentype::BOOLEAN as c_int,
                                yyvs[yyvsp - 2].Node,
                                fits_parser_yytokentype::LTE as c_int,
                                yyvs[yyvsp].Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        87 => {
                            /* bexpr: expr '~' expr  */
                            if (lParse.Nodes[yyvs[yyvsp - 2].Node as usize]).ntype
                                > (lParse.Nodes[yyvs[yyvsp].Node as usize]).ntype
                            {
                                yyvs[yyvsp].Node = New_Unary(
                                    lParse,
                                    ((lParse.Nodes)[yyvs[yyvsp - 2].Node as usize]).ntype,
                                    0,
                                    yyvs[yyvsp].Node,
                                );
                            } else if ((lParse.Nodes)[yyvs[yyvsp - 2].Node as usize]).ntype
                                < (lParse.Nodes[yyvs[yyvsp].Node as usize]).ntype
                            {
                                yyvs[yyvsp - 2].Node = New_Unary(
                                    lParse,
                                    (lParse.Nodes[yyvs[yyvsp].Node as usize]).ntype,
                                    0,
                                    yyvs[yyvsp - 2].Node,
                                );
                            }
                            yyval.Node = New_BinOp(
                                lParse,
                                fits_parser_yytokentype::BOOLEAN as c_int,
                                yyvs[yyvsp - 2].Node,
                                '~' as i32,
                                yyvs[yyvsp].Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        88 => {
                            /* bexpr: expr EQ expr  */
                            if (lParse.Nodes[yyvs[yyvsp - 2].Node as usize]).ntype
                                > (lParse.Nodes[yyvs[yyvsp].Node as usize]).ntype
                            {
                                yyvs[yyvsp].Node = New_Unary(
                                    lParse,
                                    ((lParse.Nodes)[yyvs[yyvsp - 2].Node as usize]).ntype,
                                    0,
                                    yyvs[yyvsp].Node,
                                );
                            } else if ((lParse.Nodes)[yyvs[yyvsp - 2].Node as usize]).ntype
                                < (lParse.Nodes[yyvs[yyvsp].Node as usize]).ntype
                            {
                                yyvs[yyvsp - 2].Node = New_Unary(
                                    lParse,
                                    (lParse.Nodes[yyvs[yyvsp].Node as usize]).ntype,
                                    0,
                                    yyvs[yyvsp - 2].Node,
                                );
                            }
                            yyval.Node = New_BinOp(
                                lParse,
                                fits_parser_yytokentype::BOOLEAN as c_int,
                                yyvs[yyvsp - 2].Node,
                                fits_parser_yytokentype::EQ as c_int,
                                yyvs[yyvsp].Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        89 => {
                            /* bexpr: expr NE expr  */
                            if (lParse.Nodes[yyvs[yyvsp - 2].Node as usize]).ntype
                                > (lParse.Nodes[yyvs[yyvsp].Node as usize]).ntype
                            {
                                yyvs[yyvsp].Node = New_Unary(
                                    lParse,
                                    ((lParse.Nodes)[yyvs[yyvsp - 2].Node as usize]).ntype,
                                    0,
                                    yyvs[yyvsp].Node,
                                );
                            } else if ((lParse.Nodes)[yyvs[yyvsp - 2].Node as usize]).ntype
                                < (lParse.Nodes[yyvs[yyvsp].Node as usize]).ntype
                            {
                                yyvs[yyvsp - 2].Node = New_Unary(
                                    lParse,
                                    (lParse.Nodes[yyvs[yyvsp].Node as usize]).ntype,
                                    0,
                                    yyvs[yyvsp - 2].Node,
                                );
                            }
                            yyval.Node = New_BinOp(
                                lParse,
                                fits_parser_yytokentype::BOOLEAN as c_int,
                                yyvs[yyvsp - 2].Node,
                                fits_parser_yytokentype::NE as c_int,
                                yyvs[yyvsp].Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        90 => {
                            /* bexpr: sexpr EQ sexpr  */
                            yyval.Node = New_BinOp(
                                lParse,
                                fits_parser_yytokentype::BOOLEAN as c_int,
                                yyvs[yyvsp - 2].Node,
                                fits_parser_yytokentype::EQ as c_int,
                                yyvs[yyvsp].Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                (((lParse.Nodes)[yyval.Node as usize]).value).nelem = 1;
                                current_block = 17353983478346836848;
                            }
                        }
                        91 => {
                            /* bexpr: sexpr NE sexpr  */
                            yyval.Node = New_BinOp(
                                lParse,
                                fits_parser_yytokentype::BOOLEAN as c_int,
                                yyvs[yyvsp - 2].Node,
                                fits_parser_yytokentype::NE as c_int,
                                yyvs[yyvsp].Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                (((lParse.Nodes)[yyval.Node as usize]).value).nelem = 1;
                                current_block = 17353983478346836848;
                            }
                        }
                        92 => {
                            /* bexpr: sexpr GT sexpr  */
                            yyval.Node = New_BinOp(
                                lParse,
                                fits_parser_yytokentype::BOOLEAN as c_int,
                                yyvs[yyvsp - 2].Node,
                                fits_parser_yytokentype::GT as c_int,
                                yyvs[yyvsp].Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                (((lParse.Nodes)[yyval.Node as usize]).value).nelem = 1;
                                current_block = 17353983478346836848;
                            }
                        }
                        93 => {
                            /* bexpr: sexpr GTE sexpr  */
                            yyval.Node = New_BinOp(
                                lParse,
                                fits_parser_yytokentype::BOOLEAN as c_int,
                                yyvs[yyvsp - 2].Node,
                                fits_parser_yytokentype::GTE as c_int,
                                yyvs[yyvsp].Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                (((lParse.Nodes)[yyval.Node as usize]).value).nelem = 1;
                                current_block = 17353983478346836848;
                            }
                        }
                        94 => {
                            /* bexpr: sexpr LT sexpr  */
                            yyval.Node = New_BinOp(
                                lParse,
                                fits_parser_yytokentype::BOOLEAN as c_int,
                                yyvs[yyvsp - 2].Node,
                                fits_parser_yytokentype::LT as c_int,
                                yyvs[yyvsp].Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                (((lParse.Nodes)[yyval.Node as usize]).value).nelem = 1;
                                current_block = 17353983478346836848;
                            }
                        }
                        95 => {
                            /* bexpr: sexpr LTE sexpr  */
                            yyval.Node = New_BinOp(
                                lParse,
                                fits_parser_yytokentype::BOOLEAN as c_int,
                                yyvs[yyvsp - 2].Node,
                                fits_parser_yytokentype::LTE as c_int,
                                yyvs[yyvsp].Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                (((lParse.Nodes)[yyval.Node as usize]).value).nelem = 1;
                                current_block = 17353983478346836848;
                            }
                        }
                        96 => {
                            /* bexpr: bexpr AND bexpr  */
                            yyval.Node = New_BinOp(
                                lParse,
                                fits_parser_yytokentype::BOOLEAN as c_int,
                                yyvs[yyvsp - 2].Node,
                                fits_parser_yytokentype::AND as c_int,
                                yyvs[yyvsp].Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        97 => {
                            /* bexpr: bexpr OR bexpr  */
                            yyval.Node = New_BinOp(
                                lParse,
                                fits_parser_yytokentype::BOOLEAN as c_int,
                                yyvs[yyvsp - 2].Node,
                                fits_parser_yytokentype::OR as c_int,
                                yyvs[yyvsp].Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        98 => {
                            /* bexpr: bexpr EQ bexpr  */
                            yyval.Node = New_BinOp(
                                lParse,
                                fits_parser_yytokentype::BOOLEAN as c_int,
                                yyvs[yyvsp - 2].Node,
                                fits_parser_yytokentype::EQ as c_int,
                                yyvs[yyvsp].Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        99 => {
                            /* bexpr: bexpr NE bexpr  */
                            yyval.Node = New_BinOp(
                                lParse,
                                fits_parser_yytokentype::BOOLEAN as c_int,
                                yyvs[yyvsp - 2].Node,
                                fits_parser_yytokentype::NE as c_int,
                                yyvs[yyvsp].Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        100 => {
                            /* bexpr: expr '=' expr ':' expr  */
                            if ((lParse.Nodes)[yyvs[yyvsp - 4].Node as usize]).ntype
                                > (lParse.Nodes[yyvs[yyvsp - 2].Node as usize]).ntype
                            {
                                yyvs[yyvsp - 2].Node = New_Unary(
                                    lParse,
                                    ((lParse.Nodes)[yyvs[yyvsp - 4].Node as usize]).ntype,
                                    0,
                                    yyvs[yyvsp - 2].Node,
                                );
                            } else if ((lParse.Nodes)[yyvs[yyvsp - 4].Node as usize]).ntype
                                < (lParse.Nodes[yyvs[yyvsp - 2].Node as usize]).ntype
                            {
                                yyvs[yyvsp - 4].Node = New_Unary(
                                    lParse,
                                    ((lParse.Nodes)[yyvs[yyvsp - 2].Node as usize]).ntype,
                                    0,
                                    yyvs[yyvsp - 4].Node,
                                );
                            }
                            if ((lParse.Nodes)[yyvs[yyvsp - 4].Node as usize]).ntype
                                > (lParse.Nodes[yyvs[yyvsp].Node as usize]).ntype
                            {
                                yyvs[yyvsp].Node = New_Unary(
                                    lParse,
                                    ((lParse.Nodes)[yyvs[yyvsp - 4].Node as usize]).ntype,
                                    0,
                                    yyvs[yyvsp].Node,
                                );
                            } else if ((lParse.Nodes)[yyvs[yyvsp - 4].Node as usize]).ntype
                                < (lParse.Nodes[yyvs[yyvsp].Node as usize]).ntype
                            {
                                yyvs[yyvsp - 4].Node = New_Unary(
                                    lParse,
                                    (lParse.Nodes[yyvs[yyvsp].Node as usize]).ntype,
                                    0,
                                    yyvs[yyvsp - 4].Node,
                                );
                            }
                            if (lParse.Nodes[yyvs[yyvsp - 2].Node as usize]).ntype
                                > (lParse.Nodes[yyvs[yyvsp].Node as usize]).ntype
                            {
                                yyvs[yyvsp].Node = New_Unary(
                                    lParse,
                                    ((lParse.Nodes)[yyvs[yyvsp - 2].Node as usize]).ntype,
                                    0,
                                    yyvs[yyvsp].Node,
                                );
                            } else if ((lParse.Nodes)[yyvs[yyvsp - 2].Node as usize]).ntype
                                < (lParse.Nodes[yyvs[yyvsp].Node as usize]).ntype
                            {
                                yyvs[yyvsp - 2].Node = New_Unary(
                                    lParse,
                                    (lParse.Nodes[yyvs[yyvsp].Node as usize]).ntype,
                                    0,
                                    yyvs[yyvsp - 2].Node,
                                );
                            }
                            yyvs[yyvsp - 2].Node = New_BinOp(
                                lParse,
                                fits_parser_yytokentype::BOOLEAN as c_int,
                                yyvs[yyvsp - 2].Node,
                                fits_parser_yytokentype::LTE as c_int,
                                yyvs[yyvsp - 4].Node,
                            );
                            yyvs[yyvsp].Node = New_BinOp(
                                lParse,
                                fits_parser_yytokentype::BOOLEAN as c_int,
                                yyvs[yyvsp - 4].Node,
                                fits_parser_yytokentype::LTE as c_int,
                                yyvs[yyvsp].Node,
                            );
                            yyval.Node = New_BinOp(
                                lParse,
                                fits_parser_yytokentype::BOOLEAN as c_int,
                                yyvs[yyvsp - 2].Node,
                                fits_parser_yytokentype::AND as c_int,
                                yyvs[yyvsp].Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        101 => {
                            /* bexpr: bexpr '?' bexpr ':' bexpr  */
                            if Test_Dims(lParse, yyvs[yyvsp - 2].Node, yyvs[yyvsp].Node) == 0 {
                                fits_parser_yyerror(
                                    lParse,
                                    cs!(c"Incompatible dimensions in '?:' arguments"),
                                );
                                current_block = 4830776507462815627;
                            } else {
                                yyval.Node = New_Func(
                                    lParse,
                                    0,
                                    IFTHENELSE_FCT,
                                    3 as c_int,
                                    yyvs[yyvsp - 2].Node,
                                    yyvs[yyvsp].Node,
                                    yyvs[yyvsp - 4].Node,
                                    0,
                                    0,
                                    0,
                                    0,
                                );
                                if yyval.Node < 0 {
                                    current_block = 4830776507462815627;
                                } else {
                                    if ((lParse.Nodes)[yyvs[yyvsp - 2].Node as usize]).value.nelem
                                        < (lParse.Nodes[yyvs[yyvsp].Node as usize]).value.nelem
                                    {
                                        Copy_Dims(lParse, yyval.Node, yyvs[yyvsp].Node);
                                    }
                                    if Test_Dims(lParse, yyvs[yyvsp - 4].Node, yyval.Node) == 0 {
                                        fits_parser_yyerror(
                                            lParse,
                                            cs!(c"Incompatible dimensions in '?:' condition"),
                                        );
                                        current_block = 4830776507462815627;
                                    } else {
                                        if (((lParse.Nodes)[yyval.Node as usize]).value).nelem
                                            < ((lParse.Nodes)[yyvs[yyvsp - 4].Node as usize])
                                                .value
                                                .nelem
                                        {
                                            Copy_Dims(lParse, yyval.Node, yyvs[yyvsp - 4].Node);
                                        }
                                        current_block = 17353983478346836848;
                                    }
                                }
                            }
                        }
                        102 => {
                            /* bexpr: BFUNCTION expr ')'  */
                            if (if c_int::from(yyvs[yyvsp - 2].astr[0])
                                < c_int::from((cs!(c"ISNULL("))[0])
                            {
                                -(1)
                            } else if c_int::from(yyvs[yyvsp - 2].astr[0])
                                > c_int::from((cs!(c"ISNULL("))[0])
                            {
                                1
                            } else {
                                strcmp((yyvs[yyvsp - 2].astr).as_mut_ptr(), c"ISNULL(".as_ptr())
                            }) == 0
                            {
                                yyval.Node = New_Func(
                                    lParse,
                                    0,
                                    ISNULL_FCT,
                                    1,
                                    yyvs[yyvsp - 1].Node,
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
                                    ((lParse.Nodes)[yyval.Node as usize]).ntype =
                                        fits_parser_yytokentype::BOOLEAN as c_int;
                                    current_block = 17353983478346836848;
                                }
                            } else {
                                fits_parser_yyerror(
                                    lParse,
                                    cs!(c"Boolean Function(expr) not supported"),
                                );
                                current_block = 4830776507462815627;
                            }
                        }
                        103 => {
                            /* bexpr: BFUNCTION bexpr ')'  */
                            if (if c_int::from(yyvs[yyvsp - 2].astr[0])
                                < c_int::from((cs!(c"ISNULL("))[0])
                            {
                                -(1)
                            } else if c_int::from(yyvs[yyvsp - 2].astr[0])
                                > c_int::from((cs!(c"ISNULL("))[0])
                            {
                                1
                            } else {
                                strcmp((yyvs[yyvsp - 2].astr).as_mut_ptr(), c"ISNULL(".as_ptr())
                            }) == 0
                            {
                                yyval.Node = New_Func(
                                    lParse,
                                    0,
                                    ISNULL_FCT,
                                    1,
                                    yyvs[yyvsp - 1].Node,
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
                                    ((lParse.Nodes)[yyval.Node as usize]).ntype =
                                        fits_parser_yytokentype::BOOLEAN as c_int;
                                    current_block = 17353983478346836848;
                                }
                            } else {
                                fits_parser_yyerror(
                                    lParse,
                                    cs!(c"Boolean Function(expr) not supported"),
                                );
                                current_block = 4830776507462815627;
                            }
                        }
                        104 => {
                            /* bexpr: BFUNCTION sexpr ')'  */
                            if (if c_int::from(yyvs[yyvsp - 2].astr[0])
                                < c_int::from((cs!(c"ISNULL("))[0])
                            {
                                -(1)
                            } else if c_int::from(yyvs[yyvsp - 2].astr[0])
                                > c_int::from((cs!(c"ISNULL("))[0])
                            {
                                1
                            } else {
                                strcmp((yyvs[yyvsp - 2].astr).as_mut_ptr(), c"ISNULL(".as_ptr())
                            }) == 0
                            {
                                yyval.Node = New_Func(
                                    lParse,
                                    fits_parser_yytokentype::BOOLEAN as c_int,
                                    ISNULL_FCT,
                                    1,
                                    yyvs[yyvsp - 1].Node,
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
                                    lParse,
                                    cs!(c"Boolean Function(expr) not supported"),
                                );
                                current_block = 4830776507462815627;
                            }
                        }
                        105 => {
                            /* bexpr: FUNCTION bexpr ',' bexpr ')'  */
                            if (if c_int::from(yyvs[yyvsp - 4].astr[0])
                                < c_int::from((cs!(c"DEFNULL("))[0])
                            {
                                -(1)
                            } else if c_int::from(yyvs[yyvsp - 4].astr[0])
                                > c_int::from((cs!(c"DEFNULL("))[0])
                            {
                                1
                            } else {
                                strcmp((yyvs[yyvsp - 4].astr).as_mut_ptr(), c"DEFNULL(".as_ptr())
                            }) == 0
                            {
                                if ((lParse.Nodes)[yyvs[yyvsp - 3].Node as usize]).value.nelem
                                    >= ((lParse.Nodes)[yyvs[yyvsp - 1].Node as usize]).value.nelem
                                    && Test_Dims(lParse, yyvs[yyvsp - 3].Node, yyvs[yyvsp - 1].Node)
                                        != 0
                                {
                                    yyval.Node = New_Func(
                                        lParse,
                                        0,
                                        DEFNULL_FCT,
                                        2,
                                        yyvs[yyvsp - 3].Node,
                                        yyvs[yyvsp - 1].Node,
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
                                        lParse,
                                        cs!(c"Dimensions of DEFNULL arguments are not compatible"),
                                    );
                                    current_block = 4830776507462815627;
                                }
                            } else {
                                fits_parser_yyerror(
                                    lParse,
                                    cs!(c"Boolean Function(expr,expr) not supported"),
                                );
                                current_block = 4830776507462815627;
                            }
                        }
                        106 => {
                            /* bexpr: BFUNCTION expr ',' expr ',' expr ')'  */
                            if ((lParse.Nodes)[yyvs[yyvsp - 5].Node as usize]).ntype
                                != fits_parser_yytokentype::DOUBLE as c_int
                            {
                                yyvs[yyvsp - 5].Node = New_Unary(
                                    lParse,
                                    fits_parser_yytokentype::DOUBLE as c_int,
                                    0,
                                    yyvs[yyvsp - 5].Node,
                                );
                            }
                            if ((lParse.Nodes)[yyvs[yyvsp - 3].Node as usize]).ntype
                                != fits_parser_yytokentype::DOUBLE as c_int
                            {
                                yyvs[yyvsp - 3].Node = New_Unary(
                                    lParse,
                                    fits_parser_yytokentype::DOUBLE as c_int,
                                    0,
                                    yyvs[yyvsp - 3].Node,
                                );
                            }
                            if ((lParse.Nodes)[yyvs[yyvsp - 1].Node as usize]).ntype
                                != fits_parser_yytokentype::DOUBLE as c_int
                            {
                                yyvs[yyvsp - 1].Node = New_Unary(
                                    lParse,
                                    fits_parser_yytokentype::DOUBLE as c_int,
                                    0,
                                    yyvs[yyvsp - 1].Node,
                                );
                            }
                            if !(Test_Dims(lParse, yyvs[yyvsp - 5].Node, yyvs[yyvsp - 3].Node) != 0
                                && Test_Dims(lParse, yyvs[yyvsp - 3].Node, yyvs[yyvsp - 1].Node)
                                    != 0)
                            {
                                fits_parser_yyerror(
                                    lParse,
                                    cs!(c"Dimensions of NEAR arguments are not compatible"),
                                );
                                current_block = 4830776507462815627;
                            } else if (if c_int::from(yyvs[yyvsp - 6].astr[0])
                                < c_int::from((cs!(c"NEAR("))[0])
                            {
                                -(1)
                            } else if c_int::from(yyvs[yyvsp - 6].astr[0])
                                > c_int::from((cs!(c"NEAR("))[0])
                            {
                                1
                            } else {
                                strcmp((yyvs[yyvsp - 6].astr).as_mut_ptr(), c"NEAR(".as_ptr())
                            }) == 0
                            {
                                yyval.Node = New_Func(
                                    lParse,
                                    fits_parser_yytokentype::BOOLEAN as c_int,
                                    NEAR_FCT,
                                    3 as c_int,
                                    yyvs[yyvsp - 5].Node,
                                    yyvs[yyvsp - 3].Node,
                                    yyvs[yyvsp - 1].Node,
                                    0,
                                    0,
                                    0,
                                    0,
                                );
                                if yyval.Node < 0 {
                                    current_block = 4830776507462815627;
                                } else {
                                    if (((lParse.Nodes)[yyval.Node as usize]).value).nelem
                                        < ((lParse.Nodes)[yyvs[yyvsp - 5].Node as usize])
                                            .value
                                            .nelem
                                    {
                                        Copy_Dims(lParse, yyval.Node, yyvs[yyvsp - 5].Node);
                                    }
                                    if ((lParse.Nodes)[yyvs[yyvsp - 5].Node as usize]).value.nelem
                                        < ((lParse.Nodes)[yyvs[yyvsp - 3].Node as usize])
                                            .value
                                            .nelem
                                    {
                                        Copy_Dims(lParse, yyval.Node, yyvs[yyvsp - 3].Node);
                                    }
                                    if ((lParse.Nodes)[yyvs[yyvsp - 3].Node as usize]).value.nelem
                                        < ((lParse.Nodes)[yyvs[yyvsp - 1].Node as usize])
                                            .value
                                            .nelem
                                    {
                                        Copy_Dims(lParse, yyval.Node, yyvs[yyvsp - 1].Node);
                                    }
                                    current_block = 17353983478346836848;
                                }
                            } else {
                                fits_parser_yyerror(lParse, cs!(c"Boolean Function not supported"));
                                current_block = 4830776507462815627;
                            }
                        }
                        107 => {
                            /* bexpr: BFUNCTION expr ',' expr ',' expr ',' expr ',' expr ')'  */
                            if ((lParse.Nodes)[yyvs[yyvsp - 9].Node as usize]).ntype
                                != fits_parser_yytokentype::DOUBLE as c_int
                            {
                                yyvs[yyvsp - 9].Node = New_Unary(
                                    lParse,
                                    fits_parser_yytokentype::DOUBLE as c_int,
                                    0,
                                    yyvs[yyvsp - 9].Node,
                                );
                            }
                            if ((lParse.Nodes)[yyvs[yyvsp - 7].Node as usize]).ntype
                                != fits_parser_yytokentype::DOUBLE as c_int
                            {
                                yyvs[yyvsp - 7].Node = New_Unary(
                                    lParse,
                                    fits_parser_yytokentype::DOUBLE as c_int,
                                    0,
                                    yyvs[yyvsp - 7].Node,
                                );
                            }
                            if ((lParse.Nodes)[yyvs[yyvsp - 5].Node as usize]).ntype
                                != fits_parser_yytokentype::DOUBLE as c_int
                            {
                                yyvs[yyvsp - 5].Node = New_Unary(
                                    lParse,
                                    fits_parser_yytokentype::DOUBLE as c_int,
                                    0,
                                    yyvs[yyvsp - 5].Node,
                                );
                            }
                            if ((lParse.Nodes)[yyvs[yyvsp - 3].Node as usize]).ntype
                                != fits_parser_yytokentype::DOUBLE as c_int
                            {
                                yyvs[yyvsp - 3].Node = New_Unary(
                                    lParse,
                                    fits_parser_yytokentype::DOUBLE as c_int,
                                    0,
                                    yyvs[yyvsp - 3].Node,
                                );
                            }
                            if ((lParse.Nodes)[yyvs[yyvsp - 1].Node as usize]).ntype
                                != fits_parser_yytokentype::DOUBLE as c_int
                            {
                                yyvs[yyvsp - 1].Node = New_Unary(
                                    lParse,
                                    fits_parser_yytokentype::DOUBLE as c_int,
                                    0,
                                    yyvs[yyvsp - 1].Node,
                                );
                            }
                            if !(Test_Dims(lParse, yyvs[yyvsp - 9].Node, yyvs[yyvsp - 7].Node) != 0
                                && Test_Dims(lParse, yyvs[yyvsp - 7].Node, yyvs[yyvsp - 5].Node)
                                    != 0
                                && Test_Dims(lParse, yyvs[yyvsp - 5].Node, yyvs[yyvsp - 3].Node)
                                    != 0
                                && Test_Dims(lParse, yyvs[yyvsp - 3].Node, yyvs[yyvsp - 1].Node)
                                    != 0)
                            {
                                fits_parser_yyerror(
                                    lParse,
                                    cs!(c"Dimensions of CIRCLE arguments are not compatible"),
                                );
                                current_block = 4830776507462815627;
                            } else if (if c_int::from(yyvs[yyvsp - 10].astr[0])
                                < c_int::from((cs!(c"CIRCLE("))[0])
                            {
                                -(1)
                            } else if c_int::from(yyvs[yyvsp - 10].astr[0])
                                > c_int::from((cs!(c"CIRCLE("))[0])
                            {
                                1
                            } else {
                                strcmp((yyvs[yyvsp - 10].astr).as_mut_ptr(), c"CIRCLE(".as_ptr())
                            }) == 0
                            {
                                yyval.Node = New_Func(
                                    lParse,
                                    fits_parser_yytokentype::BOOLEAN as c_int,
                                    CIRCLE_FCT,
                                    5 as c_int,
                                    yyvs[yyvsp - 9].Node,
                                    yyvs[yyvsp - 7].Node,
                                    yyvs[yyvsp - 5].Node,
                                    yyvs[yyvsp - 3].Node,
                                    yyvs[yyvsp - 1].Node,
                                    0,
                                    0,
                                );
                                if yyval.Node < 0 {
                                    current_block = 4830776507462815627;
                                } else {
                                    if (((lParse.Nodes)[yyval.Node as usize]).value).nelem
                                        < ((lParse.Nodes)[yyvs[yyvsp - 9].Node as usize])
                                            .value
                                            .nelem
                                    {
                                        Copy_Dims(lParse, yyval.Node, yyvs[yyvsp - 9].Node);
                                    }
                                    if ((lParse.Nodes)[yyvs[yyvsp - 9].Node as usize]).value.nelem
                                        < ((lParse.Nodes)[yyvs[yyvsp - 7].Node as usize])
                                            .value
                                            .nelem
                                    {
                                        Copy_Dims(lParse, yyval.Node, yyvs[yyvsp - 7].Node);
                                    }
                                    if ((lParse.Nodes)[yyvs[yyvsp - 7].Node as usize]).value.nelem
                                        < ((lParse.Nodes)[yyvs[yyvsp - 5].Node as usize])
                                            .value
                                            .nelem
                                    {
                                        Copy_Dims(lParse, yyval.Node, yyvs[yyvsp - 5].Node);
                                    }
                                    if ((lParse.Nodes)[yyvs[yyvsp - 5].Node as usize]).value.nelem
                                        < ((lParse.Nodes)[yyvs[yyvsp - 3].Node as usize])
                                            .value
                                            .nelem
                                    {
                                        Copy_Dims(lParse, yyval.Node, yyvs[yyvsp - 3].Node);
                                    }
                                    if ((lParse.Nodes)[yyvs[yyvsp - 3].Node as usize]).value.nelem
                                        < ((lParse.Nodes)[yyvs[yyvsp - 1].Node as usize])
                                            .value
                                            .nelem
                                    {
                                        Copy_Dims(lParse, yyval.Node, yyvs[yyvsp - 1].Node);
                                    }
                                    current_block = 17353983478346836848;
                                }
                            } else {
                                fits_parser_yyerror(lParse, cs!(c"Boolean Function not supported"));
                                current_block = 4830776507462815627;
                            }
                        }
                        108 => {
                            /* bexpr: BFUNCTION expr ',' expr ',' expr ',' expr ',' expr ',' expr ',' expr ')'  */
                            if ((lParse.Nodes)[yyvs[yyvsp - 13].Node as usize]).ntype
                                != fits_parser_yytokentype::DOUBLE as c_int
                            {
                                yyvs[yyvsp - 13].Node = New_Unary(
                                    lParse,
                                    fits_parser_yytokentype::DOUBLE as c_int,
                                    0,
                                    yyvs[yyvsp - 13].Node,
                                );
                            }
                            if ((lParse.Nodes)[yyvs[yyvsp - 11].Node as usize]).ntype
                                != fits_parser_yytokentype::DOUBLE as c_int
                            {
                                yyvs[yyvsp - 11].Node = New_Unary(
                                    lParse,
                                    fits_parser_yytokentype::DOUBLE as c_int,
                                    0,
                                    yyvs[yyvsp - 11].Node,
                                );
                            }
                            if ((lParse.Nodes)[yyvs[yyvsp - 9].Node as usize]).ntype
                                != fits_parser_yytokentype::DOUBLE as c_int
                            {
                                yyvs[yyvsp - 9].Node = New_Unary(
                                    lParse,
                                    fits_parser_yytokentype::DOUBLE as c_int,
                                    0,
                                    yyvs[yyvsp - 9].Node,
                                );
                            }
                            if ((lParse.Nodes)[yyvs[yyvsp - 7].Node as usize]).ntype
                                != fits_parser_yytokentype::DOUBLE as c_int
                            {
                                yyvs[yyvsp - 7].Node = New_Unary(
                                    lParse,
                                    fits_parser_yytokentype::DOUBLE as c_int,
                                    0,
                                    yyvs[yyvsp - 7].Node,
                                );
                            }
                            if ((lParse.Nodes)[yyvs[yyvsp - 5].Node as usize]).ntype
                                != fits_parser_yytokentype::DOUBLE as c_int
                            {
                                yyvs[yyvsp - 5].Node = New_Unary(
                                    lParse,
                                    fits_parser_yytokentype::DOUBLE as c_int,
                                    0,
                                    yyvs[yyvsp - 5].Node,
                                );
                            }
                            if ((lParse.Nodes)[yyvs[yyvsp - 3].Node as usize]).ntype
                                != fits_parser_yytokentype::DOUBLE as c_int
                            {
                                yyvs[yyvsp - 3].Node = New_Unary(
                                    lParse,
                                    fits_parser_yytokentype::DOUBLE as c_int,
                                    0,
                                    yyvs[yyvsp - 3].Node,
                                );
                            }
                            if ((lParse.Nodes)[yyvs[yyvsp - 1].Node as usize]).ntype
                                != fits_parser_yytokentype::DOUBLE as c_int
                            {
                                yyvs[yyvsp - 1].Node = New_Unary(
                                    lParse,
                                    fits_parser_yytokentype::DOUBLE as c_int,
                                    0,
                                    yyvs[yyvsp - 1].Node,
                                );
                            }
                            if !(Test_Dims(lParse, yyvs[yyvsp - 13].Node, yyvs[yyvsp - 11].Node)
                                != 0
                                && Test_Dims(lParse, yyvs[yyvsp - 11].Node, yyvs[yyvsp - 9].Node)
                                    != 0
                                && Test_Dims(lParse, yyvs[yyvsp - 9].Node, yyvs[yyvsp - 7].Node)
                                    != 0
                                && Test_Dims(lParse, yyvs[yyvsp - 7].Node, yyvs[yyvsp - 5].Node)
                                    != 0
                                && Test_Dims(lParse, yyvs[yyvsp - 5].Node, yyvs[yyvsp - 3].Node)
                                    != 0
                                && Test_Dims(lParse, yyvs[yyvsp - 3].Node, yyvs[yyvsp - 1].Node)
                                    != 0)
                            {
                                fits_parser_yyerror(lParse, cs!(c"Dimensions of BOX or ELLIPSE arguments are not compatible"));
                                current_block = 4830776507462815627;
                            } else {
                                if (if c_int::from(yyvs[yyvsp - 14].astr[0])
                                    < c_int::from((cs!(c"BOX("))[0])
                                {
                                    -(1)
                                } else if c_int::from(yyvs[yyvsp - 14].astr[0])
                                    > c_int::from((cs!(c"BOX("))[0])
                                {
                                    1
                                } else {
                                    strcmp((yyvs[yyvsp - 14].astr).as_mut_ptr(), c"BOX(".as_ptr())
                                }) == 0
                                {
                                    yyval.Node = New_Func(
                                        lParse,
                                        fits_parser_yytokentype::BOOLEAN as c_int,
                                        BOX_FCT,
                                        7 as c_int,
                                        yyvs[yyvsp - 13].Node,
                                        yyvs[yyvsp - 11].Node,
                                        yyvs[yyvsp - 9].Node,
                                        yyvs[yyvsp - 7].Node,
                                        yyvs[yyvsp - 5].Node,
                                        yyvs[yyvsp - 3].Node,
                                        yyvs[yyvsp - 1].Node,
                                    );
                                    current_block = 3023179740610631044;
                                } else if (if c_int::from(yyvs[yyvsp - 14].astr[0])
                                    < c_int::from((cs!(c"ELLIPSE"))[0])
                                {
                                    -(1)
                                } else if c_int::from(yyvs[yyvsp - 14].astr[0])
                                    > c_int::from((cs!(c"ELLIPSE"))[0])
                                {
                                    1
                                } else {
                                    strcmp(
                                        (yyvs[yyvsp - 14].astr).as_mut_ptr(),
                                        c"ELLIPSE(".as_ptr(),
                                    )
                                }) == 0
                                {
                                    yyval.Node = New_Func(
                                        lParse,
                                        fits_parser_yytokentype::BOOLEAN as c_int,
                                        ELPS_FCT,
                                        7 as c_int,
                                        yyvs[yyvsp - 13].Node,
                                        yyvs[yyvsp - 11].Node,
                                        yyvs[yyvsp - 9].Node,
                                        yyvs[yyvsp - 7].Node,
                                        yyvs[yyvsp - 5].Node,
                                        yyvs[yyvsp - 3].Node,
                                        yyvs[yyvsp - 1].Node,
                                    );
                                    current_block = 3023179740610631044;
                                } else {
                                    fits_parser_yyerror(
                                        lParse,
                                        cs!(c"SAO Image Function not supported"),
                                    );
                                    current_block = 4830776507462815627;
                                }
                                match current_block {
                                    4830776507462815627 => {}
                                    _ => {
                                        if yyval.Node < 0 {
                                            current_block = 4830776507462815627;
                                        } else {
                                            if (((lParse.Nodes)[yyval.Node as usize]).value).nelem
                                                < ((lParse.Nodes)[yyvs[yyvsp - 13].Node as usize])
                                                    .value
                                                    .nelem
                                            {
                                                Copy_Dims(
                                                    lParse,
                                                    yyval.Node,
                                                    yyvs[yyvsp - 13].Node,
                                                );
                                            }
                                            if ((lParse.Nodes)[yyvs[yyvsp - 13].Node as usize])
                                                .value
                                                .nelem
                                                < ((lParse.Nodes)[yyvs[yyvsp - 11].Node as usize])
                                                    .value
                                                    .nelem
                                            {
                                                Copy_Dims(
                                                    lParse,
                                                    yyval.Node,
                                                    yyvs[yyvsp - 11].Node,
                                                );
                                            }
                                            if ((lParse.Nodes)[yyvs[yyvsp - 11].Node as usize])
                                                .value
                                                .nelem
                                                < ((lParse.Nodes)[yyvs[yyvsp - 9].Node as usize])
                                                    .value
                                                    .nelem
                                            {
                                                Copy_Dims(lParse, yyval.Node, yyvs[yyvsp - 9].Node);
                                            }
                                            if ((lParse.Nodes)[yyvs[yyvsp - 9].Node as usize])
                                                .value
                                                .nelem
                                                < ((lParse.Nodes)[yyvs[yyvsp - 7].Node as usize])
                                                    .value
                                                    .nelem
                                            {
                                                Copy_Dims(lParse, yyval.Node, yyvs[yyvsp - 7].Node);
                                            }
                                            if ((lParse.Nodes)[yyvs[yyvsp - 7].Node as usize])
                                                .value
                                                .nelem
                                                < ((lParse.Nodes)[yyvs[yyvsp - 5].Node as usize])
                                                    .value
                                                    .nelem
                                            {
                                                Copy_Dims(lParse, yyval.Node, yyvs[yyvsp - 5].Node);
                                            }
                                            if ((lParse.Nodes)[yyvs[yyvsp - 5].Node as usize])
                                                .value
                                                .nelem
                                                < ((lParse.Nodes)[yyvs[yyvsp - 3].Node as usize])
                                                    .value
                                                    .nelem
                                            {
                                                Copy_Dims(lParse, yyval.Node, yyvs[yyvsp - 3].Node);
                                            }
                                            if ((lParse.Nodes)[yyvs[yyvsp - 3].Node as usize])
                                                .value
                                                .nelem
                                                < ((lParse.Nodes)[yyvs[yyvsp - 1].Node as usize])
                                                    .value
                                                    .nelem
                                            {
                                                Copy_Dims(lParse, yyval.Node, yyvs[yyvsp - 1].Node);
                                            }
                                            current_block = 17353983478346836848;
                                        }
                                    }
                                }
                            }
                        }
                        109 => {
                            /* bexpr: GTIFILTER ')'  */
                            /* Use defaults for all elements */
                            yyval.Node = New_GTI(
                                lParse,
                                GTIFILT_FCT,
                                (b"\0" as *const u8).cast::<c_char>() as *mut c_char,
                                -99,
                                -99,
                                c"*START*".as_ptr() as *mut c_char,
                                c"*STOP*".as_ptr() as *mut c_char,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        110 => {
                            /* bexpr: GTIFILTER STRING ')'  */
                            /* Use defaults for all except filename */
                            yyval.Node = New_GTI(
                                lParse,
                                GTIFILT_FCT,
                                (yyvs[yyvsp - 1].astr).as_mut_ptr(),
                                -99,
                                -99,
                                c"*START*".as_ptr() as *mut c_char,
                                c"*STOP*".as_ptr() as *mut c_char,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        111 => {
                            /* bexpr: GTIFILTER STRING ',' expr ')'  */
                            yyval.Node = New_GTI(
                                lParse,
                                GTIFILT_FCT,
                                (yyvs[yyvsp - 3].astr).as_mut_ptr(),
                                yyvs[yyvsp - 1].Node,
                                -99,
                                c"*START*".as_ptr() as *mut c_char,
                                c"*STOP*".as_ptr() as *mut c_char,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        112 => {
                            /* bexpr: GTIFILTER STRING ',' expr ',' STRING ',' STRING ')'  */
                            yyval.Node = New_GTI(
                                lParse,
                                GTIFILT_FCT,
                                (yyvs[yyvsp - 7].astr).as_mut_ptr(),
                                yyvs[yyvsp - 5].Node,
                                -99,
                                (yyvs[yyvsp - 3].astr).as_mut_ptr(),
                                (yyvs[yyvsp - 1].astr).as_mut_ptr(),
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        113 => {
                            /* bexpr: GTIFIND ')'  */
                            /* Use defaults for all elements */
                            yyval.Node = New_GTI(
                                lParse,
                                GTIFIND_FCT,
                                (b"\0" as *const u8).cast::<c_char>() as *mut c_char,
                                -99,
                                -99,
                                c"*START*".as_ptr() as *mut c_char,
                                c"*STOP*".as_ptr() as *mut c_char,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        114 => {
                            /* bexpr: GTIFIND STRING ')'  *//* Use defaults for all except filename */
                            yyval.Node = New_GTI(
                                lParse,
                                GTIFIND_FCT,
                                (yyvs[yyvsp - 1].astr).as_mut_ptr(),
                                -99,
                                -99,
                                c"*START*".as_ptr() as *mut c_char,
                                c"*STOP*".as_ptr() as *mut c_char,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        115 => {
                            /* bexpr: GTIFIND STRING ',' expr ')'  */
                            yyval.Node = New_GTI(
                                lParse,
                                GTIFIND_FCT,
                                (yyvs[yyvsp - 3].astr).as_mut_ptr(),
                                yyvs[yyvsp - 1].Node,
                                -99,
                                c"*START*".as_ptr() as *mut c_char,
                                c"*STOP*".as_ptr() as *mut c_char,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        116 => {
                            /* bexpr: GTIFIND STRING ',' expr ',' STRING ',' STRING ')'  */
                            yyval.Node = New_GTI(
                                lParse,
                                GTIFIND_FCT,
                                (yyvs[yyvsp - 7].astr).as_mut_ptr(),
                                yyvs[yyvsp - 5].Node,
                                -99,
                                (yyvs[yyvsp - 3].astr).as_mut_ptr(),
                                (yyvs[yyvsp - 1].astr).as_mut_ptr(),
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        117 => {
                            /* bexpr: REGFILTER STRING ')'  *//* Use defaults for all except filename */
                            let mut dummy = [0];
                            yyval.Node = New_REG(
                                lParse,
                                (yyvs[yyvsp - 1].astr).as_mut_ptr(),
                                -99,
                                -99,
                                dummy.as_mut_ptr(),
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        118 => {
                            /* bexpr: REGFILTER STRING ',' expr ',' expr ')'  */
                            let mut dummy = [0];
                            yyval.Node = New_REG(
                                lParse,
                                (yyvs[yyvsp - 5].astr).as_mut_ptr(),
                                yyvs[yyvsp - 3].Node,
                                yyvs[yyvsp - 1].Node,
                                dummy.as_mut_ptr(),
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        119 => {
                            /* bexpr: REGFILTER STRING ',' expr ',' expr ',' STRING ')'  */
                            yyval.Node = New_REG(
                                lParse,
                                (yyvs[yyvsp - 7].astr).as_mut_ptr(),
                                yyvs[yyvsp - 5].Node,
                                yyvs[yyvsp - 3].Node,
                                (yyvs[yyvsp - 1].astr).as_mut_ptr(),
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        120 => {
                            /* bexpr: bexpr '[' expr ']'  */
                            yyval.Node = New_Deref(
                                lParse,
                                yyvs[yyvsp - 3].Node,
                                1,
                                yyvs[yyvsp - 1].Node,
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
                            /* bexpr: bexpr '[' expr ',' expr ']'  */
                            yyval.Node = New_Deref(
                                lParse,
                                yyvs[yyvsp - 5].Node,
                                2,
                                yyvs[yyvsp - 3].Node,
                                yyvs[yyvsp - 1].Node,
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
                            /* bexpr: bexpr '[' expr ',' expr ',' expr ']'  */
                            yyval.Node = New_Deref(
                                lParse,
                                yyvs[yyvsp - 7].Node,
                                3 as c_int,
                                yyvs[yyvsp - 5].Node,
                                yyvs[yyvsp - 3].Node,
                                yyvs[yyvsp - 1].Node,
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
                            /* bexpr: bexpr '[' expr ',' expr ',' expr ',' expr ']'  */
                            yyval.Node = New_Deref(
                                lParse,
                                yyvs[yyvsp - 9].Node,
                                4 as c_int,
                                yyvs[yyvsp - 7].Node,
                                yyvs[yyvsp - 5].Node,
                                yyvs[yyvsp - 3].Node,
                                yyvs[yyvsp - 1].Node,
                                0,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        124 => {
                            /* bexpr: bexpr '[' expr ',' expr ',' expr ',' expr ',' expr ']'  */
                            yyval.Node = New_Deref(
                                lParse,
                                yyvs[yyvsp - 11].Node,
                                5 as c_int,
                                yyvs[yyvsp - 9].Node,
                                yyvs[yyvsp - 7].Node,
                                yyvs[yyvsp - 5].Node,
                                yyvs[yyvsp - 3].Node,
                                yyvs[yyvsp - 1].Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        125 => {
                            /* bexpr: NOT bexpr  */
                            yyval.Node = New_Unary(
                                lParse,
                                fits_parser_yytokentype::BOOLEAN as c_int,
                                fits_parser_yytokentype::NOT as c_int,
                                yyvs[yyvsp].Node,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        126 => {
                            /* bexpr: '(' bexpr ')'  */
                            yyval.Node = yyvs[yyvsp - 1].Node;
                            current_block = 17353983478346836848;
                        }
                        127 => {
                            /* sexpr: STRING  */
                            yyval.Node = New_Const(
                                lParse,
                                fits_parser_yytokentype::STRING as c_int,
                                (yyvs[yyvsp].astr).as_mut_ptr().cast::<c_void>(),
                                (strlen((yyvs[yyvsp].astr).as_mut_ptr()))
                                    .wrapping_add((1).try_into().unwrap())
                                    as c_long,
                            );
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                (((lParse.Nodes)[yyval.Node as usize]).value).nelem =
                                    strlen((yyvs[yyvsp].astr).as_mut_ptr()) as c_long;
                                current_block = 17353983478346836848;
                            }
                        }
                        128 => {
                            /* sexpr: SCOLUMN  */
                            yyval.Node = New_Column(lParse, yyvs[yyvsp].lng as c_int);
                            if yyval.Node < 0 {
                                current_block = 4830776507462815627;
                            } else {
                                current_block = 17353983478346836848;
                            }
                        }
                        129 => {
                            /* sexpr: SCOLUMN '{' expr '}'  */
                            if ((lParse.Nodes)[yyvs[yyvsp - 1].Node as usize]).ntype
                                != fits_parser_yytokentype::LONG as c_int
                                || ((lParse.Nodes)[yyvs[yyvsp - 1].Node as usize]).operation
                                    != CONST_OP
                            {
                                fits_parser_yyerror(
                                    lParse,
                                    cs!(c"Offset argument must be a constant integer"),
                                );
                                current_block = 4830776507462815627;
                            } else {
                                yyval.Node = New_Offset(
                                    lParse,
                                    yyvs[yyvsp - 3].lng as c_int,
                                    yyvs[yyvsp - 1].Node,
                                );
                                if yyval.Node < 0 {
                                    current_block = 4830776507462815627;
                                } else {
                                    current_block = 17353983478346836848;
                                }
                            }
                        }
                        130 => {
                            /* sexpr: SNULLREF  */
                            yyval.Node = New_Func(
                                lParse,
                                fits_parser_yytokentype::STRING as c_int,
                                NULL_FCT,
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
                            /* sexpr: '(' sexpr ')'  */
                            yyval.Node = yyvs[yyvsp - 1].Node;
                            current_block = 17353983478346836848;
                        }
                        132 => {
                            /* sexpr: sexpr '+' sexpr  */
                            if (lParse.Nodes[yyvs[yyvsp - 2].Node as usize]).value.nelem
                                + (lParse.Nodes[yyvs[yyvsp].Node as usize]).value.nelem
                                >= 256 as c_long
                            {
                                fits_parser_yyerror(
                                    lParse,
                                    cs!(c"Combined string size exceeds 255 characters"),
                                );
                                current_block = 4830776507462815627;
                            } else {
                                yyval.Node = New_BinOp(
                                    lParse,
                                    fits_parser_yytokentype::STRING as c_int,
                                    yyvs[yyvsp - 2].Node,
                                    '+' as i32,
                                    yyvs[yyvsp].Node,
                                );
                                if yyval.Node < 0 {
                                    current_block = 4830776507462815627;
                                } else {
                                    (((lParse.Nodes)[yyval.Node as usize]).value).nelem =
                                        ((lParse.Nodes)[yyvs[yyvsp - 2].Node as usize]).value.nelem
                                            + ((lParse.Nodes)[yyvs[yyvsp].Node as usize])
                                                .value
                                                .nelem;
                                    current_block = 17353983478346836848;
                                }
                            }
                        }
                        133 => {
                            /* sexpr: bexpr '?' sexpr ':' sexpr  */
                            let mut outSize: c_int = 0;
                            if ((lParse.Nodes)[yyvs[yyvsp - 4].Node as usize]).value.nelem != 1 {
                                fits_parser_yyerror(
                                    lParse,
                                    cs!(c"Cannot have a vector string column"),
                                );
                                current_block = 4830776507462815627; // goto yyerrorlab
                            } else {
                                /* Since the output can be calculated now, as a constant
                                scalar, we must precalculate the output size, in
                                order to avoid an overflow. */

                                outSize = (((lParse.Nodes)[yyvs[yyvsp - 2].Node as usize]).value)
                                    .nelem as c_int;
                                if (lParse.Nodes[yyvs[yyvsp].Node as usize]).value.nelem
                                    > c_long::from(outSize)
                                {
                                    outSize =
                                        ((lParse.Nodes)[yyvs[yyvsp].Node as usize]).value.nelem
                                            as c_int;
                                }
                                yyval.Node = New_FuncSize(
                                    lParse,
                                    0,
                                    IFTHENELSE_FCT,
                                    3 as c_int,
                                    yyvs[yyvsp - 2].Node,
                                    yyvs[yyvsp].Node,
                                    yyvs[yyvsp - 4].Node,
                                    0,
                                    0,
                                    0,
                                    0,
                                    outSize,
                                );
                                if yyval.Node < 0 {
                                    current_block = 4830776507462815627;
                                } else {
                                    if ((lParse.Nodes)[yyvs[yyvsp - 2].Node as usize]).value.nelem
                                        < (lParse.Nodes[yyvs[yyvsp].Node as usize]).value.nelem
                                    {
                                        Copy_Dims(lParse, yyval.Node, yyvs[yyvsp].Node);
                                    }
                                    current_block = 17353983478346836848;
                                }
                            }
                        }
                        134 => {
                            /* sexpr: FUNCTION sexpr ',' sexpr ')'  */
                            if (if c_int::from(yyvs[yyvsp - 4].astr[0])
                                < c_int::from((cs!(c"DEFNULL("))[0])
                            {
                                -(1)
                            } else if c_int::from(yyvs[yyvsp - 4].astr[0])
                                > c_int::from((cs!(c"DEFNULL("))[0])
                            {
                                1
                            } else {
                                strcmp((yyvs[yyvsp - 4].astr).as_mut_ptr(), c"DEFNULL(".as_ptr())
                            }) == 0
                            {
                                /* Since the output can be calculated now, as a constant
                                scalar, we must precalculate the output size, in
                                order to avoid an overflow. */

                                let mut outSize_0: c_int = 0;
                                outSize_0 = ((lParse.Nodes)[yyvs[yyvsp - 3].Node as usize])
                                    .value
                                    .nelem as c_int;
                                if ((lParse.Nodes)[yyvs[yyvsp - 1].Node as usize]).value.nelem
                                    > c_long::from(outSize_0)
                                {
                                    outSize_0 =
                                        ((lParse.Nodes)[yyvs[yyvsp - 1].Node as usize]).value.nelem
                                            as c_int;
                                }
                                yyval.Node = New_FuncSize(
                                    lParse,
                                    0,
                                    DEFNULL_FCT,
                                    2,
                                    yyvs[yyvsp - 3].Node,
                                    yyvs[yyvsp - 1].Node,
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
                                    if ((lParse.Nodes)[yyvs[yyvsp - 1].Node as usize]).value.nelem
                                        > ((lParse.Nodes)[yyvs[yyvsp - 3].Node as usize])
                                            .value
                                            .nelem
                                    {
                                        (((lParse.Nodes)[yyval.Node as usize]).value).nelem =
                                            ((lParse.Nodes)[yyvs[yyvsp - 1].Node as usize])
                                                .value
                                                .nelem;
                                    }
                                    current_block = 17353983478346836848;
                                }
                            } else {
                                fits_parser_yyerror(
                                    lParse,
                                    cs!(c"Function(string,string) not supported"),
                                );
                                current_block = 4830776507462815627;
                            }
                        }
                        135 => {
                            /* sexpr: FUNCTION sexpr ',' expr ',' expr ')'  */
                            if (if c_int::from(yyvs[yyvsp - 6].astr[0])
                                < c_int::from((cs!(c"STRMID("))[0])
                            {
                                -(1)
                            } else if c_int::from(yyvs[yyvsp - 6].astr[0])
                                > c_int::from((cs!(c"STRMID("))[0])
                            {
                                1
                            } else {
                                strcmp((yyvs[yyvsp - 6].astr).as_mut_ptr(), c"STRMID(".as_ptr())
                            }) == 0
                            {
                                let mut len: c_int = 0;
                                if ((lParse.Nodes)[yyvs[yyvsp - 3].Node as usize]).ntype
                                    != fits_parser_yytokentype::LONG as c_int
                                    || ((lParse.Nodes)[yyvs[yyvsp - 3].Node as usize]).value.nelem
                                        != 1
                                    || ((lParse.Nodes)[yyvs[yyvsp - 1].Node as usize]).ntype
                                        != fits_parser_yytokentype::LONG as c_int
                                    || ((lParse.Nodes)[yyvs[yyvsp - 1].Node as usize]).value.nelem
                                        != 1
                                {
                                    fits_parser_yyerror(lParse, cs!(c"When using STRMID(S,P,N), P and N must be integers (and not vector columns)"));
                                    current_block = 4830776507462815627;
                                } else {
                                    if ((lParse.Nodes)[yyvs[yyvsp - 1].Node as usize]).operation
                                        == CONST_OP
                                    {
                                        /* Constant value: use that directly */
                                        len = ((lParse.Nodes)[yyvs[yyvsp - 1].Node as usize])
                                            .value
                                            .data
                                            .lng
                                            as c_int;
                                    } else {
                                        /* Variable value: use the maximum possible (from $2) */
                                        len = ((lParse.Nodes)[yyvs[yyvsp - 5].Node as usize])
                                            .value
                                            .nelem
                                            as c_int;
                                    }
                                    if len <= 0 || len >= 256 as c_int {
                                        fits_parser_yyerror(
                                            lParse,
                                            cs!(c"STRMID(S,P,N), N must be 1-255"),
                                        );
                                        current_block = 4830776507462815627;
                                    } else {
                                        yyval.Node = New_FuncSize(
                                            lParse,
                                            0,
                                            STRMID_FCT,
                                            3 as c_int,
                                            yyvs[yyvsp - 5].Node,
                                            yyvs[yyvsp - 3].Node,
                                            yyvs[yyvsp - 1].Node,
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
                                    lParse,
                                    cs!(c"Function(string,expr,expr) not supported"),
                                );
                                current_block = 4830776507462815627;
                            }
                        }
                        _ => {
                            current_block = 17353983478346836848;
                        }
                    }

                    /* User semantic actions sometimes alter yychar, and that requires
                    that yytoken be updated with the new translation.  We take the
                    approach of translating immediately before every use of yytoken.
                    One alternative is translating here after every semantic action,
                    but that translation would be missed if the semantic action invokes
                    YYABORT, YYACCEPT, or YYERROR immediately after altering yychar or
                    if it invokes YYBACKUP.  In the case of YYABORT or YYACCEPT, an
                    incorrect destructor might then be invoked immediately.  In the
                    case of YYERROR or YYBACKUP, subsequent parser actions might lead
                    to an incorrect destructor call or verbose syntax error message
                    before the lookahead is translated.  */

                    YY_SYMBOL_PRINT!(
                        "-> $$ =",
                        yysymbol_kind_t::from(YYR1[yyn as usize]),
                        &yyval,
                        scanner,
                        lParse
                    );

                    match current_block {
                        4830776507462815627 => {
                            fits_parser_yynerrs += 1;
                            yyvsp -= yylen as usize;
                            yyssp -= yylen as usize;
                            yylen = 0;
                            yystate = yy_state_fast_t::from(yyss[yyssp]);
                            current_block = 1774893048582444437;
                        }
                        _ => {
                            yyvsp -= yylen as usize;
                            yyssp -= yylen as usize;
                            yylen = 0;
                            yyvsp += 1;
                            yyvs[yyvsp] = yyval;

                            /* Now 'shift' the result of the reduction.  Determine what state
                            that goes to, based on the state we popped back to and the rule
                            number reduced by.  */

                            let yylhs: c_int = c_int::from(YYR1[yyn as usize]) - YYNTOKENS as c_int;
                            let yyi: c_int =
                                c_int::from(YYPGOTO[yylhs as usize]) + c_int::from(yyss[yyssp]);
                            yystate = if 0 <= yyi
                                && yyi <= YYLAST as c_int
                                && c_int::from(YYCHECK[yyi as usize]) == c_int::from(yyss[yyssp])
                            {
                                c_int::from(YYTABLE[yyi as usize])
                            } else {
                                c_int::from(YYDEFGOTO[yylhs as usize])
                            };
                            current_block = 7872030484262409139; // goto yynewstate;
                        }
                    }
                }

                if current_block == 1774893048582444437 {
                    yyerrstatus = 3 as c_int;
                    loop {
                        yyn = c_int::from(YYPACT[yystate as usize]);
                        if yyn != -41 {
                            yyn += YYSYMBOL_YYERROR as c_int;
                            if 0 <= yyn
                                && yyn <= 1776 as c_int
                                && c_int::from(YYCHECK[yyn as usize]) == YYSYMBOL_YYERROR as c_int
                            {
                                yyn = c_int::from(YYTABLE[yyn as usize]);
                                if (0 as c_int) < yyn {
                                    break;
                                }
                            }
                        }
                        if yyssp == 0 {
                            current_block = 3964311021479492664;
                            break 's_54;
                        }
                        yydestruct(
                            c"Error: popping".as_ptr(),
                            yysymbol_kind_t::from(YYSTOS[yystate as usize]),
                            &mut yyvs[yyvsp],
                            scanner,
                            lParse,
                        );
                        yyvsp -= 1;
                        yyssp -= 1;
                        yystate = yy_state_fast_t::from(yyss[yyssp]);
                    }
                    yyvsp += 1;
                    yyvs[yyvsp] = yylval;
                    yystate = yyn;
                }

                /*------------------------------------------------------------.
                | yynewstate -- push a new state, which is found in yystate.  |
                `------------------------------------------------------------*/
                yyssp += 1;
                /* In all cases, when you get here, the value and location stacks
                have just been pushed.  So pushing a state here evens the stacks.  */
            }
        }
        match current_block {
            11794367917084412820 => {
                fits_parser_yyerror(lParse, cs!(c"memory exhausted"));
                yyresult = 2;
            }
            3964311021479492664 => {
                yyresult = 1;
            }
            _ => {}
        }
        if yychar != fits_parser_yytokentype::FITS_PARSER_YYEMPTY as c_int {
            yytoken = (if 0 <= yychar && yychar <= 292 as c_int {
                yysymbol_kind_t::from(YYTRANSLATE[yychar as usize]) as c_int
            } else {
                YYSYMBOL_YYUNDEF as c_int
            }) as yysymbol_kind_t;
            yydestruct(
                c"Cleanup: discarding lookahead".as_ptr(),
                yytoken,
                &mut yylval,
                scanner,
                lParse,
            );
        }
        yyvsp -= yylen as usize;
        yyssp -= yylen as usize;
        while yyssp != 0 {
            yydestruct(
                c"Cleanup: popping".as_ptr(),
                yysymbol_kind_t::from(YYSTOS[yyss[yyssp] as usize]),
                &mut yyvs[yyvsp],
                scanner,
                lParse,
            );
            yyvsp -= 1;
            yyssp -= 1;
        }

        let yyss_tmp_ptr = yyss.as_ptr();

        if yyss_tmp_ptr != yyssa.as_ptr() {
            free(yyss_tmp_ptr.cast::<c_void>() as *mut c_void);
        }
        yyresult
    }
}

fn New_Deref(
    lParse: &mut ParseData,
    Var: c_int,
    nDim: c_int,
    Dim1: c_int,
    Dim2: c_int,
    Dim3: c_int,
    Dim4: c_int,
    Dim5: c_int,
) -> c_int {
    let mut n: c_int = 0;
    let mut idx: c_int = 0;
    let mut constant: c_int = 0;
    let mut elem: c_long = 0;
    let this_node_idx: usize;
    let mut theDim: [usize; 5] = [0; 5];

    if Var < 0 || Dim1 < 0 || Dim2 < 0 || Dim3 < 0 || Dim4 < 0 || Dim5 < 0 {
        return -(1);
    }

    let Var: usize = Var as usize;
    let Dim1: usize = Dim1 as usize;
    let Dim2: usize = Dim2 as usize;
    let Dim3: usize = Dim3 as usize;
    let Dim4: usize = Dim4 as usize;
    let Dim5: usize = Dim5 as usize;

    let theVar = Var;

    if (lParse.Nodes[theVar]).operation == CONST_OP || (lParse.Nodes[theVar]).value.nelem == 1 {
        fits_parser_yyerror(lParse, cs!(c"Cannot index a scalar value"));
        return -(1);
    }

    n = Alloc_Node(lParse);

    if n >= 0 {
        this_node_idx = n as usize;
        (lParse.Nodes[this_node_idx]).nSubNodes = nDim + 1;
        (lParse.Nodes[this_node_idx]).SubNodes[0] = Var;

        let theVar = (lParse.Nodes[this_node_idx]).SubNodes[0];

        (lParse.Nodes[this_node_idx]).SubNodes[1] = Dim1;
        (lParse.Nodes[this_node_idx]).SubNodes[2] = Dim2;
        (lParse.Nodes[this_node_idx]).SubNodes[3] = Dim3;
        (lParse.Nodes[this_node_idx]).SubNodes[4] = Dim4;
        (lParse.Nodes[this_node_idx]).SubNodes[5] = Dim5;

        theDim[0] = Dim1;
        theDim[1] = Dim2;
        theDim[2] = Dim3;
        theDim[3] = Dim4;
        theDim[4] = Dim5;

        constant = c_int::from((lParse.Nodes[theVar]).operation == CONST_OP);
        idx = 0;
        while idx < nDim {
            constant = c_int::from(
                constant != 0 && (lParse.Nodes[theDim[idx as usize]]).operation == CONST_OP,
            );
            idx += 1;
        }
        idx = 0;
        while idx < nDim {
            if (lParse.Nodes[theDim[idx as usize]]).value.nelem > 1 {
                Free_Last_Node(lParse);
                fits_parser_yyerror(lParse, cs!(c"Cannot use an array as an index value"));
                return -(1);
            } else if (lParse.Nodes[theDim[idx as usize]]).ntype
                != fits_parser_yytokentype::LONG as c_int
            {
                Free_Last_Node(lParse);
                fits_parser_yyerror(lParse, cs!(c"Index value must be an integer type"));
                return -(1);
            }
            idx += 1;
        }

        (lParse.Nodes[this_node_idx]).operation = '[' as i32;
        (lParse.Nodes[this_node_idx]).DoOp = Some(Do_Deref);
        (lParse.Nodes[this_node_idx]).ntype = (lParse.Nodes[theVar]).ntype;

        if (lParse.Nodes[theVar]).value.naxis == nDim {
            /* All dimensions specified */
            (lParse.Nodes[this_node_idx]).value.nelem = 1;
            (lParse.Nodes[this_node_idx]).value.naxis = 1;
            (lParse.Nodes[this_node_idx]).value.naxes[0] = 1;
        } else if nDim == 1 {
            /* Dereference only one dimension */

            elem = 1;
            (lParse.Nodes[this_node_idx]).value.naxis = (lParse.Nodes[theVar]).value.naxis - 1;
            idx = 0;
            while idx < (lParse.Nodes[this_node_idx]).value.naxis {
                (lParse.Nodes[this_node_idx]).value.naxes[idx as usize] =
                    (lParse.Nodes[theVar]).value.naxes[idx as usize];
                elem *= (lParse.Nodes[this_node_idx]).value.naxes[idx as usize];
                idx += 1;
            }
            (lParse.Nodes[this_node_idx]).value.nelem = elem;
        } else {
            Free_Last_Node(lParse);
            fits_parser_yyerror(
                lParse,
                cs!(c"Must specify just one or all indices for vector"),
            );
            return -(1);
        }
        if constant != 0 {
            ((lParse.Nodes[this_node_idx]).DoOp).expect("non-null function pointer")(
                lParse,
                this_node_idx,
            );
        }
    }
    n
}

fn New_GTI(
    lParse: &mut ParseData,
    Op: funcOp,
    mut fname: *mut c_char,
    mut Node1: c_int,
    mut Node2: c_int,
    start: *mut c_char,
    stop: *mut c_char,
) -> c_int {
    unsafe {
        let mut fptr: *mut fitsfile = core::ptr::null_mut();
        let this_node_idx: usize;
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
        let mut dt: c_double = 0.0;
        let mut timeSpan: c_double = 0.0;
        let mut xcol: [c_char; 20] = [0; 20];
        let mut xexpr: [c_char; 20] = [0; 20];
        let mut colVal: FITS_PARSER_YYSTYPE = FITS_PARSER_YYSTYPE { Node: 0 };

        if (Op as c_uint == GTIFILT_FCT as c_int as c_uint
            || Op as c_uint == GTIFIND_FCT as c_int as c_uint)
            && Node1 == -99
        {
            type_0 = fits_parser_yyGetVariable(lParse, cs!(c"TIME"), &mut colVal);
            if type_0 == fits_parser_yytokentype::COLUMN as c_int {
                Node1 = New_Column(lParse, colVal.lng as c_int);
            } else {
                fits_parser_yyerror(
                    lParse,
                    cs!(c"Could not build TIME column for GTIFILTER/GTIFIND"),
                );
                return -(1);
            }
        }

        if Op as c_uint == GTIOVER_FCT as c_int as c_uint {
            if Node1 == -99 || Node2 == -99 {
                fits_parser_yyerror(
                    lParse,
                    cs!(c"startExpr and stopExpr values must be defined for GTIOVERLAP"),
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

        fptr = lParse.def_fptr;
        let fptr = fptr.as_mut().unwrap();

        ffghdn_safe(fptr, &mut evthdu);
        tstat = 0;
        if ffgkyd_safe(fptr, cs!(c"TIMEZERO"), &mut timeZeroI[0], None, &mut tstat) != 0 {
            tstat = 0;
            if ffgkyd_safe(fptr, cs!(c"TIMEZEROI"), &mut timeZeroI[0], None, &mut tstat) != 0 {
                timeZeroF[0] = 0.0;
                timeZeroI[0] = timeZeroF[0];
            } else if ffgkyd_safe(fptr, cs!(c"TIMEZEROF"), &mut timeZeroF[0], None, &mut tstat) != 0
            {
                timeZeroF[0] = 0.0;
            }
        } else {
            timeZeroF[0] = 0.0;
        }

        /*  Resolve filename parameter  */

        match *fname.offset(0) as u8 {
            0 => {
                samefile = 1;
                hdunum = 1;
            }
            b'[' => {
                samefile = 1;
                i = 1;
                while c_int::from(*fname.offset(i as isize)) != 0_i32
                    && c_int::from(*fname.offset(i as isize)) != ']' as i32
                {
                    i += 1;
                }
                if *fname.offset(i as isize) != 0 {
                    *fname.offset(i as isize) = 0;
                    fname = fname.offset(1);
                    let fname_str = CStr::from_ptr(fname).to_bytes_with_nul();
                    ffexts_safe(
                        core::slice::from_raw_parts(
                            fname_str.as_ptr().cast::<c_char>(),
                            fname_str.len(),
                        ),
                        &mut hdunum,
                        &mut extname,
                        &mut extvers,
                        &mut movetotype,
                        &mut xcol,
                        &mut xexpr,
                        &mut lParse.status,
                    );
                    if *extname.as_mut_ptr() != 0 {
                        ffmnhd_safe(fptr, movetotype, &extname, extvers, &mut lParse.status);
                        ffghdn_safe(fptr, &mut hdunum);
                    } else if hdunum != 0 {
                        hdunum += 1;
                        ffmahd_safe(fptr, hdunum, Some(&mut hdutype), &mut lParse.status);
                    } else if lParse.status == 0 {
                        fits_parser_yyerror(
                            lParse,
                            cs!(c"Cannot use primary array for GTI filter"),
                        );
                        return -(1);
                    }
                } else {
                    fits_parser_yyerror(lParse, cs!(c"File extension specifier lacks closing ']'"));
                    return -(1);
                }
            }
            b'+' => {
                samefile = 1;
                hdunum =
                    atoi(core::ffi::CStr::from_ptr(fname).to_str().unwrap_or("0")).unwrap_or(0) + 1;
                if hdunum > 1 {
                    ffmahd_safe(fptr, hdunum, Some(&mut hdutype), &mut lParse.status);
                } else {
                    fits_parser_yyerror(
                        lParse,
                        cs!(c"Cannot use primary array for GTI filter / GTIFIND"),
                    );
                    return -(1);
                }
            }
            _ => {
                samefile = 0;
                let mut fptr_tmp: Option<Box<fitsfile>> = None;
                let fname_str = CStr::from_ptr(fname).to_bytes_with_nul();
                if ffopen_safe(
                    &mut fptr_tmp,
                    core::slice::from_raw_parts(
                        fname_str.as_ptr().cast::<c_char>(),
                        fname_str.len(),
                    ),
                    0,
                    &mut lParse.status,
                ) == 0
                {
                    ffghdn_safe(fptr, &mut hdunum);
                }
            }
        }

        if lParse.status != 0 {
            return -(1);
        }

        /*  If at primary, search for GTI extension  */

        if hdunum == 1 {
            loop {
                hdunum += 1;
                if ffmahd_safe(fptr, hdunum, Some(&mut hdutype), &mut lParse.status) != 0 {
                    break;
                }
                if hdutype == 0 {
                    continue;
                }
                tstat = 0;
                if ffgkys_safe(fptr, cs!(c"EXTNAME"), &mut extname, None, &mut tstat) != 0 {
                    continue;
                }
                ffupch_safe(&mut extname);
                if !(strstr(extname.as_mut_ptr(), c"GTI".as_ptr())).is_null() {
                    break;
                }
            }
            if lParse.status != 0 {
                if lParse.status == 107 as c_int {
                    fits_parser_yyerror(lParse, cs!(c"GTI extension not found in this file"));
                }
                return -(1);
            }
        }

        /*  Locate START/STOP Columns  */
        let start_str = CStr::from_ptr(start).to_bytes_with_nul();
        let stop_str = CStr::from_ptr(stop).to_bytes_with_nul();
        ffgcno_safe(
            fptr,
            0,
            core::slice::from_raw_parts(start_str.as_ptr().cast::<c_char>(), start_str.len()),
            &mut startCol,
            &mut lParse.status,
        );
        ffgcno_safe(
            fptr,
            0,
            core::slice::from_raw_parts(stop_str.as_ptr().cast::<c_char>(), stop_str.len()),
            &mut stopCol,
            &mut lParse.status,
        );

        if lParse.status != 0 {
            return -(1);
        }

        /*  Look for TIMEZERO keywords in GTI extension  */
        tstat = 0;
        if ffgkyd_safe(fptr, cs!(c"TIMEZERO"), &mut timeZeroI[1], None, &mut tstat) != 0 {
            tstat = 0;
            if ffgkyd_safe(fptr, cs!(c"TIMEZEROI"), &mut timeZeroI[1], None, &mut tstat) != 0 {
                timeZeroF[1] = 0.0;
                timeZeroI[1] = timeZeroF[1];
            } else if ffgkyd_safe(fptr, cs!(c"TIMEZERF"), &mut timeZeroF[1], None, &mut tstat) != 0
            {
                timeZeroF[1] = 0.0;
            }
        } else {
            timeZeroF[1] = 0.0;
        }
        n = Alloc_Node(lParse);
        if n >= 0 {
            this_node_idx = n as usize;
            (lParse.Nodes[this_node_idx]).SubNodes[1] = Node1.try_into().unwrap();
            (lParse.Nodes[this_node_idx]).operation = Op as c_int;
            if Op as c_uint == GTIFILT_FCT as c_int as c_uint {
                (lParse.Nodes[this_node_idx]).nSubNodes = 2;
                (lParse.Nodes[this_node_idx]).DoOp = Some(Do_GTI);
                (lParse.Nodes[this_node_idx]).ntype = fits_parser_yytokentype::BOOLEAN as c_int;
            } else if Op as c_uint == GTIFIND_FCT as c_int as c_uint {
                (lParse.Nodes[this_node_idx]).nSubNodes = 2;
                (lParse.Nodes[this_node_idx]).DoOp = Some(Do_GTI);
                (lParse.Nodes[this_node_idx]).ntype = fits_parser_yytokentype::LONG as c_int;
            } else {
                (lParse.Nodes[this_node_idx]).nSubNodes = 3 as c_int;
                (lParse.Nodes[this_node_idx]).DoOp = Some(Do_GTI_Over);
                (lParse.Nodes[this_node_idx]).ntype = fits_parser_yytokentype::DOUBLE as c_int;
            }

            let that1_idx = Node1 as usize;
            (lParse.Nodes[this_node_idx]).value.nelem = (lParse.Nodes[that1_idx]).value.nelem;
            (lParse.Nodes[this_node_idx]).value.naxis = (lParse.Nodes[that1_idx]).value.naxis;
            i = 0;
            while i < (lParse.Nodes[that1_idx]).value.naxis {
                (lParse.Nodes[this_node_idx]).value.naxes[i as usize] =
                    (lParse.Nodes[that1_idx]).value.naxes[i as usize];
                i += 1;
            }
            if Op as c_uint == GTIOVER_FCT as c_int as c_uint {
                (lParse.Nodes[this_node_idx]).SubNodes[2] = Node2.try_into().unwrap();
                let that2_idx = Node2 as usize;
                if (lParse.Nodes[that1_idx]).value.nelem != (lParse.Nodes[that2_idx]).value.nelem {
                    fits_parser_yyerror(
                        lParse,
                        cs!(c"Dimensions of TIME and TIME_STOP must match for GTIOVERLAP"),
                    );
                    return -(1);
                }
            }

            /* Init START/STOP node to be treated as a "constant" */

            (lParse.Nodes[this_node_idx]).SubNodes[0] = Node0.try_into().unwrap();
            let that0_idx = Node0 as usize;
            (lParse.Nodes[that0_idx]).operation = CONST_OP;
            (lParse.Nodes[that0_idx]).DoOp = None;
            (lParse.Nodes[that0_idx]).value.data.ptr = core::ptr::null_mut::<c_void>();

            /*  Read in START/STOP times  */

            if ffgkyj_safe(fptr, cs!(c"NAXIS2"), &mut nrows, None, &mut lParse.status) != 0 {
                return -(1);
            }

            (lParse.Nodes[that0_idx]).value.nelem = nrows;
            if nrows != 0 {
                let mut startptr: *mut c_double = core::ptr::null_mut::<c_double>();
                let mut stopptr: *mut c_double = core::ptr::null_mut::<c_double>();

                /* We are allocating storage for both START and STOP with one pointer
                and stop is stored at dblptr+nrows, we will use aliases below to
                make this easier to read */
                (lParse.Nodes[that0_idx]).value.data.dblptr = malloc(
                    ((2 as c_long * nrows) as c_ulong)
                        .wrapping_mul(::core::mem::size_of::<c_double>() as c_ulong)
                        .try_into()
                        .unwrap(),
                )
                .cast::<c_double>();
                if ((lParse.Nodes[that0_idx]).value.data.dblptr).is_null() {
                    lParse.status = 113 as c_int;
                    return -(1);
                }
                startptr = (lParse.Nodes[that0_idx]).value.data.dblptr;
                stopptr = ((lParse.Nodes[that0_idx]).value.data.dblptr).offset(nrows as isize);
                ffgcvd_safe(
                    fptr,
                    startCol,
                    1,
                    1,
                    nrows as LONGLONG,
                    0.0,
                    core::slice::from_raw_parts_mut(startptr, nrows as usize),
                    Some(&mut i),
                    &mut lParse.status,
                );
                ffgcvd_safe(
                    fptr,
                    stopCol,
                    1,
                    1,
                    nrows as LONGLONG,
                    0.0,
                    core::slice::from_raw_parts_mut(stopptr, nrows as usize),
                    Some(&mut i),
                    &mut lParse.status,
                );
                if lParse.status != 0 {
                    free((lParse.Nodes[that0_idx]).value.data.dblptr.cast::<c_void>());
                    return -(1);
                }

                /*  Test for fully time-ordered GTI... both START && STOP  */
                (lParse.Nodes[that0_idx]).ntype = 1; /*  Assume yes  */
                i = nrows as c_int;
                loop {
                    /* the following are failure conditions for GTI ordering */

                    /* C: while( --j ) - the body never runs with j == 0, so
                    stopptr[j-1] is always in bounds */
                    i -= 1;
                    if i == 0 {
                        break;
                    }

                    if !(*startptr.offset(i as isize) > *stopptr.offset(i as isize) /* START{j} > STOP{j} */
                        || *startptr.offset(i as isize) < *stopptr.offset((i - 1) as isize))
                    /* START{j} < STOP{j-1} */
                    {
                        continue;
                    }
                    (lParse.Nodes[that0_idx]).ntype = 0;
                    break;
                }

                /* GTIOVERLAP() requires ordered GTI */
                if (lParse.Nodes[that0_idx]).ntype != 1
                    && Op as c_uint == GTIOVER_FCT as c_int as c_uint
                {
                    let mut errmsg: [c_char; 120] = [0; 120];
                    int_snprintf!(
                        &mut errmsg,
                        120,
                        "Input GTI must be time-ordered for GTIOVERLAP (row {})",
                        i + 1,
                    );
                    fits_parser_yyerror(lParse, &errmsg);
                    return -(1);
                }

                /*  Handle TIMEZERO offset, if any  */
                dt = timeZeroI[1] - timeZeroI[0] + (timeZeroF[1] - timeZeroF[0]);
                timeSpan = *stopptr.offset((nrows - 1) as isize) - *startptr.offset(0);
                if timeSpan == 0.0 {
                    timeSpan = 1.0;
                }
                if fabs(dt / timeSpan) > 1e-12f64 {
                    i = 0;
                    while c_long::from(i) < nrows {
                        *startptr.offset(i as isize) += dt;
                        *stopptr.offset(i as isize) += dt;
                        i += 1;
                    }
                }
            }

            /* If Node1 is constant (gtifilt_fct) or
            Node1 and Node2 are constant (gtiover_fct), then evaluate now */
            if ((lParse.Nodes)[Node1 as usize]).operation == CONST_OP
                && (Op as c_uint == GTIFILT_FCT as c_int as c_uint
                    || ((lParse.Nodes)[Node2 as usize]).operation == CONST_OP)
            {
                ((lParse.Nodes[this_node_idx]).DoOp).expect("non-null function pointer")(
                    lParse,
                    this_node_idx,
                );
            }
        }

        if samefile != 0 {
            ffmahd_safe(fptr, evthdu, Some(&mut hdutype), &mut lParse.status);
        } else {
            ffclos_safe(Box::from_raw(fptr), &mut lParse.status);
        }
        n
    }
}

fn New_REG(
    lParse: &mut ParseData,
    fname: *mut c_char,
    mut NodeX: c_int,
    mut NodeY: c_int,
    mut colNames: *mut c_char,
) -> c_int {
    unsafe {
        let this_node_idx: usize;
        let mut that0: &mut Node;
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
        let mut cX: *mut c_char = ptr::null_mut();
        let mut cY: *mut c_char = ptr::null_mut();
        let mut colVal: FITS_PARSER_YYSTYPE = FITS_PARSER_YYSTYPE { Node: 0 };
        if NodeX == -99 {
            type_0 = fits_parser_yyGetVariable(lParse, cs!(c"X"), &mut colVal);
            if type_0 == fits_parser_yytokentype::COLUMN as c_int {
                NodeX = New_Column(lParse, colVal.lng as c_int);
            } else {
                fits_parser_yyerror(lParse, cs!(c"Could not build X column for REGFILTER"));
                return -(1);
            }
        }
        if NodeY == -99 {
            type_0 = fits_parser_yyGetVariable(lParse, cs!(c"Y"), &mut colVal);
            if type_0 == fits_parser_yytokentype::COLUMN as c_int {
                NodeY = New_Column(lParse, colVal.lng as c_int);
            } else {
                fits_parser_yyerror(lParse, cs!(c"Could not build Y column for REGFILTER"));
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
                lParse,
                cs!(c"Dimensions of REGFILTER arguments are not compatible"),
            );
            return -(1);
        }
        n = Alloc_Node(lParse);
        if n >= 0 {
            this_node_idx = n as usize;
            (lParse.Nodes[this_node_idx]).nSubNodes = 3 as c_int;
            (lParse.Nodes[this_node_idx]).SubNodes[0] = Node0.try_into().unwrap();
            (lParse.Nodes[this_node_idx]).SubNodes[1] = NodeX.try_into().unwrap();
            (lParse.Nodes[this_node_idx]).SubNodes[2] = NodeY.try_into().unwrap();
            (lParse.Nodes[this_node_idx]).operation = REGFILT_FCT as c_int;
            (lParse.Nodes[this_node_idx]).DoOp = Some(Do_REG);
            (lParse.Nodes[this_node_idx]).ntype = fits_parser_yytokentype::BOOLEAN as c_int;
            (lParse.Nodes[this_node_idx]).value.nelem = 1;
            (lParse.Nodes[this_node_idx]).value.naxis = 1;
            (lParse.Nodes[this_node_idx]).value.naxes[0] = 1;
            Copy_Dims(lParse, n, NodeX);
            if (((lParse.Nodes)[NodeX as usize]).value).nelem
                < (((lParse.Nodes)[NodeY as usize]).value).nelem
            {
                Copy_Dims(lParse, n, NodeY);
            }
            let that0_idx = Node0 as usize;
            (lParse.Nodes[that0_idx]).operation = CONST_OP;
            (lParse.Nodes[that0_idx]).DoOp = None;
            Ycol = 0;
            Xcol = Ycol;
            if *colNames != 0 {
                while c_int::from(*colNames) == ' ' as i32 {
                    colNames = colNames.offset(1);
                }

                cY = colNames;
                cX = cY;
                while c_int::from(*cY) != 0
                    && c_int::from(*cY) != ' ' as i32
                    && c_int::from(*cY) != ',' as i32
                {
                    cY = cY.offset(1);
                }
                if *cY != 0 {
                    let fresh10 = cY;
                    cY = cY.offset(1);
                    *fresh10 = 0;
                }
                while c_int::from(*cY) == ' ' as i32 {
                    cY = cY.offset(1);
                }
                if *cY == 0 {
                    fits_parser_yyerror(
                        lParse,
                        cs!(c"Could not extract valid pair of column names from REGFILTER"),
                    );
                    Free_Last_Node(lParse);
                    return -(1);
                }

                let fptr = lParse.def_fptr.as_mut().unwrap();
                let cX_str = CStr::from_ptr(cX).to_bytes_with_nul();
                let cY_str = CStr::from_ptr(cY).to_bytes_with_nul();
                ffgcno_safe(
                    fptr,
                    0,
                    core::slice::from_raw_parts(cX_str.as_ptr().cast::<c_char>(), cX_str.len()),
                    &mut Xcol,
                    &mut lParse.status,
                );
                ffgcno_safe(
                    fptr,
                    0,
                    core::slice::from_raw_parts(cY_str.as_ptr().cast::<c_char>(), cY_str.len()),
                    &mut Ycol,
                    &mut lParse.status,
                );
                if lParse.status != 0 {
                    fits_parser_yyerror(
                        lParse,
                        cs!(c"Could not locate columns indicated for WCS info"),
                    );
                    Free_Last_Node(lParse);
                    return -(1);
                }
            } else {
                Xcol = Locate_Col(lParse, &(lParse.Nodes)[NodeX as usize]);
                Ycol = Locate_Col(lParse, &(lParse.Nodes)[NodeY as usize]);
                if Xcol < 0 || Ycol < 0 {
                    fits_parser_yyerror(
                        lParse,
                        cs!(c"Found multiple X/Y column references in REGFILTER"),
                    );
                    Free_Last_Node(lParse);
                    return -(1);
                }
            }
            wcs.exists = false;
            if Xcol > 0 && Ycol > 0 {
                tstat = 0;
                let fptr = lParse.def_fptr.as_mut().unwrap();
                ffgtcs_safe(
                    fptr,
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
                    lParse.status = tstat;
                    Free_Last_Node(lParse);
                    return -(1);
                } else {
                    wcs.exists = true;
                }
            }

            /* the region is allocated by fits_read_rgnfile; hand ownership of
            it to the node, which frees it again in ffcprs */
            let fname_slice = CStr::from_ptr(fname);
            let mut rgn: Option<Box<SAORegion>> = None;
            fits_read_rgnfile(
                cast_slice(fname_slice.to_bytes_with_nul()),
                &mut wcs,
                &mut rgn,
                &mut lParse.status,
            );
            if lParse.status != 0 {
                Free_Last_Node(lParse);
                return -(1);
            }
            let Rgn: *mut SAORegion = match rgn {
                Some(rgn) => Box::into_raw(rgn),
                None => ptr::null_mut(),
            };
            (lParse.Nodes[that0_idx]).value.data.ptr = Rgn.cast::<c_void>();
            if ((lParse.Nodes)[NodeX as usize]).operation == CONST_OP
                && ((lParse.Nodes)[NodeY as usize]).operation == CONST_OP
            {
                ((lParse.Nodes[this_node_idx]).DoOp).expect("non-null function pointer")(
                    lParse,
                    this_node_idx,
                );
            }
        }
        n
    }
}

fn New_Vector(lParse: &mut ParseData, subNode: c_int) -> c_int {
    let mut n: c_int = 0;
    n = Alloc_Node(lParse);
    if n >= 0 {
        let this_node_idx = n as usize;
        let that = &mut (lParse.Nodes)[subNode as usize];
        (lParse.Nodes[this_node_idx]).ntype = that.ntype;
        (lParse.Nodes[this_node_idx]).nSubNodes = 1;
        (lParse.Nodes[this_node_idx]).SubNodes[0] = subNode.try_into().unwrap();
        (lParse.Nodes[this_node_idx]).operation = '{' as i32;
        (lParse.Nodes[this_node_idx]).DoOp = Some(Do_Vector);
    }
    n
}

fn Close_Vec(lParse: &mut ParseData, vecNode: c_int) -> c_int {
    let mut n: c_int = 0;
    let mut nelem: c_int = 0;

    let mut this_node_idx: usize = vecNode as usize;
    n = 0;
    while n < (lParse.Nodes[this_node_idx]).nSubNodes {
        let mut subnode = lParse.Nodes[this_node_idx].SubNodes[n as usize];
        if ((lParse.Nodes)[subnode as usize]).ntype != (lParse.Nodes[this_node_idx]).ntype {
            /* New_Unary may change the lParse->Nodes pointer if
            it performs a realloc. Therefore reset 'this' just in case. */

            subnode = New_Unary(
                lParse,
                (lParse.Nodes[this_node_idx]).ntype,
                0,
                (lParse.Nodes[this_node_idx]).SubNodes[n as usize] as c_int,
            )
            .try_into()
            .unwrap();

            if subnode < 0 {
                return -(1);
            }

            this_node_idx = vecNode as usize;
            lParse.Nodes[this_node_idx].SubNodes[n as usize] = subnode;
        }
        nelem = (c_long::from(nelem)
            + ((lParse.Nodes)[(lParse.Nodes[this_node_idx]).SubNodes[n as usize] as usize])
                .value
                .nelem) as c_int;
        n += 1;
    }
    (lParse.Nodes[this_node_idx]).value.naxis = 1;
    (lParse.Nodes[this_node_idx]).value.nelem = c_long::from(nelem);
    (lParse.Nodes[this_node_idx]).value.naxes[0] = c_long::from(nelem);
    vecNode
}

fn New_Array(lParse: &mut ParseData, valueNode: c_int, mut dimNode: c_int) -> c_int {
    unsafe {
        let mut naxis: c_long = 0;
        let mut nelem: c_long = 0;
        let mut naxes: [c_long; 5] = [0; 5];
        let this_node_idx: usize;
        let mut n: c_int = 0;
        let mut i: c_int = 0;
        if valueNode < 0 || dimNode < 0 {
            return -(1);
        }

        /* Check that dimensions are {a,b,c,d}
             - vector
         - every element is constant integer
         - 5 or fewer dimensions
        */

        let dims: usize = dimNode as usize;
        i = 0;
        while i < MAXDIMS as c_int {
            naxes[i as usize] = 1;
            i += 1;
        }

        if ((lParse.Nodes)[dimNode as usize]).operation == CONST_OP {
            /* ARRAY(V,n) is a constant integer */

            if ((lParse.Nodes)[dimNode as usize]).ntype != fits_parser_yytokentype::LONG as c_int {
                dimNode = New_Unary(lParse, fits_parser_yytokentype::LONG as c_int, 0, dimNode);
            }
            if dimNode < 0 {
                return -(1);
            }
            naxis = 1;
            naxes[0] = (((lParse.Nodes)[dimNode as usize]).value).data.lng;
        } else if ((lParse.Nodes)[dimNode as usize]).operation == '{' as i32 {
            /* ARRAY(V,{a,b,c,d,e}) up to 5 dimensions */

            if (lParse.Nodes[dims]).nSubNodes > 5 as c_int {
                fits_parser_yyerror(
                    lParse,
                    cs!(c"ARRAY(V,{...}) number of dimensions must not exceed 5"),
                );
                return -(1);
            }
            naxis = c_long::from((lParse.Nodes[dims]).nSubNodes);
            i = 0;
            while i < (lParse.Nodes[dims]).nSubNodes {
                if ((lParse.Nodes)[(lParse.Nodes[dims]).SubNodes[i as usize] as usize]).ntype
                    != fits_parser_yytokentype::LONG as c_int
                {
                    (lParse.Nodes[dims]).SubNodes[i as usize] = New_Unary(
                        lParse,
                        fits_parser_yytokentype::LONG as c_int,
                        0,
                        (lParse.Nodes[dims]).SubNodes[i as usize] as c_int,
                    )
                    .try_into()
                    .unwrap();

                    /*
                    if (lParse.Nodes[dims]).SubNodes[i as usize] < 0 {
                        return -(1);
                    }
                    */
                }
                naxes[i as usize] = ((lParse.Nodes)
                    [(lParse.Nodes[dims]).SubNodes[i as usize] as usize])
                    .value
                    .data
                    .lng;
                i += 1;
            }
        } else {
            fits_parser_yyerror(
                lParse,
                cs!(c"ARRAY(V,dims) dims must be either integer or const vector"),
            );
            return -(1);
        }

        nelem = 1;
        i = 0;
        while c_long::from(i) < naxis {
            if naxes[i as usize] <= 0 {
                fits_parser_yyerror(lParse, cs!(c"ARRAY(V,dims) must have positive dimensions"));
                return -(1);
            }
            nelem *= naxes[i as usize];
            i += 1;
        }

        if !((((lParse.Nodes)[valueNode as usize]).value).nelem == nelem && nelem > 1) {
            if (((lParse.Nodes)[valueNode as usize]).value).nelem > 1 && nelem > 1 {
                fits_parser_yyerror(
                    lParse,
                    cs!(c"ARRAY(V,d) mismatch between number of elements in V and d"),
                );
                return -(1);
            } else if (((lParse.Nodes)[valueNode as usize]).value).nelem > 1 {
                fits_parser_yyerror(
                    lParse,
                    cs!(c"ARRAY(V,n) value V must have vector dimension of 1"),
                );
                return -(1);
            }
        }
        n = Alloc_Node(lParse);
        if n >= 0 {
            this_node_idx = n as usize;
            (lParse.Nodes[this_node_idx]).operation = ARRAY_FCT as c_int;
            (lParse.Nodes[this_node_idx]).nSubNodes = 1;
            (lParse.Nodes[this_node_idx]).SubNodes[0] = valueNode.try_into().unwrap();
            (lParse.Nodes[this_node_idx]).ntype = ((lParse.Nodes)[valueNode as usize]).ntype;
            (lParse.Nodes[this_node_idx]).value.nelem = nelem;
            (lParse.Nodes[this_node_idx]).value.naxis = naxis as c_int;
            i = 0;
            while c_long::from(i) < naxis {
                (lParse.Nodes[this_node_idx]).value.naxes[i as usize] = naxes[i as usize];
                i += 1;
            }
            (lParse.Nodes[this_node_idx]).DoOp = Some(Do_Array);
        }
        n
    }
}

/*  Locate the TABLE column number of any columns in "this" calculation.  */
/*  Return ZERO if none found, or negative if more than 1 found.          */
fn Locate_Col(lParse: &ParseData, this: &Node) -> c_int {
    let mut i: c_int = 0;
    let mut col: c_int = 0;
    let mut newCol: c_int = 0;
    let mut nfound: c_int = 0;

    if this.nSubNodes == 0 && this.operation <= 0 && this.operation != CONST_OP {
        return ((lParse.colData)[(-this.operation) as usize]).colnum;
    }

    i = 0;
    while i < this.nSubNodes {
        let that = &(lParse.Nodes)[this.SubNodes[i as usize] as usize];
        if that.operation > 0 {
            newCol = Locate_Col(lParse, that);
            if newCol <= 0 {
                nfound += -newCol;
            } else if nfound == 0 {
                col = newCol;
                nfound += 1;
            } else if col != newCol {
                nfound += 1;
            }
        } else if that.operation != CONST_OP {
            newCol = ((lParse.colData)[(-that.operation) as usize]).colnum;
            if nfound == 0 {
                col = newCol;
                nfound += 1;
            } else if col != newCol {
                nfound += 1;
            }
        }
        i += 1;
    }
    if nfound != 1 { -nfound } else { col }
}

fn Test_Dims(lParse: &mut ParseData, Node1: c_int, Node2: c_int) -> c_int {
    let mut valid: c_int = 0;
    let mut i: c_int = 0;

    if Node1 < 0 || Node2 < 0 {
        return 0;
    }

    let that1_idx = Node1 as usize;
    let that2_idx = Node2 as usize;

    if (lParse.Nodes[that1_idx]).value.nelem == 1 || (lParse.Nodes[that2_idx]).value.nelem == 1 {
        valid = 1;
    } else if (lParse.Nodes[that1_idx]).ntype == (lParse.Nodes[that2_idx]).ntype
        && (lParse.Nodes[that1_idx]).value.nelem == (lParse.Nodes[that2_idx]).value.nelem
        && (lParse.Nodes[that1_idx]).value.naxis == (lParse.Nodes[that2_idx]).value.naxis
    {
        valid = 1;
        i = 0;
        while i < (lParse.Nodes[that1_idx]).value.naxis {
            if (lParse.Nodes[that1_idx]).value.naxes[i as usize]
                != (lParse.Nodes[that2_idx]).value.naxes[i as usize]
            {
                valid = 0;
            }
            i += 1;
        }
    } else {
        valid = 0;
    }
    valid
}

fn Copy_Dims(lParse: &mut ParseData, Node1: c_int, Node2: c_int) {
    let mut i: c_int = 0;

    if Node1 < 0 || Node2 < 0 {
        return;
    }

    let that1_idx = Node1 as usize;
    let that2_idx = Node2 as usize;

    (lParse.Nodes[that1_idx]).value.nelem = (lParse.Nodes[that2_idx]).value.nelem;
    (lParse.Nodes[that1_idx]).value.naxis = (lParse.Nodes[that2_idx]).value.naxis;
    i = 0;
    while i < (lParse.Nodes[that2_idx]).value.naxis {
        (lParse.Nodes[that1_idx]).value.naxes[i as usize] =
            (lParse.Nodes[that2_idx]).value.naxes[i as usize];
        i += 1;
    }
}

pub(crate) fn Evaluate_Parser(lParse: &mut ParseData, firstRow: c_long, nRows: c_long) {
    unsafe {
        let mut i: c_int = 0;
        let mut column: c_int = 0;
        let mut offset: c_long = 0;
        let mut rowOffset: c_long = 0;
        static mut RAND_INITIALIZED: c_int = 0;
        if RAND_INITIALIZED == 0 {
            simplerng_srand(time(core::ptr::null_mut::<time_t>()) as c_uint);
            RAND_INITIALIZED = 1;
        }
        lParse.firstRow = firstRow;
        lParse.nRows = nRows;
        rowOffset = firstRow - lParse.firstDataRow;
        i = 0;
        while i < lParse.nNodes {
            if !(((lParse.Nodes)[i as usize]).operation > 0
                || ((lParse.Nodes)[i as usize]).operation == CONST_OP)
            {
                column = -((lParse.Nodes)[i as usize]).operation;
                offset = ((lParse.varData)[column as usize]).nelem * rowOffset;

                (((lParse.Nodes)[i as usize]).value).undef =
                    match (((lParse.varData)[column as usize]).undef).as_deref_mut() {
                        Some(ud) => ud[(offset as usize)..].as_mut_ptr(),
                        None => ptr::null_mut(),
                    };

                match ((lParse.Nodes)[i as usize]).ntype.into() {
                    fits_parser_yytokentype::BITSTR => {
                        let fresh12 = &mut (((lParse.Nodes)[i as usize]).value).data.strptr;
                        *fresh12 = ((lParse.varData)[column as usize])
                            .data
                            .cast::<*mut c_char>()
                            .offset(rowOffset as isize);
                        let fresh13 = &mut (((lParse.Nodes)[i as usize]).value).undef;
                        *fresh13 = ptr::null_mut();
                    }
                    fits_parser_yytokentype::STRING => {
                        let fresh14 = &mut (((lParse.Nodes)[i as usize]).value).data.strptr;
                        *fresh14 = ((lParse.varData)[column as usize])
                            .data
                            .cast::<*mut c_char>()
                            .offset(rowOffset as isize);
                        let fresh15 = &mut (((lParse.Nodes)[i as usize]).value).undef;
                        *fresh15 = (((lParse.varData)[column as usize]).undef)
                            .as_deref_mut()
                            .unwrap()[(rowOffset as usize)..]
                            .as_mut_ptr();
                    }
                    fits_parser_yytokentype::BOOLEAN => {
                        let fresh16 = &mut (((lParse.Nodes)[i as usize]).value).data.logptr;
                        *fresh16 = (((lParse.varData)[column as usize]).data as *const _
                            as *mut c_char)
                            .offset(offset as isize);
                    }
                    fits_parser_yytokentype::LONG => {
                        let fresh17 = &mut (((lParse.Nodes)[i as usize]).value).data.lngptr;
                        *fresh17 = (((lParse.varData)[column as usize]).data as *const _
                            as *mut c_long)
                            .offset(offset as isize);
                    }
                    fits_parser_yytokentype::DOUBLE => {
                        let fresh18 = &mut (((lParse.Nodes)[i as usize]).value).data.dblptr;
                        *fresh18 = (((lParse.varData)[column as usize]).data as *const _
                            as *mut c_double)
                            .offset(offset as isize);
                    }
                    _ => {}
                }
            }
            i += 1;
        }
        Evaluate_Node(lParse, lParse.resultNode);
    }
}

/**********************************************************************/
/*  Recursively evaluate thisNode's subNodes, then call one of the    */
/*  Do_<Action> functions pointed to by thisNode's DoOp element.      */
/**********************************************************************/
fn Evaluate_Node(lParse: &mut ParseData, thisNode: c_int) {
    let mut this_node_idx: usize;
    let mut i: c_int = 0;
    if lParse.status != 0 {
        return;
    }

    if (lParse.Nodes)[thisNode as usize].operation > 0 {
        /* <=0 indicate constants and columns */
        i = ((lParse.Nodes)[thisNode as usize]).nSubNodes;
        loop {
            if i == 0 {
                i -= 1;
                break;
            }

            i -= 1;

            Evaluate_Node(
                lParse,
                (lParse.Nodes)[thisNode as usize].SubNodes[i as usize] as c_int,
            );

            if lParse.status != 0 {
                return;
            }
        }
        ((lParse.Nodes[thisNode as usize]).DoOp).expect("non-null function pointer")(
            lParse,
            thisNode as usize,
        );
    }
}

fn Allocate_Ptrs(lParse: &mut ParseData, this_node_idx: usize) {
    unsafe {
        let mut elem: c_long = 0;
        let mut row: c_long = 0;
        let mut size: c_long = 0;
        if (lParse.Nodes[this_node_idx]).ntype == fits_parser_yytokentype::BITSTR as c_int
            || (lParse.Nodes[this_node_idx]).ntype == fits_parser_yytokentype::STRING as c_int
        {
            (lParse.Nodes[this_node_idx]).value.data.strptr = malloc(
                (lParse.nRows as c_ulong)
                    .wrapping_mul(::core::mem::size_of::<*mut c_char>() as c_ulong)
                    .try_into()
                    .unwrap(),
            )
            .cast::<*mut c_char>();
            if !((lParse.Nodes[this_node_idx]).value.data.strptr).is_null() {
                let fresh20 = &mut *((lParse.Nodes[this_node_idx]).value.data.strptr).offset(0);
                *fresh20 = malloc(
                    ((lParse.nRows * ((lParse.Nodes[this_node_idx]).value.nelem + 2 as c_long))
                        as c_ulong)
                        .wrapping_mul(::core::mem::size_of::<c_char>() as c_ulong)
                        .try_into()
                        .unwrap(),
                )
                .cast::<c_char>();
                if !(*((lParse.Nodes[this_node_idx]).value.data.strptr).offset(0)).is_null() {
                    row = 0;
                    loop {
                        row += 1;
                        if row >= lParse.nRows {
                            break;
                        }
                        let fresh21 = &mut *((lParse.Nodes[this_node_idx]).value.data.strptr)
                            .offset(row as isize);
                        *fresh21 = (*((lParse.Nodes[this_node_idx]).value.data.strptr)
                            .offset((row - 1) as isize))
                        .offset((lParse.Nodes[this_node_idx]).value.nelem as isize)
                        .offset(1);
                    }
                    if (lParse.Nodes[this_node_idx]).ntype
                        == fits_parser_yytokentype::STRING as c_int
                    {
                        (lParse.Nodes[this_node_idx]).value.undef =
                            (*((lParse.Nodes[this_node_idx]).value.data.strptr)
                                .offset((row - 1) as isize))
                            .offset((lParse.Nodes[this_node_idx]).value.nelem as isize)
                            .offset(1);
                    } else {
                        (lParse.Nodes[this_node_idx]).value.undef = ptr::null_mut(); /* BITSTRs don't use undef array */
                    }
                } else {
                    lParse.status = MEMORY_ALLOCATION;
                    free(
                        (lParse.Nodes[this_node_idx])
                            .value
                            .data
                            .strptr
                            .cast::<c_void>(),
                    );
                }
            } else {
                lParse.status = MEMORY_ALLOCATION;
            }
        } else {
            elem = (lParse.Nodes[this_node_idx]).value.nelem * lParse.nRows;
            match (lParse.Nodes[this_node_idx]).ntype.into() {
                fits_parser_yytokentype::DOUBLE => {
                    size = ::core::mem::size_of::<c_double>() as c_ulong as c_long;
                }
                fits_parser_yytokentype::LONG => {
                    size = ::core::mem::size_of::<c_long>() as c_ulong as c_long;
                }
                fits_parser_yytokentype::BOOLEAN => {
                    size = ::core::mem::size_of::<c_char>() as c_ulong as c_long;
                }
                _ => {
                    size = 1;
                }
            }
            (lParse.Nodes[this_node_idx]).value.data.ptr = calloc(
                ((size + 1) as c_ulong).try_into().unwrap(),
                (elem as c_ulong).try_into().unwrap(),
            );
            if ((lParse.Nodes[this_node_idx]).value.data.ptr).is_null() {
                lParse.status = MEMORY_ALLOCATION;
            } else {
                (lParse.Nodes[this_node_idx]).value.undef = (lParse.Nodes[this_node_idx])
                    .value
                    .data
                    .ptr
                    .cast::<c_char>()
                    .offset((elem * size) as isize);
            }
        };
    }
}

unsafe fn free_node_buffer(node: &mut Node) {
    if (node.ntype == fits_parser_yytokentype::BITSTR as c_int
        || node.ntype == fits_parser_yytokentype::STRING as c_int)
    {
        if (!node.value.data.strptr.is_null()) {
            if (!(*node.value.data.strptr).is_null()) {
                free((*node.value.data.strptr) as *mut c_void);
            }
            free(node.value.data.strptr as *mut c_void);
            node.value.data.strptr = ptr::null_mut();
        }
    } else if (!node.value.data.ptr.is_null()) {
        free(node.value.data.ptr);
        node.value.data.ptr = ptr::null_mut();
    }

    node.value.undef = ptr::null_mut();
}

fn Do_Unary(lParse: &mut ParseData, this_node_idx: usize) {
    unsafe {
        let mut that: &mut Node;
        let mut elem: c_long = 0;

        let that_idx = (lParse.Nodes[this_node_idx]).SubNodes[0];
        let (this_node, that_node) =
            get_this_that_nodes(&mut lParse.Nodes, this_node_idx, that_idx);

        if that_node.operation == CONST_OP {
            /* Operating on a constant! */
            match (this_node).operation {
                x if x == fits_parser_yytokentype::DOUBLE as c_int
                    || x == fits_parser_yytokentype::FLTCAST as c_int =>
                {
                    if that_node.ntype == fits_parser_yytokentype::LONG as c_int {
                        (this_node).value.data.dbl = that_node.value.data.lng as c_double;
                    } else if that_node.ntype == fits_parser_yytokentype::BOOLEAN as c_int {
                        (this_node).value.data.dbl = if c_int::from(that_node.value.data.log) != 0 {
                            1.0
                        } else {
                            0.0
                        };
                    }
                }
                x if x == fits_parser_yytokentype::LONG as c_int
                    || x == fits_parser_yytokentype::INTCAST as c_int =>
                {
                    if that_node.ntype == fits_parser_yytokentype::DOUBLE as c_int {
                        (this_node).value.data.lng = that_node.value.data.dbl as c_long;
                    } else if that_node.ntype == fits_parser_yytokentype::BOOLEAN as c_int {
                        (this_node).value.data.lng = if c_int::from(that_node.value.data.log) != 0 {
                            1
                        } else {
                            0
                        };
                    }
                }
                x if x == fits_parser_yytokentype::BOOLEAN as c_int => {
                    if that_node.ntype == fits_parser_yytokentype::DOUBLE as c_int {
                        (this_node).value.data.log = if that_node.value.data.dbl != 0.0 {
                            1
                        } else {
                            0
                        };
                    } else if that_node.ntype == fits_parser_yytokentype::LONG as c_int {
                        (this_node).value.data.log =
                            if that_node.value.data.lng != 0 { 1 } else { 0 };
                    }
                }
                x if x == fits_parser_yytokentype::UMINUS as c_int => {
                    if that_node.ntype == fits_parser_yytokentype::DOUBLE as c_int {
                        (this_node).value.data.dbl = -that_node.value.data.dbl;
                    } else if that_node.ntype == fits_parser_yytokentype::LONG as c_int {
                        (this_node).value.data.lng = -that_node.value.data.lng;
                    }
                }
                x if x == fits_parser_yytokentype::NOT as c_int => {
                    if that_node.ntype == fits_parser_yytokentype::BOOLEAN as c_int {
                        (this_node).value.data.log =
                            if that_node.value.data.log == 0 { 1 } else { 0 };
                    } else if that_node.ntype == fits_parser_yytokentype::BITSTR as c_int {
                        bitnot(
                            ((this_node).value.data.astr).as_mut_ptr(),
                            (that_node.value.data.astr).as_mut_ptr(),
                        );
                    }
                }
                _ => {}
            }
            (this_node).operation = CONST_OP;
        } else {
            Allocate_Ptrs(lParse, this_node_idx);
            let (this_node, that_node) =
                get_this_that_nodes(&mut lParse.Nodes, this_node_idx, that_idx);

            if lParse.status == 0 {
                if (this_node).ntype != fits_parser_yytokentype::BITSTR as c_int {
                    elem = lParse.nRows;
                    if (this_node).ntype != fits_parser_yytokentype::STRING as c_int {
                        elem *= (this_node).value.nelem;
                    }
                    loop {
                        let fresh22 = elem;
                        elem -= 1;
                        if fresh22 == 0 {
                            break;
                        }
                        *((this_node).value.undef).offset(elem as isize) =
                            *(that_node.value.undef).offset(elem as isize);
                    }
                }
                elem = lParse.nRows * (this_node).value.nelem;
                match (this_node).operation.into() {
                    fits_parser_yytokentype::BOOLEAN => {
                        if that_node.ntype == fits_parser_yytokentype::DOUBLE as c_int {
                            loop {
                                let fresh23 = elem;
                                elem -= 1;
                                if fresh23 == 0 {
                                    break;
                                }
                                *((this_node).value.data.logptr).offset(elem as isize) =
                                    if *(that_node.value.data.dblptr).offset(elem as isize) != 0.0 {
                                        1
                                    } else {
                                        0
                                    };
                            }
                        } else if that_node.ntype == fits_parser_yytokentype::LONG as c_int {
                            loop {
                                let fresh24 = elem;
                                elem -= 1;
                                if fresh24 == 0 {
                                    break;
                                }
                                *((this_node).value.data.logptr).offset(elem as isize) =
                                    if *(that_node.value.data.lngptr).offset(elem as isize) != 0 {
                                        1
                                    } else {
                                        0
                                    };
                            }
                        }
                    }
                    fits_parser_yytokentype::DOUBLE | fits_parser_yytokentype::FLTCAST => {
                        if that_node.ntype == fits_parser_yytokentype::LONG as c_int {
                            loop {
                                let fresh25 = elem;
                                elem -= 1;
                                if fresh25 == 0 {
                                    break;
                                }
                                *((this_node).value.data.dblptr).offset(elem as isize) =
                                    *(that_node.value.data.lngptr).offset(elem as isize)
                                        as c_double;
                            }
                        } else if that_node.ntype == fits_parser_yytokentype::BOOLEAN as c_int {
                            loop {
                                let fresh26 = elem;
                                elem -= 1;
                                if fresh26 == 0 {
                                    break;
                                }
                                *((this_node).value.data.dblptr).offset(elem as isize) =
                                    if c_int::from(
                                        *(that_node.value.data.logptr).offset(elem as isize),
                                    ) != 0
                                    {
                                        1.0
                                    } else {
                                        0.0
                                    };
                            }
                        }
                    }
                    fits_parser_yytokentype::LONG | fits_parser_yytokentype::INTCAST => {
                        if that_node.ntype == fits_parser_yytokentype::DOUBLE as c_int {
                            loop {
                                let fresh27 = elem;
                                elem -= 1;
                                if fresh27 == 0 {
                                    break;
                                }
                                *((this_node).value.data.lngptr).offset(elem as isize) =
                                    *(that_node.value.data.dblptr).offset(elem as isize) as c_long;
                            }
                        } else if that_node.ntype == fits_parser_yytokentype::BOOLEAN as c_int {
                            loop {
                                let fresh28 = elem;
                                elem -= 1;
                                if fresh28 == 0 {
                                    break;
                                }
                                *((this_node).value.data.lngptr).offset(elem as isize) =
                                    if c_int::from(
                                        *(that_node.value.data.logptr).offset(elem as isize),
                                    ) != 0
                                    {
                                        1
                                    } else {
                                        0
                                    };
                            }
                        }
                    }
                    fits_parser_yytokentype::UMINUS => {
                        if that_node.ntype == fits_parser_yytokentype::DOUBLE as c_int {
                            loop {
                                let fresh29 = elem;
                                elem -= 1;
                                if fresh29 == 0 {
                                    break;
                                }
                                *((this_node).value.data.dblptr).offset(elem as isize) =
                                    -*(that_node.value.data.dblptr).offset(elem as isize);
                            }
                        } else if that_node.ntype == fits_parser_yytokentype::LONG as c_int {
                            loop {
                                let fresh30 = elem;
                                elem -= 1;
                                if fresh30 == 0 {
                                    break;
                                }
                                *((this_node).value.data.lngptr).offset(elem as isize) =
                                    -*(that_node.value.data.lngptr).offset(elem as isize);
                            }
                        }
                    }
                    fits_parser_yytokentype::NOT => {
                        if that_node.ntype == fits_parser_yytokentype::BOOLEAN as c_int {
                            loop {
                                let fresh31 = elem;
                                elem -= 1;
                                if fresh31 == 0 {
                                    break;
                                }
                                *((this_node).value.data.logptr).offset(elem as isize) = c_int::from(
                                    *(that_node.value.data.logptr).offset(elem as isize) == 0,
                                )
                                    as c_char;
                            }
                        } else if that_node.ntype == fits_parser_yytokentype::BITSTR as c_int {
                            elem = lParse.nRows;
                            loop {
                                let fresh32 = elem;
                                elem -= 1;
                                if fresh32 == 0 {
                                    break;
                                }
                                bitnot(
                                    *((this_node).value.data.strptr).offset(elem as isize),
                                    *(that_node.value.data.strptr).offset(elem as isize),
                                );
                            }
                        }
                    }
                    _ => {}
                }
            }
        }

        let that_idx = (lParse.Nodes[this_node_idx]).SubNodes[0];

        if (lParse.Nodes[that_idx]).operation > 0 {
            free((lParse.Nodes[that_idx]).value.data.ptr);
        }
    }
}

fn Do_Offset(lParse: &mut ParseData, this_node_idx: usize) {
    unsafe {
        let mut fRow: c_long = 0;
        let mut nRowOverlap: c_long = 0;
        let mut nRowReload: c_long = 0;
        let rowOffset: c_long = 0;
        let mut nelem: c_long = 0;
        let mut elem: c_long = 0;
        let mut offset: c_long = 0;
        let mut nRealElem: c_long = 0;
        let mut status: c_int = 0;

        let col_idx = (lParse.Nodes[this_node_idx]).SubNodes[0];
        let rowOffset = ((lParse.Nodes)[(lParse.Nodes[this_node_idx]).SubNodes[1]])
            .value
            .data
            .lng;

        Allocate_Ptrs(lParse, this_node_idx);

        fRow = lParse.firstRow + rowOffset;
        if (lParse.Nodes[this_node_idx]).ntype == fits_parser_yytokentype::STRING as c_int
            || (lParse.Nodes[this_node_idx]).ntype == fits_parser_yytokentype::BITSTR as c_int
        {
            nRealElem = 1;
        } else {
            nRealElem = (lParse.Nodes[this_node_idx]).value.nelem;
        }

        nelem = nRealElem;

        if ((fRow >= 0 && (LONG_MAX - lParse.nRows < fRow))
            || (fRow < 0 && (LONG_MIN + lParse.firstDataRow + 1 > fRow)))
        {
            fits_parser_yyerror(
                lParse,
                cs!(c"numerical underflow or overflow for row offset value"),
            );
            if (lParse.status == 0) {
                lParse.status = PARSE_SYNTAX_ERR;
            }
            free_node_buffer(&mut (lParse.Nodes[this_node_idx]));
            return;
        }

        if fRow < lParse.firstDataRow {
            nRowReload = lParse.firstDataRow - fRow;
            if nRowReload > lParse.nRows {
                nRowReload = lParse.nRows;
            }
            nRowOverlap = lParse.nRows - nRowReload;
            offset = 0;
            while fRow < 1 && nRowReload > 0 {
                if (lParse.Nodes[this_node_idx]).ntype == fits_parser_yytokentype::BITSTR as c_int {
                    nelem = (lParse.Nodes[this_node_idx]).value.nelem;
                    *(*((lParse.Nodes[this_node_idx]).value.data.strptr).offset(offset as isize))
                        .offset(nelem as isize) = 0;
                    loop {
                        let fresh33 = nelem;
                        nelem -= 1;
                        if fresh33 == 0 {
                            break;
                        }
                        *(*((lParse.Nodes[this_node_idx]).value.data.strptr)
                            .offset(offset as isize))
                        .offset(nelem as isize) = b'0' as c_char;
                    }
                    offset += 1;
                } else {
                    loop {
                        let fresh34 = nelem;
                        nelem -= 1;
                        if fresh34 == 0 {
                            break;
                        }
                        let fresh35 = offset;
                        offset += 1;
                        *((lParse.Nodes[this_node_idx]).value.undef).offset(fresh35 as isize) = 1;
                    }
                }
                nelem = nRealElem;
                fRow += 1;
                nRowReload -= 1;
            }
        } else if fRow + lParse.nRows > lParse.firstDataRow + lParse.nDataRows {
            nRowReload = fRow + lParse.nRows - (lParse.firstDataRow + lParse.nDataRows);
            if nRowReload > lParse.nRows {
                nRowReload = lParse.nRows;
            } else {
                fRow = lParse.firstDataRow + lParse.nDataRows;
            }

            if (rowOffset > 0) {
                if (rowOffset > LONG_MAX / nelem) {
                    fits_parser_yyerror(
                        lParse,
                        cs!(c"numerical overflow for row offset * nelem value"),
                    );
                    if (lParse.status == 0) {
                        lParse.status = PARSE_SYNTAX_ERR;
                    }
                    free_node_buffer(&mut (lParse.Nodes[this_node_idx]));
                    return;
                }
            } else if (rowOffset < 0) {
                if (rowOffset < LONG_MIN / nelem) {
                    fits_parser_yyerror(
                        lParse,
                        cs!(c"numerical underflow for row offset * nelem value"),
                    );
                    if (lParse.status == 0) {
                        lParse.status = PARSE_SYNTAX_ERR;
                    }
                    free_node_buffer(&mut (lParse.Nodes[this_node_idx]));
                    return;
                }
            }

            nRowOverlap = lParse.nRows - nRowReload;
            offset = nRowOverlap * nelem;
            elem = lParse.nRows * nelem;
            while fRow + nRowReload > lParse.totalRows && nRowReload > 0 {
                if (lParse.Nodes[this_node_idx]).ntype == fits_parser_yytokentype::BITSTR as c_int {
                    nelem = (lParse.Nodes[this_node_idx]).value.nelem;
                    elem -= 1;
                    *(*((lParse.Nodes[this_node_idx]).value.data.strptr).offset(elem as isize))
                        .offset(nelem as isize) = 0;
                    loop {
                        let fresh36 = nelem;
                        nelem -= 1;
                        if fresh36 == 0 {
                            break;
                        }
                        *(*((lParse.Nodes[this_node_idx]).value.data.strptr)
                            .offset(elem as isize))
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
                        *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize) = 1;
                    }
                }
                nelem = nRealElem;
                nRowReload -= 1;
            }
        } else {
            nRowReload = 0;
            nRowOverlap = lParse.nRows;
            offset = 0;
        }
        if nRowReload > 0 {
            match (lParse.Nodes[this_node_idx]).ntype.into() {
                fits_parser_yytokentype::BITSTR | fits_parser_yytokentype::STRING => {
                    status = (lParse.loadData).expect("non-null function pointer")(
                        lParse,
                        -(lParse.Nodes[col_idx]).operation,
                        fRow,
                        nRowReload,
                        ((lParse.Nodes[this_node_idx]).value.data.strptr)
                            .offset(offset as isize)
                            .cast::<c_void>(),
                        ((lParse.Nodes[this_node_idx]).value.undef).offset(offset as isize),
                    );
                }
                fits_parser_yytokentype::BOOLEAN => {
                    status = (lParse.loadData).expect("non-null function pointer")(
                        lParse,
                        -(lParse.Nodes[col_idx]).operation,
                        fRow,
                        nRowReload,
                        ((lParse.Nodes[this_node_idx]).value.data.logptr)
                            .offset(offset as isize)
                            .cast::<c_void>(),
                        ((lParse.Nodes[this_node_idx]).value.undef).offset(offset as isize),
                    );
                }
                fits_parser_yytokentype::LONG => {
                    status = (lParse.loadData).expect("non-null function pointer")(
                        lParse,
                        -(lParse.Nodes[col_idx]).operation,
                        fRow,
                        nRowReload,
                        ((lParse.Nodes[this_node_idx]).value.data.lngptr)
                            .offset(offset as isize)
                            .cast::<c_void>(),
                        ((lParse.Nodes[this_node_idx]).value.undef).offset(offset as isize),
                    );
                }
                fits_parser_yytokentype::DOUBLE => {
                    status = (lParse.loadData).expect("non-null function pointer")(
                        lParse,
                        -(lParse.Nodes[col_idx]).operation,
                        fRow,
                        nRowReload,
                        ((lParse.Nodes[this_node_idx]).value.data.dblptr)
                            .offset(offset as isize)
                            .cast::<c_void>(),
                        ((lParse.Nodes[this_node_idx]).value.undef).offset(offset as isize),
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
            elem = lParse.nRows * nelem;
        }
        offset = nelem * rowOffset;
        loop {
            let fresh38 = nRowOverlap;
            nRowOverlap -= 1;
            if !(fresh38 != 0 && lParse.status == 0) {
                break;
            }
            loop {
                let fresh39 = nelem;
                nelem -= 1;
                if !(fresh39 != 0 && lParse.status == 0) {
                    break;
                }
                elem -= 1;

                if (lParse.Nodes[this_node_idx]).ntype != fits_parser_yytokentype::BITSTR as c_int {
                    *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize) =
                        *((lParse.Nodes[col_idx]).value.undef).offset((elem + offset) as isize);
                }
                match (lParse.Nodes[this_node_idx]).ntype.into() {
                    fits_parser_yytokentype::BITSTR => {
                        strcpy(
                            *((lParse.Nodes[this_node_idx]).value.data.strptr)
                                .offset(elem as isize),
                            *((lParse.Nodes[col_idx]).value.data.strptr)
                                .offset((elem + offset) as isize),
                        );
                    }
                    fits_parser_yytokentype::STRING => {
                        strcpy(
                            *((lParse.Nodes[this_node_idx]).value.data.strptr)
                                .offset(elem as isize),
                            *((lParse.Nodes[col_idx]).value.data.strptr)
                                .offset((elem + offset) as isize),
                        );
                    }
                    fits_parser_yytokentype::BOOLEAN => {
                        *((lParse.Nodes[this_node_idx]).value.data.logptr).offset(elem as isize) =
                            *((lParse.Nodes[col_idx]).value.data.logptr)
                                .offset((elem + offset) as isize);
                    }
                    fits_parser_yytokentype::LONG => {
                        *((lParse.Nodes[this_node_idx]).value.data.lngptr).offset(elem as isize) =
                            *((lParse.Nodes[col_idx]).value.data.lngptr)
                                .offset((elem + offset) as isize);
                    }
                    fits_parser_yytokentype::DOUBLE => {
                        *((lParse.Nodes[this_node_idx]).value.data.dblptr).offset(elem as isize) =
                            *((lParse.Nodes[col_idx]).value.data.dblptr)
                                .offset((elem + offset) as isize);
                    }
                    _ => {}
                }
            }
            nelem = nRealElem;
        }
    }
}
fn Do_BinOp_bit(lParse: &mut ParseData, this_node_idx: usize) {
    unsafe {
        let mut that1: &mut Node;
        let mut that2: &mut Node;
        let mut sptr1: *mut c_char = ptr::null_mut();
        let mut sptr2: *mut c_char = ptr::null_mut();
        let mut const1: c_int = 0;
        let mut const2: c_int = 0;
        let mut rows: c_long = 0;

        let that1_idx = (lParse.Nodes[this_node_idx]).SubNodes[0];
        let that2_idx = (lParse.Nodes[this_node_idx]).SubNodes[1];

        const1 = c_int::from((lParse.Nodes[that1_idx]).operation == CONST_OP);
        const2 = c_int::from((lParse.Nodes[that2_idx]).operation == CONST_OP);
        sptr1 = if const1 != 0 {
            ((lParse.Nodes[that1_idx]).value.data.astr).as_mut_ptr()
        } else {
            core::ptr::null_mut::<c_char>()
        };
        sptr2 = if const2 != 0 {
            ((lParse.Nodes[that2_idx]).value.data.astr).as_mut_ptr()
        } else {
            core::ptr::null_mut::<c_char>()
        };
        if const1 != 0 && const2 != 0 {
            match (lParse.Nodes[this_node_idx]).operation {
                280 => {
                    (lParse.Nodes[this_node_idx]).value.data.log =
                        if bitcmp(sptr1, sptr2) == 0 { 1 } else { 0 };
                }
                279 => {
                    (lParse.Nodes[this_node_idx]).value.data.log = bitcmp(sptr1, sptr2);
                }
                281..=284 => {
                    (lParse.Nodes[this_node_idx]).value.data.log =
                        bitlgte(sptr1, (lParse.Nodes[this_node_idx]).operation, sptr2);
                }
                124 => {
                    bitor(
                        ((lParse.Nodes[this_node_idx]).value.data.astr).as_mut_ptr(),
                        sptr1,
                        sptr2,
                    );
                }
                38 => {
                    bitand(
                        ((lParse.Nodes[this_node_idx]).value.data.astr).as_mut_ptr(),
                        sptr1,
                        sptr2,
                    );
                }
                43 => {
                    strcpy(
                        ((lParse.Nodes[this_node_idx]).value.data.astr).as_mut_ptr(),
                        sptr1,
                    );
                    strcat(
                        ((lParse.Nodes[this_node_idx]).value.data.astr).as_mut_ptr(),
                        sptr2,
                    );
                }
                291 => {
                    (lParse.Nodes[this_node_idx]).value.data.lng = 0;
                    while *sptr1 != 0 {
                        if c_int::from(*sptr1) == '1' as i32 {
                            (lParse.Nodes[this_node_idx]).value.data.lng += 1;
                            (lParse.Nodes[this_node_idx]).value.data.lng;
                        }
                        sptr1 = sptr1.offset(1);
                    }
                }
                _ => {}
            }
            (lParse.Nodes[this_node_idx]).operation = CONST_OP;
        } else {
            Allocate_Ptrs(lParse, this_node_idx);
            if lParse.status == 0 {
                rows = lParse.nRows;
                match (lParse.Nodes[this_node_idx]).operation {
                    279..=284 => loop {
                        let fresh40 = rows;
                        rows -= 1;
                        if fresh40 == 0 {
                            break;
                        }
                        if const1 == 0 {
                            sptr1 = *((lParse.Nodes[that1_idx]).value.data.strptr)
                                .offset(rows as isize);
                        }
                        if const2 == 0 {
                            sptr2 = *((lParse.Nodes[that2_idx]).value.data.strptr)
                                .offset(rows as isize);
                        }
                        match (lParse.Nodes[this_node_idx]).operation {
                            280 => {
                                *((lParse.Nodes[this_node_idx]).value.data.logptr)
                                    .offset(rows as isize) =
                                    if bitcmp(sptr1, sptr2) == 0 { 1 } else { 0 };
                            }
                            279 => {
                                *((lParse.Nodes[this_node_idx]).value.data.logptr)
                                    .offset(rows as isize) = bitcmp(sptr1, sptr2);
                            }
                            281..=284 => {
                                *((lParse.Nodes[this_node_idx]).value.data.logptr)
                                    .offset(rows as isize) =
                                    bitlgte(sptr1, (lParse.Nodes[this_node_idx]).operation, sptr2);
                            }
                            _ => {}
                        }
                        *((lParse.Nodes[this_node_idx]).value.undef).offset(rows as isize) = 0;
                    },
                    124 | 38 | 43 => loop {
                        let fresh41 = rows;
                        rows -= 1;
                        if fresh41 == 0 {
                            break;
                        }
                        if const1 == 0 {
                            sptr1 = *((lParse.Nodes[that1_idx]).value.data.strptr)
                                .offset(rows as isize);
                        }
                        if const2 == 0 {
                            sptr2 = *((lParse.Nodes[that2_idx]).value.data.strptr)
                                .offset(rows as isize);
                        }
                        if (lParse.Nodes[this_node_idx]).operation == '|' as i32 {
                            bitor(
                                *((lParse.Nodes[this_node_idx]).value.data.strptr)
                                    .offset(rows as isize),
                                sptr1,
                                sptr2,
                            );
                        } else if (lParse.Nodes[this_node_idx]).operation == '&' as i32 {
                            bitand(
                                *((lParse.Nodes[this_node_idx]).value.data.strptr)
                                    .offset(rows as isize),
                                sptr1,
                                sptr2,
                            );
                        } else {
                            strcpy(
                                *((lParse.Nodes[this_node_idx]).value.data.strptr)
                                    .offset(rows as isize),
                                sptr1,
                            );
                            strcat(
                                *((lParse.Nodes[this_node_idx]).value.data.strptr)
                                    .offset(rows as isize),
                                sptr2,
                            );
                        }
                    },
                    291 => {
                        let mut i: c_long = 0;
                        let mut previous: c_long = 0;
                        let mut curr: c_long = 0;
                        previous = (lParse.Nodes[that2_idx]).value.data.lng;
                        i = 0;
                        while i < rows {
                            sptr1 =
                                *((lParse.Nodes[that1_idx]).value.data.strptr).offset(i as isize);
                            curr = 0;
                            while *sptr1 != 0 {
                                if c_int::from(*sptr1) == '1' as i32 {
                                    curr += 1;
                                }
                                sptr1 = sptr1.offset(1);
                            }
                            previous += curr;
                            *((lParse.Nodes[this_node_idx]).value.data.lngptr).offset(i as isize) =
                                previous;
                            *((lParse.Nodes[this_node_idx]).value.undef).offset(i as isize) = 0;
                            i += 1;
                        }
                        (lParse.Nodes[that2_idx]).value.data.lng = previous;
                    }
                    _ => {}
                }
            }
        }
        if (lParse.Nodes[that1_idx]).operation > 0 {
            free((*((lParse.Nodes[that1_idx]).value.data.strptr).offset(0)).cast::<c_void>());
            free((lParse.Nodes[that1_idx]).value.data.strptr.cast::<c_void>());
        }
        if (lParse.Nodes[that2_idx]).operation > 0 {
            free((*((lParse.Nodes[that2_idx]).value.data.strptr).offset(0)).cast::<c_void>());
            free((lParse.Nodes[that2_idx]).value.data.strptr.cast::<c_void>());
        }
    }
}

fn Do_BinOp_str(lParse: &mut ParseData, this_node_idx: usize) {
    unsafe {
        let mut that1: &mut Node;
        let mut that2: &mut Node;
        let mut sptr1: *mut c_char = ptr::null_mut();
        let mut sptr2: *mut c_char = ptr::null_mut();
        let mut null1: c_char = 0;
        let mut null2: c_char = 0;
        let mut const1: c_int = 0;
        let mut const2: c_int = 0;
        let mut val: c_int = 0;
        let mut rows: c_long = 0;

        let that1_idx = (lParse.Nodes[this_node_idx]).SubNodes[0];
        let that2_idx = (lParse.Nodes[this_node_idx]).SubNodes[1];

        const1 = c_int::from((lParse.Nodes[that1_idx]).operation == CONST_OP);
        const2 = c_int::from((lParse.Nodes[that2_idx]).operation == CONST_OP);
        sptr1 = if const1 != 0 {
            ((lParse.Nodes[that1_idx]).value.data.astr).as_mut_ptr()
        } else {
            core::ptr::null_mut::<c_char>()
        };
        sptr2 = if const2 != 0 {
            ((lParse.Nodes[that2_idx]).value.data.astr).as_mut_ptr()
        } else {
            core::ptr::null_mut::<c_char>()
        };
        if const1 != 0 && const2 != 0 {
            match (lParse.Nodes[this_node_idx]).operation {
                280 | 279 => {
                    val = c_int::from(
                        (if c_int::from(*sptr1.offset(0)) < c_int::from(*sptr2.offset(0)) {
                            -(1)
                        } else if c_int::from(*sptr1.offset(0)) > c_int::from(*sptr2.offset(0)) {
                            1
                        } else {
                            strcmp(sptr1, sptr2)
                        }) == 0,
                    );
                    (lParse.Nodes[this_node_idx]).value.data.log = (if (lParse.Nodes[this_node_idx])
                        .operation
                        == fits_parser_yytokentype::EQ as c_int
                    {
                        val
                    } else {
                        c_int::from(val == 0)
                    }) as c_char;
                }
                281 => {
                    (lParse.Nodes[this_node_idx]).value.data.log =
                        if (if c_int::from(*sptr1.offset(0)) < c_int::from(*sptr2.offset(0)) {
                            -(1)
                        } else if c_int::from(*sptr1.offset(0)) > c_int::from(*sptr2.offset(0)) {
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
                    (lParse.Nodes[this_node_idx]).value.data.log =
                        if (if c_int::from(*sptr1.offset(0)) < c_int::from(*sptr2.offset(0)) {
                            -(1)
                        } else if c_int::from(*sptr1.offset(0)) > c_int::from(*sptr2.offset(0)) {
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
                    (lParse.Nodes[this_node_idx]).value.data.log =
                        if (if c_int::from(*sptr1.offset(0)) < c_int::from(*sptr2.offset(0)) {
                            -(1)
                        } else if c_int::from(*sptr1.offset(0)) > c_int::from(*sptr2.offset(0)) {
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
                    (lParse.Nodes[this_node_idx]).value.data.log =
                        if (if c_int::from(*sptr1.offset(0)) < c_int::from(*sptr2.offset(0)) {
                            -(1)
                        } else if c_int::from(*sptr1.offset(0)) > c_int::from(*sptr2.offset(0)) {
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
                    strcpy(
                        ((lParse.Nodes[this_node_idx]).value.data.astr).as_mut_ptr(),
                        sptr1,
                    );
                    strcat(
                        ((lParse.Nodes[this_node_idx]).value.data.astr).as_mut_ptr(),
                        sptr2,
                    );
                }
                _ => {}
            }
            (lParse.Nodes[this_node_idx]).operation = CONST_OP;
        } else {
            Allocate_Ptrs(lParse, this_node_idx);
            if lParse.status == 0 {
                rows = lParse.nRows;
                match (lParse.Nodes[this_node_idx]).operation {
                    280 | 279 => loop {
                        let fresh42 = rows;
                        rows -= 1;
                        if fresh42 == 0 {
                            break;
                        }
                        if const1 == 0 {
                            null1 = *((lParse.Nodes[that1_idx]).value.undef).offset(rows as isize);
                        }
                        if const2 == 0 {
                            null2 = *((lParse.Nodes[that2_idx]).value.undef).offset(rows as isize);
                        }
                        *((lParse.Nodes[this_node_idx]).value.undef).offset(rows as isize) =
                            if c_int::from(null1) != 0 || c_int::from(null2) != 0 {
                                1
                            } else {
                                0
                            };
                        if *((lParse.Nodes[this_node_idx]).value.undef).offset(rows as isize) == 0 {
                            if const1 == 0 {
                                sptr1 = *((lParse.Nodes[that1_idx]).value.data.strptr)
                                    .offset(rows as isize);
                            }
                            if const2 == 0 {
                                sptr2 = *((lParse.Nodes[that2_idx]).value.data.strptr)
                                    .offset(rows as isize);
                            }
                            val = c_int::from(
                                (if c_int::from(*sptr1.offset(0)) < c_int::from(*sptr2.offset(0)) {
                                    -(1)
                                } else if c_int::from(*sptr1.offset(0))
                                    > c_int::from(*sptr2.offset(0))
                                {
                                    1
                                } else {
                                    strcmp(sptr1, sptr2)
                                }) == 0,
                            );
                            *((lParse.Nodes[this_node_idx]).value.data.logptr)
                                .offset(rows as isize) = (if (lParse.Nodes[this_node_idx]).operation
                                == fits_parser_yytokentype::EQ as c_int
                            {
                                val
                            } else {
                                c_int::from(val == 0)
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
                            null1 = *((lParse.Nodes[that1_idx]).value.undef).offset(rows as isize);
                        }
                        if const2 == 0 {
                            null2 = *((lParse.Nodes[that2_idx]).value.undef).offset(rows as isize);
                        }
                        *((lParse.Nodes[this_node_idx]).value.undef).offset(rows as isize) =
                            if c_int::from(null1) != 0 || c_int::from(null2) != 0 {
                                1
                            } else {
                                0
                            };
                        if *((lParse.Nodes[this_node_idx]).value.undef).offset(rows as isize) == 0 {
                            if const1 == 0 {
                                sptr1 = *((lParse.Nodes[that1_idx]).value.data.strptr)
                                    .offset(rows as isize);
                            }
                            if const2 == 0 {
                                sptr2 = *((lParse.Nodes[that2_idx]).value.data.strptr)
                                    .offset(rows as isize);
                            }
                            val = if c_int::from(*sptr1.offset(0)) < c_int::from(*sptr2.offset(0)) {
                                -(1)
                            } else if c_int::from(*sptr1.offset(0)) > c_int::from(*sptr2.offset(0))
                            {
                                1
                            } else {
                                strcmp(sptr1, sptr2)
                            };
                            *((lParse.Nodes[this_node_idx]).value.data.logptr)
                                .offset(rows as isize) = (if (lParse.Nodes[this_node_idx]).operation
                                == fits_parser_yytokentype::GT as c_int
                            {
                                c_int::from(val > 0)
                            } else {
                                c_int::from(val < 0)
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
                            null1 = *((lParse.Nodes[that1_idx]).value.undef).offset(rows as isize);
                        }
                        if const2 == 0 {
                            null2 = *((lParse.Nodes[that2_idx]).value.undef).offset(rows as isize);
                        }
                        *((lParse.Nodes[this_node_idx]).value.undef).offset(rows as isize) =
                            if c_int::from(null1) != 0 || c_int::from(null2) != 0 {
                                1
                            } else {
                                0
                            };
                        if *((lParse.Nodes[this_node_idx]).value.undef).offset(rows as isize) == 0 {
                            if const1 == 0 {
                                sptr1 = *((lParse.Nodes[that1_idx]).value.data.strptr)
                                    .offset(rows as isize);
                            }
                            if const2 == 0 {
                                sptr2 = *((lParse.Nodes[that2_idx]).value.data.strptr)
                                    .offset(rows as isize);
                            }
                            val = if c_int::from(*sptr1.offset(0)) < c_int::from(*sptr2.offset(0)) {
                                -(1)
                            } else if c_int::from(*sptr1.offset(0)) > c_int::from(*sptr2.offset(0))
                            {
                                1
                            } else {
                                strcmp(sptr1, sptr2)
                            };
                            *((lParse.Nodes[this_node_idx]).value.data.logptr)
                                .offset(rows as isize) = (if (lParse.Nodes[this_node_idx]).operation
                                == fits_parser_yytokentype::GTE as c_int
                            {
                                c_int::from(val >= 0)
                            } else {
                                c_int::from(val <= 0)
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
                            null1 = *((lParse.Nodes[that1_idx]).value.undef).offset(rows as isize);
                        }
                        if const2 == 0 {
                            null2 = *((lParse.Nodes[that2_idx]).value.undef).offset(rows as isize);
                        }
                        *((lParse.Nodes[this_node_idx]).value.undef).offset(rows as isize) =
                            if c_int::from(null1) != 0 || c_int::from(null2) != 0 {
                                1
                            } else {
                                0
                            };
                        if *((lParse.Nodes[this_node_idx]).value.undef).offset(rows as isize) == 0 {
                            if const1 == 0 {
                                sptr1 = *((lParse.Nodes[that1_idx]).value.data.strptr)
                                    .offset(rows as isize);
                            }
                            if const2 == 0 {
                                sptr2 = *((lParse.Nodes[that2_idx]).value.data.strptr)
                                    .offset(rows as isize);
                            }
                            strcpy(
                                *((lParse.Nodes[this_node_idx]).value.data.strptr)
                                    .offset(rows as isize),
                                sptr1,
                            );
                            strcat(
                                *((lParse.Nodes[this_node_idx]).value.data.strptr)
                                    .offset(rows as isize),
                                sptr2,
                            );
                        }
                    },
                    _ => {}
                }
            }
        }
        if (lParse.Nodes[that1_idx]).operation > 0 {
            free((*((lParse.Nodes[that1_idx]).value.data.strptr).offset(0)).cast::<c_void>());
            free((lParse.Nodes[that1_idx]).value.data.strptr.cast::<c_void>());
        }
        if (lParse.Nodes[that2_idx]).operation > 0 {
            free((*((lParse.Nodes[that2_idx]).value.data.strptr).offset(0)).cast::<c_void>());
            free((lParse.Nodes[that2_idx]).value.data.strptr.cast::<c_void>());
        }
    }
}
fn Do_BinOp_log(lParse: &mut ParseData, this_node_idx: usize) {
    unsafe {
        let mut that1: &mut Node;
        let mut that2: &mut Node;
        let mut vector1: c_int = 0;
        let mut vector2: c_int = 0;
        let mut val1: c_char = 0;
        let mut val2: c_char = 0;
        let mut null1: c_char = 0;
        let mut null2: c_char = 0;
        let mut rows: c_long = 0;
        let mut nelem: c_long = 0;
        let mut elem: c_long = 0;
        let that1_idx = (lParse.Nodes[this_node_idx]).SubNodes[0];
        let that2_idx = (lParse.Nodes[this_node_idx]).SubNodes[1];

        vector1 = c_int::from((lParse.Nodes[that1_idx]).operation != CONST_OP);
        if vector1 != 0 {
            vector1 = (lParse.Nodes[that1_idx]).value.nelem as c_int;
        } else {
            val1 = (lParse.Nodes[that1_idx]).value.data.log;
        }
        vector2 = c_int::from((lParse.Nodes[that2_idx]).operation != CONST_OP);
        if vector2 != 0 {
            vector2 = (lParse.Nodes[that2_idx]).value.nelem as c_int;
        } else {
            val2 = (lParse.Nodes[that2_idx]).value.data.log;
        }
        if vector1 == 0 && vector2 == 0 {
            match (lParse.Nodes[this_node_idx]).operation {
                277 => {
                    (lParse.Nodes[this_node_idx]).value.data.log =
                        if c_int::from(val1) != 0 || c_int::from(val2) != 0 {
                            1
                        } else {
                            0
                        };
                }
                278 => {
                    (lParse.Nodes[this_node_idx]).value.data.log =
                        if c_int::from(val1) != 0 && c_int::from(val2) != 0 {
                            1
                        } else {
                            0
                        };
                }
                279 => {
                    (lParse.Nodes[this_node_idx]).value.data.log = c_int::from(
                        c_int::from(val1) != 0 && c_int::from(val2) != 0 || val1 == 0 && val2 == 0,
                    ) as c_char;
                }
                280 => {
                    (lParse.Nodes[this_node_idx]).value.data.log = if c_int::from(val1) != 0
                        && val2 == 0
                        || val1 == 0 && c_int::from(val2) != 0
                    {
                        1
                    } else {
                        0
                    };
                }
                291 => {
                    (lParse.Nodes[this_node_idx]).value.data.lng = c_long::from(val1);
                }
                _ => {}
            }
            (lParse.Nodes[this_node_idx]).operation = CONST_OP;
        } else if (lParse.Nodes[this_node_idx]).operation == fits_parser_yytokentype::ACCUM as c_int
        {
            let mut i: c_long = 0;
            let mut previous: c_long = 0;
            let mut curr: c_long = 0;
            rows = lParse.nRows;
            nelem = (lParse.Nodes[this_node_idx]).value.nelem;
            elem = (lParse.Nodes[this_node_idx]).value.nelem * rows;
            Allocate_Ptrs(lParse, this_node_idx);
            if lParse.status == 0 {
                previous = (lParse.Nodes[that2_idx]).value.data.lng;
                i = 0;
                while i < elem {
                    if *((lParse.Nodes[that1_idx]).value.undef).offset(i as isize) == 0 {
                        curr = c_long::from(
                            *((lParse.Nodes[that1_idx]).value.data.logptr).offset(i as isize),
                        );
                        previous += curr;
                    }
                    *((lParse.Nodes[this_node_idx]).value.data.lngptr).offset(i as isize) =
                        previous;
                    *((lParse.Nodes[this_node_idx]).value.undef).offset(i as isize) = 0;
                    i += 1;
                }
                (lParse.Nodes[that2_idx]).value.data.lng = previous;
            }
        } else {
            rows = lParse.nRows;
            nelem = (lParse.Nodes[this_node_idx]).value.nelem;
            elem = (lParse.Nodes[this_node_idx]).value.nelem * rows;
            Allocate_Ptrs(lParse, this_node_idx);
            if lParse.status == 0 {
                if (lParse.Nodes[this_node_idx]).operation
                    == fits_parser_yytokentype::ACCUM as c_int
                {
                    let mut i_0: c_long = 0;
                    let mut previous_0: c_long = 0;
                    let mut curr_0: c_long = 0;
                    previous_0 = (lParse.Nodes[that2_idx]).value.data.lng;
                    i_0 = 0;
                    while i_0 < elem {
                        if *((lParse.Nodes[that1_idx]).value.undef).offset(i_0 as isize) == 0 {
                            curr_0 = c_long::from(
                                *((lParse.Nodes[that1_idx]).value.data.logptr).offset(i_0 as isize),
                            );
                            previous_0 += curr_0;
                        }
                        *((lParse.Nodes[this_node_idx]).value.data.lngptr).offset(i_0 as isize) =
                            previous_0;
                        *((lParse.Nodes[this_node_idx]).value.undef).offset(i_0 as isize) = 0;
                        i_0 += 1;
                    }
                    (lParse.Nodes[that2_idx]).value.data.lng = previous_0;
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
                        if vector1 > 1 {
                            val1 = *((lParse.Nodes[that1_idx]).value.data.logptr)
                                .offset(elem as isize);
                            null1 = *((lParse.Nodes[that1_idx]).value.undef).offset(elem as isize);
                        } else if vector1 != 0 {
                            val1 = *((lParse.Nodes[that1_idx]).value.data.logptr)
                                .offset(rows as isize);
                            null1 = *((lParse.Nodes[that1_idx]).value.undef).offset(rows as isize);
                        }
                        if vector2 > 1 {
                            val2 = *((lParse.Nodes[that2_idx]).value.data.logptr)
                                .offset(elem as isize);
                            null2 = *((lParse.Nodes[that2_idx]).value.undef).offset(elem as isize);
                        } else if vector2 != 0 {
                            val2 = *((lParse.Nodes[that2_idx]).value.data.logptr)
                                .offset(rows as isize);
                            null2 = *((lParse.Nodes[that2_idx]).value.undef).offset(rows as isize);
                        }
                        *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize) =
                            if c_int::from(null1) != 0 || c_int::from(null2) != 0 {
                                1
                            } else {
                                0
                            };
                        match (lParse.Nodes[this_node_idx]).operation {
                            277 => {
                                if null1 == 0 && null2 == 0 {
                                    *((lParse.Nodes[this_node_idx]).value.data.logptr)
                                        .offset(elem as isize) =
                                        if c_int::from(val1) != 0 || c_int::from(val2) != 0 {
                                            1
                                        } else {
                                            0
                                        };
                                } else if c_int::from(null1) != 0
                                    && null2 == 0
                                    && c_int::from(val2) != 0
                                    || null1 == 0
                                        && c_int::from(null2) != 0
                                        && c_int::from(val1) != 0
                                {
                                    *((lParse.Nodes[this_node_idx]).value.data.logptr)
                                        .offset(elem as isize) = 1;
                                    *((lParse.Nodes[this_node_idx]).value.undef)
                                        .offset(elem as isize) = 0;
                                }
                            }
                            278 => {
                                if null1 == 0 && null2 == 0 {
                                    *((lParse.Nodes[this_node_idx]).value.data.logptr)
                                        .offset(elem as isize) =
                                        if c_int::from(val1) != 0 && c_int::from(val2) != 0 {
                                            1
                                        } else {
                                            0
                                        };
                                } else if c_int::from(null1) != 0 && null2 == 0 && val2 == 0
                                    || null1 == 0 && c_int::from(null2) != 0 && val1 == 0
                                {
                                    *((lParse.Nodes[this_node_idx]).value.data.logptr)
                                        .offset(elem as isize) = 0;
                                    *((lParse.Nodes[this_node_idx]).value.undef)
                                        .offset(elem as isize) = 0;
                                }
                            }
                            279 => {
                                *((lParse.Nodes[this_node_idx]).value.data.logptr)
                                    .offset(elem as isize) = c_int::from(
                                    c_int::from(val1) != 0 && c_int::from(val2) != 0
                                        || val1 == 0 && val2 == 0,
                                )
                                    as c_char;
                            }
                            280 => {
                                *((lParse.Nodes[this_node_idx]).value.data.logptr)
                                    .offset(elem as isize) = c_int::from(
                                    c_int::from(val1) != 0 && val2 == 0
                                        || val1 == 0 && c_int::from(val2) != 0,
                                )
                                    as c_char;
                            }
                            _ => {}
                        }
                    }
                    nelem = (lParse.Nodes[this_node_idx]).value.nelem;
                }
            }
        }
        if (lParse.Nodes[that1_idx]).operation > 0 {
            free((lParse.Nodes[that1_idx]).value.data.ptr);
        }
        if (lParse.Nodes[that2_idx]).operation > 0 {
            free((lParse.Nodes[that2_idx]).value.data.ptr);
        }
    }
}

fn Do_BinOp_lng(lParse: &mut ParseData, this_node_idx: usize) {
    unsafe {
        let mut vector1: c_int = 0;
        let mut vector2: c_int = 0;
        let mut val1: c_long = 0;
        let mut val2: c_long = 0;
        let mut null1: c_char = 0;
        let mut null2: c_char = 0;
        let mut rows: c_long = 0;
        let mut nelem: c_long = 0;
        let mut elem: c_long = 0;

        let that1_idx = (lParse.Nodes[this_node_idx]).SubNodes[0];
        let that2_idx = (lParse.Nodes[this_node_idx]).SubNodes[1];
        vector1 = c_int::from((lParse.Nodes[that1_idx]).operation != CONST_OP);

        if vector1 != 0 {
            vector1 = (lParse.Nodes[that1_idx]).value.nelem as c_int;
        } else {
            val1 = (lParse.Nodes[that1_idx]).value.data.lng;
        }

        vector2 = c_int::from((lParse.Nodes[that2_idx]).operation != CONST_OP);
        if vector2 != 0 {
            vector2 = (lParse.Nodes[that2_idx]).value.nelem as c_int;
        } else {
            val2 = (lParse.Nodes[that2_idx]).value.data.lng;
        }

        if vector1 == 0 && vector2 == 0 {
            /*  Result is a constant  */

            match (lParse.Nodes[this_node_idx]).operation {
                126 | 279 => {
                    (lParse.Nodes[this_node_idx]).value.data.log = if val1 == val2 { 1 } else { 0 };
                }
                280 => {
                    (lParse.Nodes[this_node_idx]).value.data.log = if val1 != val2 { 1 } else { 0 };
                }
                281 => {
                    (lParse.Nodes[this_node_idx]).value.data.log = if val1 > val2 { 1 } else { 0 };
                }
                282 => {
                    (lParse.Nodes[this_node_idx]).value.data.log = if val1 < val2 { 1 } else { 0 };
                }
                283 => {
                    (lParse.Nodes[this_node_idx]).value.data.log = if val1 <= val2 { 1 } else { 0 };
                }
                284 => {
                    (lParse.Nodes[this_node_idx]).value.data.log = if val1 >= val2 { 1 } else { 0 };
                }
                43 => {
                    (lParse.Nodes[this_node_idx]).value.data.lng = val1 + val2;
                }
                45 => {
                    (lParse.Nodes[this_node_idx]).value.data.lng = val1 - val2;
                }
                42 => {
                    (lParse.Nodes[this_node_idx]).value.data.lng = val1 * val2;
                }
                38 => {
                    (lParse.Nodes[this_node_idx]).value.data.lng = val1 & val2;
                }
                124 => {
                    (lParse.Nodes[this_node_idx]).value.data.lng = val1 | val2;
                }
                94 => {
                    (lParse.Nodes[this_node_idx]).value.data.lng = val1 ^ val2;
                }
                37 => {
                    if val2 != 0 {
                        if val1 == LONG_MIN && val2 == -1 {
                            *(lParse.Nodes[this_node_idx])
                                .value
                                .data
                                .lngptr
                                .offset(elem as isize) = 0;
                        } else {
                            (lParse.Nodes[this_node_idx]).value.data.lng = val1 % val2;
                        }
                    } else {
                        fits_parser_yyerror(lParse, cs!(c"Divide by Zero"));
                    }
                }
                47 => {
                    if val2 != 0 {
                        if val1 == LONG_MIN && val2 == -1 {
                            *(lParse.Nodes[this_node_idx])
                                .value
                                .data
                                .lngptr
                                .offset(elem as isize) = LONG_MAX;
                        } else {
                            (lParse.Nodes[this_node_idx]).value.data.lng = val1 / val2;
                        }
                    } else {
                        fits_parser_yyerror(lParse, cs!(c"Divide by Zero"));
                    }
                }
                286 => {
                    (lParse.Nodes[this_node_idx]).value.data.lng =
                        pow(val1 as c_double, val2 as c_double) as c_long;
                }
                291 => {
                    (lParse.Nodes[this_node_idx]).value.data.lng = val1;
                }
                292 => {
                    (lParse.Nodes[this_node_idx]).value.data.lng = 0;
                }
                _ => {}
            }
            (lParse.Nodes[this_node_idx]).operation = CONST_OP;
        } else if (lParse.Nodes[this_node_idx]).operation == fits_parser_yytokentype::ACCUM as c_int
            || (lParse.Nodes[this_node_idx]).operation == fits_parser_yytokentype::DIFF as c_int
        {
            let mut i: c_long = 0;
            let mut previous: c_long = 0;
            let mut curr: c_long = 0;
            let mut undef: c_long = 0;
            rows = lParse.nRows;
            nelem = (lParse.Nodes[this_node_idx]).value.nelem;
            elem = (lParse.Nodes[this_node_idx]).value.nelem * rows;
            Allocate_Ptrs(lParse, this_node_idx);
            if lParse.status == 0 {
                previous = (lParse.Nodes[that2_idx]).value.data.lng;
                /* the C stores this flag *in* the undef pointer field of the
                constant node ("XXX evil, but no harm here"), so read the
                pointer value back rather than dereferencing it */
                undef = (lParse.Nodes[that2_idx]).value.undef as c_long;
                if (lParse.Nodes[this_node_idx]).operation
                    == fits_parser_yytokentype::ACCUM as c_int
                {
                    i = 0;
                    while i < elem {
                        if *((lParse.Nodes[that1_idx]).value.undef).offset(i as isize) == 0 {
                            curr =
                                *((lParse.Nodes[that1_idx]).value.data.lngptr).offset(i as isize);
                            previous += curr;
                        }
                        *((lParse.Nodes[this_node_idx]).value.data.lngptr).offset(i as isize) =
                            previous;
                        *((lParse.Nodes[this_node_idx]).value.undef).offset(i as isize) = 0;
                        i += 1;
                    }
                } else {
                    i = 0;
                    while i < elem {
                        curr = *((lParse.Nodes[that1_idx]).value.data.lngptr).offset(i as isize);
                        if c_int::from(*((lParse.Nodes[that1_idx]).value.undef).offset(i as isize))
                            != 0
                            || undef != 0
                        {
                            *((lParse.Nodes[this_node_idx]).value.data.lngptr).offset(i as isize) =
                                0;
                            *((lParse.Nodes[this_node_idx]).value.undef).offset(i as isize) = 1;
                        } else {
                            *((lParse.Nodes[this_node_idx]).value.data.lngptr).offset(i as isize) =
                                curr - previous;
                            *((lParse.Nodes[this_node_idx]).value.undef).offset(i as isize) = 0;
                        }
                        previous = curr;
                        undef = c_long::from(
                            *((lParse.Nodes[that1_idx]).value.undef).offset(i as isize),
                        );
                        i += 1;
                    }
                }
                (lParse.Nodes[that2_idx]).value.data.lng = previous;
                (lParse.Nodes[that2_idx]).value.undef = undef as *mut c_char;
            }
        } else {
            rows = lParse.nRows;
            nelem = (lParse.Nodes[this_node_idx]).value.nelem;
            elem = (lParse.Nodes[this_node_idx]).value.nelem * rows;
            Allocate_Ptrs(lParse, this_node_idx);
            loop {
                let fresh48 = rows;
                rows -= 1;
                if !(fresh48 != 0 && lParse.status == 0) {
                    break;
                }
                loop {
                    let fresh49 = nelem;
                    nelem -= 1;
                    if !(fresh49 != 0 && lParse.status == 0) {
                        break;
                    }
                    elem -= 1;

                    if vector1 > 1 {
                        val1 = *((lParse.Nodes[that1_idx]).value.data.lngptr).offset(elem as isize);
                        null1 = *((lParse.Nodes[that1_idx]).value.undef).offset(elem as isize);
                    } else if vector1 != 0 {
                        val1 = *((lParse.Nodes[that1_idx]).value.data.lngptr).offset(rows as isize);
                        null1 = *((lParse.Nodes[that1_idx]).value.undef).offset(rows as isize);
                    }
                    if vector2 > 1 {
                        val2 = *((lParse.Nodes[that2_idx]).value.data.lngptr).offset(elem as isize);
                        null2 = *((lParse.Nodes[that2_idx]).value.undef).offset(elem as isize);
                    } else if vector2 != 0 {
                        val2 = *((lParse.Nodes[that2_idx]).value.data.lngptr).offset(rows as isize);
                        null2 = *((lParse.Nodes[that2_idx]).value.undef).offset(rows as isize);
                    }
                    *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize) =
                        if c_int::from(null1) != 0 || c_int::from(null2) != 0 {
                            1
                        } else {
                            0
                        };
                    match (lParse.Nodes[this_node_idx]).operation {
                        126 | 279 => {
                            *((lParse.Nodes[this_node_idx]).value.data.logptr)
                                .offset(elem as isize) = if val1 == val2 { 1 } else { 0 };
                        }
                        280 => {
                            *((lParse.Nodes[this_node_idx]).value.data.logptr)
                                .offset(elem as isize) = if val1 != val2 { 1 } else { 0 };
                        }
                        281 => {
                            *((lParse.Nodes[this_node_idx]).value.data.logptr)
                                .offset(elem as isize) = if val1 > val2 { 1 } else { 0 };
                        }
                        282 => {
                            *((lParse.Nodes[this_node_idx]).value.data.logptr)
                                .offset(elem as isize) = if val1 < val2 { 1 } else { 0 };
                        }
                        283 => {
                            *((lParse.Nodes[this_node_idx]).value.data.logptr)
                                .offset(elem as isize) = if val1 <= val2 { 1 } else { 0 };
                        }
                        284 => {
                            *((lParse.Nodes[this_node_idx]).value.data.logptr)
                                .offset(elem as isize) = if val1 >= val2 { 1 } else { 0 };
                        }
                        43 => {
                            *((lParse.Nodes[this_node_idx]).value.data.lngptr)
                                .offset(elem as isize) = val1 + val2;
                        }
                        45 => {
                            *((lParse.Nodes[this_node_idx]).value.data.lngptr)
                                .offset(elem as isize) = val1 - val2;
                        }
                        42 => {
                            *((lParse.Nodes[this_node_idx]).value.data.lngptr)
                                .offset(elem as isize) = val1 * val2;
                        }
                        38 => {
                            *((lParse.Nodes[this_node_idx]).value.data.lngptr)
                                .offset(elem as isize) = val1 & val2;
                        }
                        124 => {
                            *((lParse.Nodes[this_node_idx]).value.data.lngptr)
                                .offset(elem as isize) = val1 | val2;
                        }
                        94 => {
                            *((lParse.Nodes[this_node_idx]).value.data.lngptr)
                                .offset(elem as isize) = val1 ^ val2;
                        }
                        37 => {
                            if val2 != 0 {
                                *((lParse.Nodes[this_node_idx]).value.data.lngptr)
                                    .offset(elem as isize) = val1 % val2;
                            } else {
                                *((lParse.Nodes[this_node_idx]).value.data.lngptr)
                                    .offset(elem as isize) = 0;
                                *((lParse.Nodes[this_node_idx]).value.undef)
                                    .offset(elem as isize) = 1;
                            }
                        }
                        47 => {
                            if val2 != 0 {
                                *((lParse.Nodes[this_node_idx]).value.data.lngptr)
                                    .offset(elem as isize) = val1 / val2;
                            } else {
                                *((lParse.Nodes[this_node_idx]).value.data.lngptr)
                                    .offset(elem as isize) = 0;
                                *((lParse.Nodes[this_node_idx]).value.undef)
                                    .offset(elem as isize) = 1;
                            }
                        }
                        286 => {
                            *((lParse.Nodes[this_node_idx]).value.data.lngptr)
                                .offset(elem as isize) =
                                pow(val1 as c_double, val2 as c_double) as c_long;
                        }
                        _ => {}
                    }
                }
                nelem = (lParse.Nodes[this_node_idx]).value.nelem;
            }
        }
        if (lParse.Nodes[that1_idx]).operation > 0 {
            free((lParse.Nodes[that1_idx]).value.data.ptr);
        }
        if (lParse.Nodes[that2_idx]).operation > 0 {
            free((lParse.Nodes[that2_idx]).value.data.ptr);
        }
    }
}

fn validate_double_vector(lParse: &mut ParseData, node_idx: usize) -> c_int {
    let data = unsafe { (lParse.Nodes[node_idx]).value.data.dblptr };
    let undef = unsafe { (lParse.Nodes[node_idx]).value.undef };

    if data.is_null()
        || (data.addr()) < PARSER_VECTOR_MIN_ADDR
        || undef.is_null()
        || (undef.addr()) < PARSER_VECTOR_MIN_ADDR
    {
        fits_parser_yyerror(lParse, cs!(c"parser column data unavailable"));
        if lParse.status == 0 {
            lParse.status = PARSE_SYNTAX_ERR;
        }
        return 0;
    }

    return 1;
}

fn Do_BinOp_dbl(lParse: &mut ParseData, this_node_idx: usize) {
    unsafe {
        let mut that1: &mut Node;
        let mut that2: &mut Node;
        let mut vector1: c_int = 0;
        let mut vector2: c_int = 0;
        let mut val1: c_double = 0.0;
        let mut val2: c_double = 0.0;
        let mut null1: c_char = 0;
        let mut null2: c_char = 0;
        let mut rows: c_long = 0;
        let mut nelem: c_long = 0;
        let mut elem: c_long = 0;

        let that1_idx = (lParse.Nodes[this_node_idx]).SubNodes[0];
        let that2_idx = (lParse.Nodes[this_node_idx]).SubNodes[1];

        vector1 = c_int::from((lParse.Nodes[that1_idx]).operation != CONST_OP);
        if vector1 != 0 {
            vector1 = (lParse.Nodes[that1_idx]).value.nelem as c_int;
        } else {
            val1 = (lParse.Nodes[that1_idx]).value.data.dbl;
        }

        vector2 = c_int::from((lParse.Nodes[that2_idx]).operation != CONST_OP);
        if vector2 != 0 {
            vector2 = (lParse.Nodes[that2_idx]).value.nelem as c_int;
        } else {
            val2 = (lParse.Nodes[that2_idx]).value.data.dbl;
        }

        if vector1 != 0 && validate_double_vector(lParse, that1_idx) == 0 {
            return;
        }

        if vector2 != 0 && validate_double_vector(lParse, that2_idx) == 0 {
            return;
        }

        if vector1 == 0 && vector2 == 0 {
            /*  Result is a constant  */
            match (lParse.Nodes[this_node_idx]).operation {
                126 => {
                    (lParse.Nodes[this_node_idx]).value.data.log =
                        if fabs(val1 - val2) < APPROX { 1 } else { 0 };
                }
                279 => {
                    (lParse.Nodes[this_node_idx]).value.data.log = if val1 == val2 { 1 } else { 0 };
                }
                280 => {
                    (lParse.Nodes[this_node_idx]).value.data.log = if val1 != val2 { 1 } else { 0 };
                }
                281 => {
                    (lParse.Nodes[this_node_idx]).value.data.log = if val1 > val2 { 1 } else { 0 };
                }
                282 => {
                    (lParse.Nodes[this_node_idx]).value.data.log = if val1 < val2 { 1 } else { 0 };
                }
                283 => {
                    (lParse.Nodes[this_node_idx]).value.data.log = if val1 <= val2 { 1 } else { 0 };
                }
                284 => {
                    (lParse.Nodes[this_node_idx]).value.data.log = if val1 >= val2 { 1 } else { 0 };
                }
                43 => {
                    (lParse.Nodes[this_node_idx]).value.data.dbl = val1 + val2;
                }
                45 => {
                    (lParse.Nodes[this_node_idx]).value.data.dbl = val1 - val2;
                }
                42 => {
                    (lParse.Nodes[this_node_idx]).value.data.dbl = val1 * val2;
                }
                37 => {
                    if val2 != 0. {
                        (lParse.Nodes[this_node_idx]).value.data.dbl =
                            val1 - val2 * c_double::from((val1 / val2) as c_int);
                    } else {
                        fits_parser_yyerror(lParse, cs!(c"Divide by Zero"));
                    }
                }
                47 => {
                    if val2 != 0. {
                        (lParse.Nodes[this_node_idx]).value.data.dbl = val1 / val2;
                    } else {
                        fits_parser_yyerror(lParse, cs!(c"Divide by Zero"));
                    }
                }
                286 => {
                    (lParse.Nodes[this_node_idx]).value.data.dbl = pow(val1, val2);
                }
                291 => {
                    (lParse.Nodes[this_node_idx]).value.data.dbl = val1;
                }
                292 => {
                    (lParse.Nodes[this_node_idx]).value.data.dbl = 0.0;
                }
                _ => {}
            }
            (lParse.Nodes[this_node_idx]).operation = CONST_OP;
        } else if (lParse.Nodes[this_node_idx]).operation == fits_parser_yytokentype::ACCUM as c_int
            || (lParse.Nodes[this_node_idx]).operation == fits_parser_yytokentype::DIFF as c_int
        {
            let mut i: c_long = 0;
            let mut undef: c_long = 0;
            let mut previous: c_double = 0.0;
            let mut curr: c_double = 0.0;

            rows = lParse.nRows;
            nelem = (lParse.Nodes[this_node_idx]).value.nelem;
            elem = (lParse.Nodes[this_node_idx]).value.nelem * rows;

            Allocate_Ptrs(lParse, this_node_idx);

            if lParse.status == 0 {
                previous = (lParse.Nodes[that2_idx]).value.data.dbl;
                undef = (lParse.Nodes[that2_idx]).value.undef as c_long;
                if (lParse.Nodes[this_node_idx]).operation
                    == fits_parser_yytokentype::ACCUM as c_int
                {
                    /* Cumulative sum of this chunk */
                    for i in 0..elem {
                        if *((lParse.Nodes[that1_idx]).value.undef).offset(i as isize) == 0 {
                            curr =
                                *((lParse.Nodes[that1_idx]).value.data.dblptr).offset(i as isize);
                            previous += curr;
                        }
                        *((lParse.Nodes[this_node_idx]).value.data.dblptr).offset(i as isize) =
                            previous;
                        *((lParse.Nodes[this_node_idx]).value.undef).offset(i as isize) = 0;
                    }
                } else {
                    i = 0;
                    while i < elem {
                        curr = *((lParse.Nodes[that1_idx]).value.data.dblptr).offset(i as isize);
                        if c_int::from(*((lParse.Nodes[that1_idx]).value.undef).offset(i as isize))
                            != 0
                            || undef != 0
                        {
                            *((lParse.Nodes[this_node_idx]).value.data.dblptr).offset(i as isize) =
                                0.0;
                            *((lParse.Nodes[this_node_idx]).value.undef).offset(i as isize) = 1;
                        } else {
                            *((lParse.Nodes[this_node_idx]).value.data.dblptr).offset(i as isize) =
                                curr - previous;
                            *((lParse.Nodes[this_node_idx]).value.undef).offset(i as isize) = 0;
                        }
                        previous = curr;
                        undef = c_long::from(
                            *((lParse.Nodes[that1_idx]).value.undef).offset(i as isize),
                        );
                        i += 1;
                    }
                }
                (lParse.Nodes[that2_idx]).value.data.dbl = previous;
                (lParse.Nodes[that2_idx]).value.undef = undef as *mut c_char;
            }
        } else {
            rows = lParse.nRows;
            nelem = (lParse.Nodes[this_node_idx]).value.nelem;
            elem = (lParse.Nodes[this_node_idx]).value.nelem * rows;
            Allocate_Ptrs(lParse, this_node_idx);
            loop {
                let fresh50 = rows;
                rows -= 1;
                if !(fresh50 != 0 && lParse.status == 0) {
                    break;
                }
                loop {
                    let fresh51 = nelem;
                    nelem -= 1;
                    if !(fresh51 != 0 && lParse.status == 0) {
                        break;
                    }
                    elem -= 1;

                    if vector1 > 1 {
                        val1 = *((lParse.Nodes[that1_idx]).value.data.dblptr).offset(elem as isize);
                        null1 = *((lParse.Nodes[that1_idx]).value.undef).offset(elem as isize);
                    } else if vector1 != 0 {
                        val1 = *((lParse.Nodes[that1_idx]).value.data.dblptr).offset(rows as isize);
                        null1 = *((lParse.Nodes[that1_idx]).value.undef).offset(rows as isize);
                    }
                    if vector2 > 1 {
                        val2 = *((lParse.Nodes[that2_idx]).value.data.dblptr).offset(elem as isize);
                        null2 = *((lParse.Nodes[that2_idx]).value.undef).offset(elem as isize);
                    } else if vector2 != 0 {
                        val2 = *((lParse.Nodes[that2_idx]).value.data.dblptr).offset(rows as isize);
                        null2 = *((lParse.Nodes[that2_idx]).value.undef).offset(rows as isize);
                    }
                    *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize) =
                        if c_int::from(null1) != 0 || c_int::from(null2) != 0 {
                            1
                        } else {
                            0
                        };
                    match (lParse.Nodes[this_node_idx]).operation {
                        126 => {
                            *((lParse.Nodes[this_node_idx]).value.data.logptr)
                                .offset(elem as isize) =
                                if fabs(val1 - val2) < APPROX { 1 } else { 0 };
                        }
                        279 => {
                            *((lParse.Nodes[this_node_idx]).value.data.logptr)
                                .offset(elem as isize) = if val1 == val2 { 1 } else { 0 };
                        }
                        280 => {
                            *((lParse.Nodes[this_node_idx]).value.data.logptr)
                                .offset(elem as isize) = if val1 != val2 { 1 } else { 0 };
                        }
                        281 => {
                            *((lParse.Nodes[this_node_idx]).value.data.logptr)
                                .offset(elem as isize) = if val1 > val2 { 1 } else { 0 };
                        }
                        282 => {
                            *((lParse.Nodes[this_node_idx]).value.data.logptr)
                                .offset(elem as isize) = if val1 < val2 { 1 } else { 0 };
                        }
                        283 => {
                            *((lParse.Nodes[this_node_idx]).value.data.logptr)
                                .offset(elem as isize) = if val1 <= val2 { 1 } else { 0 };
                        }
                        284 => {
                            *((lParse.Nodes[this_node_idx]).value.data.logptr)
                                .offset(elem as isize) = if val1 >= val2 { 1 } else { 0 };
                        }
                        43 => {
                            *((lParse.Nodes[this_node_idx]).value.data.dblptr)
                                .offset(elem as isize) = val1 + val2;
                        }
                        45 => {
                            *((lParse.Nodes[this_node_idx]).value.data.dblptr)
                                .offset(elem as isize) = val1 - val2;
                        }
                        42 => {
                            *((lParse.Nodes[this_node_idx]).value.data.dblptr)
                                .offset(elem as isize) = val1 * val2;
                        }
                        37 => {
                            if val2 != 0. {
                                *((lParse.Nodes[this_node_idx]).value.data.dblptr)
                                    .offset(elem as isize) =
                                    val1 - val2 * c_double::from((val1 / val2) as c_int);
                            } else {
                                *((lParse.Nodes[this_node_idx]).value.data.dblptr)
                                    .offset(elem as isize) = 0.0;
                                *((lParse.Nodes[this_node_idx]).value.undef)
                                    .offset(elem as isize) = 1;
                            }
                        }
                        47 => {
                            if val2 != 0. {
                                *((lParse.Nodes[this_node_idx]).value.data.dblptr)
                                    .offset(elem as isize) = val1 / val2;
                            } else {
                                *((lParse.Nodes[this_node_idx]).value.data.dblptr)
                                    .offset(elem as isize) = 0.0;
                                *((lParse.Nodes[this_node_idx]).value.undef)
                                    .offset(elem as isize) = 1;
                            }
                        }
                        286 => {
                            *((lParse.Nodes[this_node_idx]).value.data.dblptr)
                                .offset(elem as isize) = pow(val1, val2);
                        }
                        _ => {}
                    }
                }
                nelem = (lParse.Nodes[this_node_idx]).value.nelem;
            }
        }
        if (lParse.Nodes[that1_idx]).operation > 0 {
            free((lParse.Nodes[that1_idx]).value.data.ptr);
        }
        if (lParse.Nodes[that2_idx]).operation > 0 {
            free((lParse.Nodes[that2_idx]).value.data.ptr);
        }
    }
}

pub(crate) fn qselect_median_lng(arr: *mut c_long, n: c_int) -> c_long {
    unsafe {
        let mut low: c_int = 0;
        let mut high: c_int = 0;
        let mut median: c_int = 0;
        let mut middle: c_int = 0;
        let mut ll: c_int = 0;
        let mut hh: c_int = 0;
        low = 0;
        high = n - 1;
        median = (low + high) / 2;
        loop {
            if high <= low {
                return *arr.offset(median as isize);
            }
            if high == low + 1 {
                if *arr.offset(low as isize) > *arr.offset(high as isize) {
                    let t: c_long = *arr.offset(low as isize);
                    *arr.offset(low as isize) = *arr.offset(high as isize);
                    *arr.offset(high as isize) = t;
                }
                return *arr.offset(median as isize);
            }
            middle = (low + high) / 2;
            if *arr.offset(middle as isize) > *arr.offset(high as isize) {
                let t_0: c_long = *arr.offset(middle as isize);
                *arr.offset(middle as isize) = *arr.offset(high as isize);
                *arr.offset(high as isize) = t_0;
            }
            if *arr.offset(low as isize) > *arr.offset(high as isize) {
                let t_1: c_long = *arr.offset(low as isize);
                *arr.offset(low as isize) = *arr.offset(high as isize);
                *arr.offset(high as isize) = t_1;
            }
            if *arr.offset(middle as isize) > *arr.offset(low as isize) {
                let t_2: c_long = *arr.offset(middle as isize);
                *arr.offset(middle as isize) = *arr.offset(low as isize);
                *arr.offset(low as isize) = t_2;
            }
            let t_3: c_long = *arr.offset(middle as isize);
            *arr.offset(middle as isize) = *arr.offset((low + 1) as isize);
            *arr.offset((low + 1) as isize) = t_3;
            ll = low + 1;
            hh = high;
            loop {
                loop {
                    ll += 1;
                    if *arr.offset(low as isize) <= *arr.offset(ll as isize) {
                        break;
                    }
                }
                loop {
                    hh -= 1;
                    if *arr.offset(hh as isize) <= *arr.offset(low as isize) {
                        break;
                    }
                }
                if hh < ll {
                    break;
                }
                let t_4: c_long = *arr.offset(ll as isize);
                *arr.offset(ll as isize) = *arr.offset(hh as isize);
                *arr.offset(hh as isize) = t_4;
            }
            let t_5: c_long = *arr.offset(low as isize);
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

pub(crate) fn qselect_median_dbl(arr: *mut c_double, n: c_int) -> c_double {
    unsafe {
        let mut low: c_int = 0;
        let mut high: c_int = 0;
        let mut median: c_int = 0;
        let mut middle: c_int = 0;
        let mut ll: c_int = 0;
        let mut hh: c_int = 0;
        low = 0;
        high = n - 1;
        median = (low + high) / 2;
        loop {
            if high <= low {
                return *arr.offset(median as isize);
            }
            if high == low + 1 {
                if *arr.offset(low as isize) > *arr.offset(high as isize) {
                    let t: c_double = *arr.offset(low as isize);
                    *arr.offset(low as isize) = *arr.offset(high as isize);
                    *arr.offset(high as isize) = t;
                }
                return *arr.offset(median as isize);
            }
            middle = (low + high) / 2;
            if *arr.offset(middle as isize) > *arr.offset(high as isize) {
                let t_0: c_double = *arr.offset(middle as isize);
                *arr.offset(middle as isize) = *arr.offset(high as isize);
                *arr.offset(high as isize) = t_0;
            }
            if *arr.offset(low as isize) > *arr.offset(high as isize) {
                let t_1: c_double = *arr.offset(low as isize);
                *arr.offset(low as isize) = *arr.offset(high as isize);
                *arr.offset(high as isize) = t_1;
            }
            if *arr.offset(middle as isize) > *arr.offset(low as isize) {
                let t_2: c_double = *arr.offset(middle as isize);
                *arr.offset(middle as isize) = *arr.offset(low as isize);
                *arr.offset(low as isize) = t_2;
            }
            let t_3: c_double = *arr.offset(middle as isize);
            *arr.offset(middle as isize) = *arr.offset((low + 1) as isize);
            *arr.offset((low + 1) as isize) = t_3;
            ll = low + 1;
            hh = high;
            loop {
                loop {
                    ll += 1;
                    if !(*arr.offset(low as isize) > *arr.offset(ll as isize)) {
                        break;
                    }
                }
                loop {
                    hh -= 1;
                    if !(*arr.offset(hh as isize) > *arr.offset(low as isize)) {
                        break;
                    }
                }
                if hh < ll {
                    break;
                }
                let t_4: c_double = *arr.offset(ll as isize);
                *arr.offset(ll as isize) = *arr.offset(hh as isize);
                *arr.offset(hh as isize) = t_4;
            }
            let t_5: c_double = *arr.offset(low as isize);
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

pub(crate) fn angsep_calc(
    ra1: c_double,
    dec1: c_double,
    ra2: c_double,
    dec2: c_double,
) -> c_double {
    let DEG: c_double = 4.0 * (1.0_f64).atan() / 180.0;
    let mut a: c_double = 0.0;
    let mut sdec: c_double = 0.0;
    let mut sra: c_double = 0.0;

    /* The algorithm is the law of Haversines.  This algorithm is
    stable even when the points are close together.  The normal
    Law of Cosines fails for angles around 0.1 arcsec. */

    sra = ((ra2 - ra1) * DEG / 2.0).sin();
    sdec = ((dec2 - dec1) * DEG / 2.0).sin();
    a = sdec * sdec + (dec1 * DEG).cos() * (dec2 * DEG).cos() * sra * sra;

    /* Sanity checking to avoid a range error in the sqrt()'s below */
    if a < 0.0 {
        a = 0.0;
    }
    if a > 1.0 {
        a = 1.0;
    }

    2.0 * atan2((a).sqrt(), (1.0 - a).sqrt()) / DEG
}

fn Do_Func(lParse: &mut ParseData, this_node_idx: usize) {
    unsafe {
        let mut theParams: [usize; 10] = [0; 10];
        let mut vector: [c_int; 10] = [0; 10];
        let mut allConst: c_int = 0;
        let mut pVals: [lval; 10] = [lval {
            nelem: 0,
            naxis: 0,
            naxes: [0; 5],
            undef: core::ptr::null_mut::<c_char>(),
            data: data_union { dbl: 0. },
        }; 10];
        let mut pNull: [c_char; 10] = [0; 10];
        let mut ival: c_long = 0;
        let mut dval: c_double = 0.0;
        let mut i: c_int = 0;
        let mut valInit: c_int = 0;
        let mut row: c_long = 0;
        let mut elem: c_long = 0;
        let mut nelem: c_long = 0;
        i = (lParse.Nodes[this_node_idx]).nSubNodes;
        allConst = 1;
        loop {
            let fresh52 = i;
            i -= 1;
            if fresh52 == 0 {
                break;
            }
            theParams[i as usize] = (lParse.Nodes[this_node_idx]).SubNodes[i as usize] as usize;
            vector[i as usize] =
                c_int::from((lParse.Nodes[theParams[i as usize]]).operation != CONST_OP);
            if vector[i as usize] != 0 {
                allConst = 0;
                vector[i as usize] = (lParse.Nodes[theParams[i as usize]]).value.nelem as c_int;
            } else {
                if (lParse.Nodes[theParams[i as usize]]).ntype
                    == fits_parser_yytokentype::DOUBLE as c_int
                {
                    pVals[i as usize].data.dbl =
                        (lParse.Nodes[theParams[i as usize]]).value.data.dbl;
                } else if (lParse.Nodes[theParams[i as usize]]).ntype
                    == fits_parser_yytokentype::LONG as c_int
                {
                    pVals[i as usize].data.lng =
                        (lParse.Nodes[theParams[i as usize]]).value.data.lng;
                } else if (lParse.Nodes[theParams[i as usize]]).ntype
                    == fits_parser_yytokentype::BOOLEAN as c_int
                {
                    pVals[i as usize].data.log =
                        (lParse.Nodes[theParams[i as usize]]).value.data.log;
                } else {
                    strcpy(
                        (pVals[i as usize].data.astr).as_mut_ptr(),
                        ((lParse.Nodes[theParams[i as usize]]).value.data.astr).as_mut_ptr(),
                    );
                }
                pNull[i as usize] = 0;
            }
        }
        if (lParse.Nodes[this_node_idx]).nSubNodes == 0 {
            allConst = 0;
        }
        if (lParse.Nodes[this_node_idx]).operation == POIRND_FCT as c_int {
            allConst = 0;
        }
        if (lParse.Nodes[this_node_idx]).operation == GASRND_FCT as c_int {
            allConst = 0;
        }
        if (lParse.Nodes[this_node_idx]).operation == RND_FCT as c_int {
            allConst = 0;
        }
        if allConst != 0 {
            let current_block_139: u64;
            match (lParse.Nodes[this_node_idx]).operation as u32 {
                1002 => {
                    if (lParse.Nodes[theParams[0]]).ntype
                        == fits_parser_yytokentype::BOOLEAN as c_int
                    {
                        (lParse.Nodes[this_node_idx]).value.data.lng =
                            c_long::from(if c_int::from(pVals[0].data.log) != 0 {
                                1
                            } else {
                                0
                            });
                    } else if (lParse.Nodes[theParams[0]]).ntype
                        == fits_parser_yytokentype::LONG as c_int
                    {
                        (lParse.Nodes[this_node_idx]).value.data.lng = pVals[0].data.lng;
                    } else if (lParse.Nodes[theParams[0]]).ntype
                        == fits_parser_yytokentype::DOUBLE as c_int
                    {
                        (lParse.Nodes[this_node_idx]).value.data.dbl = pVals[0].data.dbl;
                    } else if (lParse.Nodes[theParams[0]]).ntype
                        == fits_parser_yytokentype::BITSTR as c_int
                    {
                        strcpy(
                            ((lParse.Nodes[this_node_idx]).value.data.astr).as_mut_ptr(),
                            (pVals[0].data.astr).as_mut_ptr(),
                        );
                    }
                    current_block_139 = 7627602990488000394;
                }
                1038 => {
                    if (lParse.Nodes[theParams[0]]).ntype == fits_parser_yytokentype::LONG as c_int
                    {
                        (lParse.Nodes[this_node_idx]).value.data.dbl =
                            pVals[0].data.lng as c_double;
                    } else if (lParse.Nodes[theParams[0]]).ntype
                        == fits_parser_yytokentype::DOUBLE as c_int
                    {
                        (lParse.Nodes[this_node_idx]).value.data.dbl = pVals[0].data.dbl;
                    }
                    current_block_139 = 7627602990488000394;
                }
                1039 => {
                    (lParse.Nodes[this_node_idx]).value.data.dbl = 0.0;
                    current_block_139 = 7627602990488000394;
                }
                1037 => {
                    if (lParse.Nodes[theParams[0]]).ntype
                        == fits_parser_yytokentype::BOOLEAN as c_int
                    {
                        (lParse.Nodes[this_node_idx]).value.data.lng =
                            c_long::from(if c_int::from(pVals[0].data.log) != 0 {
                                1
                            } else {
                                0
                            });
                    } else if (lParse.Nodes[theParams[0]]).ntype
                        == fits_parser_yytokentype::LONG as c_int
                    {
                        (lParse.Nodes[this_node_idx]).value.data.lng = pVals[0].data.lng;
                    } else {
                        (lParse.Nodes[this_node_idx]).value.data.dbl = pVals[0].data.dbl;
                    }
                    current_block_139 = 7627602990488000394;
                }
                1043 => {
                    if (lParse.Nodes[theParams[0]]).ntype
                        == fits_parser_yytokentype::DOUBLE as c_int
                    {
                        (lParse.Nodes[this_node_idx]).value.data.lng =
                            c_long::from(simplerng_getpoisson(pVals[0].data.dbl));
                    } else {
                        (lParse.Nodes[this_node_idx]).value.data.lng =
                            c_long::from(simplerng_getpoisson(pVals[0].data.lng as c_double));
                    }
                    current_block_139 = 7627602990488000394;
                }
                1017 => {
                    if (lParse.Nodes[theParams[0]]).ntype
                        == fits_parser_yytokentype::DOUBLE as c_int
                    {
                        dval = pVals[0].data.dbl;
                        (lParse.Nodes[this_node_idx]).value.data.dbl =
                            if dval > 0.0 { dval } else { -dval };
                    } else {
                        ival = pVals[0].data.lng;
                        (lParse.Nodes[this_node_idx]).value.data.lng =
                            if ival > 0 { ival } else { -ival };
                    }
                    current_block_139 = 7627602990488000394;
                }
                1040 => {
                    (lParse.Nodes[this_node_idx]).value.data.lng = 1;
                    current_block_139 = 7627602990488000394;
                }
                1030 => {
                    (lParse.Nodes[this_node_idx]).value.data.log = 0;
                    current_block_139 = 7627602990488000394;
                }
                1031 => {
                    if (lParse.Nodes[this_node_idx]).ntype
                        == fits_parser_yytokentype::BOOLEAN as c_int
                    {
                        (lParse.Nodes[this_node_idx]).value.data.log = pVals[0].data.log;
                    } else if (lParse.Nodes[this_node_idx]).ntype
                        == fits_parser_yytokentype::LONG as c_int
                    {
                        (lParse.Nodes[this_node_idx]).value.data.lng = pVals[0].data.lng;
                    } else if (lParse.Nodes[this_node_idx]).ntype
                        == fits_parser_yytokentype::DOUBLE as c_int
                    {
                        (lParse.Nodes[this_node_idx]).value.data.dbl = pVals[0].data.dbl;
                    } else if (lParse.Nodes[this_node_idx]).ntype
                        == fits_parser_yytokentype::STRING as c_int
                    {
                        strcpy(
                            ((lParse.Nodes[this_node_idx]).value.data.astr).as_mut_ptr(),
                            (pVals[0].data.astr).as_mut_ptr(),
                        );
                    }
                    current_block_139 = 7627602990488000394;
                }
                1046 => {
                    if (lParse.Nodes[this_node_idx]).ntype == fits_parser_yytokentype::LONG as c_int
                    {
                        (lParse.Nodes[this_node_idx]).value.data.lng = pVals[0].data.lng;
                    } else if (lParse.Nodes[this_node_idx]).ntype
                        == fits_parser_yytokentype::DOUBLE as c_int
                    {
                        (lParse.Nodes[this_node_idx]).value.data.dbl = pVals[0].data.dbl;
                    }
                    current_block_139 = 7627602990488000394;
                }
                1004 => {
                    (lParse.Nodes[this_node_idx]).value.data.dbl = sin(pVals[0].data.dbl);
                    current_block_139 = 7627602990488000394;
                }
                1005 => {
                    (lParse.Nodes[this_node_idx]).value.data.dbl = cos(pVals[0].data.dbl);
                    current_block_139 = 7627602990488000394;
                }
                1006 => {
                    (lParse.Nodes[this_node_idx]).value.data.dbl = tan(pVals[0].data.dbl);
                    current_block_139 = 7627602990488000394;
                }
                1007 => {
                    dval = pVals[0].data.dbl;
                    if dval < -1.0 || dval > 1.0 {
                        fits_parser_yyerror(lParse, cs!(c"Out of range argument to arcsin"));
                    } else {
                        (lParse.Nodes[this_node_idx]).value.data.dbl = asin(dval);
                    }
                    current_block_139 = 7627602990488000394;
                }
                1008 => {
                    dval = pVals[0].data.dbl;
                    if dval < -1.0 || dval > 1.0 {
                        fits_parser_yyerror(lParse, cs!(c"Out of range argument to arccos"));
                    } else {
                        (lParse.Nodes[this_node_idx]).value.data.dbl = acos(dval);
                    }
                    current_block_139 = 7627602990488000394;
                }
                1009 => {
                    (lParse.Nodes[this_node_idx]).value.data.dbl = atan(pVals[0].data.dbl);
                    current_block_139 = 7627602990488000394;
                }
                1010 => {
                    (lParse.Nodes[this_node_idx]).value.data.dbl = sinh(pVals[0].data.dbl);
                    current_block_139 = 7627602990488000394;
                }
                1011 => {
                    (lParse.Nodes[this_node_idx]).value.data.dbl = cosh(pVals[0].data.dbl);
                    current_block_139 = 7627602990488000394;
                }
                1012 => {
                    (lParse.Nodes[this_node_idx]).value.data.dbl = tanh(pVals[0].data.dbl);
                    current_block_139 = 7627602990488000394;
                }
                1013 => {
                    (lParse.Nodes[this_node_idx]).value.data.dbl = exp(pVals[0].data.dbl);
                    current_block_139 = 7627602990488000394;
                }
                1014 => {
                    dval = pVals[0].data.dbl;
                    if dval <= 0.0 {
                        fits_parser_yyerror(lParse, cs!(c"Out of range argument to log"));
                    } else {
                        (lParse.Nodes[this_node_idx]).value.data.dbl = log(dval);
                    }
                    current_block_139 = 7627602990488000394;
                }
                1015 => {
                    dval = pVals[0].data.dbl;
                    if dval <= 0.0 {
                        fits_parser_yyerror(lParse, cs!(c"Out of range argument to log10"));
                    } else {
                        (lParse.Nodes[this_node_idx]).value.data.dbl = log10(dval);
                    }
                    current_block_139 = 7627602990488000394;
                }
                1016 => {
                    dval = pVals[0].data.dbl;
                    if dval < 0.0 {
                        fits_parser_yyerror(lParse, cs!(c"Out of range argument to sqrt"));
                    } else {
                        (lParse.Nodes[this_node_idx]).value.data.dbl = sqrt(dval);
                    }
                    current_block_139 = 7627602990488000394;
                }
                1019 => {
                    (lParse.Nodes[this_node_idx]).value.data.dbl = ceil(pVals[0].data.dbl);
                    current_block_139 = 7627602990488000394;
                }
                1020 => {
                    (lParse.Nodes[this_node_idx]).value.data.dbl = floor(pVals[0].data.dbl);
                    current_block_139 = 7627602990488000394;
                }
                1021 => {
                    (lParse.Nodes[this_node_idx]).value.data.dbl = floor(pVals[0].data.dbl + 0.5);
                    current_block_139 = 7627602990488000394;
                }
                1018 => {
                    (lParse.Nodes[this_node_idx]).value.data.dbl =
                        atan2(pVals[0].data.dbl, pVals[1].data.dbl);
                    current_block_139 = 7627602990488000394;
                }
                1041 => {
                    (lParse.Nodes[this_node_idx]).value.data.dbl = angsep_calc(
                        pVals[0].data.dbl,
                        pVals[1].data.dbl,
                        pVals[2].data.dbl,
                        pVals[3].data.dbl,
                    );
                    /* DEVIATION from CFITSIO 4.7.0: "case angsep_fct:" there
                    is missing its "break;" and falls through into
                    "case min1_fct:", which overwrites the separation just
                    computed with pVals[0] - so a constant-folded ANGSEP
                    returns its first argument.  Fix submitted upstream; the
                    non-constant path in Do_Func has always been correct. */
                    current_block_139 = 7627602990488000394;
                }
                1022 => {
                    current_block_139 = 15934000668868306918;
                }
                1023 => {
                    if (lParse.Nodes[this_node_idx]).ntype
                        == fits_parser_yytokentype::DOUBLE as c_int
                    {
                        (lParse.Nodes[this_node_idx]).value.data.dbl =
                            if pVals[0].data.dbl < pVals[1].data.dbl {
                                pVals[0].data.dbl
                            } else {
                                pVals[1].data.dbl
                            };
                    } else if (lParse.Nodes[this_node_idx]).ntype
                        == fits_parser_yytokentype::LONG as c_int
                    {
                        (lParse.Nodes[this_node_idx]).value.data.lng =
                            if pVals[0].data.lng < pVals[1].data.lng {
                                pVals[0].data.lng
                            } else {
                                pVals[1].data.lng
                            };
                    }
                    current_block_139 = 7627602990488000394;
                }
                1024 => {
                    if (lParse.Nodes[this_node_idx]).ntype
                        == fits_parser_yytokentype::DOUBLE as c_int
                    {
                        (lParse.Nodes[this_node_idx]).value.data.dbl = pVals[0].data.dbl;
                    } else if (lParse.Nodes[this_node_idx]).ntype
                        == fits_parser_yytokentype::LONG as c_int
                    {
                        (lParse.Nodes[this_node_idx]).value.data.lng = pVals[0].data.lng;
                    } else if (lParse.Nodes[this_node_idx]).ntype
                        == fits_parser_yytokentype::BITSTR as c_int
                    {
                        strcpy(
                            ((lParse.Nodes[this_node_idx]).value.data.astr).as_mut_ptr(),
                            (pVals[0].data.astr).as_mut_ptr(),
                        );
                    }
                    current_block_139 = 7627602990488000394;
                }
                1025 => {
                    if (lParse.Nodes[this_node_idx]).ntype
                        == fits_parser_yytokentype::DOUBLE as c_int
                    {
                        (lParse.Nodes[this_node_idx]).value.data.dbl =
                            if pVals[0].data.dbl > pVals[1].data.dbl {
                                pVals[0].data.dbl
                            } else {
                                pVals[1].data.dbl
                            };
                    } else if (lParse.Nodes[this_node_idx]).ntype
                        == fits_parser_yytokentype::LONG as c_int
                    {
                        (lParse.Nodes[this_node_idx]).value.data.lng =
                            if pVals[0].data.lng > pVals[1].data.lng {
                                pVals[0].data.lng
                            } else {
                                pVals[1].data.lng
                            };
                    }
                    current_block_139 = 7627602990488000394;
                }
                1026 => {
                    (lParse.Nodes[this_node_idx]).value.data.log =
                        bnear(pVals[0].data.dbl, pVals[1].data.dbl, pVals[2].data.dbl);
                    current_block_139 = 7627602990488000394;
                }
                1027 => {
                    (lParse.Nodes[this_node_idx]).value.data.log = circle(
                        pVals[0].data.dbl,
                        pVals[1].data.dbl,
                        pVals[2].data.dbl,
                        pVals[3].data.dbl,
                        pVals[4].data.dbl,
                    );
                    current_block_139 = 7627602990488000394;
                }
                1028 => {
                    (lParse.Nodes[this_node_idx]).value.data.log = saobox(
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
                    (lParse.Nodes[this_node_idx]).value.data.log = ellipse(
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
                    match (lParse.Nodes[this_node_idx]).ntype.into() {
                        fits_parser_yytokentype::BOOLEAN => {
                            (lParse.Nodes[this_node_idx]).value.data.log =
                                (if c_int::from(pVals[2].data.log) != 0 {
                                    c_int::from(pVals[0].data.log)
                                } else {
                                    c_int::from(pVals[1].data.log)
                                }) as c_char;
                        }
                        fits_parser_yytokentype::LONG => {
                            (lParse.Nodes[this_node_idx]).value.data.lng =
                                if c_int::from(pVals[2].data.log) != 0 {
                                    pVals[0].data.lng
                                } else {
                                    pVals[1].data.lng
                                };
                        }
                        fits_parser_yytokentype::DOUBLE => {
                            (lParse.Nodes[this_node_idx]).value.data.dbl =
                                if c_int::from(pVals[2].data.log) != 0 {
                                    pVals[0].data.dbl
                                } else {
                                    pVals[1].data.dbl
                                };
                        }
                        fits_parser_yytokentype::STRING => {
                            strcpy(
                                ((lParse.Nodes[this_node_idx]).value.data.astr).as_mut_ptr(),
                                if c_int::from(pVals[2].data.log) != 0 {
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
                STRMID_FCT => {
                    let dest_str = ((lParse.Nodes[this_node_idx]).value.data.astr).as_mut_ptr();
                    let dest_len = (lParse.Nodes[this_node_idx]).value.nelem as c_int;
                    cstrmid(
                        lParse,
                        dest_str,
                        dest_len,
                        (pVals[0].data.astr).as_mut_ptr(),
                        pVals[0].nelem as c_int,
                        pVals[1].data.lng as c_int,
                    );
                    current_block_139 = 7627602990488000394;
                }
                1045 => {
                    let res: *mut c_char = strstr(
                        (pVals[0].data.astr).as_mut_ptr(),
                        (pVals[1].data.astr).as_mut_ptr(),
                    );
                    if res.is_null() {
                        (lParse.Nodes[this_node_idx]).value.data.lng = 0;
                    } else {
                        (lParse.Nodes[this_node_idx]).value.data.lng =
                            res.offset_from((pVals[0].data.astr).as_mut_ptr()) as c_long + 1;
                    }
                    current_block_139 = 7627602990488000394;
                }
                _ => {
                    current_block_139 = 7627602990488000394;
                }
            }
            if current_block_139 == 15934000668868306918 {
                if (lParse.Nodes[this_node_idx]).ntype == fits_parser_yytokentype::DOUBLE as c_int {
                    (lParse.Nodes[this_node_idx]).value.data.dbl = pVals[0].data.dbl;
                } else if (lParse.Nodes[this_node_idx]).ntype
                    == fits_parser_yytokentype::LONG as c_int
                {
                    (lParse.Nodes[this_node_idx]).value.data.lng = pVals[0].data.lng;
                } else if (lParse.Nodes[this_node_idx]).ntype
                    == fits_parser_yytokentype::BITSTR as c_int
                {
                    strcpy(
                        ((lParse.Nodes[this_node_idx]).value.data.astr).as_mut_ptr(),
                        (pVals[0].data.astr).as_mut_ptr(),
                    );
                }
            }
            (lParse.Nodes[this_node_idx]).operation = CONST_OP;
        } else {
            Allocate_Ptrs(lParse, this_node_idx);
            row = lParse.nRows;
            elem = row * (lParse.Nodes[this_node_idx]).value.nelem;
            if lParse.status == 0 {
                match (lParse.Nodes[this_node_idx]).operation as u32 {
                    1035 => loop {
                        let fresh53 = row;
                        row -= 1;
                        if fresh53 == 0 {
                            break;
                        }
                        *((lParse.Nodes[this_node_idx]).value.data.lngptr).offset(row as isize) =
                            lParse.firstRow + row;
                        *((lParse.Nodes[this_node_idx]).value.undef).offset(row as isize) = 0;
                    },
                    1036 => {
                        if (lParse.Nodes[this_node_idx]).ntype
                            == fits_parser_yytokentype::LONG as c_int
                        {
                            loop {
                                let fresh54 = row;
                                row -= 1;
                                if fresh54 == 0 {
                                    break;
                                }
                                *((lParse.Nodes[this_node_idx]).value.data.lngptr)
                                    .offset(row as isize) = 0;
                                *((lParse.Nodes[this_node_idx]).value.undef).offset(row as isize) =
                                    1;
                            }
                        } else if (lParse.Nodes[this_node_idx]).ntype
                            == fits_parser_yytokentype::STRING as c_int
                        {
                            loop {
                                let fresh55 = row;
                                row -= 1;
                                if fresh55 == 0 {
                                    break;
                                }
                                *(*((lParse.Nodes[this_node_idx]).value.data.strptr)
                                    .offset(row as isize))
                                .offset(0) = 0;
                                *((lParse.Nodes[this_node_idx]).value.undef).offset(row as isize) =
                                    1;
                            }
                        }
                    }
                    1050 => {
                        let mut ielem: c_long = 0;
                        let mut iaxis: [c_long; 5] = [1, 1, 1, 1, 1];
                        let ipos: c_long = pVals[1].data.lng - 1;
                        let naxis: c_int = (lParse.Nodes[this_node_idx]).value.naxis;
                        let mut j: c_int = 0;
                        if ipos < 0 || ipos >= 5 as c_long {
                            fits_parser_yyerror(
                                lParse,
                                cs!(c"AXISELEM(V,n) n value exceeded maximum dimension"),
                            );
                            free((lParse.Nodes[this_node_idx]).value.data.ptr);
                        } else {
                            ielem = 0;
                            while ielem < elem {
                                *((lParse.Nodes[this_node_idx]).value.data.lngptr)
                                    .offset(ielem as isize) = iaxis[ipos as usize];
                                *((lParse.Nodes[this_node_idx]).value.undef)
                                    .offset(ielem as isize) = 0;
                                iaxis[0] += 1;
                                iaxis[0];
                                j = 0;
                                while j < naxis {
                                    if iaxis[j as usize]
                                        <= (lParse.Nodes[this_node_idx]).value.naxes[j as usize]
                                    {
                                        break;
                                    }
                                    iaxis[j as usize] = 1;
                                    if j < naxis - 1 {
                                        iaxis[(j + 1) as usize] += 1;
                                        iaxis[(j + 1) as usize];
                                    }
                                    j += 1;
                                }
                                ielem += 1;
                            }
                        }
                    }
                    1049 => {
                        let mut ielem_0: c_long = 0;
                        let mut elemnum: c_long = 1;
                        let j_0: c_int = 0;
                        ielem_0 = 0;
                        while ielem_0 < elem {
                            *((lParse.Nodes[this_node_idx]).value.data.lngptr)
                                .offset(ielem_0 as isize) = elemnum;
                            *((lParse.Nodes[this_node_idx]).value.undef).offset(ielem_0 as isize) =
                                0;
                            elemnum += 1;
                            if elemnum > (lParse.Nodes[this_node_idx]).value.nelem {
                                elemnum = 1;
                            }
                            ielem_0 += 1;
                        }
                    }
                    1001 => loop {
                        let fresh56 = elem;
                        elem -= 1;
                        if fresh56 == 0 {
                            break;
                        }
                        *((lParse.Nodes[this_node_idx]).value.data.dblptr).offset(elem as isize) =
                            simplerng_getuniform();
                        *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize) = 0;
                    },
                    1042 => loop {
                        let fresh57 = elem;
                        elem -= 1;
                        if fresh57 == 0 {
                            break;
                        }
                        *((lParse.Nodes[this_node_idx]).value.data.dblptr).offset(elem as isize) =
                            simplerng_getnorm();
                        *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize) = 0;
                    },
                    1043 => {
                        if (lParse.Nodes[theParams[0]]).ntype
                            == fits_parser_yytokentype::DOUBLE as c_int
                        {
                            if (lParse.Nodes[theParams[0]]).operation == CONST_OP {
                                loop {
                                    let fresh58 = elem;
                                    elem -= 1;
                                    if fresh58 == 0 {
                                        break;
                                    }
                                    *((lParse.Nodes[this_node_idx]).value.undef)
                                        .offset(elem as isize) =
                                        if pVals[0].data.dbl < 0.0 { 1 } else { 0 };
                                    if *((lParse.Nodes[this_node_idx]).value.undef)
                                        .offset(elem as isize)
                                        == 0
                                    {
                                        *((lParse.Nodes[this_node_idx]).value.data.lngptr)
                                            .offset(elem as isize) =
                                            c_long::from(simplerng_getpoisson(pVals[0].data.dbl));
                                    }
                                }
                            } else {
                                loop {
                                    let fresh59 = elem;
                                    elem -= 1;
                                    if fresh59 == 0 {
                                        break;
                                    }
                                    *((lParse.Nodes[this_node_idx]).value.undef)
                                        .offset(elem as isize) =
                                        *((lParse.Nodes[theParams[0]]).value.undef)
                                            .offset(elem as isize);
                                    if *((lParse.Nodes[theParams[0]]).value.data.dblptr)
                                        .offset(elem as isize)
                                        < 0.0
                                    {
                                        *((lParse.Nodes[this_node_idx]).value.undef)
                                            .offset(elem as isize) = 1;
                                    }
                                    if *((lParse.Nodes[this_node_idx]).value.undef)
                                        .offset(elem as isize)
                                        == 0
                                    {
                                        *((lParse.Nodes[this_node_idx]).value.data.lngptr)
                                            .offset(elem as isize) =
                                            c_long::from(simplerng_getpoisson(
                                                *((lParse.Nodes[theParams[0]]).value.data.dblptr)
                                                    .offset(elem as isize),
                                            ));
                                    }
                                }
                            }
                        } else if (lParse.Nodes[theParams[0]]).operation == CONST_OP {
                            loop {
                                let fresh60 = elem;
                                elem -= 1;
                                if fresh60 == 0 {
                                    break;
                                }
                                *((lParse.Nodes[this_node_idx]).value.undef)
                                    .offset(elem as isize) =
                                    if pVals[0].data.lng < 0 { 1 } else { 0 };
                                if *((lParse.Nodes[this_node_idx]).value.undef)
                                    .offset(elem as isize)
                                    == 0
                                {
                                    *((lParse.Nodes[this_node_idx]).value.data.lngptr)
                                        .offset(elem as isize) = c_long::from(
                                        simplerng_getpoisson(pVals[0].data.lng as c_double),
                                    );
                                }
                            }
                        } else {
                            loop {
                                let fresh61 = elem;
                                elem -= 1;
                                if fresh61 == 0 {
                                    break;
                                }
                                *((lParse.Nodes[this_node_idx]).value.undef)
                                    .offset(elem as isize) =
                                    *((lParse.Nodes[theParams[0]]).value.undef)
                                        .offset(elem as isize);
                                if *((lParse.Nodes[theParams[0]]).value.data.lngptr)
                                    .offset(elem as isize)
                                    < 0
                                {
                                    *((lParse.Nodes[this_node_idx]).value.undef)
                                        .offset(elem as isize) = 1;
                                }
                                if *((lParse.Nodes[this_node_idx]).value.undef)
                                    .offset(elem as isize)
                                    == 0
                                {
                                    *((lParse.Nodes[this_node_idx]).value.data.lngptr)
                                        .offset(elem as isize) =
                                        c_long::from(simplerng_getpoisson(
                                            *((lParse.Nodes[theParams[0]]).value.data.lngptr)
                                                .offset(elem as isize)
                                                as c_double,
                                        ));
                                }
                            }
                        }
                    }
                    1002 => {
                        elem = row * (lParse.Nodes[theParams[0]]).value.nelem;
                        if (lParse.Nodes[theParams[0]]).ntype
                            == fits_parser_yytokentype::BOOLEAN as c_int
                        {
                            loop {
                                let fresh62 = row;
                                row -= 1;
                                if fresh62 == 0 {
                                    break;
                                }
                                *((lParse.Nodes[this_node_idx]).value.data.lngptr)
                                    .offset(row as isize) = 0;
                                *((lParse.Nodes[this_node_idx]).value.undef).offset(row as isize) =
                                    1;
                                nelem = (lParse.Nodes[theParams[0]]).value.nelem;
                                loop {
                                    let fresh63 = nelem;
                                    nelem -= 1;
                                    if fresh63 == 0 {
                                        break;
                                    }
                                    elem -= 1;

                                    if *((lParse.Nodes[theParams[0]]).value.undef)
                                        .offset(elem as isize)
                                        == 0
                                    {
                                        *((lParse.Nodes[this_node_idx]).value.data.lngptr)
                                            .offset(row as isize) += c_long::from(
                                            if c_int::from(
                                                *((lParse.Nodes[theParams[0]]).value.data.logptr)
                                                    .offset(elem as isize),
                                            ) != 0
                                            {
                                                1
                                            } else {
                                                0
                                            },
                                        );
                                        *((lParse.Nodes[this_node_idx]).value.undef)
                                            .offset(row as isize) = 0;
                                    }
                                }
                            }
                        } else if (lParse.Nodes[theParams[0]]).ntype
                            == fits_parser_yytokentype::LONG as c_int
                        {
                            loop {
                                let fresh64 = row;
                                row -= 1;
                                if fresh64 == 0 {
                                    break;
                                }
                                *((lParse.Nodes[this_node_idx]).value.data.lngptr)
                                    .offset(row as isize) = 0;
                                *((lParse.Nodes[this_node_idx]).value.undef).offset(row as isize) =
                                    1;
                                nelem = (lParse.Nodes[theParams[0]]).value.nelem;
                                loop {
                                    let fresh65 = nelem;
                                    nelem -= 1;
                                    if fresh65 == 0 {
                                        break;
                                    }
                                    elem -= 1;

                                    if *((lParse.Nodes[theParams[0]]).value.undef)
                                        .offset(elem as isize)
                                        == 0
                                    {
                                        *((lParse.Nodes[this_node_idx]).value.data.lngptr)
                                            .offset(row as isize) +=
                                            *((lParse.Nodes[theParams[0]]).value.data.lngptr)
                                                .offset(elem as isize);
                                        *((lParse.Nodes[this_node_idx]).value.undef)
                                            .offset(row as isize) = 0;
                                    }
                                }
                            }
                        } else if (lParse.Nodes[theParams[0]]).ntype
                            == fits_parser_yytokentype::DOUBLE as c_int
                        {
                            loop {
                                let fresh66 = row;
                                row -= 1;
                                if fresh66 == 0 {
                                    break;
                                }
                                *((lParse.Nodes[this_node_idx]).value.data.dblptr)
                                    .offset(row as isize) = 0.0;
                                *((lParse.Nodes[this_node_idx]).value.undef).offset(row as isize) =
                                    1;
                                nelem = (lParse.Nodes[theParams[0]]).value.nelem;
                                loop {
                                    let fresh67 = nelem;
                                    nelem -= 1;
                                    if fresh67 == 0 {
                                        break;
                                    }
                                    elem -= 1;

                                    if *((lParse.Nodes[theParams[0]]).value.undef)
                                        .offset(elem as isize)
                                        == 0
                                    {
                                        *((lParse.Nodes[this_node_idx]).value.data.dblptr)
                                            .offset(row as isize) +=
                                            *((lParse.Nodes[theParams[0]]).value.data.dblptr)
                                                .offset(elem as isize);
                                        *((lParse.Nodes[this_node_idx]).value.undef)
                                            .offset(row as isize) = 0;
                                    }
                                }
                            }
                        } else {
                            nelem = (lParse.Nodes[theParams[0]]).value.nelem;
                            loop {
                                let fresh68 = row;
                                row -= 1;
                                if fresh68 == 0 {
                                    break;
                                }
                                let mut sptr1: *mut c_char =
                                    *((lParse.Nodes[theParams[0]]).value.data.strptr)
                                        .offset(row as isize);
                                *((lParse.Nodes[this_node_idx]).value.data.lngptr)
                                    .offset(row as isize) = 0;
                                *((lParse.Nodes[this_node_idx]).value.undef).offset(row as isize) =
                                    0;
                                while *sptr1 != 0 {
                                    if c_int::from(*sptr1) == '1' as i32 {
                                        let fresh69 =
                                            &mut *((lParse.Nodes[this_node_idx]).value.data.lngptr)
                                                .offset(row as isize);
                                        *fresh69 += 1;
                                    }
                                    sptr1 = sptr1.offset(1);
                                }
                            }
                        }
                    }
                    1038 => {
                        elem = row * (lParse.Nodes[theParams[0]]).value.nelem;
                        if (lParse.Nodes[theParams[0]]).ntype
                            == fits_parser_yytokentype::LONG as c_int
                        {
                            loop {
                                let fresh70 = row;
                                row -= 1;
                                if fresh70 == 0 {
                                    break;
                                }
                                let mut count: c_int = 0;
                                *((lParse.Nodes[this_node_idx]).value.data.dblptr)
                                    .offset(row as isize) = 0.0;
                                nelem = (lParse.Nodes[theParams[0]]).value.nelem;
                                loop {
                                    let fresh71 = nelem;
                                    nelem -= 1;
                                    if fresh71 == 0 {
                                        break;
                                    }
                                    elem -= 1;

                                    if c_int::from(
                                        *((lParse.Nodes[theParams[0]]).value.undef)
                                            .offset(elem as isize),
                                    ) == 0
                                    {
                                        *((lParse.Nodes[this_node_idx]).value.data.dblptr)
                                            .offset(row as isize) +=
                                            *((lParse.Nodes[theParams[0]]).value.data.lngptr)
                                                .offset(elem as isize)
                                                as c_double;
                                        count += 1;
                                    }
                                }
                                if count == 0 {
                                    *((lParse.Nodes[this_node_idx]).value.undef)
                                        .offset(row as isize) = 1;
                                } else {
                                    *((lParse.Nodes[this_node_idx]).value.undef)
                                        .offset(row as isize) = 0;
                                    *((lParse.Nodes[this_node_idx]).value.data.dblptr)
                                        .offset(row as isize) /= c_double::from(count);
                                }
                            }
                        } else if (lParse.Nodes[theParams[0]]).ntype
                            == fits_parser_yytokentype::DOUBLE as c_int
                        {
                            loop {
                                let fresh72 = row;
                                row -= 1;
                                if fresh72 == 0 {
                                    break;
                                }
                                let mut count_0: c_int = 0;
                                *((lParse.Nodes[this_node_idx]).value.data.dblptr)
                                    .offset(row as isize) = 0.0;
                                nelem = (lParse.Nodes[theParams[0]]).value.nelem;
                                loop {
                                    let fresh73 = nelem;
                                    nelem -= 1;
                                    if fresh73 == 0 {
                                        break;
                                    }
                                    elem -= 1;

                                    if c_int::from(
                                        *((lParse.Nodes[theParams[0]]).value.undef)
                                            .offset(elem as isize),
                                    ) == 0
                                    {
                                        *((lParse.Nodes[this_node_idx]).value.data.dblptr)
                                            .offset(row as isize) +=
                                            *((lParse.Nodes[theParams[0]]).value.data.dblptr)
                                                .offset(elem as isize);
                                        count_0 += 1;
                                    }
                                }
                                if count_0 == 0 {
                                    *((lParse.Nodes[this_node_idx]).value.undef)
                                        .offset(row as isize) = 1;
                                } else {
                                    *((lParse.Nodes[this_node_idx]).value.undef)
                                        .offset(row as isize) = 0;
                                    *((lParse.Nodes[this_node_idx]).value.data.dblptr)
                                        .offset(row as isize) /= c_double::from(count_0);
                                }
                            }
                        }
                    }
                    1039 => {
                        elem = row * (lParse.Nodes[theParams[0]]).value.nelem;
                        if (lParse.Nodes[theParams[0]]).ntype
                            == fits_parser_yytokentype::LONG as c_int
                        {
                            loop {
                                let fresh74 = row;
                                row -= 1;
                                if fresh74 == 0 {
                                    break;
                                }
                                let mut count_1: c_int = 0;
                                let mut sum: c_double = 0.0;
                                let mut sum2: c_double = 0.0;
                                nelem = (lParse.Nodes[theParams[0]]).value.nelem;
                                loop {
                                    let fresh75 = nelem;
                                    nelem -= 1;
                                    if fresh75 == 0 {
                                        break;
                                    }
                                    elem -= 1;

                                    if c_int::from(
                                        *((lParse.Nodes[theParams[0]]).value.undef)
                                            .offset(elem as isize),
                                    ) == 0
                                    {
                                        sum += *((lParse.Nodes[theParams[0]]).value.data.lngptr)
                                            .offset(elem as isize)
                                            as c_double;
                                        count_1 += 1;
                                    }
                                }
                                if count_1 > 1 {
                                    sum /= c_double::from(count_1);
                                    nelem = (lParse.Nodes[theParams[0]]).value.nelem;
                                    elem += nelem;
                                    loop {
                                        let fresh76 = nelem;
                                        nelem -= 1;
                                        if fresh76 == 0 {
                                            break;
                                        }
                                        elem -= 1;

                                        if c_int::from(
                                            *((lParse.Nodes[theParams[0]]).value.undef)
                                                .offset(elem as isize),
                                        ) == 0
                                        {
                                            let dx: c_double =
                                                *((lParse.Nodes[theParams[0]]).value.data.lngptr)
                                                    .offset(elem as isize)
                                                    as c_double
                                                    - sum;
                                            sum2 += dx * dx;
                                        }
                                    }
                                    sum2 /= c_double::from(count_1) - c_double::from(1);
                                    *((lParse.Nodes[this_node_idx]).value.undef)
                                        .offset(row as isize) = 0;
                                    *((lParse.Nodes[this_node_idx]).value.data.dblptr)
                                        .offset(row as isize) = sqrt(sum2);
                                } else {
                                    *((lParse.Nodes[this_node_idx]).value.undef)
                                        .offset(row as isize) = 0;
                                    *((lParse.Nodes[this_node_idx]).value.data.dblptr)
                                        .offset(row as isize) = 0.0;
                                }
                            }
                        } else if (lParse.Nodes[theParams[0]]).ntype
                            == fits_parser_yytokentype::DOUBLE as c_int
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
                                nelem = (lParse.Nodes[theParams[0]]).value.nelem;
                                loop {
                                    let fresh78 = nelem;
                                    nelem -= 1;
                                    if fresh78 == 0 {
                                        break;
                                    }
                                    elem -= 1;

                                    if c_int::from(
                                        *((lParse.Nodes[theParams[0]]).value.undef)
                                            .offset(elem as isize),
                                    ) == 0
                                    {
                                        sum_0 += *((lParse.Nodes[theParams[0]]).value.data.dblptr)
                                            .offset(elem as isize);
                                        count_2 += 1;
                                    }
                                }
                                if count_2 > 1 {
                                    sum_0 /= c_double::from(count_2);
                                    nelem = (lParse.Nodes[theParams[0]]).value.nelem;
                                    elem += nelem;
                                    loop {
                                        let fresh79 = nelem;
                                        nelem -= 1;
                                        if fresh79 == 0 {
                                            break;
                                        }
                                        elem -= 1;

                                        if c_int::from(
                                            *((lParse.Nodes[theParams[0]]).value.undef)
                                                .offset(elem as isize),
                                        ) == 0
                                        {
                                            let dx_0: c_double =
                                                *((lParse.Nodes[theParams[0]]).value.data.dblptr)
                                                    .offset(elem as isize)
                                                    - sum_0;
                                            sum2_0 += dx_0 * dx_0;
                                        }
                                    }
                                    sum2_0 /= c_double::from(count_2) - c_double::from(1);
                                    *((lParse.Nodes[this_node_idx]).value.undef)
                                        .offset(row as isize) = 0;
                                    *((lParse.Nodes[this_node_idx]).value.data.dblptr)
                                        .offset(row as isize) = sqrt(sum2_0);
                                } else {
                                    *((lParse.Nodes[this_node_idx]).value.undef)
                                        .offset(row as isize) = 0;
                                    *((lParse.Nodes[this_node_idx]).value.data.dblptr)
                                        .offset(row as isize) = 0.0;
                                }
                            }
                        }
                    }
                    1037 => {
                        elem = row * (lParse.Nodes[theParams[0]]).value.nelem;
                        nelem = (lParse.Nodes[theParams[0]]).value.nelem;
                        if (lParse.Nodes[theParams[0]]).ntype
                            == fits_parser_yytokentype::LONG as c_int
                        {
                            let mut dptr: *mut c_long =
                                (lParse.Nodes[theParams[0]]).value.data.lngptr;
                            let mut uptr: *mut c_char = (lParse.Nodes[theParams[0]]).value.undef;
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
                                free((lParse.Nodes[this_node_idx]).value.data.ptr);
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
                                        *((lParse.Nodes[this_node_idx]).value.undef)
                                            .offset(irow as isize) = 0;
                                        *((lParse.Nodes[this_node_idx]).value.data.lngptr)
                                            .offset(irow as isize) =
                                            qselect_median_lng(mptr, nelem1);
                                    } else {
                                        *((lParse.Nodes[this_node_idx]).value.undef)
                                            .offset(irow as isize) = 1;
                                        *((lParse.Nodes[this_node_idx]).value.data.lngptr)
                                            .offset(irow as isize) = 0;
                                    }
                                    irow += 1;
                                }
                                free(mptr.cast::<c_void>());
                            }
                        } else {
                            let mut dptr_0: *mut c_double =
                                (lParse.Nodes[theParams[0]]).value.data.dblptr;
                            let mut uptr_0: *mut c_char = (lParse.Nodes[theParams[0]]).value.undef;
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
                                free((lParse.Nodes[this_node_idx]).value.data.ptr);
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
                                        *((lParse.Nodes[this_node_idx]).value.undef)
                                            .offset(irow_0 as isize) = 0;
                                        *((lParse.Nodes[this_node_idx]).value.data.dblptr)
                                            .offset(irow_0 as isize) =
                                            qselect_median_dbl(mptr_0, nelem1_0);
                                    } else {
                                        *((lParse.Nodes[this_node_idx]).value.undef)
                                            .offset(irow_0 as isize) = 1;
                                        *((lParse.Nodes[this_node_idx]).value.data.dblptr)
                                            .offset(irow_0 as isize) = 0.0;
                                    }
                                    irow_0 += 1;
                                }
                                free(mptr_0.cast::<c_void>());
                            }
                        }
                    }
                    1017 => {
                        if (lParse.Nodes[theParams[0]]).ntype
                            == fits_parser_yytokentype::DOUBLE as c_int
                        {
                            loop {
                                let fresh84 = elem;
                                elem -= 1;
                                if fresh84 == 0 {
                                    break;
                                }
                                dval = *((lParse.Nodes[theParams[0]]).value.data.dblptr)
                                    .offset(elem as isize);
                                *((lParse.Nodes[this_node_idx]).value.data.dblptr)
                                    .offset(elem as isize) = if dval > 0.0 { dval } else { -dval };
                                *((lParse.Nodes[this_node_idx]).value.undef)
                                    .offset(elem as isize) =
                                    *((lParse.Nodes[theParams[0]]).value.undef)
                                        .offset(elem as isize);
                            }
                        } else {
                            loop {
                                let fresh85 = elem;
                                elem -= 1;
                                if fresh85 == 0 {
                                    break;
                                }
                                ival = *((lParse.Nodes[theParams[0]]).value.data.lngptr)
                                    .offset(elem as isize);
                                *((lParse.Nodes[this_node_idx]).value.data.lngptr)
                                    .offset(elem as isize) = if ival > 0 { ival } else { -ival };
                                *((lParse.Nodes[this_node_idx]).value.undef)
                                    .offset(elem as isize) =
                                    *((lParse.Nodes[theParams[0]]).value.undef)
                                        .offset(elem as isize);
                            }
                        }
                    }
                    1040 => {
                        nelem = (lParse.Nodes[theParams[0]]).value.nelem;
                        if (lParse.Nodes[theParams[0]]).ntype
                            == fits_parser_yytokentype::STRING as c_int
                        {
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
                            *((lParse.Nodes[this_node_idx]).value.data.lngptr)
                                .offset(row as isize) = 0;
                            loop {
                                let fresh87 = nelem1_1;
                                nelem1_1 -= 1;
                                if fresh87 == 0 {
                                    break;
                                }
                                elem -= 1;

                                if c_int::from(
                                    *((lParse.Nodes[theParams[0]]).value.undef)
                                        .offset(elem as isize),
                                ) == 0
                                {
                                    let fresh88 =
                                        &mut *((lParse.Nodes[this_node_idx]).value.data.lngptr)
                                            .offset(row as isize);
                                    *fresh88 += 1;
                                }
                            }
                        }
                    }
                    1030 => {
                        if (lParse.Nodes[theParams[0]]).ntype
                            == fits_parser_yytokentype::STRING as c_int
                        {
                            elem = row;
                        }
                        loop {
                            let fresh89 = elem;
                            elem -= 1;
                            if fresh89 == 0 {
                                break;
                            }
                            *((lParse.Nodes[this_node_idx]).value.data.logptr)
                                .offset(elem as isize) =
                                *((lParse.Nodes[theParams[0]]).value.undef).offset(elem as isize);
                            *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize) = 0;
                        }
                    }
                    1031 => match (lParse.Nodes[this_node_idx]).ntype.into() {
                        fits_parser_yytokentype::BOOLEAN => loop {
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
                                    if vector[i as usize] > 1 {
                                        pNull[i as usize] =
                                            *((lParse.Nodes[theParams[i as usize]]).value.undef)
                                                .offset(elem as isize);
                                        pVals[i as usize].data.log = *((lParse.Nodes
                                            [theParams[i as usize]])
                                            .value
                                            .data
                                            .logptr)
                                            .offset(elem as isize);
                                    } else if vector[i as usize] != 0 {
                                        pNull[i as usize] =
                                            *((lParse.Nodes[theParams[i as usize]]).value.undef)
                                                .offset(row as isize);
                                        pVals[i as usize].data.log = *((lParse.Nodes
                                            [theParams[i as usize]])
                                            .value
                                            .data
                                            .logptr)
                                            .offset(row as isize);
                                    }
                                }
                                if pNull[0] != 0 {
                                    *((lParse.Nodes[this_node_idx]).value.undef)
                                        .offset(elem as isize) = pNull[1];
                                    *((lParse.Nodes[this_node_idx]).value.data.logptr)
                                        .offset(elem as isize) = pVals[1].data.log;
                                } else {
                                    *((lParse.Nodes[this_node_idx]).value.undef)
                                        .offset(elem as isize) = 0;
                                    *((lParse.Nodes[this_node_idx]).value.data.logptr)
                                        .offset(elem as isize) = pVals[0].data.log;
                                }
                            }
                        },
                        fits_parser_yytokentype::LONG => loop {
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
                                    if vector[i as usize] > 1 {
                                        pNull[i as usize] =
                                            *((lParse.Nodes[theParams[i as usize]]).value.undef)
                                                .offset(elem as isize);
                                        pVals[i as usize].data.lng = *((lParse.Nodes
                                            [theParams[i as usize]])
                                            .value
                                            .data
                                            .lngptr)
                                            .offset(elem as isize);
                                    } else if vector[i as usize] != 0 {
                                        pNull[i as usize] =
                                            *((lParse.Nodes[theParams[i as usize]]).value.undef)
                                                .offset(row as isize);
                                        pVals[i as usize].data.lng = *((lParse.Nodes
                                            [theParams[i as usize]])
                                            .value
                                            .data
                                            .lngptr)
                                            .offset(row as isize);
                                    }
                                }
                                if pNull[0] != 0 {
                                    *((lParse.Nodes[this_node_idx]).value.undef)
                                        .offset(elem as isize) = pNull[1];
                                    *((lParse.Nodes[this_node_idx]).value.data.lngptr)
                                        .offset(elem as isize) = pVals[1].data.lng;
                                } else {
                                    *((lParse.Nodes[this_node_idx]).value.undef)
                                        .offset(elem as isize) = 0;
                                    *((lParse.Nodes[this_node_idx]).value.data.lngptr)
                                        .offset(elem as isize) = pVals[0].data.lng;
                                }
                            }
                        },
                        fits_parser_yytokentype::DOUBLE => loop {
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
                                    if vector[i as usize] > 1 {
                                        pNull[i as usize] =
                                            *((lParse.Nodes[theParams[i as usize]]).value.undef)
                                                .offset(elem as isize);
                                        pVals[i as usize].data.dbl = *((lParse.Nodes
                                            [theParams[i as usize]])
                                            .value
                                            .data
                                            .dblptr)
                                            .offset(elem as isize);
                                    } else if vector[i as usize] != 0 {
                                        pNull[i as usize] =
                                            *((lParse.Nodes[theParams[i as usize]]).value.undef)
                                                .offset(row as isize);
                                        pVals[i as usize].data.dbl = *((lParse.Nodes
                                            [theParams[i as usize]])
                                            .value
                                            .data
                                            .dblptr)
                                            .offset(row as isize);
                                    }
                                }
                                if pNull[0] != 0 {
                                    *((lParse.Nodes[this_node_idx]).value.undef)
                                        .offset(elem as isize) = pNull[1];
                                    *((lParse.Nodes[this_node_idx]).value.data.dblptr)
                                        .offset(elem as isize) = pVals[1].data.dbl;
                                } else {
                                    *((lParse.Nodes[this_node_idx]).value.undef)
                                        .offset(elem as isize) = 0;
                                    *((lParse.Nodes[this_node_idx]).value.data.dblptr)
                                        .offset(elem as isize) = pVals[0].data.dbl;
                                }
                            }
                        },
                        fits_parser_yytokentype::STRING => loop {
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
                                if vector[i as usize] != 0 {
                                    pNull[i as usize] =
                                        *((lParse.Nodes[theParams[i as usize]]).value.undef)
                                            .offset(row as isize);
                                    strcpy(
                                        (pVals[i as usize].data.astr).as_mut_ptr(),
                                        *((lParse.Nodes[theParams[i as usize]]).value.data.strptr)
                                            .offset(row as isize),
                                    );
                                }
                            }
                            if pNull[0] != 0 {
                                *((lParse.Nodes[this_node_idx]).value.undef).offset(row as isize) =
                                    pNull[1];
                                strcpy(
                                    *((lParse.Nodes[this_node_idx]).value.data.strptr)
                                        .offset(row as isize),
                                    (pVals[1].data.astr).as_mut_ptr(),
                                );
                            } else {
                                *((lParse.Nodes[this_node_idx]).value.undef)
                                    .offset(elem as isize) = 0;
                                strcpy(
                                    *((lParse.Nodes[this_node_idx]).value.data.strptr)
                                        .offset(row as isize),
                                    (pVals[0].data.astr).as_mut_ptr(),
                                );
                            }
                        },
                        _ => {}
                    },
                    1046 => match (lParse.Nodes[this_node_idx]).ntype.into() {
                        fits_parser_yytokentype::LONG => loop {
                            let fresh101 = elem;
                            elem -= 1;
                            if fresh101 == 0 {
                                break;
                            }
                            if (lParse.Nodes[theParams[1]]).value.data.lng
                                == *((lParse.Nodes[theParams[0]]).value.data.lngptr)
                                    .offset(elem as isize)
                            {
                                *((lParse.Nodes[this_node_idx]).value.data.lngptr)
                                    .offset(elem as isize) = 0;
                                *((lParse.Nodes[this_node_idx]).value.undef)
                                    .offset(elem as isize) = 1;
                            } else {
                                *((lParse.Nodes[this_node_idx]).value.data.lngptr)
                                    .offset(elem as isize) =
                                    *((lParse.Nodes[theParams[0]]).value.data.lngptr)
                                        .offset(elem as isize);
                                *((lParse.Nodes[this_node_idx]).value.undef)
                                    .offset(elem as isize) =
                                    *((lParse.Nodes[theParams[0]]).value.undef)
                                        .offset(elem as isize);
                            }
                        },
                        fits_parser_yytokentype::DOUBLE => loop {
                            let fresh102 = elem;
                            elem -= 1;
                            if fresh102 == 0 {
                                break;
                            }
                            if (lParse.Nodes[theParams[1]]).value.data.dbl
                                == *((lParse.Nodes[theParams[0]]).value.data.dblptr)
                                    .offset(elem as isize)
                            {
                                *((lParse.Nodes[this_node_idx]).value.data.dblptr)
                                    .offset(elem as isize) = 0.0;
                                *((lParse.Nodes[this_node_idx]).value.undef)
                                    .offset(elem as isize) = 1;
                            } else {
                                *((lParse.Nodes[this_node_idx]).value.data.dblptr)
                                    .offset(elem as isize) =
                                    *((lParse.Nodes[theParams[0]]).value.data.dblptr)
                                        .offset(elem as isize);
                                *((lParse.Nodes[this_node_idx]).value.undef)
                                    .offset(elem as isize) =
                                    *((lParse.Nodes[theParams[0]]).value.undef)
                                        .offset(elem as isize);
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
                        let fresh104 =
                            &mut *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize);
                        *fresh104 =
                            *((lParse.Nodes[theParams[0]]).value.undef).offset(elem as isize);
                        if *fresh104 == 0 {
                            *((lParse.Nodes[this_node_idx]).value.data.dblptr)
                                .offset(elem as isize) =
                                sin(*((lParse.Nodes[theParams[0]]).value.data.dblptr)
                                    .offset(elem as isize));
                        }
                    },
                    1005 => loop {
                        let fresh105 = elem;
                        elem -= 1;
                        if fresh105 == 0 {
                            break;
                        }
                        let fresh106 =
                            &mut *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize);
                        *fresh106 =
                            *((lParse.Nodes[theParams[0]]).value.undef).offset(elem as isize);
                        if *fresh106 == 0 {
                            *((lParse.Nodes[this_node_idx]).value.data.dblptr)
                                .offset(elem as isize) =
                                cos(*((lParse.Nodes[theParams[0]]).value.data.dblptr)
                                    .offset(elem as isize));
                        }
                    },
                    1006 => loop {
                        let fresh107 = elem;
                        elem -= 1;
                        if fresh107 == 0 {
                            break;
                        }
                        let fresh108 =
                            &mut *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize);
                        *fresh108 =
                            *((lParse.Nodes[theParams[0]]).value.undef).offset(elem as isize);
                        if *fresh108 == 0 {
                            *((lParse.Nodes[this_node_idx]).value.data.dblptr)
                                .offset(elem as isize) =
                                tan(*((lParse.Nodes[theParams[0]]).value.data.dblptr)
                                    .offset(elem as isize));
                        }
                    },
                    1007 => loop {
                        let fresh109 = elem;
                        elem -= 1;
                        if fresh109 == 0 {
                            break;
                        }
                        let fresh110 =
                            &mut *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize);
                        *fresh110 =
                            *((lParse.Nodes[theParams[0]]).value.undef).offset(elem as isize);
                        if *fresh110 == 0 {
                            dval = *((lParse.Nodes[theParams[0]]).value.data.dblptr)
                                .offset(elem as isize);
                            if dval < -1.0 || dval > 1.0 {
                                *((lParse.Nodes[this_node_idx]).value.data.dblptr)
                                    .offset(elem as isize) = 0.0;
                                *((lParse.Nodes[this_node_idx]).value.undef)
                                    .offset(elem as isize) = 1;
                            } else {
                                *((lParse.Nodes[this_node_idx]).value.data.dblptr)
                                    .offset(elem as isize) = asin(dval);
                            }
                        }
                    },
                    1008 => loop {
                        let fresh111 = elem;
                        elem -= 1;
                        if fresh111 == 0 {
                            break;
                        }
                        let fresh112 =
                            &mut *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize);
                        *fresh112 =
                            *((lParse.Nodes[theParams[0]]).value.undef).offset(elem as isize);
                        if *fresh112 == 0 {
                            dval = *((lParse.Nodes[theParams[0]]).value.data.dblptr)
                                .offset(elem as isize);
                            if dval < -1.0 || dval > 1.0 {
                                *((lParse.Nodes[this_node_idx]).value.data.dblptr)
                                    .offset(elem as isize) = 0.0;
                                *((lParse.Nodes[this_node_idx]).value.undef)
                                    .offset(elem as isize) = 1;
                            } else {
                                *((lParse.Nodes[this_node_idx]).value.data.dblptr)
                                    .offset(elem as isize) = acos(dval);
                            }
                        }
                    },
                    1009 => loop {
                        let fresh113 = elem;
                        elem -= 1;
                        if fresh113 == 0 {
                            break;
                        }
                        let fresh114 =
                            &mut *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize);
                        *fresh114 =
                            *((lParse.Nodes[theParams[0]]).value.undef).offset(elem as isize);
                        if *fresh114 == 0 {
                            dval = *((lParse.Nodes[theParams[0]]).value.data.dblptr)
                                .offset(elem as isize);
                            *((lParse.Nodes[this_node_idx]).value.data.dblptr)
                                .offset(elem as isize) = atan(dval);
                        }
                    },
                    1010 => loop {
                        let fresh115 = elem;
                        elem -= 1;
                        if fresh115 == 0 {
                            break;
                        }
                        let fresh116 =
                            &mut *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize);
                        *fresh116 =
                            *((lParse.Nodes[theParams[0]]).value.undef).offset(elem as isize);
                        if *fresh116 == 0 {
                            *((lParse.Nodes[this_node_idx]).value.data.dblptr)
                                .offset(elem as isize) = sinh(
                                *((lParse.Nodes[theParams[0]]).value.data.dblptr)
                                    .offset(elem as isize),
                            );
                        }
                    },
                    1011 => loop {
                        let fresh117 = elem;
                        elem -= 1;
                        if fresh117 == 0 {
                            break;
                        }
                        let fresh118 =
                            &mut *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize);
                        *fresh118 =
                            *((lParse.Nodes[theParams[0]]).value.undef).offset(elem as isize);
                        if *fresh118 == 0 {
                            *((lParse.Nodes[this_node_idx]).value.data.dblptr)
                                .offset(elem as isize) = cosh(
                                *((lParse.Nodes[theParams[0]]).value.data.dblptr)
                                    .offset(elem as isize),
                            );
                        }
                    },
                    TANH_FCT => loop {
                        let fresh119 = elem;
                        elem -= 1;
                        if fresh119 == 0 {
                            break;
                        }
                        let fresh120 =
                            &mut *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize);
                        *fresh120 =
                            *((lParse.Nodes[theParams[0]]).value.undef).offset(elem as isize);
                        if *fresh120 == 0 {
                            *((lParse.Nodes[this_node_idx]).value.data.dblptr)
                                .offset(elem as isize) = tanh(
                                *((lParse.Nodes[theParams[0]]).value.data.dblptr)
                                    .offset(elem as isize),
                            );
                        }
                    },
                    1013 => loop {
                        let fresh121 = elem;
                        elem -= 1;
                        if fresh121 == 0 {
                            break;
                        }
                        let fresh122 =
                            &mut *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize);
                        *fresh122 =
                            *((lParse.Nodes[theParams[0]]).value.undef).offset(elem as isize);
                        if *fresh122 == 0 {
                            dval = *((lParse.Nodes[theParams[0]]).value.data.dblptr)
                                .offset(elem as isize);
                            *((lParse.Nodes[this_node_idx]).value.data.dblptr)
                                .offset(elem as isize) = exp(dval);
                        }
                    },
                    1014 => loop {
                        let fresh123 = elem;
                        elem -= 1;
                        if fresh123 == 0 {
                            break;
                        }
                        let fresh124 =
                            &mut *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize);
                        *fresh124 =
                            *((lParse.Nodes[theParams[0]]).value.undef).offset(elem as isize);
                        if *fresh124 == 0 {
                            dval = *((lParse.Nodes[theParams[0]]).value.data.dblptr)
                                .offset(elem as isize);
                            if dval <= 0.0 {
                                *((lParse.Nodes[this_node_idx]).value.data.dblptr)
                                    .offset(elem as isize) = 0.0;
                                *((lParse.Nodes[this_node_idx]).value.undef)
                                    .offset(elem as isize) = 1;
                            } else {
                                *((lParse.Nodes[this_node_idx]).value.data.dblptr)
                                    .offset(elem as isize) = log(dval);
                            }
                        }
                    },
                    1015 => loop {
                        let fresh125 = elem;
                        elem -= 1;
                        if fresh125 == 0 {
                            break;
                        }
                        let fresh126 =
                            &mut *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize);
                        *fresh126 =
                            *((lParse.Nodes[theParams[0]]).value.undef).offset(elem as isize);
                        if *fresh126 == 0 {
                            dval = *((lParse.Nodes[theParams[0]]).value.data.dblptr)
                                .offset(elem as isize);
                            if dval <= 0.0 {
                                *((lParse.Nodes[this_node_idx]).value.data.dblptr)
                                    .offset(elem as isize) = 0.0;
                                *((lParse.Nodes[this_node_idx]).value.undef)
                                    .offset(elem as isize) = 1;
                            } else {
                                *((lParse.Nodes[this_node_idx]).value.data.dblptr)
                                    .offset(elem as isize) = log10(dval);
                            }
                        }
                    },
                    1016 => loop {
                        let fresh127 = elem;
                        elem -= 1;
                        if fresh127 == 0 {
                            break;
                        }
                        let fresh128 =
                            &mut *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize);
                        *fresh128 =
                            *((lParse.Nodes[theParams[0]]).value.undef).offset(elem as isize);
                        if *fresh128 == 0 {
                            dval = *((lParse.Nodes[theParams[0]]).value.data.dblptr)
                                .offset(elem as isize);
                            if dval < 0.0 {
                                *((lParse.Nodes[this_node_idx]).value.data.dblptr)
                                    .offset(elem as isize) = 0.0;
                                *((lParse.Nodes[this_node_idx]).value.undef)
                                    .offset(elem as isize) = 1;
                            } else {
                                *((lParse.Nodes[this_node_idx]).value.data.dblptr)
                                    .offset(elem as isize) = sqrt(dval);
                            }
                        }
                    },
                    1019 => loop {
                        let fresh129 = elem;
                        elem -= 1;
                        if fresh129 == 0 {
                            break;
                        }
                        let fresh130 =
                            &mut *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize);
                        *fresh130 =
                            *((lParse.Nodes[theParams[0]]).value.undef).offset(elem as isize);
                        if *fresh130 == 0 {
                            *((lParse.Nodes[this_node_idx]).value.data.dblptr)
                                .offset(elem as isize) = ceil(
                                *((lParse.Nodes[theParams[0]]).value.data.dblptr)
                                    .offset(elem as isize),
                            );
                        }
                    },
                    1020 => loop {
                        let fresh131 = elem;
                        elem -= 1;
                        if fresh131 == 0 {
                            break;
                        }
                        let fresh132 =
                            &mut *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize);
                        *fresh132 =
                            *((lParse.Nodes[theParams[0]]).value.undef).offset(elem as isize);
                        if *fresh132 == 0 {
                            *((lParse.Nodes[this_node_idx]).value.data.dblptr)
                                .offset(elem as isize) = floor(
                                *((lParse.Nodes[theParams[0]]).value.data.dblptr)
                                    .offset(elem as isize),
                            );
                        }
                    },
                    1021 => loop {
                        let fresh133 = elem;
                        elem -= 1;
                        if fresh133 == 0 {
                            break;
                        }
                        let fresh134 =
                            &mut *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize);
                        *fresh134 =
                            *((lParse.Nodes[theParams[0]]).value.undef).offset(elem as isize);
                        if *fresh134 == 0 {
                            *((lParse.Nodes[this_node_idx]).value.data.dblptr)
                                .offset(elem as isize) = floor(
                                *((lParse.Nodes[theParams[0]]).value.data.dblptr)
                                    .offset(elem as isize)
                                    + 0.5,
                            );
                        }
                    },
                    1018 => loop {
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
                                if vector[i as usize] > 1 {
                                    pVals[i as usize].data.dbl =
                                        *((lParse.Nodes[theParams[i as usize]]).value.data.dblptr)
                                            .offset(elem as isize);
                                    pNull[i as usize] =
                                        *((lParse.Nodes[theParams[i as usize]]).value.undef)
                                            .offset(elem as isize);
                                } else if vector[i as usize] != 0 {
                                    pVals[i as usize].data.dbl =
                                        *((lParse.Nodes[theParams[i as usize]]).value.data.dblptr)
                                            .offset(row as isize);
                                    pNull[i as usize] =
                                        *((lParse.Nodes[theParams[i as usize]]).value.undef)
                                            .offset(row as isize);
                                }
                            }
                            let fresh138 = &mut *((lParse.Nodes[this_node_idx]).value.undef)
                                .offset(elem as isize);
                            *fresh138 = if c_int::from(pNull[0]) != 0 || c_int::from(pNull[1]) != 0
                            {
                                1
                            } else {
                                0
                            };
                            if *fresh138 == 0 {
                                *((lParse.Nodes[this_node_idx]).value.data.dblptr)
                                    .offset(elem as isize) =
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
                                if vector[i as usize] > 1 {
                                    pVals[i as usize].data.dbl =
                                        *((lParse.Nodes[theParams[i as usize]]).value.data.dblptr)
                                            .offset(elem as isize);
                                    pNull[i as usize] =
                                        *((lParse.Nodes[theParams[i as usize]]).value.undef)
                                            .offset(elem as isize);
                                } else if vector[i as usize] != 0 {
                                    pVals[i as usize].data.dbl =
                                        *((lParse.Nodes[theParams[i as usize]]).value.data.dblptr)
                                            .offset(row as isize);
                                    pNull[i as usize] =
                                        *((lParse.Nodes[theParams[i as usize]]).value.undef)
                                            .offset(row as isize);
                                }
                            }
                            let fresh142 = &mut *((lParse.Nodes[this_node_idx]).value.undef)
                                .offset(elem as isize);
                            *fresh142 = c_int::from(
                                c_int::from(pNull[0]) != 0
                                    || c_int::from(pNull[1]) != 0
                                    || c_int::from(pNull[2]) != 0
                                    || c_int::from(pNull[3]) != 0,
                            ) as c_char;
                            if *fresh142 == 0 {
                                *((lParse.Nodes[this_node_idx]).value.data.dblptr)
                                    .offset(elem as isize) = angsep_calc(
                                    pVals[0].data.dbl,
                                    pVals[1].data.dbl,
                                    pVals[2].data.dbl,
                                    pVals[3].data.dbl,
                                );
                            }
                        }
                    },
                    1022 => {
                        elem = row * (lParse.Nodes[theParams[0]]).value.nelem;
                        if (lParse.Nodes[this_node_idx]).ntype
                            == fits_parser_yytokentype::LONG as c_int
                        {
                            let mut minVal: c_long = 0;
                            loop {
                                let fresh143 = row;
                                row -= 1;
                                if fresh143 == 0 {
                                    break;
                                }
                                valInit = 1;
                                *((lParse.Nodes[this_node_idx]).value.undef).offset(row as isize) =
                                    1;
                                nelem = (lParse.Nodes[theParams[0]]).value.nelem;
                                loop {
                                    let fresh144 = nelem;
                                    nelem -= 1;
                                    if fresh144 == 0 {
                                        break;
                                    }
                                    elem -= 1;
                                    if *((lParse.Nodes[theParams[0]]).value.undef)
                                        .offset(elem as isize)
                                        == 0
                                    {
                                        if valInit != 0 {
                                            valInit = 0;
                                            minVal =
                                                *((lParse.Nodes[theParams[0]]).value.data.lngptr)
                                                    .offset(elem as isize);
                                        } else {
                                            minVal = if minVal
                                                < *((lParse.Nodes[theParams[0]]).value.data.lngptr)
                                                    .offset(elem as isize)
                                            {
                                                minVal
                                            } else {
                                                *((lParse.Nodes[theParams[0]]).value.data.lngptr)
                                                    .offset(elem as isize)
                                            };
                                        }
                                        *((lParse.Nodes[this_node_idx]).value.undef)
                                            .offset(row as isize) = 0;
                                    }
                                }
                                *((lParse.Nodes[this_node_idx]).value.data.lngptr)
                                    .offset(row as isize) = minVal;
                            }
                        } else if (lParse.Nodes[this_node_idx]).ntype
                            == fits_parser_yytokentype::DOUBLE as c_int
                        {
                            let mut minVal_0: c_double = 0.0;
                            loop {
                                let fresh145 = row;
                                row -= 1;
                                if fresh145 == 0 {
                                    break;
                                }
                                valInit = 1;
                                *((lParse.Nodes[this_node_idx]).value.undef).offset(row as isize) =
                                    1;
                                nelem = (lParse.Nodes[theParams[0]]).value.nelem;
                                loop {
                                    let fresh146 = nelem;
                                    nelem -= 1;
                                    if fresh146 == 0 {
                                        break;
                                    }
                                    elem -= 1;
                                    if *((lParse.Nodes[theParams[0]]).value.undef)
                                        .offset(elem as isize)
                                        == 0
                                    {
                                        if valInit != 0 {
                                            valInit = 0;
                                            minVal_0 =
                                                *((lParse.Nodes[theParams[0]]).value.data.dblptr)
                                                    .offset(elem as isize);
                                        } else {
                                            minVal_0 = if minVal_0
                                                < *((lParse.Nodes[theParams[0]]).value.data.dblptr)
                                                    .offset(elem as isize)
                                            {
                                                minVal_0
                                            } else {
                                                *((lParse.Nodes[theParams[0]]).value.data.dblptr)
                                                    .offset(elem as isize)
                                            };
                                        }
                                        *((lParse.Nodes[this_node_idx]).value.undef)
                                            .offset(row as isize) = 0;
                                    }
                                }
                                *((lParse.Nodes[this_node_idx]).value.data.dblptr)
                                    .offset(row as isize) = minVal_0;
                            }
                        } else if (lParse.Nodes[this_node_idx]).ntype
                            == fits_parser_yytokentype::BITSTR as c_int
                        {
                            let mut minVal_1: c_char = 0;
                            loop {
                                let fresh147 = row;
                                row -= 1;
                                if fresh147 == 0 {
                                    break;
                                }
                                let mut sptr1_0: *mut c_char =
                                    *((lParse.Nodes[theParams[0]]).value.data.strptr)
                                        .offset(row as isize);
                                minVal_1 = b'1' as c_char;
                                while *sptr1_0 != 0 {
                                    if c_int::from(*sptr1_0) == '0' as i32 {
                                        minVal_1 = b'0' as c_char;
                                    }
                                    sptr1_0 = sptr1_0.offset(1);
                                }
                                *(*((lParse.Nodes[this_node_idx]).value.data.strptr)
                                    .offset(row as isize))
                                .offset(0) = minVal_1;
                                *(*((lParse.Nodes[this_node_idx]).value.data.strptr)
                                    .offset(row as isize))
                                .offset(1) = 0;
                            }
                        }
                    }
                    1023 => {
                        if (lParse.Nodes[this_node_idx]).ntype
                            == fits_parser_yytokentype::LONG as c_int
                        {
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
                                        if vector[i as usize] > 1 {
                                            pVals[i as usize].data.lng = *((lParse.Nodes
                                                [theParams[i as usize]])
                                                .value
                                                .data
                                                .lngptr)
                                                .offset(elem as isize);
                                            pNull[i as usize] = *((lParse.Nodes
                                                [theParams[i as usize]])
                                                .value
                                                .undef)
                                                .offset(elem as isize);
                                        } else if vector[i as usize] != 0 {
                                            pVals[i as usize].data.lng = *((lParse.Nodes
                                                [theParams[i as usize]])
                                                .value
                                                .data
                                                .lngptr)
                                                .offset(row as isize);
                                            pNull[i as usize] = *((lParse.Nodes
                                                [theParams[i as usize]])
                                                .value
                                                .undef)
                                                .offset(row as isize);
                                        }
                                    }
                                    if c_int::from(pNull[0]) != 0 && c_int::from(pNull[1]) != 0 {
                                        *((lParse.Nodes[this_node_idx]).value.undef)
                                            .offset(elem as isize) = 1;
                                        *((lParse.Nodes[this_node_idx]).value.data.lngptr)
                                            .offset(elem as isize) = 0;
                                    } else if pNull[0] != 0 {
                                        *((lParse.Nodes[this_node_idx]).value.undef)
                                            .offset(elem as isize) = 0;
                                        *((lParse.Nodes[this_node_idx]).value.data.lngptr)
                                            .offset(elem as isize) = pVals[1].data.lng;
                                    } else if pNull[1] != 0 {
                                        *((lParse.Nodes[this_node_idx]).value.undef)
                                            .offset(elem as isize) = 0;
                                        *((lParse.Nodes[this_node_idx]).value.data.lngptr)
                                            .offset(elem as isize) = pVals[0].data.lng;
                                    } else {
                                        *((lParse.Nodes[this_node_idx]).value.undef)
                                            .offset(elem as isize) = 0;
                                        *((lParse.Nodes[this_node_idx]).value.data.lngptr)
                                            .offset(elem as isize) =
                                            if pVals[0].data.lng < pVals[1].data.lng {
                                                pVals[0].data.lng
                                            } else {
                                                pVals[1].data.lng
                                            };
                                    }
                                }
                            }
                        } else if (lParse.Nodes[this_node_idx]).ntype
                            == fits_parser_yytokentype::DOUBLE as c_int
                        {
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
                                        if vector[i as usize] > 1 {
                                            pVals[i as usize].data.dbl = *((lParse.Nodes
                                                [theParams[i as usize]])
                                                .value
                                                .data
                                                .dblptr)
                                                .offset(elem as isize);
                                            pNull[i as usize] = *((lParse.Nodes
                                                [theParams[i as usize]])
                                                .value
                                                .undef)
                                                .offset(elem as isize);
                                        } else if vector[i as usize] != 0 {
                                            pVals[i as usize].data.dbl = *((lParse.Nodes
                                                [theParams[i as usize]])
                                                .value
                                                .data
                                                .dblptr)
                                                .offset(row as isize);
                                            pNull[i as usize] = *((lParse.Nodes
                                                [theParams[i as usize]])
                                                .value
                                                .undef)
                                                .offset(row as isize);
                                        }
                                    }
                                    if c_int::from(pNull[0]) != 0 && c_int::from(pNull[1]) != 0 {
                                        *((lParse.Nodes[this_node_idx]).value.undef)
                                            .offset(elem as isize) = 1;
                                        *((lParse.Nodes[this_node_idx]).value.data.dblptr)
                                            .offset(elem as isize) = 0.0;
                                    } else if pNull[0] != 0 {
                                        *((lParse.Nodes[this_node_idx]).value.undef)
                                            .offset(elem as isize) = 0;
                                        *((lParse.Nodes[this_node_idx]).value.data.dblptr)
                                            .offset(elem as isize) = pVals[1].data.dbl;
                                    } else if pNull[1] != 0 {
                                        *((lParse.Nodes[this_node_idx]).value.undef)
                                            .offset(elem as isize) = 0;
                                        *((lParse.Nodes[this_node_idx]).value.data.dblptr)
                                            .offset(elem as isize) = pVals[0].data.dbl;
                                    } else {
                                        *((lParse.Nodes[this_node_idx]).value.undef)
                                            .offset(elem as isize) = 0;
                                        *((lParse.Nodes[this_node_idx]).value.data.dblptr)
                                            .offset(elem as isize) =
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
                        elem = row * (lParse.Nodes[theParams[0]]).value.nelem;
                        if (lParse.Nodes[this_node_idx]).ntype
                            == fits_parser_yytokentype::LONG as c_int
                        {
                            let mut maxVal: c_long = 0;
                            loop {
                                let fresh154 = row;
                                row -= 1;
                                if fresh154 == 0 {
                                    break;
                                }
                                valInit = 1;
                                *((lParse.Nodes[this_node_idx]).value.undef).offset(row as isize) =
                                    1;
                                nelem = (lParse.Nodes[theParams[0]]).value.nelem;
                                loop {
                                    let fresh155 = nelem;
                                    nelem -= 1;
                                    if fresh155 == 0 {
                                        break;
                                    }
                                    elem -= 1;
                                    if *((lParse.Nodes[theParams[0]]).value.undef)
                                        .offset(elem as isize)
                                        == 0
                                    {
                                        if valInit != 0 {
                                            valInit = 0;
                                            maxVal =
                                                *((lParse.Nodes[theParams[0]]).value.data.lngptr)
                                                    .offset(elem as isize);
                                        } else {
                                            maxVal = if maxVal
                                                > *((lParse.Nodes[theParams[0]]).value.data.lngptr)
                                                    .offset(elem as isize)
                                            {
                                                maxVal
                                            } else {
                                                *((lParse.Nodes[theParams[0]]).value.data.lngptr)
                                                    .offset(elem as isize)
                                            };
                                        }
                                        *((lParse.Nodes[this_node_idx]).value.undef)
                                            .offset(row as isize) = 0;
                                    }
                                }
                                *((lParse.Nodes[this_node_idx]).value.data.lngptr)
                                    .offset(row as isize) = maxVal;
                            }
                        } else if (lParse.Nodes[this_node_idx]).ntype
                            == fits_parser_yytokentype::DOUBLE as c_int
                        {
                            let mut maxVal_0: c_double = 0.0;
                            loop {
                                let fresh156 = row;
                                row -= 1;
                                if fresh156 == 0 {
                                    break;
                                }
                                valInit = 1;
                                *((lParse.Nodes[this_node_idx]).value.undef).offset(row as isize) =
                                    1;
                                nelem = (lParse.Nodes[theParams[0]]).value.nelem;
                                loop {
                                    let fresh157 = nelem;
                                    nelem -= 1;
                                    if fresh157 == 0 {
                                        break;
                                    }
                                    elem -= 1;
                                    if *((lParse.Nodes[theParams[0]]).value.undef)
                                        .offset(elem as isize)
                                        == 0
                                    {
                                        if valInit != 0 {
                                            valInit = 0;
                                            maxVal_0 =
                                                *((lParse.Nodes[theParams[0]]).value.data.dblptr)
                                                    .offset(elem as isize);
                                        } else {
                                            maxVal_0 = if maxVal_0
                                                > *((lParse.Nodes[theParams[0]]).value.data.dblptr)
                                                    .offset(elem as isize)
                                            {
                                                maxVal_0
                                            } else {
                                                *((lParse.Nodes[theParams[0]]).value.data.dblptr)
                                                    .offset(elem as isize)
                                            };
                                        }
                                        *((lParse.Nodes[this_node_idx]).value.undef)
                                            .offset(row as isize) = 0;
                                    }
                                }
                                *((lParse.Nodes[this_node_idx]).value.data.dblptr)
                                    .offset(row as isize) = maxVal_0;
                            }
                        } else if (lParse.Nodes[this_node_idx]).ntype
                            == fits_parser_yytokentype::BITSTR as c_int
                        {
                            let mut maxVal_1: c_char = 0;
                            loop {
                                let fresh158 = row;
                                row -= 1;
                                if fresh158 == 0 {
                                    break;
                                }
                                let mut sptr1_1: *mut c_char =
                                    *((lParse.Nodes[theParams[0]]).value.data.strptr)
                                        .offset(row as isize);
                                maxVal_1 = b'0' as c_char;
                                while *sptr1_1 != 0 {
                                    if c_int::from(*sptr1_1) == '1' as i32 {
                                        maxVal_1 = b'1' as c_char;
                                    }
                                    sptr1_1 = sptr1_1.offset(1);
                                }
                                *(*((lParse.Nodes[this_node_idx]).value.data.strptr)
                                    .offset(row as isize))
                                .offset(0) = maxVal_1;
                                *(*((lParse.Nodes[this_node_idx]).value.data.strptr)
                                    .offset(row as isize))
                                .offset(1) = 0;
                            }
                        }
                    }
                    1025 => {
                        if (lParse.Nodes[this_node_idx]).ntype
                            == fits_parser_yytokentype::LONG as c_int
                        {
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
                                        if vector[i as usize] > 1 {
                                            pVals[i as usize].data.lng = *((lParse.Nodes
                                                [theParams[i as usize]])
                                                .value
                                                .data
                                                .lngptr)
                                                .offset(elem as isize);
                                            pNull[i as usize] = *((lParse.Nodes
                                                [theParams[i as usize]])
                                                .value
                                                .undef)
                                                .offset(elem as isize);
                                        } else if vector[i as usize] != 0 {
                                            pVals[i as usize].data.lng = *((lParse.Nodes
                                                [theParams[i as usize]])
                                                .value
                                                .data
                                                .lngptr)
                                                .offset(row as isize);
                                            pNull[i as usize] = *((lParse.Nodes
                                                [theParams[i as usize]])
                                                .value
                                                .undef)
                                                .offset(row as isize);
                                        }
                                    }
                                    if c_int::from(pNull[0]) != 0 && c_int::from(pNull[1]) != 0 {
                                        *((lParse.Nodes[this_node_idx]).value.undef)
                                            .offset(elem as isize) = 1;
                                        *((lParse.Nodes[this_node_idx]).value.data.lngptr)
                                            .offset(elem as isize) = 0;
                                    } else if pNull[0] != 0 {
                                        *((lParse.Nodes[this_node_idx]).value.undef)
                                            .offset(elem as isize) = 0;
                                        *((lParse.Nodes[this_node_idx]).value.data.lngptr)
                                            .offset(elem as isize) = pVals[1].data.lng;
                                    } else if pNull[1] != 0 {
                                        *((lParse.Nodes[this_node_idx]).value.undef)
                                            .offset(elem as isize) = 0;
                                        *((lParse.Nodes[this_node_idx]).value.data.lngptr)
                                            .offset(elem as isize) = pVals[0].data.lng;
                                    } else {
                                        *((lParse.Nodes[this_node_idx]).value.undef)
                                            .offset(elem as isize) = 0;
                                        *((lParse.Nodes[this_node_idx]).value.data.lngptr)
                                            .offset(elem as isize) =
                                            if pVals[0].data.lng > pVals[1].data.lng {
                                                pVals[0].data.lng
                                            } else {
                                                pVals[1].data.lng
                                            };
                                    }
                                }
                            }
                        } else if (lParse.Nodes[this_node_idx]).ntype
                            == fits_parser_yytokentype::DOUBLE as c_int
                        {
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
                                        if vector[i as usize] > 1 {
                                            pVals[i as usize].data.dbl = *((lParse.Nodes
                                                [theParams[i as usize]])
                                                .value
                                                .data
                                                .dblptr)
                                                .offset(elem as isize);
                                            pNull[i as usize] = *((lParse.Nodes
                                                [theParams[i as usize]])
                                                .value
                                                .undef)
                                                .offset(elem as isize);
                                        } else if vector[i as usize] != 0 {
                                            pVals[i as usize].data.dbl = *((lParse.Nodes
                                                [theParams[i as usize]])
                                                .value
                                                .data
                                                .dblptr)
                                                .offset(row as isize);
                                            pNull[i as usize] = *((lParse.Nodes
                                                [theParams[i as usize]])
                                                .value
                                                .undef)
                                                .offset(row as isize);
                                        }
                                    }
                                    if c_int::from(pNull[0]) != 0 && c_int::from(pNull[1]) != 0 {
                                        *((lParse.Nodes[this_node_idx]).value.undef)
                                            .offset(elem as isize) = 1;
                                        *((lParse.Nodes[this_node_idx]).value.data.dblptr)
                                            .offset(elem as isize) = 0.0;
                                    } else if pNull[0] != 0 {
                                        *((lParse.Nodes[this_node_idx]).value.undef)
                                            .offset(elem as isize) = 0;
                                        *((lParse.Nodes[this_node_idx]).value.data.dblptr)
                                            .offset(elem as isize) = pVals[1].data.dbl;
                                    } else if pNull[1] != 0 {
                                        *((lParse.Nodes[this_node_idx]).value.undef)
                                            .offset(elem as isize) = 0;
                                        *((lParse.Nodes[this_node_idx]).value.data.dblptr)
                                            .offset(elem as isize) = pVals[0].data.dbl;
                                    } else {
                                        *((lParse.Nodes[this_node_idx]).value.undef)
                                            .offset(elem as isize) = 0;
                                        *((lParse.Nodes[this_node_idx]).value.data.dblptr)
                                            .offset(elem as isize) =
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
                                if vector[i as usize] > 1 {
                                    pVals[i as usize].data.dbl =
                                        *((lParse.Nodes[theParams[i as usize]]).value.data.dblptr)
                                            .offset(elem as isize);
                                    pNull[i as usize] =
                                        *((lParse.Nodes[theParams[i as usize]]).value.undef)
                                            .offset(elem as isize);
                                } else if vector[i as usize] != 0 {
                                    pVals[i as usize].data.dbl =
                                        *((lParse.Nodes[theParams[i as usize]]).value.data.dblptr)
                                            .offset(row as isize);
                                    pNull[i as usize] =
                                        *((lParse.Nodes[theParams[i as usize]]).value.undef)
                                            .offset(row as isize);
                                }
                            }
                            let fresh168 = &mut *((lParse.Nodes[this_node_idx]).value.undef)
                                .offset(elem as isize);
                            *fresh168 = c_int::from(
                                c_int::from(pNull[0]) != 0
                                    || c_int::from(pNull[1]) != 0
                                    || c_int::from(pNull[2]) != 0,
                            ) as c_char;
                            if *fresh168 == 0 {
                                *((lParse.Nodes[this_node_idx]).value.data.logptr)
                                    .offset(elem as isize) =
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
                                if vector[i as usize] > 1 {
                                    pVals[i as usize].data.dbl =
                                        *((lParse.Nodes[theParams[i as usize]]).value.data.dblptr)
                                            .offset(elem as isize);
                                    pNull[i as usize] =
                                        *((lParse.Nodes[theParams[i as usize]]).value.undef)
                                            .offset(elem as isize);
                                } else if vector[i as usize] != 0 {
                                    pVals[i as usize].data.dbl =
                                        *((lParse.Nodes[theParams[i as usize]]).value.data.dblptr)
                                            .offset(row as isize);
                                    pNull[i as usize] =
                                        *((lParse.Nodes[theParams[i as usize]]).value.undef)
                                            .offset(row as isize);
                                }
                            }
                            let fresh172 = &mut *((lParse.Nodes[this_node_idx]).value.undef)
                                .offset(elem as isize);
                            *fresh172 = c_int::from(
                                c_int::from(pNull[0]) != 0
                                    || c_int::from(pNull[1]) != 0
                                    || c_int::from(pNull[2]) != 0
                                    || c_int::from(pNull[3]) != 0
                                    || c_int::from(pNull[4]) != 0,
                            ) as c_char;
                            if *fresh172 == 0 {
                                *((lParse.Nodes[this_node_idx]).value.data.logptr)
                                    .offset(elem as isize) = circle(
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
                                if vector[i as usize] > 1 {
                                    pVals[i as usize].data.dbl =
                                        *((lParse.Nodes[theParams[i as usize]]).value.data.dblptr)
                                            .offset(elem as isize);
                                    pNull[i as usize] =
                                        *((lParse.Nodes[theParams[i as usize]]).value.undef)
                                            .offset(elem as isize);
                                } else if vector[i as usize] != 0 {
                                    pVals[i as usize].data.dbl =
                                        *((lParse.Nodes[theParams[i as usize]]).value.data.dblptr)
                                            .offset(row as isize);
                                    pNull[i as usize] =
                                        *((lParse.Nodes[theParams[i as usize]]).value.undef)
                                            .offset(row as isize);
                                }
                            }
                            let fresh176 = &mut *((lParse.Nodes[this_node_idx]).value.undef)
                                .offset(elem as isize);
                            *fresh176 = c_int::from(
                                c_int::from(pNull[0]) != 0
                                    || c_int::from(pNull[1]) != 0
                                    || c_int::from(pNull[2]) != 0
                                    || c_int::from(pNull[3]) != 0
                                    || c_int::from(pNull[4]) != 0
                                    || c_int::from(pNull[5]) != 0
                                    || c_int::from(pNull[6]) != 0,
                            ) as c_char;
                            if *fresh176 == 0 {
                                *((lParse.Nodes[this_node_idx]).value.data.logptr)
                                    .offset(elem as isize) = saobox(
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
                                if vector[i as usize] > 1 {
                                    pVals[i as usize].data.dbl =
                                        *((lParse.Nodes[theParams[i as usize]]).value.data.dblptr)
                                            .offset(elem as isize);
                                    pNull[i as usize] =
                                        *((lParse.Nodes[theParams[i as usize]]).value.undef)
                                            .offset(elem as isize);
                                } else if vector[i as usize] != 0 {
                                    pVals[i as usize].data.dbl =
                                        *((lParse.Nodes[theParams[i as usize]]).value.data.dblptr)
                                            .offset(row as isize);
                                    pNull[i as usize] =
                                        *((lParse.Nodes[theParams[i as usize]]).value.undef)
                                            .offset(row as isize);
                                }
                            }
                            let fresh180 = &mut *((lParse.Nodes[this_node_idx]).value.undef)
                                .offset(elem as isize);
                            *fresh180 = c_int::from(
                                c_int::from(pNull[0]) != 0
                                    || c_int::from(pNull[1]) != 0
                                    || c_int::from(pNull[2]) != 0
                                    || c_int::from(pNull[3]) != 0
                                    || c_int::from(pNull[4]) != 0
                                    || c_int::from(pNull[5]) != 0
                                    || c_int::from(pNull[6]) != 0,
                            ) as c_char;
                            if *fresh180 == 0 {
                                *((lParse.Nodes[this_node_idx]).value.data.logptr)
                                    .offset(elem as isize) = ellipse(
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
                    1034 => match (lParse.Nodes[this_node_idx]).ntype.into() {
                        fits_parser_yytokentype::BOOLEAN => loop {
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
                                if vector[2] > 1 {
                                    pVals[2].data.log =
                                        *((lParse.Nodes[theParams[2]]).value.data.logptr)
                                            .offset(elem as isize);
                                    pNull[2] = *((lParse.Nodes[theParams[2]]).value.undef)
                                        .offset(elem as isize);
                                } else if vector[2] != 0 {
                                    pVals[2].data.log =
                                        *((lParse.Nodes[theParams[2]]).value.data.logptr)
                                            .offset(row as isize);
                                    pNull[2] = *((lParse.Nodes[theParams[2]]).value.undef)
                                        .offset(row as isize);
                                }
                                i = 2;
                                loop {
                                    let fresh183 = i;
                                    i -= 1;
                                    if fresh183 == 0 {
                                        break;
                                    }
                                    if vector[i as usize] > 1 {
                                        pVals[i as usize].data.log = *((lParse.Nodes
                                            [theParams[i as usize]])
                                            .value
                                            .data
                                            .logptr)
                                            .offset(elem as isize);
                                        pNull[i as usize] =
                                            *((lParse.Nodes[theParams[i as usize]]).value.undef)
                                                .offset(elem as isize);
                                    } else if vector[i as usize] != 0 {
                                        pVals[i as usize].data.log = *((lParse.Nodes
                                            [theParams[i as usize]])
                                            .value
                                            .data
                                            .logptr)
                                            .offset(row as isize);
                                        pNull[i as usize] =
                                            *((lParse.Nodes[theParams[i as usize]]).value.undef)
                                                .offset(row as isize);
                                    }
                                }
                                let fresh184 = &mut *((lParse.Nodes[this_node_idx]).value.undef)
                                    .offset(elem as isize);
                                *fresh184 = pNull[2];
                                if *fresh184 == 0 {
                                    if pVals[2].data.log != 0 {
                                        *((lParse.Nodes[this_node_idx]).value.data.logptr)
                                            .offset(elem as isize) = pVals[0].data.log;
                                        *((lParse.Nodes[this_node_idx]).value.undef)
                                            .offset(elem as isize) = pNull[0];
                                    } else {
                                        *((lParse.Nodes[this_node_idx]).value.data.logptr)
                                            .offset(elem as isize) = pVals[1].data.log;
                                        *((lParse.Nodes[this_node_idx]).value.undef)
                                            .offset(elem as isize) = pNull[1];
                                    }
                                }
                            }
                        },
                        fits_parser_yytokentype::LONG => loop {
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
                                if vector[2] > 1 {
                                    pVals[2].data.log =
                                        *((lParse.Nodes[theParams[2]]).value.data.logptr)
                                            .offset(elem as isize);
                                    pNull[2] = *((lParse.Nodes[theParams[2]]).value.undef)
                                        .offset(elem as isize);
                                } else if vector[2] != 0 {
                                    pVals[2].data.log =
                                        *((lParse.Nodes[theParams[2]]).value.data.logptr)
                                            .offset(row as isize);
                                    pNull[2] = *((lParse.Nodes[theParams[2]]).value.undef)
                                        .offset(row as isize);
                                }
                                i = 2;
                                loop {
                                    let fresh187 = i;
                                    i -= 1;
                                    if fresh187 == 0 {
                                        break;
                                    }
                                    if vector[i as usize] > 1 {
                                        pVals[i as usize].data.lng = *((lParse.Nodes
                                            [theParams[i as usize]])
                                            .value
                                            .data
                                            .lngptr)
                                            .offset(elem as isize);
                                        pNull[i as usize] =
                                            *((lParse.Nodes[theParams[i as usize]]).value.undef)
                                                .offset(elem as isize);
                                    } else if vector[i as usize] != 0 {
                                        pVals[i as usize].data.lng = *((lParse.Nodes
                                            [theParams[i as usize]])
                                            .value
                                            .data
                                            .lngptr)
                                            .offset(row as isize);
                                        pNull[i as usize] =
                                            *((lParse.Nodes[theParams[i as usize]]).value.undef)
                                                .offset(row as isize);
                                    }
                                }
                                let fresh188 = &mut *((lParse.Nodes[this_node_idx]).value.undef)
                                    .offset(elem as isize);
                                *fresh188 = pNull[2];
                                if *fresh188 == 0 {
                                    if pVals[2].data.log != 0 {
                                        *((lParse.Nodes[this_node_idx]).value.data.lngptr)
                                            .offset(elem as isize) = pVals[0].data.lng;
                                        *((lParse.Nodes[this_node_idx]).value.undef)
                                            .offset(elem as isize) = pNull[0];
                                    } else {
                                        *((lParse.Nodes[this_node_idx]).value.data.lngptr)
                                            .offset(elem as isize) = pVals[1].data.lng;
                                        *((lParse.Nodes[this_node_idx]).value.undef)
                                            .offset(elem as isize) = pNull[1];
                                    }
                                }
                            }
                        },
                        fits_parser_yytokentype::DOUBLE => loop {
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
                                if vector[2] > 1 {
                                    pVals[2].data.log =
                                        *((lParse.Nodes[theParams[2]]).value.data.logptr)
                                            .offset(elem as isize);
                                    pNull[2] = *((lParse.Nodes[theParams[2]]).value.undef)
                                        .offset(elem as isize);
                                } else if vector[2] != 0 {
                                    pVals[2].data.log =
                                        *((lParse.Nodes[theParams[2]]).value.data.logptr)
                                            .offset(row as isize);
                                    pNull[2] = *((lParse.Nodes[theParams[2]]).value.undef)
                                        .offset(row as isize);
                                }
                                i = 2;
                                loop {
                                    let fresh191 = i;
                                    i -= 1;
                                    if fresh191 == 0 {
                                        break;
                                    }
                                    if vector[i as usize] > 1 {
                                        pVals[i as usize].data.dbl = *((lParse.Nodes
                                            [theParams[i as usize]])
                                            .value
                                            .data
                                            .dblptr)
                                            .offset(elem as isize);
                                        pNull[i as usize] =
                                            *((lParse.Nodes[theParams[i as usize]]).value.undef)
                                                .offset(elem as isize);
                                    } else if vector[i as usize] != 0 {
                                        pVals[i as usize].data.dbl = *((lParse.Nodes
                                            [theParams[i as usize]])
                                            .value
                                            .data
                                            .dblptr)
                                            .offset(row as isize);
                                        pNull[i as usize] =
                                            *((lParse.Nodes[theParams[i as usize]]).value.undef)
                                                .offset(row as isize);
                                    }
                                }
                                let fresh192 = &mut *((lParse.Nodes[this_node_idx]).value.undef)
                                    .offset(elem as isize);
                                *fresh192 = pNull[2];
                                if *fresh192 == 0 {
                                    if pVals[2].data.log != 0 {
                                        *((lParse.Nodes[this_node_idx]).value.data.dblptr)
                                            .offset(elem as isize) = pVals[0].data.dbl;
                                        *((lParse.Nodes[this_node_idx]).value.undef)
                                            .offset(elem as isize) = pNull[0];
                                    } else {
                                        *((lParse.Nodes[this_node_idx]).value.data.dblptr)
                                            .offset(elem as isize) = pVals[1].data.dbl;
                                        *((lParse.Nodes[this_node_idx]).value.undef)
                                            .offset(elem as isize) = pNull[1];
                                    }
                                }
                            }
                        },
                        fits_parser_yytokentype::STRING => loop {
                            let fresh193 = row;
                            row -= 1;
                            if fresh193 == 0 {
                                break;
                            }
                            if vector[2] != 0 {
                                pVals[2].data.log =
                                    *((lParse.Nodes[theParams[2]]).value.data.logptr)
                                        .offset(row as isize);
                                pNull[2] = *((lParse.Nodes[theParams[2]]).value.undef)
                                    .offset(row as isize);
                            }
                            i = 2;
                            loop {
                                let fresh194 = i;
                                i -= 1;
                                if fresh194 == 0 {
                                    break;
                                }
                                if vector[i as usize] != 0 {
                                    strcpy(
                                        (pVals[i as usize].data.astr).as_mut_ptr(),
                                        *((lParse.Nodes[theParams[i as usize]]).value.data.strptr)
                                            .offset(row as isize),
                                    );
                                    pNull[i as usize] =
                                        *((lParse.Nodes[theParams[i as usize]]).value.undef)
                                            .offset(row as isize);
                                }
                            }
                            let fresh195 = &mut *((lParse.Nodes[this_node_idx]).value.undef)
                                .offset(row as isize);
                            *fresh195 = pNull[2];
                            if *fresh195 == 0 {
                                if pVals[2].data.log != 0 {
                                    strcpy(
                                        *((lParse.Nodes[this_node_idx]).value.data.strptr)
                                            .offset(row as isize),
                                        (pVals[0].data.astr).as_mut_ptr(),
                                    );
                                    *((lParse.Nodes[this_node_idx]).value.undef)
                                        .offset(row as isize) = pNull[0];
                                } else {
                                    strcpy(
                                        *((lParse.Nodes[this_node_idx]).value.data.strptr)
                                            .offset(row as isize),
                                        (pVals[1].data.astr).as_mut_ptr(),
                                    );
                                    *((lParse.Nodes[this_node_idx]).value.undef)
                                        .offset(row as isize) = pNull[1];
                                }
                            } else {
                                *(*((lParse.Nodes[this_node_idx]).value.data.strptr)
                                    .offset(row as isize))
                                .offset(0) = 0;
                            }
                        },
                        _ => {}
                    },
                    STRMID_FCT => {
                        let strconst: c_int =
                            c_int::from((lParse.Nodes[theParams[0]]).operation == CONST_OP);
                        let posconst: c_int =
                            c_int::from((lParse.Nodes[theParams[1]]).operation == CONST_OP);
                        let lenconst: c_int =
                            c_int::from((lParse.Nodes[theParams[2]]).operation == CONST_OP);
                        let dest_len: c_int = (lParse.Nodes[this_node_idx]).value.nelem as c_int;
                        let mut src_len: c_int = (lParse.Nodes[theParams[0]]).value.nelem as c_int;
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
                                pos = (lParse.Nodes[theParams[1]]).value.data.lng as c_int;
                            } else {
                                pos = *((lParse.Nodes[theParams[1]]).value.data.lngptr)
                                    .offset(row as isize)
                                    as c_int;
                                if *((lParse.Nodes[theParams[1]]).value.undef).offset(row as isize)
                                    != 0
                                {
                                    undef = 1;
                                }
                            }
                            if strconst != 0 {
                                str = ((lParse.Nodes[theParams[0]]).value.data.astr).as_mut_ptr();
                                if src_len == 0 {
                                    src_len = strlen(str) as c_int;
                                }
                            } else {
                                str = *((lParse.Nodes[theParams[0]]).value.data.strptr)
                                    .offset(row as isize);
                                if *((lParse.Nodes[theParams[0]]).value.undef).offset(row as isize)
                                    != 0
                                {
                                    undef = 1;
                                }
                            }
                            if lenconst != 0 {
                                len = dest_len;
                            } else {
                                len = *((lParse.Nodes[theParams[2]]).value.data.lngptr)
                                    .offset(row as isize)
                                    as c_int;
                                if *((lParse.Nodes[theParams[2]]).value.undef).offset(row as isize)
                                    != 0
                                {
                                    undef = 1;
                                }
                            }
                            *(*((lParse.Nodes[this_node_idx]).value.data.strptr)
                                .offset(row as isize))
                            .offset(0) = 0;
                            if pos == 0 {
                                undef = 1;
                            }
                            if undef == 0
                                && cstrmid(
                                    lParse,
                                    *((lParse.Nodes[this_node_idx]).value.data.strptr)
                                        .offset(row as isize),
                                    len,
                                    str,
                                    src_len,
                                    pos,
                                ) < 0
                            {
                                break;
                            }
                            *((lParse.Nodes[this_node_idx]).value.undef).offset(row as isize) =
                                undef as c_char;
                        }
                    }
                    1045 => {
                        let const1: c_int =
                            c_int::from((lParse.Nodes[theParams[0]]).operation == CONST_OP);
                        let const2: c_int =
                            c_int::from((lParse.Nodes[theParams[1]]).operation == CONST_OP);
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
                                str1 = ((lParse.Nodes[theParams[0]]).value.data.astr).as_mut_ptr();
                            } else {
                                str1 = *((lParse.Nodes[theParams[0]]).value.data.strptr)
                                    .offset(row as isize);
                                if *((lParse.Nodes[theParams[0]]).value.undef).offset(row as isize)
                                    != 0
                                {
                                    undef_0 = 1;
                                }
                            }
                            if const2 != 0 {
                                str2 = ((lParse.Nodes[theParams[1]]).value.data.astr).as_mut_ptr();
                            } else {
                                str2 = *((lParse.Nodes[theParams[1]]).value.data.strptr)
                                    .offset(row as isize);
                                if *((lParse.Nodes[theParams[1]]).value.undef).offset(row as isize)
                                    != 0
                                {
                                    undef_0 = 1;
                                }
                            }
                            *((lParse.Nodes[this_node_idx]).value.data.lngptr)
                                .offset(row as isize) = 0;
                            if undef_0 == 0 {
                                let res_0: *mut c_char = strstr(str1, str2);
                                if res_0.is_null() {
                                    undef_0 = 1;
                                    *((lParse.Nodes[this_node_idx]).value.data.lngptr)
                                        .offset(row as isize) = 0;
                                } else {
                                    *((lParse.Nodes[this_node_idx]).value.data.lngptr)
                                        .offset(row as isize) =
                                        res_0.offset_from(str1) as c_long + 1;
                                }
                            }
                            *((lParse.Nodes[this_node_idx]).value.undef).offset(row as isize) =
                                undef_0 as c_char;
                        }
                    }
                    _ => {}
                }
            }
        }
        i = (lParse.Nodes[this_node_idx]).nSubNodes;
        loop {
            let fresh198 = i;
            i -= 1;
            if fresh198 == 0 {
                break;
            }
            if (lParse.Nodes[theParams[i as usize]]).operation > 0 {
                free((lParse.Nodes[theParams[i as usize]]).value.data.ptr);
            }
        }
    }
}

fn Do_Deref(lParse: &mut ParseData, this_node_idx: usize) {
    unsafe {
        let mut theVar: &mut Node;
        let mut theDims: [usize; MAXDIMS as usize] = [0; MAXDIMS as usize];
        let mut isConst: [c_int; MAXDIMS as usize] = [0; MAXDIMS as usize];
        let mut allConst: c_int = 0;
        let mut dimVals: [c_long; MAXDIMS as usize] = [0; MAXDIMS as usize];
        let mut i: c_int = 0;
        let mut nDims: c_int = 0;
        let mut row: c_long = 0;
        let mut elem: c_long = 0;
        let mut dsize: c_long = 0;

        let theVar = (lParse.Nodes[this_node_idx]).SubNodes[0];
        nDims = (lParse.Nodes[this_node_idx]).nSubNodes - 1;
        i = nDims;
        allConst = 1;

        loop {
            if i == 0 {
                i -= 1;
                break;
            }
            i -= 1;

            theDims[i as usize] = (lParse.Nodes[this_node_idx]).SubNodes[(i + 1) as usize] as usize;
            isConst[i as usize] =
                c_int::from((lParse.Nodes[theDims[i as usize]]).operation == CONST_OP);
            if isConst[i as usize] != 0 {
                dimVals[i as usize] = (lParse.Nodes[theDims[i as usize]]).value.data.lng;
            } else {
                allConst = 0;
            }
        }

        if (lParse.Nodes[this_node_idx]).ntype == fits_parser_yytokentype::DOUBLE as c_int {
            dsize = ::core::mem::size_of::<c_double>() as c_ulong as c_long;
        } else if (lParse.Nodes[this_node_idx]).ntype == fits_parser_yytokentype::LONG as c_int {
            dsize = ::core::mem::size_of::<c_long>() as c_ulong as c_long;
        } else if (lParse.Nodes[this_node_idx]).ntype == fits_parser_yytokentype::BOOLEAN as c_int {
            dsize = ::core::mem::size_of::<c_char>() as c_ulong as c_long;
        } else {
            dsize = 0;
        }

        Allocate_Ptrs(lParse, this_node_idx);

        if lParse.status == 0 {
            if allConst != 0 && (lParse.Nodes[theVar]).value.naxis == nDims {
                /* Dereference completely using constant indices */
                elem = 0;
                i = nDims;
                loop {
                    if i == 0 {
                        i -= 1;
                        break;
                    }
                    i -= 1;

                    if dimVals[i as usize] < 1
                        || dimVals[i as usize] > (lParse.Nodes[theVar]).value.naxes[i as usize]
                    {
                        break;
                    }

                    elem = (lParse.Nodes[theVar]).value.naxes[i as usize] * elem
                        + dimVals[i as usize]
                        - 1;
                }

                if i < 0 {
                    row = 0;
                    while row < lParse.nRows {
                        if (lParse.Nodes[this_node_idx]).ntype
                            == fits_parser_yytokentype::STRING as c_int
                        {
                            *((lParse.Nodes[this_node_idx]).value.undef).offset(row as isize) =
                                *((lParse.Nodes[theVar]).value.undef).offset(row as isize);
                        } else if (lParse.Nodes[this_node_idx]).ntype
                            != fits_parser_yytokentype::BITSTR as c_int
                        {
                            *((lParse.Nodes[this_node_idx]).value.undef).offset(row as isize) =
                                *((lParse.Nodes[theVar]).value.undef).offset(elem as isize);
                            /* Dummy - BITSTRs do not have undefs */
                        }
                        if (lParse.Nodes[this_node_idx]).ntype
                            == fits_parser_yytokentype::DOUBLE as c_int
                        {
                            *((lParse.Nodes[this_node_idx]).value.data.dblptr)
                                .offset(row as isize) =
                                *((lParse.Nodes[theVar]).value.data.dblptr).offset(elem as isize);
                        } else if (lParse.Nodes[this_node_idx]).ntype
                            == fits_parser_yytokentype::LONG as c_int
                        {
                            *((lParse.Nodes[this_node_idx]).value.data.lngptr)
                                .offset(row as isize) =
                                *((lParse.Nodes[theVar]).value.data.lngptr).offset(elem as isize);
                        } else if (lParse.Nodes[this_node_idx]).ntype
                            == fits_parser_yytokentype::BOOLEAN as c_int
                        {
                            *((lParse.Nodes[this_node_idx]).value.data.logptr)
                                .offset(row as isize) =
                                *((lParse.Nodes[theVar]).value.data.logptr).offset(elem as isize);
                        } else {
                            /* XXX Note, the below expression uses knowledge of
                            the layout of the string format, namely (nelem+1)
                            characters per string, followed by (nelem+1)
                            "undef" values. */

                            *(*((lParse.Nodes[this_node_idx]).value.data.strptr)
                                .offset(row as isize))
                            .offset(0) = *(*((lParse.Nodes[theVar]).value.data.strptr).offset(0))
                                .offset((elem + row) as isize);
                            *(*((lParse.Nodes[this_node_idx]).value.data.strptr)
                                .offset(row as isize))
                            .offset(1) = 0; /* Null terminate */
                        }
                        elem += (lParse.Nodes[theVar]).value.nelem;
                        row += 1;
                    }
                } else {
                    fits_parser_yyerror(lParse, cs!(c"Index out of range"));
                    free((lParse.Nodes[this_node_idx]).value.data.ptr);
                    (lParse.Nodes[this_node_idx]).value.data.ptr = ptr::null_mut();
                }
            } else if allConst != 0 && nDims == 1 {
                /* Reduce dimensions by 1, using a constant index */

                if dimVals[0] < 1
                    || dimVals[0]
                        > (lParse.Nodes[theVar]).value.naxes
                            [((lParse.Nodes[theVar]).value.naxis - 1) as usize]
                {
                    fits_parser_yyerror(lParse, cs!(c"Index out of range"));
                    free((lParse.Nodes[this_node_idx]).value.data.ptr);
                    (lParse.Nodes[this_node_idx]).value.data.ptr = ptr::null_mut();
                } else if (lParse.Nodes[this_node_idx]).ntype
                    == fits_parser_yytokentype::BITSTR as c_int
                    || (lParse.Nodes[this_node_idx]).ntype
                        == fits_parser_yytokentype::STRING as c_int
                {
                    elem = (lParse.Nodes[this_node_idx]).value.nelem * (dimVals[0] - 1);
                    row = 0;
                    while row < lParse.nRows {
                        if !((lParse.Nodes[this_node_idx]).value.undef).is_null() {
                            *((lParse.Nodes[this_node_idx]).value.undef).offset(row as isize) =
                                *((lParse.Nodes[theVar]).value.undef).offset(row as isize);
                        }
                        memcpy(
                            (*((lParse.Nodes[this_node_idx]).value.data.strptr).offset(0))
                                .offset(
                                    (row as c_ulong)
                                        .wrapping_mul(::core::mem::size_of::<c_char>() as c_ulong)
                                        .wrapping_mul(
                                            ((lParse.Nodes[this_node_idx]).value.nelem + 1)
                                                as c_ulong,
                                        ) as isize,
                                )
                                .cast::<c_void>(),
                            (*((lParse.Nodes[theVar]).value.data.strptr).offset(0)).offset(
                                (elem as c_ulong)
                                    .wrapping_mul(::core::mem::size_of::<c_char>() as c_ulong)
                                    as isize,
                            ) as *const c_void,
                            ((lParse.Nodes[this_node_idx]).value.nelem as c_ulong)
                                .wrapping_mul(::core::mem::size_of::<c_char>() as c_ulong)
                                .try_into()
                                .unwrap(),
                        );
                        *(*((lParse.Nodes[this_node_idx]).value.data.strptr)
                            .offset(row as isize))
                        .offset((lParse.Nodes[this_node_idx]).value.nelem as isize) = 0; /* Null terminate */
                        elem += (lParse.Nodes[theVar]).value.nelem + 1;
                        row += 1;
                    }
                } else {
                    elem = (lParse.Nodes[this_node_idx]).value.nelem * (dimVals[0] - 1);
                    row = 0;
                    while row < lParse.nRows {
                        memcpy(
                            ((lParse.Nodes[this_node_idx]).value.undef)
                                .offset((row * (lParse.Nodes[this_node_idx]).value.nelem) as isize)
                                .cast::<c_void>(),
                            ((lParse.Nodes[theVar]).value.undef).offset(elem as isize)
                                as *const c_void,
                            ((lParse.Nodes[this_node_idx]).value.nelem as c_ulong)
                                .wrapping_mul(::core::mem::size_of::<c_char>() as c_ulong)
                                .try_into()
                                .unwrap(),
                        );
                        memcpy(
                            (lParse.Nodes[this_node_idx])
                                .value
                                .data
                                .ptr
                                .cast::<c_char>()
                                .offset(
                                    (row * dsize * (lParse.Nodes[this_node_idx]).value.nelem)
                                        as isize,
                                )
                                .cast::<c_void>(),
                            (lParse.Nodes[theVar])
                                .value
                                .data
                                .ptr
                                .cast::<c_char>()
                                .offset((elem * dsize) as isize)
                                as *const c_void,
                            (((lParse.Nodes[this_node_idx]).value.nelem * dsize) as c_ulong)
                                .try_into()
                                .unwrap(),
                        );
                        elem += (lParse.Nodes[theVar]).value.nelem;
                        row += 1;
                    }
                }
            } else if (lParse.Nodes[theVar]).value.naxis == nDims {
                /* Dereference completely using an expression for the indices */
                row = 0;
                while row < lParse.nRows {
                    i = 0;
                    while i < nDims {
                        if isConst[i as usize] == 0 {
                            if *((lParse.Nodes[theDims[i as usize]]).value.undef)
                                .offset(row as isize)
                                != 0
                            {
                                fits_parser_yyerror(
                                    lParse,
                                    cs!(c"Null encountered as vector index"),
                                );
                                free((lParse.Nodes[this_node_idx]).value.data.ptr);
                                (lParse.Nodes[this_node_idx]).value.data.ptr = ptr::null_mut();
                                break;
                            } else {
                                dimVals[i as usize] =
                                    *((lParse.Nodes[theDims[i as usize]]).value.data.lngptr)
                                        .offset(row as isize);
                            }
                        }
                        i += 1;
                    }
                    if lParse.status != 0 {
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
                            || dimVals[i as usize] > (lParse.Nodes[theVar]).value.naxes[i as usize]
                        {
                            break;
                        }
                        elem = (lParse.Nodes[theVar]).value.naxes[i as usize] * elem
                            + dimVals[i as usize]
                            - 1;
                    }
                    if i < 0 {
                        elem += row * (lParse.Nodes[theVar]).value.nelem;
                        if (lParse.Nodes[this_node_idx]).ntype
                            == fits_parser_yytokentype::STRING as c_int
                        {
                            *((lParse.Nodes[this_node_idx]).value.undef).offset(row as isize) =
                                *((lParse.Nodes[theVar]).value.undef).offset(row as isize);
                        } else if (lParse.Nodes[this_node_idx]).ntype
                            != fits_parser_yytokentype::BITSTR as c_int
                        {
                            /* Dummy - BITSTRs do not have undefs */
                            *((lParse.Nodes[this_node_idx]).value.undef).offset(row as isize) =
                                *((lParse.Nodes[theVar]).value.undef).offset(elem as isize);
                        }
                        if (lParse.Nodes[this_node_idx]).ntype
                            == fits_parser_yytokentype::DOUBLE as c_int
                        {
                            *((lParse.Nodes[this_node_idx]).value.data.dblptr)
                                .offset(row as isize) =
                                *((lParse.Nodes[theVar]).value.data.dblptr).offset(elem as isize);
                        } else if (lParse.Nodes[this_node_idx]).ntype
                            == fits_parser_yytokentype::LONG as c_int
                        {
                            *((lParse.Nodes[this_node_idx]).value.data.lngptr)
                                .offset(row as isize) =
                                *((lParse.Nodes[theVar]).value.data.lngptr).offset(elem as isize);
                        } else if (lParse.Nodes[this_node_idx]).ntype
                            == fits_parser_yytokentype::BOOLEAN as c_int
                        {
                            *((lParse.Nodes[this_node_idx]).value.data.logptr)
                                .offset(row as isize) =
                                *((lParse.Nodes[theVar]).value.data.logptr).offset(elem as isize);
                        } else {
                            /* XXX Note, the below expression uses knowledge of
                            the layout of the string format, namely (nelem+1)
                            characters per string, followed by (nelem+1)
                            "undef" values. */
                            *(*((lParse.Nodes[this_node_idx]).value.data.strptr)
                                .offset(row as isize))
                            .offset(0) = *(*((lParse.Nodes[theVar]).value.data.strptr).offset(0))
                                .offset((elem + row) as isize);
                            *(*((lParse.Nodes[this_node_idx]).value.data.strptr)
                                .offset(row as isize))
                            .offset(1) = 0; /* Null terminate */
                        }
                    } else {
                        fits_parser_yyerror(lParse, cs!(c"Index out of range"));
                        free((lParse.Nodes[this_node_idx]).value.data.ptr);
                        (lParse.Nodes[this_node_idx]).value.data.ptr = ptr::null_mut();
                    }
                    row += 1;
                }
            } else {
                /* Reduce dimensions by 1, using a nonconstant expression */
                row = 0;
                while row < lParse.nRows {
                    /* Index cannot be a constant */
                    if *((lParse.Nodes[theDims[0]]).value.undef).offset(row as isize) != 0 {
                        fits_parser_yyerror(lParse, cs!(c"Null encountered as vector index"));
                        free((lParse.Nodes[this_node_idx]).value.data.ptr);
                        (lParse.Nodes[this_node_idx]).value.data.ptr = ptr::null_mut();
                        break;
                    } else {
                        dimVals[0] =
                            *((lParse.Nodes[theDims[0]]).value.data.lngptr).offset(row as isize);
                        if dimVals[0] < 1
                            || dimVals[0]
                                > (lParse.Nodes[theVar]).value.naxes
                                    [((lParse.Nodes[theVar]).value.naxis - 1) as usize]
                        {
                            fits_parser_yyerror(lParse, cs!(c"Index out of range"));
                            free((lParse.Nodes[this_node_idx]).value.data.ptr);
                            (lParse.Nodes[this_node_idx]).value.data.ptr = ptr::null_mut();
                        } else if (lParse.Nodes[this_node_idx]).ntype
                            == fits_parser_yytokentype::BITSTR as c_int
                            || (lParse.Nodes[this_node_idx]).ntype
                                == fits_parser_yytokentype::STRING as c_int
                        {
                            elem = (lParse.Nodes[this_node_idx]).value.nelem * (dimVals[0] - 1);
                            elem += row * ((lParse.Nodes[theVar]).value.nelem + 1);
                            if !((lParse.Nodes[this_node_idx]).value.undef).is_null() {
                                *((lParse.Nodes[this_node_idx]).value.undef).offset(row as isize) =
                                    *((lParse.Nodes[theVar]).value.undef).offset(row as isize);
                            }
                            memcpy(
                                (*((lParse.Nodes[this_node_idx]).value.data.strptr).offset(0))
                                    .offset(
                                        (row as c_ulong)
                                            .wrapping_mul(
                                                ::core::mem::size_of::<c_char>() as c_ulong
                                            )
                                            .wrapping_mul(
                                                ((lParse.Nodes[this_node_idx]).value.nelem + 1)
                                                    as c_ulong,
                                            ) as isize,
                                    )
                                    .cast::<c_void>(),
                                (*((lParse.Nodes[theVar]).value.data.strptr).offset(0)).offset(
                                    (elem as c_ulong)
                                        .wrapping_mul(::core::mem::size_of::<c_char>() as c_ulong)
                                        as isize,
                                ) as *const c_void,
                                ((lParse.Nodes[this_node_idx]).value.nelem as c_ulong)
                                    .wrapping_mul(::core::mem::size_of::<c_char>() as c_ulong)
                                    .try_into()
                                    .unwrap(),
                            );
                            *(*((lParse.Nodes[this_node_idx]).value.data.strptr)
                                .offset(row as isize))
                            .offset((lParse.Nodes[this_node_idx]).value.nelem as isize) = 0; /* Null terminate */
                        } else {
                            elem = (lParse.Nodes[this_node_idx]).value.nelem * (dimVals[0] - 1);
                            elem += row * (lParse.Nodes[theVar]).value.nelem;
                            memcpy(
                                ((lParse.Nodes[this_node_idx]).value.undef)
                                    .offset(
                                        (row * (lParse.Nodes[this_node_idx]).value.nelem) as isize,
                                    )
                                    .cast::<c_void>(),
                                ((lParse.Nodes[theVar]).value.undef).offset(elem as isize)
                                    as *const c_void,
                                ((lParse.Nodes[this_node_idx]).value.nelem as c_ulong)
                                    .wrapping_mul(::core::mem::size_of::<c_char>() as c_ulong)
                                    .try_into()
                                    .unwrap(),
                            );
                            memcpy(
                                (lParse.Nodes[this_node_idx])
                                    .value
                                    .data
                                    .ptr
                                    .cast::<c_char>()
                                    .offset(
                                        (row * dsize * (lParse.Nodes[this_node_idx]).value.nelem)
                                            as isize,
                                    )
                                    .cast::<c_void>(),
                                (lParse.Nodes[theVar])
                                    .value
                                    .data
                                    .ptr
                                    .cast::<c_char>()
                                    .offset((elem * dsize) as isize)
                                    as *const c_void,
                                (((lParse.Nodes[this_node_idx]).value.nelem * dsize) as c_ulong)
                                    .try_into()
                                    .unwrap(),
                            );
                        }
                        row += 1;
                    }
                }
            }
        }
        if (lParse.Nodes[theVar]).operation > 0 {
            if (lParse.Nodes[theVar]).ntype == fits_parser_yytokentype::STRING as c_int
                || (lParse.Nodes[theVar]).ntype == fits_parser_yytokentype::BITSTR as c_int
            {
                free((*((lParse.Nodes[theVar]).value.data.strptr).offset(0)).cast::<c_void>());
            } else {
                free((lParse.Nodes[theVar]).value.data.ptr);
            }
        }
        i = 0;
        while i < nDims {
            if (lParse.Nodes[theDims[i as usize]]).operation > 0 {
                free((lParse.Nodes[theDims[i as usize]]).value.data.ptr);
            }
            i += 1;
        }
    }
}

fn Do_GTI(lParse: &mut ParseData, this_node_idx: usize) {
    unsafe {
        let mut theExpr: &mut Node;
        let mut theTimes: &mut Node;
        let mut start: *mut c_double = core::ptr::null_mut::<c_double>();
        let mut stop: *mut c_double = core::ptr::null_mut::<c_double>();
        let mut times: *mut c_double = core::ptr::null_mut::<c_double>();
        let mut elem: c_long = 0;
        let mut nGTI: c_long = 0;
        let mut gti: c_long = 0;
        let mut ordered: c_int = 0;
        let dorow: c_int =
            c_int::from((lParse.Nodes[this_node_idx]).operation == GTIFIND_FCT as c_int);

        let theTimes = (lParse.Nodes[this_node_idx]).SubNodes[0];
        let theExpr = (lParse.Nodes[this_node_idx]).SubNodes[1];

        nGTI = (lParse.Nodes[theTimes]).value.nelem;
        start = (lParse.Nodes[theTimes]).value.data.dblptr;
        stop = ((lParse.Nodes[theTimes]).value.data.dblptr).offset(nGTI as isize);
        ordered = (lParse.Nodes[theTimes]).ntype;
        if (lParse.Nodes[theExpr]).operation == CONST_OP {
            gti = Search_GTI(
                (lParse.Nodes[theExpr]).value.data.dbl,
                nGTI,
                start,
                stop,
                ordered,
                core::ptr::null_mut::<c_long>(),
            );
            if dorow != 0 {
                (lParse.Nodes[this_node_idx]).value.data.lng =
                    if gti >= 0 { gti + 1 } else { -(1) as c_long };
            } else {
                (lParse.Nodes[this_node_idx]).value.data.log = if gti >= 0 { 1 } else { 0 };
            }
            (lParse.Nodes[this_node_idx]).operation = CONST_OP;
        } else {
            Allocate_Ptrs(lParse, this_node_idx);
            times = (lParse.Nodes[theExpr]).value.data.dblptr;
            if lParse.status == 0 {
                elem = lParse.nRows * (lParse.Nodes[this_node_idx]).value.nelem;
                if nGTI != 0 {
                    gti = -1;
                    loop {
                        let fresh202 = elem;
                        elem -= 1;
                        if fresh202 == 0 {
                            break;
                        }
                        let fresh203 =
                            &mut *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize);
                        *fresh203 = *((lParse.Nodes[theExpr]).value.undef).offset(elem as isize);
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
                                core::ptr::null_mut::<c_long>(),
                            );
                        }
                        if dorow != 0 {
                            *((lParse.Nodes[this_node_idx]).value.data.lngptr)
                                .offset(elem as isize) =
                                if gti >= 0 { gti + 1 } else { -(1) as c_long };
                            *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize) =
                                (if gti >= 0 { 0 } else { 1 }) as c_char;
                        } else {
                            *((lParse.Nodes[this_node_idx]).value.data.logptr)
                                .offset(elem as isize) = if gti >= 0 { 1 } else { 0 };
                        }
                    }
                } else if dorow != 0 {
                    loop {
                        let fresh204 = elem;
                        elem -= 1;
                        if fresh204 == 0 {
                            break;
                        }
                        *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize) = 1;
                    }
                } else {
                    loop {
                        let fresh205 = elem;
                        elem -= 1;
                        if fresh205 == 0 {
                            break;
                        }
                        *((lParse.Nodes[this_node_idx]).value.data.logptr).offset(elem as isize) =
                            0;
                        *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize) = 0;
                    }
                }
            }
        }
        if (lParse.Nodes[theExpr]).operation > 0 {
            free((lParse.Nodes[theExpr]).value.data.ptr);
        }
    }
}

fn Do_GTI_Over(lParse: &mut ParseData, this_node_idx: usize) {
    unsafe {
        let mut theTimes: &mut Node;
        let mut theStart: &mut Node;
        let mut theStop: &mut Node;
        let mut gtiStart: *mut c_double = core::ptr::null_mut::<c_double>();
        let mut gtiStop: *mut c_double = core::ptr::null_mut::<c_double>();
        let mut evtStart: *mut c_double = core::ptr::null_mut::<c_double>();
        let mut evtStop: *mut c_double = core::ptr::null_mut::<c_double>();
        let mut elem: c_long = 0;
        let mut nGTI: c_long = 0;
        let mut gti: c_long = 0;
        let nextGTI: c_long = 0;
        let ordered: c_int = 0;

        let theTimes = (lParse.Nodes[this_node_idx]).SubNodes[0];
        let theStop = (lParse.Nodes[this_node_idx]).SubNodes[2];
        let theStart = (lParse.Nodes[this_node_idx]).SubNodes[1];

        nGTI = (lParse.Nodes[theTimes]).value.nelem;
        gtiStart = (lParse.Nodes[theTimes]).value.data.dblptr;
        gtiStop = ((lParse.Nodes[theTimes]).value.data.dblptr).offset(nGTI as isize);
        if (lParse.Nodes[theStart]).operation == CONST_OP
            && (lParse.Nodes[theStop]).operation == CONST_OP
        {
            (lParse.Nodes[this_node_idx]).value.data.dbl = GTI_Over(
                (lParse.Nodes[theStart]).value.data.dbl,
                (lParse.Nodes[theStop]).value.data.dbl,
                nGTI,
                gtiStart,
                gtiStop,
                &mut gti,
            );
            (lParse.Nodes[this_node_idx]).operation = CONST_OP;
        } else {
            let mut undefStart: c_char = 0;
            let mut undefStop: c_char = 0;
            let mut uStart: c_double = 0.0;
            let mut uStop: c_double = 0.0;
            if (lParse.Nodes[theStart]).operation == CONST_OP {
                uStart = (lParse.Nodes[theStart]).value.data.dbl;
            }
            if (lParse.Nodes[theStop]).operation == CONST_OP {
                uStop = (lParse.Nodes[theStop]).value.data.dbl;
            }

            Allocate_Ptrs(lParse, this_node_idx);

            evtStart = (lParse.Nodes[theStart]).value.data.dblptr;
            evtStop = (lParse.Nodes[theStop]).value.data.dblptr;
            if lParse.status == 0 {
                elem = lParse.nRows * (lParse.Nodes[this_node_idx]).value.nelem;
                if nGTI != 0 {
                    let mut toverlap: c_double = 0.0;
                    gti = -1;
                    loop {
                        let fresh206 = elem;
                        elem -= 1;
                        if fresh206 == 0 {
                            break;
                        }
                        if (lParse.Nodes[theStart]).operation != CONST_OP {
                            undefStart =
                                *((lParse.Nodes[theStart]).value.undef).offset(elem as isize);
                            uStart = *evtStart.offset(elem as isize);
                        }
                        if (lParse.Nodes[theStop]).operation != CONST_OP {
                            undefStop =
                                *((lParse.Nodes[theStop]).value.undef).offset(elem as isize);
                            uStop = *evtStop.offset(elem as isize);
                        }
                        let fresh207 =
                            &mut *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize);
                        *fresh207 = if c_int::from(undefStart) != 0 || c_int::from(undefStop) != 0 {
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
                        *((lParse.Nodes[this_node_idx]).value.data.dblptr).offset(elem as isize) =
                            toverlap;
                    }
                } else {
                    loop {
                        let fresh208 = elem;
                        elem -= 1;
                        if fresh208 == 0 {
                            break;
                        }
                        *((lParse.Nodes[this_node_idx]).value.data.dblptr).offset(elem as isize) =
                            0.0;
                        *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize) = 0;
                    }
                }
            }
        }
        if (lParse.Nodes[theStart]).operation > 0 {
            free((lParse.Nodes[theStart]).value.data.ptr);
        }
        if (lParse.Nodes[theStop]).operation > 0 {
            free((lParse.Nodes[theStop]).value.data.ptr);
        }
    }
}
fn GTI_Over(
    evtStart: c_double,
    evtStop: c_double,
    nGTI: c_long,
    start: *mut c_double,
    stop: *mut c_double,
    gtiout: *mut c_long,
) -> c_double {
    unsafe {
        let mut gti1: c_long = 0;
        let mut gti2: c_long = 0;
        let mut nextGTI1: c_long = 0;
        let mut nextGTI2: c_long = 0;
        let mut gti: c_long = 0;
        let mut nMax: c_long = 0;
        let mut overlap: c_double = 0.0;
        *gtiout = -(1 as c_long);
        if evtStop <= evtStart {
            return 0.0;
        }
        gti1 = Search_GTI(evtStart, nGTI, start, stop, 1, &mut nextGTI1);
        gti2 = Search_GTI(evtStop, nGTI, start, stop, 1, &mut nextGTI2);
        if gti1 >= 0 {
            *gtiout = gti1;
        }
        if nextGTI1 < 0 && nextGTI2 < 0 {
            return 0.0;
        }
        if gti1 < 0 && gti2 < 0 && nextGTI1 == nextGTI2 {
            return 0.0;
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
        }
        overlap
    }
}
fn Search_GTI(
    evtTime: c_double,
    nGTI: c_long,
    start: *mut c_double,
    stop: *mut c_double,
    ordered: c_int,
    nextGTI0: *mut c_long,
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
            nextGTI = -1;
        }
        if !nextGTI0.is_null() {
            *nextGTI0 = nextGTI;
        }
        gti
    }
}
fn Do_REG(lParse: &mut ParseData, this_node_idx: usize) {
    unsafe {
        let mut theRegion: &mut Node;
        let mut theX: &mut Node;
        let mut theY: &mut Node;
        let mut Xval: c_double = 0.0;
        let mut Yval: c_double = 0.0;
        let mut Xnull: c_char = 0;
        let mut Ynull: c_char = 0;
        let mut Xvector: c_int = 0;
        let mut Yvector: c_int = 0;
        let mut nelem: c_long = 0;
        let mut elem: c_long = 0;
        let mut rows: c_long = 0;

        let theRegion = (lParse.Nodes[this_node_idx]).SubNodes[0];
        let theX = (lParse.Nodes[this_node_idx]).SubNodes[1];
        let theY = (lParse.Nodes[this_node_idx]).SubNodes[2];

        Xvector = c_int::from((lParse.Nodes[theX]).operation != CONST_OP);
        if Xvector != 0 {
            Xvector = (lParse.Nodes[theX]).value.nelem as c_int;
        } else {
            Xval = (lParse.Nodes[theX]).value.data.dbl;
        }
        Yvector = c_int::from((lParse.Nodes[theY]).operation != CONST_OP);
        if Yvector != 0 {
            Yvector = (lParse.Nodes[theY]).value.nelem as c_int;
        } else {
            Yval = (lParse.Nodes[theY]).value.data.dbl;
        }
        if Xvector == 0 && Yvector == 0 {
            (lParse.Nodes[this_node_idx]).value.data.log = if fits_in_region(
                Xval,
                Yval,
                &mut *(lParse.Nodes[theRegion]).value.data.ptr.cast::<SAORegion>(),
            ) != 0
            {
                1
            } else {
                0
            };
            (lParse.Nodes[this_node_idx]).operation = CONST_OP;
        } else {
            Allocate_Ptrs(lParse, this_node_idx);
            if lParse.status == 0 {
                rows = lParse.nRows;
                nelem = (lParse.Nodes[this_node_idx]).value.nelem;
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
                        if Xvector > 1 {
                            Xval = *((lParse.Nodes[theX]).value.data.dblptr).offset(elem as isize);
                            Xnull = *((lParse.Nodes[theX]).value.undef).offset(elem as isize);
                        } else if Xvector != 0 {
                            Xval = *((lParse.Nodes[theX]).value.data.dblptr).offset(rows as isize);
                            Xnull = *((lParse.Nodes[theX]).value.undef).offset(rows as isize);
                        }
                        if Yvector > 1 {
                            Yval = *((lParse.Nodes[theY]).value.data.dblptr).offset(elem as isize);
                            Ynull = *((lParse.Nodes[theY]).value.undef).offset(elem as isize);
                        } else if Yvector != 0 {
                            Yval = *((lParse.Nodes[theY]).value.data.dblptr).offset(rows as isize);
                            Ynull = *((lParse.Nodes[theY]).value.undef).offset(rows as isize);
                        }
                        *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize) =
                            if c_int::from(Xnull) != 0 || c_int::from(Ynull) != 0 {
                                1
                            } else {
                                0
                            };
                        if *((lParse.Nodes[this_node_idx]).value.undef).offset(elem as isize) != 0 {
                            continue;
                        }
                        *((lParse.Nodes[this_node_idx]).value.data.logptr).offset(elem as isize) =
                            if fits_in_region(
                                Xval,
                                Yval,
                                &mut *(lParse.Nodes[theRegion]).value.data.ptr.cast::<SAORegion>(),
                            ) != 0
                            {
                                1
                            } else {
                                0
                            };
                    }
                    nelem = (lParse.Nodes[this_node_idx]).value.nelem;
                }
            }
        }
        if (lParse.Nodes[theX]).operation > 0 {
            free((lParse.Nodes[theX]).value.data.ptr);
        }
        if (lParse.Nodes[theY]).operation > 0 {
            free((lParse.Nodes[theY]).value.data.ptr);
        }
    }
}

fn get_this_that_nodes<T>(nodes: &mut [T], this_idx: usize, that_idx: usize) -> (&mut T, &mut T) {
    // Example: nodes length = 5
    // this_idx = 4 (i.e the last item)
    // that_idx = 2
    // left will contain items 0, 1, 2, 3 and right will contain items 4

    let split_idx = cmp::max(this_idx, that_idx);
    let (left, right) = nodes.split_at_mut(split_idx);
    if this_idx < that_idx {
        let that_node = &mut right[0];
        let this_node = &mut left[this_idx];
        (this_node, that_node)
    } else {
        let this_node = &mut right[0];
        let that_node = &mut left[that_idx];
        (this_node, that_node)
    }
}

/*
fn get_this_that1_that2_nodes<T>(
    nodes: &mut [T],
    this_idx: usize,
    that1_idx: usize,
    that2_idx: usize,
) -> (&mut T, &mut T, &mut T) {
    // Example: nodes length = 5
    // this_idx = 4 (i.e the last item)
    // that1_idx = 2
    // that2_idx = 1
    // `right` will contain items 4
    // `mid` will contain items 2, 3
    // `left` will contain items 0, 1

    assert_ne!(this_idx, that1_idx);
    assert_ne!(this_idx, that2_idx);
    assert_ne!(that1_idx, that2_idx);

    let mut indices = [this_idx, that1_idx, that2_idx];
    indices.sort();

    // Right will now contain one of the nodes at the 0 index
    let (tmp, right) = nodes.split_at_mut(indices[2]);

    // Mid will now contain one of the nodes at the 0 index
    let (left, mid) = tmp.split_at_mut(indices[1]);

    let a = &mut left[indices[0]];
    let b = &mut mid[0];
    let c = &mut right[0];

    if this_idx < that1_idx && that1_idx < that2_idx {
        return (a, b, c);
    }

    if this_idx < that1_idx && that1_idx > that2_idx {
        return (a, c, b);
    }

    if that1_idx < this_idx && this_idx < that2_idx {
        return (b, a, c);
    }

    if that1_idx < this_idx && this_idx > that2_idx {
        return (c, a, b);
    }

    if that2_idx < this_idx && this_idx < that1_idx {
        return (b, c, a);
    }

    if that2_idx < this_idx && this_idx > that1_idx {
        return (c, b, a);
    }

    unreachable!()
}
*/

fn Do_Vector(lParse: &mut ParseData, this_node_idx: usize) {
    unsafe {
        let mut that: &mut Node;
        let mut row: c_long = 0;
        let mut elem: c_long = 0;
        let mut idx: c_long = 0;
        let mut jdx: c_long = 0;
        let mut offset: c_long = 0;
        let mut node: c_int = 0;

        Allocate_Ptrs(lParse, this_node_idx);

        if lParse.status == 0 {
            node = 0;
            while node < (lParse.Nodes[this_node_idx]).nSubNodes {
                let that_node_idx = (lParse.Nodes[this_node_idx]).SubNodes[node as usize];
                let (this_node, that_node) =
                    get_this_that_nodes(&mut lParse.Nodes, this_node_idx, that_node_idx);

                if that_node.operation == CONST_OP {
                    idx = lParse.nRows * (this_node).value.nelem + offset;
                    loop {
                        idx -= (this_node).value.nelem;
                        if idx < 0 {
                            break;
                        }
                        *((this_node).value.undef).offset(idx as isize) = 0;
                        match (this_node).ntype.into() {
                            fits_parser_yytokentype::BOOLEAN => {
                                *((this_node).value.data.logptr).offset(idx as isize) =
                                    that_node.value.data.log;
                            }
                            fits_parser_yytokentype::LONG => {
                                *((this_node).value.data.lngptr).offset(idx as isize) =
                                    that_node.value.data.lng;
                            }
                            fits_parser_yytokentype::DOUBLE => {
                                *((this_node).value.data.dblptr).offset(idx as isize) =
                                    that_node.value.data.dbl;
                            }
                            _ => {}
                        }
                    }
                } else {
                    row = lParse.nRows;
                    idx = row * that_node.value.nelem;
                    loop {
                        let fresh212 = row;
                        row -= 1;
                        if fresh212 == 0 {
                            break;
                        }
                        elem = that_node.value.nelem;
                        jdx = row * (this_node).value.nelem + offset;
                        loop {
                            let fresh213 = elem;
                            elem -= 1;
                            if fresh213 == 0 {
                                break;
                            }
                            idx -= 1;
                            *((this_node).value.undef).offset((jdx + elem) as isize) =
                                *(that_node.value.undef).offset(idx as isize);
                            match (this_node).ntype.into() {
                                fits_parser_yytokentype::BOOLEAN => {
                                    *((this_node).value.data.logptr)
                                        .offset((jdx + elem) as isize) =
                                        *(that_node.value.data.logptr).offset(idx as isize);
                                }
                                fits_parser_yytokentype::LONG => {
                                    *((this_node).value.data.lngptr)
                                        .offset((jdx + elem) as isize) =
                                        *(that_node.value.data.lngptr).offset(idx as isize);
                                }
                                fits_parser_yytokentype::DOUBLE => {
                                    *((this_node).value.data.dblptr)
                                        .offset((jdx + elem) as isize) =
                                        *(that_node.value.data.dblptr).offset(idx as isize);
                                }
                                _ => {}
                            }
                        }
                    }
                }
                offset += that_node.value.nelem;
                node += 1;
            }
        }

        node = 0;
        while node < (lParse.Nodes[this_node_idx]).nSubNodes {
            if ((lParse.Nodes)[(lParse.Nodes[this_node_idx]).SubNodes[node as usize] as usize])
                .operation
                > 0
            {
                free(
                    ((lParse.Nodes)
                        [(lParse.Nodes[this_node_idx]).SubNodes[node as usize] as usize])
                        .value
                        .data
                        .ptr,
                );
            }
            node += 1;
        }
    }
}

fn Do_Array(lParse: &mut ParseData, this_node_idx: usize) {
    unsafe {
        let mut that: &mut Node;
        let mut row: c_long = 0;
        let mut elem: c_long = 0;
        let mut idx: c_long = 0;
        let jdx: c_long = 0;
        let offset: c_long = 0;
        let node: c_int = 0;

        Allocate_Ptrs(lParse, this_node_idx);

        let that_idx = (lParse.Nodes[this_node_idx]).SubNodes[0];
        let (this_node, that_node) =
            get_this_that_nodes(&mut lParse.Nodes, this_node_idx, that_idx);

        if lParse.status == 0 {
            if (that_node).operation == CONST_OP {
                idx = lParse.nRows * this_node.value.nelem + offset;
                loop {
                    let fresh214 = idx;
                    idx -= 1;
                    if fresh214 == 0 {
                        break;
                    }
                    *(this_node.value.undef).offset(idx as isize) = 0;
                    match this_node.ntype.into() {
                        fits_parser_yytokentype::BOOLEAN => {
                            *(this_node.value.data.logptr).offset(idx as isize) =
                                (that_node).value.data.log;
                        }
                        fits_parser_yytokentype::LONG => {
                            *(this_node.value.data.lngptr).offset(idx as isize) =
                                (that_node).value.data.lng;
                        }
                        fits_parser_yytokentype::DOUBLE => {
                            *(this_node.value.data.dblptr).offset(idx as isize) =
                                (that_node).value.data.dbl;
                        }
                        _ => {}
                    }
                }
            } else if (that_node).value.nelem > 1 {
                idx = lParse.nRows * this_node.value.nelem;
                loop {
                    let fresh215 = idx;
                    idx -= 1;
                    if fresh215 == 0 {
                        break;
                    }
                    *(this_node.value.undef).offset(idx as isize) =
                        *((that_node).value.undef).offset(idx as isize);
                    match this_node.ntype.into() {
                        fits_parser_yytokentype::BOOLEAN => {
                            *(this_node.value.data.logptr).offset(idx as isize) =
                                *((that_node).value.data.logptr).offset(idx as isize);
                        }
                        fits_parser_yytokentype::LONG => {
                            *(this_node.value.data.lngptr).offset(idx as isize) =
                                *((that_node).value.data.lngptr).offset(idx as isize);
                        }
                        fits_parser_yytokentype::DOUBLE => {
                            *(this_node.value.data.dblptr).offset(idx as isize) =
                                *((that_node).value.data.dblptr).offset(idx as isize);
                        }
                        _ => {}
                    }
                }
            } else {
                row = lParse.nRows;
                idx = row * this_node.value.nelem - 1;
                loop {
                    let fresh216 = row;
                    row -= 1;
                    if fresh216 == 0 {
                        break;
                    }
                    elem = this_node.value.nelem;
                    loop {
                        let fresh217 = elem;
                        elem -= 1;
                        if fresh217 == 0 {
                            break;
                        }
                        *(this_node.value.undef).offset(idx as isize) =
                            *((that_node).value.undef).offset(row as isize);
                        match this_node.ntype.into() {
                            fits_parser_yytokentype::BOOLEAN => {
                                *(this_node.value.data.logptr).offset(idx as isize) =
                                    *((that_node).value.data.logptr).offset(row as isize);
                            }
                            fits_parser_yytokentype::LONG => {
                                *(this_node.value.data.lngptr).offset(idx as isize) =
                                    *((that_node).value.data.lngptr).offset(row as isize);
                            }
                            fits_parser_yytokentype::DOUBLE => {
                                *(this_node.value.data.dblptr).offset(idx as isize) =
                                    *((that_node).value.data.dblptr).offset(row as isize);
                            }
                            _ => {}
                        }
                        idx -= 1;
                    }
                }
            }

            if ((lParse.Nodes)[lParse.Nodes[this_node_idx].SubNodes[0]]).operation > 0 {
                free(
                    ((lParse.Nodes)[lParse.Nodes[this_node_idx].SubNodes[0]])
                        .value
                        .data
                        .ptr,
                );
            }
        }
    }
}

fn bitlgte(mut bits1: *mut c_char, oper: c_int, mut bits2: *mut c_char) -> c_char {
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
        )
        .cast::<c_char>();
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
            if c_int::from(chr1) != 'x' as i32
                && c_int::from(chr1) != 'X' as i32
                && c_int::from(chr2) != 'x' as i32
                && c_int::from(chr2) != 'X' as i32
            {
                if c_int::from(chr1) == '1' as i32 {
                    val1 += nextbit;
                }
                if c_int::from(chr2) == '1' as i32 {
                    val2 += nextbit;
                }
                nextbit *= 2;
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
        free(stream.cast::<c_void>());
        result
    }
}
fn bitand(mut result: *mut c_char, mut bitstrm1: *mut c_char, mut bitstrm2: *mut c_char) {
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
        )
        .cast::<c_char>();
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
            if c_int::from(chr1) == 'x' as i32 || c_int::from(chr2) == 'x' as i32 {
                *result = b'x' as c_char;
            } else if c_int::from(chr1) == '1' as i32 && c_int::from(chr2) == '1' as i32 {
                *result = b'1' as c_char;
            } else {
                *result = b'0' as c_char;
            }
            result = result.offset(1);
        }
        free(stream.cast::<c_void>());
        *result = 0;
    }
}
fn bitor(mut result: *mut c_char, mut bitstrm1: *mut c_char, mut bitstrm2: *mut c_char) {
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
        )
        .cast::<c_char>();
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
            if c_int::from(chr1) == '1' as i32 || c_int::from(chr2) == '1' as i32 {
                *result = b'1' as c_char;
            } else if c_int::from(chr1) == '0' as i32 || c_int::from(chr2) == '0' as i32 {
                *result = b'0' as c_char;
            } else {
                *result = b'x' as c_char;
            }
            result = result.offset(1);
        }
        free(stream.cast::<c_void>());
        *result = 0;
    }
}
fn bitnot(mut result: *mut c_char, mut bits: *mut c_char) {
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
            *fresh255 = (if c_int::from(chr) == '1' as i32 {
                '0' as i32
            } else if c_int::from(chr) == '0' as i32 {
                '1' as i32
            } else {
                c_int::from(chr)
            }) as c_char;
        }
        *result = 0;
    }
}

fn bitcmp(bitstrm1: *mut c_char, bitstrm2: *mut c_char) -> c_char {
    unsafe {
        let mut i: c_int = 0;
        let mut ldiff: c_int = 0;
        let mut largestStream: c_int = 0;
        let mut chr1: c_char = 0;
        let mut chr2: c_char = 0;

        let mut l1 = strlen(bitstrm1) as c_int;
        let mut l2 = strlen(bitstrm2) as c_int;

        let mut bitstrm1 = core::slice::from_raw_parts_mut(bitstrm1, l1 as usize + 1);
        let mut bitstrm2 = core::slice::from_raw_parts_mut(bitstrm2, l2 as usize + 1);

        largestStream = cmp::max(l1, l2);

        let mut stream_vec: Vec<c_char> = vec![0; (largestStream + 1) as usize];
        let stream = &mut stream_vec[..];

        if l1 < l2 {
            ldiff = l2 - l1;
            i = 0;

            while ldiff > 0 {
                stream[i as usize] = b'0' as c_char;
                i += 1;
                ldiff -= 1;
            }

            while l1 > 0 {
                stream[i as usize] = bitstrm1[0];
                bitstrm1 = &mut bitstrm1[1..];
                i += 1;
                l1 -= 1;
            }

            stream[i as usize] = 0;
            bitstrm1 = stream;
        } else if l2 < l1 {
            ldiff = l1 - l2;
            i = 0;

            while ldiff > 0 {
                stream[i as usize] = b'0' as c_char;
                i += 1;
                ldiff -= 1;
            }

            while l2 > 0 {
                stream[i as usize] = bitstrm2[0];
                bitstrm2 = &mut bitstrm2[1..];
                i += 1;
                l2 -= 1;
            }

            stream[i as usize] = 0;
            bitstrm2 = stream;
        }

        loop {
            chr1 = bitstrm1[0];
            if chr1 == 0 {
                break;
            }
            bitstrm1 = &mut bitstrm1[1..];

            chr2 = bitstrm2[0];
            bitstrm2 = &mut bitstrm2[1..];

            if (chr1 == b'0' as c_char && chr2 == b'1' as c_char)
                || (chr1 == b'1' as c_char && chr2 == b'0' as c_char)
            {
                return 0;
            }
        }
        1
    }
}

fn bnear(x: c_double, y: c_double, tolerance: c_double) -> c_char {
    if (x - y).abs() < tolerance { 1 } else { 0 }
}

fn saobox(
    xcen: c_double,
    ycen: c_double,
    xwid: c_double,
    ywid: c_double,
    rot: c_double,
    xcol: c_double,
    ycol: c_double,
) -> c_char {
    let mut x: c_double = 0.0;
    let mut y: c_double = 0.0;
    let mut xprime: c_double = 0.0;
    let mut yprime: c_double = 0.0;
    let mut xmin: c_double = 0.0;
    let mut xmax: c_double = 0.0;
    let mut ymin: c_double = 0.0;
    let mut ymax: c_double = 0.0;
    let mut theta: c_double = 0.0;
    theta = rot / 180.0 * MY_PI;
    xprime = xcol - xcen;
    yprime = ycol - ycen;
    x = xprime * (theta.cos()) + yprime * (theta.sin());
    y = -xprime * (theta.sin()) + yprime * (theta.cos());
    xmin = -0.5 * xwid;
    xmax = 0.5 * xwid;
    ymin = -0.5 * ywid;
    ymax = 0.5 * ywid;

    if x >= xmin && x <= xmax && y >= ymin && y <= ymax {
        1
    } else {
        0
    }
}

fn circle(xcen: c_double, ycen: c_double, rad: c_double, xcol: c_double, ycol: c_double) -> c_char {
    let mut r2: c_double = 0.0;
    let mut dx: c_double = 0.0;
    let mut dy: c_double = 0.0;
    let mut dlen: c_double = 0.0;
    dx = xcol - xcen;
    dy = ycol - ycen;
    dx *= dx;
    dy *= dy;
    dlen = dx + dy;
    r2 = rad * rad;
    if dlen <= r2 { 1 } else { 0 }
}

fn ellipse(
    xcen: c_double,
    ycen: c_double,
    xrad: c_double,
    yrad: c_double,
    rot: c_double,
    xcol: c_double,
    ycol: c_double,
) -> c_char {
    let mut x: c_double = 0.0;
    let mut y: c_double = 0.0;
    let mut xprime: c_double = 0.0;
    let mut yprime: c_double = 0.0;
    let mut dx: c_double = 0.0;
    let mut dy: c_double = 0.0;
    let mut dlen: c_double = 0.0;
    let mut theta: c_double = 0.0;

    theta = rot / 180.0 * MY_PI;
    xprime = xcol - xcen;
    yprime = ycol - ycen;
    x = xprime * (theta.cos()) + yprime * (theta.sin());
    y = -xprime * (theta.sin()) + yprime * (theta.cos());
    dx = x / xrad;
    dy = y / yrad;
    dx *= dx;
    dy *= dy;
    dlen = dx + dy;
    if dlen <= 1.0 { 1 } else { 0 }
}

/// Extract substring
fn cstrmid(
    lParse: &mut ParseData,
    dest_str: *mut c_char,
    dest_len: c_int,
    src_str: *const c_char,
    src_len: c_int,
    pos: c_int,
) -> c_int {
    let dest_str = unsafe { core::slice::from_raw_parts_mut(dest_str, (dest_len) as usize + 1) }; // +1 for null terminate
    let dest_len = dest_len as usize;
    let mut src_len = src_len as usize;

    let mut endpos: usize = 0;

    if src_len == 0 {
        src_len = unsafe { strlen(src_str) };
    } /* .. if constant */

    let src_str = unsafe { core::slice::from_raw_parts(src_str, src_len + 1) };

    if pos < 0 {
        fits_parser_yyerror(lParse, cs!(c"STRMID(S,P,N) P must be 0 or greater"));
        return -(1);
    }

    let pos = pos as usize; // Already checked if < 0

    if pos > src_len || pos == 0 {
        /* pos==0: null string requested */
        endpos = 0;
    } else if pos - 1 + dest_len > src_len {
        /* Copy only to end of src_str */
        let nsub = src_len - pos + 1;
        endpos = nsub;
        dest_str[..nsub].copy_from_slice(&src_str[(pos - 1)..(pos - 1 + nsub)]);
    } else {
        /* Full string copy */
        endpos = dest_len;
        dest_str[..dest_len].copy_from_slice(&src_str[(pos - 1)..(pos - 1 + dest_len)]);
    }

    /* Null-terminate */
    dest_str[endpos] = 0;

    0
}

fn fits_parser_yyerror(lParse: &mut ParseData, s: &[c_char]) {
    let mut msg: [c_char; 80] = [0; 80];
    if lParse.status == 0 {
        lParse.status = PARSE_SYNTAX_ERR;
    }

    strncpy_safe(&mut msg, s, 80);
    msg[79] = 0;

    ffpmsg_slice(&msg);
}
