/************************************************************************/
/*                                                                      */
/*                       CFITSIO Lexical Parser                         */
/*                                                                      */
/* This specifies a thread-safe reentrant version of lex functions */
/* This specifies CFITSIO-unique names for lexer functions */
/* This facilitates calling between the Bison parser and this lexer */
/* This file is one of 3 files containing code which parses an          */
/* arithmetic expression and evaluates it in the context of an input    */
/* FITS file table extension.  The CFITSIO lexical parser is divided    */
/* into the following 3 parts/files: the CFITSIO "front-end",           */
/* eval_f.c, contains the interface between the user/CFITSIO and the    */
/* real core of the parser; the FLEX interpreter, eval_l.c, takes the   */
/* input string and parses it into tokens and identifies the FITS       */
/* information required to evaluate the expression (ie, keywords and    */
/* columns); and, the BISON grammar and evaluation routines, eval_y.c,  */
/* receives the FLEX output and determines and performs the actual      */
/* operations.  The files eval_l.c and eval_y.c are produced from       */
/* running flex and bison on the files eval.l and eval.y, respectively. */
/* (flex and bison are available from any GNU archive: see www.gnu.org) */
/*                                                                      */
/* The grammar rules, rather than evaluating the expression in situ,    */
/* builds a tree, or Nodal, structure mapping out the order of          */
/* operations and expression dependencies.  This "compilation" process  */
/* allows for much faster processing of multiple rows.  This technique  */
/* was developed by Uwe Lammers of the XMM Science Analysis System,     */
/* although the CFITSIO implementation is entirely code original.       */
/*                                                                      */
/*                                                                      */
/* Modification History:                                                */
/*                                                                      */
/*   Kent Blackburn      c1992  Original parser code developed for the  */
/*                              FTOOLS software package, in particular, */
/*                              the fselect task.                       */
/*   Kent Blackburn      c1995  BIT column support added                */
/*   Peter D Wilson   Feb 1998  Vector column support added             */
/*   Peter D Wilson   May 1998  Ported to CFITSIO library.  User        */
/*                              interface routines written, in essence  */
/*                              making fselect, fcalc, and maketime     */
/*                              capabilities available to all tools     */
/*                              via single function calls.              */
/*   Peter D Wilson   Jun 1998  Major rewrite of parser core, so as to  */
/*                              create a run-time evaluation tree,      */
/*                              inspired by the work of Uwe Lammers,    */
/*                              resulting in a speed increase of        */
/*                              10-100 times.                           */
/*   Peter D Wilson   Jul 1998  gtifilter(a,b,c,d) function added       */
/*   Peter D Wilson   Aug 1998  regfilter(a,b,c,d) function added       */
/*   Peter D Wilson   Jul 1999  Make parser fitsfile-independent,       */
/*                              allowing a purely vector-based usage    */
/*                                                                      */
/************************************************************************/
#![deny(dead_code)]

use errno::{Errno, errno, set_errno};

use core::ffi::CStr;
use core::ptr;
use core::slice;
use std::process::exit;

use bytemuck::cast_slice;
use libc::{ENOMEM, FILE, atof, atol, fileno, free, fwrite, isatty, size_t};

use crate::c_types::{c_char, c_int, c_long, c_short, c_uchar, c_uint, c_void};
use crate::eval_defs::{MAX_STRLEN, MAXVARNAME, P_ERROR, ParseData, ValueSort};
use crate::fitsio::PARSE_SYNTAX_ERR;
use crate::fitsio2::FSTRCMP;
use crate::helpers::boxed::box_try_new;
use crate::helpers::vec_raw_parts::vec_into_raw_parts;
use crate::wrappers::{
    isdigit_safe, strcat_safe, strcpy_safe, strlen_safe, strncat_safe, strncpy_safe,
};
use crate::{STDIN, STDOUT, cs, eval_tab::*};
use crate::{
    fitscore::{ffpmsg_slice, fits_strcasecmp, fits_strncasecmp},
    wrappers::{tolower, toupper},
};

/*****  Definitions  *****/

const OCT_0: &str = "000";
const OCT_1: &str = "001";
const OCT_2: &str = "010";
const OCT_3: &str = "011";
const OCT_4: &str = "100";
const OCT_5: &str = "101";
const OCT_6: &str = "110";
const OCT_7: &str = "111";
const OCT_X: &str = "xxx";

const HEX_0: &str = "0000";
const HEX_1: &str = "0001";
const HEX_2: &str = "0010";
const HEX_3: &str = "0011";
const HEX_4: &str = "0100";
const HEX_5: &str = "0101";
const HEX_6: &str = "0110";
const HEX_7: &str = "0111";
const HEX_8: &str = "1000";
const HEX_9: &str = "1001";
const HEX_A: &str = "1010";
const HEX_B: &str = "1011";
const HEX_C: &str = "1100";
const HEX_D: &str = "1101";
const HEX_E: &str = "1110";
const HEX_F: &str = "1111";
const HEX_X: &str = "xxxx";

pub type __uint8_t = c_uchar;
pub type __int16_t = c_short;
pub type __off_t = c_long;
pub type __off64_t = c_long;

pub type _IO_lock_t = ();

pub type int16_t = __int16_t;
pub type uint8_t = __uint8_t;
pub type flex_uint8_t = uint8_t;
pub type flex_int16_t = int16_t;

#[derive(Default)]
pub struct yy_buffer_state {
    pub yy_input_file: *mut FILE,
    pub yy_ch_buf: Option<Box<[c_char]>>,
    pub yy_buf_pos: usize,
    pub yy_buf_size: c_int,
    pub yy_n_chars: c_int,
    pub yy_is_our_buffer: bool,
    pub yy_is_interactive: c_int,
    pub yy_at_bol: c_int,
    pub yy_bs_lineno: c_int,
    pub yy_bs_column: c_int,
    pub yy_fill_buffer: c_int,
    pub yy_buffer_status: c_int,
}

/* Return values from yy_get_next_buffer(). */
const EOB_ACT_CONTINUE_SCAN: c_int = 0;
const EOB_ACT_END_OF_FILE: c_int = 1;
const EOB_ACT_LAST_MATCH: c_int = 2;

/* yy_buffer_status values. */
const YY_BUFFER_NEW: c_int = 0;
const YY_BUFFER_NORMAL: c_int = 1;
/* When an EOF's been seen but there's still some text to process
 * then we mark the buffer as YY_EOF_PENDING, to indicate that we
 * shouldn't try reading from the input source any more.  We might
 * still have a bunch of tokens to match, though, because of
 * possible backing-up.
 *
 * When we actually see the EOF, we change the status to "new"
 * (via yyrestart()), so that the user can continue scanning by
 * just pointing yyin at a new input file.
 */
const YY_BUFFER_EOF_PENDING: c_int = 2;

const YY_END_OF_BUFFER_CHAR: c_char = 0;

const YY_BUF_SIZE: usize = 16384;

/* Size of default input buffer. */
const YY_READ_BUF_SIZE: c_int = 8192;

/* The action number that stands for "end of buffer".
 * YY_STATE_EOF(state) is YY_END_OF_BUFFER + state + 1; INITIAL is start
 * condition 0, so its EOF action is YY_END_OF_BUFFER + 1.
 */
const YY_END_OF_BUFFER: c_int = 31;
const INITIAL: c_int = 0;

/* yyterminate() returns YY_NULL. */
const YY_NULL: c_int = 0;

pub type YY_BUFFER_STATE = *mut yy_buffer_state;

pub type yy_size_t = size_t;

/* Holds the entire state of the reentrant scanner. */
#[derive(Default, Copy, Clone)]
pub(crate) struct yyguts_t {
    /* User-defined. Not touched by flex. */
    /* This is a shorthand accessor to get at the "extra" data inside the
    lexer, which in our case is the lParse (ParseData) structure */
    /* How ParseData is accessed from the lexer, i.e. by yyextra */
    pub(crate) yyextra_r: *mut ParseData,

    /* The rest are the same as the globals declared in the non-reentrant scanner. */
    pub(crate) yyin_r: *mut FILE,
    pub(crate) yyout_r: *mut FILE,
    pub(crate) yy_buffer_stack_top: size_t,
    /**< index of top of stack. */
    pub(crate) yy_buffer_stack_max: size_t,
    /**< capacity of stack. */
    pub(crate) yy_buffer_stack: *mut Option<Box<yy_buffer_state>>,
    /**< Stack as an array. */
    pub(crate) yy_hold_char: c_char,
    pub(crate) yy_n_chars: c_int,
    pub(crate) yyleng_r: c_int,
    pub(crate) yy_c_buf_p: *mut c_char,
    pub(crate) yy_init: c_int,
    pub(crate) yy_start: c_int,
    pub(crate) yy_did_buffer_switch_on_eof: c_int,
    pub(crate) yy_last_accepting_state: yy_state_type,
    pub(crate) yy_last_accepting_cpos: *mut c_char,
    pub(crate) yytext_r: *mut c_char,
    pub(crate) yylval_r: *mut FITS_PARSER_YYSTYPE,
}

impl yyguts_t {
    /// The current token (flex's `yytext`) as a NUL-terminated slice.
    ///
    /// The scanner sets `yytext_r`/`yyleng_r` together and writes the
    /// terminator at `yytext_r[yyleng_r]`, so the token plus its NUL is always
    /// `yyleng_r + 1` characters. Rules use this to reach the token with the
    /// `*_safe` string helpers instead of raw pointer arithmetic.
    ///
    /// SAFETY: only valid inside a scanner action, where `yytext_r` points into
    /// the live input buffer. Returns an empty string before the first match.
    pub(crate) unsafe fn yytext(&self) -> &[c_char] {
        if self.yytext_r.is_null() {
            return &[0];
        }
        unsafe { slice::from_raw_parts(self.yytext_r, self.yyleng_r as usize + 1) }
    }
}

pub type yy_state_type = c_int;

pub type YY_CHAR = flex_uint8_t;
pub type uint = c_uint;

pub(crate) fn fits_parser_yywrap() -> c_int {
    /* MJT -- 13 June 1996
       Supplied for compatibility with
       pre-2.5.1 versions of flex which
       do not recognize %option noyywrap
    */
    1
}

/*
   expr_read is lifted from old ftools.skel.
   Now we can use any version of flex with
   no .skel file necessary! MJT - 13 June 1996

   keep a memory of how many bytes have been
   read previously, so that an unlimited-sized
   buffer can be supported. PDW - 28 Feb 1998
*/
fn expr_read(lParse: &mut ParseData, buf: &mut [c_char], nbytes: c_int) -> c_int {
    let mut n: c_int = 0;
    if lParse.is_eobuf == 0
        && let Some(ref expr_data) = lParse.expr
    {
        loop {
            let fresh0 = lParse.index;
            lParse.index += 1;
            let fresh1 = n;
            n += 1;

            // Safely access the byte array
            if (fresh0 as usize) < expr_data.len() {
                buf[fresh1 as usize] = expr_data[fresh0 as usize] as c_char;
            } else {
                break;
            }

            // Check if we've reached the end
            if !(n < nbytes
                && (lParse.index as usize) < expr_data.len()
                && expr_data[lParse.index as usize] != 0)
            {
                break;
            }
        }

        // Check if we've reached null terminator
        if (lParse.index as usize) < expr_data.len() && expr_data[lParse.index as usize] == 0 {
            lParse.is_eobuf = 1;
        }
    }
    buf[n as usize] = 0;
    n
}

/*****  Internal functions  *****/

pub(crate) fn fits_parser_yyGetVariable(
    lParse: &mut ParseData,
    varName: &[c_char],
    thelval: &mut FITS_PARSER_YYSTYPE,
) -> c_int {
    let varNum: c_int = find_variable(lParse, varName);
    let mut dtype: c_int = 0;
    let mut errMsg: [c_char; 105] = [0; 105];

    if varNum < 0 {
        if (lParse.getData).is_some() {
            dtype = (lParse.getData).expect("non-null function pointer")(lParse, varName, thelval);
        } else {
            dtype = -1;
            lParse.status = PARSE_SYNTAX_ERR;
            strcpy_safe(&mut errMsg, cs!(c"Unable to find data: "));
            strncat_safe(&mut errMsg, varName, MAXVARNAME);
            ffpmsg_slice(&errMsg);
        }
    } else {
        /*  Convert variable type into expression type  */
        /* ValueSort has no other variants, so this arm cannot be reached; it is
        kept because it is the C's error path and the field may widen again. */
        #[allow(unreachable_patterns)]
        match ((lParse.varData)[varNum as usize]).dtype {
            ValueSort::Long | ValueSort::Double => {
                dtype = fits_parser_yytokentype::COLUMN as c_int;
            }
            ValueSort::Boolean => {
                dtype = fits_parser_yytokentype::BCOLUMN as c_int;
            }
            ValueSort::String => {
                dtype = fits_parser_yytokentype::SCOLUMN as c_int;
            }
            ValueSort::Bits => {
                dtype = fits_parser_yytokentype::BITCOL as c_int;
            }
            _ => {
                dtype = P_ERROR;
                lParse.status = PARSE_SYNTAX_ERR;
                strcpy_safe(&mut errMsg, cs!(c"Bad datatype for data: "));
                strncat_safe(&mut errMsg, varName, MAXVARNAME);
                ffpmsg_slice(&errMsg);
            }
        }

        *thelval = FITS_PARSER_YYSTYPE::Long(c_long::from(varNum));
    }
    dtype
}

fn find_variable(lParse: &mut ParseData, varName: &[c_char]) -> c_int {
    if lParse.nCols != 0 {
        for i in 0..lParse.nCols {
            if fits_strncasecmp(&((lParse.varData)[i as usize]).name, varName, MAXVARNAME) == 0 {
                return i;
            }
        }
    }
    -1
}

static YY_ACCEPT: [flex_int16_t; 174] = [
    0, 0, 0, 31, 29, 1, 28, 18, 29, 29, 29, 29, 29, 29, 29, 29, 8, 8, 24, 29, 23, 13, 13, 13, 13,
    9, 13, 13, 13, 13, 13, 17, 13, 13, 13, 13, 13, 13, 13, 29, 1, 22, 0, 12, 0, 11, 0, 13, 20, 0,
    0, 0, 0, 0, 0, 0, 17, 0, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10, 8, 0, 0, 0, 0, 26,
    21, 25, 13, 13, 13, 2, 13, 13, 13, 4, 13, 13, 13, 13, 3, 13, 27, 13, 13, 13, 13, 13, 13, 13,
    13, 13, 19, 0, 11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 10, 5, 6, 7, 14, 13, 23, 24, 13, 13, 13, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 18, 0, 0,
    15, 0, 0, 0, 0, 0, 0, 0, 16, 0, 0,
];
static YY_EC: [YY_CHAR; 256] = [
    0, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    2, 4, 5, 6, 7, 1, 8, 9, 10, 11, 12, 13, 1, 13, 14, 1, 15, 16, 17, 17, 17, 17, 17, 17, 18, 18,
    1, 1, 19, 20, 21, 1, 1, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 31, 32, 31, 33, 34, 31, 35, 36,
    31, 37, 38, 31, 31, 39, 31, 31, 1, 1, 1, 40, 41, 1, 42, 43, 24, 44, 45, 46, 47, 29, 48, 31, 31,
    49, 31, 50, 51, 31, 52, 53, 31, 54, 55, 31, 31, 56, 31, 31, 1, 57, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
];
static YY_META: [YY_CHAR; 58] = [
    0, 1, 1, 2, 1, 1, 1, 3, 1, 1, 1, 1, 1, 1, 1, 4, 4, 4, 4, 1, 1, 1, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 1, 5, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 1,
];
static YY_BASE: [flex_int16_t; 182] = [
    0, 0, 0, 415, 416, 412, 416, 393, 407, 404, 403, 401, 399, 33, 395, 68, 111, 16, 386, 45, 385,
    28, 58, 363, 27, 362, 50, 153, 52, 54, 105, 362, 0, 62, 57, 89, 94, 106, 148, 344, 398, 416,
    394, 416, 391, 390, 389, 416, 416, 386, 360, 361, 359, 340, 341, 336, 416, 167, 190, 344, 337,
    79, 105, 106, 328, 327, 309, 306, 104, 105, 115, 303, 304, 196, 0, 204, 55, 157, 0, 416, 416,
    416, 65, 208, 0, 213, 216, 217, 222, 344, 227, 229, 230, 214, 249, 232, 416, 235, 246, 267,
    268, 270, 275, 281, 252, 282, 416, 346, 344, 312, 315, 311, 290, 293, 289, 302, 317, 327, 326,
    325, 324, 323, 322, 297, 318, 288, 271, 296, 293, 280, 272, 269, 261, 220, 258, 214, 282, 286,
    137, 297, 0, 416, 311, 416, 416, 316, 319, 321, 238, 231, 239, 205, 205, 227, 211, 189, 188,
    186, 179, 177, 416, 158, 151, 416, 137, 91, 112, 122, 78, 100, 97, 416, 93, 416, 362, 365, 370,
    375, 377, 379, 384, 88,
];
static YY_DEF: [flex_int16_t; 182] = [
    0, 173, 1, 173, 173, 173, 173, 173, 174, 175, 176, 173, 177, 173, 173, 173, 173, 16, 173, 173,
    173, 178, 178, 178, 178, 178, 178, 178, 178, 178, 178, 173, 179, 178, 178, 178, 178, 178, 178,
    173, 173, 173, 174, 173, 180, 175, 176, 173, 173, 177, 173, 173, 173, 173, 173, 173, 173, 173,
    173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 17, 173, 173,
    173, 181, 173, 173, 173, 178, 178, 179, 178, 178, 178, 178, 27, 178, 178, 178, 178, 178, 178,
    173, 178, 178, 178, 178, 178, 178, 178, 178, 178, 173, 180, 180, 173, 173, 173, 173, 173, 173,
    173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173,
    173, 173, 173, 173, 173, 173, 181, 173, 178, 173, 173, 178, 178, 178, 173, 173, 173, 173, 173,
    173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173,
    173, 0, 173, 173, 173, 173, 173, 173, 173, 173,
];
static YY_NXT: [flex_int16_t; 474] = [
    0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 4, 14, 4, 15, 16, 17, 17, 17, 18, 19, 20, 21, 22, 23, 23,
    24, 25, 26, 27, 23, 23, 28, 29, 30, 23, 23, 25, 23, 23, 31, 32, 33, 22, 23, 34, 25, 35, 23, 36,
    37, 38, 23, 23, 25, 23, 23, 39, 50, 173, 51, 83, 86, 52, 79, 80, 81, 173, 84, 84, 138, 138,
    173, 85, 85, 141, 87, 53, 90, 54, 92, 55, 57, 58, 58, 58, 58, 88, 93, 91, 59, 84, 140, 84, 60,
    84, 61, 85, 84, 84, 62, 63, 64, 84, 171, 118, 84, 65, 171, 98, 66, 171, 97, 67, 85, 68, 119,
    69, 70, 71, 94, 94, 94, 172, 72, 73, 74, 74, 74, 74, 84, 120, 122, 171, 99, 84, 75, 75, 170,
    101, 123, 95, 121, 100, 94, 169, 84, 84, 102, 128, 130, 103, 138, 138, 76, 75, 75, 104, 129,
    131, 132, 94, 77, 94, 94, 94, 133, 78, 89, 89, 89, 89, 139, 139, 139, 89, 89, 89, 89, 89, 89,
    57, 115, 115, 115, 115, 168, 94, 167, 84, 166, 96, 89, 160, 84, 89, 89, 89, 89, 89, 48, 105,
    96, 160, 94, 58, 58, 58, 58, 89, 57, 58, 58, 58, 58, 75, 75, 136, 141, 137, 137, 137, 137, 141,
    141, 48, 141, 141, 85, 85, 80, 81, 141, 142, 75, 75, 143, 141, 163, 141, 141, 79, 141, 144, 41,
    141, 106, 165, 164, 84, 163, 145, 85, 162, 84, 84, 141, 84, 84, 141, 80, 161, 141, 84, 94, 94,
    94, 159, 84, 85, 84, 84, 106, 84, 158, 41, 84, 141, 141, 146, 141, 81, 143, 144, 79, 141, 79,
    84, 94, 144, 84, 141, 141, 84, 143, 41, 106, 137, 137, 137, 137, 137, 137, 137, 137, 94, 147,
    81, 84, 84, 80, 84, 139, 139, 139, 157, 84, 115, 115, 115, 115, 141, 84, 84, 156, 48, 141, 75,
    75, 141, 160, 141, 106, 48, 155, 160, 41, 144, 79, 143, 81, 80, 154, 153, 152, 151, 75, 75,
    150, 149, 148, 108, 84, 108, 141, 135, 134, 84, 127, 126, 84, 125, 84, 42, 124, 42, 42, 42, 45,
    45, 45, 46, 117, 46, 46, 46, 49, 116, 49, 49, 49, 82, 82, 84, 84, 107, 114, 107, 107, 107, 113,
    112, 111, 110, 109, 43, 47, 173, 108, 43, 40, 106, 96, 84, 84, 81, 79, 56, 43, 48, 47, 44, 43,
    41, 40, 173, 3, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173,
    173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173,
    173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173,
    173, 173, 173,
];
static YY_CHK: [flex_int16_t; 474] = [
    0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 13, 17, 13, 21,
    24, 13, 19, 19, 19, 17, 24, 21, 76, 76, 17, 22, 22, 82, 26, 13, 28, 13, 29, 13, 15, 15, 15, 15,
    15, 26, 29, 28, 15, 26, 181, 28, 15, 29, 15, 22, 34, 22, 15, 15, 15, 33, 172, 61, 82, 15, 170,
    34, 15, 169, 33, 15, 22, 15, 61, 15, 15, 15, 30, 30, 30, 168, 15, 16, 16, 16, 16, 16, 35, 62,
    63, 167, 35, 36, 16, 16, 166, 36, 63, 30, 62, 35, 30, 165, 30, 37, 36, 68, 69, 37, 138, 138,
    16, 16, 16, 37, 68, 69, 70, 30, 16, 38, 38, 38, 70, 16, 27, 27, 27, 27, 77, 77, 77, 27, 27, 27,
    27, 27, 27, 57, 57, 57, 57, 57, 164, 38, 162, 38, 161, 159, 27, 158, 27, 27, 27, 27, 27, 27,
    157, 38, 156, 155, 38, 58, 58, 58, 58, 27, 73, 73, 73, 73, 73, 58, 58, 75, 83, 75, 75, 75, 75,
    85, 93, 154, 86, 87, 85, 85, 86, 87, 88, 83, 58, 58, 88, 90, 153, 91, 92, 90, 95, 91, 92, 97,
    95, 152, 151, 83, 150, 93, 85, 149, 85, 93, 98, 86, 87, 94, 98, 148, 104, 88, 94, 94, 94, 135,
    90, 85, 91, 92, 134, 95, 133, 132, 97, 99, 100, 97, 101, 99, 100, 131, 101, 102, 130, 98, 94,
    102, 94, 103, 105, 104, 129, 103, 105, 136, 136, 136, 136, 137, 137, 137, 137, 94, 104, 128,
    99, 100, 127, 101, 139, 139, 139, 126, 102, 115, 115, 115, 115, 142, 103, 105, 125, 142, 145,
    115, 115, 146, 145, 147, 124, 146, 123, 147, 122, 121, 120, 119, 118, 117, 116, 114, 113, 112,
    115, 115, 111, 110, 109, 108, 142, 107, 89, 72, 71, 145, 67, 66, 146, 65, 147, 174, 64, 174,
    174, 174, 175, 175, 175, 176, 60, 176, 176, 176, 177, 59, 177, 177, 177, 178, 178, 179, 179,
    180, 55, 180, 180, 180, 54, 53, 52, 51, 50, 49, 46, 45, 44, 42, 40, 39, 31, 25, 23, 20, 18, 14,
    12, 11, 10, 9, 8, 7, 5, 3, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173,
    173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173,
    173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173,
    173, 173, 173, 173, 173, 173, 173,
];

pub(crate) fn fits_parser_yylex(
    yylval_param: &mut FITS_PARSER_YYSTYPE,
    yyscanner: &mut yyguts_t,
) -> c_int {
    unsafe {
        let mut yy_amount_of_matched_text: c_int = 0;
        let mut yy_next_state: yy_state_type = 0;
        let mut current_block: u64;
        let mut yy_current_state: yy_state_type = 0;
        let mut yy_cp: *mut c_char = core::ptr::null_mut::<c_char>();
        let mut yy_bp: *mut c_char = core::ptr::null_mut::<c_char>();
        let mut yy_act: c_int = 0;

        yyscanner.yylval_r = yylval_param;
        if yyscanner.yy_init == 0 {
            yyscanner.yy_init = 1;
            if yyscanner.yy_start == 0 {
                yyscanner.yy_start = 1;
            }

            if (yyscanner.yyin_r).is_null() {
                yyscanner.yyin_r = STDIN!();
            }

            if (yyscanner.yyout_r).is_null() {
                yyscanner.yyout_r = STDOUT!();
            }

            if (yyscanner.yy_buffer_stack).is_null()
                || (*(yyscanner.yy_buffer_stack).add(yyscanner.yy_buffer_stack_top)).is_none()
            {
                fits_parser_yyensure_buffer_stack(yyscanner);

                let b =
                    fits_parser_yy_create_buffer(yyscanner.yyin_r, YY_BUF_SIZE as c_int, yyscanner);
                *(yyscanner.yy_buffer_stack).add(yyscanner.yy_buffer_stack_top) = Some(b);
            }
            fits_parser_yy_load_buffer_state(yyscanner);
        }

        loop {
            yy_cp = yyscanner.yy_c_buf_p;
            *yy_cp = yyscanner.yy_hold_char;
            yy_bp = yy_cp;
            yy_current_state = yyscanner.yy_start;
            '_yy_match: loop {
                loop {
                    let mut yy_c: YY_CHAR = YY_EC[*yy_cp as YY_CHAR as usize];
                    if YY_ACCEPT[yy_current_state as usize] != 0 {
                        yyscanner.yy_last_accepting_state = yy_current_state;
                        yyscanner.yy_last_accepting_cpos = yy_cp;
                    }
                    while c_int::from(
                        YY_CHK[(c_int::from(YY_BASE[yy_current_state as usize]) + c_int::from(yy_c))
                            as usize],
                    ) != yy_current_state
                    {
                        yy_current_state = c_int::from(YY_DEF[yy_current_state as usize]);
                        if yy_current_state >= 174 as c_int {
                            yy_c = YY_META[yy_c as usize];
                        }
                    }
                    yy_current_state = yy_state_type::from(
                        YY_NXT[(c_int::from(YY_BASE[yy_current_state as usize]) + c_int::from(yy_c))
                            as usize],
                    );
                    yy_cp = yy_cp.offset(1);
                    if c_int::from(YY_BASE[yy_current_state as usize]) == 416 as c_int {
                        break;
                    }
                }
                '_yy_find_action: loop {
                    yy_act = c_int::from(YY_ACCEPT[yy_current_state as usize]);
                    if yy_act == 0 {
                        yy_cp = yyscanner.yy_last_accepting_cpos;
                        yy_current_state = yyscanner.yy_last_accepting_state;
                        yy_act = c_int::from(YY_ACCEPT[yy_current_state as usize]);
                    }
                    yyscanner.yytext_r = yy_bp;
                    yyscanner.yyleng_r = yy_cp.offset_from(yy_bp) as c_long as c_int;
                    yyscanner.yy_hold_char = *yy_cp;
                    *yy_cp = 0;
                    yyscanner.yy_c_buf_p = yy_cp;
                    loop {
                        match yy_act {
                            0 => {
                                *yy_cp = yyscanner.yy_hold_char;
                                yy_cp = yyscanner.yy_last_accepting_cpos;
                                yy_current_state = yyscanner.yy_last_accepting_state;
                                continue '_yy_find_action;
                            }
                            1 => {
                                break '_yy_match;
                            }
                            2 => {
                                let mut len: c_int = 0;
                                len = strlen_safe(yyscanner.yytext()) as c_int;
                                len -= 1;
                                while c_int::from(*(yyscanner.yytext_r).offset(len as isize))
                                    == ' ' as i32
                                {
                                    len -= 1;
                                }
                                if len >= MAX_STRLEN {
                                    let mut errMsg: [c_char; 100] = [0; 100];
                                    (*yyscanner.yyextra_r).status = PARSE_SYNTAX_ERR;
                                    strcpy_safe(
                                        &mut errMsg,
                                        cs!(c"Bit string exceeds maximum length: '"),
                                    );
                                    strncat_safe(&mut errMsg, yyscanner.yytext(), 20);
                                    strcat_safe(&mut errMsg, cs!(c"...'"));
                                    ffpmsg_slice(&errMsg);
                                    (*yyscanner.yylval_r).text_mut()[0] = 0;
                                } else {
                                    strncpy_safe(
                                        (*yyscanner.yylval_r).text_mut(),
                                        &yyscanner.yytext()[1..],
                                        len as usize,
                                    );
                                    (*yyscanner.yylval_r).text_mut()[len as usize] = 0;
                                }
                                return fits_parser_yytokentype::BITSTR as c_int;
                            }
                            3 => {
                                let mut len_0: c_int = 0;
                                let mut tmpstring: [c_char; 256] = [0; 256];
                                let mut bitstring: [c_char; 256] = [0; 256];
                                len_0 = strlen_safe(yyscanner.yytext()) as c_int;
                                if len_0 >= 256 as c_int {
                                    let mut errMsg: [c_char; 100] = [0; 100];
                                    (*yyscanner.yyextra_r).status = PARSE_SYNTAX_ERR;
                                    strcpy_safe(
                                        &mut errMsg,
                                        cs!(c"Bit string exceeds maximum length: '"),
                                    );
                                    strncat_safe(&mut errMsg, yyscanner.yytext(), 20);
                                    strcat_safe(&mut errMsg, cs!(c"...'"));
                                    ffpmsg_slice(&errMsg);
                                    len_0 = 0;
                                } else {
                                    len_0 -= 1;
                                    while c_int::from(*(yyscanner.yytext_r).offset(len_0 as isize))
                                        == ' ' as i32
                                    {
                                        len_0 -= 1;
                                    }
                                    strncpy_safe(
                                        &mut tmpstring,
                                        &yyscanner.yytext()[1..],
                                        len_0 as usize,
                                    );
                                }
                                tmpstring[len_0 as usize] = 0;
                                bitstring[0] = 0;
                                len_0 = 0;

                                let mut bitlen: usize = 0;
                                let bitcap: usize = core::mem::size_of_val(&bitstring) - 1;
                                let mut overflow: i32 = 0;

                                while c_int::from(tmpstring[len_0 as usize]) != 0 {
                                    let mut chunk: &[c_char] = &[];
                                    let mut chunk_len: usize = 0;
                                    match c_int::from(tmpstring[len_0 as usize]) {
                                        48 => {
                                            chunk = cast_slice(OCT_0.as_bytes());
                                            chunk_len = core::mem::size_of_val(chunk);
                                        }
                                        49 => {
                                            chunk = cast_slice(OCT_1.as_bytes());
                                            chunk_len = core::mem::size_of_val(chunk);
                                        }
                                        50 => {
                                            chunk = cast_slice(OCT_2.as_bytes());
                                            chunk_len = core::mem::size_of_val(chunk);
                                        }
                                        51 => {
                                            chunk = cast_slice(OCT_3.as_bytes());
                                            chunk_len = core::mem::size_of_val(chunk);
                                        }
                                        52 => {
                                            chunk = cast_slice(OCT_4.as_bytes());
                                            chunk_len = core::mem::size_of_val(chunk);
                                        }
                                        53 => {
                                            chunk = cast_slice(OCT_5.as_bytes());
                                            chunk_len = core::mem::size_of_val(chunk);
                                        }
                                        54 => {
                                            chunk = cast_slice(OCT_6.as_bytes());
                                            chunk_len = core::mem::size_of_val(chunk);
                                        }
                                        55 => {
                                            chunk = cast_slice(OCT_7.as_bytes());
                                            chunk_len = core::mem::size_of_val(chunk);
                                        }
                                        120 | 88 => {
                                            chunk = cast_slice(OCT_X.as_bytes());
                                            chunk_len = core::mem::size_of_val(chunk);
                                        }
                                        _ => {}
                                    }

                                    if !chunk.is_empty() {
                                        if bitlen + chunk_len > bitcap {
                                            let mut errMsg: [c_char; 100] = [0; 100];
                                            (*yyscanner.yyextra_r).status = PARSE_SYNTAX_ERR;
                                            strcpy_safe(
                                                &mut errMsg,
                                                cs!(c"Bit string exceeds maximum length: '"),
                                            );
                                            strncat_safe(&mut errMsg, yyscanner.yytext(), 20);
                                            strcat_safe(&mut errMsg, cs!(c"...'"));
                                            ffpmsg_slice(&errMsg);
                                            overflow = 1;
                                            break;
                                        }
                                        bitstring[bitlen..(bitlen + chunk_len)]
                                            .copy_from_slice(chunk);
                                        bitlen += chunk_len;
                                        bitstring[bitlen] = 0;
                                    }

                                    len_0 += 1;
                                }
                                if overflow != 0 {
                                    (*yyscanner.yylval_r).text_mut()[0] = 0;
                                } else {
                                    strcpy_safe((*yyscanner.yylval_r).text_mut(), &bitstring);
                                }
                                return fits_parser_yytokentype::BITSTR as c_int;
                            }
                            4 => {
                                let mut len_1: c_int = 0;
                                let mut tmpstring_0: [c_char; 256] = [0; 256];
                                let mut bitstring_0: [c_char; 256] = [0; 256];
                                len_1 = strlen_safe(yyscanner.yytext()) as c_int;
                                if len_1 >= 256 as c_int {
                                    let mut errMsg_0: [c_char; 100] = [0; 100];
                                    (*yyscanner.yyextra_r).status = PARSE_SYNTAX_ERR;
                                    strcpy_safe(
                                        &mut errMsg_0,
                                        cs!(c"Hex string exceeds maximum length: '"),
                                    );
                                    strncat_safe(&mut errMsg_0, yyscanner.yytext(), 20);
                                    strcat_safe(&mut errMsg_0, cs!(c"...'"));
                                    ffpmsg_slice(&errMsg_0);
                                    len_1 = 0;
                                } else {
                                    len_1 -= 1;
                                    while c_int::from(*(yyscanner.yytext_r).offset(len_1 as isize))
                                        == ' ' as i32
                                    {
                                        len_1 -= 1;
                                    }
                                    strncpy_safe(
                                        &mut tmpstring_0,
                                        &yyscanner.yytext()[1..],
                                        len_1 as usize,
                                    );
                                }
                                tmpstring_0[len_1 as usize] = 0;
                                bitstring_0[0] = 0;
                                len_1 = 0;

                                let mut bitlen: usize = 0;
                                let bitcap: usize = core::mem::size_of_val(&bitstring_0) - 1;
                                let mut overflow: i32 = 0;
                                while c_int::from(tmpstring_0[len_1 as usize]) != 0 {
                                    let mut chunk: &[c_char] = &[];
                                    let mut chunk_len: usize = 0;
                                    match c_int::from(tmpstring_0[len_1 as usize]) {
                                        48 => {
                                            chunk = cast_slice(HEX_0.as_bytes());
                                            chunk_len = core::mem::size_of_val(chunk);
                                        }
                                        49 => {
                                            chunk = cast_slice(HEX_1.as_bytes());
                                            chunk_len = core::mem::size_of_val(chunk);
                                        }
                                        50 => {
                                            chunk = cast_slice(HEX_2.as_bytes());
                                            chunk_len = core::mem::size_of_val(chunk);
                                        }
                                        51 => {
                                            chunk = cast_slice(HEX_3.as_bytes());
                                            chunk_len = core::mem::size_of_val(chunk);
                                        }
                                        52 => {
                                            chunk = cast_slice(HEX_4.as_bytes());
                                            chunk_len = core::mem::size_of_val(chunk);
                                        }
                                        53 => {
                                            chunk = cast_slice(HEX_5.as_bytes());
                                            chunk_len = core::mem::size_of_val(chunk);
                                        }
                                        54 => {
                                            chunk = cast_slice(HEX_6.as_bytes());
                                            chunk_len = core::mem::size_of_val(chunk);
                                        }
                                        55 => {
                                            chunk = cast_slice(HEX_7.as_bytes());
                                            chunk_len = core::mem::size_of_val(chunk);
                                        }
                                        56 => {
                                            chunk = cast_slice(HEX_8.as_bytes());
                                            chunk_len = core::mem::size_of_val(chunk);
                                        }
                                        57 => {
                                            chunk = cast_slice(HEX_9.as_bytes());
                                            chunk_len = core::mem::size_of_val(chunk);
                                        }
                                        97 | 65 => {
                                            chunk = cast_slice(HEX_A.as_bytes());
                                            chunk_len = core::mem::size_of_val(chunk);
                                        }
                                        98 | 66 => {
                                            chunk = cast_slice(HEX_B.as_bytes());
                                            chunk_len = core::mem::size_of_val(chunk);
                                        }
                                        99 | 67 => {
                                            chunk = cast_slice(HEX_C.as_bytes());
                                            chunk_len = core::mem::size_of_val(chunk);
                                        }
                                        100 | 68 => {
                                            chunk = cast_slice(HEX_D.as_bytes());
                                            chunk_len = core::mem::size_of_val(chunk);
                                        }
                                        101 | 69 => {
                                            chunk = cast_slice(HEX_E.as_bytes());
                                            chunk_len = core::mem::size_of_val(chunk);
                                        }
                                        102 | 70 => {
                                            chunk = cast_slice(HEX_F.as_bytes());
                                            chunk_len = core::mem::size_of_val(chunk);
                                        }
                                        120 | 88 => {
                                            chunk = cast_slice(HEX_X.as_bytes());
                                            chunk_len = core::mem::size_of_val(chunk);
                                        }
                                        _ => {}
                                    }

                                    if !chunk.is_empty() {
                                        if bitlen + chunk_len > bitcap {
                                            let mut errMsg: [c_char; 100] = [0; 100];
                                            (*yyscanner.yyextra_r).status = PARSE_SYNTAX_ERR;
                                            strcpy_safe(
                                                &mut errMsg,
                                                cs!(c"Hex string exceeds maximum length: '"),
                                            );
                                            strncat_safe(&mut errMsg, yyscanner.yytext(), 20);
                                            strcat_safe(&mut errMsg, cs!(c"...'"));
                                            ffpmsg_slice(&errMsg);
                                            overflow = 1;
                                            break;
                                        }
                                        bitstring_0[bitlen..(bitlen + chunk_len)]
                                            .copy_from_slice(chunk);
                                        bitlen += chunk_len;
                                        bitstring_0[bitlen] = 0;
                                    }

                                    len_1 += 1;
                                }
                                if overflow != 0 {
                                    (*yyscanner.yylval_r).text_mut()[0] = 0;
                                } else {
                                    strcpy_safe((*yyscanner.yylval_r).text_mut(), &bitstring_0);
                                }
                                return fits_parser_yytokentype::BITSTR as c_int;
                            }
                            5 => {
                                let mut constval: c_long = 0;
                                let mut p: *mut c_char = core::ptr::null_mut::<c_char>();
                                p = &mut *(yyscanner.yytext_r).offset(2 as c_int as isize)
                                    as *mut c_char;
                                while *p != 0 {
                                    constval = constval << 1
                                        | c_int::from(c_int::from(*p) == '1' as i32) as c_long;
                                    p = p.offset(1);
                                }
                                (*yyscanner.yylval_r) = FITS_PARSER_YYSTYPE::Long(constval);
                                return fits_parser_yytokentype::LONG as c_int;
                            }
                            6 => {
                                let mut constval_0: c_long = 0;
                                let mut p_0: *mut c_char = core::ptr::null_mut::<c_char>();
                                p_0 = &mut *(yyscanner.yytext_r).offset(2 as c_int as isize)
                                    as *mut c_char;
                                while *p_0 != 0 {
                                    constval_0 = constval_0 << 3 as c_int
                                        | c_long::from(c_int::from(*p_0) - '0' as i32);
                                    p_0 = p_0.offset(1);
                                }
                                (*yyscanner.yylval_r) = FITS_PARSER_YYSTYPE::Long(constval_0);
                                return fits_parser_yytokentype::LONG as c_int;
                            }
                            7 => {
                                let mut constval_1: c_long = 0;
                                let mut p_1: *mut c_char = core::ptr::null_mut::<c_char>();
                                p_1 = &mut *(yyscanner.yytext_r).offset(2 as c_int as isize)
                                    as *mut c_char;
                                while *p_1 != 0 {
                                    let v: c_int = if isdigit_safe(*p_1) {
                                        c_int::from(*p_1) - '0' as i32
                                    } else {
                                        c_int::from(tolower(*p_1)) - 'a' as i32 + 10 as c_int
                                    };
                                    constval_1 = constval_1 << 4 as c_int | c_long::from(v);
                                    p_1 = p_1.offset(1);
                                }
                                (*yyscanner.yylval_r) = FITS_PARSER_YYSTYPE::Long(constval_1);
                                return fits_parser_yytokentype::LONG as c_int;
                            }
                            8 => {
                                (*yyscanner.yylval_r) =
                                    FITS_PARSER_YYSTYPE::Long(atol(yyscanner.yytext_r));
                                return fits_parser_yytokentype::LONG as c_int;
                            }
                            9 => {
                                if c_int::from(*(yyscanner.yytext_r).offset(0)) == 't' as i32
                                    || c_int::from(*(yyscanner.yytext_r).offset(0)) == 'T' as i32
                                {
                                    (*yyscanner.yylval_r) = FITS_PARSER_YYSTYPE::Logical(1);
                                } else {
                                    (*yyscanner.yylval_r) = FITS_PARSER_YYSTYPE::Logical(0);
                                }
                                return fits_parser_yytokentype::BOOLEAN as c_int;
                            }
                            10 => {
                                (*yyscanner.yylval_r) =
                                    FITS_PARSER_YYSTYPE::Double(atof(yyscanner.yytext_r));
                                return fits_parser_yytokentype::DOUBLE as c_int;
                            }
                            11 => {
                                if {
                                    let s1 = yyscanner.yytext();
                                    fits_strcasecmp(s1, cs!(c"#PI"))
                                } == 0
                                {
                                    (*yyscanner.yylval_r) =
                                        FITS_PARSER_YYSTYPE::Double(4.0 * (1.0_f64).atan());
                                    return fits_parser_yytokentype::DOUBLE as c_int;
                                } else if {
                                    let s1 = yyscanner.yytext();
                                    fits_strcasecmp(s1, cs!(c"#E"))
                                } == 0
                                {
                                    (*yyscanner.yylval_r) =
                                        FITS_PARSER_YYSTYPE::Double((1.0_f64).exp());
                                    return fits_parser_yytokentype::DOUBLE as c_int;
                                } else if {
                                    let s1 = yyscanner.yytext();
                                    fits_strcasecmp(s1, cs!(c"#DEG"))
                                } == 0
                                {
                                    (*yyscanner.yylval_r) =
                                        FITS_PARSER_YYSTYPE::Double(4.0 * (1.0_f64).atan() / 180.0);
                                    return fits_parser_yytokentype::DOUBLE as c_int;
                                } else if {
                                    let s1 = yyscanner.yytext();
                                    fits_strcasecmp(s1, cs!(c"#ROW"))
                                } == 0
                                {
                                    return fits_parser_yytokentype::ROWREF as c_int;
                                } else if {
                                    let s1 = yyscanner.yytext();
                                    fits_strcasecmp(s1, cs!(c"#NULL"))
                                } == 0
                                {
                                    return fits_parser_yytokentype::NULLREF as c_int;
                                } else if {
                                    let s1 = yyscanner.yytext();
                                    fits_strcasecmp(s1, cs!(c"#SNULL"))
                                } == 0
                                {
                                    return fits_parser_yytokentype::SNULLREF as c_int;
                                } else {
                                    let mut len_2: c_int = 0;
                                    let mut result: c_int = 0;
                                    if c_int::from(*(yyscanner.yytext_r).offset(1)) == '$' as i32 {
                                        len_2 = (strlen_safe(yyscanner.yytext())).wrapping_sub(3)
                                            as c_int;
                                        if len_2 >= MAX_STRLEN - 1 {
                                            let mut errMsg: [c_char; 100] = [0; 100];
                                            (*yyscanner.yyextra_r).status = PARSE_SYNTAX_ERR;
                                            strcpy_safe(
                                                &mut errMsg,
                                                cs!(c"Keyword string exceeds maximum length: '"),
                                            );
                                            strncat_safe(&mut errMsg, yyscanner.yytext(), 20);
                                            strcat_safe(&mut errMsg, cs!(c"...'"));
                                            ffpmsg_slice(&errMsg);
                                            (*yyscanner.yylval_r).text_mut()[0] = 0;
                                        } else {
                                            (*yyscanner.yylval_r).text_mut()[0] =
                                                '#' as i32 as c_char;
                                            strncpy_safe(
                                                &mut (*yyscanner.yylval_r).text_mut()[1..],
                                                &yyscanner.yytext()[2..],
                                                len_2 as usize,
                                            );
                                            (*yyscanner.yylval_r).text_mut()
                                                [(len_2 + 1) as usize] = 0;
                                        }
                                        yyscanner.yytext_r = (*yyscanner.yylval_r).text_mut_ptr();
                                    }

                                    let yytext_r_slice = cast_slice(
                                        CStr::from_ptr(yyscanner.yytext_r).to_bytes_with_nul(),
                                    );

                                    result = ((*yyscanner.yyextra_r).getData)
                                        .expect("non-null function pointer")(
                                        &mut *yyscanner.yyextra_r,
                                        yytext_r_slice,
                                        yyscanner
                                            .yylval_r
                                            .cast::<FITS_PARSER_YYSTYPE>()
                                            .as_mut()
                                            .unwrap(),
                                    );
                                    return result;
                                }
                            }
                            12 => {
                                let mut len_3: c_int = 0;
                                len_3 = (strlen_safe(yyscanner.yytext())).wrapping_sub(2) as c_int;
                                if len_3 >= 256 as c_int {
                                    let mut errMsg_1: [c_char; 100] = [0; 100];
                                    (*yyscanner.yyextra_r).status = PARSE_SYNTAX_ERR;
                                    strcpy_safe(
                                        &mut errMsg_1,
                                        cs!(c"String exceeds maximum length: '"),
                                    );
                                    strncat_safe(&mut errMsg_1, &yyscanner.yytext()[1..], 20);
                                    strcat_safe(&mut errMsg_1, cs!(c"...'"));
                                    ffpmsg_slice(&errMsg_1);
                                    len_3 = 0;
                                } else {
                                    strncpy_safe(
                                        (*yyscanner.yylval_r).text_mut(),
                                        &yyscanner.yytext()[1..],
                                        len_3 as usize,
                                    );
                                }
                                (*yyscanner.yylval_r).text_mut()[len_3 as usize] = 0;
                                return fits_parser_yytokentype::STRING as c_int;
                            }
                            13 => {
                                let mut len_4: c_int = 0;
                                let mut dtype: c_int = 0;
                                len_4 = strlen_safe(yyscanner.yytext()) as c_int;
                                if len_4 >= MAX_STRLEN {
                                    let mut errMsg: [c_char; 100] = [0; 100];
                                    (*yyscanner.yyextra_r).status = PARSE_SYNTAX_ERR;
                                    strcpy_safe(
                                        &mut errMsg,
                                        cs!(c"Variable exceeds maximum length: '"),
                                    );
                                    strncat_safe(&mut errMsg, yyscanner.yytext(), 20);
                                    strcat_safe(&mut errMsg, cs!(c"...'"));
                                    ffpmsg_slice(&errMsg);
                                    (*yyscanner.yylval_r).text_mut()[0] = 0;
                                    yyscanner.yytext_r = (*yyscanner.yylval_r).text_mut_ptr();
                                } else if c_int::from(*(yyscanner.yytext_r).offset(0)) == '$' as i32
                                {
                                    len_4 =
                                        (strlen_safe(yyscanner.yytext())).wrapping_sub(2) as c_int;
                                    strncpy_safe(
                                        (*yyscanner.yylval_r).text_mut(),
                                        &yyscanner.yytext()[1..],
                                        len_4 as usize,
                                    );
                                    (*yyscanner.yylval_r).text_mut()[len_4 as usize] = 0;
                                    yyscanner.yytext_r = (*yyscanner.yylval_r).text_mut_ptr();
                                }

                                let yytext_r_slice = cast_slice(
                                    CStr::from_ptr(yyscanner.yytext_r).to_bytes_with_nul(),
                                );

                                dtype = fits_parser_yyGetVariable(
                                    &mut *yyscanner.yyextra_r,
                                    yytext_r_slice,
                                    yyscanner.yylval_r.as_mut().unwrap(),
                                );
                                return dtype;
                            }
                            14 => {
                                let mut len: usize = strlen_safe(yyscanner.yytext());
                                let fname = (*yyscanner.yylval_r).text_mut();

                                if len >= MAX_STRLEN as usize {
                                    let mut errMsg: [c_char; 100] = [0; 100];
                                    (*yyscanner.yyextra_r).status = PARSE_SYNTAX_ERR;
                                    strcpy_safe(
                                        &mut errMsg,
                                        cs!(c"Function exceeds maximum length: '"),
                                    );
                                    strncat_safe(&mut errMsg, yyscanner.yytext(), 20);
                                    strcat_safe(&mut errMsg, cs!(c"...'"));
                                    ffpmsg_slice(&errMsg);
                                    len = 0;
                                    fname[0] = 0;
                                } else {
                                    len = 0;
                                    loop {
                                        fname[len] = toupper(*yyscanner.yytext_r.add(len));
                                        if fname[len] == 0 {
                                            break;
                                        }
                                        len += 1;
                                    }
                                }
                                if FSTRCMP(fname, cs!(c"BOX(")) == 0
                                    || FSTRCMP(fname, cs!(c"CIRCLE(")) == 0
                                    || FSTRCMP(fname, cs!(c"ELLIPSE(")) == 0
                                    || FSTRCMP(fname, cs!(c"NEAR(")) == 0
                                    || FSTRCMP(fname, cs!(c"ISNULL(")) == 0
                                {
                                    /* Return type is always boolean  */
                                    return fits_parser_yytokentype::BFUNCTION as c_int;
                                } else if FSTRCMP(fname, cs!(c"GTIFILTER(")) == 0 {
                                    return fits_parser_yytokentype::GTIFILTER as c_int;
                                } else if FSTRCMP(fname, cs!(c"GTIOVERLAP(")) == 0 {
                                    return fits_parser_yytokentype::GTIOVERLAP as c_int;
                                } else if FSTRCMP(fname, cs!(c"GTIFIND(")) == 0 {
                                    return fits_parser_yytokentype::GTIFIND as c_int;
                                } else if FSTRCMP(fname, cs!(c"REGFILTER(")) == 0 {
                                    return fits_parser_yytokentype::REGFILTER as c_int;
                                } else if FSTRCMP(fname, cs!(c"STRSTR(")) == 0 {
                                    return fits_parser_yytokentype::IFUNCTION as c_int; /* Returns integer */
                                } else {
                                    return fits_parser_yytokentype::FUNCTION as c_int;
                                }
                            }
                            15 => {
                                return fits_parser_yytokentype::INTCAST as c_int;
                            }
                            16 => {
                                return fits_parser_yytokentype::FLTCAST as c_int;
                            }
                            17 => {
                                return fits_parser_yytokentype::POWER as c_int;
                            }
                            18 => return fits_parser_yytokentype::NOT as c_int,
                            19 => return fits_parser_yytokentype::OR as c_int,
                            20 => return fits_parser_yytokentype::AND as c_int,
                            21 => return fits_parser_yytokentype::EQ as c_int,
                            22 => return fits_parser_yytokentype::NE as c_int,
                            23 => return fits_parser_yytokentype::GT as c_int,
                            24 => return fits_parser_yytokentype::LT as c_int,
                            25 => return fits_parser_yytokentype::GTE as c_int,
                            26 => return fits_parser_yytokentype::LTE as c_int,
                            27 => return fits_parser_yytokentype::XOR as c_int,
                            28 => return '\n' as i32,
                            29 => {
                                return c_int::from(*(yyscanner.yytext_r).offset(0));
                            }
                            30 => {
                                fwrite(
                                    yyscanner.yytext_r as *const c_void,
                                    yyscanner.yyleng_r as usize,
                                    1,
                                    yyscanner.yyout_r,
                                );
                                break '_yy_match;
                            }
                            x if x == YY_END_OF_BUFFER + INITIAL + 1 => return YY_NULL,
                            x if x == YY_END_OF_BUFFER => {
                                /* Amount of text matched not including the EOB char. */
                                yy_amount_of_matched_text =
                                    yy_cp.offset_from(yyscanner.yytext_r) as c_long as c_int - 1;

                                /* Undo the effects of YY_DO_BEFORE_ACTION. */
                                *yy_cp = yyscanner.yy_hold_char;

                                let top_state = (*(yyscanner.yy_buffer_stack)
                                    .add(yyscanner.yy_buffer_stack_top))
                                .as_deref_mut()
                                .unwrap();
                                /* We're scanning a new file or input source.  It's
                                 * possible that this happened because the user
                                 * just pointed yyin at a new source and called
                                 * yylex().  If so, then we have to assure
                                 * consistency between YY_CURRENT_BUFFER and our
                                 * globals.  Here is the right place to do so, because
                                 * this is the first action (other than possibly a
                                 * back-up) that will match for the new input source.
                                 */
                                if top_state.yy_buffer_status == YY_BUFFER_NEW {
                                    yyscanner.yy_n_chars = top_state.yy_n_chars;
                                    let fresh4 = &mut top_state.yy_input_file;
                                    *fresh4 = yyscanner.yyin_r;
                                    top_state.yy_buffer_status = YY_BUFFER_NORMAL;
                                }
                                /* Note that here we test for yy_c_buf_p "<=" to the position
                                 * of the first EOB in the buffer, since yy_c_buf_p will
                                 * already have been incremented past the NUL character
                                 * (since all states make transitions on EOB to the
                                 * end-of-buffer state).  Contrast this with the test
                                 * in input().
                                 */
                                if yyscanner.yy_c_buf_p
                                    <= (top_state.yy_ch_buf).as_deref_mut().unwrap()
                                        [yyscanner.yy_n_chars as usize..]
                                        .as_mut_ptr()
                                {
                                    /* This was really a NUL. */
                                    yy_next_state = 0;
                                    yyscanner.yy_c_buf_p = (yyscanner.yytext_r)
                                        .offset(yy_amount_of_matched_text as isize);
                                    yy_current_state = yy_get_previous_state(yyscanner);

                                    /* Okay, we're now positioned to make the NUL
                                     * transition.  We couldn't have
                                     * yy_get_previous_state() go ahead and do it
                                     * for us because it doesn't know how to deal
                                     * with the possibility of jamming (and we don't
                                     * want to build jamming into it because then it
                                     * will run more slowly).
                                     */
                                    yy_next_state = yy_try_NUL_trans(yy_current_state, yyscanner);
                                    yy_bp = (yyscanner.yytext_r).offset(0);
                                    if yy_next_state != 0 {
                                        current_block = 2606663910910355487;
                                        break;
                                    } else {
                                        current_block = 7986280648684494582;
                                        break;
                                    }
                                } else {
                                    match yy_get_next_buffer(yyscanner) {
                                        EOB_ACT_END_OF_FILE => {
                                            yyscanner.yy_did_buffer_switch_on_eof = 0;
                                            if fits_parser_yywrap() != 0 {
                                                /* Note: because we've taken care in
                                                 * yy_get_next_buffer() to have set up
                                                 * yytext, we can now set up
                                                 * yy_c_buf_p so that if some total
                                                 * hoser (like flex itself) wants to
                                                 * call the scanner after we return the
                                                 * YY_NULL, it'll still work - another
                                                 * YY_NULL will get returned.
                                                 */
                                                yyscanner.yy_c_buf_p =
                                                    (yyscanner.yytext_r).offset(0);
                                                /* YY_STATE_EOF(YY_START) */
                                                yy_act = YY_END_OF_BUFFER
                                                    + (yyscanner.yy_start - 1) / 2
                                                    + 1;
                                            } else {
                                                if yyscanner.yy_did_buffer_switch_on_eof == 0 {
                                                    fits_parser_yyrestart(
                                                        yyscanner.yyin_r,
                                                        yyscanner,
                                                    );
                                                }
                                                break '_yy_match;
                                            }
                                        }
                                        EOB_ACT_CONTINUE_SCAN => {
                                            yyscanner.yy_c_buf_p = (yyscanner.yytext_r)
                                                .offset(yy_amount_of_matched_text as isize);
                                            yy_current_state = yy_get_previous_state(yyscanner);
                                            yy_cp = yyscanner.yy_c_buf_p;
                                            yy_bp = (yyscanner.yytext_r).offset(0);
                                            break '_yy_find_action;
                                        }
                                        EOB_ACT_LAST_MATCH => {
                                            let top_state = (*(yyscanner.yy_buffer_stack)
                                                .add(yyscanner.yy_buffer_stack_top))
                                            .as_deref_mut()
                                            .unwrap();

                                            yyscanner.yy_c_buf_p =
                                                (top_state.yy_ch_buf).as_deref_mut().unwrap()
                                                    [yyscanner.yy_n_chars as usize..]
                                                    .as_mut_ptr();
                                            yy_current_state = yy_get_previous_state(yyscanner);
                                            yy_cp = yyscanner.yy_c_buf_p;
                                            yy_bp = (yyscanner.yytext_r).offset(0);
                                            continue '_yy_find_action;
                                        }
                                        _ => {
                                            break '_yy_match;
                                        }
                                    }
                                }
                            }
                            _ => {
                                yy_fatal_error(
                                    "fatal flex scanner internal error--no action found",
                                );
                            }
                        }
                    }
                    match current_block {
                        7986280648684494582 => {
                            yy_cp = yyscanner.yy_c_buf_p;
                        }
                        _ => {
                            yyscanner.yy_c_buf_p = (yyscanner.yy_c_buf_p).offset(1);
                            yy_cp = yyscanner.yy_c_buf_p;
                            yy_current_state = yy_next_state;
                            break;
                        }
                    }
                }
            }
        }
    }
}

/* yy_get_next_buffer - try to read in a new buffer
 *
 * Returns a code representing an action:
 *	EOB_ACT_LAST_MATCH -
 *	EOB_ACT_CONTINUE_SCAN - continue scanning from current position
 *	EOB_ACT_END_OF_FILE - end of file
 */
fn yy_get_next_buffer(yyscanner: &mut yyguts_t) -> c_int {
    unsafe {
        let mut source: *mut c_char = yyscanner.yytext_r;
        let mut number_to_move: c_int = 0;
        let mut i: c_int = 0;
        let mut ret_val: c_int = 0;

        let top_state = (*(yyscanner.yy_buffer_stack).add(yyscanner.yy_buffer_stack_top))
            .as_deref_mut()
            .unwrap();
        let mut dest: *mut c_char = top_state.yy_ch_buf.as_deref_mut().unwrap().as_mut_ptr();

        if yyscanner.yy_c_buf_p
            > (top_state.yy_ch_buf).as_deref_mut().unwrap()[(yyscanner.yy_n_chars + 1) as usize..]
                .as_mut_ptr()
        {
            yy_fatal_error("fatal flex scanner internal error--end of buffer missed");
        }

        if top_state.yy_fill_buffer == 0 {
            /* Don't try to fill the buffer, so this is an EOF. */
            if ((yyscanner.yy_c_buf_p).offset_from(yyscanner.yytext_r) as c_long) == 1 {
                /* We matched a single character, the EOB, so
                 * treat this as a final EOF.
                 */
                return EOB_ACT_END_OF_FILE;
            } else {
                /* We matched some text prior to the EOB, first
                 * process it.
                 */
                return EOB_ACT_LAST_MATCH;
            }
        }

        /* Try to read more data. */

        /* First move last chars to start of buffer. */
        number_to_move =
            ((yyscanner.yy_c_buf_p).offset_from(yyscanner.yytext_r) as c_long - 1) as c_int;
        i = 0;
        while i < number_to_move {
            let fresh5 = source;
            source = source.offset(1);
            let fresh6 = dest;
            dest = dest.offset(1);
            *fresh6 = *fresh5;
            i += 1;
        }
        if top_state.yy_buffer_status == YY_BUFFER_EOF_PENDING {
            /* don't do the read, it's not guaranteed to return an EOF,
             * just force an EOF
             */
            yyscanner.yy_n_chars = 0;
            top_state.yy_n_chars = yyscanner.yy_n_chars;
        } else {
            let mut num_to_read: c_int = top_state.yy_buf_size - number_to_move - 1;
            while num_to_read <= 0 {
                /* Not enough room in the buffer - grow it. */
                let yy_c_buf_p_offset: c_int = (yyscanner.yy_c_buf_p)
                    .offset_from((top_state).yy_ch_buf.as_deref().unwrap().as_ptr())
                    as c_long as c_int;
                if (top_state).yy_is_our_buffer {
                    let new_size: c_int = (top_state).yy_buf_size * 2;
                    if new_size <= 0 {
                        (top_state).yy_buf_size += (top_state).yy_buf_size / 8 as c_int;
                    } else {
                        (top_state).yy_buf_size *= 2;
                    }

                    /* Include room in for 2 EOB chars. */
                    let mut v = Vec::new();
                    let size = new_size as usize;
                    if v.try_reserve_exact(size).is_err() {
                        yy_fatal_error("out of dynamic memory in yy_get_next_buffer()");
                    } else {
                        v.resize(size, 0);
                    }

                    for i in (top_state).yy_ch_buf.as_deref().unwrap().iter() {
                        v.push(*i);
                    }

                    (top_state).yy_ch_buf = Some(v.into_boxed_slice());
                } else {
                    /* Can't grow it, we don't own it. */
                    (top_state).yy_ch_buf = None;
                }

                if ((top_state).yy_ch_buf).is_none() {
                    yy_fatal_error("fatal error - scanner input buffer overflow");
                }

                yyscanner.yy_c_buf_p = ((top_state).yy_ch_buf).as_deref_mut().unwrap()
                    [yy_c_buf_p_offset as usize..]
                    .as_mut_ptr();
                num_to_read = top_state.yy_buf_size - number_to_move - 1;
            }
            if num_to_read > YY_READ_BUF_SIZE {
                num_to_read = YY_READ_BUF_SIZE;
            }
            /* Read in more data. */
            /*
               MJT - 13 June 1996
               read from buffer instead of stdin
               (as per old ftools.skel)
            */
            yyscanner.yy_n_chars = expr_read(
                &mut *yyscanner.yyextra_r,
                &mut (top_state.yy_ch_buf).as_deref_mut().unwrap()[number_to_move as usize..],
                num_to_read,
            );
            if yyscanner.yy_n_chars < 0 {
                yy_fatal_error("read() in flex scanner failed");
            }
            top_state.yy_n_chars = yyscanner.yy_n_chars;
        }

        if yyscanner.yy_n_chars == 0 {
            if number_to_move == 0 {
                ret_val = EOB_ACT_END_OF_FILE;
                fits_parser_yyrestart(yyscanner.yyin_r, yyscanner);
            } else {
                ret_val = EOB_ACT_LAST_MATCH;
                top_state.yy_buffer_status = YY_BUFFER_EOF_PENDING;
            }
        } else {
            ret_val = EOB_ACT_CONTINUE_SCAN;
        }

        let top_state = (*(yyscanner.yy_buffer_stack).add(yyscanner.yy_buffer_stack_top))
            .as_deref_mut()
            .unwrap();

        if yyscanner.yy_n_chars + number_to_move > top_state.yy_buf_size {
            /* Extend the array by 50%, plus the number we really need. */
            let new_size_0: c_int =
                yyscanner.yy_n_chars + number_to_move + (yyscanner.yy_n_chars >> 1);

            let mut v = Vec::new();
            let size = new_size_0 as usize;
            if v.try_reserve_exact(size).is_err() {
                yy_fatal_error("out of dynamic memory in yy_get_next_buffer()");
            } else {
                v.resize(size, 0);
            }

            for i in top_state.yy_ch_buf.as_deref().unwrap().iter() {
                v.push(*i);
            }

            top_state.yy_ch_buf = Some(v.into_boxed_slice());

            /* "- 2" to take care of EOB's */
            top_state.yy_buf_size = new_size_0 - 2;
        }

        yyscanner.yy_n_chars += number_to_move;

        (top_state.yy_ch_buf).as_deref_mut().unwrap()[yyscanner.yy_n_chars as usize] =
            YY_END_OF_BUFFER_CHAR;
        (top_state.yy_ch_buf).as_deref_mut().unwrap()[(yyscanner.yy_n_chars + 1) as usize] =
            YY_END_OF_BUFFER_CHAR;
        yyscanner.yytext_r = &mut *(top_state.yy_ch_buf).as_deref_mut().unwrap().as_mut_ptr();
        ret_val
    }
}

/* yy_get_previous_state - get the state just before the EOB char was reached */
fn yy_get_previous_state(yyscanner: &mut yyguts_t) -> yy_state_type {
    unsafe {
        let mut yy_current_state: yy_state_type = 0;
        let mut yy_cp: *mut c_char = core::ptr::null_mut::<c_char>();

        yy_current_state = yyscanner.yy_start;
        yy_cp = yyscanner.yytext_r;
        while yy_cp < yyscanner.yy_c_buf_p {
            let mut yy_c: YY_CHAR = (if c_int::from(*yy_cp) != 0 {
                c_int::from(YY_EC[*yy_cp as YY_CHAR as usize])
            } else {
                1
            }) as YY_CHAR;
            if YY_ACCEPT[yy_current_state as usize] != 0 {
                yyscanner.yy_last_accepting_state = yy_current_state;
                yyscanner.yy_last_accepting_cpos = yy_cp;
            }
            while c_int::from(
                YY_CHK[(c_int::from(YY_BASE[yy_current_state as usize]) + c_int::from(yy_c))
                    as usize],
            ) != yy_current_state
            {
                yy_current_state = c_int::from(YY_DEF[yy_current_state as usize]);
                if yy_current_state >= 174 as c_int {
                    yy_c = YY_META[yy_c as usize];
                }
            }
            yy_current_state = yy_state_type::from(
                YY_NXT[(c_int::from(YY_BASE[yy_current_state as usize]) + c_int::from(yy_c))
                    as usize],
            );
            yy_cp = yy_cp.offset(1);
        }
        yy_current_state
    }
}

/* yy_try_NUL_trans - try to make a transition on the NUL character
 *
 * synopsis
 *	next_state = yy_try_NUL_trans( current_state );
 */
fn yy_try_NUL_trans(
    mut yy_current_state: yy_state_type,
    yyscanner: &mut yyguts_t,
) -> yy_state_type {
    let mut yy_is_jam: c_int = 0;

    let yy_cp: *mut c_char = yyscanner.yy_c_buf_p;
    let mut yy_c: YY_CHAR = 1;
    if YY_ACCEPT[yy_current_state as usize] != 0 {
        yyscanner.yy_last_accepting_state = yy_current_state;
        yyscanner.yy_last_accepting_cpos = yy_cp;
    }
    while c_int::from(
        YY_CHK[(c_int::from(YY_BASE[yy_current_state as usize]) + c_int::from(yy_c)) as usize],
    ) != yy_current_state
    {
        yy_current_state = c_int::from(YY_DEF[yy_current_state as usize]);
        if yy_current_state >= 174 as c_int {
            yy_c = YY_META[yy_c as usize];
        }
    }
    yy_current_state = yy_state_type::from(
        YY_NXT[(c_int::from(YY_BASE[yy_current_state as usize]) + c_int::from(yy_c)) as usize],
    );
    yy_is_jam = c_int::from(yy_current_state == 173 as c_int);
    if yy_is_jam != 0 { 0 } else { yy_current_state }
}

/// Immediately switch to a different input stream.
/// @param input_file A readable stream.
/// @param yyscanner The scanner object.
/// @note This function does not reset the start condition to @c INITIAL .
pub(crate) fn fits_parser_yyrestart(input_file: *mut FILE, yyscanner: &mut yyguts_t) {
    unsafe {
        if yyscanner.yy_buffer_stack.is_null()
            || (*(yyscanner.yy_buffer_stack).add(yyscanner.yy_buffer_stack_top)).is_none()
        {
            fits_parser_yyensure_buffer_stack(yyscanner);

            let fresh8 = &mut *(yyscanner.yy_buffer_stack).add(yyscanner.yy_buffer_stack_top);
            *fresh8 = Some(fits_parser_yy_create_buffer(
                yyscanner.yyin_r,
                YY_BUF_SIZE as c_int,
                yyscanner,
            ));
        }

        if !(yyscanner.yy_buffer_stack).is_null()
            && (*(yyscanner.yy_buffer_stack).add(yyscanner.yy_buffer_stack_top))
                .as_deref_mut()
                .is_some()
        {
            fits_parser_yy_init_buffer(
                (*(yyscanner.yy_buffer_stack).add(yyscanner.yy_buffer_stack_top))
                    .as_deref_mut()
                    .unwrap(),
                input_file,
                yyscanner,
            );
        }

        fits_parser_yy_load_buffer_state(yyscanner);
    }
}

fn fits_parser_yy_load_buffer_state(yyscanner: &mut yyguts_t) {
    unsafe {
        let top_state = (*(yyscanner.yy_buffer_stack).add(yyscanner.yy_buffer_stack_top))
            .as_deref_mut()
            .unwrap();

        yyscanner.yy_n_chars = top_state.yy_n_chars;

        // Convert index to pointer
        let buf_pos_index = top_state.yy_buf_pos;
        yyscanner.yy_c_buf_p = if let Some(ch_buf) = top_state.yy_ch_buf.as_deref_mut() {
            ch_buf[buf_pos_index..].as_mut_ptr()
        } else {
            ptr::null_mut()
        };

        yyscanner.yytext_r = yyscanner.yy_c_buf_p;
        yyscanner.yyin_r = top_state.yy_input_file;
        if !yyscanner.yy_c_buf_p.is_null() {
            yyscanner.yy_hold_char = *yyscanner.yy_c_buf_p;
        }
    }
}

/// Allocate and initialize an input buffer state.
/// @param file A readable stream.
/// @param size The character buffer size in bytes. When in doubt, use @c YY_BUF_SIZE.
/// @param yyscanner The scanner object.
/// @return the allocated buffer state.
pub(crate) fn fits_parser_yy_create_buffer(
    file: *mut FILE,
    size: c_int,
    yyscanner: &mut yyguts_t,
) -> Box<yy_buffer_state> {
    let b = box_try_new(yy_buffer_state::default());

    if b.is_err() {
        yy_fatal_error("out of dynamic memory in yy_create_buffer()");
    }

    let mut b = b.unwrap();
    b.yy_buf_size = size;

    /* yy_ch_buf has to be 2 characters longer than the size given because
     * we need to put in 2 end-of-buffer characters.
     */

    let mut v = Vec::new();
    let size = (b.yy_buf_size + 2) as usize;
    if v.try_reserve_exact(size).is_err() {
        yy_fatal_error("out of dynamic memory in yy_create_buffer()");
    } else {
        v.resize(size, 0);
    }

    b.yy_ch_buf = None;
    b.yy_ch_buf = Some(v.into_boxed_slice());
    b.yy_is_our_buffer = true;

    fits_parser_yy_init_buffer(&mut (b), file, yyscanner);

    b
}

/// Destroy the buffer.
/// @param b a buffer created with yy_create_buffer()
/// @param yyscanner The scanner object.
pub(crate) fn fits_parser_yy_delete_buffer(
    b: Option<Box<yy_buffer_state>>,
    yyscanner: &mut yyguts_t,
) {
    unsafe {
        if b.is_none() {
            return;
        }

        let mut b = b.unwrap();
        let top_stack = match (*(yyscanner.yy_buffer_stack).add(yyscanner.yy_buffer_stack_top))
            .as_deref_mut()
        {
            Some(x) => x,
            None => ptr::null_mut(),
        };

        let mut z = None;

        if ptr::addr_eq(b.as_mut(), top_stack)
        /* Not sure if we should pop here. */
        {
            z = (*(yyscanner.yy_buffer_stack).add(yyscanner.yy_buffer_stack_top)).take();
        }

        if b.yy_is_our_buffer {
            eprintln!()
            // (*b).yy_ch_buf = None;
        }

        // drop b
        drop(b);
    }
}

/* Initializes or reinitializes a buffer.
 * This function is sometimes called more than once on the same buffer,
 * such as during a yyrestart() or at EOF.
 */
fn fits_parser_yy_init_buffer(b: &mut yy_buffer_state, file: *mut FILE, yyscanner: &mut yyguts_t) {
    unsafe {
        let oerrno = errno().0;

        fits_parser_yy_flush_buffer(b, yyscanner);
        b.yy_input_file = file;
        b.yy_fill_buffer = 1;

        let top_state = match (*(yyscanner.yy_buffer_stack).add(yyscanner.yy_buffer_stack_top))
            .as_deref_mut()
        {
            Some(x) => x as *mut yy_buffer_state,
            None => ptr::null_mut(),
        };

        /* If b is the current buffer, then yy_init_buffer was _probably_
         * called from yyrestart() or through yy_get_next_buffer.
         * In that case, we don't want to reset the lineno or column.
         */
        if !ptr::addr_eq(b, top_state) {
            b.yy_bs_lineno = 1;
            b.yy_bs_column = 0;
        }
        b.yy_is_interactive = if !file.is_null() {
            c_int::from(isatty(fileno(file)) > 0)
        } else {
            0
        };
        set_errno(Errno(oerrno));
    }
}

/// Discard all buffered characters. On the next scan, YY_INPUT will be called.
/// @param b the buffer state to be flushed, usually @c YY_CURRENT_BUFFER.
/// @param yyscanner The scanner object.
pub(crate) fn fits_parser_yy_flush_buffer(b: &mut yy_buffer_state, yyscanner: &mut yyguts_t) {
    unsafe {
        b.yy_n_chars = 0;

        /* We always need two end-of-buffer characters.  The first causes
         * a transition to the end-of-buffer state.  The second causes
         * a jam in that state.
         */
        (b.yy_ch_buf).as_deref_mut().unwrap()[0] = YY_END_OF_BUFFER_CHAR;
        (b.yy_ch_buf).as_deref_mut().unwrap()[1] = YY_END_OF_BUFFER_CHAR;
        b.yy_buf_pos = 0; // Set index to start of buffer
        b.yy_at_bol = 1;
        b.yy_buffer_status = YY_BUFFER_NEW;

        let top_state = match (*(yyscanner.yy_buffer_stack).add(yyscanner.yy_buffer_stack_top))
            .as_deref_mut()
        {
            Some(x) => x as *mut yy_buffer_state,
            None => ptr::null_mut(),
        };

        if ptr::addr_eq(b, top_state) {
            fits_parser_yy_load_buffer_state(yyscanner);
        }
    }
}

/// Removes and deletes the top of the stack, if present.
///  The next element becomes the new top.
///  @param yyscanner The scanner object.
pub(crate) fn fits_parser_yypop_buffer_state(yyscanner: &mut yyguts_t) {
    unsafe {
        if (yyscanner.yy_buffer_stack).is_null()
            || (*(yyscanner.yy_buffer_stack).add(yyscanner.yy_buffer_stack_top)).is_none()
        {
            return;
        }

        let z = (*(yyscanner.yy_buffer_stack).add(yyscanner.yy_buffer_stack_top)).take();
        fits_parser_yy_delete_buffer(z, yyscanner);

        (*(yyscanner.yy_buffer_stack).add(yyscanner.yy_buffer_stack_top)) = None;
        if yyscanner.yy_buffer_stack_top > 0 {
            yyscanner.yy_buffer_stack_top = (yyscanner.yy_buffer_stack_top).wrapping_sub(1);
        }

        if (*(yyscanner.yy_buffer_stack).add(yyscanner.yy_buffer_stack_top)).is_some() {
            fits_parser_yy_load_buffer_state(yyscanner);
            yyscanner.yy_did_buffer_switch_on_eof = 1;
        }
    }
}

/* Allocates the stack if it does not exist.
 *  Guarantees space for at least one push.
 */
fn fits_parser_yyensure_buffer_stack(yyscanner: &mut yyguts_t) {
    unsafe {
        let mut num_to_alloc: yy_size_t = 0;

        if (yyscanner.yy_buffer_stack).is_null() {
            /* First allocation is just for 2 elements, since we don't know if this
             * scanner will even need a stack. We use 2 instead of 1 to avoid an
             * immediate realloc on the next call.
             */
            num_to_alloc = 1; /* After all that talk, this was set to 1 anyways... */

            let mut v = Vec::new();
            if v.try_reserve_exact(num_to_alloc).is_err() {
                yy_fatal_error("out of dynamic memory in yyensure_buffer_stack()");
            } else {
                v.resize_with(num_to_alloc, || None);
            }

            let (ptr, len, capacity) = vec_into_raw_parts(v);
            yyscanner.yy_buffer_stack = ptr;
            yyscanner.yy_buffer_stack_max = num_to_alloc;
            yyscanner.yy_buffer_stack_top = 0;
            return;
        }

        if yyscanner.yy_buffer_stack_top
            >= (yyscanner.yy_buffer_stack_max).wrapping_sub(1 as c_int as size_t)
        {
            /* Increase the buffer to prepare for a possible push. */
            let grow_size: yy_size_t = 8 /* arbitrary grow size */;

            num_to_alloc = (yyscanner.yy_buffer_stack_max) + (grow_size);

            let mut v = Vec::from_raw_parts(
                yyscanner.yy_buffer_stack,
                yyscanner.yy_buffer_stack_max,
                yyscanner.yy_buffer_stack_max,
            );

            if v.try_reserve_exact(grow_size).is_err() {
                yy_fatal_error("out of dynamic memory in yyensure_buffer_stack()");
            } else {
                v.resize_with(num_to_alloc, || None);
            }

            let (ptr, len, capacity) = vec_into_raw_parts(v);
            yyscanner.yy_buffer_stack = ptr;
            yyscanner.yy_buffer_stack_max = num_to_alloc;
        }
    }
}

const YY_EXIT_FAILURE: c_int = 2;

fn yy_fatal_error(msg: &str) -> ! {
    eprintln!("{}", msg);
    exit(YY_EXIT_FAILURE);
}

/** Set the user-defined data. This data is never touched by the scanner.
 * @param user_defined The data to be associated with this scanner.
 * @param yyscanner The scanner object.
 */
pub(crate) fn fits_parser_yyset_extra(user_defined: &mut ParseData, yyscanner: &mut yyguts_t) {
    yyscanner.yyextra_r = user_defined;
}

/* yylex_init_extra has the same functionality as yylex_init, but follows the
 * convention of taking the scanner as the last argument. Note however, that
 * this is a *pointer* to a scanner, as it will be allocated by this call (and
 * is the reason, too, why this function also must handle its own declaration).
 * The user defined value in the first argument will be available to yyalloc in
 * the yyextra field.
 */
pub(crate) fn fits_parser_yylex_init_extra(
    yy_user_defined: &mut ParseData,
    ptr_yy_globals: &mut Option<Box<yyguts_t>>,
) -> c_int {
    let mut dummy_yyguts = yyguts_t::default();

    fits_parser_yyset_extra(yy_user_defined, &mut dummy_yyguts);

    let b = box_try_new(yyguts_t::default());

    if b.is_err() {
        set_errno(Errno(ENOMEM));
        return 1;
    }

    let mut b = b.unwrap();

    fits_parser_yyset_extra(yy_user_defined, &mut b);

    *ptr_yy_globals = Some(b);

    yy_init_globals(ptr_yy_globals.as_deref_mut().unwrap())
}

fn yy_init_globals(yyscanner: &mut yyguts_t) -> c_int {
    /* Initialization is the same as for the non-reentrant scanner.
     * This function is called from yylex_destroy(), so don't allocate here.
     */
    yyscanner.yy_buffer_stack = core::ptr::null_mut();
    yyscanner.yy_buffer_stack_top = 0;
    yyscanner.yy_buffer_stack_max = 0;
    yyscanner.yy_c_buf_p = core::ptr::null_mut::<c_char>();
    yyscanner.yy_init = 0;
    yyscanner.yy_start = 0;
    /* Defined in main.c */
    yyscanner.yyin_r = core::ptr::null_mut::<FILE>();
    yyscanner.yyout_r = core::ptr::null_mut::<FILE>();

    /* For future reference: Set errno on error, since we are called by
     * yylex_init()
     */
    0
}

/* yylex_destroy is for both reentrant and non-reentrant scanners. */
pub(crate) fn fits_parser_yylex_destroy(mut yyscanner: Box<yyguts_t>) -> c_int {
    unsafe {
        let yyg = yyscanner.as_mut();

        /* Pop the buffer stack, destroying each element. */
        while !(yyscanner.yy_buffer_stack).is_null()
            && (*(yyscanner.yy_buffer_stack).add(yyscanner.yy_buffer_stack_top)).is_some()
        {
            let z = (*(yyscanner.yy_buffer_stack).add(yyscanner.yy_buffer_stack_top)).take();
            fits_parser_yy_delete_buffer(z, &mut yyscanner);

            (*(yyscanner.yy_buffer_stack).add(yyscanner.yy_buffer_stack_top)) = None;
            fits_parser_yypop_buffer_state(&mut yyscanner);
        }
        /* Destroy the stack itself. */
        free(yyscanner.yy_buffer_stack.cast::<c_void>());
        yyscanner.yy_buffer_stack = core::ptr::null_mut();

        /* Reset the globals. This is important in a non-reentrant scanner so the next time
         * yylex() is called, initialization will occur. */
        yy_init_globals(&mut yyscanner);

        0
    }
}
