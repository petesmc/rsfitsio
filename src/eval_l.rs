#![deny(dead_code)]

use errno::{Errno, errno, set_errno};

use std::ffi::CStr;
use std::process::exit;
use std::ptr;

use bytemuck::cast_slice;
use libc::{ENOMEM, FILE, atof, atol, fileno, free, fwrite, isatty, size_t};

use crate::c_types::{c_char, c_int, c_long, c_short, c_uchar, c_uint, c_void};
use crate::eval_defs::{MAXVARNAME, P_ERROR, ParseData};
use crate::fitsio::PARSE_SYNTAX_ERR;
use crate::fitsio2::FSTRCMP;
use crate::helpers::boxed::box_try_new;
use crate::helpers::vec_raw_parts::vec_into_raw_parts;
use crate::wrappers::{isdigit_safe, strcat_safe, strcpy_safe, strncat_safe};
use crate::{STDIN, STDOUT, cs, eval_tab::*};
use crate::{
    fitscore::{ffpmsg_slice, fits_strcasecmp, fits_strncasecmp},
    wrappers::{strcat, strcpy, strlen, strncat, strncpy, toupper},
};

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

const YY_BUF_SIZE: usize = 16384;

pub type YY_BUFFER_STATE = *mut yy_buffer_state;

pub type yy_size_t = size_t;

/* Holds the entire state of the reentrant scanner. */
#[derive(Default, Copy, Clone)]
pub(crate) struct yyguts_t {
    /* User-defined. Not touched by flex. */
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

pub type yy_state_type = c_int;

pub type YY_CHAR = flex_uint8_t;
pub type uint = c_uint;

pub(crate) fn fits_parser_yywrap() -> c_int {
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
            lParse.status = 431 as c_int;
            strcpy_safe(&mut errMsg, cs!(c"Unable to find data: "));
            strncat_safe(&mut errMsg, varName, MAXVARNAME);
            ffpmsg_slice(&errMsg);
        }
    } else {
        /*  Convert variable type into expression type  */
        match ((lParse.varData)[varNum as usize]).dtype.into() {
            fits_parser_yytokentype::LONG | fits_parser_yytokentype::DOUBLE => {
                dtype = fits_parser_yytokentype::COLUMN as c_int;
            }
            fits_parser_yytokentype::BOOLEAN => {
                dtype = fits_parser_yytokentype::BCOLUMN as c_int;
            }
            fits_parser_yytokentype::STRING => {
                dtype = fits_parser_yytokentype::SCOLUMN as c_int;
            }
            fits_parser_yytokentype::BITSTR => {
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

        thelval.lng = c_long::from(varNum);
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
    0, 0, 0, 31, 29, 1, 28, 18, 29, 29, 29, 29, 29, 29, 29, 10, 8, 8, 24, 29, 23, 13, 13, 13, 13,
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
    31, 37, 38, 31, 31, 39, 31, 31, 1, 1, 40, 41, 42, 1, 43, 44, 24, 45, 46, 47, 48, 29, 49, 31,
    31, 50, 31, 51, 52, 31, 53, 54, 31, 55, 56, 31, 31, 57, 31, 31, 1, 58, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
];
static YY_META: [YY_CHAR; 59] = [
    0, 1, 1, 2, 1, 1, 1, 3, 1, 1, 1, 1, 1, 1, 1, 4, 4, 4, 4, 1, 1, 1, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 1, 1, 5, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 1,
];
static YY_BASE: [flex_int16_t; 182] = [
    0, 0, 0, 412, 413, 409, 413, 390, 404, 401, 400, 398, 396, 34, 392, 70, 114, 16, 383, 46, 382,
    29, 84, 359, 28, 358, 52, 157, 64, 91, 128, 358, 0, 40, 27, 69, 92, 100, 171, 340, 395, 413,
    391, 413, 388, 387, 386, 413, 413, 383, 357, 358, 356, 336, 337, 335, 413, 139, 190, 352, 349,
    71, 111, 135, 347, 348, 330, 327, 59, 64, 116, 325, 323, 175, 0, 59, 120, 326, 0, 413, 413,
    413, 153, 184, 0, 202, 209, 210, 219, 351, 220, 228, 229, 211, 230, 240, 413, 221, 246, 254,
    263, 264, 265, 266, 239, 275, 413, 346, 342, 310, 313, 309, 289, 292, 288, 275, 317, 327, 326,
    325, 324, 323, 322, 298, 320, 297, 287, 317, 315, 314, 312, 311, 310, 249, 289, 243, 294, 298,
    134, 246, 0, 413, 285, 413, 413, 288, 308, 309, 261, 261, 256, 221, 215, 246, 241, 223, 218,
    213, 208, 197, 413, 166, 160, 413, 128, 122, 150, 154, 105, 101, 96, 413, 84, 413, 351, 354,
    359, 364, 366, 368, 373, 89,
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
static YY_NXT: [flex_int16_t; 472] = [
    0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 4, 14, 4, 15, 16, 17, 17, 17, 18, 19, 20, 21, 22, 23, 23,
    24, 25, 26, 27, 23, 23, 28, 29, 30, 23, 23, 25, 23, 23, 4, 31, 32, 33, 22, 23, 34, 25, 35, 23,
    36, 37, 38, 23, 23, 25, 23, 23, 39, 50, 173, 51, 83, 86, 52, 79, 80, 81, 173, 84, 84, 84, 136,
    173, 137, 137, 137, 137, 87, 53, 98, 54, 84, 55, 57, 58, 58, 58, 58, 88, 90, 97, 59, 140, 84,
    171, 60, 118, 61, 85, 85, 91, 62, 63, 64, 128, 84, 171, 119, 65, 130, 84, 171, 66, 129, 99, 67,
    92, 68, 131, 69, 70, 71, 85, 100, 93, 84, 72, 73, 74, 74, 74, 74, 84, 84, 138, 138, 120, 101,
    75, 75, 85, 84, 94, 94, 94, 103, 102, 121, 138, 138, 172, 104, 57, 115, 115, 115, 115, 76, 75,
    75, 122, 132, 141, 95, 171, 77, 94, 133, 123, 84, 78, 89, 89, 89, 89, 170, 169, 168, 89, 89,
    89, 89, 89, 89, 94, 94, 94, 94, 57, 58, 58, 58, 58, 141, 84, 89, 167, 166, 84, 89, 89, 89, 89,
    89, 58, 58, 58, 58, 142, 94, 96, 141, 84, 89, 75, 75, 85, 85, 141, 141, 141, 160, 80, 81, 105,
    84, 48, 94, 141, 141, 141, 96, 143, 79, 75, 75, 160, 141, 141, 141, 85, 144, 41, 84, 94, 94,
    94, 145, 141, 141, 84, 84, 84, 106, 48, 141, 163, 165, 85, 80, 84, 84, 84, 141, 164, 146, 163,
    81, 94, 84, 84, 84, 141, 141, 141, 141, 143, 79, 144, 41, 84, 84, 162, 161, 141, 139, 94, 84,
    106, 115, 115, 115, 115, 147, 141, 84, 159, 141, 48, 75, 75, 160, 106, 158, 84, 84, 84, 84,
    137, 137, 137, 137, 137, 137, 137, 137, 84, 141, 141, 75, 75, 48, 160, 41, 144, 79, 84, 143,
    81, 84, 80, 157, 156, 106, 155, 41, 144, 79, 143, 81, 80, 154, 153, 152, 151, 150, 149, 148,
    108, 84, 84, 42, 108, 42, 42, 42, 45, 45, 45, 46, 141, 46, 46, 46, 49, 139, 49, 49, 49, 82, 82,
    84, 84, 107, 135, 107, 107, 107, 134, 127, 126, 125, 124, 117, 116, 114, 113, 112, 111, 110,
    109, 43, 47, 173, 108, 43, 40, 106, 96, 84, 84, 81, 79, 56, 43, 48, 47, 44, 43, 41, 40, 173, 3,
    173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173,
    173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173,
    173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173,
    173,
];
static YY_CHK: [flex_int16_t; 472] = [
    0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 13, 17, 13,
    21, 24, 13, 19, 19, 19, 17, 34, 24, 21, 75, 17, 75, 75, 75, 75, 26, 13, 34, 13, 33, 13, 15, 15,
    15, 15, 15, 26, 28, 33, 15, 181, 26, 172, 15, 61, 15, 22, 22, 28, 15, 15, 15, 68, 28, 170, 61,
    15, 69, 35, 169, 15, 68, 35, 15, 29, 15, 69, 15, 15, 15, 22, 35, 29, 22, 15, 16, 16, 16, 16,
    16, 29, 36, 76, 76, 62, 36, 16, 16, 22, 37, 30, 30, 30, 37, 36, 62, 138, 138, 168, 37, 57, 57,
    57, 57, 57, 16, 16, 16, 63, 70, 82, 30, 167, 16, 30, 70, 63, 30, 16, 27, 27, 27, 27, 166, 165,
    164, 27, 27, 27, 27, 27, 27, 30, 38, 38, 38, 73, 73, 73, 73, 73, 83, 82, 27, 162, 161, 27, 27,
    27, 27, 27, 27, 58, 58, 58, 58, 83, 38, 159, 85, 38, 27, 58, 58, 85, 85, 86, 87, 93, 158, 86,
    87, 38, 83, 157, 38, 88, 90, 97, 156, 88, 90, 58, 58, 155, 91, 92, 94, 85, 91, 92, 85, 94, 94,
    94, 93, 104, 95, 86, 87, 93, 95, 154, 98, 153, 152, 85, 98, 88, 90, 97, 99, 151, 97, 150, 99,
    94, 91, 92, 94, 100, 101, 102, 103, 100, 101, 102, 103, 104, 95, 149, 148, 105, 139, 94, 98,
    105, 115, 115, 115, 115, 104, 142, 99, 135, 145, 142, 115, 115, 145, 134, 133, 100, 101, 102,
    103, 136, 136, 136, 136, 137, 137, 137, 137, 105, 146, 147, 115, 115, 146, 147, 132, 131, 130,
    142, 129, 128, 145, 127, 126, 125, 124, 123, 122, 121, 120, 119, 118, 117, 116, 114, 113, 112,
    111, 110, 109, 108, 146, 147, 174, 107, 174, 174, 174, 175, 175, 175, 176, 89, 176, 176, 176,
    177, 77, 177, 177, 177, 178, 178, 179, 179, 180, 72, 180, 180, 180, 71, 67, 66, 65, 64, 60, 59,
    55, 54, 53, 52, 51, 50, 49, 46, 45, 44, 42, 40, 39, 31, 25, 23, 20, 18, 14, 12, 11, 10, 9, 8,
    7, 5, 3, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173,
    173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173,
    173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173,
    173, 173, 173, 173,
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
        let mut yy_cp: *mut c_char = std::ptr::null_mut::<c_char>();
        let mut yy_bp: *mut c_char = std::ptr::null_mut::<c_char>();
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
                    while YY_CHK[(c_int::from(YY_BASE[yy_current_state as usize])
                        + c_int::from(yy_c)) as usize] as c_int
                        != yy_current_state
                    {
                        yy_current_state = c_int::from(YY_DEF[yy_current_state as usize]);
                        if yy_current_state >= 174 as c_int {
                            yy_c = YY_META[yy_c as usize];
                        }
                    }
                    yy_current_state = YY_NXT[(c_int::from(YY_BASE[yy_current_state as usize])
                        + c_int::from(yy_c)) as usize]
                        as yy_state_type;
                    yy_cp = yy_cp.offset(1);
                    if c_int::from(YY_BASE[yy_current_state as usize]) == 413 as c_int {
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
                                len = strlen(yyscanner.yytext_r) as c_int;
                                while c_int::from(*(yyscanner.yytext_r).offset(len as isize))
                                    == ' ' as i32
                                {
                                    len -= 1;
                                }
                                len -= 1;
                                strncpy(
                                    ((*yyscanner.yylval_r).astr).as_mut_ptr(),
                                    &*(yyscanner.yytext_r).offset(1),
                                    len as usize,
                                );
                                (*yyscanner.yylval_r).astr[len as usize] = 0;
                                return fits_parser_yytokentype::BITSTR as c_int;
                            }
                            3 => {
                                let mut len_0: c_int = 0;
                                let mut tmpstring: [c_char; 256] = [0; 256];
                                let mut bitstring: [c_char; 256] = [0; 256];
                                len_0 = strlen(yyscanner.yytext_r) as c_int;
                                if len_0 >= 256 as c_int {
                                    let mut errMsg: [c_char; 100] = [0; 100];
                                    (*yyscanner.yyextra_r).status = 431 as c_int;
                                    strcpy_safe(
                                        &mut errMsg,
                                        cs!(c"Bit string exceeds maximum length: '"),
                                    );
                                    strncat(
                                        errMsg.as_mut_ptr(),
                                        &*(yyscanner.yytext_r).offset(0),
                                        20,
                                    );
                                    strcat_safe(&mut errMsg, cs!(c"...'"));
                                    ffpmsg_slice(&errMsg);
                                    len_0 = 0;
                                } else {
                                    while c_int::from(*(yyscanner.yytext_r).offset(len_0 as isize))
                                        == ' ' as i32
                                    {
                                        len_0 -= 1;
                                    }
                                    len_0 -= 1;
                                    strncpy(
                                        tmpstring.as_mut_ptr(),
                                        &*(yyscanner.yytext_r).offset(1),
                                        len_0 as usize,
                                    );
                                }
                                tmpstring[len_0 as usize] = 0;
                                bitstring[0] = 0;
                                len_0 = 0;
                                while c_int::from(tmpstring[len_0 as usize]) != 0 {
                                    match c_int::from(tmpstring[len_0 as usize]) {
                                        48 => {
                                            strcat(bitstring.as_mut_ptr(), c"000".as_ptr());
                                        }
                                        49 => {
                                            strcat(bitstring.as_mut_ptr(), c"001".as_ptr());
                                        }
                                        50 => {
                                            strcat(bitstring.as_mut_ptr(), c"010".as_ptr());
                                        }
                                        51 => {
                                            strcat(bitstring.as_mut_ptr(), c"011".as_ptr());
                                        }
                                        52 => {
                                            strcat(bitstring.as_mut_ptr(), c"100".as_ptr());
                                        }
                                        53 => {
                                            strcat(bitstring.as_mut_ptr(), c"101".as_ptr());
                                        }
                                        54 => {
                                            strcat(bitstring.as_mut_ptr(), c"110".as_ptr());
                                        }
                                        55 => {
                                            strcat(bitstring.as_mut_ptr(), c"111".as_ptr());
                                        }
                                        120 | 88 => {
                                            strcat(bitstring.as_mut_ptr(), c"xxx".as_ptr());
                                        }
                                        _ => {}
                                    }
                                    len_0 += 1;
                                }
                                strcpy(
                                    ((*yyscanner.yylval_r).astr).as_mut_ptr(),
                                    bitstring.as_mut_ptr(),
                                );
                                return fits_parser_yytokentype::BITSTR as c_int;
                            }
                            4 => {
                                let mut len_1: c_int = 0;
                                let mut tmpstring_0: [c_char; 256] = [0; 256];
                                let mut bitstring_0: [c_char; 256] = [0; 256];
                                len_1 = strlen(yyscanner.yytext_r) as c_int;
                                if len_1 >= 256 as c_int {
                                    let mut errMsg_0: [c_char; 100] = [0; 100];
                                    (*yyscanner.yyextra_r).status = 431 as c_int;
                                    strcpy_safe(
                                        &mut errMsg_0,
                                        cs!(c"Hex string exceeds maximum length: '"),
                                    );
                                    strncat(
                                        errMsg_0.as_mut_ptr(),
                                        &*(yyscanner.yytext_r).offset(0),
                                        20,
                                    );
                                    strcat_safe(&mut errMsg_0, cs!(c"...'"));
                                    ffpmsg_slice(&errMsg_0);
                                    len_1 = 0;
                                } else {
                                    while c_int::from(*(yyscanner.yytext_r).offset(len_1 as isize))
                                        == ' ' as i32
                                    {
                                        len_1 -= 1;
                                    }
                                    len_1 -= 1;
                                    strncpy(
                                        tmpstring_0.as_mut_ptr(),
                                        &*(yyscanner.yytext_r).offset(1),
                                        len_1 as usize,
                                    );
                                }
                                tmpstring_0[len_1 as usize] = 0;
                                bitstring_0[0] = 0;
                                len_1 = 0;
                                while c_int::from(tmpstring_0[len_1 as usize]) != 0 {
                                    match c_int::from(tmpstring_0[len_1 as usize]) {
                                        48 => {
                                            strcat_safe(&mut bitstring_0, cs!(c"0000"));
                                        }
                                        49 => {
                                            strcat_safe(&mut bitstring_0, cs!(c"0001"));
                                        }
                                        50 => {
                                            strcat_safe(&mut bitstring_0, cs!(c"0010"));
                                        }
                                        51 => {
                                            strcat_safe(&mut bitstring_0, cs!(c"0011"));
                                        }
                                        52 => {
                                            strcat_safe(&mut bitstring_0, cs!(c"0100"));
                                        }
                                        53 => {
                                            strcat_safe(&mut bitstring_0, cs!(c"0101"));
                                        }
                                        54 => {
                                            strcat_safe(&mut bitstring_0, cs!(c"0110"));
                                        }
                                        55 => {
                                            strcat_safe(&mut bitstring_0, cs!(c"0111"));
                                        }
                                        56 => {
                                            strcat_safe(&mut bitstring_0, cs!(c"1000"));
                                        }
                                        57 => {
                                            strcat_safe(&mut bitstring_0, cs!(c"1001"));
                                        }
                                        97 | 65 => {
                                            strcat_safe(&mut bitstring_0, cs!(c"1010"));
                                        }
                                        98 | 66 => {
                                            strcat_safe(&mut bitstring_0, cs!(c"1011"));
                                        }
                                        99 | 67 => {
                                            strcat_safe(&mut bitstring_0, cs!(c"1100"));
                                        }
                                        100 | 68 => {
                                            strcat_safe(&mut bitstring_0, cs!(c"1101"));
                                        }
                                        101 | 69 => {
                                            strcat_safe(&mut bitstring_0, cs!(c"1110"));
                                        }
                                        102 | 70 => {
                                            strcat_safe(&mut bitstring_0, cs!(c"1111"));
                                        }
                                        120 | 88 => {
                                            strcat_safe(&mut bitstring_0, cs!(c"xxxx"));
                                        }
                                        _ => {}
                                    }
                                    len_1 += 1;
                                }
                                strcpy_safe(&mut ((*yyscanner.yylval_r).astr), &bitstring_0);
                                return fits_parser_yytokentype::BITSTR as c_int;
                            }
                            5 => {
                                let mut constval: c_long = 0;
                                let mut p: *mut c_char = std::ptr::null_mut::<c_char>();
                                p = &mut *(yyscanner.yytext_r).offset(2 as c_int as isize)
                                    as *mut c_char;
                                while *p != 0 {
                                    constval = constval << 1
                                        | (c_int::from(*p) == '1' as i32) as c_int as c_long;
                                    p = p.offset(1);
                                }
                                (*yyscanner.yylval_r).lng = constval;
                                return fits_parser_yytokentype::LONG as c_int;
                            }
                            6 => {
                                let mut constval_0: c_long = 0;
                                let mut p_0: *mut c_char = std::ptr::null_mut::<c_char>();
                                p_0 = &mut *(yyscanner.yytext_r).offset(2 as c_int as isize)
                                    as *mut c_char;
                                while *p_0 != 0 {
                                    constval_0 = constval_0 << 3 as c_int
                                        | (c_int::from(*p_0) - '0' as i32) as c_long;
                                    p_0 = p_0.offset(1);
                                }
                                (*yyscanner.yylval_r).lng = constval_0;
                                return fits_parser_yytokentype::LONG as c_int;
                            }
                            7 => {
                                let mut constval_1: c_long = 0;
                                let mut p_1: *mut c_char = std::ptr::null_mut::<c_char>();
                                p_1 = &mut *(yyscanner.yytext_r).offset(2 as c_int as isize)
                                    as *mut c_char;
                                while *p_1 != 0 {
                                    let v: c_int = if isdigit_safe(*p_1) {
                                        c_int::from(*p_1) - '0' as i32
                                    } else {
                                        c_int::from(*p_1) - 'a' as i32 + 10 as c_int
                                    };
                                    constval_1 = constval_1 << 4 as c_int | c_long::from(v);
                                    p_1 = p_1.offset(1);
                                }
                                (*yyscanner.yylval_r).lng = constval_1;
                                return fits_parser_yytokentype::LONG as c_int;
                            }
                            8 => {
                                (*yyscanner.yylval_r).lng = atol(yyscanner.yytext_r);
                                return fits_parser_yytokentype::LONG as c_int;
                            }
                            9 => {
                                if c_int::from(*(yyscanner.yytext_r).offset(0)) == 't' as i32
                                    || c_int::from(*(yyscanner.yytext_r).offset(0)) == 'T' as i32
                                {
                                    (*yyscanner.yylval_r).log = 1;
                                } else {
                                    (*yyscanner.yylval_r).log = 0;
                                }
                                return fits_parser_yytokentype::BOOLEAN as c_int;
                            }
                            10 => {
                                (*yyscanner.yylval_r).dbl = atof(yyscanner.yytext_r);
                                return fits_parser_yytokentype::DOUBLE as c_int;
                            }
                            11 => {
                                if {
                                    let s1 = std::slice::from_raw_parts(
                                        yyscanner.yytext_r,
                                        strlen(yyscanner.yytext_r) as usize + 1,
                                    );
                                    fits_strcasecmp(s1, cs!(c"#PI"))
                                } == 0
                                {
                                    (*yyscanner.yylval_r).dbl = 4.0 * (1.0_f64).atan();
                                    return fits_parser_yytokentype::DOUBLE as c_int;
                                } else if {
                                    let s1 = std::slice::from_raw_parts(
                                        yyscanner.yytext_r,
                                        strlen(yyscanner.yytext_r) as usize + 1,
                                    );
                                    fits_strcasecmp(s1, cs!(c"#E"))
                                } == 0
                                {
                                    (*yyscanner.yylval_r).dbl = (1.0_f64).exp();
                                    return fits_parser_yytokentype::DOUBLE as c_int;
                                } else if {
                                    let s1 = std::slice::from_raw_parts(
                                        yyscanner.yytext_r,
                                        strlen(yyscanner.yytext_r) as usize + 1,
                                    );
                                    fits_strcasecmp(s1, cs!(c"#DEG"))
                                } == 0
                                {
                                    (*yyscanner.yylval_r).dbl = 4.0 * (1.0_f64).atan() / 180.0;
                                    return fits_parser_yytokentype::DOUBLE as c_int;
                                } else if {
                                    let s1 = std::slice::from_raw_parts(
                                        yyscanner.yytext_r,
                                        strlen(yyscanner.yytext_r) as usize + 1,
                                    );
                                    fits_strcasecmp(s1, cs!(c"#ROW"))
                                } == 0
                                {
                                    return fits_parser_yytokentype::ROWREF as c_int;
                                } else if {
                                    let s1 = std::slice::from_raw_parts(
                                        yyscanner.yytext_r,
                                        strlen(yyscanner.yytext_r) as usize + 1,
                                    );
                                    fits_strcasecmp(s1, cs!(c"#NULL"))
                                } == 0
                                {
                                    return fits_parser_yytokentype::NULLREF as c_int;
                                } else if {
                                    let s1 = std::slice::from_raw_parts(
                                        yyscanner.yytext_r,
                                        strlen(yyscanner.yytext_r) as usize + 1,
                                    );
                                    fits_strcasecmp(s1, cs!(c"#SNULL"))
                                } == 0
                                {
                                    return fits_parser_yytokentype::SNULLREF as c_int;
                                } else {
                                    let mut len_2: c_int = 0;
                                    let mut result: c_int = 0;
                                    if c_int::from(*(yyscanner.yytext_r).offset(1)) == '$' as i32 {
                                        len_2 =
                                            (strlen(yyscanner.yytext_r)).wrapping_sub(3) as c_int;
                                        (*yyscanner.yylval_r).astr[0] = '#' as i32 as c_char;
                                        strncpy(
                                            ((*yyscanner.yylval_r).astr).as_mut_ptr().offset(1),
                                            &*(yyscanner.yytext_r).offset(2 as c_int as isize),
                                            len_2 as usize,
                                        );
                                        (*yyscanner.yylval_r).astr[(len_2 + 1) as usize] = 0;
                                        yyscanner.yytext_r =
                                            ((*yyscanner.yylval_r).astr).as_mut_ptr();
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
                                len_3 = (strlen(yyscanner.yytext_r)).wrapping_sub(2) as c_int;
                                if len_3 >= 256 as c_int {
                                    let mut errMsg_1: [c_char; 100] = [0; 100];
                                    (*yyscanner.yyextra_r).status = 431 as c_int;
                                    strcpy_safe(
                                        &mut errMsg_1,
                                        cs!(c"String exceeds maximum length: '"),
                                    );
                                    strncat(
                                        errMsg_1.as_mut_ptr(),
                                        &*(yyscanner.yytext_r).offset(1),
                                        20,
                                    );
                                    strcat_safe(&mut errMsg_1, cs!(c"...'"));
                                    ffpmsg_slice(&errMsg_1);
                                    len_3 = 0;
                                } else {
                                    strncpy(
                                        ((*yyscanner.yylval_r).astr).as_mut_ptr(),
                                        &*(yyscanner.yytext_r).offset(1),
                                        len_3 as usize,
                                    );
                                }
                                (*yyscanner.yylval_r).astr[len_3 as usize] = 0;
                                return fits_parser_yytokentype::STRING as c_int;
                            }
                            13 => {
                                let mut len_4: c_int = 0;
                                let mut dtype: c_int = 0;
                                if c_int::from(*(yyscanner.yytext_r).offset(0)) == '$' as i32 {
                                    len_4 = (strlen(yyscanner.yytext_r)).wrapping_sub(2) as c_int;
                                    strncpy(
                                        ((*yyscanner.yylval_r).astr).as_mut_ptr(),
                                        &*(yyscanner.yytext_r).offset(1),
                                        len_4 as usize,
                                    );
                                    (*yyscanner.yylval_r).astr[len_4 as usize] = 0;
                                    yyscanner.yytext_r = ((*yyscanner.yylval_r).astr).as_mut_ptr();
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
                                let mut len: usize = 0;
                                let fname = &mut ((*yyscanner.yylval_r).astr);

                                loop {
                                    fname[len] = toupper(*yyscanner.yytext_r.add(len));
                                    if fname[len] == 0 {
                                        break;
                                    }
                                    len += 1;
                                }
                                if FSTRCMP(fname, cs!(c"BOX(")) == 0
                                    || FSTRCMP(fname, cs!(c"CIRCLE(")) == 0
                                    || FSTRCMP(fname, cs!(c"ELLIPSE(")) == 0
                                    || FSTRCMP(fname, cs!(c"NEAR(")) == 0
                                    || FSTRCMP(fname, cs!(c"ISNULL(")) == 0
                                {
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
                                0;
                                break '_yy_match;
                            }
                            32 => return 0,
                            31 => {
                                yy_amount_of_matched_text =
                                    yy_cp.offset_from(yyscanner.yytext_r) as c_long as c_int - 1;
                                *yy_cp = yyscanner.yy_hold_char;

                                let top_state = (*(yyscanner.yy_buffer_stack)
                                    .add(yyscanner.yy_buffer_stack_top))
                                .as_deref_mut()
                                .unwrap();
                                if top_state.yy_buffer_status == 0 {
                                    yyscanner.yy_n_chars = top_state.yy_n_chars;
                                    let fresh4 = &mut top_state.yy_input_file;
                                    *fresh4 = yyscanner.yyin_r;
                                    top_state.yy_buffer_status = 1;
                                }
                                if yyscanner.yy_c_buf_p
                                    <= (top_state.yy_ch_buf).as_deref_mut().unwrap()
                                        [yyscanner.yy_n_chars as usize..]
                                        .as_mut_ptr()
                                {
                                    yy_next_state = 0;
                                    yyscanner.yy_c_buf_p = (yyscanner.yytext_r)
                                        .offset(yy_amount_of_matched_text as isize);
                                    yy_current_state = yy_get_previous_state(yyscanner);
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
                                        1 => {
                                            yyscanner.yy_did_buffer_switch_on_eof = 0;
                                            if fits_parser_yywrap() != 0 {
                                                yyscanner.yy_c_buf_p =
                                                    (yyscanner.yytext_r).offset(0);
                                                yy_act =
                                                    31 as c_int + (yyscanner.yy_start - 1) / 2 + 1;
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
                                        0 => {
                                            yyscanner.yy_c_buf_p = (yyscanner.yytext_r)
                                                .offset(yy_amount_of_matched_text as isize);
                                            yy_current_state = yy_get_previous_state(yyscanner);
                                            yy_cp = yyscanner.yy_c_buf_p;
                                            yy_bp = (yyscanner.yytext_r).offset(0);
                                            break '_yy_find_action;
                                        }
                                        2 => {
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
            if ((yyscanner.yy_c_buf_p).offset_from(yyscanner.yytext_r) as c_long) == 1 {
                return 1;
            } else {
                return 2;
            }
        }
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
        if top_state.yy_buffer_status == 2 {
            yyscanner.yy_n_chars = 0;
            top_state.yy_n_chars = yyscanner.yy_n_chars;
        } else {
            let mut num_to_read: c_int = top_state.yy_buf_size - number_to_move - 1;
            while num_to_read <= 0 {
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
            if num_to_read > 8192 as c_int {
                num_to_read = 8192 as c_int;
            }
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
                ret_val = 1;
                fits_parser_yyrestart(yyscanner.yyin_r, yyscanner);
            } else {
                ret_val = 2;
                top_state.yy_buffer_status = 2;
            }
        } else {
            ret_val = 0;
        }

        let top_state = (*(yyscanner.yy_buffer_stack).add(yyscanner.yy_buffer_stack_top))
            .as_deref_mut()
            .unwrap();

        if yyscanner.yy_n_chars + number_to_move > top_state.yy_buf_size {
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

            top_state.yy_buf_size = new_size_0 - 2;
        }

        yyscanner.yy_n_chars += number_to_move;

        (top_state.yy_ch_buf).as_deref_mut().unwrap()[yyscanner.yy_n_chars as usize] = 0;
        (top_state.yy_ch_buf).as_deref_mut().unwrap()[(yyscanner.yy_n_chars + 1) as usize] = 0;
        yyscanner.yytext_r = &mut *(top_state.yy_ch_buf).as_deref_mut().unwrap().as_mut_ptr();
        ret_val
    }
}

fn yy_get_previous_state(yyscanner: &mut yyguts_t) -> yy_state_type {
    unsafe {
        let mut yy_current_state: yy_state_type = 0;
        let mut yy_cp: *mut c_char = std::ptr::null_mut::<c_char>();

        yy_current_state = yyscanner.yy_start;
        yy_cp = (yyscanner.yytext_r).offset(0);
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
            while YY_CHK
                [(c_int::from(YY_BASE[yy_current_state as usize]) + c_int::from(yy_c)) as usize]
                as c_int
                != yy_current_state
            {
                yy_current_state = c_int::from(YY_DEF[yy_current_state as usize]);
                if yy_current_state >= 174 as c_int {
                    yy_c = YY_META[yy_c as usize];
                }
            }
            yy_current_state = YY_NXT
                [(c_int::from(YY_BASE[yy_current_state as usize]) + c_int::from(yy_c)) as usize]
                as yy_state_type;
            yy_cp = yy_cp.offset(1);
        }
        yy_current_state
    }
}

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
    while YY_CHK[(c_int::from(YY_BASE[yy_current_state as usize]) + c_int::from(yy_c)) as usize]
        as c_int
        != yy_current_state
    {
        yy_current_state = c_int::from(YY_DEF[yy_current_state as usize]);
        if yy_current_state >= 174 as c_int {
            yy_c = YY_META[yy_c as usize];
        }
    }
    yy_current_state = YY_NXT
        [(c_int::from(YY_BASE[yy_current_state as usize]) + c_int::from(yy_c)) as usize]
        as yy_state_type;
    yy_is_jam = c_int::from(yy_current_state == 173 as c_int);
    if yy_is_jam != 0 { 0 } else { yy_current_state }
}

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

        if ptr::addr_eq(b.as_mut(), top_stack) {
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

pub(crate) fn fits_parser_yy_flush_buffer(b: &mut yy_buffer_state, yyscanner: &mut yyguts_t) {
    unsafe {
        b.yy_n_chars = 0;
        (b.yy_ch_buf).as_deref_mut().unwrap()[0] = 0;
        (b.yy_ch_buf).as_deref_mut().unwrap()[1] = 0;
        b.yy_buf_pos = 0; // Set index to start of buffer
        b.yy_at_bol = 1;
        b.yy_buffer_status = 0;

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
            yyscanner.yy_buffer_stack_top;
        }

        if (*(yyscanner.yy_buffer_stack).add(yyscanner.yy_buffer_stack_top)).is_some() {
            fits_parser_yy_load_buffer_state(yyscanner);
            yyscanner.yy_did_buffer_switch_on_eof = 1;
        }
    }
}

fn fits_parser_yyensure_buffer_stack(yyscanner: &mut yyguts_t) {
    unsafe {
        let mut num_to_alloc: yy_size_t = 0;

        if (yyscanner.yy_buffer_stack).is_null() {
            num_to_alloc = 1;

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
            let grow_size: yy_size_t = 8;
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
    yyscanner.yy_buffer_stack = std::ptr::null_mut();
    yyscanner.yy_buffer_stack_top = 0;
    yyscanner.yy_buffer_stack_max = 0;
    yyscanner.yy_c_buf_p = std::ptr::null_mut::<c_char>();
    yyscanner.yy_init = 0;
    yyscanner.yy_start = 0;
    yyscanner.yyin_r = std::ptr::null_mut::<FILE>();
    yyscanner.yyout_r = std::ptr::null_mut::<FILE>();
    0
}

pub(crate) fn fits_parser_yylex_destroy(mut yyscanner: Box<yyguts_t>) -> c_int {
    unsafe {
        let yyg = yyscanner.as_mut();
        while !(yyscanner.yy_buffer_stack).is_null()
            && (*(yyscanner.yy_buffer_stack).add(yyscanner.yy_buffer_stack_top)).is_some()
        {
            let z = (*(yyscanner.yy_buffer_stack).add(yyscanner.yy_buffer_stack_top)).take();
            fits_parser_yy_delete_buffer(z, &mut yyscanner);

            (*(yyscanner.yy_buffer_stack).add(yyscanner.yy_buffer_stack_top)) = None;
            fits_parser_yypop_buffer_state(&mut yyscanner);
        }
        free(yyscanner.yy_buffer_stack.cast::<c_void>());
        yyscanner.yy_buffer_stack = std::ptr::null_mut();
        yy_init_globals(&mut yyscanner);

        0
    }
}
