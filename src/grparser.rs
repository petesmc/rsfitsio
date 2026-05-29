use core::slice;
use std::env;
use std::ffi::{CStr, c_void};
use std::fs::File;
use std::io::{BufReader, Read};
use std::path::{Path, PathBuf};
use std::ptr;

use crate::aliases::c_api::*;
use crate::c_types::{c_char, c_int, c_long};
use crate::fitscore::{fits_strcasecmp, fits_strncasecmp};
use crate::fitsio::{
    ASCII_TBL, BINARY_TBL, FLEN_KEYWORD, GT_ID_ALL_URI, NGP_BAD_ARG, NGP_EMPTY_CURLINE, NGP_EOF,
    NGP_ERR_FOPEN, NGP_INC_NESTING, NGP_NO_MEMORY, NGP_NUL_PTR, NGP_OK, NGP_READ_ERR,
    NGP_TOKEN_NOT_EXPECT, NGP_UNREAD_QUEUE_FULL, NULL_MSG, OPT_RM_GPT, TDBLCOMPLEX, TDOUBLE,
    TLOGICAL, TLONG, TSTRING, VALUE_UNDEFINED, fitsfile,
};
use crate::relibc::header::stdio::{sscanf_d_c, sscanf_d_n, sscanf_lg_lg_n, sscanf_lg_n};
use crate::wrappers::{strcmp, strcmp_safe, strcpy, strlen, strlen_safe, strncmp, strncpy};
use crate::{FFLOCK, bb, cs, int_snprintf};
use bytemuck::{cast_slice, cast_slice_mut};

use libc::{c_uint, free, malloc, memcmp, memcpy, realloc, size_t};

const NGP_ALLOCCHUNK: c_int = 1000;
const NGP_MAX_INCLUDE: usize = 10; /* include file nesting limit */
const NGP_MAX_COMMENT: usize = 80; /* max size for comment */
const NGP_MAX_NAME: usize = FLEN_KEYWORD; /* max size for KEYWORD (FITS limits it to 8 chars) */
/* except HIERARCH can have longer effective keyword names */
const NGP_MAX_STRING: usize = 80; /* max size for various strings */
const NGP_MAX_ARRAY_DIM: usize = 999; /* max. number of dimensions in array */
const NGP_MAX_FNAME: usize = 1000; /* max size of combined path+fname */
const NGP_MAX_ENVFILES: usize = 10000; /* max size of CFITSIO_INCLUDE_FILES env. variable */

const NGP_TOKEN_UNKNOWN: c_int = -1; /* token type unknown */
const NGP_TOKEN_INCLUDE: c_int = 0; /* \INCLUDE token */
const NGP_TOKEN_GROUP: c_int = 1; /* \GROUP token */
const NGP_TOKEN_END: c_int = 2; /* \END token */
const NGP_TOKEN_XTENSION: c_int = 3; /* XTENSION token */
const NGP_TOKEN_SIMPLE: c_int = 4; /* SIMPLE token */
const NGP_TOKEN_EOF: c_int = 5; /* End Of File pseudo token */

const NGP_TTYPE_UNKNOWN: c_int = 0; /* undef (yet) token type - invalid to print/write to disk */
const NGP_TTYPE_BOOL: c_int = 1; /* boolean, it is 'T' or 'F' */
const NGP_TTYPE_STRING: c_int = 2; /* something withing "" or starting with letter */
const NGP_TTYPE_INT: c_int = 3; /* starting with digit and not with '.' */
const NGP_TTYPE_REAL: c_int = 4; /* digits + '.' */
const NGP_TTYPE_COMPLEX: c_int = 5; /* 2 reals, separated with ',' */
const NGP_TTYPE_NULL: c_int = 6; /* NULL token, format is : NAME = / comment */
const NGP_TTYPE_RAW: c_int = 7; /* HISTORY/COMMENT/8SPACES + comment string without / */

const NGP_FOUND_EQUAL_SIGN: c_int = 1; /* line contains '=' after keyword name */

const NGP_FORMAT_OK: c_int = 0; /* line format OK */
const NGP_FORMAT_ERROR: c_int = 1; /* line format error */

const NGP_NODE_INVALID: c_int = 0; /* default node type - invalid (to catch errors) */
const NGP_NODE_IMAGE: c_int = 1; /* IMAGE type */
const NGP_NODE_ATABLE: c_int = 2; /* ASCII table type */
const NGP_NODE_BTABLE: c_int = 3; /* BINARY table type */

const NGP_NON_SYSTEM_ONLY: c_int = 0; /* save all keywords except NAXIS,BITPIX,etc.. */
const NGP_REALLY_ALL: c_int = 1; /* save really all keywords */

const NGP_XTENSION_SIMPLE: c_int = 1; /* HDU defined with SIMPLE T */
const NGP_XTENSION_FIRST: c_int = 2; /* this is first extension in template */

const NGP_LINE_REREAD: c_int = 1; /* reread line */

const NGP_BITPIX_INVALID: c_int = -12345; /* default BITPIX (to catch errors) */

unsafe fn ngp_alloc(x: size_t) -> *mut c_void {
    unsafe { malloc(x) }
}
unsafe fn ngp_free(x: *mut c_void) {
    unsafe { free(x) }
}
unsafe fn ngp_realloc(x: *mut c_void, y: size_t) -> *mut c_void {
    unsafe { realloc(x, y) }
}

/// Safe allocation function that returns a boxed slice
/// Uses try_reserve_exact for proper error handling
fn ngp_alloc_boxed(size: usize) -> Result<Box<[c_char]>, c_int> {
    let mut vec = Vec::new();
    if vec.try_reserve_exact(size).is_err() {
        return Err(NGP_NO_MEMORY);
    }
    vec.resize(size, 0);
    Ok(vec.into_boxed_slice())
}

/* type definitions */

#[repr(C)]
#[derive(Clone, Default)]
pub struct NgpRawLine {
    pub line: Box<[c_char]>,
    pub name_idx: Option<usize>,
    pub value_idx: Option<usize>,
    pub type_: c_int,
    pub comment_idx: Option<usize>,
    pub format: c_int,
    pub flags: c_int,
}

#[repr(C)]
#[derive(Copy, Clone)]
pub struct NgpComplex {
    pub re: f64,
    pub im: f64,
}

#[repr(C)]
#[derive(Copy, Clone)]
pub union NgpTokval {
    pub s: *mut c_char,
    pub b: c_char,
    pub i: c_int,
    pub l: c_long,
    pub d: f64,
    pub c: NgpComplex,
}

#[repr(C)]
#[derive(Copy, Clone)]
pub struct NgpToken {
    pub type_: c_int,
    pub name: [c_char; NGP_MAX_NAME],
    pub value: NgpTokval,
    pub comment: [c_char; NGP_MAX_COMMENT],
}

impl Default for NgpToken {
    fn default() -> Self {
        NgpToken {
            type_: NGP_TTYPE_UNKNOWN,
            name: [0; NGP_MAX_NAME],
            value: NgpTokval { i: 0 },
            comment: [0; NGP_MAX_COMMENT],
        }
    }
}

#[repr(C)]
pub struct NgpHdu {
    pub tokcnt: c_int,
    pub tok: Vec<NgpToken>,
}

#[repr(C)]
pub struct NgpTkdef {
    pub name: &'static CStr,
    pub code: c_int,
}

#[repr(C)]
pub struct NgpExtverTab {
    pub extname: Box<[c_char]>,
    pub version: c_int,
}

/// Token definitions - constant lookup table for recognized tokens
pub const NGP_TKDEF: [NgpTkdef; 6] = [
    NgpTkdef {
        name: c"\\INCLUDE",
        code: NGP_TOKEN_INCLUDE,
    },
    NgpTkdef {
        name: c"\\GROUP",
        code: NGP_TOKEN_GROUP,
    },
    NgpTkdef {
        name: c"\\END",
        code: NGP_TOKEN_END,
    },
    NgpTkdef {
        name: c"XTENSION",
        code: NGP_TOKEN_XTENSION,
    },
    NgpTkdef {
        name: c"SIMPLE",
        code: NGP_TOKEN_SIMPLE,
    },
    NgpTkdef {
        name: c"",
        code: NGP_TOKEN_UNKNOWN,
    },
];

struct GRParseState {
    /// Current line being parsed
    NGP_CURLINE: NgpRawLine,

    /// Previous line that was parsed
    NGP_PREVLINE: NgpRawLine,

    /// Extension version table
    NGP_EXTVER_TAB: Vec<NgpExtverTab>,

    // More global variables from grparser.c
    /// Include file nesting level (1 means main file)
    NGP_INCLEVEL: c_int,

    /// Group nesting level (0 means no grouping)
    NGP_GRPLEVEL: c_int,

    /// Stack of included file handles
    NGP_FP: [Option<BufReader<File>>; NGP_MAX_INCLUDE + 1],

    /// Index of token in current line
    NGP_KEYIDX: c_int,

    /// Keyword after line analyze
    NGP_LINKEY: NgpToken,

    /// Directory of top level include file
    NGP_MASTER_DIR: PathBuf,

    /// Current unnamed group in object
    MASTER_GRP_IDX: c_int,
}

impl Default for GRParseState {
    fn default() -> Self {
        Self {
            NGP_CURLINE: Default::default(),
            NGP_PREVLINE: Default::default(),
            NGP_EXTVER_TAB: Vec::new(),
            NGP_INCLEVEL: Default::default(),
            NGP_GRPLEVEL: Default::default(),
            NGP_FP: Default::default(),
            NGP_KEYIDX: NGP_TOKEN_UNKNOWN,
            NGP_LINKEY: Default::default(),
            NGP_MASTER_DIR: PathBuf::new(),
            MASTER_GRP_IDX: 1,
        }
    }
}

fn ngp_get_extver(
    parser_state: &mut GRParseState,
    extname: &[c_char],
    version: &mut c_int,
) -> c_int {
    if extname.is_empty() {
        return NGP_BAD_ARG;
    }

    // Search for existing entry
    for entry in &mut parser_state.NGP_EXTVER_TAB {
        if 0 == strcmp_safe(extname, &entry.extname) {
            entry.version += 1;
            *version = entry.version;
            return NGP_OK;
        }
    }

    // Not found, create new entry
    // Calculate length including null terminator
    let len = strlen_safe(extname) + 1;

    // Allocate Vec with error checking
    let mut vec = Vec::new();
    if vec.try_reserve_exact(len).is_err() {
        return NGP_NO_MEMORY;
    }

    // Copy the C string into the Vec
    vec.extend_from_slice(&extname[..len]);

    let new_entry = NgpExtverTab {
        extname: vec.into_boxed_slice(),
        version: 1,
    };

    // Check if we can allocate space for the new entry
    if parser_state.NGP_EXTVER_TAB.try_reserve_exact(1).is_err() {
        return NGP_NO_MEMORY;
    }

    parser_state.NGP_EXTVER_TAB.push(new_entry);
    *version = 1;

    NGP_OK
}

fn ngp_set_extver(parser_state: &mut GRParseState, extname: &[c_char], version: c_int) -> c_int {
    if extname.is_empty() {
        return NGP_BAD_ARG;
    }

    // Search for existing entry
    for entry in &mut parser_state.NGP_EXTVER_TAB {
        if 0 == strcmp_safe(extname, &entry.extname) {
            if version > entry.version {
                entry.version = version;
            }
            return NGP_OK;
        }
    }

    // Not found, create new entry
    // Calculate length including null terminator
    let len = strlen_safe(extname) + 1;

    // Allocate Vec with error checking
    let mut vec = Vec::new();
    if vec.try_reserve_exact(len).is_err() {
        return NGP_NO_MEMORY;
    }

    // Copy the C string into the Vec
    vec.extend_from_slice(&extname[..len]);

    let new_entry = NgpExtverTab {
        extname: vec.into_boxed_slice(),
        version,
    };

    // Check if we can allocate space for the new entry
    if parser_state.NGP_EXTVER_TAB.try_reserve_exact(1).is_err() {
        return NGP_NO_MEMORY;
    }

    parser_state.NGP_EXTVER_TAB.push(new_entry);

    NGP_OK
}

fn ngp_delete_extver_tab(parser_state: &mut GRParseState) -> c_int {
    // Simply clear the vector - Box will handle cleanup automatically
    parser_state.NGP_EXTVER_TAB.clear();
    NGP_OK
}

/// read one line from file
unsafe fn ngp_line_from_file(reader: &mut BufReader<File>, line_out: &mut Box<[c_char]>) -> c_int {
    let mut r = NGP_OK; /* initialize stuff, reset err code */
    let mut llen: usize = 0; /* 0 characters read so far */
    let mut allocsize: usize = 1; /* preallocate 1 byte */

    let mut vec = Vec::new();
    if vec.try_reserve_exact(1).is_err() {
        return NGP_NO_MEMORY;
    }
    vec.resize(1, 0);

    let mut byte_buf = [0u8; 1];
    loop {
        // Read one byte at a time
        let bytes_read = match reader.read(&mut byte_buf) {
            Ok(0) => {
                // EOF
                if llen == 0 {
                    return NGP_EOF; /* signal EOF only if 0 characters read so far */
                }
                break;
            }
            Ok(n) => n,
            Err(_) => {
                r = NGP_READ_ERR;
                if llen == 0 {
                    return NGP_EOF;
                }
                break;
            }
        };

        let c = byte_buf[0] as c_char;

        if c == bb(b'\r') as c_char {
            continue; /* carriage return character ?  Just ignore it */
        }
        if c == bb(b'\n') as c_char {
            break; /* end of line character ? */
        }

        llen += 1; /* we have new character, make room for it */
        let alen =
            ((llen + NGP_ALLOCCHUNK as usize) / NGP_ALLOCCHUNK as usize) * NGP_ALLOCCHUNK as usize;
        if alen > allocsize {
            /* realloc buffer if needed - grow in chunks */
            let additional = alen - vec.len();
            if vec.try_reserve_exact(additional).is_err() {
                r = NGP_NO_MEMORY;
                break;
            }
            vec.resize(alen, 0);
            allocsize = alen;
        }
        vec[llen - 1] = c; /* copy character to buffer */
    }

    if (NGP_EOF != r) && (NGP_OK != r) {
        /* in case of errors, return empty box */
        *line_out = Box::default();
        return r;
    }

    llen += 1; /* place for terminating \0 */
    if llen != allocsize {
        /* shrink to exact size */
        vec.truncate(llen);
        vec.shrink_to_fit();
    }
    if llen > 0 {
        vec[llen - 1] = 0; /* add terminating \0 */
    }

    *line_out = vec.into_boxed_slice();

    r /* return  status code */
}

/// free current line structure
fn ngp_free_line(parser_state: &mut GRParseState) -> c_int {
    if !parser_state.NGP_CURLINE.line.is_empty() {
        parser_state.NGP_CURLINE.line = Box::default();
        parser_state.NGP_CURLINE.name_idx = None;
        parser_state.NGP_CURLINE.value_idx = None;
        parser_state.NGP_CURLINE.comment_idx = None;
        parser_state.NGP_CURLINE.type_ = NGP_TTYPE_UNKNOWN;
        parser_state.NGP_CURLINE.format = NGP_FORMAT_OK;
        parser_state.NGP_CURLINE.flags = 0;
    }
    NGP_OK
}

/// free cached line structure
fn ngp_free_prevline(parser_state: &mut GRParseState) -> c_int {
    if !parser_state.NGP_PREVLINE.line.is_empty() {
        parser_state.NGP_PREVLINE.line = Box::default();
        parser_state.NGP_PREVLINE.name_idx = None;
        parser_state.NGP_PREVLINE.value_idx = None;
        parser_state.NGP_PREVLINE.comment_idx = None;
        parser_state.NGP_PREVLINE.type_ = NGP_TTYPE_UNKNOWN;
        parser_state.NGP_PREVLINE.format = NGP_FORMAT_OK;
        parser_state.NGP_PREVLINE.flags = 0;
    }
    NGP_OK
}

/// read one line
unsafe fn ngp_read_line_buffered(
    parser_state: &mut GRParseState,
    reader: &mut BufReader<File>,
) -> c_int {
    ngp_free_line(parser_state); /* first free current line (if any) */

    if !parser_state.NGP_PREVLINE.line.is_empty()
    /* if cached, return cached line */
    {
        parser_state.NGP_CURLINE = parser_state.NGP_PREVLINE.clone();
        parser_state.NGP_PREVLINE = NgpRawLine::default();
        parser_state.NGP_CURLINE.flags = NGP_LINE_REREAD;
        return NGP_OK;
    }

    parser_state.NGP_CURLINE.flags = 0; /* if not cached really read line from file */
    unsafe { ngp_line_from_file(reader, &mut parser_state.NGP_CURLINE.line) }
}

/// unread line
fn ngp_unread_line(parser_state: &mut GRParseState) -> c_int {
    if parser_state.NGP_CURLINE.line.is_empty() {
        /* nothing to unread */
        return NGP_EMPTY_CURLINE;
    }

    if !parser_state.NGP_PREVLINE.line.is_empty() {
        /* we cannot unread line twice */
        return NGP_UNREAD_QUEUE_FULL;
    }

    parser_state.NGP_PREVLINE = parser_state.NGP_CURLINE.clone();
    parser_state.NGP_CURLINE = NgpRawLine::default();
    NGP_OK
}

/// a first guess line decomposition
fn ngp_extract_tokens(cl: &mut NgpRawLine) -> c_int {
    let mut p_idx: usize; // Index into line buffer
    let mut s_idx: usize;
    let mut cl_flags: c_int;
    let mut i: c_int;

    if cl.line.is_empty() {
        return NGP_NUL_PTR;
    }

    cl.comment_idx = None;
    cl.name_idx = None;
    cl.value_idx = None;
    cl.type_ = NGP_TTYPE_UNKNOWN;
    cl.format = NGP_FORMAT_OK;

    cl_flags = 0;
    p_idx = 0; /* start from beginning of line */

    i = 0;
    loop {
        /* if 8 spaces at beginning then line is comment */

        if (0 == cl.line[p_idx]) || (bb(b'\n') == cl.line[p_idx]) {
            /* if line has only blanks -> write blank keyword */
            cl.line[0] = 0; /* create empty name (0 length string) */
            cl.name_idx = Some(0);
            cl.comment_idx = Some(0);
            cl.type_ = NGP_TTYPE_RAW; /* signal write unformatted to FITS file */
            return NGP_OK;
        }
        if (bb(b' ') != cl.line[p_idx]) && (bb(b'\t') != cl.line[p_idx]) {
            break;
        }
        if i >= 7 {
            cl.comment_idx = Some(p_idx + 1);
            s_idx = cl.comment_idx.unwrap();
            loop {
                /* filter out any EOS characters in comment */
                if bb(b'\n') == cl.line[s_idx] {
                    cl.line[s_idx] = 0;
                }
                if 0 == cl.line[s_idx] {
                    break;
                }
                s_idx += 1;
            }
            cl.line[0] = 0; /* create empty name (0 length string) */
            cl.name_idx = Some(0);
            cl.type_ = NGP_TTYPE_RAW;
            return NGP_OK;
        }
        p_idx += 1;
        i += 1;
    }

    cl.name_idx = Some(p_idx);

    loop
    /* we need to find 1st whitespace */
    {
        if (0 == cl.line[p_idx]) || (bb(b'\n') == cl.line[p_idx]) {
            cl.line[p_idx] = 0;
            break;
        }

        /*
          from Richard Mathar, 2002-05-03, add 10 lines:
          if upper/lowercase HIERARCH followed also by an equal sign...
        */
        if p_idx + 9 <= cl.line.len()
            && fits_strncasecmp(cs!(c"HIERARCH"), &cl.line[p_idx..p_idx + 9], 8) == 0
        {
            // Find '=' sign starting from current position
            let mut found_eq = None;
            for idx in p_idx..cl.line.len() {
                if cl.line[idx] == bb(b'=') {
                    found_eq = Some(idx);
                    break;
                }
                if cl.line[idx] == 0 {
                    break;
                }
            }
            if let Some(eq_idx) = found_eq {
                cl_flags |= NGP_FOUND_EQUAL_SIGN;
                p_idx = eq_idx;
                break;
            }
        }

        if (bb(b' ') == cl.line[p_idx]) || (bb(b'\t') == cl.line[p_idx]) {
            break;
        }
        if bb(b'=') == cl.line[p_idx] {
            cl_flags |= NGP_FOUND_EQUAL_SIGN;
            break;
        }

        p_idx += 1;
    }

    if 0 != cl.line[p_idx] {
        cl.line[p_idx] = 0; /* found end of keyname so terminate string with zero */
        p_idx += 1;
    }

    // Get name as a slice for comparison
    let name_idx = cl.name_idx.unwrap();
    let name_slice = &cl.line[name_idx..];

    if (0
        == fits_strcasecmp(
            cast_slice::<u8, c_char>(c"HISTORY".to_bytes_with_nul()),
            name_slice,
        ))
        || (0
            == fits_strcasecmp(
                cast_slice::<u8, c_char>(c"COMMENT".to_bytes_with_nul()),
                name_slice,
            ))
        || (0
            == fits_strcasecmp(
                cast_slice::<u8, c_char>(c"CONTINUE".to_bytes_with_nul()),
                name_slice,
            ))
    {
        cl.comment_idx = Some(p_idx);
        s_idx = cl.comment_idx.unwrap();
        loop {
            /* filter out any EOS characters in comment */

            if bb(b'\n') == cl.line[s_idx] {
                cl.line[s_idx] = 0;
            }
            if 0 == cl.line[s_idx] {
                break;
            }
            s_idx += 1;
        }
        cl.type_ = NGP_TTYPE_RAW;
        return NGP_OK;
    }

    if 0 == fits_strcasecmp(
        cast_slice::<u8, c_char>(c"\\INCLUDE".to_bytes_with_nul()),
        name_slice,
    ) {
        loop {
            if (bb(b' ') != cl.line[p_idx]) && (bb(b'\t') != cl.line[p_idx]) {
                break; /* skip whitespace */
            }
            p_idx += 1;
        }

        cl.value_idx = Some(p_idx);
        s_idx = cl.value_idx.unwrap();
        loop {
            /* filter out any EOS characters */
            if bb(b'\n') == cl.line[s_idx] {
                cl.line[s_idx] = 0;
            }
            if 0 == cl.line[s_idx] {
                break;
            }
            s_idx += 1;
        }
        cl.type_ = NGP_TTYPE_UNKNOWN;
        return NGP_OK;
    }

    loop {
        if (0 == cl.line[p_idx]) || (bb(b'\n') == cl.line[p_idx]) {
            return NGP_OK; /* test if at end of string */
        }
        if (bb(b' ') == cl.line[p_idx]) || (bb(b'\t') == cl.line[p_idx]) {
            p_idx += 1;
            continue; /* skip whitespace */
        }
        if 0 != (cl_flags & NGP_FOUND_EQUAL_SIGN) {
            break;
        }
        if bb(b'=') != cl.line[p_idx] {
            break; /* ignore initial equal sign */
        }
        cl_flags |= NGP_FOUND_EQUAL_SIGN;
        p_idx += 1;
    }

    if bb(b'/') == cl.line[p_idx]
    /* no value specified, comment only */
    {
        p_idx += 1;
        if (bb(b' ') == cl.line[p_idx]) || (bb(b'\t') == cl.line[p_idx]) {
            p_idx += 1;
        }
        cl.comment_idx = Some(p_idx);
        s_idx = cl.comment_idx.unwrap();
        loop {
            /* filter out any EOS characters in comment */
            if bb(b'\n') == cl.line[s_idx] {
                cl.line[s_idx] = 0;
            }
            if 0 == cl.line[s_idx] {
                break;
            }
            s_idx += 1;
        }
        return NGP_OK;
    }

    if bb(b'\'') == cl.line[p_idx]
    /* we have found string within quotes */
    {
        p_idx += 1;
        s_idx = p_idx;
        cl.value_idx = Some(s_idx); /* set pointer to beginning of that string */
        cl.type_ = NGP_TTYPE_STRING; /* signal that it is of string type */

        loop {
            /* analyze it */

            if (0 == cl.line[p_idx]) || (bb(b'\n') == cl.line[p_idx])
            /* end of line -> end of string */
            {
                cl.line[s_idx] = 0;
                return NGP_OK;
            }

            if bb(b'\'') == cl.line[p_idx]
            /* we have found doublequote */
            {
                if (0 == cl.line[p_idx + 1]) || (bb(b'\n') == cl.line[p_idx + 1])
                /* doublequote is the last character in line */
                {
                    cl.line[s_idx] = 0;
                    return NGP_OK;
                }
                if (bb(b'\t') == cl.line[p_idx + 1]) || (bb(b' ') == cl.line[p_idx + 1])
                /* duoblequote was string terminator */
                {
                    cl.line[s_idx] = 0;
                    p_idx += 1;
                    break;
                }
                if bb(b'\'') == cl.line[p_idx + 1] {
                    p_idx += 1; /* doublequote is inside string, convert "" -> " */
                }
            }

            cl.line[s_idx] = cl.line[p_idx]; /* compact string in place, necess. by "" -> " conversion */
            s_idx += 1;
            p_idx += 1;
        }
    } else
    /* regular token */
    {
        cl.value_idx = Some(p_idx); /* set pointer to token */
        cl.type_ = NGP_TTYPE_UNKNOWN; /* we dont know type at the moment */

        loop {
            /* we need to find 1st whitespace */

            if (0 == cl.line[p_idx]) || (bb(b'\n') == cl.line[p_idx]) {
                cl.line[p_idx] = 0;
                return NGP_OK;
            }
            if (bb(b' ') == cl.line[p_idx]) || (bb(b'\t') == cl.line[p_idx]) {
                break;
            }
            p_idx += 1;
        }
        if cl.line[p_idx] != 0 {
            cl.line[p_idx] = 0; /* found so terminate string with zero */
            p_idx += 1;
        }
    }

    loop {
        if (0 == cl.line[p_idx]) || (bb(b'\n') == cl.line[p_idx]) {
            return NGP_OK; /* test if at end of string */
        }
        if (bb(b' ') != cl.line[p_idx]) && (bb(b'\t') != cl.line[p_idx]) {
            break; /* skip whitespace */
        }
        p_idx += 1;
    }

    if bb(b'/') == cl.line[p_idx]
    /* no value specified, comment only */
    {
        p_idx += 1;
        if (bb(b' ') == cl.line[p_idx]) || (bb(b'\t') == cl.line[p_idx]) {
            p_idx += 1;
        }
        cl.comment_idx = Some(p_idx);
        s_idx = cl.comment_idx.unwrap();
        loop {
            /* filter out any EOS characters in comment */
            if bb(b'\n') == cl.line[s_idx] {
                cl.line[s_idx] = 0;
            }
            if 0 == cl.line[s_idx] {
                break;
            }
            s_idx += 1;
        }
        return NGP_OK;
    }

    cl.format = NGP_FORMAT_ERROR;
    NGP_OK /* too many tokens ... */
}

/// Helper function to search for an include file in multiple locations.
/// Returns the full path to the file if found.
///
/// Search order:
/// 1. Direct path (as provided)
/// 2. Directories in CFITSIO_INCLUDE_FILES environment variable
/// 3. Relative to master_dir
unsafe fn find_include_file(fname: &[c_char], master_dir: &Path) -> Option<PathBuf> {
    unsafe {
        // Convert C string to Rust string
        let fname_str = match CStr::from_ptr(fname.as_ptr()).to_str() {
            Ok(s) => s,
            Err(_) => return None,
        };

        // 1. Try direct path first
        let direct_path = Path::new(fname_str);
        if direct_path.exists() {
            return Some(direct_path.to_path_buf());
        }

        // 2. Try directories from CFITSIO_INCLUDE_FILES environment variable
        if let Ok(env_paths) = env::var("CFITSIO_INCLUDE_FILES") {
            // Split on ':' for Unix or ';' for Windows
            let separator = if cfg!(target_os = "windows") {
                ';'
            } else {
                ':'
            };

            for dir_str in env_paths.split(separator) {
                let dir_path = Path::new(dir_str);
                let full_path = dir_path.join(fname_str);
                if full_path.exists() {
                    return Some(full_path);
                }
            }
        }

        // 3. Try relative to master directory
        // Only if fname is not an absolute path and master_dir is set
        if !direct_path.is_absolute() && master_dir != Path::new("") {
            let master_path = master_dir.join(fname_str);
            if master_path.exists() {
                return Some(master_path);
            }
        }

        None
    }
}

/// try to open include file.
/// If open fails and fname
/// does not specify absolute pathname, try to open fname
/// in any directory specified in CFITSIO_INCLUDE_FILES
/// environment variable. Finally try to open fname
/// relative to NGP_MASTER_DIR, which is directory of top
/// level include file
unsafe fn ngp_include_file(parser_state: &mut GRParseState, fname: &[c_char]) -> c_int /* try to open include file */
{
    unsafe {
        // Check nesting limit
        if parser_state.NGP_INCLEVEL >= NGP_MAX_INCLUDE as c_int {
            return NGP_INC_NESTING;
        }

        // Use the helper to find the file
        let file_path = match find_include_file(fname, &parser_state.NGP_MASTER_DIR) {
            Some(path) => path,
            None => return NGP_ERR_FOPEN,
        };

        // Try to open the file
        let file = match File::open(&file_path) {
            Ok(f) => f,
            Err(_) => return NGP_ERR_FOPEN,
        };

        // Create buffered reader and store it
        let reader = BufReader::new(file);
        parser_state.NGP_FP[parser_state.NGP_INCLEVEL as usize] = Some(reader);

        parser_state.NGP_INCLEVEL += 1;
        NGP_OK
    }
}

/// read line in the intelligent way.
/// All \INCLUDE directives are handled,
/// empty and comment line skipped. If this function returns NGP_OK, than
/// decomposed line (name, type, value in proper type and comment) are
/// stored in NGP_LINKEY structure. ignore_blank_lines parameter is zero
/// when parser is inside GROUP or HDU definition. Nonzero otherwise.
unsafe fn ngp_read_line(parser_state: &mut GRParseState, ignore_blank_lines: c_int) -> c_int {
    unsafe {
        let mut r: c_int;
        let mut nc: c_int = 0;
        let savec: c_int;
        let mut k: c_uint;

        if parser_state.NGP_INCLEVEL <= 0
        /* do some sanity checking first */
        {
            parser_state.NGP_KEYIDX = NGP_TOKEN_EOF; /* no parents, so report error */
            return NGP_OK;
        }
        if parser_state.NGP_INCLEVEL > NGP_MAX_INCLUDE as c_int {
            return NGP_INC_NESTING;
        }
        if parser_state.NGP_FP[(parser_state.NGP_INCLEVEL - 1) as usize].is_none() {
            return NGP_NUL_PTR;
        }

        loop {
            // Temporarily take the reader out of the array to avoid borrowing issues
            let idx = (parser_state.NGP_INCLEVEL - 1) as usize;
            let mut reader = parser_state.NGP_FP[idx].take().unwrap();

            r = ngp_read_line_buffered(parser_state, &mut reader);

            // Put the reader back
            parser_state.NGP_FP[idx] = Some(reader);

            match r {
                NGP_EOF => {
                    parser_state.NGP_INCLEVEL -= 1; /* end of file, revert to parent */

                    // Close file by dropping the BufReader (RAII - automatic cleanup)
                    parser_state.NGP_FP[parser_state.NGP_INCLEVEL as usize] = None;

                    if parser_state.NGP_INCLEVEL <= 0 {
                        parser_state.NGP_KEYIDX = NGP_TOKEN_EOF; /* no parents, so report error */
                        return NGP_OK;
                    }
                    continue;
                }
                NGP_OK => {
                    if (parser_state.NGP_CURLINE.flags & NGP_LINE_REREAD) != 0 {
                        return r;
                    }
                }
                _ => {
                    return r;
                }
            }

            match parser_state.NGP_CURLINE.line.get(0).copied().unwrap_or(0) as u8 {
                0 => {
                    if 0 == ignore_blank_lines {
                        // break;
                        /* ignore empty lines if told so */
                    }
                }
                b'#' => {
                    continue; /* ignore comment lines */
                }
                _ => {}
            }

            r = ngp_extract_tokens(&mut (parser_state.NGP_CURLINE)); /* analyse line, extract tokens and comment */
            if NGP_OK != r {
                return r;
            }

            if parser_state.NGP_CURLINE.name_idx.is_none() {
                continue; /* skip lines consisting only of whitespaces */
            }

            let name_idx = parser_state.NGP_CURLINE.name_idx.unwrap();
            for k in 0..8 {
                if name_idx + k >= parser_state.NGP_CURLINE.line.len() {
                    break;
                }
                if parser_state.NGP_CURLINE.line[name_idx + k] == 0 {
                    break;
                }
                if (parser_state.NGP_CURLINE.line[name_idx + k] >= bb(b'a'))
                    && (parser_state.NGP_CURLINE.line[name_idx + k] <= bb(b'z'))
                {
                    parser_state.NGP_CURLINE.line[name_idx + k] += bb(b'A') - bb(b'a'); /* force keyword to be upper case */
                }
            }

            k = 0;
            loop {
                /* find index of keyword in keyword table */
                if NGP_TOKEN_UNKNOWN == NGP_TKDEF[k as usize].code {
                    break;
                }
                let name_slice = &parser_state.NGP_CURLINE.line[name_idx..];
                if 0 == strcmp_safe(name_slice, cast_slice(NGP_TKDEF[k as usize].name.to_bytes_with_nul())) {
                    break;
                }
                k += 1;
            }

            parser_state.NGP_KEYIDX = NGP_TKDEF[k as usize].code; /* save this index, grammar parser will need this */

            if NGP_TOKEN_INCLUDE == parser_state.NGP_KEYIDX
            /* if this is \INCLUDE keyword, try to include file */
            {
                // Copy the filename to avoid borrow checker issues
                let filename_vec: Vec<c_char> = parser_state
                    .NGP_CURLINE
                    .value_idx
                    .map(|idx| parser_state.NGP_CURLINE.line[idx..].to_vec())
                    .unwrap_or_else(|| vec![0]);
                r = ngp_include_file(parser_state, &filename_vec);
                if NGP_OK != (r) {
                    return r;
                }
                continue; /* and read next line */
            }

            parser_state.NGP_LINKEY.type_ = NGP_TTYPE_UNKNOWN; /* now, get the keyword type, it's a long story ... */

            if let Some(value_idx) = parser_state.NGP_CURLINE.value_idx
            /* if no value given signal it */
            {
                let value_slice = &parser_state.NGP_CURLINE.line[value_idx..];
                let value_ptr = value_slice.as_ptr() as *mut c_char;

                if NGP_TTYPE_STRING == parser_state.NGP_CURLINE.type_
                /* string type test */
                {
                    parser_state.NGP_LINKEY.type_ = NGP_TTYPE_STRING;
                    parser_state.NGP_LINKEY.value.s = value_ptr;
                }
                if NGP_TTYPE_UNKNOWN == parser_state.NGP_LINKEY.type_
                /* bool type test */
                    && ((0
                        == fits_strcasecmp(
                            cast_slice::<u8, c_char>(c"T".to_bytes_with_nul()),
                            value_slice,
                        ))
                        || (0
                            == fits_strcasecmp(
                                cast_slice::<u8, c_char>(c"F".to_bytes_with_nul()),
                                value_slice,
                            )))
                {
                    parser_state.NGP_LINKEY.type_ = NGP_TTYPE_BOOL;
                    parser_state.NGP_LINKEY.value.b = if fits_strcasecmp(
                        cast_slice::<u8, c_char>(c"T".to_bytes_with_nul()),
                        value_slice,
                    ) != 0
                    {
                        0
                    } else {
                        1
                    };
                }

                let cstr_tmp = CStr::from_ptr(value_ptr);
                let cstr_len = cstr_tmp.to_bytes_with_nul().len();
                let str_slice =
                    slice::from_raw_parts_mut(cstr_tmp.as_ptr() as *mut c_char, cstr_len);

                if NGP_TTYPE_UNKNOWN == parser_state.NGP_LINKEY.type_
                /* complex type test */
                    && 2
                        == sscanf_lg_lg_n(
                            str_slice,
                            cs!(c"(%lg,%lg)%n"),
                            &mut (parser_state.NGP_LINKEY.value.c.re),
                            &mut (parser_state.NGP_LINKEY.value.c.im),
                            &mut nc,
                        )
                        && ((bb(b' ') == parser_state.NGP_CURLINE.line[value_idx + nc as usize])
                            || (bb(b'\t') == parser_state.NGP_CURLINE.line[value_idx + nc as usize])
                            || (bb(b'\n') == parser_state.NGP_CURLINE.line[value_idx + nc as usize])
                            || (0 == parser_state.NGP_CURLINE.line[value_idx + nc as usize]))
                {
                    parser_state.NGP_LINKEY.type_ = NGP_TTYPE_COMPLEX;
                }

                let cstr_tmp = CStr::from_ptr(value_ptr);
                let cstr_len = cstr_tmp.to_bytes_with_nul().len();
                let str_slice =
                    slice::from_raw_parts_mut(cstr_tmp.as_ptr() as *mut c_char, cstr_len);

                // Check for decimal point in value
                let has_decimal = value_slice.iter().any(|&ch| ch == bb(b'.'));

                if NGP_TTYPE_UNKNOWN == parser_state.NGP_LINKEY.type_
                /* real type test */
                    && has_decimal
                        && (1
                            == sscanf_lg_n(
                                str_slice,
                                cs!(c"%lg%n"),
                                &mut (parser_state.NGP_LINKEY.value.d),
                                &mut nc,
                            ))
                {
                    if bb(b'D') == parser_state.NGP_CURLINE.line[value_idx + nc as usize] {
                        /* test if template used a 'D' rather than an 'E' as the exponent character (added by WDP in 12/2010) */
                        savec = nc;
                        parser_state.NGP_CURLINE.line[value_idx + nc as usize] = bb(b'E');

                        let cstr_tmp = CStr::from_ptr(value_ptr);
                        let cstr_len = cstr_tmp.to_bytes_with_nul().len();
                        let str_slice =
                            slice::from_raw_parts_mut(cstr_tmp.as_ptr() as *mut c_char, cstr_len);

                        sscanf_lg_n(
                            str_slice,
                            cs!(c"%lg%n"),
                            &mut (parser_state.NGP_LINKEY.value.d),
                            &mut nc,
                        );
                        if (bb(b' ') == parser_state.NGP_CURLINE.line[value_idx + nc as usize])
                            || (bb(b'\t') == parser_state.NGP_CURLINE.line[value_idx + nc as usize])
                            || (bb(b'\n') == parser_state.NGP_CURLINE.line[value_idx + nc as usize])
                            || (0 == parser_state.NGP_CURLINE.line[value_idx + nc as usize])
                        {
                            parser_state.NGP_LINKEY.type_ = NGP_TTYPE_REAL;
                        } else {
                            /* no, this is not a real value */
                            parser_state.NGP_CURLINE.line[value_idx + savec as usize] = bb(b'D'); /* restore the original D character */
                        }
                    } else if (bb(b' ') == parser_state.NGP_CURLINE.line[value_idx + nc as usize])
                        || (bb(b'\t') == parser_state.NGP_CURLINE.line[value_idx + nc as usize])
                        || (bb(b'\n') == parser_state.NGP_CURLINE.line[value_idx + nc as usize])
                        || (0 == parser_state.NGP_CURLINE.line[value_idx + nc as usize])
                    {
                        parser_state.NGP_LINKEY.type_ = NGP_TTYPE_REAL;
                    }
                }

                let cstr_tmp = CStr::from_ptr(value_ptr);
                let cstr_len = cstr_tmp.to_bytes_with_nul().len();
                let str_slice =
                    slice::from_raw_parts_mut(cstr_tmp.as_ptr() as *mut c_char, cstr_len);

                if NGP_TTYPE_UNKNOWN == parser_state.NGP_LINKEY.type_
                /* integer type test */
                    && 1
                        == sscanf_d_n(
                            str_slice,
                            cs!(c"%d%n"),
                            &mut parser_state.NGP_LINKEY.value.i,
                            &mut nc,
                        )
                        && ((bb(b' ') == parser_state.NGP_CURLINE.line[value_idx + nc as usize])
                            || (bb(b'\t') == parser_state.NGP_CURLINE.line[value_idx + nc as usize])
                            || (bb(b'\n') == parser_state.NGP_CURLINE.line[value_idx + nc as usize])
                            || (0 == parser_state.NGP_CURLINE.line[value_idx + nc as usize]))
                {
                    parser_state.NGP_LINKEY.type_ = NGP_TTYPE_INT;
                }
                if NGP_TTYPE_UNKNOWN == parser_state.NGP_LINKEY.type_
                /* force string type */
                {
                    parser_state.NGP_LINKEY.type_ = NGP_TTYPE_STRING;
                    parser_state.NGP_LINKEY.value.s = value_ptr;
                }
            } else if NGP_TTYPE_RAW == parser_state.NGP_CURLINE.type_ {
                parser_state.NGP_LINKEY.type_ = NGP_TTYPE_RAW;
            } else {
                parser_state.NGP_LINKEY.type_ = NGP_TTYPE_NULL;
            }

            if let Some(comment_idx) = parser_state.NGP_CURLINE.comment_idx {
                let comment_ptr = parser_state.NGP_CURLINE.line[comment_idx..].as_ptr();
                strncpy(
                    std::ptr::addr_of_mut!(parser_state.NGP_LINKEY.comment).cast(),
                    comment_ptr,
                    NGP_MAX_COMMENT,
                ); /* store comment */
                parser_state.NGP_LINKEY.comment[NGP_MAX_COMMENT - 1] = 0;
            } else {
                parser_state.NGP_LINKEY.comment[0] = 0;
            }

            let name_ptr = parser_state.NGP_CURLINE.line[name_idx..].as_ptr();
            strncpy(
                std::ptr::addr_of_mut!(parser_state.NGP_LINKEY.name).cast(),
                name_ptr,
                NGP_MAX_NAME,
            ); /* and keyword's name */
            parser_state.NGP_LINKEY.name[NGP_MAX_NAME - 1] = 0;

            if strlen(std::ptr::addr_of!(parser_state.NGP_LINKEY.name).cast()) > FLEN_KEYWORD
            /* WDP: 20-Jun-2002:  mod to support HIERARCH */
            {
                return NGP_BAD_ARG; /* cfitsio does not allow names > 8 chars */
            }

            return NGP_OK; /* we have valid non empty line, so return success */
        }
    }
}

/// check whether keyword can be written as is
unsafe fn ngp_keyword_is_write(ngp_tok: *const NgpToken) -> c_int {
    unsafe {
        let mut i: c_int;
        let mut j: c_int;
        let mut l: c_int;
        let mut spc: c_int;

        /* indexed variables not to write */
        let nm: [*const c_char; 4] = [
            c"NAXIS".as_ptr(),
            c"TFORM".as_ptr(),
            c"TTYPE".as_ptr(),
            ptr::null(),
        ];

        /* non indexed variables not allowed to write */
        let nmni: [*const c_char; 11] = [
            c"SIMPLE".as_ptr(),
            c"XTENSION".as_ptr(),
            c"BITPIX".as_ptr(),
            c"NAXIS".as_ptr(),
            c"PCOUNT".as_ptr(),
            c"GCOUNT".as_ptr(),
            c"TFIELDS".as_ptr(),
            c"THEAP".as_ptr(),
            c"EXTEND".as_ptr(),
            c"EXTVER".as_ptr(),
            ptr::null(),
        ];

        if ngp_tok.is_null() {
            return NGP_NUL_PTR;
        }

        j = 0;

        /* first check non indexed */
        loop {
            if nmni[j as usize].is_null() {
                break;
            }
            if 0 == strcmp(nmni[j as usize], (*ngp_tok).name.as_ptr()) {
                return NGP_BAD_ARG;
            }
            j += 1;
        }

        j = 0;
        /* now check indexed */
        loop {
            if nm[j as usize].is_null() {
                return NGP_OK;
            }
            l = strlen(nm[j as usize]) as c_int;
            if (l < 1) || (l > 5) {
                j += 1;
                continue;
            }
            if 0 == strncmp(nm[j as usize], (*ngp_tok).name.as_ptr(), l as usize) {
                break;
            }
            j += 1;
        }

        if ((*ngp_tok).name[l as usize] < bb(b'1')) || ((*ngp_tok).name[l as usize] > bb(b'9')) {
            return NGP_OK;
        }
        spc = 0;
        for i in (l + 1)..8 {
            if spc != 0 {
                if bb(b' ') != (*ngp_tok).name[i as usize] {
                    return NGP_OK;
                }
            } else {
                if ((*ngp_tok).name[i as usize] >= bb(b'0'))
                    && ((*ngp_tok).name[i as usize] <= bb(b'9'))
                {
                    continue;
                }
                if bb(b' ') == (*ngp_tok).name[i as usize] {
                    spc = 1;
                    continue;
                }
                if 0 == (*ngp_tok).name[i as usize] {
                    break;
                }
                return NGP_OK;
            }
        }
        NGP_BAD_ARG
    }
}

/// write (almost) all keywords from given HDU to disk
unsafe fn ngp_keyword_all_write(ngph: &mut NgpHdu, ffp: &mut fitsfile, mode: c_int) -> c_int {
    unsafe {
        let mut i: c_int;
        let mut r: c_int;
        let mut ib: c_int;
        let mut buf: [c_char; 200] = [0; 200];
        let mut l: c_long;

        r = NGP_OK;

        for i in 0..(*ngph).tokcnt {
            let token = &(&*ngph).tok[i as usize];
            r = ngp_keyword_is_write(token);
            if (NGP_REALLY_ALL & mode) != 0 || (NGP_OK == r) {
                match token.type_ {
                    NGP_TTYPE_BOOL => {
                        ib = c_int::from(token.value.b);
                        fits_write_key(
                            ffp,
                            TLOGICAL,
                            token.name.as_ptr(),
                            (&ib as *const c_int).cast::<c_void>(),
                            token.comment.as_ptr(),
                            &mut r,
                        );
                    }
                    NGP_TTYPE_STRING => {
                        fits_write_key_longstr(
                            ffp,
                            token.name.as_ptr(),
                            token.value.s,
                            token.comment.as_ptr(),
                            &mut r,
                        );
                    }
                    NGP_TTYPE_INT => {
                        l = c_long::from(token.value.i); /* bugfix - 22-Jan-99, BO - nonalignment of OSF/Alpha */
                        fits_write_key(
                            ffp,
                            TLONG,
                            token.name.as_ptr(),
                            (&l as *const c_long).cast::<c_void>(),
                            token.comment.as_ptr(),
                            &mut r,
                        );
                    }
                    NGP_TTYPE_REAL => {
                        fits_write_key(
                            ffp,
                            TDOUBLE,
                            token.name.as_ptr(),
                            (&token.value.d as *const f64).cast::<c_void>(),
                            token.comment.as_ptr(),
                            &mut r,
                        );
                    }
                    NGP_TTYPE_COMPLEX => {
                        fits_write_key(
                            ffp,
                            TDBLCOMPLEX,
                            token.name.as_ptr(),
                            (&token.value.c as *const NgpComplex).cast::<c_void>(),
                            token.comment.as_ptr(),
                            &mut r,
                        );
                    }
                    NGP_TTYPE_NULL => {
                        fits_write_key_null(
                            ffp,
                            token.name.as_ptr(),
                            token.comment.as_ptr(),
                            &mut r,
                        );
                    }
                    NGP_TTYPE_RAW => {
                        let mut skip = false;
                        if 0 == strcmp(c"HISTORY".as_ptr(), token.name.as_ptr()) {
                            fits_write_history(ffp, token.comment.as_ptr(), &mut r);
                            skip = true;
                        }
                        if !skip && (0 == strcmp(c"COMMENT".as_ptr(), token.name.as_ptr())) {
                            fits_write_comment(ffp, token.comment.as_ptr(), &mut r);
                            skip = true;
                        }
                        if !skip {
                            int_snprintf!(
                                &mut buf,
                                200,
                                "{:-8.8}{}",
                                CStr::from_bytes_with_nul(cast_slice(&token.name))
                                    .unwrap()
                                    .to_str()
                                    .unwrap(),
                                CStr::from_bytes_with_nul(cast_slice(&token.comment))
                                    .unwrap()
                                    .to_str()
                                    .unwrap(),
                            );
                            fits_write_record(ffp, buf.as_mut_ptr(), &mut r);
                        }
                    }
                    _ => {}
                }
            } else if NGP_BAD_ARG == r
            /* enhancement 10 dec 2003, James Peachey: template comments replace defaults */
            {
                r = NGP_OK; /* update comments of special keywords like TFORM */
                if 0 != token.comment[0]
                /* do not update with a blank comment */
                {
                    fits_modify_comment(ffp, token.name.as_ptr(), token.comment.as_ptr(), &mut r);
                }
            } else
            /* other problem, typically a blank token */
            {
                r = NGP_OK; /* skip this token, but continue */
            }
            if r != 0 {
                return r;
            }
        }

        fits_set_hdustruc(ffp, &mut r); /* resync cfitsio */
        r
    }
}

/// init HDU structure
fn ngp_hdu_init(ngph: &mut NgpHdu) -> c_int {
    ngph.tok = Vec::new();
    ngph.tokcnt = 0;
    NGP_OK
}

/// clear HDU structure
unsafe fn ngp_hdu_clear(ngph: &mut NgpHdu) -> c_int {
    unsafe {
        for i in 0..(*ngph).tokcnt {
            let token = &mut (&mut *ngph).tok[i as usize];
            if NGP_TTYPE_STRING == token.type_ && !token.value.s.is_null() {
                // Free by reconstructing the Box
                let str_len = 1 + strlen(token.value.s);
                drop(Box::from_raw(slice::from_raw_parts_mut(
                    token.value.s,
                    str_len,
                )));
                token.value.s = ptr::null_mut();
            }
        }

        (*ngph).tok.clear();
        (*ngph).tokcnt = 0;

        NGP_OK
    }
}

/// insert new token to HDU structure
fn ngp_hdu_insert_token(ngph: &mut NgpHdu, newtok: &mut NgpToken) -> c_int {
    unsafe {
        (&mut *ngph).tok.push(*newtok);

        if NGP_TTYPE_STRING == (*newtok).type_ && !(*newtok).value.s.is_null() {
            let last_idx = (&*ngph).tok.len() - 1;
            let str_len = 1 + strlen((*newtok).value.s);
            let str_boxed = match ngp_alloc_boxed(str_len) {
                Ok(b) => b,
                Err(e) => return e,
            };
            (&mut *ngph).tok[last_idx].value.s = Box::into_raw(str_boxed).cast::<c_char>();
            strcpy((&*ngph).tok[last_idx].value.s, (*newtok).value.s);
        }

        (*ngph).tokcnt += 1;
        NGP_OK
    }
}

fn ngp_append_columns(ff: &mut fitsfile, ngph: &mut NgpHdu, aftercol: c_int) -> c_int {
    unsafe {
        let mut r: c_int;
        let mut i: c_int;
        let j: c_int = 0;
        let mut exitflg: c_int;
        let mut ngph_i: c_int = 0;
        let mut my_tform: *mut c_char;
        let mut my_ttype: *mut c_char;
        let mut ngph_ctmp: c_char = 0;

        if 0 == (*ngph).tokcnt {
            return NGP_OK; /* nothing to do ! */
        }

        r = NGP_OK;
        exitflg = 0;

        for j in aftercol..NGP_MAX_ARRAY_DIM as c_int {
            /* 0 for table, 6 for group */

            my_tform = ptr::null_mut();
            my_ttype = c"".as_ptr() as *mut c_char;

            i = 0;
            loop {
                let token = &(&*ngph).tok[i as usize];
                if 1 == sscanf_d_c(
                    &mut token.name.clone(),
                    cs!(c"TFORM%d%c"),
                    &mut ngph_i,
                    &mut ngph_ctmp,
                ) {
                    if (NGP_TTYPE_STRING == token.type_) && (ngph_i == (j + 1)) {
                        my_tform = token.value.s;
                    }
                } else if 1
                    == sscanf_d_c(
                        &mut token.name.clone(),
                        cs!(c"TTYPE%d%c"),
                        &mut ngph_i,
                        &mut ngph_ctmp,
                    )
                    && (NGP_TTYPE_STRING == token.type_)
                    && (ngph_i == (j + 1))
                {
                    my_ttype = token.value.s;
                }

                if !my_tform.is_null() && (*my_ttype.offset(0)) != 0 {
                    break;
                }

                if i < ((*ngph).tokcnt - 1) {
                    i += 1;
                    continue;
                }
                exitflg = 1;
                break;
            }
            if (NGP_OK == r) && !my_tform.is_null() {
                fits_insert_col(ff, j + 1, my_ttype, my_tform, &mut r);
            }

            if (NGP_OK != r) || exitflg != 0 {
                break;
            }
        }
        r
    }
}

/// read complete HDU
unsafe fn ngp_read_xtension(
    parser_state: &mut GRParseState,
    ff: &mut fitsfile,
    parent_hn: c_int,
    simple_mode: c_int,
) -> c_int {
    unsafe {
        let mut r: c_int;
        let mut exflg: c_int;
        let mut l: c_int;
        let mut my_hn: c_int = 0;
        let mut tmp0: c_int = 0;
        let mut incrementor_index: c_int;
        let mut i: c_int;
        let mut j: c_int = 0;
        let mut ngph_dim: c_int;
        let mut ngph_bitpix: c_int;
        let mut ngph_node_type: c_int;
        let mut my_version: c_int = 0;
        let mut incrementor_name: [c_char; NGP_MAX_STRING] = [0; NGP_MAX_STRING];
        let mut ngph_ctmp: c_char = 0;
        let mut ngph_extname: *mut c_char = ptr::null_mut();
        let mut ngph_size: [c_long; NGP_MAX_ARRAY_DIM] = [0; NGP_MAX_ARRAY_DIM];
        let mut ngph: NgpHdu = NgpHdu {
            tokcnt: 0,
            tok: Vec::new(),
        };
        let lv: c_long;

        incrementor_name[0] = 0; /* signal no keyword+'#' found yet */
        incrementor_index = 0;

        r = ngp_hdu_init(&mut ngph);
        if NGP_OK != (r) {
            return r;
        }
        r = ngp_read_line(parser_state, 0);
        if NGP_OK != (r) {
            return r; /* EOF always means error here */
        }
        match NGP_XTENSION_SIMPLE & simple_mode {
            0 => {
                if NGP_TOKEN_XTENSION != parser_state.NGP_KEYIDX {
                    return NGP_TOKEN_NOT_EXPECT;
                }
            }
            _ => {
                if NGP_TOKEN_SIMPLE != parser_state.NGP_KEYIDX {
                    return NGP_TOKEN_NOT_EXPECT;
                }
            }
        }

        r = ngp_hdu_insert_token(&mut ngph, &mut (parser_state.NGP_LINKEY));
        if NGP_OK != (r) {
            return r;
        }

        loop {
            r = ngp_read_line(parser_state, 0);
            if NGP_OK != (r) {
                return r; /* EOF always means error here */
            }
            exflg = 0;
            match parser_state.NGP_KEYIDX {
                NGP_TOKEN_SIMPLE => {
                    r = NGP_TOKEN_NOT_EXPECT;
                }
                NGP_TOKEN_END | NGP_TOKEN_XTENSION | NGP_TOKEN_GROUP => {
                    r = ngp_unread_line(parser_state); /* WARNING - not break here .... */
                    exflg = 1;
                }
                NGP_TOKEN_EOF => {
                    exflg = 1;
                }
                _ => {
                    // default case - process normal keyword
                    l = strlen(std::ptr::addr_of!(parser_state.NGP_LINKEY.name).cast()) as c_int;
                    if (l >= 2)
                        && (l <= 6)
                        && bb(b'#') == parser_state.NGP_LINKEY.name[(l - 1) as usize]
                    {
                        if 0 == incrementor_name[0] {
                            memcpy(
                                incrementor_name.as_mut_ptr().cast::<c_void>(),
                                std::ptr::addr_of!(parser_state.NGP_LINKEY.name).cast::<c_void>(),
                                (l - 1) as size_t,
                            );
                            incrementor_name[(l - 1) as usize] = 0;
                        }
                        if ((l - 1) == strlen(incrementor_name.as_ptr()) as c_int)
                            && (0
                                == memcmp(
                                    incrementor_name.as_ptr().cast::<c_void>(),
                                    std::ptr::addr_of!(parser_state.NGP_LINKEY.name)
                                        .cast::<c_void>(),
                                    (l - 1) as size_t,
                                ))
                        {
                            incrementor_index += 1;
                        }
                        int_snprintf!(
                            &mut (parser_state.NGP_LINKEY.name[(l - 1) as usize..]),
                            ((NGP_MAX_NAME as c_int) - l + 1) as size_t,
                            "{}",
                            incrementor_index,
                        );
                    }
                    r = ngp_hdu_insert_token(&mut ngph, &mut (parser_state.NGP_LINKEY));
                    // NO break here - continue the loop to process next keyword
                }
            }
            if (NGP_OK != r) || exflg != 0 {
                break;
            }
        }

        if NGP_OK == r {
            /* we should scan keywords, and calculate HDU's */
            /* structure ourselves .... */

            ngph_node_type = NGP_NODE_INVALID; /* init variables */
            ngph_bitpix = 0;
            ngph_extname = ptr::null_mut();
            for i in 0..NGP_MAX_ARRAY_DIM {
                ngph_size[i] = 0;
            }
            ngph_dim = 0;
            for i in 0..ngph.tokcnt as usize {
                let token = &ngph.tok[i];
                if strcmp(c"XTENSION".as_ptr(), token.name.as_ptr()) == 0 {
                    if NGP_TTYPE_STRING == token.type_ {
                        if fits_strncasecmp(
                            cast_slice::<u8, c_char>(c"BINTABLE".to_bytes_with_nul()),
                            std::slice::from_raw_parts(token.value.s as *const c_char, 9),
                            8,
                        ) == 0
                        {
                            ngph_node_type = NGP_NODE_BTABLE;
                        }
                        if fits_strncasecmp(
                            cast_slice::<u8, c_char>(c"TABLE".to_bytes_with_nul()),
                            std::slice::from_raw_parts(token.value.s as *const c_char, 6),
                            5,
                        ) == 0
                        {
                            ngph_node_type = NGP_NODE_ATABLE;
                        }
                        if fits_strncasecmp(
                            cast_slice::<u8, c_char>(c"IMAGE".to_bytes_with_nul()),
                            std::slice::from_raw_parts(token.value.s as *const c_char, 6),
                            5,
                        ) == 0
                        {
                            ngph_node_type = NGP_NODE_IMAGE;
                        }
                    }
                } else if strcmp(c"SIMPLE".as_ptr(), token.name.as_ptr()) == 0 {
                    if NGP_TTYPE_BOOL == token.type_ && token.value.b != 0 {
                        ngph_node_type = NGP_NODE_IMAGE;
                    }
                } else if strcmp(c"BITPIX".as_ptr(), token.name.as_ptr()) == 0 {
                    if NGP_TTYPE_INT == token.type_ {
                        ngph_bitpix = token.value.i;
                    }
                } else if strcmp(c"NAXIS".as_ptr(), token.name.as_ptr()) == 0 {
                    if NGP_TTYPE_INT == token.type_ {
                        ngph_dim = token.value.i;
                    }
                } else if strcmp(c"EXTNAME".as_ptr(), token.name.as_ptr()) == 0
                /* assign EXTNAME, I hope struct does not move */
                {
                    if NGP_TTYPE_STRING == token.type_ {
                        ngph_extname = token.value.s;
                    }
                } else if 1
                    == sscanf_d_c(
                        &mut token.name.clone(),
                        cs!(c"NAXIS%d%c"),
                        &mut j,
                        &mut ngph_ctmp,
                    )
                    && NGP_TTYPE_INT == token.type_
                    && (j >= 1)
                    && (j <= NGP_MAX_ARRAY_DIM as c_int)
                {
                    ngph_size[(j - 1) as usize] = c_long::from(token.value.i);
                }
            }

            match ngph_node_type {
                NGP_NODE_IMAGE => {
                    if NGP_XTENSION_FIRST
                        == ((NGP_XTENSION_FIRST | NGP_XTENSION_SIMPLE) & simple_mode)
                    {
                        /* if caller signals that this is 1st HDU in file */
                        /* and it is IMAGE defined with XTENSION, then we */
                        /* need create dummy Primary HDU */
                        fits_create_img(ff, 16, 0, ptr::null_mut(), &mut r);
                    }
                    /* create image */
                    fits_create_img(ff, ngph_bitpix, ngph_dim, ngph_size.as_ptr(), &mut r);

                    /* update keywords */
                    if NGP_OK == r {
                        r = ngp_keyword_all_write(&mut ngph, ff, NGP_NON_SYSTEM_ONLY);
                    }
                }
                NGP_NODE_ATABLE | NGP_NODE_BTABLE => {
                    /* create table, 0 rows and 0 columns for the moment */
                    fits_create_tbl(
                        ff,
                        if NGP_NODE_ATABLE == ngph_node_type {
                            ASCII_TBL
                        } else {
                            BINARY_TBL
                        },
                        0,
                        0,
                        ptr::null_mut(),
                        ptr::null_mut(),
                        ptr::null_mut(),
                        ptr::null_mut(),
                        &mut r,
                    );

                    let continue_match = true;

                    while continue_match {
                        if NGP_OK != r {
                            break;
                        }

                        /* add columns ... */
                        r = ngp_append_columns(ff, &mut ngph, 0);
                        if NGP_OK != r {
                            break;
                        }

                        /* add remaining keywords */
                        r = ngp_keyword_all_write(&mut ngph, ff, NGP_NON_SYSTEM_ONLY);
                        if NGP_OK != r {
                            break;
                        }

                        /* if requested add rows */
                        if ngph_size[1] > 0 {
                            fits_insert_rows(ff, 0, ngph_size[1].into(), &mut r);
                        }
                        break;
                    }
                }
                _ => {
                    r = NGP_BAD_ARG;
                }
            }
        }

        if (NGP_OK == r) && !ngph_extname.is_null() {
            let extname_slice = slice::from_raw_parts(ngph_extname, strlen(ngph_extname) + 1);
            r = ngp_get_extver(parser_state, extname_slice, &mut my_version); /* write correct ext version number */
            lv = c_long::from(my_version); /* bugfix - 22-Jan-99, BO - nonalignment of OSF/Alpha */
            fits_write_key(
                ff,
                TLONG,
                c"EXTVER".as_ptr(),
                (&lv as *const c_long).cast::<c_void>(),
                c"auto assigned by template parser".as_ptr(),
                &mut r,
            );
        }

        if NGP_OK == r && parent_hn > 0 {
            fits_get_hdu_num(ff, &mut my_hn);
            fits_movabs_hdu(ff, parent_hn, &mut tmp0, &mut r); /* link us to parent */
            fits_add_group_member(ff, ptr::null_mut(), my_hn, &mut r);
            fits_movabs_hdu(ff, my_hn, &mut tmp0, &mut r);
            if NGP_OK != r {
                return r;
            }
        }

        if NGP_OK != r
        /* in case of error - delete hdu */
        {
            tmp0 = 0;
            fits_delete_hdu(ff, ptr::null_mut(), &mut tmp0);
        }

        ngp_hdu_clear(&mut ngph);
        r
    }
}

/// read complete GROUP
unsafe fn ngp_read_group(
    parser_state: &mut GRParseState,
    ff: &mut fitsfile,
    grpname: &[c_char],
    parent_hn: c_int,
) -> c_int {
    unsafe {
        let mut r: c_int = 0;
        let mut exitflg: c_int = 0;
        let mut l: c_int = 0;
        let mut my_hn: c_int = 0;
        let mut tmp0: c_int = 0;
        let mut incrementor_index: c_int = 0;
        let mut grnm: [c_char; NGP_MAX_STRING] = [0; NGP_MAX_STRING]; /* keyword holding group name */
        let mut incrementor_name: [c_char; NGP_MAX_STRING] = [0; NGP_MAX_STRING];
        let mut ngph: NgpHdu = NgpHdu {
            tokcnt: 0,
            tok: Vec::new(),
        };

        incrementor_name[0] = 0; /* signal no keyword+'#' found yet */
        incrementor_index = 6; /* first 6 cols are used by group */

        parser_state.NGP_GRPLEVEL += 1;
        r = ngp_hdu_init(&mut ngph);
        if NGP_OK != (r) {
            return r;
        }

        r = NGP_OK;
        r = fits_create_group(ff, grpname.as_ptr(), GT_ID_ALL_URI as c_int, &mut r);
        if NGP_OK != (r) {
            return r;
        }
        fits_get_hdu_num(ff, &mut my_hn);
        if parent_hn > 0 {
            fits_movabs_hdu(ff, parent_hn, &mut tmp0, &mut r); /* link us to parent */
            fits_add_group_member(ff, ptr::null_mut(), my_hn, &mut r);
            fits_movabs_hdu(ff, my_hn, &mut tmp0, &mut r);
            if NGP_OK != r {
                return r;
            }
        }

        exitflg = 0;
        while 0 == exitflg {
            r = ngp_read_line(parser_state, 0);
            if NGP_OK != (r) {
                break; /* EOF always means error here */
            }
            match parser_state.NGP_KEYIDX {
                NGP_TOKEN_SIMPLE | NGP_TOKEN_EOF => {
                    r = NGP_TOKEN_NOT_EXPECT;
                }
                NGP_TOKEN_END => {
                    parser_state.NGP_GRPLEVEL -= 1;
                    exitflg = 1;
                }

                NGP_TOKEN_GROUP => {
                    if NGP_TTYPE_STRING == parser_state.NGP_LINKEY.type_ {
                        strncpy(
                            grnm.as_mut_ptr(),
                            parser_state.NGP_LINKEY.value.s,
                            NGP_MAX_STRING,
                        );
                    } else {
                        int_snprintf!(
                            &mut grnm,
                            NGP_MAX_STRING as size_t,
                            "DEFAULT_GROUP_{}",
                            parser_state.MASTER_GRP_IDX,
                        );
                        parser_state.MASTER_GRP_IDX += 1;
                    }
                    grnm[NGP_MAX_STRING - 1] = 0;
                    r = ngp_read_group(parser_state, ff, &grnm, my_hn);
                    /* we can have many subsequent GROUP defs */
                }
                NGP_TOKEN_XTENSION => {
                    r = ngp_unread_line(parser_state);
                    if NGP_OK == r {
                        r = ngp_read_xtension(parser_state, ff, my_hn, 0);
                    }
                    /* we can have many subsequent HDU defs */
                }
                _ => {
                    l = strlen(std::ptr::addr_of!(parser_state.NGP_LINKEY.name).cast()) as c_int;
                    if (l >= 2)
                        && (l <= 6)
                        && bb(b'#') == parser_state.NGP_LINKEY.name[(l - 1) as usize]
                    {
                        if 0 == incrementor_name[0] {
                            memcpy(
                                incrementor_name.as_mut_ptr().cast::<c_void>(),
                                std::ptr::addr_of!(parser_state.NGP_LINKEY.name).cast::<c_void>(),
                                (l - 1) as size_t,
                            );
                            incrementor_name[(l - 1) as usize] = 0;
                        }
                        if ((l - 1) == strlen(incrementor_name.as_ptr()) as c_int)
                            && (0
                                == memcmp(
                                    incrementor_name.as_ptr().cast::<c_void>(),
                                    std::ptr::addr_of!(parser_state.NGP_LINKEY.name)
                                        .cast::<c_void>(),
                                    (l - 1) as size_t,
                                ))
                        {
                            incrementor_index += 1;
                        }
                        int_snprintf!(
                            parser_state.NGP_LINKEY.name[(l - 1) as usize..],
                            ((NGP_MAX_NAME as c_int) - l + 1) as usize,
                            "{}",
                            incrementor_index
                        );
                    }
                    r = ngp_hdu_insert_token(&mut ngph, &mut (parser_state.NGP_LINKEY));
                    /* here we can add keyword */
                }
            }
            if NGP_OK != r {
                break;
            }
        }

        fits_movabs_hdu(ff, my_hn, &mut tmp0, &mut r); /* back to our HDU */

        if NGP_OK == r {
            /* create additional columns, if requested */
            r = ngp_append_columns(ff, &mut ngph, 6);
        }

        if NGP_OK == r {
            /* and write keywords */
            r = ngp_keyword_all_write(&mut ngph, ff, NGP_NON_SYSTEM_ONLY);
        }

        if NGP_OK != r
        /* delete group in case of error */
        {
            tmp0 = 0;
            fits_remove_group(ff, OPT_RM_GPT as c_int, &mut tmp0);
        }

        ngp_hdu_clear(&mut ngph); /* we are done with this HDU, so delete it */
        r
    }
}

/* top level API functions */

/*--------------------------------------------------------------------------*/
/// Read whole template
/// ff should point to the opened empty fits file.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn fits_execute_template(
    ff: *mut fitsfile,           /* I - FITS file pointer */
    ngp_template: *const c_char, /* I - template string */
    status: *mut c_int,          /* IO - error status */
) -> c_int {
    unsafe {
        let mut r: c_int = 0;
        let mut exit_flg: c_int = 0;
        let mut first_extension: c_int = 0;
        let mut i: c_int = 0;
        let mut my_hn: c_int = 0;
        let mut tmp0: c_int = 0;
        let mut keys_exist: c_int = 0;
        let mut more_keys: c_int = 0;
        let mut used_ver: c_int = 0;
        let mut grnm: [c_char; NGP_MAX_STRING] = [0; NGP_MAX_STRING];
        let mut used_name: [c_char; NGP_MAX_STRING] = [0; NGP_MAX_STRING];
        let mut luv: c_long = 0;

        let mut parser_state = GRParseState::default();

        if status.is_null() {
            return NGP_NUL_PTR;
        }

        let status = status.as_mut().expect(NULL_MSG);

        if NGP_OK != *status {
            return *status;
        }

        /* This function uses many global variables (local to this file) and therefore is not thread-safe. */
        let _lock = FFLOCK();

        if ff.is_null() || ngp_template.is_null() {
            *status = NGP_NUL_PTR;
            return *status;
        }

        let ff = ff.as_mut().expect(NULL_MSG);

        parser_state.NGP_INCLEVEL = 0; /* initialize things, not all should be zero */
        parser_state.NGP_GRPLEVEL = 0;
        parser_state.MASTER_GRP_IDX = 1;
        exit_flg = 0;
        parser_state.NGP_MASTER_DIR = PathBuf::new(); /* this should be before 1st call to ngp_include_file */
        first_extension = 1; /* we need to create PHDU */

        r = ngp_delete_extver_tab(&mut parser_state);
        if NGP_OK != r {
            *status = r;
            return r;
        }

        fits_get_hdu_num(ff, &mut my_hn); /* our HDU position */
        if my_hn <= 1
        /* check whether we really need to create PHDU */
        {
            fits_movabs_hdu(ff, 1, &mut tmp0, status);
            fits_get_hdrspace(ff, &mut keys_exist, &mut more_keys, status);
            fits_movabs_hdu(ff, my_hn, &mut tmp0, status);
            if NGP_OK != *status
            /* error here means file is corrupted */
            {
                return *status;
            }
            if keys_exist > 0 {
                first_extension = 0; /* if keywords exist assume PHDU already exist */
            }
        } else {
            first_extension = 0; /* PHDU (followed by 1+ extensions) exist */

            i = 2;
            while i <= my_hn {
                *status = NGP_OK;
                fits_movabs_hdu(ff, 1, &mut tmp0, status);
                if NGP_OK != *status {
                    break;
                }

                fits_read_key(
                    ff,
                    TSTRING,
                    c"EXTNAME".as_ptr(),
                    used_name.as_mut_ptr().cast::<c_void>(),
                    ptr::null_mut(),
                    status,
                );
                if NGP_OK != *status {
                    continue;
                }

                fits_read_key(
                    ff,
                    TLONG,
                    c"EXTVER".as_ptr(),
                    (&mut luv as *mut c_long).cast::<c_void>(),
                    ptr::null_mut(),
                    status,
                );
                used_ver = luv as c_int; /* bugfix - 22-Jan-99, BO - nonalignment of OSF/Alpha */
                if VALUE_UNDEFINED == *status {
                    used_ver = 1;
                    *status = NGP_OK;
                }

                if NGP_OK == *status {
                    *status = ngp_set_extver(&mut parser_state, &used_name, used_ver);
                }
                i += 1;
            }

            fits_movabs_hdu(ff, my_hn, &mut tmp0, status);
        }

        if NGP_OK != *status {
            return *status;
        }

        // Convert C string pointer to slice for ngp_include_file
        let template_cstr = CStr::from_ptr(ngp_template);
        let template_slice = cast_slice(template_cstr.to_bytes_with_nul());
        *status = ngp_include_file(&mut parser_state, template_slice);
        if NGP_OK != (*status) {
            return *status;
        }

        // Extract directory from template path
        if let Ok(template_str) = CStr::from_ptr(ngp_template).to_str() {
            let template_path = Path::new(template_str);
            if let Some(parent) = template_path.parent() {
                if parent != Path::new("") {
                    parser_state.NGP_MASTER_DIR = parent.to_path_buf();
                }
            }
        }

        loop {
            r = ngp_read_line(&mut parser_state, 1);
            if NGP_OK != r {
                break; /* EOF always means error here */
            }
            match parser_state.NGP_KEYIDX {
                NGP_TOKEN_SIMPLE => {
                    let mut skip = false;
                    if 0 == first_extension
                    /* simple only allowed in first HDU */
                    {
                        r = NGP_TOKEN_NOT_EXPECT;
                        skip = true;
                    }

                    if !skip {
                        r = ngp_unread_line(&mut parser_state);
                        if NGP_OK == r {
                            r = ngp_read_xtension(
                                &mut parser_state,
                                ff,
                                0,
                                NGP_XTENSION_SIMPLE | NGP_XTENSION_FIRST,
                            );
                            first_extension = 0;
                        }
                    }
                }
                NGP_TOKEN_XTENSION => {
                    r = ngp_unread_line(&mut parser_state);
                    let mut skip = false;
                    if NGP_OK != (r) {
                        skip = true;
                    }
                    if !skip {
                        r = ngp_read_xtension(
                            &mut parser_state,
                            ff,
                            0,
                            if first_extension != 0 {
                                NGP_XTENSION_FIRST
                            } else {
                                0
                            },
                        );
                        first_extension = 0;
                    }
                }
                NGP_TOKEN_GROUP => {
                    if NGP_TTYPE_STRING == parser_state.NGP_LINKEY.type_ {
                        strncpy(
                            grnm.as_mut_ptr(),
                            parser_state.NGP_LINKEY.value.s,
                            NGP_MAX_STRING,
                        );
                    } else {
                        int_snprintf!(
                            &mut grnm,
                            NGP_MAX_STRING as size_t,
                            "DEFAULT_GROUP_{}",
                            parser_state.MASTER_GRP_IDX,
                        );
                        parser_state.MASTER_GRP_IDX += 1;
                    }
                    grnm[NGP_MAX_STRING - 1] = 0;
                    r = ngp_read_group(&mut parser_state, ff, &grnm, 0);
                    first_extension = 0;
                }
                NGP_TOKEN_EOF => {
                    exit_flg = 1;
                }
                _ => {
                    r = NGP_TOKEN_NOT_EXPECT;
                }
            }
            if exit_flg != 0 || (NGP_OK != r) {
                break;
            }
        }

        /* all top level HDUs up to faulty one are left intact in case of i/o error. It is up
        to the caller to call fits_close_file or fits_delete_file when this function returns
        error. */

        ngp_free_line(&mut parser_state); /* deallocate last line (if any) */
        ngp_free_prevline(&mut parser_state); /* deallocate cached line (if any) */
        ngp_delete_extver_tab(&mut parser_state); /* delete extver table (if present), error ignored */

        *status = r;
        r
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::ffi::CString;

    // Helper function to create a C string
    fn to_cstring(s: &str) -> CString {
        CString::new(s).unwrap()
    }

    // Helper function to create a test token
    unsafe fn create_test_token(
        type_: c_int,
        name: &str,
        value_int: c_int,
        comment: &str,
    ) -> NgpToken {
        let mut token = NgpToken::default();
        token.type_ = type_;

        // Copy name
        let name_bytes = name.as_bytes();
        let name_len = name_bytes.len().min(NGP_MAX_NAME - 1);
        for i in 0..name_len {
            token.name[i] = name_bytes[i] as c_char;
        }
        token.name[name_len] = 0;

        // Set value based on type
        if type_ == NGP_TTYPE_STRING {
            // For string tokens, the value_int parameter is ignored
            // and we use the name as the string value for testing
            let str_value = to_cstring(name);
            let str_bytes = str_value.as_bytes_with_nul();
            unsafe {
                let str_boxed = match ngp_alloc_boxed(str_bytes.len()) {
                    Ok(b) => b,
                    Err(_) => {
                        token.value = NgpTokval { s: ptr::null_mut() };
                        return token;
                    }
                };
                let ptr = Box::into_raw(str_boxed).cast::<c_char>();
                std::ptr::copy_nonoverlapping(
                    str_bytes.as_ptr().cast::<c_char>(),
                    ptr,
                    str_bytes.len(),
                );
                token.value = NgpTokval { s: ptr };
            }
        } else {
            token.value = NgpTokval { i: value_int };
        }

        // Copy comment
        let comment_bytes = comment.as_bytes();
        let comment_len = comment_bytes.len().min(NGP_MAX_COMMENT - 1);
        for i in 0..comment_len {
            token.comment[i] = comment_bytes[i] as c_char;
        }
        token.comment[comment_len] = 0;

        token
    }

    // Helper function to free memory allocated in test tokens
    unsafe fn free_test_token(token: &mut NgpToken) {
        unsafe {
            if token.type_ == NGP_TTYPE_STRING && !token.value.s.is_null() {
                // Free by reconstructing the Box
                let str_len = 1 + strlen(token.value.s);
                drop(Box::from_raw(slice::from_raw_parts_mut(
                    token.value.s,
                    str_len,
                )));
                token.value.s = ptr::null_mut();
            }
        }
    }

    //
    // Extension Version Management Tests
    //

    #[test]
    fn test_ngp_get_extver_new_name() {
        unsafe {
            let mut parser_state = GRParseState::default();

            // Clean slate
            ngp_delete_extver_tab(&mut parser_state);

            let extname = to_cstring("TEST_EXT");
            let mut version: c_int = 0;

            let status =
                ngp_get_extver(&mut parser_state, cast_slice(extname.as_bytes_with_nul()), &mut version);

            assert_eq!(status, NGP_OK);
            assert_eq!(version, 1);

            // Cleanup
            ngp_delete_extver_tab(&mut parser_state);
        }
    }

    #[test]
    fn test_ngp_get_extver_increment() {
        unsafe {
            let mut parser_state = GRParseState::default();

            let extname = to_cstring("TEST_EXT");
            let mut version1: c_int = 0;
            let mut version2: c_int = 0;
            let mut version3: c_int = 0;

            ngp_get_extver(
                &mut parser_state,
                cast_slice(extname.as_bytes_with_nul()),
                &mut version1,
            );
            ngp_get_extver(
                &mut parser_state,
                cast_slice(extname.as_bytes_with_nul()),
                &mut version2,
            );
            ngp_get_extver(
                &mut parser_state,
                cast_slice(extname.as_bytes_with_nul()),
                &mut version3,
            );

            assert_eq!(version1, 1);
            assert_eq!(version2, 2);
            assert_eq!(version3, 3);

            // Cleanup
            ngp_delete_extver_tab(&mut parser_state);
        }
    }

    #[test]
    fn test_ngp_get_extver_multiple_names() {
        unsafe {
            let mut parser_state = GRParseState::default();

            ngp_delete_extver_tab(&mut parser_state);

            let ext1 = to_cstring("EXT1");
            let ext2 = to_cstring("EXT2");
            let mut ver_ext1: c_int = 0;
            let mut ver_ext2: c_int = 0;
            let mut ver_ext1_again: c_int = 0;

            ngp_get_extver(&mut parser_state, cast_slice(ext1.as_bytes_with_nul()), &mut ver_ext1);
            ngp_get_extver(&mut parser_state, cast_slice(ext2.as_bytes_with_nul()), &mut ver_ext2);
            ngp_get_extver(
                &mut parser_state,
                cast_slice(ext1.as_bytes_with_nul()),
                &mut ver_ext1_again,
            );

            assert_eq!(ver_ext1, 1);
            assert_eq!(ver_ext2, 1);
            assert_eq!(ver_ext1_again, 2);

            ngp_delete_extver_tab(&mut parser_state);
        }
    }

    #[test]
    fn test_ngp_get_extver_null_pointer() {
        unsafe {
            let mut parser_state = GRParseState::default();
            ngp_delete_extver_tab(&mut parser_state);

            let extname = to_cstring("TEST");
            let mut version: c_int = 0;

            // Null extname
            let status = ngp_get_extver(&mut parser_state, &[], &mut version);
            assert_eq!(status, NGP_BAD_ARG, "Should fail with null extname");

            ngp_delete_extver_tab(&mut parser_state);
        }
    }

    #[test]
    fn test_ngp_set_extver_new_name() {
        let mut parser_state = GRParseState::default();
        ngp_delete_extver_tab(&mut parser_state);

        let extname = to_cstring("TEST_EXT");

        let status = ngp_set_extver(&mut parser_state, cast_slice(extname.as_bytes_with_nul()), 5);
        assert_eq!(status, NGP_OK);

        ngp_delete_extver_tab(&mut parser_state);
    }

    #[test]
    fn test_ngp_set_extver_null_pointer() {
        let mut parser_state = GRParseState::default();
        let status = ngp_set_extver(&mut parser_state, &[], 1);
        assert_eq!(status, NGP_BAD_ARG);
    }

    #[test]
    fn test_ngp_set_extver_update_higher() {
        unsafe {
            let mut parser_state = GRParseState::default();
            ngp_delete_extver_tab(&mut parser_state);

            let extname = to_cstring("TEST_EXT");

            // Set to 3 first
            let status = ngp_set_extver(&mut parser_state, cast_slice(extname.as_bytes_with_nul()), 3);
            assert_eq!(status, NGP_OK);

            // Set to 5 - should update to higher value
            let status = ngp_set_extver(&mut parser_state, cast_slice(extname.as_bytes_with_nul()), 5);
            assert_eq!(status, NGP_OK);

            // Next get should return 6 (5 + 1)
            let mut version: c_int = 0;
            ngp_get_extver(&mut parser_state, cast_slice(extname.as_bytes_with_nul()), &mut version);
            assert_eq!(
                version, 6,
                "Should update to higher version and return 6 on next get"
            );

            ngp_delete_extver_tab(&mut parser_state);
        }
    }

    #[test]
    fn test_ngp_set_extver_keep_higher() {
        unsafe {
            let mut parser_state = GRParseState::default();
            ngp_delete_extver_tab(&mut parser_state);

            let extname = to_cstring("TEST_EXT");

            // Set to 5 first
            let status = ngp_set_extver(&mut parser_state, cast_slice(extname.as_bytes_with_nul()), 5);
            assert_eq!(status, NGP_OK);

            // Set to 3 - should keep 5 (the higher value)
            let status = ngp_set_extver(&mut parser_state, cast_slice(extname.as_bytes_with_nul()), 3);
            assert_eq!(status, NGP_OK);

            // Next get should return 6 (5 + 1, not 4)
            let mut version: c_int = 0;
            ngp_get_extver(&mut parser_state, cast_slice(extname.as_bytes_with_nul()), &mut version);
            assert_eq!(
                version, 6,
                "Should keep higher version 5 and return 6 on next get"
            );

            ngp_delete_extver_tab(&mut parser_state);
        }
    }

    #[test]
    fn test_ngp_delete_extver_tab_empty() {
        let mut parser_state = GRParseState::default();
        ngp_delete_extver_tab(&mut parser_state);

        // Deleting empty table should succeed
        let status = ngp_delete_extver_tab(&mut parser_state);
        assert_eq!(status, NGP_OK);
    }

    #[test]
    fn test_ngp_delete_extver_tab_populated() {
        unsafe {
            let mut parser_state = GRParseState::default();
            ngp_delete_extver_tab(&mut parser_state);

            let ext1 = to_cstring("EXT1");
            let ext2 = to_cstring("EXT2");
            let mut ver1: c_int = 0;
            let mut ver2: c_int = 0;

            ngp_get_extver(&mut parser_state, cast_slice(ext1.as_bytes_with_nul()), &mut ver1);
            ngp_get_extver(&mut parser_state, cast_slice(ext2.as_bytes_with_nul()), &mut ver2);

            let status = ngp_delete_extver_tab(&mut parser_state);
            assert_eq!(status, NGP_OK);
        }
    }

    //
    // HDU Structure Management Tests
    //

    #[test]
    fn test_ngp_hdu_init_basic() {
        let mut hdu = NgpHdu {
            tokcnt: 999,
            tok: Vec::new(), // Will be initialized properly
        };

        let status = ngp_hdu_init(&mut hdu);

        assert_eq!(status, NGP_OK);
        assert_eq!(hdu.tokcnt, 0);
        assert!(hdu.tok.is_empty());
    }

    #[test]
    fn test_ngp_hdu_insert_token_first() {
        unsafe {
            let mut hdu = NgpHdu {
                tokcnt: 0,
                tok: Vec::new(),
            };
            ngp_hdu_init(&mut hdu);

            let mut token = create_test_token(NGP_TTYPE_INT, "NAXIS", 2, "Number of axes");

            let status = ngp_hdu_insert_token(&mut hdu, &mut token);

            assert_eq!(status, NGP_OK);
            assert_eq!(hdu.tokcnt, 1);
            assert!(!hdu.tok.is_empty());

            // Cleanup
            ngp_hdu_clear(&mut hdu);
        }
    }

    #[test]
    fn test_ngp_hdu_insert_token_string() {
        unsafe {
            let mut hdu = NgpHdu {
                tokcnt: 0,
                tok: Vec::new(),
            };
            ngp_hdu_init(&mut hdu);

            let mut token = create_test_token(NGP_TTYPE_STRING, "EXTNAME", 0, "Extension name");

            let status = ngp_hdu_insert_token(&mut hdu, &mut token);

            assert_eq!(status, NGP_OK);
            assert_eq!(hdu.tokcnt, 1);

            // Cleanup
            free_test_token(&mut token);
            ngp_hdu_clear(&mut hdu);
        }
    }

    #[test]
    fn test_ngp_hdu_insert_token_multiple() {
        unsafe {
            let mut hdu = NgpHdu {
                tokcnt: 0,
                tok: Vec::new(),
            };
            ngp_hdu_init(&mut hdu);

            // Insert first token
            let mut token1 = create_test_token(NGP_TTYPE_INT, "NAXIS", 2, "Number of axes");
            ngp_hdu_insert_token(&mut hdu, &mut token1);

            // Insert second token
            let mut token2 = create_test_token(NGP_TTYPE_INT, "BITPIX", 16, "Bits per pixel");
            ngp_hdu_insert_token(&mut hdu, &mut token2);

            // Insert third token
            let mut token3 = create_test_token(NGP_TTYPE_BOOL, "SIMPLE", 1, "Standard FITS");
            ngp_hdu_insert_token(&mut hdu, &mut token3);

            assert_eq!(hdu.tokcnt, 3, "Should have 3 tokens");

            // Cleanup
            ngp_hdu_clear(&mut hdu);
        }
    }

    #[test]
    fn test_ngp_hdu_clear_with_data() {
        unsafe {
            let mut hdu = NgpHdu {
                tokcnt: 0,
                tok: Vec::new(),
            };
            ngp_hdu_init(&mut hdu);

            // Add some tokens
            let mut token1 = create_test_token(NGP_TTYPE_INT, "NAXIS", 2, "");
            let mut token2 = create_test_token(NGP_TTYPE_INT, "BITPIX", 16, "");
            ngp_hdu_insert_token(&mut hdu, &mut token1);
            ngp_hdu_insert_token(&mut hdu, &mut token2);

            assert_eq!(hdu.tokcnt, 2);
            assert!(!hdu.tok.is_empty());

            // Clear
            let status = ngp_hdu_clear(&mut hdu);

            assert_eq!(status, NGP_OK);
            assert_eq!(hdu.tokcnt, 0);
            assert!(hdu.tok.is_empty());
        }
    }

    #[test]
    fn test_ngp_hdu_clear_empty() {
        unsafe {
            let mut hdu = NgpHdu {
                tokcnt: 0,
                tok: Vec::new(),
            };
            ngp_hdu_init(&mut hdu);

            // Clear an already empty HDU
            let status = ngp_hdu_clear(&mut hdu);
            assert_eq!(status, NGP_OK);
        }
    }

    //
    // Line Management Tests
    //

    #[test]
    fn test_ngp_free_line_empty() {
        // Should not crash when freeing empty line
        let mut parser_state = GRParseState::default();
        let status = ngp_free_line(&mut parser_state);
        assert_eq!(status, NGP_OK);
    }

    #[test]
    fn test_ngp_free_prevline_empty() {
        // Should not crash when freeing empty prevline
        let mut parser_state = GRParseState::default();
        let status = ngp_free_prevline(&mut parser_state);
        assert_eq!(status, NGP_OK);
    }

    //
    // Keyword Validation Tests
    //

    #[test]
    fn test_ngp_keyword_is_write_user_keyword() {
        unsafe {
            let name = to_cstring("MYKEY");
            let tok = create_test_token(NGP_TTYPE_INT, "MYKEY", 42, "User keyword");

            let status = ngp_keyword_is_write(&tok);
            assert_eq!(status, NGP_OK, "User keyword should be allowed");
        }
    }

    #[test]
    fn test_ngp_keyword_is_write_simple() {
        unsafe {
            let tok = create_test_token(NGP_TTYPE_BOOL, "SIMPLE", 1, "Simple keyword");

            let status = ngp_keyword_is_write(&tok);
            assert_eq!(status, NGP_BAD_ARG, "SIMPLE keyword should not be writable");
        }
    }

    #[test]
    fn test_ngp_keyword_is_write_xtension() {
        unsafe {
            let tok = create_test_token(NGP_TTYPE_INT, "XTENSION", 0, "Extension keyword");

            let status = ngp_keyword_is_write(&tok);
            assert_eq!(
                status, NGP_BAD_ARG,
                "XTENSION keyword should not be writable"
            );
        }
    }

    #[test]
    fn test_ngp_keyword_is_write_bitpix() {
        unsafe {
            let tok = create_test_token(NGP_TTYPE_INT, "BITPIX", 16, "Bitpix keyword");

            let status = ngp_keyword_is_write(&tok);
            assert_eq!(status, NGP_BAD_ARG, "BITPIX keyword should not be writable");
        }
    }

    #[test]
    fn test_ngp_keyword_is_write_naxis() {
        unsafe {
            let tok = create_test_token(NGP_TTYPE_INT, "NAXIS", 2, "Naxis keyword");

            let status = ngp_keyword_is_write(&tok);
            assert_eq!(status, NGP_BAD_ARG, "NAXIS keyword should not be writable");
        }
    }

    #[test]
    fn test_ngp_keyword_is_write_naxis1() {
        unsafe {
            let tok = create_test_token(NGP_TTYPE_INT, "NAXIS1", 100, "Naxis1 keyword");

            let status = ngp_keyword_is_write(&tok);
            assert_eq!(status, NGP_BAD_ARG, "NAXIS1 keyword should not be writable");
        }
    }

    #[test]
    fn test_ngp_keyword_is_write_tform1() {
        unsafe {
            let tok = create_test_token(NGP_TTYPE_INT, "TFORM1", 0, "Tform1 keyword");

            let status = ngp_keyword_is_write(&tok);
            assert_eq!(status, NGP_BAD_ARG, "TFORM1 keyword should not be writable");
        }
    }

    #[test]
    fn test_ngp_keyword_is_write_ttype5() {
        unsafe {
            let tok = create_test_token(NGP_TTYPE_INT, "TTYPE5", 0, "Ttype5 keyword");

            let status = ngp_keyword_is_write(&tok);
            assert_eq!(status, NGP_BAD_ARG, "TTYPE5 keyword should not be writable");
        }
    }

    #[test]
    fn test_ngp_keyword_is_write_extver() {
        unsafe {
            let tok = create_test_token(NGP_TTYPE_INT, "EXTVER", 1, "Extver keyword");

            let status = ngp_keyword_is_write(&tok);
            assert_eq!(status, NGP_BAD_ARG, "EXTVER keyword should not be writable");
        }
    }

    #[test]
    fn test_ngp_keyword_is_write_null_pointer() {
        unsafe {
            let status = ngp_keyword_is_write(ptr::null());
            assert_eq!(
                status, NGP_NUL_PTR,
                "Null pointer should return NGP_NUL_PTR"
            );
        }
    }

    //
    // Token Extraction Tests
    //

    // Helper to compare C strings
    unsafe fn cstr_eq(ptr: *const c_char, expected: &str) -> bool {
        unsafe {
            if ptr.is_null() {
                return false;
            }
            let c_str = CStr::from_ptr(ptr);
            c_str.to_str().unwrap() == expected
        }
    }

    // Helper to create NgpRawLine from string for testing
    fn create_test_line(s: &str) -> NgpRawLine {
        let line_str = to_cstring(s);
        let bytes = cast_slice(line_str.as_bytes_with_nul());
        let mut vec: Vec<c_char> = Vec::with_capacity(bytes.len());
        vec.extend_from_slice(bytes);
        NgpRawLine {
            line: vec.into_boxed_slice(),
            name_idx: None,
            value_idx: None,
            type_: 0,
            comment_idx: None,
            format: 0,
            flags: 0,
        }
    }

    #[test]
    fn test_ngp_extract_tokens_simple_keyword() {
        unsafe {
            let mut line = create_test_line("KEYWORD");

            let status = ngp_extract_tokens(&mut line);
            assert_eq!(status, NGP_OK);
            assert!(line.name_idx.is_some(), "Name should be set");
            let name_slice = &line.line[line.name_idx.unwrap()..];
            assert!(
                cstr_eq(name_slice.as_ptr(), "KEYWORD"),
                "Name should be KEYWORD"
            );
            assert!(line.value_idx.is_none(), "Value should be null");
            assert!(line.comment_idx.is_none(), "Comment should be null");
        }
    }

    #[test]
    fn test_ngp_extract_tokens_keyword_with_int() {
        unsafe {
            let mut line = create_test_line("BITPIX = 16");

            let status = ngp_extract_tokens(&mut line);
            assert_eq!(status, NGP_OK);
            assert!(
                cstr_eq(line.line[line.name_idx.unwrap()..].as_ptr(), "BITPIX"),
                "Name should be BITPIX"
            );
            assert!(
                cstr_eq(line.line[line.value_idx.unwrap()..].as_ptr(), "16"),
                "Value should be 16"
            );
        }
    }

    #[test]
    fn test_ngp_extract_tokens_keyword_with_string() {
        unsafe {
            let mut line = create_test_line("EXTNAME = 'test value'");

            let status = ngp_extract_tokens(&mut line);
            assert_eq!(status, NGP_OK);
            assert!(
                cstr_eq(line.line[line.name_idx.unwrap()..].as_ptr(), "EXTNAME"),
                "Name should be EXTNAME"
            );
            assert!(
                cstr_eq(line.line[line.value_idx.unwrap()..].as_ptr(), "test value"),
                "Value should be 'test value'"
            );
            assert_eq!(line.type_, NGP_TTYPE_STRING, "Type should be STRING");
        }
    }

    #[test]
    fn test_ngp_extract_tokens_keyword_with_comment() {
        unsafe {
            let mut line = create_test_line("BITPIX = 16 / bits per pixel");

            let status = ngp_extract_tokens(&mut line);
            assert_eq!(status, NGP_OK);
            assert!(
                cstr_eq(line.line[line.name_idx.unwrap()..].as_ptr(), "BITPIX"),
                "Name should be BITPIX"
            );
            assert!(
                cstr_eq(line.line[line.value_idx.unwrap()..].as_ptr(), "16"),
                "Value should be 16"
            );
            assert!(
                cstr_eq(
                    line.line[line.comment_idx.unwrap()..].as_ptr(),
                    "bits per pixel"
                ),
                "Comment should be 'bits per pixel'"
            );
        }
    }

    #[test]
    fn test_ngp_extract_tokens_comment_only() {
        unsafe {
            let mut line = create_test_line("KEYWORD = / comment only");

            let status = ngp_extract_tokens(&mut line);
            assert_eq!(status, NGP_OK);
            assert!(
                cstr_eq(line.line[line.name_idx.unwrap()..].as_ptr(), "KEYWORD"),
                "Name should be KEYWORD"
            );
            assert!(line.value_idx.is_none(), "Value should be null");
            assert!(
                cstr_eq(
                    line.line[line.comment_idx.unwrap()..].as_ptr(),
                    "comment only"
                ),
                "Comment should be 'comment only'"
            );
        }
    }

    #[test]
    fn test_ngp_extract_tokens_history() {
        unsafe {
            let mut line = create_test_line("HISTORY This is a history record");

            let status = ngp_extract_tokens(&mut line);
            assert_eq!(status, NGP_OK);
            assert!(
                cstr_eq(line.line[line.name_idx.unwrap()..].as_ptr(), "HISTORY"),
                "Name should be HISTORY"
            );
            assert!(
                cstr_eq(
                    line.line[line.comment_idx.unwrap()..].as_ptr(),
                    "This is a history record"
                ),
                "Comment should be 'This is a history record'"
            );
            assert_eq!(line.type_, NGP_TTYPE_RAW, "Type should be RAW");
        }
    }

    #[test]
    fn test_ngp_extract_tokens_comment_keyword() {
        unsafe {
            let mut line = create_test_line("COMMENT This is a comment");

            let status = ngp_extract_tokens(&mut line);
            assert_eq!(status, NGP_OK);
            assert!(
                cstr_eq(line.line[line.name_idx.unwrap()..].as_ptr(), "COMMENT"),
                "Name should be COMMENT"
            );
            assert!(
                cstr_eq(
                    line.line[line.comment_idx.unwrap()..].as_ptr(),
                    "This is a comment"
                ),
                "Comment should be 'This is a comment'"
            );
            assert_eq!(line.type_, NGP_TTYPE_RAW, "Type should be RAW");
        }
    }

    #[test]
    fn test_ngp_extract_tokens_null_pointer() {
        let mut line = NgpRawLine {
            line: Box::new([]),
            name_idx: None,
            value_idx: None,
            type_: 0,
            comment_idx: None,
            format: 0,
            flags: 0,
        };

        let status = ngp_extract_tokens(&mut line);
        assert_eq!(status, NGP_NUL_PTR, "Empty line should return NGP_NUL_PTR");
    }

    #[test]
    fn test_ngp_extract_tokens_string_with_double_quotes() {
        unsafe {
            let mut line = create_test_line("TEST = 'it''s a test'");

            let status = ngp_extract_tokens(&mut line);
            assert_eq!(status, NGP_OK);
            assert!(
                cstr_eq(line.line[line.value_idx.unwrap()..].as_ptr(), "it's a test"),
                "Value should convert '' to '"
            );
        }
    }

    #[test]
    fn test_ngp_extract_tokens_blank_line_8spaces() {
        unsafe {
            let mut line = create_test_line("        This is a blank keyword");

            let status = ngp_extract_tokens(&mut line);
            assert_eq!(status, NGP_OK);
            assert!(
                cstr_eq(line.line[line.name_idx.unwrap()..].as_ptr(), ""),
                "Name should be empty for 8-space line"
            );
            assert!(
                cstr_eq(
                    line.line[line.comment_idx.unwrap()..].as_ptr(),
                    "This is a blank keyword"
                ),
                "Comment should be preserved"
            );
            assert_eq!(line.type_, NGP_TTYPE_RAW, "Type should be RAW");
        }
    }

    #[test]
    fn test_ngp_extract_tokens_empty_line() {
        let mut line = create_test_line("");

        let status = ngp_extract_tokens(&mut line);
        assert_eq!(status, NGP_OK);
        assert_eq!(line.type_, NGP_TTYPE_RAW, "Empty line should have RAW type");
    }

    #[test]
    fn test_ngp_extract_tokens_hierarch() {
        unsafe {
            let mut line = create_test_line("HIERARCH LongKeywordName = 'value'");

            let status = ngp_extract_tokens(&mut line);
            assert_eq!(status, NGP_OK);

            // Name should contain HIERARCH
            if let Some(name_idx) = line.name_idx {
                let name_str = CStr::from_ptr(line.line[name_idx..].as_ptr())
                    .to_str()
                    .unwrap();
                assert!(
                    name_str.contains("HIERARCH"),
                    "Name should contain HIERARCH"
                );
            }
            assert!(
                cstr_eq(line.line[line.value_idx.unwrap()..].as_ptr(), "value"),
                "Value should be 'value'"
            );
        }
    }
}
