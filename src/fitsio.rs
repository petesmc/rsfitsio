use crate::c_types::c_ulonglong;
use bytemuck::{cast_slice, cast_slice_mut};
use std::{
    cmp,
    ffi::{CStr, c_char, c_int, c_long, c_longlong},
    os::raw::c_void,
};

use crate::{fitscore::ffpmsg_slice, int_snprintf, slice_to_str, wrappers::strlen_safe};

use crate::cfileio::fitsdriver;

pub const HAVE_FTRUNCATE: bool = true;

#[macro_export]
macro_rules! BL {
    () => {
        2880
    };
}

pub const NULL_MSG: &str = "Null Pointer";

pub const BLOCK_LEN: usize = 2880;

pub const CFITSIO_VERSION: [u8; 6] = *b"4.6.2\0";

/* Minor and micro numbers must not exceed 99 under current method of version representataion in ffvers(). */
pub const CFITSIO_MICRO: u64 = 2;
pub const CFITSIO_MINOR: u64 = 6;
pub const CFITSIO_MAJOR: u64 = 4;
pub const CFITSIO_SONAME: u64 = 10;

/* the SONAME is incremented in a new release if the binary shared */
/* library (on linux and Mac systems) is not backward compatible */
/* with the previous release of CFITSIO */

pub const USE_LL_SUFFIX: u64 = 1;

pub type LONGLONG = c_longlong;
pub type ULONGLONG = c_ulonglong;

/*  define a default value, even if it is never used */
pub const LONGLONG_MAX: c_longlong = c_longlong::MAX;
pub const LONGLONG_MIN: c_longlong = c_longlong::MIN;

pub const LONG_MAX: c_long = c_long::MAX;
pub const LONG_MIN: c_long = c_long::MIN;

pub const NIOBUF: u64 = 40; /* number of IO buffers to create (default = 40) */
/* !! Significantly increasing NIOBUF may degrade performance !! */

pub const IOBUFLEN: i64 = 2880; /* size in bytes of each IO buffer (DONT CHANGE!) */

/* global variables */

pub const FLEN_FILENAME: usize = 1025; /* max length of a filename  */
pub const FLEN_KEYWORD: usize = 75; /* max length of a keyword (HIERARCH convention) */
pub const FLEN_CARD: usize = 81; /* length of a FITS header card */
pub const FLEN_VALUE: usize = 71; /* max length of a keyword value string */
pub const FLEN_COMMENT: usize = 73; /* max length of a keyword comment string */
pub const FLEN_ERRMSG: usize = 81; /* max length of a FITSIO error message */
pub const FLEN_STATUS: usize = 31; /* max length of a FITSIO status text string */

pub const TBIT: c_int = 1; /* codes for FITS table data types */
pub const TBYTE: c_int = 11;
pub const TSBYTE: c_int = 12;
pub const TLOGICAL: c_int = 14;
pub const TSTRING: c_int = 16;
pub const TUSHORT: c_int = 20;
pub const TSHORT: c_int = 21;
pub const TUINT: c_int = 30;
pub const TINT: c_int = 31;
pub const TULONG: c_int = 40;
pub const TLONG: c_int = 41;
pub const TINT32BIT: c_int = 41; /* used when returning datatype of a column */
pub const TFLOAT: c_int = 42;
pub const TULONGLONG: c_int = 80;
pub const TLONGLONG: c_int = 81;
pub const TDOUBLE: c_int = 82;
pub const TCOMPLEX: c_int = 83;
pub const TDBLCOMPLEX: c_int = 163;

pub const TYP_STRUC_KEY: c_int = 10;
pub const TYP_CMPRS_KEY: c_int = 20;
pub const TYP_SCAL_KEY: c_int = 30;
pub const TYP_NULL_KEY: c_int = 40;
pub const TYP_DIM_KEY: c_int = 50;
pub const TYP_RANG_KEY: c_int = 60;
pub const TYP_UNIT_KEY: c_int = 70;
pub const TYP_DISP_KEY: c_int = 80;
pub const TYP_HDUID_KEY: c_int = 90;
pub const TYP_CKSUM_KEY: c_int = 100;
pub const TYP_WCS_KEY: c_int = 110;
pub const TYP_REFSYS_KEY: c_int = 120;
pub const TYP_COMM_KEY: c_int = 130;
pub const TYP_CONT_KEY: c_int = 140;
pub const TYP_USER_KEY: c_int = 150;

pub type INT32BIT = c_int; /* 32-bit integer datatype.  Currently this       */
/* datatype is an 'int' on all useful platforms   */
/* however, it is possible that that are cases    */
/* where 'int' is a 2-byte integer, in which case */
/* INT32BIT would need to be defined as 'long'.   */

pub const BYTE_IMG: c_int = 8; /* BITPIX code values for FITS image types */
pub const SHORT_IMG: c_int = 16;
pub const LONG_IMG: c_int = 32;
pub const LONGLONG_IMG: c_int = 64;
pub const FLOAT_IMG: c_int = -32;
pub const DOUBLE_IMG: c_int = -64;
/* The following 2 codes are not true FITS         */
/* datatypes; these codes are only used internally */
/* within cfitsio to make it easier for users      */
/* to deal with unsigned integers.                 */
pub const SBYTE_IMG: c_int = 10;
pub const USHORT_IMG: c_int = 20;
pub const ULONG_IMG: c_int = 40;
pub const ULONGLONG_IMG: c_int = 80;

pub const IMAGE_HDU: c_int = 0; /* Primary Array or IMAGE HDU */
pub const ASCII_TBL: c_int = 1; /* ASCII table HDU  */
pub const BINARY_TBL: c_int = 2; /* Binary table HDU */
pub const ANY_HDU: c_int = -1; /* matches any HDU type */

pub const READONLY: c_int = 0; /* options when opening a file */
pub const READWRITE: c_int = 1;

/* adopt a hopefully obscure number to use as a null value flag */
pub const FLOATNULLVALUE: f32 = -9.11912E-36;
pub const DOUBLENULLVALUE: f64 = -9.1191291391491E-36;

/* compression algorithm codes */
pub const NO_DITHER: c_int = -1;
pub const SUBTRACTIVE_DITHER_1: c_int = 1;
pub const SUBTRACTIVE_DITHER_2: c_int = 2;

pub const MAX_COMPRESS_DIM: usize = 6;
pub const RICE_1: c_int = 11;
pub const GZIP_1: c_int = 21;
pub const GZIP_2: c_int = 22;
pub const PLIO_1: c_int = 31;
pub const HCOMPRESS_1: c_int = 41;
pub const BZIP2_1: c_int = 51; /* not publicly supported; only for test purposes */
pub const NOCOMPRESS: c_int = -1;

pub const TRUE: u64 = 1;

pub const FALSE: u64 = 0;

pub const CASESEN: u64 = 1; /* do case-sensitive string match */
pub const CASEINSEN: u64 = 0; /* do case-insensitive string match */

pub const GT_ID_ALL_URI: u64 = 0; /* hierarchical grouping parameters */
pub const GT_ID_REF: u64 = 1;
pub const GT_ID_POS: u64 = 2;
pub const GT_ID_ALL: u64 = 3;
pub const GT_ID_REF_URI: u64 = 11;
pub const GT_ID_POS_URI: u64 = 12;

pub const OPT_RM_GPT: u64 = 0;
pub const OPT_RM_ENTRY: u64 = 1;
pub const OPT_RM_MBR: u64 = 2;
pub const OPT_RM_ALL: u64 = 3;

pub const OPT_GCP_GPT: u64 = 0;
pub const OPT_GCP_MBR: u64 = 1;
pub const OPT_GCP_ALL: u64 = 2;

pub const OPT_MCP_ADD: u64 = 0;
pub const OPT_MCP_NADD: u64 = 1;
pub const OPT_MCP_REPL: u64 = 2;
pub const OPT_MCP_MOV: u64 = 3;

pub const OPT_MRG_COPY: u64 = 0;
pub const OPT_MRG_MOV: u64 = 1;

pub const OPT_CMT_MBR: u64 = 1;
pub const OPT_CMT_MBR_DEL: u64 = 11;

pub const VALIDSTRUC: c_int = 555; /* magic value used to identify if structure is valid */

pub const INPUT_COL: u64 = 0; /* flag for input only iterator column       */
pub const INPUT_OUTPUT_COL: u64 = 1; /* flag for input and output iterator column */
pub const OUTPUT_COL: u64 = 2; /* flag for output only iterator column      */
pub const TEMPORARY_COL: u64 = 3; /* flag for temporary iterator column INTERNAL */

/* error status codes */

pub const CREATE_DISK_FILE: c_int = -106; /* create disk file, without extended filename syntax */
pub const OPEN_DISK_FILE: c_int = -105; /* open disk file, without extended filename syntax */
pub const SKIP_TABLE: c_int = -104; /* move to 1st image when opening file */
pub const SKIP_IMAGE: c_int = -103; /* move to 1st table when opening file */
pub const SKIP_NULL_PRIMARY: c_int = -102; /* skip null primary array when opening file */
pub const USE_MEM_BUFF: c_int = -101; /* use memory buffer when opening file */
pub const OVERFLOW_ERR: c_int = -11; /* overflow during datatype conversion */
pub const PREPEND_PRIMARY: c_int = -9; /* used in ffiimg to insert new primary array */
pub const SAME_FILE: c_int = 101; /* input and output files are the same */
pub const TOO_MANY_FILES: c_int = 103; /* tried to open too many FITS files */
pub const FILE_NOT_OPENED: c_int = 104; /* could not open the named file */
pub const FILE_NOT_CREATED: c_int = 105; /* could not create the named file */
pub const WRITE_ERROR: c_int = 106; /* error writing to FITS file */
pub const END_OF_FILE: c_int = 107; /* tried to move past end of file */
pub const READ_ERROR: c_int = 108; /* error reading from FITS file */
pub const FILE_NOT_CLOSED: c_int = 110; /* could not close the file */
pub const ARRAY_TOO_BIG: c_int = 111; /* array dimensions exceed internal limit */
pub const READONLY_FILE: c_int = 112; /* Cannot write to readonly file */
pub const MEMORY_ALLOCATION: c_int = 113; /* Could not allocate memory */
pub const BAD_FILEPTR: c_int = 114; /* invalid fitsfile pointer */
pub const NULL_INPUT_PTR: c_int = 115; /* NULL input pointer to routine */
pub const SEEK_ERROR: c_int = 116; /* error seeking position in file */
pub const BAD_NETTIMEOUT: c_int = 117; /* bad value for file download timeout setting */

pub const BAD_URL_PREFIX: c_int = 121; /* invalid URL prefix on file name */
pub const TOO_MANY_DRIVERS: c_int = 122; /* tried to register too many IO drivers */
pub const DRIVER_INIT_FAILED: c_int = 123; /* driver initialization failed */
pub const NO_MATCHING_DRIVER: c_int = 124; /* matching driver is not registered */
pub const URL_PARSE_ERROR: c_int = 125; /* failed to parse input file URL */
pub const RANGE_PARSE_ERROR: c_int = 126; /* failed to parse input file URL */

pub const SHARED_ERRBASE: c_int = 150;
pub const SHARED_BADARG: c_int = SHARED_ERRBASE + 1;
pub const SHARED_NULPTR: c_int = SHARED_ERRBASE + 2;
pub const SHARED_TABFULL: c_int = SHARED_ERRBASE + 3;
pub const SHARED_NOTINIT: c_int = SHARED_ERRBASE + 4;
pub const SHARED_IPCERR: c_int = SHARED_ERRBASE + 5;
pub const SHARED_NOMEM: c_int = SHARED_ERRBASE + 6;
pub const SHARED_AGAIN: c_int = SHARED_ERRBASE + 7;
pub const SHARED_NOFILE: c_int = SHARED_ERRBASE + 8;
pub const SHARED_NORESIZE: c_int = SHARED_ERRBASE + 9;

pub const HEADER_NOT_EMPTY: c_int = 201; /* header already contains keywords */
pub const KEY_NO_EXIST: c_int = 202; /* keyword not found in header */
pub const KEY_OUT_BOUNDS: c_int = 203; /* keyword record number is out of bounds */
pub const VALUE_UNDEFINED: c_int = 204; /* keyword value field is blank */
pub const NO_QUOTE: c_int = 205; /* string is missing the closing quote */
pub const BAD_INDEX_KEY: c_int = 206; /* illegal indexed keyword name */
pub const BAD_KEYCHAR: c_int = 207; /* illegal character in keyword name or card */
pub const BAD_ORDER: c_int = 208; /* required keywords out of order */
pub const NOT_POS_INT: c_int = 209; /* keyword value is not a positive integer */
pub const NO_END: c_int = 210; /* couldn't find END keyword */
pub const BAD_BITPIX: c_int = 211; /* illegal BITPIX keyword value*/
pub const BAD_NAXIS: c_int = 212; /* illegal NAXIS keyword value */
pub const BAD_NAXES: c_int = 213; /* illegal NAXISn keyword value */
pub const BAD_PCOUNT: c_int = 214; /* illegal PCOUNT keyword value */
pub const BAD_GCOUNT: c_int = 215; /* illegal GCOUNT keyword value */
pub const BAD_TFIELDS: c_int = 216; /* illegal TFIELDS keyword value */
pub const NEG_WIDTH: c_int = 217; /* negative table row size */
pub const NEG_ROWS: c_int = 218; /* negative number of rows in table */
pub const COL_NOT_FOUND: c_int = 219; /* column with this name not found in table */
pub const BAD_SIMPLE: c_int = 220; /* illegal value of SIMPLE keyword  */
pub const NO_SIMPLE: c_int = 221; /* Primary array doesn't start with SIMPLE */
pub const NO_BITPIX: c_int = 222; /* Second keyword not BITPIX */
pub const NO_NAXIS: c_int = 223; /* Third keyword not NAXIS */
pub const NO_NAXES: c_int = 224; /* Couldn't find all the NAXISn keywords */
pub const NO_XTENSION: c_int = 225; /* HDU doesn't start with XTENSION keyword */
pub const NOT_ATABLE: c_int = 226; /* the CHDU is not an ASCII table extension */
pub const NOT_BTABLE: c_int = 227; /* the CHDU is not a binary table extension */
pub const NO_PCOUNT: c_int = 228; /* couldn't find PCOUNT keyword */
pub const NO_GCOUNT: c_int = 229; /* couldn't find GCOUNT keyword */
pub const NO_TFIELDS: c_int = 230; /* couldn't find TFIELDS keyword */
pub const NO_TBCOL: c_int = 231; /* couldn't find TBCOLn keyword */
pub const NO_TFORM: c_int = 232; /* couldn't find TFORMn keyword */
pub const NOT_IMAGE: c_int = 233; /* the CHDU is not an IMAGE extension */
pub const BAD_TBCOL: c_int = 234; /* TBCOLn keyword value < 0 or > rowlength */
pub const NOT_TABLE: c_int = 235; /* the CHDU is not a table */
pub const COL_TOO_WIDE: c_int = 236; /* column is too wide to fit in table */
pub const COL_NOT_UNIQUE: c_int = 237; /* more than 1 column name matches template */
pub const BAD_ROW_WIDTH: c_int = 241; /* sum of column widths not = NAXIS1 */
pub const UNKNOWN_EXT: c_int = 251; /* unrecognizable FITS extension type */
pub const UNKNOWN_REC: c_int = 252; /* unrecognizable FITS record */
pub const END_JUNK: c_int = 253; /* END keyword is not blank */
pub const BAD_HEADER_FILL: c_int = 254; /* Header fill area not blank */
pub const BAD_DATA_FILL: c_int = 255; /* Data fill area not blank or zero */
pub const BAD_TFORM: c_int = 261; /* illegal TFORM format code */
pub const BAD_TFORM_DTYPE: c_int = 262; /* unrecognizable TFORM datatype code */
pub const BAD_TDIM: c_int = 263; /* illegal TDIMn keyword value */
pub const BAD_HEAP_PTR: c_int = 264; /* invalid BINTABLE heap address */

pub const BAD_HDU_NUM: c_int = 301; /* HDU number < 1 or > MAXHDU */
pub const BAD_COL_NUM: c_int = 302; /* column number < 1 or > tfields */
pub const NEG_FILE_POS: c_int = 304; /* tried to move before beginning of file  */
pub const NEG_BYTES: c_int = 306; /* tried to read or write negative bytes */
pub const BAD_ROW_NUM: c_int = 307; /* illegal starting row number in table */
pub const BAD_ELEM_NUM: c_int = 308; /* illegal starting element number in vector */
pub const NOT_ASCII_COL: c_int = 309; /* this is not an ASCII string column */
pub const NOT_LOGICAL_COL: c_int = 310; /* this is not a logical datatype column */
pub const BAD_ATABLE_FORMAT: c_int = 311; /* ASCII table column has wrong format */
pub const BAD_BTABLE_FORMAT: c_int = 312; /* Binary table column has wrong format */
pub const NO_NULL: c_int = 314; /* null value has not been defined */
pub const NOT_VARI_LEN: c_int = 317; /* this is not a variable length column */
pub const BAD_DIMEN: c_int = 320; /* illegal number of dimensions in array */
pub const BAD_PIX_NUM: c_int = 321; /* first pixel number greater than last pixel */
pub const ZERO_SCALE: c_int = 322; /* illegal BSCALE or TSCALn keyword = 0 */
pub const NEG_AXIS: c_int = 323; /* illegal axis length < 1 */

pub const NOT_GROUP_TABLE: c_int = 340;
pub const HDU_ALREADY_MEMBER: c_int = 341;
pub const MEMBER_NOT_FOUND: c_int = 342;
pub const GROUP_NOT_FOUND: c_int = 343;
pub const BAD_GROUP_ID: c_int = 344;
pub const TOO_MANY_HDUS_TRACKED: c_int = 345;
pub const HDU_ALREADY_TRACKED: c_int = 346;
pub const BAD_OPTION: c_int = 347;
pub const IDENTICAL_POINTERS: c_int = 348;
pub const BAD_GROUP_ATTACH: c_int = 349;
pub const BAD_GROUP_DETACH: c_int = 350;

pub const BAD_I2C: c_int = 401; /* bad int to formatted string conversion */
pub const BAD_F2C: c_int = 402; /* bad float to formatted string conversion */
pub const BAD_INTKEY: c_int = 403; /* can't interprete keyword value as integer */
pub const BAD_LOGICALKEY: c_int = 404; /* can't interprete keyword value as logical */
pub const BAD_FLOATKEY: c_int = 405; /* can't interprete keyword value as float */
pub const BAD_DOUBLEKEY: c_int = 406; /* can't interprete keyword value as double */
pub const BAD_C2I: c_int = 407; /* bad formatted string to int conversion */
pub const BAD_C2F: c_int = 408; /* bad formatted string to float conversion */
pub const BAD_C2D: c_int = 409; /* bad formatted string to double conversion */
pub const BAD_DATATYPE: c_int = 410; /* bad keyword datatype code */
pub const BAD_DECIM: c_int = 411; /* bad number of decimal places specified */
pub const NUM_OVERFLOW: c_int = 412; /* overflow during datatype conversion */

pub const DATA_COMPRESSION_ERR: c_int = 413; /* error in imcompress routines */
pub const DATA_DECOMPRESSION_ERR: c_int = 414; /* error in imcompress routines */
pub const NO_COMPRESSED_TILE: c_int = 415; /* compressed tile doesn't exist */

pub const BAD_DATE: c_int = 420; /* error in date or time conversion */

pub const PARSE_SYNTAX_ERR: c_int = 431; /* syntax error in parser expression */
pub const PARSE_BAD_TYPE: c_int = 432; /* expression did not evaluate to desired type */
pub const PARSE_LRG_VECTOR: c_int = 433; /* vector result too large to return in array */
pub const PARSE_NO_OUTPUT: c_int = 434; /* data parser failed not sent an out column */
pub const PARSE_BAD_COL: c_int = 435; /* bad data encounter while parsing column */
pub const PARSE_BAD_OUTPUT: c_int = 436; /* Output file not of proper type          */

pub const ANGLE_TOO_BIG: c_int = 501; /* celestial angle too large for projection */
pub const BAD_WCS_VAL: c_int = 502; /* bad celestial coordinate or pixel value */
pub const WCS_ERROR: c_int = 503; /* error in celestial coordinate calculation */
pub const BAD_WCS_PROJ: c_int = 504; /* unsupported type of celestial projection */
pub const NO_WCS_KEY: c_int = 505; /* celestial coordinate keywords not found */
pub const APPROX_WCS_KEY: c_int = 506; /* approximate WCS keywords were calculated */

pub const NO_CLOSE_ERROR: c_int = 999; /* special value used internally to switch off */
/* the error message from ffclos and ffchdu */

/*------- following error codes are used in the grparser.c file -----------*/
pub const NGP_ERRBASE: c_int = 360; /* base chosen so not to interfere with CFITSIO */
pub const NGP_OK: c_int = 0;
pub const NGP_NO_MEMORY: c_int = NGP_ERRBASE; /* malloc failed */
pub const NGP_READ_ERR: c_int = NGP_ERRBASE + 1; /* read error from file */
pub const NGP_NUL_PTR: c_int = NGP_ERRBASE + 2; /* null pointer passed as argument */
pub const NGP_EMPTY_CURLINE: c_int = NGP_ERRBASE + 3; /* line read seems to be empty */
pub const NGP_UNREAD_QUEUE_FULL: c_int = NGP_ERRBASE + 4; /* cannot unread more then 1 line (or single line twice) */
pub const NGP_INC_NESTING: c_int = NGP_ERRBASE + 5; /* too deep include file nesting (inf. loop ?) */
pub const NGP_ERR_FOPEN: c_int = NGP_ERRBASE + 6; /* fopen() failed, cannot open file */
pub const NGP_EOF: c_int = NGP_ERRBASE + 7; /* end of file encountered */
pub const NGP_BAD_ARG: c_int = NGP_ERRBASE + 8; /* bad arguments passed */
pub const NGP_TOKEN_NOT_EXPECT: c_int = NGP_ERRBASE + 9; /* token not expected here */

#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct tcolumn {
    /* structure used to store table column information */
    pub ttype: [c_char; 70],
    /* column name = FITS TTYPEn keyword; */
    pub tbcol: LONGLONG,
    /* offset in row to first byte of each column */
    pub tdatatype: c_int,
    /* datatype code of each column */
    pub trepeat: LONGLONG,
    /* repeat count of column; number of elements */
    pub tscale: f64,
    /* FITS TSCALn linear scaling factor */
    pub tzero: f64,
    /* FITS TZEROn linear scaling zero point */
    pub tnull: LONGLONG,
    /* FITS null value for int image or binary table cols */
    pub strnull: [c_char; 20],
    /* FITS null value string for ASCII table columns */
    pub tform: [c_char; 10],
    /* FITS tform keyword value  */
    pub twidth: c_long, /* width of each ASCII table column */
}

impl Default for tcolumn {
    fn default() -> Self {
        Self {
            ttype: [0; 70],
            tbcol: Default::default(),
            tdatatype: Default::default(),
            trepeat: Default::default(),
            tscale: Default::default(),
            tzero: Default::default(),
            tnull: Default::default(),
            strnull: Default::default(),
            tform: Default::default(),
            twidth: Default::default(),
        }
    }
}

/// Structure used to store basic FITS file information.
#[repr(C)]
pub struct FITSfile {
    /// Handle returned by the file open function.
    pub filehandle: c_int,
    /// Defines which set of I/O drivers should be used.
    pub driver: c_int,
    /// Number of opened 'fitsfiles' using this structure.
    pub open_count: c_int,
    /// File name.
    pub filename: *mut c_char,
    /// Magic value used to verify that structure is valid.
    pub validcode: c_int,
    /// Flag meaning only copy the specified extension.
    pub only_one: c_int,
    /// Flag for file opened with request to ignore extended syntax.
    pub noextsyntax: c_int,
    /// Current size of the physical disk file in bytes.
    pub filesize: LONGLONG,
    /// Logical size of file, including unflushed buffers.
    pub logfilesize: LONGLONG,
    /// Is this the last HDU in the file? 0 = no, else yes.
    pub lasthdu: c_int,
    /// Current logical I/O pointer position in file.
    pub bytepos: LONGLONG,
    /// Current I/O pointer position in the physical file.
    pub io_pos: LONGLONG,
    /// Number of I/O buffer currently in use.
    pub curbuf: c_int,
    /// Current HDU number; 0 = primary array.
    pub curhdu: c_int,
    /// 0 = primary array, 1 = ASCII table, 2 = binary table.
    pub hdutype: c_int,
    /// 0 = readonly, 1 = readwrite.
    pub writemode: c_int,
    /// Highest numbered HDU known to exist in the file.
    pub maxhdu: c_int,
    /// Dynamically allocated dimension of headstart array.
    pub MAXHDU: c_int,
    /// Byte offset in file to start of each HDU.
    pub headstart: *mut LONGLONG,
    /// Byte offset in file to end of the current HDU header.
    pub headend: LONGLONG,
    /// Byte offest to where the END keyword was last written.
    pub ENDpos: LONGLONG,
    /// Byte offset in file to beginning of next keyword.
    pub nextkey: LONGLONG,
    /// Byte offset in file to start of the current data unit.
    pub datastart: LONGLONG,
    /// Dimension of image; cached for fast access.
    pub imgdim: c_int,
    /// Length of each axis; cached for fast access.
    pub imgnaxis: [LONGLONG; 99],
    /// Number of fields in the table (primary array has 2).
    pub tfield: c_int,
    /// Used by ffgcnn to record starting column number.
    pub startcol: c_int,
    /// Original number of rows (value of NAXIS2 keyword).
    pub origrows: LONGLONG,
    /// Number of rows in the table (dynamically updated).
    pub numrows: LONGLONG,
    /// Length of a table row or image size (bytes).
    pub rowlength: LONGLONG,
    /// Pointer to the table structure.
    pub tableptr: *mut tcolumn,
    /// Heap start byte relative to start of data unit.
    pub heapstart: LONGLONG,
    /// Size of the heap, in bytes.
    pub heapsize: LONGLONG,

    // The following elements are related to compressed images

    // These record the 'requested' options to be used when the image is compressed
    /// Requested image compression algorithm.
    pub request_compress_type: c_int,
    /// Requested tiling size.
    pub request_tilesize: [c_long; 6],
    /// Requested quantize level.
    pub request_quantize_level: f32,
    /// Requested quantizing method.
    pub request_quantize_method: c_int,
    /// Starting offset into the array of random dithering.
    pub request_dither_seed: c_int,
    /// Lossy compress integer image as if float image?
    pub request_lossy_int_compress: c_int,
    /// Use '1Q' rather then '1P' variable length arrays.
    pub request_huge_hdu: c_int,
    /// Requested HCOMPRESS scale factor.
    pub request_hcomp_scale: f32,
    /// Requested HCOMPRESS smooth parameter.
    pub request_hcomp_smooth: c_int,

    // These record the actual options that were used when the image was compressed
    /// Type of compression algorithm.
    pub compress_type: c_int,
    /// Size of compression tiles.
    pub tilesize: [c_long; 6],
    /// Floating point quantization level.
    pub quantize_level: f32,
    /// Floating point pixel quantization algorithm.
    pub quantize_method: c_int,
    /// Starting offset into the array of random dithering.
    pub dither_seed: c_int,

    // Other compression parameters
    /// 1 if HDU contains a compressed image, else 0.
    pub compressimg: c_int,
    /// Compression type string.
    pub zcmptype: [c_char; 12],
    /// FITS data type of image (BITPIX).
    pub zbitpix: c_int,
    /// Dimension of image.
    pub zndim: c_int,
    /// Length of each axis.
    pub znaxis: [c_long; 6],
    /// Max number of pixels in each image tile.
    pub maxtilelen: c_long,
    /// Maximum byte length of tile compressed arrays.
    pub maxelem: c_long,

    /// Column number for COMPRESSED_DATA column.
    pub cn_compressed: c_int,
    /// Column number for UNCOMPRESSED_DATA column.
    pub cn_uncompressed: c_int,
    /// Column number for GZIP2 lossless compressed data.
    pub cn_gzip_data: c_int,
    /// Column number for ZSCALE column.
    pub cn_zscale: c_int,
    /// Column number for ZZERO column.
    pub cn_zzero: c_int,
    /// Column number for the ZBLANK column.
    pub cn_zblank: c_int,

    /// Scaling value, if same for all tiles.
    pub zscale: f64,
    /// Zero pt, if same for all tiles.
    pub zzero: f64,
    /// Value of the BSCALE keyword in header.
    pub cn_bscale: f64,
    /// Value of the BZERO keyword (may be reset).
    pub cn_bzero: f64,
    /// Actual value of the BZERO keyword.
    pub cn_actual_bzero: f64,
    /// Value for null pixels, if not a column.
    pub zblank: c_int,

    /// First compression parameter: Rice pixels/block.
    pub rice_blocksize: c_int,
    /// Second compression parameter: Rice bytes/pixel.
    pub rice_bytepix: c_int,
    /// First hcompress compression parameter.
    pub hcomp_scale: f32,
    /// Second hcompress compression parameter.
    pub hcomp_smooth: c_int,

    /// Row number of the array of uncompressed tiledata.
    pub tilerow: *mut c_int,
    /// Length of the array of tile data in bytes.
    pub tiledatasize: *mut c_long,
    /// Datatype of the array of tile (TINT, TSHORT, etc).
    pub tiletype: *mut c_int,
    /// Array of uncompressed tile of data, for row *tilerow.
    pub tiledata: *mut *mut c_void,
    /// Array of optional array of null value flags.
    pub tilenullarray: *mut *mut c_void,
    /// Anynulls in the array of tile?
    pub tileanynull: *mut c_int,

    /// Pointer to FITS file I/O buffers.
    pub iobuffer: Box<[c_char; NIOBUF as usize * IOBUFLEN as usize]>,
    /// File record number of each of the buffers.
    pub bufrecnum: [c_long; 40],
    /// Has the corresponding buffer been modified?
    pub dirty: [c_int; 40],
    /// Relative age of each buffer.
    pub ageindex: [c_int; 40],
}

impl Default for FITSfile {
    fn default() -> Self {
        Self {
            filehandle: Default::default(),
            driver: Default::default(),
            open_count: Default::default(),
            filename: std::ptr::null_mut(),
            validcode: Default::default(),
            only_one: Default::default(),
            noextsyntax: Default::default(),
            filesize: Default::default(),
            logfilesize: Default::default(),
            lasthdu: Default::default(),
            bytepos: Default::default(),
            io_pos: Default::default(),
            curbuf: Default::default(),
            curhdu: Default::default(),
            hdutype: Default::default(),
            writemode: Default::default(),
            maxhdu: Default::default(),
            MAXHDU: Default::default(),
            headstart: std::ptr::null_mut(),
            headend: Default::default(),
            ENDpos: Default::default(),
            nextkey: Default::default(),
            datastart: Default::default(),
            imgdim: Default::default(),
            imgnaxis: [0; 99],
            tfield: Default::default(),
            startcol: Default::default(),
            origrows: Default::default(),
            numrows: Default::default(),
            rowlength: Default::default(),
            tableptr: std::ptr::null_mut(),
            heapstart: Default::default(),
            heapsize: Default::default(),
            request_compress_type: Default::default(),
            request_tilesize: Default::default(),
            request_quantize_level: Default::default(),
            request_quantize_method: Default::default(),
            request_dither_seed: Default::default(),
            request_lossy_int_compress: Default::default(),
            request_huge_hdu: Default::default(),
            request_hcomp_scale: Default::default(),
            request_hcomp_smooth: Default::default(),
            compress_type: Default::default(),
            tilesize: Default::default(),
            quantize_level: Default::default(),
            quantize_method: Default::default(),
            dither_seed: Default::default(),
            compressimg: Default::default(),
            zcmptype: Default::default(),
            zbitpix: Default::default(),
            zndim: Default::default(),
            znaxis: Default::default(),
            maxtilelen: Default::default(),
            maxelem: Default::default(),
            cn_compressed: Default::default(),
            cn_uncompressed: Default::default(),
            cn_gzip_data: Default::default(),
            cn_zscale: Default::default(),
            cn_zzero: Default::default(),
            cn_zblank: Default::default(),
            zscale: Default::default(),
            zzero: Default::default(),
            cn_bscale: Default::default(),
            cn_bzero: Default::default(),
            cn_actual_bzero: Default::default(),
            zblank: Default::default(),
            rice_blocksize: Default::default(),
            rice_bytepix: Default::default(),
            hcomp_scale: Default::default(),
            hcomp_smooth: Default::default(),
            tilerow: std::ptr::null_mut(),
            tiledatasize: std::ptr::null_mut(),
            tiletype: std::ptr::null_mut(),
            tiledata: std::ptr::null_mut(),
            tilenullarray: std::ptr::null_mut(),
            tileanynull: std::ptr::null_mut(),
            iobuffer: Box::new([0; (NIOBUF as usize * IOBUFLEN as usize)]),
            bufrecnum: [0; 40],
            dirty: [0; 40],
            ageindex: [0; 40],
        }
    }
}

impl FITSfile {
    pub fn new(
        driver: &fitsdriver,
        handle: c_int,
        url: &[c_char],
        caller: &[c_char],
        status: &mut c_int,
    ) -> Result<Box<Self>, c_int> {
        let mut errmsg: [c_char; FLEN_ERRMSG] = [0; FLEN_ERRMSG];

        let mut slen = strlen_safe(url) + 1;
        slen = cmp::max(slen, 32); /* reserve at least 32 chars */

        // HEAP ALLOCATION
        /* mem for file name */
        let mut f_filename = Vec::new();
        if f_filename.try_reserve_exact(slen).is_err() {
            (driver.close)(handle); /* close the file */
            int_snprintf!(
                &mut errmsg,
                FLEN_ERRMSG,
                "failed to allocate memory for filename: ({}),",
                slice_to_str!(caller),
            );
            ffpmsg_slice(&errmsg);
            ffpmsg_slice(url);
            *status = MEMORY_ALLOCATION;
            return Err(*status);
        } else {
            f_filename.resize(slen, 0);
        }

        // HEAP ALLOCATION
        /* mem for headstart array */
        let mut f_headstart = Vec::new();
        if f_headstart.try_reserve_exact(1001).is_err() {
            (driver.close)(handle); /* close the file */
            int_snprintf!(
                &mut errmsg,
                FLEN_ERRMSG,
                "failed to allocate memory for headstart array: ({}),",
                slice_to_str!(caller),
            );
            ffpmsg_slice(&errmsg);
            ffpmsg_slice(url);
            *status = MEMORY_ALLOCATION;
            return Err(*status);
        } else {
            f_headstart.resize(1001, 0);
        }

        // HEAP ALLOCATION
        /* mem for file I/O buffers */
        let mut f_iobuffer = Vec::new();
        let iosize = NIOBUF as usize * IOBUFLEN as usize;
        if f_iobuffer.try_reserve_exact(iosize).is_err() {
            (driver.close)(handle);
            int_snprintf!(
                &mut errmsg,
                FLEN_ERRMSG,
                "failed to allocate memory for iobuffer array: ({}),",
                slice_to_str!(caller),
            );
            ffpmsg_slice(&errmsg);
            ffpmsg_slice(url);
            *status = MEMORY_ALLOCATION;
            return Err(*status);
        } else {
            f_iobuffer.resize(iosize, 0);
        }

        let f_iobuffer: Box<[c_char; NIOBUF as usize * IOBUFLEN as usize]> =
            f_iobuffer.into_boxed_slice().try_into().unwrap();

        let (f_headstart, _, _) = f_headstart.into_raw_parts();

        let (f_filename, _, _) = f_filename.into_raw_parts();

        // HEAP ALLOCATION
        /* allocate FITSfile structure and initialize = 0 */
        let b = Box::try_new(FITSfile {
            iobuffer: f_iobuffer,
            headstart: f_headstart,
            filename: f_filename,
            ..Default::default()
        });
        if b.is_err() {
            (driver.close)(handle); /* close the file */
            int_snprintf!(
                &mut errmsg,
                FLEN_ERRMSG,
                "failed to allocate structure for following file: ({}),",
                slice_to_str!(caller),
            );
            ffpmsg_slice(&errmsg);
            ffpmsg_slice(url);
            *status = MEMORY_ALLOCATION;
            return Err(*status);
        }
        Ok(b.unwrap())
    }

    pub fn get_headstart_as_slice(&self) -> &[LONGLONG] {
        unsafe { std::slice::from_raw_parts(self.headstart, self.MAXHDU as usize + 1) }
    }

    pub fn get_headstart_as_mut_slice(&mut self) -> &mut [LONGLONG] {
        unsafe { std::slice::from_raw_parts_mut(self.headstart, self.MAXHDU as usize + 1) }
    }

    pub fn get_iobuffer(iobuffer: &mut *mut c_char) -> &[c_char] {
        unsafe { std::slice::from_raw_parts(*iobuffer, NIOBUF as usize * IOBUFLEN as usize) }
    }

    pub fn get_iobuffer_mut(iobuffer: &mut *mut c_char) -> &mut [c_char] {
        unsafe { std::slice::from_raw_parts_mut(*iobuffer, NIOBUF as usize * IOBUFLEN as usize) }
    }

    pub fn get_tableptr_as_slice(&self) -> &[tcolumn] {
        unsafe { std::slice::from_raw_parts(self.tableptr, self.tfield as usize) }
    }

    pub fn get_tableptr_as_slice_unlinked<'a>(
        tableptr: &*mut tcolumn,
        len: usize,
    ) -> &'a [tcolumn] {
        unsafe { std::slice::from_raw_parts(*tableptr, len) }
    }

    pub fn get_tableptr_as_mut_slice(&mut self) -> &mut [tcolumn] {
        unsafe { std::slice::from_raw_parts_mut(self.tableptr, self.tfield as usize) }
    }

    pub fn get_tile_alloc_len(&self) -> usize {
        ((((self).znaxis[0] - 1) / ((self).tilesize[0])) + 1) as usize
    }

    pub fn get_tilerow_as_slice(&self) -> &[c_int] {
        unsafe { std::slice::from_raw_parts(self.tilerow, self.get_tile_alloc_len()) }
    }

    pub fn get_tilerow_as_mut_slice(&mut self) -> &mut [c_int] {
        unsafe { std::slice::from_raw_parts_mut(self.tilerow, self.get_tile_alloc_len()) }
    }

    pub fn get_tiledatasize_as_slice(&self) -> &[c_long] {
        unsafe { std::slice::from_raw_parts(self.tiledatasize, self.get_tile_alloc_len()) }
    }

    pub fn get_tiledatasize_as_mut_slice(&mut self) -> &mut [c_long] {
        unsafe { std::slice::from_raw_parts_mut(self.tiledatasize, self.get_tile_alloc_len()) }
    }

    pub fn get_tiletype_as_slice(&self) -> &[c_int] {
        unsafe { std::slice::from_raw_parts(self.tiletype, self.get_tile_alloc_len()) }
    }

    pub fn get_tiletype_as_mut_slice(&mut self) -> &mut [c_int] {
        unsafe { std::slice::from_raw_parts_mut(self.tiletype, self.get_tile_alloc_len()) }
    }

    pub fn get_tiledata_as_slice(&self) -> &[*mut c_void] {
        unsafe { std::slice::from_raw_parts(self.tiledata, self.get_tile_alloc_len()) }
    }

    pub fn get_tiledata_as_mut_slice(&mut self) -> &mut [*mut c_void] {
        unsafe { std::slice::from_raw_parts_mut(self.tiledata, self.get_tile_alloc_len()) }
    }

    pub fn get_tilenullarray_as_slice(&self) -> &[*mut c_void] {
        unsafe { std::slice::from_raw_parts(self.tilenullarray, self.get_tile_alloc_len()) }
    }

    pub fn get_tilenullarray_as_mut_slice(&mut self) -> &mut [*mut c_void] {
        unsafe { std::slice::from_raw_parts_mut(self.tilenullarray, self.get_tile_alloc_len()) }
    }

    /*

    pub fn get_tiledata_as_slice(&self) -> &[&[c_void]] {
        let len = self.get_tile_alloc_len();
        let mut v = Vec::with_capacity(len);
        for i in 0..len {
            let start = i * self.tiledatasize[i] as usize;
            let end = start + self.tiledatasize[i] as usize;
            v.push(&self.tiledata[start..end]);
        }
        v
    }


    pub fn get_tiledata_as_mut_slice(&mut self) -> &mut [&mut [c_void]] {
        let len = self.get_tile_alloc_len();
        let mut v = Vec::with_capacity(len);
        for i in 0..len {
            let start = i * self.tiledatasize[i] as usize;
            let end = start + self.tiledatasize[i] as usize;
            v.push(&mut self.tiledata[start..end]);
        }
        v
    }

    pub fn get_tilenullarray_as_slice(&self) -> &[&[c_void]] {
        let len = self.get_tile_alloc_len();
        let mut v = Vec::with_capacity(len);
        for i in 0..len {
            let start = i * self.tiledatasize[i] as usize;
            let end = start + self.tiledatasize[i] as usize;
            v.push(&self.tilenullarray[start..end]);
        }
        v
    }

    pub fn get_tilenullarray_as_mut_slice(&mut self) -> &mut [&mut [c_void]] {
        let len = self.get_tile_alloc_len();
        let mut v = Vec::with_capacity(len);
        for i in 0..len {
            let start = i * self.tiledatasize[i] as usize;
            let end = start + self.tiledatasize[i] as usize;
            v.push(&mut self.tilenullarray[start..end]);
        }
        v
    }
    */

    pub fn get_tileanynull_as_slice(&self) -> &[c_int] {
        unsafe { std::slice::from_raw_parts(self.tileanynull, self.get_tile_alloc_len()) }
    }

    pub fn get_tileanynull_as_mut_slice(&mut self) -> &mut [c_int] {
        unsafe { std::slice::from_raw_parts_mut(self.tileanynull, self.get_tile_alloc_len()) }
    }
}

impl Drop for FITSfile {
    fn drop(&mut self) {
        unsafe {
            if !self.filename.is_null() {
                // Can't use CString::from_raw here because there is a minimum of 32 bytes reserved
                let c = CStr::from_ptr(self.filename);
                let mut slen = c.to_bytes_with_nul().len();
                slen = cmp::max(slen, 32);

                let _ = Vec::from_raw_parts(self.filename, slen, slen); //HEAP DEALLOCATION
            }

            if !self.headstart.is_null() {
                let _ = Vec::from_raw_parts(
                    self.headstart,
                    (self.MAXHDU as usize) + 1,
                    (self.MAXHDU as usize) + 1,
                ); //HEAP DEALLOCATION
            }

            if !self.tableptr.is_null() {
                let _ =
                    Vec::from_raw_parts(self.tableptr, self.tfield as usize, self.tfield as usize);
                //HEAP DEALLOCATION
            }

            if !self.tilerow.is_null() {
                let _ = Vec::from_raw_parts(
                    self.tilerow,
                    self.get_tile_alloc_len(),
                    self.get_tile_alloc_len(),
                ); //HEAP DEALLOCATION
            }

            if !self.tiledatasize.is_null() {
                let _ = Vec::from_raw_parts(
                    self.tiledatasize,
                    self.get_tile_alloc_len(),
                    self.get_tile_alloc_len(),
                ); //HEAP DEALLOCATION
            }

            if !self.tiletype.is_null() {
                let _ = Vec::from_raw_parts(
                    self.tiletype,
                    self.get_tile_alloc_len(),
                    self.get_tile_alloc_len(),
                ); //HEAP DEALLOCATION
            }

            if !self.tiledata.is_null() {
                let _ = Vec::from_raw_parts(
                    self.tiledata,
                    self.get_tile_alloc_len(),
                    self.get_tile_alloc_len(),
                ); //HEAP DEALLOCATION
            }

            if !self.tilenullarray.is_null() {
                let _ = Vec::from_raw_parts(
                    self.tilenullarray,
                    self.get_tile_alloc_len(),
                    self.get_tile_alloc_len(),
                ); //HEAP DEALLOCATION
            }

            if !self.tileanynull.is_null() {
                let _ = Vec::from_raw_parts(
                    self.tileanynull,
                    self.get_tile_alloc_len(),
                    self.get_tile_alloc_len(),
                ); //HEAP DEALLOCATION
            }
        }
    }
}

#[repr(C)]
/// Structure used to store basic HDU information
pub struct fitsfile {
    /// HDU position in file; 0 = first HDU
    pub HDUposition: c_int,

    /// Pointer to FITS file structure
    pub Fptr: Box<FITSfile>,
}

#[repr(C)]
pub struct iteratorCol {
    /* structure for the iterator function column information */
    /* elements required as input to fits_iterate_data: */
    pub fptr: *mut fitsfile,
    /* pointer to the HDU containing the column */
    pub colnum: c_int,
    /* column number in the table (use name if < 1) */
    pub colname: [c_char; 70],
    /* name (= TTYPEn value) of the column (optional) */
    pub datatype: c_int,
    /* output datatype (converted if necessary  */
    pub iotype: c_int,
    /* = InputCol, InputOutputCol, or OutputCol */
    /* output elements that may be useful for the work function: */
    pub array: *mut c_void,
    /* pointer to the array (and the null value) */
    pub repeat: c_long,
    /* binary table vector repeat value */
    pub tlmin: c_long,
    /* legal minimum data value */
    pub tlmax: c_long,
    /* legal maximum data value */
    pub tunit: [c_char; 70],
    /* physical unit string */
    pub tdisp: [c_char; 70], /* suggested display format */
}

/*=============================================================================
*
*       The following wtbarr typedef is used in the fits_read_wcstab() routine,
*       which is intended for use with the WCSLIB library written by Mark
*       Calabretta, http://www.atnf.csiro.au/~mcalabre/index.html
*
*       In order to maintain WCSLIB and CFITSIO as independent libraries it
*       was not permissible for any CFITSIO library code to include WCSLIB
*       header files, or vice versa.  However, the CFITSIO function
*       fits_read_wcstab() accepts an array of structs defined by wcs.h within
*       WCSLIB.  The problem then was to define this struct within fitsio.h
*       without including wcs.h, especially noting that wcs.h will often (but
*       not always) be included together with fitsio.h in an applications
*       program that uses fits_read_wcstab().
*
*       Of the various possibilities, the solution adopted was for WCSLIB to
*       define "struct wtbarr" while fitsio.h defines "typedef wtbarr", a
*       untagged struct with identical members.  This allows both wcs.h and
*       fitsio.h to define a wtbarr data type without conflict by virtue of
*       the fact that structure tags and typedef names share different
*       namespaces in C. Therefore, declarations within WCSLIB look like
*
*          struct wtbarr *w;
*
*       while within CFITSIO they are simply
*
*          wtbarr *w;
*
*       but as suggested by the commonality of the names, these are really the
*       same aggregate data type.  However, in passing a (struct wtbarr *) to
*       fits_read_wcstab() a cast to (wtbarr *) is formally required.
*===========================================================================*/
#[repr(C)]
pub struct wtbarr {
    pub i: c_int,             /* Image axis number.                       */
    pub m: c_int,             /* Array axis number for index vectors.     */
    pub kind: c_int,          /* Array type, 'c' (coord) or 'i' (index).  */
    pub extnam: [c_char; 72], /* EXTNAME of binary table extension.       */
    pub extver: c_int,        /* EXTVER  of binary table extension.       */
    pub extlev: c_int,        /* EXTLEV  of binary table extension.       */
    pub ttype: [c_char; 72],  /* TTYPEn of column containing the array.   */
    pub row: c_long,          /* Table row number.                        */
    pub ndim: c_int,          /* Expected array dimensionality.           */
    pub dimlen: *mut c_int,   /* Where to write the array axis lengths.   */
    pub arrayp: *mut *mut f64, /* Where to write the address of the array  */
                              /* allocated to store the array.            */
}

/// We need to implement the Drop trait for wtbarr to ensure that the memory
/// allocated for the dimlen and arrayp fields is properly freed when the
/// wtbarr struct is dropped. This is necessary because these fields are
/// pointers to dynamically allocated memory, and Rust's default behavior
/// does not automatically free memory allocated with malloc or similar
/// functions.
impl Drop for wtbarr {
    fn drop(&mut self) {
        unsafe {
            if !self.dimlen.is_null() {
                let _ = Vec::from_raw_parts(self.dimlen, self.ndim as usize, self.ndim as usize);
            }
            if !self.arrayp.is_null() {
                let _ = Vec::from_raw_parts(*self.arrayp, self.ndim as usize, self.ndim as usize);
            }
        }
    }
}

#[repr(C)]
pub struct PixelFilter {
    /* input(s) */
    pub count: c_int,
    pub path: *mut *mut c_char,
    pub tag: *mut *mut c_char,
    pub ifptr: *mut *mut fitsfile,
    pub expression: *mut c_char,
    /* output control */
    pub bitpix: c_int,
    pub blank: c_long,
    pub ofptr: *mut fitsfile,
    pub keyword: [c_char; FLEN_KEYWORD],
    pub comment: [c_char; FLEN_COMMENT],
}
