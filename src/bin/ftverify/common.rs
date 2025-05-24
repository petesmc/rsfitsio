use rsfitsio::c_types::{c_char, c_int, c_ulong};
use rsfitsio::fitsio::{FLEN_FILENAME, FLEN_KEYWORD, FLEN_VALUE, LONGLONG};

pub(crate) const MAXERRORS: usize = 200;
pub(crate) const MAXWRNS: usize = 200;

static errmes: [c_char; 256] = [0; 256]; /* error message buffer */
static comm: [c_char; FLEN_FILENAME + 6] = [0; FLEN_FILENAME + 6]; /* comment buffer */

/********************************
*				*
*       Keywords 		*
*				*
********************************/

enum kwdtyp {
    STR_KEY, /* string   key */
    LOG_KEY, /* Logical key */
    INT_KEY, /* Integer key */
    FLT_KEY, /* Float key   */
    CMI_KEY, /* Complex integer key */
    CMF_KEY, /* Complex float key */
    COM_KEY, /* history, comment, "  ", and end */
    UNKNOWN, /* Unknown types */
}

/* error number masks of  the keyword test */
const BAD_STR: c_ulong = 0x0001;
const NO_TRAIL_QUOTE: c_ulong = 0x0002;
const BAD_NUM: c_ulong = 0x0004;
const LOWCASE_EXPO: c_ulong = 0x0008;
const NO_TRAIL_PAREN: c_ulong = 0x0010;
const NO_COMMA: c_ulong = 0x0020;
const TOO_MANY_COMMA: c_ulong = 0x0040;
const BAD_REAL: c_ulong = 0x0080;
const BAD_IMG: c_ulong = 0x0100;
const BAD_LOGICAL: c_ulong = 0x0200;
const NO_START_SLASH: c_ulong = 0x0400;
const BAD_COMMENT: c_ulong = 0x0800;
const UNKNOWN_TYPE: c_ulong = 0x1000;

/* Number of possible WCS descriptions to check.*/
/* 1 for the primary + 26 for [A-Z] suffix. */
const NWCSDESCR: usize = 27;

/* keyword structure */
struct FitsKey {
    kname: [char; FLEN_KEYWORD], /* fits keyword name */
    ktype: kwdtyp,               /* fits keyword type */
    kvalue: [char; FLEN_VALUE],  /* fits keyword name */
    kindex: c_int,               /* position at the header */
    goodkey: c_int,              /* good keyword flag (=1 good)*/
}

/********************************
*				*
*       Headers  		*
*				*
********************************/
struct FitsHdu {
    hdutype: c_int,                /* hdutype */
    hdunum: c_int,                 /* hdunum  */
    isgroup: c_int,                /* random group flag */
    istilecompressed: c_int,       /* tile compressed image */
    gcount: c_int,                 /* gcount  */
    pcount: LONGLONG,              /* pcount  */
    bitpix: c_int,                 /* pix number */
    naxis: c_int,                  /* number of the axis,used for image array*/
    naxes: *mut LONGLONG,          /* dimension of each axis,used for image array*/
    ncols: c_int,                  /* number of the columns, used for image only*/
    extname: [c_char; FLEN_VALUE], /* EXTENSION NAME */
    extver: c_int,                 /* extension version */
    datamax: *mut *mut c_char,     /* strings for the maximum of the data in a column */
    datamin: *mut *mut c_char,     /* strings for the minimum of the data in a column */
    tnull: *mut *mut c_char,       /* number of NULL values */
    nkeys: c_int,                  /* number of keys */
    tkeys: c_int,                  /* total of the keys tested*/
    heap: c_int,                   /* heap */
    kwds: *mut *mut FitsKey,       /* keywords list starting from the
                                   last NAXISn keyword. The array
                                   is sorted in the ascending alphabetical
                                   order of keyword names. The last keyword END
                                   and commentary keywords are  excluded.
                                   The total number of element, tkey, is
                                   nkeys - 4 - naxis - ncomm. */
    use_longstr: c_int, /* flag indicates that the long string
                        convention is used */
}

struct ColName {
    name: *mut c_char, /* column name */
    index: c_int,      /* column index */
}

struct HduName {
    hdutype: c_int,                /* hdutype */
    hdunum: c_int,                 /* hdunum  */
    extname: [c_char; FLEN_VALUE], /* extension name, used for extension*/
    extver: c_int,                 /* extension version, used for extension */
    errnum: c_int,                 /* number of errors in this hdu */
    wrnno: c_int,                  /* number of warnning in this hdu */
}
