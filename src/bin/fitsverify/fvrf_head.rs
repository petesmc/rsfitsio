/* Transpiled from cfitsio/utilities/fvrf_head.c

   The `char **' working arrays (cards, tmpkwds, ttype, tform, tunit) are
   file-statics in the C; here they are thread-locals holding owned copies.
   `temp' / `ptemp' were only ever used to hand a pattern to key_match(), so
   the pattern is passed directly instead.  */

use std::cell::{Cell, RefCell};
use std::cmp::Ordering;

use bytemuck::cast_slice;
use rsfitsio::aliases::rust_api::{
    fits_close_file, fits_decode_tdim, fits_get_coltype, fits_get_hdrspace, fits_get_num_cols,
    fits_get_num_hdus, fits_movabs_hdu, fits_null_check, fits_open_diskfile, fits_read_key,
    fits_read_record, fits_str2time,
};
use rsfitsio::c_types::{c_char, c_double, c_int, c_long};
use rsfitsio::cs;
use rsfitsio::fitscore::ffchfl_safe;
use rsfitsio::fitsio::{
    ASCII_TBL, BINARY_TBL, BYTE_IMG, DOUBLE_IMG, FLEN_CARD, FLEN_FILENAME, FLEN_KEYWORD,
    FLEN_VALUE, FLOAT_IMG, IMAGE_HDU, LONGLONG, LONGLONG_IMG, LONG_IMG, READONLY, SHORT_IMG,
    ULONG_IMG, USHORT_IMG, fitsfile,
};
use rsfitsio::{KeywordDatatypeMut, NullValue};

use crate::common::*;
use crate::fvrf_data::test_data;
use crate::fvrf_file::*;
use crate::fvrf_key::*;
use crate::fvrf_misc::*;
use crate::{scat, spf};

thread_local! {
    /* array to store the keywords  */
    static CARDS: RefCell<Vec<[c_char; FLEN_CARD]>> = const { RefCell::new(Vec::new()) };
    /* String array holding the keyword name.  It is sorted in alphabetical
       ascending order and does not include the keywords before the first
       non-reserved keyword and END keyword. */
    static TMPKWDS: RefCell<Vec<[c_char; FLEN_KEYWORD]>> = const { RefCell::new(Vec::new()) };

    static TTYPE: RefCell<Vec<Vec<c_char>>> = const { RefCell::new(Vec::new()) };
    static TFORM: RefCell<Vec<Vec<c_char>>> = const { RefCell::new(Vec::new()) };
    static TUNIT: RefCell<Vec<Vec<c_char>>> = const { RefCell::new(Vec::new()) };
}

/* total number of the keywords */
thread_local! { static NCARDS: Cell<c_int> = const { Cell::new(0) }; }
thread_local! { static CURHDU: Cell<c_int> = const { Cell::new(0) }; } /* current HDU index */
thread_local! { static CURTYPE: Cell<c_int> = const { Cell::new(0) }; } /* current HDU type  */
/* print_title's `static int oldhdu' */
thread_local! { static OLDHDU: Cell<c_int> = const { Cell::new(0) }; }

fn cards_len() -> usize {
    CARDS.with(|c| c.borrow().len())
}

fn card_at(i: usize) -> [c_char; FLEN_CARD] {
    CARDS.with(|c| c.borrow()[i])
}

fn tmpkwds_get() -> Vec<[c_char; FLEN_KEYWORD]> {
    TMPKWDS.with(|t| t.borrow().clone())
}

fn tform_at(i: usize) -> Vec<c_char> {
    TFORM.with(|t| {
        let t = t.borrow();
        t.get(i).cloned().unwrap_or_else(|| vec![0])
    })
}

fn ttype_at(i: usize) -> Vec<c_char> {
    TTYPE.with(|t| {
        let t = t.borrow();
        t.get(i).cloned().unwrap_or_else(|| vec![0])
    })
}

fn tunit_at(i: usize) -> Vec<c_char> {
    TUNIT.with(|t| {
        let t = t.borrow();
        t.get(i).cloned().unwrap_or_else(|| vec![0])
    })
}

/******************************************************************************
* Function
*      verify_fits
*
* DESCRIPTION:
*      Verify individual fits file.
*
*******************************************************************************/
/* routine to verify individual fitsfile */
/* NB: like the C, this strips trailing whitespace in the caller's buffer -- 
   ftverify_work() relies on that when it echoes the name fgets() read. */
pub(crate) fn verify_fits(infile: &mut [c_char], out: Out) -> c_int {
    let rootnam: [c_char; FLEN_FILENAME] = [0; FLEN_FILENAME]; /* Input Fits file root name */
    let mut infits: Option<Box<fitsfile>> = None; /* input fits file pointer */
    let mut fitshdu = FitsHdu::default(); /* hdu information */
    let mut hdutype: c_int;
    let mut status: c_int = 0;
    let mut i: c_int;
    let len: usize;
    let mut p: usize;
    let pfile: usize;
    let mut xtension: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut comm: [c_char; COMM_LEN] = [0; COMM_LEN];

    /* take out the leading and trailing space and skip the empty line*/
    p = 0;
    while isspace_c(infile[p]) {
        p += 1;
    }
    len = cstrlen(&infile[p..]);
    pfile = p;
    if len > 0 {
        let mut q = p + len - 1;
        let mut ii = len as isize - 1;
        while ii >= 0 && isspace_c(infile[q]) {
            infile[q] = 0;
            if q == 0 {
                break;
            }
            q -= 1;
            ii -= 1;
        }
    }
    if cstrlen(&infile[pfile..]) == 0 {
        return status;
    }

    /* #ifndef WEBTOOL */
    wrtout_str(out, " ");
    spf!(comm; "File: ", CS(&infile[pfile..]));
    wrtout(out, &comm);

    set_totalhdu(0);

    /* #else  (STANDALONE) */
    if fits_open_diskfile(&mut infits, &infile[pfile..], READONLY, &mut status) != 0 {
        wrtserr_str(out, "", &mut status, 2);
        leave_early(out);
        status = 1;
        return status;
    }

    let mut infits = infits.unwrap();

    /* get the total hdus */
    let mut nhdu: c_int = 0;
    if fits_get_num_hdus(&mut infits, &mut nhdu, &mut status) != 0 {
        set_totalhdu(nhdu);
        wrtserr_str(out, "", &mut status, 2);
        leave_early(out);
        status = 1;
        return status;
    }
    set_totalhdu(nhdu);

    /* initialize the report */
    init_report(out, &rootnam);
    /*------------------  Hdu Loop --------------------------------*/
    i = 1;
    while i <= totalhdu() {
        /* move to the right hdu and do the CFITSIO test */
        hdutype = -1;
        if fits_movabs_hdu(&mut infits, i, Some(&mut hdutype), &mut status) != 0 {
            print_title(out, i, hdutype);
            wrtferr_str(out, "", &mut status, 2);
            set_hdubasic(i, hdutype);
            break;
        }

        if i != 1 && hdutype == IMAGE_HDU {
            /* test if this is a tile compressed image in a binary table */
            fits_read_key(
                &mut infits,
                KeywordDatatypeMut::TSTRING(&mut xtension),
                cs!(c"XTENSION"),
                None,
                &mut status,
            );
            if ceq(&xtension, b"BINTABLE") {
                print_title(out, i, BINARY_TBL);
            } else {
                print_title(out, i, hdutype);
            }
        } else {
            print_title(out, i, hdutype);
        }

        init_hdu(&mut infits, out, i, hdutype, &mut fitshdu); /* initialize fitshdu  */

        test_hdu(&mut infits, out, &mut fitshdu); /* test hdu header */

        if testdata() != 0 {
            test_data(&mut infits, out, &mut fitshdu);
        }

        close_err(out); /* end of error report */

        if prhead() != 0 {
            print_header(out);
        }
        if prstat() != 0 {
            print_summary(&mut infits, out, &fitshdu);
        }
        close_hdu(&mut fitshdu); /* clear the fitshdu  */

        i += 1;
    }
    /* test the end of file  */
    test_end(&mut infits, out);

    /*------------------ Closing  --------------------------------*/
    /* closing the report*/
    close_report(out);

    /* close the input fitsfile  */
    fits_close_file(infits, &mut status);

    status
}

pub(crate) fn leave_early(out: Out) {
    let mut comm: [c_char; COMM_LEN] = [0; COMM_LEN];
    spf!(comm; "**** Abort Verification: Fatal Error. ****");
    wrtout(out, &comm);

    /* write the total number of errors and warnings to parfile*/
    crate::main_ftverify::update_parfile(1, 0);
}

fn close_err(out: Out) {
    let mut merr = 0;
    let mut mwrn = 0;
    num_err_wrn(&mut merr, &mut mwrn);
    if merr != 0 || mwrn != 0 {
        wrtout_str(out, " ");
    }
}

/*************************************************************
*
*      init_hdu
*
*   Initialize the FitsHdu, HduName and ttype, tform, tunit if
* the hdu is a table.
*
*
*************************************************************/
fn init_hdu(
    infits: &mut fitsfile, /* input fits file   */
    out: Out,              /* output ascii file */
    hdunum: c_int,         /* hdu index 	     */
    hdutype: c_int,        /* hdutype	     */
    hduptr: &mut FitsHdu,
) {
    let mut morekeys: c_int;
    let mut i: usize;
    let mut j: usize;
    let k: usize;
    let m: usize;
    let n: usize;
    let mut status: c_int = 0;
    let mut p: usize;
    let numusrkey: usize;
    let mut lv: LONGLONG;
    let mut lu: LONGLONG = 0;
    let mut errmes: [c_char; ERRMES_LEN] = [0; ERRMES_LEN];
    let mut comm: [c_char; COMM_LEN] = [0; COMM_LEN];

    let mut tmpkey = FitsKey::default();

    hduptr.hdunum = hdunum;
    hduptr.hdutype = hdutype;

    /* curhdu and curtype are shared with print_title */
    CURHDU.with(|c| c.set(hdunum)); /* set the current hdu number */
    CURTYPE.with(|c| c.set(hdutype)); /* set the current hdu number */

    /* check the null character in the header.(only the first one will
       be recorded */
    let lvi = fits_null_check(infits, &mut status) as LONGLONG;
    if lvi > 0 {
        let mm = (lvi - 1) / 80 + 1;
        let nn = lvi - (mm - 1) * 80;
        spf!(errmes; "Byte #", nn as i64, " in Card#", mm as i64, " is a null(\\0).");
        wrterr(out, &errmes, 1);
        status = 0;
    } else if status != 0 {
        wrtserr_str(out, "", &mut status, 1);
        status = 0;
    }

    /* get the total number of keywords */
    hduptr.nkeys = 0;
    morekeys = 0;
    if fits_get_hdrspace(
        infits,
        Some(&mut hduptr.nkeys),
        Some(&mut morekeys),
        &mut status,
    ) != 0
    {
        wrtferr_str(out, "", &mut status, 1);
    }
    hduptr.nkeys += 1; /* include END keyword */

    /* read all the keywords  */
    let ncards = hduptr.nkeys as usize;
    NCARDS.with(|c| c.set(ncards as c_int));
    CARDS.with(|c| {
        let mut c = c.borrow_mut();
        c.clear();
        c.resize(ncards, [0; FLEN_CARD]);
    });

    for i in 1..=ncards {
        let mut card: [c_char; FLEN_CARD] = [0; FLEN_CARD];
        if fits_read_record(infits, i as c_int, Some(&mut card), &mut status) != 0 {
            wrtferr_str(out, "", &mut status, 1);
        }
        CARDS.with(|c| c.borrow_mut()[i - 1] = card);
    }

    /* if there were blank cards prior to the END card, then
       make a fake END card, because CFITSIO blocks us from reading
       the real END card */

    if !cstarts(&card_at(ncards - 1), b"END     ") {
        let mut c: [c_char; FLEN_CARD] = [0; FLEN_CARD];
        strcpy_c(&mut c, cs!(c"END     "));
        CARDS.with(|cc| cc.borrow_mut()[ncards - 1] = c);
    }

    /* Parse the XTENSION/SIMPLEX  keyword */
    let mut card0 = card_at(0);
    fits_parse_card(
        out,
        1,
        &mut card0,
        &mut tmpkey.kname,
        &mut tmpkey.ktype,
        &mut tmpkey.kvalue,
        &mut comm,
    );
    if tmpkey.kvalue[0] as u8 == b' ' {
        spf!(errmes;
            "Keyword #1, ", CS(&tmpkey.kname), " \"", CS(&tmpkey.kvalue),
            "\" should not have leading space.");
        wrterr(out, &errmes, 1);
    }
    if hdunum == 1 {
        /* SIMPLE should be logical T */
        if !ceq(&tmpkey.kname, b"SIMPLE") {
            wrterr_str(out, "The 1st keyword of a primary array is not SIMPLE.", 1);
        }
        if check_log(&tmpkey, out) == 0 || !ceq(&tmpkey.kvalue, b"T") {
            wrtwrn_str(
                out,
                "SIMPLE != T indicates file may not conform to the FITS Standard.",
                0,
            );
        }

        check_fixed_log(&card_at(0), out);
    } else {
        if !ceq(&tmpkey.kname, b"XTENSION") {
            wrterr_str(out, "The 1st keyword of a extension is not XTENSION.", 1);
        }
        check_str(&tmpkey, out);

        check_fixed_str(&card_at(0), out);

        /* Get the original string */
        let c0 = card_at(0);
        p = 10;
        while c0[p] as u8 == b' ' {
            p += 1;
        }
        p += 1; /* skip the  quote */
        if !cstarts(&c0[p..], b"TABLE   ")
            && !cstarts(&c0[p..], b"BINTABLE")
            && !cstarts(&c0[p..], b"A3DTABLE")
            && !cstarts(&c0[p..], b"IUEIMAGE")
            && !cstarts(&c0[p..], b"FOREIGN ")
            && !cstarts(&c0[p..], b"DUMP    ")
            && !cstarts(&c0[p..], b"IMAGE   ")
        {
            spf!(errmes;
                "Unregistered XTENSION value \"", CSW(&c0[p..], 0, Some(8)), "\".");
            wrterr(out, &errmes, 1);
        } else if c0[p + 8] as u8 != b'\'' {
            spf!(errmes;
                "Extra '", CHR(c0[p + 8]), "' follows the XTENSION value \"",
                CSW(&c0[p..], 0, Some(8)), "\".");
            wrterr(out, &errmes, 1);
        }

        /* test if this is a tile compressed image, stored in a binary table */
        /* If so then test the extension as binary table rather than an image */

        if cstarts(&c0[p..], b"BINTABLE") && hduptr.hdutype == IMAGE_HDU {
            hduptr.hdutype = BINARY_TBL;
            hduptr.istilecompressed = 1;
        } else {
            hduptr.istilecompressed = 0;
        }
    }

    /* read the BITPIX keywords */
    if fits_read_key(
        infits,
        KeywordDatatypeMut::TINT(&mut hduptr.bitpix),
        cs!(c"BITPIX"),
        None,
        &mut status,
    ) != 0
    {
        wrtferr_str(out, "", &mut status, 2);
    }
    check_fixed_int(&card_at(1), out);

    /* Read and Parse the NAXIS */
    hduptr.naxis = 0;
    if fits_read_key(
        infits,
        KeywordDatatypeMut::TINT(&mut hduptr.naxis),
        cs!(c"NAXIS"),
        None,
        &mut status,
    ) != 0
    {
        wrtferr_str(out, "", &mut status, 2);
    }
    check_fixed_int(&card_at(2), out);

    hduptr.naxes = if hduptr.naxis != 0 {
        vec![-1; hduptr.naxis as usize]
    } else {
        Vec::new()
    };

    /* Parse the keywords NAXISn */
    j = 3;
    while j < 3 + hduptr.naxis as usize {
        let mut cardj = card_at(j);
        fits_parse_card(
            out,
            1 + j as c_int,
            &mut cardj,
            &mut tmpkey.kname,
            &mut tmpkey.ktype,
            &mut tmpkey.kvalue,
            &mut comm,
        );
        p = 5;
        if !isdigit_c(tmpkey.kname[p]) {
            j += 1;
            continue;
        }
        if check_int(&tmpkey, out) != 0 {
            lu = strtoll_c(&tmpkey.kvalue);
        }
        lv = strtol_c(&tmpkey.kname[p..]).0 as LONGLONG;
        if lv > hduptr.naxis as LONGLONG && lv <= 0 {
            spf!(errmes;
                "Keyword #", tmpkey.kindex, ", ", CS(&tmpkey.kname),
                " is not allowed (with n > NAXIS = ", hduptr.naxis, ").");
            wrterr(out, &errmes, 1);
        } else if lv >= 1 && (lv as usize) <= hduptr.naxes.len() {
            if hduptr.naxes[(lv - 1) as usize] == -1 {
                hduptr.naxes[(lv - 1) as usize] = lu;
            } else {
                spf!(errmes;
                    "Keyword #", tmpkey.kindex, ", ", CS(&tmpkey.kname), " is duplicated.");
                wrterr(out, &errmes, 1);
            }
        }

        check_fixed_int(&card_at(j), out);
        j += 1;
    }

    /* check all the NAXISn are there */
    for j in 0..hduptr.naxis as usize {
        if hduptr.naxes[j] == -1 {
            spf!(errmes; "Keyword NAXIS", j + 1, " is not present or is out of order.");
            wrterr(out, &errmes, 2);
        }
    }

    /* get the column number */
    hduptr.ncols = 1;
    if hduptr.hdutype == ASCII_TBL || hduptr.hdutype == BINARY_TBL {
        /* get the total number of columns  */
        if fits_get_num_cols(infits, &mut hduptr.ncols, &mut status) != 0 {
            wrtferr_str(out, "", &mut status, 2);
        }
    }

    /* parse the keywords after NAXISn and prepare the array for
       sorting. We only check the keywords after the NAXISn */
    n = (hduptr.nkeys - 4 - hduptr.naxis).max(0) as usize; /* excluding the SIMPLE/XTENSION,
                                                           BITPIX, NAXIS, NAXISn
                                                           and END */
    hduptr.kwds = vec![FitsKey::default(); n];
    k = 3 + hduptr.naxis as usize; /* index of first keyword following NAXISn. */
    m = (hduptr.nkeys - 1) as usize; /* last key  */
    i = 0;
    hduptr.use_longstr = 0;
    j = k;
    while j < m {
        let mut cardj = card_at(j);
        hduptr.kwds[i].kindex = j as c_int + 1; /* record number */
        hduptr.kwds[i].goodkey = 1;
        let (mut kn, mut kt, mut kv) = (
            hduptr.kwds[i].kname,
            hduptr.kwds[i].ktype,
            hduptr.kwds[i].kvalue,
        );
        if fits_parse_card(out, 1 + j as c_int, &mut cardj, &mut kn, &mut kt, &mut kv, &mut comm)
            != 0
        {
            hduptr.kwds[i].goodkey = 0;
        }
        hduptr.kwds[i].kname = kn;
        hduptr.kwds[i].ktype = kt;
        hduptr.kwds[i].kvalue = kv;

        if hduptr.kwds[i].ktype == kwdtyp::UNKNOWN && hduptr.kwds[i].kvalue[0] == 0 {
            spf!(errmes;
                "Keyword #", j + 1, ", ", CS(&hduptr.kwds[i].kname), " has a null value.");
            wrtwrn(out, &errmes, 0);
        }

        /* only count the non-commentary keywords */
        if ceq(&hduptr.kwds[i].kname, b"CONTINUE") {
            hduptr.use_longstr = 1;
        }
        if !ceq(&hduptr.kwds[i].kname, b"COMMENT")
            && !ceq(&hduptr.kwds[i].kname, b"HISTORY")
            && (!ceq(&hduptr.kwds[i].kname, b"HIERARCH") || testhierarch() != 0)
            && !ceq(&hduptr.kwds[i].kname, b"CONTINUE")
            && !ceq(&hduptr.kwds[i].kname, b"")
        {
            i += 1;
        }
        j += 1;
    }
    numusrkey = i;
    hduptr.tkeys = i as c_int;

    /* parse the END key */
    let mut cardend = card_at(hduptr.nkeys as usize - 1);
    fits_parse_card(
        out,
        m as c_int + 1,
        &mut cardend,
        &mut tmpkey.kname,
        &mut tmpkey.ktype,
        &mut tmpkey.kvalue,
        &mut comm,
    );

    /* sort the keyword in the ascending order of kname field*/
    /* glibc's qsort is a stable merge sort here, so use a stable sort. */
    hduptr.kwds[..numusrkey].sort_by(compkey);

    /* store addresses of sorted keyword names in a working array */
    TMPKWDS.with(|t| {
        let mut t = t.borrow_mut();
        t.clear();
        for i in 0..numusrkey {
            t.push(hduptr.kwds[i].kname);
        }
    });

    /* Initialize  the PCOUNT, GCOUNT and heap values */
    hduptr.pcount = -99;
    hduptr.gcount = -99;
    hduptr.heap = -99;

    /* set the random group flag (will be determined later) */
    hduptr.isgroup = 0;

    /* allocate memory for datamax and datamin (will determined later)*/
    hduptr.datamax = Vec::new();
    hduptr.datamin = Vec::new();
    hduptr.tnull = Vec::new();
    if hduptr.ncols > 0 {
        for _ in 0..hduptr.ncols {
            hduptr.datamax.push(vec![0; 13]);
            hduptr.datamin.push(vec![0; 13]);
            hduptr.tnull.push(vec![0; 12]);
        }
    }

    /* initialize  the extension  name and version */
    hduptr.extname[0] = 0;
    hduptr.extver = -999;
}

/*************************************************************
*
*      test_hdu
*
*   Test the  HDU header
*    This includes many tests of WCS keywords
*
*************************************************************/
fn test_hdu(
    infits: &mut fitsfile, /* input fits file   */
    out: Out,              /* output ascii file */
    hduptr: &mut FitsHdu,
) {
    let mut status: c_int = 0;
    let numusrkey: c_int;
    let hdunum: c_int;
    let mut pname: [Option<[c_char; FLEN_KEYWORD]>; NWCSDESCR] = [None; NWCSDESCR];
    let mut p: usize;
    let mut j: c_int;
    let mut k: c_int = 0;
    let mut m: c_int;
    let mut n: c_int = 0;
    let mut wcsaxesmax: c_int = 0;
    let mut taxes: c_int;
    let mut iaxis: c_int = 0;
    let mut ialt: c_int = 0;
    let mut wcsaxesExists: c_int = 0;
    let mut wcsaxesvalue: c_int = 0;
    let mut wcsaxespos: [c_int; NWCSDESCR] = [0; NWCSDESCR];
    let mut wcskeypos: [c_int; NWCSDESCR] = [0; NWCSDESCR];
    let mut crota2_exists: c_int = 0;
    let mut matrix_exists: [c_int; 2] = [0, 0];
    let mut dvalue: c_double;
    let mut primary_naxis: c_int = 0;
    let mut errmes: [c_char; ERRMES_LEN] = [0; ERRMES_LEN];

    /* floating WCS keywords  */
    let cfltkeys: [&std::ffi::CStr; 7] = [
        c"CRPIX", c"CRVAL", c"CDELT", c"CROTA", c"CRDER", c"CSYER", c"PV",
    ];
    let ncfltkeys = 7;
    let mut keynum: [c_int; 7] = [0; 7];
    let mut nmax: c_int = 0;

    /* floating non-indexed WCS keywords  */
    let cfltnkeys: [&std::ffi::CStr; 11] = [
        c"RESTFRQ", c"RESTFREQ", c"RESTWAV", c"OBSGEO-X", c"OBSGEO-Y", c"OBSGEO-Z", c"VELOSYS",
        c"ZSOURCE", c"VELANGL", c"LONPOLE", c"LATPOLE",
    ];
    let ncfltnkeys = 11;

    /* floating WCS keywords w/ underscore  */
    let cflt_keys: [&std::ffi::CStr; 2] = [c"PC", c"CD"];
    let ncflt_keys = 2;

    /* string WCS keywords  */
    let cstrkeys: [&std::ffi::CStr; 4] = [c"CTYPE", c"CUNIT", c"PS", c"CNAME"];
    let ncstrkeys = 4;

    /* string RADESYS keywords with list of allowed values  */
    let rastrkeys: [&std::ffi::CStr; 2] = [c"RADESYS", c"RADECSYS"];
    let nrastrkeys = 2;

    /* string spectral ref frame keywords with list of allowed values  */
    let specstrkeys: [&std::ffi::CStr; 3] = [c"SPECSYS", c"SSYSOBS", c"SSYSSRC"];
    let nspecstrkeys = 3;

    for i in 0..NWCSDESCR {
        wcsaxespos[i] = -1;
        wcskeypos[i] = 1000000000;
        pname[i] = None;
    }

    numusrkey = hduptr.tkeys;
    let tmpkwds = tmpkwds_get();

    /* find the extension  name and version */
    key_match(&tmpkwds, numusrkey, cs!(c"EXTNAME"), 1, &mut k, &mut n);
    if k > -1 && hduptr.kwds[k as usize].ktype == kwdtyp::STR_KEY {
        strcpy_c(&mut hduptr.extname, &hduptr.kwds[k as usize].kvalue.clone());
    }

    key_match(&tmpkwds, numusrkey, cs!(c"EXTVER"), 1, &mut k, &mut n);
    if k > -1 && hduptr.kwds[k as usize].ktype == kwdtyp::INT_KEY {
        hduptr.extver = strtol_c(&hduptr.kwds[k as usize].kvalue).0 as c_int;
    }

    /* set the HduName structure */
    hdunum = hduptr.hdunum;
    let extname = hduptr.extname;
    set_hduname(hdunum, hduptr.hdutype, Some(&extname), hduptr.extver);

    if hduptr.hdunum == 1 {
        test_prm(infits, out, hduptr);
        primary_naxis = hduptr.naxis;
    } else {
        /* test the keywords specific to the hdutype*/
        match hduptr.hdutype {
            IMAGE_HDU => test_img_ext(infits, out, hduptr),
            ASCII_TBL => test_asc_ext(infits, out, hduptr),
            BINARY_TBL => test_bin_ext(infits, out, hduptr),
            _ => {}
        }
    }
    /* test the general keywords */
    test_header(infits, out, hduptr);

    /* Check INHERIT keyword; must not be used if primary contains data */
    key_match(&tmpkwds, numusrkey, cs!(c"INHERIT"), 1, &mut k, &mut n);
    if k > -1 {
        if primary_naxis != 0 {
            spf!(errmes;
                "Keyword #", hduptr.kwds[k as usize].kindex, ", ",
                CS(&hduptr.kwds[k as usize].kname),
                " cannot be used if the primary array contains data (NAXIS != 0).");
            wrtwrn(out, &errmes, 0);
        }
        check_log(&hduptr.kwds[k as usize], out);
    }

    /* test if CROTA2 exists; if so, then PCi_j must not exist */
    key_match(&tmpkwds, numusrkey, cs!(c"CROTA2"), 1, &mut k, &mut n);
    if n == 1 {
        crota2_exists = hduptr.kwds[k as usize].kindex;
    }

    /* first find the primary WCSAXES value, if it exists */
    key_match(&tmpkwds, numusrkey, cs!(c"WCSAXES"), 1, &mut k, &mut n);
    if k >= 0 {
        j = k;
        if check_int(&hduptr.kwds[j as usize], out) != 0 {
            wcsaxesvalue = strtol_c(&hduptr.kwds[j as usize].kvalue).0 as c_int;
            nmax = wcsaxesvalue;
            if wcsaxesvalue > wcsaxesmax {
                wcsaxesmax = wcsaxesvalue;
            }
            wcsaxesExists = 1;

            /* store index of the wcsaxes keyword */
            /*  (it must appear before other WCS keywords) */
            if hduptr.kwds[j as usize].kindex > wcsaxespos[0] {
                wcsaxespos[0] = hduptr.kwds[j as usize].kindex;
            }
        }
    }

    /* Check and find max value of the WCSAXESa keywords */
    /* Use the max value when checking the range of the indexed WCS keywords. */
    /* This is a less rigorous test than if one were to test the range of the */
    /* keywords for each of the alternate WCS systems (A - Z) against the */
    /* corresponding WCSAXESa keyword.  */

    key_match(&tmpkwds, numusrkey, cs!(c"WCSAXES"), 0, &mut k, &mut n);

    j = k;
    while j < n + k {
        if check_int(&hduptr.kwds[j as usize], out) != 0 {
            taxes = strtol_c(&hduptr.kwds[j as usize].kvalue).0 as c_int;
            if taxes > wcsaxesmax {
                wcsaxesmax = taxes;
            }
            wcsaxesExists = 1;

            /*  Removed this check on 6/28/2012.  See discussion on FITSBITS related
                to this requirement.  */

            parse_wcskey_suffix(
                &hduptr.kwds[j as usize].kname,
                cs!(c"WCSAXES"),
                &mut iaxis,
                &mut ialt,
            );
            wcsaxespos[ialt as usize] = hduptr.kwds[j as usize].kindex;
        }
        j += 1;
    }

    /* test datatype of reserved indexed floating point WCS keywords */
    for i in 0..ncfltkeys {
        let temp: &[c_char] = cast_slice(cfltkeys[i].to_bytes_with_nul());
        key_match(&tmpkwds, numusrkey, temp, 0, &mut k, &mut n);
        if k < 0 {
            continue;
        }

        j = k;
        while j < k + n {
            let kindex = hduptr.kwds[j as usize].kindex;
            let kname = hduptr.kwds[j as usize].kname;
            let kvalue = hduptr.kwds[j as usize].kvalue;

            p = cstrlen(temp);
            if !isdigit_c(kname[p]) {
                j += 1;
                continue;
            }

            if check_flt(&hduptr.kwds[j as usize], out) == 0 {
                j += 1;
                continue;
            }

            if i == 2 {
                /* test that CDELTi != 0 */
                dvalue = strtod_c(&kvalue);
                if dvalue == 0. {
                    spf!(errmes;
                        "Keyword #", kindex, ", ", CS(&kname), ": must have non-zero value.");
                    wrterr(out, &errmes, 1);
                }
            }

            if i == 4 || i == 5 {
                /* test that CRDERi and CSYSERi are non-negative */
                dvalue = strtod_c(&kvalue);
                if dvalue < 0. {
                    spf!(errmes;
                        "Keyword #", kindex, ", ", CS(&kname),
                        ": must have non-negative value: ", CS(&kvalue));
                    wrterr(out, &errmes, 1);
                }
            }

            parse_wcskey_suffix(&kname, temp, &mut iaxis, &mut ialt);
            if wcsaxesExists != 0 {
                /* WCSAXES keyword exists */

                if iaxis < 1 || iaxis > wcsaxesmax {
                    spf!(errmes;
                        "Keyword #", kindex, ", ", CS(&kname), ": index ", iaxis,
                        " is not in range 1-", wcsaxesmax, " (WCSAXES).");
                    wrterr(out, &errmes, 1);
                }
            } else if iaxis < 1 || iaxis > hduptr.naxis {
                spf!(errmes;
                    "Keyword #", kindex, ", ", CS(&kname), ": index ", iaxis,
                    " is not in range 1-", hduptr.naxis, " (NAXIS).");
                wrtwrn(out, &errmes, 0);
            }

            /* count the number of each keyword */
            if ialt == 0 {
                /* only test the primary set of WCS keywords */
                keynum[i] += 1;
                if iaxis > nmax {
                    nmax = iaxis;
                }
            }

            /* store lowest index of any wcs keyword */
            if kindex < wcskeypos[ialt as usize] {
                wcskeypos[ialt as usize] = kindex;
                pname[ialt as usize] = Some(kname);
            }
            j += 1;
        }
    }

    if wcsaxesvalue == 0 {
        /* limit value of nmax to the legal maximum */
        if nmax > hduptr.naxis {
            nmax = hduptr.naxis;
        }
    } else if nmax > wcsaxesvalue {
        nmax = wcsaxesvalue;
    }

    if keynum[0] < nmax {
        /* test number of CRPIXi keywords */
        spf!(errmes; "Some CRPIXi keywords appear to be missing; expected ", nmax, ".");
        wrtwrn(out, &errmes, 0);
    }
    if keynum[1] < nmax {
        /* test number of CRVALi keywords */
        spf!(errmes; "Some CRVALi keywords appear to be missing; expected ", nmax, ".");
        wrtwrn(out, &errmes, 0);
    }

    /* test datatype of reserved non-indexed floating point WCS keywords */
    for i in 0..ncfltnkeys {
        let temp: &[c_char] = cast_slice(cfltnkeys[i].to_bytes_with_nul());
        key_match(&tmpkwds, numusrkey, temp, 0, &mut k, &mut n);

        if k < 0 {
            continue;
        }

        j = k;
        while j < k + n {
            if check_flt(&hduptr.kwds[j as usize], out) == 0 {
                j += 1;
                continue;
            }
            j += 1;
        }
    }

    /* test datatype of reserved indexed floating point WCS keywords with "_" */
    for i in 0..ncflt_keys {
        let temp: &[c_char] = cast_slice(cflt_keys[i].to_bytes_with_nul());
        key_match(&tmpkwds, numusrkey, temp, 0, &mut k, &mut n);
        if k < 0 {
            continue;
        }

        j = k;
        while j < k + n {
            let kindex = hduptr.kwds[j as usize].kindex;
            let kname = hduptr.kwds[j as usize].kname;

            p = cstrlen(temp);
            if !isdigit_c(kname[p]) {
                j += 1;
                continue;
            }

            /* 2 digits must be separated by a '_' */
            let p2 = match kname[p..].iter().position(|&c| c as u8 == b'_') {
                Some(o) => p + o,
                None => {
                    j += 1;
                    continue;
                }
            };

            if check_flt(&hduptr.kwds[j as usize], out) == 0 {
                j += 1;
                continue;
            }

            /* test the first digit (strtol stops at the '_' anyway, so the C's
               temporary NUL there is not needed) */
            m = strtol_c(&kname[p..]).0 as c_int;

            if wcsaxesExists != 0 {
                /* WCSAXES keyword exists */

                if m < 1 || m > wcsaxesmax {
                    spf!(errmes;
                        "Keyword #", kindex, ", ", CS(&kname), ": 1st index ", m,
                        " is not in range 1-", wcsaxesmax, " (WCSAXES).");
                    wrterr(out, &errmes, 1);
                }
            } else if m < 1 || m > hduptr.naxis {
                spf!(errmes;
                    "Keyword #", kindex, ", ", CS(&kname), ": 1st index ", m,
                    " is not in range 1-", hduptr.naxis, " (NAXIS).");
                wrtwrn(out, &errmes, 0);
            }

            /* test the second digit */
            let pp = p2 + 1;
            let (mv, endoff) = strtol_c(&kname[pp..]);
            m = mv as c_int;
            let p2 = pp + endoff;

            if wcsaxesExists != 0 {
                /* WCSAXES keyword exists */

                if m < 1 || m > wcsaxesmax {
                    spf!(errmes;
                        "Keyword #", kindex, ", ", CS(&kname), ": 2nd index ", m,
                        " is not in range 1-", wcsaxesmax, " (WCSAXES).");
                    wrterr(out, &errmes, 1);
                }
            } else if m < 1 || m > hduptr.naxis {
                spf!(errmes;
                    "Keyword #", kindex, ", ", CS(&kname), ": 2nd index ", m,
                    " is not in range 1-", hduptr.naxis, " (NAXIS).");
                wrtwrn(out, &errmes, 0);
            }

            ialt = 0;
            if kname[p2] == 0 {
                /* no alternate suffix on the PC or CD name */
                matrix_exists[i] = kindex;
            } else if cstrlen(&kname[p2..]) == 1 {
                let c = kname[p2] as u8;
                if c.is_ascii_uppercase() {
                    ialt = c as c_int - 64;
                }
            }

            /* store lowest index of any wcs keyword */
            if kindex < wcskeypos[ialt as usize] {
                wcskeypos[ialt as usize] = kindex;
                pname[ialt as usize] = Some(kname);
            }
            j += 1;
        }
    }

    if matrix_exists[0] > 0 && matrix_exists[1] > 0 {
        spf!(errmes;
            "Keywords PCi_j (#", matrix_exists[0], ") and CDi_j (#", matrix_exists[1],
            ") are mutually exclusive.");
        wrterr(out, &errmes, 1);
    }

    if matrix_exists[0] > 0 && crota2_exists > 0 {
        spf!(errmes;
            "Keywords PCi_j (#", matrix_exists[0], ") and CROTA2 (#", crota2_exists,
            ") are mutually exclusive.");
        wrterr(out, &errmes, 1);
    }

    /* test datatype of reserved indexed string WCS keywords */
    for i in 0..ncstrkeys {
        let temp: &[c_char] = cast_slice(cstrkeys[i].to_bytes_with_nul());
        keynum[i] = 0;
        key_match(&tmpkwds, numusrkey, temp, 0, &mut k, &mut n);

        if k < 0 {
            continue;
        }

        j = k;
        while j < k + n {
            let kindex = hduptr.kwds[j as usize].kindex;
            let kname = hduptr.kwds[j as usize].kname;

            p = cstrlen(temp);
            if !isdigit_c(kname[p]) {
                j += 1;
                continue;
            }

            if check_str(&hduptr.kwds[j as usize], out) == 0 {
                j += 1;
                continue;
            }

            parse_wcskey_suffix(&kname, temp, &mut iaxis, &mut ialt);

            if wcsaxesExists != 0 {
                /* WCSAXES keyword exists */

                if iaxis < 1 || iaxis > wcsaxesmax {
                    spf!(errmes;
                        "Keyword #", kindex, ", ", CS(&kname), ": index ", iaxis,
                        " is not in range 1-", wcsaxesmax, " (WCSAXES).");
                    wrterr(out, &errmes, 1);
                }
            } else if iaxis < 1 || iaxis > hduptr.naxis {
                spf!(errmes;
                    "Keyword #", kindex, ", ", CS(&kname), ": index ", iaxis,
                    " is not in range 1-", hduptr.naxis, " (NAXIS).");
                wrtwrn(out, &errmes, 0);
            }

            if ialt == 0 {
                /* only test the primary set of WCS keywords */
                keynum[i] += 1;
            }

            /* store lowest index of any wcs keyword */
            if kindex < wcskeypos[ialt as usize] {
                wcskeypos[ialt as usize] = kindex;
                pname[ialt as usize] = Some(kname);
            }
            j += 1;
        }
    }

    if keynum[0] < nmax {
        spf!(errmes; "Some CTYPEi keywords appear to be missing; expected ", nmax, ".");
        wrtwrn(out, &errmes, 0);
    }

    for i in 0..NWCSDESCR {
        if wcskeypos[i] < wcsaxespos[i] {
            let nm = pname[i].unwrap_or([0; FLEN_KEYWORD]);
            spf!(errmes;
                "WCSAXES keyword #", wcsaxespos[i], " appears after other WCS keyword ",
                CS(&nm), " #", wcskeypos[i]);
            wrterr(out, &errmes, 1);
        }
    }

    /* test datatype and value of reserved RADECSYS WCS keywords */
    for i in 0..nrastrkeys {
        let temp: &[c_char] = cast_slice(rastrkeys[i].to_bytes_with_nul());
        keynum[i] = 0;
        key_match(&tmpkwds, numusrkey, temp, 0, &mut k, &mut n);

        if k < 0 {
            continue;
        }

        j = k;
        while j < k + n {
            let kindex = hduptr.kwds[j as usize].kindex;
            let kname = hduptr.kwds[j as usize].kname;
            let kvalue = hduptr.kwds[j as usize].kvalue;

            if check_str(&hduptr.kwds[j as usize], out) == 0 {
                j += 1;
                continue;
            }

            if !ceq(&kvalue, b"ICRS")
                && !ceq(&kvalue, b"FK5")
                && !ceq(&kvalue, b"FK4")
                && !ceq(&kvalue, b"FK4-NO-E")
                && !ceq(&kvalue, b"GAPPT")
            {
                spf!(errmes;
                    "Keyword #", kindex, ", ", CS(&kname), " has non-allowed value: ",
                    CS(&kvalue));
                wrtwrn(out, &errmes, 0);
            }
            j += 1;
        }
    }

    /* test datatype and value of reserved spectral ref frame WCS keywords */
    for i in 0..nspecstrkeys {
        let temp: &[c_char] = cast_slice(specstrkeys[i].to_bytes_with_nul());
        keynum[i] = 0;
        key_match(&tmpkwds, numusrkey, temp, 0, &mut k, &mut n);

        if k < 0 {
            continue;
        }

        j = k;
        while j < k + n {
            let kindex = hduptr.kwds[j as usize].kindex;
            let kname = hduptr.kwds[j as usize].kname;
            let kvalue = hduptr.kwds[j as usize].kvalue;

            if check_str(&hduptr.kwds[j as usize], out) == 0 {
                j += 1;
                continue;
            }

            if !ceq(&kvalue, b"TOPOCENT")
                && !ceq(&kvalue, b"GEOCENTR")
                && !ceq(&kvalue, b"BARYCENT")
                && !ceq(&kvalue, b"HELIOCEN")
                && !ceq(&kvalue, b"LSRK")
                && !ceq(&kvalue, b"LSRD")
                && !ceq(&kvalue, b"GALACTOC")
                && !ceq(&kvalue, b"LOCALGRP")
                && !ceq(&kvalue, b"CMBDIPOL")
                && !ceq(&kvalue, b"SOURCE")
            {
                spf!(errmes;
                    "Keyword #", kindex, ", ", CS(&kname), " has non-allowed value: ",
                    CS(&kvalue));
                wrtwrn(out, &errmes, 0);
            }
            j += 1;
        }
    }

    /* test the fill area */
    if testfill() != 0 && ffchfl_safe(infits, &mut status) != 0 {
        wrterr_str(
            out,
            "The header fill area is not totally filled with blanks.",
            1,
        );
    }
}

/*************************************************************
*
*      test_prm
*
*   Test the primary array header
*
*************************************************************/
fn test_prm(
    infits: &mut fitsfile, /* input fits file   */
    out: Out,              /* output ascii file */
    hduptr: &mut FitsHdu,  /* hdu information structure */
) {
    let mut i: c_int;
    let mut j: c_int;
    let mut k: c_int = 0;
    let mut n: c_int = 0;
    let mut p: usize;
    let mut errmes: [c_char; ERRMES_LEN] = [0; ERRMES_LEN];

    let exlkey: [&std::ffi::CStr; 2] = [c"XTENSION", c"INHERIT"];
    let nexlkey = 1;

    let numusrkey = hduptr.tkeys;
    let tmpkwds = tmpkwds_get();

    /* The SIMPLE, BITPIX, NAXIS, and NAXISn keywords  have been
       checked in CFITSIO */

    /* excluded keywords cannot be used. */
    for i in 0..nexlkey {
        let temp: &[c_char] = cast_slice(exlkey[i].to_bytes_with_nul());
        key_match(&tmpkwds, numusrkey, temp, 1, &mut k, &mut n);
        if n > 0 {
            spf!(errmes;
                "Keyword #", hduptr.kwds[k as usize].kindex, ", ", CS(temp),
                " is not allowed in a primary array.");
            wrterr(out, &errmes, 1);
        }
    }

    /* Check if Random Groups file */
    key_match(&tmpkwds, numusrkey, cs!(c"GROUPS"), 1, &mut k, &mut n);
    if k > -1 {
        let kindex = hduptr.kwds[k as usize].kindex;
        if hduptr.kwds[k as usize].kvalue[0] as u8 == b'T'
            && hduptr.naxis > 0
            && hduptr.naxes[0] == 0
        {
            hduptr.isgroup = 1;

            check_fixed_log(&card_at(kindex as usize - 1), out);
        }
    }

    /* check the position of the EXTEND  */

    /*  the EXTEND keyword is no longer required if the file contains extensions */

    if hduptr.isgroup == 0 {
        key_match(&tmpkwds, numusrkey, cs!(c"EXTEND"), 1, &mut k, &mut n);
        if k > 0 {
            if check_log(&hduptr.kwds[k as usize], out) != 0
                && hduptr.kwds[k as usize].kvalue[0] as u8 != b'T'
                && totalhdu() > 1
            {
                spf!(errmes; "There are extensions but EXTEND = F.");
                wrterr(out, &errmes, 1);
            }
        }
    }

    /* Check PCOUNT and GCOUNT  keyword */
    key_match(&tmpkwds, numusrkey, cs!(c"PCOUNT"), 1, &mut k, &mut n);
    if k > -1 {
        let kindex = hduptr.kwds[k as usize].kindex;
        let kname = hduptr.kwds[k as usize].kname;
        /* Primary array cannot have PCOUNT */
        if hduptr.isgroup == 0 {
            spf!(errmes;
                " Keyword #", kindex, ", ", CS(&kname), " is not allowed in a primary array.");
            wrterr(out, &errmes, 1);
        } else {
            if check_int(&hduptr.kwds[k as usize], out) != 0 {
                hduptr.pcount = strtod_c(&hduptr.kwds[k as usize].kvalue) as LONGLONG;
            }

            check_fixed_int(&card_at(kindex as usize - 1), out);
        }
    }

    key_match(&tmpkwds, numusrkey, cs!(c"GCOUNT"), 1, &mut k, &mut n);
    if k > -1 {
        let kindex = hduptr.kwds[k as usize].kindex;
        let kname = hduptr.kwds[k as usize].kname;
        /* Primary array cannot have GCOUNT */
        if hduptr.isgroup == 0 {
            spf!(errmes;
                " Keyword #", kindex, ", ", CS(&kname), " is not allowed in a primary array.");
            wrterr(out, &errmes, 1);
        } else {
            if check_int(&hduptr.kwds[k as usize], out) != 0 {
                hduptr.gcount = strtol_c(&hduptr.kwds[k as usize].kvalue).0 as c_int;
            }

            check_fixed_int(&card_at(kindex as usize - 1), out);
        }
    }

    key_match(&tmpkwds, numusrkey, cs!(c"BLOCKED"), 1, &mut k, &mut n);
    if k > -1 {
        spf!(errmes;
            "Keyword #", hduptr.kwds[k as usize].kindex, ", ",
            CS(&hduptr.kwds[k as usize].kname), " is deprecated.");
        wrtwrn(out, &errmes, 0);
        check_log(&hduptr.kwds[k as usize], out);
    }

    /*  Check PSCALn keywords (only in Random Groups) */
    key_match(&tmpkwds, numusrkey, cs!(c"PSCAL"), 0, &mut k, &mut n);
    j = k;
    while j < k + n {
        let kindex = hduptr.kwds[j as usize].kindex;
        let kname = hduptr.kwds[j as usize].kname;
        p = 5;
        if !isdigit_c(kname[p]) {
            j += 1;
            continue;
        }

        if hduptr.isgroup == 0 {
            spf!(errmes; "Keyword #", kindex, ", ", CS(&kname), " ");
            scat!(errmes; "is only allowed in Random Groups structures.");
            wrterr(out, &errmes, 1);
            j += 1;
            continue;
        }

        if check_flt(&hduptr.kwds[j as usize], out) != 0
            && strtod_c(&hduptr.kwds[j as usize].kvalue) == 0.0
        {
            spf!(errmes; "Keyword #", kindex, ", ", CS(&kname), ": ");
            scat!(errmes; "The scaling factor is zero.");
            wrtwrn(out, &errmes, 0);
        }

        i = strtol_c(&kname[p..]).0 as c_int - 1;
        if i < 0 || i >= hduptr.gcount {
            spf!(errmes;
                "Keyword #", kindex, ", ", CS(&kname), ": invalid index ", i + 1,
                " (> GCOUNT = ", hduptr.gcount, ").");
            wrterr(out, &errmes, 1);
            j += 1;
            continue;
        }
        j += 1;
    }

    /*  Check PZEROn keywords (only in Random Groups) */
    key_match(&tmpkwds, numusrkey, cs!(c"PZERO"), 0, &mut k, &mut n);
    j = k;
    while j < k + n {
        let kindex = hduptr.kwds[j as usize].kindex;
        let kname = hduptr.kwds[j as usize].kname;
        p = 5;
        if !isdigit_c(kname[p]) {
            j += 1;
            continue;
        }

        if hduptr.isgroup == 0 {
            spf!(errmes; "Keyword #", kindex, ", ", CS(&kname), " ");
            scat!(errmes; "is only allowed in Random Groups structures.");
            wrterr(out, &errmes, 1);
            j += 1;
            continue;
        }

        check_flt(&hduptr.kwds[j as usize], out);
        i = strtol_c(&kname[p..]).0 as c_int - 1;
        if i < 0 || i >= hduptr.gcount {
            spf!(errmes;
                "Keyword #", kindex, ", ", CS(&kname), ": invalid index ", i + 1,
                " (> GCOUNT = ", hduptr.gcount, ").");
            wrterr(out, &errmes, 1);
            j += 1;
            continue;
        }
        j += 1;
    }

    /*  Check PTYPEn keywords (only in Random Groups) */
    key_match(&tmpkwds, numusrkey, cs!(c"PTYPE"), 0, &mut k, &mut n);
    j = k;
    while j < k + n {
        let kindex = hduptr.kwds[j as usize].kindex;
        let kname = hduptr.kwds[j as usize].kname;
        p = 5;
        if !isdigit_c(kname[p]) {
            j += 1;
            continue;
        }

        if hduptr.isgroup == 0 {
            spf!(errmes; "Keyword #", kindex, ", ", CS(&kname), " ");
            scat!(errmes; "is only allowed in Random Groups structures.");
            wrterr(out, &errmes, 1);
            j += 1;
            continue;
        }

        check_str(&hduptr.kwds[j as usize], out);
        i = strtol_c(&kname[p..]).0 as c_int - 1;
        if i < 0 || i >= hduptr.gcount {
            spf!(errmes;
                "Keyword #", kindex, ", ", CS(&kname), ": invalid index ", i + 1,
                " (> GCOUNT = ", hduptr.gcount, ").");
            wrterr(out, &errmes, 1);
            j += 1;
            continue;
        }
        j += 1;
    }
    test_array(infits, out, hduptr);
}

/*************************************************************
*
*      test_ext
*
*   Test the extension header
*
*************************************************************/
fn test_ext(
    _infits: &mut fitsfile, /* input fits file   */
    out: Out,               /* output ascii file */
    hduptr: &mut FitsHdu,   /* information about header */
) {
    let mut i: c_int;
    let mut j: c_int;
    let mut k: c_int = 0;
    let mut n: c_int = 0;
    let mut p: usize;
    let exlkey: [&std::ffi::CStr; 3] = [c"SIMPLE", c"EXTEND", c"BLOCKED"];
    let nexlkey = 3;
    let exlnkey: [&std::ffi::CStr; 4] = [c"PTYPE", c"PSCAL", c"PZERO", c"GROUPS"];
    let nexlnkey = 4;
    let hdunum: c_int;
    let mut errmes: [c_char; ERRMES_LEN] = [0; ERRMES_LEN];
    let mut comm: [c_char; COMM_LEN] = [0; COMM_LEN];

    let numusrkey = hduptr.tkeys;
    let tmpkwds = tmpkwds_get();
    hdunum = hduptr.hdunum;

    /* check the duplicate extensions */
    i = hdunum - 1;
    while i > 0 {
        if test_hduname(hdunum, i) != 0 {
            spf!(comm; "The HDU ", hdunum, " and ", i, " have identical type/name/version");
            wrtwrn(out, &comm, 0);
        }
        i -= 1;
    }

    /* check the position of the PCOUNT  */
    key_match(&tmpkwds, numusrkey, cs!(c"PCOUNT"), 1, &mut k, &mut n);
    if k < 0 {
        spf!(errmes; "cannot find the PCOUNT keyword.");
        wrterr(out, &errmes, 1);
    } else {
        let kindex = hduptr.kwds[k as usize].kindex;
        if check_int(&hduptr.kwds[k as usize], out) != 0 {
            hduptr.pcount = strtod_c(&hduptr.kwds[k as usize].kvalue) as LONGLONG;
        }
        if kindex != 4 + hduptr.naxis {
            spf!(errmes; "PCOUNT is not in record ", hduptr.naxis + 4, " of the header.");
            wrterr(out, &errmes, 1);
        }

        check_fixed_int(&card_at(kindex as usize - 1), out);
    }

    /* check the position of the GCOUNT */
    key_match(&tmpkwds, numusrkey, cs!(c"GCOUNT"), 1, &mut k, &mut n);
    if k < 0 {
        spf!(errmes; "cannot find the GCOUNT keyword.");
        wrterr(out, &errmes, 1);
    } else {
        let kindex = hduptr.kwds[k as usize].kindex;
        if check_int(&hduptr.kwds[k as usize], out) != 0 {
            hduptr.gcount = strtol_c(&hduptr.kwds[k as usize].kvalue).0 as c_int;
        }
        if kindex != 5 + hduptr.naxis {
            spf!(errmes; "GCOUNT is not in record ", hduptr.naxis + 5, " of the header.");
            wrterr(out, &errmes, 1);
        }

        check_fixed_int(&card_at(kindex as usize - 1), out);
    }

    for i in 0..nexlkey {
        let temp: &[c_char] = cast_slice(exlkey[i].to_bytes_with_nul());
        key_match(&tmpkwds, numusrkey, temp, 1, &mut k, &mut n);
        if k > -1 {
            spf!(errmes;
                "Keyword #", hduptr.kwds[k as usize].kindex, ", ",
                CS(&hduptr.kwds[k as usize].kname), " is not allowed in extensions.");
            wrterr(out, &errmes, 1);
        }
    }

    for i in 0..nexlnkey {
        let temp: &[c_char] = cast_slice(exlnkey[i].to_bytes_with_nul());
        key_match(&tmpkwds, numusrkey, temp, 0, &mut k, &mut n);
        if k > -1 {
            j = k;
            while j < k + n {
                let kname = hduptr.kwds[j as usize].kname;
                p = 5;
                if !isdigit_c(kname[p]) {
                    j += 1;
                    continue;
                }

                spf!(errmes;
                    "Keyword #", hduptr.kwds[j as usize].kindex, ", ", CS(&kname),
                    " is only allowed in Random Groups structures.");
                wrterr(out, &errmes, 1);
                j += 1;
            }
        }
    }
}

/*************************************************************
*
*      test_img_ext
*
*   Test the image extension header
*
*************************************************************/
fn test_img_ext(
    infits: &mut fitsfile, /* input fits file   */
    out: Out,              /* output ascii file */
    hduptr: &mut FitsHdu,  /* information about header */
) {
    let mut errmes: [c_char; ERRMES_LEN] = [0; ERRMES_LEN];

    test_ext(infits, out, hduptr);

    /* The XTENSION, BITPIX, NAXIS, and NAXISn keywords  have been
       checked in CFITSIO */

    if hduptr.pcount != 0 && hduptr.pcount != -99 {
        spf!(errmes; "Illegal pcount value ", hduptr.pcount as i64, " for image ext.");
        wrterr(out, &errmes, 1);
    }

    if hduptr.gcount != 1 && hduptr.gcount != -99 {
        spf!(errmes; "Illegal gcount value ", hduptr.gcount, " for image ext.");
        wrterr(out, &errmes, 1);
    }

    test_array(infits, out, hduptr);
}

/*************************************************************
*
*      test_array
*
*   Test the keywords which are used by both the primary array
* and image Extension.
*
*************************************************************/
fn test_array(
    _infits: &mut fitsfile, /* input fits file   */
    out: Out,               /* output ascii file */
    hduptr: &mut FitsHdu,   /* information about header */
) {
    let mut p: usize;
    let mut j: c_int;
    let mut k: c_int = 0;
    let mut n: c_int = 0;
    let mut errmes: [c_char; ERRMES_LEN] = [0; ERRMES_LEN];

    /* excluded non-indexed keywords  */
    let exlkeys: [&std::ffi::CStr; 2] = [c"TFIELDS", c"THEAP"];
    let nexlkeys = 2;

    /* excluded indexed keywords */
    let exlnkeys: [&std::ffi::CStr; 15] = [
        c"TBCOL", c"TFORM", c"TSCAL", c"TZERO", c"TNULL", c"TTYPE", c"TUNIT", c"TDISP", c"TDIM",
        c"TCTYP", c"TCUNI", c"TCRVL", c"TCDLT", c"TCRPX", c"TCROT",
    ];
    let nexlnkeys = 15;

    /* non-indexed floating keywords  (excluding BSCALE) */
    let fltkeys: [&std::ffi::CStr; 3] = [c"BZERO", c"DATAMAX", c"DATAMIN"];
    let nfltkeys = 3;

    /* non-indexed string keywords */
    let strkeys: [&std::ffi::CStr; 1] = [c"BUNIT"];
    let nstrkeys = 1;

    let numusrkey = hduptr.tkeys;
    let tmpkwds = tmpkwds_get();

    /*  Check BLANK, BSCALE keywords */
    key_match(&tmpkwds, numusrkey, cs!(c"BLANK"), 1, &mut k, &mut n);
    if k >= 0 {
        check_int(&hduptr.kwds[k as usize], out);
        if hduptr.bitpix < 0 {
            spf!(errmes;
                "Keyword #", hduptr.kwds[k as usize].kindex, ", ",
                CS(&hduptr.kwds[k as usize].kname),
                " must not be used with floating point data (BITPIX = ", hduptr.bitpix, ").");
            wrterr(out, &errmes, 2);
        }
    }

    key_match(&tmpkwds, numusrkey, cs!(c"BSCALE"), 1, &mut k, &mut n);
    if k >= 0
        && check_flt(&hduptr.kwds[k as usize], out) != 0
        && strtod_c(&hduptr.kwds[k as usize].kvalue) == 0.0
    {
        spf!(errmes;
            "Keyword #", hduptr.kwds[k as usize].kindex, ", ",
            CS(&hduptr.kwds[k as usize].kname), ": The scaling factor is 0.");
        wrtwrn(out, &errmes, 0);
    }

    /* search for excluded, non-indexed keywords */
    for i in 0..nexlkeys {
        let temp: &[c_char] = cast_slice(exlkeys[i].to_bytes_with_nul());
        key_match(&tmpkwds, numusrkey, temp, 1, &mut k, &mut n);
        if k < 0 {
            continue;
        }
        j = k;
        while j < k + n {
            spf!(errmes;
                "Keyword #", hduptr.kwds[j as usize].kindex, ", ",
                CS(&hduptr.kwds[j as usize].kname), " is not allowed in the array HDU.");
            wrterr(out, &errmes, 1);
            j += 1;
        }
    }

    /* search for excluded, indexed keywords */
    for i in 0..nexlnkeys {
        let temp: &[c_char] = cast_slice(exlnkeys[i].to_bytes_with_nul());
        key_match(&tmpkwds, numusrkey, temp, 0, &mut k, &mut n);
        if k < 0 {
            continue;
        }
        j = k;
        while j < k + n {
            let kname = hduptr.kwds[j as usize].kname;
            p = cstrlen(temp);
            if !isdigit_c(kname[p]) {
                j += 1;
                continue;
            }

            spf!(errmes;
                "Keyword #", hduptr.kwds[j as usize].kindex, ", ", CS(&kname),
                " is not allowed in the array HDU.");
            wrterr(out, &errmes, 1);
            j += 1;
        }
    }

    /* test datatype of reserved non-indexed floating point keywords */
    for i in 0..nfltkeys {
        let temp: &[c_char] = cast_slice(fltkeys[i].to_bytes_with_nul());
        key_match(&tmpkwds, numusrkey, temp, 1, &mut k, &mut n);
        if k < 0 {
            continue;
        }
        j = k;
        while j < k + n {
            if check_flt(&hduptr.kwds[j as usize], out) == 0 {
                j += 1;
                continue;
            }
            j += 1;
        }
    }

    /* test datatype of reserved non-indexed string keywords */
    for i in 0..nstrkeys {
        let temp: &[c_char] = cast_slice(strkeys[i].to_bytes_with_nul());
        key_match(&tmpkwds, numusrkey, temp, 1, &mut k, &mut n);
        if k < 0 {
            continue;
        }
        j = k;
        while j < k + n {
            check_str(&hduptr.kwds[j as usize], out);
            j += 1;
        }
    }
}

/*************************************************************
*
*      test_tbl
*
*   Test the table extension header and fill the tform, ttype,
*   tunit.
*
*************************************************************/
fn test_tbl(
    _infits: &mut fitsfile, /* input fits file   */
    out: Out,               /* output ascii file */
    hduptr: &mut FitsHdu,   /* information about header */
) {
    let mut p: usize;
    let mut m: c_int;
    let mut n: c_int = 0;
    let mut i: c_int;
    let mut j: c_int;
    let mut k: c_int = 0;
    let mut w: c_long;
    let mut d: c_long;
    let mut e: c_long;
    let mut lm: c_long;
    let mcol: c_int;
    let mut errmes: [c_char; ERRMES_LEN] = [0; ERRMES_LEN];

    /* excluded, non-index keywords (allowed in tile-compressed images) */
    let exlkey: [&std::ffi::CStr; 6] = [
        c"BSCALE", c"BZERO", c"BUNIT", c"BLANK", c"DATAMAX", c"DATAMIN",
    ];
    let nexlkey = 6;

    /* floating WCS keywords  */
    let cfltkeys: [&std::ffi::CStr; 4] = [c"TCRVL", c"TCDLT", c"TCRPX", c"TCROT"];
    let ncfltkeys = 4;

    /* string WCS keywords  */
    let cstrkeys: [&std::ffi::CStr; 2] = [c"TCTYP", c"TCUNI"];
    let ncstrkeys = 2;

    let numusrkey = hduptr.tkeys;
    mcol = hduptr.ncols;
    let tmpkwds = tmpkwds_get();

    'otherkey: {
        if mcol <= 0 {
            break 'otherkey;
        }
        /* set the ttype, ttform, tunit for tables */
        TTYPE.with(|t| {
            let mut t = t.borrow_mut();
            t.clear();
            t.resize(mcol as usize, vec![0]);
        });
        TFORM.with(|t| {
            let mut t = t.borrow_mut();
            t.clear();
            t.resize(mcol as usize, vec![0]);
        });
        TUNIT.with(|t| {
            let mut t = t.borrow_mut();
            t.clear();
            t.resize(mcol as usize, vec![0]);
        });

        key_match(&tmpkwds, numusrkey, cs!(c"TFIELDS"), 1, &mut k, &mut n);
        if k >= 0 {
            check_fixed_int(&card_at(hduptr.kwds[k as usize].kindex as usize - 1), out);
        }

        key_match(&tmpkwds, numusrkey, cs!(c"TTYPE"), 0, &mut k, &mut n);
        j = k;
        while j < k + n {
            let kindex = hduptr.kwds[j as usize].kindex;
            let kname = hduptr.kwds[j as usize].kname;
            let kvalue = hduptr.kwds[j as usize].kvalue;
            p = 5;
            if !isdigit_c(kname[p]) {
                j += 1;
                continue;
            }

            check_str(&hduptr.kwds[j as usize], out);
            i = strtol_c(&kname[p..]).0 as c_int - 1;
            if i >= 0 && i < mcol {
                TTYPE.with(|t| t.borrow_mut()[i as usize] = cvec(&kvalue));
            } else {
                spf!(errmes;
                    "Keyword #", kindex, ", ", CS(&kname), ": invalid index ", i + 1,
                    " (> TFIELD = ", mcol, ").");
                wrterr(out, &errmes, 2);
            }
            j += 1;
        }

        key_match(&tmpkwds, numusrkey, cs!(c"TFORM"), 0, &mut k, &mut n);
        j = k;
        while j < k + n {
            let kindex = hduptr.kwds[j as usize].kindex;
            let kname = hduptr.kwds[j as usize].kname;
            let kvalue = hduptr.kwds[j as usize].kvalue;
            p = 5;
            if !isdigit_c(kname[p]) {
                j += 1;
                continue;
            }

            check_str(&hduptr.kwds[j as usize], out);

            /*  TFORMn keyword no longer required to be padded to at least 8 characters
                    check_fixed_str(cards[pkey->kindex - 1], out);
            */

            if kvalue[0] as u8 == b' ' {
                spf!(errmes;
                    "Keyword #", kindex, ", ", CS(&kname), ": TFORM=\"", CS(&kvalue), "\" ");
                scat!(errmes; "should not have leading space.");
                wrterr(out, &errmes, 1);
            }

            i = strtol_c(&kname[p..]).0 as c_int - 1;
            if i >= 0 && i < mcol {
                TFORM.with(|t| t.borrow_mut()[i as usize] = cvec(&kvalue));
            } else {
                spf!(errmes;
                    "Keyword #", kindex, ", ", CS(&kname), ": invalid index ", i + 1,
                    " (> TFIELD = ", mcol, ").");
                wrterr(out, &errmes, 2);
            }

            p = 0;
            while kvalue[p] as u8 != b' ' && kvalue[p] != 0 {
                let c = kvalue[p] as u8;
                if !c.is_ascii_digit()
                    && !isupper_c(kvalue[p])
                    && c != b'.'
                    && c != b')'
                    && c != b'('
                {
                    spf!(errmes;
                        "Keyword #", kindex, ", ", CS(&kname), ": The value ", CS(&kvalue),
                        " has character ", CHR(kvalue[p]),
                        " which is not uppercase letter.");
                    wrterr(out, &errmes, 1);
                }

                p += 1;
            }
            j += 1;
        }

        key_match(&tmpkwds, numusrkey, cs!(c"TUNIT"), 0, &mut k, &mut n);
        j = k;
        while j < k + n {
            let kindex = hduptr.kwds[j as usize].kindex;
            let kname = hduptr.kwds[j as usize].kname;
            let kvalue = hduptr.kwds[j as usize].kvalue;
            p = 5;
            if !isdigit_c(kname[p]) {
                j += 1;
                continue;
            }

            check_str(&hduptr.kwds[j as usize], out);
            i = strtol_c(&kname[p..]).0 as c_int - 1;
            if i >= 0 && i < mcol {
                TUNIT.with(|t| t.borrow_mut()[i as usize] = cvec(&kvalue));
            } else {
                spf!(errmes;
                    "Keyword #", kindex, ", ", CS(&kname), ": invalid index ", i + 1,
                    " (> TFIELD = ", mcol, ").");
                wrterr(out, &errmes, 1);
            }
            j += 1;
        }

        /*  Check TDISPn keywords */
        key_match(&tmpkwds, numusrkey, cs!(c"TDISP"), 0, &mut k, &mut n);
        j = k;
        while j < k + n {
            let kindex = hduptr.kwds[j as usize].kindex;
            let kname = hduptr.kwds[j as usize].kname;
            let kvalue = hduptr.kwds[j as usize].kvalue;
            p = 5;
            if !isdigit_c(kname[p]) {
                j += 1;
                continue;
            }

            if kvalue[0] == 0 {
                j += 1;
                continue;
            } /* ignore blank string */
            check_str(&hduptr.kwds[j as usize], out);
            if kvalue[0] as u8 == b' ' {
                spf!(errmes;
                    "Keyword #", kindex, ", ", CS(&kname), ": TDISP=\"", CS(&kvalue), "\" ");
                scat!(errmes; "should not have leading space.");
                wrterr(out, &errmes, 1);
            }

            i = strtol_c(&kname[p..]).0 as c_int - 1;
            if i < 0 || i >= mcol {
                spf!(errmes;
                    "Keyword #", kindex, ", ", CS(&kname), ": invalid index ", i + 1,
                    " (> TFIELD = ", mcol, ").");
                wrterr(out, &errmes, 1);
                j += 1;
                continue;
            }
            let tf = tform_at(i as usize);
            p = 0;
            match kvalue[p] as u8 {
                b'A' => {
                    p += 1;
                    w = strtol_c(&kvalue[p..]).0;
                    if w == 0 || w == c_long::MAX || w == c_long::MIN {
                        spf!(errmes;
                            "Keyword #", kindex, ", ", CS(&kname), ": invalid format \"",
                            CS(&kvalue), "\".");
                        wrterr(out, &errmes, 1);
                    }
                    if !chr_in(&tf, b'A') {
                        spf!(errmes;
                            "Keyword #", kindex, ", ", CS(&kname), ":  Format \"", CS(&kvalue),
                            "\" cannot be used for TFORM \"", CS(&tf), "\".");
                        wrterr(out, &errmes, 1);
                    }
                }
                b'L' => {
                    p += 1;
                    w = strtol_c(&kvalue[p..]).0;
                    if w == 0 || w == c_long::MAX || w == c_long::MIN {
                        spf!(errmes;
                            "Keyword #", kindex, ", ", CS(&kname), ": invalid format \"",
                            CS(&kvalue), "\".");
                        wrterr(out, &errmes, 1);
                    }
                    if !chr_in(&tf, b'L') {
                        spf!(errmes;
                            "Keyword #", kindex, ", ", CS(&kname), ":  Format ", CS(&kvalue),
                            " cannot be used for TFORM \"", CS(&tf), "\".");
                        wrterr(out, &errmes, 1);
                    }
                }
                b'I' | b'B' | b'O' | b'Z' => {
                    p += 1;
                    w = strtol_c(&kvalue[p..]).0;
                    if let Some(q) = chr_pos(&kvalue[p..], b'.') {
                        p += q;
                        p += 1;
                        lm = strtol_c(&kvalue[p..]).0;
                    } else {
                        lm = -1; /* no minimum digit field */
                    }
                    if w == 0
                        || w == c_long::MAX
                        || w == c_long::MIN
                        || lm == c_long::MAX
                        || lm == c_long::MIN
                        || w < lm
                    {
                        spf!(errmes;
                            "Keyword #", kindex, ", ", CS(&kname), ": invalid format \"",
                            CS(&kvalue), "\".");
                        wrterr(out, &errmes, 1);
                    }
                    if !chr_in(&tf, b'I')
                        && !chr_in(&tf, b'J')
                        && !chr_in(&tf, b'K')
                        && !chr_in(&tf, b'B')
                        && !chr_in(&tf, b'X')
                    {
                        spf!(errmes;
                            "Keyword #", kindex, ", ", CS(&kname), ":  Format \"", CS(&kvalue),
                            "\" cannot be used for TFORM \"", CS(&tf), "\".");
                        wrterr(out, &errmes, 1);
                    }
                }
                b'F' => {
                    p += 1;
                    d = -1;
                    w = strtol_c(&kvalue[p..]).0;
                    if let Some(q) = chr_pos(&kvalue[p..], b'.') {
                        p += q;
                        p += 1;
                        d = strtol_c(&kvalue[p..]).0;
                    }
                    if w == 0
                        || w == c_long::MAX
                        || w == c_long::MIN
                        || d == -1
                        || d == c_long::MAX
                        || d == c_long::MIN
                        || w < d + 1
                    {
                        spf!(errmes;
                            "Keyword #", kindex, ", ", CS(&kname), ": invalid format \"",
                            CS(&kvalue), "\".");
                        wrterr(out, &errmes, 1);
                    }
                    if !chr_in(&tf, b'E')
                        && !chr_in(&tf, b'F')
                        && !chr_in(&tf, b'C')
                        && !chr_in(&tf, b'D')
                        && !chr_in(&tf, b'M')
                        && !chr_in(&tf, b'I')
                        && !chr_in(&tf, b'J')
                        && !chr_in(&tf, b'K')
                        && !chr_in(&tf, b'B')
                        && !chr_in(&tf, b'X')
                    {
                        spf!(errmes;
                            "Keyword #", kindex, ", ", CS(&kname), ":  Format \"", CS(&kvalue),
                            "\" cannot be used for TFORM \"", CS(&tf), "\".");
                        wrterr(out, &errmes, 1);
                    }
                }
                b'E' | b'D' => {
                    p += 1;
                    d = 0;
                    if kvalue[p] as u8 == b'N' || kvalue[p] as u8 == b'S' {
                        p += 1;
                    }
                    w = strtol_c(&kvalue[p..]).0;
                    if let Some(q) = chr_pos(&kvalue[p..], b'.') {
                        p += q;
                        p += 1;
                        d = strtol_c(&kvalue[p..]).0;
                    }
                    if let Some(q) = chr_pos(&kvalue[p..], b'E') {
                        p += q;
                        p += 1;
                        e = strtol_c(&kvalue[p..]).0;
                    } else {
                        e = 2;
                    }
                    if w == 0
                        || w == c_long::MAX
                        || w == c_long::MIN
                        || d == 0
                        || d == c_long::MAX
                        || d == c_long::MIN
                        || e == 0
                        || e == c_long::MAX
                        || e == c_long::MIN
                        || w < d + e + 3
                    {
                        spf!(errmes;
                            "Keyword #", kindex, ", ", CS(&kname), ": invalid format \"",
                            CS(&kvalue), "\".");
                        wrterr(out, &errmes, 1);
                    }
                    if !chr_in(&tf, b'E')
                        && !chr_in(&tf, b'F')
                        && !chr_in(&tf, b'C')
                        && !chr_in(&tf, b'D')
                        && !chr_in(&tf, b'M')
                        && !chr_in(&tf, b'I')
                        && !chr_in(&tf, b'J')
                        && !chr_in(&tf, b'K')
                        && !chr_in(&tf, b'B')
                        && !chr_in(&tf, b'X')
                    {
                        spf!(errmes;
                            "Keyword #", kindex, ", ", CS(&kname), ":  Format \"", CS(&kvalue),
                            "\" cannot be used for TFORM \"", CS(&tf), "\".");
                        wrterr(out, &errmes, 1);
                    }
                }
                b'G' => {
                    p += 1;
                    d = 0;
                    w = strtol_c(&kvalue[p..]).0;
                    if let Some(q) = chr_pos(&kvalue[p..], b'.') {
                        p += q;
                        p += 1;
                        d = strtol_c(&kvalue[p..]).0;
                    }
                    if let Some(q) = chr_pos(&kvalue[p..], b'E') {
                        p += q;
                        p += 1;
                        e = strtol_c(&kvalue[p..]).0;
                    } else {
                        e = 2;
                    }
                    if w == 0
                        || w == c_long::MAX
                        || w == c_long::MIN
                        || d == 0
                        || d == c_long::MAX
                        || d == c_long::MIN
                        || e == 0
                        || e == c_long::MAX
                        || e == c_long::MIN
                    {
                        spf!(errmes;
                            "Keyword #", kindex, ", ", CS(&kname), ": invalid format \"",
                            CS(&kvalue), "\".");
                        wrterr(out, &errmes, 1);
                    }
                }
                _ => {
                    spf!(errmes;
                        "Keyword #", kindex, ", ", CS(&kname), ": invalid format \"",
                        CS(&kvalue), "\".");
                    wrterr(out, &errmes, 1);
                }
            }
            j += 1;
        }
    }

    /* OTHERKEY: */
    if hduptr.istilecompressed == 0 {
        /* tile compressed images can have these keywords */
        for i in 0..nexlkey {
            let temp: &[c_char] = cast_slice(exlkey[i].to_bytes_with_nul());
            key_match(&tmpkwds, numusrkey, temp, 1, &mut k, &mut n);
            if k > -1 {
                spf!(errmes;
                    "Keyword #", hduptr.kwds[k as usize].kindex, ", ",
                    CS(&hduptr.kwds[k as usize].kname),
                    " is not allowed in the Bin/ASCII table.");
                wrterr(out, &errmes, 1);
            }
        }

        /* these WCS keywords are all allowed (changed July 2010) */
    }

    /* test datatype of reserved indexed floating point WCS keywords */
    for i in 0..ncfltkeys {
        let temp: &[c_char] = cast_slice(cfltkeys[i].to_bytes_with_nul());
        key_match(&tmpkwds, numusrkey, temp, 0, &mut k, &mut n);
        if k < 0 {
            continue;
        }

        j = k;
        while j < k + n {
            let kindex = hduptr.kwds[j as usize].kindex;
            let kname = hduptr.kwds[j as usize].kname;

            p = cstrlen(temp);
            if !isdigit_c(kname[p]) {
                j += 1;
                continue;
            }

            if check_flt(&hduptr.kwds[j as usize], out) == 0 {
                j += 1;
                continue;
            }

            m = strtol_c(&kname[p..]).0 as c_int;
            if m < 1 || m > mcol {
                spf!(errmes;
                    "Keyword #", kindex, ", ", CS(&kname), ": index ", m,
                    " is not in range 1-", mcol, " (TFIELD).");
                wrterr(out, &errmes, 1);
            }
            j += 1;
        }
    }

    /* test datatype of reserved indexed string WCS keywords */
    for i in 0..ncstrkeys {
        let temp: &[c_char] = cast_slice(cstrkeys[i].to_bytes_with_nul());
        key_match(&tmpkwds, numusrkey, temp, 0, &mut k, &mut n);
        if k < 0 {
            continue;
        }

        j = k;
        while j < k + n {
            let kindex = hduptr.kwds[j as usize].kindex;
            let kname = hduptr.kwds[j as usize].kname;

            p = cstrlen(temp);
            if !isdigit_c(kname[p]) {
                j += 1;
                continue;
            }

            if check_str(&hduptr.kwds[j as usize], out) == 0 {
                j += 1;
                continue;
            }

            m = strtol_c(&kname[p..]).0 as c_int;
            if m < 1 || m > mcol {
                spf!(errmes;
                    "Keyword #", kindex, ", ", CS(&kname), ": index ", m,
                    " is not in range 1-", mcol, " (TFIELD).");
                wrterr(out, &errmes, 1);
            }
            j += 1;
        }
    }
}

/*************************************************************
*
*      test_asc_ext
*
*   Test the ascii table extension header
*
*************************************************************/
fn test_asc_ext(
    infits: &mut fitsfile, /* input fits file   */
    out: Out,              /* output ascii file */
    hduptr: &mut FitsHdu,  /* information about header */
) {
    let mut p: usize;
    let mut i: c_int;
    let mut j: c_int;
    let mut k: c_int = 0;
    let mut n: c_int = 0;
    let mcol: c_int;
    let mut errmes: [c_char; ERRMES_LEN] = [0; ERRMES_LEN];

    let numusrkey = hduptr.tkeys;
    mcol = hduptr.ncols;

    /* The XTENSION, BITPIX, NAXIS, NAXISn, TFIELDS, PCOUNT, GCOUNT, TFORMn,
       TBCOLn, TTYPEn keywords  have been checked in CFITSIO */

    /* General extension */
    test_ext(infits, out, hduptr);

    /* general table */
    test_tbl(infits, out, hduptr);

    let tmpkwds = tmpkwds_get();

    /* Check TBCOLn */
    key_match(&tmpkwds, numusrkey, cs!(c"TBCOL"), 0, &mut k, &mut n);
    j = k;
    while j < k + n {
        let kindex = hduptr.kwds[j as usize].kindex;
        let kname = hduptr.kwds[j as usize].kname;
        p = 5;
        if !isdigit_c(kname[p]) {
            j += 1;
            continue;
        }

        check_int(&hduptr.kwds[j as usize], out);

        i = strtol_c(&kname[p..]).0 as c_int;
        if i < 0 || i > mcol {
            spf!(errmes;
                "Keyword #", kindex, ", ", CS(&kname), ": invalid index ", i,
                " (> TFIELD = ", mcol, ").");
            wrterr(out, &errmes, 1);
        } else {
            check_fixed_int(&card_at(kindex as usize - 1), out);
        }
        j += 1;
    }

    /*  Check TNULLn, TSCALn, and TZEORn keywords */
    key_match(&tmpkwds, numusrkey, cs!(c"TNULL"), 0, &mut k, &mut n);
    j = k;
    while j < k + n {
        let kindex = hduptr.kwds[j as usize].kindex;
        let kname = hduptr.kwds[j as usize].kname;
        p = 5;
        if !isdigit_c(kname[p]) {
            j += 1;
            continue;
        }
        i = strtol_c(&kname[p..]).0 as c_int - 1;
        if i < 0 || i >= mcol {
            spf!(errmes;
                "Keyword #", kindex, ", ", CS(&kname), ": invalid index ", i + 1,
                " (> TFIELD = ", mcol, ").");
            wrterr(out, &errmes, 1);
        }
        check_str(&hduptr.kwds[j as usize], out);
        j += 1;
    }

    key_match(&tmpkwds, numusrkey, cs!(c"TSCAL"), 0, &mut k, &mut n);
    j = k;
    while j < k + n {
        let kindex = hduptr.kwds[j as usize].kindex;
        let kname = hduptr.kwds[j as usize].kname;
        p = 5;
        if !isdigit_c(kname[p]) {
            j += 1;
            continue;
        }
        i = strtol_c(&kname[p..]).0 as c_int - 1;
        if check_flt(&hduptr.kwds[j as usize], out) != 0
            && strtod_c(&hduptr.kwds[j as usize].kvalue) == 0.0
        {
            spf!(errmes; "Keyword #", kindex, ", ", CS(&kname), ": Scaling factor is zero.");
            wrtwrn(out, &errmes, 0);
        }
        if i < 0 || i >= mcol {
            spf!(errmes;
                "Keyword #", kindex, ", ", CS(&kname), ": invalid index ", i + 1,
                " (> TFIELD = ", mcol, ").");
            wrterr(out, &errmes, 1);
            j += 1;
            continue;
        }
        if chr_in(&tform_at(i as usize), b'A') {
            spf!(errmes;
                "Keyword #", kindex, ", ", CS(&kname),
                " may not be used for the A-format fields.");
            wrterr(out, &errmes, 1);
        }
        j += 1;
    }

    key_match(&tmpkwds, numusrkey, cs!(c"TZERO"), 0, &mut k, &mut n);
    j = k;
    while j < k + n {
        let kindex = hduptr.kwds[j as usize].kindex;
        let kname = hduptr.kwds[j as usize].kname;
        p = 5;
        if !isdigit_c(kname[p]) {
            j += 1;
            continue;
        }
        check_flt(&hduptr.kwds[j as usize], out);
        i = strtol_c(&kname[p..]).0 as c_int - 1;
        if i < 0 || i >= mcol {
            spf!(errmes;
                "Keyword #", kindex, ", ", CS(&kname), ": invalid index ", i + 1,
                " (> TFIELD = ", mcol, ").");
            wrterr(out, &errmes, 1);
            j += 1;
            continue;
        }
        if chr_in(&tform_at(i as usize), b'A') {
            spf!(errmes;
                "Keyword #", kindex, ", ", CS(&kname),
                " may not be used for the A-format fields.");
            wrterr(out, &errmes, 1);
        }
        j += 1;
    }

    key_match(&tmpkwds, numusrkey, cs!(c"TDIM"), 0, &mut k, &mut n);
    j = k;
    while j < k + n {
        let kname = hduptr.kwds[j as usize].kname;
        p = 4;
        if !isdigit_c(kname[p]) {
            j += 1;
            continue;
        }

        spf!(errmes;
            "Keyword #", hduptr.kwds[j as usize].kindex, ", ", CS(&kname),
            " is not allowed in the ASCII table.");
        wrterr(out, &errmes, 1);
        j += 1;
    }

    key_match(&tmpkwds, numusrkey, cs!(c"THEAP"), 1, &mut k, &mut n);
    if k > -1 {
        spf!(errmes;
            "Keyword #", hduptr.kwds[k as usize].kindex, ", ",
            CS(&hduptr.kwds[k as usize].kname), " is not allowed in the ASCII table.");
        wrterr(out, &errmes, 1);
    }

    /* check whether the column name is unique  */
    test_colnam(out, hduptr);
}

/*************************************************************
*
*      test_bin_ext
*
*   Test the binary table extension header
*
*************************************************************/
fn test_bin_ext(
    infits: &mut fitsfile, /* input fits file   */
    out: Out,              /* output ascii file */
    hduptr: &mut FitsHdu,  /* information about header */
) {
    let mut i: c_int;
    let mut j: c_int;
    let mut k: c_int = 0;
    let mut n: c_int = 0;
    let mut l: c_long;
    let mut status: c_int = 0;
    let mut p: usize;

    let mut ntdim: c_int = 0;
    let mut tdim: [c_long; 10] = [0; 10];
    let mut repeat: c_int;
    let mut width: c_int;
    let mcol: c_int;
    let mut vla: c_int;
    let mut datatype: c_int = 0;
    let mut errmes: [c_char; ERRMES_LEN] = [0; ERRMES_LEN];

    /* The indexed keywords excluded from ascii table */
    let exlkeys: [&std::ffi::CStr; 1] = [c"TBCOL"];
    let nexlkeys = 1;

    let numusrkey = hduptr.tkeys;
    mcol = hduptr.ncols;

    /* General extension */
    test_ext(infits, out, hduptr);

    /* General table */
    test_tbl(infits, out, hduptr);

    let tmpkwds = tmpkwds_get();

    /* The XTENSION, BITPIX, NAXIS, NAXISn, TFIELDS, PCOUNT, GCOUNT, TFORMn,
       TTYPEn keywords  have been checked in CFITSIO */

    /*  Check TNULLn keywords */
    key_match(&tmpkwds, numusrkey, cs!(c"TNULL"), 0, &mut k, &mut n);
    j = k;
    while j < k + n {
        let kindex = hduptr.kwds[j as usize].kindex;
        let kname = hduptr.kwds[j as usize].kname;
        let kvalue = hduptr.kwds[j as usize].kvalue;
        p = 5;
        if !isdigit_c(kname[p]) {
            j += 1;
            continue;
        }
        check_int(&hduptr.kwds[j as usize], out);
        i = strtol_c(&kname[p..]).0 as c_int - 1;
        if i < 0 || i >= mcol {
            spf!(errmes;
                "Keyword #", kindex, ", ", CS(&kname), ": invalid index ", i + 1,
                " (> TFIELD = ", mcol, ").");
            wrterr(out, &errmes, 1);
            j += 1;
            continue;
        }
        let tf = tform_at(i as usize);
        if !chr_in(&tf, b'B') && !chr_in(&tf, b'I') && !chr_in(&tf, b'J') && !chr_in(&tf, b'K') {
            spf!(errmes;
                "Keyword #", kindex, ", ", CS(&kname),
                " is used for the column with format \"", CS(&tf), " \".");
            wrterr(out, &errmes, 2);
        }
        l = strtol_c(&kvalue).0;
        if chr_in(&tf, b'B') && (l < 0 || l > 255) {
            spf!(errmes; "Keyword #", kindex, ", ", CS(&kname), ": The value ", l as i64);
            scat!(errmes; " is not in the range of datatype B.");
            wrtwrn(out, &errmes, 0);
        }
        l = strtol_c(&kvalue).0;
        if chr_in(&tf, b'I') && (l < -32768 || l > 32767) {
            spf!(errmes; "Keyword #", kindex, ", ", CS(&kname), ": The value ", l as i64);
            scat!(errmes; " is not in the range of datatype I ");
            wrtwrn(out, &errmes, 0);
        }
        j += 1;
    }

    /*  Check TSCALn keywords */
    key_match(&tmpkwds, numusrkey, cs!(c"TSCAL"), 0, &mut k, &mut n);
    j = k;
    while j < k + n {
        let kindex = hduptr.kwds[j as usize].kindex;
        let kname = hduptr.kwds[j as usize].kname;
        p = 5;
        if !isdigit_c(kname[p]) {
            j += 1;
            continue;
        }
        if check_flt(&hduptr.kwds[j as usize], out) != 0
            && strtod_c(&hduptr.kwds[j as usize].kvalue) == 0.0
        {
            spf!(errmes; "Keyword #", kindex, ", ", CS(&kname), ":");
            scat!(errmes; "The scaling factor is zero.");
            wrtwrn(out, &errmes, 0);
        }
        i = strtol_c(&kname[p..]).0 as c_int - 1;
        if i < 0 || i >= mcol {
            spf!(errmes;
                "Keyword #", kindex, ", ", CS(&kname), ": invalid index ", i + 1,
                " (> TFIELD = ", mcol, ").");
            wrterr(out, &errmes, 1);
            j += 1;
            continue;
        }
        let tf = tform_at(i as usize);
        if chr_in(&tf, b'A') || chr_in(&tf, b'L') || chr_in(&tf, b'X') {
            spf!(errmes;
                "Keyword #", kindex, ", ", CS(&kname), " is used in A, L, or X column. ");
            wrterr(out, &errmes, 1);
        }
        j += 1;
    }

    /*  Check TZEROn keywords */
    key_match(&tmpkwds, numusrkey, cs!(c"TZERO"), 0, &mut k, &mut n);
    j = k;
    while j < k + n {
        let kindex = hduptr.kwds[j as usize].kindex;
        let kname = hduptr.kwds[j as usize].kname;
        p = 5;
        if !isdigit_c(kname[p]) {
            j += 1;
            continue;
        }
        check_flt(&hduptr.kwds[j as usize], out);
        i = strtol_c(&kname[p..]).0 as c_int - 1;
        if i < 0 || i >= mcol {
            spf!(errmes;
                "Keyword #", kindex, ", ", CS(&kname), ": invalid index ", i + 1,
                " (> TFIELD = ", mcol, ").");
            wrterr(out, &errmes, 1);
            j += 1;
            continue;
        }
        let tf = tform_at(i as usize);
        if chr_in(&tf, b'A') && chr_in(&tf, b'L') && chr_in(&tf, b'X') {
            spf!(errmes;
                "Keyword #", kindex, ", ", CS(&kname), " is used in A, L, or X column. ");
            wrterr(out, &errmes, 1);
        }
        j += 1;
    }

    /* Check THEAP keyword */
    /* The C indexes naxes[0] and naxes[1] unconditionally; a corrupt header
       with NAXIS < 2 makes that an out-of-bounds read, so missing axes read
       as 0 here rather than as whatever happened to be on the stack. */
    let naxes0 = hduptr.naxes.first().copied().unwrap_or(0);
    let naxes1 = hduptr.naxes.get(1).copied().unwrap_or(0);
    hduptr.heap = (naxes0 * naxes1) as c_int;
    key_match(&tmpkwds, numusrkey, cs!(c"THEAP"), 1, &mut k, &mut n);
    if k > -1 {
        if check_int(&hduptr.kwds[k as usize], out) != 0 {
            hduptr.heap = strtol_c(&hduptr.kwds[k as usize].kvalue).0 as c_int;
        }
        if hduptr.pcount == 0 {
            spf!(errmes;
                "Pcount is zero, but keyword THEAP is present at record #",
                hduptr.kwds[k as usize].kindex, "). ");
            wrterr(out, &errmes, 1);
        }
    }

    /* if PCOUNT != 0, test that there is at least 1 variable length array column */
    vla = 0;
    if hduptr.pcount != 0 {
        for i in 0..mcol {
            if fits_get_coltype(infits, i + 1, Some(&mut datatype), None, None, &mut status) != 0 {
                spf!(errmes; "Column #", i, ": ");
                wrtferr(out, &errmes, &mut status, 2);
            }
            if datatype < 0 {
                vla = 1;
                break;
            }
        }

        if vla == 0 {
            spf!(errmes;
                "PCOUNT = ", hduptr.pcount as i64,
                ", but there are no variable-length array columns.");
            wrtwrn(out, &errmes, 0);
        }
    }

    /* Check TDIMn  keywords */
    key_match(&tmpkwds, numusrkey, cs!(c"TDIM"), 0, &mut k, &mut n);
    j = k;
    while j < k + n {
        let kindex = hduptr.kwds[j as usize].kindex;
        let kname = hduptr.kwds[j as usize].kname;
        let kvalue = hduptr.kwds[j as usize].kvalue;
        p = 4;
        if !isdigit_c(kname[p]) {
            j += 1;
            continue;
        }
        check_str(&hduptr.kwds[j as usize], out);
        if kvalue[0] as u8 == b' ' {
            spf!(errmes;
                "Keyword #", kindex, ", ", CS(&kname), ": TDIM=\"", CS(&kvalue), "\" ");
            scat!(errmes; "should not have leading space.");
            wrterr(out, &errmes, 1);
            j += 1;
            continue;
        }
        i = strtol_c(&kname[p..]).0 as c_int - 1;
        if i < 0 || i >= mcol {
            spf!(errmes;
                "Keyword #", kindex, ", ", CS(&kname), ": invalid index ", i + 1,
                " (> TFIELD = ", mcol, ").");
            wrterr(out, &errmes, 1);
            j += 1;
            continue;
        }
        if fits_decode_tdim(
            infits,
            &kvalue,
            i + 1,
            10,
            &mut ntdim,
            &mut tdim,
            &mut status,
        ) != 0
        {
            spf!(errmes; "Keyword #", kindex, ", ", CS(&kname), ": ");
            wrtferr(out, &errmes, &mut status, 1);
        }
        j += 1;
    }

    /* check the local convension "rAw"*/
    for i in 0..hduptr.ncols {
        let tf = tform_at(i as usize);
        let ap = match chr_pos(&tf, b'A') {
            Some(o) => o,
            None => continue,
        };
        repeat = strtol_c(&tf).0 as c_int;
        let p = ap + 1;
        if !isdigit_c(tf[p]) {
            continue;
        }
        width = strtol_c(&tf[p..]).0 as c_int;
        if width != 0 && repeat % width != 0 {
            spf!(errmes;
                "TFORM ", CS(&tf), " of column ", i + 1, ": repeat ", repeat,
                " is not the multiple of the width ", width);
            wrtwrn(out, &errmes, 0);
        }
    }

    for i in 0..nexlkeys {
        let temp: &[c_char] = cast_slice(exlkeys[i].to_bytes_with_nul());
        key_match(&tmpkwds, numusrkey, temp, 0, &mut k, &mut n);
        if k < 0 {
            continue;
        }
        j = k;
        while j < k + n {
            let kname = hduptr.kwds[j as usize].kname;
            p = cstrlen(temp);
            if !isdigit_c(kname[p]) {
                j += 1;
                continue;
            }

            spf!(errmes;
                "Keyword #", hduptr.kwds[j as usize].kindex, ", ", CS(&kname),
                " is not allowed in the Binary table.");
            wrterr(out, &errmes, 1);
            j += 1;
        }
    }

    /* check whether the column name is unique */
    test_colnam(out, hduptr);
}

/*************************************************************
*
*      test_header
*
*   Test the general keywords that can be in any header
*
*************************************************************/
fn test_header(
    _infits: &mut fitsfile, /* input fits file   */
    out: Out,               /* output ascii file */
    hduptr: &mut FitsHdu,   /* information about header  */
) {
    /* common mandatory  keywords */
    let mandkey: [&std::ffi::CStr; 5] = [c"SIMPLE", c"BITPIX", c"NAXIS", c"XTENSION", c"END"];
    let nmandkey = 5;

    /* string keywords */
    let strkey: [&std::ffi::CStr; 9] = [
        c"EXTNAME", c"ORIGIN", c"AUTHOR", c"CREATOR", c"REFERENC", c"TELESCOP", c"INSTRUME",
        c"OBSERVER", c"OBJECT",
    ];
    let nstrkey = 9;

    /* int keywords  */
    let intkey: [&std::ffi::CStr; 2] = [c"EXTVER", c"EXTLEVEL"];
    let nintkey = 2;

    /* floating keywords  */
    let fltkey: [&std::ffi::CStr; 3] = [c"EQUINOX", c"MJD-OBS", c"MJD-AVG"];
    let nfltkey = 3;

    let mut j: c_int;
    let mut k: c_int = 0;
    let mut n: c_int = 0;
    let mut lv: c_long;
    let mut stat: rsfitsio::c_types::c_ulong = 0;
    let mut vtemp: [c_char; 72] = [0; 72];
    let mut status: c_int = 0;
    let mut yr: c_int = 0;
    let mut mn: c_int = 0;
    let mut dy: c_int = 0;
    let mut hr: c_int = 0;
    let mut min: c_int = 0; /* time */
    let mut sec: c_double = 0.0;
    let mut yy: c_int;
    let mut ktype = kwdtyp::UNKNOWN;
    let mut errmes: [c_char; ERRMES_LEN] = [0; ERRMES_LEN];

    let numusrkey = hduptr.tkeys;
    let tmpkwds = tmpkwds_get();

    /* Check the mandatory keywords */
    for i in 0..nmandkey {
        let pp: &[c_char] = cast_slice(mandkey[i].to_bytes_with_nul());
        key_match(&tmpkwds, numusrkey, pp, 1, &mut k, &mut n);
        if k > -1 {
            j = k;
            while j < k + n {
                spf!(errmes;
                    "Keyword #", hduptr.kwds[j as usize].kindex, ", ",
                    CS(&hduptr.kwds[j as usize].kname), " is duplicated or out of order.");
                wrterr(out, &errmes, 1);
                j += 1;
            }
        }
    }

    /* check the NAXIS index keyword */
    key_match(&tmpkwds, numusrkey, cs!(c"NAXIS"), 0, &mut k, &mut n);
    j = k;
    while j < k + n {
        let kindex = hduptr.kwds[j as usize].kindex;
        let kname = hduptr.kwds[j as usize].kname;
        lv = strtol_c(&kname[5..]).0;
        if lv > 0 {
            if kindex as c_long != 3 + lv {
                spf!(errmes;
                    "Keyword #", kindex, ", ", CS(&kname), " is duplicated or out of order.");
                wrterr(out, &errmes, 1);
            }
            if lv > hduptr.naxis as c_long {
                spf!(errmes;
                    "Keyword #", kindex, ", ", CS(&kname), " is not allowed (with n > NAXIS =",
                    hduptr.naxis, ").");
                wrterr(out, &errmes, 1);
            }
        }
        j += 1;
    }

    /* Check the deprecated keywords */
    key_match(&tmpkwds, numusrkey, cs!(c"EPOCH"), 1, &mut k, &mut n);
    if k > -1 {
        spf!(errmes;
            "Keyword #", hduptr.kwds[k as usize].kindex, ", ",
            CS(&hduptr.kwds[k as usize].kname), " is deprecated. Use EQUINOX instead.");
        wrtwrn(out, &errmes, 0);
        check_flt(&hduptr.kwds[k as usize], out);
    }

    /* Check the DATExxxx keyword */
    key_match(&tmpkwds, numusrkey, cs!(c"DATE"), 0, &mut k, &mut n);
    j = k;
    while j < n + k {
        let kindex = hduptr.kwds[j as usize].kindex;
        let kname = hduptr.kwds[j as usize].kname;
        let kvalue = hduptr.kwds[j as usize].kvalue;
        check_str(&hduptr.kwds[j as usize], out);
        if fits_str2time(
            Some(&kvalue),
            Some(&mut yr),
            Some(&mut mn),
            Some(&mut dy),
            Some(&mut hr),
            Some(&mut min),
            Some(&mut sec),
            &mut status,
        ) != 0
        {
            spf!(errmes; "Keyword #", kindex, ", ", CS(&kname), ": ");
            wrtserr(out, &errmes, &mut status, 1);
        }
        if let Some(o) = chr_pos(&kvalue, b'/') {
            let pt = o + 4;
            yy = strtol_c(&kvalue[pt..]).0 as c_int;
            if (0..=10).contains(&yy) {
                spf!(errmes;
                    "Keyword #", kindex, ", ", CS(&kname), " ", CS(&kvalue),
                    " intends to mean year 20", SW(&format!("{yy:02}"), -2, None), "?");
                wrtwrn(out, &errmes, 0);
            }
        }
        j += 1;
    }

    /* Check the reserved string keywords */
    for i in 0..nstrkey {
        let temp: &[c_char] = cast_slice(strkey[i].to_bytes_with_nul());
        key_match(&tmpkwds, numusrkey, temp, 1, &mut k, &mut n);
        if k > -1 {
            check_str(&hduptr.kwds[k as usize], out);
        }
    }

    /* Check the reserved int keywords */
    for i in 0..nintkey {
        let temp: &[c_char] = cast_slice(intkey[i].to_bytes_with_nul());
        key_match(&tmpkwds, numusrkey, temp, 1, &mut k, &mut n);
        if k > -1 {
            check_int(&hduptr.kwds[k as usize], out);
        }
    }

    /* Check  reserved floating  keywords */
    for i in 0..nfltkey {
        let temp: &[c_char] = cast_slice(fltkey[i].to_bytes_with_nul());
        key_match(&tmpkwds, numusrkey, temp, 1, &mut k, &mut n);
        if k > -1 {
            check_flt(&hduptr.kwds[k as usize], out);
        }
    }

    /* Check the duplication of the keywords */
    for i in 0..(numusrkey - 1).max(0) as usize {
        if cbytes(&tmpkwds[i]) == cbytes(&tmpkwds[i + 1])
            && !ceq(&hduptr.kwds[i].kname, b"HIERARCH")
        {
            spf!(errmes;
                "Keyword ", CS(&hduptr.kwds[i].kname), " is duplicated in card #",
                hduptr.kwds[i].kindex, " and card #", hduptr.kwds[i + 1].kindex, ".");
            wrtwrn(out, &errmes, 0);
        }
    }

    /* check the long string convention */
    if hduptr.use_longstr == 1 {
        key_match(&tmpkwds, numusrkey, cs!(c"LONGSTRN"), 1, &mut k, &mut n);
        if k <= -1 {
            spf!(errmes;
                "The OGIP long string keyword convention is used without the recommended LONGSTRN keyword. ");
            wrtwrn(out, &errmes, 1);
        }
    }

    /* Check the HIERARCH keywords */
    if testhierarch() != 0 {
        key_match(&tmpkwds, numusrkey, cs!(c"HIERARCH"), 1, &mut k, &mut n);
        j = k;
        while j < n + k {
            let i = (hduptr.kwds[j as usize].kindex - 1) as usize; /* index number of the keyword */
            let cardi = card_at(i);

            /* Must have a space character following "HIERARCH" */
            if cardi[8] as u8 != b' ' {
                spf!(errmes;
                    "Keyword #", i + 1,
                    ": does not have a space character in byte 9:               ",
                    CSW(&cardi, 66, None));
                wrterr(out, &errmes, 1);
            }

            /* Whether the characters in HIERARCH token names are valid */
            let mut pt = 0usize;
            while cardi[pt] != 0 {
                if cardi[pt] as u8 == b'=' {
                    /* look for the required "=" sign */
                    break;
                }

                let c = cardi[pt] as u8;
                if !c.is_ascii_uppercase()
                    && !c.is_ascii_digit()
                    && c != b'-'
                    && c != b'_'
                    && c != b' '
                {
                    spf!(errmes;
                        "Keyword #", i + 1, ": token contains illegal char \"", CHR(cardi[pt]),
                        "\" (only A-Z,0-9,-,_): ", CSW(&cardi, 66, None));
                    wrterr(out, &errmes, 1);
                    break;
                }

                pt += 1;
            }

            if cardi[pt] == 0 {
                /* the "=" is not present on the card */
                spf!(errmes;
                    "Keyword #", i + 1,
                    ": does not contain the \"=\" value indicator char:      ",
                    CSW(&cardi, 66, None));
                wrterr(out, &errmes, 1);
            } else if cardi[pt] as u8 == b'=' {
                /* now check that the keyword has a legal value */
                pt += 1;
                while isspace_c(cardi[pt]) && cardi[pt] != 0 {
                    pt += 1;
                }
                match cardi[pt] as u8 {
                    b'\'' => get_str(&cardi, &mut pt, &mut vtemp, &mut stat), /* string */
                    b'T' | b'F' => get_log(&cardi, &mut pt, &mut vtemp, &mut stat), /*logical */
                    /* number */
                    b'+' | b'-' | b'.' | b'0'..=b'9' => {
                        get_num(&cardi, &mut pt, &mut vtemp, &mut ktype, &mut stat)
                    }
                    b'(' => get_cmp(&cardi, &mut pt, &mut vtemp, &mut ktype, &mut stat), /* complex */
                    b'/' => {} /* comment */
                    _ => get_unknown(&cardi, &mut pt, &mut vtemp, &mut ktype, &mut stat),
                }

                pr_kval_err(out, i as c_int + 1, cs!(c"HIERARCH"), &vtemp, stat);
                stat = 0; /* reset error status for next time */
            } /* end of keyword value test */
            j += 1;
        } /* end of test of individual HIERARCH keywords */

        /* now test for any duplicate HIERARCH keywords */
        j = k;
        while j < n + k - 1 {
            /* loop over all the HIERARCH keywords except the last */

            let i = (hduptr.kwds[j as usize].kindex - 1) as usize;
            let cardi = card_at(i);
            if chr_pos(&cardi, b'=').is_some() {
                let mut jj = j + 1;
                while jj < n + k {
                    /* loop over any other HIERARCH keywords */
                    let ii = (hduptr.kwds[jj as usize].kindex - 1) as usize;
                    let cardii = card_at(ii);
                    if chr_pos(&cardi, b'=').is_some() {
                        /* compare names char by char, ignoring extra spaces */
                        let mut p1 = 8usize; /*start at the end of the HIERARCH name */
                        let mut p2 = 8usize;

                        while cardi[p1] == cardii[p2] {
                            /* chars are the same in both */
                            if cardi[p1] as u8 == b' ' {
                                /* this is a space char */

                                /* skip over non-significant spaces in both name,
                                then continue testing next chars */
                                while cardi[p1] as u8 == b' ' {
                                    p1 += 1;
                                }
                                while cardii[p2] as u8 == b' ' {
                                    p2 += 1;
                                }
                            } else if cardi[p1] as u8 == b'=' {
                                /* found '=' in both keywords */
                                spf!(errmes;
                                    "HIERARCH keyword name is duplicated in cards #", i + 1,
                                    " and card #", ii + 1, ":   ", CSW(&cardi, 66, None));
                                wrtwrn(out, &errmes, 0);
                                p1 += 1; /* do this to prevent duplicate warning message, below */
                                break;
                            } else {
                                p1 += 1;
                                p2 += 1;
                            }
                        } /* end of identical chars test */

                        /* chars are not identical */
                        /* test for special case where one is a '=' and other
                        non-significant spaces followed by '=' */
                        /* first, skip over spaces in either name */
                        while cardi[p1] as u8 == b' ' {
                            p1 += 1;
                        }
                        while cardii[p2] as u8 == b' ' {
                            p2 += 1;
                        }

                        if cardi[p1] as u8 == b'=' && cardii[p2] as u8 == b'=' {
                            /* if we got here, then the names are the same except for
                            non-significant spacing differences */
                            spf!(errmes;
                                "HIERARCH keyword name is duplicated in cards #", i + 1,
                                " and card #", ii + 1, ":   ", CSW(&cardi, 66, None));
                            wrtwrn(out, &errmes, 0);
                            break;
                        } else {
                            /* these HIERARCH keywords do not have idential tokens */
                            break;
                        }
                    } /* end of if second HIERARCH keyword has a '=' */
                    #[allow(unreachable_code)]
                    {
                        jj += 1;
                    }
                } /* end of loop over other HIERARCH keywords */
            } /* end of if first HIERARCH keyword has a '=' */
            j += 1;
        } /* end of loop over all HIERARCH keywords */
    }

    /*  disabled this routine because it doesn't perform an useful tests
        test_img_wcs(infits, out, hduptr);
    */
}

/*************************************************************
*
*      key_match
*
*   find the keywords whose name match the pattern. The keywords
*   name is stored in a sorted array.
*
*************************************************************/
fn key_match(
    strs: &[[c_char; FLEN_KEYWORD]], /* fits keyname  array */
    nstr: c_int,                     /* total number of keys */
    pattern: &[c_char],              /* wanted pattern  */
    exact: c_int,                    /* exact matching or pattern matching
                                     exact = 1: exact matching.
                                     exact = 0: pattern matching.
                                     Any keywords with "patten"* is included
                                     */
    ikey: &mut c_int, /* The element number of first key
                      Return -99 if not found */
    mkey: &mut c_int, /* total number of key matched
                      return -999 if not found */
) {
    let mut i: c_int;
    let fnpt: fn(&[c_char], &[c_char]) -> Ordering = if exact != 0 { compstre } else { compstrp };
    *mkey = -999;
    *ikey = -99;

    /* bsearch(): glibc's binary search, so that the element it lands on for a
       run of equal keys matches the C exactly. */
    let mut l: c_int = 0;
    let mut u: c_int = nstr;
    let mut found: Option<c_int> = None;
    while l < u {
        let idx = (l + u) / 2;
        match fnpt(pattern, &strs[idx as usize]) {
            Ordering::Less => u = idx,
            Ordering::Greater => l = idx + 1,
            Ordering::Equal => {
                found = Some(idx);
                break;
            }
        }
    }

    if let Some(pos) = found {
        *mkey = 1;
        *ikey = pos;
        i = *ikey - 1;
        let mut p = pos - 1;
        /* NB: `i > 0', not `i >= 0' -- the C never revisits index 0, so a run
           of matches starting there is only partly reported.  Kept as-is. */
        while i > 0 && fnpt(pattern, &strs[p as usize]) == Ordering::Equal {
            *mkey += 1;
            *ikey = i;
            i -= 1;
            p -= 1;
        }
        i = *ikey + *mkey;
        let mut p = pos + 1;
        while i < nstr && fnpt(pattern, &strs[p as usize]) == Ordering::Equal {
            *mkey += 1;
            i += 1;
            p += 1;
        }
    }
}

/*************************************************************
*
*      test_colnam
*
*   Test the whether the column name is unique.
*
*************************************************************/
fn test_colnam(out: Out, hduptr: &FitsHdu) {
    let n: c_int;
    let mut errmes: [c_char; ERRMES_LEN] = [0; ERRMES_LEN];

    n = hduptr.ncols;

    if n <= 0 {
        return;
    }
    /* make a local working copy of ttype */
    let ttype: Vec<Vec<c_char>> = (0..n as usize).map(ttype_at).collect();
    let mut ttypecopy: Vec<Vec<c_char>> = Vec::with_capacity(n as usize);
    for i in 0..n as usize {
        let mut b: Vec<c_char> = vec![0; FLEN_VALUE];
        strcpy_c(&mut b, &ttype[i]);
        ttypecopy.push(b);
    }

    /* check whether there are any other non ASCII-text characters
      (FITS standard R14). Also "uppercase" the working copies. */
    for i in 0..n as usize {
        if cstrlen(&ttype[i]) == 0 {
            spf!(errmes;
                "Column #", i + 1, " has no name (No TTYPE", i + 1, " keyword).");
            wrtwrn(out, &errmes, 0);
            continue;
        }

        /*      disable this check (it was only a warning) */

        let mut p = 0usize;
        while ttype[i][p] != 0 {
            let c = ttype[i][p] as u8;
            if !c.is_ascii_lowercase()
                && !c.is_ascii_uppercase()
                && !c.is_ascii_digit()
                && c != b'_'
            {
                if c == b'&' {
                    spf!(errmes;
                        "Column #", i + 1, ": Reserved column name keyword (TTYPE", i + 1,
                        ") may use an illegal CONTINUE ('", CHR(ttype[i][p]), "')");
                    wrtwrn(out, &errmes, 0);
                } else {
                    spf!(errmes;
                        "Column #", i + 1, ": Name \"", CS(&ttype[i]), "\" contains character '",
                        CHR(ttype[i][p]), "' other than letters, digits, and \"_\".");
                    wrtwrn(out, &errmes, 0);
                }
            }
            ttypecopy[i][p] = toupper_c(ttype[i][p]);
            p += 1;
        }
    }

    let mut cols: Vec<ColName> = Vec::with_capacity(n as usize);
    for i in 0..n as usize {
        cols.push(ColName {
            name: ttypecopy[i].clone(),
            index: i as c_int + 1,
        });
    }

    /* sort the column name in the ascending order of name field*/
    cols.sort_by(compcol);

    /* Check the duplication of the column name */
    for i in 0..(n - 1).max(0) as usize {
        if cstrlen(&cols[i].name) == 0 {
            continue;
        }

        if cbytes(&cols[i].name) == cbytes(&cols[i + 1].name) {
            spf!(errmes;
                "Columns #", cols[i].index, ", ", CS(&ttype[(cols[i].index - 1) as usize]),
                " and #", cols[i + 1].index, ", ",
                CS(&ttype[(cols[i + 1].index - 1) as usize]),
                " are not unique (case insensitive).");
            wrtwrn(out, &errmes, 0);
        }
    }
}

/*************************************************************
*
*     parse_vtform
*
*   Parse the tform of the variable length vector.
*
*************************************************************/
pub(crate) fn parse_vtform(
    infits: &mut fitsfile,
    out: Out,
    _hduptr: &FitsHdu,
    colnum: c_int,        /* column number */
    datacode: &mut c_int, /* data code */
    maxlen: &mut c_long,  /* maximum length of the vector */
    isQFormat: &mut c_int, /* true if var col is 'Q' format */
) {
    let mut i: c_int = 0;
    let mut status: c_int = 0;
    let mut p: usize;
    let mut errmes: [c_char; ERRMES_LEN] = [0; ERRMES_LEN];

    *maxlen = -1;
    let temp = tform_at((colnum - 1) as usize);
    p = 0;

    if isdigit_c(temp[p]) {
        i = strtol_c(&temp).0 as c_int; /* sscanf(ptemp, "%d", &i) */
    }
    if i > 1 {
        spf!(errmes;
            "Illegal repeat value for value ", CS(&temp), " of TFORM", colnum, ".");
        wrterr(out, &errmes, 1);
    }
    while isdigit_c(temp[p]) {
        p += 1;
    }

    if temp[p] as u8 != b'P' && temp[p] as u8 != b'Q' {
        spf!(errmes;
            "TFORM", colnum, " is not for the variable length array: ", CS(&temp), ".");
        wrterr(out, &errmes, 1);
    }
    *isQFormat = if temp[p] as u8 == b'Q' { 1 } else { 0 };

    fits_get_coltype(infits, colnum, Some(datacode), None, None, &mut status);
    status = 0;
    let _ = status;
    p += 2;
    if p >= temp.len() || temp[p] as u8 != b'(' {
        return;
    }
    p += 1;
    if !isdigit_c(temp[p]) {
        spf!(errmes; "Bad value of TFORM", colnum, ": ", CS(&temp), ".");
        wrterr(out, &errmes, 1);
    }
    *maxlen = strtol_c(&temp[p..]).0; /* sscanf(p,"%ld",maxlen) */
    while isdigit_c(temp[p]) {
        p += 1;
    }
    if temp[p] as u8 != b')' {
        spf!(errmes; "Bad value of TFORM", colnum, ": ", CS(&temp), ".");
        wrterr(out, &errmes, 1);
    }
}

/*************************************************************
*
*     parse_wcskey_suffix
*
*   Retrieve the axis number and alternative coordinate suffix
*    from WCS keywords.  Return -1 if keyword does not match the
*    expected format.
*
*************************************************************/
fn parse_wcskey_suffix(
    fullname: &[c_char],
    rootname: &[c_char],
    axis: &mut c_int,
    alt: &mut c_int,
) -> c_int {
    let mut status: c_int = 0;
    let testAxis: c_int;

    *axis = 0;
    *alt = 0;
    if cstrlen(fullname) < cstrlen(rootname) {
        status = -1;
    } else {
        let suffix = cstrlen(rootname);
        if cstrlen(&fullname[suffix..]) != 0 {
            let (v, end) = strtol_c(&fullname[suffix..]);
            testAxis = v as c_int;
            let testAlt = suffix + end;
            if testAxis != 0 {
                *axis = testAxis;
            }
            if cstrlen(&fullname[testAlt..]) == 1 {
                let c = fullname[testAlt] as u8;
                if c.is_ascii_uppercase() {
                    *alt = c as c_int - 64;
                } else {
                    status = -1;
                }
            } else if cstrlen(&fullname[testAlt..]) > 1 {
                status = -1;
            }
        }
    }

    status
}

/*************************************************************
*
*      print_title
*
*  Print the title of the HDU.
*  when verbose < 2, called by wrterr and wrtwrn.
*
*************************************************************/
fn print_title(out: Out, hdunum: c_int, hdutype: c_int) {
    let mut hdutitle: [c_char; 64] = [0; 64];

    /* print out the title */
    CURHDU.with(|c| c.set(hdunum));
    CURTYPE.with(|c| c.set(hdutype));
    let curhdu = hdunum;
    let curtype = hdutype;

    if OLDHDU.with(Cell::get) == curhdu {
        return;
    } /* Do not print it twice */
    if curhdu == 1 {
        spf!(hdutitle; " HDU ", curhdu, ": Primary Array ");
    } else {
        match curtype {
            IMAGE_HDU => spf!(hdutitle; " HDU ", curhdu, ": Image Exten. "),
            ASCII_TBL => spf!(hdutitle; " HDU ", curhdu, ": ASCII Table "),
            BINARY_TBL => spf!(hdutitle; " HDU ", curhdu, ": BINARY Table "),
            _ => spf!(hdutitle; " HDU ", curhdu, ": Unknown Ext. "),
        }
    }
    wrtsep(out, b'=' as c_char, &hdutitle, 60);
    wrtout_str(out, " ");
    OLDHDU.with(|c| c.set(curhdu));
    if curhdu == totalhdu() {
        OLDHDU.with(|c| c.set(0)); /* reset the old hdu at the last hdu */
    }
}

/*************************************************************
*
*      print_header
*
*  Print the header of the HDU.
*
*************************************************************/
fn print_header(out: Out) {
    let mut htemp: [c_char; 100] = [0; 100];
    let ncards = NCARDS.with(Cell::get);
    for i in 1..=ncards {
        spf!(htemp; DW(i, 4), " | ", CS(&card_at(i as usize - 1)));
        wrtout(out, &htemp);
    }
    wrtout_str(out, " ");
}

/*************************************************************
*
*      print_summary
*
*  Print out the summary of this hdu.
*
**************************************************************/
fn print_summary(
    _infits: &mut fitsfile, /* input fits file   */
    out: Out,               /* output ascii file */
    hduptr: &FitsHdu,
) {
    let mut extver: [c_char; 10] = [0; 10];
    let mut extnv: [c_char; 2 * FLEN_VALUE + 4] = [0; 2 * FLEN_VALUE + 4];
    let mut npix: c_long;
    let hdutype: c_int;
    let mut temp: [c_char; 80] = [0; 80];
    let mut comm: [c_char; COMM_LEN] = [0; COMM_LEN];

    /* get the error number and wrn number */
    set_hduerr(hduptr.hdunum);

    hdutype = hduptr.hdutype;
    spf!(comm; " ", hduptr.nkeys, " header keywords");
    wrtout(out, &comm);
    wrtout_str(out, " ");
    if hdutype == ASCII_TBL || hdutype == BINARY_TBL {
        spf!(extnv; CS(&hduptr.extname));
        if hduptr.extver != -999 {
            spf!(extver; "(", hduptr.extver, ")");
            scat!(extnv; CS(&extver));
        }

        spf!(comm;
            " ", CS(&extnv), "  (", hduptr.ncols, " columns x ",
            hduptr.naxes.get(1).copied().unwrap_or(0) as i64, " rows)");
        wrtout(out, &comm);
        if hduptr.ncols != 0 {
            wrtout_str(out, " ");
            spf!(comm; " Col# Name (Units)       Format");
            wrtout(out, &comm);
        }
        for i in 0..hduptr.ncols as usize {
            let tt = ttype_at(i);
            let tu = tunit_at(i);
            let tf = tform_at(i);
            if cstrlen(&tu) != 0 {
                spf!(extnv; CS(&tt), " (", CS(&tu), ")");
            } else {
                spf!(extnv; CS(&tt));
            }
            spf!(comm;
                " ", DW(i + 1, 3), " ", CSW(&extnv, -20, Some(20)), " ",
                CSW(&tf, -10, Some(10)));
            wrtout(out, &comm);
        }
    } else if hdutype == IMAGE_HDU && hduptr.isgroup != 0 {
        spf!(comm; " ", hduptr.gcount, " Random Groups, ");

        bitpix_text(hduptr.bitpix, &mut temp);
        scat!(comm; CS(&temp));

        spf!(temp; " ", hduptr.naxis, " axes ");
        scat!(comm; CS(&temp));

        spf!(temp; "(", hduptr.naxes[0] as i64);
        scat!(comm; CS(&temp));

        npix = hduptr.naxes[0] as c_long;
        for i in 1..hduptr.naxis as usize {
            npix *= hduptr.naxes[i] as c_long;
            spf!(temp; " x ", hduptr.naxes[i] as i64);
            scat!(comm; CS(&temp));
        }
        let _ = npix;
        scat!(comm; "), ");
        wrtout(out, &comm);
    } else if hdutype == IMAGE_HDU {
        if hduptr.naxis > 0 {
            if hduptr.hdunum == 1 {
                extnv[0] = 0;
            } else {
                spf!(extnv; CS(&hduptr.extname));
                if hduptr.extver != -999 {
                    spf!(extver; " (", hduptr.extver, ")");
                    scat!(extnv; CS(&extver));
                }
            }
            spf!(comm; CS(&extnv));

            bitpix_text(hduptr.bitpix, &mut temp);
            scat!(comm; CS(&temp));

            spf!(temp; " ", hduptr.naxis, " axes ");
            scat!(comm; CS(&temp));

            spf!(temp; "(", hduptr.naxes[0] as i64);
            scat!(comm; CS(&temp));

            npix = hduptr.naxes[0] as c_long;
            for i in 1..hduptr.naxis as usize {
                npix *= hduptr.naxes[i] as c_long;
                spf!(temp; " x ", hduptr.naxes[i] as i64);
                scat!(comm; CS(&temp));
            }
            let _ = npix;
            scat!(comm; "), ");
            wrtout(out, &comm);
        } else {
            spf!(comm; " Null data array; NAXIS = 0 ");
            wrtout(out, &comm);
        }
    }
    wrtout_str(out, " ");
}

fn bitpix_text(bitpix: c_int, temp: &mut [c_char]) {
    match bitpix {
        BYTE_IMG => spf!(temp; " 8-bit integer pixels, "),
        SHORT_IMG => spf!(temp; " 16-bit integer pixels, "),
        USHORT_IMG => spf!(temp; " 16-bit unsigned integer pixels, "),
        LONG_IMG => spf!(temp; " 32-bit integer pixels, "),
        LONGLONG_IMG => spf!(temp; " 64-bit long integer pixels, "),
        ULONG_IMG => spf!(temp; " 32-bit unsigned integer pixels, "),
        FLOAT_IMG => spf!(temp; " 32-bit floating point pixels, "),
        DOUBLE_IMG => spf!(temp; " 64-bit double precision pixels, "),
        _ => spf!(temp; " unknown datatype, "),
    }
}

/*************************************************************
*
*      close_hdu
*
*  Free the memory allocated to the FitsHdu structure and
*  other temporary  spaces.
*
**************************************************************/
fn close_hdu(hduptr: &mut FitsHdu) {
    /* The C free()s cards, hduptr->kwds, datamin/datamax/tnull, naxes and
       tmpkwds here; RAII does that for us.  (Note the C's
       `if(hdutype == ASCII_TBL && hdutype == BINARY_TBL)' guard around the
       ttype/tform/tunit frees can never be true, so those are leaked there
       and simply left in place here.) */
    hduptr.kwds.clear();
    hduptr.datamin.clear();
    hduptr.datamax.clear();
    hduptr.tnull.clear();
    hduptr.naxes.clear();
    CARDS.with(|c| c.borrow_mut().clear());
    TMPKWDS.with(|t| t.borrow_mut().clear());
}

/*===========================================================================
 *  local helpers
 *==========================================================================*/


#[allow(dead_code)]
fn _unused_api() {
    let _ = NullValue::Int(0);
    let _ = cards_len();
}

#[cfg(test)]
mod tests {
    use super::*;

    fn kw(s: &str) -> [c_char; FLEN_KEYWORD] {
        let mut b = [0 as c_char; FLEN_KEYWORD];
        for (i, &c) in s.as_bytes().iter().enumerate() {
            b[i] = c as c_char;
        }
        b
    }

    #[test]
    fn test_key_match_exact() {
        let strs = [kw("BITPIX"), kw("CRPIX1"), kw("CRPIX2"), kw("CRVAL1"), kw("NAXIS")];
        let (mut k, mut n) = (0, 0);

        key_match(&strs, 5, cs!(c"NAXIS"), 1, &mut k, &mut n);
        assert_eq!((k, n), (4, 1));

        key_match(&strs, 5, cs!(c"BITPIX"), 1, &mut k, &mut n);
        assert_eq!((k, n), (0, 1));

        /* not found: the C's sentinels, which the `for (j = k; j < k+n; j++)'
           loops rely on to iterate zero times */
        key_match(&strs, 5, cs!(c"THEAP"), 1, &mut k, &mut n);
        assert_eq!((k, n), (-99, -999));
        assert!(k + n < k);

        /* exact matching must not match a longer keyword */
        key_match(&strs, 5, cs!(c"CRPIX"), 1, &mut k, &mut n);
        assert_eq!((k, n), (-99, -999));
    }

    #[test]
    fn test_key_match_prefix() {
        let strs = [kw("BITPIX"), kw("CRPIX1"), kw("CRPIX2"), kw("CRVAL1"), kw("NAXIS")];
        let (mut k, mut n) = (0, 0);

        key_match(&strs, 5, cs!(c"CRPIX"), 0, &mut k, &mut n);
        assert_eq!((k, n), (1, 2));

        key_match(&strs, 5, cs!(c"CR"), 0, &mut k, &mut n);
        assert_eq!((k, n), (1, 3));

        key_match(&strs, 5, cs!(c"TC"), 0, &mut k, &mut n);
        assert_eq!((k, n), (-99, -999));
    }

    /* The C's backward scan is guarded by `i > 0', so a run of matches that
       begins at index 0 is only partly reported when bsearch lands above it.
       Locked in here because the port has to reproduce it. */
    #[test]
    fn test_key_match_index_zero_quirk() {
        let strs = [kw("CRPIX1"), kw("CRPIX2"), kw("NAXIS")];
        let (mut k, mut n) = (0, 0);
        key_match(&strs, 3, cs!(c"CRPIX"), 0, &mut k, &mut n);
        /* bsearch lands on index 1; index 0 is never revisited */
        assert_eq!((k, n), (1, 1));
    }

    #[test]
    fn test_parse_wcskey_suffix() {
        let (mut axis, mut alt) = (0, 0);

        assert_eq!(parse_wcskey_suffix(&kw("CRPIX1"), cs!(c"CRPIX"), &mut axis, &mut alt), 0);
        assert_eq!((axis, alt), (1, 0));

        assert_eq!(parse_wcskey_suffix(&kw("CRPIX12"), cs!(c"CRPIX"), &mut axis, &mut alt), 0);
        assert_eq!((axis, alt), (12, 0));

        /* alternate WCS description: trailing A-Z maps to 1-26 */
        assert_eq!(parse_wcskey_suffix(&kw("CRPIX1A"), cs!(c"CRPIX"), &mut axis, &mut alt), 0);
        assert_eq!((axis, alt), (1, 1));
        assert_eq!(parse_wcskey_suffix(&kw("CRPIX2Z"), cs!(c"CRPIX"), &mut axis, &mut alt), 0);
        assert_eq!((axis, alt), (2, 26));

        /* more than one trailing char, or a non A-Z one, is rejected */
        assert_eq!(parse_wcskey_suffix(&kw("CRPIX1AB"), cs!(c"CRPIX"), &mut axis, &mut alt), -1);
        assert_eq!(parse_wcskey_suffix(&kw("CRPIX1a"), cs!(c"CRPIX"), &mut axis, &mut alt), -1);

        /* name shorter than the root */
        assert_eq!(parse_wcskey_suffix(&kw("CR"), cs!(c"CRPIX"), &mut axis, &mut alt), -1);

        /* no suffix at all leaves both at 0 */
        assert_eq!(parse_wcskey_suffix(&kw("CRPIX"), cs!(c"CRPIX"), &mut axis, &mut alt), 0);
        assert_eq!((axis, alt), (0, 0));
    }
}

/* File-level tests.

   Each fixture is a small hand-built FITS file that exercises one report path.
   The expected (errors, warnings) pairs were read off the C fitsverify acting
   as the oracle -- they are not recordings of what this implementation happens
   to produce.  Every fixture here was also checked byte-for-byte (stdout,
   stderr and exit status) against the C across all seven flag combinations. */
#[cfg(test)]
mod fits_tests {
    use super::*;
    use crate::fvrf_file::{get_total_err, get_total_warn};

    fn card(s: &str) -> String {
        format!("{s:<80.80}")
    }

    fn block(cards: &[&str]) -> String {
        let mut h: String = cards.iter().map(|c| card(c)).collect();
        h.push_str(&card("END"));
        while h.len() % 2880 != 0 {
            h.push(' ');
        }
        h
    }

    const PRIM: [&str; 3] = [
        "SIMPLE  =                    T",
        "BITPIX  =                    8",
        "NAXIS   =                    0",
    ];
    const IMG: [&str; 5] = [
        "SIMPLE  =                    T",
        "BITPIX  =                  -32",
        "NAXIS   =                    2",
        "NAXIS1  =                    4",
        "NAXIS2  =                    4",
    ];
    const PRIM_EXT: [&str; 4] = [
        "SIMPLE  =                    T",
        "BITPIX  =                    8",
        "NAXIS   =                    0",
        "EXTEND  =                    T",
    ];

    /* a primary header with no data, plus `extra' */
    fn pw(extra: &[&str]) -> Vec<u8> {
        let mut c: Vec<&str> = PRIM.to_vec();
        c.extend_from_slice(extra);
        block(&c).into_bytes()
    }

    /* a 4x4 float image, plus `extra' */
    fn img(extra: &[&str]) -> Vec<u8> {
        let mut c: Vec<&str> = IMG.to_vec();
        c.extend_from_slice(extra);
        let mut v = block(&c).into_bytes();
        v.extend(std::iter::repeat_n(0u8, 2880));
        v
    }

    /* a table extension carrying `data', plus `extra' header cards */
    fn tbl(kind: &str, nax1: usize, nax2: usize, extra: &[&str], data: &[u8]) -> Vec<u8> {
        let n1 = format!("NAXIS1  = {nax1:>20}");
        let n2 = format!("NAXIS2  = {nax2:>20}");
        let x = format!("XTENSION= '{kind}'");
        let mut c: Vec<&str> = vec![
            &x,
            "BITPIX  =                    8",
            "NAXIS   =                    2",
            &n1,
            &n2,
            "PCOUNT  =                    0",
            "GCOUNT  =                    1",
        ];
        c.extend_from_slice(extra);
        let mut v = block(&PRIM_EXT).into_bytes();
        v.extend(block(&c).into_bytes());
        let mut d = data.to_vec();
        d.resize(2880, 0);
        v.extend(d);
        v
    }

    /* Runs verify_fits over `bytes' and returns (errors, warnings).  Out::Null
       suppresses all reporting, so only the counters are under test.  The
       fitsverify counters are thread-locals, so these run in parallel safely. */
    fn counts(bytes: &[u8]) -> (c_int, c_int) {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("t.fits");
        std::fs::write(&path, bytes).unwrap();

        let mut buf = vec![0 as c_char; FLEN_FILENAME];
        set_cstr(&mut buf, path.to_str().unwrap().as_bytes());
        verify_fits(&mut buf, Out::Null);
        (get_total_err(), get_total_warn())
    }

    /*--- clean input reports nothing -------------------------------------*/

    #[test]
    fn test_clean_primary_array() {
        assert_eq!(counts(&pw(&[])), (0, 0));
    }

    /*--- warnings --------------------------------------------------------*/

    #[test]
    fn test_simple_false_warns() {
        let mut c: Vec<&str> = PRIM.to_vec();
        c[0] = "SIMPLE  =                    F";
        assert_eq!(counts(&block(&c).into_bytes()), (0, 1));
    }

    #[test]
    fn test_deprecated_keywords_warn() {
        assert_eq!(counts(&pw(&["EPOCH   =               1950.0"])), (0, 1));
        assert_eq!(counts(&pw(&["BLOCKED =                    T"])), (0, 1));
    }

    #[test]
    fn test_duplicate_keyword_warns() {
        assert_eq!(
            counts(&pw(&[
                "DUPKEY  =                    1",
                "DUPKEY  =                    2",
            ])),
            (0, 1)
        );
    }

    #[test]
    fn test_y2k_date_warns() {
        assert_eq!(counts(&pw(&["DATE    = '01/02/09'"])), (0, 1));
        /* a four-digit year is fine */
        assert_eq!(counts(&pw(&["DATE    = '2009-02-01'"])), (0, 0));
    }

    #[test]
    fn test_long_string_convention_needs_longstrn() {
        assert_eq!(
            counts(&pw(&["LONGKEY = 'abcdefghij&'", "CONTINUE  'klmnop'"])),
            (0, 1)
        );
        /* declaring the convention silences it */
        assert_eq!(
            counts(&pw(&[
                "LONGSTRN= 'OGIP 1.0'",
                "LONGKEY = 'abcdefghij&'",
                "CONTINUE  'klmnop'",
            ])),
            (0, 0)
        );
    }

    #[test]
    fn test_null_keyword_value_warns() {
        assert_eq!(counts(&pw(&["NULLKEY ="])), (0, 1));
    }

    #[test]
    fn test_bscale_zero_warns() {
        assert_eq!(counts(&img(&["BSCALE  =                  0.0"])), (0, 1));
    }

    /*--- keyword syntax errors -------------------------------------------*/

    #[test]
    fn test_lower_case_keyword_name_errors() {
        assert_eq!(counts(&pw(&["lowerkey=                    1"])), (1, 0));
    }

    #[test]
    fn test_unterminated_string_errors() {
        assert_eq!(counts(&pw(&["OBJECT  = 'unterminated"])), (1, 0));
    }

    #[test]
    fn test_lower_case_exponent_errors() {
        assert_eq!(counts(&pw(&["VALUE   =              1.25e+02"])), (1, 0));
        assert_eq!(counts(&pw(&["VALUE   =              1.25E+02"])), (0, 0));
    }

    /* The C dereferences a NULL pr_end here and dies with SIGSEGV; this port
       reports the malformed value instead.  A deliberate divergence. */
    #[test]
    fn test_complex_without_comma_is_reported_not_fatal() {
        let (e, _w) = counts(&pw(&["CPLX    = (1 2)"]));
        assert!(e > 0, "a complex value with no comma must be reported");
    }

    /*--- keywords in the wrong kind of HDU -------------------------------*/

    #[test]
    fn test_xtension_in_primary_errors() {
        assert_eq!(counts(&pw(&["XTENSION= 'IMAGE   '"])), (2, 0));
    }

    #[test]
    fn test_pscal_outside_random_groups_errors() {
        assert_eq!(counts(&pw(&["PSCAL1  =                  1.0"])), (1, 0));
    }

    #[test]
    fn test_blank_with_floating_point_data_errors() {
        assert_eq!(counts(&img(&["BLANK   =                   -1"])), (1, 0));
    }

    /*--- WCS -------------------------------------------------------------*/

    #[test]
    fn test_wcs_cdelt_zero_warns() {
        assert_eq!(
            counts(&img(&[
                "CDELT1  =                  0.0",
                "CDELT2  =                  1.0",
            ])),
            (0, 3)
        );
    }

    #[test]
    fn test_wcs_negative_crder_errors() {
        assert_eq!(counts(&img(&["CRDER1  =                 -1.0"])), (1, 3));
    }

    #[test]
    fn test_wcs_pc_and_cd_are_mutually_exclusive() {
        assert_eq!(
            counts(&img(&[
                "PC1_1   =                  1.0",
                "CD1_1   =                  1.0",
            ])),
            (1, 0)
        );
    }

    #[test]
    fn test_wcs_pc_and_crota2_are_mutually_exclusive() {
        assert_eq!(
            counts(&img(&[
                "PC1_1   =                  1.0",
                "CROTA2  =                  0.0",
            ])),
            (1, 3)
        );
    }

    #[test]
    fn test_wcs_non_allowed_frame_values_warn() {
        assert_eq!(counts(&img(&["RADESYS = 'BOGUS   '"])), (0, 1));
        assert_eq!(counts(&img(&["SPECSYS = 'NOPE    '"])), (0, 1));
    }

    /*--- tables ----------------------------------------------------------*/

    #[test]
    fn test_tdim_not_allowed_in_ascii_table() {
        assert_eq!(
            counts(&tbl(
                "TABLE   ",
                4,
                1,
                &[
                    "TFIELDS =                    1",
                    "TBCOL1  =                    1",
                    "TFORM1  = 'I4      '",
                    "TTYPE1  = 'C1      '",
                    "TDIM1   = '(4)     '",
                ],
                b"1234",
            )),
            (2, 0)
        );
    }

    #[test]
    fn test_tbcol_not_allowed_in_binary_table() {
        assert_eq!(
            counts(&tbl(
                "BINTABLE",
                4,
                1,
                &[
                    "TFIELDS =                    1",
                    "TFORM1  = '4A      '",
                    "TTYPE1  = 'C1      '",
                    "TBCOL1  =                    1",
                ],
                b"abcd",
            )),
            (1, 0)
        );
    }

    #[test]
    fn test_missing_column_name_warns() {
        assert_eq!(
            counts(&tbl(
                "BINTABLE",
                4,
                1,
                &["TFIELDS =                    1", "TFORM1  = '4A      '"],
                b"abcd",
            )),
            (0, 1)
        );
    }

    #[test]
    fn test_illegal_character_in_column_name_warns() {
        assert_eq!(
            counts(&tbl(
                "BINTABLE",
                4,
                1,
                &[
                    "TFIELDS =                    1",
                    "TFORM1  = '4A      '",
                    "TTYPE1  = 'a b     '",
                ],
                b"abcd",
            )),
            (0, 1)
        );
    }

    #[test]
    fn test_duplicate_column_names_warn() {
        /* the comparison is case insensitive */
        assert_eq!(
            counts(&tbl(
                "BINTABLE",
                4,
                1,
                &[
                    "TFIELDS =                    2",
                    "TFORM1  = '2A      '",
                    "TFORM2  = '2A      '",
                    "TTYPE1  = 'SAME    '",
                    "TTYPE2  = 'same    '",
                ],
                b"abcd",
            )),
            (0, 1)
        );
    }

    #[test]
    fn test_tform_leading_space_errors() {
        assert_eq!(
            counts(&tbl(
                "BINTABLE",
                4,
                1,
                &[
                    "TFIELDS =                    1",
                    "TFORM1  = ' 4A     '",
                    "TTYPE1  = 'C1      '",
                ],
                b"abcd",
            )),
            (1, 0)
        );
    }

    /* Every offending row is reported, not just the first: the C's `break'
       leaves the inner character scan only, and the row loop carries on. */
    #[test]
    fn test_non_ascii_reported_for_every_bad_row() {
        assert_eq!(
            counts(&tbl(
                "BINTABLE",
                4,
                3,
                &[
                    "TFIELDS =                    1",
                    "TFORM1  = '4A      '",
                    "TTYPE1  = 'TXT     '",
                ],
                b"a\xffbcabcd\xfeefg", /* rows 1 and 3 are bad, row 2 is clean */
            )),
            (2, 0)
        );
    }

    /*--- file structure --------------------------------------------------*/

    #[test]
    fn test_extraneous_bytes_after_last_hdu_error() {
        let mut v = pw(&[]);
        v.extend_from_slice(b"junkjunk");
        assert_eq!(counts(&v), (1, 0));
    }
}
