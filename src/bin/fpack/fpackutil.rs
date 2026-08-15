/* Transpiled from cfitsio/utilities/fpackutil.c
 *
 * FPACK utility routines
 * R. Seaman, NOAO & W. Pence, NASA/GSFC
 *
 * Shared by both binaries, as the C links fpackutil.o into both.  Nothing here
 * calls fp_usage/fp_help or fu_usage/fu_help, so funpack/main.rs can
 * #[path]-include it unchanged.
 *
 * DEVIATIONS from the C, each also marked at its site:
 *
 *   exit()      -> Err(FpExit(n)); see fpack_h.rs.
 *   temp files  tempfile::TempPath instead of the abort_fpack signal handler
 *               and its three globals, which are not ported.
 *   system()    funpack -Z streams through flate2 rather than `gzip -1'.
 *   bug fixes   Only where the C is undefined and Rust cannot reproduce it.
 *               Everything else upstream gets wrong is reproduced faithfully,
 *               marked NOTE (upstream bug N) at its site.
 */

use std::cell::{Cell, RefCell};
use std::fs::File;
use std::io::{BufWriter, Write as _};
use std::path::PathBuf;
use std::time::Instant;

use bytemuck::{cast_slice, cast_slice_mut};
use tempfile::TempPath;

use rsfitsio::aliases::rust_api::{
    fits_close_file, fits_compress_table, fits_copy_hdu, fits_copy_header, fits_create_file,
    fits_create_img, fits_delete_file, fits_delete_hdu, fits_file_name, fits_get_chksum,
    fits_get_compression_type, fits_get_hdu_num, fits_get_hdu_type, fits_get_hduaddr,
    fits_get_hduaddrll, fits_get_img_param, fits_get_num_cols, fits_get_num_hdus,
    fits_get_num_rowsll, fits_get_version, fits_img_compress, fits_img_decompress,
    fits_img_stats_float, fits_img_stats_int, fits_img_stats_short, fits_is_compressed_image,
    fits_movabs_hdu, fits_movnam_hdu, fits_movrel_hdu, fits_open_file, fits_read_img_int,
    fits_read_img_sht, fits_read_key, fits_read_keyword, fits_read_pix, fits_read_subset,
    fits_read_subset_flt, fits_read_subset_int, fits_read_subset_sht, fits_report_error,
    fits_set_bscale, fits_set_compression_type, fits_set_dither_offset, fits_set_hcomp_scale,
    fits_set_hcomp_smooth, fits_set_hdustruc, fits_set_lossy_int, fits_set_quantize_level,
    fits_set_quantize_method, fits_set_tile_dim, fits_uncompress_table, fits_update_key,
    fits_write_chksum, fits_write_img_int, fits_write_img_sht,
};
use rsfitsio::c_types::{c_char, c_int, c_long, c_short, c_ulong};
use rsfitsio::fitscore::ffrhdu_safe;
use rsfitsio::fitsio::{
    ASCII_TBL, BINARY_TBL, BYTE_IMG, DATA_COMPRESSION_ERR, DOUBLE_IMG, END_OF_FILE, FLEN_COMMENT,
    FLEN_FILENAME, FLEN_VALUE, FLOAT_IMG, FLOATNULLVALUE, GZIP_1, GZIP_2, HCOMPRESS_1, IMAGE_HDU,
    KEY_NO_EXIST, LONG_IMG, LONGLONG, LONGLONG_IMG, MAX_COMPRESS_DIM, MEMORY_ALLOCATION,
    NOCOMPRESS, NOT_IMAGE, PLIO_1, READONLY, RICE_1, SHORT_IMG, TBYTE, TDOUBLE, TFLOAT, TINT,
    TSHORT, fitsfile,
};
use rsfitsio::iraffits::fits_delete_iraf_file_safe as fits_delete_iraf_file;
use rsfitsio::wrappers::{
    strcat_safe, strchr_safe, strcmp_safe, strcpy_safe, strlen_safe, strncmp_safe,
};
use rsfitsio::{KeywordDatatype, KeywordDatatypeMut, NullValue, STDERR, bb, cs};

use crate::cfmt::{dbl, not_int};
use crate::fpack_h::*;

/* nearest integer function.  DEVIATION: C's conversion of an out-of-range
double is undefined; Rust's `as' saturates. */
#[allow(non_snake_case)]
fn NINT(x: f64) -> c_int {
    if x >= 0. {
        (x + 0.5) as c_int
    } else {
        (x - 0.5) as c_int
    }
}
#[allow(non_snake_case)]
fn NSHRT(x: f64) -> c_short {
    if x >= 0. {
        (x + 0.5) as c_short
    } else {
        (x - 0.5) as c_short
    }
}

thread_local! {
    /* FILE *outreport; */
    static OUTREPORT: RefCell<Option<BufWriter<File>>> = const { RefCell::new(None) };

    /* dimension of central image area to be sampled for test statistics */
    static XSAMPLE: Cell<c_int> = const { Cell::new(4100) };
    static YSAMPLE: Cell<c_int> = const { Cell::new(4100) };

    /* elapsed and CPU time; the C keeps scpu/ecpu plus startsec/startmilli */
    static SCPU: Cell<f64> = const { Cell::new(0.0) };
    static START: Cell<Option<Instant>> = const { Cell::new(None) };
}

/// C: `clock() / CLOCKTICKS`, in seconds.  Not a runtime `cfg!()` because libc
/// exposes `clock()` on windows but not on linux-gnu.
#[cfg(unix)]
fn cpu_seconds() -> f64 {
    let mut ts = libc::timespec {
        tv_sec: 0,
        tv_nsec: 0,
    };
    // SAFETY: `ts` is a fully initialised timespec that outlives the call.
    unsafe {
        libc::clock_gettime(libc::CLOCK_PROCESS_CPUTIME_ID, &mut ts);
    }
    ts.tv_sec as f64 + ts.tv_nsec as f64 / 1e9
}

#[cfg(windows)]
fn cpu_seconds() -> f64 {
    /* CLOCKS_PER_SEC is 1000 in the Microsoft CRT and libc does not export it */
    // SAFETY: clock() reads a process-global counter and has no preconditions.
    (unsafe { libc::clock() } as f64) / 1000.0
}

#[cfg(not(any(unix, windows)))]
fn cpu_seconds() -> f64 {
    0.0
}

/* ------------------------------------------------------------------------ *
 * small helpers the C gets from <stdio.h>/<string.h>
 * ------------------------------------------------------------------------ */

/// The bytes of a NUL-terminated buffer, without the NUL.
pub(crate) fn cbytes(s: &[c_char]) -> &[u8] {
    &cast_slice(s)[..strlen_safe(s)]
}

/// A `[c_char]` buffer as a filesystem path.
fn cpath(s: &[c_char]) -> PathBuf {
    PathBuf::from(String::from_utf8_lossy(cbytes(s)).into_owned())
}

/// C: `strncpy(dst, s, N-1); dst[N-1] = 0;`
fn strncpy_trunc(dst: &mut [c_char], s: &[c_char]) {
    let n = strlen_safe(s).min(dst.len() - 1);
    dst[..n].copy_from_slice(&s[..n]);
    dst[n] = 0;
}

/* argv as NUL-terminated buffers, so argv[iarg][0] transpiles unchanged */
pub(crate) type Argv = Vec<Vec<c_char>>;

pub(crate) fn c_argv() -> Argv {
    std::env::args()
        .map(|a| {
            let mut v: Vec<c_char> = cast_slice(a.as_bytes()).to_vec();
            v.push(0);
            v
        })
        .collect()
}

/*--------------------------------------------------------------------------*/
pub(crate) fn fp_msg(msg: &[c_char]) -> c_int {
    fp_msg_bytes(cbytes(msg))
}

/// `fp_msg` of a C string literal.  C: `fp_msg ("Error: ...\n");`
pub(crate) fn fp_msg_str(msg: &str) -> c_int {
    fp_msg_bytes(msg.as_bytes())
}

fn fp_msg_bytes(msg: &[u8]) -> c_int {
    /* printf ("%s", msg); */
    let stdout = std::io::stdout();
    let mut lock = stdout.lock();
    let _ = lock.write_all(msg);
    0
}

/// `printf(...)` of already-formatted text.
fn pf(text: &str) {
    fp_msg_bytes(text.as_bytes());
}

/// C: `fits_report_error (stderr, stat)`.
fn report_error(stat: c_int) {
    // SAFETY: `stderr` is libc's process-global FILE*.  It is only read here
    // and handed straight back to CFITSIO, which writes the error stack to it.
    fits_report_error(unsafe { STDERR!() }, stat);
}

/// `fprintf(outreport, ...)`
fn report(text: &str) {
    OUTREPORT.with_borrow_mut(|r| {
        if let Some(f) = r.as_mut() {
            let _ = f.write_all(text.as_bytes());
        }
    });
}

/*--------------------------------------------------------------------------*/
pub(crate) fn fp_noop() -> c_int {
    fp_msg_str("Input and output files are unchanged.\n");
    0
}
/*--------------------------------------------------------------------------*/
/// C: `fp_abort_output(...)`, which ends in `exit (stat)` -- hence the FpExit
/// return.  Both files are taken by value: one is closed, the other deleted.
pub(crate) fn fp_abort_output(
    infptr: Option<Box<fitsfile>>,
    mut outfptr: Option<Box<fitsfile>>,
    stat: c_int,
) -> FpExit {
    let mut status = 0;
    let mut hdunum = 0;

    if let Some(mut infptr) = infptr {
        /* DEVIATION (upstream bugs 2 and 3): the C copies into the 513-byte
        global `tempfilename', which is both too small and the name its signal
        handler needed. */
        let mut filename: [c_char; FLEN_FILENAME] = [0; FLEN_FILENAME];
        fits_file_name(&mut infptr, &mut filename, &mut status);
        fits_get_hdu_num(&mut infptr, &mut hdunum);

        fits_close_file(infptr, &mut status);

        pf(&format!(
            "Error processing file: {}\n",
            String::from_utf8_lossy(cbytes(&filename))
        ));
        pf(&format!("  in HDU number {hdunum}\n"));
    } else {
        fp_msg_str("Error: Unable to process input file\n");
    }
    report_error(stat);

    if outfptr.is_some() {
        fits_delete_file(&mut outfptr, &mut status);
        fp_msg_str("Input file is unchanged.\n");
    }

    FpExit(stat)
}

/*--------------------------------------------------------------------------*/
pub(crate) fn fp_version() -> c_int {
    let mut version: f32 = 0.0;

    fp_msg_str(FPACK_VERSION);
    fits_get_version(&mut version);
    /* snprintf(cfitsioversion, 40," CFITSIO version %5.3f", version); */
    pf(&format!(
        " CFITSIO version {}",
        dbl(c"%5.3f", f64::from(version))
    ));
    fp_msg_str("\n");
    0
}
/*--------------------------------------------------------------------------*/
pub(crate) fn fp_access(filename: &[c_char]) -> c_int {
    /* test if a file exists */

    match File::open(cpath(filename)) {
        Ok(_) => 0,
        Err(_) => -1,
    }
}
/*--------------------------------------------------------------------------*/
/// C: `fp_tmpnam(suffix, rootname, tmpnam)`, transpiled verbatim -- including
/// its check-then-use race -- because its names appear in the C's messages.
/// Callers wrap the result in a `TempPath`.
pub(crate) fn fp_tmpnam(suffix: &[c_char], rootname: &[c_char], tmpnam: &mut [c_char]) -> FpResult {
    /* create temporary file name */

    let maxtry: usize;
    let mut ii: usize;

    if strlen_safe(suffix) + strlen_safe(rootname) > SZ_STR - 5 {
        fp_msg_str("Error: filename is too long to create temporary file\n");
        return Err(FpExit(-1));
    }

    strcpy_safe(tmpnam, rootname); /* start with rootname */
    strcat_safe(tmpnam, suffix); /* append the suffix */

    maxtry = SZ_STR - strlen_safe(tmpnam) - 1;

    ii = 0;
    while ii < maxtry {
        if fp_access(tmpnam) != 0 {
            break; /* good, the file does not exist */
        }
        if strlen_safe(tmpnam) > SZ_STR - 2 {
            fp_msg_str("\nCould not create temporary file name:\n");
            fp_msg(tmpnam);
            fp_msg_str("\n");
            return Err(FpExit(-1));
        }
        strcat_safe(tmpnam, cs!(c"x")); /* append an x to the name, and try again */
        ii += 1;
    }

    if ii == maxtry {
        fp_msg_str("\nCould not create temporary file name:\n");
        fp_msg(tmpnam);
        fp_msg_str("\n");
        return Err(FpExit(-1));
    }

    Ok(())
}

/// `fp_tmpnam` plus the guard that replaces `abort_fpack`.
fn fp_tmpnam_guard(
    suffix: &[c_char],
    rootname: &[c_char],
    tmpnam: &mut [c_char],
) -> FpResult<TempPath> {
    fp_tmpnam(suffix, rootname, tmpnam)?;
    /* `try_from_path' only rejects a path containing an interior NUL, which
    a NUL-terminated buffer cannot produce */
    Ok(TempPath::try_from_path(cpath(tmpnam)).expect("temp file name has no interior NUL"))
}

/*--------------------------------------------------------------------------*/
pub(crate) fn fp_init(fpptr: &mut fpstate) -> c_int {
    let mut ii: usize;

    fpptr.comptype = RICE_1;
    fpptr.quantize_level = DEF_QLEVEL;
    fpptr.no_dither = 0;
    fpptr.dither_method = 1;
    fpptr.dither_offset = 0;
    fpptr.int_to_float = 0;

    /* thresholds when using the -i2f flag */
    fpptr.n3ratio = 2.0; /* minimum ratio of image noise sigma / q */
    fpptr.n3min = 6.; /* minimum noise sigma. */

    fpptr.scale = DEF_HCOMP_SCALE;
    fpptr.smooth = DEF_HCOMP_SMOOTH;
    fpptr.rescale_noise = DEF_RESCALE_NOISE;
    fpptr.ntile[0] = -1; /* -1 means extent of axis */

    ii = 1;
    while ii < MAX_COMPRESS_DIM {
        fpptr.ntile[ii] = 1;
        ii += 1;
    }

    fpptr.to_stdout = 0;
    fpptr.listonly = 0;
    fpptr.clobber = 0;
    fpptr.delete_input = 0;
    fpptr.do_not_prompt = 0;
    fpptr.do_checksums = 1;
    fpptr.do_gzip_file = 0;
    fpptr.do_tables = 0; /* this is intended for testing purposes  */
    fpptr.do_images = 1; /* can be turned off with -tableonly switch */
    fpptr.test_all = 0;
    fpptr.verbose = 0;

    fpptr.prefix[0] = 0;
    fpptr.extname[0] = 0;
    fpptr.delete_suffix = 0;
    fpptr.outfile[0] = 0;

    fpptr.firstfile = 1;

    /* magic number for initialization check, boolean for preflight
     */
    fpptr.initialized = FP_INIT_MAGIC;
    fpptr.preflight_checked = 0;
    0
}
/*--------------------------------------------------------------------------*/
pub(crate) fn fp_list(argc: c_int, argv: &Argv, fpvar: fpstate) -> FpResult<c_int> {
    let mut infits: [c_char; SZ_STR] = [0; SZ_STR];
    let mut hdunum: c_int = 0;
    let mut iarg: c_int;
    let mut stat: c_int = 0;
    let mut sizell: LONGLONG = 0;

    if fpvar.initialized != FP_INIT_MAGIC {
        fp_msg_str("Error: internal initialization error\n");
        return Err(FpExit(-1));
    }

    iarg = fpvar.firstfile;
    while iarg < argc {
        strncpy_trunc(&mut infits, &argv[iarg as usize]);

        if strchr_safe(&infits, bb(b'[')).is_some() || strchr_safe(&infits, bb(b']')).is_some() {
            fp_msg_str("Error: section/extension notation not supported: ");
            fp_msg(&infits);
            fp_msg_str("\n");
            return Err(FpExit(-1));
        }

        if fp_access(&infits) != 0 {
            fp_msg_str("Error: can't find or read input file ");
            fp_msg(&infits);
            fp_msg_str("\n");
            fp_noop();
            return Err(FpExit(-1));
        }

        let mut infptr: Option<Box<fitsfile>> = None;
        fits_open_file(&mut infptr, &infits, READONLY, &mut stat);
        if stat != 0 {
            report_error(stat);
            return Err(FpExit(stat));
        }

        /* move to the end of file, to get the total size in bytes */
        {
            let f = infptr.as_mut().unwrap();
            fits_get_num_hdus(f, &mut hdunum, &mut stat);
            fits_movabs_hdu(f, hdunum, None, &mut stat);
            fits_get_hduaddrll(f, None, None, Some(&mut sizell), &mut stat);
        }

        if stat != 0 {
            return Err(fp_abort_output(infptr, None, stat));
        }

        pf(&format!("# {} (", String::from_utf8_lossy(cbytes(&infits))));
        /* the C selects "%I64d"/"%lld"/"%ld" by platform; all three print the
        same digits for a LONGLONG */
        pf(&format!("{sizell} bytes)\n"));

        fp_info_hdu(&mut infptr)?;

        fits_close_file(infptr.take().unwrap(), &mut stat);
        if stat != 0 {
            report_error(stat);
            return Err(FpExit(stat));
        }

        iarg += 1;
    }
    Ok(0)
}
/*--------------------------------------------------------------------------*/
/// Takes the owning `Option` because `fp_abort_output` closes the file.
pub(crate) fn fp_info_hdu(infptr: &mut Option<Box<fitsfile>>) -> FpResult<c_int> {
    /* NOTE (upstream bug 5): declared once but reread every pass, so a
    low-dimension HDU inherits the previous one's axes.  Reproduced. */
    let mut naxes: [c_long; 9] = [1, 1, 1, 1, 1, 1, 1, 1, 1];
    let mut val: [c_char; SZ_CARD] = [0; SZ_CARD];
    /* DEVIATION: FLEN_COMMENT rather than the C's SZ_CARD; never read. */
    let mut com: [c_char; FLEN_COMMENT] = [0; FLEN_COMMENT];
    let mut naxis: c_int = 0;
    let mut hdutype: c_int = 0;
    let mut bitpix: c_int = 0;
    let mut hdupos: c_int;
    let mut stat: c_int = 0;
    let mut ii: c_int;
    let mut datasum: c_ulong = 0;
    let mut hdusum: c_ulong = 0;

    fits_movabs_hdu(infptr.as_mut().unwrap(), 1, None, &mut stat);
    if stat != 0 {
        return Err(fp_abort_output(infptr.take(), None, stat));
    }

    hdupos = 1;
    while stat == 0 {
        fits_get_hdu_type(infptr.as_mut().unwrap(), &mut hdutype, &mut stat);
        if stat != 0 {
            return Err(fp_abort_output(infptr.take(), None, stat));
        }

        /* fits_get_hdu_type calls unknown extensions "IMAGE_HDU"
         * so consult XTENSION keyword itself
         */
        fits_read_keyword(
            infptr.as_mut().unwrap(),
            cs!(c"XTENSION"),
            &mut val,
            Some(&mut com),
            &mut stat,
        );
        if stat == KEY_NO_EXIST {
            /* in primary HDU which by definition is an "image" */
            stat = 0; /* clear for later error handling */
        } else if stat != 0 {
            return Err(fp_abort_output(infptr.take(), None, stat));
        } else if hdutype == IMAGE_HDU {
            /* that is, if XTENSION != "IMAGE" AND != "BINTABLE" */
            /* NOTE (upstream bug 9): only 5 characters of "BINTABLE" */
            if strncmp_safe(&val[1..], cs!(c"IMAGE"), 5) != 0
                && strncmp_safe(&val[1..], cs!(c"BINTABLE"), 5) != 0
            {
                /* assign something other than any of these */
                hdutype = IMAGE_HDU + ASCII_TBL + BINARY_TBL;
            }
        }

        fits_get_chksum(
            infptr.as_mut().unwrap(),
            &mut datasum,
            &mut hdusum,
            &mut stat,
        );

        if hdutype == IMAGE_HDU {
            pf(&format!("  {hdupos} IMAGE"));
            pf(&format!(" SUMS={}/{}", not_int(hdusum), datasum));

            fits_get_img_param(
                infptr.as_mut().unwrap(),
                9,
                Some(&mut bitpix),
                Some(&mut naxis),
                Some(&mut naxes),
                &mut stat,
            );

            pf(&format!(" BITPIX={bitpix}"));

            if naxis == 0 {
                pf(" [no_pixels]");
            } else if naxis == 1 {
                /* NOTE (upstream bug 4): should be naxes[0] */
                pf(&format!(" [{}]", naxes[1]));
            } else {
                pf(&format!(" [{}", naxes[0]));
                ii = 1;
                while ii < naxis {
                    pf(&format!("x{}", naxes[ii as usize]));
                    ii += 1;
                }
                fp_msg_str("]");
            }

            if fits_is_compressed_image(infptr.as_mut().unwrap(), &mut stat) != 0 {
                fits_read_keyword(
                    infptr.as_mut().unwrap(),
                    cs!(c"ZCMPTYPE"),
                    &mut val,
                    Some(&mut com),
                    &mut stat,
                );

                /* allow for quote in keyword value */
                if strncmp_safe(&val[1..], cs!(c"RICE_1"), 6) == 0 {
                    fp_msg_str(" tiled_rice\n");
                } else if strncmp_safe(&val[1..], cs!(c"GZIP_1"), 6) == 0 {
                    fp_msg_str(" tiled_gzip_1\n");
                } else if strncmp_safe(&val[1..], cs!(c"GZIP_2"), 6) == 0 {
                    fp_msg_str(" tiled_gzip_2\n");
                } else if strncmp_safe(&val[1..], cs!(c"PLIO_1"), 6) == 0 {
                    fp_msg_str(" tiled_plio\n");
                } else if strncmp_safe(&val[1..], cs!(c"HCOMPRESS_1"), 11) == 0 {
                    fp_msg_str(" tiled_hcompress\n");
                } else {
                    fp_msg_str(" unknown\n");
                }
            } else {
                fp_msg_str(" not_tiled\n");
            }
        } else if hdutype == ASCII_TBL {
            pf(&format!("  {hdupos} ASCII_TBL"));
            pf(&format!(" SUMS={}/{}\n", not_int(hdusum), datasum));
        } else if hdutype == BINARY_TBL {
            pf(&format!("  {hdupos} BINARY_TBL"));
            pf(&format!(" SUMS={}/{}\n", not_int(hdusum), datasum));
        } else {
            pf(&format!("  {hdupos} OTHER"));
            pf(&format!(" SUMS={}/{}", not_int(hdusum), datasum));
            pf(&format!(" {}\n", String::from_utf8_lossy(cbytes(&val))));
        }

        fits_movrel_hdu(infptr.as_mut().unwrap(), 1, None, &mut stat);
        hdupos += 1;
    }
    Ok(0)
}

/*--------------------------------------------------------------------------*/
pub(crate) fn fp_preflight(
    argc: c_int,
    argv: &Argv,
    unpack: c_int,
    fpptr: &mut fpstate,
) -> FpResult<c_int> {
    let mut infits: [c_char; SZ_STR] = [0; SZ_STR];
    let mut outfits: [c_char; SZ_STR] = [0; SZ_STR];
    let mut iarg: c_int;
    let mut namelen: usize;
    let mut nfiles: c_int = 0;

    if fpptr.initialized != FP_INIT_MAGIC {
        fp_msg_str("Error: internal initialization error\n");
        return Err(FpExit(-1));
    }

    iarg = fpptr.firstfile;
    'files: while iarg < argc {
        let arg = &argv[iarg as usize];

        outfits[0] = b'\0' as c_char;

        if strlen_safe(arg) > SZ_STR - 4 {
            /* allow for .fz or .gz suffix */
            fp_msg_str("Error: input file name\n   ");
            fp_msg(arg);
            fp_msg_str("\n   is too long\n");
            fp_noop();
            return Err(FpExit(-1));
        }

        strncpy_trunc(&mut infits, arg);
        if infits[0] == bb(b'-') && infits[1] != 0 {
            /* don't interpret this as intending to read input file from stdin */
            fp_msg_str("Error: invalid input file name\n   ");
            fp_msg(arg);
            fp_msg_str("\n");
            fp_noop();
            return Err(FpExit(-1));
        }

        if strchr_safe(&infits, bb(b'[')).is_some() || strchr_safe(&infits, bb(b']')).is_some() {
            fp_msg_str("Error: section/extension notation not supported: ");
            fp_msg(&infits);
            fp_msg_str("\n");
            fp_noop();
            return Err(FpExit(-1));
        }

        if unpack != 0 {
            /* ********** This section applies to funpack ************ */

            /* check that input file  exists */
            if infits[0] != bb(b'-') {
                /* if not reading from stdin stream */
                if fp_access(&infits) != 0 {
                    /* if not, then check if */
                    strcat_safe(&mut infits, cs!(c".fz")); /* a .fz version exsits */
                    if fp_access(&infits) != 0 {
                        namelen = strlen_safe(&infits);
                        infits[namelen - 3] = 0; /* remove the .fz suffix */
                        fp_msg_str("Error: can't find or read input file ");
                        fp_msg(&infits);
                        fp_msg_str("\n");
                        fp_noop();
                        return Err(FpExit(-1));
                    }
                } else {
                    /* make sure a .fz version of the same file doesn't exist */
                    namelen = strlen_safe(&infits);
                    strcat_safe(&mut infits, cs!(c".fz"));
                    if fp_access(&infits) == 0 {
                        infits[namelen] = 0; /* remove the .fz suffix */
                        fp_msg_str(
                            "Error: ambiguous input file name.  Which file should be unpacked?:\n  ",
                        );
                        fp_msg(&infits);
                        fp_msg_str("\n  ");
                        fp_msg(&infits);
                        fp_msg_str(".fz\n");
                        fp_noop();
                        return Err(FpExit(-1));
                    } else {
                        infits[namelen] = 0; /* remove the .fz suffix */
                    }
                }
            }

            /* if writing to stdout, then we are all done */
            if fpptr.to_stdout != 0 {
                iarg += 1;
                continue 'files;
            }

            if fpptr.outfile[0] != 0 {
                /* user specified output file name */
                nfiles += 1;
                if nfiles > 1 {
                    fp_msg_str("Error: cannot use same output file name for multiple files:\n   ");
                    fp_msg(&fpptr.outfile);
                    fp_msg_str("\n");
                    fp_noop();
                    return Err(FpExit(-1));
                }

                /* check that output file doesn't exist */
                if fp_access(&fpptr.outfile) == 0 {
                    fp_msg_str("Error: output file already exists:\n ");
                    fp_msg(&fpptr.outfile);
                    fp_msg_str("\n ");
                    fp_noop();
                    return Err(FpExit(-1));
                }
                iarg += 1;
                continue 'files;
            }

            /* construct output file name to test */
            if fpptr.prefix[0] != 0 {
                if strlen_safe(&fpptr.prefix) + strlen_safe(&infits) > SZ_STR - 1 {
                    fp_msg_str("Error: output file name for\n   ");
                    fp_msg(&infits);
                    fp_msg_str("\n   is too long with the prefix\n");
                    fp_noop();
                    return Err(FpExit(-1));
                }
                strcpy_safe(&mut outfits, &fpptr.prefix);
            }

            /* construct output file name */
            /* NOTE (upstream bug 12): discards the prefix copied in above */
            if infits[0] == bb(b'-') {
                strcpy_safe(&mut outfits, cs!(c"output.fits"));
            } else {
                strcat_safe(&mut outfits, &infits);
            }

            /* remove .gz or .bz2 suffix, if present (output is not gzipped) */
            namelen = strlen_safe(&outfits);
            if namelen >= 3 && strcmp_safe(cs!(c".gz"), &outfits[namelen - 3..]) == 0 {
                outfits[namelen - 3] = 0;
            } else if namelen >= 4 && strcmp_safe(cs!(c".bz2"), &outfits[namelen - 4..]) == 0 {
                outfits[namelen - 4] = 0;
            }

            /* check for .fz suffix that is sometimes required */
            /* and remove it if present */
            if infits[0] != bb(b'-') {
                /* if not reading from stdin stream */
                namelen = strlen_safe(&outfits);
                if namelen >= 3 && strcmp_safe(cs!(c".fz"), &outfits[namelen - 3..]) == 0 {
                    /* suffix is present */
                    outfits[namelen - 3] = 0;
                } else if fpptr.delete_suffix != 0 {
                    /* required suffix is missing */
                    fp_msg_str("Error: input compressed file ");
                    fp_msg(&infits);
                    fp_msg_str("\n does not have the default .fz suffix.\n");
                    fp_noop();
                    return Err(FpExit(-1));
                }
            }

            /* if infits != outfits, make sure outfits doesn't already exist */
            if strcmp_safe(&infits, &outfits) != 0 && fp_access(&outfits) == 0 {
                fp_msg_str("Error: output file already exists:\n ");
                fp_msg(&outfits);
                fp_msg_str("\n ");
                fp_noop();
                return Err(FpExit(-1));
            }

            /* if gzipping the output, make sure .gz file doesn't exist */
            if fpptr.do_gzip_file != 0 {
                if strlen_safe(&outfits) + 3 > SZ_STR - 1 {
                    fp_msg_str("Error: output file name too long:\n ");
                    fp_msg(&outfits);
                    fp_msg_str("\n ");
                    fp_noop();
                    return Err(FpExit(-1));
                }
                strcat_safe(&mut outfits, cs!(c".gz"));
                if fp_access(&outfits) == 0 {
                    fp_msg_str("Error: output file already exists:\n ");
                    fp_msg(&outfits);
                    fp_msg_str("\n ");
                    fp_noop();
                    return Err(FpExit(-1));
                }
                namelen = strlen_safe(&outfits);
                outfits[namelen - 3] = 0; /* remove the .gz suffix again */
            }
        } else {
            /* ********** This section applies to fpack ************ */

            /* check that input file  exists */
            if infits[0] != bb(b'-') {
                /* if not reading from stdin stream */
                if fp_access(&infits) != 0 {
                    /* if not, then check if */
                    if strlen_safe(&infits) + 3 > SZ_STR - 1 {
                        fp_msg_str("Error: input file name too long:\n ");
                        fp_msg(&infits);
                        fp_msg_str("\n ");
                        fp_noop();
                        return Err(FpExit(-1));
                    }
                    strcat_safe(&mut infits, cs!(c".gz")); /* a gzipped version exsits */
                    if fp_access(&infits) != 0 {
                        namelen = strlen_safe(&infits);
                        infits[namelen - 3] = 0; /* remove the .gz suffix */
                        fp_msg_str("Error: can't find or read input file ");
                        fp_msg(&infits);
                        fp_msg_str("\n");
                        fp_noop();
                        return Err(FpExit(-1));
                    }
                }
            }

            /* make sure the file to pack does not already have a .fz suffix */
            namelen = strlen_safe(&infits);
            if namelen >= 3 && strcmp_safe(cs!(c".fz"), &infits[namelen - 3..]) == 0 {
                fp_msg_str("Error: fpack input file already has '.fz' suffix\n");
                fp_msg(&infits);
                fp_msg_str("\n");
                fp_noop();
                return Err(FpExit(-1));
            }

            /* if writing to stdout, or just testing the files, then we are all done */
            if fpptr.to_stdout != 0 || fpptr.test_all != 0 {
                iarg += 1;
                continue 'files;
            }

            if fpptr.outfile[0] != 0 {
                /* user specified output file name */
                nfiles += 1;
                if nfiles > 1 {
                    fp_msg_str("Error: cannot use same output file name for multiple files:\n   ");
                    fp_msg(&fpptr.outfile);
                    fp_msg_str("\n");
                    fp_noop();
                    return Err(FpExit(-1));
                }

                /* check that output file doesn't exist */
                if fp_access(&fpptr.outfile) == 0 {
                    fp_msg_str("Error: output file already exists:\n ");
                    fp_msg(&fpptr.outfile);
                    fp_msg_str("\n ");
                    fp_noop();
                    return Err(FpExit(-1));
                }
                iarg += 1;
                continue 'files;
            }

            /* construct output file name */
            if infits[0] == bb(b'-') {
                strcpy_safe(&mut outfits, cs!(c"input.fits"));
            } else {
                strcpy_safe(&mut outfits, &infits);
            }

            /* remove .gz suffix, if present (output is not gzipped) */
            /* do the same if compression suffix is bz2 */
            namelen = strlen_safe(&outfits);
            if namelen >= 3 && strcmp_safe(cs!(c".gz"), &outfits[namelen - 3..]) == 0 {
                outfits[namelen - 3] = 0;
            } else if namelen >= 4 && strcmp_safe(cs!(c".bz2"), &outfits[namelen - 4..]) == 0 {
                outfits[namelen - 4] = 0;
            }

            /* remove .imh suffix (IRAF format image), and replace with .fits */
            namelen = strlen_safe(&outfits);
            if namelen >= 4 && strcmp_safe(cs!(c".imh"), &outfits[namelen - 4..]) == 0 {
                outfits[namelen - 4] = 0;
                if strlen_safe(&outfits) == SZ_STR - 5 {
                    strcat_safe(&mut outfits, cs!(c".fit"));
                } else {
                    strcat_safe(&mut outfits, cs!(c".fits"));
                }
            }

            /* If not clobbering the input file, add .fz suffix to output name */
            if fpptr.clobber == 0 {
                if strlen_safe(&outfits) > SZ_STR - 4 {
                    fp_msg_str("Error: output file name too long:\n ");
                    fp_msg(&outfits);
                    fp_msg_str("\n ");
                    fp_noop();
                    return Err(FpExit(-1));
                } else {
                    strcat_safe(&mut outfits, cs!(c".fz"));
                }
            }

            /* if infits != outfits, make sure outfits doesn't already exist */
            if strcmp_safe(&infits, &outfits) != 0 && fp_access(&outfits) == 0 {
                fp_msg_str("Error: output file already exists:\n ");
                fp_msg(&outfits);
                fp_msg_str("\n ");
                fp_noop();
                return Err(FpExit(-1));
            }
        } /* end of fpack section */

        iarg += 1;
    }

    fpptr.preflight_checked += 1;
    Ok(0)
}

/*--------------------------------------------------------------------------*/
/* must run fp_preflight() before fp_loop()
 */
pub(crate) fn fp_loop(argc: c_int, argv: &Argv, unpack: c_int, fpvar: fpstate) -> FpResult<c_int> {
    let mut infits: [c_char; SZ_STR] = [0; SZ_STR];
    let mut outfits: [c_char; SZ_STR] = [0; SZ_STR];
    let mut temp: [c_char; SZ_STR] = [0; SZ_STR];
    let mut iarg: c_int;
    let mut islossless: c_int;
    let mut namelen: usize;
    let mut iraf_infile: c_int = 0;
    let mut status: c_int = 0;

    if fpvar.initialized != FP_INIT_MAGIC {
        fp_msg_str("Error: internal initialization error\n");
        return Err(FpExit(-1));
    } else if fpvar.preflight_checked == 0 {
        fp_msg_str("Error: internal preflight error\n");
        return Err(FpExit(-1));
    }

    if fpvar.test_all != 0 && fpvar.outfile[0] != 0 {
        match File::create(cpath(&fpvar.outfile)) {
            Ok(f) => OUTREPORT.with_borrow_mut(|r| *r = Some(BufWriter::new(f))),
            Err(_) => {
                /* the C never checks fopen; a NULL stream makes every
                fprintf a no-op, which None reproduces */
            }
        }
        report(
            " Filename Extension BITPIX NAXIS1 NAXIS2 Size N_nulls Minval Maxval Mean Sigm Noise1 Noise2 Noise3 Noise5 T_whole T_rowbyrow ",
        );
        report(
            "[Comp_ratio, Pack_cpu, Unpack_cpu, Lossless readtimes] (repeated for Rice, Hcompress, and GZIP)\n",
        );
    }

    /* DEVIATION: the C installs abort_fpack here to remove its three temp
    files; each is a TempPath below instead. */

    iarg = fpvar.firstfile;
    'files: while iarg < argc {
        /* dropping these unlinks them, as abort_fpack did */
        let tempfile1: TempPath;
        let tempfile2: TempPath;
        #[allow(unused_assignments)]
        let mut clobber_tmp: Option<TempPath> = None;

        temp[0] = 0;
        outfits[0] = 0;
        islossless = 1;

        strncpy_trunc(&mut infits, &argv[iarg as usize]);

        if unpack != 0 {
            /* ********** This section applies to funpack ************ */

            /* find input file */
            if infits[0] != bb(b'-') {
                /* if not reading from stdin stream */
                if fp_access(&infits) != 0 {
                    /* if not, then */
                    strcat_safe(&mut infits, cs!(c".fz")); /* a .fz version must exsit */
                    /* fp_preflight already checked for enough size to add '.fz' */
                }
            }

            if fpvar.to_stdout != 0 {
                strcpy_safe(&mut outfits, cs!(c"-"));
            } else if fpvar.outfile[0] != 0 {
                /* user specified output file name */
                strcpy_safe(&mut outfits, &fpvar.outfile);
            } else {
                /* construct output file name */
                if fpvar.prefix[0] != 0 {
                    /* fp_preflight already checked this */
                    strcpy_safe(&mut outfits, &fpvar.prefix);
                }

                /* construct output file name */
                if infits[0] == bb(b'-') {
                    strcpy_safe(&mut outfits, cs!(c"output.fits"));
                } else {
                    strcat_safe(&mut outfits, &infits);
                }

                /* remove .gz suffix, if present (output is not gzipped) */
                namelen = strlen_safe(&outfits);
                if namelen >= 3 && strcmp_safe(cs!(c".gz"), &outfits[namelen - 3..]) == 0 {
                    outfits[namelen - 3] = 0;
                } else if namelen >= 4 && strcmp_safe(cs!(c".bz2"), &outfits[namelen - 4..]) == 0 {
                    outfits[namelen - 4] = 0;
                }

                /* check for .fz suffix that is sometimes required */
                /* and remove it if present */
                namelen = strlen_safe(&outfits);
                if namelen >= 3 && strcmp_safe(cs!(c".fz"), &outfits[namelen - 3..]) == 0 {
                    /* suffix is present */
                    outfits[namelen - 3] = 0;
                }
            }
        } else {
            /* ********** This section applies to fpack ************ */

            if fpvar.to_stdout != 0 {
                strcpy_safe(&mut outfits, cs!(c"-"));
            } else if fpvar.test_all == 0 {
                if fpvar.outfile[0] != 0 {
                    /* user specified output file name */
                    strcpy_safe(&mut outfits, &fpvar.outfile);
                } else {
                    /* construct output file name */
                    if infits[0] == bb(b'-') {
                        strcpy_safe(&mut outfits, cs!(c"input.fits"));
                    } else {
                        strcpy_safe(&mut outfits, &infits);
                    }
                    /* Remove .gz suffix, if present (output is not gzipped).
                    Do the same for .bz2 */
                    namelen = strlen_safe(&outfits);
                    if namelen >= 3 && strcmp_safe(cs!(c".gz"), &outfits[namelen - 3..]) == 0 {
                        outfits[namelen - 3] = 0;
                    } else if namelen >= 4
                        && strcmp_safe(cs!(c".bz2"), &outfits[namelen - 4..]) == 0
                    {
                        outfits[namelen - 4] = 0;
                    }

                    /* remove .imh suffix (IRAF format image), and replace with .fits */
                    namelen = strlen_safe(&outfits);
                    if namelen >= 4 && strcmp_safe(cs!(c".imh"), &outfits[namelen - 4..]) == 0 {
                        outfits[namelen - 4] = 0;
                        if strlen_safe(&outfits) == SZ_STR - 5 {
                            strcat_safe(&mut outfits, cs!(c".fit"));
                        } else {
                            strcat_safe(&mut outfits, cs!(c".fits"));
                        }
                        iraf_infile = 1; /* this is an IRAF format input file */
                        /* change the output name to "NAME.fits.fz" */
                    }

                    /* If not clobbering the input file, add .fz suffix to output name */
                    if fpvar.clobber == 0 {
                        strcat_safe(&mut outfits, cs!(c".fz"));
                    }
                }
            }
        }

        strncpy_trunc(&mut temp, &outfits);

        if infits[0] != bb(b'-') {
            /* if not reading from stdin stream */
            if strcmp_safe(&infits, &outfits) == 0 {
                /* are input and output names the same? */

                /* clobber the input file with the output file with the same name */
                if fpvar.clobber == 0 {
                    fp_msg_str("\nError: must use -F flag to clobber input file.\n");
                    return Err(FpExit(-1));
                }

                /* create temporary file name in the output directory (same as input directory)*/
                clobber_tmp = Some(fp_tmpnam_guard(cs!(c"Tmp1"), &infits, &mut outfits)?);
            }
        }

        /* *************** now do the real work ********************* */

        if fpvar.verbose != 0 && fpvar.to_stdout == 0 {
            pf(&format!("{} ", String::from_utf8_lossy(cbytes(&infits))));
        }

        if fpvar.test_all != 0 {
            /* compare all the algorithms */

            /* create 2 temporary file names, in the CWD */
            let mut tf1: [c_char; SZ_STR] = [0; SZ_STR];
            let mut tf2: [c_char; SZ_STR] = [0; SZ_STR];
            tempfile1 = fp_tmpnam_guard(cs!(c"Tmpfile1"), cs!(c""), &mut tf1)?;
            tempfile2 = fp_tmpnam_guard(cs!(c"Tmpfile2"), cs!(c""), &mut tf2)?;

            fp_test(&infits, &tf1, &tf2, fpvar)?;

            /* the C removes both here; TempPath does it on drop */
            drop(tempfile1);
            drop(tempfile2);
            iarg += 1;
            continue 'files;
        } else if unpack != 0 {
            if fpvar.to_stdout != 0 {
                /* unpack the input file to the stdout stream */
                fp_unpack(&infits, &outfits, fpvar)?;
            } else {
                /* unpack to temporary file, so other tasks can't open it until it is renamed */

                /* create  temporary file name, in the output directory */
                let mut tf2: [c_char; SZ_STR] = [0; SZ_STR];
                tempfile2 = fp_tmpnam_guard(cs!(c"Tmp2"), &outfits, &mut tf2)?;

                /* unpack the input file to the temporary file */
                fp_unpack(&infits, &tf2, fpvar)?;

                /* rename the temporary file to it's real name */
                match std::fs::rename(cpath(&tf2), cpath(&outfits)) {
                    Err(_) => {
                        fp_msg_str("Failed to rename temporary file name:\n  ");
                        fp_msg(&tf2);
                        fp_msg_str(" -> ");
                        fp_msg(&outfits);
                        fp_msg_str("\n");
                        return Err(FpExit(-1));
                    }
                    Ok(()) => {
                        /* the file has moved, so disarm the guard */
                        let _ = tempfile2.keep();
                    }
                }
            }
        } else {
            fp_pack(&infits, &outfits, fpvar, &mut islossless)?;
        }

        if fpvar.to_stdout != 0 {
            iarg += 1;
            continue 'files;
        }

        /* ********** clobber and/or delete files, if needed ************** */

        if strcmp_safe(&infits, &temp) == 0 && fpvar.clobber != 0 {
            if islossless == 0 && fpvar.do_not_prompt == 0 {
                fp_msg_str("\nFile ");
                fp_msg(&infits);
                fp_msg_str("\nwas compressed with a LOSSY method.  Overwrite the\n");
                fp_msg_str("original file with the compressed version? (Y/N) ");
                let answer = read_answer();
                if answer != b'Y' && answer != b'y' {
                    fp_msg_str("\noriginal file NOT overwritten!\n");
                    let _ = std::fs::remove_file(cpath(&outfits));
                    iarg += 1;
                    continue 'files;
                }
            }

            if iraf_infile != 0 {
                /* special case of deleting an IRAF format header and pixel file */
                if fits_delete_iraf_file(&infits, &mut status) != 0 {
                    fp_msg_str("\nError deleting IRAF .imh and .pix files.\n");
                    fp_msg(&infits);
                    fp_msg_str("\n");
                    return Err(FpExit(-1));
                }
            }

            if cfg!(unix) {
                /* rename clobbers input on Unix platforms */
                if std::fs::rename(cpath(&outfits), cpath(&temp)).is_err() {
                    fp_msg_str("\nError renaming tmp file to ");
                    fp_msg(&temp);
                    fp_msg_str("\n");
                    return Err(FpExit(-1));
                }
            } else {
                /* rename DOES NOT clobber existing files on Windows platforms */
                /* so explicitly remove any existing file before renaming the file */
                let _ = std::fs::remove_file(cpath(&temp));
                if std::fs::rename(cpath(&outfits), cpath(&temp)).is_err() {
                    fp_msg_str("\nError renaming tmp file to ");
                    fp_msg(&temp);
                    fp_msg_str("\n");
                    return Err(FpExit(-1));
                }
            }

            /* the temp file has moved, so disarm its guard */
            if let Some(g) = clobber_tmp.take() {
                let _ = g.keep();
            }
            strcpy_safe(&mut outfits, &temp);
        } else if fpvar.clobber != 0 || fpvar.delete_input != 0 {
            /* delete the input file */
            if islossless == 0 && fpvar.do_not_prompt == 0 {
                /* user did not turn off delete prompt */
                fp_msg_str("\nFile ");
                fp_msg(&infits);
                fp_msg_str("\nwas compressed with a LOSSY method.  \n");
                fp_msg_str("Delete the original file? (Y/N) ");
                let answer = read_answer();
                if answer != b'Y' && answer != b'y' {
                    /* user abort */
                    fp_msg_str("\noriginal file NOT deleted!\n");
                } else if iraf_infile != 0 {
                    /* special case of deleting an IRAF format header and pixel file */
                    if fits_delete_iraf_file(&infits, &mut status) != 0 {
                        fp_msg_str("\nError deleting IRAF .imh and .pix files.\n");
                        fp_msg(&infits);
                        fp_msg_str("\n");
                        return Err(FpExit(-1));
                    }
                } else if std::fs::remove_file(cpath(&infits)).is_err() {
                    /* normal case of deleting input FITS file */
                    fp_msg_str("\nError deleting input file ");
                    fp_msg(&infits);
                    fp_msg_str("\n");
                    return Err(FpExit(-1));
                }
            } else {
                /* user said don't prompt, so just delete the input file */
                if iraf_infile != 0 {
                    /* special case of deleting an IRAF format header and pixel file */
                    if fits_delete_iraf_file(&infits, &mut status) != 0 {
                        fp_msg_str("\nError deleting IRAF .imh and .pix files.\n");
                        fp_msg(&infits);
                        fp_msg_str("\n");
                        return Err(FpExit(-1));
                    }
                } else if std::fs::remove_file(cpath(&infits)).is_err() {
                    /* normal case of deleting input FITS file */
                    fp_msg_str("\nError deleting input file ");
                    fp_msg(&infits);
                    fp_msg_str("\n");
                    return Err(FpExit(-1));
                }
            }
        }
        iraf_infile = 0;

        if fpvar.do_gzip_file != 0 {
            /* gzip the output file.  DEVIATION: the C hands "gzip -1 <file>"
            to system() behind a filename whitelist; there is no shell here, so
            the whitelist goes with it. */
            if let Err(e) = gzip_file(&outfits) {
                fp_msg_str("\nError gzipping output file ");
                fp_msg(&outfits);
                pf(&format!("\n  {e}\n"));
                return Err(FpExit(-1));
            }
            strcat_safe(&mut outfits, cs!(c".gz")); /* only possibible with funpack */
        }

        if fpvar.verbose != 0 && fpvar.to_stdout == 0 {
            pf(&format!(
                "-> {}\n",
                String::from_utf8_lossy(cbytes(&outfits))
            ));
        }

        iarg += 1;
    }

    if fpvar.test_all != 0 && fpvar.outfile[0] != 0 {
        OUTREPORT.with_borrow_mut(|r| *r = None); /* fclose(outreport) */
    }
    Ok(0)
}

/// C: `fgets(answer, 29, stdin); if (answer[0] != 'Y' && ...)`.
fn read_answer() -> u8 {
    use std::io::BufRead;
    let mut line = String::new();
    let _ = std::io::stdin().lock().read_line(&mut line);
    line.as_bytes().first().copied().unwrap_or(0)
}

/// Replacement for `system("gzip -1 <file>")`, streamed rather than buffered.
fn gzip_file(file: &[c_char]) -> std::io::Result<()> {
    use flate2::Compression;
    use flate2::write::GzEncoder;

    let src = cpath(file);
    /* from the path, not by strcat into the caller's buffer */
    let mut dst = src.clone().into_os_string();
    dst.push(".gz");
    let dst = PathBuf::from(dst);

    let mut input = File::open(&src)?;
    let mut encoder = GzEncoder::new(BufWriter::new(File::create(&dst)?), Compression::new(1));
    std::io::copy(&mut input, &mut encoder)?;
    encoder.finish()?.flush()?;
    drop(input);

    std::fs::remove_file(&src)
}

/*--------------------------------------------------------------------------*/
/* fp_pack assumes the output file does not exist (checked by preflight)
 */
pub(crate) fn fp_pack(
    infits: &[c_char],
    outfits: &[c_char],
    fpvar: fpstate,
    islossless: &mut c_int,
) -> FpResult<c_int> {
    let mut infptr: Option<Box<fitsfile>> = None;
    let mut outfptr: Option<Box<fitsfile>> = None;
    let mut stat: c_int = 0;

    fits_open_file(&mut infptr, infits, READONLY, &mut stat);
    if stat != 0 {
        report_error(stat);
        return Err(FpExit(stat));
    }

    fits_create_file(&mut outfptr, outfits, &mut stat);
    if stat != 0 {
        return Err(fp_abort_output(infptr, None, stat));
    }

    while stat == 0 {
        /*  LOOP OVER EACH HDU */

        let inf = infptr.as_mut().unwrap();
        let outf = outfptr.as_mut().unwrap();

        fits_set_lossy_int(outf, fpvar.int_to_float, &mut stat);
        fits_set_compression_type(outf, fpvar.comptype, &mut stat);
        fits_set_tile_dim(outf, 6, &fpvar.ntile, &mut stat);

        if fpvar.no_dither != 0 {
            fits_set_quantize_method(outf, -1, &mut stat);
        } else {
            fits_set_quantize_method(outf, fpvar.dither_method, &mut stat);
        }

        fits_set_quantize_level(outf, fpvar.quantize_level, &mut stat);
        fits_set_dither_offset(outf, fpvar.dither_offset, &mut stat);
        fits_set_hcomp_scale(outf, fpvar.scale, &mut stat);
        fits_set_hcomp_smooth(outf, fpvar.smooth, &mut stat);

        fp_pack_hdu(inf, outf, fpvar, islossless, &mut stat);

        if fpvar.do_checksums != 0 {
            fits_write_chksum(outf, &mut stat);
        }

        fits_movrel_hdu(inf, 1, None, &mut stat);
    }

    if stat == END_OF_FILE {
        stat = 0;
    }

    /* set checksum for case of newly created primary HDU	 */

    if fpvar.do_checksums != 0 {
        let outf = outfptr.as_mut().unwrap();
        fits_movabs_hdu(outf, 1, None, &mut stat);
        fits_write_chksum(outf, &mut stat);
    }

    if stat != 0 {
        return Err(fp_abort_output(infptr, outfptr, stat));
    }

    fits_close_file(outfptr.take().unwrap(), &mut stat);
    fits_close_file(infptr.take().unwrap(), &mut stat);

    Ok(0)
}

/*--------------------------------------------------------------------------*/
/* fp_unpack assumes the output file does not exist
 */
pub(crate) fn fp_unpack(infits: &[c_char], outfits: &[c_char], fpvar: fpstate) -> FpResult<c_int> {
    let mut infptr: Option<Box<fitsfile>> = None;
    let mut outfptr: Option<Box<fitsfile>> = None;
    let mut stat: c_int = 0;
    let mut hdutype: c_int = 0;
    let mut extnum: c_int;
    let mut single = 0;
    let mut hduname: [c_char; SZ_STR] = [0; SZ_STR];

    /* the C walks fpvar.extname with char pointers, writing NULs over the
    commas; here that walk is an index into a local copy */
    let mut extname: [c_char; SZ_STR] = fpvar.extname;
    let mut hduloc: usize = 0;

    fits_open_file(&mut infptr, infits, READONLY, &mut stat);
    fits_create_file(&mut outfptr, outfits, &mut stat);

    if stat != 0 {
        return Err(fp_abort_output(infptr, outfptr, stat));
    }

    if fpvar.extname[0] != 0 {
        /* unpack a list of HDUs? */

        /* move to the first HDU in the list */
        let loc = strchr_safe(&extname[hduloc..], bb(b',')).map(|i| hduloc + i);

        if let Some(loc) = loc {
            extname[loc] = 0; /* terminate the first name in the string */
        }

        strcpy_safe(&mut hduname, &extname[hduloc..]); /* copy the first name into temporary string */

        if let Some(loc) = loc {
            hduloc = loc + 1; /* advance to the beginning of the next name, if any */
        } else {
            hduloc += strlen_safe(&hduname); /* end of the list */
            single = 1; /* only 1 HDU is being unpacked */
        }

        if (hduname[0] as u8).is_ascii_digit() {
            /* read the string as an integer */
            let (val, end) = strtol_c(&hduname);
            extnum = val as c_int;

            /* check for junk following the integer */
            if hduname[end] == 0 {
                /* no junk, so move to this HDU number (+1) */
                fits_movabs_hdu(
                    infptr.as_mut().unwrap(),
                    extnum + 1,
                    Some(&mut hdutype),
                    &mut stat,
                ); /* move to HDU number */
                if hdutype != IMAGE_HDU {
                    stat = NOT_IMAGE;
                }
            } else {
                /* the string is not an integer, so must be the column name */
                hdutype = IMAGE_HDU;
                fits_movnam_hdu(infptr.as_mut().unwrap(), hdutype, &hduname, 0, &mut stat);
            }
        } else {
            /* move to the named image extension */
            hdutype = IMAGE_HDU;
            fits_movnam_hdu(infptr.as_mut().unwrap(), hdutype, &hduname, 0, &mut stat);
        }
    }

    if stat != 0 {
        fp_msg_str("Unable to find and move to extension '");
        fp_msg(&hduname);
        fp_msg_str("'\n");
        return Err(fp_abort_output(infptr, outfptr, stat));
    }

    while stat == 0 {
        if single != 0 {
            stat = -1; /* special status flag to force output primary array */
        }

        fp_unpack_hdu(
            infptr.as_mut().unwrap(),
            outfptr.as_mut().unwrap(),
            fpvar,
            &mut stat,
        );

        if fpvar.do_checksums != 0 {
            fits_write_chksum(outfptr.as_mut().unwrap(), &mut stat);
        }

        /* move to the next HDU */
        if fpvar.extname[0] != 0 {
            /* unpack a list of HDUs? */

            if extname[hduloc] == 0 {
                stat = END_OF_FILE; /* we reached the end of the list */
            } else {
                /* parse the next HDU name and move to it */
                let loc = strchr_safe(&extname[hduloc..], bb(b',')).map(|i| hduloc + i);

                if let Some(loc) = loc {
                    /* look for 'comma' delimiter between names */
                    extname[loc] = 0; /* terminate the first name in the string */
                }

                strcpy_safe(&mut hduname, &extname[hduloc..]); /* copy the next name into temporary string */

                if let Some(loc) = loc {
                    hduloc = loc + 1; /* advance to the beginning of the next name, if any */
                } else {
                    extname[hduloc] = 0; /* end of the list */
                }

                if (hduname[0] as u8).is_ascii_digit() {
                    let (val, end) = strtol_c(&hduname); /* read the string as an integer */
                    extnum = val as c_int;

                    /* check for junk following the integer */
                    if hduname[end] == 0 {
                        /* no junk, so move to this HDU number (+1) */
                        fits_movabs_hdu(
                            infptr.as_mut().unwrap(),
                            extnum + 1,
                            Some(&mut hdutype),
                            &mut stat,
                        ); /* move to HDU number */
                        if hdutype != IMAGE_HDU {
                            stat = NOT_IMAGE;
                        }
                    } else {
                        /* the string is not an integer, so must be the column name */
                        hdutype = IMAGE_HDU;
                        fits_movnam_hdu(infptr.as_mut().unwrap(), hdutype, &hduname, 0, &mut stat);
                    }
                } else {
                    /* move to the named image extension */
                    hdutype = IMAGE_HDU;
                    fits_movnam_hdu(infptr.as_mut().unwrap(), hdutype, &hduname, 0, &mut stat);
                }

                if stat != 0 {
                    fp_msg_str("Unable to find and move to extension '");
                    fp_msg(&hduname);
                    fp_msg_str("'\n");
                }
            }
        } else {
            /* increment to the next HDU */
            fits_movrel_hdu(infptr.as_mut().unwrap(), 1, None, &mut stat);
        }
    }

    if stat == END_OF_FILE {
        stat = 0;
    }

    /* set checksum for case of newly created primary HDU
     */
    if fpvar.do_checksums != 0 {
        let outf = outfptr.as_mut().unwrap();
        fits_movabs_hdu(outf, 1, None, &mut stat);
        fits_write_chksum(outf, &mut stat);
    }

    if stat != 0 {
        return Err(fp_abort_output(infptr, outfptr, stat));
    }

    fits_close_file(outfptr.take().unwrap(), &mut stat);
    fits_close_file(infptr.take().unwrap(), &mut stat);

    Ok(0)
}

/// C: `atol(s)` -- the leading integer, or 0 if there is not one.
pub(crate) fn atol_c(s: &[c_char]) -> c_long {
    strtol_c(s).0
}

/// C: `atoi(s)`.
pub(crate) fn atoi_c(s: &[c_char]) -> c_int {
    atol_c(s) as c_int
}

/// C: `atof(s)` -- the longest leading prefix that parses as a double, or 0.0.
pub(crate) fn atof_c(s: &[c_char]) -> f64 {
    let b = cbytes(s);
    let mut i = 0;
    while i < b.len() && b[i].is_ascii_whitespace() {
        i += 1;
    }
    let start = i;
    if i < b.len() && (b[i] == b'+' || b[i] == b'-') {
        i += 1;
    }
    while i < b.len() && b[i].is_ascii_digit() {
        i += 1;
    }
    if i < b.len() && b[i] == b'.' {
        i += 1;
        while i < b.len() && b[i].is_ascii_digit() {
            i += 1;
        }
    }
    /* only take an exponent if it is complete, as strtod does */
    if i < b.len() && (b[i] == b'e' || b[i] == b'E') {
        let mut j = i + 1;
        if j < b.len() && (b[j] == b'+' || b[j] == b'-') {
            j += 1;
        }
        if j < b.len() && b[j].is_ascii_digit() {
            while j < b.len() && b[j].is_ascii_digit() {
                j += 1;
            }
            i = j;
        }
    }
    std::str::from_utf8(&b[start..i])
        .ok()
        .and_then(|t| t.parse().ok())
        .unwrap_or(0.0)
}

/// C: `strtol(hduname, &loc, 10)`, returning the value and the index `loc`
/// ended at so the caller can check for junk.
fn strtol_c(s: &[c_char]) -> (c_long, usize) {
    let b = cbytes(s);
    let mut i = 0;
    while i < b.len() && (b[i] as char).is_whitespace() {
        i += 1;
    }
    let start = i;
    if i < b.len() && (b[i] == b'+' || b[i] == b'-') {
        i += 1;
    }
    let dstart = i;
    while i < b.len() && b[i].is_ascii_digit() {
        i += 1;
    }
    if i == dstart {
        return (0, 0); /* no conversion: strtol leaves endptr at the start */
    }
    let val: c_long = std::str::from_utf8(&b[start..i])
        .ok()
        .and_then(|t| t.parse().ok())
        .unwrap_or(if b[start] == b'-' {
            c_long::MIN
        } else {
            c_long::MAX
        });
    (val, i)
}

/*--------------------------------------------------------------------------*/
/* fp_test assumes the output files do not exist
 */
pub(crate) fn fp_test(
    infits: &[c_char],
    outfits: &[c_char],
    outfits2: &[c_char],
    fpvar: fpstate,
) -> FpResult<c_int> {
    let mut inputfptr: Option<Box<fitsfile>> = None;
    let mut outfptr: Option<Box<fitsfile>> = None;
    let mut outfptr2: Option<Box<fitsfile>> = None;
    let mut tempfile: Option<Box<fitsfile>> = None;
    /* the rescale temp file's name guard, if one is created */
    let mut tempfile3: Option<TempPath> = None;

    /* NOTE (upstream bug 5): declared once, outside the per-HDU loop, and only
    the first `naxis' entries are refilled.  Reproduced. */
    let mut naxes: [c_long; 9] = [1, 1, 1, 1, 1, 1, 1, 1, 1];
    let mut stat: c_int = 0;
    let mut totpix: c_int = 0;
    let mut naxis: c_int = 0;
    let mut ii: c_int;
    let mut hdutype: c_int = 0;
    let mut bitpix: c_int = 0;
    let mut extnum: c_int = 0;
    let mut len: usize;
    let mut tstatus: c_int;
    let mut hdunum: c_int = 0;
    let mut rescale_flag: bool;
    let mut bpix: c_int = 0;
    let mut ncols: c_int = 0;
    let mut dtype: [c_char; 8] = [0; 8];
    let mut dimen: [c_char; 100] = [0; 100];
    let mut bscale: f64 = 0.0;
    let mut rescale: f64;
    let mut noisemin: f64;
    let mut headstart: c_long = 0;
    let mut datastart: c_long = 0;
    let mut dataend: c_long = 0;
    let mut origdata: f32;
    let mut whole_cpu: f32 = 0.0;
    let mut whole_elapse: f32 = 0.0;
    let mut row_elapse: f32 = 0.0;
    let mut row_cpu: f32 = 0.0;
    let mut xbits: f32;

    let mut nrows: LONGLONG = 0;
    /* structure to hold image statistics (defined in fpack.h)
     *
     * DEVIATION (upstream bug 7): the C leaves this uninitialised and reads it
     * anyway for a BITPIX the statistics branch does not cover (LONGLONG_IMG,
     * and any first HDU that is not an image).  Rust has no uninitialised
     * read; zero is the closest defined equivalent. */
    let mut imagestats = imgstats::default();

    fits_open_file(&mut inputfptr, infits, READONLY, &mut stat);
    fits_create_file(&mut outfptr, outfits, &mut stat);
    fits_create_file(&mut outfptr2, outfits2, &mut stat);

    if stat != 0 {
        report_error(stat);
        return Err(FpExit(stat));
    }

    while stat == 0 {
        /*  LOOP OVER EACH HDU */
        rescale_flag = false;
        fits_get_hdu_type(inputfptr.as_mut().unwrap(), &mut hdutype, &mut stat);

        if hdutype == IMAGE_HDU {
            fits_get_img_param(
                inputfptr.as_mut().unwrap(),
                9,
                Some(&mut bitpix),
                Some(&mut naxis),
                Some(&mut naxes),
                &mut stat,
            );
            /* NOTE (upstream bug 8): totpix is an int and this overflows for
            images above 2^31 pixels.  Reproduced with wrapping arithmetic,
            which is what every compiler CFITSIO supports actually emits. */
            totpix = 1;
            ii = 0;
            while ii < 9 {
                totpix = totpix.wrapping_mul(naxes[ii as usize] as c_int);
                ii += 1;
            }
        }

        if fits_is_compressed_image(inputfptr.as_mut().unwrap(), &mut stat) == 0
            && hdutype == IMAGE_HDU
            && naxis != 0
            && totpix != 0
            && fpvar.do_images != 0
        {
            /* rescale a scaled integer image to reduce noise? */
            if fpvar.rescale_noise != 0. && bitpix > 0 && bitpix < LONGLONG_IMG {
                tstatus = 0;
                fits_read_key(
                    inputfptr.as_mut().unwrap(),
                    KeywordDatatypeMut::TDOUBLE(&mut bscale),
                    cs!(c"BSCALE"),
                    None,
                    &mut tstatus,
                );

                if tstatus == 0 && bscale != 1.0 {
                    /* image must be scaled */

                    if bitpix == LONG_IMG {
                        fp_i4stat(
                            inputfptr.as_mut().unwrap(),
                            naxis,
                            &naxes,
                            &mut imagestats,
                            &mut stat,
                        );
                    } else {
                        fp_i2stat(
                            inputfptr.as_mut().unwrap(),
                            naxis,
                            &naxes,
                            &mut imagestats,
                            &mut stat,
                        );
                    }

                    /* use the minimum of the MAD 2nd, 3rd, and 5th order noise estimates */
                    noisemin = imagestats.noise3;
                    if imagestats.noise2 != 0. && imagestats.noise2 < noisemin {
                        noisemin = imagestats.noise2;
                    }
                    if imagestats.noise5 != 0. && imagestats.noise5 < noisemin {
                        noisemin = imagestats.noise5;
                    }

                    rescale = noisemin / f64::from(fpvar.rescale_noise);
                    if rescale > 1.0 {
                        /* all the criteria are met, so create a temporary file that */
                        /* contains a rescaled version of the image, in CWD */

                        /* create temporary file name */
                        let mut tf3: [c_char; SZ_STR] = [0; SZ_STR];
                        tempfile3 = Some(fp_tmpnam_guard(cs!(c"Tmpfile3"), cs!(c""), &mut tf3)?);

                        fits_create_file(&mut tempfile, &tf3, &mut stat);

                        fits_get_hdu_num(inputfptr.as_mut().unwrap(), &mut hdunum);
                        if hdunum != 1 {
                            /* the input hdu is an image extension, so create dummy primary */
                            fits_create_img(tempfile.as_mut().unwrap(), 8, 0, &naxes, &mut stat);
                        }

                        fits_copy_header(
                            inputfptr.as_mut().unwrap(),
                            tempfile.as_mut().unwrap(),
                            &mut stat,
                        ); /* copy the header */

                        /* rescale the data, so that it will compress more efficiently */
                        if bitpix == LONG_IMG {
                            fp_i4rescale(
                                inputfptr.as_mut().unwrap(),
                                naxis,
                                &naxes,
                                rescale,
                                tempfile.as_mut().unwrap(),
                                &mut stat,
                            );
                        } else {
                            fp_i2rescale(
                                inputfptr.as_mut().unwrap(),
                                naxis,
                                &naxes,
                                rescale,
                                tempfile.as_mut().unwrap(),
                                &mut stat,
                            );
                        }

                        /* scale the BSCALE keyword by the inverse factor */

                        bscale *= rescale;
                        fits_update_key(
                            tempfile.as_mut().unwrap(),
                            KeywordDatatype::TDOUBLE(&bscale),
                            cs!(c"BSCALE"),
                            None,
                            &mut stat,
                        );

                        /* rescan the header, to reset the actual scaling parameters */
                        fits_set_hdustruc(tempfile.as_mut().unwrap(), &mut stat);

                        rescale_flag = true;
                    }
                }
            }

            /* C: infptr aliases either inputfptr or tempfile.  Two owning
            handles cannot share one &mut, so it is re-derived at each use. */

            /* compute basic statistics about the input image */
            {
                let infptr = if rescale_flag {
                    tempfile.as_deref_mut().unwrap()
                } else {
                    inputfptr.as_deref_mut().unwrap()
                };
                if bitpix == BYTE_IMG {
                    bpix = 8;
                    strcpy_safe(&mut dtype, cs!(c"8  "));
                    fp_i2stat(infptr, naxis, &naxes, &mut imagestats, &mut stat);
                } else if bitpix == SHORT_IMG {
                    bpix = 16;
                    strcpy_safe(&mut dtype, cs!(c"16 "));
                    fp_i2stat(infptr, naxis, &naxes, &mut imagestats, &mut stat);
                } else if bitpix == LONG_IMG {
                    bpix = 32;
                    strcpy_safe(&mut dtype, cs!(c"32 "));
                    fp_i4stat(infptr, naxis, &naxes, &mut imagestats, &mut stat);
                } else if bitpix == LONGLONG_IMG {
                    bpix = 64;
                    strcpy_safe(&mut dtype, cs!(c"64 "));
                } else if bitpix == FLOAT_IMG {
                    bpix = 32;
                    strcpy_safe(&mut dtype, cs!(c"-32"));
                    fp_r4stat(infptr, naxis, &naxes, &mut imagestats, &mut stat);
                } else if bitpix == DOUBLE_IMG {
                    bpix = 64;
                    strcpy_safe(&mut dtype, cs!(c"-64"));
                    fp_r4stat(infptr, naxis, &naxes, &mut imagestats, &mut stat);
                }
            }

            /* use the minimum of the MAD 2nd, 3rd, and 5th order noise estimates */
            noisemin = imagestats.noise3;
            if imagestats.noise2 != 0. && imagestats.noise2 < noisemin {
                noisemin = imagestats.noise2;
            }
            if imagestats.noise5 != 0. && imagestats.noise5 < noisemin {
                noisemin = imagestats.noise5;
            }

            xbits = (noisemin.log10() / 0.301 + 1.792) as f32;

            pf(&format!(
                "\n File: {}\n",
                String::from_utf8_lossy(cbytes(infits))
            ));
            pf(
                "  Ext BITPIX Dimens.   Nulls    Min    Max     Mean    Sigma  Noise2  Noise3  Noise5  Nbits   MaxR\n",
            );

            pf(&format!(
                "  {:3}  {}",
                extnum,
                String::from_utf8_lossy(cbytes(&dtype))
            ));
            snprintf_into(&mut dimen, &format!(" ({}", naxes[0]));
            len = strlen_safe(&dimen);
            ii = 1;
            while ii < naxis {
                if len < 99 {
                    let text = format!(",{}", naxes[ii as usize]);
                    snprintf_into(&mut dimen[len..100.min(len + (100 - len))], &text);
                }
                len = strlen_safe(&dimen);
                ii += 1;
            }
            if strlen_safe(&dimen) < 99 {
                strcat_safe(&mut dimen, cs!(c")"));
            }
            pf(&format!("{:<12}", String::from_utf8_lossy(cbytes(&dimen))));

            fits_get_hduaddr(
                inputfptr.as_mut().unwrap(),
                Some(&mut headstart),
                Some(&mut datastart),
                Some(&mut dataend),
                &mut stat,
            );
            origdata = ((dataend - datastart) as f64 / 1000000.) as f32;

            /* get elapsed and cpu times need to read the uncompressed image */
            {
                let infptr = if rescale_flag {
                    tempfile.as_deref_mut().unwrap()
                } else {
                    inputfptr.as_deref_mut().unwrap()
                };
                fits_read_image_speed(
                    infptr,
                    &mut whole_elapse,
                    &mut whole_cpu,
                    &mut row_elapse,
                    &mut row_cpu,
                    &mut stat,
                );
            }

            pf(&format!(
                " {:5} {} {} {} {} {} {} {} {} {}\n",
                imagestats.n_nulls,
                dbl(c"%6.0f", imagestats.minval),
                dbl(c"%6.0f", imagestats.maxval),
                dbl(c"%8.1f", imagestats.mean),
                dbl(c"%#8.2g", imagestats.sigma),
                dbl(c"%#7.3g", imagestats.noise2),
                dbl(c"%#7.3g", imagestats.noise3),
                dbl(c"%#7.3g", imagestats.noise5),
                dbl(c"%#5.1f", f64::from(xbits)),
                dbl(c"%#6.2f", f64::from(bpix as f32 / xbits)),
            ));

            pf(
                "\n       Type   Ratio       Size (MB)     Pk (Sec) UnPk Exact ElpN CPUN  Elp1  CPU1\n",
            );

            pf(&format!(
                "       Native                                             {} {} {} {}\n",
                dbl(c"%5.3f", f64::from(whole_elapse)),
                dbl(c"%5.3f", f64::from(whole_cpu)),
                dbl(c"%5.3f", f64::from(row_elapse)),
                dbl(c"%5.3f", f64::from(row_cpu)),
            ));

            if fpvar.outfile[0] != 0 {
                report(&format!(
                    " {}  {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}",
                    String::from_utf8_lossy(cbytes(infits)),
                    extnum,
                    bitpix,
                    naxes[0],
                    naxes[1],
                    dbl(c"%#10.4g", f64::from(origdata)),
                    imagestats.n_nulls,
                    dbl(c"%#10.4g", imagestats.minval),
                    dbl(c"%#10.4g", imagestats.maxval),
                    dbl(c"%#10.4g", imagestats.mean),
                    dbl(c"%#10.4g", imagestats.sigma),
                    dbl(c"%#10.4g", imagestats.noise1),
                    dbl(c"%#10.4g", imagestats.noise2),
                    dbl(c"%#10.4g", imagestats.noise3),
                    dbl(c"%#10.4g", imagestats.noise5),
                    dbl(c"%#10.4g", f64::from(whole_elapse)),
                    dbl(c"%#10.4g", f64::from(whole_cpu)),
                    dbl(c"%#10.4g", f64::from(row_elapse)),
                    dbl(c"%#10.4g", f64::from(row_cpu)),
                ));
            }

            fits_set_lossy_int(outfptr.as_mut().unwrap(), fpvar.int_to_float, &mut stat);
            if bitpix > 0 && fpvar.int_to_float != 0 {
                if noisemin < f64::from(fpvar.n3ratio * fpvar.quantize_level)
                    || noisemin < f64::from(fpvar.n3min)
                {
                    /* image contains too little noise to quantize effectively */
                    fits_set_lossy_int(outfptr.as_mut().unwrap(), 0, &mut stat);
                    let infptr = if rescale_flag {
                        tempfile.as_deref_mut().unwrap()
                    } else {
                        inputfptr.as_deref_mut().unwrap()
                    };
                    fits_get_hdu_num(infptr, &mut hdunum);

                    pf(&format!(
                        "    HDU {hdunum} does not meet noise criteria to be quantized, so losslessly compressed.\n"
                    ));
                }
            }

            /* test compression ratio and speed for each algorithm */

            if fpvar.quantize_level != 0. {
                let outf = outfptr.as_mut().unwrap();
                fits_set_compression_type(outf, RICE_1, &mut stat);
                fits_set_tile_dim(outf, 6, &fpvar.ntile, &mut stat);
                if fpvar.no_dither != 0 {
                    fits_set_quantize_method(outf, -1, &mut stat);
                } else {
                    fits_set_quantize_method(outf, fpvar.dither_method, &mut stat);
                }

                fits_set_quantize_level(outf, fpvar.quantize_level, &mut stat);
                fits_set_dither_offset(outf, fpvar.dither_offset, &mut stat);
                fits_set_hcomp_scale(outf, fpvar.scale, &mut stat);
                fits_set_hcomp_smooth(outf, fpvar.smooth, &mut stat);

                let infptr = if rescale_flag {
                    tempfile.as_deref_mut().unwrap()
                } else {
                    inputfptr.as_deref_mut().unwrap()
                };
                fp_test_hdu(
                    infptr,
                    outfptr.as_deref_mut().unwrap(),
                    outfptr2.as_deref_mut().unwrap(),
                    fpvar,
                    &mut stat,
                );
            }

            if fpvar.quantize_level != 0. {
                let outf = outfptr.as_mut().unwrap();
                fits_set_compression_type(outf, HCOMPRESS_1, &mut stat);
                fits_set_tile_dim(outf, 6, &fpvar.ntile, &mut stat);

                if fpvar.no_dither != 0 {
                    fits_set_quantize_method(outf, -1, &mut stat);
                } else {
                    fits_set_quantize_method(outf, fpvar.dither_method, &mut stat);
                }

                fits_set_quantize_level(outf, fpvar.quantize_level, &mut stat);
                fits_set_dither_offset(outf, fpvar.dither_offset, &mut stat);
                fits_set_hcomp_scale(outf, fpvar.scale, &mut stat);
                fits_set_hcomp_smooth(outf, fpvar.smooth, &mut stat);

                let infptr = if rescale_flag {
                    tempfile.as_deref_mut().unwrap()
                } else {
                    inputfptr.as_deref_mut().unwrap()
                };
                fp_test_hdu(
                    infptr,
                    outfptr.as_deref_mut().unwrap(),
                    outfptr2.as_deref_mut().unwrap(),
                    fpvar,
                    &mut stat,
                );
            }

            {
                let outf = outfptr.as_mut().unwrap();
                if fpvar.comptype == GZIP_2 {
                    fits_set_compression_type(outf, GZIP_2, &mut stat);
                } else {
                    fits_set_compression_type(outf, GZIP_1, &mut stat);
                }

                fits_set_tile_dim(outf, 6, &fpvar.ntile, &mut stat);

                if fpvar.no_dither != 0 {
                    fits_set_quantize_method(outf, -1, &mut stat);
                } else {
                    fits_set_quantize_method(outf, fpvar.dither_method, &mut stat);
                }

                fits_set_quantize_level(outf, fpvar.quantize_level, &mut stat);
                fits_set_dither_offset(outf, fpvar.dither_offset, &mut stat);
                fits_set_hcomp_scale(outf, fpvar.scale, &mut stat);
                fits_set_hcomp_smooth(outf, fpvar.smooth, &mut stat);
            }

            {
                let infptr = if rescale_flag {
                    tempfile.as_deref_mut().unwrap()
                } else {
                    inputfptr.as_deref_mut().unwrap()
                };
                fp_test_hdu(
                    infptr,
                    outfptr.as_deref_mut().unwrap(),
                    outfptr2.as_deref_mut().unwrap(),
                    fpvar,
                    &mut stat,
                );
            }

            /* the C also has BZIP2_1, PLIO_1 and NOCOMPRESS passes here,
            commented out upstream; not ported */

            if fpvar.outfile[0] != 0 {
                report("\n");
            }

            /* delete the temporary file */
            if rescale_flag {
                fits_delete_file(&mut tempfile, &mut stat);
                if let Some(g) = tempfile3.take() {
                    /* fits_delete_file already removed it */
                    let _ = g.keep();
                }
            }
        } else if hdutype == BINARY_TBL && fpvar.do_tables != 0 {
            fits_get_num_rowsll(inputfptr.as_mut().unwrap(), &mut nrows, &mut stat);
            fits_get_num_cols(inputfptr.as_mut().unwrap(), &mut ncols, &mut stat);
            pf(&format!(
                "\n File: {}, HDU {},  {} cols X {} rows\n",
                String::from_utf8_lossy(cbytes(infits)),
                extnum,
                ncols,
                nrows
            ));
            fp_test_table(
                inputfptr.as_deref_mut().unwrap(),
                outfptr.as_deref_mut().unwrap(),
                outfptr2.as_deref_mut().unwrap(),
                fpvar,
                &mut stat,
            );
        } else {
            fits_copy_hdu(
                inputfptr.as_deref_mut().unwrap(),
                outfptr.as_deref_mut().unwrap(),
                0,
                &mut stat,
            );
            fits_copy_hdu(
                inputfptr.as_deref_mut().unwrap(),
                outfptr2.as_deref_mut().unwrap(),
                0,
                &mut stat,
            );
        }

        fits_movrel_hdu(inputfptr.as_mut().unwrap(), 1, None, &mut stat);
        extnum += 1;
    }

    if stat == END_OF_FILE {
        stat = 0;
    }

    fits_close_file(outfptr2.take().unwrap(), &mut stat);
    fits_close_file(outfptr.take().unwrap(), &mut stat);
    fits_close_file(inputfptr.take().unwrap(), &mut stat);

    if stat != 0 {
        report_error(stat);
    }
    Ok(0)
}

/// C: `snprintf(buf, sizeof buf, ...)`
fn snprintf_into(buf: &mut [c_char], text: &str) {
    let src: &[c_char] = cast_slice(text.as_bytes());
    let n = src.len().min(buf.len() - 1);
    buf[..n].copy_from_slice(&src[..n]);
    buf[n] = 0;
}

/*--------------------------------------------------------------------------*/
pub(crate) fn fp_pack_hdu(
    infptr: &mut fitsfile,
    outfptr: &mut fitsfile,
    fpvar: fpstate,
    islossless: &mut c_int,
    status: &mut c_int,
) -> c_int {
    let mut naxes: [c_long; 9] = [1, 1, 1, 1, 1, 1, 1, 1, 1];
    let mut stat: c_int = 0;
    let mut totpix: c_int = 0;
    let mut naxis: c_int = 0;
    let mut ii: c_int;
    let mut hdutype: c_int = 0;
    let mut bitpix: c_int = 0;
    let mut tstatus: c_int;
    let mut hdunum: c_int = 0;
    let mut bscale: f64 = 0.0;
    let rescale: f64;

    let mut fzalgor: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut headstart: c_long = 0;
    let mut datastart: c_long = 0;
    let mut dataend: c_long = 0;
    let datasize: c_long;
    let mut noisemin: f64;
    /* structure to hold image statistics (defined in fpack.h) */
    let mut imagestats = imgstats::default();

    if *status != 0 {
        return 0;
    }

    fits_get_hdu_type(infptr, &mut hdutype, &mut stat);

    if hdutype == IMAGE_HDU {
        fits_get_img_param(
            infptr,
            9,
            Some(&mut bitpix),
            Some(&mut naxis),
            Some(&mut naxes),
            &mut stat,
        );
        totpix = 1;
        ii = 0;
        while ii < 9 {
            totpix = totpix.wrapping_mul(naxes[ii as usize] as c_int);
            ii += 1;
        }
    }

    /* check directive keyword to see if this HDU should not be compressed */
    tstatus = 0;
    if fits_read_key(
        infptr,
        KeywordDatatypeMut::TSTRING(&mut fzalgor),
        cs!(c"FZALGOR"),
        None,
        &mut tstatus,
    ) == 0
        && (strcmp_safe(&fzalgor, cs!(c"NONE")) == 0 || strcmp_safe(&fzalgor, cs!(c"none")) == 0)
    {
        fits_copy_hdu(infptr, outfptr, 0, &mut stat);

        *status = stat;
        return 0;
    }

    /* =============================================================== */
    /* This block is only for  binary table compression */
    if hdutype == BINARY_TBL && fpvar.do_tables != 0 {
        fits_get_hduaddr(
            infptr,
            Some(&mut headstart),
            Some(&mut datastart),
            Some(&mut dataend),
            status,
        );
        datasize = dataend - datastart;

        if datasize <= 2880 {
            /* data is less than 1 FITS block in size, so don't compress */
            fits_copy_hdu(infptr, outfptr, 0, &mut stat);
        } else {
            fits_compress_table(infptr, outfptr, &mut stat);
        }

        *status = stat;
        return 0;
    }
    /* =============================================================== */

    /* If this is not a non-null image HDU, just copy it verbatim */
    if fits_is_compressed_image(infptr, &mut stat) != 0
        || hdutype != IMAGE_HDU
        || naxis == 0
        || totpix == 0
        || fpvar.do_images == 0
    {
        fits_copy_hdu(infptr, outfptr, 0, &mut stat);
    } else {
        /* remaining code deals only with IMAGE HDUs */

        /* special case: rescale a scaled integer image to reduce noise? */
        if fpvar.rescale_noise != 0. && bitpix > 0 && bitpix < LONGLONG_IMG {
            tstatus = 0;
            fits_read_key(
                infptr,
                KeywordDatatypeMut::TDOUBLE(&mut bscale),
                cs!(c"BSCALE"),
                None,
                &mut tstatus,
            );
            if tstatus == 0 && bscale != 1.0 {
                /* image must be scaled */

                if bitpix == LONG_IMG {
                    fp_i4stat(infptr, naxis, &naxes, &mut imagestats, &mut stat);
                } else {
                    fp_i2stat(infptr, naxis, &naxes, &mut imagestats, &mut stat);
                }

                /* use the minimum of the MAD 2nd, 3rd, and 5th order noise estimates */
                noisemin = imagestats.noise3;
                if imagestats.noise2 != 0. && imagestats.noise2 < noisemin {
                    noisemin = imagestats.noise2;
                }
                if imagestats.noise5 != 0. && imagestats.noise5 < noisemin {
                    noisemin = imagestats.noise5;
                }

                rescale = noisemin / f64::from(fpvar.rescale_noise);
                if rescale > 1.0 {
                    /* all the criteria are met, so create a temporary file that */
                    /* contains a rescaled version of the image, in output directory */

                    /* create temporary file name
                     *
                     * DEVIATION (upstream bug 2): the C hands fits_file_name a
                     * 513-byte buffer for a name that may be FLEN_FILENAME
                     * long. */
                    let mut outfits: [c_char; FLEN_FILENAME] = [0; FLEN_FILENAME];
                    fits_file_name(outfptr, &mut outfits, &mut stat); /* get the output file name */
                    let mut tf3: [c_char; SZ_STR] = [0; SZ_STR];
                    let guard = match fp_tmpnam_guard(cs!(c"Tmp3"), &outfits[..SZ_STR], &mut tf3) {
                        Ok(g) => g,
                        Err(FpExit(e)) => {
                            *status = e;
                            return 0;
                        }
                    };

                    let mut tempfile: Option<Box<fitsfile>> = None;
                    fits_create_file(&mut tempfile, &tf3, &mut stat);

                    fits_get_hdu_num(infptr, &mut hdunum);
                    if hdunum != 1 {
                        /* the input hdu is an image extension, so create dummy primary */
                        fits_create_img(tempfile.as_mut().unwrap(), 8, 0, &naxes, &mut stat);
                    }

                    fits_copy_header(infptr, tempfile.as_mut().unwrap(), &mut stat); /* copy the header */

                    /* rescale the data, so that it will compress more efficiently */
                    if bitpix == LONG_IMG {
                        fp_i4rescale(
                            infptr,
                            naxis,
                            &naxes,
                            rescale,
                            tempfile.as_mut().unwrap(),
                            &mut stat,
                        );
                    } else {
                        fp_i2rescale(
                            infptr,
                            naxis,
                            &naxes,
                            rescale,
                            tempfile.as_mut().unwrap(),
                            &mut stat,
                        );
                    }

                    /* scale the BSCALE keyword by the inverse factor */

                    bscale *= rescale;
                    fits_update_key(
                        tempfile.as_mut().unwrap(),
                        KeywordDatatype::TDOUBLE(&bscale),
                        cs!(c"BSCALE"),
                        None,
                        &mut stat,
                    );

                    /* rescan the header, to reset the actual scaling parameters */
                    fits_set_hdustruc(tempfile.as_mut().unwrap(), &mut stat);

                    fits_img_compress(tempfile.as_mut().unwrap(), outfptr, &mut stat);
                    fits_delete_file(&mut tempfile, &mut stat);
                    let _ = guard.keep(); /* fits_delete_file removed it */
                    *islossless = 0; /* used a lossy compression method */

                    *status = stat;
                    return 0;
                }
            }
        }

        /* if requested to do lossy compression of integer images (by */
        /* converting to float), then check if this HDU qualifies */
        if bitpix > 0 && fpvar.int_to_float != 0 {
            if bitpix >= LONG_IMG {
                fp_i4stat(infptr, naxis, &naxes, &mut imagestats, &mut stat);
            } else {
                fp_i2stat(infptr, naxis, &naxes, &mut imagestats, &mut stat);
            }

            /* rescan the image header to reset scaling values (changed by fp_iNstat) */
            ffrhdu_safe(infptr, Some(&mut hdutype), &mut stat);

            /* use the minimum of the MAD 2nd, 3rd, and 5th order noise estimates */
            noisemin = imagestats.noise3;
            if imagestats.noise2 != 0. && imagestats.noise2 < noisemin {
                noisemin = imagestats.noise2;
            }
            if imagestats.noise5 != 0. && imagestats.noise5 < noisemin {
                noisemin = imagestats.noise5;
            }

            if noisemin < f64::from(fpvar.n3ratio * fpvar.quantize_level)
                || imagestats.noise3 < f64::from(fpvar.n3min)
            {
                /* image contains too little noise to quantize effectively */
                fits_set_lossy_int(outfptr, 0, &mut stat);

                fits_get_hdu_num(infptr, &mut hdunum);

                pf(&format!(
                    "    HDU {hdunum} does not meet noise criteria to be quantized, so losslessly compressed.\n"
                ));
            } else {
                /* compressed image is not identical to original */
                *islossless = 0;
            }
        }

        /* finally, do the actual image compression */
        fits_img_compress(infptr, outfptr, &mut stat);

        if bitpix < 0 || (fpvar.comptype == HCOMPRESS_1 && fpvar.scale != 0.) {
            /* compressed image is not identical to original */
            *islossless = 0;
        }
    }

    *status = stat;
    0
}

/*--------------------------------------------------------------------------*/
pub(crate) fn fp_unpack_hdu(
    infptr: &mut fitsfile,
    outfptr: &mut fitsfile,
    _fpvar: fpstate,
    status: &mut c_int,
) -> c_int {
    let mut hdutype: c_int = 0;
    let mut lval: c_int = 0;

    if *status > 0 {
        return 0;
    }

    fits_get_hdu_type(infptr, &mut hdutype, status);

    /* =============================================================== */
    /* This block is only for beta testing of binary table compression */
    if hdutype == BINARY_TBL {
        fits_read_key(
            infptr,
            KeywordDatatypeMut::TLOGICAL(&mut lval),
            cs!(c"ZTABLE"),
            None,
            status,
        );

        if *status == 0 && lval != 0 {
            /*  uncompress the table */
            fits_uncompress_table(infptr, outfptr, status);
        } else {
            if *status == KEY_NO_EXIST {
                /* table is not compressed */
                *status = 0;
            }
            fits_copy_hdu(infptr, outfptr, 0, status);
        }

        return 0;
    /* =============================================================== */
    } else if fits_is_compressed_image(infptr, status) != 0 {
        /* uncompress the compressed image HDU */
        fits_img_decompress(infptr, outfptr, status);
    } else {
        /* not a compressed image HDU, so just copy it to the output */
        fits_copy_hdu(infptr, outfptr, 0, status);
    }

    0
}
/*--------------------------------------------------------------------------*/
/// DEVIATION: the C guards each output pointer against NULL; no call site
/// passes one, so they are plain `&mut` here.
pub(crate) fn fits_read_image_speed(
    infptr: &mut fitsfile,
    whole_elapse: &mut f32,
    whole_cpu: &mut f32,
    row_elapse: &mut f32,
    row_cpu: &mut f32,
    status: &mut c_int,
) -> c_int {
    let cnull: u8 = 0;
    let snull: c_short = 0;
    let mut bitpix: c_int = 0;
    let mut naxis: c_int = 0;
    let mut anynull: c_int = 0;
    let inull: c_int = 0;
    let mut ii: c_long;
    let mut naxes: [c_long; 9] = [0; 9];
    let mut fpixel: [c_long; 9] = [1, 1, 1, 1, 1, 1, 1, 1, 1];
    let mut lpixel: [c_long; 9] = [1, 1, 1, 1, 1, 1, 1, 1, 1];
    let inc: [c_long; 9] = [1, 1, 1, 1, 1, 1, 1, 1, 1];
    let enull: f32 = 0.0;
    let filesize: f32;
    let dnull: f64 = 0.0;

    if *status != 0 {
        return *status;
    }

    fits_get_img_param(
        infptr,
        9,
        Some(&mut bitpix),
        Some(&mut naxis),
        Some(&mut naxes),
        status,
    );

    if naxis != 2 {
        return *status;
    }

    lpixel[0] = naxes[0];
    lpixel[1] = naxes[1];

    /* filesize in MB
     *
     * NOTE (upstream bug 11): a small enough image rounds this to 0.0 and the
     * four divisions at the bottom then yield inf/nan.  Reproduced. */
    filesize = (naxes[0] as f64 * f64::from(bitpix.abs()) / 8000000. * naxes[1] as f64) as f32;

    /* measure time required to read the raw image */
    fits_set_bscale(infptr, 1.0, 0.0, status);
    *whole_elapse = 0.;
    *whole_cpu = 0.;

    let npix = (naxes[0] * naxes[1]) as usize;

    if bitpix == BYTE_IMG {
        let mut carray: Vec<u8> = vec![0; npix];

        marktime(status);
        fits_read_subset(
            infptr,
            TBYTE,
            &fpixel,
            &lpixel,
            &inc,
            Some(NullValue::UByte(cnull)),
            &mut carray,
            Some(&mut anynull),
            status,
        );

        /* get elapsped times */
        gettime(whole_elapse, whole_cpu, status);

        /* now read the image again, row by row */
        marktime(status);
        ii = 0;
        while ii < naxes[1] {
            fpixel[1] = ii + 1;
            fits_read_pix(
                infptr,
                TBYTE,
                &fpixel,
                naxes[0] as LONGLONG,
                Some(NullValue::UByte(cnull)),
                &mut carray,
                Some(&mut anynull),
                status,
            );
            ii += 1;
        }
        /* get elapsped times */
        gettime(row_elapse, row_cpu, status);
    } else if bitpix == SHORT_IMG {
        let mut sarray: Vec<c_short> = vec![0; npix];

        marktime(status);
        fits_read_subset(
            infptr,
            TSHORT,
            &fpixel,
            &lpixel,
            &inc,
            Some(NullValue::Short(snull)),
            cast_slice_mut(&mut sarray),
            Some(&mut anynull),
            status,
        );

        gettime(whole_elapse, whole_cpu, status); /* get elapsped times */

        /* now read the image again, row by row */
        marktime(status);
        ii = 0;
        while ii < naxes[1] {
            fpixel[1] = ii + 1;
            fits_read_pix(
                infptr,
                TSHORT,
                &fpixel,
                naxes[0] as LONGLONG,
                Some(NullValue::Short(snull)),
                cast_slice_mut(&mut sarray),
                Some(&mut anynull),
                status,
            );
            ii += 1;
        }
        /* get elapsped times */
        gettime(row_elapse, row_cpu, status);
    } else if bitpix == LONG_IMG {
        let mut iarray: Vec<c_int> = vec![0; npix];

        marktime(status);

        fits_read_subset(
            infptr,
            TINT,
            &fpixel,
            &lpixel,
            &inc,
            Some(NullValue::Int(inull)),
            cast_slice_mut(&mut iarray),
            Some(&mut anynull),
            status,
        );

        /* get elapsped times */
        gettime(whole_elapse, whole_cpu, status);

        /* now read the image again, row by row */
        marktime(status);
        ii = 0;
        while ii < naxes[1] {
            fpixel[1] = ii + 1;
            fits_read_pix(
                infptr,
                TINT,
                &fpixel,
                naxes[0] as LONGLONG,
                Some(NullValue::Int(inull)),
                cast_slice_mut(&mut iarray),
                Some(&mut anynull),
                status,
            );
            ii += 1;
        }
        /* get elapsped times */
        gettime(row_elapse, row_cpu, status);
    } else if bitpix == FLOAT_IMG {
        let mut earray: Vec<f32> = vec![0.0; npix];

        marktime(status);

        fits_read_subset(
            infptr,
            TFLOAT,
            &fpixel,
            &lpixel,
            &inc,
            Some(NullValue::Float(enull)),
            cast_slice_mut(&mut earray),
            Some(&mut anynull),
            status,
        );

        /* get elapsped times */
        gettime(whole_elapse, whole_cpu, status);

        /* now read the image again, row by row */
        marktime(status);
        ii = 0;
        while ii < naxes[1] {
            fpixel[1] = ii + 1;
            fits_read_pix(
                infptr,
                TFLOAT,
                &fpixel,
                naxes[0] as LONGLONG,
                Some(NullValue::Float(enull)),
                cast_slice_mut(&mut earray),
                Some(&mut anynull),
                status,
            );
            ii += 1;
        }
        /* get elapsped times */
        gettime(row_elapse, row_cpu, status);
    } else if bitpix == DOUBLE_IMG {
        let mut darray: Vec<f64> = vec![0.0; npix];

        marktime(status);

        fits_read_subset(
            infptr,
            TDOUBLE,
            &fpixel,
            &lpixel,
            &inc,
            Some(NullValue::Double(dnull)),
            cast_slice_mut(&mut darray),
            Some(&mut anynull),
            status,
        );

        /* get elapsped times */
        gettime(whole_elapse, whole_cpu, status);

        /* now read the image again, row by row */
        marktime(status);
        ii = 0;
        while ii < naxes[1] {
            fpixel[1] = ii + 1;
            fits_read_pix(
                infptr,
                TDOUBLE,
                &fpixel,
                naxes[0] as LONGLONG,
                Some(NullValue::Double(dnull)),
                cast_slice_mut(&mut darray),
                Some(&mut anynull),
                status,
            );
            ii += 1;
        }
        /* get elapsped times */
        gettime(row_elapse, row_cpu, status);
    }

    *whole_elapse /= filesize;
    *row_elapse /= filesize;
    *whole_cpu /= filesize;
    *row_cpu /= filesize;

    *status
}
/*--------------------------------------------------------------------------*/
pub(crate) fn fp_test_hdu(
    infptr: &mut fitsfile,
    outfptr: &mut fitsfile,
    outfptr2: &mut fitsfile,
    fpvar: fpstate,
    status: &mut c_int,
) -> c_int {
    /*   This routine is only used for performance testing of image HDUs. */
    /*   Use fp_test_table for testing binary table HDUs.    */

    let mut stat: c_int = 0;
    let mut hdutype: c_int = 0;
    let mut comptype: c_int = 0;
    let mut ctype: [c_char; 20] = [0; 20];
    let mut lossless: [c_char; 4] = [0; 4];
    let mut headstart: c_long = 0;
    let mut datastart: c_long = 0;
    let mut dataend: c_long = 0;
    let mut origdata: f32;
    let mut compressdata: f32;
    let mut compratio: f32 = 0.;
    let mut packcpu: f32 = 0.;
    let mut unpackcpu: f32 = 0.;
    let mut elapse: f32 = 0.;
    let mut whole_elapse: f32 = 0.;
    let mut row_elapse: f32 = 0.;
    let mut whole_cpu: f32 = 0.;
    let mut row_cpu: f32 = 0.;
    let mut datasum1: c_ulong = 0;
    let mut datasum2: c_ulong = 0;
    let mut hdusum: c_ulong = 0;

    if *status != 0 {
        return 0;
    }

    origdata = 0.;
    compressdata = 0.;
    compratio = 0.;
    lossless[0] = 0;

    fits_get_compression_type(outfptr, &mut comptype, &mut stat);
    if comptype == RICE_1 {
        strcpy_safe(&mut ctype, cs!(c"RICE"));
    } else if comptype == GZIP_1 {
        strcpy_safe(&mut ctype, cs!(c"GZIP1"));
    } else if comptype == GZIP_2 {
        strcpy_safe(&mut ctype, cs!(c"GZIP2"));
    } else if comptype == PLIO_1 {
        strcpy_safe(&mut ctype, cs!(c"PLIO"));
    } else if comptype == HCOMPRESS_1 {
        strcpy_safe(&mut ctype, cs!(c"HCOMP"));
    } else if comptype == NOCOMPRESS {
        strcpy_safe(&mut ctype, cs!(c"NONE"));
    } else {
        fp_msg_str("Error: unsupported image compression type ");
        *status = DATA_COMPRESSION_ERR;
        return 0;
    }

    /* -------------- COMPRESS the image ------------------ */

    marktime(&mut stat);

    fits_img_compress(infptr, outfptr, &mut stat);

    /* get elapsped times */
    gettime(&mut elapse, &mut packcpu, &mut stat);

    /* get elapsed and cpu times need to read the compressed image */
    fits_read_image_speed(
        outfptr,
        &mut whole_elapse,
        &mut whole_cpu,
        &mut row_elapse,
        &mut row_cpu,
        &mut stat,
    );

    if stat == 0 {
        /* -------------- UNCOMPRESS the image ------------------ */

        marktime(&mut stat);

        fits_img_decompress(outfptr, outfptr2, &mut stat);

        /* get elapsped times */
        gettime(&mut elapse, &mut unpackcpu, &mut stat);

        /* ----------------------------------------------------- */

        /* get sizes of original and compressed images */

        fits_get_hduaddr(
            infptr,
            Some(&mut headstart),
            Some(&mut datastart),
            Some(&mut dataend),
            &mut stat,
        );
        origdata = ((dataend - datastart) as f64 / 1000000.) as f32;

        fits_get_hduaddr(
            outfptr,
            Some(&mut headstart),
            Some(&mut datastart),
            Some(&mut dataend),
            &mut stat,
        );
        compressdata = ((dataend - datastart) as f64 / 1000000.) as f32;

        if compressdata != 0. {
            compratio = origdata / compressdata;
        }

        /* is this uncompressed image identical to the original? */

        fits_get_chksum(infptr, &mut datasum1, &mut hdusum, &mut stat);
        fits_get_chksum(outfptr2, &mut datasum2, &mut hdusum, &mut stat);

        if datasum1 == datasum2 {
            strcpy_safe(&mut lossless, cs!(c"Yes"));
        } else {
            strcpy_safe(&mut lossless, cs!(c"No"));
        }

        pf(&format!(
            "       {:<5} {} {} ->{} {} {} {} {} {} {} {}\n",
            String::from_utf8_lossy(cbytes(&ctype)),
            dbl(c"%6.2f", f64::from(compratio)),
            dbl(c"%7.2f", f64::from(origdata)),
            dbl(c"%7.2f", f64::from(compressdata)),
            dbl(c"%7.2f", f64::from(packcpu)),
            dbl(c"%7.2f", f64::from(unpackcpu)),
            String::from_utf8_lossy(cbytes(&lossless)),
            dbl(c"%5.3f", f64::from(whole_elapse)),
            dbl(c"%5.3f", f64::from(whole_cpu)),
            dbl(c"%5.3f", f64::from(row_elapse)),
            dbl(c"%5.3f", f64::from(row_cpu)),
        ));

        if fpvar.outfile[0] != 0 {
            report(&format!(
                " {} {} {} {} {} {} {} {}",
                dbl(c"%6.3f", f64::from(compratio)),
                dbl(c"%5.2f", f64::from(packcpu)),
                dbl(c"%5.2f", f64::from(unpackcpu)),
                String::from_utf8_lossy(cbytes(&lossless)),
                dbl(c"%7.3f", f64::from(whole_elapse)),
                dbl(c"%7.3f", f64::from(whole_cpu)),
                dbl(c"%7.3f", f64::from(row_elapse)),
                dbl(c"%7.3f", f64::from(row_cpu)),
            ));
        }

        /* delete the output HDUs to concerve disk space */

        fits_delete_hdu(outfptr, Some(&mut hdutype), &mut stat);
        fits_delete_hdu(outfptr2, Some(&mut hdutype), &mut stat);
    } else {
        pf(&format!(
            "       {:<5}     (unable to compress image)\n",
            String::from_utf8_lossy(cbytes(&ctype))
        ));
    }

    /* try to recover from any compression errors */
    if stat == DATA_COMPRESSION_ERR {
        stat = 0;
    }

    *status = stat;
    0
}
/*--------------------------------------------------------------------------*/
pub(crate) fn fp_test_table(
    infptr: &mut fitsfile,
    outfptr: &mut fitsfile,
    _outfptr2: &mut fitsfile,
    _fpvar: fpstate,
    status: &mut c_int,
) -> c_int {
    /* this routine is for performance testing of the table compression methods */

    let mut stat: c_int = 0;
    let mut hdutype: c_int = 0;
    let mut tstatus: c_int = 0;
    let mut fzalgor: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut headstart: LONGLONG = 0;
    let mut datastart: LONGLONG = 0;
    let mut dataend: LONGLONG = 0;
    let mut elapse: f32 = 0.;
    let mut cpu: f32 = 0.;

    if *status != 0 {
        return 0;
    }

    /* check directive keyword to see if this HDU should not be compressed */
    if fits_read_key(
        infptr,
        KeywordDatatypeMut::TSTRING(&mut fzalgor),
        cs!(c"FZALGOR"),
        None,
        &mut tstatus,
    ) == 0
        && (strcmp_safe(&fzalgor, cs!(c"NONE")) == 0 || strcmp_safe(&fzalgor, cs!(c"none")) == 0)
    {
        return 0;
    }

    fits_get_hduaddrll(
        infptr,
        Some(&mut headstart),
        Some(&mut datastart),
        Some(&mut dataend),
        status,
    );

    /* can't compress small tables with less than 2880 bytes of data */
    if dataend - datastart <= 2880 {
        return 0;
    }

    marktime(&mut stat);
    stat = -999; /* set special flag value */
    fits_compress_table(infptr, outfptr, &mut stat);

    /* get elapsped times */
    gettime(&mut elapse, &mut cpu, &mut stat);

    fits_delete_hdu(outfptr, Some(&mut hdutype), &mut stat);

    pf(&format!(
        "\nElapsed time = {}, cpu = {}\n",
        dbl(c"%f", f64::from(elapse)),
        dbl(c"%f", f64::from(cpu))
    ));

    report_error(stat);

    0
}
/*--------------------------------------------------------------------------*/
pub(crate) fn marktime(status: &mut c_int) -> c_int {
    /* the C reads gettimeofday on unix and gives up on elapsed time elsewhere;
    Instant is monotonic and available on both, so only the reporting in
    gettime() below keeps the platform split */
    START.set(Some(Instant::now()));
    SCPU.set(cpu_seconds());

    *status
}
/*--------------------------------------------------------------------------*/
pub(crate) fn gettime(elapse: &mut f32, elapscpu: &mut f32, status: &mut c_int) -> c_int {
    let ecpu = cpu_seconds();
    let scpu = SCPU.get();

    *elapscpu = (ecpu - scpu) as f32;

    if cfg!(unix) {
        *elapse = START
            .get()
            .map(|s| s.elapsed().as_secs_f64())
            .unwrap_or(0.0) as f32;
    } else {
        /* the C does not support high timing precision on Windows machines and
        sets the elapsed time the same as the CPU time */
        *elapse = *elapscpu;
    }

    *status
}
/*--------------------------------------------------------------------------*/
pub(crate) fn fp_i2stat(
    infptr: &mut fitsfile,
    naxis: c_int,
    naxes: &[c_long],
    imagestats: &mut imgstats,
    status: &mut c_int,
) -> c_int {
    /*
        read the central XSAMPLE by YSAMPLE region of pixels in the int*2 image,
        and then compute basic statistics: min, max, mean, sigma, mean diff, etc.
    */

    let mut fpixel: [c_long; 9] = [1, 1, 1, 1, 1, 1, 1, 1, 1];
    let mut lpixel: [c_long; 9] = [1, 1, 1, 1, 1, 1, 1, 1, 1];
    let inc: [c_long; 9] = [1, 1, 1, 1, 1, 1, 1, 1, 1];
    let (mut i1, mut i2): (c_long, c_long);
    let npix: c_long;
    let mut ngood: c_long = 0;
    let nx: c_long;
    let ny: c_long;
    let mut minvalue: c_short = 0;
    let mut maxvalue: c_short = 0;
    let mut nullvalue: c_short = 0;
    let mut anynul: c_int = 0;
    let mut tstatus: c_int;
    let mut checknull = true;
    let mut mean: f64 = 0.;
    let mut sigma: f64 = 0.;
    let mut noise1: f64 = 0.;
    let mut noise2: f64 = 0.;
    let mut noise3: f64 = 0.;
    let mut noise5: f64 = 0.;

    let xsample = c_long::from(XSAMPLE.get());
    let ysample = c_long::from(YSAMPLE.get());

    /* select the middle XSAMPLE by YSAMPLE area of the image */
    i1 = naxes[0] / 2 - (xsample / 2 - 1);
    i2 = naxes[0] / 2 + (xsample / 2);
    if i1 < 1 {
        i1 = 1;
    }
    if i2 > naxes[0] {
        i2 = naxes[0];
    }
    fpixel[0] = i1;
    lpixel[0] = i2;
    nx = i2 - i1 + 1;

    if naxis > 1 {
        i1 = naxes[1] / 2 - (ysample / 2 - 1);
        i2 = naxes[1] / 2 + (ysample / 2);
        if i1 < 1 {
            i1 = 1;
        }
        if i2 > naxes[1] {
            i2 = naxes[1];
        }
        fpixel[1] = i1;
        lpixel[1] = i2;
    }
    /* NOTE (upstream bug 6): for a 1-D image the block above never runs, so
    ny comes out equal to nx and the statistics run over nx*nx elements of a
    buffer only nx of which was read.  Reproduced. */
    ny = i2 - i1 + 1;

    npix = nx * ny;

    /* if there are higher dimensions, read the middle plane of the cube */
    if naxis > 2 {
        fpixel[2] = naxes[2] / 2 + 1;
        lpixel[2] = naxes[2] / 2 + 1;
    }

    if npix < 0 {
        *status = MEMORY_ALLOCATION;
        return *status;
    }
    let mut intarray: Vec<c_short> = vec![0; npix as usize];

    /* turn off any scaling of the integer pixel values */
    fits_set_bscale(infptr, 1.0, 0.0, status);

    fits_read_subset_sht(
        infptr,
        0,
        naxis,
        naxes,
        &fpixel,
        &lpixel,
        &inc,
        0,
        &mut intarray,
        Some(&mut anynul),
        status,
    );

    /* read the null value keyword (BLANK) if present */
    tstatus = 0;
    fits_read_key(
        infptr,
        KeywordDatatypeMut::TSHORT(&mut nullvalue),
        cs!(c"BLANK"),
        None,
        &mut tstatus,
    );
    if tstatus != 0 {
        nullvalue = 0;
        checknull = false;
    }

    /* compute statistics of the image */

    fits_img_stats_short(
        &intarray,
        nx,
        ny,
        checknull,
        nullvalue,
        Some(&mut ngood),
        Some(&mut minvalue),
        Some(&mut maxvalue),
        Some(&mut mean),
        Some(&mut sigma),
        Some(&mut noise1),
        Some(&mut noise2),
        Some(&mut noise3),
        Some(&mut noise5),
        status,
    );

    imagestats.n_nulls = (npix - ngood) as c_int;
    imagestats.minval = f64::from(minvalue);
    imagestats.maxval = f64::from(maxvalue);
    imagestats.mean = mean;
    imagestats.sigma = sigma;
    imagestats.noise1 = noise1;
    imagestats.noise2 = noise2;
    imagestats.noise3 = noise3;
    imagestats.noise5 = noise5;

    *status
}
/*--------------------------------------------------------------------------*/
pub(crate) fn fp_i4stat(
    infptr: &mut fitsfile,
    naxis: c_int,
    naxes: &[c_long],
    imagestats: &mut imgstats,
    status: &mut c_int,
) -> c_int {
    /*
        read the central XSAMPLE by YSAMPLE region of pixels in the int*2 image,
        and then compute basic statistics: min, max, mean, sigma, mean diff, etc.
    */

    let mut fpixel: [c_long; 9] = [1, 1, 1, 1, 1, 1, 1, 1, 1];
    let mut lpixel: [c_long; 9] = [1, 1, 1, 1, 1, 1, 1, 1, 1];
    let inc: [c_long; 9] = [1, 1, 1, 1, 1, 1, 1, 1, 1];
    let (mut i1, mut i2): (c_long, c_long);
    let npix: c_long;
    let mut ngood: c_long = 0;
    let nx: c_long;
    let ny: c_long;
    let mut minvalue: c_int = 0;
    let mut maxvalue: c_int = 0;
    let mut nullvalue: c_int = 0;
    let mut anynul: c_int = 0;
    let mut tstatus: c_int;
    let mut checknull = true;
    let mut mean: f64 = 0.;
    let mut sigma: f64 = 0.;
    let mut noise1: f64 = 0.;
    let mut noise2: f64 = 0.;
    let mut noise3: f64 = 0.;
    let mut noise5: f64 = 0.;

    let xsample = c_long::from(XSAMPLE.get());
    let ysample = c_long::from(YSAMPLE.get());

    /* select the middle XSAMPLE by YSAMPLE area of the image */
    i1 = naxes[0] / 2 - (xsample / 2 - 1);
    i2 = naxes[0] / 2 + (xsample / 2);
    if i1 < 1 {
        i1 = 1;
    }
    if i2 > naxes[0] {
        i2 = naxes[0];
    }
    fpixel[0] = i1;
    lpixel[0] = i2;
    nx = i2 - i1 + 1;

    if naxis > 1 {
        i1 = naxes[1] / 2 - (ysample / 2 - 1);
        i2 = naxes[1] / 2 + (ysample / 2);
        if i1 < 1 {
            i1 = 1;
        }
        if i2 > naxes[1] {
            i2 = naxes[1];
        }
        fpixel[1] = i1;
        lpixel[1] = i2;
    }
    /* see the note in fp_i2stat (upstream bug 6) */
    ny = i2 - i1 + 1;

    npix = nx * ny;

    /* if there are higher dimensions, read the middle plane of the cube */
    if naxis > 2 {
        fpixel[2] = naxes[2] / 2 + 1;
        lpixel[2] = naxes[2] / 2 + 1;
    }

    if npix < 0 {
        *status = MEMORY_ALLOCATION;
        return *status;
    }
    let mut intarray: Vec<c_int> = vec![0; npix as usize];

    /* turn off any scaling of the integer pixel values */
    fits_set_bscale(infptr, 1.0, 0.0, status);

    fits_read_subset_int(
        infptr,
        0,
        naxis,
        naxes,
        &fpixel,
        &lpixel,
        &inc,
        0,
        &mut intarray,
        Some(&mut anynul),
        status,
    );

    /* read the null value keyword (BLANK) if present */
    tstatus = 0;
    fits_read_key(
        infptr,
        KeywordDatatypeMut::TINT(&mut nullvalue),
        cs!(c"BLANK"),
        None,
        &mut tstatus,
    );
    if tstatus != 0 {
        nullvalue = 0;
        checknull = false;
    }

    /* compute statistics of the image */

    fits_img_stats_int(
        &intarray,
        nx,
        ny,
        checknull,
        nullvalue,
        Some(&mut ngood),
        Some(&mut minvalue),
        Some(&mut maxvalue),
        Some(&mut mean),
        Some(&mut sigma),
        Some(&mut noise1),
        Some(&mut noise2),
        Some(&mut noise3),
        Some(&mut noise5),
        status,
    );

    imagestats.n_nulls = (npix - ngood) as c_int;
    imagestats.minval = f64::from(minvalue);
    imagestats.maxval = f64::from(maxvalue);
    imagestats.mean = mean;
    imagestats.sigma = sigma;
    imagestats.noise1 = noise1;
    imagestats.noise2 = noise2;
    imagestats.noise3 = noise3;
    imagestats.noise5 = noise5;

    *status
}
/*--------------------------------------------------------------------------*/
pub(crate) fn fp_r4stat(
    infptr: &mut fitsfile,
    naxis: c_int,
    naxes: &[c_long],
    imagestats: &mut imgstats,
    status: &mut c_int,
) -> c_int {
    /*
        read the central XSAMPLE by YSAMPLE region of pixels in the int*2 image,
        and then compute basic statistics: min, max, mean, sigma, mean diff, etc.
    */

    let mut fpixel: [c_long; 9] = [1, 1, 1, 1, 1, 1, 1, 1, 1];
    let mut lpixel: [c_long; 9] = [1, 1, 1, 1, 1, 1, 1, 1, 1];
    let inc: [c_long; 9] = [1, 1, 1, 1, 1, 1, 1, 1, 1];
    let (mut i1, mut i2): (c_long, c_long);
    let npix: c_long;
    let mut ngood: c_long = 0;
    let nx: c_long;
    let ny: c_long;
    let mut minvalue: f32 = 0.;
    let mut maxvalue: f32 = 0.;
    let mut nullvalue: f32 = FLOATNULLVALUE;
    let mut anynul: c_int = 0;
    let mut checknull = true;
    let mut mean: f64 = 0.;
    let mut sigma: f64 = 0.;
    let mut noise1: f64 = 0.;
    let mut noise2: f64 = 0.;
    let mut noise3: f64 = 0.;
    let mut noise5: f64 = 0.;

    let xsample = c_long::from(XSAMPLE.get());
    let ysample = c_long::from(YSAMPLE.get());

    /* select the middle XSAMPLE by YSAMPLE area of the image */
    i1 = naxes[0] / 2 - (xsample / 2 - 1);
    i2 = naxes[0] / 2 + (xsample / 2);
    if i1 < 1 {
        i1 = 1;
    }
    if i2 > naxes[0] {
        i2 = naxes[0];
    }
    fpixel[0] = i1;
    lpixel[0] = i2;
    nx = i2 - i1 + 1;

    if naxis > 1 {
        i1 = naxes[1] / 2 - (ysample / 2 - 1);
        i2 = naxes[1] / 2 + (ysample / 2);
        if i1 < 1 {
            i1 = 1;
        }
        if i2 > naxes[1] {
            i2 = naxes[1];
        }
        fpixel[1] = i1;
        lpixel[1] = i2;
    }
    /* see the note in fp_i2stat (upstream bug 6) */
    ny = i2 - i1 + 1;

    npix = nx * ny;

    /* if there are higher dimensions, read the middle plane of the cube */
    if naxis > 2 {
        fpixel[2] = naxes[2] / 2 + 1;
        lpixel[2] = naxes[2] / 2 + 1;
    }

    if npix < 0 {
        *status = MEMORY_ALLOCATION;
        return *status;
    }
    let mut array: Vec<f32> = vec![0.0; npix as usize];

    fits_read_subset_flt(
        infptr,
        0,
        naxis,
        naxes,
        &fpixel,
        &lpixel,
        &inc,
        nullvalue,
        &mut array,
        Some(&mut anynul),
        status,
    );

    /* are there any null values in the array? */
    if anynul == 0 {
        nullvalue = 0.;
        checknull = false;
    }

    /* compute statistics of the image */

    fits_img_stats_float(
        &array,
        nx,
        ny,
        checknull,
        nullvalue,
        Some(&mut ngood),
        Some(&mut minvalue),
        Some(&mut maxvalue),
        Some(&mut mean),
        Some(&mut sigma),
        Some(&mut noise1),
        Some(&mut noise2),
        Some(&mut noise3),
        Some(&mut noise5),
        status,
    );

    imagestats.n_nulls = (npix - ngood) as c_int;
    imagestats.minval = f64::from(minvalue);
    imagestats.maxval = f64::from(maxvalue);
    imagestats.mean = mean;
    imagestats.sigma = sigma;
    imagestats.noise1 = noise1;
    imagestats.noise2 = noise2;
    imagestats.noise3 = noise3;
    imagestats.noise5 = noise5;

    *status
}
/*--------------------------------------------------------------------------*/
pub(crate) fn fp_i2rescale(
    infptr: &mut fitsfile,
    naxis: c_int,
    naxes: &[c_long],
    rescale: f64,
    outfptr: &mut fitsfile,
    status: &mut c_int,
) -> c_int {
    /*
        divide the integer pixel values in the input file by rescale,
        and write back out to the output file..
    */

    let mut ii: c_long;
    let mut jj: usize;
    let mut nelem: LONGLONG = 1;
    let nx: c_long;
    let mut ny: c_long;
    let mut nullvalue: c_short = 0;
    let mut anynul: c_int = 0;
    let mut tstatus: c_int;
    let mut checknull = true;

    nx = naxes[0];
    ny = 1;

    ii = 1;
    while ii < c_long::from(naxis) {
        ny *= naxes[ii as usize];
        ii += 1;
    }

    let mut intarray: Vec<c_short> = vec![0; nx as usize];

    /* read the null value keyword (BLANK) if present */
    tstatus = 0;
    fits_read_key(
        infptr,
        KeywordDatatypeMut::TSHORT(&mut nullvalue),
        cs!(c"BLANK"),
        None,
        &mut tstatus,
    );
    if tstatus != 0 {
        checknull = false;
    }

    /* turn off any scaling of the integer pixel values */
    fits_set_bscale(infptr, 1.0, 0.0, status);
    fits_set_bscale(outfptr, 1.0, 0.0, status);

    ii = 0;
    while ii < ny {
        fits_read_img_sht(
            infptr,
            1,
            nelem,
            nx as LONGLONG,
            0,
            &mut intarray,
            Some(&mut anynul),
            status,
        );

        if checknull {
            jj = 0;
            while jj < nx as usize {
                if intarray[jj] != nullvalue {
                    intarray[jj] = NSHRT(f64::from(intarray[jj]) / rescale);
                }
                jj += 1;
            }
        } else {
            jj = 0;
            while jj < nx as usize {
                intarray[jj] = NSHRT(f64::from(intarray[jj]) / rescale);
                jj += 1;
            }
        }

        fits_write_img_sht(outfptr, 1, nelem, nx as LONGLONG, &intarray, status);

        nelem += nx as LONGLONG;
        ii += 1;
    }

    *status
}
/*--------------------------------------------------------------------------*/
pub(crate) fn fp_i4rescale(
    infptr: &mut fitsfile,
    naxis: c_int,
    naxes: &[c_long],
    rescale: f64,
    outfptr: &mut fitsfile,
    status: &mut c_int,
) -> c_int {
    /*
        divide the integer pixel values in the input file by rescale,
        and write back out to the output file..
    */

    let mut ii: c_long;
    let mut jj: usize;
    let mut nelem: LONGLONG = 1;
    let nx: c_long;
    let mut ny: c_long;
    let mut nullvalue: c_int = 0;
    let mut anynul: c_int = 0;
    let mut tstatus: c_int;
    let mut checknull = true;

    nx = naxes[0];
    ny = 1;

    ii = 1;
    while ii < c_long::from(naxis) {
        ny *= naxes[ii as usize];
        ii += 1;
    }

    let mut intarray: Vec<c_int> = vec![0; nx as usize];

    /* read the null value keyword (BLANK) if present */
    tstatus = 0;
    fits_read_key(
        infptr,
        KeywordDatatypeMut::TINT(&mut nullvalue),
        cs!(c"BLANK"),
        None,
        &mut tstatus,
    );
    if tstatus != 0 {
        checknull = false;
    }

    /* turn off any scaling of the integer pixel values */
    fits_set_bscale(infptr, 1.0, 0.0, status);
    fits_set_bscale(outfptr, 1.0, 0.0, status);

    ii = 0;
    while ii < ny {
        fits_read_img_int(
            infptr,
            1,
            nelem,
            nx as LONGLONG,
            0,
            &mut intarray,
            Some(&mut anynul),
            status,
        );

        if checknull {
            jj = 0;
            while jj < nx as usize {
                if intarray[jj] != nullvalue {
                    intarray[jj] = NINT(f64::from(intarray[jj]) / rescale);
                }
                jj += 1;
            }
        } else {
            jj = 0;
            while jj < nx as usize {
                intarray[jj] = NINT(f64::from(intarray[jj]) / rescale);
                jj += 1;
            }
        }

        fits_write_img_int(outfptr, 1, nelem, nx as LONGLONG, &intarray, status);

        nelem += nx as LONGLONG;
        ii += 1;
    }

    *status
}

/* DEVIATION: `abort_fpack', the signal handler that removed the three
temp-file globals, is not ported; TempPath covers both paths. */

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use std::io::Read as _;

    /// `&str` -> the NUL-terminated `[c_char]` buffer the C works in.
    pub(crate) fn cbuf(s: &str) -> Vec<c_char> {
        let mut v: Vec<c_char> = cast_slice(s.as_bytes()).to_vec();
        v.push(0);
        v
    }

    /// The C's `argv`, for the *_get_param tests.
    pub(crate) fn argv_of(args: &[&str]) -> Argv {
        args.iter().map(|a| cbuf(a)).collect()
    }

    fn s(b: &[c_char]) -> String {
        String::from_utf8_lossy(cbytes(b)).into_owned()
    }

    /* ---------------- fp_init ------------------------------------------- */

    /// Field for field against fpackutil.c:159-206.
    #[test]
    fn test_fp_init_defaults() {
        let mut fpvar = fpstate::default();
        fp_init(&mut fpvar);

        assert_eq!(fpvar.comptype, RICE_1);
        assert_eq!(fpvar.quantize_level, 4.0);
        assert_eq!(fpvar.no_dither, 0);
        assert_eq!(fpvar.dither_method, 1);
        assert_eq!(fpvar.dither_offset, 0);
        assert_eq!(fpvar.int_to_float, 0);
        assert_eq!(fpvar.n3ratio, 2.0);
        assert_eq!(fpvar.n3min, 6.0);
        assert_eq!(fpvar.scale, 0.0);
        assert_eq!(fpvar.smooth, 0);
        assert_eq!(fpvar.rescale_noise, 0.0);
        /* ntile[0] is -1, "extent of axis"; the rest are 1 */
        assert_eq!(fpvar.ntile, [-1, 1, 1, 1, 1, 1]);
        assert_eq!(fpvar.do_checksums, 1);
        assert_eq!(fpvar.do_images, 1);
        assert_eq!(fpvar.do_tables, 0);
        assert_eq!(fpvar.delete_suffix, 0);
        assert_eq!(fpvar.firstfile, 1);
        assert_eq!(fpvar.initialized, FP_INIT_MAGIC);
        assert_eq!(fpvar.preflight_checked, 0);
        assert_eq!(fpvar.prefix[0], 0);
        assert_eq!(fpvar.extname[0], 0);
        assert_eq!(fpvar.outfile[0], 0);
    }

    /// Every entry point checks the magic first; without fp_init they exit(-1).
    #[test]
    fn test_uninitialized_fpstate_is_rejected() {
        let mut fpvar = fpstate::default();
        let argv = argv_of(&["fpack", "x.fits"]);
        assert_eq!(fp_list(2, &argv, fpvar), Err(FpExit(-1)));
        assert_eq!(fp_preflight(2, &argv, FPACK, &mut fpvar), Err(FpExit(-1)));
        assert_eq!(fp_loop(2, &argv, FPACK, fpvar), Err(FpExit(-1)));
    }

    /// fp_loop also refuses to run before fp_preflight has.
    #[test]
    fn test_fp_loop_requires_preflight() {
        let mut fpvar = fpstate::default();
        fp_init(&mut fpvar);
        let argv = argv_of(&["fpack", "x.fits"]);
        assert_eq!(fp_loop(2, &argv, FPACK, fpvar), Err(FpExit(-1)));
    }

    /* ---------------- the C stdlib shims -------------------------------- */

    #[test]
    fn test_atol_c_and_strtol_c() {
        /* (input, value, index strtol's endptr lands on) */
        assert_eq!(strtol_c(&cbuf("12")), (12, 2));
        assert_eq!(strtol_c(&cbuf("12abc")), (12, 2));
        assert_eq!(strtol_c(&cbuf("  -7")), (-7, 4));
        assert_eq!(strtol_c(&cbuf("+3")), (3, 2));
        /* no conversion: strtol leaves endptr at the start, so fp_unpack sees
        junk and falls through to the movnam_hdu branch */
        assert_eq!(strtol_c(&cbuf("abc")), (0, 0));
        assert_eq!(strtol_c(&cbuf("")), (0, 0));

        assert_eq!(atol_c(&cbuf("100")), 100);
        assert_eq!(atol_c(&cbuf("100,200")), 100);
        assert_eq!(atoi_c(&cbuf("0")), 0);
        assert_eq!(atoi_c(&cbuf("10000")), 10000);
        assert_eq!(atoi_c(&cbuf("x")), 0);
    }

    #[test]
    fn test_atof_c() {
        assert_eq!(atof_c(&cbuf("4")), 4.0);
        assert_eq!(atof_c(&cbuf("4.5")), 4.5);
        assert_eq!(atof_c(&cbuf("-2.25")), -2.25);
        assert_eq!(atof_c(&cbuf("1e3")), 1000.0);
        assert_eq!(atof_c(&cbuf("1.5e-2")), 0.015);
        assert_eq!(atof_c(&cbuf("  0.5xyz")), 0.5);
        /* atof stops at the first character that cannot extend the number,
        and yields 0.0 when that is the first one */
        assert_eq!(atof_c(&cbuf("abc")), 0.0);
        assert_eq!(atof_c(&cbuf("")), 0.0);
        /* an incomplete exponent is not consumed, as strtod does */
        assert_eq!(atof_c(&cbuf("2e")), 2.0);
        assert_eq!(atof_c(&cbuf("2e+")), 2.0);
    }

    #[test]
    fn test_strncpy_trunc_terminates() {
        let mut dst: [c_char; 5] = [9; 5];
        strncpy_trunc(&mut dst, &cbuf("abcdefgh"));
        assert_eq!(s(&dst), "abcd");
        strncpy_trunc(&mut dst, &cbuf("xy"));
        assert_eq!(s(&dst), "xy");
    }

    /* ---------------- fp_tmpnam ----------------------------------------- */

    #[test]
    fn test_fp_tmpnam_builds_rootname_plus_suffix() {
        let dir = tempfile::tempdir().unwrap();
        let root = format!("{}/out.fits", dir.path().display());

        let mut name: [c_char; SZ_STR] = [0; SZ_STR];
        fp_tmpnam(&cbuf("Tmp1"), &cbuf(&root), &mut name).unwrap();
        assert_eq!(s(&name), format!("{root}Tmp1"));
    }

    /// When the name is taken the C appends an 'x' and tries again.
    #[test]
    fn test_fp_tmpnam_appends_x_when_taken() {
        let dir = tempfile::tempdir().unwrap();
        let root = format!("{}/out.fits", dir.path().display());
        std::fs::write(format!("{root}Tmp1"), b"").unwrap();

        let mut name: [c_char; SZ_STR] = [0; SZ_STR];
        fp_tmpnam(&cbuf("Tmp1"), &cbuf(&root), &mut name).unwrap();
        assert_eq!(s(&name), format!("{root}Tmp1x"));
    }

    #[test]
    fn test_fp_tmpnam_rejects_overlong_name() {
        let root = "a".repeat(SZ_STR - 5);
        let mut name: [c_char; SZ_STR] = [0; SZ_STR];
        assert_eq!(
            fp_tmpnam(&cbuf("Tmp1"), &cbuf(&root), &mut name),
            Err(FpExit(-1))
        );
    }

    /// The guard is what replaces abort_fpack: the file goes away when the
    /// TempPath is dropped, however the caller leaves scope.
    #[test]
    fn test_fp_tmpnam_guard_removes_the_file() {
        let dir = tempfile::tempdir().unwrap();
        let root = format!("{}/out.fits", dir.path().display());
        let mut name: [c_char; SZ_STR] = [0; SZ_STR];

        let path = {
            let guard = fp_tmpnam_guard(&cbuf("Tmp2"), &cbuf(&root), &mut name).unwrap();
            let p = cpath(&name);
            std::fs::write(&p, b"partial output").unwrap();
            assert!(p.exists());
            drop(guard);
            p
        };
        assert!(!path.exists(), "the temp file outlived its guard");
    }

    /* ---------------- fp_access ----------------------------------------- */

    #[test]
    fn test_fp_access() {
        let dir = tempfile::tempdir().unwrap();
        let there = format!("{}/there", dir.path().display());
        std::fs::write(&there, b"x").unwrap();

        assert_eq!(fp_access(&cbuf(&there)), 0);
        assert_eq!(
            fp_access(&cbuf(&format!("{}/missing", dir.path().display()))),
            -1
        );
    }

    /* ---------------- fp_preflight -------------------------------------- */

    fn touch(dir: &std::path::Path, name: &str) -> String {
        let p = format!("{}/{}", dir.display(), name);
        std::fs::write(&p, b"").unwrap();
        p
    }

    fn preflight(unpack: c_int, files: &[&str], edit: impl Fn(&mut fpstate)) -> FpResult<c_int> {
        let mut fpvar = fpstate::default();
        fp_init(&mut fpvar);
        edit(&mut fpvar);
        let mut args = vec!["fpack"];
        args.extend_from_slice(files);
        let argv = argv_of(&args);
        fpvar.firstfile = 1;
        fp_preflight(argv.len() as c_int, &argv, unpack, &mut fpvar)
    }

    #[test]
    fn test_preflight_rejects_extension_notation() {
        assert_eq!(
            preflight(FPACK, &["a.fits[1]"], |_| {}),
            Err(FpExit(-1)),
            "section/extension notation is not supported"
        );
    }

    #[test]
    fn test_preflight_rejects_leading_dash_that_is_not_stdin() {
        /* a bare "-" means stdin and is fine; "-x" is a rejected file name */
        assert_eq!(preflight(FPACK, &["-x"], |_| {}), Err(FpExit(-1)));
    }

    #[test]
    fn test_preflight_rejects_overlong_name() {
        let long = "a".repeat(SZ_STR - 3);
        assert_eq!(preflight(FPACK, &[&long], |_| {}), Err(FpExit(-1)));
    }

    #[test]
    fn test_preflight_fpack_missing_input() {
        let dir = tempfile::tempdir().unwrap();
        let missing = format!("{}/nope.fits", dir.path().display());
        assert_eq!(preflight(FPACK, &[&missing], |_| {}), Err(FpExit(-1)));
    }

    #[test]
    fn test_preflight_fpack_rejects_fz_input() {
        let dir = tempfile::tempdir().unwrap();
        let f = touch(dir.path(), "a.fits.fz");
        assert_eq!(
            preflight(FPACK, &[&f], |_| {}),
            Err(FpExit(-1)),
            "fpack refuses an input that already has a .fz suffix"
        );
    }

    #[test]
    fn test_preflight_fpack_rejects_existing_output() {
        let dir = tempfile::tempdir().unwrap();
        let f = touch(dir.path(), "a.fits");
        touch(dir.path(), "a.fits.fz"); /* the name fpack would write */
        assert_eq!(preflight(FPACK, &[&f], |_| {}), Err(FpExit(-1)));
    }

    #[test]
    fn test_preflight_fpack_accepts_plain_input() {
        let dir = tempfile::tempdir().unwrap();
        let f = touch(dir.path(), "a.fits");
        assert_eq!(preflight(FPACK, &[&f], |_| {}), Ok(0));
    }

    /// -S and -T skip every output-name check.
    #[test]
    fn test_preflight_fpack_stdout_and_test_all_skip_output_checks() {
        let dir = tempfile::tempdir().unwrap();
        let f = touch(dir.path(), "a.fits");
        touch(dir.path(), "a.fits.fz");
        assert_eq!(preflight(FPACK, &[&f], |v| v.to_stdout = 1), Ok(0));
        assert_eq!(preflight(FPACK, &[&f], |v| v.test_all = 1), Ok(0));
    }

    #[test]
    fn test_preflight_fpack_one_outfile_per_run() {
        let dir = tempfile::tempdir().unwrap();
        let a = touch(dir.path(), "a.fits");
        let b = touch(dir.path(), "b.fits");
        let out = format!("{}/out.fz", dir.path().display());
        assert_eq!(
            preflight(FPACK, &[&a], |v| strcpy_safe(&mut v.outfile, &cbuf(&out))),
            Ok(0)
        );
        assert_eq!(
            preflight(FPACK, &[&a, &b], |v| strcpy_safe(
                &mut v.outfile,
                &cbuf(&out)
            )),
            Err(FpExit(-1)),
            "-O may not be reused across several inputs"
        );
    }

    /// funpack finds `name` when given `name` and `name.fz` when given either.
    #[test]
    fn test_preflight_funpack_finds_fz() {
        let dir = tempfile::tempdir().unwrap();
        let base = format!("{}/a.fits", dir.path().display());
        std::fs::write(format!("{base}.fz"), b"").unwrap();
        /* given the bare name, funpack appends .fz */
        assert_eq!(preflight(FUNPACK, &[&base], |_| {}), Ok(0));
        /* and given the full name it strips it again */
        assert_eq!(preflight(FUNPACK, &[&format!("{base}.fz")], |_| {}), Ok(0));
    }

    /// Both `a.fits` and `a.fits.fz` present: funpack cannot tell which to use.
    #[test]
    fn test_preflight_funpack_ambiguous_input() {
        let dir = tempfile::tempdir().unwrap();
        let base = touch(dir.path(), "a.fits");
        std::fs::write(format!("{base}.fz"), b"").unwrap();
        assert_eq!(preflight(FUNPACK, &[&base], |_| {}), Err(FpExit(-1)));
    }

    /// delete_suffix is set by fu_get_param; without a .fz suffix to remove
    /// there would be no output name, so funpack refuses.
    #[test]
    fn test_preflight_funpack_requires_fz_suffix_when_deleting_it() {
        let dir = tempfile::tempdir().unwrap();
        let f = touch(dir.path(), "a.fits");
        assert_eq!(
            preflight(FUNPACK, &[&f], |v| v.delete_suffix = 1),
            Err(FpExit(-1))
        );
        /* -F clears delete_suffix, and then the same input is accepted */
        assert_eq!(preflight(FUNPACK, &[&f], |v| v.delete_suffix = 0), Ok(0));
    }

    #[test]
    fn test_preflight_funpack_rejects_existing_output() {
        let dir = tempfile::tempdir().unwrap();
        let base = format!("{}/a.fits", dir.path().display());
        std::fs::write(format!("{base}.fz"), b"").unwrap();
        std::fs::write(&base, b"").unwrap();
        /* now a.fits exists as the would-be output *and* as an ambiguous input */
        assert_eq!(preflight(FUNPACK, &[&base], |_| {}), Err(FpExit(-1)));
    }

    /* ---------------- gzip_file ----------------------------------------- */

    /// funpack -Z.  The C shells out to `gzip -1`; this checks the replacement
    /// writes a real gzip stream and removes the original, which is the whole
    /// of the contract fp_loop depends on.
    #[test]
    fn test_gzip_file_replaces_the_input() {
        let dir = tempfile::tempdir().unwrap();
        let p = format!("{}/a.fits", dir.path().display());
        let payload: Vec<u8> = (0..10000u32).map(|i| (i % 251) as u8).collect();
        std::fs::write(&p, &payload).unwrap();

        gzip_file(&cbuf(&p)).unwrap();

        assert!(!std::path::Path::new(&p).exists(), "original not removed");
        let gz = std::fs::read(format!("{p}.gz")).unwrap();
        assert_eq!(&gz[..2], &[0x1f, 0x8b], "not a gzip stream");
        assert!(gz.len() < payload.len(), "did not actually compress");

        /* and it round-trips through the crate's own inflate */
        let mut out = Vec::new();
        flate2::read::GzDecoder::new(&gz[..])
            .read_to_end(&mut out)
            .unwrap();
        assert_eq!(out, payload);
    }
}
