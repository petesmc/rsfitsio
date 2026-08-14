/* Transpiled from cfitsio/utilities/fpack.h
 *
 * used by FPACK and FUNPACK
 * R. Seaman, NOAO
 * W. Pence, NASA/GSFC
 *
 * The C header is just constants, two structs and the prototypes; the
 * prototypes are carried by Rust's module system instead.
 *
 * DEVIATION: each `exit(n)' becomes `return Err(FpExit(n))', propagated to
 * main().  This keeps every error path reachable from a unit test instead of
 * killing the test runner, and is the only structural change to the control
 * flow.
 */

use rsfitsio::c_types::{c_char, c_float, c_int, c_long};
use rsfitsio::fitsio::MAX_COMPRESS_DIM;

pub(crate) const FPACK_VERSION: &str = "1.7.0 (Dec 2013)";

pub(crate) const FP_INIT_MAGIC: c_int = 42;
pub(crate) const FPACK: c_int = 0;
pub(crate) const FUNPACK: c_int = 1;

/* changed from 16 in Jan. 2010 */
pub(crate) const DEF_QLEVEL: c_float = 4.;

pub(crate) const DEF_HCOMP_SCALE: c_float = 0.;
pub(crate) const DEF_HCOMP_SMOOTH: c_int = 0;
pub(crate) const DEF_RESCALE_NOISE: c_float = 0.;

pub(crate) const SZ_STR: usize = 513;
pub(crate) const SZ_CARD: usize = 81;

#[derive(Clone, Copy)]
pub(crate) struct fpstate {
    pub comptype: c_int,
    pub quantize_level: c_float,
    pub no_dither: c_int,
    pub dither_offset: c_int,
    pub dither_method: c_int,
    pub scale: c_float,
    pub rescale_noise: c_float,
    pub smooth: c_int,
    pub int_to_float: c_int,
    pub n3ratio: c_float,
    pub n3min: c_float,
    pub ntile: [c_long; MAX_COMPRESS_DIM],

    pub to_stdout: c_int,
    pub listonly: c_int,
    pub clobber: c_int,
    pub delete_input: c_int,
    pub do_not_prompt: c_int,
    pub do_checksums: c_int,
    pub do_gzip_file: c_int,
    pub do_images: c_int,
    pub do_tables: c_int,
    pub test_all: c_int,
    pub verbose: c_int,

    pub prefix: [c_char; SZ_STR],
    pub extname: [c_char; SZ_STR],
    pub delete_suffix: c_int,
    pub outfile: [c_char; SZ_STR],
    pub firstfile: c_int,

    pub initialized: c_int,
    pub preflight_checked: c_int,
}

impl Default for fpstate {
    /* the C leaves fpstate uninitialised until fp_init(), which is the only
    constructor either program uses and checks the magic */
    fn default() -> Self {
        fpstate {
            comptype: 0,
            quantize_level: 0.,
            no_dither: 0,
            dither_offset: 0,
            dither_method: 0,
            scale: 0.,
            rescale_noise: 0.,
            smooth: 0,
            int_to_float: 0,
            n3ratio: 0.,
            n3min: 0.,
            ntile: [0; MAX_COMPRESS_DIM],
            to_stdout: 0,
            listonly: 0,
            clobber: 0,
            delete_input: 0,
            do_not_prompt: 0,
            do_checksums: 0,
            do_gzip_file: 0,
            do_images: 0,
            do_tables: 0,
            test_all: 0,
            verbose: 0,
            prefix: [0; SZ_STR],
            extname: [0; SZ_STR],
            delete_suffix: 0,
            outfile: [0; SZ_STR],
            firstfile: 0,
            initialized: 0,
            preflight_checked: 0,
        }
    }
}

#[derive(Clone, Copy, Default)]
pub(crate) struct imgstats {
    pub n_nulls: c_int,
    pub minval: f64,
    pub maxval: f64,
    pub mean: f64,
    pub sigma: f64,
    pub noise1: f64,
    pub noise2: f64,
    pub noise3: f64,
    pub noise5: f64,
}

/* ---------------------------------------------------------------------- */
/* C: exit(n).  See the DEVIATION note at the top of this file.            */

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct FpExit(pub c_int);

pub(crate) type FpResult<T = ()> = Result<T, FpExit>;
