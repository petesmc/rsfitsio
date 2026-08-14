/* Transpiled from cfitsio/utilities/fpack.c
 *
 * FPACK
 * R. Seaman, NOAO, with a few enhancements by W. Pence, HEASARC
 *
 * Calls fits_img_compress in the CFITSIO library by W. Pence, HEASARC
 */

// The C's names are its own: fpstate/imgstats are lower-case types, and fpvar
// and fpptr are the locals it uses everywhere.
#![allow(non_camel_case_types, non_snake_case, non_upper_case_globals)]
// The C is full of dead stores -- `origdata = 0;` immediately before the real
// assignment, `compratio = 0.;` twice in a row -- which are kept so the
// transpile stays line-for-line with fpackutil.c.
#![allow(unused_assignments)]
// fpack_h.rs and fpackutil.rs are shared by both binaries, so each one
// carries items only the other uses (FPACK/FUNPACK, the atoi/atof shims).
#![allow(dead_code)]

mod cfmt;
mod fpack_h;
mod fpackutil;

use std::process::ExitCode;

use bytemuck::cast_slice;
use rsfitsio::c_types::{c_char, c_int, c_long};
use rsfitsio::fitsio::{GZIP_1, GZIP_2, HCOMPRESS_1, MAX_COMPRESS_DIM, NOCOMPRESS, PLIO_1, RICE_1};
use rsfitsio::wrappers::{strcmp_safe, strlen_safe, strncmp_safe};
use rsfitsio::{bb, cs};

use crate::fpack_h::*;
use crate::fpackutil::{
    Argv, atof_c, atoi_c, atol_c, c_argv, fp_init, fp_list, fp_loop, fp_msg, fp_msg_str,
    fp_preflight, fp_version,
};

/* ================================================================== */
pub fn main() -> ExitCode {
    /* C: `exit (n)' from anywhere in the call tree; see fpack_h.rs.  main()
    turns the propagated code back into the process exit status -- exit(-1) and
    ExitCode::from(255) are the same 255 to the shell. */
    match run() {
        Ok(()) => ExitCode::from(0),
        Err(FpExit(n)) => ExitCode::from(n as u8),
    }
}

fn run() -> FpResult {
    let argv: Argv = c_argv();
    let argc = argv.len() as c_int;

    let mut fpvar = fpstate::default();

    if argc <= 1 {
        fp_usage();
        fp_hint();
        return Err(FpExit(-1));
    }

    fp_init(&mut fpvar);
    fp_get_param(argc, &argv, &mut fpvar)?;

    if fpvar.listonly != 0 {
        fp_list(argc, &argv, fpvar)?;
    } else {
        fp_preflight(argc, &argv, FPACK, &mut fpvar)?;
        fp_loop(argc, &argv, FPACK, fpvar)?;
    }

    Ok(())
}
/* ================================================================== */
pub(crate) fn fp_get_param(argc: c_int, argv: &Argv, fpptr: &mut fpstate) -> FpResult<c_int> {
    let mut gottype = 0;
    let mut gottile = 0;
    let mut wholetile = 0;
    let mut iarg: c_int;
    let len: usize;
    let mut ndim: usize;
    let mut ii: usize;
    let mut doffset: c_int;
    let mut gotR = 0;
    let mut gotO = 0;
    let mut tile: [c_char; SZ_STR] = [0; SZ_STR];

    if fpptr.initialized != FP_INIT_MAGIC {
        fp_msg_str("Error: internal initialization error\n");
        return Err(FpExit(-1));
    }

    tile[0] = 0;

    /* flags must come first and be separately specified
     */
    iarg = 1;
    while iarg < argc {
        let arg = &argv[iarg as usize];

        if (arg[0] == bb(b'-') && strlen_safe(arg) == 2)
            || strncmp_safe(arg, cs!(c"-q"), 2) == 0
            || strncmp_safe(arg, cs!(c"-qz"), 3) == 0
            || strncmp_safe(arg, cs!(c"-g1"), 3) == 0
            || strncmp_safe(arg, cs!(c"-g2"), 3) == 0
            || strncmp_safe(arg, cs!(c"-i2f"), 4) == 0
            || strncmp_safe(arg, cs!(c"-n3ratio"), 8) == 0
            || strncmp_safe(arg, cs!(c"-n3min"), 6) == 0
            || strncmp_safe(arg, cs!(c"-tableonly"), 10) == 0
            || strncmp_safe(arg, cs!(c"-table"), 6) == 0
        {
            /* Rice is the default, so -r is superfluous  */
            if arg[1] == bb(b'r') {
                fpptr.comptype = RICE_1;
                if gottype != 0 {
                    fp_msg_str("Error: multiple compression flags\n");
                    fp_usage();
                    return Err(FpExit(-1));
                } else {
                    gottype += 1;
                }
            } else if arg[1] == bb(b'p') {
                fpptr.comptype = PLIO_1;
                if gottype != 0 {
                    fp_msg_str("Error: multiple compression flags\n");
                    fp_usage();
                    return Err(FpExit(-1));
                } else {
                    gottype += 1;
                }
            } else if arg[1] == bb(b'g') {
                /* test for modifiers following the 'g' */
                if arg[2] == bb(b'2') {
                    fpptr.comptype = GZIP_2;
                } else {
                    fpptr.comptype = GZIP_1;
                }

                if gottype != 0 {
                    fp_msg_str("Error: multiple compression flags\n");
                    fp_usage();
                    return Err(FpExit(-1));
                } else {
                    gottype += 1;
                }
            /* the C also has a -b/BZIP2_1 branch here, commented out
            upstream; not ported */
            } else if arg[1] == bb(b'h') {
                fpptr.comptype = HCOMPRESS_1;
                if gottype != 0 {
                    fp_msg_str("Error: multiple compression flags\n");
                    fp_usage();
                    return Err(FpExit(-1));
                } else {
                    gottype += 1;
                }
            } else if arg[1] == bb(b'd') {
                fpptr.comptype = NOCOMPRESS;
                if gottype != 0 {
                    fp_msg_str("Error: multiple compression flags\n");
                    fp_usage();
                    return Err(FpExit(-1));
                } else {
                    gottype += 1;
                }
            } else if strcmp_safe(arg, cs!(c"-i2f")) == 0 {
                /* this means convert integer images to float, and then */
                /* quantize and compress the float image.  This lossy */
                /* compression method may give higher compression than the */
                /* lossless compression method that is usually applied to */
                /* integer images. */

                fpptr.int_to_float = 1;
            } else if strcmp_safe(arg, cs!(c"-n3ratio")) == 0 {
                /* this is the minimum ratio between the MAD noise sigma */
                /* and the q parameter value in the case where the integer */
                /* image is quantized and compressed like a float image. */
                iarg += 1;
                if iarg >= argc {
                    fp_usage();
                    return Err(FpExit(-1));
                } else {
                    fpptr.n3ratio = atof_c(&argv[iarg as usize]) as f32;
                }
            } else if strcmp_safe(arg, cs!(c"-n3min")) == 0 {
                /* this is the minimum  MAD noise sigma in the case where the */
                /* integer image is quantized and compressed like a float image. */
                iarg += 1;
                if iarg >= argc {
                    fp_usage();
                    return Err(FpExit(-1));
                } else {
                    fpptr.n3min = atof_c(&argv[iarg as usize]) as f32;
                }
            } else if arg[1] == bb(b'q') {
                /* test for modifiers following the 'q' */

                if arg[2] == bb(b'z') {
                    fpptr.dither_method = 2; /* preserve zero pixels */

                    if arg[3] == bb(b't') {
                        fpptr.dither_offset = -1; /* dither based on tile checksum */
                    } else if (arg[3] as u8).is_ascii_digit() {
                        /* is a number appended to q? */
                        doffset = atoi_c(&arg[3..]);

                        if doffset == 0 {
                            fpptr.no_dither = 1; /* don't dither the quantized values */
                        } else if doffset > 0 && doffset <= 10000 {
                            fpptr.dither_offset = doffset;
                        } else {
                            fp_msg_str("Error: invalid q suffix\n");
                            fp_usage();
                            return Err(FpExit(-1));
                        }
                    }
                } else if arg[2] == bb(b't') {
                    fpptr.dither_offset = -1; /* dither based on tile checksum */
                } else if (arg[2] as u8).is_ascii_digit() {
                    /* is a number appended to q? */
                    doffset = atoi_c(&arg[2..]);

                    if doffset == 0 {
                        fpptr.no_dither = 1; /* don't dither the quantized values */
                    } else if doffset > 0 && doffset <= 10000 {
                        fpptr.dither_offset = doffset;
                    } else {
                        fp_msg_str("Error: invalid q suffix\n");
                        fp_usage();
                        return Err(FpExit(-1));
                    }
                }

                iarg += 1;
                if iarg >= argc {
                    fp_usage();
                    return Err(FpExit(-1));
                } else {
                    fpptr.quantize_level = atof_c(&argv[iarg as usize]) as f32;
                }
            } else if arg[1] == bb(b'n') {
                iarg += 1;
                if iarg >= argc {
                    fp_usage();
                    return Err(FpExit(-1));
                } else {
                    fpptr.rescale_noise = atof_c(&argv[iarg as usize]) as f32;
                }
            } else if arg[1] == bb(b's') {
                iarg += 1;
                if iarg >= argc {
                    fp_usage();
                    return Err(FpExit(-1));
                } else {
                    fpptr.scale = atof_c(&argv[iarg as usize]) as f32;
                }
            } else if strcmp_safe(arg, cs!(c"-tableonly")) == 0 {
                fpptr.do_tables = 1;
                fpptr.do_images = 0;
                /* Do not write this to stdout via fp_msg.  Otherwise it will be placed at start of piped FITS
                file, which will then be corrupted. */
                eprint!("Note: The table compression method used by fpack has been\n");
                eprint!(" officially approved as part of FITS format standard since 2016.\n");
                eprint!(" However users should be aware that the compressed table files may\n");
                eprint!(
                    " only be readable by a limited number of applications (including fpack).\n"
                );
            } else if strcmp_safe(arg, cs!(c"-table")) == 0 {
                fpptr.do_tables = 1;
                eprint!("Note: The table compression method used by fpack has been\n");
                eprint!(" officially approved as part of FITS format standard since 2016.\n");
                eprint!(" However users should be aware that the compressed table files may\n");
                eprint!(
                    " only be readable by a limited number of applications (including fpack).\n"
                );
            } else if arg[1] == bb(b't') {
                if gottile != 0 {
                    fp_msg_str("Error: multiple tile specifications\n");
                    fp_usage();
                    return Err(FpExit(-1));
                } else {
                    gottile += 1;
                }

                iarg += 1;
                if iarg >= argc {
                    fp_usage();
                    return Err(FpExit(-1));
                } else {
                    /* strncpy (tile, argv[iarg], SZ_STR-1); tile[SZ_STR-1]=0; */
                    let src = &argv[iarg as usize];
                    let n = strlen_safe(src).min(SZ_STR - 1);
                    tile[..n].copy_from_slice(&src[..n]);
                    tile[n] = 0;
                }
            } else if arg[1] == bb(b'v') {
                fpptr.verbose = 1;
            } else if arg[1] == bb(b'w') {
                wholetile += 1;
                if gottile != 0 {
                    fp_msg_str("Error: multiple tile specifications\n");
                    fp_usage();
                    return Err(FpExit(-1));
                } else {
                    gottile += 1;
                }
            } else if arg[1] == bb(b'F') {
                fpptr.clobber += 1; /* overwrite existing file */
            } else if arg[1] == bb(b'D') {
                fpptr.delete_input += 1;
            } else if arg[1] == bb(b'Y') {
                fpptr.do_not_prompt += 1;
            } else if arg[1] == bb(b'S') {
                fpptr.to_stdout += 1;
            } else if arg[1] == bb(b'L') {
                fpptr.listonly += 1;
            } else if arg[1] == bb(b'C') {
                fpptr.do_checksums = 0;
            } else if arg[1] == bb(b'T') {
                fpptr.test_all = 1;
            } else if arg[1] == bb(b'R') {
                if gotO != 0 {
                    fp_msg_str("Error: -R option is not allowed with -O\n");
                    return Err(FpExit(-1));
                }
                iarg += 1;
                if iarg >= argc {
                    fp_usage();
                    fp_hint();
                    return Err(FpExit(-1));
                } else {
                    let src = &argv[iarg as usize];
                    let n = strlen_safe(src).min(SZ_STR - 1);
                    fpptr.outfile[..n].copy_from_slice(&src[..n]);
                    fpptr.outfile[n] = 0;
                    gotR = 1;
                }
            } else if arg[1] == bb(b'H') {
                fp_help();
                return Err(FpExit(0));
            } else if arg[1] == bb(b'V') {
                fp_version();
                return Err(FpExit(0));
            } else if arg[1] == bb(b'O') {
                if gotR != 0 {
                    fp_msg_str("Error: -O option is not allowed with -R\n");
                    return Err(FpExit(-1));
                }
                iarg += 1;
                if iarg >= argc {
                    fp_usage();
                    fp_hint();
                    return Err(FpExit(-1));
                } else {
                    let src = &argv[iarg as usize];
                    let n = strlen_safe(src).min(SZ_STR - 1);
                    fpptr.outfile[..n].copy_from_slice(&src[..n]);
                    fpptr.outfile[n] = 0;
                    gotO = 1;
                }
            } else {
                fp_msg_str("Error: unknown command line flag `");
                fp_msg(arg);
                fp_msg_str("'\n");
                fp_usage();
                fp_hint();
                return Err(FpExit(-1));
            }
        } else {
            break;
        }

        iarg += 1;
    }

    /* In earlier loop, already made sure both -O and -R are not being used.
    This is essential, as each must store info in the same 'outfile' array.
    Now do additional tests of -O and -R with other flags. */

    if gotR != 0 && fpptr.test_all == 0 {
        fp_msg_str("Error: -R option may only be used with -T\n");
        return Err(FpExit(-1));
    }

    if gotO != 0 && (fpptr.test_all != 0 || fpptr.to_stdout != 0) {
        fp_msg_str("Error: -O option may not be used with -S or -T\n");
        return Err(FpExit(-1));
    }

    if fpptr.scale != 0. && fpptr.comptype != HCOMPRESS_1 && fpptr.test_all != 1 {
        fp_msg_str("Error: `-s' requires `-h or -T'\n");
        return Err(FpExit(-1));
    }

    if fpptr.quantize_level == 0. {
        if fpptr.comptype != GZIP_1 && fpptr.comptype != GZIP_2 {
            fp_msg_str("Error: `-q 0' only allowed with GZIP\n");
            return Err(FpExit(-1));
        }

        if fpptr.int_to_float == 1 {
            fp_msg_str("Error: `-q 0' not allowed with -i2f\n");
            return Err(FpExit(-1));
        }
    }

    if wholetile != 0 {
        ndim = 0;
        while ndim < MAX_COMPRESS_DIM {
            fpptr.ntile[ndim] = -1;
            ndim += 1;
        }
    } else if gottile != 0 {
        len = strlen_safe(&tile);
        ii = 0;
        ndim = 0;
        while ii < len {
            if !((tile[ii] as u8).is_ascii_digit() || tile[ii] == bb(b',')) {
                fp_msg_str("Error: `-t' requires comma separated tile dims, ");
                fp_msg_str("e.g., `-t 100,100'\n");
                return Err(FpExit(-1));
            }

            if tile[ii] == bb(b',') {
                ii += 1;
                continue;
            }

            /* DEVIATION (upstream bug 1): the C stores fpptr->ntile[ndim]
            *before* the `if (++ndim > MAX_COMPRESS_DIM)' check, so a seventh
            dimension writes one `long' past the end of the array.  The bounds
            check is hoisted above the store; it still fires on the same
            (seventh) dimension, with the same message and exit status. */
            if ndim >= MAX_COMPRESS_DIM {
                fp_msg_str("Error: too many dimensions for `-t', max=");
                fp_msg_str(&format!("{MAX_COMPRESS_DIM}\n"));
                return Err(FpExit(-1));
            }

            fpptr.ntile[ndim] = atol_c(&tile[ii..]) as c_long;
            while ii < len && (tile[ii] as u8).is_ascii_digit() {
                ii += 1;
            }

            ndim += 1;
        }
    }

    if iarg >= argc {
        fp_msg_str("Error: no FITS files to compress\n");
        fp_usage();
        return Err(FpExit(-1));
    } else {
        fpptr.firstfile = iarg;
    }

    Ok(0)
}

/* ================================================================== */
pub(crate) fn fp_usage() -> c_int {
    fp_msg_str("usage: fpack ");
    fp_msg_str("[-r|-h|-g|-p] [-w|-t <axes>] [-q <level>] [-s <scale>] [-n <noise>] -v <FITS>\n");
    fp_msg_str("more:   [-T] [-R] [-F] [-D] [-Y] [-O <file>] [-S] [-L] [-C] [-H] [-V] [-i2f]\n");
    0
}

/* ================================================================== */
pub(crate) fn fp_hint() -> c_int {
    fp_msg_str("      `fpack -H' for help\n");
    0
}

/* ================================================================== */
pub(crate) fn fp_help() -> c_int {
    fp_msg_str("fpack, a FITS image compression program.  Version ");
    fp_version();
    fp_usage();
    fp_msg_str("\n");
    fp_msg_str("NOTE: the compression parameters specified on the fpack command line may\n");
    fp_msg_str("be over-ridden by compression directive keywords in the header of each HDU\n");
    fp_msg_str("of the input file(s). See the fpack User's Guide for more details\n");
    fp_msg_str("\n");

    fp_msg_str("Flags must be separate and appear before filenames:\n");
    fp_msg_str(" -r          Rice compression [default], or\n");
    fp_msg_str(" -h          Hcompress compression, or\n");
    fp_msg_str(" -g  or -g1  GZIP_1 (per-tile) compression, or\n");
    fp_msg_str(" -g2         GZIP_2 (per-tile) compression (with byte shuffling), or\n");
    /* the C's -b/BZIP2 help line is commented out upstream; not ported */
    fp_msg_str(" -p          PLIO compression (only for positive 8 or 16-bit integer images).\n");
    fp_msg_str(" -d          Tile the image without compression (debugging mode).\n");

    fp_msg_str(" -w          Compress the whole image as a single large tile.\n");
    fp_msg_str(" -t <axes>   Comma separated list of tile dimensions [default is row by row].\n");

    fp_msg_str(" -q <level>  Quantized level spacing when converting floating point images to\n");
    fp_msg_str("             scaled integers. (+value relative to sigma of background noise;\n");
    fp_msg_str(
        "             -value is absolute). Default q value of 4 gives a compression ratio\n",
    );
    fp_msg_str("             of about 6 with very high fidelity (only 0.26% increase in noise).\n");
    fp_msg_str("             Using q values of  2, or 1 will give compression ratios of\n");
    fp_msg_str("             about 8, or 10, respectively (with 1.0% or 4.1% noise increase).\n");
    fp_msg_str("             The scaled quantized values are randomly dithered using a seed \n");
    fp_msg_str("             value determined from the system clock at run time.\n");
    fp_msg_str("             Use -q0 instead of -q to suppress random dithering.\n");
    fp_msg_str("             Use -qz instead of -q to not dither zero-valued pixels.\n");
    fp_msg_str(
        "             Use -qt or -qzt to compute random dithering seed from first tile checksum.\n",
    );
    fp_msg_str(
        "             Use -qN or -qzN, (N in range 1 to 10000) to use a specific dithering seed.\n",
    );
    fp_msg_str("             Floating-point images can be losslessly compressed by selecting\n");
    fp_msg_str(
        "             the GZIP algorithm and specifying -q 0, but this is slower and often\n",
    );
    fp_msg_str(
        "             produces much less compression than the default quantization method.\n",
    );
    fp_msg_str(
        " -i2f        Convert integer images to floating point, then quantize and compress\n",
    );
    fp_msg_str("             using the specified q level.  When used appropriately, this lossy\n");
    fp_msg_str(
        "             compression method can give much better compression than the normal\n",
    );
    fp_msg_str(
        "             lossless compression methods without significant loss of information.\n",
    );
    fp_msg_str(
        "             The -n3ratio and -n3min flags control the minimum noise thresholds;\n",
    );
    fp_msg_str("             Images below these thresholds will be losslessly compressed.\n");
    fp_msg_str(
        " -n3ratio    Minimum ratio of background noise sigma divided by q.  Default = 2.0.\n",
    );
    fp_msg_str(
        " -n3min      Minimum background noise sigma. Default = 6. The -i2f flag will be ignored\n",
    );
    fp_msg_str("             if the noise level in the image does not exceed both thresholds.\n");
    fp_msg_str(" -s <scale>  Scale factor for lossy Hcompress [default = 0 = lossless]\n");
    fp_msg_str("             (+values relative to RMS noise; -value is absolute)\n");
    fp_msg_str(
        " -n <noise>  Rescale scaled-integer images to reduce noise and improve compression.\n",
    );
    fp_msg_str(" -v          Verbose mode; list each file as it is processed.\n");
    fp_msg_str(
        " -T          Show compression algorithm comparison test statistics; files unchanged.\n",
    );
    fp_msg_str(" -R <file>   Write the comparison test report (above) to a text file.\n");
    fp_msg_str(" -table      Compress FITS binary tables as well as compress any image HDUs.\n");
    fp_msg_str(" -tableonly  Compress only FITS binary tables; do not compress any image HDUs.\n");
    fp_msg_str("             \n");

    fp_msg_str("\nkeywords shared with funpack:\n");

    fp_msg_str(" -F          Overwrite input file by output file with same name.\n");
    fp_msg_str(" -D          Delete input file after writing output.\n");
    fp_msg_str(" -Y          Suppress prompts to confirm -F or -D options.\n");

    fp_msg_str(" -O <file>   Specify full output file name. This may be used only when fpack\n");
    fp_msg_str("               is run on a single input file.\n");
    fp_msg_str(" -S          Output compressed FITS files to STDOUT.\n");
    fp_msg_str(" -L          List contents; files unchanged.\n");

    fp_msg_str(" -C          Don't update FITS checksum keywords.\n");

    fp_msg_str(" -H          Show this message.\n");
    fp_msg_str(" -V          Show version number.\n");

    fp_msg_str(
        "\n <FITS>      FITS files to pack; enter '-' (a hyphen) to read input from stdin stream.\n",
    );
    fp_msg_str(" Refer to the fpack User's Guide for more extensive help.\n");
    0
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fpackutil::fp_init;
    use crate::fpackutil::tests::argv_of;

    /// Run fp_get_param over `args` (which must start with the program name)
    /// on a freshly fp_init'ed state, as main() does.
    fn parse(args: &[&str]) -> (FpResult<c_int>, fpstate) {
        let mut fpvar = fpstate::default();
        fp_init(&mut fpvar);
        let argv = argv_of(args);
        let r = fp_get_param(argv.len() as c_int, &argv, &mut fpvar);
        (r, fpvar)
    }

    fn ok(args: &[&str]) -> fpstate {
        let (r, v) = parse(args);
        assert_eq!(r, Ok(0), "expected {args:?} to be accepted");
        v
    }

    fn err(args: &[&str], code: c_int) {
        let (r, _) = parse(args);
        assert_eq!(r, Err(FpExit(code)), "expected {args:?} to exit({code})");
    }

    /* ---------------- compression flags --------------------------------- */

    #[test]
    fn test_compression_flags() {
        assert_eq!(ok(&["fpack", "a.fits"]).comptype, RICE_1, "Rice is default");
        assert_eq!(ok(&["fpack", "-r", "a.fits"]).comptype, RICE_1);
        assert_eq!(ok(&["fpack", "-p", "a.fits"]).comptype, PLIO_1);
        assert_eq!(ok(&["fpack", "-h", "a.fits"]).comptype, HCOMPRESS_1);
        assert_eq!(ok(&["fpack", "-d", "a.fits"]).comptype, NOCOMPRESS);
        /* -g and -g1 are GZIP_1; only the '2' modifier selects GZIP_2 */
        assert_eq!(ok(&["fpack", "-g", "a.fits"]).comptype, GZIP_1);
        assert_eq!(ok(&["fpack", "-g1", "a.fits"]).comptype, GZIP_1);
        assert_eq!(ok(&["fpack", "-g2", "a.fits"]).comptype, GZIP_2);
    }

    #[test]
    fn test_multiple_compression_flags_rejected() {
        err(&["fpack", "-r", "-g", "a.fits"], -1);
        err(&["fpack", "-h", "-p", "a.fits"], -1);
        err(&["fpack", "-g", "-g2", "a.fits"], -1);
    }

    /* ---------------- the -q suffix matrix ------------------------------ */

    /* fpack.c:136-179.  The suffix selects the dithering; the *next*
    argument is always the quantization level. */
    #[test]
    fn test_q_suffixes() {
        let v = ok(&["fpack", "-q", "2", "a.fits"]);
        assert_eq!(v.quantize_level, 2.0);
        assert_eq!(v.no_dither, 0);
        assert_eq!(v.dither_offset, 0);
        assert_eq!(v.dither_method, 1);

        /* -q0 means "do not dither"; the 0 is a suffix, not the level */
        let v = ok(&["fpack", "-q0", "4", "a.fits"]);
        assert_eq!(v.no_dither, 1);
        assert_eq!(v.quantize_level, 4.0);

        /* -qN picks a specific seed */
        assert_eq!(ok(&["fpack", "-q5", "4", "a.fits"]).dither_offset, 5);
        assert_eq!(
            ok(&["fpack", "-q10000", "4", "a.fits"]).dither_offset,
            10000
        );

        /* -qt seeds from the first tile checksum */
        assert_eq!(ok(&["fpack", "-qt", "4", "a.fits"]).dither_offset, -1);

        /* -qz preserves zero pixels, and composes with t and N */
        let v = ok(&["fpack", "-qz", "4", "a.fits"]);
        assert_eq!(v.dither_method, 2);
        assert_eq!(v.dither_offset, 0);
        assert_eq!(v.quantize_level, 4.0);

        let v = ok(&["fpack", "-qzt", "4", "a.fits"]);
        assert_eq!(v.dither_method, 2);
        assert_eq!(v.dither_offset, -1);

        let v = ok(&["fpack", "-qz3", "4", "a.fits"]);
        assert_eq!(v.dither_method, 2);
        assert_eq!(v.dither_offset, 3);

        let v = ok(&["fpack", "-qz0", "4", "a.fits"]);
        assert_eq!(v.dither_method, 2);
        assert_eq!(v.no_dither, 1);
    }

    #[test]
    fn test_q_suffix_out_of_range_rejected() {
        err(&["fpack", "-q10001", "4", "a.fits"], -1);
        err(&["fpack", "-qz10001", "4", "a.fits"], -1);
    }

    /// `-q 0` (a *level* of zero) means lossless, which only GZIP can do.
    #[test]
    fn test_q_zero_level_only_with_gzip() {
        err(&["fpack", "-q", "0", "a.fits"], -1);
        err(&["fpack", "-r", "-q", "0", "a.fits"], -1);
        assert_eq!(
            ok(&["fpack", "-g", "-q", "0", "a.fits"]).quantize_level,
            0.0
        );
        assert_eq!(
            ok(&["fpack", "-g2", "-q", "0", "a.fits"]).quantize_level,
            0.0
        );
        /* ...and never together with -i2f */
        err(&["fpack", "-g", "-i2f", "-q", "0", "a.fits"], -1);
    }

    /* ---------------- tile dimensions ----------------------------------- */

    #[test]
    fn test_tile_dims() {
        assert_eq!(
            ok(&["fpack", "-t", "100,100", "a.fits"]).ntile,
            [100, 100, 1, 1, 1, 1],
            "unspecified dimensions keep the fp_init default of 1"
        );
        assert_eq!(
            ok(&["fpack", "-t", "1,2,3,4,5,6", "a.fits"]).ntile,
            [1, 2, 3, 4, 5, 6],
            "MAX_COMPRESS_DIM dimensions is the most that fits"
        );
        /* -w is "one tile for the whole image": -1 in every dimension */
        assert_eq!(ok(&["fpack", "-w", "a.fits"]).ntile, [-1; MAX_COMPRESS_DIM]);
    }

    /// Upstream bug 1: the C stores ntile[6] -- one past the end -- before
    /// noticing there are too many dimensions.  The bounds check is hoisted,
    /// so the seventh dimension is refused instead of corrupting the struct,
    /// with the same exit status the C reaches a moment later.
    #[test]
    fn test_seventh_tile_dimension_is_refused_not_written() {
        err(&["fpack", "-t", "1,2,3,4,5,6,7", "a.fits"], -1);
    }

    #[test]
    fn test_tile_dims_must_be_digits_and_commas() {
        err(&["fpack", "-t", "100x100", "a.fits"], -1);
        err(&["fpack", "-t", "-5", "a.fits"], -1);
    }

    #[test]
    fn test_multiple_tile_specifications_rejected() {
        err(&["fpack", "-w", "-t", "10,10", "a.fits"], -1);
        err(&["fpack", "-t", "10,10", "-w", "a.fits"], -1);
        err(&["fpack", "-t", "10,10", "-t", "20,20", "a.fits"], -1);
    }

    /* ---------------- -s, -n, -i2f, the noise thresholds ---------------- */

    #[test]
    fn test_scale_requires_hcompress_or_test_all() {
        err(&["fpack", "-s", "1", "a.fits"], -1);
        assert_eq!(ok(&["fpack", "-h", "-s", "1", "a.fits"]).scale, 1.0);
        assert_eq!(ok(&["fpack", "-T", "-s", "1", "a.fits"]).scale, 1.0);
        /* a scale of 0 is the default and needs no -h */
        assert_eq!(ok(&["fpack", "-s", "0", "a.fits"]).scale, 0.0);
    }

    #[test]
    fn test_noise_and_i2f_options() {
        assert_eq!(ok(&["fpack", "-n", "3.5", "a.fits"]).rescale_noise, 3.5);
        assert_eq!(ok(&["fpack", "-i2f", "a.fits"]).int_to_float, 1);
        assert_eq!(ok(&["fpack", "-n3ratio", "1.5", "a.fits"]).n3ratio, 1.5);
        assert_eq!(ok(&["fpack", "-n3min", "8", "a.fits"]).n3min, 8.0);
    }

    /* ---------------- -O and -R exclusions ------------------------------ */

    #[test]
    fn test_output_file_options() {
        /* -R is the -T report file, and is only meaningful with -T */
        err(&["fpack", "-R", "rep.txt", "a.fits"], -1);
        assert!(ok(&["fpack", "-T", "-R", "rep.txt", "a.fits"]).outfile[0] != 0);

        /* -O names the output file, and cannot combine with -S or -T */
        assert!(ok(&["fpack", "-O", "out.fz", "a.fits"]).outfile[0] != 0);
        err(&["fpack", "-O", "out.fz", "-S", "a.fits"], -1);
        err(&["fpack", "-O", "out.fz", "-T", "a.fits"], -1);

        /* the two share `outfile', so they can never be used together */
        err(&["fpack", "-R", "rep.txt", "-O", "out.fz", "a.fits"], -1);
        err(&["fpack", "-O", "out.fz", "-R", "rep.txt", "a.fits"], -1);
    }

    /* ---------------- the plain switches -------------------------------- */

    #[test]
    fn test_simple_switches() {
        assert_eq!(ok(&["fpack", "-v", "a.fits"]).verbose, 1);
        assert_eq!(ok(&["fpack", "-F", "a.fits"]).clobber, 1);
        assert_eq!(ok(&["fpack", "-D", "a.fits"]).delete_input, 1);
        assert_eq!(ok(&["fpack", "-Y", "a.fits"]).do_not_prompt, 1);
        assert_eq!(ok(&["fpack", "-S", "a.fits"]).to_stdout, 1);
        assert_eq!(ok(&["fpack", "-L", "a.fits"]).listonly, 1);
        assert_eq!(ok(&["fpack", "-C", "a.fits"]).do_checksums, 0);
        assert_eq!(ok(&["fpack", "-T", "a.fits"]).test_all, 1);

        let v = ok(&["fpack", "-table", "a.fits"]);
        assert_eq!((v.do_tables, v.do_images), (1, 1));
        let v = ok(&["fpack", "-tableonly", "a.fits"]);
        assert_eq!((v.do_tables, v.do_images), (1, 0));
    }

    /* ---------------- where the flags stop ------------------------------ */

    #[test]
    fn test_firstfile_marks_the_end_of_the_flags() {
        assert_eq!(ok(&["fpack", "a.fits"]).firstfile, 1);
        assert_eq!(ok(&["fpack", "-v", "a.fits"]).firstfile, 2);
        assert_eq!(
            ok(&["fpack", "-q", "4", "-v", "a.fits", "b.fits"]).firstfile,
            4
        );
    }

    #[test]
    fn test_no_files_is_an_error() {
        err(&["fpack", "-v"], -1);
        err(&["fpack", "-r"], -1);
    }

    #[test]
    fn test_flag_needing_a_value_at_the_end_of_argv() {
        err(&["fpack", "-q"], -1);
        err(&["fpack", "-t"], -1);
        err(&["fpack", "-s"], -1);
        err(&["fpack", "-n"], -1);
        err(&["fpack", "-O"], -1);
        err(&["fpack", "-R"], -1);
        err(&["fpack", "-n3min"], -1);
        err(&["fpack", "-n3ratio"], -1);
    }

    #[test]
    fn test_unknown_flag() {
        err(&["fpack", "-Q", "a.fits"], -1);
        err(&["fpack", "-Z", "a.fits"], -1); /* -Z is funpack's, not fpack's */
    }

    /// -H and -V print and exit(0); as an Err they still stop the run.
    #[test]
    fn test_help_and_version_exit_zero() {
        err(&["fpack", "-H"], 0);
        err(&["fpack", "-V"], 0);
    }

    #[test]
    fn test_uninitialized_state_is_rejected() {
        let mut fpvar = fpstate::default();
        let argv = argv_of(&["fpack", "a.fits"]);
        assert_eq!(
            fp_get_param(2, &argv, &mut fpvar),
            Err(FpExit(-1)),
            "fp_get_param checks the FP_INIT_MAGIC first"
        );
    }
}
