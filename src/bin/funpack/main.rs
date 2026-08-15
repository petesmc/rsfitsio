/* Transpiled from cfitsio/utilities/funpack.c
 *
 * FUNPACK
 * R. Seaman, NOAO
 * uses fits_img_compress by W. Pence, HEASARC
 */

#![allow(non_camel_case_types, non_snake_case, non_upper_case_globals)]
#![allow(unused_assignments)]
// items only the *other* binary uses out of the two shared modules
#![allow(dead_code)]

/* The C links funpack.o against fpackutil.o; Rust binaries cannot share a
module tree, so the shared files are pulled in by path. */
#[path = "../fpack/cfmt.rs"]
mod cfmt;
#[path = "../fpack/fpack_h.rs"]
mod fpack_h;
#[path = "../fpack/fpackutil.rs"]
mod fpackutil;

use std::process::ExitCode;

use rsfitsio::bb;
use rsfitsio::c_types::{c_char, c_int};
use rsfitsio::wrappers::strlen_safe;

use crate::fpack_h::*;
use crate::fpackutil::{
    Argv, c_argv, fp_init, fp_list, fp_loop, fp_msg, fp_msg_str, fp_preflight, fp_version,
};

pub fn main() -> ExitCode {
    /* C: `exit (n)'; see fpack_h.rs. */
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
        fu_usage();
        fu_hint();
        return Err(FpExit(-1));
    }

    fp_init(&mut fpvar);
    fu_get_param(argc, &argv, &mut fpvar)?;

    if fpvar.listonly != 0 {
        fp_list(argc, &argv, fpvar)?;
    } else {
        fp_preflight(argc, &argv, FUNPACK, &mut fpvar)?;
        fp_loop(argc, &argv, FUNPACK, fpvar)?;
    }

    Ok(())
}

pub(crate) fn fu_get_param(argc: c_int, argv: &Argv, fpptr: &mut fpstate) -> FpResult<c_int> {
    let mut iarg: c_int;
    let mut tile: [c_char; SZ_STR] = [0; SZ_STR];

    if fpptr.initialized != FP_INIT_MAGIC {
        fp_msg_str("Error: internal initialization error\n");
        return Err(FpExit(-1));
    }

    tile[0] = 0;

    /* by default, .fz suffix characters to be deleted from compressed file */
    fpptr.delete_suffix = 1;

    /* flags must come first and be separately specified
     */
    iarg = 1;
    while iarg < argc {
        let arg = &argv[iarg as usize];

        if arg[0] == bb(b'-') && strlen_safe(arg) == 2 {
            if arg[1] == bb(b'F') {
                fpptr.clobber += 1;
                fpptr.delete_suffix = 0; /* no suffix in this case */
            } else if arg[1] == bb(b'D') {
                fpptr.delete_input += 1;
            } else if arg[1] == bb(b'P') {
                iarg += 1;
                if iarg >= argc {
                    fu_usage();
                    fu_hint();
                    return Err(FpExit(-1));
                } else {
                    copy_arg(&mut fpptr.prefix, &argv[iarg as usize]);
                }
            } else if arg[1] == bb(b'E') {
                iarg += 1;
                if iarg >= argc {
                    fu_usage();
                    fu_hint();
                    return Err(FpExit(-1));
                } else {
                    copy_arg(&mut fpptr.extname, &argv[iarg as usize]);
                }
            } else if arg[1] == bb(b'S') {
                fpptr.to_stdout += 1;
            } else if arg[1] == bb(b'L') {
                fpptr.listonly += 1;
            } else if arg[1] == bb(b'C') {
                fpptr.do_checksums = 0;
            } else if arg[1] == bb(b'H') {
                fu_help();
                return Err(FpExit(0));
            } else if arg[1] == bb(b'V') {
                fp_version();
                return Err(FpExit(0));
            } else if arg[1] == bb(b'Z') {
                fpptr.do_gzip_file += 1;
            } else if arg[1] == bb(b'v') {
                fpptr.verbose = 1;
            } else if arg[1] == bb(b'O') {
                iarg += 1;
                if iarg >= argc {
                    fu_usage();
                    fu_hint();
                    return Err(FpExit(-1));
                } else {
                    copy_arg(&mut fpptr.outfile, &argv[iarg as usize]);
                }
            } else {
                fp_msg_str("Error: unknown command line flag `");
                fp_msg(arg);
                fp_msg_str("'\n");
                fu_usage();
                fu_hint();
                return Err(FpExit(-1));
            }
        } else {
            break;
        }

        iarg += 1;
    }

    if fpptr.extname[0] != 0 && (fpptr.clobber != 0 || fpptr.delete_input != 0) {
        fp_msg_str("Error: -E option may not be used with -F or -D\n");
        fu_usage();
        return Err(FpExit(-1));
    }

    if fpptr.to_stdout != 0 && (fpptr.outfile[0] != 0 || fpptr.prefix[0] != 0) {
        fp_msg_str("Error: -S option may not be used with -P or -O\n");
        fu_usage();
        return Err(FpExit(-1));
    }

    if fpptr.outfile[0] != 0 && fpptr.prefix[0] != 0 {
        fp_msg_str("Error: -P and -O options may not be used together\n");
        fu_usage();
        return Err(FpExit(-1));
    }

    if iarg >= argc {
        fp_msg_str("Error: no FITS files to uncompress\n");
        fu_usage();
        return Err(FpExit(-1));
    } else {
        fpptr.firstfile = iarg;
    }

    Ok(0)
}

/// C: `strncpy (dst, argv[iarg], SZ_STR-1); dst[SZ_STR-1] = 0;`
fn copy_arg(dst: &mut [c_char; SZ_STR], src: &[c_char]) {
    let n = strlen_safe(src).min(SZ_STR - 1);
    dst[..n].copy_from_slice(&src[..n]);
    dst[n] = 0;
}

pub(crate) fn fu_usage() -> c_int {
    fp_msg_str("usage: funpack [-E <HDUlist>] [-P <pre>] [-O <name>] [-Z] -v <FITS>\n");
    fp_msg_str("more:   [-F] [-D] [-S] [-L] [-C] [-H] [-V] \n");
    0
}

pub(crate) fn fu_hint() -> c_int {
    fp_msg_str("      `funpack -H' for help\n");
    0
}

pub(crate) fn fu_help() -> c_int {
    fp_msg_str("funpack, decompress fpacked files.  Version ");
    fp_version();
    fu_usage();
    fp_msg_str("\n");

    fp_msg_str("Flags must be separate and appear before filenames:\n");
    fp_msg_str(" -E <HDUlist> Unpack only the list of HDU names or numbers in the file.\n");
    fp_msg_str(" -P <pre>    Prepend <pre> to create new output filenames.\n");
    fp_msg_str(" -O <name>   Specify full output file name.\n");
    fp_msg_str(" -Z          Recompress the output file with host GZIP program.\n");
    fp_msg_str(" -F          Overwrite input file by output file with same name.\n");
    fp_msg_str(" -D          Delete input file after writing output.\n");
    fp_msg_str(" -S          Output uncompressed file to STDOUT file stream.\n");
    fp_msg_str(" -L          List contents, files unchanged.\n");

    fp_msg_str(" -C          Don't update FITS checksum keywords.\n");

    fp_msg_str(" -v          Verbose mode; list each file as it is processed.\n");
    fp_msg_str(" -H          Show this message.\n");
    fp_msg_str(" -V          Show version number.\n");

    fp_msg_str(" \n<FITS>       FITS files to unpack; enter '-' (a hyphen) to read from stdin.\n");
    fp_msg_str(" Refer to the fpack User's Guide for more extensive help.\n");
    0
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fpackutil::fp_init;
    use crate::fpackutil::tests::{argv_of, cbuf};
    use rsfitsio::wrappers::strcmp_safe;

    fn parse(args: &[&str]) -> (FpResult<c_int>, fpstate) {
        let mut fpvar = fpstate::default();
        fp_init(&mut fpvar);
        let argv = argv_of(args);
        let r = fu_get_param(argv.len() as c_int, &argv, &mut fpvar);
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

    /// funpack.c:41 sets delete_suffix before reading any flag, and -F clears
    /// it again because the output then has the input's name.
    #[test]
    fn test_delete_suffix_default_and_clobber() {
        assert_eq!(ok(&["funpack", "a.fits.fz"]).delete_suffix, 1);
        let v = ok(&["funpack", "-F", "a.fits.fz"]);
        assert_eq!(v.clobber, 1);
        assert_eq!(v.delete_suffix, 0);
    }

    #[test]
    fn test_simple_switches() {
        assert_eq!(ok(&["funpack", "-D", "a.fz"]).delete_input, 1);
        assert_eq!(ok(&["funpack", "-S", "a.fz"]).to_stdout, 1);
        assert_eq!(ok(&["funpack", "-L", "a.fz"]).listonly, 1);
        assert_eq!(ok(&["funpack", "-C", "a.fz"]).do_checksums, 0);
        assert_eq!(ok(&["funpack", "-Z", "a.fz"]).do_gzip_file, 1);
        assert_eq!(ok(&["funpack", "-v", "a.fz"]).verbose, 1);
    }

    #[test]
    fn test_string_options() {
        let v = ok(&["funpack", "-P", "new_", "a.fz"]);
        assert_eq!(strcmp_safe(&v.prefix, &cbuf("new_")), 0);

        let v = ok(&["funpack", "-E", "1,SCI", "a.fz"]);
        assert_eq!(strcmp_safe(&v.extname, &cbuf("1,SCI")), 0);

        let v = ok(&["funpack", "-O", "out.fits", "a.fz"]);
        assert_eq!(strcmp_safe(&v.outfile, &cbuf("out.fits")), 0);
    }

    /* funpack.c:110-124 -- the three mutually exclusive combinations */

    #[test]
    fn test_extname_excludes_clobber_and_delete() {
        err(&["funpack", "-E", "SCI", "-F", "a.fz"], -1);
        err(&["funpack", "-E", "SCI", "-D", "a.fz"], -1);
        err(&["funpack", "-F", "-E", "SCI", "a.fz"], -1);
    }

    #[test]
    fn test_stdout_excludes_prefix_and_outfile() {
        err(&["funpack", "-S", "-P", "new_", "a.fz"], -1);
        err(&["funpack", "-S", "-O", "out.fits", "a.fz"], -1);
    }

    #[test]
    fn test_prefix_excludes_outfile() {
        err(&["funpack", "-P", "new_", "-O", "out.fits", "a.fz"], -1);
        err(&["funpack", "-O", "out.fits", "-P", "new_", "a.fz"], -1);
    }

    #[test]
    fn test_firstfile_marks_the_end_of_the_flags() {
        assert_eq!(ok(&["funpack", "a.fz"]).firstfile, 1);
        assert_eq!(ok(&["funpack", "-v", "a.fz"]).firstfile, 2);
        assert_eq!(
            ok(&["funpack", "-P", "p_", "-v", "a.fz", "b.fz"]).firstfile,
            4
        );
    }

    #[test]
    fn test_no_files_is_an_error() {
        err(&["funpack", "-v"], -1);
    }

    #[test]
    fn test_flag_needing_a_value_at_the_end_of_argv() {
        err(&["funpack", "-P"], -1);
        err(&["funpack", "-E"], -1);
        err(&["funpack", "-O"], -1);
    }

    #[test]
    fn test_unknown_flag() {
        err(&["funpack", "-Q", "a.fz"], -1);
    }

    /// funpack recognises only two-character flags (funpack.c:46), so fpack's
    /// long options and -q suffixes do not reach the unknown-flag branch at
    /// all: the loop breaks and they become the first *file name*.  They are
    /// rejected a step later, by fp_preflight's "invalid input file name"
    /// check on a leading '-'.  Verified against the C binary.
    #[test]
    fn test_long_options_are_taken_as_filenames_not_flags() {
        let v = ok(&["funpack", "-q0", "a.fz"]);
        assert_eq!(v.firstfile, 1, "-q0 is argv[1], the first file");
        let v = ok(&["funpack", "-tableonly", "a.fz"]);
        assert_eq!(v.firstfile, 1);
    }

    #[test]
    fn test_help_and_version_exit_zero() {
        err(&["funpack", "-H"], 0);
        err(&["funpack", "-V"], 0);
    }

    #[test]
    fn test_uninitialized_state_is_rejected() {
        let mut fpvar = fpstate::default();
        let argv = argv_of(&["funpack", "a.fz"]);
        assert_eq!(fu_get_param(2, &argv, &mut fpvar), Err(FpExit(-1)));
    }
}
