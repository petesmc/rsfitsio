//! End-to-end tests for the `fpack` and `funpack` binaries, driving the real
//! executables so exit status and both streams are covered.  Everything below
//! the CLI is covered by the `#[cfg(test)]` modules inside `src/bin/fpack/`,
//! and `scripts/fpack_diff.sh` compares against the C over a wider matrix.

mod common;

use std::path::Path;
use std::process::{Command, Output};

use libc::c_int;
use rsfitsio::aliases::rust_api::*;
use rsfitsio::fitsio::{
    DOUBLE_IMG, FLOAT_IMG, LONG_IMG, LONGLONG, SHORT_IMG, USHORT_IMG, fitsfile,
};
use tempfile::TempDir;

const FPACK: &str = env!("CARGO_BIN_EXE_fpack");
const FUNPACK: &str = env!("CARGO_BIN_EXE_funpack");

const NX: LONGLONG = 60;
const NY: LONGLONG = 40;
const NELEM: usize = (NX * NY) as usize;

/* ---------------------------------------------------------------------- */
/* harness                                                                 */
/* ---------------------------------------------------------------------- */

fn run(exe: &str, dir: &Path, args: &[&str]) -> Output {
    Command::new(exe)
        .args(args)
        .current_dir(dir)
        .output()
        .unwrap_or_else(|e| panic!("failed to run {exe}: {e}"))
}

fn stdout(o: &Output) -> String {
    String::from_utf8_lossy(&o.stdout).into_owned()
}

fn status(o: &Output) -> i32 {
    o.status.code().expect("process was killed by a signal")
}

/// C's `exit(-1)`, as the shell sees it.
const EXIT_FAILURE: i32 = 255;

fn pixels(seed: u32) -> Vec<i32> {
    let mut x = seed.wrapping_mul(2_654_435_761).wrapping_add(1);
    (0..NELEM)
        .map(|_| {
            x ^= x << 13;
            x ^= x >> 17;
            x ^= x << 5;
            (x % 4000) as i32
        })
        .collect()
}

/// Write a 2-D image of `bitpix` into `path` and return its byte length.
fn make_image(path: &str, bitpix: c_int) {
    use bytemuck::cast_slice;
    use std::ffi::CString;

    let cpath = CString::new(path).unwrap();
    let name = cast_slice(cpath.to_bytes_with_nul());
    let mut status: c_int = 0;
    let mut fptr: Option<Box<fitsfile>> = None;

    fits_create_diskfile(&mut fptr, name, &mut status);
    let f = fptr.as_mut().unwrap();
    fits_create_imgll(f, bitpix, 2, &[NX, NY], &mut status);

    let v = pixels(bitpix.unsigned_abs());
    match bitpix {
        SHORT_IMG => {
            let d: Vec<i16> = v.iter().map(|&x| (x - 2000) as i16).collect();
            fits_write_img_sht(f, 1, 1, NELEM as LONGLONG, &d, &mut status);
        }
        USHORT_IMG => {
            let d: Vec<u16> = v.iter().map(|&x| x as u16).collect();
            fits_write_img_usht(f, 1, 1, NELEM as LONGLONG, &d, &mut status);
        }
        LONG_IMG => {
            fits_write_img_int(f, 1, 1, NELEM as LONGLONG, &v, &mut status);
        }
        FLOAT_IMG => {
            let d: Vec<f32> = v.iter().map(|&x| x as f32 * 0.25).collect();
            fits_write_img_flt(f, 1, 1, NELEM as LONGLONG, &d, &mut status);
        }
        DOUBLE_IMG => {
            let d: Vec<f64> = v.iter().map(|&x| f64::from(x) * 0.125).collect();
            fits_write_img_dbl(f, 1, 1, NELEM as LONGLONG, &d, &mut status);
        }
        other => panic!("no fixture for BITPIX {other}"),
    }
    fits_close_file(fptr.take().unwrap(), &mut status);
    assert_eq!(status, 0, "failed to build the {bitpix} fixture");
}

/// A non-negative `SHORT_IMG` fixture, for the PLIO cases: PLIO rejects
/// negatives, and `USHORT_IMG` would drag in the BZERO scaling that
/// `test_c_packed_ushort_is_unreadable` is about.
fn nonneg_fixture() -> TempDir {
    use bytemuck::cast_slice;
    use std::ffi::CString;

    let dir = tempfile::Builder::new()
        .prefix("fpack-nn-")
        .tempdir()
        .unwrap();
    let path = dir.path().join("a.fits");
    let cpath = CString::new(path.to_str().unwrap()).unwrap();
    let name = cast_slice(cpath.to_bytes_with_nul());
    let mut status: c_int = 0;
    let mut fptr: Option<Box<fitsfile>> = None;

    fits_create_diskfile(&mut fptr, name, &mut status);
    let f = fptr.as_mut().unwrap();
    fits_create_imgll(f, SHORT_IMG, 2, &[NX, NY], &mut status);
    let d: Vec<i16> = pixels(3).iter().map(|&x| x as i16).collect();
    fits_write_img_sht(f, 1, 1, NELEM as LONGLONG, &d, &mut status);
    fits_close_file(fptr.take().unwrap(), &mut status);
    assert_eq!(status, 0, "could not build the non-negative fixture");
    dir
}

/// A tempdir holding `a.fits` of the given BITPIX.
fn fixture(bitpix: c_int) -> TempDir {
    let dir = tempfile::Builder::new().prefix("fpack-").tempdir().unwrap();
    make_image(dir.path().join("a.fits").to_str().unwrap(), bitpix);
    dir
}

/* ---------------------------------------------------------------------- */
/* usage, help, version, bad input                                         */
/* ---------------------------------------------------------------------- */

#[test]
fn test_no_arguments_prints_usage_and_fails() {
    let dir = tempfile::tempdir().unwrap();

    let o = run(FPACK, dir.path(), &[]);
    assert_eq!(status(&o), EXIT_FAILURE);
    assert!(stdout(&o).starts_with("usage: fpack "), "{}", stdout(&o));
    assert!(stdout(&o).contains("`fpack -H' for help"));

    let o = run(FUNPACK, dir.path(), &[]);
    assert_eq!(status(&o), EXIT_FAILURE);
    assert!(stdout(&o).starts_with("usage: funpack "), "{}", stdout(&o));
    assert!(stdout(&o).contains("`funpack -H' for help"));
}

#[test]
fn test_help_and_version_succeed() {
    let dir = tempfile::tempdir().unwrap();

    for (exe, banner) in [
        (FPACK, "fpack, a FITS image compression program."),
        (FUNPACK, "funpack, decompress fpacked files."),
    ] {
        let o = run(exe, dir.path(), &["-H"]);
        assert_eq!(status(&o), 0, "-H exits 0");
        assert!(stdout(&o).starts_with(banner), "{}", stdout(&o));

        let o = run(exe, dir.path(), &["-V"]);
        assert_eq!(status(&o), 0, "-V exits 0");
        /* FPACK_VERSION, then the CFITSIO version fp_version appends */
        assert!(stdout(&o).starts_with("1.7.0 (Dec 2013) CFITSIO version "));
    }
}

#[test]
fn test_unknown_flag_is_reported() {
    let dir = tempfile::tempdir().unwrap();

    let o = run(FPACK, dir.path(), &["-Q", "a.fits"]);
    assert_eq!(status(&o), EXIT_FAILURE);
    assert!(stdout(&o).contains("Error: unknown command line flag `-Q'"));

    let o = run(FUNPACK, dir.path(), &["-Q", "a.fits"]);
    assert_eq!(status(&o), EXIT_FAILURE);
    assert!(stdout(&o).contains("Error: unknown command line flag `-Q'"));
}

#[test]
fn test_missing_input_file() {
    let dir = tempfile::tempdir().unwrap();
    let o = run(FPACK, dir.path(), &["nope.fits"]);
    assert_eq!(status(&o), EXIT_FAILURE);
    assert!(stdout(&o).contains("Error: can't find or read input file nope.fits"));
    assert!(stdout(&o).contains("Input and output files are unchanged."));
}

#[test]
fn test_extension_notation_is_refused() {
    let dir = fixture(SHORT_IMG);
    let o = run(FPACK, dir.path(), &["a.fits[1]"]);
    assert_eq!(status(&o), EXIT_FAILURE);
    assert!(stdout(&o).contains("section/extension notation not supported"));
}

#[test]
fn test_fpack_refuses_an_input_that_is_already_packed() {
    let dir = fixture(SHORT_IMG);
    std::fs::rename(dir.path().join("a.fits"), dir.path().join("a.fits.fz")).unwrap();

    let o = run(FPACK, dir.path(), &["a.fits.fz"]);
    assert_eq!(status(&o), EXIT_FAILURE);
    assert!(stdout(&o).contains("already has '.fz' suffix"));
}

#[test]
fn test_fpack_refuses_to_overwrite_an_existing_output() {
    let dir = fixture(SHORT_IMG);
    std::fs::write(dir.path().join("a.fits.fz"), b"in the way").unwrap();

    let o = run(FPACK, dir.path(), &["a.fits"]);
    assert_eq!(status(&o), EXIT_FAILURE);
    assert!(stdout(&o).contains("Error: output file already exists"));
}

/// Regression for upstream bug 1: the C stores `ntile[6]` -- past the end of
/// the array inside `fpstate` -- before deciding there are too many
/// dimensions.  Seven dimensions must be refused, not corrupt the struct.
#[test]
fn test_seventh_tile_dimension_is_refused() {
    let dir = fixture(SHORT_IMG);

    let o = run(FPACK, dir.path(), &["-t", "1,2,3,4,5,6", "a.fits"]);
    assert_eq!(status(&o), 0, "six dimensions is the most that fits");

    let dir = fixture(SHORT_IMG);
    let o = run(FPACK, dir.path(), &["-t", "1,2,3,4,5,6,7", "a.fits"]);
    assert_eq!(status(&o), EXIT_FAILURE);
    assert!(stdout(&o).contains("Error: too many dimensions for `-t', max=6"));
    assert!(!dir.path().join("a.fits.fz").exists());
}

/* ---------------------------------------------------------------------- */
/* round trips                                                             */
/* ---------------------------------------------------------------------- */

/// fpack then funpack, and the bytes must come back unchanged.  `-C` drops the
/// timestamped checksum cards; `flags` must select a lossless combination.
fn assert_lossless_round_trip(bitpix: c_int, flags: &[&str]) {
    let dir = fixture(bitpix);
    let original = std::fs::read(dir.path().join("a.fits")).unwrap();

    let mut pack: Vec<&str> = vec!["-C"];
    pack.extend_from_slice(flags);
    pack.push("a.fits");

    let o = run(FPACK, dir.path(), &pack);
    assert_eq!(
        status(&o),
        0,
        "fpack {flags:?} BITPIX {bitpix}: {}",
        stdout(&o)
    );
    assert!(
        dir.path().join("a.fits.fz").exists(),
        "fpack {flags:?} produced no output"
    );
    assert!(
        dir.path().join("a.fits").exists(),
        "fpack should leave the input alone without -F or -D"
    );

    /* funpack refuses to overwrite an existing output */
    std::fs::remove_file(dir.path().join("a.fits")).unwrap();

    let o = run(FUNPACK, dir.path(), &["-C", "a.fits.fz"]);
    assert_eq!(
        status(&o),
        0,
        "funpack {flags:?} BITPIX {bitpix}: {}",
        stdout(&o)
    );

    let back = std::fs::read(dir.path().join("a.fits")).unwrap();
    assert_eq!(
        back, original,
        "round trip of BITPIX {bitpix} through {flags:?} was not lossless"
    );
}

macro_rules! round_trip_tests {
    ($($name:ident: $bitpix:expr, [$($flag:literal),*];)*) => { $(
        #[test]
        fn $name() { assert_lossless_round_trip($bitpix, &[$($flag),*]); }
    )* };
}

/* Integer images are lossless under every algorithm; `-q0 4` only fixes the
dithering so that successive runs agree. */
round_trip_tests! {
    test_round_trip_i16_rice:    SHORT_IMG,  ["-q0", "4", "-r"];
    test_round_trip_i16_gzip1:   SHORT_IMG,  ["-q0", "4", "-g"];
    test_round_trip_i16_gzip2:   SHORT_IMG,  ["-q0", "4", "-g2"];
    test_round_trip_i16_hcomp:   SHORT_IMG,  ["-q0", "4", "-h"];
    test_round_trip_i16_none:    SHORT_IMG,  ["-q0", "4", "-d"];
    test_round_trip_i16_default: SHORT_IMG,  ["-q0", "4"];
    test_round_trip_i16_whole:   SHORT_IMG,  ["-q0", "4", "-w"];
    test_round_trip_i16_tiled:   SHORT_IMG,  ["-q0", "4", "-t", "20,10"];
    test_round_trip_u16_plio:    USHORT_IMG, ["-q0", "4", "-p"];
    test_round_trip_u16_rice:    USHORT_IMG, ["-q0", "4", "-r"];
    test_round_trip_u16_hcomp:   USHORT_IMG, ["-q0", "4", "-h"];
    test_round_trip_i32_rice:    LONG_IMG,   ["-q0", "4", "-r"];
    test_round_trip_i32_gzip2:   LONG_IMG,   ["-q0", "4", "-g2"];
    /* floats only survive intact with quantization switched off entirely */
    test_round_trip_f32_gzip1:   FLOAT_IMG,  ["-g", "-q", "0"];
    test_round_trip_f32_gzip2:   FLOAT_IMG,  ["-g2", "-q", "0"];
    test_round_trip_f64_gzip1:   DOUBLE_IMG, ["-g", "-q", "0"];
    test_round_trip_f64_gzip2:   DOUBLE_IMG, ["-g2", "-q", "0"];
}

/// PLIO only handles positive 8- and 16-bit values, so a signed image with
/// negative pixels must fail -- in both implementations alike.
#[test]
fn test_plio_refuses_negative_pixels() {
    let dir = fixture(SHORT_IMG);
    let o = run(FPACK, dir.path(), &["-C", "-p", "a.fits"]);
    assert_ne!(status(&o), 0, "PLIO should refuse negative pixel values");
}

/// Quantizing a float image is lossy on purpose: `-q 4` is a quarter of the
/// *estimated background noise*, so the error is small only for data with a
/// smooth background -- hence the gradient rather than the noise the lossless
/// fixtures use.
#[test]
fn test_quantized_float_round_trip_is_close() {
    use bytemuck::cast_slice;
    use std::ffi::CString;

    let dir = tempfile::Builder::new()
        .prefix("fpack-q-")
        .tempdir()
        .unwrap();
    let path = dir.path().join("a.fits");

    let smooth: Vec<f32> = (0..NELEM)
        .map(|i| {
            let (x, y) = ((i as LONGLONG % NX) as f32, (i as LONGLONG / NX) as f32);
            1000.0 + x + 0.5 * y + ((x * 0.7).sin() + (y * 0.3).cos()) * 0.5
        })
        .collect();

    {
        let cpath = CString::new(path.to_str().unwrap()).unwrap();
        let name = cast_slice(cpath.to_bytes_with_nul());
        let mut status: c_int = 0;
        let mut fptr: Option<Box<fitsfile>> = None;
        fits_create_diskfile(&mut fptr, name, &mut status);
        let f = fptr.as_mut().unwrap();
        fits_create_imgll(f, FLOAT_IMG, 2, &[NX, NY], &mut status);
        fits_write_img_flt(f, 1, 1, NELEM as LONGLONG, &smooth, &mut status);
        fits_close_file(fptr.take().unwrap(), &mut status);
        assert_eq!(status, 0, "could not build the smooth fixture");
    }

    let o = run(FPACK, dir.path(), &["-C", "-q0", "4", "a.fits"]);
    assert_eq!(status(&o), 0, "{}", stdout(&o));
    std::fs::remove_file(&path).unwrap();

    let o = run(FUNPACK, dir.path(), &["-C", "a.fits.fz"]);
    assert_eq!(status(&o), 0, "{}", stdout(&o));

    let back = read_float_pixels(path.to_str().unwrap());
    assert_eq!(back.len(), smooth.len());

    let worst = smooth
        .iter()
        .zip(&back)
        .map(|(a, b)| (a - b).abs())
        .fold(0.0f32, f32::max);
    assert!(
        worst < 1.0,
        "quantized round trip of a smooth image drifted by {worst}"
    );
    assert_ne!(back, smooth, "q = 4 is meant to be lossy");
}

/// Read a 2-D float image back as pixels.
fn read_float_pixels(path: &str) -> Vec<f32> {
    use bytemuck::cast_slice;
    use rsfitsio::fitsio::READONLY;
    use std::ffi::CString;

    let cpath = CString::new(path).unwrap();
    let name = cast_slice(cpath.to_bytes_with_nul());
    let mut status: c_int = 0;
    let mut fptr: Option<Box<fitsfile>> = None;
    fits_open_diskfile(&mut fptr, name, READONLY, &mut status);
    let f = fptr.as_mut().unwrap();
    let mut out = vec![0f32; NELEM];
    fits_read_img_flt(f, 1, 1, NELEM as LONGLONG, 0.0, &mut out, None, &mut status);
    fits_close_file(fptr.take().unwrap(), &mut status);
    assert_eq!(status, 0, "could not read {path} back");
    out
}

/* ---------------------------------------------------------------------- */
/* the file-handling options                                               */
/* ---------------------------------------------------------------------- */

#[test]
fn test_listing_leaves_the_file_alone() {
    let dir = fixture(SHORT_IMG);
    run(FPACK, dir.path(), &["-C", "-q0", "4", "a.fits"]);
    let before = std::fs::read(dir.path().join("a.fits.fz")).unwrap();

    let o = run(FPACK, dir.path(), &["-L", "a.fits.fz"]);
    assert_eq!(status(&o), 0, "{}", stdout(&o));

    let out = stdout(&o);
    assert!(out.starts_with("# a.fits.fz ("), "{out}");
    assert!(out.contains(" bytes)"), "{out}");
    assert!(out.contains("  1 IMAGE"), "{out}");
    assert!(out.contains("tiled_rice"), "{out}");

    let after = std::fs::read(dir.path().join("a.fits.fz")).unwrap();
    assert_eq!(before, after, "-L must not modify the file");
}

#[test]
fn test_verbose_names_the_files() {
    let dir = fixture(SHORT_IMG);
    let o = run(FPACK, dir.path(), &["-C", "-q0", "4", "-v", "a.fits"]);
    assert_eq!(status(&o), 0);
    assert_eq!(stdout(&o), "a.fits -> a.fits.fz\n");
}

#[test]
fn test_clobber_replaces_the_input() {
    let dir = fixture(SHORT_IMG);
    let original = std::fs::read(dir.path().join("a.fits")).unwrap();

    let o = run(FPACK, dir.path(), &["-C", "-q0", "4", "-F", "a.fits"]);
    assert_eq!(status(&o), 0, "{}", stdout(&o));
    assert!(!dir.path().join("a.fits.fz").exists(), "-F writes in place");

    let packed = std::fs::read(dir.path().join("a.fits")).unwrap();
    assert_ne!(packed, original, "a.fits should now be the compressed file");

    /* and funpack -F puts it back */
    let o = run(FUNPACK, dir.path(), &["-C", "-F", "a.fits"]);
    assert_eq!(status(&o), 0, "{}", stdout(&o));
    assert_eq!(std::fs::read(dir.path().join("a.fits")).unwrap(), original);
}

#[test]
fn test_delete_input_removes_it() {
    let dir = fixture(SHORT_IMG);
    let o = run(FPACK, dir.path(), &["-C", "-q0", "4", "-D", "-Y", "a.fits"]);
    assert_eq!(status(&o), 0, "{}", stdout(&o));
    assert!(!dir.path().join("a.fits").exists(), "-D should delete it");
    assert!(dir.path().join("a.fits.fz").exists());
}

#[test]
fn test_named_output_file() {
    let dir = fixture(SHORT_IMG);
    let o = run(
        FPACK,
        dir.path(),
        &["-C", "-q0", "4", "-O", "packed.fz", "a.fits"],
    );
    assert_eq!(status(&o), 0, "{}", stdout(&o));
    assert!(dir.path().join("packed.fz").exists());
    assert!(!dir.path().join("a.fits.fz").exists());

    let o = run(FUNPACK, dir.path(), &["-C", "-O", "back.fits", "packed.fz"]);
    assert_eq!(status(&o), 0, "{}", stdout(&o));
    assert_eq!(
        std::fs::read(dir.path().join("back.fits")).unwrap(),
        std::fs::read(dir.path().join("a.fits")).unwrap()
    );
}

#[test]
fn test_funpack_prefix() {
    let dir = fixture(SHORT_IMG);
    run(FPACK, dir.path(), &["-C", "-q0", "4", "a.fits"]);
    std::fs::remove_file(dir.path().join("a.fits")).unwrap();

    let o = run(FUNPACK, dir.path(), &["-C", "-P", "new_", "a.fits.fz"]);
    assert_eq!(status(&o), 0, "{}", stdout(&o));
    assert!(
        dir.path().join("new_a.fits").exists(),
        "-P should prepend the prefix"
    );
}

#[test]
fn test_funpack_requires_the_fz_suffix() {
    let dir = fixture(SHORT_IMG);
    run(FPACK, dir.path(), &["-C", "-q0", "4", "a.fits"]);
    std::fs::remove_file(dir.path().join("a.fits")).unwrap();
    std::fs::rename(dir.path().join("a.fits.fz"), dir.path().join("packed.dat")).unwrap();

    let o = run(FUNPACK, dir.path(), &["packed.dat"]);
    assert_eq!(status(&o), EXIT_FAILURE);
    assert!(stdout(&o).contains("does not have the default .fz suffix"));
}

#[test]
fn test_funpack_ambiguous_input() {
    let dir = fixture(SHORT_IMG);
    run(FPACK, dir.path(), &["-C", "-q0", "4", "a.fits"]);
    /* both a.fits and a.fits.fz now exist */
    let o = run(FUNPACK, dir.path(), &["a.fits"]);
    assert_eq!(status(&o), EXIT_FAILURE);
    assert!(stdout(&o).contains("Error: ambiguous input file name"));
}

/// `funpack -Z` gzips its output; the C shells out to `gzip -1`.
#[test]
fn test_funpack_gzips_its_output() {
    let dir = fixture(SHORT_IMG);
    run(FPACK, dir.path(), &["-C", "-q0", "4", "a.fits"]);
    let original = std::fs::read(dir.path().join("a.fits")).unwrap();
    std::fs::remove_file(dir.path().join("a.fits")).unwrap();

    let o = run(FUNPACK, dir.path(), &["-C", "-Z", "a.fits.fz"]);
    assert_eq!(status(&o), 0, "{}", stdout(&o));
    assert!(
        !dir.path().join("a.fits").exists(),
        "-Z leaves only the .gz"
    );

    let gz = std::fs::read(dir.path().join("a.fits.gz")).unwrap();
    assert_eq!(&gz[..2], &[0x1f, 0x8b], "not a gzip stream");

    let mut back = Vec::new();
    std::io::Read::read_to_end(&mut flate2::read::GzDecoder::new(&gz[..]), &mut back).unwrap();
    assert_eq!(
        back, original,
        "the gzipped output is not the unpacked file"
    );
}

/// `-S` writes the FITS stream to stdout.
///
/// The stream must be byte-identical to the file `fpack` writes without `-S`.
/// That is what catches a stdout path that mangles the bytes: writing through
/// the C runtime on windows, where its stdout is in text mode, turns every
/// 0x0A into 0x0D 0x0A and silently corrupts the file.
#[test]
fn test_stdout_stream_matches_the_file() {
    let dir = fixture(SHORT_IMG);
    let original = std::fs::read(dir.path().join("a.fits")).unwrap();

    /* the same input packed to a file, for comparison */
    std::fs::copy(dir.path().join("a.fits"), dir.path().join("b.fits")).unwrap();
    let o = run(FPACK, dir.path(), &["-C", "-q0", "4", "b.fits"]);
    assert_eq!(status(&o), 0, "{}", stdout(&o));
    let on_disk = std::fs::read(dir.path().join("b.fits.fz")).unwrap();

    let o = run(FPACK, dir.path(), &["-C", "-q0", "4", "-S", "a.fits"]);
    assert_eq!(status(&o), 0, "{}", String::from_utf8_lossy(&o.stderr));
    assert!(
        !dir.path().join("a.fits.fz").exists(),
        "-S must not also write a file"
    );
    assert_eq!(
        o.stdout.len(),
        on_disk.len(),
        "the -S stream is {} bytes where the file is {}",
        o.stdout.len(),
        on_disk.len()
    );
    assert_eq!(o.stdout, on_disk, "the -S stream differs from the file");

    /* and it is readable: feed it back through funpack */
    std::fs::write(dir.path().join("c.fits.fz"), &o.stdout).unwrap();
    let o = run(FUNPACK, dir.path(), &["-C", "c.fits.fz"]);
    assert_eq!(status(&o), 0, "{}", stdout(&o));
    assert_eq!(std::fs::read(dir.path().join("c.fits")).unwrap(), original);
}

/// `-E` unpacks only the named or numbered HDUs.
#[test]
fn test_funpack_extension_list() {
    let dir = fixture(SHORT_IMG);
    run(FPACK, dir.path(), &["-C", "-q0", "4", "a.fits"]);
    std::fs::remove_file(dir.path().join("a.fits")).unwrap();

    let o = run(FUNPACK, dir.path(), &["-C", "-E", "1", "a.fits.fz"]);
    assert_eq!(status(&o), 0, "{}", stdout(&o));
    assert!(dir.path().join("a.fits").exists());

    /* a name that is not in the file is reported, not ignored */
    let dir2 = fixture(SHORT_IMG);
    run(FPACK, dir2.path(), &["-C", "-q0", "4", "a.fits"]);
    std::fs::remove_file(dir2.path().join("a.fits")).unwrap();
    let o = run(
        FUNPACK,
        dir2.path(),
        &["-C", "-E", "NOSUCHHDU", "a.fits.fz"],
    );
    assert_ne!(status(&o), 0);
    assert!(stdout(&o).contains("Unable to find and move to extension 'NOSUCHHDU'"));
}

/* ---------------------------------------------------------------------- */
/* -T, the algorithm comparison report                                     */
/* ---------------------------------------------------------------------- */

/// `-T` benchmarks every algorithm and leaves the input untouched.  The
/// numbers are timings, so only the structure is asserted here.
#[test]
fn test_report_leaves_files_unchanged() {
    let dir = fixture(SHORT_IMG);
    let before = std::fs::read(dir.path().join("a.fits")).unwrap();

    let o = run(FPACK, dir.path(), &["-T", "a.fits"]);
    assert_eq!(status(&o), 0, "{}", String::from_utf8_lossy(&o.stderr));

    let out = stdout(&o);
    assert!(out.contains(" File: a.fits"), "{out}");
    assert!(out.contains("Ext BITPIX Dimens."), "{out}");
    assert!(out.contains("Type   Ratio"), "{out}");
    assert!(out.contains("       Native"), "{out}");
    assert!(out.contains("RICE"), "{out}");
    assert!(out.contains("GZIP1"), "{out}");
    assert!(out.contains("HCOMP"), "{out}");

    assert_eq!(
        std::fs::read(dir.path().join("a.fits")).unwrap(),
        before,
        "-T must not modify the input"
    );
    assert!(
        !dir.path().join("a.fits.fz").exists(),
        "-T must not write an output file"
    );
    /* and the two scratch files it makes in the CWD are gone again */
    let leftovers: Vec<_> = std::fs::read_dir(dir.path())
        .unwrap()
        .map(|e| e.unwrap().file_name().to_string_lossy().into_owned())
        .filter(|n| n != "a.fits")
        .collect();
    assert!(
        leftovers.is_empty(),
        "temporary files left behind: {leftovers:?}"
    );
}

/// `-R` writes the same report to a text file, and only makes sense with -T.
#[test]
fn test_report_file() {
    let dir = fixture(SHORT_IMG);

    let o = run(FPACK, dir.path(), &["-R", "report.txt", "a.fits"]);
    assert_eq!(status(&o), EXIT_FAILURE, "-R alone is refused");
    assert!(stdout(&o).contains("Error: -R option may only be used with -T"));

    let o = run(FPACK, dir.path(), &["-T", "-R", "report.txt", "a.fits"]);
    assert_eq!(status(&o), 0, "{}", String::from_utf8_lossy(&o.stderr));
    let report = std::fs::read_to_string(dir.path().join("report.txt")).unwrap();
    assert!(report.starts_with(" Filename Extension BITPIX"), "{report}");
    assert!(report.contains("a.fits"), "{report}");
}

/* ---------------------------------------------------------------------- */
/* interoperability with the C implementation                              */
/* ---------------------------------------------------------------------- */

/// A file packed here must be readable by the C `funpack`, and vice versa, for
/// every algorithm.  PLIO was the one exception until pliocomp 0.5.0 fixed its
/// `LL_LEN` off-by-one, so the empty `expected` list is what keeps that fixed.
/// Needs `C_FPACK`/`C_FUNPACK`; run with `--ignored`.
#[test]
#[ignore = "needs the C fpack/funpack; set C_FPACK and C_FUNPACK"]
fn test_interoperates_with_the_c_implementation() {
    let (c_fpack, c_funpack) = match (std::env::var("C_FPACK"), std::env::var("C_FUNPACK")) {
        (Ok(a), Ok(b)) => (a, b),
        _ => panic!("set C_FPACK and C_FUNPACK to the C binaries"),
    };

    let mut failures = Vec::new();
    for alg in ["-r", "-g", "-g2", "-h", "-d", "-p"] {
        for (packer, unpacker, label) in [
            (FPACK, c_funpack.as_str(), "rust -> C"),
            (c_fpack.as_str(), FUNPACK, "C -> rust"),
        ] {
            /* PLIO needs non-negative pixels */
            let dir = if alg == "-p" {
                nonneg_fixture()
            } else {
                fixture(SHORT_IMG)
            };
            let original = std::fs::read(dir.path().join("a.fits")).unwrap();

            let o = run(packer, dir.path(), &["-C", "-q0", "4", alg, "a.fits"]);
            if status(&o) != 0 {
                failures.push(format!("{alg} {label}: pack failed"));
                continue;
            }
            std::fs::remove_file(dir.path().join("a.fits")).unwrap();

            let o = run(unpacker, dir.path(), &["-C", "a.fits.fz"]);
            if status(&o) != 0 {
                failures.push(format!("{alg} {label}: unpack failed"));
                continue;
            }
            if std::fs::read(dir.path().join("a.fits")).unwrap() != original {
                failures.push(format!("{alg} {label}: data changed"));
            }
        }
    }

    let expected: [&str; 0] = [];
    let unexpected: Vec<_> = failures
        .iter()
        .filter(|f| !expected.contains(&f.as_str()))
        .collect();
    assert!(unexpected.is_empty(), "{unexpected:#?}");
}

/// A `USHORT_IMG` (BZERO = 32768) file compressed by the *C* fpack cannot be
/// uncompressed by this port: every algorithm fails alike with
/// `NUM_OVERFLOW` (412), "Numerical overflow during type conversion while
/// writing FITS data", so the fault is in the unsigned-16-bit scaling on the
/// decompression path and not in any codec.
///
/// Rust -> Rust round trips fine (see `test_round_trip_u16_rice`), which is
/// why this only shows up against the C.  Ignored by default because it needs
/// the C binaries; run with `C_FPACK` set and `--ignored`.
#[test]
#[ignore = "known failure: USHORT decompression overflow; needs C_FPACK"]
fn test_c_packed_ushort_is_unreadable() {
    let c_fpack = std::env::var("C_FPACK").expect("set C_FPACK to the C binary");

    for alg in ["-r", "-g", "-p", "-d"] {
        let dir = fixture(USHORT_IMG);
        let original = std::fs::read(dir.path().join("a.fits")).unwrap();

        let o = run(&c_fpack, dir.path(), &["-C", "-q0", "4", alg, "a.fits"]);
        assert_eq!(status(&o), 0, "the C should pack {alg} happily");
        std::fs::remove_file(dir.path().join("a.fits")).unwrap();

        let o = run(FUNPACK, dir.path(), &["-C", "a.fits.fz"]);
        assert_eq!(status(&o), 0, "{alg}: funpack failed: {}", stdout(&o));
        assert_eq!(
            std::fs::read(dir.path().join("a.fits")).unwrap(),
            original,
            "{alg}: USHORT round trip through the C changed the data"
        );
    }
}

/// Compressed binary tables, over the column types the transposer branches on.
#[test]
fn test_compressed_table_round_trip() {
    for flag in ["-tableonly", "-table"] {
        let dir = tempfile::Builder::new()
            .prefix("fpack-tbl-")
            .tempdir()
            .unwrap();
        let path = dir.path().join("a.fits");
        make_binary_table(path.to_str().unwrap());
        let original = std::fs::read(&path).unwrap();

        let o = run(FPACK, dir.path(), &["-C", flag, "a.fits"]);
        assert_eq!(status(&o), 0, "fpack {flag}: {}", stdout(&o));
        assert!(
            std::fs::metadata(dir.path().join("a.fits.fz"))
                .unwrap()
                .len()
                < std::fs::metadata(&path).unwrap().len(),
            "{flag} did not actually compress the table"
        );
        std::fs::remove_file(&path).unwrap();

        let o = run(FUNPACK, dir.path(), &["-C", "a.fits.fz"]);
        assert_eq!(status(&o), 0, "funpack after {flag}: {}", stdout(&o));
        assert_eq!(
            std::fs::read(&path).unwrap(),
            original,
            "compressed-table round trip through {flag} changed the data"
        );
    }
}

/// A variable-length string column too large for the transposed buffer must be
/// refused, not written past the end of it.  The C aborts with a heap overflow
/// here (heasarc/cfitsio#134).
#[test]
fn test_oversized_vla_column_is_refused_not_overrun() {
    let c_fpack = match std::env::var("C_FPACK") {
        Ok(p) => p,
        Err(_) => return, // needs a file the C can produce; skipped without it
    };
    let dir = tempfile::Builder::new()
        .prefix("fpack-vla-")
        .tempdir()
        .unwrap();
    let src = Path::new("testprog.fit");
    if !src.exists() {
        return;
    }
    std::fs::copy(src, dir.path().join("a.fits")).unwrap();

    let o = run(
        &c_fpack,
        dir.path(),
        &["-C", "-q0", "4", "-tableonly", "a.fits"],
    );
    assert_eq!(status(&o), 0, "the C should pack testprog.fit");
    std::fs::remove_file(dir.path().join("a.fits")).unwrap();

    let o = run(FUNPACK, dir.path(), &["-C", "a.fits.fz"]);
    /* a normal exit, not a signal -- the C dies here with SIGABRT */
    assert!(
        o.status.code().is_some(),
        "funpack was killed by a signal: {:?}",
        o.status
    );
    assert_ne!(status(&o), 0, "an oversized VLA column must be refused");
    /* the CFITSIO error stack goes to stderr */
    let stderr = String::from_utf8_lossy(&o.stderr);
    assert!(
        stderr.contains("cfitsio#134"),
        "expected the labelled limitation, got stdout [{}] stderr [{stderr}]",
        stdout(&o)
    );
}

/// A binary table over one FITS block, so fpack compresses rather than copies.
fn make_binary_table(path: &str) {
    use bytemuck::cast_slice;
    use rsfitsio::c_types::c_char;
    use rsfitsio::fitsio::{BINARY_TBL, TDOUBLE, TFLOAT, TINT, TSHORT};
    use std::ffi::CString;

    let cpath = CString::new(path).unwrap();
    let name = cast_slice(cpath.to_bytes_with_nul());
    let mut status: c_int = 0;
    let mut fptr: Option<Box<fitsfile>> = None;

    /* one column per branch of the transpose */
    let nrows: LONGLONG = 400;
    let ttype: [&[c_char]; 5] = [
        cast_slice(c"IVAL".to_bytes_with_nul()),
        cast_slice(c"DVAL".to_bytes_with_nul()),
        cast_slice(c"SVAL".to_bytes_with_nul()),
        cast_slice(c"EVAL".to_bytes_with_nul()),
        cast_slice(c"BVAL".to_bytes_with_nul()),
    ];
    let tform: [&[c_char]; 5] = [
        cast_slice(c"1J".to_bytes_with_nul()),
        cast_slice(c"1D".to_bytes_with_nul()),
        cast_slice(c"16A".to_bytes_with_nul()),
        cast_slice(c"1E".to_bytes_with_nul()),
        cast_slice(c"1I".to_bytes_with_nul()),
    ];

    fits_create_diskfile(&mut fptr, name, &mut status);
    let f = fptr.as_mut().unwrap();
    fits_create_tbl(
        f,
        BINARY_TBL,
        0,
        5,
        &ttype.map(Some),
        &tform,
        None,
        None,
        &mut status,
    );
    let n = nrows as usize;
    let ivals: Vec<c_int> = (0..n).map(|i| i as c_int * 7919).collect();
    let dvals: Vec<f64> = (0..n).map(|i| i as f64 * 1.5).collect();
    let evals: Vec<f32> = (0..n).map(|i| i as f32 * 0.25).collect();
    let bvals: Vec<i16> = (0..n).map(|i| (i % 30000) as i16).collect();
    fits_write_col(f, TINT, 1, 1, 1, nrows, cast_slice(&ivals), &mut status);
    fits_write_col(f, TDOUBLE, 2, 1, 1, nrows, cast_slice(&dvals), &mut status);
    for (row, sval) in (0..n).enumerate() {
        let text = CString::new(format!("name{sval:011}")).unwrap();
        let cell: [&[c_char]; 1] = [cast_slice(text.to_bytes_with_nul())];
        fits_write_col_str(f, 3, row as LONGLONG + 1, 1, 1, &cell, &mut status);
    }
    fits_write_col(f, TFLOAT, 4, 1, 1, nrows, cast_slice(&evals), &mut status);
    fits_write_col(f, TSHORT, 5, 1, 1, nrows, cast_slice(&bvals), &mut status);
    fits_close_file(fptr.take().unwrap(), &mut status);
    assert_eq!(status, 0, "could not build the binary-table fixture");
}
