/* H-compress routines */

use bytemuck::cast_slice_mut;
use hcompress::write::HCEncoder;

use crate::c_types::{c_char, c_int, c_long};

use crate::fitsio::{DATA_COMPRESSION_ERR, LONGLONG, NULL_MSG};

/* ---------------------------------------------------------------------- */
/// Compress the input image using the H-compress algorithm
///
/// a  - input image array
/// nx - size of X axis of image
/// ny - size of Y axis of image
/// scale - quantization scale factor. Larger values results in more (lossy) compression
/// scale = 0 does lossless compression
/// output - pre-allocated array to hold the output compressed stream of bytes
/// nbytes  - input value = size of the output buffer;
/// returned value = size of the compressed byte stream, in bytes
///
/// NOTE: the nx and ny dimensions as defined within this code are reversed from
/// the usual FITS notation.  ny is the fastest varying dimension, which is
/// usually considered the X axis in the FITS image display
///
/// NOTE: `a` is modified in place (the H-transform and digitize steps overwrite
/// it), which is why it is a `*mut` here and a `&mut [c_int]` in the safe form.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn fits_hcompress(
    a: *mut c_int,
    ny: c_int,
    nx: c_int,
    scale: c_int,
    output: *mut c_char,
    nbytes: *mut c_long,
    status: *mut c_int,
) -> c_int {
    // FFI WRAPPER
    unsafe {
        let nbytes = nbytes.as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);

        let a = core::slice::from_raw_parts_mut(a, nx as usize * ny as usize);
        let output = core::slice::from_raw_parts_mut(output, *nbytes as usize);

        fits_hcompress_safe(a, ny, nx, scale, output, nbytes, status)
    }
}

/* ---------------------------------------------------------------------- */
/// Compress the input image using the H-compress algorithm
///
/// a  - input image array
/// nx - size of X axis of image
/// ny - size of Y axis of image
/// scale - quantization scale factor. Larger values results in more (lossy) compression
/// scale = 0 does lossless compression
/// output - pre-allocated array to hold the output compressed stream of bytes
/// nbytes  - input value = size of the output buffer;
/// returned value = size of the compressed byte stream, in bytes
///
/// NOTE: the nx and ny dimensions as defined within this code are reversed from
/// the usual FITS notation.  ny is the fastest varying dimension, which is
/// usually considered the X axis in the FITS image display
pub fn fits_hcompress_safe(
    a: &mut [c_int],
    ny: c_int,
    nx: c_int,
    scale: c_int,
    output: &mut [c_char],
    nbytes: &mut c_long,
    status: &mut c_int,
) -> c_int {
    if *status > 0 {
        return *status;
    }

    /* The C's htrans() / digitize() / encode() are not in this crate -- they
    live in the external `hcompress` crate, where HCEncoder::write() is exactly
    that same sequence.  The C's FFLOCK/FFUNLOCK around encode() guarded the
    file-scope `noutmax`/bit-buffer globals; the crate keeps that state inside
    the HCEncoder, so no lock is needed here. */

    /* noutmax = *nbytes;  input value is the allocated size of the array */
    let noutmax = (*nbytes as usize).min(output.len());
    *nbytes = 0; /* reset */

    /* DEVIATION: the C's encode() compares its running byte count against
    `noutmax` and fails with "encode: output buffer too small".  The crate's
    encoder discards its sink's io::Errors (`let _ = write_all(..)`), so writing
    straight into `output` would truncate silently.  Encode into a Vec instead
    and enforce `noutmax` here, so the caller still sees DATA_COMPRESSION_ERR. */
    let mut buf: Vec<u8> = Vec::with_capacity(noutmax);

    /* H-transform, digitize, then encode and write to the output array */
    let stat = HCEncoder::new(&mut buf).write(a, ny as usize, nx as usize, scale);

    if stat.is_err() || buf.len() > noutmax {
        *status = DATA_COMPRESSION_ERR;
        return *status;
    }

    cast_slice_mut(&mut output[..buf.len()]).copy_from_slice(&buf);
    *nbytes = buf.len() as c_long;

    *status
}

/* ---------------------------------------------------------------------- */
/// Compress the input image using the H-compress algorithm
///   
/// a  - input image array
/// nx - size of X axis of image
/// ny - size of Y axis of image
/// scale - quantization scale factor. Larger values results in more (lossy) compression
///         scale = 0 does lossless compression
/// output - pre-allocated array to hold the output compressed stream of bytes
/// nbyts  - size of the compressed byte stream, in bytes
///
/// NOTE: the nx and ny dimensions as defined within this code are reversed from
/// the usual FITS notation.  ny is the fastest varying dimension, which is
/// usually considered the X axis in the FITS image display
///
/// NOTE: `a` is modified in place, as in the 32-bit routine above.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn fits_hcompress64(
    a: *mut LONGLONG,
    ny: c_int,
    nx: c_int,
    scale: c_int,
    output: *mut c_char,
    nbytes: *mut c_long,
    status: *mut c_int,
) -> c_int {
    // FFI WRAPPER
    unsafe {
        let nbytes = nbytes.as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);

        let a = core::slice::from_raw_parts_mut(a, nx as usize * ny as usize);
        let output = core::slice::from_raw_parts_mut(output, *nbytes as usize);

        fits_hcompress64_safe(a, ny, nx, scale, output, nbytes, status)
    }
}

/* ---------------------------------------------------------------------- */
/// Compress the input image using the H-compress algorithm
///   
/// a  - input image array
/// nx - size of X axis of image
/// ny - size of Y axis of image
/// scale - quantization scale factor. Larger values results in more (lossy) compression
///         scale = 0 does lossless compression
/// output - pre-allocated array to hold the output compressed stream of bytes
/// nbyts  - size of the compressed byte stream, in bytes
///
/// NOTE: the nx and ny dimensions as defined within this code are reversed from
/// the usual FITS notation.  ny is the fastest varying dimension, which is
/// usually considered the X axis in the FITS image display
pub fn fits_hcompress64_safe(
    a: &mut [LONGLONG],
    ny: c_int,
    nx: c_int,
    scale: c_int,
    output: &mut [c_char],
    nbytes: &mut c_long,
    status: &mut c_int,
) -> c_int {
    if *status > 0 {
        return *status;
    }

    /* As above: HCEncoder::write64() is the C's htrans64() + digitize64() +
    encode64(), and it carries its own state so FFLOCK/FFUNLOCK is moot. */

    /* noutmax = *nbytes;  input value is the allocated size of the array */
    let noutmax = (*nbytes as usize).min(output.len());
    *nbytes = 0; /* reset */

    /* DEVIATION: see fits_hcompress_safe -- the crate cannot report a full
    sink, so the C's "output buffer too small" check is done here instead. */
    let mut buf: Vec<u8> = Vec::with_capacity(noutmax);

    /* H-transform, digitize, then encode and write to the output array */
    let stat = HCEncoder::new(&mut buf).write64(a, ny as usize, nx as usize, scale);

    if stat.is_err() || buf.len() > noutmax {
        *status = DATA_COMPRESSION_ERR;
        return *status;
    }

    cast_slice_mut(&mut output[..buf.len()]).copy_from_slice(&buf);
    *nbytes = buf.len() as c_long;

    *status
}

/// Tests ported from cfitsio's test_hcompress.c
///
/// rsfitsio delegates H-compress to the external `hcompress` crate
/// (`HCEncoder`/`HCDecoder`), which is what `imcompress` drives internally and
/// what `fits_hcompress_safe`/`fits_hdecompress_safe` are written against. Most
/// of these tests drive the crate API directly, mirroring the C test's
/// round-trips; `test_hcompress_safe_*` go through the `_safe` wrappers instead.
///
/// NOTE: like the C API, `HCEncoder::write` takes the dimensions in `(ny, nx)`
/// order and modifies its input array in place, and `HCDecoder::read` returns
/// `(nx, ny, scale)`.
#[cfg(test)]
mod tests {
    use hcompress::read::HCDecoder;
    use hcompress::write::HCEncoder;

    use super::{fits_hcompress_safe, fits_hcompress64_safe};
    use crate::c_types::{c_char, c_int, c_long};
    use crate::fits_hdecompress::{fits_hdecompress_safe, fits_hdecompress64_safe};
    use crate::fitsio::{DATA_COMPRESSION_ERR, LONGLONG};

    /// Compress `original` (32-bit) losslessly-or-lossily, then decompress.
    /// Returns the decompressed image and the `(nx, ny)` reported by the decoder.
    fn roundtrip32(
        original: &[i32],
        nx: usize,
        ny: usize,
        scale: i32,
        smooth: i32,
    ) -> (Vec<i32>, usize, usize) {
        let mut a = original.to_vec(); // write() modifies its input in place
        let mut compressed = Vec::new();
        HCEncoder::new(&mut compressed)
            .write(&mut a, ny, nx, scale)
            .expect("hcompress write failed");
        assert!(!compressed.is_empty(), "compressed stream is empty");

        let mut decompressed = vec![0i32; nx * ny];
        let (rnx, rny, _scale) = HCDecoder::new()
            .read(&compressed, smooth, &mut decompressed)
            .expect("hcompress read failed");
        (decompressed, rnx, rny)
    }

    /// 64-bit variant. Returns the `(nx, ny)` reported by the decoder.
    fn roundtrip64(
        original: &[LONGLONG],
        nx: usize,
        ny: usize,
        scale: i32,
        smooth: i32,
    ) -> (usize, usize) {
        let mut a = original.to_vec();
        let mut compressed = Vec::new();
        HCEncoder::new(&mut compressed)
            .write64(&mut a, ny, nx, scale)
            .expect("hcompress write64 failed");
        assert!(!compressed.is_empty());

        let mut decompressed = vec![0 as LONGLONG; nx * ny];
        let (rnx, rny, _scale) = HCDecoder::new()
            .read64(&compressed, smooth, &mut decompressed)
            .expect("hcompress read64 failed");
        (rnx, rny)
    }

    #[test]
    fn test_hcompress_lossless_roundtrip() {
        let nx = 16;
        let ny = 16;
        let original: Vec<i32> = (0..(nx * ny) as i32).map(|i| i * 10).collect();
        let (decompressed, rnx, rny) = roundtrip32(&original, nx, ny, 0, 0);
        assert_eq!(rnx, nx);
        assert_eq!(rny, ny);
        assert_eq!(decompressed, original);
    }

    #[test]
    fn test_hcompress_lossy() {
        let nx = 16;
        let ny = 16;
        let original: Vec<i32> = (0..(nx * ny) as i32).map(|i| i * 100).collect();
        let (decompressed, rnx, rny) = roundtrip32(&original, nx, ny, 4, 0);
        assert_eq!(rnx, nx);
        assert_eq!(rny, ny);

        // With scale=4 the maximum absolute error should be bounded.
        let max_diff = original
            .iter()
            .zip(&decompressed)
            .map(|(o, d)| (o - d).abs())
            .max()
            .unwrap();
        assert!(max_diff <= 100, "max_diff {max_diff} exceeded bound");
    }

    #[test]
    fn test_hcompress_smooth() {
        let nx = 16;
        let ny = 16;
        let original: Vec<i32> = (0..(nx * ny) as i32).map(|i| i * 50).collect();
        let (_decompressed, rnx, rny) = roundtrip32(&original, nx, ny, 2, 1);
        assert_eq!(rnx, nx);
        assert_eq!(rny, ny);
    }

    #[test]
    fn test_hcompress_uniform() {
        let nx = 16;
        let ny = 16;
        let original = vec![42i32; nx * ny];

        let mut a = original.clone();
        let mut compressed = Vec::new();
        HCEncoder::new(&mut compressed)
            .write(&mut a, ny, nx, 0)
            .unwrap();
        // Uniform data compresses very well.
        assert!(compressed.len() <= 100, "uniform image did not compress");

        let mut decompressed = vec![0i32; nx * ny];
        HCDecoder::new()
            .read(&compressed, 0, &mut decompressed)
            .unwrap();
        assert!(decompressed.iter().all(|&v| v == 42));
    }

    #[test]
    fn test_hcompress_large() {
        let nx = 64;
        let ny = 64;
        let original: Vec<i32> = (0..(nx * ny))
            .map(|i| (i % nx) as i32 + (i / nx) as i32 * 100)
            .collect();
        let (decompressed, _rnx, _rny) = roundtrip32(&original, nx, ny, 0, 0);
        assert_eq!(decompressed, original);
    }

    #[test]
    fn test_hcompress_odd_dimensions() {
        let nx = 17;
        let ny = 19;
        let original: Vec<i32> = (0..(nx * ny) as i32).map(|i| i % 1000).collect();
        let (decompressed, rnx, rny) = roundtrip32(&original, nx, ny, 0, 0);
        assert_eq!(rnx, nx);
        assert_eq!(rny, ny);
        assert_eq!(decompressed, original);
    }

    #[test]
    fn test_hcompress_odd_x() {
        let nx = 15;
        let ny = 16;
        let original: Vec<i32> = (0..(nx * ny) as i32).map(|i| (i * 7) % 500).collect();
        let (decompressed, rnx, rny) = roundtrip32(&original, nx, ny, 0, 0);
        assert_eq!(rnx, nx);
        assert_eq!(rny, ny);
        assert_eq!(decompressed, original);
    }

    #[test]
    fn test_hcompress_odd_y() {
        let nx = 16;
        let ny = 15;
        let original: Vec<i32> = (0..(nx * ny) as i32).map(|i| (i * 13) % 800).collect();
        let (decompressed, rnx, rny) = roundtrip32(&original, nx, ny, 0, 0);
        assert_eq!(rnx, nx);
        assert_eq!(rny, ny);
        assert_eq!(decompressed, original);
    }

    #[test]
    fn test_hcompress_non_power2() {
        let nx = 20;
        let ny = 24;
        let original: Vec<i32> = (0..(nx * ny) as i32).map(|i| i * 3).collect();
        let (decompressed, _rnx, _rny) = roundtrip32(&original, nx, ny, 0, 0);
        assert_eq!(decompressed, original);
    }

    #[test]
    fn test_hcompress_buffer_too_small() {
        // A fixed 50-byte sink is too small for the compressed stream; the
        // encoder must fail cleanly (return Err) rather than panic.
        let nx = 16;
        let ny = 16;
        let mut a = vec![100i32; nx * ny];
        let mut buf = [0u8; 50];
        let res = HCEncoder::new(&mut buf[..]).write(&mut a, ny, nx, 0);
        // Either outcome is acceptable; the point is that it does not crash.
        let _ = res;
    }

    #[test]
    fn test_hcompress_minimal() {
        let nx = 4;
        let ny = 4;
        let original: Vec<i32> = (0..(nx * ny) as i32).collect();
        let (decompressed, _rnx, _rny) = roundtrip32(&original, nx, ny, 0, 0);
        assert_eq!(decompressed, original);
    }

    #[test]
    fn test_hcompress_prime_dims() {
        let nx = 31;
        let ny = 37;
        let original: Vec<i32> = (0..(nx * ny) as i32).map(|i| (i * 17) % 1000).collect();
        let (decompressed, rnx, rny) = roundtrip32(&original, nx, ny, 0, 0);
        assert_eq!(rnx, nx);
        assert_eq!(rny, ny);
        assert_eq!(decompressed, original);
    }

    #[test]
    fn test_hcompress_odd_lossy_smooth() {
        let nx = 23;
        let ny = 29;
        let original: Vec<i32> = (0..(nx * ny) as i32).map(|i| (i * 11) % 500).collect();
        // Lossy compression followed by smoothing decompression: just succeed.
        let (_decompressed, _rnx, _rny) = roundtrip32(&original, nx, ny, 4, 1);
    }

    #[test]
    fn test_hcompress_high_scale() {
        let nx = 32;
        let ny = 32;
        let original: Vec<i32> = (0..(nx * ny) as i32).map(|i| (i * 1000) % 100000).collect();
        let (_decompressed, _rnx, _rny) = roundtrip32(&original, nx, ny, 16, 1);
    }

    #[test]
    fn test_hcompress64_odd() {
        let nx = 17;
        let ny = 19;
        let original: Vec<LONGLONG> = (0..(nx * ny) as i64).map(|i| i * 13 + 100).collect();
        let (rnx, rny) = roundtrip64(&original, nx, ny, 0, 0);
        assert_eq!(rnx, nx);
        assert_eq!(rny, ny);
    }

    #[test]
    fn test_hcompress64_smooth() {
        let nx = 32;
        let ny = 32;
        let original: Vec<LONGLONG> = (0..(nx * ny) as i64).map(|i| i * 100).collect();
        let (rnx, rny) = roundtrip64(&original, nx, ny, 4, 1);
        assert_eq!(rnx, nx);
        assert_eq!(rny, ny);
    }

    #[test]
    fn test_hcompress64_large_smooth() {
        let nx = 64;
        let ny = 64;
        let original: Vec<LONGLONG> = (0..(nx * ny))
            .map(|i| ((i % nx) * 100 + (i / nx) * 10) as i64)
            .collect();
        let (_rnx, _rny) = roundtrip64(&original, nx, ny, 8, 1);
    }

    /* ---- round-trips through the `_safe` wrappers themselves ---- */

    /// `fits_hcompress_safe` -> `fits_hdecompress_safe`, lossless.
    #[test]
    fn test_hcompress_safe_roundtrip() {
        let nx: c_int = 17;
        let ny: c_int = 19;
        let original: Vec<c_int> = (0..(nx * ny)).map(|i| (i * 7) % 500).collect();

        let mut a = original.clone(); // compression works in place
        let mut output = vec![0 as c_char; 4 * original.len() + 1024];
        let mut nbytes = output.len() as c_long;
        let mut status: c_int = 0;

        fits_hcompress_safe(&mut a, ny, nx, 0, &mut output, &mut nbytes, &mut status);
        assert_eq!(status, 0);
        assert!(nbytes > 0 && (nbytes as usize) <= output.len());

        let input: Vec<u8> = output[..nbytes as usize].iter().map(|&c| c as u8).collect();
        let mut decompressed = vec![0 as c_int; original.len()];
        let (mut rny, mut rnx, mut rscale) = (0, 0, 0);

        fits_hdecompress_safe(
            &input,
            0,
            &mut decompressed,
            original.len() as c_int,
            &mut rny,
            &mut rnx,
            &mut rscale,
            &mut status,
        );
        assert_eq!(status, 0);
        assert_eq!(rnx, nx);
        assert_eq!(rny, ny);
        assert_eq!(rscale, 0);
        assert_eq!(decompressed, original);
    }

    /// The 64-bit pair. `fits_hdecompress64_safe` packs the I*8 result back
    /// into an I*4 array in place, as the C does.
    #[test]
    fn test_hcompress64_safe_roundtrip() {
        let nx: c_int = 16;
        let ny: c_int = 16;
        let nval = (nx * ny) as usize;
        let original: Vec<LONGLONG> = (0..nval as LONGLONG).map(|i| i * 13 + 100).collect();

        let mut a = original.clone();
        let mut output = vec![0 as c_char; 8 * nval + 1024];
        let mut nbytes = output.len() as c_long;
        let mut status: c_int = 0;

        fits_hcompress64_safe(&mut a, ny, nx, 0, &mut output, &mut nbytes, &mut status);
        assert_eq!(status, 0);
        assert!(nbytes > 0);

        let input: Vec<u8> = output[..nbytes as usize].iter().map(|&c| c as u8).collect();
        let mut decompressed = vec![0 as LONGLONG; nval];
        let (mut rny, mut rnx, mut rscale) = (0, 0, 0);

        fits_hdecompress64_safe(
            &input,
            0,
            &mut decompressed,
            nval as c_int,
            &mut rny,
            &mut rnx,
            &mut rscale,
            &mut status,
        );
        assert_eq!(status, 0);
        assert_eq!(rnx, nx);
        assert_eq!(rny, ny);

        /* the values come back packed as I*4 in the first half of the array */
        let packed: &[c_int] = bytemuck::cast_slice(&decompressed);
        let got: Vec<LONGLONG> = packed[..nval].iter().map(|&v| v as LONGLONG).collect();
        assert_eq!(got, original);
    }

    /// A non-zero input status is passed straight through, untouched.
    #[test]
    fn test_hcompress_safe_status_passthrough() {
        let mut a = vec![0 as c_int; 16];
        let mut output = vec![0 as c_char; 256];
        let mut nbytes = output.len() as c_long;
        let mut status: c_int = 104;

        assert_eq!(
            fits_hcompress_safe(&mut a, 4, 4, 0, &mut output, &mut nbytes, &mut status),
            104
        );
        assert_eq!(nbytes, 256); /* not reset */

        let mut decompressed = vec![0 as c_int; 16];
        let (mut rny, mut rnx, mut rscale) = (0, 0, 0);
        assert_eq!(
            fits_hdecompress_safe(
                &[],
                0,
                &mut decompressed,
                16,
                &mut rny,
                &mut rnx,
                &mut rscale,
                &mut status
            ),
            104
        );
    }

    /// An output buffer too small for the compressed stream reports
    /// DATA_COMPRESSION_ERR, as the C's "encode: output buffer too small" does.
    #[test]
    fn test_hcompress_safe_buffer_too_small() {
        let nx: c_int = 16;
        let ny: c_int = 16;
        let mut a: Vec<c_int> = (0..(nx * ny)).map(|i| (i * 977) % 65536).collect();
        let mut output = vec![0 as c_char; 20];
        let mut nbytes = output.len() as c_long;
        let mut status: c_int = 0;

        fits_hcompress_safe(&mut a, ny, nx, 0, &mut output, &mut nbytes, &mut status);
        assert_eq!(status, DATA_COMPRESSION_ERR);
    }

    /// A truncated/garbage stream is rejected rather than panicking.
    #[test]
    fn test_hdecompress_safe_bad_stream() {
        let mut decompressed = vec![0 as c_int; 256];
        let (mut rny, mut rnx, mut rscale) = (0, 0, 0);
        let mut status: c_int = 0;

        fits_hdecompress_safe(
            &[0xde, 0xad, 0xbe, 0xef],
            0,
            &mut decompressed,
            256,
            &mut rny,
            &mut rnx,
            &mut rscale,
            &mut status,
        );
        assert_ne!(status, 0);
    }
}
