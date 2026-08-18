/* H-compress routines */

use crate::c_types::{c_char, c_int, c_long};

use crate::fitsio::{LONGLONG, NULL_MSG};

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
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn fits_hcompress(
    a: *const c_int,
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

        let a = core::slice::from_raw_parts(a, nx as usize * ny as usize);
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
    _a: &[c_int],
    _ny: c_int,
    _nx: c_int,
    _scale: c_int,
    _output: &mut [c_char],
    _nbytes: &mut c_long,
    _status: &mut c_int,
) -> c_int {
    todo!()
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
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn fits_hcompress64(
    a: *const LONGLONG,
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

        let a = core::slice::from_raw_parts(a, nx as usize * ny as usize);
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
    _a: &[LONGLONG],
    _ny: c_int,
    _nx: c_int,
    _scale: c_int,
    _output: &mut [c_char],
    _nbytes: &mut c_long,
    _status: &mut c_int,
) -> c_int {
    todo!();
}

/// Tests ported from cfitsio's test_hcompress.c
///
/// rsfitsio delegates H-compress to the external `hcompress` crate
/// (`HCEncoder`/`HCDecoder`), which is what `imcompress` drives internally (the
/// `fits_hcompress_safe`/`fits_hdecompress_safe` shims in this file are unused
/// stubs). The original C test exercised `fits_hcompress`/`fits_hdecompress`
/// directly; here we drive the same round-trips through the crate API.
///
/// NOTE: like the C API, `HCEncoder::write` takes the dimensions in `(ny, nx)`
/// order and modifies its input array in place, and `HCDecoder::read` returns
/// `(nx, ny, scale)`.
#[cfg(test)]
mod tests {
    use hcompress::read::HCDecoder;
    use hcompress::write::HCEncoder;

    use crate::fitsio::LONGLONG;

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
}
