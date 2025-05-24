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
    unsafe {
        let nbytes = nbytes.as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);

        let a = std::slice::from_raw_parts(a, nx as usize * ny as usize);
        let output = std::slice::from_raw_parts_mut(output, *nbytes as usize);

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
pub(crate) fn fits_hcompress_safe(
    a: &[c_int],
    ny: c_int,
    nx: c_int,
    scale: c_int,
    output: &mut [c_char],
    nbytes: &mut c_long,
    status: &mut c_int,
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
    unsafe {
        let nbytes = nbytes.as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);

        let a = std::slice::from_raw_parts(a, nx as usize * ny as usize);
        let output = std::slice::from_raw_parts_mut(output, *nbytes as usize);

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
pub(crate) fn fits_hcompress64_safe(
    a: &[LONGLONG],
    ny: c_int,
    nx: c_int,
    scale: c_int,
    output: &mut [c_char],
    nbytes: &mut c_long,
    status: &mut c_int,
) -> c_int {
    todo!();
}
