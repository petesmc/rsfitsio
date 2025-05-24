use crate::c_types::{c_int, c_uchar};

use crate::fitsio::{LONGLONG, NULL_MSG};

/* ---------------------------------------------------------------------- */
/// Decompress the input byte stream using the H-compress algorithm
///
/// input  - input array of compressed bytes
/// a - pre-allocated array to hold the output uncompressed image
/// nx - returned X axis size
/// ny - returned Y axis size
///
/// NOTE: the nx and ny dimensions as defined within this code are reversed from
/// the usual FITS notation.  ny is the fastest varying dimension, which is
/// usually considered the X axis in the FITS image display
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn fits_hdecompress(
    input: *const c_uchar,
    smooth: c_int,
    a: *mut c_int,
    ny: *mut c_int,
    nx: *mut c_int,
    scale: *mut c_int,
    status: *mut c_int,
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let nx = nx.as_mut().expect(NULL_MSG);
        let ny = ny.as_mut().expect(NULL_MSG);
        let scale = scale.as_mut().expect(NULL_MSG);

        let a = std::slice::from_raw_parts_mut(a, 42);
        let input = std::slice::from_raw_parts(input, 42);

        fits_hdecompress_safe(input, smooth, a, ny, nx, scale, status)
    }
}

/* ---------------------------------------------------------------------- */
/// Decompress the input byte stream using the H-compress algorithm
///
/// input  - input array of compressed bytes
/// a - pre-allocated array to hold the output uncompressed image
/// nx - returned X axis size
/// ny - returned Y axis size
///
/// NOTE: the nx and ny dimensions as defined within this code are reversed from
/// the usual FITS notation.  ny is the fastest varying dimension, which is
/// usually considered the X axis in the FITS image display
pub(crate) fn fits_hdecompress_safe(
    input: &[c_uchar],
    smooth: c_int,
    a: &mut [c_int],
    ny: &mut c_int,
    nx: &mut c_int,
    scale: &mut c_int,
    status: &mut c_int,
) -> c_int {
    todo!();
}

/* ---------------------------------------------------------------------- */
/// Decompress the input byte stream using the H-compress algorithm
///
/// input  - input array of compressed bytes
/// a - pre-allocated array to hold the output uncompressed image
/// nx - returned X axis size
/// ny - returned Y axis size
///
/// NOTE: the nx and ny dimensions as defined within this code are reversed from
/// the usual FITS notation.  ny is the fastest varying dimension, which is
/// usually considered the X axis in the FITS image display
///
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn fits_hdecompress64(
    input: *const c_uchar,
    smooth: c_int,
    a: *mut LONGLONG,
    ny: *mut c_int,
    nx: *mut c_int,
    scale: *mut c_int,
    status: *mut c_int,
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let nx = nx.as_mut().expect(NULL_MSG);
        let ny = ny.as_mut().expect(NULL_MSG);
        let scale = scale.as_mut().expect(NULL_MSG);

        let a = std::slice::from_raw_parts_mut(a, 42);
        let input = std::slice::from_raw_parts(input, 42);

        fits_hdecompress64_safe(input, smooth, a, ny, nx, scale, status)
    }
}

/* ---------------------------------------------------------------------- */
/// Decompress the input byte stream using the H-compress algorithm
///
/// input  - input array of compressed bytes
/// a - pre-allocated array to hold the output uncompressed image
/// nx - returned X axis size
/// ny - returned Y axis size
///
/// NOTE: the nx and ny dimensions as defined within this code are reversed from
/// the usual FITS notation.  ny is the fastest varying dimension, which is
/// usually considered the X axis in the FITS image display
pub(crate) fn fits_hdecompress64_safe(
    input: &[c_uchar],
    smooth: c_int,
    a: &mut [LONGLONG],
    ny: &mut c_int,
    nx: &mut c_int,
    scale: &mut c_int,
    status: &mut c_int,
) -> c_int {
    todo!();
}
