/* H-decompress routines */

use hcompress::read::HCDecoder;

use crate::c_types::{c_int, c_uchar};

use crate::fitsio::{DATA_DECOMPRESSION_ERR, LONGLONG, NULL_MSG};

/* ---------------------------------------------------------------------- */
/// Decompress the input byte stream using the H-compress algorithm
///
/// input  - input array of compressed bytes
/// a - pre-allocated array to hold the output uncompressed image
/// na - number of integers allocated for a (eg, sizeof a / sizeof *a)
/// nx - returned X axis size
/// ny - returned Y axis size
///
/// NOTE: the nx and ny dimensions as defined within this code are reversed from
/// the usual FITS notation.  ny is the fastest varying dimension, which is
/// usually considered the X axis in the FITS image display
///
/// # Safety
///
/// The C API gives no length for `input`, so this wrapper cannot derive one:
/// the caller must guarantee that `input` points at a readable buffer of at
/// least `na * 4 + HDRSIZE` bytes.  That is the largest stream a valid
/// H-compressed `na`-pixel image can occupy, so a stream produced by
/// `fits_hcompress` always fits.  Prefer [`fits_hdecompress_safe`], which takes
/// the compressed stream as a slice and needs no such assumption.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn fits_hdecompress(
    input: *const c_uchar,
    smooth: c_int,
    a: *mut c_int,
    na: c_int,
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

        let a = core::slice::from_raw_parts_mut(a, na as usize);
        // SAFETY: bound documented above -- the C signature carries no input length.
        let input = core::slice::from_raw_parts(input, na as usize * 4 + HDRSIZE);

        fits_hdecompress_safe(input, smooth, a, na, ny, nx, scale, status)
    }
}

/// Size of the H-compress stream header (magic + nx + ny + scale), plus slack
/// for the trailing quadtree/bit-plane bookkeeping.
const HDRSIZE: usize = 1024;

/* ---------------------------------------------------------------------- */
/// Decompress the input byte stream using the H-compress algorithm
///
/// input  - input array of compressed bytes
/// a - pre-allocated array to hold the output uncompressed image
/// na - number of integers allocated for a (eg, sizeof a / sizeof *a)
/// nx - returned X axis size
/// ny - returned Y axis size
///
/// NOTE: the nx and ny dimensions as defined within this code are reversed from
/// the usual FITS notation.  ny is the fastest varying dimension, which is
/// usually considered the X axis in the FITS image display
pub fn fits_hdecompress_safe(
    input: &[c_uchar],
    smooth: c_int,
    a: &mut [c_int],
    na: c_int,
    ny: &mut c_int,
    nx: &mut c_int,
    scale: &mut c_int,
    status: &mut c_int,
) -> c_int {
    if *status > 0 {
        return *status;
    }

    /* decode the input array */

    /* The C's decode() / undigitize() / hinv() are not in this crate -- they
    live in the external `hcompress` crate, where HCDecoder::read() is exactly
    that same sequence.  The C took FFLOCK around decode() because decode used
    the `nextchar` global; the crate reads through a Cursor it owns, so no lock
    is needed here. */
    let na = (na as usize).min(a.len());

    match HCDecoder::new().read(input, smooth, &mut a[..na]) {
        Ok((rnx, rny, rscale)) => {
            /* decode() returns the dimensions and the digitization scale */
            *nx = rnx as c_int;
            *ny = rny as c_int;
            *scale = rscale;
        }
        Err(_) => {
            /* every failure path in the C's decode()/hinv() returns
            DATA_DECOMPRESSION_ERR */
            *status = DATA_DECOMPRESSION_ERR;
        }
    }

    *status
}

/* ---------------------------------------------------------------------- */
/// Decompress the input byte stream using the H-compress algorithm
///
/// input  - input array of compressed bytes
/// a - pre-allocated array to hold the output uncompressed image
/// na - number of integers allocated for a (eg, sizeof a / sizeof *a)
/// nx - returned X axis size
/// ny - returned Y axis size
///
/// NOTE: the nx and ny dimensions as defined within this code are reversed from
/// the usual FITS notation.  ny is the fastest varying dimension, which is
/// usually considered the X axis in the FITS image display
///
/// # Safety
///
/// As for [`fits_hdecompress`], the caller must guarantee `input` points at a
/// readable buffer of at least `na * 8 + HDRSIZE` bytes.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn fits_hdecompress64(
    input: *const c_uchar,
    smooth: c_int,
    a: *mut LONGLONG,
    na: c_int,
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

        let a = core::slice::from_raw_parts_mut(a, na as usize);
        // SAFETY: bound documented above -- the C signature carries no input length.
        let input = core::slice::from_raw_parts(input, na as usize * 8 + HDRSIZE);

        fits_hdecompress64_safe(input, smooth, a, na, ny, nx, scale, status)
    }
}

/* ---------------------------------------------------------------------- */
/// Decompress the input byte stream using the H-compress algorithm
///
/// input  - input array of compressed bytes
/// a - pre-allocated array to hold the output uncompressed image
/// na - number of integers allocated for a (eg, sizeof a / sizeof *a)
/// nx - returned X axis size
/// ny - returned Y axis size
///
/// NOTE: the nx and ny dimensions as defined within this code are reversed from
/// the usual FITS notation.  ny is the fastest varying dimension, which is
/// usually considered the X axis in the FITS image display
///
pub fn fits_hdecompress64_safe(
    input: &[c_uchar],
    smooth: c_int,
    a: &mut [LONGLONG],
    na: c_int,
    ny: &mut c_int,
    nx: &mut c_int,
    scale: &mut c_int,
    status: &mut c_int,
) -> c_int {
    if *status > 0 {
        return *status;
    }

    /* decode the input array */

    /* HCDecoder::read64() is the C's decode64() + undigitize64() + hinv64(),
    and it also performs the C's trailing "pack the I*8 values back into an I*4
    array" loop, so there is nothing left to do here.  As above, the C's FFLOCK
    guarded decode64's `nextchar` global, which the crate does not have. */
    let na = (na as usize).min(a.len());

    match HCDecoder::new().read64(input, smooth, &mut a[..na]) {
        Ok((rnx, rny, rscale)) => {
            *nx = rnx as c_int;
            *ny = rny as c_int;
            *scale = rscale;
        }
        Err(_) => {
            *status = DATA_DECOMPRESSION_ERR;
        }
    }

    *status
}
