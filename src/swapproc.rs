//! Byte-swapping utilities used by the other FITSIO routines.
//!
//! FITS stores numbers most significant byte first, so on a little-endian
//! machine ([`BYTESWAPPED`](crate::fitsio2::BYTESWAPPED)) every numeric value
//! is swapped on the way in and out.
//!
//! Ported from CFITSIO's `swapproc.c`, written by William Pence at the High
//! Energy Astrophysics Science Archive Research Center (HEASARC), NASA Goddard
//! Space Flight Center.
#![warn(missing_docs)]

use crate::c_types::c_long;

use crate::fitsio::INT32BIT;

/// Swaps the bytes of `nvals` shorts, one at a time.
fn ffswap2_slow(svalues: &mut [i16], nvals: c_long) {
    svalues
        .iter_mut()
        .take(nvals as usize)
        .for_each(|x| *x = (*x).swap_bytes());
}

/// Swap the bytes in the input short integers: ( 0 1 -> 1 0 )
///
/// # Parameters
///
/// * `svalues` — (IO) pointer to shorts to be swapped
/// * `nvals`   — (I) number of shorts to be swapped
pub(crate) fn ffswap2(svalues: &mut [i16], nvals: c_long) {
    ffswap2_slow(svalues, nvals);
}

/// Swaps the bytes of `nvals` 4-byte integers, one at a time.
fn ffswap4_slow(svalues: &mut [INT32BIT], nvals: c_long) {
    //let mut v = std::slice::from_raw_parts_mut(svalues, nvals as usize);
    svalues
        .iter_mut()
        .take(nvals as usize)
        .for_each(|x| *x = (*x).swap_bytes());
}

/// Swap the bytes in the input 4-byte integers: ( 0 1 2 3 -> 3 2 1 0 )
///
/// # Parameters
///
/// * `ivalues` — (IO) pointer to INT*4 to be swapped
/// * `nvals`   — (I) number of floats to be swapped
pub(crate) fn ffswap4(ivalues: &mut [INT32BIT], nvals: c_long) {
    ffswap4_slow(ivalues, nvals);
}

/// Swaps the bytes of `nvals` 8-byte values, one at a time.
fn ffswap8_slow(svalues: &mut [i64], nvals: c_long) {
    svalues
        .iter_mut()
        .take(nvals as usize)
        .for_each(|x| *x = (*x).swap_bytes());
}

/// Swap the bytes in the input doubles: ( 01234567  -> 76543210 )
///
/// # Parameters
///
/// * `ivalues` — (IO) pointer to doubles to be swapped
/// * `nvals`   — (I) number of doubles to be swapped
pub(crate) fn ffswap8(ivalues: &mut [i64], nvals: c_long) {
    ffswap8_slow(ivalues, nvals);
}
