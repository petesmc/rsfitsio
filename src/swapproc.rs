/*  This file, swapproc.c, contains general utility routines that are      */
/*  used by other FITSIO routines to swap bytes.                           */

/*  The FITSIO software was written by William Pence at the High Energy    */
/*  Astrophysic Science Archive Research Center (HEASARC) at the NASA      */
/*  Goddard Space Flight Center.                                           */

use crate::c_types::c_long;

use crate::fitsio::INT32BIT;

fn ffswap2_slow(svalues: &mut [i16], nvals: c_long) {
    svalues
        .iter_mut()
        .take(nvals as usize)
        .for_each(|x| *x = (*x).swap_bytes());
}

/// Swap the bytes in the input short integers: ( 0 1 -> 1 0 )
pub(crate) fn ffswap2(
    svalues: &mut [i16], /* IO - pointer to shorts to be swapped    */
    nvals: c_long,       /* I  - number of shorts to be swapped     */
) {
    ffswap2_slow(svalues, nvals);
}

fn ffswap4_slow(svalues: &mut [INT32BIT], nvals: c_long) {
    //let mut v = std::slice::from_raw_parts_mut(svalues, nvals as usize);
    svalues
        .iter_mut()
        .take(nvals as usize)
        .for_each(|x| *x = (*x).swap_bytes());
}

/// swap the bytes in the input 4-byte integer: ( 0 1 2 3 -> 3 2 1 0 )
pub(crate) fn ffswap4(
    ivalues: &mut [INT32BIT], /* IO - pointer to INT*4 to be swapped    */
    nvals: c_long,            /* I  - number of floats to be swapped     */
) {
    ffswap4_slow(ivalues, nvals);
}

fn ffswap8_slow(svalues: &mut [i64], nvals: c_long) {
    svalues
        .iter_mut()
        .take(nvals as usize)
        .for_each(|x| *x = (*x).swap_bytes());
}

/// Swap the bytes in the input doubles: ( 01234567  -> 76543210 )
pub(crate) fn ffswap8(
    ivalues: &mut [i64], /* IO - pointer to doubles to be swapped     */
    nvals: c_long,       /* I  - number of doubles to be swapped      */
) {
    ffswap8_slow(ivalues, nvals);
}
