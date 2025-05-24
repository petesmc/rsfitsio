use pliocomp;

/// Convert a pixel array to a line list.
///
/// Arguments
///
/// * `pxsrc` - input pixel array
/// * `xs` - starting index in pxsrc (?)
/// * `lldst` - encoded line list
/// * `npix` - number of pixels to convert
///
/// Returns
///
/// * The length of the list is returned as the function value
pub fn pl_p2li(pxsrc: &[i32], xs: i32, lldst: &mut [i16], npix: usize) -> usize {
    pliocomp::pl_p2li(pxsrc, xs, lldst, npix)
}

/// Translate a PLIO line list into an integer pixel array.
///
/// Arguments
///
/// * `ll_src` - encoded line list
/// * `xs` - starting index in ll_src
/// * `px_dst` - output pixel array
/// * `npix` - number of pixels to convert
///
/// Returns
///
/// * The number of pixels output (always npix) is returned as the function value.
pub fn pl_l2pi(ll_src: &[i16], xs: i32, px_dst: &mut [i32], npix: usize) -> usize {
    pliocomp::pl_l2pi(ll_src, xs, px_dst, npix)
}