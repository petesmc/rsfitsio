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
/// * The length of the list, or `None` if `lldst` is too small to hold it.
///   Size it with [`pl_p2li_max_len`] and it cannot be.
pub fn pl_p2li(pxsrc: &[i32], xs: i32, lldst: &mut [i16], npix: usize) -> Option<usize> {
    pliocomp::pl_p2li(pxsrc, xs, lldst, npix)
}

/// The worst-case length, in 16-bit words, of the line list `pl_p2li` produces
/// for `npix` pixels: the counterpart of `imcomp_calc_max_elem`.  An encode
/// buffer sized with this cannot overflow.
pub fn pl_p2li_max_len(npix: usize) -> usize {
    pliocomp::pl_p2li_max_len(npix)
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

/// Tests ported from cfitsio's test_pliocomp.c
///
/// NOTE: The C `pl_p2li`/`pl_l2pi` use a 1-based `xs` starting index (the
/// Fortran-derived original applies an internal `--pxsrc`/`--ll_src`). The Rust
/// `pliocomp` crate used here is 0-based instead — which is how `imcompress`
/// drives it (it always passes `xs = 0`). The C tests below are therefore ported
/// with `xs = 0` to start at the first pixel, and offsets translated by -1.
#[cfg(test)]
mod tests {
    use super::*;

    /// Encode `input`, decode the full length back, and check exact roundtrip.
    fn roundtrip(input: &[i32]) {
        let siz = input.len();
        let mut linelist = [0i16; 100];

        let nbytes = pl_p2li(input, 0, &mut linelist, siz).expect("buffer was sized with pl_p2li_max_len");
        assert!(nbytes > 0, "pl_p2li returned zero-length list");

        let mut output = vec![0i32; siz];
        let ndecoded = pl_l2pi(&linelist, 0, &mut output, siz);
        assert_eq!(ndecoded, siz, "pl_l2pi did not decode all pixels");

        assert_eq!(&output[..], input, "roundtrip mismatch");
    }

    #[test]
    fn test_empty_input() {
        let pixels = [0i32];
        let mut linelist = [0i16; 100];

        // Encoding zero pixels produces an empty (zero-length) line list.
        assert_eq!(pl_p2li(&pixels, 0, &mut linelist, 0), Some(0));
    }

    #[test]
    fn test_data() {
        roundtrip(&[0, 0, 0, 0, 0, 0, 0, 0, 0, 0]);
        roundtrip(&[42, 42, 42, 42, 42]);
        roundtrip(&[0, 1, 2, 3, 4, 5, 6, 7]);
        roundtrip(&[0, 0, 5, 5, 0, 0, 0, 10, 0, 0]);
        roundtrip(&[0, 100000, 0]);
        roundtrip(&[1, 10000, 1, 10000]);
        roundtrip(&[99]);
        roundtrip(&[1, 0, 1, 0, 1, 0]);
        roundtrip(&[100, 80, 60, 40, 20]);
        roundtrip(&[10000, 1, 10000]);
    }

    #[test]
    fn test_negative_clamp() {
        // Negative pixel values are clamped to zero by the encoder.
        let pixels = [-5i32, 10, -3];
        let mut linelist = [0i16; 100];
        let mut output = [0i32; 3];

        assert!(pl_p2li(&pixels, 0, &mut linelist, 3).unwrap() > 0);
        assert_eq!(pl_l2pi(&linelist, 0, &mut output, 3), 3);

        assert_eq!(output[0], 0);
        assert_eq!(output[1], 10);
        assert_eq!(output[2], 0);
    }

    #[test]
    fn test_partial_decode() {
        // Encode 10 pixels but only decode the first 5.
        let pixels = [1i32, 2, 3, 4, 5, 6, 7, 8, 9, 10];
        let mut linelist = [0i16; 100];
        let mut output = [0i32; 5];

        assert!(pl_p2li(&pixels, 0, &mut linelist, 10).unwrap() > 0);
        assert_eq!(pl_l2pi(&linelist, 0, &mut output, 5), 5);

        for i in 0..5 {
            assert_eq!(output[i], pixels[i]);
        }
    }

    #[test]
    fn test_l2pi_empty() {
        // Decoding zero pixels returns zero, regardless of the list contents.
        let mut linelist = [0i16; 10];
        linelist[3] = 0;
        let mut output = [0i32; 10];

        assert_eq!(pl_l2pi(&linelist, 0, &mut output, 0), 0);
    }

    #[test]
    fn test_decode_more_than_encoded() {
        // Requesting more pixels than were encoded zero-fills the remainder.
        let pixels = [1i32, 2, 3, 4, 5];
        let mut linelist = [0i16; 100];
        let mut output = [0i32; 10];

        assert!(pl_p2li(&pixels, 0, &mut linelist, 5).unwrap() > 0);
        assert_eq!(pl_l2pi(&linelist, 0, &mut output, 10), 10);

        for i in 0..5 {
            assert_eq!(output[i], pixels[i]);
        }
        for i in 5..10 {
            assert_eq!(output[i], 0);
        }
    }
}
