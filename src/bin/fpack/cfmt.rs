/* C printf conversions that Rust's format!() cannot express: %g in any form,
 * and the `#' alternate flag.  Everything format!() can do is left to it at
 * the call sites; the double conversions come through here so that there is
 * one rounding and one specials-spelling engine, and it is the C one.
 */

use std::ffi::CStr;

use bytemuck::cast_slice;
use rsfitsio::c_types::{c_char, c_long, c_ulong};
use rsfitsio::relibc::header::stdio::printf::{CustomVaList, VaArg};
use rsfitsio::relibc::header::stdio::snprintf_va;

/* long enough for %f of any finite f64 */
const CFMT_BUF: usize = 512;

/// One C floating-point conversion.  `spec` is the whole conversion exactly as
/// it appears in the C format string, e.g. `dbl(c"%#8.2g", x)`.
pub(crate) fn dbl(spec: &CStr, val: f64) -> String {
    let mut buf: [c_char; CFMT_BUF] = [0; CFMT_BUF];
    let mut args = CustomVaList::new();
    args.push(VaArg::c_double(val));

    snprintf_va(
        &mut buf,
        CFMT_BUF,
        cast_slice(spec.to_bytes_with_nul()),
        args,
    );

    /* the conversions used here are all ASCII, so this cannot fail */
    CStr::from_bytes_until_nul(cast_slice(&buf))
        .unwrap()
        .to_string_lossy()
        .into_owned()
}

/// C: `(unsigned long) (~((int) hdusum))`, from fp_info_hdu.  The narrowing to
/// a signed int is what makes the result sign-extend.
pub(crate) fn not_int(hdusum: c_ulong) -> c_ulong {
    (!(hdusum as u32 as i32)) as c_long as c_ulong
}

#[cfg(test)]
// 3.14159265 below is a printf input paired with the exact digits C prints for
// it, not a use of pi as a constant; substituting f64::consts::PI would change
// the value being formatted.
#[allow(clippy::approx_constant)]
mod tests {
    use super::*;

    /* Every expected string below was produced by C printf on this machine
    (gcc/glibc), not by this implementation.  See the `golden' program in the
    port notes: the inputs are the ones fpack's statistics actually produce --
    zero, unity, sub-normal noise estimates, megapixel counts and the NaN/inf
    that a division by a zero file size yields. */

    #[test]
    fn test_g_conversions_match_c() {
        // (value, %#8.2g, %#7.3g, %#10.4g)
        let cases: &[(f64, &str, &str, &str)] = &[
            (0.0, "     0.0", "   0.00", "     0.000"),
            (1.0, "     1.0", "   1.00", "     1.000"),
            (-1.0, "    -1.0", "  -1.00", "    -1.000"),
            (0.5, "    0.50", "  0.500", "    0.5000"),
            (1234.5, " 1.2e+03", "1.23e+03", "     1234."),
            (1e-5, " 1.0e-05", "1.00e-05", " 1.000e-05"),
            (1e5, " 1.0e+05", "1.00e+05", " 1.000e+05"),
            (123456.0, " 1.2e+05", "1.23e+05", " 1.235e+05"),
            (9.999e-5, " 0.00010", "0.000100", " 9.999e-05"),
            (0.000123456, " 0.00012", "0.000123", " 0.0001235"),
            (3.14159265, "     3.1", "   3.14", "     3.142"),
            (100.0, " 1.0e+02", "   100.", "     100.0"),
            (1e10, " 1.0e+10", "1.00e+10", " 1.000e+10"),
            (-0.0025, " -0.0025", "-0.00250", " -0.002500"),
            (2.0, "     2.0", "   2.00", "     2.000"),
            (1e-300, "1.0e-300", "1.00e-300", "1.000e-300"),
        ];

        for &(v, g82, g73, g104) in cases {
            assert_eq!(dbl(c"%#8.2g", v), g82, "%#8.2g of {v}");
            assert_eq!(dbl(c"%#7.3g", v), g73, "%#7.3g of {v}");
            assert_eq!(dbl(c"%#10.4g", v), g104, "%#10.4g of {v}");
        }
    }

    #[test]
    fn test_f_conversions_match_c() {
        // (value, %#5.1f, %#6.2f, %8.1f, %6.0f, %5.3f, %7.2f)
        let cases: &[(f64, &str, &str, &str, &str, &str, &str)] = &[
            (
                0.0, "  0.0", "  0.00", "     0.0", "     0", "0.000", "   0.00",
            ),
            (
                1.0, "  1.0", "  1.00", "     1.0", "     1", "1.000", "   1.00",
            ),
            (
                -1.0, " -1.0", " -1.00", "    -1.0", "    -1", "-1.000", "  -1.00",
            ),
            (
                0.5, "  0.5", "  0.50", "     0.5", "     0", "0.500", "   0.50",
            ),
            (
                1234.5, "1234.5", "1234.50", "  1234.5", "  1234", "1234.500", "1234.50",
            ),
            (
                3.14159265, "  3.1", "  3.14", "     3.1", "     3", "3.142", "   3.14",
            ),
            (
                -0.0025, " -0.0", " -0.00", "    -0.0", "    -0", "-0.003", "  -0.00",
            ),
            (
                123456.0,
                "123456.0",
                "123456.00",
                "123456.0",
                "123456",
                "123456.000",
                "123456.00",
            ),
        ];

        for &(v, f51, f62, f81, f60, f53, f72) in cases {
            assert_eq!(dbl(c"%#5.1f", v), f51, "%#5.1f of {v}");
            assert_eq!(dbl(c"%#6.2f", v), f62, "%#6.2f of {v}");
            assert_eq!(dbl(c"%8.1f", v), f81, "%8.1f of {v}");
            assert_eq!(dbl(c"%6.0f", v), f60, "%6.0f of {v}");
            assert_eq!(dbl(c"%5.3f", v), f53, "%5.3f of {v}");
            assert_eq!(dbl(c"%7.2f", v), f72, "%7.2f of {v}");
        }
    }

    /* fits_read_image_speed divides its four timings by a file size that
    rounds to 0.0 for a small enough image (upstream bug 11), so the report
    really can be handed a NaN or an infinity.  C spells them lower-case. */
    #[test]
    fn test_specials_match_c() {
        assert_eq!(dbl(c"%#8.2g", f64::NAN), "     nan");
        assert_eq!(dbl(c"%#7.3g", f64::NAN), "    nan");
        assert_eq!(dbl(c"%#5.1f", f64::NAN), "  nan");
        assert_eq!(dbl(c"%#8.2g", f64::INFINITY), "     inf");
        assert_eq!(dbl(c"%#7.3g", f64::INFINITY), "    inf");
        assert_eq!(dbl(c"%#5.1f", f64::INFINITY), "  inf");
    }

    #[test]
    fn test_not_int_matches_c_cast_chain() {
        /* C: printf("%lu", (unsigned long)(~((int) hdusum))) */
        let all_ones = c_ulong::MAX; /* 64-bit on unix, 32-bit on windows */
        assert_eq!(not_int(0), all_ones);
        assert_eq!(not_int(1), all_ones - 1);
        assert_eq!(not_int(0x8000_0000), 0x7fff_ffff);
        assert_eq!(not_int(0xffff_ffff), 0);
        /* bits above 32 are discarded by the (int) cast.  Only reachable
        where `unsigned long' is wider than 32 bits, which is not the case on
        windows even at 64-bit -- hence the runtime check rather than a
        target_pointer_width cfg, and the cast rather than a literal that
        would be out of range for a 32-bit c_ulong. */
        if size_of::<c_ulong>() > 4 {
            assert_eq!(not_int(0x1_ffff_ffff_u64 as c_ulong), 0);
        }
    }
}
