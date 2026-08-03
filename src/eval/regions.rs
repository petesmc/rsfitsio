//! The geometric predicates and the angular separation.
//!
//! Transcriptions of `bnear`, `circle`, `saobox`, `ellipse` and `angsep_calc`
//! in `eval_y.rs`. All five are pure functions of their arguments, so the only
//! thing the kernels around them add is the rule those functions share: a row
//! is undefined if *any* argument is undefined there.

use crate::c_types::c_double;

/// `NEAR(x, y, tol)`: whether two values are within a tolerance.
pub(crate) fn near(x: c_double, y: c_double, tolerance: c_double) -> bool {
    (x - y).abs() < tolerance
}

/// `CIRCLE(xcen, ycen, rad, x, y)`: whether a point is inside a circle.
pub(crate) fn circle(
    xcen: c_double,
    ycen: c_double,
    rad: c_double,
    xcol: c_double,
    ycol: c_double,
) -> bool {
    let (dx, dy) = (xcol - xcen, ycol - ycen);
    dx * dx + dy * dy <= rad * rad
}

/// The point rotated into the shape's own frame, which `saobox` and `ellipse`
/// both start from. `rot` is in degrees.
fn rotated(
    xcen: c_double,
    ycen: c_double,
    rot: c_double,
    xcol: c_double,
    ycol: c_double,
) -> (c_double, c_double) {
    let theta = rot / 180.0 * core::f64::consts::PI;
    let (xprime, yprime) = (xcol - xcen, ycol - ycen);
    (
        xprime * theta.cos() + yprime * theta.sin(),
        -xprime * theta.sin() + yprime * theta.cos(),
    )
}

/// `BOX(xcen, ycen, xwid, ywid, rot, x, y)`: whether a point is inside a
/// rotated rectangle. The widths are full widths, so the box spans half of
/// each either side of the centre.
#[allow(clippy::too_many_arguments)]
pub(crate) fn saobox(
    xcen: c_double,
    ycen: c_double,
    xwid: c_double,
    ywid: c_double,
    rot: c_double,
    xcol: c_double,
    ycol: c_double,
) -> bool {
    let (x, y) = rotated(xcen, ycen, rot, xcol, ycol);
    x >= -0.5 * xwid && x <= 0.5 * xwid && y >= -0.5 * ywid && y <= 0.5 * ywid
}

/// `ELLIPSE(xcen, ycen, xrad, yrad, rot, x, y)`: whether a point is inside a
/// rotated ellipse. Here the radii are semi-axes, unlike `BOX`'s widths.
#[allow(clippy::too_many_arguments)]
pub(crate) fn ellipse(
    xcen: c_double,
    ycen: c_double,
    xrad: c_double,
    yrad: c_double,
    rot: c_double,
    xcol: c_double,
    ycol: c_double,
) -> bool {
    let (x, y) = rotated(xcen, ycen, rot, xcol, ycol);
    let (dx, dy) = (x / xrad, y / yrad);
    dx * dx + dy * dy <= 1.0
}

/// `ANGSEP(ra1, dec1, ra2, dec2)`: the angle between two sky positions, in
/// degrees.
///
/// This is the law of Haversines rather than the law of Cosines, which is
/// unstable for angles around 0.1 arcsec, and `a` is clamped before the square
/// roots so that rounding cannot put it outside `[0, 1]`.
pub(crate) fn angsep(ra1: c_double, dec1: c_double, ra2: c_double, dec2: c_double) -> c_double {
    let deg = 4.0 * (1.0_f64).atan() / 180.0;
    let sra = ((ra2 - ra1) * deg / 2.0).sin();
    let sdec = ((dec2 - dec1) * deg / 2.0).sin();
    let a = sdec * sdec + (dec1 * deg).cos() * (dec2 * deg).cos() * sra * sra;
    let a = a.clamp(0.0, 1.0);
    2.0 * a.sqrt().atan2((1.0 - a).sqrt()) / deg
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn near_compares_against_the_tolerance() {
        assert!(near(1.0, 1.05, 0.1));
        assert!(!near(1.0, 1.2, 0.1));
        /* the bound is strict */
        assert!(!near(1.0, 1.1, 0.1));
    }

    #[test]
    fn a_point_on_the_circle_is_inside_it() {
        assert!(circle(0.0, 0.0, 5.0, 3.0, 4.0), "on the boundary counts");
        assert!(circle(0.0, 0.0, 5.0, 0.0, 0.0));
        assert!(!circle(0.0, 0.0, 5.0, 4.0, 4.0));
    }

    #[test]
    fn a_box_spans_half_its_width_either_side() {
        assert!(saobox(0.0, 0.0, 4.0, 2.0, 0.0, 2.0, 1.0), "corner counts");
        assert!(!saobox(0.0, 0.0, 4.0, 2.0, 0.0, 2.1, 0.0));
        assert!(!saobox(0.0, 0.0, 4.0, 2.0, 0.0, 0.0, 1.1));
    }

    #[test]
    fn rotating_a_box_moves_what_is_inside_it() {
        /* a tall thin box, and a point that only falls inside once it turns */
        assert!(!saobox(0.0, 0.0, 1.0, 8.0, 0.0, 3.0, 0.0));
        assert!(saobox(0.0, 0.0, 1.0, 8.0, 90.0, 3.0, 0.0));
    }

    #[test]
    fn an_ellipse_measures_in_semi_axes() {
        /* unlike BOX, the radii reach the boundary rather than half of it */
        assert!(ellipse(0.0, 0.0, 4.0, 2.0, 0.0, 4.0, 0.0));
        assert!(ellipse(0.0, 0.0, 4.0, 2.0, 0.0, 0.0, 2.0));
        assert!(!ellipse(0.0, 0.0, 4.0, 2.0, 0.0, 4.1, 0.0));
        assert!(!ellipse(0.0, 0.0, 4.0, 2.0, 0.0, 3.0, 1.9));
    }

    #[test]
    fn angsep_measures_along_the_sphere() {
        /* the pole is 90 degrees from anywhere on the equator */
        assert!((angsep(0.0, 0.0, 0.0, 90.0) - 90.0).abs() < 1e-9);
        assert!((angsep(0.0, 0.0, 180.0, 0.0) - 180.0).abs() < 1e-9);
        assert!(angsep(12.0, 34.0, 12.0, 34.0).abs() < 1e-12);
    }

    #[test]
    fn angsep_stays_stable_for_tiny_separations() {
        /* the reason for the Haversine form: a Law of Cosines version loses
        all its precision here */
        let tiny = 1.0 / 3_600_000.0; /* one milliarcsecond, in degrees */
        let got = angsep(10.0, 20.0, 10.0 + tiny, 20.0);
        let want = tiny * (20.0_f64).to_radians().cos();
        assert!((got - want).abs() < want * 1e-6, "got {got}, want ~{want}");
    }

    #[test]
    fn angsep_is_symmetric() {
        let a = angsep(10.0, -20.0, 250.0, 65.0);
        let b = angsep(250.0, 65.0, 10.0, -20.0);
        assert!((a - b).abs() < 1e-12);
    }
}
