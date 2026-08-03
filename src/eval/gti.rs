//! Good-time-interval search.
//!
//! A transcription of `Search_GTI` in `eval_y.rs`. The intervals themselves are
//! read from the file while the arena is built and handed here by the lowering
//! -- see `parser::lower::GtiData` -- so this module is only the per-row
//! lookup.

use crate::c_types::{c_double, c_long};

/// The index of the interval containing `time`, or `-1`.
///
/// Two searches, as the C has: a bisection when the intervals are known to be
/// ordered and there are enough of them to pay for it, and a linear scan from
/// the last interval backwards otherwise. The bisection is transcribed rather
/// than replaced with a standard binary search, because its stepping is not a
/// textbook one -- it halves a fixed step and leans on the intervals being
/// disjoint -- and the two only agree on well-formed input.
pub(crate) fn search(
    time: c_double,
    start: &[c_double],
    stop: &[c_double],
    ordered: bool,
) -> c_long {
    let n = start.len().min(stop.len()) as c_long;
    if n == 0 {
        return -1;
    }
    if ordered && n > 15 {
        bisect(time, start, stop, n)
    } else {
        scan(time, start, stop, n)
    }
}

fn bisect(time: c_double, start: &[c_double], stop: &[c_double], n: c_long) -> c_long {
    let (lo, hi) = (start[0], stop[(n - 1) as usize]);
    if time < lo || time > hi {
        return -1;
    }
    let mut step = n >> 1;
    let mut gti = step;
    loop {
        if step > 1 {
            step >>= 1;
        }
        let g = gti as usize;
        if time > stop[g] {
            /* past this interval: move on only if the next one can hold it */
            if gti + 1 < n && time >= start[(gti + 1) as usize] {
                gti += step;
            } else {
                return -1;
            }
        } else if time < start[g] {
            if gti >= 1 && time <= stop[(gti - 1) as usize] {
                gti -= step;
            } else {
                return -1;
            }
        } else {
            return gti;
        }
        if gti < 0 || gti >= n {
            return -1;
        }
    }
}

fn scan(time: c_double, start: &[c_double], stop: &[c_double], n: c_long) -> c_long {
    for gti in (0..n).rev() {
        let g = gti as usize;
        if time >= start[g] && time <= stop[g] {
            return gti;
        }
    }
    -1
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Small enough to take the linear path.
    const S: [c_double; 3] = [10.0, 20.0, 30.0];
    const E: [c_double; 3] = [12.0, 22.0, 32.0];

    #[test]
    fn a_time_inside_an_interval_gives_its_index() {
        assert_eq!(search(10.0, &S, &E, false), 0);
        assert_eq!(search(21.0, &S, &E, false), 1);
        assert_eq!(search(32.0, &S, &E, false), 2);
    }

    #[test]
    fn the_bounds_are_inclusive_at_both_ends() {
        assert_eq!(search(12.0, &S, &E, false), 0);
        assert_eq!(search(20.0, &S, &E, false), 1);
    }

    #[test]
    fn a_time_in_no_interval_gives_minus_one() {
        for t in [9.0, 13.0, 25.0, 33.0] {
            assert_eq!(search(t, &S, &E, false), -1, "t={t}");
        }
    }

    #[test]
    fn an_empty_gti_admits_nothing() {
        assert_eq!(search(10.0, &[], &[], true), -1);
        assert_eq!(search(10.0, &[], &[], false), -1);
    }

    /// Above fifteen intervals the ordered path bisects; it must agree with
    /// the scan everywhere, which is the only reason it is safe to have two.
    #[test]
    fn the_two_searches_agree_wherever_both_apply() {
        let start: Vec<c_double> = (0..40).map(|i| f64::from(i) * 10.0).collect();
        let stop: Vec<c_double> = start.iter().map(|s| s + 4.0).collect();
        let mut checked = 0;
        for t in 0..4000 {
            let t = f64::from(t) * 0.1;
            let a = search(t, &start, &stop, true);
            let b = scan(t, &start, &stop, start.len() as c_long);
            assert_eq!(a, b, "disagree at t={t}");
            checked += 1;
        }
        assert!(checked > 0);
    }

    #[test]
    fn ordering_is_only_claimed_when_it_holds() {
        /* out of order, so the caller must not pass ordered=true; the scan
        still finds what is there */
        let start = [30.0, 10.0, 20.0];
        let stop = [32.0, 12.0, 22.0];
        assert_eq!(search(21.0, &start, &stop, false), 2);
        assert_eq!(search(31.0, &start, &stop, false), 0);
    }
}
