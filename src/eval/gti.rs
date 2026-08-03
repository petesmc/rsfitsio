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
    search_with_next(time, start, stop, ordered).0
}

/// The containing interval and the *next* one at or after `time`, both `-1`
/// when there is none.
///
/// `Search_GTI` reports the second through an out-parameter; `GTI_Over` needs
/// it to know which intervals to sum across.
pub(crate) fn search_with_next(
    time: c_double,
    start: &[c_double],
    stop: &[c_double],
    ordered: bool,
) -> (c_long, c_long) {
    let n = start.len().min(stop.len()) as c_long;
    if n == 0 {
        return (-1, -1);
    }
    let (gti, mut next) = if ordered && n > 15 {
        bisect(time, start, stop, n)
    } else {
        scan(time, start, stop, n)
    };
    if next >= n {
        next = -1;
    }
    (gti, next)
}

fn bisect(time: c_double, start: &[c_double], stop: &[c_double], n: c_long) -> (c_long, c_long) {
    if time < start[0] || time > stop[(n - 1) as usize] {
        /* outside every interval; the next one is the first if we are before
        it, and there is none if we are past the end */
        return (-1, if start[0] > time { 0 } else { -1 });
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
                return (-1, gti + 1);
            }
        } else if time < start[g] {
            if gti >= 1 && time <= stop[(gti - 1) as usize] {
                gti -= step;
            } else {
                return (-1, gti);
            }
        } else {
            return (gti, gti);
        }
        if gti < 0 || gti >= n {
            return (-1, -1);
        }
    }
}

fn scan(time: c_double, start: &[c_double], stop: &[c_double], n: c_long) -> (c_long, c_long) {
    let mut next = -1;
    for gti in (0..n).rev() {
        let g = gti as usize;
        /* the scan runs downwards, so the last interval that still ends at or
        after `time` is the earliest one that does */
        if stop[g] >= time {
            next = gti;
        }
        if time >= start[g] && time <= stop[g] {
            return (gti, next);
        }
    }
    (-1, next)
}

/// `GTIOVERLAP(file, start, stop)`: how much of `[evt_start, evt_stop]` falls
/// inside the good-time intervals.
///
/// A transcription of `GTI_Over`. It always searches as if the intervals were
/// ordered, which the C does too -- `New_GTI` refuses an unordered GTI for this
/// function rather than letting the search go wrong.
pub(crate) fn overlap(
    evt_start: c_double,
    evt_stop: c_double,
    start: &[c_double],
    stop: &[c_double],
) -> c_double {
    let n = start.len().min(stop.len()) as c_long;
    if n == 0 || evt_stop <= evt_start {
        return 0.0;
    }
    let (gti1, next1) = search_with_next(evt_start, start, stop, true);
    let (gti2, next2) = search_with_next(evt_stop, start, stop, true);

    /* nothing at or after either end, or both fall in the same gap */
    if (next1 < 0 && next2 < 0) || (gti1 < 0 && gti2 < 0 && next1 == next2) {
        return 0.0;
    }
    /* wholly inside one interval */
    if gti1 >= 0 && gti1 == gti2 {
        return evt_stop - evt_start;
    }

    let n_max = if next2 < 0 {
        n - 1
    } else if gti2 >= 0 {
        next2
    } else {
        next2 - 1
    };
    /* the guards above mean `next1` is a real index by here; clamping keeps
    that an assumption the compiler can see rather than one it cannot */
    let mut overlap = 0.0;
    let mut gti = next1.max(0);
    while gti <= n_max && gti < n {
        let g = gti as usize;
        let lo = start[g].max(evt_start);
        let hi = stop[g].min(evt_stop);
        overlap += hi - lo;
        gti += 1;
    }
    overlap
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
        let n = start.len() as c_long;
        let mut checked = 0;
        for t in 0..4000 {
            let t = f64::from(t) * 0.1;
            let a = search_with_next(t, &start, &stop, true);
            /* the scan's own answer, put through the same clamp */
            let (g, mut nx) = scan(t, &start, &stop, n);
            if nx >= n {
                nx = -1;
            }
            assert_eq!(a, (g, nx), "disagree at t={t}");
            checked += 1;
        }
        assert!(checked > 0);
    }

    #[test]
    fn an_event_inside_one_interval_overlaps_its_whole_span() {
        assert_eq!(overlap(10.5, 11.5, &S, &E), 1.0);
        /* exactly the interval */
        assert_eq!(overlap(10.0, 12.0, &S, &E), 2.0);
    }

    #[test]
    fn an_event_in_a_gap_overlaps_nothing() {
        assert_eq!(overlap(13.0, 19.0, &S, &E), 0.0);
        /* and before the first interval, and after the last */
        assert_eq!(overlap(0.0, 5.0, &S, &E), 0.0);
        assert_eq!(overlap(40.0, 50.0, &S, &E), 0.0);
    }

    #[test]
    fn an_event_spanning_intervals_sums_only_the_good_parts() {
        /* 10..12 and 20..22 are good; 12..20 is not */
        assert_eq!(overlap(10.0, 22.0, &S, &E), 4.0);
        /* clipped at both ends: 11..12 plus 20..21 */
        assert_eq!(overlap(11.0, 21.0, &S, &E), 2.0);
        /* across all three */
        assert_eq!(overlap(0.0, 100.0, &S, &E), 6.0);
    }

    #[test]
    fn an_empty_or_backwards_event_overlaps_nothing() {
        assert_eq!(overlap(11.0, 11.0, &S, &E), 0.0);
        assert_eq!(overlap(12.0, 10.0, &S, &E), 0.0);
        assert_eq!(overlap(10.0, 20.0, &[], &[]), 0.0);
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
