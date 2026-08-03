//! String kernels.
//!
//! Transcriptions of `Do_BinOp_str`, `cstrmid` and `str_pos_*` in `eval_y.rs`.
//! A string value holds one entry per row, and the entry is the text up to the
//! NUL rather than the column's full declared width -- that width is carried
//! separately because `STRMID` measures against it, not against the text.

use core::cmp::Ordering;

/// `strcmp`'s ordering.
///
/// `Do_BinOp_str` shortcuts on the first character, comparing it as a *signed*
/// `char` before falling back to `strcmp`, which compares unsigned. The two
/// agree on ASCII and disagree only when the first bytes differ and one is
/// 0x80 or above; the shortcut is reproduced rather than simplified away so
/// that a high-byte column orders the same way it does today.
pub(crate) fn compare(a: &[u8], b: &[u8]) -> Ordering {
    let first = |s: &[u8]| s.first().copied().unwrap_or(0) as i8;
    match first(a).cmp(&first(b)) {
        Ordering::Equal => a.cmp(b),
        other => other,
    }
}

/// `a + b`: strings concatenate.
pub(crate) fn concat(a: &[u8], b: &[u8]) -> Vec<u8> {
    let mut out = a.to_vec();
    out.extend_from_slice(b);
    out
}

/// `STRSTR(a, b)`: the 1-based position of `b` in `a`, or `None` when it does
/// not occur -- which the engine reports as a *null*, not as 0.
pub(crate) fn find(a: &[u8], b: &[u8]) -> Option<i64> {
    if b.is_empty() {
        return Some(1);
    }
    a.windows(b.len())
        .position(|w| w == b)
        .map(|i| i as i64 + 1)
}

/// `STRMID(s, pos, len)`, per `cstrmid`.
///
/// `pos` is 1-based and measured against `src_len` -- the column's declared
/// width, not the text's length -- so a position inside the declared width but
/// past the text yields the NUL padding there, which reads as empty. A `pos`
/// past the width gives an empty string, and `pos == 0` is undefined (the
/// caller checks that before calling). A negative `pos` is a parse-time error,
/// so it cannot reach here.
pub(crate) fn mid(s: &[u8], src_len: usize, pos: usize, len: usize) -> Vec<u8> {
    if pos > src_len || pos == 0 {
        return Vec::new();
    }
    /* the source is NUL-padded out to its declared width */
    let take = if pos - 1 + len > src_len {
        src_len - pos + 1
    } else {
        len
    };
    let mut out = Vec::with_capacity(take);
    for i in 0..take {
        match s.get(pos - 1 + i) {
            /* padding, and the copy stops at the first NUL as the C's
            null-terminate does */
            None => break,
            Some(&c) => out.push(c),
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn comparison_orders_bytewise() {
        assert_eq!(compare(b"abc", b"abc"), Ordering::Equal);
        assert_eq!(compare(b"abc", b"abd"), Ordering::Less);
        assert_eq!(compare(b"beta", b"alpha"), Ordering::Greater);
        /* a prefix sorts before the longer string */
        assert_eq!(compare(b"ab", b"abc"), Ordering::Less);
        assert_eq!(compare(b"", b"a"), Ordering::Less);
    }

    #[test]
    fn concatenation_joins_the_text() {
        assert_eq!(concat(b"alpha", b"!"), b"alpha!".to_vec());
    }

    #[test]
    fn find_reports_a_one_based_position() {
        assert_eq!(find(b"abc", b"b"), Some(2));
        assert_eq!(find(b"beta", b"ta"), Some(3));
        assert_eq!(find(b"abc", b"abc"), Some(1));
    }

    #[test]
    fn find_reports_absence_rather_than_zero() {
        /* the engine marks the row undefined, so this must be distinguishable
        from a real position */
        assert_eq!(find(b"abc", b"z"), None);
        assert_eq!(find(b"a", b"alpha"), None);
    }

    #[test]
    fn mid_takes_a_one_based_substring() {
        assert_eq!(mid(b"alpha", 10, 1, 3), b"alp".to_vec());
        assert_eq!(mid(b"alpha", 10, 2, 3), b"lph".to_vec());
        assert_eq!(mid(b"abcdef", 6, 2, 3), b"bcd".to_vec());
    }

    #[test]
    fn mid_stops_at_the_end_of_the_text() {
        /* 'alpha' in a width-10 column: asking past the text gives what is
        there, which is nothing */
        assert_eq!(mid(b"alpha", 10, 3, 5), b"pha".to_vec());
        assert_eq!(mid(b"alpha", 10, 4, 5), b"ha".to_vec());
        assert_eq!(mid(b"alpha", 10, 8, 2), Vec::<u8>::new());
    }

    #[test]
    fn mid_outside_the_declared_width_is_empty() {
        assert_eq!(mid(b"alpha", 10, 11, 3), Vec::<u8>::new());
        assert_eq!(mid(b"alpha", 10, 0, 3), Vec::<u8>::new());
    }

    #[test]
    fn mid_truncates_to_the_declared_width() {
        assert_eq!(mid(b"abcdef", 6, 5, 10), b"ef".to_vec());
    }
}
