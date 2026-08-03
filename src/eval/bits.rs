//! Bit-string kernels.
//!
//! A bit string is held as one character per bit -- `'0'`, `'1'` or `'x'` for
//! an undefined bit -- exactly as the engine's `char*` bit streams are, so
//! these are transcriptions of `bitand`, `bitor`, `bitnot`, `bitcmp` and
//! `bitlgte` in `eval_y.rs` rather than reimplementations over packed words.
//!
//! Every binary operation left-pads the shorter operand with `'0'` to the
//! longer one's width first, which is what those functions do with their
//! scratch buffer.

use crate::eval_defs::OpCode;

/// Left-pad `bits` with `'0'` to `width`.
fn padded(bits: &[u8], width: usize) -> Vec<u8> {
    let mut out = vec![b'0'; width.saturating_sub(bits.len())];
    out.extend_from_slice(bits);
    out
}

/// Pad both operands to a common width.
fn align(a: &[u8], b: &[u8]) -> (Vec<u8>, Vec<u8>) {
    let w = a.len().max(b.len());
    (padded(a, w), padded(b, w))
}

/// `a & b`. An undefined bit in *either* operand wins, so `x & 0` is `x`
/// rather than `0` -- the `'x'` test comes first in `bitand`.
pub(crate) fn and(a: &[u8], b: &[u8]) -> Vec<u8> {
    let (a, b) = align(a, b);
    a.iter()
        .zip(&b)
        .map(|(&x, &y)| {
            if x == b'x' || y == b'x' {
                b'x'
            } else if x == b'1' && y == b'1' {
                b'1'
            } else {
                b'0'
            }
        })
        .collect()
}

/// `a | b`. A set bit wins, then a clear one, so `0 | x` is `0` and only
/// `x | x` is undefined -- the order of the tests in `bitor`.
pub(crate) fn or(a: &[u8], b: &[u8]) -> Vec<u8> {
    let (a, b) = align(a, b);
    a.iter()
        .zip(&b)
        .map(|(&x, &y)| {
            if x == b'1' || y == b'1' {
                b'1'
            } else if x == b'0' || y == b'0' {
                b'0'
            } else {
                b'x'
            }
        })
        .collect()
}

/// `!a`, leaving undefined bits undefined.
pub(crate) fn not(a: &[u8]) -> Vec<u8> {
    a.iter()
        .map(|&c| match c {
            b'1' => b'0',
            b'0' => b'1',
            other => other,
        })
        .collect()
}

/// `a + b`: bit strings concatenate.
pub(crate) fn concat(a: &[u8], b: &[u8]) -> Vec<u8> {
    let mut out = a.to_vec();
    out.extend_from_slice(b);
    out
}

/// Equality. Two bits differ only when one is `'0'` and the other `'1'`, so an
/// undefined bit matches anything -- which is why `BITS == b1111xxxx` holds.
pub(crate) fn cmp_eq(a: &[u8], b: &[u8]) -> bool {
    let (a, b) = align(a, b);
    !a.iter()
        .zip(&b)
        .any(|(&x, &y)| (x == b'0' && y == b'1') || (x == b'1' && y == b'0'))
}

/// The pair of values `bitlgte` compares.
///
/// A position where *either* stream has a wildcard takes no part: it is
/// skipped without advancing the bit weight, so the surviving bits compact
/// down. That is what makes `BITS .gt. bxxx100xx` false for `11110000` --
/// only the four positions the mask leaves defined are compared, and there
/// the two are equal.
///
/// The sum accumulates into a C `int`, so a stream wider than that wraps;
/// `wrapping_*` keeps that rather than panicking in debug.
fn values(a: &[u8], b: &[u8]) -> (i32, i32) {
    let (mut v1, mut v2) = (0i32, 0i32);
    let mut nextbit: i32 = 1;
    for (&x, &y) in a.iter().rev().zip(b.iter().rev()) {
        if matches!(x, b'x' | b'X') || matches!(y, b'x' | b'X') {
            continue;
        }
        if x == b'1' {
            v1 = v1.wrapping_add(nextbit);
        }
        if y == b'1' {
            v2 = v2.wrapping_add(nextbit);
        }
        nextbit = nextbit.wrapping_mul(2);
    }
    (v1, v2)
}

/// Ordering comparisons, which `bitlgte` does on the two streams' values.
pub(crate) fn cmp_ord(a: &[u8], op: OpCode, b: &[u8]) -> bool {
    let (a, b) = align(a, b);
    let (v1, v2) = values(&a, &b);
    match op {
        OpCode::Lt => v1 < v2,
        OpCode::Lte => v1 <= v2,
        OpCode::Gt => v1 > v2,
        OpCode::Gte => v1 >= v2,
        _ => false,
    }
}

/// The number of set bits, which `ACCUM` over a bit string counts.
pub(crate) fn count_ones(bits: &[u8]) -> i64 {
    bits.iter().filter(|&&c| c == b'1').count() as i64
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn operands_of_different_widths_are_left_padded() {
        assert_eq!(and(b"1", b"0001"), b"0001".to_vec());
        assert_eq!(and(b"1", b"0000"), b"0000".to_vec());
        assert_eq!(or(b"1", b"0000"), b"0001".to_vec());
    }

    #[test]
    fn an_undefined_bit_follows_the_engines_test_order() {
        /* x wins in AND, but in OR a defined bit wins over x */
        assert_eq!(and(b"x0", b"00"), b"x0".to_vec());
        assert_eq!(or(b"x0", b"00"), b"00".to_vec());
        assert_eq!(or(b"xx", b"xx"), b"xx".to_vec());
        assert_eq!(or(b"x1", b"00"), b"01".to_vec());
    }

    #[test]
    fn not_leaves_undefined_bits_alone() {
        assert_eq!(not(b"10x"), b"01x".to_vec());
    }

    #[test]
    fn an_undefined_bit_matches_either_value() {
        assert!(cmp_eq(b"11110000", b"1111xxxx"));
        assert!(cmp_eq(b"11110000", b"xxxx0000"));
        assert!(!cmp_eq(b"11110000", b"00001111"));
    }

    #[test]
    fn equality_pads_before_comparing() {
        assert!(cmp_eq(b"1", b"0001"));
        assert!(!cmp_eq(b"1", b"0000"));
    }

    #[test]
    fn ordering_compares_the_streams_as_numbers() {
        assert!(cmp_ord(b"11110000", OpCode::Gt, b"00000000"));
        assert!(cmp_ord(b"11110000", OpCode::Gte, b"11110000"));
        assert!(cmp_ord(b"11110000", OpCode::Lt, b"11111111"));
        assert!(!cmp_ord(b"11111111", OpCode::Lt, b"11110000"));
    }

    #[test]
    fn a_wildcard_excludes_that_position_from_an_ordering_test() {
        /* only the four positions bxxx100xx leaves defined are compared, and
        11110000 matches the mask there, so it is not greater */
        assert!(!cmp_ord(b"11110000", OpCode::Gt, b"xxx100xx"));
        assert!(cmp_ord(b"11110000", OpCode::Gte, b"xxx100xx"));
        assert!(cmp_ord(b"11110000", OpCode::Lte, b"1xxxxxxx"));
        /* skipping must not advance the weight, or the one surviving bit
        would keep its original place value and 1xxxxxxx would read as 128 */
        assert!(cmp_ord(b"11110000", OpCode::Gte, b"1xxxxxxx"));
    }

    #[test]
    fn an_uppercase_wildcard_counts_too() {
        assert!(!cmp_ord(b"11110000", OpCode::Gt, b"XXX100XX"));
    }

    #[test]
    fn concatenation_joins_the_streams() {
        assert_eq!(concat(b"101", b"110"), b"101110".to_vec());
    }

    #[test]
    fn accumulating_counts_the_set_bits() {
        assert_eq!(count_ones(b"1111x000"), 4);
    }
}
