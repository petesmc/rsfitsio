//! Per-row values of the bit-string operators in `src/eval_y.rs`.
//!
//! `Do_BinOp_bit`, `Do_Unary` and the `bitand`/`bitor`/`bitnot`/`bitcmp`/
//! `bitlgte` helpers were converted from raw `*mut c_char` walks to slices.
//! The `eval_corpus` golden covers the *shapes* of these expressions but not
//! their values -- its probe skips the value of a `TSTRING`-typed result, and a
//! bit-string result is typed `TSTRING` -- so the arithmetic itself is pinned
//! here, one boolean per row.
//!
//! Expected values are derived by hand from the three BITS rows, `0xF0`,
//! `0x0F` and `0xAA`, i.e. `11110000`, `00001111` and `10101010`.

mod common;

#[cfg(test)]
mod tests {
    use rsfitsio::aliases::rust_api::*;
    use rsfitsio::c_types::{c_char, c_long};
    use rsfitsio::fitsio::{BINARY_TBL, BYTE_IMG, LONGLONG, fitsfile};

    const NROWS: c_long = 3;

    fn cc(s: &str) -> Vec<c_char> {
        let mut v: Vec<c_char> = s.bytes().map(|b| b as c_char).collect();
        v.push(0);
        v
    }

    /// A table with one 8-bit column, the same three rows the corpus uses.
    fn bits_table() -> Box<fitsfile> {
        let mut status = 0;
        let mut f: Option<Box<fitsfile>> = None;
        fits_create_file(&mut f, &cc("mem://bits.fits"), &mut status);
        fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);

        let tt = [cc("BITS"), cc("OTHER")];
        let tf = [cc("8X"), cc("8X")];
        let tt_ref: Vec<Option<&[c_char]>> = tt.iter().map(|v| Some(v.as_slice())).collect();
        let tf_ref: Vec<&[c_char]> = tf.iter().map(|v| v.as_slice()).collect();
        fits_create_tbl(
            f.as_deref_mut().unwrap(),
            BINARY_TBL,
            NROWS as LONGLONG,
            2,
            &tt_ref,
            &tf_ref,
            None,
            None,
            &mut status,
        );

        let fp = f.as_deref_mut().unwrap();
        /* 11110000, 00001111, 10101010 */
        fits_write_col_byt(fp, 1, 1, 1, 3, &[0xF0, 0x0F, 0xAA], &mut status);
        /* 10000001, 11111111, 00000000 */
        fits_write_col_byt(fp, 2, 1, 1, 3, &[0x81, 0xFF, 0x00], &mut status);

        assert_eq!(status, 0, "bits table setup failed");
        f.unwrap()
    }

    /// Evaluate a boolean expression over every row.
    fn rows(f: &mut fitsfile, expr: &str) -> Vec<c_char> {
        let mut out = [0 as c_char; NROWS as usize];
        let mut n_good: c_long = 0;
        let mut status = 0;

        fits_find_rows(f, &cc(expr), 1, NROWS, &mut n_good, &mut out, &mut status);
        assert_eq!(status, 0, "evaluating `{expr}` failed");
        out.to_vec()
    }

    fn t() -> c_char {
        1
    }
    fn f_() -> c_char {
        0
    }

    #[test]
    fn bitand_over_rows() {
        let mut f = bits_table();
        /* 11110000 & 10000001 = 10000000 -> only row 1 matches 10000000 */
        assert_eq!(
            rows(&mut f, "(BITS & OTHER) == b10000000"),
            vec![t(), f_(), f_()]
        );
        /* 00001111 & 11111111 = 00001111 */
        assert_eq!(
            rows(&mut f, "(BITS & OTHER) == b00001111"),
            vec![f_(), t(), f_()]
        );
        /* 10101010 & 00000000 = 00000000 */
        assert_eq!(
            rows(&mut f, "(BITS & OTHER) == b00000000"),
            vec![f_(), f_(), t()]
        );
    }

    #[test]
    fn bitor_over_rows() {
        let mut f = bits_table();
        /* 11110000 | 10000001 = 11110001 */
        assert_eq!(
            rows(&mut f, "(BITS | OTHER) == b11110001"),
            vec![t(), f_(), f_()]
        );
        /* 00001111 | 11111111 = 11111111 */
        assert_eq!(
            rows(&mut f, "(BITS | OTHER) == b11111111"),
            vec![f_(), t(), f_()]
        );
        /* 10101010 | 00000000 = 10101010 */
        assert_eq!(
            rows(&mut f, "(BITS | OTHER) == b10101010"),
            vec![f_(), f_(), t()]
        );
    }

    #[test]
    fn bitnot_over_rows() {
        /* Do_Unary's Bits arm -- the one bit operator the corpus does not
        reach, in both its row and its constant form. */
        let mut f = bits_table();
        assert_eq!(rows(&mut f, "(!BITS) == b00001111"), vec![t(), f_(), f_()]);
        assert_eq!(rows(&mut f, "(!BITS) == b11110000"), vec![f_(), t(), f_()]);
        assert_eq!(rows(&mut f, "(!BITS) == b01010101"), vec![f_(), f_(), t()]);
        /* constant-folded form */
        assert_eq!(
            rows(&mut f, "(!b1010) == b0101"),
            vec![t(), t(), t()],
            "constant bitnot"
        );
    }

    #[test]
    fn bit_concatenation_over_rows() {
        let mut f = bits_table();
        /* '+' concatenates: 11110000 ++ 10000001 */
        assert_eq!(
            rows(&mut f, "(BITS + OTHER) == b1111000010000001"),
            vec![t(), f_(), f_()]
        );
        assert_eq!(
            rows(&mut f, "(BITS + BITS) == b1111000011110000"),
            vec![t(), f_(), f_()]
        );
    }

    #[test]
    fn bit_comparisons_over_rows() {
        let mut f = bits_table();
        /* bitcmp: equality, with 'x' matching anything */
        assert_eq!(rows(&mut f, "BITS == b1111xxxx"), vec![t(), f_(), f_()]);
        assert_eq!(rows(&mut f, "BITS == bxxxx0000"), vec![t(), f_(), f_()]);
        assert_eq!(rows(&mut f, "BITS != b11110000"), vec![f_(), t(), t()]);

        /* bitlgte: 0xF0=240, 0x0F=15, 0xAA=170 */
        assert_eq!(rows(&mut f, "BITS > b10101010"), vec![t(), f_(), f_()]);
        assert_eq!(rows(&mut f, "BITS >= b10101010"), vec![t(), f_(), t()]);
        assert_eq!(rows(&mut f, "BITS < b10101010"), vec![f_(), t(), f_()]);
        assert_eq!(rows(&mut f, "BITS <= b10101010"), vec![f_(), t(), t()]);
    }

    #[test]
    fn unequal_length_operands_are_left_padded() {
        /* bit_align: the shorter operand is padded with '0' on the left, so a
        4-bit literal compares as its 8-bit zero-extension. */
        let mut f = bits_table();
        assert_eq!(rows(&mut f, "BITS == b0000"), vec![f_(), f_(), f_()]);
        assert_eq!(rows(&mut f, "OTHER == b0"), vec![f_(), f_(), t()]);
        /* 00001111 vs 1111 -> equal once padded */
        assert_eq!(rows(&mut f, "BITS == b1111"), vec![f_(), t(), f_()]);
        /* and the padding also applies through & / | */
        assert_eq!(
            rows(&mut f, "(BITS & b1111) == b00000000"),
            vec![t(), f_(), f_()]
        );
    }

    #[test]
    fn constant_folding_matches_the_row_path() {
        /* Both operands constant takes Do_BinOp_bit's other branch, which
        writes the scalar Text arm rather than a row buffer. */
        let mut f = bits_table();
        assert_eq!(rows(&mut f, "(b101 & b110) == b100"), vec![t(), t(), t()]);
        assert_eq!(rows(&mut f, "(b101 | b110) == b111"), vec![t(), t(), t()]);
        assert_eq!(
            rows(&mut f, "(b101 + b110) == b101110"),
            vec![t(), t(), t()]
        );
        assert_eq!(rows(&mut f, "b101 > b100"), vec![t(), t(), t()]);
        assert_eq!(rows(&mut f, "b101 < b110"), vec![t(), t(), t()]);
    }

    #[test]
    fn row_offset_copies_bit_strings() {
        /* Do_Offset's Bits arm -- the strcpy that moves a whole row buffer. */
        let mut f = bits_table();
        /* row i reads row i+1: rows 1,2 see 00001111 and 10101010 */
        let got = rows(&mut f, "BITS{1} == b00001111");
        assert_eq!(got[0], t(), "row 1 offset by 1 should see row 2");
        let got = rows(&mut f, "BITS{1} == b10101010");
        assert_eq!(got[1], t(), "row 2 offset by 1 should see row 3");
    }
}
