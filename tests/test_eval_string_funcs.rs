//! Per-row values of `STRMID` and `STRSTR`, and of string concatenation and
//! comparison, in `src/eval_y.rs`.
//!
//! `cstrmid` now takes slices instead of raw pointers, and the `STRMID_FCT` /
//! `STRPOS_FCT` row loops reach their subject through `lval::str_row` rather
//! than a `char **` row-pointer array. The `eval_corpus` golden covers the
//! shapes of these expressions, but its probe skips the *value* of a
//! `TSTRING`-typed result, so `STRMID`'s output and the string operators are
//! pinned here through boolean comparisons, one per row.
//!
//! The rows are `SA = "abcdefg", "xy", "mnop"` in an `8A` column and
//! `SB = "cd", "xy", "zz"` in a `4A` one. Two things about `STRMID` follow from
//! the engine rather than from the strings, and the expectations below depend
//! on both:
//!
//!  * the source length is the *column width*, not `strlen` of the row, and
//!    the reader zero-fills the rest of the field -- so a window past the end
//!    of a short row yields the empty string rather than an error;
//!  * a window running off the end copies only as far as the field goes.
//!
//! Column names avoid `T` and `F`, which the grammar reads as the FITS logical
//! constants rather than as column references.

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

    /// Two string columns, three rows.
    fn str_table() -> Box<fitsfile> {
        let mut status = 0;
        let mut f: Option<Box<fitsfile>> = None;
        fits_create_file(&mut f, &cc("mem://strs.fits"), &mut status);
        fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);

        let tt = [cc("SA"), cc("SB")];
        let tf = [cc("8A"), cc("4A")];
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
        for (i, s) in ["abcdefg", "xy", "mnop"].iter().enumerate() {
            let sv = cc(s);
            let arr: [&[c_char]; 1] = [sv.as_slice()];
            fits_write_col_str(fp, 1, (i + 1) as LONGLONG, 1, 1, &arr, &mut status);
        }
        for (i, s) in ["cd", "xy", "zz"].iter().enumerate() {
            let sv = cc(s);
            let arr: [&[c_char]; 1] = [sv.as_slice()];
            fits_write_col_str(fp, 2, (i + 1) as LONGLONG, 1, 1, &arr, &mut status);
        }

        assert_eq!(status, 0, "string table setup failed");
        f.unwrap()
    }

    fn rows(f: &mut fitsfile, expr: &str) -> Vec<c_char> {
        let mut out = [0 as c_char; NROWS as usize];
        let mut n_good: c_long = 0;
        let mut status = 0;

        fits_find_rows(f, &cc(expr), 1, NROWS, &mut n_good, &mut out, &mut status);
        assert_eq!(status, 0, "evaluating `{expr}` failed");
        out.to_vec()
    }

    const T: c_char = 1;
    const F: c_char = 0;

    #[test]
    fn strmid_over_row_buffers() {
        /* cstrmid with a per-row subject. */
        let mut f = str_table();
        assert_eq!(rows(&mut f, "STRMID(SA,1,3) == 'abc'"), vec![T, F, F]);
        assert_eq!(rows(&mut f, "STRMID(SA,1,3) == 'xy'"), vec![F, T, F]);
        assert_eq!(rows(&mut f, "STRMID(SA,1,3) == 'mno'"), vec![F, F, T]);

        /* A window past the end of a short row lands in the field's zero fill,
        so rows 2 ("xy") and 3 ("mnop") both come back empty. */
        assert_eq!(rows(&mut f, "STRMID(SA,5,2) == 'ef'"), vec![T, F, F]);
        assert_eq!(rows(&mut f, "STRMID(SA,5,2) == ''"), vec![F, T, T]);

        /* A window running off the end copies only to the end of the field. */
        assert_eq!(rows(&mut f, "STRMID(SA,3,9) == 'cdefg'"), vec![T, F, F]);
        assert_eq!(rows(&mut f, "STRMID(SA,3,9) == 'op'"), vec![F, F, T]);
    }

    #[test]
    fn strmid_over_a_constant_subject() {
        /* The strconst arm, where the subject is copied out of the node once
        before the row loop. */
        let mut f = str_table();
        assert_eq!(rows(&mut f, "STRMID('abcdef',2,3) == 'bcd'"), vec![T, T, T]);
        assert_eq!(
            rows(&mut f, "STRMID('abcdef',1,6) == 'abcdef'"),
            vec![T, T, T]
        );
        /* pos == 0 asks for a null string. The constant path returns it as an
        empty string; only the per-row path additionally flags it undefined. */
        assert_eq!(rows(&mut f, "STRMID('abcdef',0,3) == ''"), vec![T, T, T]);
    }

    #[test]
    fn strstr_over_row_buffers() {
        /* STRPOS_FCT: 1-based index. Per-row subject, constant needle. */
        let mut f = str_table();
        assert_eq!(rows(&mut f, "STRSTR(SA,'cd') == 3"), vec![T, F, F]);
        assert_eq!(rows(&mut f, "STRSTR(SA,'y') == 2"), vec![F, T, F]);
        assert_eq!(rows(&mut f, "STRSTR(SA,'op') == 3"), vec![F, F, T]);
        /* An absent needle is *undefined* on the per-row path. */
        assert_eq!(
            rows(&mut f, "DEFNULL(STRSTR(SA,'zz'),-1) == -1"),
            vec![T, T, T]
        );
    }

    #[test]
    fn strstr_with_both_operands_per_row() {
        /* Subject and needle both come from row buffers: "cd" in "abcdefg" at
        3, "xy" in "xy" at 1, and "zz" absent from "mnop". */
        let mut f = str_table();
        assert_eq!(rows(&mut f, "STRSTR(SA,SB) == 3"), vec![T, F, F]);
        assert_eq!(rows(&mut f, "STRSTR(SA,SB) == 1"), vec![F, T, F]);
        assert_eq!(
            rows(&mut f, "DEFNULL(STRSTR(SA,SB),-1) == -1"),
            vec![F, F, T]
        );
    }

    #[test]
    fn strstr_with_both_operands_constant() {
        /* The constant-folded branch. Note the asymmetry with the row path
        above: when the needle is absent this returns 0 rather than marking the
        result undefined, so DEFNULL has nothing to replace. That is the C's
        behaviour and is preserved here deliberately. */
        let mut f = str_table();
        assert_eq!(rows(&mut f, "STRSTR('abc','b') == 2"), vec![T, T, T]);
        assert_eq!(rows(&mut f, "STRSTR('abc','a') == 1"), vec![T, T, T]);
        assert_eq!(rows(&mut f, "STRSTR('abc','z') == 0"), vec![T, T, T]);
        assert_eq!(
            rows(&mut f, "DEFNULL(STRSTR('abc','z'),-1) == -1"),
            vec![F, F, F]
        );
    }

    #[test]
    fn string_concatenation_and_comparison() {
        /* Do_BinOp_str: '+' concatenates, and the comparisons go through the
        first-character fast path and then strcmp_safe. */
        let mut f = str_table();
        assert_eq!(rows(&mut f, "(SA + SB) == 'abcdefgcd'"), vec![T, F, F]);
        assert_eq!(rows(&mut f, "(SA + SB) == 'xyxy'"), vec![F, T, F]);
        assert_eq!(rows(&mut f, "(SB + SB) == 'zzzz'"), vec![F, F, T]);

        assert_eq!(rows(&mut f, "SA == SB"), vec![F, T, F]);
        assert_eq!(rows(&mut f, "SA != SB"), vec![T, F, T]);
        assert_eq!(rows(&mut f, "SA < SB"), vec![T, F, T]);
        assert_eq!(rows(&mut f, "SA > SB"), vec![F, F, F]);
        assert_eq!(rows(&mut f, "SA <= SB"), vec![T, T, T]);
    }

    #[test]
    fn row_offset_copies_strings() {
        /* Do_Offset's String arm -- the strcpy that moves a whole row. */
        let mut f = str_table();
        let got = rows(&mut f, "SA{1} == 'xy'");
        assert_eq!(got[0], T, "row 1 offset by 1 should see row 2");
        let got = rows(&mut f, "SA{1} == 'mnop'");
        assert_eq!(got[1], T, "row 2 offset by 1 should see row 3");
    }
}
