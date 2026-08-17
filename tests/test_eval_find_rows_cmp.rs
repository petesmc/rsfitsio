//! `fits_find_rows_cmp` (`fffrwc`) over a compressed-housekeeping table.
//!
//! The table is one row per (time, parameter, value) sample; the call collapses
//! it to one flag per distinct time. `ffiprs` is told the file is compressed, so
//! every identifier in the expression is looked up as a *parameter name* in the
//! parameter column rather than as a table column.
//!
//! This function had no coverage at all, which is how the two defects below
//! survived. Expected values are hand-derived from the rows written below.

mod common;

#[cfg(test)]
mod tests {
    use rsfitsio::aliases::rust_api::*;
    use rsfitsio::c_types::{c_char, c_int, c_long};
    use rsfitsio::fitsio::{BAD_COL_NUM, BINARY_TBL, BYTE_IMG, LONGLONG, fitsfile};

    fn cc(s: &str) -> Vec<c_char> {
        let mut v: Vec<c_char> = s.bytes().map(|b| b as c_char).collect();
        v.push(0);
        v
    }

    /// Four samples over two distinct times, two parameters each.
    ///
    /// | TIME | PARAM | VALUE |
    /// |------|-------|-------|
    /// | 1.0  | TEMP  |  3.0  |
    /// | 1.0  | PRES  | 10.0  |
    /// | 2.0  | TEMP  |  7.0  |
    /// | 2.0  | PRES  | 20.0  |
    fn hk_table() -> Box<fitsfile> {
        let mut status = 0;
        let mut f: Option<Box<fitsfile>> = None;
        fits_create_file(&mut f, &cc("mem://hk.fits"), &mut status);
        fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);

        let tt = [cc("TIME"), cc("PARAM"), cc("VALUE")];
        let tf = [cc("1D"), cc("8A"), cc("1D")];
        let ttr: Vec<Option<&[c_char]>> = tt.iter().map(|v| Some(v.as_slice())).collect();
        let tfr: Vec<&[c_char]> = tf.iter().map(|v| v.as_slice()).collect();
        fits_create_tbl(
            f.as_deref_mut().unwrap(),
            BINARY_TBL,
            4,
            3,
            &ttr,
            &tfr,
            None,
            None,
            &mut status,
        );

        let fp = f.as_deref_mut().unwrap();
        fits_write_col_dbl(fp, 1, 1, 1, 4, &[1.0, 1.0, 2.0, 2.0], &mut status);
        for (i, s) in ["TEMP", "PRES", "TEMP", "PRES"].iter().enumerate() {
            let sv = cc(s);
            let arr: [&[c_char]; 1] = [sv.as_slice()];
            fits_write_col_str(fp, 2, (i + 1) as LONGLONG, 1, 1, &arr, &mut status);
        }
        fits_write_col_dbl(fp, 3, 1, 1, 4, &[3.0, 10.0, 7.0, 20.0], &mut status);

        assert_eq!(status, 0, "hk table setup failed");
        f.unwrap()
    }

    /// Returns (status, times, flags).
    fn find_rows(f: &mut fitsfile, expr: &str) -> (c_int, Vec<f64>, Vec<c_char>) {
        let mut times = [0.0f64; 2];
        let mut flags = [0 as c_char; 2];
        let mut status = 0;

        fits_find_rows_cmp(
            f,
            &cc(expr),
            &cc("TIME"),
            &cc("PARAM"),
            &cc("VALUE"),
            2 as c_long,
            &mut times,
            &mut flags,
            &mut status,
        );
        if status != 0 {
            fits_clear_errmsg();
        }
        (status, times.to_vec(), flags.to_vec())
    }

    #[test]
    fn constant_expression_flags_every_time() {
        /* A constant expression takes the `nelem < 0` path, which zeroes nCols
        and so allocates no parameter arrays at all. The distinct times are
        still collected. */
        let mut f = hk_table();

        let (st, times, flags) = find_rows(&mut f, "1 > 0");
        assert_eq!(st, 0);
        assert_eq!(times, vec![1.0, 2.0]);
        assert_eq!(flags, vec![1, 1]);

        let (st, times, flags) = find_rows(&mut f, "0 > 1");
        assert_eq!(st, 0);
        assert_eq!(times, vec![1.0, 2.0]);
        assert_eq!(flags, vec![0, 0]);
    }

    #[test]
    fn distinct_times_are_collected_not_the_raw_rows() {
        /* Four rows, two distinct times: the output is per time, not per row. */
        let mut f = hk_table();
        let (st, times, _) = find_rows(&mut f, "1 > 0");
        assert_eq!(st, 0);
        assert_eq!(times.len(), 2);
        assert_eq!(times, vec![1.0, 2.0]);
    }

    #[test]
    fn parameter_reference_is_rejected_upstream_bug() {
        /* NOTE (upstream bug): an expression that names a parameter cannot
        work. `ffiprs` parses the expression *before* `fits_get_colnum` fills in
        lParse.valCol, and in compressed mode the parser resolves every
        identifier to that column number -- still 0 at parse time -- so
        fits_get_coltype rejects it with BAD_COL_NUM.

        cfitsio's eval_f.c has the calls in exactly this order
        (`ffiprs(...)` then the three `fits_get_colnum` calls, with
        `ParseData lParse = {0}`), so the C fails the same way. Reproduced here
        rather than fixed, and written up in notes/CFITSIO_BUGS_EVAL.md.

        The practical consequence is that only the constant path above is
        reachable, which is why the parameter-array allocation underneath was
        able to stay broken unnoticed. */
        let mut f = hk_table();
        for expr in ["TEMP > 5", "PRES > 15", "TEMP + PRES > 0"] {
            let (st, _, _) = find_rows(&mut f, expr);
            assert_eq!(st, BAD_COL_NUM, "`{expr}` should hit the upstream defect");
        }
    }

    #[test]
    fn a_missing_column_name_is_reported() {
        let mut f = hk_table();
        let mut times = [0.0f64; 2];
        let mut flags = [0 as c_char; 2];
        let mut status = 0;
        fits_find_rows_cmp(
            &mut f,
            &cc("1 > 0"),
            &cc("NOSUCH"),
            &cc("PARAM"),
            &cc("VALUE"),
            2 as c_long,
            &mut times,
            &mut flags,
            &mut status,
        );
        assert_ne!(status, 0, "an unknown time column should be an error");
        fits_clear_errmsg();
    }

    #[test]
    fn a_non_logical_expression_is_rejected() {
        let mut f = hk_table();
        let (st, _, _) = find_rows(&mut f, "1 + 1");
        assert_ne!(st, 0, "a non-logical expression should be rejected");
    }
}
