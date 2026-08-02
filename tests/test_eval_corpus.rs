//! Differential corpus for the expression parser.
//!
//! Each line of `fixtures/eval_corpus.txt` is parsed with `fits_test_expr` and,
//! when it parses and is deterministic, evaluated with `fits_calc_rows`. The
//! resulting one-line summary is compared against `fixtures/eval_corpus.golden`,
//! which was generated from the flex/bison parser before the `nom` rewrite.
//!
//! This is the acceptance test for the parser migration: the new front end must
//! reproduce the old one's accept/reject decision, inferred type and shape, and
//! evaluated values, for every line.
//!
//! Regenerate the golden file after an *intentional* language change with:
//!
//! ```text
//! UPDATE_EVAL_GOLDEN=1 cargo test --test test_eval_corpus
//! ```
//!
//! and then check the result against the real CFITSIO parser — see
//! `tests/oracle/README.md`.  To find which expression is responsible when the
//! evaluation engine crashes:
//!
//! ```text
//! CORPUS_TRACE=1 cargo test --test test_eval_corpus -- --nocapture
//! ```

#[cfg(test)]
mod tests {
    use rsfitsio::aliases::rust_api::*;
    use rsfitsio::c_types::{c_char, c_int, c_long};
    use rsfitsio::fitsio::{
        BINARY_TBL, BYTE_IMG, LONGLONG, TBYTE, TDOUBLE, TINT, TLOGICAL, TLONG, TLONGLONG, TSHORT,
        TSTRING, fitsfile,
    };
    use std::fmt::Write as _;

    const CORPUS: &str = include_str!("fixtures/eval_corpus.txt");
    const GOLDEN: &str = include_str!("fixtures/eval_corpus.golden");

    /// Number of rows in the corpus table.
    const NROWS: c_long = 3;
    /// Cap on how many vector elements are evaluated for one expression.
    const MAX_ELEM: c_long = 8;

    fn cc(s: &str) -> Vec<c_char> {
        let mut v: Vec<c_char> = s.bytes().map(|b| b as c_char).collect();
        v.push(0);
        v
    }

    fn as_bytes_mut<T>(s: &mut [T]) -> &mut [u8] {
        unsafe {
            core::slice::from_raw_parts_mut(s.as_mut_ptr().cast::<u8>(), core::mem::size_of_val(s))
        }
    }

    /// Build the corpus table: one column per value sort, plus a 2-D vector.
    ///
    /// Values are small, exact in binary floating point, and distinct across
    /// rows so that ordering mistakes show up.
    fn create_corpus_table() -> Box<fitsfile> {
        let mut status = 0;
        let mut f: Option<Box<fitsfile>> = None;
        fits_create_file(&mut f, &cc("mem://corpus.fits"), &mut status);
        fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);

        let ttype = [
            "INTCOL", "FLOATCOL", "STRCOL", "BOOLCOL", "DBLCOL", "VECCOL", "BITS", "MATRIX",
        ];
        let tform = ["1J", "1E", "10A", "1L", "1D", "3E", "8X", "6E"];
        let tt: Vec<Vec<c_char>> = ttype.iter().map(|s| cc(s)).collect();
        let tf: Vec<Vec<c_char>> = tform.iter().map(|s| cc(s)).collect();
        let tt_ref: Vec<Option<&[c_char]>> = tt.iter().map(|v| Some(v.as_slice())).collect();
        let tf_ref: Vec<&[c_char]> = tf.iter().map(|v| v.as_slice()).collect();
        fits_create_tbl(
            f.as_deref_mut().unwrap(),
            BINARY_TBL,
            NROWS as LONGLONG,
            ttype.len() as c_int,
            &tt_ref,
            &tf_ref,
            None,
            None,
            &mut status,
        );

        let fp = f.as_deref_mut().unwrap();
        fits_write_key_str(fp, &cc("TDIM8"), &cc("(2,3)"), None, &mut status);
        fits_write_key_lng(fp, &cc("INTKEY"), 42, None, &mut status);
        fits_write_key_dbl(fp, &cc("DBLKEY"), 2.5, 4, None, &mut status);
        fits_write_key_log(fp, &cc("LOGKEY"), 1, None, &mut status);
        fits_write_key_str(fp, &cc("STRKEY"), &cc("hello"), None, &mut status);

        fits_write_col_lng(fp, 1, 1, 1, 3, &[7, -3, 10], &mut status);
        fits_write_col_flt(fp, 2, 1, 1, 3, &[2.5, 4.0, 0.5], &mut status);
        for (i, s) in ["abc", "de", "fghij"].iter().enumerate() {
            let sv = cc(s);
            let arr: [&[c_char]; 1] = [sv.as_slice()];
            fits_write_col_str(fp, 3, (i + 1) as LONGLONG, 1, 1, &arr, &mut status);
        }
        fits_write_col_log(fp, 4, 1, 1, 3, &[1, 0, 1], &mut status);
        fits_write_col_dbl(fp, 5, 1, 1, 3, &[1.25, -2.5, 8.0], &mut status);
        fits_write_col_flt(
            fp,
            6,
            1,
            1,
            9,
            &[1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0],
            &mut status,
        );
        /* BITS: 0xF0, 0x0F, 0xAA across the three rows */
        fits_write_col_byt(fp, 7, 1, 1, 3, &[0xF0, 0x0F, 0xAA], &mut status);
        fits_write_col_flt(
            fp,
            8,
            1,
            1,
            18,
            &[
                1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 21.0, 22.0, 23.0,
                24.0, 25.0, 26.0,
            ],
            &mut status,
        );

        assert_eq!(status, 0, "corpus table setup failed");
        f.unwrap()
    }

    /// Expressions whose value is not reproducible run parse-only.
    fn is_nondeterministic(expr: &str) -> bool {
        expr.to_ascii_uppercase().contains("RANDOM")
    }

    /// Render one corpus line's behaviour as a single deterministic string.
    fn probe(f: &mut fitsfile, expr: &str) -> String {
        /* CORPUS_TRACE=1 ... -- --nocapture names each expression before it is
        evaluated, so a crash inside the engine can be pinned to a line. */
        if std::env::var_os("CORPUS_TRACE").is_some() {
            eprintln!("PROBE {expr}");
        }
        let mut datatype = 0;
        let mut nelem: c_long = 0;
        let mut naxis = 0;
        let mut naxes = [0 as c_long; 5];
        let mut status = 0;

        fits_test_expr(
            f,
            &cc(expr),
            5,
            &mut datatype,
            &mut nelem,
            &mut naxis,
            &mut naxes,
            &mut status,
        );
        if status != 0 {
            fits_clear_errmsg();
            return format!("ERR {status}");
        }

        let mut out = String::new();
        let _ = write!(out, "OK dt={datatype} nelem={nelem} naxis={naxis}");
        if naxis > 0 && (naxis as usize) <= naxes.len() {
            let _ = write!(out, " naxes={:?}", &naxes[..naxis as usize]);
        }

        /* String results carry no numeric value; the corpus exercises string
        semantics through the boolean comparison lines instead. */
        if datatype == TSTRING || is_nondeterministic(expr) {
            return out;
        }

        let n = nelem.clamp(1, MAX_ELEM);
        let mut anynul = 0;
        let mut st = 0;

        /* Integer-valued expressions are read back as LONGLONG so that large
        magnitudes stay exact; f64 cannot distinguish i64::MAX from 2^63. */
        let integral = matches!(datatype, TLOGICAL | TBYTE | TSHORT | TINT | TLONG | TLONGLONG);
        let rendered = if integral {
            let mut results = vec![0 as LONGLONG; (n * NROWS) as usize];
            fits_calc_rows(
                f,
                TLONGLONG,
                &cc(expr),
                1,
                n * NROWS,
                core::ptr::null(),
                as_bytes_mut(&mut results),
                &mut anynul,
                &mut st,
            );
            results
                .iter()
                .map(|v| v.to_string())
                .collect::<Vec<_>>()
                .join(", ")
        } else {
            let mut results = vec![0.0f64; (n * NROWS) as usize];
            fits_calc_rows(
                f,
                TDOUBLE,
                &cc(expr),
                1,
                n * NROWS,
                core::ptr::null(),
                as_bytes_mut(&mut results),
                &mut anynul,
                &mut st,
            );
            results
                .iter()
                .map(|v| format!("{v:.6}"))
                .collect::<Vec<_>>()
                .join(", ")
        };

        if st != 0 {
            fits_clear_errmsg();
            let _ = write!(out, " | EVALERR {st}");
            return out;
        }
        let _ = write!(out, " | [{rendered}]");
        if anynul != 0 {
            out.push_str(" (null)");
        }
        out
    }

    /// A corpus line is a comment when it is `#` alone or starts with `# `.
    /// Bare `#NAME` lines are expressions (the `#KEYWORD` constant form).
    fn is_comment(line: &str) -> bool {
        line == "#" || line.starts_with("# ")
    }

    fn corpus_lines() -> Vec<&'static str> {
        CORPUS
            .lines()
            .map(str::trim_end)
            .filter(|l| !l.is_empty() && !is_comment(l))
            .collect()
    }

    #[test]
    fn eval_corpus_matches_golden() {
        let mut status = 0;
        let mut f = create_corpus_table();
        let mut actual = String::new();
        for expr in corpus_lines() {
            let _ = writeln!(actual, "{expr}\t{}", probe(&mut f, expr));
        }
        fits_close_file(f, &mut status);
        assert_eq!(status, 0);

        if std::env::var_os("UPDATE_EVAL_GOLDEN").is_some() {
            let path = concat!(
                env!("CARGO_MANIFEST_DIR"),
                "/tests/fixtures/eval_corpus.golden"
            );
            std::fs::write(path, &actual).expect("write golden");
            eprintln!("wrote {path}");
            return;
        }

        if actual != GOLDEN {
            let a: Vec<&str> = actual.lines().collect();
            let g: Vec<&str> = GOLDEN.lines().collect();
            let mut diffs = 0;
            let mut msg = String::new();
            for i in 0..a.len().max(g.len()) {
                let x = a.get(i).copied().unwrap_or("<missing>");
                let y = g.get(i).copied().unwrap_or("<missing>");
                if x != y {
                    diffs += 1;
                    if diffs <= 40 {
                        let _ = writeln!(msg, "  line {}:\n    golden: {y}\n    actual: {x}", i + 1);
                    }
                }
            }
            panic!("{diffs} corpus line(s) diverge from the golden file:\n{msg}");
        }
    }
}
