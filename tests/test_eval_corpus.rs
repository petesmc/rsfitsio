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

    /// Whether a golden line records an answer a narrower C `long` cannot hold.
    ///
    /// Integer results travel through `c_long`, which is 64 bits on LP64 and 32
    /// on LLP64 (Windows). The golden was captured on LP64, so any line whose
    /// recorded value is outside the 32-bit range describes something Windows
    /// cannot reproduce -- correctly so, since C's `long` is genuinely narrower
    /// there. Those lines are skipped where a long is narrow and checked
    /// everywhere else.
    ///
    /// The test is on the *recorded value*, not on the expression, because the
    /// dependency does not have to come from a literal. `INTCOL ** INTCOL` is
    /// `10000000000` with no large literal anywhere in it, and
    /// `0o777777777777777777777` is a literal but not a decimal one -- an
    /// earlier version of this guard scanned the source text and missed both.
    fn depends_on_long_width(golden_line: &str) -> bool {
        size_of::<c_long>() < 8 && records_value_outside(golden_line, 32)
    }

    /// Whether the values a golden line records include an integer that does
    /// not fit in a signed integer of `bits` bits.
    ///
    /// Only the value section after `|` is considered: the metadata before it
    /// carries datatype and shape numbers that are the same on every platform.
    /// Anything with a decimal point is a double and travels in an `f64`
    /// regardless, so only bare integers count.
    fn records_value_outside(golden_line: &str, bits: u32) -> bool {
        let Some(values) = golden_line.split_once('|').map(|(_, v)| v) else {
            return false;
        };
        let lo = -(1i128 << (bits - 1));
        let hi = (1i128 << (bits - 1)) - 1;

        let b: Vec<char> = values.chars().collect();
        let mut i = 0;
        while i < b.len() {
            if !b[i].is_ascii_digit() {
                i += 1;
                continue;
            }
            let start = i;
            while i < b.len() && b[i].is_ascii_digit() {
                i += 1;
            }
            /* a decimal point on either side makes this a double's digits */
            let is_float = b.get(i) == Some(&'.') || (start > 0 && b[start - 1] == '.');
            let negative = start > 0 && b[start - 1] == '-';
            if is_float {
                continue;
            }
            let text: String = b[start..i].iter().collect();
            /* parse into i128 so a value far past i64 cannot wrap the check */
            match text.parse::<i128>() {
                Ok(v) => {
                    let v = if negative { -v } else { v };
                    if v < lo || v > hi {
                        return true;
                    }
                }
                /* longer than i128 can hold, so certainly outside */
                Err(_) => return true,
            }
        }
        false
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
        let integral = matches!(
            datatype,
            TLOGICAL | TBYTE | TSHORT | TINT | TLONG | TLONGLONG
        );
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

    /// Floor on the corpus size.
    ///
    /// The test reports a diff per line, so a *truncated* corpus would still be
    /// caught by the golden being longer — but if both files shrank together
    /// (a bad merge, a stray editor save) it would pass while testing far less
    /// than it claims. Raise this when the corpus grows meaningfully.
    const MIN_CORPUS: usize = 1800;

    /// The 32-bit-long skip list, checked here because CI's LLP64 machines are
    /// the ones that need it and this machine cannot run them.
    ///
    /// The positive cases are the golden lines Windows actually reported, not
    /// invented ones: two of them have no large literal in the expression at
    /// all, which is why the guard reads the recorded value instead.
    #[test]
    fn the_long_width_guard_picks_out_the_right_lines() {
        for g in [
            /* a literal past the range, saturated by atol */
            "9223372036854775807\tOK dt=41 nelem=-1 naxis=1 naxes=[1] | [9223372036854775807]",
            /* an octal literal, truncated by the cast in radix_const */
            "0o777777777777777777777\tOK dt=41 | [9223372036854775807, 9223372036854775807]",
            /* and a computed result, with no large literal anywhere */
            "INTCOL ** INTCOL\tOK dt=41 nelem=1 naxis=1 naxes=[1] | [823543, 0, 10000000000]",
            "INTCOL ^ INTCOL\tOK dt=41 nelem=1 naxis=1 naxes=[1] | [823543, 0, 10000000000]",
            "-9223372036854775807\tOK dt=41 | [-9223372036854775807]",
            /* longer than any integer type, so certainly outside */
            "big\tOK dt=41 | [123456789012345678901234567890123456789012]",
        ] {
            assert!(records_value_outside(g, 32), "should skip: {g}");
        }

        for g in [
            /* the boundary itself fits */
            "2147483647\tOK dt=41 nelem=-1 naxis=1 naxes=[1] | [2147483647, 2147483647]",
            "-2147483648\tOK dt=41 | [-2147483648]",
            "INTCOL\tOK dt=41 nelem=1 naxis=1 naxes=[1] | [7, -3, 10]",
            /* doubles travel in an f64 whatever a long is */
            "1e308\tOK dt=82 nelem=-1 naxis=1 naxes=[1] | [100000000000000000000.000000]",
            "FLOATCOL\tOK dt=82 nelem=1 naxis=1 naxes=[1] | [2.500000, 4.000000, 0.500000]",
            /* a parse error records no values at all */
            "INTCOL{'a'}\tERR 431",
            /* and a big number in the expression does not matter if the
            recorded answer is small */
            "9223372036854775807 > 0\tOK dt=14 nelem=-1 naxis=1 naxes=[1] | [1, 1, 1]",
        ] {
            assert!(!records_value_outside(g, 32), "should not skip: {g}");
        }

        /* on this platform the guard is inert whatever the line says */
        if size_of::<c_long>() >= 8 {
            assert!(!depends_on_long_width(
                "x\tOK dt=41 | [9223372036854775807]"
            ));
        }
    }

    #[test]
    fn eval_corpus_matches_golden() {
        let exprs = corpus_lines();
        assert!(
            exprs.len() >= MIN_CORPUS,
            "corpus has shrunk to {} expressions, expected at least {MIN_CORPUS}",
            exprs.len()
        );

        let mut status = 0;
        let mut f = create_corpus_table();
        let mut actual = String::new();
        for expr in exprs {
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
            let mut skipped = 0;
            let mut msg = String::new();
            for i in 0..a.len().max(g.len()) {
                let x = a.get(i).copied().unwrap_or("<missing>");
                let y = g.get(i).copied().unwrap_or("<missing>");
                if depends_on_long_width(y) {
                    skipped += 1;
                    continue;
                }
                if x != y {
                    diffs += 1;
                    if diffs <= 40 {
                        let _ =
                            writeln!(msg, "  line {}:\n    golden: {y}\n    actual: {x}", i + 1);
                    }
                }
            }
            if diffs == 0 {
                /* only the lines whose answer depends on the width of a long,
                which this platform cannot be held to */
                eprintln!("{skipped} corpus line(s) skipped: a C long is 32 bits here");
                return;
            }
            panic!("{diffs} corpus line(s) diverge from the golden file:\n{msg}");
        }
    }
}
