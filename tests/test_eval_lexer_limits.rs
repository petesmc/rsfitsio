//! Over-long-token paths in the expression lexer (`src/eval_l.rs`).
//!
//! Every lexer rule that builds a bit string, hex string, quoted string,
//! `#KEYWORD` reference or function name has an "exceeds maximum length" arm
//! that truncates the offending token into a 100-byte error card. Those arms
//! are the ones the `eval_corpus` golden cannot reach -- it only carries
//! well-formed expressions -- so they are pinned here instead.
//!
//! Each case asserts the parse is rejected with `PARSE_SYNTAX_ERR` rather than
//! overrunning a buffer, which is what the raw `strncat`/`strncpy` calls these
//! rules used to make possible.

mod common;

#[cfg(test)]
mod tests {
    use crate::common::with_temp_file;
    use rsfitsio::aliases::rust_api::*;
    use rsfitsio::c_types::{c_char, c_int, c_long};
    use rsfitsio::fitsio::{BINARY_TBL, BYTE_IMG, PARSE_SYNTAX_ERR, fitsfile};

    /// `MAX_STRLEN` from `eval_defs.rs` -- the lexer's per-token buffer size.
    const MAX_STRLEN: usize = 256;

    fn cc(s: &str) -> Vec<c_char> {
        let mut v: Vec<c_char> = s.bytes().map(|b| b as c_char).collect();
        v.push(0);
        v
    }

    /// Parse `expr` against a one-column table and return the status.
    fn parse_status(f: &mut fitsfile, expr: &str) -> c_int {
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
        }
        status
    }

    fn with_table<F: Fn(&mut fitsfile)>(body: F) {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = cc(filename);
            let ttype = [Some(&cc("X")[..])];
            let tform = [&cc("1J")[..]];

            let mut fptr: Option<Box<fitsfile>> = None;
            fits_create_file(&mut fptr, &name, &mut status);
            fits_write_imghdr(fptr.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
            fits_create_tbl(
                fptr.as_deref_mut().unwrap(),
                BINARY_TBL,
                3,
                1,
                &ttype,
                &tform,
                None,
                None,
                &mut status,
            );
            assert_eq!(status, 0, "table setup failed");

            body(fptr.as_deref_mut().unwrap());

            fits_close_file(fptr.take().unwrap(), &mut status);
        });
    }

    /// The bit buffer the octal and hex rules expand into holds 255 bits.
    /// (`PARSER_SPEC.md` §3, rules 3 and 4.)
    const BIT_CAP: usize = 255;

    #[test]
    fn short_tokens_still_parse() {
        /* Guards the tests below: the same shapes, under the limit, are fine. */
        with_table(|f| {
            assert_eq!(parse_status(f, "b1010"), 0, "short bit string");
            assert_eq!(parse_status(f, "o77"), 0, "short octal bit string");
            assert_eq!(parse_status(f, "h1F"), 0, "short hex bit string");
            assert_eq!(parse_status(f, "'hello'"), 0, "short quoted string");
            assert_eq!(parse_status(f, "abs(X)"), 0, "short function name");
        });
    }

    #[test]
    fn overlong_bit_string_is_rejected() {
        /* Rule 2: len >= MAX_STRLEN takes the truncating error arm. */
        with_table(|f| {
            assert_eq!(
                parse_status(f, &format!("b{}", "1".repeat(MAX_STRLEN - 1))),
                0,
                "one under the limit still parses"
            );
            let expr = format!("b{}", "1".repeat(MAX_STRLEN));
            assert_eq!(parse_status(f, &expr), PARSE_SYNTAX_ERR);
        });
    }

    #[test]
    fn overlong_hex_bit_string_is_rejected() {
        /* Rule 4: `h…`, each hex digit expanding to four bits. A token far
        under MAX_STRLEN characters still overflows the 255-bit buffer. */
        with_table(|f| {
            let ok = BIT_CAP / 4; // 63 digits -> 252 bits
            assert_eq!(parse_status(f, &format!("h{}", "f".repeat(ok))), 0);
            let expr = format!("h{}", "f".repeat(ok + 1)); // 256 bits
            assert_eq!(parse_status(f, &expr), PARSE_SYNTAX_ERR);
        });
    }

    #[test]
    fn overlong_octal_bit_string_is_rejected() {
        /* Rule 3: `o…`, each octal digit expanding to three bits. */
        with_table(|f| {
            let ok = BIT_CAP / 3; // 85 digits -> 255 bits, exactly the cap
            assert_eq!(parse_status(f, &format!("o{}", "7".repeat(ok))), 0);
            let expr = format!("o{}", "7".repeat(ok + 1)); // 258 bits
            assert_eq!(parse_status(f, &expr), PARSE_SYNTAX_ERR);
        });
    }

    #[test]
    fn overlong_quoted_string_is_rejected() {
        /* Rule 12: len-2 >= MAX_STRLEN, the arm that truncates to 20 chars. */
        with_table(|f| {
            let expr = format!("'{}'", "a".repeat(MAX_STRLEN + 10));
            assert_eq!(parse_status(f, &expr), PARSE_SYNTAX_ERR);
        });
    }

    #[test]
    fn overlong_keyword_reference_is_rejected() {
        /* The `#KEYWORD` rule, which copies from yytext+2 after writing '#'.
        The truncated name is then looked up in the header, so this comes back
        KEY_NO_EXIST rather than a syntax error -- either way it is rejected
        without overrunning the 256-byte token buffer. */
        with_table(|f| {
            let expr = format!("#{}", "A".repeat(MAX_STRLEN + 10));
            assert_ne!(parse_status(f, &expr), 0);
        });
    }

    #[test]
    fn overlong_function_name_is_rejected() {
        /* Rule 14: the function-name rule, which reads yytext by length. */
        with_table(|f| {
            let expr = format!("{}(X)", "z".repeat(MAX_STRLEN + 10));
            assert_eq!(parse_status(f, &expr), PARSE_SYNTAX_ERR);
        });
    }

    #[test]
    fn token_at_exactly_the_limit_does_not_overrun() {
        /* The boundary between the accepting and the truncating arm -- the
        index the raw strncpy calls used to write one past. */
        with_table(|f| {
            for n in [MAX_STRLEN - 2, MAX_STRLEN - 1, MAX_STRLEN, MAX_STRLEN + 1] {
                let expr = format!("'{}'", "a".repeat(n));
                /* Either outcome is acceptable; not panicking is the point. */
                let st = parse_status(f, &expr);
                assert!(st == 0 || st == PARSE_SYNTAX_ERR, "n = {n}, status = {st}");
            }
        });
    }
}
