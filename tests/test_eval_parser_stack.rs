//! Growth of the bison parser's state and value stacks (`fits_parser_yyparse`).
//!
//! The two stacks were a pair of `YYINITDEPTH` C arrays that, on overflow, were
//! replaced by a single `alloca`-turned-`malloc` block carved in two with
//! alignment arithmetic. They are now owned `Vec`s that `resize`.
//!
//! `YYINITDEPTH` is 100 and `YYMAXDEPTH` is 10000, while the deepest expression
//! in the `eval_corpus` golden nests 20 parentheses -- so nothing there reaches
//! the growth path at all. These do.

mod common;

#[cfg(test)]
mod tests {
    use rsfitsio::aliases::rust_api::*;
    use rsfitsio::c_types::{c_char, c_int, c_long};
    use rsfitsio::fitsio::{BINARY_TBL, BYTE_IMG, fitsfile};

    /// `YYINITDEPTH` from eval_y.rs: the depth the stacks start at.
    const YYINITDEPTH: usize = 100;

    fn cc(s: &str) -> Vec<c_char> {
        let mut v: Vec<c_char> = s.bytes().map(|b| b as c_char).collect();
        v.push(0);
        v
    }

    fn table() -> Box<fitsfile> {
        let mut status = 0;
        let mut f: Option<Box<fitsfile>> = None;
        fits_create_file(&mut f, &cc("mem://stack.fits"), &mut status);
        fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);

        let tt = [cc("X")];
        let tf = [cc("1J")];
        let tt_ref: Vec<Option<&[c_char]>> = tt.iter().map(|v| Some(v.as_slice())).collect();
        let tf_ref: Vec<&[c_char]> = tf.iter().map(|v| v.as_slice()).collect();
        fits_create_tbl(
            f.as_deref_mut().unwrap(),
            BINARY_TBL,
            2,
            1,
            &tt_ref,
            &tf_ref,
            None,
            None,
            &mut status,
        );
        fits_write_col_lng(f.as_deref_mut().unwrap(), 1, 1, 1, 2, &[3, 4], &mut status);
        assert_eq!(status, 0, "table setup failed");
        f.unwrap()
    }

    /// Parse `expr` and return (status, evaluated rows).
    fn eval(f: &mut fitsfile, expr: &str) -> (c_int, Vec<c_char>) {
        let mut out = [0 as c_char; 2];
        let mut n_good: c_long = 0;
        let mut status = 0;
        fits_find_rows(f, &cc(expr), 1, 2, &mut n_good, &mut out, &mut status);
        if status != 0 {
            fits_clear_errmsg();
        }
        (status, out.to_vec())
    }

    /// `((((…X > 0…))))` with `n` parentheses on each side.
    fn nested(n: usize) -> String {
        format!("{}X > 0{}", "(".repeat(n), ")".repeat(n))
    }

    #[test]
    fn shallow_nesting_stays_on_the_initial_stack() {
        let mut f = table();
        let (st, rows) = eval(&mut f, &nested(10));
        assert_eq!(st, 0);
        assert_eq!(rows, vec![1, 1]);
    }

    #[test]
    fn nesting_past_the_initial_depth_grows_the_stacks() {
        /* Each parenthesis pushes, so this is comfortably past YYINITDEPTH and
        forces at least one resize of both stacks. */
        let mut f = table();
        for depth in [YYINITDEPTH + 5, YYINITDEPTH * 2, YYINITDEPTH * 4] {
            let (st, rows) = eval(&mut f, &nested(depth));
            assert_eq!(st, 0, "depth {depth} failed to parse");
            assert_eq!(rows, vec![1, 1], "depth {depth} gave the wrong answer");
        }
    }

    #[test]
    fn values_survive_the_stack_growth() {
        /* The resize has to preserve what is already on the value stack -- the
        C did that with a memcpy. If it did not, the operands of the outer
        expression would be lost, so check a value that depends on tokens
        pushed before the growth and reduced after it. */
        let mut f = table();
        let deep = nested(YYINITDEPTH + 20);
        assert_eq!(eval(&mut f, &format!("({deep}) && X == 3")).1, vec![1, 0]);
        assert_eq!(eval(&mut f, &format!("({deep}) && X == 4")).1, vec![0, 1]);
        /* a long left-associative chain also walks the value stack deep */
        let sum = std::iter::repeat_n("X", 200).collect::<Vec<_>>().join("+");
        assert_eq!(eval(&mut f, &format!("({sum}) == X * 200")).1, vec![1, 1]);
    }

    #[test]
    fn nesting_past_the_maximum_is_rejected_not_crashed() {
        /* Beyond YYMAXDEPTH the parser takes yyexhaustedlab. The point is that
        it reports an error rather than overrunning or aborting. */
        let mut f = table();
        let (st, _) = eval(&mut f, &nested(20_000));
        assert_ne!(st, 0, "an over-deep expression should be rejected");
    }
}
