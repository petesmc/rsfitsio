# Pending rsfitsio fixes: expression lexer defects

> **DONE.** Both fixes landed in rsfitsio in commit `7a97251`, together with
> the differential corpus. `src/eval_l.rs` has since been deleted entirely by
> the `nom` migration (see [`PARSER_MIGRATION.md`](PARSER_MIGRATION.md)), so the
> file and line references below are historical. Kept as a record of the
> analysis; nothing here is outstanding.

Two defects in the transpiled expression lexer (`src/eval_l.rs`), inherited
verbatim from CFITSIO's `eval.l`. Both were fixed in `~/code/cfitsio` commit
`47359ca` and submitted upstream as [cfitsio#153](https://github.com/HEASARC/cfitsio/pull/153), which also adds C
regression tests to `tests/test_eval.c`. This document ported those fixes to
rsfitsio.

Background and verified behaviour: [`PARSER_SPEC.md`](PARSER_SPEC.md) §6.1, §6.2.

Both edits are small and localised. Line numbers are deliberately omitted —
use the quoted search anchors, since `src/eval_l.rs` may have moved.

---

## Bug 1 — Uppercase hex literals evaluate to garbage

### Symptom

| Expression | rsfitsio today | Correct |
|---|---|---|
| `0x1f` | 31 | 31 |
| `0x1F` | **−1** | 31 |
| `0xff` | 255 | 255 |
| `0xFF` | **−1** | 255 |
| `0xA` | **−22** | 10 |
| `0xABCDEF` | garbage | 11259375 |

### Cause

The `{hexconst}` action converts non-digits with `*p - 'a' + 10` and never
folds case, so `'F'` (70) yields `70 − 97 + 10 = −17`. The negative value is
then OR-ed into the accumulator, corrupting all higher bits.

### Fix — step 1: add `tolower` to `src/wrappers.rs`

There is a `toupper`/`islower` pair but no `tolower`/`isupper`. Add them
immediately after `islower`, matching the existing style:

```rust
pub(crate) fn tolower(c: c_char) -> c_char {
    if isupper(c) {
        return c | 0x20;
    }
    c
}

pub(crate) fn isupper(c: c_char) -> bool {
    (c >= 65) && (c <= 90)
}
```

### Fix — step 2: use it in `src/eval_l.rs`

`src/eval_l.rs` has two `use crate::wrappers::{…}` groups near the top: one
supplying `isdigit_safe`, and one supplying `toupper`. Add `tolower` to the
latter.

Then find the `{hexconst}` action, which is the match arm `7 => {` containing
`constval_1`, and change the digit conversion:

```rust
                                while *p_1 != 0 {
                                    let v: c_int = if isdigit_safe(*p_1) {
                                        c_int::from(*p_1) - '0' as i32
                                    } else {
-                                       c_int::from(*p_1) - 'a' as i32 + 10 as c_int
+                                       c_int::from(tolower(*p_1)) - 'a' as i32 + 10 as c_int
                                    };
                                    constval_1 = constval_1 << 4 as c_int | c_long::from(v);
                                    p_1 = p_1.offset(1);
                                }
```

This mirrors the upstream C exactly:

```c
int c = (unsigned char)*p;
int v = (isdigit(c) ? (c - '0') : (tolower(c) - 'a' + 10));
```

> The upstream cast to `unsigned char` guards against `char` being signed for
> bytes ≥ 128. It is not needed here: the `{hexconst}` pattern only ever
> matches `[0-9a-fA-F]`, all of which are < 128, so `c_char` (i8 on Linux)
> is never negative in this loop. Keep the `c_int::from(...)` conversions as
> they are.

### Do not touch

The `{bitconst}` (`5 =>`, `constval`) and `{octconst}` (`6 =>`, `constval_0`)
arms are correct — their patterns admit no letters. Leave them alone.

---

## Bug 2 — A bare `.` lexes as the double `0.0`

### Symptom

`.` and `..` are accepted as valid expressions evaluating to `0.0` /
two zero constants, instead of being syntax errors. `INTCOL + .` parses.

### Cause

The `{real}` pattern's third alternative made the integer part optional:

```
real  ([0-9]*"."[0-9]+)|([0-9]*"."*[0-9]+[eEdD][+-]?[0-9]+)|([0-9]*".")
                                                              ^^^^^ zero digits allowed
```

Upstream changed the last alternative to `([0-9]+".")`.

### Fix — one array element

rsfitsio does not contain the regex; it contains the compiled flex DFA. The
good news: **regenerating `eval_l.c` from the corrected `eval.l` with flex
2.6.4 changes exactly one table entry and nothing else.** Verified by diffing
the regenerated file — `yy_base`, `yy_def`, `yy_nxt`, `yy_chk`, `yy_ec` and
`yy_meta` are all byte-identical; only `yy_accept[15]` moves.

In `src/eval_l.rs`, in `static YY_ACCEPT: [flex_int16_t; 174]`, change
**element index 15** from `10` to `29`:

```rust
 static YY_ACCEPT: [flex_int16_t; 174] = [
-    0, 0, 0, 31, 29, 1, 28, 18, 29, 29, 29, 29, 29, 29, 29, 10, 8, 8, 24, 29, 23, 13, 13, 13, 13,
+    0, 0, 0, 31, 29, 1, 28, 18, 29, 29, 29, 29, 29, 29, 29, 29, 8, 8, 24, 29, 23, 13, 13, 13, 13,
     9, 13, 13, 13, 13, 13, 17, 13, 13, 13, 13, 13, 13, 13, 29, 1, 22, 0, 12, 0, 11, 0, 13, 20, 0,
     ...
```

Counting from `[0] = 0`, the elements are:

```
[0]=0 [1]=0 [2]=0 [3]=31 [4]=29 [5]=1 [6]=28 [7]=18 [8]=29 [9]=29 [10]=29
[11]=29 [12]=29 [13]=29 [14]=29 [15]=10→29 [16]=8 [17]=8 [18]=24 [19]=29 [20]=23
```

It is the 16th value — the only `10` on the first wrapped line of the literal,
sitting between a run of `29`s and the pair of `8`s.

### What this means

DFA state 15 is the state reached after consuming a single `.`. It previously
reported rule **10** (`{real}`) as an accepting rule; it now reports rule
**29**, the `.` catch-all that returns the raw character. The parser then sees
a bare `'.'` token, which no production accepts, and reports
`PARSE_SYNTAX_ERR`.

### Do not touch

`YY_ACCEPT` contains `10` at three other indices — **58, 73 and 137**. Those
are the genuine `{real}` accepting states (`1.5`, `.5`, `1.5e3`, `12.`) and
must stay `10`. Only index 15 changes. If a diff shows more than one changed
value, it is wrong.

### If you would rather regenerate than hand-edit

```sh
cd ~/code/cfitsio && make eval    # flex -o eval_l.c eval.l  (flex 2.6.4)
```
then diff `eval_l.c`'s `yy_accept` against `src/eval_l.rs`'s `YY_ACCEPT`. The
`eval.l` in that repo already carries both fixes.

---

## Tests to add

Add to the `#[cfg(test)] mod tests` block in `src/eval_f.rs`, mirroring the C
tests added in `~/code/cfitsio/tests/test_eval.c`
(`test_ffcrow_hex_constant_case`, `test_fftexp_bare_dot_rejected`).

`PARSE_SYNTAX_ERR` is **not** currently in the module's `use crate::fitsio::{…}`
list — add it.

There is an existing `test_ffcrow_integer_constant_bases` that already covers
`0x12f3` (lowercase, so it passes today). Extend it rather than duplicating,
then add a separate test for the bare dot.

```rust
    #[test]
    fn test_ffcrow_hex_constant_case() {
        /* 0x literals must be case-insensitive.  The lexer used to convert
        non-digits with (*p - 'a' + 10) without folding case, so uppercase
        A-F produced negative digit values: 0x1F evaluated to -1. */
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));

            for (expr, want) in [
                ("0x1f", 31),
                ("0x1F", 31),
                ("0xff", 255),
                ("0xFF", 255),
                ("0xFf", 255),
                ("0xa", 10),
                ("0xA", 10),
                ("0xabcdef", 11259375),
                ("0xABCDEF", 11259375),
                ("0x10", 16),
                ("0x0", 0),
            ] {
                assert_eq!(eval_lng::<1>(&mut f, expr), [want], "expr: {expr}");
            }

            /* the other radix prefixes must keep working */
            assert_eq!(eval_lng::<1>(&mut f, "0b1011"), [11]);
            assert_eq!(eval_lng::<1>(&mut f, "0o17"), [15]);

            fits_close_file(f, &mut status);
        });
    }

    #[test]
    fn test_fftexp_bare_dot_rejected() {
        /* A lone '.' must not lex as the double 0.0.  The {real} pattern's
        third alternative used to allow zero digits before the point. */
        with_temp_file(|filename| {
            let mut status = 0;
            let mut f = create_test_table(&to_buf(filename));

            for expr in [".", "..", "INTCOL + .", ". + 1", ".e5"] {
                let mut datatype = 0;
                let mut nelem: c_long = 0;
                let mut naxis = 0;
                let mut naxes = [0 as c_long; 5];
                let mut st = 0;
                fits_test_expr(
                    &mut f,
                    &cc(expr),
                    5,
                    &mut datatype,
                    &mut nelem,
                    &mut naxis,
                    &mut naxes,
                    &mut st,
                );
                assert_eq!(st, PARSE_SYNTAX_ERR, "expr should be rejected: {expr}");
            }

            /* valid reals are unaffected */
            for (expr, want) in [
                ("1.", 1.0),
                ("12.", 12.0),
                (".5", 0.5),
                ("1.5", 1.5),
                ("1.5e3", 1500.0),
                ("1.5E3", 1500.0),
                ("0.", 0.0),
            ] {
                assert_eq!(eval_dbl::<1>(&mut f, expr), [want], "expr: {expr}");
            }

            fits_close_file(f, &mut status);
        });
    }
```

Helpers used (`cc`, `to_buf`, `with_temp_file`, `create_test_table`,
`eval_lng`, `eval_dbl`) all already exist in that test module.

`fits_test_expr` comes from `crate::aliases::rust_api` (already glob-imported)
and resolves to `fftexp_safe`.

---

## Verification

Run each fix's test **before** applying that fix and confirm it fails — the two
bugs are independent, and it is worth proving each test catches its own:

```sh
cargo test --lib eval_f::tests::test_ffcrow_hex_constant_case
cargo test --lib eval_f::tests::test_fftexp_bare_dot_rejected
cargo test --lib
cargo test
cargo clippy --all-targets
```

The full suite must stay green (`CLAUDE.md`). Upstream's equivalent change
passed all 54 C tests with no regressions, so no fallout is expected — but
`{real}` is a lexer-wide pattern, so run the whole suite, not just `eval_f`.

Sanity check against the C oracle described in
[`PARSER_MIGRATION.md`](PARSER_MIGRATION.md) Appendix A, built against the
**fixed** `~/code/cfitsio` — it should now agree with rsfitsio on all of
`0x1F`, `0xFF`, `0xA`, `.`, `..`, `1.`, `.5`.

---

## Documentation to update after the fix

1. **`PARSER_SPEC.md` §6.1** — retitle from "Known defects in the current
   implementation" framing; note both are fixed in rsfitsio and upstream
   `47359ca`, and that `0x1F` → 31. Keep the historical description, since the
   spec documents a language whose older releases behave differently.
2. **`PARSER_SPEC.md` §6.2** — same; `.` is now `PARSE_SYNTAX_ERR`.
3. **`PARSER_SPEC.md` §2.6** — the `.`-related rows of the lexer quirk table
   (`.`, `..`, `.e5`, `1.e5`) need revising: `.` no longer produces `DOUBLE`.
   `1.e5` still splits as `DOUBLE(1.)` + `VAR(e5)`; `.e5` now errors.
4. **`PARSER_MIGRATION.md` §4.1** — the "Decision required — the hex bug"
   paragraph is resolved: option (b), fix it. Update the text.
5. **`PARSER_MIGRATION.md` §6, risk R7** — bug-compat decisions are now made;
   downgrade or close.
6. Add both expressions to the phase-0 differential corpus so the nom
   migration cannot silently reintroduce either.

---

## Coordination note

Written 2026-08-02 while another agent held `src/eval_l.rs`. The two source
edits are:

- `src/wrappers.rs` — add `tolower` + `isupper` after `islower`.
- `src/eval_l.rs` — one `use` addition, one expression in the `7 =>` arm, one
  array element in `YY_ACCEPT`.
- `src/eval_f.rs` — one `use` addition (`PARSE_SYNTAX_ERR`), two new tests.

None of them overlap the evaluation engine (`src/eval_y.rs`) or the public API,
so they should rebase cleanly over unrelated work.
