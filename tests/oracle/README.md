# Expression-parser differential oracle

Tools for checking rsfitsio's expression parser (`src/parser/`) against the
real CFITSIO one. Use these whenever you change the parser, the lexer, or the
node builders — or when you want to know what upstream actually does before
deciding whether a difference is a bug.

Nothing here is built or run by `cargo test`. It needs a C compiler, `flex`,
and a **built** CFITSIO source tree, none of which CI is guaranteed to have.

### Seeing CFITSIO's error messages

```text
ORACLE_MSGS=1 ./oracle < ../fixtures/eval_corpus.txt
```

appends the text CFITSIO pushed onto its error stack -- what a library consumer
reads back with `fits_read_errmsg` -- to each `ERR` line. Off by default so the
output still matches the golden file.

Worth knowing before changing any parser error message: **CFITSIO answers
"syntax error" to 840 of the corpus's rejections.** Its grammar is
type-stratified, so `1 && 2`, `'ab' > 1`, `STRCOL[1]` and
`T ? b1010 : b0101` are not type errors it diagnoses -- they simply match no
production, and bison supplies its own text. The specific messages
(`Function(bool) not supported`, `Index value must be an integer type`,
`Array sizes/dims do not match for binary operator`) come from semantic actions
and are a much smaller set.

---

## Why this exists

`src/parser/` was written to replace a flex scanner and a bison LALR parser
that had been transpiled from CFITSIO's generated C
(see [`../../PARSER_MIGRATION.md`](../../PARSER_MIGRATION.md)). The language it
implements is specified in [`../../PARSER_SPEC.md`](../../PARSER_SPEC.md), and
almost every behavioural claim in that document was established by *running*
these two programs, not by reading the grammar. The language has a lot of
surprising corners — `-2**2` is `+4`, `1|2*3` is `9`, `!X>5` parses as
`!(X>5)` but `!B||B` parses as `(!B)||B` — and the only reliable way to settle
a question about any of them is to ask the C.

## Prerequisites

A CFITSIO checkout that has been configured and built, so that
`.libs/libcfitsio.so*` exists:

```sh
git clone https://github.com/HEASARC/cfitsio ~/code/cfitsio
cd ~/code/cfitsio && ./configure && make
```

Every command below takes `CFITSIO=<path>`; it defaults to `~/code/cfitsio`.

---

## 1. `oracle` — expression-level

Runs `../fixtures/eval_corpus.txt` through CFITSIO's parser and prints the
**same one-line summary** that `../test_eval_corpus.rs` prints, so the output
diffs directly against our golden file.

```sh
make CFITSIO=~/code/cfitsio check
```

`check` prints `identical to ../fixtures/eval_corpus.golden`, or a unified diff.
See "Known divergences" below before treating a diff as a bug.

To ask about specific expressions instead of the whole corpus:

```sh
make CFITSIO=~/code/cfitsio oracle
printf '%s\n' '-2**2' '1|2*3' '!X>5' 'ANGSEP(0,0,0,90)' | ./oracle
```

The table it builds has one column per value sort and three rows:

| Column | TFORM | Row 1 | Row 2 | Row 3 |
|---|---|---|---|---|
| `INTCOL` | `1J` | 7 | −3 | 10 |
| `FLOATCOL` | `1E` | 2.5 | 4.0 | 0.5 |
| `STRCOL` | `10A` | `abc` | `de` | `fghij` |
| `BOOLCOL` | `1L` | T | F | T |
| `DBLCOL` | `1D` | 1.25 | −2.5 | 8.0 |
| `VECCOL` | `3E` | {1,2,3} | {4,5,6} | {7,8,9} |
| `BITS` | `8X` | `0xF0` | `0x0F` | `0xAA` |
| `MATRIX` | `6E`, `TDIM8='(2,3)'` | 1..6 | 11..16 | 21..26 |

plus keywords `INTKEY=42`, `DBLKEY=2.5`, `LOGKEY=T`, `STRKEY='hello'`.

> **This must stay in step with `create_corpus_table` and `probe` in
> `../test_eval_corpus.rs`.** If you change the schema, the values, or the
> output format on either side, change both — otherwise `make check` compares
> two different things and quietly reports garbage.

## 2. `lexdump` — token-level

The pattern definitions of `eval.l` with every action replaced by a print, so
you can compare token streams against `src/parser/lexer.rs` without needing a
FITS file. This is the tool that produced §2.6 of the spec.

```sh
make lexdump
printf '%s\n' 'b101' 'b1x0y' 'T' 'TT' '1.e5' '.' '0x1F' 'ab(' 'a (' | ./lexdump
```

```
BIT<b101>
VAR<b1x0y>          <- longest match beats rule order
BOOLEAN<T>
VAR<TT>             <- 'boolean' is single letters only
DOUBLE<1.>  VAR<e5> <- the exponent branch needs a digit before [eEdD]
CH<.>               <- a lone '.' is not a number (fixed; see spec 6.2)
HEXCONST<0x1F>
FUNC<ab(>
VAR<a>  CH<(>       <- no space allowed before the parenthesis
```

`lexdump.l` embeds a copy of `eval.l`'s definitions block. If upstream changes
those patterns, regenerate it — the header comment in `lexdump.l` says how.

---

## Regenerating the golden file

`../fixtures/eval_corpus.golden` is the contract for the Rust parser. It was
captured from rsfitsio's own flex/bison parser immediately before that parser
was deleted, so it encodes the behaviour we promised to preserve.

**Do not regenerate it to make a failing test pass.** Regenerate only when you
have *deliberately* changed the language, and say so in the commit message:

```sh
UPDATE_EVAL_GOLDEN=1 cargo test --test test_eval_corpus
git diff tests/fixtures/eval_corpus.golden   # review every line
```

To add coverage, append expressions to `../fixtures/eval_corpus.txt`, then
regenerate and check the new lines against `./oracle` output by hand.

---

## Known divergences from CFITSIO

`make check` currently reports two, both **pre-existing and deliberate** —
they predate the `nom` migration and live in the evaluation engine and the
keyword lookup, not the parser.

### `ANGSEP` with four constant arguments

```
ANGSEP(0,0,0,90)    C: 0.0        rsfitsio: 90.0
```

When all four arguments are constants the call is folded at parse time. In
CFITSIO 4.7.0 `case angsep_fct:` in `New_Func` is missing its `break;` and
falls through into `case min1_fct:`, so the folded form returns its first
argument instead of the separation. With a non-constant argument
(`ANGSEP(0.0*INTCOL,0.0,0.0,90.0)`) both agree on 90.

rsfitsio fixes this; see `test_ffcrow_angsep_constant_arguments` in
`src/eval_f.rs`. Fix submitted upstream.

### `#STRKEY` — a string-valued header keyword

```
#STRKEY             C: OK, TSTRING     rsfitsio: ERR 407 (BAD_C2I)
```

CFITSIO resolves a `#KEYWORD` reference whose value is a quoted string to a
string constant. rsfitsio's `find_keywd` (in `src/eval_f.rs`) reports
`BAD_C2I`. **This one looks like an rsfitsio bug rather than a deliberate
deviation** — it is recorded here because the golden file captured it, so the
Rust parser reproduces it faithfully and no test flags it. Worth investigating:
the divergence is in `find_keywd`'s use of `fits_get_keytype`, which the parser
migration did not touch.

---

## If you are here because a corpus test failed

1. `make CFITSIO=... check` — does the C agree with the golden? If yes, the
   regression is in `src/parser/`.
2. If the C also disagrees, the golden encodes an rsfitsio-specific behaviour.
   Check the two entries above, then `git log -p tests/fixtures/eval_corpus.golden`.
3. For a lexing question, reach for `lexdump` before `oracle`: it isolates
   tokenization from everything else, and most surprises live there.
