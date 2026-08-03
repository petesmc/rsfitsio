# Migrating the CFITSIO Expression Parser off flex/bison

**Goal:** replace `src/eval_l.rs` (transpiled flex DFA) and the LALR half of
`src/eval_y.rs` (transpiled bison tables + driver) with hand-written safe Rust
built on `nom`, while keeping the `Node`-arena evaluation engine and the entire
public API bit-for-bit compatible.

Do not make yourself as a co-author on any commits.

Read [`PARSER_SPEC.md`](PARSER_SPEC.md) first — this plan assumes its
terminology and its list of verified quirks.

---

## 1. Why this is worth doing

| Property | Today | After |
|---|---|---|
| Lines of transpiled generated code | ~8,000 | 0 |
| DFA / LALR tables | 4,700 integer literals across 11 arrays | none |
| `unsafe` blocks in the parse path | pervasive (`malloc`, raw `*mut`, `memcpy`, `slice::from_raw_parts_mut`, `static mut`) | none |
| Error messages | `"syntax error"` + a byte offset nobody surfaces | span-anchored, with the offending token |
| Ability to change the language | regenerate with flex+bison, re-transpile ~9k lines by hand | edit a Rust function |
| Re-entrancy | `yyguts_t` scanner object threaded through raw pointers | plain `&mut` |
| Stack behaviour | `malloc`-grown parse stack, `alloca`→`malloc` leak on error paths (§6.6) | bounded recursion with an explicit depth guard |

The `Do_*` evaluation engine (~7,300 lines) is **not** touched. Neither is
`eval_f.rs`. This is a surgical replacement of the front half.

---

## 2. Library selection

### 2.1 Recommendation

| Layer | Choice | Version |
|---|---|---|
| Tokenizer | **`nom`** | `8.0` |
| Parser | **hand-written Pratt / precedence-climbing over the token slice** | — |
| Errors | `nom`'s `Err::Error` + a crate-local `ParseError { span, kind }` | — |

**Do not try to express the grammar in nom combinators.** Reasons:

1. Every binary operator in `eval.y` is **left-recursive** (`expr '+' expr`).
   nom is PEG-flavoured and left recursion is an infinite loop; the standard
   workaround (`many0` of `(op, operand)` folded left) has to be repeated for
   all 15 precedence levels and becomes unreadable and slow.
2. nom 8 moved precedence handling out of core into the separate
   [`nom-language`](https://crates.io/crates/nom-language) crate (`0.1`).
   Its `precedence` combinator handles uniform operator tables — it does *not*
   handle the **type-directed** productions this grammar needs (`NOT` reaching
   through numeric prefixes, `{…}` bool-vs-numeric vectors, function overload
   on argument sort). You would end up writing the special cases by hand
   anyway, on top of a combinator you are fighting.
3. A Pratt parser for 15 levels is ~150 lines, is the *exact* dual of a bison
   precedence table (so review against `eval.y` is line-by-line), and handles
   right-associativity and the `?:`/range ternaries naturally.

nom earns its place on the **lexer**, where it is genuinely good: 29 rules of
byte-level pattern matching with longest-match/priority semantics, zero-copy
`&[u8]` slicing, and clean error positions.

### 2.2 Alternatives considered

| Option | Verdict |
|---|---|
| **`logos` 0.16** for the lexer | *Arguably the better technical fit.* It is a derive-based lexer generator whose semantics — longest match, then declaration priority — map 1:1 onto flex, so §2.6 of the spec transfers for free. Downsides: it is a proc-macro (a code generator, which is what we are moving away from), it compiles to a DFA you can't read, and the token *values* here need arbitrary Rust (bit-string expansion, `getData` callbacks) which means callbacks anyway. **Use it only if hand-written nom tokenization turns out to be a measurable hot spot** — it will not be; parsing happens once per query, not once per row. |
| **`chumsky` 1.0** | Best-in-class error recovery and a built-in `pratt` extension. Still `1.0.0-alpha.8` — not appropriate for a library at `0.464` that other crates depend on. Revisit if it stabilises. |
| **`winnow` 1.0** | nom fork, better ergonomics, actively developed. A legitimate substitute for nom here; the plan is identical. Pick nom if you want the larger ecosystem, winnow if you want the nicer API. Not both. |
| **`lalrpop` 0.23** | The *lowest-risk* option: `eval.y` translates almost mechanically into a `.lalrpop` file, preserving the type-stratified nonterminals and therefore every quirk in §3.4/§3.5 **for free**. But it is still a build-time code generator producing tables you can't read — the exact objection to bison, just cargo-integrated. Recommend only as a **fallback** if §4.4 (type-directed `NOT`) proves harder than expected. |
| **`pest`** | PEG + external `.pest` grammar file. Same code-generator objection, worse: PEG's ordered choice cannot express the type-directed productions at all. Reject. |
| Add `nom-language` for `precedence` | Optional. Try it for the pure-numeric levels if you like, but budget for replacing it. Not on the critical path. |

### 2.3 Dependency delta

```toml
[dependencies]
nom = { version = "8.0", default-features = false, features = ["alloc"] }
```

One new dependency, `no_std`-compatible, no proc macros, no build script.

---

## 3. Target architecture

```
src/parser/
    mod.rs          pub(crate) fn parse_expression(&mut ParseData, &[u8]) -> Result<i32, ParseError>
    token.rs        Token enum, Spanned<Token>, TokenKind
    lexer.rs        nom-based tokenizer  (replaces eval_l.rs)
    ast.rs          Ast enum — untyped syntax tree
    grammar.rs      Pratt parser: &[Spanned<Token>] -> Ast    (replaces yyparse's driver)
    lower.rs        Ast + ParseData -> Node arena             (replaces yyparse's 135 actions)
    error.rs        ParseError, diagnostics -> ffpmsg
```

### 3.1 The key design decision: split the type stratification

`eval.y` fuses three jobs into one LALR machine:

1. **syntax** — precedence, associativity, bracket matching;
2. **sort checking** — which of `expr`/`bexpr`/`sexpr`/`bits` an operand is;
3. **lowering** — building `Node`s, promotions, dimension checks.

The migration separates them:

- `grammar.rs` does (1) only, over **one untyped `Ast`**. `a + b` is a single
  `Ast::Binary` node regardless of operand sorts.
- `lower.rs` does (2) and (3): it walks the `Ast`, resolves identifiers, learns
  each subtree's sort, and calls the *existing, unchanged* `New_Const`,
  `New_Column`, `New_Offset`, `New_Unary`, `New_BinOp`, `New_Func`,
  `New_FuncSize`, `New_Deref`, `New_GTI`, `New_REG`, `New_Vector`,
  `Close_Vec`, `New_Array`, `Test_Dims`, `Copy_Dims`.

This is what collapses 135 rules into roughly 30 `Ast` variants: the four
`NELEM(` implementations become one, the six `?:` rules become one, the
twenty-five deref rules become one.

**The cost of this decision** is that the syntactic constraints listed in
spec §3.4 — which the LALR grammar enforced *for free* by not having a rule —
now have to be enforced explicitly in `lower.rs` as errors. There are nine of
them and they are enumerated; §4.5 covers the mechanism.

### 3.2 `Ast` shape (sketch)

```rust
pub(crate) enum Ast {
    Long(c_long),
    Double(f64),
    Boolean(bool),
    Str(Box<[u8]>),          // quotes already stripped
    BitStr(Box<[u8]>),       // already expanded to '0'/'1'/'x'
    RowRef, NullRef, SNullRef,
    KeywordConst(Box<[u8]>), // '#'-prefixed, resolved during lowering
    Ident(Box<[u8]>),        // resolved during lowering -> COLUMN/BCOLUMN/...
    Offset { col: Box<Ast>, off: Box<Ast> },    // col must be Ast::Ident
    Unary  { op: UnOp, arg: Box<Ast> },
    Binary { op: BinOp, lhs: Box<Ast>, rhs: Box<Ast> },
    Ternary{ cond: Box<Ast>, then: Box<Ast>, els: Box<Ast> },
    Range  { val: Box<Ast>, lo: Box<Ast>, hi: Box<Ast> },   // x = lo : hi
    Deref  { base: Box<Ast>, idx: Vec<Ast> },               // 1..=5
    Vector (Vec<Ast>),                                      // { … }
    Call   { kind: CallKind, name: Box<[u8]>, args: Vec<Ast> },
}

pub(crate) enum CallKind { Function, BFunction, IFunction,
                           GtiFilter, GtiOverlap, GtiFind, RegFilter }
```

`Ast` is `Box`-heavy and allocated once per `ffiprs` call — irrelevant next to
the per-row evaluation cost.

### 3.3 Public surface: unchanged

`ffiprs` currently does:

```rust
fits_parser_yylex_init_extra(lParse, &mut yylex_scanner);
fits_parser_yyrestart(ptr::null_mut(), yylex_scanner.as_deref_mut().unwrap());
*status = fits_parser_yyparse(yylex_scanner.as_deref_mut().unwrap(), lParse);
fits_parser_yylex_destroy(yylex_scanner.unwrap());
```

becomes

```rust
*status = crate::parser::parse_expression(lParse);
```

`ParseData.expr` / `.index` / `.is_eobuf` stay for now (other code reads them);
`.is_eobuf` and `.index` become dead and can be removed in a later cleanup.
`lParse.resultNode`, `lParse.status`, `lParse.nNodes` are set exactly as
before. Nothing in `eval_f.rs` above `ffiprs` changes.

---

## 4. Implementation detail, by phase

### Phase 0 — Build the differential oracle *(prerequisite, ~1 day)*

**This is the highest-value single step and must come first.** It has already
been prototyped and works.

Two binaries, both driven from a shared corpus file, one line per expression:

- **C oracle** — links `~/code/cfitsio/.libs/libcfitsio.so`, creates a
  `mem://` binary table with one column per sort (`1J`, `1D`, `10A`, `1L`,
  `3E`, `8X`), and for each input line calls `fits_test_expr` then
  `fits_calc_rows` for rows 1..N, printing `(status, datatype, nelem, naxis,
  values…)`.
- **Rust oracle** — same table via `new_mem_file`, same calls through
  `crate::aliases::rust_api`.

A working prototype of the C side is in Appendix A — it was used to produce
every **[verified]** claim in `PARSER_SPEC.md`. Build with:

```sh
gcc -I ~/code/cfitsio -o oracle oracle.c \
    ~/code/cfitsio/.libs/libcfitsio.so.10 -lm \
    -Wl,-rpath,$HOME/code/cfitsio/.libs
```

Put both into `tests/oracle/` (gated behind a `cfitsio-oracle` feature so CI
without a cfitsio build still passes). Diffing the two outputs over a corpus is
the acceptance test for every later phase.

**Seed corpus** (~600 lines, assembled in this order):
1. Every expression string appearing in the 144 `#[test]`s in `src/eval_f.rs`.
2. Every expression in `tests/*.rs` and `docs/` filter examples.
3. One line per grammar production (135) — mechanically derivable from the
   bison `-v` output (`eval_y.output`, already generated).
4. Every row of spec §2.6 (lexer quirks) and §3.4/§3.5 (structural
   constraints), which is where a rewrite will actually break.
5. Operator precedence probes: all 15·15 adjacent-level pairs.

### Phase 1 — Tokenizer *(replaces `src/eval_l.rs`, ~700 lines)*

`lexer.rs` exposes:

```rust
pub(crate) fn tokenize(input: &[u8]) -> Result<Vec<Spanned<Token>>, ParseError>
```

**Do not attempt to emulate flex's DFA.** Emulate its *semantics*: at each
position, try every rule, keep the longest match, break ties by declaration
order. With 29 rules and short inputs this is trivially fast and — crucially —
directly auditable against spec §2.4.

```rust
fn next_token(i: &[u8]) -> Option<(usize, Token)> {
    // rules in eval.l declaration order
    const RULES: &[fn(&[u8]) -> IResult<&[u8], Token>] = &[
        lex_bit, lex_oct, lex_hex, lex_bitconst, lex_octconst, lex_hexconst,
        lex_integer, lex_boolean, lex_real, lex_constant, lex_string,
        lex_variable, lex_function, lex_intcast, lex_fltcast, lex_power,
        lex_not, lex_or, lex_and, lex_eq, lex_ne, lex_gt, lex_lt,
        lex_gte, lex_lte, lex_xor, lex_newline,
    ];
    let mut best: Option<(usize, usize, Token)> = None;   // (len, rule_idx, tok)
    for (idx, rule) in RULES.iter().enumerate() {
        if let Ok((rest, tok)) = rule(i) {
            let len = i.len() - rest.len();
            if len > 0 && best.as_ref().is_none_or(|(bl, _, _)| len > *bl) {
                best = Some((len, idx, tok));
            }
        }
    }
    best.map(|(len, _, tok)| (len, tok))
        .or_else(|| i.first().map(|&c| (1, Token::Char(c))))   // flex's `.` rule
}
```

Each `lex_*` is a small nom combinator; e.g.

```rust
// bit  ([bB][01xX]+)
fn lex_bit(i: &[u8]) -> IResult<&[u8], Token> {
    let (i, _)    = one_of("bB")(i)?;
    let (i, body) = take_while1(|c| matches!(c, b'0'|b'1'|b'x'|b'X'))(i)?;
    Ok((i, Token::BitStr(body.to_vec().into())))
}
```

**Ordering hazards to get right** (all in spec §2.6, all in the corpus):
`b101` vs `b1x0y`; `T` vs `TT`; `oct123`; `.` as `DOUBLE`; `1.e5` splitting;
`0x`/`0b`/`0o` with no digits; `X(` not being a function; `#$a$b$` greedy vs
`$a$b$` non-greedy.

**Identifier resolution moves out of the lexer.** `lex_variable` emits
`Token::Ident(bytes)`; `lex_constant` emits `Token::KeywordConst(bytes)` for
non-builtin `#…`. The `getData`/`find_variable` callback fires in `lower.rs`.
This removes the `ParseData` back-pointer from the lexer entirely — the lexer
becomes a pure function, which is why it can be unit-tested in isolation
against the C.

> **Behavioural risk (must be tested):** today a lookup failure aborts the
> *lex*, so the parse dies at the first unknown column. Deferred resolution
> means the whole expression is tokenized and parsed first. For a *valid*
> parse the resolution order is unchanged (left-to-right during lowering, so
> `varData[]`/`colData[]` registration order is preserved — this matters,
> because `ffiprs` reports columns in registration order). For an *invalid*
> parse, a syntax error may now be reported where a "column not found" was
> reported before. Corpus entries: `NOSUCHCOL +`, `NOSUCHCOL[`,
> `1 + NOSUCHCOL +`. Decide and document; matching the old ordering exactly is
> possible but requires a pre-pass and is probably not worth it.

**Bug-compat: resolved.** The hex-case bug (spec §6.1) and the bare-`.` bug
(spec §6.2) were **fixed upstream** in `~/code/cfitsio` commit `47359ca`, so
there is no bug-for-bug question left: the new lexer parses `0x` literals
case-insensitively and rejects a lone `.`. rsfitsio itself has not been
updated yet — port the fixes per [`BUGS_TO_FIX.md`](BUGS_TO_FIX.md) *before*
phase 0, so the oracle and the corpus agree from the start. Add corpus lines
for `0x1F`, `0xFF`, `0xABCDEF`, `.`, `..` and `.e5` either way.

### Phase 2 — Pratt parser *(replaces the yyparse driver, ~400 lines)*

`grammar.rs` operates on `&[Spanned<Token>]` with a cursor. Structure:

```rust
fn parse_expr(&mut self, min_bp: u8) -> Result<Ast, ParseError> {
    let mut lhs = self.parse_prefix()?;
    loop {
        let (lbp, rbp) = match self.peek_infix() { Some(bp) => bp, None => break };
        if lbp < min_bp { break }
        self.advance();
        lhs = self.parse_infix_rest(lhs, op, rbp)?;
    }
    Ok(lhs)
}
```

**Binding powers, transcribed from spec §3.1.** Left-assoc level *n* →
`(2n, 2n+1)`; right-assoc → `(2n+1, 2n)`:

| Level | Assoc | Tokens | (lbp, rbp) |
|---:|---|---|---|
| 1 | left | `,` `=` `:` `{` `}` | (2, 3) |
| 2 | right | `?` | (5, 4) |
| 3 | left | `OR` | (6, 7) |
| 4 | left | `AND` | (8, 9) |
| 5 | left | `EQ` `NE` `~` | (10, 11) |
| 6 | left | `GT` `LT` `LTE` `GTE` | (12, 13) |
| 7 | left | `+` `-` `%` | (14, 15) |
| 8 | left | `*` `/` | (16, 17) |
| 9 | left | `\|` `&` `XOR` | (18, 19) |
| 10 | right | `POWER` | (21, 20) |
| 11 | prefix | `NOT` | rbp 22 — **but see §4.4** |
| 12 | prefix | `INTCAST` `FLTCAST` | rbp 24 |
| 13 | prefix | unary `+`/`-` (`UMINUS`) | rbp 26 |
| 14 | postfix | `[` | lbp 28 |

The two counter-intuitive levels (bitwise above `*`/`/`; `UMINUS` above
`POWER`) fall straight out of this table — do not "fix" them.

Prefix parsing handles: literals, `(`…`)`, `{`…`}` vector literals, `NOT`,
`INTCAST`/`FLTCAST`, unary `±`, and the seven call kinds (a
`Token::Function(name)` already carries its `(`, so a call is
`Function` → arg list → `)`).

Postfix handles `[` index lists (1–5, comma-separated, `]`) and `{`…`}`
offsets. The ternaries are infix: `?` consumes `then`, expects `:`, consumes
`else` at rbp 4; `=` consumes `lo`, expects `:`, consumes `hi` and builds
`Ast::Range`.

### 4.4 The `NOT` problem — the one place a naive rewrite breaks

Spec §3.5: `!X > 5` must parse as `!(X > 5)`, but `!B || B` must parse as
`(!B) || B`. A prefix operator with a single binding power cannot do both,
because the correct extent depends on the *sort* of what follows.

Three viable implementations, in order of preference:

**(a) Parse wide, re-associate during lowering.** Parse `NOT` with rbp = 6
(just above `OR`), i.e. let it swallow up to but not including `||`. Then in
`lower.rs`, when lowering `Unary{NOT, arg}`:
- lower `arg`;
- if `arg`'s sort is `BOOLEAN` or `BITSTR` → emit `New_Unary(NOT)`, done;
- if `arg` is a `Binary` whose *result* sort is boolean → correct already;
- the failing shape is `!B || B` parsed as `!(B||B)` — so **exclude `OR`, `AND`
  and `EQ`/`NE`-on-booleans from NOT's reach** by giving NOT rbp 10 (above
  `AND`, below `EQ NE ~`). Check: `!X>5` → `>` at lbp 12 > 10 → swallowed →
  `!(X>5)` ✓. `!B||B` → `||` at lbp 6 < 10 → stops → `(!B)||B` ✓.
  `!B&&1>2` → `&&` at lbp 8 < 10 → `(!B)&&(1>2)` = 0 ✓
  (**[verified]**: 0, whereas `!(B&&1>2)` = 1 — this pair discriminates,
  `!B&&B` does not since both readings give 0).
  `!X+1>5` → `+` (14) and `>` (12) both > 10 → `!((X+1)>5)` ✓.

  **rbp = 10 reproduces all six discriminating cases.** This is the
  recommendation. The one known divergence is `!B == B` (level 5, lbp 10 ≥ 10
  → swallowed → `!(B==B)`), where LALR gives `(!B)==B`. For booleans
  `!(a==b) ≡ (!a)==b`, so the value is identical (**[verified]**: both 0) and
  only the node tree differs. `~` is not defined on booleans at all
  (**[verified]**: `!B~B` is a syntax error), so it poses no risk. Add corpus
  lines for all of these.

**(b) Two-pass: parse `NOT` greedily to the tightest operand, and if lowering
finds the operand is numeric, re-parse.** Correct but ugly.

**(c) Sort-aware Pratt: thread an expected-sort parameter through
`parse_expr`.** This re-introduces the stratification the design just removed.
Only do this if (a) fails the corpus.

If none of these survive the corpus, that is the signal to fall back to
`lalrpop` (§2.2), which preserves the behaviour by construction.

### 4.5 Phase 3 — Lowering *(replaces the 135 reduce actions, ~1,200 lines)*

`lower.rs` is a recursive walk returning `(node_index, sort)`. It is the
largest new file and it is where all the semantics live.

```rust
struct Lowerer<'a> { p: &'a mut ParseData }

impl Lowerer<'_> {
    fn lower(&mut self, ast: &Ast) -> Result<i32, ParseError>;
    fn sort_of(&self, node: i32) -> Sort;   // reads Nodes[node].type
}
```

Work items, each a direct transcription of a group of `eval.y` actions:

1. **Identifier resolution** — `Ast::Ident` → `fits_parser_yyGetVariable`
   equivalent → `New_Column`. `Ast::KeywordConst` → `find_keywd` → literal.
2. **Binary dispatch** — one `match (op, lhs_sort, rhs_sort)` reproducing
   `expr+expr` / `sexpr+sexpr` / `bits+bits`, the `PROMOTE` macro, the
   `expr*bexpr` / `bexpr*expr` coercions, the LONG-only check for `& | XOR`
   on `expr`, and the `SIZE()` adjustments (bit/string concat sums, comparison
   results forced to `nelem = 1`).
3. **The nine structural constraints (spec §3.4)** — these no longer come for
   free. Each needs an explicit rejection with the *same* `PARSE_SYNTAX_ERR`
   status:
   - offset only on `Ast::Ident`;
   - no deref on `sexpr`;
   - `NOT` only on boolean/bits;
   - GTI/REG filenames must be `Ast::Str`, not a general string expression;
   - no `?:` over bits;
   - no `+` between numeric and boolean;
   - no `ISNULL(bits)`;
   - no `DEFNULL(bits,bits)`;
   - no range chaining.

   **Write these nine as a checklist and a matching corpus block. This is the
   single most likely source of silent behaviour drift.**
4. **Vector literals** — collect elements, boolean iff *all* elements are
   boolean, else numeric with `type = max(element types)`; chunk at
   `MAXSUBS = 10` via `Close_Vec` + `New_Vector` exactly as
   `bvector`/`vector` do.
5. **Function dispatch** — one table keyed by
   `(CallKind, upper_name, arity, arg_sorts)` → builder closure, transcribed
   from spec §4. Roughly 60 entries. Unknown name/arity → the corresponding
   `"Function(...) not supported"` message. This replaces ~2,000 lines of
   `FSTRCMP` chains with a table, and is where most of the readability win is.
6. **Constant folding** — `NELEM`, `NAXIS`, `NAXES`, and the `OPER(x) ==
   CONST_OP` fast paths for `AXISELEM` / `ARRAY` / `SETNULL` / `STRMID`
   must fold at lowering time, because downstream code inspects
   `Nodes[n].value.data.lng` directly.
7. **Dimension checks** — `Test_Dims` / `Copy_Dims` calls at every site
   `eval.y` has them (`?:`, `ARCTAN2`, `MIN`/`MAX`/`DEFNULL` 2-arg, `ANGSEP`,
   `NEAR`, `CIRCLE`, `BOX`, `ELLIPSE`).

### 4.6 Phase 4 — Robustness

- **Recursion depth guard.** Recursive descent on `((((…))))` will blow the
  native stack where bison grew a heap array. Add a `depth` counter to the
  parser, limit ~200 (bison's `YYINITDEPTH` is 100, `YYMAXDEPTH` far higher),
  return `PARSE_SYNTAX_ERR` past it. Corpus: 1,000 nested parens, 1,000
  nested `!`, 1,000-element vector literal.
- **Fuzz target.** `cargo-fuzz` over `parse_expression` with the corpus as the
  seed set. The old code has `malloc`/`memcpy`/raw-pointer paths; the new code
  should be panic-free by construction, so any panic is a bug. Run the same
  inputs through the C oracle to catch divergence, not just crashes.
- **Error messages.** Keep `ffpmsg` output identical for the messages the
  tests assert on. Everything else can improve — attach the byte span and
  render a caret line. `ParseData.status` and the returned `*status` must stay
  exactly as before.

### 4.7 Phase 5 — Deletion

Once the corpus is green:

- delete `src/eval_l.rs` entirely;
- delete `src/eval_y.rs` lines 1–525 (shims, symbol constants, LALR tables),
  834–1029 (bison debug scaffolding), 1120–6964 (`fits_parser_yyparse`);
- delete `src/eval_tab.rs`'s `yysymbol_kind_t` half, keep the token ids that
  double as type tags (`BOOLEAN`/`LONG`/`DOUBLE`/`STRING`/`BITSTR` are stored
  in `Node.type` and read all over the evaluation engine — **do not renumber
  them**);
- remove `ParseData.index` / `.is_eobuf`;
- rename what remains of `eval_y.rs` to `src/eval_engine.rs` (optional but
  clarifying — it is no longer "the y file").

---

## 5. Sequencing and effort

| Phase | Deliverable | Est. | Gate |
|---|---|---|---|
| 0 | C + Rust oracles, 600-line corpus, diff harness | 1–2 d | oracles agree on today's code |
| 1 | `token.rs` + `lexer.rs` + unit tests | 2–3 d | token-stream diff vs. instrumented `eval_l.rs` on the whole corpus |
| 2 | `ast.rs` + `grammar.rs` | 2–3 d | every corpus line parses or errors; AST snapshot tests |
| 3 | `lower.rs` | 5–8 d | **full oracle diff green**; `cargo test` green |
| 4 | depth guard, fuzzing, error messages | 2–3 d | 24 h fuzz, no panics, no divergence |
| 5 | delete dead code, docs | 1 d | line count drops ~8,000 |

**≈ 3 weeks of focused work.** Phases 1 and 2 are independently testable and
can land behind a `#[cfg(feature = "new-parser")]` switch so the old path
stays live until phase 3 completes.

### 5.1 Landing strategy

Add a temporary feature flag:

```toml
[features]
new-parser = []
```

```rust
// ffiprs
#[cfg(feature = "new-parser")]
{ *status = crate::parser::parse_expression(lParse); }
#[cfg(not(feature = "new-parser"))]
{ /* existing flex/bison path */ }
```

Run the full test suite both ways in CI until phase 5. This keeps every
intermediate commit shippable, which matters for a crate at `0.464` with a
`cdylib` ABI surface.

---

## 6. Risk register

| # | Risk | Likelihood | Mitigation |
|---|---|---|---|
| R1 | Type-directed `NOT` (§4.4) can't be matched by a fixed binding power | **medium** | rbp = 10 reproduces all five verified cases; corpus block; fallback = `lalrpop` |
| R2 | One of the nine structural constraints (§3.4) silently becomes *accepted* | **high** if not tested | explicit checklist in `lower.rs`, one corpus line each, `#[test]` per constraint |
| R3 | Deferred identifier resolution changes `varData[]`/`colData[]` registration order | medium | lower left-to-right; assert column order in `ffiprs` tests |
| R4 | Constant-folding sites missed → downstream reads uninitialised `value.data.lng` | medium | grep `CONST_OP` in `eval.y`, one corpus line per site |
| R5 | Recursion depth → stack overflow on adversarial input | **high** without a guard | explicit depth limit + fuzzing (§4.6) |
| R6 | Lexer longest-match ties resolved differently | medium | rule order is explicit in the `RULES` array; §2.6 table is the test |
| R7 | ~~Bug-compat decisions (hex §6.1, bare `.` §6.2)~~ — **closed**, fixed upstream in `47359ca`; port to rsfitsio via `BUGS_TO_FIX.md` before phase 0 | — | corpus lines for `0x1F`, `0xFF`, `.`, `..` |
| R8 | `Node.type` tag renumbering breaks the evaluation engine | low but fatal | do not touch `eval_tab.rs` numeric values; add a `const _: () = assert!(LONG == 259)` guard |
| R9 | Scope creep into the `Do_*` engine | medium | explicitly out of scope; the engine is not touched in any phase |

---

## 7. What "done" looks like

1. `cargo test --lib && cargo test` green.
2. The oracle diff over the ≥600-line corpus is empty, including every row of
   spec §2.6, §3.4, §3.5 and the precedence probes.
3. `src/eval_l.rs` gone; `src/eval_y.rs` down from 15,305 to ~7,900 lines with
   zero LALR tables and zero `malloc` in the parse path.
4. `rg 'unsafe' src/parser/` returns nothing.
5. `TODO.md` items relating to the parser's `libc`/`malloc` usage closed.

---

## 8. Follow-on opportunities (explicitly out of scope)

Once the front end is Rust, these become cheap — but none of them should be
attempted during the migration:

- Surface real diagnostics (`expected ')', found '\n'` with a caret) through
  `ffpmsg`, instead of the current bare `"syntax error"`.
- Cache parsed `Ast`s keyed by expression string — `ffcalc` on a large table
  re-parses per call in some paths.
- Constant-fold the whole tree at lowering (`2+3*4` currently allocates three
  `Node`s and evaluates per row).
- Replace the `Nodes: Vec<Node>` arena + `i32` indices with a typed tree —
  large, invasive, and touches the entire `Do_*` engine.
- Apply the same treatment to `src/grparser.rs` (the group-template parser),
  which is hand-written C-style recursive descent rather than generated, and
  is a separate, smaller project.

---

## Appendix A — C differential oracle (working prototype)

Reads one expression per line from stdin; prints parse status, inferred
datatype/shape, and the evaluated values for row 1. Columns: `X` int,
`Y` double, `S` string, `B` logical, `V` 3-element float vector, `BIT` 8-bit.

```c
#include <stdio.h>
#include <string.h>
#include "fitsio.h"

int main(void) {
  fitsfile *fptr; int status = 0;
  char *ttype[] = {"X","Y","S","B","V","BIT"};
  char *tform[] = {"1J","1D","10A","1L","3E","8X"};
  long   xv[3] = {7, -3, 10};
  double yv[3] = {2.5, 4.0, 0.5};
  char   bv[3] = {1, 0, 1};
  float  vv[9] = {1,2,3,4,5,6,7,8,9};
  char  *sv[3] = {"abc","de","fghij"};

  fits_create_file(&fptr, "mem://t.fits", &status);
  fits_create_tbl(fptr, BINARY_TBL, 3, 6, ttype, tform, NULL, "T", &status);
  fits_write_col(fptr, TLONG,    1, 1, 1, 3, xv, &status);
  fits_write_col(fptr, TDOUBLE,  2, 1, 1, 3, yv, &status);
  fits_write_col(fptr, TSTRING,  3, 1, 1, 3, sv, &status);
  fits_write_col(fptr, TLOGICAL, 4, 1, 1, 3, bv, &status);
  fits_write_col(fptr, TFLOAT,   5, 1, 1, 9, vv, &status);
  if (status) { fprintf(stderr, "setup %d\n", status); return 1; }

  char line[4096];
  while (fgets(line, sizeof line, stdin)) {
    line[strcspn(line, "\n")] = 0;
    if (!line[0]) continue;
    int dt = 0, naxis = 0, st = 0, anyn = 0;
    long nelem = 0, naxes[5] = {0};

    fits_test_expr(fptr, line, 5, &dt, &nelem, &naxis, naxes, &st);
    if (st) { printf("%-34s ERR %d\n", line, st); fits_clear_errmsg(); continue; }

    long n = (nelem > 0 ? nelem : 1); if (n > 8) n = 8;
    double out[8] = {0}, nul = -999;
    fits_calc_rows(fptr, TDOUBLE, line, 1, n, &nul, out, &anyn, &st);
    if (st) { printf("%-34s EVALERR %d\n", line, st); fits_clear_errmsg(); continue; }

    printf("%-34s dt=%d [", line, dt);
    for (long i = 0; i < n; i++) printf("%s%g", i ? ", " : "", out[i]);
    printf("]%s\n", anyn ? " (null)" : "");
  }
  return 0;
}
```

**Known gap:** the `TDOUBLE` output path cannot render `sexpr` or `bits`
results (they report `EVALERR 433`). Add a `dt == TSTRING` branch using
`fits_calc_rows(..., TSTRING, ...)` with a pre-allocated `char*[]`, and a
`dt == TBIT` branch, before using the corpus for string/bit coverage.

## Appendix B — Lexer oracle

For token-level differential testing (phase 1), build a standalone flex
scanner from the `%{ … %}`-free portion of `eval.l`: copy the pattern
definitions block verbatim, replace every action with
`printf("TOKENNAME<%s> ", yytext)`, add `%option noyywrap` and a
`main(){ yylex(); }`. This reproduces flex's exact longest-match/rule-order
behaviour with no `ParseData` dependency, and is what produced the §2.6 table
in `PARSER_SPEC.md`.

---

## 9. Outcome

The migration is **complete**. What was actually built, and where it differed
from this plan:

| Plan | Outcome |
|---|---|
| `nom` for the lexer, hand-written Pratt parser | As planned. `nom` 8.0, one new dependency, no proc macros or build script. |
| ~8,000 lines of generated code removed | `src/eval_l.rs` deleted (1,688 lines); `src/eval_y.rs` 15,305 → ~8,950. All 11 DFA/LALR tables gone. |
| New front end ~2,500 lines | 3,180 lines across `src/parser/` (7 files), including ~700 lines of unit tests. |
| Phase 0 oracle + corpus | `tests/test_eval_corpus.rs` + `tests/fixtures/eval_corpus.{txt,golden}`. Started at 503 expressions (golden captured from the flex/bison parser before removal), later expanded to 1,852 and validated line-by-line against CFITSIO. The oracle itself is checked in under `tests/oracle/`. |
| `NOT` binding power (§4.4 option (a)) | Worked. `BP_NOT_OPERAND = 10` reproduces every discriminating case in the corpus. `lalrpop` was not needed. |
| Feature-flagged landing | Skipped. The corpus made the cutover verifiable in one step, so the flag would have been dead weight. |

### 9.1 What the plan underestimated

Three behaviours cost more than the estimate, all of them consequences of
moving name resolution out of the lexer (risk **R3**, which was real):

1. **Status codes, not just accept/reject.** `find_column` set `lParse.status`
   and returned `pERROR` (-1); bison treats any token `<= YYEOF` as
   end-of-input, so the parse *succeeded* on the truncated prefix and `ffiprs`
   reported the resolver's status. `NOSUCHCOL` is therefore a `COL_NOT_FOUND`
   (202), not a syntax error. 26 corpus lines depended on this.

2. **Truncation semantics.** The same mechanism means a failed lookup silently
   cuts the token stream — including the terminating newline. Since
   `line: expr '\n'` requires that newline, `1.e5` (which lexes as `1.` then
   the unresolvable name `e5`) is a *syntax* error, while a bare `NOSUCHCOL`
   is not. This is reproduced by `src/parser/resolve.rs`, a pre-pass that
   resolves every name in source order and truncates at the first failure —
   which also preserves column registration order for free.

3. **`BOOLEAN` is not `expr`.** A boolean-valued node is a `bexpr`, a different
   nonterminal, so most arithmetic and comparison productions reject it:
   `INTCOL + BOOLCOL` and `1<2==1` are syntax errors, while `INTCOL * BOOLCOL`
   is fine because `expr '*' bexpr` exists. The lowering pass distinguishes
   `is_expr` (LONG/DOUBLE) from `is_numeric` (adds BOOLEAN) for exactly this.

### 9.2 One nom trap worth recording

`digit0` followed by `digit1` does not backtrack. The `{real}` exponent branch
`[0-9]*"."*[0-9]+[eEdD]…` fails on `1e308` under naive translation, because
`digit0` eats the `1` and `digit1` then has nothing. A regex engine backtracks;
nom does not. `m_real::expo` splits the mantissa by hand. This is the only
place the impedance mismatch bit, but it would have been silent without the
corpus.

### 9.3 Deferred

* **`ACCUM(BOOLCOL)` / `ACCUM(BITS)` crash.** The grammar lowers these to a
  `LONG` binop over a char-valued operand, so `Do_BinOp_lng` reads a `char`
  buffer through a `*long`. CFITSIO returns garbage (`ACCUM(bool)` → 65537);
  the Rust debug build traps on the misaligned read. This is a defect in the
  *evaluation engine*, which this migration deliberately did not touch, so the
  two expressions are excluded from the corpus with a comment. Worth fixing
  separately, upstream and here.
* **Fuzzing** (§4.6). The corpus and the depth guard are in; a `cargo-fuzz`
  target over `parse_expression` is still worth adding.
* **Error-message quality** (§8). The parser now carries byte spans, but
  `ffpmsg` output is still one line. The machinery to do better is in place.

### 9.4 What the corpus expansion found

The corpus was later grown from 503 to 1,852 expressions, mostly as systematic
grids: every binary operator against every pair of operand sorts, every unary
operator against every sort, every `?:` branch combination, function arities 0
to 5, the `MAXSUBS` vector-chunking boundary, and a block of malformed input.
Every line was then checked against CFITSIO with `tests/oracle`.

That found **four real bugs in the new parser**, three of which were memory
unsafety rather than a wrong answer:

1. `ABS(BOOLCOL)` — the `BOOLEAN` arm of `func1` fell through to the numeric
   function list. `FUNCTION bexpr ')'` supports only `SUM`, `NELEM` and
   `ACCUM`; the others handed `Do_Func` a `char` buffer to read through a
   `*long`. **Segfault.**
2. `ARCTAN2(INTCOL,BITS)`, `MIN(INTCOL,STRCOL)`, `SETNULL(1,STRCOL)` … — the
   two-argument functions did not check operand sorts at all. Same crash.
   Fixed with a `require_expr` / `require_numeric` guard applied uniformly to
   the 2-, 4- and region-function paths and to the GTI/region node arguments.
3. `ANGSEP` and the region functions had the same gap.
4. `+BOOLCOL` — unary `+` is `$$ = $2` in `eval.y`, so it was implemented as a
   pass-through with no node and therefore no sort check. But the production
   is `'+' expr`, so `+BOOLCOL`, `+STRCOL` and `+BITS` are syntax errors. Now
   modelled as `UnOp::Plus`, which checks and then returns its operand.

and **two pre-existing bugs elsewhere in the crate**, both reachable from a
filter expression:

5. `ffdtyp` (`fits_get_keytype`) tested `cval[0]` against a backslash where the
   C tests against a single quote — `'\''` had been transpiled as `'\\'`. Every
   quoted keyword value was classified as an integer, so `#STRKEY` reported
   `BAD_C2I`, as did `ffc2x` / `ffc2r` / `ffc2d` for any string. The extern
   wrapper had a second copy of the same body; it now delegates to
   `ffdtyp_safe` per the crate convention.
6. `ffgcrd_safe` computed `namelen - 1` on a `usize`. The C does it in `int`,
   where an empty name gives `max(-1, 1) = 1`; in Rust it underflowed and
   **panicked**. Reachable from `$$`, which lexes as a variable with an empty
   name.

The lesson is that the systematic grids paid for themselves: none of these
were found by the 503-line corpus, the 2,300-test unit suite, or code review,
and three of them were crashes on inputs a user can type.

### 9.5 Upstream status

Two of the CFITSIO defects this work turned up are pending upstream:

| PR | Defect | Where it shows here |
|---|---|---|
| [cfitsio#152](https://github.com/HEASARC/cfitsio/pull/152) | `angsep_fct` falls through to `min1_fct` in `New_Func` for lack of a `break`, so a fully-constant `ANGSEP` returns its first argument | the only 5 lines where `tests/oracle` still differs |
| [cfitsio#153](https://github.com/HEASARC/cfitsio/pull/153) | uppercase `0x` literals miscomputed; a bare `.` lexes as `0.0` | §6.1 and §6.2 of `PARSER_SPEC.md`, fixed here first |

When #152 merges and CFITSIO is rebuilt, `make -C tests/oracle check` should
report no differences at all, and the "Known divergences" section of
`tests/oracle/README.md` can be deleted.

The other defects found here are rsfitsio-only transpilation errors
(`ffdtyp`'s quote test, `ffgcrd_safe`'s underflow) and do not apply upstream.
The `ACCUM(BOOLCOL)` type confusion in §9.3 *does* apply upstream and has not
been reported.

---

## 10. Replacing the evaluation engine

The parser migration left an `Ast` between the front end and the `Node` arena.
That seam is what makes the engine replaceable without touching the parser.

### 10.1 What is wrong with the arena

Measured, not asserted:

| | |
|---|---|
| `lParse.Nodes[…]` lookups in `eval_y.rs` | 1,371 |
| scalar-versus-vector special cases | 57 |
| `DoOp` | `Option<fn(&mut ParseData, usize)>` |

`DoOp` is a hand-rolled vtable whose `self` is an *index*, so no operation can
hold a reference to its own node — hence the 1,371 lookups, and hence the raw
pointers, because `&mut Nodes[this]` together with `&Nodes[that]` is exactly
what the borrow checker forbids.

**The root cause is that nodes own their output buffers.** Every alternative
fixes it the same way: separate value storage from expression structure.

Two more consequences worth naming:

* the 57 scalar-versus-vector sites are one concept — "is this operand a single
  value or a column" — written out 57 times;
* `undef` is a null mask, but it is a `*mut c_char` aliasing into the *tail of
  the data allocation*, so the two cannot be separated without changing the
  allocation layout.

### 10.2 The target

What this has always been is a **vectorized expression evaluator over column
batches with null masks**, row-chunked by `ffiter`. The shape is right; the
type erasure and the manual memory are not. So:

* `ColumnarValue` — `Scalar | Null | Array`, where the `Scalar` arm *is* the 57
  special cases;
* `Array` — owned `Vec` per sort, owning its own null mask;
* kernels — `fn(op, &ColumnarValue, &ColumnarValue) -> Result<ColumnarValue>`,
  sort dispatch done once instead of once per `Do_BinOp_*` variant;
* `Expr` — an owned tree evaluated recursively, so children return values and
  the aliasing problem never arises.

### 10.3 Status

The new evaluator is **wired in behind `--features new-eval`, off by default**.
`src/eval/` holds the value model, the kernels, the `Expr` tree, the
`Ast -> Expr` lowering and the bridge that writes a result back into the arena's
result node. With the feature on, **1,353 of the corpus's 1,852 expressions go
through it** and all 1,852 still match the golden file byte for byte; the rest
hit a construct the lowering does not cover yet and fall back.

The fallback is per expression and explicit: `eval::lower::lower` returns
`Unsupported("function call")`, `Unsupported("bit-string column")` and so on,
so what is missing is greppable rather than mysterious. Counting the corpus by
reason gives the order to work in:

| remaining | reason |
|---:|---|
| 132 | function call — the region and GTI tests, the shape functions (`NELEM`, `NAXIS`, `NAXES`, `ARRAY`, `AXISELEM`, `ELEMENTNUM`), `ACCUM`/`SEQDIFF`, the random generators |
| 90 | bit-string **result** (see below) |
| 6 | row offset of a string column |
| 2 | subscript slice (see below) |
| 2 | row offset written as an expression (`INTCOL{1+1}`) |
| 2 | `#SNULL` |
| 1 | string keyword |

Counting fallbacks by reason overstates what is left, because many of those
expressions do not evaluate under the old engine either. Of the 128 that once
fell back for a bit string, only 29 produce a retrievable answer; the rest are
parse errors reaching an operand, or bit-valued results the caller can never
read. **A bit-valued result is always an error** — 432 for a row-varying one,
433 for a constant — so those are deliberately left with the arena rather than
lowered just to reproduce the error, which would move the counter without
moving any real work. Bit-valued *subexpressions* are evaluated normally;
`parser::mod` checks only the top-level sort.

The same caution applies to reading the string count, for the opposite reason:
the corpus harness prints only *metadata* for a `dt=16` result, never the text,
so no corpus line arbitrates a string-valued answer. Those are covered by the
`ffcalc`/`fffrow` tests in `eval_f.rs` instead, which compare through a boolean
(`STRMID(STRCOL,1,3) == 'alp'`) or read the written column back.

What is left that does real work: the individual functions.

One structural item remains:

* **The shape functions are folded before the new lowering sees them.** The
  arena lowering turns `NELEM(V)` into a constant from the node's `nelem`;
  `eval::lower` works from the `Ast`, where that size is not known. They need
  the sizes threading through, or the fold moving earlier.

The other structural item — subscripting — is **done**. `Array` now carries
`naxes`, and `parser::mod` hands the per-column shapes to the lowering so it
can tell an element from a slice at lowering time: `M[1,1]` on a 2×3 column
lowers, `M[1]` returns `Unsupported("subscript slice")` and falls back, because
selecting a slice needs a shape-aware gather rather than one element per row.
Out-of-range subscripts are a range error at evaluation, matching `Do_Deref`.

Done: the sixteen transcendentals, `ABS`, `ARCTAN2`, two-argument `MIN`/`MAX`,
`ISNULL`, `DEFNULL`, `SETNULL`, the bitwise operators, the reductions
`SUM`, `AVERAGE`, `STDDEV`, `MEDIAN`, `NVALID` and one-argument `MIN`/`MAX`,
vector literals, fully-indexed subscripts, row offsets, the bit strings — `&`, `|`, `!`, `+`
and the six comparisons, plus `SUM` and `NVALID` over one — and the strings:
the six comparisons, `+`, `STRSTR`, `STRMID`, `NVALID`, and a conditional whose
branches are strings.

### 10.4 Row offsets, and declining a batch

`COL{k}` reads a column `k` rows from the row being evaluated, so the batch now
gathers the **whole loaded chunk** rather than just the rows being evaluated,
and a column reference is the chunk sliced at a shift of zero. Three cases:

* the row is outside the **table** — undefined, as `Do_Offset` reports it;
* the row is inside the chunk — read at the shift;
* the row is inside the table but outside the **chunk** — the file has to be
  read again, which only the engine's reload path does.

That last case is why `bridge::evaluate` now returns a bool. It declines the
batch before writing anything to the result node, and `Evaluate_Parser` walks
the arena for it instead. The fallback is therefore per *batch* as well as per
expression.

Worth knowing when changing this: **the corpus cannot reach the reload case.**
Its table is 3 rows in a single chunk, so every offset either lands in the
chunk or falls off the table. The reload path is covered by unit tests on
`Batch::shifted` alone. The same test set pins the case that made this worth
guarding — `totalRows` is 0 for an image with no `NAXIS2`, and since a plain
column reference now goes through the same path, bounding against it
unconditionally would have made every column read null.

### 10.4 What the corpus caught

Wiring it up produced four wrong answers, and only one was a slip. The other
three were design errors in the value model, which is the useful part:

1. **A value is not just scalar-or-array.** It carries *elements per row*. A
   batch is `n_rows * nelem` elements, so a scalar column and a 3-vector column
   over the same rows have lengths 3 and 9, and combining them broadcasts each
   scalar element across its row. The old engine spelled this out in every loop
   as `vector1 > 1 ? buf[elem] : buf[row]`; the first cut of `ColumnarValue`
   dropped the *branch* and the *concept* with it, so `BOOLCOL ? VECCOL : 1`
   was a length mismatch. `Array` now carries `nelem`.

2. **Integer arithmetic must not detour through `f64`.** Everything went
   through `f64`, which silently loses every bit past 2^53:
   `-9223372036854775807` came back as `-9223372036854775808`. There is now a
   separate `i64` path, with `**` the deliberate exception — the C computes it
   with `pow()` and truncates, so `(-3)**(-3)` is 0 rather than 1.

3. **`#NULL` is a rows kernel, not a folded constant.** The engine fills a
   buffer and sets every undef flag, so representing it as a scalar meant
   `anynul` never reached the caller.

4. The slip: `Expr::Column` did not copy the column's `nelem` into the array it
   built, so nothing downstream could broadcast correctly even after (1).

None of these would have been found by the unit tests, which is the argument
for wiring an incomplete evaluator in early rather than finishing it first.

### 10.5 More of what the corpus caught

Porting the elementwise functions produced four more divergences, all of them
places where the engine does something other than the obvious:

* **A domain error is a null, not a NaN.** `SQRT(-3)` sets the undef flag; the
  same goes for `LOG` and `LOG10` at or below zero and the inverse
  trigonometric functions outside [-1, 1].
* **`ROUND` is `floor(x + 0.5)`**, not Rust's `round`, which rounds half away
  from zero — they differ at `-2.5`.
* **`SETNULL`'s arguments are reordered.** The parser writes
  `SETNULL(sentinel, value)` and `New_Func` swaps them, so the value is the
  second argument by the time a kernel sees it.
* **`SETNULL` ignores its sentinel when both arguments are constant.**
  `set_null_const` copies the value and never compares, so `SETNULL(1,1)` is 1
  rather than null. CFITSIO agrees, so this is the contract rather than a
  transpilation slip — but it is almost certainly an upstream bug.
