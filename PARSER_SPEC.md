# CFITSIO Row-Filter / Calculator Expression Language — Formal Specification

This document specifies the row-filter / calculator expression language that
CFITSIO defines in `~/code/cfitsio/eval.l` and `eval.y`, and that rsfitsio
implements in `src/parser/`.

It is written to be complete enough to re-implement the parser from scratch
without reading the flex/bison sources. Every behavioural claim marked
**[verified]** was checked against the real `libcfitsio.so` via an oracle
harness (see `PARSER_MIGRATION.md` §7), not inferred from the grammar.

> **Status:** rsfitsio no longer uses flex or bison. The front end is
> hand-written safe Rust in `src/parser/` — see
> [`PARSER_MIGRATION.md`](PARSER_MIGRATION.md). **The language described here is
> unchanged**, and `tests/test_eval_corpus.rs` holds 503 expressions checked
> against output captured from the flex/bison parser before it was removed.
> Section 1 below describes the pre-migration layout and is kept as a record of
> what the generated code looked like; everything else describes the language
> as it is today.
>
> The two lexer defects in §6.1 and §6.2 are fixed in rsfitsio and submitted
> upstream as [cfitsio#153](https://github.com/HEASARC/cfitsio/pull/153).

---

## 1. Architecture

### 1.0 Today

| Layer | Rust file | Lines | Role |
|---|---|---|---|
| Public API | `src/eval_f.rs` | ~10,100 | `ffcalc`, `fffrow`, `ffsrow`, `ffcrow`, `fftexp`, `ffffrw`, `fffrwc`, `fits_pixel_filter`; `ffiprs` sets up `ParseData`, `ffcprs` tears down; column and keyword lookup callbacks. |
| Tokenizer | `src/parser/lexer.rs` | ~800 | `nom` patterns plus flex's longest-match / declaration-order rule. |
| Name resolution | `src/parser/resolve.rs` | ~100 | Resolves every name token in source order, truncating on failure. |
| Parser | `src/parser/grammar.rs` | ~640 | Pratt / precedence climbing over an untyped `Ast`. |
| Lowering | `src/parser/lower.rs` | ~1,330 | `Ast` → `Node` arena: sorts, promotions, dimension checks, function dispatch. |
| Evaluation engine | `src/eval_y.rs` | ~8,900 | `New_*` node builders and the `Do_*` per-row evaluators. Unchanged from the C. |
| Shared types | `src/eval_defs.rs`, `src/eval_tab.rs` | 245 + 103 | `ParseData`, `Node`, `lval`, `funcOp`, token ids. |

Data flow:

```
user expr (char*)
  → ffiprs: copy into ParseData.expr, append '\n'
  → parser::parse_expression
        ↳ lexer::tokenize          bytes  -> [Spanned<Tok>]
        ↳ resolve::resolve_names   names  -> (token kind, value), in source order
        ↳ grammar::parse           tokens -> Ast
        ↳ lower::Lowerer::lower    Ast    -> Nodes[], resultNode
  → Evaluate_Parser(lParse, firstRow, nRows) walks Nodes bottom-up per row-chunk
```

### 1.1 Before the migration (historical)

| Layer | C source | Rust file | Lines | Role |
|---|---|---|---|---|
| Front end | `eval_f.c` | `src/eval_f.rs` | 10,160 | Public API (`ffcalc`, `fffrow`, `ffsrow`, `ffcrow`, `fftexp`, `ffffrw`, `fffrwc`, `fits_pixel_filter`); `ffiprs` sets up `ParseData`, `ffcprs` tears down; column/keyword lookup callbacks. |
| Lexer | `eval.l` → `eval_l.c` | `src/eval_l.rs` | 1,688 | Flex DFA + buffer machinery. |
| Parser + evaluator | `eval.y` → `eval_y.c` | `src/eval_y.rs` | 15,305 | Bison LALR tables + driver + semantic actions (build `Node` tree), then the whole run-time evaluation engine. |
| Shared types | `eval_defs.h`, `eval_tab.h` | `src/eval_defs.rs`, `src/eval_tab.rs` | 248 + 103 | `ParseData`, `Node`, `lval`, `funcOp`, token ids. |

#### Breakdown of the old `src/eval_y.rs`

| Line range | Contents | Fate under migration |
|---|---|---|
| 1–257 | libm shims, `yysymbol_kind_t` constants | delete |
| 258–525 | `YYTRANSLATE`, `YYPACT`, `YYDEFACT`, `YYPGOTO`, `YYDEFGOTO`, `YYTABLE` (1777), `YYCHECK` (1777), `YYSTOS`, `YYR1` | **delete** (LALR tables) |
| 526–833 | `Alloc_Node`, `Free_Last_Node`, `New_Const`, `New_Column`, `New_Offset`, `New_Unary`, `New_BinOp`, `New_Func` | **keep** |
| 834–1029 | `yysymbol_name`, `yy_symbol_print`, `yy_stackprint`, `yy_reduce_print`, `yydestruct` | **delete** (bison debug scaffolding) |
| 1030–1119 | `New_FuncSize` | **keep** |
| 1120–6964 | `fits_parser_yyparse` — LALR driver + 135 inlined reduce actions (~5,845 lines) | **replace** |
| 6965–7989 | `New_Deref`, `New_GTI`, `New_REG`, `New_Vector`, `Close_Vec`, `New_Array`, `Locate_Col`, `Test_Dims`, `Copy_Dims` | **keep** |
| 7990–15305 | `Evaluate_Parser`, `Evaluate_Node`, `Allocate_Ptrs`, `Do_Unary`, `Do_Offset`, `Do_BinOp_*`, `Do_Func`, `Do_Deref`, `Do_GTI`, `Do_REG`, `Do_Vector`, `Do_Array`, `bitand/bitor/bitnot/bitcmp/bitlgte`, `saobox`, `ellipse`, `circle`, `bnear`, `cstrmid`, GTI search | **keep** — the evaluation engine, ~7,300 lines |

**Replaceable surface = 1,688 (all of `eval_l.rs`) + ~6,300 (`eval_y.rs` tables +
driver + debug) ≈ 8,000 lines.** Everything else stays.

#### Old data flow

```
user expr (char*)
  → ffiprs: copy into ParseData.expr, append '\n', index = 0
  → fits_parser_yylex_init_extra(lParse) / yyrestart
  → fits_parser_yyparse(scanner, lParse)
        ↳ fits_parser_yylex()  ──── YY_INPUT → expr_read() reads ParseData.expr byte-wise
              ↳ fits_parser_yyGetVariable() → find_variable() → lParse->getData
                 (= find_column in eval_f.rs, which also handles '#' keywords)
        ↳ reduce actions → New_*() → lParse.Nodes[] (flat arena, i32 indices)
  → lParse.resultNode = index of root
  → Evaluate_Parser(lParse, firstRow, nRows) walks Nodes bottom-up per row-chunk
```

The lexer is **not** context-free with respect to the host FITS file: token
classification of an identifier requires a column/keyword lookup. This is the
single most important coupling in the design.

---

## 2. Lexical specification

### 2.1 Input model

`expr_read()` feeds the flex buffer from `ParseData.expr` one chunk at a time,
stopping at the first NUL. `ffiprs` appends a `'\n'` before the NUL, so every
expression is terminated by exactly one newline, which is the grammar's
statement terminator.

> Note: `expr_read` writes `buf[n] = '\0'` where `n` may equal `nbytes`,
> i.e. one byte past the requested size. Flex over-allocates its buffer by 2,
> so this is benign in practice, but it is a real out-of-bounds-by-design in
> the C and should not be reproduced.

### 2.2 Matching discipline

Standard flex: **longest match wins; ties broken by rule order in the file.**
Both halves matter and produce several surprising results (§2.6).

### 2.3 Pattern definitions (verbatim from `eval.l`)

```
bit         ([bB][01xX]+)
oct         ([oO][01234567xX]+)
hex         ([hH][0123456789aAbBcCdDeEfFxX]+)
bitconst    (0b[01]+)
hexconst    (0x[0123456789aAbBcCdDeEfF]+)
octconst    (0o[01234567]+)
integer     [0-9]+
boolean     (t|f|T|F)
real        ([0-9]*"."[0-9]+)|([0-9]*"."*[0-9]+[eEdD][+-]?[0-9]+)|([0-9]+".")
constant    ("#"[a-zA-Z0-9_]+)|("#""$"[^\n]*"$")
string      ([\"][^\"\n]*[\"])|([\'][^\'\n]*[\'])
variable    ([a-zA-Z_][a-zA-Z0-9_]*)|("$"[^$\n]*"$")
function    [a-zA-Z][a-zA-Z0-9]+"("
intcast     ("(int)"|"(INT)")
fltcast     ("(float)"|"(FLOAT)"|"(double)"|"(DOUBLE)")
power       ("^"|"**")
not         ("!"|".not."|".NOT."|"not."|"NOT.")
or          ("||"|".or."|".OR."|"or."|"OR.")
and         ("&&"|".and."|".AND."|"and."|"AND.")
equal       ("=="|".eq."|".EQ."|"eq."|"EQ.")
not_equal   ("!="|".ne."|".NE."|"ne."|"NE.")
greater     (">"|".gt."|".GT."|"gt."|"GT.")
lesser      ("<"|".lt."|".LT."|"lt."|"LT.")
greater_eq  (">="|"=>"|".ge."|".GE."|"ge."|"GE.")
lesser_eq   ("<="|"=<"|".le."|".LE."|"le."|"LE.")
xor         ("^^"|".xor."|".XOR.")
nl          \n
```

### 2.4 Rules, in file order (order is semantically significant)

| # | Pattern | Token emitted | Semantic value / action |
|---|---|---|---|
| 1 | `[ \t]+` | *(skipped)* | — |
| 2 | `{bit}` | `BITSTR` | `str` = text after leading `b`/`B`. Errors `PARSE_SYNTAX_ERR` if length ≥ `MAX_STRLEN` (256). |
| 3 | `{oct}` | `BITSTR` | each octal digit expanded to 3 bits (`x`→`xxx`); overflow past 255 bits → error. |
| 4 | `{hex}` | `BITSTR` | each hex digit expanded to 4 bits (`x`→`xxxx`); overflow past 255 bits → error. |
| 5 | `{bitconst}` | `LONG` | `0b…` parsed base-2 into `lng`. |
| 6 | `{octconst}` | `LONG` | `0o…` parsed base-8 into `lng`. |
| 7 | `{hexconst}` | `LONG` | `0x…` parsed base-16, case-insensitively. *(Was miscomputed for uppercase A–F before the fix — see §6.1.)* |
| 8 | `{integer}` | `LONG` | `atol`. No overflow check. |
| 9 | `{boolean}` | `BOOLEAN` | `log` = 1 for `t`/`T`, 0 for `f`/`F`. |
| 10 | `{real}` | `DOUBLE` | `atof` (accepts Fortran `d`/`D` exponent because `atof` sees… see §6.2). |
| 11 | `{constant}` | varies | `#PI`→`DOUBLE` 4·atan(1); `#E`→`DOUBLE` e; `#DEG`→`DOUBLE` π/180; `#ROW`→`ROWREF`; `#NULL`→`NULLREF`; `#SNULL`→`SNULLREF`; otherwise `#$name$`/`#name` is rewritten to `#name` and passed to `getData` (keyword lookup). Case-insensitive compare. |
| 12 | `{string}` | `STRING` | quotes stripped; length ≥ 256 → error, value becomes `""`. |
| 13 | `{variable}` | *(from lookup)* | `$…$` form has the `$` delimiters stripped. Then `fits_parser_yyGetVariable` → §2.5. |
| 14 | `{function}` | `BFUNCTION` / `GTIFILTER` / `GTIOVERLAP` / `GTIFIND` / `REGFILTER` / `IFUNCTION` / `FUNCTION` | Name **upper-cased** (including the trailing `(`) into `str`. Mapping in §2.7. |
| 15 | `{intcast}` | `INTCAST` | |
| 16 | `{fltcast}` | `FLTCAST` | |
| 17 | `{power}` | `POWER` | |
| 18 | `{not}` | `NOT` | |
| 19 | `{or}` | `OR` | |
| 20 | `{and}` | `AND` | |
| 21 | `{equal}` | `EQ` | |
| 22 | `{not_equal}` | `NE` | |
| 23 | `{greater}` | `GT` | |
| 24 | `{lesser}` | `LT` | |
| 25 | `{greater_eq}` | `GTE` | |
| 26 | `{lesser_eq}` | `LTE` | |
| 27 | `{xor}` | `XOR` | |
| 28 | `{nl}` | `'\n'` (literal 10) | statement terminator |
| 29 | `.` | the character itself | covers `( ) [ ] { } , : ? + - * / % & | ~ =` and anything else |

Tokens `ACCUM` (291) and `DIFF` (292) are declared but **never produced by the
lexer and never appear in any grammar rule** — they exist only as opcode
constants passed to `New_BinOp` by the `ACCUM(`/`SEQDIFF(` function actions.
Bison reports them as "Terminals unused in grammar". **[verified]**

### 2.5 Identifier resolution (the semantic callback)

```
fits_parser_yyGetVariable(lParse, name, &lval):
    varNum = find_variable(lParse, name)      # case-insensitive over first 80 chars
    if varNum < 0:
        return lParse->getData(lParse, name, &lval)   # = find_column()
    else:
        lval.lng = varNum
        switch varData[varNum].type:
            LONG | DOUBLE → COLUMN
            BOOLEAN       → BCOLUMN
            STRING        → SCOLUMN
            BITSTR        → BITCOL
            default       → pERROR (-1), status = PARSE_SYNTAX_ERR
```

`find_column` (in `eval_f.rs`) additionally:
- routes names starting with `#` to `find_keywd` (header keyword lookup, which
  can yield a literal `LONG` / `DOUBLE` / `BOOLEAN` / `STRING` token);
- for `IMAGE_HDU` inputs, resolves against `ParseData.pixFilter` tags;
- registers the column in `varData[]` and `colData[]` on first use.

So the token returned for an identifier is one of:
`COLUMN`, `BCOLUMN`, `SCOLUMN`, `BITCOL`, `LONG`, `DOUBLE`, `BOOLEAN`,
`STRING`, or `pERROR`.

### 2.6 Verified lexical quirks

All **[verified]** against a standalone flex build of the exact `eval.l`
patterns.

| Input | Tokens | Why |
|---|---|---|
| `b101` | `BITSTR("101")` | `bit` matches 4 chars; `variable` also 4; `bit` rule comes first. |
| `b1x0y` | `VAR("b1x0y")` | `bit` matches only 4 (`b1x0`), `variable` matches 5 → longest wins. |
| `T` | `BOOLEAN(1)` | tie at length 1, `boolean` rule precedes `variable`. |
| `TT` | `VAR("TT")` | `variable` is longer. |
| `false` | `VAR("false")` | not a boolean literal — only single letters are. |
| `.` | `'.'` (catch-all) → syntax error | `real` alt 3 is `[0-9]+"."`, so a bare dot is not a number. Was `DOUBLE(0.0)` before the fix — see §6.2. |
| `..` | `'.'` `'.'` → syntax error | |
| `1.2.3` | `DOUBLE(1.2)` `DOUBLE(.3)` | |
| `1..2` | `DOUBLE(1.)` `DOUBLE(.2)` | |
| `.5.5` | `DOUBLE(.5)` `DOUBLE(.5)` | alt 1 keeps the leading-dot form. |
| `1.e5` | `DOUBLE(1.)` `VAR("e5")` | exponent alt requires a digit immediately before `[eEdD]`. |
| `.e5` | `'.'` `VAR("e5")` → syntax error | was `DOUBLE(.)` `VAR("e5")` before the fix. |
| `1.5e` | `DOUBLE(1.5)` `VAR("e")` | |
| `0x` / `0b` / `0o` | `LONG(0)` `VAR("x"/"b"/"o")` | the `…const` patterns need ≥1 digit. |
| `12abc` | `LONG(12)` `VAR("abc")` | |
| `oct123` | `VAR("oct123")` | `oct` can't extend past `o`, so `variable` wins. |
| `h1F` | `BITSTR("00011111")` | hex bit-string. |
| `X(1)` | `VAR("X")` `'('` `LONG(1)` `')'` | `function` needs **≥2** name chars. |
| `ab(` | `FUNCTION("AB(")` | |
| `a (` | `VAR("a")` `'('` | no whitespace allowed before `(`. |
| `(int )x` | `'('` `VAR("int")` `')'` `VAR("x")` | cast spelling is exact. |
| `a.Ne.b` | `VAR("a")` `'.'` `VAR("Ne")` `'.'` `VAR("b")` | mixed case not accepted — only all-lower or all-upper. |
| `a***b` | `VAR` `POWER` `'*'` `VAR` | |
| `A=>B` | `A` `GTE` `B` | `=>` is an alias for `>=`. |
| `A=<B` | `A` `LTE` `B` | |
| `#$a$b$` | `CONST("#$a$b$")` | greedy `[^\n]*"$"`. |
| `$a$b$` | `VAR("$a$")` `VAR("b")` `'$'` | non-greedy via `[^$\n]*`. |
| `"unterminated` | `'"'` `VAR("unterminated")` | unterminated strings degrade silently. |
| `not .a` | `VAR("not")` `'.'` `VAR("a")` | no space allowed inside `.not.`. |

---

## 3. Grammar

Bison reports **0 shift/reduce and 0 reduce/reduce conflicts**, 322 LALR
states, 135 rules. The grammar is unambiguous LALR(1). **[verified]**

### 3.1 Precedence, lowest → highest

| Level | Assoc | Tokens |
|---:|---|---|
| 1 | left | `,` `=` `:` `{` `}` |
| 2 | right | `?` |
| 3 | left | `OR` |
| 4 | left | `AND` |
| 5 | left | `EQ` `NE` `~` |
| 6 | left | `GT` `LT` `LTE` `GTE` |
| 7 | left | `+` `-` `%` |
| 8 | left | `*` `/` |
| 9 | left | `\|` `&` `XOR` |
| 10 | right | `POWER` |
| 11 | left | `NOT` |
| 12 | left | `INTCAST` `FLTCAST` |
| 13 | left | `UMINUS` |
| 14 | left | `[` |
| 15 | right | `ACCUM` `DIFF` *(unused)* |

**Two levels differ sharply from C/Fortran intuition and must be preserved:**

1. **Bitwise `| & XOR` bind *tighter* than `* /`.**
   `1|2*3` = `(1|2)*3` = 9 **[verified]**; `1+2|4` = `1+(2|4)` = 7 **[verified]**.
2. **Unary minus binds *tighter* than `**`/`^`.**
   `-2**2` = `(-2)**2` = **4**, not −4 **[verified]**.
   `2**3**2` = 512 (right-associative) **[verified]**.

Casts bind tighter than `*`: `(int)2.7*2` = 4, `(int)(2.7*2)` = 5 **[verified]**.

### 3.2 Nonterminals

The grammar is **type-stratified**: four mutually-recursive nonterminals, one
per value sort.

| Nonterminal | Value sort | `lval` payload |
|---|---|---|
| `expr` | numeric (`LONG` or `DOUBLE`) and numeric vectors | `Node` index |
| `bexpr` | boolean scalar/vector | `Node` index |
| `sexpr` | string | `Node` index |
| `bits` | bit-string | `Node` index |
| `vector` / `bvector` | partially-built `{…}` literals | `Node` index |

Stratification is what makes the grammar LALR(1) without conflicts: the *type*
of an operand selects the production, so `a + b` has three distinct rules
(`expr+expr`, `sexpr+sexpr`, `bits+bits`) instead of one ambiguous one.

### 3.3 Productions

```
lines   : ε | lines line

line    : '\n'
        | expr  '\n'      → resultNode = $1
        | bexpr '\n'      → resultNode = $1
        | sexpr '\n'      → resultNode = $1
        | bits  '\n'      → resultNode = $1
        | error '\n'      → yyerrok

bvector : '{' bexpr | bvector ',' bexpr
vector  : '{' expr  | vector ',' expr | vector ',' bexpr | bvector ',' expr

expr    : vector  '}'
bexpr   : bvector '}'
```

A `{…}` literal is a **boolean** vector iff *every* element is boolean;
the first non-boolean element promotes the whole literal to numeric.
**[verified]** `{B,B}`→bool, `{1,B}`→long, `{B,1}`→long.
Elements are chunked at `MAXSUBS` = 10: on the 11th element the current vector
node is closed with `Close_Vec` and nested inside a fresh `New_Vector`.

```
bits    : BITSTR                                       (size = strlen)
        | BITCOL
        | BITCOL '{' expr '}'                          (offset; expr must be constant LONG)
        | bits '&' bits                                (size = max)
        | bits '|' bits                                (size = max)
        | bits '+' bits                                (concat; size = sum, must be < 256)
        | bits '[' expr (',' expr){0,4} ']'            (1–5 dimensional deref)
        | NOT bits
        | '(' bits ')'

sexpr   : STRING                                       (size = strlen)
        | SCOLUMN
        | SCOLUMN '{' expr '}'
        | SNULLREF
        | '(' sexpr ')'
        | sexpr '+' sexpr                              (concat; size = sum, must be < 256)
        | bexpr '?' sexpr ':' sexpr                    (condition must be scalar)
        | FUNCTION sexpr ',' sexpr ')'                 (DEFNULL)
        | FUNCTION sexpr ',' expr ',' expr ')'         (STRMID)

expr    : LONG | DOUBLE | COLUMN | ROWREF | NULLREF
        | COLUMN '{' expr '}'
        | expr ('%'|'+'|'-'|'*'|'/') expr              (with numeric promotion)
        | expr ('&'|'|'|XOR) expr                      (both operands must be LONG)
        | expr POWER expr
        | ('+'|'-') expr            %prec UMINUS
        | '(' expr ')'
        | expr '*' bexpr | bexpr '*' expr              (bool coerced to numeric)
        | bexpr '?' expr  ':' expr
        | bexpr '?' bexpr ':' expr
        | bexpr '?' expr  ':' bexpr
        | expr '[' expr (',' expr){0,4} ']'
        | INTCAST expr | INTCAST bexpr
        | FLTCAST expr | FLTCAST bexpr
        | FUNCTION ')'                                 (see §4)
        | FUNCTION (expr|bexpr|sexpr|bits) ')'
        | FUNCTION expr ',' expr ')'
        | FUNCTION bexpr ',' expr ')'
        | FUNCTION expr ',' expr ',' expr ',' expr ')'
        | IFUNCTION sexpr ',' sexpr ')'
        | GTIOVERLAP STRING ',' expr ',' expr ')'
        | GTIOVERLAP STRING ',' expr ',' expr ',' STRING ',' STRING ')'

bexpr   : BOOLEAN | BCOLUMN
        | BCOLUMN '{' expr '}'
        | bits  (EQ|NE|LT|LTE|GT|GTE) bits             (result size forced to 1)
        | expr  (GT|LT|GTE|LTE|'~'|EQ|NE) expr
        | sexpr (EQ|NE|GT|GTE|LT|LTE) sexpr            (result size forced to 1)
        | bexpr (AND|OR|EQ|NE) bexpr
        | expr '=' expr ':' expr                       (range test, see below)
        | bexpr '?' bexpr ':' bexpr
        | BFUNCTION (expr|bexpr|sexpr) ')'             (ISNULL)
        | FUNCTION bexpr ',' bexpr ')'                 (DEFNULL)
        | BFUNCTION expr (',' expr){2}  ')'            (NEAR)
        | BFUNCTION expr (',' expr){4}  ')'            (CIRCLE)
        | BFUNCTION expr (',' expr){6}  ')'            (BOX, ELLIPSE)
        | GTIFILTER ')' | GTIFILTER STRING ')'
        | GTIFILTER STRING ',' expr ')'
        | GTIFILTER STRING ',' expr ',' STRING ',' STRING ')'
        | GTIFIND  ')' | GTIFIND STRING ')'
        | GTIFIND  STRING ',' expr ')'
        | GTIFIND  STRING ',' expr ',' STRING ',' STRING ')'
        | REGFILTER STRING ')'
        | REGFILTER STRING ',' expr ',' expr ')'
        | REGFILTER STRING ',' expr ',' expr ',' STRING ')'
        | bexpr '[' expr (',' expr){0,4} ']'
        | NOT bexpr
        | '(' bexpr ')'
```

### 3.4 Structural constraints that are *syntactic*, not semantic

These fall out of the stratified grammar and must be reproduced by any
replacement, because they determine *accept vs. reject*, not just diagnostics:

1. **`{…}` offset applies only to a bare column token**, never to a
   parenthesised or computed expression.
   `(X){1}` and `(X+1){1}` are syntax errors. **[verified]**
2. **Deref `[…]` is defined for `expr`, `bexpr`, `bits` — not `sexpr`.**
   `S[1]` is a syntax error. **[verified]**
3. **`NOT` applies to `bexpr` and `bits` only** — there is no `NOT expr`.
4. **GTI/REG filter filenames must be literal `STRING` tokens**, not string
   expressions. `GTIFILTER(S)` where `S` is a string column is a syntax error.
   **[verified]**
5. **There is no `bexpr '?' bits ':' bits`.** `B?BIT:BIT` is a syntax error.
   **[verified]**
6. **There is no `expr '+' bexpr`.** `X+B` and `B+B` are syntax errors, but
   `X*B` and `B*X` are accepted (dedicated rules). **[verified]**
7. **`ISNULL` is not defined for `bits`.** `ISNULL(BIT)` is a syntax error.
   **[verified]**
8. **`DEFNULL(bits,bits)` does not exist.** **[verified]**
9. **Range tests do not chain.** `X=1:2` is fine; `X=1:2:3` is a syntax error.
   **[verified]** A bare `1:2` is also a syntax error. **[verified]**

### 3.5 Type-directed parsing of `NOT` — the one genuinely hard case

`NOT` sits at precedence 11 (above everything except casts, `UMINUS`, `[`),
but its only rules are `NOT bexpr` and `NOT bits`. When the token after `!`
starts a *numeric* expression the parser **cannot** reduce, so it keeps
shifting until a boolean has been formed. Consequences (X = 7):

| Input | Parse | Value |
|---|---|---|
| `!B \|\| B` | `(!B) \|\| B` | 1 **[verified]** (`!(B\|\|B)` would be 0) |
| `!B && 1>2` | `(!B) && (1>2)` | 0 **[verified]** (`!(B&&1>2)` is 1) |
| `!X > 5` | `!(X > 5)` | 0 **[verified]** |
| `!X > 10` | `!(X > 10)` | 1 **[verified]** |
| `!X + 1 > 5` | `!((X+1) > 5)` | 0 **[verified]** |
| `!X > 5 \|\| B` | `(!(X>5)) \|\| B` | 1 **[verified]** |

So `NOT` binds at level 11 *over boolean operands*, but "reaches through" any
numeric prefix until the first sub-expression that is boolean-valued. A naive
Pratt parser that binds `!` to the tightest operand would reject `!X>5`.
This is the single behaviour most likely to be broken by a straightforward
rewrite; it needs an explicit rule (see `PARSER_MIGRATION.md` §4.4).

### 3.6 `?:` and the range operator `=`/`:`

`=` and `:` sit at level 1, *below* `?` at level 2. The range form
`expr '=' expr ':' expr` desugars to `(lo <= x) AND (x <= hi)`.

Interaction is exercised and unambiguous:
- `B?1=2:3:4` → `B ? (1=2:3) : 4` → LONG **[verified]**
- `B?X=1:2:Y` → `B ? (X=1:2) : Y` → DOUBLE **[verified]**
- `B?B?1:2:3` → `B ? (B?1:2) : 3` → 1 **[verified]**

---

## 4. Function catalogue

The lexer classifies the *upper-cased* function name into one of six tokens:

| Lexer test (on upper-cased name incl. `(`) | Token |
|---|---|
| `BOX(`, `CIRCLE(`, `ELLIPSE(`, `NEAR(`, `ISNULL(` | `BFUNCTION` |
| `GTIFILTER(` | `GTIFILTER` |
| `GTIOVERLAP(` | `GTIOVERLAP` |
| `GTIFIND(` | `GTIFIND` |
| `REGFILTER(` | `REGFILTER` |
| `STRSTR(` | `IFUNCTION` |
| anything else | `FUNCTION` |

Name matching is therefore **case-insensitive** everywhere. Dispatch on the
specific name happens in the *grammar action*, by `FSTRCMP` against the
upper-cased text — so an unknown name is a **parse-time** error
("Function(expr) not supported"), not a lexer error.

### 4.1 `FUNCTION` — dispatch by (arity, argument sorts)

| Arity / arg sorts | Names | Result |
|---|---|---|
| `()` | `RANDOM`, `RANDOMN` | `DOUBLE` scalar |
| `(bexpr)` | `SUM`→LONG, `NELEM`→const LONG, `ACCUM`→LONG | |
| `(sexpr)` | `NELEM`→const LONG, `NVALID`→LONG | |
| `(bits)` | `NELEM`, `NVALID` (= size), `SUM`, `MIN`, `MAX` (size forced to 1), `ACCUM` | |
| `(expr)` | `SUM`, `AVERAGE`→DOUBLE, `STDDEV`→DOUBLE, `MEDIAN`, `NELEM`, `NVALID`→LONG, `ACCUM`, `SEQDIFF`, `ABS`, `MIN`, `MAX`, `RANDOM`→DOUBLE, `RANDOMN`→DOUBLE, `ELEMENTNUM`→LONG, `NAXIS`→const LONG | |
| `(expr)` coerced to DOUBLE first | `SIN` `COS` `TAN` `ARCSIN`/`ASIN` `ARCCOS`/`ACOS` `ARCTAN`/`ATAN` `SINH` `COSH` `TANH` `EXP` `LOG` `LOG10` `SQRT` `ROUND` `FLOOR` `CEIL` `RANDOMP`(→LONG) | `DOUBLE` |
| `(bexpr, expr)` | `AXISELEM`, `NAXES`, `ARRAY` | |
| `(expr, expr)` | `DEFNULL`, `ARCTAN2`, `MIN`, `MAX`, `SETNULL`, `AXISELEM`, `NAXES`, `ARRAY` | |
| `(expr, expr, expr, expr)` | `ANGSEP` | `DOUBLE` |
| `(bexpr, bexpr)` | `DEFNULL` | `BOOLEAN` |
| `(sexpr, sexpr)` | `DEFNULL` | `STRING` |
| `(sexpr, expr, expr)` | `STRMID` | `STRING` |

`NAXIS`, `NAXES`, `NELEM` are **constant-folded at parse time** from the
argument node's `naxis` / `naxes[]` / `nelem` metadata. `AXISELEM`/`NAXES`
require their second argument to be a scalar constant.

### 4.2 `BFUNCTION`

| Signature | Name | Notes |
|---|---|---|
| `(expr)` / `(bexpr)` / `(sexpr)` | `ISNULL` | result `BOOLEAN`, size = argument size |
| `(expr, expr, expr)` | `NEAR(x, y, tol)` | all coerced to DOUBLE, dims must match |
| `(expr ×5)` | `CIRCLE(xcen, ycen, rad, xcol, ycol)` | |
| `(expr ×7)` | `BOX(xcen, ycen, xwid, ywid, rot, xcol, ycol)`, `ELLIPSE(xcen, ycen, xrad, yrad, rot, xcol, ycol)` | |

### 4.3 `IFUNCTION`

`STRSTR(sexpr, sexpr)` → `LONG` (1-based position, 0 if not found).

### 4.4 GTI / region filters

| Form | Meaning |
|---|---|
| `GTIFILTER()` | defaults: file `""`, time expr default, `*START*`/`*STOP*` columns |
| `GTIFILTER(f)` | `f` must be a literal `STRING` |
| `GTIFILTER(f, t)` | |
| `GTIFILTER(f, t, start, stop)` | `start`, `stop` literal `STRING`s |
| `GTIFIND(…)` | identical four shapes; returns the GTI index |
| `GTIOVERLAP(f, tstart, tstop)` | |
| `GTIOVERLAP(f, tstart, tstop, start, stop)` | |
| `REGFILTER(f)` | |
| `REGFILTER(f, x, y)` | |
| `REGFILTER(f, x, y, colnames)` | `colnames` literal `STRING` |

Sentinel `-99` is passed for omitted node arguments.

---

## 5. Type system

### 5.1 Sorts and promotion

Token ids double as type tags, ordered for promotion:

```
BOOLEAN (258) < LONG (259) < DOUBLE (260)      # numeric promotion ladder
STRING  (261), BITSTR (262)                    # not on the ladder
```

`PROMOTE(a,b)` inserts a `New_Unary(type, 0, …)` conversion on whichever side
has the lower tag. This is why the `%token` declaration order in `eval.y`
carries the comment *"First 3 must be in order of increasing promotion"*.

### 5.2 Dimensions

Every `Node` carries `value.nelem`, `value.naxis`, `value.naxes[MAXDIMS=5]`.
`Test_Dims(a,b)` accepts when either side is a scalar (`nelem == 1`) or the
shapes match exactly; `Copy_Dims` propagates the larger shape to the result.
Comparison operators on `bits` and `sexpr` force result `nelem = 1`.

### 5.3 Sizes

`MAX_STRLEN` = 256 (255 usable). String concat, bit concat and `STRMID`
lengths are all checked against it at parse time.

---

## 6. Defects

§6.1 and §6.2 are **fixed** in rsfitsio and submitted upstream as
[cfitsio#153](https://github.com/HEASARC/cfitsio/pull/153), which also adds `test_ffcrow_hex_constant_case` and
`test_fftexp_bare_dot_rejected` to `tests/test_eval.c`. §6.3–§6.7 described the *generated* code and no longer
apply to rsfitsio, which has none of it; they are retained because they still
describe upstream CFITSIO.

### 6.1 Uppercase hex literals were miscomputed — **fixed** **[verified]**

`eval.l` rule 7 read:
```c
int v = (isdigit(*p) ? (*p - '0') : (*p - 'a' + 10));
```
There was no `tolower`, so `'F'` yielded `70 - 97 + 10 = -17`, and the
negative digit corrupted every higher bit of the accumulator.

| Input | Before | After / correct |
|---|---|---|
| `0x1f` | 31 | 31 |
| `0x1F` | **−1** | 31 |
| `0xff` | 255 | 255 |
| `0xFF` | **−1** | 255 |
| `0xA`  | **−22** | 10 |
| `0xABCDEF` | garbage | 11259375 |

Now reads:
```c
int c = (unsigned char)*p;
int v = (isdigit(c) ? (c - '0') : (tolower(c) - 'a' + 10));
```

In rsfitsio this is `Rule::HexConst` in `src/parser/lexer.rs`, which folds
case via `to_digit(16)`.

### 6.2 A bare `.` lexed as `DOUBLE 0.0` — **fixed** **[verified]**

`real` alternative 3 was `[0-9]*"."`, with the integer part optional, so `.`
matched and `..` yielded two `DOUBLE` tokens. It is now `[0-9]+"."`, so a lone
dot falls through to the `.` catch-all rule and the parser rejects it with
`PARSE_SYNTAX_ERR`. `1.`, `.5`, `12.`, `1.5e3` are unaffected (§2.6).

In rsfitsio this is the `trailing_dot` branch of `m_real` in
`src/parser/lexer.rs`, which requires `digit1` before the point.

### 6.3 Dead trailing-space stripping in bit/oct/hex rules

Rules 2–4 contain `while (yytext[len] == ' ') len--;`, but none of the patterns
can match a space, so the loop never executes. Harmless, but it is a decrement
loop with no lower bound — a latent underflow if the patterns ever change.

### 6.4 Fortran `d`/`D` exponents rely on `atof`

`real` accepts `1.5d-3`, but `atof` does not understand `d`. Whether this
works is libc-dependent; the Rust transpile must be checked against it.
*(Not yet verified in isolation — flagged for the migration test matrix.)*

### 6.5 `expr_read` writes one byte past the requested size

See §2.1.

### 6.6 `alloca` is `#define`d to `malloc` in `eval.y`

Line 92: `#define alloca malloc` — every "stack" allocation in the bison
skeleton is a heap allocation that is never freed on the error paths. The
Rust transpile inherits `malloc` calls in `fits_parser_yyparse`
(`src/eval_y.rs:~1180`). Removing the LALR driver removes this entirely.

### 6.7 `#define YYINITDEPTH 100` with doubling `malloc` growth

Deeply nested expressions reallocate the parse stack via raw `malloc`/`memcpy`
in the transpiled code, with `YYMAXDEPTH` as the ceiling. A recursive-descent
replacement needs its own depth limit to avoid stack overflow (§`PARSER_MIGRATION.md` §4.6).

### 6.8 `STRSTR` reports a miss two different ways

`str_pos_rows` marks the row **undefined** when the needle does not occur, but
`str_pos_const` — the folded form used when both arguments are constants —
stores a plain **0** instead. So `STRSTR(STRCOL,'z')` is null while
`STRSTR('abc','z')` is `0`, and `ISNULL(STRSTR(...))` answers differently
depending only on whether the arguments happen to be literals.

Found while porting the string kernels; not reported upstream. The port
reproduces both behaviours.

---

## 7. Reference: token id table

From `src/eval_tab.rs` (`fits_parser_yytokentype`):

```
BOOLEAN 258   LONG 259    DOUBLE 260   STRING 261   BITSTR 262
FUNCTION 263  BFUNCTION 264  IFUNCTION 265
GTIFILTER 266 GTIOVERLAP 267 GTIFIND 268 REGFILTER 269
COLUMN 270    BCOLUMN 271 SCOLUMN 272  BITCOL 273
ROWREF 274    NULLREF 275 SNULLREF 276
OR 277  AND 278  EQ 279  NE 280  GT 281  LT 282  LTE 283  GTE 284
XOR 285  POWER 286  NOT 287  INTCAST 288  FLTCAST 289  UMINUS 290
ACCUM 291  DIFF 292        # never lexed, never in a rule
```

`funcOp` opcodes (`src/eval_defs.rs`) start at `rnd_fct = 1001` and are the
values stored in `Node.operation` for function nodes; binary/unary operators
store the token id or the raw character.
