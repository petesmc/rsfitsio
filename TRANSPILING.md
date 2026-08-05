# Transpiling CFITSIO C into safe Rust

This crate is a port of CFITSIO. The standing task is to replace `todo!()` /
raw-pointer stubs with safe Rust transpiled from the original C. This guide
captures the conventions so a transpile lands consistently.

## The workflow

A request usually looks like: *"I've pasted the C body of `ffxyz` into `ffxyz_safe`,
transpile it."* The pasted C does not compile, so:

1. **Read** the pasted C and the function's current signature.
2. **Find the safe signatures** of every function it calls (see "Finding dependencies").
3. **Replace the body.** The C often contains tabs; matching it exactly with `Edit`
   is fragile. Prefer deleting the body by line range, then inserting Rust:
   ```bash
   # find the `) -> c_int {` line (S) and the closing `}` line (E)
   sed -i 'S+1,E-1d' src/<file>.rs        # delete just the body
   ```
   then `Edit` the now-empty `sig { }` into `sig { <rust body> }`. The signature
   lines are space-indented and match reliably. When deleting several bodies, work
   **bottom-to-top** so earlier line numbers stay valid.
4. **Build and test** (`cargo build --lib`, `cargo test --lib`). Fix callers whose
   signatures you changed. Run the full `cargo test` before declaring done.

### Rules of engagement

- **Keep the structure and the original `/* ... */` comments.** Adding clarifying
  comments is encouraged; deleting the C's comments is not.
- **Keep the C variable names verbatim** (`xtensionCol`, `ttypeBuff`, …). The crate
  has `#![allow(non_snake_case)]`, so camelCase locals are fine.
- **Use the C integer types** (`c_int`, `c_long`, `c_uchar`, `LONGLONG`), not `i32`/`i64`.
- A transpiled function may still depend on `todo!()` helpers; it will `panic!` at
  runtime until they're implemented. That's expected and matches its siblings — the
  goal is a faithful, compiling transpile with a coherent safe signature.

## Naming & the wrapper convention

Every public CFITSIO function `ffxyz` exists in three forms:

- `ffxyz` — `#[no_mangle] pub unsafe extern "C" fn`, the C ABI entry point. It only
  null-checks the raw pointers and delegates.
- `ffxyz_safe` — the safe implementation (`pub fn`, references/slices/`Option`s). This
  is what you transpile into. (`_safer` is an older name for one that still carries a
  raw pointer in its signature; prefer making it fully `_safe`.)
- A `rust_api` / `c_api` alias in `aliases.rs` (e.g. `fits_open_file`).

In `group.rs` and similar files, `use crate::aliases::rust_api::*;` is in scope, so
call dependencies by their `fits_*` names — which conveniently match the C.

## Finding dependencies

For each `fits_foo(...)` / `ffxyz(...)` call in the C:

```bash
grep -nE 'as fits_foo\b' src/aliases.rs          # rust_api line -> the _safe target
grep -nE 'pub(\(crate\))? fn ffxyz_safe' src/*.rs # then read its real signature
```

Read the safe signature before writing the call — parameter shapes vary (see below),
and some stubs have **wrong** provisional signatures that you must fix.

## Pointer / type translations

| C | Rust (safe) |
|---|---|
| `fitsfile *fptr` (in/out) | `&mut fitsfile` |
| `int *status`, `long *n` (out scalar) | `&mut c_int`, `&mut c_long` |
| nullable out scalar (`int *x`, may be `NULL`) | `Option<&mut c_int>` |
| `char *buf` (string buffer, out) | `&mut [c_char]` |
| `const char *s` / `char buf[N]` (in) | `&[c_char]` |
| nullable `char *` (in / out) | `Option<&[c_char]>` / `Option<&mut [c_char]>` |
| `char *arr[N]` (array of string ptrs) | `&[Option<&[c_char]>]` / `&mut [&mut [c_char]]` |
| string literal `"FOO"` | `cs!(c"FOO")`  → `&[c_char]` |
| char literal `'x'` | `bb(b'x')` → `c_char` |
| an opened file (`fitsfile *` from `ffopen`) | `Option<Box<fitsfile>>` (see below) |
| `datatype, void *value` pair (`ffpky`) | `KeywordDatatype::TXxx(&value)` |

### Open files are `Option<Box<fitsfile>>`

The canonical pattern (`ffopen_safe` / `ffclos_safe`): an opened file is an owned
`Box<fitsfile>`, and `fits_close_file` **consumes** it.

```rust
let mut fptr: Option<Box<fitsfile>> = None;
fits_open_member(gfptr, i, &mut fptr, status);          // output: &mut Option<Box<fitsfile>>
fits_add_group_member(out, fptr.as_deref_mut(), 0, st); // borrow: Option<&mut fitsfile>
if let Some(f) = fptr.take() { fits_close_file(f, st); }// move out to close
```

If a dependency still outputs a raw `*mut fitsfile` (e.g. `ffreopen_safer`), bridge it
in a small `unsafe` block with a `// SAFETY:` note:
```rust
let mut raw: *mut fitsfile = std::ptr::null_mut();
fits_reopen_file(mfptr, &mut raw, status);
*gfptr = if raw.is_null() { None } else { Some(unsafe { Box::from_raw(raw) }) };
```

## Control flow

- `do { ... } while(0);` → `loop { ...; break; }`. A C `continue` inside it →
  `break`. **Nested** do/whiles → labeled `'outer:` / `'inner:` loops; pick `break
  'outer` vs `break 'inner` to match which `while(0)` the `continue` targeted.
- `for (i = a; cond; ++i) { ... continue ... }`: if the body has `continue`, bump the
  counter at the **top** so `continue` still advances it:
  ```rust
  i = a;                 // C: for(i=a; ...; ++i)
  while cond {
      i += 1;            // ... ++i moved up so `continue` works
      /* body with continue */
  }
  ```
  But if `i` is **read after the loop** (e.g. recording a failed index), keep the
  increment at the **bottom** to preserve C's exit value.
- `switch (x) { case FOO: ... }` → `match x as u64 { FOO => ..., _ => ... }`. The
  CFITSIO option constants (`OPT_*`, `GT_ID_*`) are `u64`, so match on `x as u64`.
  C `case` `break`s are switch-breaks → just separate match arms (no `break`).

## String / buffer helpers (src/wrappers.rs)

Walk strings with **index variables**, not pointers. Use the safe wrappers instead of
libc:

`strcpy_safe(dst, src)`, `strcat_safe(dst, src)`, `strlen_safe(s)`,
`strncpy_safe`, `strcmp_safe`, `strncmp_safe`, `strchr_safe(s, c) -> Option<usize>`,
`strrchr_safe(s, c) -> Option<usize>`, `strstr_safe(s, t) -> Option<usize>`,
`fits_strcasecmp(s, t)`, `fits_strncasecmp(s, t, n)`.

Pointer arithmetic into a buffer becomes index math: `tmpStr1 = strstr(p, q)` and
later `p + i` → track `let idx = base + strstr_safe(&buf[base..], q)?;`. If a helper
you need is missing (a `strrchr_safe` was added this way), add it next to its
`strchr_safe` sibling rather than writing a local closure.

Other substitutions:

- `snprintf(buf, n, "fmt", a)` → `int_snprintf!(&mut buf, n, "fmt{}", a)`.
- `ffpmsg("literal")` → `ffpmsg_str("literal")`; `ffpmsg(buf)` → `ffpmsg_slice(&buf)`.
- `malloc`/`free` → own the memory (`Vec`, `Box`) and let RAII free it; **drop the
  `free` loop** with a note. A few APIs hand back a raw heap string tracked in the
  `ALLOCATIONS` map — copy it out, then release via that map under `unsafe`.
- Fixed column buffers (`char *ttype[6]` into `char buf[102]`) → a 2-D array
  `[[c_char; 17]; 6]`, viewed with `buf.each_mut().map(|c| c.as_mut_slice())` (for an
  output `&mut [&mut [c_char]]`) and `buf.each_ref().map(|c| c.as_slice())` (immutable).
  Avoid `Vec` for these — the C uses fixed stack arrays.

## Fixing a dependency's wrong stub signature

Stubs were sometimes generated with provisional/incorrect signatures (a `char*` typed
as `&mut c_char`; an output array typed immutably; a `fitsfile**` left raw). If the C
can't be transpiled against the current signature, **fix the stub** — it's a `todo!()`
with a known caller set:

1. `grep` its callers first (`grep -rn '\bffxyz_safe\b' src/`); confirm you can update
   them all.
2. Change the `_safe` signature, then update its `extern "C"` wrapper to bridge
   (`Option<&mut T>` ↔ raw, `Box::into_raw` for outputs, etc.) and every caller.

Examples done this way: `ffgtam_safe` (`&mut fitsfile` → `Option<&mut fitsfile>` for a
NULL member), `ffgmop_safe`/`ffgtop_safe` (raw `fitsfile**` → `&mut Option<Box<fitsfile>>`),
`ffgtdc`/`ffgtgc` (immutable arrays → `&mut [&mut [c_char]]`), `fits_relurl2url`
(by-value arrays → `&[c_char]`/`&mut [c_char]`).

## Aliasing & two `&mut` to the "same" object

Two `&mut fitsfile` can never alias in safe Rust, so a C check like
`if (infptr == outfptr)` is effectively always false; keep it (`std::ptr::eq`) for
fidelity with a note, but don't rely on it. When a C function takes the same pointer
twice (output + input), restructure so only one `&mut` is live at a time, or change
the dependency to `Option<&mut T>`.

## Platform `#if` blocks

Transpile **every** platform branch, selecting it at runtime with `cfg!(target_os =
...)` — both arms are compiled (and linted), only one runs. Follow the existing
`fits_path2url` / `fits_url2path`:

```rust
if cfg!(target_os = "windows") {
    /* the C  #if defined(MSDOS) || defined(WIN32)  branch */
} else if cfg!(target_os = "macos") {
    /* the C  #elif defined(macintosh)  branch */
} else {
    /* the C  #else  (default Unix) branch */
}
```

Map the C macros to targets: `MSDOS`/`__WIN32__`/`WIN32` → `target_os = "windows"`,
`macintosh` → `target_os = "macos"`, `#else` → the `else`. Branches with **no Rust
target** (`VMS`, and the older `WINNT` `//disk/...` variant that `WIN32` subsumes) are
not ported — say so in a comment, as `fits_path2url` does.

`ffstrtok(buff, "/", &saveptr)` token loops become a safe split that skips empty
tokens (strtok merges consecutive/leading delimiters):
`buff[..strlen_safe(&buff)].split(|&c| c == bb(b'/')).filter(|t| !t.is_empty())`. Note
the split tokens are **not** NUL-terminated, so append them by length
(`out[l..l+tok.len()].copy_from_slice(tok); out[l+tok.len()] = 0;`) rather than with
`strcat_safe` (which scans for a NUL).

When a transpiled function's output is platform-dependent, make its **test**
platform-aware too — pick the expected values with `if cfg!(target_os = ...)` (see
`test_fits_url2path`). Where genuine platform behavior is needed, the crate also uses
`std` (see `fits_get_cwd`).

## Unit tests

Add tests to the file's `#[cfg(test)] mod tests`. Convert between `&str` and
NUL-terminated `[c_char; N]` buffers:

```rust
fn to_buf(s: &str) -> [c_char; FLEN_FILENAME] {
    let mut b = [0 as c_char; FLEN_FILENAME];
    for (i, &c) in s.as_bytes().iter().enumerate() { b[i] = c as c_char; }
    b
}
fn from_buf(b: &[c_char]) -> &str {
    CStr::from_bytes_until_nul(cast_slice(b)).unwrap().to_str().unwrap()
}
```

Derive expected values by **hand-tracing the C**, not by trusting the new Rust. Cover
the tricky edges (truncation, trailing escapes, the empty string). Beware functions
that post-process: e.g. `fits_relurl2url` ends by calling `fits_clean_url`, which
collapses `//` and resolves `.`/`..`, so a `scheme://host/...` input gets mangled —
test with plain paths where the result is unambiguous.

Functions that need a live `fitsfile` or depend on a `todo!()` helper generally can't
be unit-tested yet; say so rather than writing a test that just panics.
