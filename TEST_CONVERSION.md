# Converting CFITSIO C unit tests into Rust unit tests

This crate ports CFITSIO. Alongside porting functions (see
[TRANSPILING.md](TRANSPILING.md)), we port the C test programs (`test_*.c`,
exercised with the macros in `test_macros.h`) into Rust `#[cfg(test)]` modules so
the safe functions get real coverage. This guide captures the conventions, using
the `test_group.c` → `src/group.rs` conversion as the worked example.

## The workflow

A request usually looks like: *"port the unit tests in `test_foo.c` into
`foo.rs`."* Each C `test_xxx(void)` function becomes one Rust `#[test] fn
test_xxx()`. To convert one:

1. **Read the C test** and the `_safe` signatures of every function it calls
   (see "Finding the functions"). Confirm none of them are still `todo!()` — if a
   dependency panics, the test can't pass yet; say so rather than writing it.
2. **Decide the file backing** (see "Creating a FITS file"): an in-memory
   `mem://` file if the test only creates/operates/closes one handle, or a real
   temp file if it re-opens the file by name (extended `[EXT]` syntax).
3. **Translate the body** call-by-call, mapping the `call_0N`/`fail_*` macros and
   the pointer/output arguments (see below).
4. **Derive the expected values by hand-tracing the C**, not by trusting the new
   Rust. Keep the C's `fail_if` conditions as `assert!`/`assert_eq!`.
5. **Run** `cargo test --lib <module>::tests`, then the full `cargo test` before
   declaring done.

Add the tests to the file's existing `#[cfg(test)] mod tests { ... }` block
(create one at the end of the file if absent, with `use super::*;`).

## The C test macros

`test_macros.h` wraps each library call so it passes `&status` last and aborts on
failure. Translate them as:

| C macro | Meaning | Rust |
|---|---|---|
| `call_0N(ff, a, …)` | `ff(a, …, &status)`; abort if return != 0 **or** `status != 0` | call `ff(a, …, &mut status)`; then `assert_eq!(status, 0, "ff failed")` |
| `fail_if(x)` | abort if `x` is true | `assert!(!x)` — usually rephrase positively, e.g. `fail_if(n != 1)` → `assert_eq!(n, 1)` |
| `fail_st(x)` | abort if `x` **or** `status != 0` | the assertion plus an `assert_eq!(status, 0)` |

You don't need to assert `status == 0` after *every* call — assert it at the
points the C macro would have caught a failure (after setup blocks and before
reading a result), and always check the specific values the test cares about.
`fail_if(strcmp(extname, "GROUPING") != 0)` → `assert_eq!(read_str_key(f,
cs!(c"EXTNAME"), &mut status), "GROUPING")`.

## Finding the functions

The C calls public names (`ffgtcr`, `ffphps`, …). In a `mod tests` with `use
super::*`, the module's `use crate::aliases::rust_api::*` is in scope, so prefer
the **`fits_*` rust_api aliases** — they resolve to the `_safe` impls and read
like the C:

```bash
grep -nE 'as fits_create_group\b' src/aliases.rs   # ffgtcr_safe -> fits_create_group
grep -nE 'pub fn ffgtcr_safe' src/group.rs         # then read its real signature
```

A few file-level functions have no `fits_*` alias (e.g. `ffopen`); import the
`_safe` form directly:

```rust
use crate::cfileio::{ffclos_safe, ffinit_safe, ffopen_safe};
```

Read the `_safe` signature before writing the call — argument shapes differ from
the C (see next).

## Argument translation

The C passes raw pointers and `&status`; the `_safe` functions take
references/slices/`Option`s. Common mappings (mirror
[TRANSPILING.md](TRANSPILING.md)):

| C call site | Rust |
|---|---|
| `&status` (always last) | `&mut status` |
| out scalar `&hdutype` | `&mut hdutype` |
| nullable out scalar `…, NULL, …` | `None` (param is `Option<&mut _>`) |
| `NULL` axes array (`ffphps(f, BITPIX, 0, NULL)`) | `&[]` |
| string literal `"GROUPING"` | `cs!(c"GROUPING")` |
| out string buffer `char buf[FLEN_VALUE]` | a `[c_char; FLEN_VALUE]`, passed `&mut buf` |
| option constant `GT_ID_ALL_URI`, `OPT_RM_GPT` | `GT_ID_ALL_URI as c_int` — these consts are `u64`, the params are `c_int` |

The `GT_ID_*` / `OPT_*` casts are the easiest thing to forget: the constants in
`fitsio.rs` are `u64`, but the grouptype/rmopt/cmopt/etc. parameters are `c_int`,
so write `... as c_int`.

## Creating a FITS file

The shared helpers live in `src/helpers/testhelpers.rs` (a normal, always-compiled
`pub` module — `with_temp_file` is already used by integration tests). Import what
you need:

```rust
use crate::helpers::testhelpers::{
    from_buf, new_mem_file, path_with_ext, read_str_key, to_buf, with_temp_file,
};
```

**In-memory (`mem://`) — the default.** Use when the test only ever touches the
handle(s) it created, never re-opening by filename. `new_mem_file` returns a file
with an empty `BYTE_IMG` primary HDU (the C's `ffinit` + `ffphps`):

```rust
let mut status = 0;
let mut fptr = new_mem_file(&mut status);
let f = fptr.as_deref_mut().unwrap();
fits_create_group(f, cs!(c"TestGroup"), GT_ID_ALL_URI as c_int, &mut status);
assert_eq!(status, 0, "ffgtcr failed");
// ... assertions ...
ffclos_safe(fptr.take().unwrap(), &mut status);
```

Two independent `mem://` files (e.g. a copy test) are just two `new_mem_file`
calls — copy handle-to-handle, no reopening needed.

**Real temp file — when the C re-opens by name.** Tests that use CFITSIO's
extended-filename syntax (`test_path "[GROUPING]"`, `"[IMAGE1]"`,
`"[GROUPING,1]"`) need a path on disk. Use `with_temp_file`, build the name with
`to_buf` / `path_with_ext`, and create with `ffinit_safe` (the temp path doesn't
exist yet, so no clobber is needed — drop the C's leading `!`):

```rust
with_temp_file(|filename| {
    let mut status = 0;
    let name = to_buf(filename);
    let mut fptr: Option<Box<fitsfile>> = None;
    ffinit_safe(&mut fptr, &name, &mut status);
    fits_write_imghdr(fptr.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);
    // ... build HDUs, then re-open a second handle: ...
    let mut gfptr: Option<Box<fitsfile>> = None;
    let gname = path_with_ext(filename, "[GROUPING]");
    ffopen_safe(&mut gfptr, &gname, READWRITE, &mut status);
    assert_eq!(status, 0, "ffopen [GROUPING] failed");
    // ...
    ffclos_safe(gfptr.take().unwrap(), &mut status);
    ffclos_safe(fptr.take().unwrap(), &mut status);
});
```

`with_temp_file` makes a fresh temp dir and cleans it up on drop, so there's no C
`remove(test_path)` to port.

## Opened files are `Option<Box<fitsfile>>`

Every file handle (created or opened) is an owned `Option<Box<fitsfile>>`, exactly
as in the ported code:

- borrow it for a call: `fptr.as_deref_mut().unwrap()` (or keep `let f =
  fptr.as_deref_mut().unwrap();` when used repeatedly).
- pass it as an `Option<&mut fitsfile>` member argument (`fits_add_group_member`):
  `Some(f)` or `mfptr.as_deref_mut()`.
- output of an opener (`ffopen_safe`, `fits_open_member`) is `&mut fptr`.
- **close consumes the box**: `ffclos_safe(fptr.take().unwrap(), &mut status)`.

Two handles on the same physical file (a grouping table plus a member) is fine —
they're distinct `Box`es, so they never alias.

## Strings, buffers, and reading keywords

- `cs!(c"LITERAL")` for an input C string (`&[c_char]`); `bb(b'x')` for a `char`.
- `to_buf(s) -> [c_char; FLEN_FILENAME]` copies a `&str` into a NUL-terminated
  buffer (filenames, keyword names).
- `from_buf(&buf) -> &str` reads a NUL-terminated buffer back.
- `read_str_key(f, cs!(c"EXTNAME"), &mut status) -> String` wraps the common
  "read a string keyword and compare" pattern; values come back already
  unquoted/trimmed, so `assert_eq!(read_str_key(...), "GROUPING")` works directly.

If a helper you need is missing from `testhelpers.rs`, add it there (next to the
existing ones) rather than redefining it per module.

## What not to port

- **Tests the C itself disabled.** `test_group.c` removed
  `test_ffgmop_no_members`, `test_ffgtam_null_member`, etc. with notes that the
  library "crashes instead of returning an error." Don't resurrect these — honor
  the C's note (port it as a comment).
- **Tests whose dependencies are still `todo!()`.** They'll panic at runtime; skip
  with a one-line note until the dependency is transpiled.
- **The core compression codecs are EXTERNAL crates, not part of this repo.**
  PLIO, Rice and Hcompress have been externalised into the published crates
  `pliocomp`, `ricecomp` and `hcompress` (by cruzzil; see `Cargo.toml`).
  `imcompress.rs` calls them directly (`use pliocomp::…`, `hcompress::…`,
  `ricecomp::…`); GZIP/GZIP_2 go through `libz-rs-sys`; BZIP2 through
  `libbz2-rs-sys` (the `bzip2` default feature). The orphaned `src/pliocomp.rs`
  and `src/ricecomp.rs` thin wrappers are NOT wired into `lib.rs`. For golden
  comparison the external crate sources can be cloned into
  `~/code/external/{pliocomp,ricecomp,hcompress}` from cruzzil's GitHub. The
  standalone codec tests (`test_pliocomp.c`, `test_ricecomp.c`,
  `test_hcompress.c`) therefore exercise external crates and are out of scope
  here; bugs in those belong upstream.
- **`test_quantize.c` and `test_imcompress.c` ARE ported** (quantization lives in
  `src/quantize.rs`, tiled-image (de)compression in `src/imcompress.rs`, both
  part of this repo). Known imcompress gaps surfaced by the tests: the
  `SUBTRACTIVE_DITHER_2` zero-value preservation is commented out in the
  byte/short unquantize paths, and the BZIP2 compress branch silently no-ops
  when the `bzip2` feature is off.
- The C `main()` driver and its `remove()` cleanup — `#[test]` discovery and
  `with_temp_file`'s drop replace both.

## Deriving expected values

Take the numbers straight from the C `fail_if` conditions (`nmembers != 1`,
`hdunum != 2`, `numhdus != 1`, `firstfailed != 0`) and phrase them positively as
`assert_eq!`. Initialise out-params to a sentinel the call must overwrite (e.g.
`let mut nmembers: c_long = -1;`) so a no-op call is caught. Add a short message
to each assertion naming the C function, so a regression points at the call.

## Checklist

- [ ] One `#[test] fn test_xxx` per C `test_xxx`, same name.
- [ ] `mem://` unless the test re-opens by filename; then `with_temp_file`.
- [ ] `&status` → `&mut status`; `NULL` out-scalars → `None`; `GT_ID_*`/`OPT_*` →
      `... as c_int`.
- [ ] Every opened/created handle closed with `ffclos_safe(fptr.take().unwrap(),
      …)`.
- [ ] Expected values hand-derived from the C, asserted with messages.
- [ ] `cargo test --lib <module>::tests` green, then full `cargo test` green.
