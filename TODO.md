- [ ] Refactor the NullCheckType and NullValue. NullValue isn't Copy and these concepts are intrinsicly linked.
- [ ] Investigate all code with 'WARNING'
- [ ] Investigate all code with 'TODO'
- [ ] Remove use a malloc and free
- [X] Remove use of libc unsafe functions. Every raw C string function is gone:
      the whole `strcpy`/`strncpy`/`strcmp`/`strncmp`/`strlen`/`strnlen`/
      `strchr`/`strstr`/`strcat`/`strncat`/`strtok_r`/`strpbrk`/`strspn`/
      `strcspn` family has been deleted from `wrappers.rs`, along with the
      `cbitset` dependency it pulled in. The expression engine reaches its
      per-row string buffers through `lval::str_row`/`str_row_mut` instead of a
      `char **` row-pointer array, and `bitand`/`bitor`/`bitnot`/`bitcmp`/
      `bitlgte`/`cstrmid` are now safe functions taking slices.
      Remaining follow-up, none of it string-related: stage 3 of
      `notes/EVAL_STRING_STORAGE.md` (own the row buffers in `ParseData` so the
      two `str_row` accessors stop being `unsafe`), then stage 4 (the numeric
      buffers and `value.undef`, ~660 sites, which would remove the last
      `malloc`/`free` and raw pointer from `NodeValue`).
- [X] Fix the Stacked Borrows violation in the lexer's buffer stack
      (`yy_buffer_stack`). yyrestart held a reference to the stack slot across
      yy_create_buffer, which re-enters the scanner; yy_init_buffer and
      yy_flush_buffer took `&mut yy_buffer_state` alongside `yyscanner`, which
      also reaches that buffer; and yy_init_buffer flushed before its own
      writes, so yy_load_buffer_state invalidated them. `b` is a raw pointer
      now and the comparisons go through yy_current_buffer_ptr.
- [X] Fix the yyextra_r aliasing in the eval engine. ParseData is passed to
      yylex as an argument and the yyextra_r field is gone, so the scanner no
      longer refers to it and cannot alias the `&mut ParseData` that ffiprs
      and yyparse hold.
- [ ] Fix the Info.parseData aliasing in the eval engine. `eval_f.rs` stores
      `Info.parseData = &mut lParse` and then invalidates it two lines later
      with `&mut lParse.colData[..]`, so the iterator work function's later
      `as_mut()` on it is UB. This is the front of the queue for the eval
      engine now that yyextra_r is gone. Same shape as the FPTR_TABLE and
      infptr/outfptr entries below: a struct holding a pointer to something
      that is also reached through a `&mut`.
- [X] The `test_ffcalc_*` tests alias `&mut *fp_self` twice to pass one file
      as both input and output. Done: every `fp_self` is gone; the tests call
      the `_inplace` forms instead. See the infptr/outfptr entry below.
- [ ] Clean up all warnings
- [ ] Remove clippy allow(unused_assignments)
- [ ] Remove clippy allow(unused_variables)
- [ ] Test IRAF code. Ported but untested.
- [ ] Fix all todo!()s
- [ ] Implement Utility programs (done: fitsverify, fpack, funpack, imcopy,
      fitscopy, speed, smem; remaining: cookbook, iter_*, testprog helpers)
- [ ] Compressed tables: a variable-length *string* column whose payload
      exceeds cm_buffer (sized for descriptors only, as in the C) is refused
      with DATA_DECOMPRESSION_ERR rather than uncompressed. This is
      heasarc/cfitsio#134, where the C overruns the buffer instead; fixing it
      properly needs the buffer sizing reworked. Pinned by
      test_oversized_vla_column_is_refused_not_overrun
- [ ] Fix USHORT decompression: a BZERO=32768 image compressed by the C fpack
      fails to uncompress with NUM_OVERFLOW under every algorithm. Pinned by
      the ignored test_c_packed_ushort_is_unreadable
- [ ] Report the CFITSIO defects the fpack/funpack port turned up. Two are
      already filed upstream (heasarc/cfitsio#134 and #136); the rest are
      marked NOTE (upstream bug N) in src/bin/fpack/
- [ ] Restructure modules, ::api ??
- [X] Every extern function should be a wrapper around a safe interface
- [ ] Keep shrinking `unsafe` outside the FFI wrappers. Measure with
      `notes/unsafe_metric.py` (non-FFI blocks and the lines inside them);
      the crate-wide `#![allow(deprecated)]` is gone from lib.rs so the
      deprecation warnings once again show which internal callers are still on
      the C ABI (`cargo check --all-targets 2>&1 | grep -c deprecated`).
      Done so far: ffopen_safe (998-line block -> none; it needed `unsafe` for
      exactly one line, a `CStr::from_ptr` on a `[*const c_char; 3]` that is now
      `[&CStr; 3]`), fits_already_open, ffsrow_safe, fffrwc_safe,
      fits_execute_template (now has a `_safe` form), the zcompress entry points
      with already-safe signatures (which deleted 15 blocks in imcompress.rs and
      one in drvrfile.rs), the `.filename` sites that bypassed
      `get_filename_as_cstr`, and `strto_float_impl`.
      Still open, in rough value order:
      * `uncompress2mem` / `uncompress2mem_from_mem` still take
        `buffptr: *mut *mut u8` plus a `mem_realloc` callback whose body is a
        `panic!("not implemented")`. Converting them to `&mut [u8]`/`&mut Vec<u8>`
        makes them safe and unblocks fits_uncompress_table_safe, whose 1012-line
        block exists *only* because it calls uncompress2mem_from_mem -- the body
        has no unsafe operation of its own.
      * ffiter_safe (putcol.rs, 1312 lines) and fits_parser_workfn_safe
        (eval_f.rs, 661 lines) have genuine pointer work; narrow, do not remove.
      * The putkey.rs ffpkn*_safe family take `&[*const c_char]` and pay ~20
        `CStr::from_ptr` blocks for it; `&[&[c_char]]` (as ffhdr2str_safe already
        uses) moves that to the extern wrapper.
      * src/bin/speed/main.rs is the only binary still on the C ABI (48 sites).
      * fits_pixel_filter_safer keeps a function-wide block on purpose:
        PixelFilter is a C struct of raw pointers, so narrowing it means changing
        that struct. Documented in place.
- [X] Fix broken testprog.out comparison
- [ ] Miri still reports Undefined Behavior. Run with
      `MIRIFLAGS="-Zmiri-disable-isolation" cargo +nightly miri nextest run
      --no-fail-fast -j 8` (cap the threads: the default fans out to one miri
      process per core at ~1.1GB each and thrashes swap). All of the findings
      below are Stacked Borrows violations, i.e. real aliasing UB that LLVM's
      `noalias` is entitled to miscompile, not miri pedantry. A full run takes
      ~3.8h and reports 2287 passed / 485 failed: 435 UB, 42 unsupported
      operations, 8 alignment panics. Counts below are failing tests per site.
      Root causes, in order of blast radius. The line numbers are from that
      run; the lexer, drvrmem, getcolui, getcol and grparser ones are fixed:
      * eval_l.rs:1571 (258) — fixed; those tests now stop on the yyextra_r
        aliasing above instead.
      * drvrmem.rs:259 (53), getcolui.rs:1380 (38), getcol.rs:2778 (8),
        grparser.rs:1046 (5) — fixed. Two uncovered a further problem behind
        them: drvrmem frees Rust-allocated buffers through libc realloc, and
        the getcol TSTRING path assumes every caller-supplied char* is at
        least FLEN_VALUE bytes.
      * cfileio.rs:1880 (20) — the FPTR_TABLE entry below.
      * The `&mut`/`&mut` aliasing sites are 1 test each, but see the first
        entry: they are a signature problem, not a local one.
      * FIXED for ffsrow/ffcalc/ffcalc_rng/ffcpcl/ffcpky. The C API allows
        `infptr == outfptr` but two `&mut` to one object can never be legal, so
        each of those now has an `_inplace` twin taking a single `&mut`, and the
        `extern "C"` wrapper compares the two raw pointers and dispatches --
        matching what ffcopy/ffcpfl/ffcphd/ffcpdt already did. The bodies are
        deliberately duplicated rather than parameterised; each pair carries a
        NOTE pointing at its twin. The `&mut *f` launders in cfileio.rs
        (ffselect_table, ffedit_columns) and the eval_f.rs/editcol.rs tests are
        all gone. Note the wrapper test is *handle* identity: two distinct
        handles sharing one FITSfile (fits_reopen_file) is a supported two-file
        case, and the FITSfile-level `ptr::eq(a.Fptr.as_ptr(), ..)` checks in
        ffcpcl_safe/ffcpdt_safe are a different test and stay.
      * FIXED: `fitsfile.Fptr` is an `FptrRef` (repr(transparent) NonNull with
        Deref/DerefMut), not a `Box`, so sharing a FITSfile between handles no
        longer means two live Boxes owning one allocation.
      * `FPTR_TABLE` now stores `Fptr.as_ptr()`, the same pointer the handles
        deref through, and fits_already_open takes a shared `&FITSfile` from it
        rather than minting a `&mut` that aliased every live handle. What is
        left is that `FPTR_TABLE` is still a `static mut`; wrapping it in a
        Mutex over a Send newtype would confine that to one `unsafe impl`.
      * eval_l.rs:1571 — the lexer's yy_buffer_stack hands out `&mut` into a
        raw buffer that is then invalidated at eval_l.rs:1752.
      * drvrmem.rs:259 — `m[ii].memaddrptr = &mut m[ii].memaddr` is a
        self-referential struct; any later access through the slice kills it.
      * cfileio.rs:8907 — returns a raw pointer into the caller's buffer which
        is used after a later write invalidates the tag.
      * getcolui.rs:1380, getcol.rs:2778, grparser.rs:1046 — each rebuilds a
        slice with slice::from_raw_parts(_mut) from a pointer that is still
        borrowed elsewhere. grparser's also makes a *mutable* slice from a
        CStr's const pointer, which miri reports as a write through a
        SharedReadOnly tag.
      Not UB, but surfaced by the same run:
      * 8 tests panic in bytemuck `cast_slice_mut` with
        TargetAlignmentGreaterAndInputNotAligned, at getcol.rs:2299 reading
        TSHORT into the tiled-image `cbuf` (imcompress.rs:7720), which is a
        `Vec<u8>` and so only 1-byte aligned. Native malloc happens to return
        over-aligned blocks, so this passes outside miri; it is a real latent
        panic that fires whenever the buffer lands on an odd address. Fixing
        it means giving cbuf an alignment suitable for the widest type read
        into it rather than relying on the allocator.
      * 42 tests fail on miri's own gaps, not on defects: `socketpair` and
        `copy_file_range` (the fpack/funpack/fitsverify tests spawn
        subprocesses). Exactly one is ours: a remaining C `fopen` call.
- [X] Fix dodgy safety code in ffedit_columns. The `same_ftpr` launder is gone;
      it calls `ffcalc_inplace_safe` instead.
- [X] Mark all extern functions as deprecated so that we can detect usage
- [X] Feature "bzip2" doesn't work
- [X] Feature "shared_mem" doesn't work
- [ ] fits_parser_yytokentype refactor to enum