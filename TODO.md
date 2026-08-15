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
- [ ] Fix the Stacked Borrows violation in the lexer's buffer stack. Miri
      rejects `eval_l.rs:1522` vs `:1703` (`yy_buffer_stack`) on any test that
      invokes the parser, so Miri cannot currently validate the eval engine.
      Pre-existing and verified against a clean worktree at fed4b4f. A second,
      separate one is in the `test_ffcalc_*` tests themselves, which alias
      `&mut *fp_self` twice to pass one file as both input and output.
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
- [X] Fix broken testprog.out comparison
- [ ] Miri still reports Undefined Behavior. Run with
      `MIRIFLAGS="-Zmiri-disable-isolation" cargo +nightly miri nextest run
      --no-fail-fast -j 8` (cap the threads: the default fans out to one miri
      process per core at ~1.1GB each and thrashes swap). All of the findings
      below are Stacked Borrows violations, i.e. real aliasing UB that LLVM's
      `noalias` is entitled to miscompile, not miri pedantry. A full run takes
      ~3.8h and reports 2287 passed / 485 failed: 435 UB, 42 unsupported
      operations, 8 alignment panics. Counts below are failing tests per site.
      Root causes, in order of blast radius:
      * eval_l.rs:1571 (258) — see the lexer entry below; by far the biggest.
      * drvrmem.rs:259 (53), getcolui.rs:1380 (38), getcol.rs:2778 (8),
        grparser.rs:1046 (5) — each already carries a TODO/SAFETY-TODO.
      * cfileio.rs:1880 (20) — the FPTR_TABLE entry below.
      * The `&mut`/`&mut` aliasing sites are 1 test each, but see the first
        entry: they are a signature problem, not a local one.
      * The ported `_safe` signatures take `infptr: &mut fitsfile` and
        `outfptr: &mut fitsfile`, but the C API explicitly allows
        `infptr == outfptr` and callers rely on it. Two `&mut` to one object
        can never be legal, so every site that reproduces the C aliasing is UB:
        `ffsrow_safe(&mut *f, &mut *f, ..)` in cfileio.rs:4893 (library code,
        with a SAFETY comment acknowledging it), plus the same trick in the
        eval_f.rs and editcol.rs tests. These signatures need to take a single
        `&mut` plus a flag, or raw pointers, before the aliasing can go away.
      * `fitsfile.Fptr` is a `Box<FITSfile>` but several FITSfiles are shared
        between handles by `Box::from_raw`, so two live Boxes own one
        allocation (cfileio.rs:1854, already marked "TODO this is very
        unsafe!"). Confirmed by miri, invalidated at getkey.rs:242.
      * `FPTR_TABLE` (cfileio.rs:1880) stores a raw pointer derived from a
        `&mut FITSfile` borrowed out of that Box. Any later write through the
        Box (e.g. fitscore.rs:5757) invalidates it, and the next
        `FPTR_TABLE[ii].as_mut()` in fits_already_open (cfileio.rs:1997) is UB.
        Reproducer: open one file, then open a second.
        Fixing this soundly means owning FITSfile through a raw pointer
        (NonNull + Deref/DerefMut wrapper, which is layout- and C-ABI-neutral)
        rather than Box, so all accesses share one provenance root.
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
- [ ] Fix dodgy safety code in ffedit_columns. WARNING / SAFETY / TODO
- [X] Mark all extern functions as deprecated so that we can detect usage
- [X] Feature "bzip2" doesn't work
- [X] Feature "shared_mem" doesn't work
- [ ] fits_parser_yytokentype refactor to enum