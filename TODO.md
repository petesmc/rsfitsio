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
- [X] Fix the Info.parseData aliasing in the eval engine. `Info.parseData` is
      a raw borrow (`&raw mut lParse`) rather than a `&mut`, and the raw pointer
      to `lParse.colData` is taken before it rather than after, so storing one
      no longer invalidates the other and the work function's `as_mut()` is
      sound.
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
- [ ] Report the zlib-rs soundness bug upstream, then revert Cargo.toml to
      `flate2 = { ..., features = ["zlib-rs"] }` so flate2 shares the
      libz-rs-sys already in the tree instead of pulling in miniz_oxide.
      `DeflateStream` holds `state: &'a mut State<'a>` pointing into the block
      that `deflate::end` then frees while that reference is still a protected
      function argument, so `GzEncoder::finish` is UB: miri reports
      "deallocating while item is strongly protected". A 12-line GzEncoder
      write-and-finish reproduces it with no rsfitsio code involved. Present in
      0.6.6, 0.6.7 and upstream main. Only zlib-rs's `stable` API is affected --
      libz-rs-sys's C `deflateEnd` builds its `&mut DeflateStream` differently
      and is clean -- so the library's own GZIP paths were never involved, only
      the two fpack/funpack gzip tests.
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
- [X] Miri is clean. The confirming end-to-end run reported 2780 passed, 0
      failed, 59 skipped in 5h12m: no UB, no leaks, no unsupported operations.
      Re-run it after any change to the aliasing-sensitive code below. Use
      `MIRIFLAGS="-Zmiri-disable-isolation" cargo +nightly miri nextest run
      --no-fail-fast -j 8` (cap the threads: the default fans out to one miri
      process per core at ~1.1GB each and thrashes swap). A full run takes
      ~5h. Do not run any other cargo command while one is in flight -- it
      corrupts the run.
      The run before the fixes reported 485 failures. What it found, and where
      each part went:
      * 435 UB, all Stacked Borrows violations -- real aliasing UB that LLVM's
        `noalias` is entitled to miscompile, not miri pedantry. The dominant
        shape, five times over, was a struct holding a pointer to something
        also reached through a `&mut`: FPTR_TABLE, the lexer's yy_buffer_stack,
        yyextra_r, Info.parseData, histData.iterCols and def_fptr. The rest
        were slices rebuilt with slice::from_raw_parts(_mut) over a pointer
        still borrowed elsewhere (eval_l, drvrmem, getcolui, getcol, grparser),
        `infptr == outfptr` aliasing (fixed with `_inplace` twins), a TSTRING
        path that assumed every caller's char* was FLEN_VALUE bytes, an
        ffgcvn_safe length multiplied twice, and ffomem declaring `void **` as
        `*const *const`.
      * 19 leaks, in the eval engine's row buffers: the BITSTR column buffers,
        const nodes' buffers, and Do_Deref's error paths releasing only one of
        an Allocate_Ptrs node's two allocations.
      * 8 alignment panics in bytemuck cast_slice_mut, from reading TSHORT into
        a `Vec<u8>` tiled-image buffer that is only 1-byte aligned. Native
        malloc happens to over-align, so this passed outside miri; it was a
        real latent panic. Fixed with helpers/aligned.rs.
      * 42 unsupported operations, all miri's own gaps rather than defects:
        `socketpair`/`copy_file_range` in tests/test_fpack_cli.rs, which spawns
        the real executables, and one libc `fopen` in a test of ffwrhdu, which
        takes a C `FILE *`. Both are now skipped under miri so that a run
        reports only things worth looking at.
      Still open, and not a miri finding as such: FPTR_TABLE is a `static mut`.
      Wrapping it in a Mutex over a Send newtype would confine that to one
      `unsafe impl`.
- [X] Fix dodgy safety code in ffedit_columns. The `same_ftpr` launder is gone;
      it calls `ffcalc_inplace_safe` instead.
- [X] Mark all extern functions as deprecated so that we can detect usage
- [X] Feature "bzip2" doesn't work
- [X] Feature "shared_mem" doesn't work
- [ ] fits_parser_yytokentype refactor to enum