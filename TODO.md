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
- [X] Miri currently fails, fix.
- [ ] Fix dodgy safety code in ffedit_columns. WARNING / SAFETY / TODO
- [X] Mark all extern functions as deprecated so that we can detect usage
- [X] Feature "bzip2" doesn't work
- [X] Feature "shared_mem" doesn't work
- [ ] fits_parser_yytokentype refactor to enum