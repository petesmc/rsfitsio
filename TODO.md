- [ ] Refactor the NullCheckType and NullValue. NullValue isn't Copy and these concepts are intrinsicly linked.
- [ ] Investigate all code with 'WARNING'
- [ ] Investigate all code with 'TODO'
- [ ] Remove use a malloc and free
- [ ] Remove use of libc unsafe functions. The raw C string functions are gone
      from every file except `eval_y.rs` (38 sites) and `eval_f.rs` (6). Those
      are one connected knot: the `*mut *mut c_char` per-row string buffers
      behind `NodeValue::Buffer`, the `sptr1`/`sptr2` aliases into them, and the
      `bit*` helpers (`bitcmp`, `bitlgte`, `bitand`, `bitor`, `bitnot`) that take
      raw pointers. Untangling needs the string-node storage reworked, not a
      call-site swap. `strtok_r`, `strpbrk`, `strspn`, `strcspn`, `strncpy` and
      `strchr` have been deleted outright as dead.
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