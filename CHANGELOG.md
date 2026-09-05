# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed
- Allocation ownership put on one ladder. `ffiter_internal` owns its per-column work
  arrays instead of `calloc`ing seventeen of them, so they are released on every path
  and a failed reservation reports `MEMORY_ALLOCATION` rather than a null pointer that
  only a later check would notice. The allocation registry moved behind
  `helpers::raw_owned`, which records the element size and alignment alongside the
  length and capacity and refuses to release a pointer it never handed out.
- No `libc` allocation is left in the library. The expression engine's node buffers
  come from the Rust allocator, and `NodeValue::Buffer` carries the `Layout` to
  release them with plus whether the node owns the block at all -- so a column node
  pointing into `varData` can no longer be freed through, and both of
  `Allocate_Ptrs`' allocations go together instead of the character block having to
  be released by hand at each call site.
- The decompression routines' `mem_realloc` callback is no longer `libc::realloc`
  where the buffer being grown is a Rust `Vec`.

### Fixed
- Compressed files were not recognised at all. Two independent defects in the same
  path: the two-byte magic numbers for GZIP, PACK, LZW and LZH were transpiled from the
  C's octal escapes as if they were decimal, so `file_is_compressed` and
  `mem_compress_open` rejected every `.gz`, `.Z` and `.bz2` file; and
  `file_is_compressed` read those two bytes with a bare `read()`, which is allowed to
  come back short, and treated a short read as "not compressed".
- Memory-file buffers allocated by `mem_createmem` (a `Vec`, i.e. the Rust allocator)
  were resized with the C `realloc` by `stdin2mem`, both compressed-open paths and the
  decompression grow callback -- a cross-allocator free. They route through
  `owned_realloc` now, and `stdin2mem` publishes the new address as it grows instead of
  leaving a stale one behind for its own error path to free.
- Two frees on a guessed layout: `fits_read_wcstab` released an `nelem`-element `f64`
  array as if it held `ndim` elements, and `FITSfile::drop` freed the six `tile*`
  pointers that `TILE_STRUCTS` owns. The latter is the double free that had
  `test_copy_within_same_file` ignored.

## [0.470.1] - 2026.08.16

### Changed
- Large reduction in raw `libc` and unsafe pointer usage: the last raw libc string functions
  are gone from the expression engine, references and constness changes now go through
  `ptr::from_ref`/`from_mut`, `cast_mut`/`cast_const`, `pointer::cast` and `&raw`, and the five
  raw-pointer cast lints are denied package-wide.
- Pedantic clippy lints cleared across the library, tests and the transpiled binaries.
- C-ABI buffer reinterpretations now assert pointer alignment, and the tiled-image scratch
  buffers are allocated with explicit alignment.
- Miri findings documented and the remaining known issues recorded.

### Fixed
- Undefined behaviour reported by Miri: aliasing in the lexer buffer stack, the
  self-referential memory driver, and pointer provenance in `grparser` and `getcol`.
- The lexer no longer defaults `yyin`/`yyout` to the C stdio externs.
- Short-read handling in the memory, file, IRAF and region readers.
- Double free when the `headstart` array fails to grow (#100).

## [0.470.0] - 2026.08.15

### Added
- `fpack` and `funpack`, transpiled from CFITSIO.
- `fitsverify`, transpiled from CFITSIO, with unit tests and coverage.
- Differential oracle test running the expression parser against a recorded corpus.

### Changed
- Synced to CFITSIO 4.7.0 (`d6d2765`).
- Expression parser rewritten around typed enums (`NodeValue`, `ValueSort`, `funcOp`, typed
  operators) instead of unions and integer tags; GTI and region loading split out of
  `New_GTI`/`New_REG`, with regions shared behind an `Rc`.
- Large reduction in `unsafe`: 68 redundant blocks removed, raw pointers dropped from several
  `_safe` functions, error stack messages moved to a fixed-size buffer.
- Clippy clean across the crate, tests and examples; VAX/VMS/Alpha conditionals removed.
- `pliocomp` 0.6, for its worst-case buffer bound.

### Fixed
- Quickselect pivot bug breaking lossy compression reproducibility against CFITSIO (#82).
- `usize` underflow when unshuffling GZIP_2 tiles (#89).
- Image and table compression defects exposed by `fpack`/`funpack`; `BITPIX` and `BLOCKED`
  restored to the img2comp keyword table.
- `printf` conformance for `%g` and the `#` flag.
- `ffdtyp` string detection, hex literal case and bare `.` lexing, constant-folded `ANGSEP`,
  four parser bugs found while testing the expression operators, and a pending syntax error
  being discarded during error recovery.
- The stdout/stdin stream driver, which now writes its memory file through Rust rather than
  the C runtime.

## [0.464.0] - 2026.07.12

### Added
- Test harness ported from CFITSIO, and a large body of unit tests including the compression
  routines.
- Implementations for the remaining `todo!()` stubs: group and hierarchy routines (`ffgt*`,
  `ffgm*`), URL helpers (`fits_get_url`, `fits_url2path`, `fits_relurl2url`, ...), `ffimem`,
  `ffextn`, `ffexist`, `ffpcln`, `ffparsecompspec` and `fits_copy_cell2image`.

### Changed
- Use `core` instead of `std` where possible.

### Fixed
- A significant number of bugs found by the new test harness.
- Compatibility problems found by integration testing against the Python `fitsio` library.
- `iraffits` handling of long Windows filenames and `/` in place of `\`.
- Incorrect `DLONG_MAX`/`DLONG_MIN` on Windows, incorrect `c_char` usage, and mutex poisoning.

## [0.463.0] - 2026.05.30

### Added
- `grparser` (template parser) implementation.
- `histo` (histogram) implementation.

### Changed
- Minimum supported Rust version is 1.96.
- `grparser` uses indices and `Box<[T]>` in place of raw pointers.

### Fixed
- Use after free bug.
- `c_char` casting issue.
- Crash caused by bit string row filters.

## [0.462.10] - 2025.09.20

### Fixed
- Fixed compile issue with smem unexpectingly being enabled on windows.

## [0.462.9] - 2025.09.20

### Fixed
- Fixed compile error

## [0.462.8] - 2025.09.20

### Added
- Comprehensive eval* tests for table where clauses and calculator expressions
- Expanded cookbook examples for binary tables
- Unit tests for eval functions and row expressions
- Python fitsio compatibility improvements

### Changed
- Rewritten error handling with new ErrorStack and ErrorMessage structs for thread safety

### Fixed
- Bzip2 feature now fully functional
- Multiple safety issues and TODO implementations across modkey.rs, putkey.rs
- Various bugs found through Python fitsio integration testing
- Overflow handling in eval functions
- Table where clause parsing and execution

## [0.462.7] - 2025-08-04

### Fixed
- Critical bug typo

## [0.462.6] - 2025-08-04

### Added
- Shared memory code and utility functionality (previously non-functional feature)
- Expanded test coverage

### Changed
- Further reduced libc dependencies from cfileio.rs and edithdu.rs modules

### Fixed
- Bugs with long keywords and null NAXES handling
- fgets usage in region.rs

## [0.462.5] - 2025-06-29

### Added
- CI support for MacOS (Intel and ARM), Windows (x86_64 and ARM), and Linux ARM64
- Iterator examples demonstrating FITS file traversal

### Changed
- ZCompress module now uses Read/Write traits for better Rust integration

### Removed
- Dependency on nightly variadic functions - now works on stable Rust

### Fixed
- Platform-specific issues for MacOS and Windows builds
- stdin/stdout/stderr linking on MacOS
- Type compatibility issues on Windows (c_long vs LONGLONG)
- Various TODOs and documentation improvements

## [0.462.4] - 2025-06-17

### Added
- Cookbook examples (`cookbook_c` and `cookbook_rust`) demonstrating FITS file operations
- Speed utility benchmark tool (`speed` binary) for performance testing
- GitHub Actions integration with dependabot for automated dependency updates

### Removed
- Removed dependency on `errno` crate

### Fixed
- Various TODO items throughout the codebase
- Minor code improvements and optimizations

## [0.462.3] - 2025-06-06

### Added
- Split aliases into separate `c_api` and `rust_api` modules for better organization
- More safe wrapper functions for improved memory safety

### Changed
- Improved visibility controls for internal functions
- Better separation between C-compatible and Rust-native APIs

### Fixed
- Fixed visibility issues with certain functions
- Various formatting and clippy warnings

## [0.462.2] - 2025-06-01

### Added
- Full Miri support - all tests now pass under Miri memory safety checker
- Codecov integration for test coverage reporting

### Changed
- Updated from Rust edition 2021 to 2024
- Minimum supported Rust version set to 1.87
- Removed dependency on nightly features (except variadics)

### Removed
- Removed various libc function dependencies in favor of safe Rust alternatives

### Fixed
- Many Miri-detected memory safety issues
- Various unsafe code patterns replaced with safe alternatives
- Improved error handling throughout the codebase

## [0.462.1] - 2024-XX-XX

### Changed
- Initial line-by-line Rust translation of cfitsio v4.6.2
- Maintains full C API compatibility

### Added
- Dual API design: C-compatible functions and safer Rust APIs
- Pluggable I/O driver system
- Compression support (rice, plio, hcompress)
- Thread-safe error handling

### Known Issues
- Features `bzip2` and `shared_mem` not fully functional
- Limited macOS testing
- Some platform-specific code uses outdated cfg values (winnt, vms)

## Notes

This project is a direct Rust translation of the cfitsio C library (v4.6.2), maintaining
compatibility while providing Rust's memory safety guarantees. The version numbering
follows cfitsio's version with an additional patch number for Rust-specific changes.

Synchronization with upstream cfitsio is tracked in `SYNCED_COMMIT.md`.