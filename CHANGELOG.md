# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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