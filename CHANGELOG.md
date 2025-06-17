# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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