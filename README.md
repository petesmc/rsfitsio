[![Crates.io](https://img.shields.io/crates/v/rsfitsio.svg)](https://crates.io/crates/rsfitsio)
[![Actions Status](https://github.com/petesmc/rsfitsio/workflows/CI/badge.svg)](https://github.com/petesmc/rsfitsio/actions)
[![Documentation](https://docs.rs/rsfitsio/badge.svg)](https://docs.rs/rsfitsio/)
[![codecov](https://codecov.io/gh/petesmc/rsfitsio/branch/master/graph/badge.svg?token=0DNCU7VRH2)](https://codecov.io/gh/petesmc/rsfitsio)
[![Dependency status](https://deps.rs/repo/github/petesmc/rsfitsio/status.svg)](https://deps.rs/repo/github/petesmc/rsfitsio)

Rust rewrite of cfitsio.

This is a line for line translation of cfitsio, attempting to keep it as compatible as possible and synced with the main repo. The repo exposes a compatible C API.

There are many safety issues to still deal with and we do not support all the targets that the original C code does.

Has not been tested on MacOS.

## C API Compatibility

This library exposes the complete cfitsio C API with one exception:

- `CUnit2FITS()` - This function is not exposed as it is primarily used for Fortran wrapper support, which is not supported in this Rust implementation.

All other 662 C API functions from cfitsio are available.

[cfitsio Repo](https://github.com/HEASARC/cfitsio)

[cfitsio Homepage](https://heasarc.gsfc.nasa.gov/fitsio/)

[Fits File Standard](https://fits.gsfc.nasa.gov/fits_standard.html)