//! Shared helpers for the integration tests.
//!
//! These are duplicated from the crate-internal `helpers::testhelpers` module
//! so that the integration tests do not force `tempfile` to be a normal
//! (non-dev) dependency of the library.  Not every test uses every helper, so
//! dead-code warnings are suppressed.
#![allow(dead_code)]

use tempfile::Builder;

/// Function to allow access to a temporary file
pub fn with_temp_file<F>(callback: F)
where
    F: for<'a> Fn(&'a str),
{
    let tdir = Builder::new().prefix("rsfitsio-").tempdir().unwrap();
    let tdir_path = tdir.path();
    let filename = tdir_path.join("test.fits");

    let filename_str = filename.to_str().expect("cannot create string filename");
    callback(filename_str);
}

/// Helper function for float comparisons
pub fn floats_close_f32(a: f32, b: f32) -> bool {
    (a - b).abs() < f32::EPSILON
}

pub fn floats_close_f64(a: f64, b: f64) -> bool {
    (a - b).abs() < f64::EPSILON
}
