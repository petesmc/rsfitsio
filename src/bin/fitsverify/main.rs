#![allow(deprecated)]
// The C-derived declarations in common.rs cover the whole of fitsverify; the
// port is still incomplete, so many are not referenced yet.
#![allow(dead_code)]
// kwdtyp, data and the *_KEY variants keep their fitsverify C spellings.
#![allow(non_camel_case_types, clippy::upper_case_acronyms)]

use std::process::ExitCode;

mod common;
mod fvrf_data;
mod fvrf_file;
mod fvrf_head;
mod fvrf_key;
mod fvrf_misc;

pub fn main() -> ExitCode {
    ExitCode::from(0)
}
