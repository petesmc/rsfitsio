#![allow(deprecated)]

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
