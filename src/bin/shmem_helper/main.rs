//! Cross-process test helper for the shared-memory driver.
//!
//! SKETCH (Phase 0 / cross-process test in HARDENING_PLAN_drvrsmem.md). This is
//! a single-operation-per-invocation helper so a parent test process can drive
//! genuinely separate processes that share one SysV segment. It honors the
//! `SHMEM_LIB_KEYBASE` / `SHMEM_LIB_MAXSEG` env vars (read by the driver at
//! init) so the parent can put each run in a private key namespace.
//!
//! Usage:
//!   shmem_helper create <hN> <seed>   create shmem FITS image, fill from seed
//!   shmem_helper read   <hN> <seed>   open + verify pixels match seed pattern
//!   shmem_helper delete <hN>          delete the segment (teardown)
//!
//! Exit code 0 = success/match; non-zero = error or mismatch.
//!
//! The segment is created PERSIST (the driver's smem_create sets
//! SHARED_RESIZE | SHARED_PERSIST), so it survives this process exiting and is
//! visible to a later `read` invocation in another process.

#[cfg(not(all(feature = "shared_mem", not(windows))))]
fn main() {
    eprintln!("shmem_helper requires --features shared_mem on a non-Windows target");
    std::process::exit(2);
}

#[cfg(all(feature = "shared_mem", not(windows)))]
fn main() -> std::process::ExitCode {
    use bytemuck::{cast_slice, cast_slice_mut};
    use rsfitsio::aliases::rust_api::{
        fits_close_file, fits_create_file, fits_delete_file, fits_open_file, fits_read_img,
        fits_write_img, fits_write_imghdr,
    };
    use rsfitsio::c_types::{c_char, c_int, c_long};
    use rsfitsio::fitsio::{READONLY, READWRITE, SHORT_IMG, TSHORT, fitsfile};
    use std::process::ExitCode;

    const NPIX: usize = 100; // 10x10 SHORT_IMG, mirrors the in-crate tests

    /// NUL-terminated Vec<c_char> from &str (same helper the unit tests use).
    fn cc(s: &str) -> Vec<c_char> {
        let mut v: Vec<c_char> = s.bytes().map(|b| b as c_char).collect();
        v.push(0);
        v
    }

    /// Deterministic pattern so writer and reader agree without shared state.
    fn pattern(seed: i32) -> [i16; NPIX] {
        let mut d = [0i16; NPIX];
        for (i, x) in d.iter_mut().enumerate() {
            *x = (seed.wrapping_add(i as i32) & 0x7fff) as i16;
        }
        d
    }

    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        eprintln!("usage: shmem_helper <create|read|delete> <hN> [seed]");
        return ExitCode::from(2);
    }
    let op = args[1].as_str();
    let url = format!("shmem://{}", args[2]); // e.g. shmem://h5
    let seed: i32 = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(0);

    let mut status: c_int = 0;
    let mut f: Option<Box<fitsfile>> = None;

    match op {
        "create" => {
            let naxes: [c_long; 2] = [10, 10];
            let data = pattern(seed);
            fits_create_file(&mut f, &cc(&url), &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), SHORT_IMG, 2, &naxes, &mut status);
            fits_write_img(
                f.as_deref_mut().unwrap(),
                TSHORT,
                1,
                NPIX as c_long,
                cast_slice(&data),
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);
            if status != 0 {
                eprintln!("create failed: status={status}");
                return ExitCode::from(1);
            }
            ExitCode::SUCCESS
        }
        "read" => {
            let want = pattern(seed);
            let mut got = [0i16; NPIX];
            let mut anynull: c_int = 0;
            fits_open_file(&mut f, &cc(&url), READONLY, &mut status);
            fits_read_img(
                f.as_deref_mut().unwrap(),
                TSHORT,
                1,
                NPIX as c_long,
                None,
                cast_slice_mut(&mut got),
                Some(&mut anynull),
                &mut status,
            );
            let _ = f.take().map(|h| fits_close_file(h, &mut status));
            if status != 0 {
                eprintln!("read failed: status={status}");
                return ExitCode::from(1);
            }
            if got != want {
                eprintln!("MISMATCH: cross-process data differs from written pattern");
                return ExitCode::from(3);
            }
            ExitCode::SUCCESS
        }
        "delete" => {
            fits_open_file(&mut f, &cc(&url), READWRITE, &mut status);
            fits_delete_file(&mut f, &mut status);
            if status != 0 {
                eprintln!("delete failed: status={status}");
                return ExitCode::from(1);
            }
            ExitCode::SUCCESS
        }
        other => {
            eprintln!("unknown op: {other}");
            ExitCode::from(2)
        }
    }
}
