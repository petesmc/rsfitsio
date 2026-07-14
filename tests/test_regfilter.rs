//! Integration test for the `regfilter()` row filter (REGFILTER expression).
//!
//! Builds a FITS table with X/Y spatial columns and a DS9-style region file
//! (a circle in pixel coordinates), then filters the rows, keeping only those
//! whose (X, Y) falls inside the region. Exercises the `New_REG` code path and
//! its region-pointer handling, which otherwise has no coverage.

use libc::{c_char, c_int, c_long};

use rsfitsio::aliases::rust_api::*;
use rsfitsio::fitsio::{BINARY_TBL, BYTE_IMG, READONLY, fitsfile};

mod common;
use common::with_temp_file;

fn cc(s: &str) -> Vec<c_char> {
    let mut v: Vec<c_char> = s.bytes().map(|b| b as c_char).collect();
    v.push(0);
    v
}

#[test]
fn test_regfilter_circle() {
    with_temp_file(|filename| {
        let mut status: c_int = 0;
        let path = cc(filename);

        // Region file: a circle centred at (100,100) with radius 30 (pixels).
        let region_path = format!("{filename}.reg");
        std::fs::write(&region_path, "circle(100,100,30)\n").unwrap();

        // Event X/Y coordinates.
        let xs: [f64; 4] = [100.0, 120.0, 100.0, 200.0];
        let ys: [f64; 4] = [100.0, 100.0, 140.0, 200.0];

        // --- Build the file: primary + EVENTS table with X and Y ---
        let mut fptr: Option<Box<fitsfile>> = None;
        fits_create_file(&mut fptr, &path, &mut status);
        assert_eq!(status, 0, "create_file");
        {
            let f = fptr.as_deref_mut().unwrap();
            let naxes: [c_long; 0] = [];
            fits_create_img(f, BYTE_IMG, 0, &naxes, &mut status);
            assert_eq!(status, 0, "primary");

            let ttype = [Some(cc("X")), Some(cc("Y"))];
            let ttype_ref: Vec<Option<&[c_char]>> = ttype.iter().map(|o| o.as_deref()).collect();
            let tform_v = [cc("1D"), cc("1D")];
            let tform_ref: Vec<&[c_char]> = tform_v.iter().map(|v| v.as_slice()).collect();
            fits_create_tbl(
                f,
                BINARY_TBL,
                xs.len() as i64,
                2,
                &ttype_ref,
                &tform_ref,
                None,
                Some(&cc("EVENTS")),
                &mut status,
            );
            assert_eq!(status, 0, "create EVENTS");
            fits_write_col_dbl(f, 1, 1, 1, xs.len() as i64, &xs, &mut status);
            fits_write_col_dbl(f, 2, 1, 1, ys.len() as i64, &ys, &mut status);
            assert_eq!(status, 0, "write X/Y");
        }
        fits_close_file(fptr.take().unwrap(), &mut status);
        assert_eq!(status, 0, "close");

        // --- Reopen, move to EVENTS, apply regfilter("<region file>") ---
        let mut rfptr: Option<Box<fitsfile>> = None;
        fits_open_file(&mut rfptr, &path, READONLY, &mut status);
        assert_eq!(status, 0, "reopen");
        let f = rfptr.as_deref_mut().unwrap();
        let mut hdutype = 0;
        fits_movabs_hdu(f, 2, Some(&mut hdutype), &mut status);
        assert_eq!(status, 0, "move to EVENTS");

        let expr = cc(&format!("regfilter(\"{region_path}\")"));
        let mut n_good: c_long = -1;
        let mut row_status: [c_char; 4] = [0; 4];
        fits_find_rows(f, &expr, 1, 4, &mut n_good, &mut row_status, &mut status);
        assert_eq!(status, 0, "regfilter evaluation (status={status})");

        // (100,100) centre and (120,100) dist 20 are inside r=30; (100,140)
        // dist 40 and (200,200) are outside.
        assert_eq!(row_status, [1, 1, 0, 0], "matched rows");
        assert_eq!(n_good, 2, "number of good rows");

        fits_close_file(rfptr.take().unwrap(), &mut status);
        let _ = std::fs::remove_file(&region_path);
    });
}
