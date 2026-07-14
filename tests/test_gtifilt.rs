//! Integration test for the `gtifilt()` row filter (GTIFILTER expression).
//!
//! Builds a FITS file with an events table (a `TIME` column) plus a Good-Time-
//! Interval extension (`START`/`STOP` columns), then filters the events, keeping
//! only rows whose `TIME` falls inside one of the intervals. This exercises the
//! `New_GTI` code path and its per-node result buffer, which otherwise has no
//! coverage.

use libc::{c_char, c_int, c_long};

use rsfitsio::aliases::rust_api::*;
use rsfitsio::fitsio::{BINARY_TBL, BYTE_IMG, READONLY, fitsfile};

mod common;
use common::with_temp_file;

/// NUL-terminated `c_char` buffer from a `str`.
fn cc(s: &str) -> Vec<c_char> {
    let mut v: Vec<c_char> = s.bytes().map(|b| b as c_char).collect();
    v.push(0);
    v
}

#[test]
fn test_gtifilt_basic() {
    with_temp_file(|filename| {
        let mut status: c_int = 0;
        let path = cc(filename);

        // TIME values of the six events.
        let times: [f64; 6] = [0.5, 1.5, 2.5, 3.5, 4.5, 5.5];
        // Good-time intervals [1, 2] and [4, 5].
        let start: [f64; 2] = [1.0, 4.0];
        let stop: [f64; 2] = [2.0, 5.0];

        // --- Build the file: primary + EVENTS table + GTI extension ---
        let mut fptr: Option<Box<fitsfile>> = None;
        fits_create_file(&mut fptr, &path, &mut status);
        assert_eq!(status, 0, "create_file");

        {
            let f = fptr.as_deref_mut().unwrap();
            let naxes: [c_long; 0] = [];
            fits_create_img(f, BYTE_IMG, 0, &naxes, &mut status);
            assert_eq!(status, 0, "create null primary");
        }

        {
            let f = fptr.as_deref_mut().unwrap();
            let ttype = [Some(cc("TIME"))];
            let ttype_ref: Vec<Option<&[c_char]>> = ttype.iter().map(|o| o.as_deref()).collect();
            let tform_v = [cc("1D")];
            let tform_ref: Vec<&[c_char]> = tform_v.iter().map(|v| v.as_slice()).collect();
            fits_create_tbl(
                f,
                BINARY_TBL,
                times.len() as i64,
                1,
                &ttype_ref,
                &tform_ref,
                None,
                Some(&cc("EVENTS")),
                &mut status,
            );
            assert_eq!(status, 0, "create EVENTS table");
            fits_write_col_dbl(f, 1, 1, 1, times.len() as i64, &times, &mut status);
            assert_eq!(status, 0, "write TIME column");
        }

        {
            let f = fptr.as_deref_mut().unwrap();
            let ttype = [Some(cc("START")), Some(cc("STOP"))];
            let ttype_ref: Vec<Option<&[c_char]>> = ttype.iter().map(|o| o.as_deref()).collect();
            let tform_v = [cc("1D"), cc("1D")];
            let tform_ref: Vec<&[c_char]> = tform_v.iter().map(|v| v.as_slice()).collect();
            fits_create_tbl(
                f,
                BINARY_TBL,
                start.len() as i64,
                2,
                &ttype_ref,
                &tform_ref,
                None,
                Some(&cc("GTI")),
                &mut status,
            );
            assert_eq!(status, 0, "create GTI table");
            fits_write_col_dbl(f, 1, 1, 1, start.len() as i64, &start, &mut status);
            fits_write_col_dbl(f, 2, 1, 1, stop.len() as i64, &stop, &mut status);
            assert_eq!(status, 0, "write START/STOP columns");
        }

        fits_close_file(fptr.take().unwrap(), &mut status);
        assert_eq!(status, 0, "close after build");

        // --- Reopen, move to EVENTS, and apply gtifilt() ---
        let mut rfptr: Option<Box<fitsfile>> = None;
        fits_open_file(&mut rfptr, &path, READONLY, &mut status);
        assert_eq!(status, 0, "reopen");

        let f = rfptr.as_deref_mut().unwrap();
        let mut hdutype = 0;
        fits_movabs_hdu(f, 2, Some(&mut hdutype), &mut status);
        assert_eq!(status, 0, "move to EVENTS");

        let expr = cc("gtifilter()");
        let mut n_good: c_long = -1;
        let mut row_status: [c_char; 6] = [0; 6];
        fits_find_rows(f, &expr, 1, 6, &mut n_good, &mut row_status, &mut status);
        assert_eq!(status, 0, "gtifilt() evaluation (status={status})");

        // TIME 1.5 is in [1, 2] and 4.5 is in [4, 5]; the rest fall outside.
        assert_eq!(row_status, [0, 1, 0, 0, 1, 0], "matched rows");
        assert_eq!(n_good, 2, "number of good rows");

        fits_close_file(rfptr.take().unwrap(), &mut status);
        assert_eq!(status, 0, "close after filter");
    });
}
