/*
 * Integration tests for histogram public API functions in histo.c
 * Tests real end-to-end histogram creation with actual FITS files
 * This version uses the SAFE Rust API (from rsfitsio::aliases::rust_api)
 */

use bytemuck::{cast_slice, cast_slice_mut};
use libc::{c_char, c_double, c_int, c_long};
use std::ffi::CString;
use tempfile::Builder;

use rsfitsio::aliases::rust_api::*;
use rsfitsio::fitsio::{
    BINARY_TBL, FLEN_VALUE, LONG_IMG, LONGLONG, READONLY, READWRITE, SHORT_IMG, TDOUBLE, TINT,
    TLONG, fitsfile,
};
use rsfitsio::helpers::testhelpers::with_temp_file;
use rsfitsio::histo::{
    ffhist2_safe, ffhist3_safe, fits_calc_binning_safe, fits_calc_binningd_safe,
    fits_make_hist_safe, fits_make_histd_safe, fits_rebin_wcs_safe, fits_rebin_wcsd_safe,
    fits_write_keys_histo_safe,
};
use rsfitsio::{KeywordDatatype, KeywordDatatypeMut};

/// Helper: Create a test FITS table with X, Y, Z, WEIGHT columns
fn create_test_table(filename: &str, nrows: i32) -> c_int {
    let mut fptr: Option<Box<fitsfile>> = None;
    let mut status: c_int = 0;
    let naxes: [c_long; 2] = [0, 0];

    let filename_with_bang = format!("!{}", filename);
    let filename_cstr = CString::new(filename_with_bang).unwrap();

    // Create new FITS file
    fits_create_file(
        &mut fptr,
        cast_slice(filename_cstr.to_bytes_with_nul()),
        &mut status,
    );
    if status != 0 {
        return status;
    }

    // Create empty primary array
    if let Some(ref mut fptr_box) = fptr {
        fits_create_img(fptr_box, SHORT_IMG, 0, &naxes, &mut status);
        if status != 0 {
            if let Some(fptr_box) = fptr {
                fits_close_file(fptr_box, &mut status);
            }
            return status;
        }

        // Create binary table extension
        let ttype_x = CString::new("X").unwrap();
        let ttype_y = CString::new("Y").unwrap();
        let ttype_z = CString::new("Z").unwrap();
        let ttype_w = CString::new("WEIGHT").unwrap();
        let ttype: Vec<&[u8]> = vec![
            ttype_x.as_bytes_with_nul(),
            ttype_y.as_bytes_with_nul(),
            ttype_z.as_bytes_with_nul(),
            ttype_w.as_bytes_with_nul(),
        ];

        let tform_d = CString::new("1D").unwrap();
        let tform: Vec<&[u8]> = vec![
            tform_d.as_bytes_with_nul(),
            tform_d.as_bytes_with_nul(),
            tform_d.as_bytes_with_nul(),
            tform_d.as_bytes_with_nul(),
        ];

        let tunit_empty = CString::new("").unwrap();
        let tunit: Vec<&[u8]> = vec![
            tunit_empty.as_bytes_with_nul(),
            tunit_empty.as_bytes_with_nul(),
            tunit_empty.as_bytes_with_nul(),
            tunit_empty.as_bytes_with_nul(),
        ];

        let extname = CString::new("TEST_DATA").unwrap();

        let ttype_casted: Vec<Option<&[c_char]>> =
            ttype.iter().map(|s| Some(cast_slice(s))).collect();
        let tform_casted: Vec<&[c_char]> = tform.iter().map(|s| cast_slice(s)).collect();
        let tunit_casted: Vec<Option<&[c_char]>> =
            tunit.iter().map(|s| Some(cast_slice(s))).collect();

        fits_create_tbl(
            fptr_box,
            BINARY_TBL,
            nrows as LONGLONG,
            4,
            &ttype_casted,
            &tform_casted,
            Some(&tunit_casted),
            Some(cast_slice(extname.as_bytes_with_nul())),
            &mut status,
        );
        if status != 0 {
            if let Some(fptr_box) = fptr {
                fits_close_file(fptr_box, &mut status);
            }
            return status;
        }

        // Fill table with test data
        let mut xdata: Vec<c_double> = Vec::with_capacity(nrows as usize);
        let mut ydata: Vec<c_double> = Vec::with_capacity(nrows as usize);
        let mut zdata: Vec<c_double> = Vec::with_capacity(nrows as usize);
        let mut wdata: Vec<c_double> = Vec::with_capacity(nrows as usize);

        for i in 0..nrows {
            xdata.push((i % 100) as c_double);
            ydata.push((i % 50) as c_double);
            zdata.push((i % 30) as c_double);
            wdata.push(1.0);
        }

        fits_write_col(
            fptr_box,
            TDOUBLE,
            1,
            1,
            1,
            nrows as LONGLONG,
            cast_slice(&xdata),
            &mut status,
        );
        fits_write_col(
            fptr_box,
            TDOUBLE,
            2,
            1,
            1,
            nrows as LONGLONG,
            cast_slice(&ydata),
            &mut status,
        );
        fits_write_col(
            fptr_box,
            TDOUBLE,
            3,
            1,
            1,
            nrows as LONGLONG,
            cast_slice(&zdata),
            &mut status,
        );
        fits_write_col(
            fptr_box,
            TDOUBLE,
            4,
            1,
            1,
            nrows as LONGLONG,
            cast_slice(&wdata),
            &mut status,
        );
    }

    if let Some(fptr_box) = fptr {
        fits_close_file(fptr_box, &mut status);
    }
    status
}

/// Helper: Read histogram values and verify sum
fn verify_histogram_sum(filename: &str, expected_sum: i64) -> bool {
    let mut fptr: Option<Box<fitsfile>> = None;
    let mut status: c_int = 0;
    let mut naxis: c_int = 0;
    let mut naxes: [c_long; 4] = [0; 4];
    let mut npixels: c_long = 1;

    let filename_cstr = CString::new(filename).unwrap();
    fits_open_file(
        &mut fptr,
        cast_slice(filename_cstr.to_bytes_with_nul()),
        READONLY,
        &mut status,
    );
    if status != 0 {
        return false;
    }

    if let Some(ref mut fptr_box) = fptr {
        // Get image dimensions
        fits_get_img_dim(fptr_box, &mut naxis, &mut status);
        if status != 0 {
            if let Some(fptr_box) = fptr {
                fits_close_file(fptr_box, &mut status);
            }
            return false;
        }

        fits_get_img_size(fptr_box, naxis, &mut naxes, &mut status);
        if status != 0 {
            if let Some(fptr_box) = fptr {
                fits_close_file(fptr_box, &mut status);
            }
            return false;
        }

        // Calculate total pixels
        for i in 0..naxis as usize {
            npixels *= naxes[i];
        }

        // Read all pixel values
        let mut pixels: Vec<c_long> = vec![0; npixels as usize];
        fits_read_img(
            fptr_box,
            TLONG,
            1,
            npixels.into(),
            None,
            cast_slice_mut(&mut pixels),
            None,
            &mut status,
        );
        if status != 0 {
            if let Some(fptr_box) = fptr {
                fits_close_file(fptr_box, &mut status);
            }
            return false;
        }

        let sum: i64 = pixels.iter().map(|&x| x as i64).sum();

        if let Some(fptr_box) = fptr {
            fits_close_file(fptr_box, &mut status);
        }

        return sum == expected_sum;
    }

    false
}

/// Helper for tests that need two temp files
fn with_two_temp_files<F>(callback: F)
where
    F: FnOnce(&str, &str),
{
    let tdir = Builder::new().prefix("rsfitsio-").tempdir().unwrap();
    let tdir_path = tdir.path();

    let filename1 = tdir_path.join("test1.fits");
    let filename2 = tdir_path.join("test2.fits");

    let filename1_str = filename1.to_str().expect("cannot create string filename");
    let filename2_str = filename2.to_str().expect("cannot create string filename");

    callback(filename1_str, filename2_str);
}

#[test]
fn test_fits_calc_binning_1d() {
    with_temp_file(|test_file| {
        let mut status: c_int = 0;
        let mut fptr: Option<Box<fitsfile>> = None;

        // Create test table
        assert_eq!(create_test_table(test_file, 100), 0);

        // Open table
        let filename = CString::new(test_file).unwrap();
        fits_open_file(
            &mut fptr,
            cast_slice(filename.as_bytes_with_nul()),
            READONLY,
            &mut status,
        );
        assert_eq!(status, 0);

        if let Some(ref mut fptr_box) = fptr {
            fits_movabs_hdu(fptr_box, 2, None, &mut status);
            assert_eq!(status, 0);

            // Set up binning parameters - need proper 2D array format
            let colname_x = CString::new("X").unwrap();
            let mut colname: [[c_char; FLEN_VALUE]; 4] = [[0; FLEN_VALUE]; 4];

            // Copy X into first array
            let x_bytes = colname_x.as_bytes_with_nul();
            for (i, &b) in x_bytes.iter().enumerate() {
                if i < FLEN_VALUE {
                    colname[0][i] = b as c_char;
                }
            }

            let mut minin: [f64; 4] = [0.0, 0.0, 0.0, 0.0];
            let mut maxin: [f64; 4] = [100.0, 0.0, 0.0, 0.0];
            let mut binsizein: [f64; 4] = [10.0, 0.0, 0.0, 0.0];

            let mut minname: [[c_char; FLEN_VALUE]; 4] = [[0; FLEN_VALUE]; 4];
            let mut maxname: [[c_char; FLEN_VALUE]; 4] = [[0; FLEN_VALUE]; 4];
            let mut binname: [[c_char; FLEN_VALUE]; 4] = [[0; FLEN_VALUE]; 4];

            let mut colnum: [c_int; 4] = [0; 4];
            let mut haxes: [c_long; 4] = [0; 4];
            let mut amin: [f32; 4] = [0.0; 4];
            let mut amax: [f32; 4] = [0.0; 4];
            let mut binsize: [f32; 4] = [0.0; 4];

            // Calculate binning
            fits_calc_binning_safe(
                fptr_box,
                1,
                &mut colname,
                &mut minin,
                &mut maxin,
                &mut binsizein,
                &mut minname,
                &mut maxname,
                &mut binname,
                &mut colnum,
                &mut haxes,
                &mut amin,
                &mut amax,
                &mut binsize,
                &mut status,
            );
            assert_eq!(status, 0);

            // Verify results
            assert_eq!(haxes[0], 10); // 100/10 = 10 bins
            assert!((amin[0] - 0.0).abs() < 0.001);
            assert!((amax[0] - 100.0).abs() < 0.001);
            assert!((binsize[0] - 10.0).abs() < 0.001);
        }

        if let Some(fptr_box) = fptr {
            fits_close_file(fptr_box, &mut status);
        }
    });
}

#[test]
fn test_fits_calc_binningd_2d() {
    with_temp_file(|test_file| {
        let mut status: c_int = 0;
        let mut fptr: Option<Box<fitsfile>> = None;

        // Create test table
        assert_eq!(create_test_table(test_file, 100), 0);

        // Open table
        let filename = CString::new(test_file).unwrap();
        fits_open_file(
            &mut fptr,
            cast_slice(filename.as_bytes_with_nul()),
            READONLY,
            &mut status,
        );
        assert_eq!(status, 0);

        if let Some(ref mut fptr_box) = fptr {
            fits_movabs_hdu(fptr_box, 2, None, &mut status);
            assert_eq!(status, 0);

            // Set up binning parameters for 2D histogram
            let colname_x = CString::new("X").unwrap();
            let colname_y = CString::new("Y").unwrap();
            let mut colname: [[c_char; FLEN_VALUE]; 4] = [[0; FLEN_VALUE]; 4];

            // Copy X and Y into arrays
            let x_bytes = colname_x.as_bytes_with_nul();
            for (i, &b) in x_bytes.iter().enumerate() {
                if i < FLEN_VALUE {
                    colname[0][i] = b as c_char;
                }
            }
            let y_bytes = colname_y.as_bytes_with_nul();
            for (i, &b) in y_bytes.iter().enumerate() {
                if i < FLEN_VALUE {
                    colname[1][i] = b as c_char;
                }
            }

            let mut minin: [f64; 4] = [0.0, 0.0, 0.0, 0.0];
            let mut maxin: [f64; 4] = [100.0, 50.0, 0.0, 0.0];
            let mut binsizein: [f64; 4] = [10.0, 10.0, 0.0, 0.0];

            let minname: [[c_char; FLEN_VALUE]; 4] = [[0; FLEN_VALUE]; 4];
            let maxname: [[c_char; FLEN_VALUE]; 4] = [[0; FLEN_VALUE]; 4];
            let binname: [[c_char; FLEN_VALUE]; 4] = [[0; FLEN_VALUE]; 4];

            let mut colnum: [c_int; 4] = [0; 4];
            let mut haxes: [c_long; 4] = [0; 4];
            let mut amin: [f64; 4] = [0.0; 4];
            let mut amax: [f64; 4] = [0.0; 4];
            let mut binsize: [f64; 4] = [0.0; 4];

            // Calculate binning using double version
            fits_calc_binningd_safe(
                fptr_box,
                2,
                &mut colname,
                &mut minin,
                &mut maxin,
                &mut binsizein,
                &minname,
                &maxname,
                &binname,
                &mut colnum,
                &mut haxes,
                &mut amin,
                &mut amax,
                &mut binsize,
                &mut status,
            );
            assert_eq!(status, 0);

            // Verify results
            assert_eq!(haxes[0], 10); // 100/10 = 10 bins
            assert_eq!(haxes[1], 5); // 50/10 = 5 bins
            assert!((amin[0] - 0.0).abs() < 0.001);
            assert!((amin[1] - 0.0).abs() < 0.001);
        }

        if let Some(fptr_box) = fptr {
            fits_close_file(fptr_box, &mut status);
        }
    });
}

#[test]
fn test_ffhist2_1d() {
    with_two_temp_files(|test_file, hist_file| {
        let mut status: c_int = 0;
        let mut fptr: Option<Box<fitsfile>> = None;

        // Create test table
        assert_eq!(create_test_table(test_file, 100), 0);

        // Open table
        let filename = CString::new(test_file).unwrap();
        fits_open_file(
            &mut fptr,
            cast_slice(filename.as_bytes_with_nul()),
            READWRITE,
            &mut status,
        );
        assert_eq!(status, 0);

        if let Some(ref mut fptr_box) = fptr {
            fits_movabs_hdu(fptr_box, 2, None, &mut status);
            assert_eq!(status, 0);

            // Set up histogram parameters
            let colname_x = CString::new("X").unwrap();
            let mut colname: [[c_char; FLEN_VALUE]; 4] = [[0; FLEN_VALUE]; 4];
            let x_bytes = colname_x.as_bytes_with_nul();
            for (i, &b) in x_bytes.iter().enumerate() {
                if i < FLEN_VALUE {
                    colname[0][i] = b as c_char;
                }
            }

            let minin: [f64; 4] = [0.0, 0.0, 0.0, 0.0];
            let maxin: [f64; 4] = [100.0, 0.0, 0.0, 0.0];
            let binsizein: [f64; 4] = [10.0, 0.0, 0.0, 0.0];
            let minname: [[c_char; FLEN_VALUE]; 4] = [[0; FLEN_VALUE]; 4];
            let maxname: [[c_char; FLEN_VALUE]; 4] = [[0; FLEN_VALUE]; 4];
            let binname: [[c_char; FLEN_VALUE]; 4] = [[0; FLEN_VALUE]; 4];
            let wtcol: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];

            let hist_filename_str = format!("!{}", hist_file);
            let hist_filename = CString::new(hist_filename_str).unwrap();

            // Call ffhist2 to create 1D histogram - fptr will be reassigned to histogram
            ffhist2_safe(
                &mut fptr,
                cast_slice(hist_filename.as_bytes_with_nul()),
                TINT,
                1,
                &colname,
                &minin,
                &maxin,
                &binsizein,
                &minname,
                &maxname,
                &binname,
                1.0,
                &wtcol,
                0,
                None,
                &mut status,
            );
            assert_eq!(status, 0);
        }

        // fptr now points to histogram
        if let Some(fptr_box) = fptr {
            fits_close_file(fptr_box, &mut status);
        }

        // Verify histogram sum
        assert!(verify_histogram_sum(hist_file, 100));
    });
}

#[test]
fn test_ffhist2_2d() {
    with_two_temp_files(|test_file, hist_file| {
        let mut status: c_int = 0;
        let mut fptr: Option<Box<fitsfile>> = None;

        // Create test table
        assert_eq!(create_test_table(test_file, 100), 0);

        // Open table
        let filename = CString::new(test_file).unwrap();
        fits_open_file(
            &mut fptr,
            cast_slice(filename.as_bytes_with_nul()),
            READWRITE,
            &mut status,
        );
        assert_eq!(status, 0);

        if let Some(ref mut fptr_box) = fptr {
            fits_movabs_hdu(fptr_box, 2, None, &mut status);
            assert_eq!(status, 0);

            // Set up 2D histogram parameters
            let colname_x = CString::new("X").unwrap();
            let colname_y = CString::new("Y").unwrap();
            let mut colname: [[c_char; FLEN_VALUE]; 4] = [[0; FLEN_VALUE]; 4];

            let x_bytes = colname_x.as_bytes_with_nul();
            for (i, &b) in x_bytes.iter().enumerate() {
                if i < FLEN_VALUE {
                    colname[0][i] = b as c_char;
                }
            }
            let y_bytes = colname_y.as_bytes_with_nul();
            for (i, &b) in y_bytes.iter().enumerate() {
                if i < FLEN_VALUE {
                    colname[1][i] = b as c_char;
                }
            }

            let minin: [f64; 4] = [0.0, 0.0, 0.0, 0.0];
            let maxin: [f64; 4] = [100.0, 50.0, 0.0, 0.0];
            let binsizein: [f64; 4] = [10.0, 10.0, 0.0, 0.0];
            let minname: [[c_char; FLEN_VALUE]; 4] = [[0; FLEN_VALUE]; 4];
            let maxname: [[c_char; FLEN_VALUE]; 4] = [[0; FLEN_VALUE]; 4];
            let binname: [[c_char; FLEN_VALUE]; 4] = [[0; FLEN_VALUE]; 4];
            let wtcol: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];

            let hist_filename_str = format!("!{}", hist_file);
            let hist_filename = CString::new(hist_filename_str).unwrap();

            // Call ffhist2 to create 2D histogram
            ffhist2_safe(
                &mut fptr,
                cast_slice(hist_filename.as_bytes_with_nul()),
                TINT,
                2,
                &colname,
                &minin,
                &maxin,
                &binsizein,
                &minname,
                &maxname,
                &binname,
                1.0,
                &wtcol,
                0,
                None,
                &mut status,
            );
            assert_eq!(status, 0);
        }

        // fptr now points to histogram
        if let Some(fptr_box) = fptr {
            fits_close_file(fptr_box, &mut status);
        }

        // Verify histogram sum
        assert!(verify_histogram_sum(hist_file, 100));
    });
}

#[test]
fn test_ffhist2_large_dataset() {
    with_two_temp_files(|test_file, hist_file| {
        let mut status: c_int = 0;
        let mut fptr: Option<Box<fitsfile>> = None;

        // Create test table with 200 rows
        assert_eq!(create_test_table(test_file, 200), 0);

        // Open table
        let filename = CString::new(test_file).unwrap();
        fits_open_file(
            &mut fptr,
            cast_slice(filename.as_bytes_with_nul()),
            READWRITE,
            &mut status,
        );
        assert_eq!(status, 0);

        if let Some(ref mut fptr_box) = fptr {
            fits_movabs_hdu(fptr_box, 2, None, &mut status);
            assert_eq!(status, 0);

            // Set up histogram parameters
            let colname_x = CString::new("X").unwrap();
            let mut colname: [[c_char; FLEN_VALUE]; 4] = [[0; FLEN_VALUE]; 4];
            let x_bytes = colname_x.as_bytes_with_nul();
            for (i, &b) in x_bytes.iter().enumerate() {
                if i < FLEN_VALUE {
                    colname[0][i] = b as c_char;
                }
            }

            let minin: [f64; 4] = [0.0, 0.0, 0.0, 0.0];
            let maxin: [f64; 4] = [100.0, 0.0, 0.0, 0.0];
            let binsizein: [f64; 4] = [10.0, 0.0, 0.0, 0.0];
            let minname: [[c_char; FLEN_VALUE]; 4] = [[0; FLEN_VALUE]; 4];
            let maxname: [[c_char; FLEN_VALUE]; 4] = [[0; FLEN_VALUE]; 4];
            let binname: [[c_char; FLEN_VALUE]; 4] = [[0; FLEN_VALUE]; 4];
            let wtcol: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];

            let hist_filename_str = format!("!{}", hist_file);
            let hist_filename = CString::new(hist_filename_str).unwrap();

            // Call ffhist2 to create 1D histogram
            ffhist2_safe(
                &mut fptr,
                cast_slice(hist_filename.as_bytes_with_nul()),
                TINT,
                1,
                &colname,
                &minin,
                &maxin,
                &binsizein,
                &minname,
                &maxname,
                &binname,
                1.0,
                &wtcol,
                0,
                None,
                &mut status,
            );
            assert_eq!(status, 0);
        }

        // fptr now points to histogram
        if let Some(fptr_box) = fptr {
            fits_close_file(fptr_box, &mut status);
        }

        // Verify histogram sum - 200 rows
        assert!(verify_histogram_sum(hist_file, 200));
    });
}

#[test]
fn test_ffhist2_with_weight() {
    with_two_temp_files(|test_file, hist_file| {
        let mut status: c_int = 0;
        let mut fptr: Option<Box<fitsfile>> = None;

        // Create test table
        assert_eq!(create_test_table(test_file, 100), 0);

        // Open table
        let filename = CString::new(test_file).unwrap();
        fits_open_file(
            &mut fptr,
            cast_slice(filename.as_bytes_with_nul()),
            READWRITE,
            &mut status,
        );
        assert_eq!(status, 0);

        if let Some(ref mut fptr_box) = fptr {
            fits_movabs_hdu(fptr_box, 2, None, &mut status);
            assert_eq!(status, 0);

            // Set up histogram parameters with weight column
            let colname_x = CString::new("X").unwrap();
            let mut colname: [[c_char; FLEN_VALUE]; 4] = [[0; FLEN_VALUE]; 4];
            let x_bytes = colname_x.as_bytes_with_nul();
            for (i, &b) in x_bytes.iter().enumerate() {
                if i < FLEN_VALUE {
                    colname[0][i] = b as c_char;
                }
            }

            let minin: [f64; 4] = [0.0, 0.0, 0.0, 0.0];
            let maxin: [f64; 4] = [100.0, 0.0, 0.0, 0.0];
            let binsizein: [f64; 4] = [10.0, 0.0, 0.0, 0.0];
            let minname: [[c_char; FLEN_VALUE]; 4] = [[0; FLEN_VALUE]; 4];
            let maxname: [[c_char; FLEN_VALUE]; 4] = [[0; FLEN_VALUE]; 4];
            let binname: [[c_char; FLEN_VALUE]; 4] = [[0; FLEN_VALUE]; 4];

            let weight_str = CString::new("WEIGHT").unwrap();
            let mut wtcol: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
            let weight_bytes = weight_str.as_bytes_with_nul();
            for (i, &b) in weight_bytes.iter().enumerate() {
                if i < FLEN_VALUE {
                    wtcol[i] = b as c_char;
                }
            }

            let hist_filename_str = format!("!{}", hist_file);
            let hist_filename = CString::new(hist_filename_str).unwrap();

            // Create weighted histogram
            ffhist2_safe(
                &mut fptr,
                cast_slice(hist_filename.as_bytes_with_nul()),
                TINT,
                1,
                &colname,
                &minin,
                &maxin,
                &binsizein,
                &minname,
                &maxname,
                &binname,
                1.0,
                &wtcol,
                0,
                None,
                &mut status,
            );
            assert_eq!(status, 0);
        }

        if let Some(fptr_box) = fptr {
            fits_close_file(fptr_box, &mut status);
        }

        // With weight=1.0 for all rows, sum should still be 100
        assert!(verify_histogram_sum(hist_file, 100));
    });
}

#[test]
fn test_ffhist3_1d() {
    with_two_temp_files(|test_file, hist_file| {
        let mut status: c_int = 0;
        let mut fptr: Option<Box<fitsfile>> = None;

        // Create test table with 150 rows
        assert_eq!(create_test_table(test_file, 150), 0);

        // Open table
        let filename = CString::new(test_file).unwrap();
        fits_open_file(
            &mut fptr,
            cast_slice(filename.as_bytes_with_nul()),
            READONLY,
            &mut status,
        );
        assert_eq!(status, 0);

        let mut histptr: Option<Box<fitsfile>> = None;

        if let Some(ref mut fptr_box) = fptr {
            fits_movabs_hdu(fptr_box, 2, None, &mut status);
            assert_eq!(status, 0);

            // Set up histogram parameters
            let colname_x = CString::new("X").unwrap();
            let mut colname: [[c_char; FLEN_VALUE]; 4] = [[0; FLEN_VALUE]; 4];
            let x_bytes = colname_x.as_bytes_with_nul();
            for (i, &b) in x_bytes.iter().enumerate() {
                if i < FLEN_VALUE {
                    colname[0][i] = b as c_char;
                }
            }

            let minin: [f64; 4] = [0.0, 0.0, 0.0, 0.0];
            let maxin: [f64; 4] = [100.0, 0.0, 0.0, 0.0];
            let binsizein: [f64; 4] = [10.0, 0.0, 0.0, 0.0];
            let minname: [[c_char; FLEN_VALUE]; 4] = [[0; FLEN_VALUE]; 4];
            let maxname: [[c_char; FLEN_VALUE]; 4] = [[0; FLEN_VALUE]; 4];
            let binname: [[c_char; FLEN_VALUE]; 4] = [[0; FLEN_VALUE]; 4];
            let wtcol: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];

            let hist_filename_str = format!("!{}", hist_file);
            let hist_filename = CString::new(hist_filename_str).unwrap();

            // Create histogram - ffhist3 returns a new fitsfile pointer
            histptr = ffhist3_safe(
                fptr_box,
                cast_slice(hist_filename.as_bytes_with_nul()),
                TINT,
                1,
                &colname,
                &minin,
                &maxin,
                &binsizein,
                &minname,
                &maxname,
                &binname,
                1.0,
                &wtcol,
                0,
                None,
                &mut status,
            );
            assert_eq!(status, 0);
            assert!(histptr.is_some());
        }

        // Close both file pointers
        if let Some(fptr_box) = fptr {
            fits_close_file(fptr_box, &mut status);
        }

        if let Some(histptr_box) = histptr {
            fits_close_file(histptr_box, &mut status);
        }

        // Verify histogram sum - 150 rows
        assert!(verify_histogram_sum(hist_file, 150));
    });
}

#[test]
fn test_fits_make_hist_1d() {
    with_two_temp_files(|test_file, hist_file| {
        let mut status: c_int = 0;
        let mut fptr: Option<Box<fitsfile>> = None;
        let mut histptr: Option<Box<fitsfile>> = None;

        // Create test table
        assert_eq!(create_test_table(test_file, 100), 0);

        // Open table
        let filename = CString::new(test_file).unwrap();
        fits_open_file(
            &mut fptr,
            cast_slice(filename.as_bytes_with_nul()),
            READONLY,
            &mut status,
        );
        assert_eq!(status, 0);

        if let Some(ref mut fptr_box) = fptr {
            fits_movabs_hdu(fptr_box, 2, None, &mut status);
            assert_eq!(status, 0);

            // Create output histogram file
            let hist_filename_str = format!("!{}", hist_file);
            let hist_filename = CString::new(hist_filename_str).unwrap();
            fits_create_file(
                &mut histptr,
                cast_slice(hist_filename.as_bytes_with_nul()),
                &mut status,
            );
            assert_eq!(status, 0);

            if let Some(ref mut histptr_box) = histptr {
                // Set up histogram parameters
                let colnum: [c_int; 4] = [1, 0, 0, 0]; // Column 1 = X
                let naxes: [c_long; 1] = [10]; // 10 bins
                let amin: [f32; 1] = [0.0];
                let amax: [f32; 1] = [100.0];
                let binsize: [f32; 1] = [10.0];

                // Create the output image HDU first
                fits_create_img(histptr_box, LONG_IMG, 1, &naxes, &mut status);
                assert_eq!(status, 0);

                // Create histogram using low-level function
                fits_make_hist_safe(
                    fptr_box,
                    histptr_box,
                    LONG_IMG,
                    1,
                    &naxes,
                    &colnum,
                    &amin,
                    &amax,
                    &binsize,
                    1.0,
                    0,
                    0,
                    None,
                    &mut status,
                );
                assert_eq!(status, 0);
            }
        }

        if let Some(fptr_box) = fptr {
            fits_close_file(fptr_box, &mut status);
        }

        if let Some(histptr_box) = histptr {
            fits_close_file(histptr_box, &mut status);
        }

        // Verify histogram sum
        assert!(verify_histogram_sum(hist_file, 100));
    });
}

#[test]
fn test_fits_make_histd_2d() {
    with_two_temp_files(|test_file, hist_file| {
        let mut status: c_int = 0;
        let mut fptr: Option<Box<fitsfile>> = None;
        let mut histptr: Option<Box<fitsfile>> = None;

        // Create test table
        assert_eq!(create_test_table(test_file, 100), 0);

        // Open table
        let filename = CString::new(test_file).unwrap();
        fits_open_file(
            &mut fptr,
            cast_slice(filename.as_bytes_with_nul()),
            READONLY,
            &mut status,
        );
        assert_eq!(status, 0);

        if let Some(ref mut fptr_box) = fptr {
            fits_movabs_hdu(fptr_box, 2, None, &mut status);
            assert_eq!(status, 0);

            // Create output histogram file
            let hist_filename_str = format!("!{}", hist_file);
            let hist_filename = CString::new(hist_filename_str).unwrap();
            fits_create_file(
                &mut histptr,
                cast_slice(hist_filename.as_bytes_with_nul()),
                &mut status,
            );
            assert_eq!(status, 0);

            if let Some(ref mut histptr_box) = histptr {
                // Set up 2D histogram parameters
                let colnum: [c_int; 4] = [1, 2, 0, 0]; // Columns 1,2 = X,Y
                let naxes: [c_long; 2] = [10, 5]; // 10x5 bins
                let amin: [f64; 2] = [0.0, 0.0];
                let amax: [f64; 2] = [100.0, 50.0];
                let binsize: [f64; 2] = [10.0, 10.0];

                // Create the output image HDU first
                fits_create_img(histptr_box, LONG_IMG, 2, &naxes, &mut status);
                assert_eq!(status, 0);

                // Create 2D histogram using double version
                fits_make_histd_safe(
                    fptr_box,
                    histptr_box,
                    LONG_IMG,
                    2,
                    &naxes,
                    &colnum,
                    &amin,
                    &amax,
                    &binsize,
                    1.0,
                    0,
                    0,
                    None,
                    &mut status,
                );
                assert_eq!(status, 0);
            }
        }

        if let Some(fptr_box) = fptr {
            fits_close_file(fptr_box, &mut status);
        }

        if let Some(histptr_box) = histptr {
            fits_close_file(histptr_box, &mut status);
        }

        // Verify histogram sum
        assert!(verify_histogram_sum(hist_file, 100));
    });
}

#[test]
fn test_fits_write_keys_histo() {
    with_two_temp_files(|test_file, hist_file| {
        let mut status: c_int = 0;
        let mut fptr: Option<Box<fitsfile>> = None;
        let mut histptr: Option<Box<fitsfile>> = None;

        // Create test table
        assert_eq!(create_test_table(test_file, 100), 0);

        // Open table
        let filename = CString::new(test_file).unwrap();
        fits_open_file(
            &mut fptr,
            cast_slice(filename.as_bytes_with_nul()),
            READONLY,
            &mut status,
        );
        assert_eq!(status, 0);

        if let Some(ref mut fptr_box) = fptr {
            fits_movabs_hdu(fptr_box, 2, None, &mut status);
            assert_eq!(status, 0);

            // Create output histogram file manually
            let hist_filename_str = format!("!{}", hist_file);
            let hist_filename = CString::new(hist_filename_str).unwrap();
            fits_create_file(
                &mut histptr,
                cast_slice(hist_filename.as_bytes_with_nul()),
                &mut status,
            );
            assert_eq!(status, 0);

            if let Some(ref mut histptr_box) = histptr {
                // Create a 1D image
                let naxes: [c_long; 1] = [10];
                fits_create_img(histptr_box, LONG_IMG, 1, &naxes, &mut status);
                assert_eq!(status, 0);

                // Write histogram keywords
                let colnum: [c_int; 4] = [1, 0, 0, 0];
                fits_write_keys_histo_safe(fptr_box, histptr_box, 1, &colnum, &mut status);
                assert_eq!(status, 0);
            }
        }

        // Success if no errors occurred
        if let Some(fptr_box) = fptr {
            fits_close_file(fptr_box, &mut status);
        }

        if let Some(histptr_box) = histptr {
            fits_close_file(histptr_box, &mut status);
        }

        assert_eq!(status, 0);
    });
}

#[test]
fn test_fits_rebin_wcsd() {
    with_temp_file(|hist_file| {
        let mut status: c_int = 0;
        let mut fptr: Option<Box<fitsfile>> = None;
        let naxes: [c_long; 2] = [10, 10];

        // Create a simple image with WCS keywords
        let filename_str = format!("!{}", hist_file);
        let filename = CString::new(filename_str).unwrap();
        fits_create_file(
            &mut fptr,
            cast_slice(filename.as_bytes_with_nul()),
            &mut status,
        );
        assert_eq!(status, 0);

        if let Some(ref mut fptr_box) = fptr {
            fits_create_img(fptr_box, LONG_IMG, 2, &naxes, &mut status);
            assert_eq!(status, 0);

            // Write some basic WCS keywords
            let crpix1: c_double = 5.0;
            let crpix2: c_double = 5.0;
            let crval1: c_double = 100.0;
            let crval2: c_double = 50.0;
            let cdelt1: c_double = 1.0;
            let cdelt2: c_double = 1.0;

            let key_crpix1 = CString::new("CRPIX1").unwrap();
            let key_crpix2 = CString::new("CRPIX2").unwrap();
            let key_crval1 = CString::new("CRVAL1").unwrap();
            let key_crval2 = CString::new("CRVAL2").unwrap();
            let key_cdelt1 = CString::new("CDELT1").unwrap();
            let key_cdelt2 = CString::new("CDELT2").unwrap();
            let comment = CString::new("Reference pixel").unwrap();
            let comment2 = CString::new("Reference value").unwrap();
            let comment3 = CString::new("Pixel scale").unwrap();

            fits_write_key(
                fptr_box,
                KeywordDatatype::TDOUBLE(&crpix1),
                cast_slice(key_crpix1.as_bytes_with_nul()),
                Some(cast_slice(comment.as_bytes_with_nul())),
                &mut status,
            );
            fits_write_key(
                fptr_box,
                KeywordDatatype::TDOUBLE(&crpix2),
                cast_slice(key_crpix2.as_bytes_with_nul()),
                Some(cast_slice(comment.as_bytes_with_nul())),
                &mut status,
            );
            fits_write_key(
                fptr_box,
                KeywordDatatype::TDOUBLE(&crval1),
                cast_slice(key_crval1.as_bytes_with_nul()),
                Some(cast_slice(comment2.as_bytes_with_nul())),
                &mut status,
            );
            fits_write_key(
                fptr_box,
                KeywordDatatype::TDOUBLE(&crval2),
                cast_slice(key_crval2.as_bytes_with_nul()),
                Some(cast_slice(comment2.as_bytes_with_nul())),
                &mut status,
            );
            fits_write_key(
                fptr_box,
                KeywordDatatype::TDOUBLE(&cdelt1),
                cast_slice(key_cdelt1.as_bytes_with_nul()),
                Some(cast_slice(comment3.as_bytes_with_nul())),
                &mut status,
            );
            fits_write_key(
                fptr_box,
                KeywordDatatype::TDOUBLE(&cdelt2),
                cast_slice(key_cdelt2.as_bytes_with_nul()),
                Some(cast_slice(comment3.as_bytes_with_nul())),
                &mut status,
            );
            assert_eq!(status, 0);

            // Rebin WCS with new binning
            let mut amin: [f64; 2] = [0.0, 0.0];
            let mut binsize: [f64; 2] = [2.0, 2.0]; // 2x larger bins

            fits_rebin_wcsd_safe(fptr_box, 2, &mut amin, &mut binsize, &mut status);
            assert_eq!(status, 0);

            // Verify CDELT was updated (should be 2x larger)
            let mut new_cdelt1: c_double = 0.0;
            let mut new_cdelt2: c_double = 0.0;
            fits_read_key(
                fptr_box,
                KeywordDatatypeMut::TDOUBLE(&mut new_cdelt1),
                cast_slice(key_cdelt1.as_bytes_with_nul()),
                None,
                &mut status,
            );
            fits_read_key(
                fptr_box,
                KeywordDatatypeMut::TDOUBLE(&mut new_cdelt2),
                cast_slice(key_cdelt2.as_bytes_with_nul()),
                None,
                &mut status,
            );
            assert_eq!(status, 0);

            // Verify cdelt values doubled
            assert!((new_cdelt1 - 2.0).abs() < 0.001);
            assert!((new_cdelt2 - 2.0).abs() < 0.001);
        }

        if let Some(fptr_box) = fptr {
            fits_close_file(fptr_box, &mut status);
        }
    });
}

#[test]
fn test_fits_rebin_wcs() {
    with_temp_file(|hist_file| {
        let mut status: c_int = 0;
        let mut fptr: Option<Box<fitsfile>> = None;
        let naxes: [c_long; 1] = [20];

        // Create a simple 1D image with WCS keywords
        let filename_str = format!("!{}", hist_file);
        let filename = CString::new(filename_str).unwrap();
        fits_create_file(
            &mut fptr,
            cast_slice(filename.as_bytes_with_nul()),
            &mut status,
        );
        assert_eq!(status, 0);

        if let Some(ref mut fptr_box) = fptr {
            fits_create_img(fptr_box, LONG_IMG, 1, &naxes, &mut status);
            assert_eq!(status, 0);

            // Write basic WCS keywords
            let crpix1: c_double = 10.0;
            let crval1: c_double = 100.0;
            let cdelt1: c_double = 1.0;

            let key_crpix1 = CString::new("CRPIX1").unwrap();
            let key_crval1 = CString::new("CRVAL1").unwrap();
            let key_cdelt1 = CString::new("CDELT1").unwrap();
            let comment = CString::new("Reference pixel").unwrap();
            let comment2 = CString::new("Reference value").unwrap();
            let comment3 = CString::new("Pixel scale").unwrap();

            fits_write_key(
                fptr_box,
                KeywordDatatype::TDOUBLE(&crpix1),
                cast_slice(key_crpix1.as_bytes_with_nul()),
                Some(cast_slice(comment.as_bytes_with_nul())),
                &mut status,
            );
            fits_write_key(
                fptr_box,
                KeywordDatatype::TDOUBLE(&crval1),
                cast_slice(key_crval1.as_bytes_with_nul()),
                Some(cast_slice(comment2.as_bytes_with_nul())),
                &mut status,
            );
            fits_write_key(
                fptr_box,
                KeywordDatatype::TDOUBLE(&cdelt1),
                cast_slice(key_cdelt1.as_bytes_with_nul()),
                Some(cast_slice(comment3.as_bytes_with_nul())),
                &mut status,
            );
            assert_eq!(status, 0);

            // Rebin WCS with new binning (float version)
            let mut amin: [f32; 1] = [0.0];
            let mut binsize: [f32; 1] = [2.5];

            fits_rebin_wcs_safe(fptr_box, 1, &mut amin, &mut binsize, &mut status);
            assert_eq!(status, 0);

            // Verify CDELT was updated
            let mut new_cdelt1: c_double = 0.0;
            fits_read_key(
                fptr_box,
                KeywordDatatypeMut::TDOUBLE(&mut new_cdelt1),
                cast_slice(key_cdelt1.as_bytes_with_nul()),
                None,
                &mut status,
            );
            assert_eq!(status, 0);

            // CDELT should be scaled by binsize
            assert!((new_cdelt1 - 2.5).abs() < 0.001);
        }

        if let Some(fptr_box) = fptr {
            fits_close_file(fptr_box, &mut status);
        }
    });
}
