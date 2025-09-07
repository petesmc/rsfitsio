use std::ffi::{CStr, CString};

use bytemuck::{cast_slice, cast_slice_mut};
use libc::{c_char, c_float, c_int, c_long, c_ushort};

use rsfitsio::aliases::rust_api::*;
use rsfitsio::fitsio::{ASCII_TBL, KEY_NO_EXIST, LONGLONG, TUSHORT, USHORT_IMG};
use rsfitsio::fitsio::{BINARY_TBL, FLEN_VALUE, READONLY, READWRITE, TFLOAT, TLONG, fitsfile};
use rsfitsio::helpers::testhelpers::with_temp_file;
use rsfitsio::{KeywordDatatype, NullValue};

fn writeimage(filename: &str) {
    let mut fptr: Option<Box<fitsfile>> = None;
    let mut status: c_int = 0;
    let fpixel: c_long = 1;

    /* initialize FITS image parameters */
    let bitpix: c_int = USHORT_IMG; /* 16-bit unsigned short pixel values       */
    let naxis: c_long = 2; /* 2-dimensional image                            */
    let naxes: [c_long; 2] = [300, 200]; /* image is 300 pixels wide by 200 rows */
    let nelements: c_long = naxes[0] * naxes[1];

    /* allocate memory for the whole image */
    let array_size = (naxes[0] * naxes[1]) as usize;
    let mut array: Vec<c_ushort> = vec![0; array_size];

    let _ = std::fs::remove_file(filename); /* Delete old file if it already exists */

    let filename_cstr = CString::new(filename).unwrap();

    unsafe {
        fits_create_file(
            &mut fptr,
            cast_slice(filename_cstr.to_bytes_with_nul()),
            &mut status,
        );
    }
    assert_eq!(status, 0);

    /* write the required keywords for the primary array image.     */
    /* Since bitpix = USHORT_IMG, this will cause cfitsio to create */
    /* a FITS image with BITPIX = 16 (signed short integers) with   */
    /* BSCALE = 1.0 and BZERO = 32768.  This is the convention that */
    /* FITS uses to store unsigned integers.  Note that the BSCALE  */
    /* and BZERO keywords will be automatically written by cfitsio  */
    /* in this case.                                                */

    if let Some(ref mut fptr_box) = fptr {
        unsafe { fits_create_img(fptr_box, bitpix, naxis as c_int, &naxes, &mut status) };
        assert_eq!(status, 0);
    }

    /* initialize the values in the image with a linear ramp function */
    for jj in 0..naxes[1] {
        for ii in 0..naxes[0] {
            let index = (jj * naxes[0] + ii) as usize;
            array[index] = (ii + jj) as c_ushort;
        }
    }

    /* write the array of unsigned integers to the FITS file */
    if let Some(ref mut fptr_box) = fptr {
        fits_write_img(
            fptr_box,
            TUSHORT,
            fpixel as LONGLONG,
            nelements as LONGLONG,
            cast_slice::<c_ushort, u8>(&array),
            &mut status,
        );
        assert_eq!(status, 0);
    }

    /* write another optional keyword to the header */
    /* Note that the ADDRESS of the value is passed in the routine */
    let exposure: c_long = 1500;
    if let Some(ref mut fptr_box) = fptr {
        let exposure_key = CString::new("EXPOSURE").unwrap();
        let exposure_comment = CString::new("Total Exposure Time").unwrap();

        fits_update_key(
            fptr_box,
            KeywordDatatype::TLONG(&exposure),
            cast_slice(exposure_key.to_bytes_with_nul()),
            Some(cast_slice(exposure_comment.to_bytes_with_nul())),
            &mut status,
        );
        assert_eq!(status, 0);
    }

    if let Some(fptr_box) = fptr {
        unsafe { fits_close_file(fptr_box, &mut status) };
        assert_eq!(status, 0);
    }
}

fn write_binary_tbl(filename: &str) {
    let mut fptr: Option<Box<fitsfile>> = None;
    let mut status: c_int = 0;
    let mut hdutype: c_int = 0;
    let firstrow: c_long = 1;
    let firstelem: c_long = 1;

    let tfields: c_int = 3; /* table will have 3 columns */
    let nrows: LONGLONG = 6; /* table will have 6 rows    */

    let extname = "PLANETS_Binary"; /* extension name */

    /* define the name, datatype, and physical units for the 3 columns */
    let ttype_strs = [
        CString::new("Planet").unwrap(),
        CString::new("Diameter").unwrap(),
        CString::new("Density").unwrap(),
    ];
    let tform_strs = [
        CString::new("8a").unwrap(),
        CString::new("1J").unwrap(),
        CString::new("1E").unwrap(),
    ];
    let tunit_strs = [
        CString::new("").unwrap(),
        CString::new("km").unwrap(),
        CString::new("g/cm").unwrap(),
    ];

    let ttype: Vec<&[c_char]> = ttype_strs
        .iter()
        .map(|s| cast_slice(s.to_bytes_with_nul()))
        .collect();
    let tform: Vec<&[c_char]> = tform_strs
        .iter()
        .map(|s| cast_slice(s.to_bytes_with_nul()))
        .collect();
    let tunit: Vec<&[c_char]> = tunit_strs
        .iter()
        .map(|s| cast_slice(s.to_bytes_with_nul()))
        .collect();

    /* define the name diameter, and density of each planet */
    let planet_strs = [
        CString::new("Mercury").unwrap(),
        CString::new("Venus").unwrap(),
        CString::new("Earth").unwrap(),
        CString::new("Mars").unwrap(),
        CString::new("Jupiter").unwrap(),
        CString::new("Saturn").unwrap(),
    ];
    let planet: Vec<&[c_char]> = planet_strs
        .iter()
        .map(|s| cast_slice(s.to_bytes_with_nul()))
        .collect();
    let diameter: [c_long; 6] = [4880, 12112, 12742, 6800, 143000, 121000];
    let density: [c_float; 6] = [5.1, 5.3, 5.52, 3.94, 1.33, 0.69];

    let filename_cstr = CString::new(filename).unwrap();
    let extname_cstr = CString::new(extname).unwrap();

    /* open the FITS file containing a primary array and an ASCII table */
    fits_open_file(
        &mut fptr,
        cast_slice(filename_cstr.to_bytes_with_nul()),
        READWRITE,
        &mut status,
    );
    assert_eq!(status, 0);

    if let Some(ref mut fptr_box) = fptr {
        fits_movabs_hdu(fptr_box, 1, Some(&mut hdutype), &mut status);
        assert_eq!(status, 0);
    }

    /* append a new empty binary table onto the FITS file */
    if let Some(ref mut fptr_box) = fptr {
        let ttype_opts: Vec<Option<&[c_char]>> = ttype.iter().map(|s| Some(*s)).collect();
        let tunit_opts: Vec<Option<&[c_char]>> = tunit.iter().map(|s| Some(*s)).collect();

        unsafe {
            fits_create_tbl(
                fptr_box,
                BINARY_TBL,
                nrows,
                tfields,
                &ttype_opts,
                &tform,
                Some(&tunit_opts),
                Some(cast_slice(extname_cstr.to_bytes_with_nul())),
                &mut status,
            );
        }
        assert_eq!(status, 0);
    }

    /* write names to the first column (character strings) */
    /* write diameters to the second column (longs) */
    /* write density to the third column (floats) */

    if let Some(ref mut fptr_box) = fptr {
        fits_write_col_str(
            fptr_box,
            1,
            firstrow as LONGLONG,
            firstelem as LONGLONG,
            planet.len() as LONGLONG,
            &planet,
            &mut status,
        );
        assert_eq!(status, 0);

        fits_write_col(
            fptr_box,
            TLONG,
            2,
            firstrow as LONGLONG,
            firstelem as LONGLONG,
            diameter.len() as LONGLONG,
            cast_slice(&diameter),
            &mut status,
        );
        assert_eq!(status, 0);

        fits_write_col(
            fptr_box,
            TFLOAT,
            3,
            firstrow as LONGLONG,
            firstelem as LONGLONG,
            density.len() as LONGLONG,
            cast_slice(&density),
            &mut status,
        );
        assert_eq!(status, 0);
    }

    if let Some(fptr_box) = fptr {
        unsafe { fits_close_file(fptr_box, &mut status) };
        assert_eq!(status, 0);
    }
}

fn readtable(filename: &str) {
    let mut fptr: Option<Box<fitsfile>> = None;
    let mut status: c_int = 0;
    let mut hdutype: c_int = 0;
    let mut nfound: c_int = 0;
    let frow: c_long = 1;
    let felem: c_long = 1;
    let nelem: c_long = 6;
    let longnull: c_long = 0;
    let mut dia: [c_long; 6] = [0; 6];
    let floatnull: c_float = 0.0;
    let mut den: [c_float; 6] = [0.0; 6];

    let filename_cstr = CString::new(filename).unwrap();

    fits_open_file(
        &mut fptr,
        cast_slice(filename_cstr.to_bytes_with_nul()),
        READONLY,
        &mut status,
    );
    assert_eq!(status, 0);

    /* allocate space for the column labels - use Vec instead of malloc */
    let mut ttype_vecs: Vec<Vec<u8>> = Vec::new();
    for _i in 0..3 {
        ttype_vecs.push(vec![0; FLEN_VALUE]);
    }

    /* allocate space for string column values - use Vec instead of malloc */
    let mut name_vecs: Vec<Vec<u8>> = Vec::new();
    for _i in 0..6 {
        name_vecs.push(vec![0; 10]);
    }

    for hdunum in 2..=2 {
        /*read ASCII, then binary table */
        /* move to the HDU */
        if let Some(ref mut fptr_box) = fptr {
            fits_movabs_hdu(fptr_box, hdunum, Some(&mut hdutype), &mut status);
            assert_eq!(status, 0);
        }

        if hdutype == ASCII_TBL {
            println!("\nReading ASCII table in HDU {hdunum}:");
        } else if hdutype == BINARY_TBL {
            println!("\nReading binary table in HDU {hdunum}:");
        } else {
            println!("Error: this HDU is not an ASCII or binary table");
            assert_eq!(1, 0);
        }

        /* read the column names from the TTYPEn keywords */
        if let Some(ref mut fptr_box) = fptr {
            let ttype_key = CString::new("TTYPE").unwrap();
            let mut ttype_refs: Vec<&mut [u8]> =
                ttype_vecs.iter_mut().map(|v| v.as_mut_slice()).collect();
            fits_read_keys_str(
                fptr_box,
                cast_slice(ttype_key.to_bytes_with_nul()),
                1,
                3,
                &mut ttype_refs
                    .iter_mut()
                    .map(|r| cast_slice_mut(r))
                    .collect::<Vec<_>>(),
                &mut nfound,
                &mut status,
            );
            assert_eq!(status, 0);
        }

        let ttype0 =
            unsafe { CStr::from_ptr(ttype_vecs[0].as_ptr() as *const c_char) }.to_string_lossy();
        let ttype1 =
            unsafe { CStr::from_ptr(ttype_vecs[1].as_ptr() as *const c_char) }.to_string_lossy();
        let ttype2 =
            unsafe { CStr::from_ptr(ttype_vecs[2].as_ptr() as *const c_char) }.to_string_lossy();
        println!(" Row  {ttype0:>10} {ttype1:>10} {ttype2:>10}");

        /*  read the columns */
        if let Some(ref mut fptr_box) = fptr {
            let mut name_refs: Vec<&mut [u8]> =
                name_vecs.iter_mut().map(|v| v.as_mut_slice()).collect();

            fits_read_col_str(
                fptr_box,
                1,
                frow as LONGLONG,
                felem as LONGLONG,
                nelem as LONGLONG,
                None,
                &mut name_refs
                    .iter_mut()
                    .map(|r| cast_slice_mut(r))
                    .collect::<Vec<_>>(),
                None,
                &mut status,
            );
            assert_eq!(status, 0);

            fits_read_col(
                fptr_box,
                TLONG,
                2,
                frow as LONGLONG,
                felem as LONGLONG,
                nelem as LONGLONG,
                Some(NullValue::Long(longnull)),
                cast_slice_mut(&mut dia),
                None,
                &mut status,
            );
            assert_eq!(status, 0);

            fits_read_col(
                fptr_box,
                TFLOAT,
                3,
                frow as LONGLONG,
                felem as LONGLONG,
                nelem as LONGLONG,
                Some(NullValue::Float(floatnull)),
                cast_slice_mut(&mut den),
                None,
                &mut status,
            );
            assert_eq!(status, 0);

            for ii in 0..6 {
                let name_str = unsafe { CStr::from_ptr(name_vecs[ii].as_ptr() as *const c_char) }
                    .to_string_lossy();
                println!(
                    "{:5} {:>10} {:>10} {:>10.2}",
                    ii + 1,
                    name_str,
                    dia[ii],
                    den[ii]
                );
            }

            let where_clause = cast_slice(c"D > 3".to_bytes_with_nul());
            let mut n_matched_rows = -1;
            let mut row_status = [0; 6];

            // Try to query an invalid row.
            fits_find_rows(
                fptr_box,
                where_clause,
                1,
                6,
                &mut n_matched_rows,
                &mut row_status,
                &mut status,
            );
            assert_eq!(status, KEY_NO_EXIST);
            assert_eq!(n_matched_rows, -1);

            // Try a valid query: Planet = Earth
            status = 0;
            let where_clause = cast_slice(c"Planet == \"Earth\"".to_bytes_with_nul());
            fits_find_rows(
                fptr_box,
                where_clause,
                1,
                6,
                &mut n_matched_rows,
                &mut row_status,
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(n_matched_rows, 1);

            let where_clause =
                cast_slice(c"(DENSITY > 3.0) && (DIAMETER > 10000)".to_bytes_with_nul());
            let mut n_matched_rows = -1;
            let mut row_status = [0; 6];
            fits_find_rows(
                fptr_box,
                where_clause,
                1,
                6,
                &mut n_matched_rows,
                &mut row_status,
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(n_matched_rows, 2);

            let where_clause = cast_slice(
                c"(DENSITY > 3.0) && (DIAMETER > 10000) && (Planet == \"Earth\")"
                    .to_bytes_with_nul(),
            );
            let mut n_matched_rows = -1;
            let mut row_status = [0; 6];
            fits_find_rows(
                fptr_box,
                where_clause,
                1,
                6,
                &mut n_matched_rows,
                &mut row_status,
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(n_matched_rows, 1);

            let where_clause = cast_slice(c"#ROW > 2".to_bytes_with_nul());
            let mut n_matched_rows = -1;
            let mut row_status = [0; 6];
            fits_find_rows(
                fptr_box,
                where_clause,
                1,
                6,
                &mut n_matched_rows,
                &mut row_status,
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(n_matched_rows, 4);
        }
    }

    if let Some(fptr_box) = fptr {
        unsafe {
            fits_close_file(fptr_box, &mut status);
        }
        assert_eq!(status, 0);
    }
}

#[test]
fn test_read_table_where() {
    with_temp_file(|filename| {
        writeimage(filename);
        write_binary_tbl(filename);
        readtable(filename);
    });
}
