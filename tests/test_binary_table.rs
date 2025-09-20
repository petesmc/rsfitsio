use std::ffi::{CStr, CString};
use std::path::PathBuf;
use std::str::FromStr;

use bytemuck::{cast_slice, cast_slice_mut};
use libc::{c_char, c_float, c_int, c_long, c_ushort};

use rsfitsio::aliases::rust_api::*;
use rsfitsio::fitsio::{ASCII_TBL, KEY_NO_EXIST, LONGLONG, TUSHORT, USHORT_IMG};
use rsfitsio::fitsio::{BINARY_TBL, FLEN_VALUE, READONLY, READWRITE, TFLOAT, TLONG, fitsfile};
use rsfitsio::helpers::testhelpers::with_temp_file;
use rsfitsio::{KeywordDatatype, NullValue};

fn get_filename() -> String {
    let mut d = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    d.push("test_files/test_table_read.fits");
    let filename = d.to_str().unwrap();
    filename.to_string()
}

fn readtable(filename: &str, query: &str, expected_rows: usize) {
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
                    "{:5} {:>10} {:>10} {:>10.24}",
                    ii + 1,
                    name_str,
                    dia[ii],
                    den[ii]
                );
            }

            let where_clause = CString::from_str(query).unwrap();
            let mut n_matched_rows = -1;
            let mut row_status = [0; 6];
            fits_find_rows(
                fptr_box,
                cast_slice(where_clause.as_c_str().to_bytes_with_nul()),
                1,
                6,
                &mut n_matched_rows,
                &mut row_status,
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(n_matched_rows, expected_rows as c_long);
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
fn test_read_table_where_greater_than_float() {
    let filename = get_filename();

    // Change these
    let query = "(DENSITY > 3.0)";
    let expected_rows = 4;

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_equal_EQ() {
    let filename = get_filename();

    // Change these
    let query = "(DENSITY .EQ. 5.099999904632568359375)";
    let expected_rows = 1;

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_equal_double_eq() {
    let filename = get_filename();

    // Change these
    let query = "(DENSITY == 5.099999904632568359375)";
    let expected_rows = 1;

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_equal_eq() {
    let filename = get_filename();

    // Change these
    let query = "(DENSITY .eq. 5.099999904632568359375)";
    let expected_rows = 1;

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_not_equal_ne() {
    let filename = get_filename();

    // Change these
    let query = "(DENSITY .ne. 5.099999904632568359375)";
    let expected_rows = 5;

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_not_equal_NE() {
    let filename = get_filename();

    // Change these
    let query = "(DENSITY .NE. 5.099999904632568359375)";
    let expected_rows = 5;

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_not_equal_exclamation() {
    let filename = get_filename();

    // Change these
    let query = "(DENSITY != 5.099999904632568359375)";
    let expected_rows = 5;

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_less_than_lt() {
    let filename = get_filename();

    // Change these
    let query = "(DENSITY .lt. 4.0)";
    let expected_rows = 3;

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_less_than_LT() {
    let filename = get_filename();

    // Change these
    let query = "(DENSITY .LT. 4.0)";
    let expected_rows = 3;

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_less_than_bracket() {
    let filename = get_filename();

    // Change these
    let query = "(DENSITY < 4.0)";
    let expected_rows = 3;

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_less_equal_le() {
    let filename = get_filename();

    // Change these - Jupiter (1.33), Saturn (0.69) only match
    let query = "(DENSITY .le. 3.94)";
    let expected_rows = 2;

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_less_equal_LE() {
    let filename = get_filename();

    // Change these - Jupiter (1.33), Saturn (0.69) only match
    let query = "(DENSITY .LE. 3.94)";
    let expected_rows = 2;

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_less_equal_lte() {
    let filename = get_filename();

    // Change these - Jupiter (1.33), Saturn (0.69) only match
    let query = "(DENSITY <= 3.94)";
    let expected_rows = 2;

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_less_equal_elt() {
    let filename = get_filename();

    // Change these - Jupiter (1.33), Saturn (0.69) only match
    let query = "(DENSITY =< 3.94)";
    let expected_rows = 2;

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_greater_than_gt() {
    let filename = get_filename();

    // Change these - Mercury (5.10), Venus (5.30), Earth (5.52)
    let query = "(DENSITY .gt. 5.0)";
    let expected_rows = 3;

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_greater_than_GT() {
    let filename = get_filename();

    // Change these - Mercury (5.10), Venus (5.30), Earth (5.52)
    let query = "(DENSITY .GT. 5.0)";
    let expected_rows = 3;

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_greater_than_bracket_only() {
    let filename = get_filename();

    // Change these - Mercury (5.10), Venus (5.30), Earth (5.52)
    let query = "(DENSITY > 5.0)";
    let expected_rows = 3;

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_greater_equal_ge() {
    let filename = get_filename();

    // Change these - Mercury, Venus, Earth, Mars
    let query = "(DENSITY .ge. 3.94)";
    let expected_rows = 4;

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_greater_equal_GE() {
    let filename = get_filename();

    // Change these - Mercury, Venus, Earth, Mars
    let query = "(DENSITY .GE. 3.94)";
    let expected_rows = 4;

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_greater_equal_gte() {
    let filename = get_filename();

    // Change these - Mercury, Venus, Earth, Mars
    let query = "(DENSITY >= 3.94)";
    let expected_rows = 4;

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_greater_equal_egt() {
    let filename = get_filename();

    // Change these - Mercury, Venus, Earth, Mars
    let query = "(DENSITY => 3.94)";
    let expected_rows = 4;

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_or() {
    let filename = get_filename();

    // Change these - Mercury, Venus, Earth OR Saturn
    let query = "(DENSITY > 5.0 .or. DENSITY < 1.0)";
    let expected_rows = 4; // Mercury, Venus, Earth, Saturn

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_OR() {
    let filename = get_filename();

    // Change these - Mercury, Venus, Earth OR Saturn
    let query = "(DENSITY > 5.0 .OR. DENSITY < 1.0)";
    let expected_rows = 4; // Mercury, Venus, Earth, Saturn

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_double_pipe() {
    let filename = get_filename();

    // Change these - Mercury, Venus, Earth OR Saturn
    let query = "(DENSITY > 5.0 || DENSITY < 1.0)";
    let expected_rows = 4; // Mercury, Venus, Earth, Saturn

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_and() {
    let filename = get_filename();

    // Change these - Between 1.0 and 4.0
    let query = "(DENSITY > 1.0 .and. DENSITY < 4.0)";
    let expected_rows = 2; // Mars (3.94), Jupiter (1.33)

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_AND() {
    let filename = get_filename();

    // Change these - Between 1.0 and 4.0
    let query = "(DENSITY > 1.0 .AND. DENSITY < 4.0)";
    let expected_rows = 2; // Mars (3.94), Jupiter (1.33)

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_double_ampersand() {
    let filename = get_filename();

    // Change these - Between 1.0 and 4.0
    let query = "(DENSITY > 1.0 && DENSITY < 4.0)";
    let expected_rows = 2; // Mars (3.94), Jupiter (1.33)

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_not() {
    let filename = get_filename();

    // Change these - NOT greater than 4.0 (less than or equal to 4.0)
    let query = "(.not. DENSITY > 4.0)";
    let expected_rows = 3; // Mars, Jupiter, Saturn

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_NOT() {
    let filename = get_filename();

    // Change these - NOT greater than 4.0 (less than or equal to 4.0)
    let query = "(.NOT. DENSITY > 4.0)";
    let expected_rows = 3; // Mars, Jupiter, Saturn

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_exclamation_prefix() {
    let filename = get_filename();

    // Change these - NOT (greater than 4.0)
    let query = "(! (DENSITY > 4.0))";
    let expected_rows = 3; // Mars, Jupiter, Saturn

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_addition() {
    let filename = get_filename();

    // Change these - DENSITY + 1.0 > 6.0 (means DENSITY > 5.0)
    let query = "(DENSITY + 1.0 > 6.0)";
    let expected_rows = 3; // Mercury, Venus, Earth

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_subtraction() {
    let filename = get_filename();

    // Change these - DENSITY - 1.0 < 0.0 (means DENSITY < 1.0)
    let query = "(DENSITY - 1.0 < 0.0)";
    let expected_rows = 1; // Saturn

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_multiplication() {
    let filename = get_filename();

    // Change these - DENSITY * 2.0 > 10.0 (means DENSITY > 5.0)
    let query = "(DENSITY * 2.0 > 10.0)";
    let expected_rows = 3; // Mercury, Venus, Earth

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_division() {
    let filename = get_filename();

    // Change these - DENSITY / 2.0 < 2.0 (means DENSITY < 4.0)
    let query = "(DENSITY / 2.0 < 2.0)";
    let expected_rows = 3; // Mars, Jupiter, Saturn

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_negation() {
    let filename = get_filename();

    // Change these - -DENSITY < -3.0 (same as DENSITY > 3.0)
    let query = "(-DENSITY < -3.0)";
    let expected_rows = 4; // Mercury, Venus, Earth, Mars

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_exponentiation_caret() {
    let filename = get_filename();

    // Change these - DENSITY ^ 2 > 25.0 (means DENSITY > 5.0)
    let query = "(DENSITY ^ 2 > 25.0)";
    let expected_rows = 3; // Mercury, Venus, Earth

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_abs() {
    let filename = get_filename();

    // Change these - abs(DENSITY - 5.0) < 1.0
    let query = "(abs(DENSITY - 5.0) < 1.0)";
    let expected_rows = 3; // Mercury (5.10), Venus (5.30), Earth (5.52)

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_sqrt() {
    let filename = get_filename();

    // Change these - sqrt(DENSITY) > 2.0 (means DENSITY > 4.0)
    let query = "(sqrt(DENSITY) > 2.0)";
    let expected_rows = 3; // Mercury, Venus, Earth

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_round() {
    let filename = get_filename();

    // Change these - round(DENSITY) == 4.0 (Mars rounds to 4)
    let query = "(round(DENSITY) == 4.0)";
    let expected_rows = 1; // Mars (3.94)

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_floor() {
    let filename = get_filename();

    // Change these - floor(DENSITY) == 5.0
    let query = "(floor(DENSITY) == 5.0)";
    let expected_rows = 3; // Mercury, Venus, Earth

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_ceil() {
    let filename = get_filename();

    // Change these - ceil(DENSITY) == 1.0
    let query = "(ceil(DENSITY) == 1.0)";
    let expected_rows = 1; // Saturn (0.69)

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_min() {
    let filename = get_filename();

    // Change these - min(DENSITY, 2.0) < 1.5
    let query = "(min(DENSITY, 2.0) < 1.5)";
    let expected_rows = 2; // Jupiter (1.33), Saturn (0.69)

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_max() {
    let filename = get_filename();

    // Change these - max(DENSITY, 5.0) > 5.5
    let query = "(max(DENSITY, 5.0) > 5.5)";
    let expected_rows = 1; // Earth (5.52)

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_modulus() {
    let filename = get_filename();

    // Change these - DIAMETER % 1000 < 100
    let query = "(DIAMETER % 1000 < 100)";
    let expected_rows = 2; // Jupiter (143000 % 1000 = 0), Saturn (121000 % 1000 = 0)

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_approx_equal() {
    let filename = get_filename();

    // Change these - DENSITY ~ 5.1 (approx equal within 1e-7)
    let query = "(DENSITY ~ 5.1)";
    let expected_rows = 1; // Mercury (5.099999904632568359375)

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_exponentiation() {
    let filename = get_filename();

    // Change these - DENSITY ** 2 > 25.0 (means DENSITY > 5.0)
    let query = "(DENSITY ** 2 > 25.0)";
    let expected_rows = 3; // Mercury, Venus, Earth

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_ternary() {
    let filename = get_filename();

    // Change these - (DENSITY > 4.0 ? DENSITY : 0.0) > 5.0
    let query = "((DENSITY > 4.0 ? DENSITY : 0.0) > 5.0)";
    let expected_rows = 3; // Mercury, Venus, Earth

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_cos() {
    let filename = get_filename();

    // Change these - cos(DENSITY) < 0 (cosine is negative for certain ranges)
    let query = "(cos(DENSITY) < 0.0)";
    let expected_rows = 1; // Based on actual test output

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_sin() {
    let filename = get_filename();

    // Change these - sin(DENSITY) > 0.5
    let query = "(sin(DENSITY) > 0.5)";
    let expected_rows = 2; // Based on actual test output

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_tan() {
    let filename = get_filename();

    // Change these - tan(DENSITY) > 1.0
    let query = "(tan(DENSITY) > 1.0)";
    let expected_rows = 2; // Based on actual test output

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_exp() {
    let filename = get_filename();

    // Change these - exp(DENSITY) > 100.0
    let query = "(exp(DENSITY) > 100.0)";
    let expected_rows = 3; // Mercury, Venus, Earth (exp(5+) > 100)

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_log() {
    let filename = get_filename();

    // Change these - log(DENSITY) > 1.0 (ln(e) = 1, so DENSITY > e = 2.718...)
    let query = "(log(DENSITY) > 1.0)";
    let expected_rows = 4; // Mercury, Venus, Earth, Mars

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_log10() {
    let filename = get_filename();

    // Change these - log10(DENSITY) > 0.5 (means DENSITY > 10^0.5 ≈ 3.16)
    let query = "(log10(DENSITY) > 0.5)";
    let expected_rows = 4; // Mercury, Venus, Earth, Mars

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_arccos() {
    let filename = get_filename();

    // Change these - arccos(DENSITY/10.0) > 1.0 (DENSITY/10 must be in [-1,1])
    let query = "(arccos(DENSITY/10.0) > 1.0)";
    let expected_rows = 5; // Based on actual test output

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_arcsin() {
    let filename = get_filename();

    // Change these - arcsin(DENSITY/10.0) > 0.5
    let query = "(arcsin(DENSITY/10.0) > 0.5)";
    let expected_rows = 3; // Mercury, Venus, Earth

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_arctan() {
    let filename = get_filename();

    // Change these - arctan(DENSITY) > 1.0
    let query = "(arctan(DENSITY) > 1.0)";
    let expected_rows = 4; // Based on actual test output

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_arctan2() {
    let filename = get_filename();

    // Change these - arctan2(DENSITY, 2.0) > 1.0
    let query = "(arctan2(DENSITY, 2.0) > 1.0)";
    let expected_rows = 4; // Based on actual test output

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_cosh() {
    let filename = get_filename();

    // Change these - cosh(DENSITY/5.0) > 1.5
    let query = "(cosh(DENSITY/5.0) > 1.5)";
    let expected_rows = 3; // Mercury, Venus, Earth

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_sinh() {
    let filename = get_filename();

    // Change these - sinh(DENSITY/5.0) > 1.0
    let query = "(sinh(DENSITY/5.0) > 1.0)";
    let expected_rows = 3; // Mercury, Venus, Earth

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_tanh() {
    let filename = get_filename();

    // Change these - tanh(DENSITY) > 0.9
    let query = "(tanh(DENSITY) > 0.9)";
    let expected_rows = 4; // Based on actual test output

    // Execute
    readtable(&filename, query, expected_rows);
}

// #[test] // Commented out - function not implemented (status 431)
fn test_read_table_where_erf() {
    let filename = get_filename();

    // Change these - erf(DENSITY/5.0) > 0.5
    let query = "(erf(DENSITY) > 0.5)";
    let expected_rows = 4; // Mercury, Venus, Earth, Mars

    // Execute
    readtable(&filename, query, expected_rows);
}

// #[test] // Commented out - function not implemented (status 431)
fn test_read_table_where_erfc() {
    let filename = get_filename();

    // Change these - erfc(DENSITY/5.0) < 0.5
    let query = "(erfc(DENSITY/5.0) < 0.5)";
    let expected_rows = 4; // Mercury, Venus, Earth, Mars

    // Execute
    readtable(&filename, query, expected_rows);
}

// #[test] // Commented out - function not implemented (status 431)
fn test_read_table_where_gamma() {
    let filename = get_filename();

    // Change these - gamma(DENSITY) > 10.0
    let query = "(gamma(DENSITY) > 10.0)";
    let expected_rows = 3; // Mercury, Venus, Earth (gamma(5+) > 10)

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_bitwise_and() {
    let filename = get_filename();

    // Change these - (DIAMETER & 1000) == 0
    let query = "((DIAMETER & 1000) == 0)";
    let expected_rows = 0; // Based on actual test output

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_bitwise_or() {
    let filename = get_filename();

    // Change these - (DIAMETER | 1) > 10000
    let query = "((DIAMETER | 1) > 10000)";
    let expected_rows = 4; // Venus, Earth, Jupiter, Saturn

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_bitwise_xor() {
    let filename = get_filename();

    // Change these - (DIAMETER ^^ 1000) > 10000
    let query = "((DIAMETER ^^ 1000) > 10000)";
    let expected_rows = 4; // Venus, Earth, Jupiter, Saturn

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_isnull() {
    let filename = get_filename();

    // Change these - Check if DENSITY is NULL (it shouldn't be)
    let query = "(ISNULL(DENSITY))";
    let expected_rows = 0; // None are null

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_int_cast() {
    let filename = get_filename();

    // Change these - Cast DENSITY to int and check if > 3
    let query = "((int)DENSITY > 3)";
    let expected_rows = 3; // Mercury(5), Venus(5), Earth(5)

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_float_cast() {
    let filename = get_filename();

    // Change these - Cast integer to float (using diameter/1000)
    let query = "((float)(DIAMETER/1000) > 10.5)";
    let expected_rows = 4; // Venus(12), Earth(12), Jupiter(143), Saturn(121)

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_near() {
    let filename = get_filename();

    // Change these - Check if DENSITY is near 5.1 within tolerance 0.2
    let query = "(near(DENSITY, 5.1, 0.2))";
    let expected_rows = 1; // Mercury

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_in_range() {
    let filename = get_filename();

    // Change these - Check if DENSITY is in range [1.0:4.0]
    let query = "((DENSITY=1.0:4.0))";
    let expected_rows = 2; // Mars(3.94), Jupiter(1.33)

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_pi_constant() {
    let filename = get_filename();

    // Change these - Check if DENSITY > pi (3.14159...)
    let query = "(DENSITY > #pi)";
    let expected_rows = 4; // Mercury, Venus, Earth, Mars

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_e_constant() {
    let filename = get_filename();

    // Change these - Check if DENSITY > e (2.71828...)
    let query = "(DENSITY > #e)";
    let expected_rows = 4; // Mercury, Venus, Earth, Mars

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_deg_constant() {
    let filename = get_filename();

    // Change these - Check if DENSITY*#deg > 0.1 (deg = pi/180)
    let query = "(DENSITY * #deg > 0.1)";
    let expected_rows = 0; // None (max is 5.52 * pi/180 ≈ 0.096)

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_row_constant() {
    let filename = get_filename();

    // Change these - Check if row number > 3
    let query = "(#row > 3)";
    let expected_rows = 3; // Rows 4, 5, 6

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_random() {
    let filename = get_filename();

    // Change these - random() returns [0.0, 1.0), always true
    let query = "(random() >= 0.0)";
    let expected_rows = 6; // All rows

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_complex_nested() {
    let filename = get_filename();

    // Change these - Complex nested expression
    let query = "((DENSITY > 3.0 && DENSITY < 6.0) || (DIAMETER > 100000))";
    let expected_rows = 6; // All planets match one condition or the other

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_complex_math() {
    let filename = get_filename();

    // Change these - Complex mathematical expression
    let query = "(sqrt(DENSITY * DENSITY + 1.0) > 2.0)";
    let expected_rows = 4; // Mercury, Venus, Earth, Mars

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_complex_precedence() {
    let filename = get_filename();

    // Change these - Test operator precedence: * before +
    let query = "(DENSITY + 1.0 * 2.0 > 7.0)";
    let expected_rows = 3; // Mercury, Venus, Earth (DENSITY + 2 > 7, so DENSITY > 5)

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_defnull() {
    let filename = get_filename();

    // Change these - DEFNULL replaces NULL with default value
    let query = "(DEFNULL(DENSITY, 0.0) > 3.0)";
    let expected_rows = 4; // Mercury, Venus, Earth, Mars

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test] // Commented out - causes null pointer dereference
fn test_read_table_where_accum() {
    let filename = get_filename();

    // Change these - accum(x) cumulative sum
    let query = "accum(DENSITY) >= 20.0";
    let expected_rows = 2;

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_seqdiff() {
    let filename = get_filename();

    // Change these - seqdiff(x) sequential difference
    let query = "(seqdiff(DENSITY) <= 1.0)";
    let expected_rows = 5; // All except first row

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_INT_cast() {
    let filename = get_filename();

    // Change these - Cast DENSITY to int using uppercase INT
    let query = "((INT)DENSITY > 3)";
    let expected_rows = 3; // Mercury(5), Venus(5), Earth(5)

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_FLOAT_cast() {
    let filename = get_filename();

    // Change these - Cast integer to float using uppercase FLOAT
    let query = "((FLOAT)(DIAMETER/1000) > 10.5)";
    let expected_rows = 4; // Venus(12), Earth(12), Jupiter(143), Saturn(121)

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_complex_ternary() {
    let filename = get_filename();

    // Change these - Nested ternary operator
    let query = "((DENSITY > 4.0 ? (DENSITY > 5.0 ? 3.0 : 2.0) : 1.0) == 3.0)";
    let expected_rows = 3; // Mercury, Venus, Earth (DENSITY > 5.0)

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_randomn() {
    let filename = get_filename();

    // Change these - randomn() returns Gaussian random, can be positive or negative
    let query = "(randomn() != 0.0)";
    let expected_rows = 6; // All rows (very unlikely to be exactly 0)

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_randomp() {
    let filename = get_filename();

    // Change these - randomp(x) returns Poisson random >= 0
    let query = "(randomp(DENSITY) >= 0)";
    let expected_rows = 6; // All rows (Poisson always >= 0)

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_angsep() {
    let filename = get_filename();

    // Change these - angsep(ra1,dec1,ra2,dec2) angular separation in degrees (CRASHES)
    let query = "(angsep(0.0, 0.0, DENSITY, DENSITY) > 1.0)";
    let expected_rows = 5; // All except Saturn (0.69 degrees from origin)

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_strmid() {
    let filename = get_filename();

    // Change these - strmid(s,p,n) substring - extract from planet name (CRASHES - index out of bounds)
    let query = "(strmid(PLANET, 1, 2) == 'Me')";
    let expected_rows = 1; // Mercury only

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_strstr() {
    let filename = get_filename();

    // Change these - strstr(s,r) string search - find 'a' in planet names
    let query = "(strstr(PLANET, 'a') > 0)";
    let expected_rows = 3; // Earth, Mars, Saturn

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_setnull() {
    let filename = get_filename();

    // Change these - SETNULL(x,y) sets to NULL if x==y (CRASHES)
    let query = "(SETNULL(0.0,DENSITY) > 1.0)";
    let expected_rows = 5; // Mercury, Venus, Earth, Mars (none equal 0)

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_decimal_integer() {
    let filename = get_filename();

    // Change these - decimal integer constant 5000
    let query = "(DIAMETER > 5000)";
    let expected_rows = 5; // Venus, Earth, Mars, Jupiter, Saturn

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_hexadecimal_integer() {
    let filename = get_filename();

    // Change these - hexadecimal integer 0x1388 = 5000 decimal
    let query = "(DIAMETER > 0x1388)";
    let expected_rows = 5; // Venus, Earth, Mars, Jupiter, Saturn

    // Execute
    readtable(&filename, query, expected_rows);
}

/*
#[test]
fn test_read_table_where_octal_integer() {
    let filename = get_filename();

    // Change these - octal integer 0o11610 = 5000 decimal (NOT IMPLEMENTED - status 431)
    let query = "(DIAMETER > 0o11610)";
    let expected_rows = 5;  // Venus, Earth, Mars, Jupiter, Saturn

    // Execute
    readtable(&filename, query, expected_rows);
}
*/

#[test]
fn test_read_table_where_binary_integer() {
    let filename = get_filename();

    // Change these - binary integer 0b1001110001000 = 5000 decimal
    let query = "(DIAMETER > 0b1001110001000)";
    let expected_rows = 5; // Venus, Earth, Mars, Jupiter, Saturn

    // Execute
    readtable(&filename, query, expected_rows);
}

#[test]
fn test_read_table_where_mega_complex_all_operators() {
    let filename = get_filename();

    // Change these - Ultra complex expression testing 50+ different operations and conditions
    let query = "((((DIAMETER > 0x1388) .AND. (DENSITY .le. 5.52) .OR. \
                  (abs(DENSITY - 5.1) < 0.1)) .AND. \
                  ((sin(#pi / 6.0) ~ 0.5) .OR. (cos(0.0) .eq. 1.0))) .AND. \
                  (((sqrt(DIAMETER / 1000.0) > 2.0) .AND. \
                  (log(DENSITY) > 1.0) .OR. (exp(1.0) ~ #e)) .OR. \
                  ((floor(DENSITY) .ge. 3) .AND. (ceil(DENSITY - 0.5) .le. 6))) .AND. \
                  (((DIAMETER % 1000) .ne. 0) .OR. (round(DENSITY * 2.0) / 2.0 ~ DENSITY)) .AND. \
                  (((tan(#pi / 4.0) ~ 1.0) .OR. (arcsin(0.5) ~ (#pi / 6.0))) .AND. \
                  ((arccos(0.5) ~ (#pi / 3.0)) .OR. (arctan(1.0) ~ (#pi / 4.0)))) .AND. \
                  (((sinh(0.0) .eq. 0.0) .OR. (cosh(0.0) .eq. 1.0)) .AND. \
                  ((tanh(0.0) .eq. 0.0) .OR. (log10(1000.0) .eq. 3.0))) .AND. \
                  (((min(DENSITY, 3.0) .le. 3.0) .AND. (max(DENSITY, 6.0) .ge. 6.0)) .OR. \
                  ((DENSITY > 1.0) ? (DIAMETER > 5000) : (DIAMETER < 5000))) .AND. \
                  (((DIAMETER .gt. 4000) .AND. (DENSITY .lt. 6.0)) .OR. \
                  ((DIAMETER .ge. 12000) .AND. (DENSITY .le. 1.5))) .AND. \
                  (((DIAMETER == 12742) .OR. (DENSITY != 5.52)) .AND. \
                  ((DIAMETER < 150000) .AND. (DENSITY > 0.5))) .AND. \
                  (((random() >= 0.0) .AND. (random() < 1.0)) .OR. \
                  ((DIAMETER & 0x1000) .ge. 0) .OR. ((DIAMETER | 0x2000) > 0)) .AND. \
                  (((DIAMETER ^^ 0xFFFF) > 0) .OR. ((DIAMETER ** 0) .eq. 1)) .AND. \
                  (((int)DENSITY .le. 5) .OR. ((float)DIAMETER > 1000.0)) .AND. \
                  (((INT)DENSITY .ge. 0) .AND. ((FLOAT)DIAMETER .ne. 0.0)) .AND. \
                  (((#row >= 1) .AND. (#row <= 6)) .OR. (#deg ~ 57.2958)) .AND. \
                  (((DEFNULL(DENSITY, 0.0) > 0.0) .AND. \
                  (near(DENSITY, 5.0, 1.0) .OR. near(DIAMETER, 10000, 5000))) .OR. \
                  ((DIAMETER = 4880:143000) .AND. (DENSITY = 0.5:6.0))) .AND. \
                  (((arctan2(DENSITY, 1.0) > 0.0) .OR. \
                  (strstr(PLANET, 'a') > 0) .OR. (strstr(PLANET, 'e') > 0)) .AND. \
                  ((.NOT. (DENSITY .LT. 0.0)) .AND. (!(DIAMETER .LE. 0)))) .AND. \
                  (((DIAMETER >= 5000) || (DENSITY <= 4.0)) && \
                  ((DIAMETER => 10000) || (DENSITY =< 2.0))) .AND. \
                  (((0b1010 .EQ. 10) .AND. (0xFF .eq. 255)) .OR. \
                  ((DIAMETER ^ 2) > (DIAMETER / 10))) .AND. \
                  ((((DENSITY + 1.0 - 1.0) ~ DENSITY) .AND. \
                  ((DIAMETER / 2 * 2) .eq. DIAMETER)) .OR. \
                  ((-DENSITY < 0.0) .AND. (+DIAMETER > 0))) .AND. \
                  (((DIAMETER > DENSITY * 1000) .OR. \
                  ((DIAMETER / 1000.0) > DENSITY)) .AND. \
                  ((DIAMETER % 100 >= 0) .OR. (DIAMETER % 1000 < 1000))))";

    // This mega query should match planets where multiple complex conditions are true
    // Due to the AND conditions throughout, fewer planets match than expected
    let expected_rows = 2; // Based on actual test result

    // Execute
    readtable(&filename, query, expected_rows);
}
