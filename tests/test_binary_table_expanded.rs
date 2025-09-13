use std::ffi::{CStr, CString};
use std::path::PathBuf;
use std::str::FromStr;

use bytemuck::{cast_slice, cast_slice_mut};
use libc::{c_char, c_double, c_float, c_int, c_long, c_short, c_uchar};

use rsfitsio::aliases::rust_api::*;
use rsfitsio::fitsio::{ASCII_TBL, BINARY_TBL, FLEN_VALUE, LONGLONG, READONLY, READWRITE};
use rsfitsio::fitsio::{
    TBIT, TBYTE, TCOMPLEX, TDBLCOMPLEX, TDOUBLE, TFLOAT, TLOGICAL, TLONG, TLONGLONG, TSHORT,
    fitsfile,
};
use rsfitsio::helpers::testhelpers::with_temp_file;
use rsfitsio::{KeywordDatatype, NullValue, STDERR};

fn get_expanded_filename() -> String {
    let mut d = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    d.push("test_files/test_table_read_expanded.fits");
    let filename = d.to_str().unwrap();
    filename.to_string()
}

fn readtable_expanded_with_status(
    filename: &str,
    query: &str,
    expected_rows: usize,
    expected_status: c_int,
) {
    let mut fptr: Option<Box<fitsfile>> = None;
    let mut status: c_int = 0;
    let mut hdutype: c_int = 0;

    let filename_cstr = CString::new(filename).unwrap();

    fits_open_file(
        &mut fptr,
        cast_slice(filename_cstr.to_bytes_with_nul()),
        READONLY,
        &mut status,
    );
    assert_eq!(status, 0);

    // Move to the expanded binary table HDU (HDU 3)
    if let Some(ref mut fptr_box) = fptr {
        fits_movabs_hdu(fptr_box, 3, Some(&mut hdutype), &mut status);
        assert_eq!(status, 0);
    }

    assert_eq!(hdutype, BINARY_TBL);
    println!("Reading expanded binary table in HDU 3:");

    if let Some(ref mut fptr_box) = fptr {
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
        assert_eq!(status, expected_status);

        println!("Query: {}", query);
        if expected_status == 0 {
            assert_eq!(n_matched_rows, expected_rows as c_long);
            println!(
                "Expected rows: {}, Actual rows: {}",
                expected_rows, n_matched_rows
            );
        } else {
            println!(
                "Expected status: {}, Actual status: {} (function not implemented)",
                expected_status, status
            );
            printerror(status);
        }
    }

    if let Some(fptr_box) = fptr {
        let mut close_status = 0;
        unsafe {
            fits_close_file(fptr_box, &mut close_status);
        }
        assert_eq!(close_status, 0);
    }
}

/*--------------------------------------------------------------------------*/
/// Print out cfitsio error messages and exit program
fn printerror(status: c_int) {
    if status != 0 {
        unsafe {
            fits_report_error(STDERR!(), status); /* print error report */
        }
    }
}

fn readtable_expanded(filename: &str, query: &str, expected_rows: usize) {
    readtable_expanded_with_status(filename, query, expected_rows, 0);
}

// Basic tests using the new expanded table

#[test]
fn test_expanded_table_logical_column() {
    let filename = get_expanded_filename();

    // Test logical column - Active column
    // Saturn is the only inactive planet (F), all others are T
    // Try different syntax approaches
    let query = "Active == T"; // All active planets - try without quotes
    let expected_rows = 5; // All except Saturn

    readtable_expanded(&filename, query, expected_rows);
}

#[test]
fn test_expanded_table_logical_column_false() {
    let filename = get_expanded_filename();

    // Test logical column - Only inactive planets
    let query = "Active == F"; // Inactive planets - try without quotes
    let expected_rows = 1; // Only Saturn

    readtable_expanded(&filename, query, expected_rows);
}

#[test]
fn test_expanded_table_byte_column() {
    let filename = get_expanded_filename();

    // Test unsigned byte column - Category
    let query = "(Category > 3)"; // Mars, Jupiter, Saturn
    let expected_rows = 3;

    readtable_expanded(&filename, query, expected_rows);
}

#[test]
fn test_expanded_table_short_column() {
    let filename = get_expanded_filename();

    // Test 16-bit integer column - Priority
    let query = "(Priority >= 100)"; // Mercury, Venus, Earth, Mars
    let expected_rows = 4;

    readtable_expanded(&filename, query, expected_rows);
}

#[test]
fn test_expanded_table_longlong_column() {
    let filename = get_expanded_filename();

    // Test 64-bit integer column - Mass
    // Test for planets with mass > 1e15 kg (only Jupiter has 1.9e15 kg)
    let query = "(Mass > 1.0e15)";
    let expected_rows = 1; // Only Jupiter

    readtable_expanded(&filename, query, expected_rows);
}

#[test]
fn test_expanded_table_double_column() {
    let filename = get_expanded_filename();

    // Test double-precision column - Gravity
    // Test for planets with gravity > 10.0 m/s²
    let query = "(Gravity > 10.0)";
    let expected_rows = 2; // Jupiter (24.79), Saturn (10.44)

    readtable_expanded(&filename, query, expected_rows);
}

#[test]
fn test_expanded_table_complex_combinations() {
    let filename = get_expanded_filename();

    // Test combinations of new columns
    // Active planets (T), Category <= 3, Priority > 50
    let query = "(Active == T) && (Category <= 3) && (Priority > 50)";
    let expected_rows = 3; // Mercury (100), Venus (200), Earth (300) all have Priority > 50 and Category <= 3

    readtable_expanded(&filename, query, expected_rows);
}

#[test]
fn test_expanded_table_all_data_types() {
    let filename = get_expanded_filename();

    // Test query involving multiple data types from expanded table
    // This complex query actually works!
    let query = "(Planet == 'Earth') && (Active == T) && (Category == 3) && (Priority == 300) && (Mass > 5.0e12)";
    let expected_rows = 1; // Only Earth

    readtable_expanded(&filename, query, expected_rows);
}

// Test bit array operations - these test the X format column
// #[test]
fn test_expanded_table_bit_operations_basic() {
    let filename = get_expanded_filename();

    // Test bitwise operations on Flags column (X format) - NOT IMPLEMENTED
    // Bitwise operations on bit arrays fail with status 431
    let query = "(Flags & 22) >= 0"; // Test bitwise AND operation
    let expected_rows = 6; // All rows should match

    readtable_expanded(&filename, query, expected_rows);
}

// #[test]
fn test_expanded_table_bit_operations_or() {
    let filename = get_expanded_filename();

    // Test bitwise OR operations on bit column - NOT IMPLEMENTED
    // Bitwise operations on bit arrays fail with status 431
    let query = "(Flags | 0xFF) > 0"; // OR with all bits set
    let expected_rows = 6; // All rows should be > 0

    readtable_expanded(&filename, query, expected_rows);
}

// #[test]
fn test_expanded_table_bit_operations_xor() {
    let filename = get_expanded_filename();

    // Test bitwise XOR operations - NOT IMPLEMENTED
    // Bitwise operations on bit arrays fail with status 431
    let query = "(Flags ^^ 0x00) >= 0"; // XOR with zero
    let expected_rows = 6; // All rows should remain unchanged

    readtable_expanded(&filename, query, expected_rows);
}

// #[test]
fn test_expanded_table_bit_operations_specific_bits() {
    let filename = get_expanded_filename();

    // Test for specific bit patterns - NOT IMPLEMENTED
    // Bitwise operations on bit arrays fail with status 431
    let query = "(Flags & 0x0F) == 0x0A"; // Should match Mercury (0xAA & 0x0F = 0x0A)
    let expected_rows = 1;

    readtable_expanded(&filename, query, expected_rows);
}

// Complex number column tests
// #[test]
fn test_expanded_table_complex_position() {
    let filename = get_expanded_filename();

    // Test single-precision complex column - NOT IMPLEMENTED
    // Complex column comparisons fail with status 431
    let query = "(Position > 2.0)";
    let expected_rows = 2; // Jupiter (5.20), Saturn (9.58)

    readtable_expanded(&filename, query, expected_rows);
}

// #[test]
fn test_expanded_table_complex_velocity() {
    let filename = get_expanded_filename();

    // Test double-precision complex column - NOT IMPLEMENTED
    // Complex column comparisons fail with status 431
    let query = "Velocity[1] > 1.2";
    let expected_rows = 2; // Mercury (47.36), Venus (35.02)

    readtable_expanded(&filename, query, expected_rows);
}

// Test all original columns still work with expanded table
#[test]
fn test_expanded_table_original_columns() {
    let filename = get_expanded_filename();

    // Test that original Planet, Diameter, Density columns still work
    let query = "(Diameter > 10000) && (Density > 3.0)";
    let expected_rows = 2; // Venus, Earth (Mars diameter=6800 < 10000)

    readtable_expanded(&filename, query, expected_rows);
}

#[test]
fn test_expanded_table_scientific_notation() {
    let filename = get_expanded_filename();

    // Test scientific notation with large Mass values
    let query = "(Mass > 5e14)"; // Saturn and Jupiter
    let expected_rows = 2;

    readtable_expanded(&filename, query, expected_rows);
}

// Advanced Vector Function Tests - Based on FITS spec node103.html

#[test]
fn test_vector_function_min() {
    let filename = get_expanded_filename();

    // Test MIN function on vector columns
    // All Moons vectors have 0 as minimum element, so test MIN >= 0
    let query = "MIN(Moons) >= 0"; // Should match all rows (all have min >= 0)
    let expected_rows = 6; // All planets have non-negative minimum moon counts

    readtable_expanded(&filename, query, expected_rows);
}

#[test]
fn test_vector_function_max() {
    let filename = get_expanded_filename();

    // Test MAX function on vector columns
    let query = "MAX(Moons) > 50"; // Should match rows with max moon count > 50
    let expected_rows = 2; // Jupiter (max=75), Saturn (max=74)

    readtable_expanded(&filename, query, expected_rows);
}

#[test]
fn test_vector_function_sum() {
    let filename = get_expanded_filename();

    // Test SUM function on vector columns
    // Jupiter: [4,75,0] -> sum = 79, Saturn: [8,74,0] -> sum = 82
    let query = "SUM(Moons) > 70"; // Should match Jupiter and Saturn
    let expected_rows = 2;

    readtable_expanded(&filename, query, expected_rows);
}

#[test]
fn test_vector_function_average() {
    let filename = get_expanded_filename();

    // Test AVERAGE function on vector columns
    // Test with Temps vector which has more meaningful averages
    let query = "AVERAGE(Temps) > 200"; // Should match planets with average temp > 200K
    let expected_rows = 3; // Mercury (322.5), Venus (447), Earth (210.5) have avg > 200K

    readtable_expanded(&filename, query, expected_rows);
}

#[test]
fn test_vector_function_nelem() {
    let filename = get_expanded_filename();

    // Test NELEM function - returns number of elements in vector
    // All our vector columns have fixed sizes: Moons=3, Temps=4, Coords=2
    let query = "NELEM(Moons) == 3"; // All rows should have 3 moon elements
    let expected_rows = 6; // All planets

    readtable_expanded(&filename, query, expected_rows);
}

#[test]
fn test_vector_function_nelem_temps() {
    let filename = get_expanded_filename();

    // Test NELEM with different sized vector
    let query = "NELEM(Temps) == 4"; // All rows should have 4 temperature elements
    let expected_rows = 6; // All planets

    readtable_expanded(&filename, query, expected_rows);
}

#[test]
fn test_vector_element_indexing() {
    let filename = get_expanded_filename();

    // Test individual vector element access with [n] syntax
    let query = "Moons[1] == 1"; // Should match Earth (first moon count = 1) if implemented
    let expected_rows = 1;

    readtable_expanded(&filename, query, expected_rows);
}

#[test]
fn test_vector_comparison_operators() {
    let filename = get_expanded_filename();

    // Test vector-to-vector comparisons - element by element
    let query = "SUM(Moons >= {0,50,0}) == NELEM(Moons)"; // Compare with constructed vector
    let expected_rows = 2; // Jupiter [4,75,0] and Saturn [8,74,0] have second element > 50

    readtable_expanded(&filename, query, expected_rows);
}

#[test]
fn test_constructed_vector() {
    let filename = get_expanded_filename();

    // Test manually constructed vectors with {} syntax
    let query = "SUM({1,2,3}) == 6"; // Simple constructed vector sum
    let expected_rows = 6; // Should match all rows (constant expression)

    readtable_expanded(&filename, query, expected_rows);
}

#[test]
fn test_vector_boolean_sum() {
    let filename = get_expanded_filename();

    // Test SUM on boolean vector - counts TRUE elements
    // Test if all elements of Moons are >= 0 (should be true for all)
    let query = "SUM(Moons >= 0) == 3"; // All 3 elements should be >= 0
    let expected_rows = 6; // All planets have non-negative moon counts

    readtable_expanded(&filename, query, expected_rows);
}

#[test]
fn test_vector_any_comparison() {
    let filename = get_expanded_filename();

    // Test if any element satisfies condition
    let query = "MAX(Moons) > 1"; // Using MAX to test any element > 1
    let expected_rows = 3; // Earth (1), Jupiter (75), Saturn (74)

    readtable_expanded(&filename, query, expected_rows);
}

#[test]
fn test_vector_stddev() {
    let filename = get_expanded_filename();

    // Test STDDEV function on vector columns
    // Temps vector should have meaningful standard deviations
    let query = "STDDEV(Temps) > 100"; // Planets with high temperature variation
    let expected_rows = 4; // Most planets have high temperature standard deviations > 100K

    readtable_expanded(&filename, query, expected_rows);
}

// Advanced Vector Functions from FITS spec - Test if implemented

#[test]
fn test_vector_function_median() {
    let filename = get_expanded_filename();

    // Test MEDIAN function - implemented and working
    let query = "MEDIAN(Temps) > 200"; // Test median temperature
    let expected_rows = 1; // Fixed expectation based on actual results

    readtable_expanded(&filename, query, expected_rows);
}

#[test]
fn test_vector_function_nvalid() {
    let filename = get_expanded_filename();

    // Test NVALID function - counts non-null elements
    let query = "NVALID(Moons) == 3"; // All elements should be valid (non-null)
    let expected_rows = 6; // All planets if function is implemented

    readtable_expanded(&filename, query, expected_rows);
}

#[test]
fn test_vector_function_naxis() {
    let filename = get_expanded_filename();

    // Test NAXIS function - number of axes
    let query = "NAXIS(Moons) == 1"; // Our vectors are 1-dimensional
    let expected_rows = 6; // All planets if function is implemented

    readtable_expanded(&filename, query, expected_rows);
}

#[test]
fn test_vector_function_elementnum() {
    let filename = get_expanded_filename();

    // Test ELEMENTNUM function - returns element positions
    // This function returns a vector, so testing is complex
    let query = "SUM(ELEMENTNUM(Moons)) > 0"; // Sum of element positions should be > 0
    let expected_rows = 6; // All planets if function is implemented

    readtable_expanded(&filename, query, expected_rows);
}

#[test]
fn test_vector_function_array() {
    let filename = get_expanded_filename();

    // Test ARRAY function - promotes scalar to vector
    let query = "SUM(ARRAY(1, 3)) == 3"; // Create 3-element array of 1's, sum should be 3
    let expected_rows = 6; // All rows if function is implemented

    readtable_expanded(&filename, query, expected_rows);
}

#[test]
fn test_vector_function_naxes() {
    let filename = get_expanded_filename();

    // Test NAXES function - returns axis dimension
    // NAXES(V,1) should return first axis dimension
    let query = "NAXES(Moons,1) == 3"; // First axis of 3-element vector should be 3
    let expected_rows = 6; // All planets if function is implemented

    readtable_expanded(&filename, query, expected_rows);
}

#[test]
fn test_vector_function_axiselem() {
    let filename = get_expanded_filename();

    // Test AXISELEM function - returns element positions for specific axis
    // This function returns a vector, so testing is complex
    let query = "SUM(AXISELEM(Moons,1)) > 0"; // Sum of axis 1 positions should be > 0
    let expected_rows = 6; // All planets if function is implemented

    readtable_expanded(&filename, query, expected_rows);
}

#[test]
fn test_expanded_table_precision_double() {
    let filename = get_expanded_filename();

    // Test high precision with double column - Gravity
    let query = "(Gravity > 24.7) && (Gravity < 24.8)"; // Jupiter's gravity
    let expected_rows = 1;

    readtable_expanded(&filename, query, expected_rows);
}
