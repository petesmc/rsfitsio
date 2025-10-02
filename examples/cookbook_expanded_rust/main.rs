#![allow(deprecated)]

use std::ffi::{CStr, CString};
use std::fs::remove_file;
use std::process::ExitCode;
use std::process::exit;

use bytemuck::{cast_slice, cast_slice_mut};
use libc::{c_char, c_double, c_float, c_int, c_long, c_short, c_uchar, c_ushort};
use rsfitsio::STDERR;
use rsfitsio::aliases::rust_api::*;
use rsfitsio::fitsio::LONGLONG;
use rsfitsio::fitsio::{
    ASCII_TBL, BINARY_TBL, CASEINSEN, END_OF_FILE, FLEN_CARD, FLEN_VALUE, READONLY, READWRITE,
    TBIT, TBYTE, TCOMPLEX, TDBLCOMPLEX, TDOUBLE, TFLOAT, TLOGICAL, TLONG, TLONGLONG, TSHORT,
    TUSHORT, USHORT_IMG, fitsfile,
};
use rsfitsio::{KeywordDatatype, NullValue};

pub fn main() -> ExitCode {
    /*************************************************************************
       This is a simple main program that calls the following routines:

        writeimage        - write a FITS primary array image
        writeascii        - write a FITS ASCII table extension
        writebintable     - write a FITS binary table extension with all data types
        copyhdu           - copy a header/data unit from one FITS file to another
        selectrows        - copy selected row from one HDU to another
        readheader        - read and print the header keywords in every extension
        readimage         - read a FITS image and compute the min and max value
        read_ascii_table  - read columns of data from ASCII table
        read_binary_table - read columns of data from binary table

    **************************************************************************/

    writeimage();
    writeascii();
    writebintable();
    copyhdu();
    selectrows();
    readheader();
    readimage();
    read_ascii_table();
    read_binary_table();

    println!("\nAll the cfitsio cookbook routines ran successfully.");

    ExitCode::from(0)
}

/*--------------------------------------------------------------------------*/
/// Create a FITS primary array containing a 2-D image
fn writeimage() {
    let mut fptr: Option<Box<fitsfile>> = None;
    let mut status: c_int = 0;
    let fpixel: c_long = 1;

    /* initialize FITS image parameters */
    let filename = "cetestfil.fit"; /* name for new FITS file */
    let bitpix: c_int = USHORT_IMG; /* 16-bit unsigned short pixel values       */
    let naxis: c_long = 2; /* 2-dimensional image                            */
    let naxes: [c_long; 2] = [300, 200]; /* image is 300 pixels wide by 200 rows */
    let nelements: c_long = naxes[0] * naxes[1];

    /* allocate memory for the whole image */
    let array_size = (naxes[0] * naxes[1]) as usize;
    let mut array: Vec<c_ushort> = vec![0; array_size];

    let _ = remove_file(filename); /* Delete old file if it already exists */

    let filename_cstr = CString::new(filename).unwrap();

    if fits_create_file(
        &mut fptr,
        cast_slice(filename_cstr.to_bytes_with_nul()),
        &mut status,
    ) != 0
    {
        printerror(status); /* call printerror if error occurs */
    }

    /* write the required keywords for the primary array image.     */
    /* Since bitpix = USHORT_IMG, this will cause cfitsio to create */
    /* a FITS image with BITPIX = 16 (signed short integers) with   */
    /* BSCALE = 1.0 and BZERO = 32768.  This is the convention that */
    /* FITS uses to store unsigned integers.  Note that the BSCALE  */
    /* and BZERO keywords will be automatically written by cfitsio  */
    /* in this case.                                                */

    if let Some(ref mut fptr_box) = fptr
        && fits_create_img(fptr_box, bitpix, naxis as c_int, &naxes, &mut status) != 0
    {
        printerror(status);
    }

    /* initialize the values in the image with a linear ramp function */
    for jj in 0..naxes[1] {
        for ii in 0..naxes[0] {
            let index = (jj * naxes[0] + ii) as usize;
            array[index] = (ii + jj) as c_ushort;
        }
    }

    /* write the array of unsigned integers to the FITS file */
    if let Some(ref mut fptr_box) = fptr
        && fits_write_img(
            fptr_box,
            TUSHORT,
            fpixel as LONGLONG,
            nelements as LONGLONG,
            cast_slice::<c_ushort, u8>(&array),
            &mut status,
        ) != 0
    {
        printerror(status);
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

        if status != 0 {
            printerror(status);
        }
    }

    if let Some(fptr_box) = fptr
        && fits_close_file(fptr_box, &mut status) != 0
    {
        /* close the file */
        printerror(status);
    }
}

/*--------------------------------------------------------------------------*/
/// Create an ASCII table extension containing 3 columns and 6 rows
fn writeascii() {
    let mut fptr: Option<Box<fitsfile>> = None;
    let mut status: c_int = 0;
    let firstrow: c_long = 1;
    let firstelem: c_long = 1;

    let tfields: c_int = 3; /* table will have 3 columns */
    let nrows: LONGLONG = 6; /* table will have 6 rows    */

    let filename = "cetestfil.fit"; /* name for new FITS file */
    let extname = "PLANETS_ASCII"; /* extension name */

    /* define the name, datatype, and physical units for the 3 columns */
    let ttype_strs = [
        CString::new("Planet").unwrap(),
        CString::new("Diameter").unwrap(),
        CString::new("Density").unwrap(),
    ];
    let tform_strs = [
        CString::new("a8").unwrap(),
        CString::new("I6").unwrap(),
        CString::new("F4.2").unwrap(),
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

    /* open with write access the FITS file containing a primary array */
    if fits_open_file(
        &mut fptr,
        cast_slice(filename_cstr.to_bytes_with_nul()),
        READWRITE,
        &mut status,
    ) != 0
    {
        printerror(status);
    }

    /* append a new empty ASCII table onto the FITS file */
    if let Some(ref mut fptr_box) = fptr {
        let ttype_opts: Vec<Option<&[c_char]>> = ttype.iter().map(|s| Some(*s)).collect();
        let tunit_opts: Vec<Option<&[c_char]>> = tunit.iter().map(|s| Some(*s)).collect();

        if fits_create_tbl(
            fptr_box,
            ASCII_TBL,
            nrows,
            tfields,
            &ttype_opts,
            &tform,
            Some(&tunit_opts),
            Some(cast_slice(extname_cstr.to_bytes_with_nul())),
            &mut status,
        ) != 0
        {
            printerror(status);
        }
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
    }

    if let Some(fptr_box) = fptr
        && fits_close_file(fptr_box, &mut status) != 0
    {
        /* close the FITS file */
        printerror(status);
    }
}

/*--------------------------------------------------------------------------*/
/// Create a binary table extension with all data types (except P and Q)
///
/// IMPORTANT NOTE ON TBIT (Bit Array) COLUMNS:
/// ============================================
/// TBIT columns in FITS binary tables require special handling that differs from other data types:
///
/// 1. DATA REPRESENTATION:
///    - TBIT columns store individual bits, not bytes
///    - Each bit is represented as a logical value in memory:
///      * 0 (zero) = FALSE = bit will be 0
///      * Non-zero = TRUE = bit will be 1
///    - For an "8X" column format, each row contains 8 bits
///
/// 2. MEMORY LAYOUT:
///    - The data is passed as an array of c_char (logical values)
///    - For 6 rows with 8 bits each, you need 48 c_char elements total
///    - NOT 6 bytes as you might expect!
///
/// 3. API FUNCTIONS:
///    - Use fits_write_col() and fits_read_col() with TBIT datatype for TBIT columns
///    - Parameters are the same as other column types, but:
///      * datatype: TBIT
///      * nelem: Total number of bits to read/write (NOT number of rows)
///      * array: Array of logical values (c_char) where 0=FALSE, non-zero=TRUE
///
/// 4. PARAMETER COUNTING:
///    - The 'nelem' parameter is the TOTAL number of bits across all rows
///    - Example: 6 rows × 8 bits/row = 48 bits total
///    - This is different from other column types where nelem = number of rows
///
/// This example demonstrates proper TBIT handling alongside all other FITS data types.
fn writebintable() {
    let mut fptr: Option<Box<fitsfile>> = None;
    let mut status: c_int = 0;
    let mut hdutype: c_int = 0;
    let firstrow: c_long = 1;
    let firstelem: c_long = 1;

    let tfields: c_int = 14; /* table will have 14 columns */
    let nrows: LONGLONG = 6; /* table will have 6 rows    */

    let filename = "cetestfil.fit"; /* name for new FITS file */
    let extname = "PLANETS_Binary"; /* extension name */

    /* define the name, datatype, and physical units for all columns */
    let ttype_strs = [
        CString::new("Planet").unwrap(),   // A - Character string
        CString::new("Active").unwrap(),   // L - Logical
        CString::new("Flags").unwrap(),    // X - Bit array
        CString::new("Category").unwrap(), // B - Unsigned byte
        CString::new("Priority").unwrap(), // I - 16-bit integer
        CString::new("Diameter").unwrap(), // J - 32-bit integer
        CString::new("Mass").unwrap(),     // K - 64-bit integer
        CString::new("Density").unwrap(),  // E - Single-precision float
        CString::new("Gravity").unwrap(),  // D - Double-precision float
        CString::new("Position").unwrap(), // C - Single-precision complex
        CString::new("Velocity").unwrap(), // M - Double-precision complex
        CString::new("Moons").unwrap(),    // 3I - Vector of 3 16-bit integers
        CString::new("Temps").unwrap(),    // 4E - Vector of 4 floats
        CString::new("Coords").unwrap(),   // 2D - Vector of 2 doubles
    ];
    let tform_strs = [
        CString::new("8A").unwrap(), // Character string
        CString::new("1L").unwrap(), // Logical
        CString::new("8X").unwrap(), // 8 bits
        CString::new("1B").unwrap(), // Unsigned byte
        CString::new("1I").unwrap(), // 16-bit integer
        CString::new("1J").unwrap(), // 32-bit integer
        CString::new("1K").unwrap(), // 64-bit integer
        CString::new("1E").unwrap(), // Single-precision float
        CString::new("1D").unwrap(), // Double-precision float
        CString::new("1C").unwrap(), // Single-precision complex
        CString::new("1M").unwrap(), // Double-precision complex
        CString::new("3I").unwrap(), // Vector of 3 16-bit integers
        CString::new("4E").unwrap(), // Vector of 4 floats
        CString::new("2D").unwrap(), // Vector of 2 doubles
    ];
    let tunit_strs = [
        CString::new("").unwrap(),
        CString::new("").unwrap(),
        CString::new("").unwrap(),
        CString::new("").unwrap(),
        CString::new("").unwrap(),
        CString::new("km").unwrap(),
        CString::new("kg").unwrap(),
        CString::new("g/cm^3").unwrap(),
        CString::new("m/s^2").unwrap(),
        CString::new("AU").unwrap(),
        CString::new("km/s").unwrap(),
        CString::new("count").unwrap(),
        CString::new("K").unwrap(),
        CString::new("deg").unwrap(),
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

    /* define the data for each column */
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

    // L - Logical values (1=T, 0=F for FITS internal representation)
    let active: [c_char; 6] = [1, 1, 1, 1, 1, 0];

    // X - TBIT columns (Bit arrays in FITS binary tables)
    // IMPORTANT: TBIT columns in FITS store bits as character arrays, not byte values!
    // Each bit is represented as a character '0' or '1' (ASCII 48 or 49).
    // For an "8X" column format, each row contains 8 bits, requiring 8 char elements per row.
    // The nelem parameter in fits_write_col/fits_read_col is the total number of bits, not rows.
    //
    // Here we define the bit patterns we want to write:
    // - flags_bytes: The actual bit patterns we want (as bytes for convenience)
    // - flags_chars: Will store the logical representation for FITS (0 or 1 values)
    let flags_bytes: [c_uchar; 6] = [
        0b10101010, // Mercury: alternating pattern
        0b11110000, // Venus: first 4 bits set
        0b00001111, // Earth: last 4 bits set
        0b11001100, // Mars: bits 0,1,4,5 set
        0b01010101, // Jupiter: alternating opposite pattern
        0b10011001, // Saturn: bits 0,3,4,7 set
    ];

    // Convert bit patterns to logical arrays for FITS TBIT format
    // For fits_write_col_bit, non-zero values are TRUE (bit set to 1)
    // and zero values are FALSE (bit set to 0)
    // Total: 48 characters (8 bits/row × 6 rows)
    let mut flags_chars: [c_char; 48] = [0; 48];
    for row in 0..6 {
        for bit in 0..8 {
            let bit_index = row * 8 + bit;
            // Check if bit is set in the byte
            // For FITS: non-zero = TRUE (1), zero = FALSE (0)
            if (flags_bytes[row] >> (7 - bit)) & 1 == 1 {
                flags_chars[bit_index] = 1; // TRUE - bit will be set to 1
            } else {
                flags_chars[bit_index] = 0; // FALSE - bit will be set to 0
            }
        }
    }

    // B - Unsigned bytes
    let category: [c_uchar; 6] = [1, 2, 3, 4, 5, 6];

    // I - 16-bit integers
    let priority: [c_short; 6] = [100, 200, 300, 150, 50, 75];

    // J - 32-bit integers
    let diameter: [c_long; 6] = [4880, 12112, 12742, 6800, 143000, 121000];

    // K - 64-bit integers
    let mass: [LONGLONG; 6] = [
        330110000000,
        4867500000000,
        5972400000000,
        641710000000,
        1898200000000000,
        568340000000000,
    ];

    // E - Single-precision floats
    let density: [c_float; 6] = [5.1, 5.3, 5.52, 3.94, 1.33, 0.69];

    // D - Double-precision floats
    let gravity: [c_double; 6] = [3.7, 8.87, 9.807, 3.71, 24.79, 10.44];

    // C - Single-precision complex (real, imaginary pairs)
    let position: [[c_float; 2]; 6] = [
        [0.39, 5.0], // Mercury
        [0.72, 4.0], // Venus
        [1.00, 3.0], // Earth
        [1.52, 2.0], // Mars
        [5.20, 1.0], // Jupiter
        [9.58, 0.0], // Saturn
    ];

    // M - Double-precision complex (real, imaginary pairs)
    let velocity: [[c_double; 2]; 6] = [
        [47.36, 0.0], // Mercury orbital velocity
        [35.02, 1.0], // Venus
        [29.78, 2.0], // Earth
        [24.07, 3.0], // Mars
        [13.07, 4.0], // Jupiter
        [9.68, 5.0],  // Saturn
    ];

    // 3I - Vector of 3 16-bit integers (number of major moons in different size categories)
    let moons: [[c_short; 3]; 6] = [
        [0, 0, 0],  // Mercury: no moons
        [0, 0, 0],  // Venus: no moons
        [1, 0, 0],  // Earth: 1 major moon
        [0, 2, 0],  // Mars: 2 small moons
        [4, 75, 0], // Jupiter: 4 Galilean + many others
        [8, 74, 0], // Saturn: many moons
    ];

    // 4E - Vector of 4 floats (min, max, mean, std dev temperatures in Kelvin)
    let temps: [[c_float; 4]; 6] = [
        [100.0, 700.0, 340.0, 150.0], // Mercury
        [228.0, 773.0, 737.0, 50.0],  // Venus
        [184.0, 330.0, 288.0, 40.0],  // Earth
        [130.0, 308.0, 210.0, 60.0],  // Mars
        [88.0, 165.0, 124.0, 20.0],   // Jupiter
        [82.0, 134.0, 95.0, 15.0],    // Saturn
    ];

    // 2D - Vector of 2 doubles (latitude, longitude of notable feature)
    let coords: [[c_double; 2]; 6] = [
        [15.0, 45.0],   // Mercury: Caloris Basin
        [-5.0, 165.0],  // Venus: Maxwell Montes
        [0.0, 0.0],     // Earth: Prime Meridian
        [-14.6, 175.5], // Mars: Olympus Mons
        [-22.0, 115.0], // Jupiter: Great Red Spot
        [0.0, 0.0],     // Saturn: Equator
    ];

    let filename_cstr = CString::new(filename).unwrap();
    let extname_cstr = CString::new(extname).unwrap();

    /* open the FITS file containing a primary array and an ASCII table */
    if fits_open_file(
        &mut fptr,
        cast_slice(filename_cstr.to_bytes_with_nul()),
        READWRITE,
        &mut status,
    ) != 0
    {
        printerror(status);
    }

    if let Some(ref mut fptr_box) = fptr
        && fits_movabs_hdu(fptr_box, 2, Some(&mut hdutype), &mut status) != 0
    {
        /* move to 2nd HDU */
        printerror(status);
    }

    /* append a new empty binary table onto the FITS file */
    if let Some(ref mut fptr_box) = fptr {
        let ttype_opts: Vec<Option<&[c_char]>> = ttype.iter().map(|s| Some(*s)).collect();
        let tunit_opts: Vec<Option<&[c_char]>> = tunit.iter().map(|s| Some(*s)).collect();

        if fits_create_tbl(
            fptr_box,
            BINARY_TBL,
            nrows,
            tfields,
            &ttype_opts,
            &tform,
            Some(&tunit_opts),
            Some(cast_slice(extname_cstr.to_bytes_with_nul())),
            &mut status,
        ) != 0
        {
            printerror(status);
        }
    }

    /* write data to all columns */
    if let Some(ref mut fptr_box) = fptr {
        // Column 1: A - Character strings
        fits_write_col_str(
            fptr_box,
            1,
            firstrow as LONGLONG,
            firstelem as LONGLONG,
            planet.len() as LONGLONG,
            &planet,
            &mut status,
        );

        // Column 2: L - Logical
        fits_write_col(
            fptr_box,
            TLOGICAL,
            2,
            firstrow as LONGLONG,
            firstelem as LONGLONG,
            active.len() as LONGLONG,
            cast_slice(&active),
            &mut status,
        );

        // Column 3: X - Bit array (TBIT)
        // For TBIT columns with fits_write_col:
        // - Use TBIT as the datatype
        // - nelem parameter is the total number of BITS, not rows
        // - Data is passed as logical values: 0=FALSE (bit 0), non-zero=TRUE (bit 1)
        // - For 6 rows × 8 bits = 48 total bits
        fits_write_col(
            fptr_box,
            TBIT,
            3,
            firstrow as LONGLONG,
            firstelem as LONGLONG,
            48 as LONGLONG, // Total number of bits (8 bits/row × 6 rows)
            cast_slice(&flags_chars),
            &mut status,
        );

        // Column 4: B - Unsigned byte
        fits_write_col(
            fptr_box,
            TBYTE,
            4,
            firstrow as LONGLONG,
            firstelem as LONGLONG,
            category.len() as LONGLONG,
            cast_slice(&category),
            &mut status,
        );

        // Column 5: I - 16-bit integer
        fits_write_col(
            fptr_box,
            TSHORT,
            5,
            firstrow as LONGLONG,
            firstelem as LONGLONG,
            priority.len() as LONGLONG,
            cast_slice(&priority),
            &mut status,
        );

        // Column 6: J - 32-bit integer
        fits_write_col(
            fptr_box,
            TLONG,
            6,
            firstrow as LONGLONG,
            firstelem as LONGLONG,
            diameter.len() as LONGLONG,
            cast_slice(&diameter),
            &mut status,
        );

        // Column 7: K - 64-bit integer
        fits_write_col(
            fptr_box,
            TLONGLONG,
            7,
            firstrow as LONGLONG,
            firstelem as LONGLONG,
            mass.len() as LONGLONG,
            cast_slice(&mass),
            &mut status,
        );

        // Column 8: E - Single-precision float
        fits_write_col(
            fptr_box,
            TFLOAT,
            8,
            firstrow as LONGLONG,
            firstelem as LONGLONG,
            density.len() as LONGLONG,
            cast_slice(&density),
            &mut status,
        );

        // Column 9: D - Double-precision float
        fits_write_col(
            fptr_box,
            TDOUBLE,
            9,
            firstrow as LONGLONG,
            firstelem as LONGLONG,
            gravity.len() as LONGLONG,
            cast_slice(&gravity),
            &mut status,
        );

        // Column 10: C - Single-precision complex
        fits_write_col(
            fptr_box,
            TCOMPLEX,
            10,
            firstrow as LONGLONG,
            firstelem as LONGLONG,
            position.len() as LONGLONG,
            cast_slice(&position),
            &mut status,
        );

        // Column 11: M - Double-precision complex
        fits_write_col(
            fptr_box,
            TDBLCOMPLEX,
            11,
            firstrow as LONGLONG,
            firstelem as LONGLONG,
            velocity.len() as LONGLONG,
            cast_slice(&velocity),
            &mut status,
        );

        // Column 12: 3I - Vector of 3 16-bit integers
        fits_write_col(
            fptr_box,
            TSHORT,
            12,
            firstrow as LONGLONG,
            firstelem as LONGLONG,
            (moons.len() * 3) as LONGLONG, // Total elements
            cast_slice(&moons),
            &mut status,
        );

        // Column 13: 4E - Vector of 4 floats
        fits_write_col(
            fptr_box,
            TFLOAT,
            13,
            firstrow as LONGLONG,
            firstelem as LONGLONG,
            (temps.len() * 4) as LONGLONG, // Total elements
            cast_slice(&temps),
            &mut status,
        );

        // Column 14: 2D - Vector of 2 doubles
        fits_write_col(
            fptr_box,
            TDOUBLE,
            14,
            firstrow as LONGLONG,
            firstelem as LONGLONG,
            (coords.len() * 2) as LONGLONG, // Total elements
            cast_slice(&coords),
            &mut status,
        );
    }

    if let Some(fptr_box) = fptr
        && fits_close_file(fptr_box, &mut status) != 0
    {
        /* close the FITS file */
        printerror(status);
    }
}

/*--------------------------------------------------------------------------*/
/// Copy the 1st and 3rd HDUs from the input file to a new FITS file
fn copyhdu() {
    let mut infptr: Option<Box<fitsfile>> = None;
    let mut outfptr: Option<Box<fitsfile>> = None;

    let infilename = "cetestfil.fit"; /* name for existing FITS file   */
    let outfilename = "dtestfil.fit"; /* name for new FITS file        */

    let mut status: c_int = 0;
    let morekeys: c_int = 0;
    let mut hdutype: c_int = 0;

    let _ = remove_file(outfilename); /* Delete old file if it already exists */

    let infilename_cstr = CString::new(infilename).unwrap();
    let outfilename_cstr = CString::new(outfilename).unwrap();

    /* open the existing FITS file */
    if fits_open_file(
        &mut infptr,
        cast_slice(infilename_cstr.to_bytes_with_nul()),
        READONLY,
        &mut status,
    ) != 0
    {
        printerror(status);
    }

    if fits_create_file(
        &mut outfptr,
        cast_slice(outfilename_cstr.to_bytes_with_nul()),
        &mut status,
    ) != 0
    {
        /*create FITS file*/
        printerror(status); /* call printerror if error occurs */
    }

    /* copy the primary array from the input file to the output file */
    if let (Some(infptr_box), Some(outfptr_box)) = (&mut infptr, &mut outfptr)
        && fits_copy_hdu(infptr_box, outfptr_box, morekeys, &mut status) != 0
    {
        printerror(status);
    }

    /* move to the 3rd HDU in the input file */
    if let Some(ref mut infptr_box) = infptr
        && fits_movabs_hdu(infptr_box, 3, Some(&mut hdutype), &mut status) != 0
    {
        printerror(status);
    }

    /* copy 3rd HDU from the input file to the output file (to 2nd HDU) */
    if let (Some(infptr_box), Some(outfptr_box)) = (&mut infptr, &mut outfptr)
        && fits_copy_hdu(infptr_box, outfptr_box, morekeys, &mut status) != 0
    {
        printerror(status);
    }

    let mut close_status = 0;
    if let Some(outfptr_box) = outfptr
        && fits_close_file(outfptr_box, &mut close_status) != 0
    {
        printerror(close_status);
    }
    if let Some(infptr_box) = infptr
        && fits_close_file(infptr_box, &mut status) != 0
    {
        /* close files */
        printerror(status);
    }
}

/*--------------------------------------------------------------------------*/
/// select rows from an input table and copy them to the output table
fn selectrows() {
    let mut infptr: Option<Box<fitsfile>> = None;
    let mut outfptr: Option<Box<fitsfile>> = None;

    let mut card = [0u8; FLEN_CARD];
    let mut status: c_int = 0;
    let mut hdutype: c_int = 0;
    let mut nkeys: c_int = 0;
    let mut keypos: c_int = 0;
    let mut nfound: c_int = 0;
    let mut colnum: c_int = 0;
    let mut anynulls: c_int = 0;
    let mut naxes: [c_long; 2] = [0; 2];
    let frow: LONGLONG = 1;
    let felem: LONGLONG = 1;
    let mut noutrows: c_long;
    let nullval: c_float = -99.0;
    let mut density: [c_float; 6] = [0.0; 6];

    let infilename = "cetestfil.fit"; /* name for existing FITS file   */
    let outfilename = "dtestfil.fit"; /* name for new FITS file        */

    let infilename_cstr = CString::new(infilename).unwrap();
    let outfilename_cstr = CString::new(outfilename).unwrap();

    /* open the existing FITS files */
    if fits_open_file(
        &mut infptr,
        cast_slice(infilename_cstr.to_bytes_with_nul()),
        READONLY,
        &mut status,
    ) != 0
        || fits_open_file(
            &mut outfptr,
            cast_slice(outfilename_cstr.to_bytes_with_nul()),
            READWRITE,
            &mut status,
        ) != 0
    {
        printerror(status);
    }

    /* move to the 3rd HDU in the input file (a binary table in this case) */
    if let Some(ref mut infptr_box) = infptr
        && fits_movabs_hdu(infptr_box, 3, Some(&mut hdutype), &mut status) != 0
    {
        printerror(status);
    }

    if hdutype != BINARY_TBL {
        println!("Error: expected to find a binary table in this HDU");
        return;
    }

    /* move to the last (2rd) HDU in the output file */
    if let Some(ref mut outfptr_box) = outfptr
        && fits_movabs_hdu(outfptr_box, 2, Some(&mut hdutype), &mut status) != 0
    {
        printerror(status);
    }

    /* create new extension in the output file */
    if let Some(ref mut outfptr_box) = outfptr
        && fits_create_hdu(outfptr_box, &mut status) != 0
    {
        printerror(status);
    }

    /* get number of keywords */
    if let Some(ref mut infptr_box) = infptr
        && fits_get_hdrpos(infptr_box, Some(&mut nkeys), Some(&mut keypos), &mut status) != 0
    {
        printerror(status);
    }

    /* copy all the keywords from the input to the output extension */
    for ii in 1..=nkeys {
        if let (Some(infptr_box), Some(outfptr_box)) = (&mut infptr, &mut outfptr) {
            fits_read_record(infptr_box, ii, Some(cast_slice_mut(&mut card)), &mut status);

            fits_write_record(outfptr_box, cast_slice(&card), &mut status);
        }
    }

    /* read the NAXIS1 and NAXIS2 keyword to get table size */
    if let Some(ref mut infptr_box) = infptr {
        let naxis_key = CString::new("NAXIS").unwrap();
        if fits_read_keys_lng(
            infptr_box,
            cast_slice(naxis_key.to_bytes_with_nul()),
            1,
            2,
            &mut naxes,
            &mut nfound,
            &mut status,
        ) != 0
        {
            printerror(status);
        }
    }

    /* find which column contains the DENSITY values */
    if let Some(ref mut infptr_box) = infptr {
        let density_key = CString::new("density").unwrap();
        if fits_get_colnum(
            infptr_box,
            CASEINSEN as c_int,
            cast_slice(density_key.to_bytes_with_nul()),
            &mut colnum,
            &mut status,
        ) != 0
        {
            printerror(status);
        }
    }

    /* read the DENSITY column values */
    if let Some(ref mut infptr_box) = infptr
        && fits_read_col(
            infptr_box,
            TFLOAT,
            colnum,
            frow as LONGLONG,
            felem as LONGLONG,
            naxes[1] as LONGLONG,
            Some(NullValue::Float(nullval)),
            cast_slice_mut(&mut density),
            Some(&mut anynulls),
            &mut status,
        ) != 0
    {
        printerror(status);
    }

    /* allocate buffer large enough for 1 row of the table - use Vec instead of malloc */
    let mut buffer: Vec<u8> = vec![0; naxes[0] as usize];

    /* If the density is less than 3.0, copy the row to the output table */
    noutrows = 0;
    for irow in 1..=naxes[1] {
        if density[(irow - 1) as usize] < 3.0 {
            noutrows += 1;
            if let (Some(infptr_box), Some(outfptr_box)) = (&mut infptr, &mut outfptr) {
                fits_read_tblbytes(
                    infptr_box,
                    irow as LONGLONG,
                    1,
                    naxes[0] as LONGLONG,
                    &mut buffer,
                    &mut status,
                );

                fits_write_tblbytes(
                    outfptr_box,
                    noutrows as LONGLONG,
                    1,
                    naxes[0] as LONGLONG,
                    &buffer,
                    &mut status,
                );
            }
        }
    }

    /* update the NAXIS2 keyword with the correct number of rows */
    if let Some(ref mut outfptr_box) = outfptr {
        let naxis2_key = CString::new("NAXIS2").unwrap();
        if fits_update_key(
            outfptr_box,
            KeywordDatatype::TLONG(&noutrows),
            cast_slice(naxis2_key.to_bytes_with_nul()),
            None,
            &mut status,
        ) != 0
        {
            printerror(status);
        }
    }

    let mut close_status = 0;
    if let Some(outfptr_box) = outfptr
        && fits_close_file(outfptr_box, &mut close_status) != 0
    {
        printerror(close_status);
    }
    if let Some(infptr_box) = infptr
        && fits_close_file(infptr_box, &mut status) != 0
    {
        printerror(status);
    }
}

/*--------------------------------------------------------------------------*/
/// Print out all the header keywords in all extensions of a FITS file
fn readheader() {
    let mut fptr: Option<Box<fitsfile>> = None;

    let mut status: c_int = 0;
    let mut nkeys: c_int = 0;
    let mut keypos: c_int = 0;
    let mut hdutype: c_int = 0;
    let mut ii: c_int;
    let filename = "cetestfil.fit"; /* name of existing FITS file   */
    let mut card = [0u8; FLEN_CARD]; /* standard string lengths defined in fitsioc.h */

    let filename_cstr = CString::new(filename).unwrap();

    if fits_open_file(
        &mut fptr,
        cast_slice(filename_cstr.to_bytes_with_nul()),
        READONLY,
        &mut status,
    ) != 0
    {
        printerror(status);
    }

    /* attempt to move to next HDU, until we get an EOF error */
    ii = 1;
    while let Some(ref mut fptr_box) = fptr {
        if fits_movabs_hdu(fptr_box, ii, Some(&mut hdutype), &mut status) != 0 {
            break;
        }

        /* get no. of keywords */
        if fits_get_hdrpos(fptr_box, Some(&mut nkeys), Some(&mut keypos), &mut status) != 0 {
            printerror(status);
        }

        println!("Header listing for HDU #{ii}:");
        for jj in 1..=nkeys {
            if fits_read_record(fptr_box, jj, Some(cast_slice_mut(&mut card)), &mut status) != 0 {
                printerror(status);
            }

            let card_str = CStr::from_bytes_until_nul(&card).unwrap();
            println!("{}", card_str.to_string_lossy()); /* print the keyword card */
        }
        println!("END\n"); /* terminate listing with END */
        ii += 1;
    }

    if status == END_OF_FILE {
        /* status values are defined in fitsioc.h */
        status = 0; /* got the expected EOF error; reset = 0  */
    } else {
        printerror(status); /* got an unexpected error                */
    }

    if let Some(fptr_box) = fptr
        && fits_close_file(fptr_box, &mut status) != 0
    {
        printerror(status);
    }
}

/*--------------------------------------------------------------------------*/
/// Read a FITS image and determine the minimum and maximum pixel values
fn readimage() {
    let mut fptr: Option<Box<fitsfile>> = None;
    let mut status: c_int = 0;
    let mut nfound: c_int = 0;
    let mut naxes: [c_long; 2] = [0; 2];
    let mut fpixel: LONGLONG;
    let mut nbuffer: LONGLONG;
    let mut npixels: LONGLONG;

    const BUFFSIZE: usize = 1000;
    let mut datamin: c_float = 1.0e30;
    let mut datamax: c_float = -1.0e30;
    let nullval: c_float = 0.0;
    let mut buffer: [c_float; BUFFSIZE] = [0.0; BUFFSIZE];
    let filename = "cetestfil.fit"; /* name of existing FITS file   */

    let filename_cstr = CString::new(filename).unwrap();

    if fits_open_file(
        &mut fptr,
        cast_slice(filename_cstr.to_bytes_with_nul()),
        READONLY,
        &mut status,
    ) != 0
    {
        printerror(status);
    }

    /* read the NAXIS1 and NAXIS2 keyword to get image size */
    if let Some(ref mut fptr_box) = fptr {
        let naxis_key = CString::new("NAXIS").unwrap();
        if fits_read_keys_lng(
            fptr_box,
            cast_slice(naxis_key.to_bytes_with_nul()),
            1,
            2,
            &mut naxes,
            &mut nfound,
            &mut status,
        ) != 0
        {
            printerror(status);
        }
    }

    npixels = (naxes[0] * naxes[1]) as LONGLONG; /* number of pixels in the image */
    fpixel = 1;

    while npixels > 0 {
        nbuffer = npixels;
        if npixels > BUFFSIZE as LONGLONG {
            nbuffer = BUFFSIZE as LONGLONG; /* read as many pixels as will fit in buffer */
        }

        /* Note that even though the FITS images contains unsigned integer */
        /* pixel values (or more accurately, signed integer pixels with    */
        /* a bias of 32768),  this routine is reading the values into a    */
        /* float array.   Cfitsio automatically performs the datatype      */
        /* conversion in cases like this.                                  */

        if let Some(ref mut fptr_box) = fptr {
            let buffer_slice = &mut buffer[0..nbuffer as usize];
            if fits_read_img(
                fptr_box,
                TFLOAT,
                fpixel,
                nbuffer as LONGLONG,
                Some(NullValue::Float(nullval)),
                cast_slice_mut(buffer_slice),
                None,
                &mut status,
            ) != 0
            {
                printerror(status);
            }
        }

        for ii in 0..nbuffer {
            if buffer[ii as usize] < datamin {
                datamin = buffer[ii as usize];
            }

            if buffer[ii as usize] > datamax {
                datamax = buffer[ii as usize];
            }
        }
        npixels -= nbuffer; /* increment remaining number of pixels */
        fpixel += nbuffer; /* next pixel to be read in image */
    }

    println!("\nMin and max image pixels =  {datamin:.0}, {datamax:.0}");

    if let Some(fptr_box) = fptr
        && fits_close_file(fptr_box, &mut status) != 0
    {
        printerror(status);
    }
}

/*--------------------------------------------------------------------------*/
/// Read and print data values from an ASCII table
fn read_ascii_table() {
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

    let filename = "cetestfil.fit"; /* name of existing FITS file   */

    let filename_cstr = CString::new(filename).unwrap();

    if fits_open_file(
        &mut fptr,
        cast_slice(filename_cstr.to_bytes_with_nul()),
        READONLY,
        &mut status,
    ) != 0
    {
        printerror(status);
    }

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

    /* move to the ASCII table HDU */
    if let Some(ref mut fptr_box) = fptr
        && fits_movabs_hdu(fptr_box, 2, Some(&mut hdutype), &mut status) != 0
    {
        printerror(status);
    }

    if hdutype != ASCII_TBL {
        println!("Error: expected to find an ASCII table in HDU 2");
        return;
    }

    println!("\nReading ASCII table in HDU 2:");

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
        if status != 0 {
            printerror(status);
        }
    }

    let ttype0 = CStr::from_bytes_until_nul(&ttype_vecs[0])
        .unwrap()
        .to_string_lossy();
    let ttype1 = CStr::from_bytes_until_nul(&ttype_vecs[1])
        .unwrap()
        .to_string_lossy();
    let ttype2 = CStr::from_bytes_until_nul(&ttype_vecs[2])
        .unwrap()
        .to_string_lossy();
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
    }

    for ii in 0..6 {
        let name_str = CStr::from_bytes_until_nul(&name_vecs[ii])
            .unwrap()
            .to_string_lossy();
        println!(
            "{:5} {:>10} {:>10} {:>10.2}",
            ii + 1,
            name_str,
            dia[ii],
            den[ii]
        );
    }

    if let Some(fptr_box) = fptr
        && fits_close_file(fptr_box, &mut status) != 0
    {
        printerror(status);
    }
}

/*--------------------------------------------------------------------------*/
/// Read and print data values from a binary table
fn read_binary_table() {
    let mut fptr: Option<Box<fitsfile>> = None;
    let mut status: c_int = 0;
    let mut hdutype: c_int = 0;
    let frow: c_long = 1;
    let felem: c_long = 1;
    let nelem: c_long = 6;
    let longnull: c_long = 0;
    let floatnull: c_float = 0.0;
    let doublenull: c_double = 0.0;

    // Arrays to store data from various columns
    let mut active: [c_char; 6] = [0; 6];
    // For TBIT columns: We need 48 logical values (8 bits/row × 6 rows)
    // Each bit is represented as a logical value: 0 = FALSE (bit is 0), non-zero = TRUE (bit is 1)
    let mut flags_chars: [c_char; 48] = [0; 48]; // X - Bit array (stored as logical values)
    let mut flags_bytes: [c_uchar; 6] = [0; 6]; // Converted byte representation for display
    let mut category: [c_uchar; 6] = [0; 6];
    let mut priority: [c_short; 6] = [0; 6];
    let mut diameter: [c_long; 6] = [0; 6];
    let mut mass: [LONGLONG; 6] = [0; 6];
    let mut density: [c_float; 6] = [0.0; 6];
    let mut gravity: [c_double; 6] = [0.0; 6];
    let mut position: [[c_float; 2]; 6] = [[0.0; 2]; 6]; // C - Single complex
    let mut velocity: [[c_double; 2]; 6] = [[0.0; 2]; 6]; // M - Double complex
    let mut moons: [[c_short; 3]; 6] = [[0; 3]; 6];
    let mut temps: [[c_float; 4]; 6] = [[0.0; 4]; 6]; // 4E - Vector of 4 floats
    let mut coords: [[c_double; 2]; 6] = [[0.0; 2]; 6]; // 2D - Vector of 2 doubles

    let filename = "cetestfil.fit"; /* name of existing FITS file   */

    let filename_cstr = CString::new(filename).unwrap();

    if fits_open_file(
        &mut fptr,
        cast_slice(filename_cstr.to_bytes_with_nul()),
        READONLY,
        &mut status,
    ) != 0
    {
        printerror(status);
    }

    /* allocate space for the column labels */
    let mut ttype_vecs: Vec<Vec<u8>> = Vec::new();
    for _i in 0..9 {
        // We'll read 9 column headers
        ttype_vecs.push(vec![0; FLEN_VALUE]);
    }

    /* allocate space for string column values */
    let mut name_vecs: Vec<Vec<u8>> = Vec::new();
    for _i in 0..6 {
        name_vecs.push(vec![0; 10]);
    }

    /* move to the binary table HDU */
    if let Some(ref mut fptr_box) = fptr
        && fits_movabs_hdu(fptr_box, 3, Some(&mut hdutype), &mut status) != 0
    {
        printerror(status);
    }

    if hdutype != BINARY_TBL {
        println!("Error: expected to find a binary table in HDU 3");
        return;
    }

    println!("\nReading binary table in HDU 3:");

    /* read the column names from specific TTYPEn keywords */
    if let Some(ref mut fptr_box) = fptr {
        // Read column names for multiple data types: Planet, Active, Category, Priority, Diameter, Mass, Density, Gravity, Moons
        let keys = [
            "TTYPE1", "TTYPE2", "TTYPE4", "TTYPE5", "TTYPE6", "TTYPE7", "TTYPE8", "TTYPE9",
            "TTYPE12",
        ];
        for (i, key) in keys.iter().enumerate() {
            let ttype_key = CString::new(*key).unwrap();
            let mut value = vec![0u8; FLEN_VALUE];
            fits_read_keyword(
                fptr_box,
                cast_slice(ttype_key.to_bytes_with_nul()),
                cast_slice_mut(&mut value),
                None,
                &mut status,
            );
            if status == 0 {
                // Extract the value between quotes
                let value_str = CStr::from_bytes_until_nul(&value)
                    .unwrap()
                    .to_string_lossy();
                let trimmed = value_str.trim().trim_matches('\'').trim();
                let trimmed_bytes = CString::new(trimmed).unwrap_or(CString::new("").unwrap());
                let bytes = trimmed_bytes.to_bytes_with_nul();
                ttype_vecs[i][..bytes.len().min(FLEN_VALUE)]
                    .copy_from_slice(&bytes[..bytes.len().min(FLEN_VALUE)]);
            }
        }
        if status != 0 {
            printerror(status);
        }
    }

    // Print header showing all columns with exact alignment to match data
    println!("\nBinary table columns and their data types:");
    println!(
        "{:>8} {:>5} {:>8} {:>5} {:>5} {:>11} {:>12} {:>10} {:>10} {:>7} {:>8} {:>11} {:>18} {:>13}",
        "Planet(A)",
        "Act(L)",
        "Flags(X)",
        "Cat(B)",
        "Pri(I)",
        "Diameter(J)",
        "Mass(K)",
        "Density(E)",
        "Gravity(D)",
        "Pos(C)",
        "Vel(M)",
        "Moons(3I)",
        "Temps(4E)",
        "Coords(2D)"
    );
    println!(
        "{:>8} {:>5} {:>8} {:>5} {:>5} {:>11} {:>12} {:>10} {:>10} {:>7} {:>8} {:>11} {:>18} {:>13}",
        "--------",
        "-----",
        "----------",
        "-----",
        "-----",
        "-----------",
        "------------",
        "----------",
        "----------",
        "-------",
        "--------",
        "-----------",
        "------------------",
        "-------------"
    );

    /* read multiple columns showcasing different data types */
    if let Some(ref mut fptr_box) = fptr {
        let mut name_refs: Vec<&mut [u8]> =
            name_vecs.iter_mut().map(|v| v.as_mut_slice()).collect();

        // Read column 1: Planet (A - Character strings)
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

        // Read column 2: Active (L - Logical)
        fits_read_col(
            fptr_box,
            TLOGICAL,
            2,
            frow as LONGLONG,
            felem as LONGLONG,
            nelem as LONGLONG,
            None,
            cast_slice_mut(&mut active),
            None,
            &mut status,
        );

        // Read column 3: Flags (X - Bit array)
        // For TBIT columns with fits_read_col:
        // - Use TBIT as the datatype
        // - nelem parameter is the total number of BITS, not rows
        // - Data is returned as logical values: 0=FALSE (bit was 0), non-zero=TRUE (bit was 1)
        fits_read_col(
            fptr_box,
            TBIT,
            3,
            frow as LONGLONG,
            felem as LONGLONG,
            48 as LONGLONG, // Total number of bits (8 bits/row × 6 rows)
            None,
            cast_slice_mut(&mut flags_chars),
            None,
            &mut status,
        );

        // Convert the logical values back to bytes for easier handling
        // Non-zero values mean the bit is set to 1
        for row in 0..6 {
            let mut byte_val: c_uchar = 0;
            for bit in 0..8 {
                let char_index = row * 8 + bit;
                if flags_chars[char_index] != 0 {
                    // Non-zero = TRUE (bit is 1)
                    byte_val |= 1 << (7 - bit);
                }
            }
            flags_bytes[row] = byte_val;
        }

        // Read column 4: Category (B - Unsigned byte)
        fits_read_col(
            fptr_box,
            TBYTE,
            4,
            frow as LONGLONG,
            felem as LONGLONG,
            nelem as LONGLONG,
            None,
            cast_slice_mut(&mut category),
            None,
            &mut status,
        );

        // Read column 5: Priority (I - 16-bit integer)
        fits_read_col(
            fptr_box,
            TSHORT,
            5,
            frow as LONGLONG,
            felem as LONGLONG,
            nelem as LONGLONG,
            None,
            cast_slice_mut(&mut priority),
            None,
            &mut status,
        );

        // Read column 6: Diameter (J - 32-bit integer)
        fits_read_col(
            fptr_box,
            TLONG,
            6,
            frow as LONGLONG,
            felem as LONGLONG,
            nelem as LONGLONG,
            Some(NullValue::Long(longnull)),
            cast_slice_mut(&mut diameter),
            None,
            &mut status,
        );

        // Read column 7: Mass (K - 64-bit integer)
        fits_read_col(
            fptr_box,
            TLONGLONG,
            7,
            frow as LONGLONG,
            felem as LONGLONG,
            nelem as LONGLONG,
            None,
            cast_slice_mut(&mut mass),
            None,
            &mut status,
        );

        // Read column 8: Density (E - Single-precision float)
        fits_read_col(
            fptr_box,
            TFLOAT,
            8,
            frow as LONGLONG,
            felem as LONGLONG,
            nelem as LONGLONG,
            Some(NullValue::Float(floatnull)),
            cast_slice_mut(&mut density),
            None,
            &mut status,
        );

        // Read column 9: Gravity (D - Double-precision float)
        fits_read_col(
            fptr_box,
            TDOUBLE,
            9,
            frow as LONGLONG,
            felem as LONGLONG,
            nelem as LONGLONG,
            Some(NullValue::Double(doublenull)),
            cast_slice_mut(&mut gravity),
            None,
            &mut status,
        );

        // Read column 10: Position (C - Single-precision complex)
        fits_read_col(
            fptr_box,
            TCOMPLEX,
            10,
            frow as LONGLONG,
            felem as LONGLONG,
            nelem as LONGLONG,
            None,
            cast_slice_mut(&mut position),
            None,
            &mut status,
        );

        // Read column 11: Velocity (M - Double-precision complex)
        fits_read_col(
            fptr_box,
            TDBLCOMPLEX,
            11,
            frow as LONGLONG,
            felem as LONGLONG,
            nelem as LONGLONG,
            None,
            cast_slice_mut(&mut velocity),
            None,
            &mut status,
        );

        // Read column 12: Moons (3I - Vector of 3 16-bit integers)
        fits_read_col(
            fptr_box,
            TSHORT,
            12,
            frow as LONGLONG,
            felem as LONGLONG,
            (nelem * 3) as LONGLONG, // 3 elements per row
            None,
            cast_slice_mut(&mut moons),
            None,
            &mut status,
        );

        // Read column 13: Temps (4E - Vector of 4 floats)
        fits_read_col(
            fptr_box,
            TFLOAT,
            13,
            frow as LONGLONG,
            felem as LONGLONG,
            (nelem * 4) as LONGLONG, // 4 elements per row
            None,
            cast_slice_mut(&mut temps),
            None,
            &mut status,
        );

        // Read column 14: Coords (2D - Vector of 2 doubles)
        fits_read_col(
            fptr_box,
            TDOUBLE,
            14,
            frow as LONGLONG,
            felem as LONGLONG,
            (nelem * 2) as LONGLONG, // 2 elements per row
            None,
            cast_slice_mut(&mut coords),
            None,
            &mut status,
        );
    }

    // Define expected values to verify against (same as written in writebintable)
    let expected_planets = ["Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn"];
    // For logical values, FITS uses 1 for T and 0 for F when stored internally
    let expected_active: [c_char; 6] = [1, 1, 1, 1, 1, 0]; // 1=T, 0=F
    // Expected flags: The varied bit patterns we wrote
    let expected_flags: [c_uchar; 6] = [
        0b10101010, // Mercury: alternating pattern
        0b11110000, // Venus: first 4 bits set
        0b00001111, // Earth: last 4 bits set
        0b11001100, // Mars: bits 0,1,4,5 set
        0b01010101, // Jupiter: alternating opposite pattern
        0b10011001, // Saturn: bits 0,3,4,7 set
    ];
    let expected_category: [c_uchar; 6] = [1, 2, 3, 4, 5, 6];
    let expected_priority: [c_short; 6] = [100, 200, 300, 150, 50, 75];
    let expected_diameter: [c_long; 6] = [4880, 12112, 12742, 6800, 143000, 121000];
    let expected_mass: [LONGLONG; 6] = [
        330110000000,
        4867500000000,
        5972400000000,
        641710000000,
        1898200000000000,
        568340000000000,
    ];
    let expected_density: [c_float; 6] = [5.1, 5.3, 5.52, 3.94, 1.33, 0.69];
    let expected_gravity: [c_double; 6] = [3.7, 8.87, 9.807, 3.71, 24.79, 10.44];
    let expected_position: [[c_float; 2]; 6] = [
        [0.39, 5.0],
        [0.72, 4.0],
        [1.00, 3.0],
        [1.52, 2.0],
        [5.20, 1.0],
        [9.58, 0.0],
    ];
    let expected_velocity: [[c_double; 2]; 6] = [
        [47.36, 0.0],
        [35.02, 1.0],
        [29.78, 2.0],
        [24.07, 3.0],
        [13.07, 4.0],
        [9.68, 5.0],
    ];
    let expected_moons: [[c_short; 3]; 6] = [
        [0, 0, 0],
        [0, 0, 0],
        [1, 0, 0],
        [0, 2, 0],
        [4, 75, 0],
        [8, 74, 0],
    ];
    let expected_temps: [[c_float; 4]; 6] = [
        [100.0, 700.0, 340.0, 150.0],
        [228.0, 773.0, 737.0, 50.0],
        [184.0, 330.0, 288.0, 40.0],
        [130.0, 308.0, 210.0, 60.0],
        [88.0, 165.0, 124.0, 20.0],
        [82.0, 134.0, 95.0, 15.0],
    ];
    let expected_coords: [[c_double; 2]; 6] = [
        [15.0, 45.0],
        [-5.0, 165.0],
        [0.0, 0.0],
        [-14.6, 175.5],
        [-22.0, 115.0],
        [0.0, 0.0],
    ];

    for ii in 0..6 {
        let name_str = CStr::from_bytes_until_nul(&name_vecs[ii])
            .unwrap()
            .to_string_lossy();

        // Assert that all read values match expected values exactly
        assert_eq!(
            name_str.trim(),
            expected_planets[ii],
            "Planet name mismatch at row {}: expected {}, got {}",
            ii + 1,
            expected_planets[ii],
            name_str.trim()
        );

        assert_eq!(
            active[ii],
            expected_active[ii],
            "Active value mismatch at row {}: expected {}, got {}",
            ii + 1,
            expected_active[ii],
            active[ii]
        );

        // Verify TBIT column data
        assert_eq!(
            flags_bytes[ii],
            expected_flags[ii],
            "Flags mismatch at row {}: expected {:#010b}, got {:#010b}",
            ii + 1,
            expected_flags[ii],
            flags_bytes[ii]
        );

        assert_eq!(
            category[ii],
            expected_category[ii],
            "Category mismatch at row {}: expected {}, got {}",
            ii + 1,
            expected_category[ii],
            category[ii]
        );

        assert_eq!(
            priority[ii],
            expected_priority[ii],
            "Priority mismatch at row {}: expected {}, got {}",
            ii + 1,
            expected_priority[ii],
            priority[ii]
        );

        assert_eq!(
            diameter[ii],
            expected_diameter[ii],
            "Diameter mismatch at row {}: expected {}, got {}",
            ii + 1,
            expected_diameter[ii],
            diameter[ii]
        );

        assert_eq!(
            mass[ii],
            expected_mass[ii],
            "Mass mismatch at row {}: expected {}, got {}",
            ii + 1,
            expected_mass[ii],
            mass[ii]
        );

        // Use approximate comparison for floating-point values
        assert!(
            (density[ii] - expected_density[ii]).abs() < 1e-6,
            "Density mismatch at row {}: expected {}, got {} (diff: {})",
            ii + 1,
            expected_density[ii],
            density[ii],
            (density[ii] - expected_density[ii]).abs()
        );

        assert!(
            (gravity[ii] - expected_gravity[ii]).abs() < 1e-10,
            "Gravity mismatch at row {}: expected {}, got {} (diff: {})",
            ii + 1,
            expected_gravity[ii],
            gravity[ii],
            (gravity[ii] - expected_gravity[ii]).abs()
        );

        // Check complex position (real and imaginary parts)
        assert!(
            (position[ii][0] - expected_position[ii][0]).abs() < 1e-6,
            "Position real part mismatch at row {}: expected {}, got {}",
            ii + 1,
            expected_position[ii][0],
            position[ii][0]
        );
        assert!(
            (position[ii][1] - expected_position[ii][1]).abs() < 1e-6,
            "Position imaginary part mismatch at row {}: expected {}, got {}",
            ii + 1,
            expected_position[ii][1],
            position[ii][1]
        );

        // Check double complex velocity (real and imaginary parts)
        assert!(
            (velocity[ii][0] - expected_velocity[ii][0]).abs() < 1e-10,
            "Velocity real part mismatch at row {}: expected {}, got {}",
            ii + 1,
            expected_velocity[ii][0],
            velocity[ii][0]
        );
        assert!(
            (velocity[ii][1] - expected_velocity[ii][1]).abs() < 1e-10,
            "Velocity imaginary part mismatch at row {}: expected {}, got {}",
            ii + 1,
            expected_velocity[ii][1],
            velocity[ii][1]
        );

        // Check all 3 elements of the moons vector
        for jj in 0..3 {
            assert_eq!(
                moons[ii][jj],
                expected_moons[ii][jj],
                "Moons vector element {} mismatch at row {}: expected {}, got {}",
                jj,
                ii + 1,
                expected_moons[ii][jj],
                moons[ii][jj]
            );
        }

        // Check all 4 elements of the temps vector
        for jj in 0..4 {
            assert!(
                (temps[ii][jj] - expected_temps[ii][jj]).abs() < 1e-6,
                "Temps vector element {} mismatch at row {}: expected {}, got {} (diff: {})",
                jj,
                ii + 1,
                expected_temps[ii][jj],
                temps[ii][jj],
                (temps[ii][jj] - expected_temps[ii][jj]).abs()
            );
        }

        // Check all 2 elements of the coords vector
        for jj in 0..2 {
            assert!(
                (coords[ii][jj] - expected_coords[ii][jj]).abs() < 1e-10,
                "Coords vector element {} mismatch at row {}: expected {}, got {} (diff: {})",
                jj,
                ii + 1,
                expected_coords[ii][jj],
                coords[ii][jj],
                (coords[ii][jj] - expected_coords[ii][jj]).abs()
            );
        }

        // Convert logical value to readable format (1=T, 0=F)
        let active_str = if active[ii] == 1 { "T" } else { "F" };

        // Format bit flags as binary (8 chars)
        let flags_str = format!("{:08b}", flags_bytes[ii]);

        // Format mass in scientific notation for readability (12 chars max)
        let mass_str = if mass[ii] > 1e12 as LONGLONG {
            format!("{:>12.1e}", mass[ii] as f64)
        } else {
            format!("{:>12}", mass[ii])
        };

        // Show position complex number (real part only for brevity)
        let pos_str = format!("{:.2}", position[ii][0]);

        // Show velocity complex number (real part only for brevity)
        let vel_str = format!("{:.1}", velocity[ii][0]);

        // Show the 3-element moons vector
        let moons_str = format!("[{},{},{}]", moons[ii][0], moons[ii][1], moons[ii][2]);

        // Show the 4-element temps vector (abbreviated)
        let temps_str = format!(
            "[{:.0},{:.0},{:.0},{:.0}]",
            temps[ii][0], temps[ii][1], temps[ii][2], temps[ii][3]
        );

        // Show the 2-element coords vector
        let coords_str = format!("[{:.1},{:.1}]", coords[ii][0], coords[ii][1]);

        println!(
            "{:>8} {:>5} {:>10} {:>5} {:>5} {:>11} {:>12} {:>10.2} {:>10.3} {:>7} {:>8} {:>11} {:>18} {:>13}",
            name_str,
            active_str,
            flags_str,
            category[ii],
            priority[ii],
            diameter[ii],
            mass_str,
            density[ii],
            gravity[ii],
            pos_str,
            vel_str,
            moons_str,
            temps_str,
            coords_str
        );
    }

    println!("\n✅ All assertion tests passed! Read values match written values exactly.");
    println!(
        "   Verified all 14 columns × 6 rows = {} individual data values across all FITS data types!",
        14 * 6
    );

    println!("\nData type examples shown and verified:");
    println!("• A (Character): Planet names");
    println!("• L (Logical): Active status (T/F)");
    println!("• X (Bit array): Flags stored as char arrays ('0'/'1' per bit)");
    println!("• B (Unsigned byte): Category numbers (1-6)");
    println!("• I (16-bit int): Priority values");
    println!("• J (32-bit int): Diameter in km");
    println!("• K (64-bit int): Mass in kg (scientific notation)");
    println!("• E (Float): Density in g/cm³");
    println!("• D (Double): Gravity in m/s²");
    println!("• C (Complex): Position AU (real and imaginary parts)");
    println!("• M (Double complex): Velocity (real and imaginary parts)");
    println!("• 3I (Vector): Moon counts [major,minor,tiny]");
    println!("• 4E (Vector): Temperature stats [min,max,mean,stddev]");
    println!("• 2D (Vector): Coordinates [lat,lon]");

    if let Some(fptr_box) = fptr
        && fits_close_file(fptr_box, &mut status) != 0
    {
        printerror(status);
    }
}

/*--------------------------------------------------------------------------*/
/// Print out cfitsio error messages and exit program
fn printerror(status: c_int) {
    if status != 0 {
        unsafe {
            fits_report_error(STDERR!(), status); /* print error report */
        }
        exit(status); /* terminate the program, returning error status */
    }
}
