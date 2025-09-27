#[cfg(test)]
mod tests {
    use bytemuck::{cast_slice, cast_slice_mut};
    use libc::{
        c_char, c_double, c_float, c_int, c_long, c_longlong, c_short, c_uchar, c_ulong,
        c_ulonglong,
    };
    use rsfitsio::aliases::rust_api::*;
    use rsfitsio::fitsio::{
        BYTE_IMG, DOUBLE_IMG, FLOAT_IMG, LONG_IMG, LONGLONG, LONGLONG_IMG, READWRITE, SBYTE_IMG,
        SHORT_IMG, TBYTE, TDOUBLE, TFLOAT, TLONG, TLONGLONG, TSBYTE, TSHORT, TULONG, TULONGLONG,
        TUSHORT, ULONG_IMG, ULONGLONG_IMG, USHORT_IMG, fitsfile,
    };
    use rsfitsio::helpers::testhelpers::{floats_close_f32, floats_close_f64, with_temp_file};
    use std::ffi::CString;

    const IMAGE_WIDTH: c_long = 50;
    const IMAGE_HEIGHT: c_long = 50;

    /// Generate SHORT_IMG test data (includes negatives)
    fn generate_short_test_data() -> Vec<c_short> {
        let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;
        let mut data = vec![0; nelements];
        for j in 0..IMAGE_HEIGHT {
            for i in 0..IMAGE_WIDTH {
                let index = (j * IMAGE_WIDTH + i) as usize;
                data[index] = ((i - 25) * (j - 25)) as c_short;
            }
        }
        data
    }

    #[test]
    fn test_short_write_short_read() {
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;
            let write_data = generate_short_test_data();

            // Write as SHORT_IMG
            let filename_cstr = CString::new(filename).unwrap();
            fits_create_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                &mut status,
            );
            fits_create_img(
                fptr.as_mut().unwrap(),
                SHORT_IMG,
                naxis as c_int,
                &naxes,
                &mut status,
            );
            fits_write_img(
                fptr.as_mut().unwrap(),
                TSHORT,
                1,
                nelements as LONGLONG,
                cast_slice(&write_data),
                &mut status,
            );
            fits_close_file(fptr.take().unwrap(), &mut status);

            assert_eq!(status, 0, "Failed to write SHORT_IMG");

            // Read back as SHORT
            fptr = None;
            fits_open_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                READWRITE,
                &mut status,
            );

            if let Some(ref mut fptr_box) = fptr {
                let mut read_data: Vec<c_short> = vec![0; nelements];
                let mut anynull: c_int = 0;

                fits_read_img(
                    fptr_box,
                    TSHORT,
                    1,
                    nelements as LONGLONG,
                    None,
                    cast_slice_mut(&mut read_data),
                    Some(&mut anynull),
                    &mut status,
                );

                assert_eq!(status, 0, "Failed to read SHORT_IMG as TSHORT");

                for i in 0..nelements {
                    assert_eq!(write_data[i], read_data[i], "Data mismatch at index {}", i);
                }
            }

            fits_close_file(fptr.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_short_write_byte_read_overflow() {
        // This test verifies that CFITSIO correctly returns NUM_OVERFLOW (412)
        // when SHORT data contains values outside BYTE range (0-255)
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;
            let write_data = generate_short_test_data(); // Contains negative values

            // Write as SHORT_IMG
            let filename_cstr = CString::new(filename).unwrap();
            fits_create_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                &mut status,
            );
            fits_create_img(
                fptr.as_mut().unwrap(),
                SHORT_IMG,
                naxis as c_int,
                &naxes,
                &mut status,
            );
            fits_write_img(
                fptr.as_mut().unwrap(),
                TSHORT,
                1,
                nelements as LONGLONG,
                cast_slice(&write_data),
                &mut status,
            );
            fits_close_file(fptr.take().unwrap(), &mut status);

            assert_eq!(status, 0, "Failed to write SHORT_IMG");

            // Attempt to read as BYTE - should fail with NUM_OVERFLOW (412)
            fptr = None;
            fits_open_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                READWRITE,
                &mut status,
            );

            if let Some(ref mut fptr_box) = fptr {
                let mut read_data: Vec<c_uchar> = vec![0; nelements];
                let mut anynull: c_int = 0;

                fits_read_img(
                    fptr_box,
                    TBYTE,
                    1,
                    nelements as LONGLONG,
                    None,
                    cast_slice_mut(&mut read_data),
                    Some(&mut anynull),
                    &mut status,
                );

                // Verify that NUM_OVERFLOW (412) is returned due to negative SHORT values
                assert_eq!(
                    status, 412,
                    "Expected NUM_OVERFLOW (412) when reading SHORT with negatives as BYTE"
                );

                // Verify that CFITSIO clamps overflowed values while reporting NUM_OVERFLOW
                for i in 0..nelements {
                    let expected = if write_data[i] < 0 {
                        0 // Negative values clamped to 0
                    } else if write_data[i] > 255 {
                        255 // Values > 255 clamped to 255 (max BYTE)
                    } else {
                        write_data[i] as c_uchar // Values in range converted normally
                    };

                    assert_eq!(
                        read_data[i], expected,
                        "Value at index {} should be clamped: SHORT {} → BYTE {} (got {})",
                        i, write_data[i], expected, read_data[i]
                    );
                }
            }

            fits_close_file(fptr.take().unwrap(), &mut 0);
        });
    }

    #[test]
    fn test_short_write_byte_read_success() {
        // Test successful SHORT → BYTE conversion with data in valid range (0-255)
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;

            // Create SHORT data that fits in BYTE range (0-255)
            let mut write_data: Vec<c_short> = vec![0; nelements];
            for j in 0..IMAGE_HEIGHT {
                for i in 0..IMAGE_WIDTH {
                    let index = (j * IMAGE_WIDTH + i) as usize;
                    write_data[index] = ((i + j) % 256) as c_short; // 0-255 range
                }
            }

            // Write as SHORT_IMG
            let filename_cstr = CString::new(filename).unwrap();
            fits_create_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                &mut status,
            );
            fits_create_img(
                fptr.as_mut().unwrap(),
                SHORT_IMG,
                naxis as c_int,
                &naxes,
                &mut status,
            );
            fits_write_img(
                fptr.as_mut().unwrap(),
                TSHORT,
                1,
                nelements as LONGLONG,
                cast_slice(&write_data),
                &mut status,
            );
            fits_close_file(fptr.take().unwrap(), &mut status);

            assert_eq!(status, 0, "Failed to write SHORT_IMG");

            // Read back as BYTE - should succeed since all values are 0-255
            fptr = None;
            fits_open_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                READWRITE,
                &mut status,
            );

            if let Some(ref mut fptr_box) = fptr {
                let mut read_data: Vec<c_uchar> = vec![0; nelements];
                let mut anynull: c_int = 0;

                fits_read_img(
                    fptr_box,
                    TBYTE,
                    1,
                    nelements as LONGLONG,
                    None,
                    cast_slice_mut(&mut read_data),
                    Some(&mut anynull),
                    &mut status,
                );

                assert_eq!(
                    status, 0,
                    "Should succeed reading SHORT in BYTE range as TBYTE"
                );

                // Verify data matches
                for i in 0..nelements {
                    assert_eq!(
                        write_data[i] as c_uchar, read_data[i],
                        "Data mismatch at index {}: SHORT {}, BYTE {}",
                        i, write_data[i], read_data[i]
                    );
                }
            }

            fits_close_file(fptr.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_short_write_byte_read_detailed_overflow() {
        // Create a mixed dataset with some values in BYTE range and some out of range
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [10, 10]; // Smaller for easier verification
            let nelements = 100;

            // Create mixed SHORT data: some in BYTE range (0-255), some negative, some > 255
            let mut write_data: Vec<c_short> = vec![0; nelements];
            write_data[0] = -1; // Underflow
            write_data[1] = -100; // Underflow  
            write_data[2] = 0; // Valid
            write_data[3] = 127; // Valid
            write_data[4] = 255; // Valid (max BYTE)
            write_data[5] = 256; // Overflow
            write_data[6] = 1000; // Overflow
            // Fill rest with valid values
            for i in 7..nelements {
                write_data[i] = (i % 200) as c_short; // 0-199 range
            }

            // Write as SHORT_IMG
            let filename_cstr = CString::new(filename).unwrap();
            fits_create_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                &mut status,
            );
            fits_create_img(
                fptr.as_mut().unwrap(),
                SHORT_IMG,
                naxis as c_int,
                &naxes,
                &mut status,
            );
            fits_write_img(
                fptr.as_mut().unwrap(),
                TSHORT,
                1,
                nelements as LONGLONG,
                cast_slice(&write_data),
                &mut status,
            );
            fits_close_file(fptr.take().unwrap(), &mut status);

            assert_eq!(status, 0, "Failed to write SHORT_IMG");

            // Attempt to read as BYTE
            fptr = None;
            fits_open_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                READWRITE,
                &mut status,
            );

            if let Some(ref mut fptr_box) = fptr {
                let mut read_data: Vec<c_uchar> = vec![99; nelements]; // Initialize to non-zero to detect changes
                let mut anynull: c_int = 0;

                fits_read_img(
                    fptr_box,
                    TBYTE,
                    1,
                    nelements as LONGLONG,
                    None,
                    cast_slice_mut(&mut read_data),
                    Some(&mut anynull),
                    &mut status,
                );

                // Should get NUM_OVERFLOW due to out-of-range values
                assert_eq!(
                    status, 412,
                    "Expected NUM_OVERFLOW (412) for mixed-range data"
                );

                // Check specific values - CFITSIO clamps overflowed values but reports NUM_OVERFLOW
                assert_eq!(
                    read_data[0], 0,
                    "Underflow value (-1) should be clamped to 0"
                );
                assert_eq!(
                    read_data[1], 0,
                    "Underflow value (-100) should be clamped to 0"
                );
                assert_eq!(
                    read_data[5], 255,
                    "Overflow value (256) should be clamped to 255"
                );
                assert_eq!(
                    read_data[6], 255,
                    "Overflow value (1000) should be clamped to 255"
                );

                // Check that valid values are properly converted
                assert_eq!(
                    read_data[2], 0,
                    "Valid value (0) should be converted correctly"
                );
                assert_eq!(
                    read_data[3], 127,
                    "Valid value (127) should be converted correctly"
                );
                assert_eq!(
                    read_data[4], 255,
                    "Valid value (255) should be converted correctly"
                );

                // Check some of the other valid values
                for i in 7..nelements {
                    let expected = write_data[i] as c_uchar; // Should be 0-199 range
                    assert_eq!(
                        read_data[i], expected,
                        "Valid value at index {} should be converted: SHORT {} → BYTE {}",
                        i, write_data[i], expected
                    );
                }

                println!(
                    "Read data verification complete - overflowed values clamped, NUM_OVERFLOW reported"
                );
            }

            fits_close_file(fptr.take().unwrap(), &mut 0);
        });
    }

    #[test]
    fn test_short_write_float_read() {
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;
            let write_data = generate_short_test_data();

            // Write as SHORT_IMG
            let filename_cstr = CString::new(filename).unwrap();
            fits_create_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                &mut status,
            );
            fits_create_img(
                fptr.as_mut().unwrap(),
                SHORT_IMG,
                naxis as c_int,
                &naxes,
                &mut status,
            );
            fits_write_img(
                fptr.as_mut().unwrap(),
                TSHORT,
                1,
                nelements as LONGLONG,
                cast_slice(&write_data),
                &mut status,
            );
            fits_close_file(fptr.take().unwrap(), &mut status);

            assert_eq!(status, 0, "Failed to write SHORT_IMG");

            // Read back as FLOAT
            fptr = None;
            fits_open_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                READWRITE,
                &mut status,
            );

            if let Some(ref mut fptr_box) = fptr {
                let mut read_data: Vec<c_float> = vec![0.0; nelements];
                let mut anynull: c_int = 0;

                fits_read_img(
                    fptr_box,
                    TFLOAT,
                    1,
                    nelements as LONGLONG,
                    None,
                    cast_slice_mut(&mut read_data),
                    Some(&mut anynull),
                    &mut status,
                );

                assert_eq!(status, 0, "Failed to read SHORT_IMG as TFLOAT");

                for i in 0..nelements {
                    let expected = write_data[i] as c_float;
                    assert!(
                        floats_close_f32(expected, read_data[i]),
                        "Data mismatch at index {}: expected {}, got {}",
                        i,
                        expected,
                        read_data[i]
                    );
                }
            }

            fits_close_file(fptr.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_short_write_ushort_read_success() {
        // SHORT → USHORT should succeed for positive values, clamp negatives
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;
            let write_data = generate_short_test_data();

            // Write as SHORT_IMG
            let filename_cstr = CString::new(filename).unwrap();
            fits_create_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                &mut status,
            );
            fits_create_img(
                fptr.as_mut().unwrap(),
                SHORT_IMG,
                naxis as c_int,
                &naxes,
                &mut status,
            );
            fits_write_img(
                fptr.as_mut().unwrap(),
                TSHORT,
                1,
                nelements as LONGLONG,
                cast_slice(&write_data),
                &mut status,
            );
            fits_close_file(fptr.take().unwrap(), &mut status);

            assert_eq!(status, 0, "Failed to write SHORT_IMG");

            // Read as USHORT - should get overflow for negative values
            fptr = None;
            fits_open_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                READWRITE,
                &mut status,
            );

            if let Some(ref mut fptr_box) = fptr {
                let mut read_data: Vec<u16> = vec![0; nelements];
                let mut anynull: c_int = 0;

                fits_read_img(
                    fptr_box,
                    TUSHORT,
                    1,
                    nelements as LONGLONG,
                    None,
                    cast_slice_mut(&mut read_data),
                    Some(&mut anynull),
                    &mut status,
                );

                assert_eq!(
                    status, 412,
                    "Expected NUM_OVERFLOW (412) for SHORT → USHORT with negatives"
                );

                // Verify clamping behavior
                for i in 0..nelements {
                    let expected = if write_data[i] < 0 {
                        0
                    } else {
                        write_data[i] as u16
                    };
                    assert_eq!(
                        read_data[i], expected,
                        "Value at index {} should be clamped: SHORT {} → USHORT {}",
                        i, write_data[i], expected
                    );
                }
            }

            fits_close_file(fptr.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_short_write_long_read_success() {
        // SHORT → LONG should always succeed (SHORT range fits in LONG)
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;
            let write_data = generate_short_test_data();

            // Write as SHORT_IMG
            let filename_cstr = CString::new(filename).unwrap();
            fits_create_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                &mut status,
            );
            fits_create_img(
                fptr.as_mut().unwrap(),
                SHORT_IMG,
                naxis as c_int,
                &naxes,
                &mut status,
            );
            fits_write_img(
                fptr.as_mut().unwrap(),
                TSHORT,
                1,
                nelements as LONGLONG,
                cast_slice(&write_data),
                &mut status,
            );
            fits_close_file(fptr.take().unwrap(), &mut status);

            assert_eq!(status, 0, "Failed to write SHORT_IMG");

            // Read as LONG - should succeed
            fptr = None;
            fits_open_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                READWRITE,
                &mut status,
            );

            if let Some(ref mut fptr_box) = fptr {
                let mut read_data: Vec<c_long> = vec![0; nelements];
                let mut anynull: c_int = 0;

                fits_read_img(
                    fptr_box,
                    TLONG,
                    1,
                    nelements as LONGLONG,
                    None,
                    cast_slice_mut(&mut read_data),
                    Some(&mut anynull),
                    &mut status,
                );

                assert_eq!(status, 0, "SHORT → LONG should succeed without overflow");

                for i in 0..nelements {
                    assert_eq!(
                        write_data[i] as c_long, read_data[i],
                        "Data mismatch at index {}: SHORT {} → LONG {}",
                        i, write_data[i], read_data[i]
                    );
                }
            }

            fits_close_file(fptr.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_short_write_longlong_read_success() {
        // SHORT → LONGLONG should always succeed (SHORT range fits in LONGLONG)
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;
            let write_data = generate_short_test_data();

            // Write as SHORT_IMG
            let filename_cstr = CString::new(filename).unwrap();
            fits_create_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                &mut status,
            );
            fits_create_img(
                fptr.as_mut().unwrap(),
                SHORT_IMG,
                naxis as c_int,
                &naxes,
                &mut status,
            );
            fits_write_img(
                fptr.as_mut().unwrap(),
                TSHORT,
                1,
                nelements as LONGLONG,
                cast_slice(&write_data),
                &mut status,
            );
            fits_close_file(fptr.take().unwrap(), &mut status);

            assert_eq!(status, 0, "Failed to write SHORT_IMG");

            // Read as LONGLONG - should succeed
            fptr = None;
            fits_open_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                READWRITE,
                &mut status,
            );

            if let Some(ref mut fptr_box) = fptr {
                let mut read_data: Vec<c_longlong> = vec![0; nelements];
                let mut anynull: c_int = 0;

                fits_read_img(
                    fptr_box,
                    TLONGLONG,
                    1,
                    nelements as LONGLONG,
                    None,
                    cast_slice_mut(&mut read_data),
                    Some(&mut anynull),
                    &mut status,
                );

                assert_eq!(
                    status, 0,
                    "SHORT → LONGLONG should succeed without overflow"
                );

                for i in 0..nelements {
                    assert_eq!(
                        write_data[i] as c_longlong, read_data[i],
                        "Data mismatch at index {}: SHORT {} → LONGLONG {}",
                        i, write_data[i], read_data[i]
                    );
                }
            }

            fits_close_file(fptr.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_short_write_double_read_success() {
        // SHORT → DOUBLE should always succeed
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;
            let write_data = generate_short_test_data();

            // Write as SHORT_IMG
            let filename_cstr = CString::new(filename).unwrap();
            fits_create_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                &mut status,
            );
            fits_create_img(
                fptr.as_mut().unwrap(),
                SHORT_IMG,
                naxis as c_int,
                &naxes,
                &mut status,
            );
            fits_write_img(
                fptr.as_mut().unwrap(),
                TSHORT,
                1,
                nelements as LONGLONG,
                cast_slice(&write_data),
                &mut status,
            );
            fits_close_file(fptr.take().unwrap(), &mut status);

            assert_eq!(status, 0, "Failed to write SHORT_IMG");

            // Read as DOUBLE - should succeed
            fptr = None;
            fits_open_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                READWRITE,
                &mut status,
            );

            if let Some(ref mut fptr_box) = fptr {
                let mut read_data: Vec<c_double> = vec![0.0; nelements];
                let mut anynull: c_int = 0;

                fits_read_img(
                    fptr_box,
                    TDOUBLE,
                    1,
                    nelements as LONGLONG,
                    None,
                    cast_slice_mut(&mut read_data),
                    Some(&mut anynull),
                    &mut status,
                );

                assert_eq!(status, 0, "SHORT → DOUBLE should succeed without overflow");

                for i in 0..nelements {
                    let expected = write_data[i] as c_double;
                    assert!(
                        floats_close_f64(expected, read_data[i]),
                        "Data mismatch at index {}: SHORT {} → DOUBLE {} (expected {})",
                        i,
                        write_data[i],
                        read_data[i],
                        expected
                    );
                }
            }

            fits_close_file(fptr.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_short_write_sbyte_read_overflow() {
        // SHORT → SBYTE should overflow for values outside -128 to 127
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;
            let write_data = generate_short_test_data(); // Range: ~-625 to ~576

            // Write as SHORT_IMG
            let filename_cstr = CString::new(filename).unwrap();
            fits_create_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                &mut status,
            );
            fits_create_img(
                fptr.as_mut().unwrap(),
                SHORT_IMG,
                naxis as c_int,
                &naxes,
                &mut status,
            );
            fits_write_img(
                fptr.as_mut().unwrap(),
                TSHORT,
                1,
                nelements as LONGLONG,
                cast_slice(&write_data),
                &mut status,
            );
            fits_close_file(fptr.take().unwrap(), &mut status);

            assert_eq!(status, 0, "Failed to write SHORT_IMG");

            // Read as SBYTE - should get overflow
            fptr = None;
            fits_open_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                READWRITE,
                &mut status,
            );

            if let Some(ref mut fptr_box) = fptr {
                let mut read_data: Vec<i8> = vec![0; nelements];
                let mut anynull: c_int = 0;

                fits_read_img(
                    fptr_box,
                    TSBYTE,
                    1,
                    nelements as LONGLONG,
                    None,
                    cast_slice_mut(&mut read_data),
                    Some(&mut anynull),
                    &mut status,
                );

                assert_eq!(status, 412, "Expected NUM_OVERFLOW (412) for SHORT → SBYTE");

                // Verify clamping behavior
                for i in 0..nelements {
                    let expected = if write_data[i] < -128 {
                        -128
                    } else if write_data[i] > 127 {
                        127
                    } else {
                        write_data[i] as i8
                    };
                    assert_eq!(
                        read_data[i], expected,
                        "Value at index {} should be clamped: SHORT {} → SBYTE {}",
                        i, write_data[i], expected
                    );
                }
            }

            fits_close_file(fptr.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_short_write_ulong_read_overflow() {
        // SHORT → ULONG should clamp negative values to 0
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;
            let write_data = generate_short_test_data();

            // Write as SHORT_IMG
            let filename_cstr = CString::new(filename).unwrap();
            fits_create_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                &mut status,
            );
            fits_create_img(
                fptr.as_mut().unwrap(),
                SHORT_IMG,
                naxis as c_int,
                &naxes,
                &mut status,
            );
            fits_write_img(
                fptr.as_mut().unwrap(),
                TSHORT,
                1,
                nelements as LONGLONG,
                cast_slice(&write_data),
                &mut status,
            );
            fits_close_file(fptr.take().unwrap(), &mut status);

            assert_eq!(status, 0, "Failed to write SHORT_IMG");

            // Read as ULONG - should get overflow for negative values
            fptr = None;
            fits_open_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                READWRITE,
                &mut status,
            );

            if let Some(ref mut fptr_box) = fptr {
                let mut read_data: Vec<c_ulong> = vec![0; nelements];
                let mut anynull: c_int = 0;

                fits_read_img(
                    fptr_box,
                    TULONG,
                    1,
                    nelements as LONGLONG,
                    None,
                    cast_slice_mut(&mut read_data),
                    Some(&mut anynull),
                    &mut status,
                );

                assert_eq!(
                    status, 412,
                    "Expected NUM_OVERFLOW (412) for SHORT → ULONG with negatives"
                );

                // Verify clamping behavior
                for i in 0..nelements {
                    let expected = if write_data[i] < 0 {
                        0
                    } else {
                        write_data[i] as c_ulong
                    };
                    assert_eq!(
                        read_data[i], expected,
                        "Value at index {} should be clamped: SHORT {} → ULONG {}",
                        i, write_data[i], expected
                    );
                }
            }

            fits_close_file(fptr.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_short_write_ulonglong_read_overflow() {
        // SHORT → ULONGLONG should clamp negative values to 0
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;
            let write_data = generate_short_test_data();

            // Write as SHORT_IMG
            let filename_cstr = CString::new(filename).unwrap();
            fits_create_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                &mut status,
            );
            fits_create_img(
                fptr.as_mut().unwrap(),
                SHORT_IMG,
                naxis as c_int,
                &naxes,
                &mut status,
            );
            fits_write_img(
                fptr.as_mut().unwrap(),
                TSHORT,
                1,
                nelements as LONGLONG,
                cast_slice(&write_data),
                &mut status,
            );
            fits_close_file(fptr.take().unwrap(), &mut status);

            assert_eq!(status, 0, "Failed to write SHORT_IMG");

            // Read as ULONGLONG - should get overflow for negative values
            fptr = None;
            fits_open_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                READWRITE,
                &mut status,
            );

            if let Some(ref mut fptr_box) = fptr {
                let mut read_data: Vec<c_ulonglong> = vec![0; nelements];
                let mut anynull: c_int = 0;

                fits_read_img(
                    fptr_box,
                    TULONGLONG,
                    1,
                    nelements as LONGLONG,
                    None,
                    cast_slice_mut(&mut read_data),
                    Some(&mut anynull),
                    &mut status,
                );

                assert_eq!(
                    status, 412,
                    "Expected NUM_OVERFLOW (412) for SHORT → ULONGLONG with negatives"
                );

                // Verify clamping behavior
                for i in 0..nelements {
                    let expected = if write_data[i] < 0 {
                        0
                    } else {
                        write_data[i] as c_ulonglong
                    };
                    assert_eq!(
                        read_data[i], expected,
                        "Value at index {} should be clamped: SHORT {} → ULONGLONG {}",
                        i, write_data[i], expected
                    );
                }
            }

            fits_close_file(fptr.take().unwrap(), &mut status);
        });
    }
}
