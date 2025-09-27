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

    /// Generate ULONG_IMG test data (0 to large positive values)
    fn generate_ulong_test_data() -> Vec<c_ulong> {
        let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;
        let mut data = vec![0; nelements];
        for j in 0..IMAGE_HEIGHT {
            for i in 0..IMAGE_WIDTH {
                let index = (j * IMAGE_WIDTH + i) as usize;
                data[index] = ((i + j) * 1000) as c_ulong; // 0 to ~98,000 range
            }
        }
        data
    }

    #[test]
    fn test_ulong_write_ulong_read() {
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;
            let write_data = generate_ulong_test_data();

            // Write as ULONG_IMG - KNOWN ISSUE: This may fail
            let filename_cstr = CString::new(filename).unwrap();
            fits_create_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                &mut status,
            );
            fits_create_img(
                fptr.as_mut().unwrap(),
                ULONG_IMG,
                naxis as c_int,
                &naxes,
                &mut status,
            );
            fits_write_img(
                fptr.as_mut().unwrap(),
                TULONG,
                1,
                nelements as LONGLONG,
                cast_slice(&write_data),
                &mut status,
            );
            fits_close_file(fptr.take().unwrap(), &mut status);

            assert_eq!(status, 0, "Failed to write ULONG_IMG");

            // Read back as ULONG
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

                assert_eq!(status, 0, "Failed to read ULONG_IMG as TULONG");

                for i in 0..nelements {
                    assert_eq!(write_data[i], read_data[i], "Data mismatch at index {}", i);
                }
            }

            fits_close_file(fptr.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_ulong_write_byte_read_overflow() {
        // ULONG data (0-98000) will overflow BYTE range (0-255)
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;
            let write_data = generate_ulong_test_data();

            // Write as ULONG_IMG
            let filename_cstr = CString::new(filename).unwrap();
            fits_create_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                &mut status,
            );
            fits_create_img(
                fptr.as_mut().unwrap(),
                ULONG_IMG,
                naxis as c_int,
                &naxes,
                &mut status,
            );
            fits_write_img(
                fptr.as_mut().unwrap(),
                TULONG,
                1,
                nelements as LONGLONG,
                cast_slice(&write_data),
                &mut status,
            );
            fits_close_file(fptr.take().unwrap(), &mut status);

            assert_eq!(status, 0, "Failed to write ULONG_IMG");

            // Read as BYTE - should get overflow for values > 255
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

                assert_eq!(status, 412, "Expected NUM_OVERFLOW (412) for ULONG → BYTE");

                // Verify clamping behavior
                for i in 0..nelements {
                    let expected = if write_data[i] > 255 {
                        255
                    } else {
                        write_data[i] as c_uchar
                    };
                    assert_eq!(
                        read_data[i], expected,
                        "Value at index {} should be clamped: ULONG {} → BYTE {}",
                        i, write_data[i], expected
                    );
                }
            }

            fits_close_file(fptr.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_ulong_range_analysis() {
        // Verify our test data generation produces expected range
        let test_data = generate_ulong_test_data();
        let min_val = *test_data.iter().min().unwrap();
        let max_val = *test_data.iter().max().unwrap();

        println!("ULONG test data range: {} to {}", min_val, max_val);
        assert_eq!(min_val, 0, "Minimum should be 0");
        assert!(max_val <= 98000, "Maximum should be around 98,000");

        // Count values that would overflow different target types
        let byte_overflows = test_data.iter().filter(|&&x| x > 255).count();
        let short_overflows = test_data.iter().filter(|&&x| x > 32767).count();
        let ushort_overflows = test_data.iter().filter(|&&x| x > 65535).count();

        println!("Values that would overflow BYTE: {}", byte_overflows);
        println!("Values that would overflow SHORT: {}", short_overflows);
        println!("Values that would overflow USHORT: {}", ushort_overflows);

        // Most values should overflow BYTE
        assert!(
            byte_overflows > 2000,
            "Most ULONG values should overflow BYTE"
        );
        // Many values should overflow SHORT
        assert!(
            short_overflows > 1500,
            "Many ULONG values should overflow SHORT"
        );
        // Some values should overflow USHORT (561 found, so adjust expectation)
        assert!(
            ushort_overflows > 500,
            "Some ULONG values should overflow USHORT"
        );
    }

    #[test]
    fn test_ulong_theoretical_clamping_behavior() {
        // Test the clamping logic we expect CFITSIO to apply
        let test_cases = vec![
            (0u32, 0u8, 0i16, 0u16, 0i32),                   // Zero
            (100u32, 100u8, 100i16, 100u16, 100i32),         // Small value
            (255u32, 255u8, 255i16, 255u16, 255i32),         // BYTE max
            (256u32, 255u8, 256i16, 256u16, 256i32),         // BYTE overflow
            (32767u32, 255u8, 32767i16, 32767u16, 32767i32), // SHORT max
            (32768u32, 255u8, 32767i16, 32768u16, 32768i32), // SHORT overflow
            (65535u32, 255u8, 32767i16, 65535u16, 65535i32), // USHORT max
            (65536u32, 255u8, 32767i16, 65535u16, 65536i32), // USHORT overflow
        ];

        for (ulong_val, expected_byte, expected_short, expected_ushort, expected_long) in test_cases
        {
            // Test BYTE clamping
            let byte_result = if ulong_val > 255 {
                255u8
            } else {
                ulong_val as u8
            };
            assert_eq!(
                byte_result, expected_byte,
                "BYTE clamping for {}",
                ulong_val
            );

            // Test SHORT clamping
            let short_result = if ulong_val > 32767 {
                32767i16
            } else {
                ulong_val as i16
            };
            assert_eq!(
                short_result, expected_short,
                "SHORT clamping for {}",
                ulong_val
            );

            // Test USHORT clamping
            let ushort_result = if ulong_val > 65535 {
                65535u16
            } else {
                ulong_val as u16
            };
            assert_eq!(
                ushort_result, expected_ushort,
                "USHORT clamping for {}",
                ulong_val
            );

            // Test LONG clamping (ULONG max > LONG max on some systems)
            let long_result = if ulong_val > 2147483647 {
                2147483647i32
            } else {
                ulong_val as i32
            };
            assert_eq!(
                long_result, expected_long,
                "LONG clamping for {}",
                ulong_val
            );
        }

        println!("ULONG clamping logic verification complete");
    }

    #[test]
    fn test_ulong_write_short_read_overflow() {
        // ULONG data (0 to ~98,000) will overflow SHORT range (-32768 to +32767)
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;
            let write_data = generate_ulong_test_data();

            // Write as ULONG_IMG
            let filename_cstr = CString::new(filename).unwrap();
            fits_create_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                &mut status,
            );
            fits_create_img(
                fptr.as_mut().unwrap(),
                ULONG_IMG,
                naxis as c_int,
                &naxes,
                &mut status,
            );
            fits_write_img(
                fptr.as_mut().unwrap(),
                TULONG,
                1,
                nelements as LONGLONG,
                cast_slice(&write_data),
                &mut status,
            );
            fits_close_file(fptr.take().unwrap(), &mut status);

            assert_eq!(status, 0, "Failed to write ULONG_IMG");

            // Read as SHORT - should get overflow
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

                assert_eq!(status, 412, "Expected NUM_OVERFLOW (412) for ULONG → SHORT");

                // Verify clamping behavior
                for i in 0..nelements {
                    let expected = if write_data[i] > 32767 {
                        32767
                    } else {
                        write_data[i] as c_short
                    };
                    assert_eq!(
                        read_data[i], expected,
                        "Value at index {} should be clamped: ULONG {} → SHORT {} (got {})",
                        i, write_data[i], expected, read_data[i]
                    );
                }
            }

            fits_close_file(fptr.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_ulong_write_ushort_read_overflow() {
        // ULONG data (0 to ~98,000) will overflow USHORT range (0 to 65535)
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;
            let write_data = generate_ulong_test_data();

            // Write as ULONG_IMG
            let filename_cstr = CString::new(filename).unwrap();
            fits_create_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                &mut status,
            );
            fits_create_img(
                fptr.as_mut().unwrap(),
                ULONG_IMG,
                naxis as c_int,
                &naxes,
                &mut status,
            );
            fits_write_img(
                fptr.as_mut().unwrap(),
                TULONG,
                1,
                nelements as LONGLONG,
                cast_slice(&write_data),
                &mut status,
            );
            fits_close_file(fptr.take().unwrap(), &mut status);

            assert_eq!(status, 0, "Failed to write ULONG_IMG");

            // Read as USHORT - should get overflow for large values
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
                    "Expected NUM_OVERFLOW (412) for ULONG → USHORT"
                );

                // Verify clamping behavior
                for i in 0..nelements {
                    let expected = if write_data[i] > 65535 {
                        65535
                    } else {
                        write_data[i] as u16
                    };
                    assert_eq!(
                        read_data[i], expected,
                        "Value at index {} should be clamped: ULONG {} → USHORT {} (got {})",
                        i, write_data[i], expected, read_data[i]
                    );
                }
            }

            fits_close_file(fptr.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_ulong_write_long_read_success() {
        // ULONG → LONG should mostly succeed since our test data is small
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;
            let write_data = generate_ulong_test_data();

            // Write as ULONG_IMG
            let filename_cstr = CString::new(filename).unwrap();
            fits_create_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                &mut status,
            );
            fits_create_img(
                fptr.as_mut().unwrap(),
                ULONG_IMG,
                naxis as c_int,
                &naxes,
                &mut status,
            );
            fits_write_img(
                fptr.as_mut().unwrap(),
                TULONG,
                1,
                nelements as LONGLONG,
                cast_slice(&write_data),
                &mut status,
            );
            fits_close_file(fptr.take().unwrap(), &mut status);

            assert_eq!(status, 0, "Failed to write ULONG_IMG");

            // Read as LONG - should succeed for our data range
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

                assert_eq!(status, 0, "ULONG → LONG should succeed for our data range");

                for i in 0..nelements {
                    let expected = write_data[i] as c_long;
                    assert_eq!(
                        read_data[i], expected,
                        "Data mismatch at index {}: ULONG {} → LONG {} (expected {})",
                        i, write_data[i], read_data[i], expected
                    );
                }
            }

            fits_close_file(fptr.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_ulong_write_longlong_read_success() {
        // ULONG → LONGLONG should always succeed
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;
            let write_data = generate_ulong_test_data();

            // Write as ULONG_IMG
            let filename_cstr = CString::new(filename).unwrap();
            fits_create_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                &mut status,
            );
            fits_create_img(
                fptr.as_mut().unwrap(),
                ULONG_IMG,
                naxis as c_int,
                &naxes,
                &mut status,
            );
            fits_write_img(
                fptr.as_mut().unwrap(),
                TULONG,
                1,
                nelements as LONGLONG,
                cast_slice(&write_data),
                &mut status,
            );
            fits_close_file(fptr.take().unwrap(), &mut status);

            assert_eq!(status, 0, "Failed to write ULONG_IMG");

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
                    "ULONG → LONGLONG should succeed without overflow"
                );

                for i in 0..nelements {
                    let expected = write_data[i] as c_longlong;
                    assert_eq!(
                        read_data[i], expected,
                        "Data mismatch at index {}: ULONG {} → LONGLONG {} (expected {})",
                        i, write_data[i], read_data[i], expected
                    );
                }
            }

            fits_close_file(fptr.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_ulong_write_float_read_success() {
        // ULONG → FLOAT should succeed (though may lose precision for very large values)
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;
            let write_data = generate_ulong_test_data();

            // Write as ULONG_IMG
            let filename_cstr = CString::new(filename).unwrap();
            fits_create_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                &mut status,
            );
            fits_create_img(
                fptr.as_mut().unwrap(),
                ULONG_IMG,
                naxis as c_int,
                &naxes,
                &mut status,
            );
            fits_write_img(
                fptr.as_mut().unwrap(),
                TULONG,
                1,
                nelements as LONGLONG,
                cast_slice(&write_data),
                &mut status,
            );
            fits_close_file(fptr.take().unwrap(), &mut status);

            assert_eq!(status, 0, "Failed to write ULONG_IMG");

            // Read as FLOAT - should succeed
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

                assert_eq!(status, 0, "ULONG → FLOAT should succeed without overflow");

                for i in 0..nelements {
                    let expected = write_data[i] as c_float;
                    assert!(
                        floats_close_f32(expected, read_data[i]),
                        "Data mismatch at index {}: ULONG {} → FLOAT {} (expected {})",
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
    fn test_ulong_write_double_read_success() {
        // ULONG → DOUBLE should always succeed
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;
            let write_data = generate_ulong_test_data();

            // Write as ULONG_IMG
            let filename_cstr = CString::new(filename).unwrap();
            fits_create_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                &mut status,
            );
            fits_create_img(
                fptr.as_mut().unwrap(),
                ULONG_IMG,
                naxis as c_int,
                &naxes,
                &mut status,
            );
            fits_write_img(
                fptr.as_mut().unwrap(),
                TULONG,
                1,
                nelements as LONGLONG,
                cast_slice(&write_data),
                &mut status,
            );
            fits_close_file(fptr.take().unwrap(), &mut status);

            assert_eq!(status, 0, "Failed to write ULONG_IMG");

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

                assert_eq!(status, 0, "ULONG → DOUBLE should succeed without overflow");

                for i in 0..nelements {
                    let expected = write_data[i] as c_double;
                    assert!(
                        floats_close_f64(expected, read_data[i]),
                        "Data mismatch at index {}: ULONG {} → DOUBLE {} (expected {})",
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
    fn test_ulong_write_sbyte_read_overflow() {
        // ULONG data (0 to ~98,000) will overflow SBYTE range (-128 to 127)
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;
            let write_data = generate_ulong_test_data();

            // Write as ULONG_IMG
            let filename_cstr = CString::new(filename).unwrap();
            fits_create_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                &mut status,
            );
            fits_create_img(
                fptr.as_mut().unwrap(),
                ULONG_IMG,
                naxis as c_int,
                &naxes,
                &mut status,
            );
            fits_write_img(
                fptr.as_mut().unwrap(),
                TULONG,
                1,
                nelements as LONGLONG,
                cast_slice(&write_data),
                &mut status,
            );
            fits_close_file(fptr.take().unwrap(), &mut status);

            assert_eq!(status, 0, "Failed to write ULONG_IMG");

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

                assert_eq!(status, 412, "Expected NUM_OVERFLOW (412) for ULONG → SBYTE");

                // Verify clamping behavior
                for i in 0..nelements {
                    let expected = if write_data[i] > 127 {
                        127
                    } else {
                        write_data[i] as i8
                    };
                    assert_eq!(
                        read_data[i], expected,
                        "Value at index {} should be clamped: ULONG {} → SBYTE {} (got {})",
                        i, write_data[i], expected, read_data[i]
                    );
                }
            }

            fits_close_file(fptr.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_ulong_write_ulonglong_read_success() {
        // ULONG → ULONGLONG should always succeed
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;
            let write_data = generate_ulong_test_data();

            // Write as ULONG_IMG
            let filename_cstr = CString::new(filename).unwrap();
            fits_create_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                &mut status,
            );
            fits_create_img(
                fptr.as_mut().unwrap(),
                ULONG_IMG,
                naxis as c_int,
                &naxes,
                &mut status,
            );
            fits_write_img(
                fptr.as_mut().unwrap(),
                TULONG,
                1,
                nelements as LONGLONG,
                cast_slice(&write_data),
                &mut status,
            );
            fits_close_file(fptr.take().unwrap(), &mut status);

            assert_eq!(status, 0, "Failed to write ULONG_IMG");

            // Read as ULONGLONG - should succeed
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
                    status, 0,
                    "ULONG → ULONGLONG should succeed without overflow"
                );

                for i in 0..nelements {
                    let expected = write_data[i] as c_ulonglong;
                    assert_eq!(
                        read_data[i], expected,
                        "Data mismatch at index {}: ULONG {} → ULONGLONG {} (expected {})",
                        i, write_data[i], read_data[i], expected
                    );
                }
            }

            fits_close_file(fptr.take().unwrap(), &mut status);
        });
    }
}
