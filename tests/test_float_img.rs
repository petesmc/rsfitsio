#[cfg(test)]
mod tests {
    use bytemuck::{cast_slice, cast_slice_mut};
    use libc::{
        c_double, c_float, c_int, c_long, c_longlong, c_short, c_uchar, c_ulong, c_ulonglong,
    };
    use rsfitsio::aliases::rust_api::*;
    use rsfitsio::fitsio::{
        FLOAT_IMG, LONGLONG, READWRITE, TBYTE, TDOUBLE, TFLOAT, TLONG, TLONGLONG, TSBYTE, TSHORT,
        TULONG, TULONGLONG, TUSHORT, fitsfile,
    };
    use rsfitsio::helpers::testhelpers::{floats_close_f32, floats_close_f64, with_temp_file};
    use std::ffi::CString;

    const IMAGE_WIDTH: c_long = 50;
    const IMAGE_HEIGHT: c_long = 50;

    /// Generate FLOAT_IMG test data (includes negatives and fractional values)
    fn generate_float_test_data() -> Vec<c_float> {
        let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;
        let mut data = vec![0.0; nelements];
        for j in 0..IMAGE_HEIGHT {
            for i in 0..IMAGE_WIDTH {
                let index = (j * IMAGE_WIDTH + i) as usize;
                data[index] = (i - 25) as c_float * 10.0 + (j - 25) as c_float * 1.0;
            }
        }
        data
    }

    #[test]
    fn test_float_write_float_read() {
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;
            let write_data = generate_float_test_data();

            // Write as FLOAT_IMG
            let filename_cstr = CString::new(filename).unwrap();
            fits_create_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                &mut status,
            );
            fits_create_img(
                fptr.as_mut().unwrap(),
                FLOAT_IMG,
                naxis as c_int,
                &naxes,
                &mut status,
            );
            fits_write_img(
                fptr.as_mut().unwrap(),
                TFLOAT,
                1,
                nelements as LONGLONG,
                cast_slice(&write_data),
                &mut status,
            );
            fits_close_file(fptr.take().unwrap(), &mut status);

            assert_eq!(status, 0, "Failed to write FLOAT_IMG");

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

                assert_eq!(status, 0, "Failed to read FLOAT_IMG as TFLOAT");

                for i in 0..nelements {
                    assert!(
                        floats_close_f32(write_data[i], read_data[i]),
                        "Data mismatch at index {}: expected {}, got {}",
                        i,
                        write_data[i],
                        read_data[i]
                    );
                }
            }

            fits_close_file(fptr.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_float_write_byte_read_overflow() {
        // FLOAT data with negatives and fractional values will overflow/truncate for BYTE
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [10, 10];
            let nelements = 100;

            // Create FLOAT data with specific test cases
            let mut write_data: Vec<c_float> = vec![0.0; nelements];
            write_data[0] = -1.5; // Negative
            write_data[1] = 0.0; // Zero
            write_data[2] = 0.7; // Fractional < 1
            write_data[3] = 1.8; // Fractional > 1
            write_data[4] = 127.9; // Truncates to 127
            write_data[5] = 255.0; // Exact max BYTE
            write_data[6] = 255.9; // Truncates to 255
            write_data[7] = 256.0; // Overflow
            write_data[8] = 1000.5; // Large overflow
            for i in 9..nelements {
                write_data[i] = (i as c_float - 50.0) * 2.5; // Mix of values
            }

            // Write as FLOAT_IMG
            let filename_cstr = CString::new(filename).unwrap();
            fits_create_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                &mut status,
            );
            fits_create_img(
                fptr.as_mut().unwrap(),
                FLOAT_IMG,
                naxis as c_int,
                &naxes,
                &mut status,
            );
            fits_write_img(
                fptr.as_mut().unwrap(),
                TFLOAT,
                1,
                nelements as LONGLONG,
                cast_slice(&write_data),
                &mut status,
            );
            fits_close_file(fptr.take().unwrap(), &mut status);

            assert_eq!(status, 0, "Failed to write FLOAT_IMG");

            // Read as BYTE - should get overflow due to negatives and > 255 values
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

                assert_eq!(status, 412, "Expected NUM_OVERFLOW (412) for FLOAT → BYTE");

                // Verify specific conversions (truncation + clamping)
                assert_eq!(read_data[0], 0, "Negative should clamp to 0");
                assert_eq!(read_data[1], 0, "Zero should convert to 0");
                assert_eq!(read_data[2], 0, "0.7 should truncate to 0");
                assert_eq!(read_data[3], 1, "1.8 should truncate to 1");
                assert_eq!(read_data[4], 127, "127.9 should truncate to 127");
                assert_eq!(read_data[5], 255, "255.0 should convert to 255");
                assert_eq!(read_data[6], 255, "255.9 should truncate to 255");
                assert_eq!(read_data[7], 255, "256.0 should clamp to 255");
                assert_eq!(read_data[8], 255, "1000.5 should clamp to 255");
            }

            fits_close_file(fptr.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_float_write_short_read_overflow() {
        // Create FLOAT data that tests SHORT range boundaries
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [10, 10];
            let nelements = 100;

            let mut write_data: Vec<c_float> = vec![0.0; nelements];
            write_data[0] = -32769.5; // Underflow SHORT
            write_data[1] = -32768.0; // Min SHORT
            write_data[2] = -100.7; // Negative with fraction
            write_data[3] = 0.0; // Zero
            write_data[4] = 100.9; // Positive with fraction
            write_data[5] = 32767.0; // Max SHORT
            write_data[6] = 32767.9; // Just under overflow (truncates to 32767)
            write_data[7] = 32768.0; // Overflow SHORT
            write_data[8] = 50000.5; // Large overflow
            for i in 9..nelements {
                write_data[i] = (i as c_float - 50.0) * 100.0;
            }

            // Write as FLOAT_IMG
            let filename_cstr = CString::new(filename).unwrap();
            fits_create_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                &mut status,
            );
            fits_create_img(
                fptr.as_mut().unwrap(),
                FLOAT_IMG,
                naxis as c_int,
                &naxes,
                &mut status,
            );
            fits_write_img(
                fptr.as_mut().unwrap(),
                TFLOAT,
                1,
                nelements as LONGLONG,
                cast_slice(&write_data),
                &mut status,
            );
            fits_close_file(fptr.take().unwrap(), &mut status);

            assert_eq!(status, 0, "Failed to write FLOAT_IMG");

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

                assert_eq!(status, 412, "Expected NUM_OVERFLOW (412) for FLOAT → SHORT");

                // Verify specific conversions (truncation + clamping)
                assert_eq!(read_data[0], -32768, "Underflow should clamp to SHORT min");
                assert_eq!(read_data[1], -32768, "Min SHORT should convert exactly");
                assert_eq!(read_data[2], -100, "Negative should truncate to -100");
                assert_eq!(read_data[3], 0, "Zero should convert exactly");
                assert_eq!(read_data[4], 100, "Positive should truncate to 100");
                assert_eq!(read_data[5], 32767, "Max SHORT should convert exactly");
                assert_eq!(read_data[6], 32767, "32767.9 should truncate to 32767");
                assert_eq!(read_data[7], 32767, "Overflow should clamp to SHORT max");
                assert_eq!(
                    read_data[8], 32767,
                    "Large overflow should clamp to SHORT max"
                );
            }

            fits_close_file(fptr.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_float_write_double_read_success() {
        // FLOAT → DOUBLE should always succeed
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;
            let write_data = generate_float_test_data();

            // Write as FLOAT_IMG
            let filename_cstr = CString::new(filename).unwrap();
            fits_create_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                &mut status,
            );
            fits_create_img(
                fptr.as_mut().unwrap(),
                FLOAT_IMG,
                naxis as c_int,
                &naxes,
                &mut status,
            );
            fits_write_img(
                fptr.as_mut().unwrap(),
                TFLOAT,
                1,
                nelements as LONGLONG,
                cast_slice(&write_data),
                &mut status,
            );
            fits_close_file(fptr.take().unwrap(), &mut status);

            assert_eq!(status, 0, "Failed to write FLOAT_IMG");

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

                assert_eq!(status, 0, "FLOAT → DOUBLE should succeed without overflow");

                for i in 0..nelements {
                    let expected = write_data[i] as c_double;
                    assert!(
                        floats_close_f64(expected, read_data[i]),
                        "Data mismatch at index {}: FLOAT {} → DOUBLE {} (expected {})",
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
    fn test_float_write_ushort_read_overflow() {
        // Create FLOAT data that tests USHORT range boundaries
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [10, 10];
            let nelements = 100;

            let mut write_data: Vec<c_float> = vec![0.0; nelements];
            write_data[0] = -1.5; // Negative (clamps to 0)
            write_data[1] = 0.0; // Zero
            write_data[2] = 100.7; // Positive with fraction
            write_data[3] = 65535.0; // Max USHORT
            write_data[4] = 65535.9; // Just under overflow (truncates to 65535)
            write_data[5] = 65536.0; // Overflow USHORT
            write_data[6] = 100000.5; // Large overflow
            for i in 7..nelements {
                write_data[i] = (i as c_float - 50.0) * 500.0;
            }

            // Write as FLOAT_IMG
            let filename_cstr = CString::new(filename).unwrap();
            fits_create_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                &mut status,
            );
            fits_create_img(
                fptr.as_mut().unwrap(),
                FLOAT_IMG,
                naxis as c_int,
                &naxes,
                &mut status,
            );
            fits_write_img(
                fptr.as_mut().unwrap(),
                TFLOAT,
                1,
                nelements as LONGLONG,
                cast_slice(&write_data),
                &mut status,
            );
            fits_close_file(fptr.take().unwrap(), &mut status);

            assert_eq!(status, 0, "Failed to write FLOAT_IMG");

            // Read as USHORT - should get overflow
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
                    "Expected NUM_OVERFLOW (412) for FLOAT → USHORT"
                );

                // Verify specific conversions (truncation + clamping)
                assert_eq!(read_data[0], 0, "Negative should clamp to 0");
                assert_eq!(read_data[1], 0, "Zero should convert to 0");
                assert_eq!(read_data[2], 100, "100.7 should truncate to 100");
                assert_eq!(read_data[3], 65535, "65535.0 should convert to 65535");
                assert_eq!(read_data[4], 65535, "65535.9 should truncate to 65535");
                assert_eq!(read_data[5], 65535, "65536.0 should clamp to 65535");
                assert_eq!(read_data[6], 65535, "100000.5 should clamp to 65535");
            }

            fits_close_file(fptr.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_float_write_long_read_overflow() {
        // Create FLOAT data that tests LONG range boundaries
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [10, 10];
            let nelements = 100;

            let mut write_data: Vec<c_float> = vec![0.0; nelements];
            write_data[0] = c_long::MIN as f32 - 1.0; // Underflow LONG
            write_data[1] = c_long::MIN as f32; // Min LONG
            write_data[2] = -100.7; // Negative with fraction
            write_data[3] = 0.0; // Zero
            write_data[4] = 100.9; // Positive with fraction
            write_data[5] = c_long::MAX as f32; // Max LONG
            write_data[6] = c_long::MAX as f32 + 1.0; // Overflow LONG
            write_data[7] = 1e20; // Large overflow
            for i in 8..nelements {
                write_data[i] = (i as c_float - 50.0) * 1000000.0;
            }

            // Write as FLOAT_IMG
            let filename_cstr = CString::new(filename).unwrap();
            fits_create_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                &mut status,
            );
            fits_create_img(
                fptr.as_mut().unwrap(),
                FLOAT_IMG,
                naxis as c_int,
                &naxes,
                &mut status,
            );
            fits_write_img(
                fptr.as_mut().unwrap(),
                TFLOAT,
                1,
                nelements as LONGLONG,
                cast_slice(&write_data),
                &mut status,
            );
            fits_close_file(fptr.take().unwrap(), &mut status);

            assert_eq!(status, 0, "Failed to write FLOAT_IMG");

            // Read as LONG - should get overflow
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

                assert_eq!(status, 412, "Expected NUM_OVERFLOW (412) for FLOAT → LONG");

                // Verify specific conversions (truncation + clamping)
                assert_eq!(
                    read_data[0],
                    c_long::MIN,
                    "Underflow should clamp to LONG min"
                );
                assert_eq!(read_data[1], c_long::MIN, "Min LONG should convert exactly");
                assert_eq!(read_data[2], -100, "Negative should truncate to -100");
                assert_eq!(read_data[3], 0, "Zero should convert exactly");
                assert_eq!(read_data[4], 100, "Positive should truncate to 100");
                assert_eq!(read_data[5], c_long::MAX, "Max LONG should convert exactly");
                assert_eq!(
                    read_data[6],
                    c_long::MAX,
                    "Overflow should clamp to LONG max"
                );
                assert_eq!(
                    read_data[7],
                    c_long::MAX,
                    "Large overflow should clamp to LONG max"
                );
            }

            fits_close_file(fptr.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_float_write_longlong_read_success() {
        // FLOAT → LONGLONG should succeed (though may lose precision)
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;
            let write_data = generate_float_test_data();

            // Write as FLOAT_IMG
            let filename_cstr = CString::new(filename).unwrap();
            fits_create_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                &mut status,
            );
            fits_create_img(
                fptr.as_mut().unwrap(),
                FLOAT_IMG,
                naxis as c_int,
                &naxes,
                &mut status,
            );
            fits_write_img(
                fptr.as_mut().unwrap(),
                TFLOAT,
                1,
                nelements as LONGLONG,
                cast_slice(&write_data),
                &mut status,
            );
            fits_close_file(fptr.take().unwrap(), &mut status);

            assert_eq!(status, 0, "Failed to write FLOAT_IMG");

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
                    "FLOAT → LONGLONG should succeed without overflow"
                );

                for i in 0..nelements {
                    let expected = write_data[i] as c_longlong;
                    assert_eq!(
                        read_data[i], expected,
                        "Data mismatch at index {}: FLOAT {} → LONGLONG {} (expected {})",
                        i, write_data[i], read_data[i], expected
                    );
                }
            }

            fits_close_file(fptr.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_float_write_sbyte_read_overflow() {
        // FLOAT data will overflow SBYTE range (-128 to 127)
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [10, 10];
            let nelements = 100;

            let mut write_data: Vec<c_float> = vec![0.0; nelements];
            write_data[0] = -129.5; // Underflow SBYTE
            write_data[1] = -128.0; // Min SBYTE
            write_data[2] = -50.7; // Negative with fraction
            write_data[3] = 0.0; // Zero
            write_data[4] = 50.9; // Positive with fraction
            write_data[5] = 127.0; // Max SBYTE
            write_data[6] = 127.9; // Just under overflow (truncates to 127)
            write_data[7] = 128.0; // Overflow SBYTE
            write_data[8] = 1000.5; // Large overflow
            for i in 9..nelements {
                write_data[i] = (i as c_float - 50.0) * 10.0;
            }

            // Write as FLOAT_IMG
            let filename_cstr = CString::new(filename).unwrap();
            fits_create_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                &mut status,
            );
            fits_create_img(
                fptr.as_mut().unwrap(),
                FLOAT_IMG,
                naxis as c_int,
                &naxes,
                &mut status,
            );
            fits_write_img(
                fptr.as_mut().unwrap(),
                TFLOAT,
                1,
                nelements as LONGLONG,
                cast_slice(&write_data),
                &mut status,
            );
            fits_close_file(fptr.take().unwrap(), &mut status);

            assert_eq!(status, 0, "Failed to write FLOAT_IMG");

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

                assert_eq!(status, 412, "Expected NUM_OVERFLOW (412) for FLOAT → SBYTE");

                // Verify specific conversions (truncation + clamping)
                assert_eq!(read_data[0], -128, "Underflow should clamp to SBYTE min");
                assert_eq!(read_data[1], -128, "Min SBYTE should convert exactly");
                assert_eq!(read_data[2], -50, "Negative should truncate to -50");
                assert_eq!(read_data[3], 0, "Zero should convert exactly");
                assert_eq!(read_data[4], 50, "Positive should truncate to 50");
                assert_eq!(read_data[5], 127, "Max SBYTE should convert exactly");
                assert_eq!(read_data[6], 127, "127.9 should truncate to 127");
                assert_eq!(read_data[7], 127, "Overflow should clamp to SBYTE max");
                assert_eq!(
                    read_data[8], 127,
                    "Large overflow should clamp to SBYTE max"
                );
            }

            fits_close_file(fptr.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_float_write_ulong_read_overflow() {
        // FLOAT with negatives will overflow ULONG (0 to ~4.3 billion)
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [10, 10];
            let nelements = 100;

            let mut write_data: Vec<c_float> = vec![0.0; nelements];
            write_data[0] = -1.5; // Negative (clamps to 0)
            write_data[1] = 0.0; // Zero
            write_data[2] = 100.7; // Positive with fraction
            write_data[3] = c_ulong::MAX as f32; // Max ULONG
            write_data[4] = (c_ulong::MAX as f32) + 1.0; // Overflow ULONG (may not be exactly representable in f32)
            write_data[5] = 1e20; // Large overflow
            for i in 6..nelements {
                write_data[i] = (i as c_float - 50.0) * 1000000.0;
            }

            // Write as FLOAT_IMG
            let filename_cstr = CString::new(filename).unwrap();
            fits_create_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                &mut status,
            );
            fits_create_img(
                fptr.as_mut().unwrap(),
                FLOAT_IMG,
                naxis as c_int,
                &naxes,
                &mut status,
            );
            fits_write_img(
                fptr.as_mut().unwrap(),
                TFLOAT,
                1,
                nelements as LONGLONG,
                cast_slice(&write_data),
                &mut status,
            );
            fits_close_file(fptr.take().unwrap(), &mut status);

            assert_eq!(status, 0, "Failed to write FLOAT_IMG");

            // Read as ULONG - should get overflow due to negatives and large values
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

                assert_eq!(status, 412, "Expected NUM_OVERFLOW (412) for FLOAT → ULONG");

                // Verify clamping: negatives → 0, large values → max
                assert_eq!(read_data[0], 0, "Negative should clamp to 0");
                assert_eq!(read_data[1], 0, "Zero should convert to 0");
                assert_eq!(read_data[2], 100, "100.7 should truncate to 100");
                assert_eq!(
                    read_data[3],
                    c_ulong::MAX,
                    "Max ULONG should convert exactly"
                );
                assert_eq!(
                    read_data[4],
                    c_ulong::MAX,
                    "Overflow should clamp to ULONG max"
                );
                assert_eq!(
                    read_data[5],
                    c_ulong::MAX,
                    "Large overflow should clamp to ULONG max"
                );
            }

            fits_close_file(fptr.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_float_write_ulonglong_read_overflow() {
        // FLOAT with negatives will overflow ULONGLONG (0 to ~18 quintillion)
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;
            let write_data = generate_float_test_data(); // Contains negatives

            // Write as FLOAT_IMG
            let filename_cstr = CString::new(filename).unwrap();
            fits_create_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                &mut status,
            );
            fits_create_img(
                fptr.as_mut().unwrap(),
                FLOAT_IMG,
                naxis as c_int,
                &naxes,
                &mut status,
            );
            fits_write_img(
                fptr.as_mut().unwrap(),
                TFLOAT,
                1,
                nelements as LONGLONG,
                cast_slice(&write_data),
                &mut status,
            );
            fits_close_file(fptr.take().unwrap(), &mut status);

            assert_eq!(status, 0, "Failed to write FLOAT_IMG");

            // Read as ULONGLONG - should get overflow due to negatives
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
                    "Expected NUM_OVERFLOW (412) for FLOAT with negatives → ULONGLONG"
                );

                // Verify clamping: negatives → 0
                for i in 0..nelements {
                    let expected = if write_data[i] < 0.0 {
                        0
                    } else {
                        write_data[i] as c_ulonglong
                    };
                    assert_eq!(
                        read_data[i], expected,
                        "Value at index {} should be clamped: FLOAT {} → ULONGLONG {} (got {})",
                        i, write_data[i], expected, read_data[i]
                    );
                }
            }

            fits_close_file(fptr.take().unwrap(), &mut status);
        });
    }
}
