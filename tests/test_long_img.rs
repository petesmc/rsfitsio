#[cfg(test)]
mod tests {
    use bytemuck::{cast_slice, cast_slice_mut};
    use libc::{
        c_double, c_float, c_int, c_long, c_longlong, c_short, c_uchar, c_ulong, c_ulonglong,
    };
    use rsfitsio::aliases::rust_api::*;
    use rsfitsio::fitsio::{
        LONG_IMG, LONGLONG, READWRITE, TBYTE, TDOUBLE, TFLOAT, TLONG, TLONGLONG, TSBYTE, TSHORT,
        TULONG, TULONGLONG, TUSHORT, fitsfile,
    };
    use rsfitsio::helpers::testhelpers::{floats_close_f32, floats_close_f64, with_temp_file};
    use std::ffi::CString;

    const IMAGE_WIDTH: c_long = 50;
    const IMAGE_HEIGHT: c_long = 50;

    /// Generate LONG_IMG test data (includes negatives)
    fn generate_long_test_data() -> Vec<c_long> {
        let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;
        let mut data = vec![0; nelements];
        for j in 0..IMAGE_HEIGHT {
            for i in 0..IMAGE_WIDTH {
                let index = (j * IMAGE_WIDTH + i) as usize;
                data[index] = ((i - 25) * 1000 + (j - 25) * 10) as c_long; // -25,250 to +24,240
            }
        }
        data
    }

    #[test]
    fn test_long_write_long_read() {
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;
            let write_data = generate_long_test_data();

            // Write as LONG_IMG
            let filename_cstr = CString::new(filename).unwrap();
            fits_create_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                &mut status,
            );
            fits_create_img(
                fptr.as_mut().unwrap(),
                LONG_IMG,
                naxis as c_int,
                &naxes,
                &mut status,
            );
            fits_write_img(
                fptr.as_mut().unwrap(),
                TLONG,
                1,
                nelements as LONGLONG,
                cast_slice(&write_data),
                &mut status,
            );
            fits_close_file(fptr.take().unwrap(), &mut status);

            assert_eq!(status, 0, "Failed to write LONG_IMG");

            // Read back as LONG
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

                assert_eq!(status, 0, "Failed to read LONG_IMG as TLONG");

                for i in 0..nelements {
                    assert_eq!(write_data[i], read_data[i], "Data mismatch at index {}", i);
                }
            }

            fits_close_file(fptr.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_long_write_byte_read_overflow() {
        // LONG data (-25250 to +24240) will overflow BYTE range (0-255)
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;
            let write_data = generate_long_test_data();

            // Write as LONG_IMG
            let filename_cstr = CString::new(filename).unwrap();
            fits_create_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                &mut status,
            );
            fits_create_img(
                fptr.as_mut().unwrap(),
                LONG_IMG,
                naxis as c_int,
                &naxes,
                &mut status,
            );
            fits_write_img(
                fptr.as_mut().unwrap(),
                TLONG,
                1,
                nelements as LONGLONG,
                cast_slice(&write_data),
                &mut status,
            );
            fits_close_file(fptr.take().unwrap(), &mut status);

            assert_eq!(status, 0, "Failed to write LONG_IMG");

            // Read as BYTE - should get overflow
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

                assert_eq!(status, 412, "Expected NUM_OVERFLOW (412) for LONG → BYTE");

                // Verify clamping behavior
                for i in 0..nelements {
                    let expected = if write_data[i] < 0 {
                        0
                    } else if write_data[i] > 255 {
                        255
                    } else {
                        write_data[i] as c_uchar
                    };
                    assert_eq!(
                        read_data[i], expected,
                        "Value at index {} should be clamped: LONG {} → BYTE {} (got {})",
                        i, write_data[i], expected, read_data[i]
                    );
                }
            }

            fits_close_file(fptr.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_long_write_short_read_overflow() {
        // LONG data (-25250 to +24240) will mostly fit in SHORT (-32768 to +32767)
        // but let's create data that specifically tests SHORT overflow
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [10, 10];
            let nelements = 100;

            // Create LONG data that tests SHORT boundaries
            let mut write_data: Vec<c_long> = vec![0; nelements];
            write_data[0] = -32769; // Underflow SHORT
            write_data[1] = -32768; // Min SHORT
            write_data[2] = 0; // Valid
            write_data[3] = 32767; // Max SHORT
            write_data[4] = 32768; // Overflow SHORT
            write_data[5] = 100000; // Large overflow
            for i in 6..nelements {
                write_data[i] = (i as c_long - 50) * 100; // Mix of values
            }

            // Write as LONG_IMG
            let filename_cstr = CString::new(filename).unwrap();
            fits_create_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                &mut status,
            );
            fits_create_img(
                fptr.as_mut().unwrap(),
                LONG_IMG,
                naxis as c_int,
                &naxes,
                &mut status,
            );
            fits_write_img(
                fptr.as_mut().unwrap(),
                TLONG,
                1,
                nelements as LONGLONG,
                cast_slice(&write_data),
                &mut status,
            );
            fits_close_file(fptr.take().unwrap(), &mut status);

            assert_eq!(status, 0, "Failed to write LONG_IMG");

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

                assert_eq!(status, 412, "Expected NUM_OVERFLOW (412) for LONG → SHORT");

                // Verify specific clamping
                assert_eq!(read_data[0], -32768, "Underflow should clamp to SHORT min");
                assert_eq!(read_data[1], -32768, "Min SHORT should convert exactly");
                assert_eq!(read_data[2], 0, "Zero should convert exactly");
                assert_eq!(read_data[3], 32767, "Max SHORT should convert exactly");
                assert_eq!(read_data[4], 32767, "Overflow should clamp to SHORT max");
                assert_eq!(
                    read_data[5], 32767,
                    "Large overflow should clamp to SHORT max"
                );
            }

            fits_close_file(fptr.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_long_write_ushort_read_overflow() {
        // LONG with negatives will overflow USHORT (0-65535)
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;
            let write_data = generate_long_test_data(); // Contains negatives

            // Write as LONG_IMG
            let filename_cstr = CString::new(filename).unwrap();
            fits_create_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                &mut status,
            );
            fits_create_img(
                fptr.as_mut().unwrap(),
                LONG_IMG,
                naxis as c_int,
                &naxes,
                &mut status,
            );
            fits_write_img(
                fptr.as_mut().unwrap(),
                TLONG,
                1,
                nelements as LONGLONG,
                cast_slice(&write_data),
                &mut status,
            );
            fits_close_file(fptr.take().unwrap(), &mut status);

            assert_eq!(status, 0, "Failed to write LONG_IMG");

            // Read as USHORT - should get overflow due to negatives
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
                    "Expected NUM_OVERFLOW (412) for LONG with negatives → USHORT"
                );

                // Verify clamping: negatives → 0, large values → 65535
                for i in 0..nelements {
                    let expected = if write_data[i] < 0 {
                        0
                    } else if write_data[i] > 65535 {
                        65535
                    } else {
                        write_data[i] as u16
                    };
                    assert_eq!(
                        read_data[i], expected,
                        "Value at index {} should be clamped: LONG {} → USHORT {} (got {})",
                        i, write_data[i], expected, read_data[i]
                    );
                }
            }

            fits_close_file(fptr.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_long_write_longlong_read_success() {
        // LONG → LONGLONG should always succeed
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;
            let write_data = generate_long_test_data();

            // Write as LONG_IMG
            let filename_cstr = CString::new(filename).unwrap();
            fits_create_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                &mut status,
            );
            fits_create_img(
                fptr.as_mut().unwrap(),
                LONG_IMG,
                naxis as c_int,
                &naxes,
                &mut status,
            );
            fits_write_img(
                fptr.as_mut().unwrap(),
                TLONG,
                1,
                nelements as LONGLONG,
                cast_slice(&write_data),
                &mut status,
            );
            fits_close_file(fptr.take().unwrap(), &mut status);

            assert_eq!(status, 0, "Failed to write LONG_IMG");

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

                assert_eq!(status, 0, "LONG → LONGLONG should succeed without overflow");

                for i in 0..nelements {
                    assert_eq!(
                        write_data[i] as c_longlong, read_data[i],
                        "Data mismatch at index {}: LONG {} → LONGLONG {}",
                        i, write_data[i], read_data[i]
                    );
                }
            }

            fits_close_file(fptr.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_long_write_float_read_success() {
        // LONG → FLOAT should succeed (though may lose precision for very large values)
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;
            let write_data = generate_long_test_data();

            // Write as LONG_IMG
            let filename_cstr = CString::new(filename).unwrap();
            fits_create_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                &mut status,
            );
            fits_create_img(
                fptr.as_mut().unwrap(),
                LONG_IMG,
                naxis as c_int,
                &naxes,
                &mut status,
            );
            fits_write_img(
                fptr.as_mut().unwrap(),
                TLONG,
                1,
                nelements as LONGLONG,
                cast_slice(&write_data),
                &mut status,
            );
            fits_close_file(fptr.take().unwrap(), &mut status);

            assert_eq!(status, 0, "Failed to write LONG_IMG");

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

                assert_eq!(status, 0, "LONG → FLOAT should succeed without overflow");

                for i in 0..nelements {
                    let expected = write_data[i] as c_float;
                    assert!(
                        floats_close_f32(expected, read_data[i]),
                        "Data mismatch at index {}: LONG {} → FLOAT {} (expected {})",
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
    fn test_long_write_double_read_success() {
        // LONG → DOUBLE should always succeed
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;
            let write_data = generate_long_test_data();

            // Write as LONG_IMG
            let filename_cstr = CString::new(filename).unwrap();
            fits_create_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                &mut status,
            );
            fits_create_img(
                fptr.as_mut().unwrap(),
                LONG_IMG,
                naxis as c_int,
                &naxes,
                &mut status,
            );
            fits_write_img(
                fptr.as_mut().unwrap(),
                TLONG,
                1,
                nelements as LONGLONG,
                cast_slice(&write_data),
                &mut status,
            );
            fits_close_file(fptr.take().unwrap(), &mut status);

            assert_eq!(status, 0, "Failed to write LONG_IMG");

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

                assert_eq!(status, 0, "LONG → DOUBLE should succeed without overflow");

                for i in 0..nelements {
                    let expected = write_data[i] as c_double;
                    assert!(
                        floats_close_f64(expected, read_data[i]),
                        "Data mismatch at index {}: LONG {} → DOUBLE {} (expected {})",
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
    fn test_long_write_sbyte_read_overflow() {
        // LONG data (-25250 to +24240) will overflow SBYTE range (-128 to 127)
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;
            let write_data = generate_long_test_data();

            // Write as LONG_IMG
            let filename_cstr = CString::new(filename).unwrap();
            fits_create_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                &mut status,
            );
            fits_create_img(
                fptr.as_mut().unwrap(),
                LONG_IMG,
                naxis as c_int,
                &naxes,
                &mut status,
            );
            fits_write_img(
                fptr.as_mut().unwrap(),
                TLONG,
                1,
                nelements as LONGLONG,
                cast_slice(&write_data),
                &mut status,
            );
            fits_close_file(fptr.take().unwrap(), &mut status);

            assert_eq!(status, 0, "Failed to write LONG_IMG");

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

                assert_eq!(status, 412, "Expected NUM_OVERFLOW (412) for LONG → SBYTE");

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
                        "Value at index {} should be clamped: LONG {} → SBYTE {} (got {})",
                        i, write_data[i], expected, read_data[i]
                    );
                }
            }

            fits_close_file(fptr.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_long_write_ulong_read_overflow() {
        // LONG with negatives will overflow ULONG (0 to ~4.3 billion)
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;
            let write_data = generate_long_test_data(); // Contains negatives

            // Write as LONG_IMG
            let filename_cstr = CString::new(filename).unwrap();
            fits_create_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                &mut status,
            );
            fits_create_img(
                fptr.as_mut().unwrap(),
                LONG_IMG,
                naxis as c_int,
                &naxes,
                &mut status,
            );
            fits_write_img(
                fptr.as_mut().unwrap(),
                TLONG,
                1,
                nelements as LONGLONG,
                cast_slice(&write_data),
                &mut status,
            );
            fits_close_file(fptr.take().unwrap(), &mut status);

            assert_eq!(status, 0, "Failed to write LONG_IMG");

            // Read as ULONG - should get overflow due to negatives
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
                    "Expected NUM_OVERFLOW (412) for LONG with negatives → ULONG"
                );

                // Verify clamping: negatives → 0
                for i in 0..nelements {
                    let expected = if write_data[i] < 0 {
                        0
                    } else {
                        write_data[i] as c_ulong
                    };
                    assert_eq!(
                        read_data[i], expected,
                        "Value at index {} should be clamped: LONG {} → ULONG {} (got {})",
                        i, write_data[i], expected, read_data[i]
                    );
                }
            }

            fits_close_file(fptr.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_long_write_ulonglong_read_overflow() {
        // LONG with negatives will overflow ULONGLONG (0 to ~18 quintillion)
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;
            let write_data = generate_long_test_data(); // Contains negatives

            // Write as LONG_IMG
            let filename_cstr = CString::new(filename).unwrap();
            fits_create_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                &mut status,
            );
            fits_create_img(
                fptr.as_mut().unwrap(),
                LONG_IMG,
                naxis as c_int,
                &naxes,
                &mut status,
            );
            fits_write_img(
                fptr.as_mut().unwrap(),
                TLONG,
                1,
                nelements as LONGLONG,
                cast_slice(&write_data),
                &mut status,
            );
            fits_close_file(fptr.take().unwrap(), &mut status);

            assert_eq!(status, 0, "Failed to write LONG_IMG");

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
                    "Expected NUM_OVERFLOW (412) for LONG with negatives → ULONGLONG"
                );

                // Verify clamping: negatives → 0
                for i in 0..nelements {
                    let expected = if write_data[i] < 0 {
                        0
                    } else {
                        write_data[i] as c_ulonglong
                    };
                    assert_eq!(
                        read_data[i], expected,
                        "Value at index {} should be clamped: LONG {} → ULONGLONG {} (got {})",
                        i, write_data[i], expected, read_data[i]
                    );
                }
            }

            fits_close_file(fptr.take().unwrap(), &mut status);
        });
    }
}
