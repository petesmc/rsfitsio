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

    /// Generate ULONGLONG_IMG test data (0 to values that exceed LONG_MAX)
    fn generate_ulonglong_test_data() -> Vec<c_ulonglong> {
        let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;
        let mut data = vec![0; nelements];
        for j in 0..IMAGE_HEIGHT {
            for i in 0..IMAGE_WIDTH {
                let index = (j * IMAGE_WIDTH + i) as usize;
                // Create values that will exceed LONG_MAX for testing overflow
                // LONG_MAX = 9,223,372,036,854,775,807
                let base = (i as c_ulonglong + j as c_ulonglong);
                if base == 0 {
                    data[index] = 0;
                } else if base < 50 {
                    // For small bases, use safe multiplication that creates values > LONG_MAX
                    data[index] = (c_long::MAX as c_ulonglong) + base * 1_000_000_000u64;
                } else {
                    // For larger bases, use direct large values
                    data[index] = c_ulonglong::MAX - base * 1000u64;
                }
            }
        }
        data
    }

    #[test]
    fn test_ulonglong_write_ulonglong_read() {
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;
            let write_data = generate_ulonglong_test_data();

            // Write as ULONGLONG_IMG
            let filename_cstr = CString::new(filename).unwrap();
            unsafe {
                fits_create_file(
                    &mut fptr,
                    cast_slice(filename_cstr.to_bytes_with_nul()),
                    &mut status,
                );
                fits_create_img(
                    fptr.as_mut().unwrap(),
                    ULONGLONG_IMG,
                    naxis as c_int,
                    &naxes,
                    &mut status,
                );
                fits_write_img(
                    fptr.as_mut().unwrap(),
                    TULONGLONG,
                    1,
                    nelements as LONGLONG,
                    cast_slice(&write_data),
                    &mut status,
                );
                fits_close_file(fptr.take().unwrap(), &mut status);
            }
            assert_eq!(status, 0, "Failed to write ULONGLONG_IMG");

            // Read back as ULONGLONG
            fptr = None;
            unsafe {
                fits_open_file(
                    &mut fptr,
                    cast_slice(filename_cstr.to_bytes_with_nul()),
                    READWRITE,
                    &mut status,
                );
            }

            if let Some(ref mut fptr_box) = fptr {
                let mut read_data: Vec<c_ulonglong> = vec![0; nelements];
                let mut anynull: c_int = 0;

                unsafe {
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
                }
                assert_eq!(status, 0, "Failed to read ULONGLONG_IMG as TULONGLONG");

                for i in 0..nelements {
                    assert_eq!(write_data[i], read_data[i], "Data mismatch at index {}", i);
                }
            }

            unsafe {
                fits_close_file(fptr.take().unwrap(), &mut status);
            }
        });
    }

    #[test]
    fn test_ulonglong_write_byte_read_overflow() {
        // ULONGLONG data (0-9800000) will overflow BYTE range (0-255)
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;
            let write_data = generate_ulonglong_test_data();

            // Write as ULONGLONG_IMG
            let filename_cstr = CString::new(filename).unwrap();
            unsafe {
                fits_create_file(
                    &mut fptr,
                    cast_slice(filename_cstr.to_bytes_with_nul()),
                    &mut status,
                );
                fits_create_img(
                    fptr.as_mut().unwrap(),
                    ULONGLONG_IMG,
                    naxis as c_int,
                    &naxes,
                    &mut status,
                );
                fits_write_img(
                    fptr.as_mut().unwrap(),
                    TULONGLONG,
                    1,
                    nelements as LONGLONG,
                    cast_slice(&write_data),
                    &mut status,
                );
                fits_close_file(fptr.take().unwrap(), &mut status);
            }
            assert_eq!(status, 0, "Failed to write ULONGLONG_IMG");

            // Read as BYTE - should get overflow for values > 255
            fptr = None;
            unsafe {
                fits_open_file(
                    &mut fptr,
                    cast_slice(filename_cstr.to_bytes_with_nul()),
                    READWRITE,
                    &mut status,
                );
            }

            if let Some(ref mut fptr_box) = fptr {
                let mut read_data: Vec<c_uchar> = vec![0; nelements];
                let mut anynull: c_int = 0;

                unsafe {
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
                }

                assert_eq!(
                    status, 412,
                    "Expected NUM_OVERFLOW (412) for ULONGLONG → BYTE"
                );

                // Verify clamping behavior - almost all values should clamp to 255
                for i in 0..nelements {
                    let expected = if write_data[i] > 255 {
                        255
                    } else {
                        write_data[i] as c_uchar
                    };
                    assert_eq!(
                        read_data[i], expected,
                        "Value at index {} should be clamped: ULONGLONG {} → BYTE {}",
                        i, write_data[i], expected
                    );
                }
            }

            unsafe {
                fits_close_file(fptr.take().unwrap(), &mut status);
            }
        });
    }

    #[test]
    fn test_ulonglong_write_longlong_read_overflow() {
        // Create ULONGLONG data that exceeds LONGLONG max
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [10, 10];
            let nelements = 100;

            // Create ULONGLONG data that tests LONGLONG boundaries
            let mut write_data: Vec<c_ulonglong> = vec![0; nelements];
            write_data[0] = 0; // Valid
            write_data[1] = 9223372036854775807; // Max LONGLONG
            write_data[2] = 9223372036854775808; // Overflow LONGLONG (LONGLONG max + 1)
            write_data[3] = 18446744073709551615; // Max ULONGLONG
            for i in 4..nelements {
                write_data[i] = (i as c_ulonglong) * 1000000000000; // Mix of large values
            }

            // Write as ULONGLONG_IMG
            let filename_cstr = CString::new(filename).unwrap();
            unsafe {
                fits_create_file(
                    &mut fptr,
                    cast_slice(filename_cstr.to_bytes_with_nul()),
                    &mut status,
                );
                fits_create_img(
                    fptr.as_mut().unwrap(),
                    ULONGLONG_IMG,
                    naxis as c_int,
                    &naxes,
                    &mut status,
                );
                fits_write_img(
                    fptr.as_mut().unwrap(),
                    TULONGLONG,
                    1,
                    nelements as LONGLONG,
                    cast_slice(&write_data),
                    &mut status,
                );
                fits_close_file(fptr.take().unwrap(), &mut status);
            }
            assert_eq!(status, 0, "Failed to write ULONGLONG_IMG");

            // Read as LONGLONG - should get overflow
            fptr = None;
            unsafe {
                fits_open_file(
                    &mut fptr,
                    cast_slice(filename_cstr.to_bytes_with_nul()),
                    READWRITE,
                    &mut status,
                );
            }

            if let Some(ref mut fptr_box) = fptr {
                let mut read_data: Vec<c_longlong> = vec![0; nelements];
                let mut anynull: c_int = 0;

                unsafe {
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
                }

                assert_eq!(
                    status, 412,
                    "Expected NUM_OVERFLOW (412) for ULONGLONG → LONGLONG"
                );

                // Verify clamping to LONGLONG range
                assert_eq!(read_data[0], 0, "Zero should convert exactly");
                assert_eq!(
                    read_data[1], 9223372036854775807,
                    "Max LONGLONG should convert exactly"
                );
                assert_eq!(
                    read_data[2], 9223372036854775807,
                    "Overflow should clamp to LONGLONG max"
                );
                assert_eq!(
                    read_data[3], 9223372036854775807,
                    "Max ULONGLONG should clamp to LONGLONG max"
                );
            }

            unsafe {
                fits_close_file(fptr.take().unwrap(), &mut status);
            }
        });
    }

    #[test]
    fn test_ulonglong_write_double_read_success() {
        // ULONGLONG → DOUBLE should succeed (though may lose precision for very large values)
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;
            let write_data = generate_ulonglong_test_data();

            // Write as ULONGLONG_IMG
            let filename_cstr = CString::new(filename).unwrap();
            unsafe {
                fits_create_file(
                    &mut fptr,
                    cast_slice(filename_cstr.to_bytes_with_nul()),
                    &mut status,
                );
                fits_create_img(
                    fptr.as_mut().unwrap(),
                    ULONGLONG_IMG,
                    naxis as c_int,
                    &naxes,
                    &mut status,
                );
                fits_write_img(
                    fptr.as_mut().unwrap(),
                    TULONGLONG,
                    1,
                    nelements as LONGLONG,
                    cast_slice(&write_data),
                    &mut status,
                );
                fits_close_file(fptr.take().unwrap(), &mut status);
            }
            assert_eq!(status, 0, "Failed to write ULONGLONG_IMG");

            // Read as DOUBLE - should succeed
            fptr = None;
            unsafe {
                fits_open_file(
                    &mut fptr,
                    cast_slice(filename_cstr.to_bytes_with_nul()),
                    READWRITE,
                    &mut status,
                );
            }

            if let Some(ref mut fptr_box) = fptr {
                let mut read_data: Vec<c_double> = vec![0.0; nelements];
                let mut anynull: c_int = 0;

                unsafe {
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
                }
                assert_eq!(
                    status, 0,
                    "ULONGLONG → DOUBLE should succeed without overflow"
                );

                for i in 0..nelements {
                    let expected = write_data[i] as c_double;
                    assert!(
                        floats_close_f64(expected, read_data[i]),
                        "Data mismatch at index {}: ULONGLONG {} → DOUBLE {} (expected {})",
                        i,
                        write_data[i],
                        read_data[i],
                        expected
                    );
                }
            }

            unsafe {
                fits_close_file(fptr.take().unwrap(), &mut status);
            }
        });
    }

    #[test]
    fn test_ulonglong_write_short_read_overflow() {
        // ULONGLONG data (0 to ~9.8 million) will overflow SHORT range (-32768 to +32767)
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;
            let write_data = generate_ulonglong_test_data();

            // Write as ULONGLONG_IMG
            let filename_cstr = CString::new(filename).unwrap();
            unsafe {
                fits_create_file(
                    &mut fptr,
                    cast_slice(filename_cstr.to_bytes_with_nul()),
                    &mut status,
                );
                fits_create_img(
                    fptr.as_mut().unwrap(),
                    ULONGLONG_IMG,
                    naxis as c_int,
                    &naxes,
                    &mut status,
                );
                fits_write_img(
                    fptr.as_mut().unwrap(),
                    TULONGLONG,
                    1,
                    nelements as LONGLONG,
                    cast_slice(&write_data),
                    &mut status,
                );
                fits_close_file(fptr.take().unwrap(), &mut status);
            }
            assert_eq!(status, 0, "Failed to write ULONGLONG_IMG");

            // Read as SHORT - should get overflow
            fptr = None;
            unsafe {
                fits_open_file(
                    &mut fptr,
                    cast_slice(filename_cstr.to_bytes_with_nul()),
                    READWRITE,
                    &mut status,
                );
            }

            if let Some(ref mut fptr_box) = fptr {
                let mut read_data: Vec<c_short> = vec![0; nelements];
                let mut anynull: c_int = 0;

                unsafe {
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
                }

                assert_eq!(
                    status, 412,
                    "Expected NUM_OVERFLOW (412) for ULONGLONG → SHORT"
                );

                // Verify clamping behavior
                for i in 0..nelements {
                    let expected = if write_data[i] > 32767 {
                        32767
                    } else {
                        write_data[i] as c_short
                    };
                    assert_eq!(
                        read_data[i], expected,
                        "Value at index {} should be clamped: ULONGLONG {} → SHORT {} (got {})",
                        i, write_data[i], expected, read_data[i]
                    );
                }
            }

            unsafe {
                fits_close_file(fptr.take().unwrap(), &mut status);
            }
        });
    }

    #[test]
    fn test_ulonglong_write_ushort_read_overflow() {
        // ULONGLONG data (0 to ~9.8 million) will overflow USHORT range (0 to 65535)
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;
            let write_data = generate_ulonglong_test_data();

            // Write as ULONGLONG_IMG
            let filename_cstr = CString::new(filename).unwrap();
            unsafe {
                fits_create_file(
                    &mut fptr,
                    cast_slice(filename_cstr.to_bytes_with_nul()),
                    &mut status,
                );
                fits_create_img(
                    fptr.as_mut().unwrap(),
                    ULONGLONG_IMG,
                    naxis as c_int,
                    &naxes,
                    &mut status,
                );
                fits_write_img(
                    fptr.as_mut().unwrap(),
                    TULONGLONG,
                    1,
                    nelements as LONGLONG,
                    cast_slice(&write_data),
                    &mut status,
                );
                fits_close_file(fptr.take().unwrap(), &mut status);
            }
            assert_eq!(status, 0, "Failed to write ULONGLONG_IMG");

            // Read as USHORT - should get overflow
            fptr = None;
            unsafe {
                fits_open_file(
                    &mut fptr,
                    cast_slice(filename_cstr.to_bytes_with_nul()),
                    READWRITE,
                    &mut status,
                );
            }

            if let Some(ref mut fptr_box) = fptr {
                let mut read_data: Vec<u16> = vec![0; nelements];
                let mut anynull: c_int = 0;

                unsafe {
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
                }

                assert_eq!(
                    status, 412,
                    "Expected NUM_OVERFLOW (412) for ULONGLONG → USHORT"
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
                        "Value at index {} should be clamped: ULONGLONG {} → USHORT {} (got {})",
                        i, write_data[i], expected, read_data[i]
                    );
                }
            }

            unsafe {
                fits_close_file(fptr.take().unwrap(), &mut status);
            }
        });
    }

    #[test]
    fn test_ulonglong_write_long_read_overflow() {
        // ULONGLONG data (0 to ~9.8 million) will overflow LONG range (-2.1B to +2.1B)
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;
            let write_data = generate_ulonglong_test_data();

            // Write as ULONGLONG_IMG
            let filename_cstr = CString::new(filename).unwrap();
            unsafe {
                fits_create_file(
                    &mut fptr,
                    cast_slice(filename_cstr.to_bytes_with_nul()),
                    &mut status,
                );
                fits_create_img(
                    fptr.as_mut().unwrap(),
                    ULONGLONG_IMG,
                    naxis as c_int,
                    &naxes,
                    &mut status,
                );
                fits_write_img(
                    fptr.as_mut().unwrap(),
                    TULONGLONG,
                    1,
                    nelements as LONGLONG,
                    cast_slice(&write_data),
                    &mut status,
                );
                fits_close_file(fptr.take().unwrap(), &mut status);
            }
            assert_eq!(status, 0, "Failed to write ULONGLONG_IMG");

            // Read as LONG - should get overflow for large values
            fptr = None;
            unsafe {
                fits_open_file(
                    &mut fptr,
                    cast_slice(filename_cstr.to_bytes_with_nul()),
                    READWRITE,
                    &mut status,
                );
            }

            if let Some(ref mut fptr_box) = fptr {
                let mut read_data: Vec<c_long> = vec![0; nelements];
                let mut anynull: c_int = 0;

                unsafe {
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
                }

                assert_eq!(
                    status, 412,
                    "Expected NUM_OVERFLOW (412) for ULONGLONG → LONG"
                );

                // Verify clamping behavior
                for i in 0..nelements {
                    let expected = if write_data[i] > c_long::MAX as c_ulonglong {
                        c_long::MAX
                    } else {
                        write_data[i] as c_long
                    };
                    assert_eq!(
                        read_data[i], expected,
                        "Value at index {} should be clamped: ULONGLONG {} → LONG {} (got {})",
                        i, write_data[i], expected, read_data[i]
                    );
                }
            }

            unsafe {
                fits_close_file(fptr.take().unwrap(), &mut status);
            }
        });
    }

    #[test]
    fn test_ulonglong_write_float_read_success() {
        // ULONGLONG → FLOAT should succeed (though may lose precision for very large values)
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;
            let write_data = generate_ulonglong_test_data();

            // Write as ULONGLONG_IMG
            let filename_cstr = CString::new(filename).unwrap();
            unsafe {
                fits_create_file(
                    &mut fptr,
                    cast_slice(filename_cstr.to_bytes_with_nul()),
                    &mut status,
                );
                fits_create_img(
                    fptr.as_mut().unwrap(),
                    ULONGLONG_IMG,
                    naxis as c_int,
                    &naxes,
                    &mut status,
                );
                fits_write_img(
                    fptr.as_mut().unwrap(),
                    TULONGLONG,
                    1,
                    nelements as LONGLONG,
                    cast_slice(&write_data),
                    &mut status,
                );
                fits_close_file(fptr.take().unwrap(), &mut status);
            }
            assert_eq!(status, 0, "Failed to write ULONGLONG_IMG");

            // Read as FLOAT - should succeed
            fptr = None;
            unsafe {
                fits_open_file(
                    &mut fptr,
                    cast_slice(filename_cstr.to_bytes_with_nul()),
                    READWRITE,
                    &mut status,
                );
            }

            if let Some(ref mut fptr_box) = fptr {
                let mut read_data: Vec<c_float> = vec![0.0; nelements];
                let mut anynull: c_int = 0;

                unsafe {
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
                }
                assert_eq!(
                    status, 0,
                    "ULONGLONG → FLOAT should succeed without overflow"
                );

                for i in 0..nelements {
                    let expected = write_data[i] as c_float;
                    assert!(
                        floats_close_f32(expected, read_data[i]),
                        "Data mismatch at index {}: ULONGLONG {} → FLOAT {} (expected {})",
                        i,
                        write_data[i],
                        read_data[i],
                        expected
                    );
                }
            }

            unsafe {
                fits_close_file(fptr.take().unwrap(), &mut status);
            }
        });
    }

    #[test]
    fn test_ulonglong_write_sbyte_read_overflow() {
        // ULONGLONG data (0 to ~9.8 million) will overflow SBYTE range (-128 to 127)
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;
            let write_data = generate_ulonglong_test_data();

            // Write as ULONGLONG_IMG
            let filename_cstr = CString::new(filename).unwrap();
            unsafe {
                fits_create_file(
                    &mut fptr,
                    cast_slice(filename_cstr.to_bytes_with_nul()),
                    &mut status,
                );
                fits_create_img(
                    fptr.as_mut().unwrap(),
                    ULONGLONG_IMG,
                    naxis as c_int,
                    &naxes,
                    &mut status,
                );
                fits_write_img(
                    fptr.as_mut().unwrap(),
                    TULONGLONG,
                    1,
                    nelements as LONGLONG,
                    cast_slice(&write_data),
                    &mut status,
                );
                fits_close_file(fptr.take().unwrap(), &mut status);
            }
            assert_eq!(status, 0, "Failed to write ULONGLONG_IMG");

            // Read as SBYTE - should get overflow
            fptr = None;
            unsafe {
                fits_open_file(
                    &mut fptr,
                    cast_slice(filename_cstr.to_bytes_with_nul()),
                    READWRITE,
                    &mut status,
                );
            }

            if let Some(ref mut fptr_box) = fptr {
                let mut read_data: Vec<i8> = vec![0; nelements];
                let mut anynull: c_int = 0;

                unsafe {
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
                }

                assert_eq!(
                    status, 412,
                    "Expected NUM_OVERFLOW (412) for ULONGLONG → SBYTE"
                );

                // Verify clamping behavior
                for i in 0..nelements {
                    let expected = if write_data[i] > 127 {
                        127
                    } else {
                        write_data[i] as i8
                    };
                    assert_eq!(
                        read_data[i], expected,
                        "Value at index {} should be clamped: ULONGLONG {} → SBYTE {} (got {})",
                        i, write_data[i], expected, read_data[i]
                    );
                }
            }

            unsafe {
                fits_close_file(fptr.take().unwrap(), &mut status);
            }
        });
    }

    #[test]
    fn test_ulonglong_write_ulong_read_overflow() {
        // ULONGLONG data (0 to ~9.8 million) may overflow ULONG (depends on system, usually 32-bit ~4.3B)
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;
            let write_data = generate_ulonglong_test_data();

            // Write as ULONGLONG_IMG
            let filename_cstr = CString::new(filename).unwrap();
            unsafe {
                fits_create_file(
                    &mut fptr,
                    cast_slice(filename_cstr.to_bytes_with_nul()),
                    &mut status,
                );
                fits_create_img(
                    fptr.as_mut().unwrap(),
                    ULONGLONG_IMG,
                    naxis as c_int,
                    &naxes,
                    &mut status,
                );
                fits_write_img(
                    fptr.as_mut().unwrap(),
                    TULONGLONG,
                    1,
                    nelements as LONGLONG,
                    cast_slice(&write_data),
                    &mut status,
                );
                fits_close_file(fptr.take().unwrap(), &mut status);
            }
            assert_eq!(status, 0, "Failed to write ULONGLONG_IMG");

            // Read as ULONG - should succeed for our data range (max ~9.8M)
            fptr = None;
            unsafe {
                fits_open_file(
                    &mut fptr,
                    cast_slice(filename_cstr.to_bytes_with_nul()),
                    READWRITE,
                    &mut status,
                );
            }

            if let Some(ref mut fptr_box) = fptr {
                let mut read_data: Vec<c_ulong> = vec![0; nelements];
                let mut anynull: c_int = 0;

                unsafe {
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
                }
                assert_eq!(
                    status, 0,
                    "ULONGLONG → ULONG should succeed for our data range"
                );

                for i in 0..nelements {
                    let expected = write_data[i] as c_ulong;
                    assert_eq!(
                        read_data[i], expected,
                        "Data mismatch at index {}: ULONGLONG {} → ULONG {} (expected {})",
                        i, write_data[i], read_data[i], expected
                    );
                }
            }

            unsafe {
                fits_close_file(fptr.take().unwrap(), &mut status);
            }
        });
    }
}
