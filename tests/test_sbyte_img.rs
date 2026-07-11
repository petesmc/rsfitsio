mod common;

#[cfg(test)]
mod tests {
    use crate::common::{floats_close_f32, floats_close_f64, with_temp_file};
    use bytemuck::{cast_slice, cast_slice_mut};
    use libc::{
        c_double, c_float, c_int, c_long, c_longlong, c_short, c_uchar, c_ulong, c_ulonglong,
    };
    use rsfitsio::aliases::rust_api::*;
    use rsfitsio::fitsio::{
        LONGLONG, READWRITE, SBYTE_IMG, TBYTE, TDOUBLE, TFLOAT, TLONG, TLONGLONG, TSBYTE, TSHORT,
        TULONG, TULONGLONG, TUSHORT, fitsfile,
    };
    use std::ffi::CString;

    const IMAGE_WIDTH: c_long = 50;
    const IMAGE_HEIGHT: c_long = 50;

    /// Generate SBYTE_IMG test data (-128 to 127 range)
    fn generate_sbyte_test_data() -> Vec<i8> {
        let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;
        let mut data = vec![0; nelements];
        for j in 0..IMAGE_HEIGHT {
            for i in 0..IMAGE_WIDTH {
                let index = (j * IMAGE_WIDTH + i) as usize;
                data[index] = ((i + j) % 255 - 127) as i8; // -127 to 127 range
            }
        }
        data
    }

    #[test]
    fn test_sbyte_write_sbyte_read() {
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;
            let write_data = generate_sbyte_test_data();

            // Write as SBYTE_IMG
            let filename_cstr = CString::new(filename).unwrap();
            fits_create_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                &mut status,
            );
            fits_create_img(
                fptr.as_mut().unwrap(),
                SBYTE_IMG,
                naxis as c_int,
                &naxes,
                &mut status,
            );
            fits_write_img(
                fptr.as_mut().unwrap(),
                TSBYTE,
                1,
                nelements as LONGLONG,
                cast_slice(&write_data),
                &mut status,
            );
            fits_close_file(fptr.take().unwrap(), &mut status);

            assert_eq!(status, 0, "Failed to write SBYTE_IMG");

            // Read back as SBYTE
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

                assert_eq!(status, 0, "Failed to read SBYTE_IMG as TSBYTE");

                for i in 0..nelements {
                    assert_eq!(write_data[i], read_data[i], "Data mismatch at index {}", i);
                }
            }

            fits_close_file(fptr.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_sbyte_write_byte_read_overflow() {
        // SBYTE with negative values will overflow BYTE (0-255)
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;
            let write_data = generate_sbyte_test_data(); // Contains negative values

            // Write as SBYTE_IMG
            let filename_cstr = CString::new(filename).unwrap();
            fits_create_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                &mut status,
            );
            fits_create_img(
                fptr.as_mut().unwrap(),
                SBYTE_IMG,
                naxis as c_int,
                &naxes,
                &mut status,
            );
            fits_write_img(
                fptr.as_mut().unwrap(),
                TSBYTE,
                1,
                nelements as LONGLONG,
                cast_slice(&write_data),
                &mut status,
            );
            fits_close_file(fptr.take().unwrap(), &mut status);

            assert_eq!(status, 0, "Failed to write SBYTE_IMG");

            // Read as BYTE - should get overflow due to negatives
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
                    status, 412,
                    "Expected NUM_OVERFLOW (412) for SBYTE with negatives → BYTE"
                );

                // Verify clamping: negatives → 0, positives unchanged
                for i in 0..nelements {
                    let expected = if write_data[i] < 0 {
                        0
                    } else {
                        write_data[i] as c_uchar
                    };
                    assert_eq!(
                        read_data[i], expected,
                        "Value at index {} should be clamped: SBYTE {} → BYTE {}",
                        i, write_data[i], expected
                    );
                }
            }

            fits_close_file(fptr.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_sbyte_write_ushort_read_overflow() {
        // SBYTE with negative values will overflow USHORT (0-65535)
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;
            let write_data = generate_sbyte_test_data(); // Contains negative values

            // Write as SBYTE_IMG
            let filename_cstr = CString::new(filename).unwrap();
            fits_create_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                &mut status,
            );
            fits_create_img(
                fptr.as_mut().unwrap(),
                SBYTE_IMG,
                naxis as c_int,
                &naxes,
                &mut status,
            );
            fits_write_img(
                fptr.as_mut().unwrap(),
                TSBYTE,
                1,
                nelements as LONGLONG,
                cast_slice(&write_data),
                &mut status,
            );
            fits_close_file(fptr.take().unwrap(), &mut status);

            assert_eq!(status, 0, "Failed to write SBYTE_IMG");

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
                    "Expected NUM_OVERFLOW (412) for SBYTE with negatives → USHORT"
                );

                // Verify clamping: negatives → 0, positives unchanged
                for i in 0..nelements {
                    let expected = if write_data[i] < 0 {
                        0
                    } else {
                        write_data[i] as u16
                    };
                    assert_eq!(
                        read_data[i], expected,
                        "Value at index {} should be clamped: SBYTE {} → USHORT {}",
                        i, write_data[i], expected
                    );
                }
            }

            fits_close_file(fptr.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_sbyte_write_short_read_success() {
        // SBYTE → SHORT should always succeed (SBYTE range -128 to 127 fits in SHORT)
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;
            let write_data = generate_sbyte_test_data();

            // Write as SBYTE_IMG
            let filename_cstr = CString::new(filename).unwrap();
            fits_create_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                &mut status,
            );
            fits_create_img(
                fptr.as_mut().unwrap(),
                SBYTE_IMG,
                naxis as c_int,
                &naxes,
                &mut status,
            );
            fits_write_img(
                fptr.as_mut().unwrap(),
                TSBYTE,
                1,
                nelements as LONGLONG,
                cast_slice(&write_data),
                &mut status,
            );
            fits_close_file(fptr.take().unwrap(), &mut status);

            assert_eq!(status, 0, "Failed to write SBYTE_IMG");

            // Read as SHORT - should succeed
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

                assert_eq!(status, 0, "SBYTE → SHORT should succeed without overflow");

                for i in 0..nelements {
                    assert_eq!(
                        write_data[i] as c_short, read_data[i],
                        "Data mismatch at index {}: SBYTE {} → SHORT {}",
                        i, write_data[i], read_data[i]
                    );
                }
            }

            fits_close_file(fptr.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_sbyte_write_float_read_success() {
        // SBYTE → FLOAT should always succeed
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;
            let write_data = generate_sbyte_test_data();

            // Write as SBYTE_IMG
            let filename_cstr = CString::new(filename).unwrap();
            fits_create_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                &mut status,
            );
            fits_create_img(
                fptr.as_mut().unwrap(),
                SBYTE_IMG,
                naxis as c_int,
                &naxes,
                &mut status,
            );
            fits_write_img(
                fptr.as_mut().unwrap(),
                TSBYTE,
                1,
                nelements as LONGLONG,
                cast_slice(&write_data),
                &mut status,
            );
            fits_close_file(fptr.take().unwrap(), &mut status);

            assert_eq!(status, 0, "Failed to write SBYTE_IMG");

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

                assert_eq!(status, 0, "SBYTE → FLOAT should succeed without overflow");

                for i in 0..nelements {
                    let expected = write_data[i] as c_float;
                    assert!(
                        floats_close_f32(expected, read_data[i]),
                        "Data mismatch at index {}: SBYTE {} → FLOAT {} (expected {})",
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
    fn test_sbyte_write_long_read_success() {
        // SBYTE → LONG should always succeed
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;
            let write_data = generate_sbyte_test_data();

            // Write as SBYTE_IMG
            let filename_cstr = CString::new(filename).unwrap();
            fits_create_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                &mut status,
            );
            fits_create_img(
                fptr.as_mut().unwrap(),
                SBYTE_IMG,
                naxis as c_int,
                &naxes,
                &mut status,
            );
            fits_write_img(
                fptr.as_mut().unwrap(),
                TSBYTE,
                1,
                nelements as LONGLONG,
                cast_slice(&write_data),
                &mut status,
            );
            fits_close_file(fptr.take().unwrap(), &mut status);

            assert_eq!(status, 0, "Failed to write SBYTE_IMG");

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

                assert_eq!(status, 0, "SBYTE → LONG should succeed without overflow");

                for i in 0..nelements {
                    let expected = write_data[i] as c_long;
                    assert_eq!(
                        read_data[i], expected,
                        "Data mismatch at index {}: SBYTE {} → LONG {} (expected {})",
                        i, write_data[i], read_data[i], expected
                    );
                }
            }

            fits_close_file(fptr.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_sbyte_write_longlong_read_success() {
        // SBYTE → LONGLONG should always succeed
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;
            let write_data = generate_sbyte_test_data();

            // Write as SBYTE_IMG
            let filename_cstr = CString::new(filename).unwrap();
            fits_create_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                &mut status,
            );
            fits_create_img(
                fptr.as_mut().unwrap(),
                SBYTE_IMG,
                naxis as c_int,
                &naxes,
                &mut status,
            );
            fits_write_img(
                fptr.as_mut().unwrap(),
                TSBYTE,
                1,
                nelements as LONGLONG,
                cast_slice(&write_data),
                &mut status,
            );
            fits_close_file(fptr.take().unwrap(), &mut status);

            assert_eq!(status, 0, "Failed to write SBYTE_IMG");

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
                    "SBYTE → LONGLONG should succeed without overflow"
                );

                for i in 0..nelements {
                    let expected = write_data[i] as c_longlong;
                    assert_eq!(
                        read_data[i], expected,
                        "Data mismatch at index {}: SBYTE {} → LONGLONG {} (expected {})",
                        i, write_data[i], read_data[i], expected
                    );
                }
            }

            fits_close_file(fptr.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_sbyte_write_double_read_success() {
        // SBYTE → DOUBLE should always succeed
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;
            let write_data = generate_sbyte_test_data();

            // Write as SBYTE_IMG
            let filename_cstr = CString::new(filename).unwrap();
            fits_create_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                &mut status,
            );
            fits_create_img(
                fptr.as_mut().unwrap(),
                SBYTE_IMG,
                naxis as c_int,
                &naxes,
                &mut status,
            );
            fits_write_img(
                fptr.as_mut().unwrap(),
                TSBYTE,
                1,
                nelements as LONGLONG,
                cast_slice(&write_data),
                &mut status,
            );
            fits_close_file(fptr.take().unwrap(), &mut status);

            assert_eq!(status, 0, "Failed to write SBYTE_IMG");

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

                assert_eq!(status, 0, "SBYTE → DOUBLE should succeed without overflow");

                for i in 0..nelements {
                    let expected = write_data[i] as c_double;
                    assert!(
                        floats_close_f64(expected, read_data[i]),
                        "Data mismatch at index {}: SBYTE {} → DOUBLE {} (expected {})",
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
    fn test_sbyte_write_ulong_read_overflow() {
        // SBYTE with negatives will overflow ULONG (0 to ~4.3 billion)
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;
            let write_data = generate_sbyte_test_data(); // Contains negatives

            // Write as SBYTE_IMG
            let filename_cstr = CString::new(filename).unwrap();
            fits_create_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                &mut status,
            );
            fits_create_img(
                fptr.as_mut().unwrap(),
                SBYTE_IMG,
                naxis as c_int,
                &naxes,
                &mut status,
            );
            fits_write_img(
                fptr.as_mut().unwrap(),
                TSBYTE,
                1,
                nelements as LONGLONG,
                cast_slice(&write_data),
                &mut status,
            );
            fits_close_file(fptr.take().unwrap(), &mut status);

            assert_eq!(status, 0, "Failed to write SBYTE_IMG");

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
                    "Expected NUM_OVERFLOW (412) for SBYTE with negatives → ULONG"
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
                        "Value at index {} should be clamped: SBYTE {} → ULONG {} (got {})",
                        i, write_data[i], expected, read_data[i]
                    );
                }
            }

            fits_close_file(fptr.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_sbyte_write_ulonglong_read_overflow() {
        // SBYTE with negatives will overflow ULONGLONG (0 to ~18 quintillion)
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;
            let write_data = generate_sbyte_test_data(); // Contains negatives

            // Write as SBYTE_IMG
            let filename_cstr = CString::new(filename).unwrap();
            fits_create_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                &mut status,
            );
            fits_create_img(
                fptr.as_mut().unwrap(),
                SBYTE_IMG,
                naxis as c_int,
                &naxes,
                &mut status,
            );
            fits_write_img(
                fptr.as_mut().unwrap(),
                TSBYTE,
                1,
                nelements as LONGLONG,
                cast_slice(&write_data),
                &mut status,
            );
            fits_close_file(fptr.take().unwrap(), &mut status);

            assert_eq!(status, 0, "Failed to write SBYTE_IMG");

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
                    "Expected NUM_OVERFLOW (412) for SBYTE with negatives → ULONGLONG"
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
                        "Value at index {} should be clamped: SBYTE {} → ULONGLONG {} (got {})",
                        i, write_data[i], expected, read_data[i]
                    );
                }
            }

            fits_close_file(fptr.take().unwrap(), &mut status);
        });
    }
}
