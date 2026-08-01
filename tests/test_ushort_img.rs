// Ported from the CFITSIO C test programs; indexed loops and the
// operator-named test fns mirror the originals.
#![allow(clippy::needless_range_loop)]

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
        LONGLONG, READWRITE, TBYTE, TDOUBLE, TFLOAT, TLONG, TLONGLONG, TSBYTE, TSHORT, TULONG,
        TULONGLONG, TUSHORT, USHORT_IMG, fitsfile,
    };
    use std::ffi::CString;

    const IMAGE_WIDTH: c_long = 50;
    const IMAGE_HEIGHT: c_long = 50;

    /// Generate USHORT_IMG test data (0-65535 range, all positive)
    fn generate_ushort_test_data() -> Vec<u16> {
        let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;
        let mut data = vec![0; nelements];
        for j in 0..IMAGE_HEIGHT {
            for i in 0..IMAGE_WIDTH {
                let index = (j * IMAGE_WIDTH + i) as usize;
                data[index] = ((i + j) * 200) as u16; // 0 to ~19,600 range
            }
        }
        data
    }

    #[test]
    fn test_ushort_write_ushort_read() {
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;
            let write_data = generate_ushort_test_data();

            // Write as USHORT_IMG
            let filename_cstr = CString::new(filename).unwrap();
            fits_create_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                &mut status,
            );
            fits_create_img(
                fptr.as_mut().unwrap(),
                USHORT_IMG,
                naxis as c_int,
                &naxes,
                &mut status,
            );
            fits_write_img(
                fptr.as_mut().unwrap(),
                TUSHORT,
                1,
                nelements as LONGLONG,
                cast_slice(&write_data),
                &mut status,
            );
            fits_close_file(fptr.take().unwrap(), &mut status);

            assert_eq!(status, 0, "Failed to write USHORT_IMG");

            // Read back as USHORT
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

                assert_eq!(status, 0, "Failed to read USHORT_IMG as TUSHORT");

                for i in 0..nelements {
                    assert_eq!(write_data[i], read_data[i], "Data mismatch at index {}", i);
                }
            }

            fits_close_file(fptr.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_ushort_write_byte_read_overflow() {
        // USHORT data (0-19600) will overflow BYTE range (0-255)
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;
            let write_data = generate_ushort_test_data();

            // Write as USHORT_IMG
            let filename_cstr = CString::new(filename).unwrap();
            fits_create_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                &mut status,
            );
            fits_create_img(
                fptr.as_mut().unwrap(),
                USHORT_IMG,
                naxis as c_int,
                &naxes,
                &mut status,
            );
            fits_write_img(
                fptr.as_mut().unwrap(),
                TUSHORT,
                1,
                nelements as LONGLONG,
                cast_slice(&write_data),
                &mut status,
            );
            fits_close_file(fptr.take().unwrap(), &mut status);

            assert_eq!(status, 0, "Failed to write USHORT_IMG");

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

                assert_eq!(
                    status, 412,
                    "Expected NUM_OVERFLOW (412) for USHORT > 255 → BYTE"
                );

                // Verify clamping behavior
                for i in 0..nelements {
                    let expected = if write_data[i] > 255 {
                        255
                    } else {
                        write_data[i] as c_uchar
                    };
                    assert_eq!(
                        read_data[i], expected,
                        "Value at index {} should be clamped: USHORT {} → BYTE {} (got {})",
                        i, write_data[i], expected, read_data[i]
                    );
                }
            }

            fits_close_file(fptr.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_ushort_write_short_read_overflow() {
        // USHORT max value (~19600) fits in SHORT, but let's test with larger values
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [10, 10];
            let nelements = 100;

            // Create USHORT data that exceeds SHORT max (32767)
            let mut write_data: Vec<u16> = vec![0; nelements];
            write_data[0] = 32766; // Valid for SHORT
            write_data[1] = 32767; // Max SHORT
            write_data[2] = 32768; // Overflow SHORT
            write_data[3] = 40000; // Overflow SHORT
            write_data[4] = 65535; // Max USHORT, overflow SHORT
            for i in 5..nelements {
                write_data[i] = (i * 100) as u16; // Mix of valid values
            }

            // Write as USHORT_IMG
            let filename_cstr = CString::new(filename).unwrap();
            fits_create_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                &mut status,
            );
            fits_create_img(
                fptr.as_mut().unwrap(),
                USHORT_IMG,
                naxis as c_int,
                &naxes,
                &mut status,
            );
            fits_write_img(
                fptr.as_mut().unwrap(),
                TUSHORT,
                1,
                nelements as LONGLONG,
                cast_slice(&write_data),
                &mut status,
            );
            fits_close_file(fptr.take().unwrap(), &mut status);

            assert_eq!(status, 0, "Failed to write USHORT_IMG");

            // Read as SHORT - should get overflow for values > 32767
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

                assert_eq!(
                    status, 412,
                    "Expected NUM_OVERFLOW (412) for USHORT > 32767 → SHORT"
                );

                // Verify specific clamping
                assert_eq!(read_data[0], 32766, "Valid value should convert exactly");
                assert_eq!(read_data[1], 32767, "Max SHORT should convert exactly");
                assert_eq!(read_data[2], 32767, "32768 should clamp to 32767");
                assert_eq!(read_data[3], 32767, "40000 should clamp to 32767");
                assert_eq!(read_data[4], 32767, "65535 should clamp to 32767");
            }

            fits_close_file(fptr.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_ushort_write_long_read_success() {
        // USHORT → LONG should always succeed (USHORT max 65535 < LONG max)
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;
            let write_data = generate_ushort_test_data();

            // Write as USHORT_IMG
            let filename_cstr = CString::new(filename).unwrap();
            fits_create_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                &mut status,
            );
            fits_create_img(
                fptr.as_mut().unwrap(),
                USHORT_IMG,
                naxis as c_int,
                &naxes,
                &mut status,
            );
            fits_write_img(
                fptr.as_mut().unwrap(),
                TUSHORT,
                1,
                nelements as LONGLONG,
                cast_slice(&write_data),
                &mut status,
            );
            fits_close_file(fptr.take().unwrap(), &mut status);

            assert_eq!(status, 0, "Failed to write USHORT_IMG");

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

                assert_eq!(status, 0, "USHORT → LONG should succeed without overflow");

                for i in 0..nelements {
                    assert_eq!(
                        write_data[i] as c_long, read_data[i],
                        "Data mismatch at index {}: USHORT {} → LONG {}",
                        i, write_data[i], read_data[i]
                    );
                }
            }

            fits_close_file(fptr.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_ushort_write_float_read_success() {
        // USHORT → FLOAT should always succeed
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;
            let write_data = generate_ushort_test_data();

            // Write as USHORT_IMG
            let filename_cstr = CString::new(filename).unwrap();
            fits_create_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                &mut status,
            );
            fits_create_img(
                fptr.as_mut().unwrap(),
                USHORT_IMG,
                naxis as c_int,
                &naxes,
                &mut status,
            );
            fits_write_img(
                fptr.as_mut().unwrap(),
                TUSHORT,
                1,
                nelements as LONGLONG,
                cast_slice(&write_data),
                &mut status,
            );
            fits_close_file(fptr.take().unwrap(), &mut status);

            assert_eq!(status, 0, "Failed to write USHORT_IMG");

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

                assert_eq!(status, 0, "USHORT → FLOAT should succeed without overflow");

                for i in 0..nelements {
                    let expected = write_data[i] as c_float;
                    assert!(
                        floats_close_f32(expected, read_data[i]),
                        "Data mismatch at index {}: USHORT {} → FLOAT {} (expected {})",
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
    fn test_ushort_write_longlong_read_success() {
        // USHORT → LONGLONG should always succeed (USHORT range fits in LONGLONG)
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;
            let write_data = generate_ushort_test_data();

            // Write as USHORT_IMG
            let filename_cstr = CString::new(filename).unwrap();
            fits_create_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                &mut status,
            );
            fits_create_img(
                fptr.as_mut().unwrap(),
                USHORT_IMG,
                naxis as c_int,
                &naxes,
                &mut status,
            );
            fits_write_img(
                fptr.as_mut().unwrap(),
                TUSHORT,
                1,
                nelements as LONGLONG,
                cast_slice(&write_data),
                &mut status,
            );
            fits_close_file(fptr.take().unwrap(), &mut status);

            assert_eq!(status, 0, "Failed to write USHORT_IMG");

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
                    "USHORT → LONGLONG should succeed without overflow"
                );

                for i in 0..nelements {
                    assert_eq!(
                        write_data[i] as c_longlong, read_data[i],
                        "Data mismatch at index {}: USHORT {} → LONGLONG {}",
                        i, write_data[i], read_data[i]
                    );
                }
            }

            fits_close_file(fptr.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_ushort_write_double_read_success() {
        // USHORT → DOUBLE should always succeed
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;
            let write_data = generate_ushort_test_data();

            // Write as USHORT_IMG
            let filename_cstr = CString::new(filename).unwrap();
            fits_create_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                &mut status,
            );
            fits_create_img(
                fptr.as_mut().unwrap(),
                USHORT_IMG,
                naxis as c_int,
                &naxes,
                &mut status,
            );
            fits_write_img(
                fptr.as_mut().unwrap(),
                TUSHORT,
                1,
                nelements as LONGLONG,
                cast_slice(&write_data),
                &mut status,
            );
            fits_close_file(fptr.take().unwrap(), &mut status);

            assert_eq!(status, 0, "Failed to write USHORT_IMG");

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

                assert_eq!(status, 0, "USHORT → DOUBLE should succeed without overflow");

                for i in 0..nelements {
                    let expected = write_data[i] as c_double;
                    assert!(
                        floats_close_f64(expected, read_data[i]),
                        "Data mismatch at index {}: USHORT {} → DOUBLE {} (expected {})",
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
    fn test_ushort_write_sbyte_read_overflow() {
        // USHORT → SBYTE should overflow for values > 127
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;
            let write_data = generate_ushort_test_data(); // Range: 0 to ~19,600

            // Write as USHORT_IMG
            let filename_cstr = CString::new(filename).unwrap();
            fits_create_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                &mut status,
            );
            fits_create_img(
                fptr.as_mut().unwrap(),
                USHORT_IMG,
                naxis as c_int,
                &naxes,
                &mut status,
            );
            fits_write_img(
                fptr.as_mut().unwrap(),
                TUSHORT,
                1,
                nelements as LONGLONG,
                cast_slice(&write_data),
                &mut status,
            );
            fits_close_file(fptr.take().unwrap(), &mut status);

            assert_eq!(status, 0, "Failed to write USHORT_IMG");

            // Read as SBYTE - should get overflow for values > 127
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

                assert_eq!(
                    status, 412,
                    "Expected NUM_OVERFLOW (412) for USHORT → SBYTE"
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
                        "Value at index {} should be clamped: USHORT {} → SBYTE {}",
                        i, write_data[i], expected
                    );
                }
            }

            fits_close_file(fptr.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_ushort_write_ulong_read_success() {
        // USHORT → ULONG should always succeed (USHORT range fits in ULONG)
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;
            let write_data = generate_ushort_test_data();

            // Write as USHORT_IMG
            let filename_cstr = CString::new(filename).unwrap();
            fits_create_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                &mut status,
            );
            fits_create_img(
                fptr.as_mut().unwrap(),
                USHORT_IMG,
                naxis as c_int,
                &naxes,
                &mut status,
            );
            fits_write_img(
                fptr.as_mut().unwrap(),
                TUSHORT,
                1,
                nelements as LONGLONG,
                cast_slice(&write_data),
                &mut status,
            );
            fits_close_file(fptr.take().unwrap(), &mut status);

            assert_eq!(status, 0, "Failed to write USHORT_IMG");

            // Read as ULONG - should succeed
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

                assert_eq!(status, 0, "USHORT → ULONG should succeed without overflow");

                for i in 0..nelements {
                    assert_eq!(
                        write_data[i] as c_ulong, read_data[i],
                        "Data mismatch at index {}: USHORT {} → ULONG {}",
                        i, write_data[i], read_data[i]
                    );
                }
            }

            fits_close_file(fptr.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_ushort_write_ulonglong_read_success() {
        // USHORT → ULONGLONG should always succeed (USHORT range fits in ULONGLONG)
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;
            let write_data = generate_ushort_test_data();

            // Write as USHORT_IMG
            let filename_cstr = CString::new(filename).unwrap();
            fits_create_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                &mut status,
            );
            fits_create_img(
                fptr.as_mut().unwrap(),
                USHORT_IMG,
                naxis as c_int,
                &naxes,
                &mut status,
            );
            fits_write_img(
                fptr.as_mut().unwrap(),
                TUSHORT,
                1,
                nelements as LONGLONG,
                cast_slice(&write_data),
                &mut status,
            );
            fits_close_file(fptr.take().unwrap(), &mut status);

            assert_eq!(status, 0, "Failed to write USHORT_IMG");

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
                    "USHORT → ULONGLONG should succeed without overflow"
                );

                for i in 0..nelements {
                    assert_eq!(
                        write_data[i] as c_ulonglong, read_data[i],
                        "Data mismatch at index {}: USHORT {} → ULONGLONG {}",
                        i, write_data[i], read_data[i]
                    );
                }
            }

            fits_close_file(fptr.take().unwrap(), &mut status);
        });
    }
}
