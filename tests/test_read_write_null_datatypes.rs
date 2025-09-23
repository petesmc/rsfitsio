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
    use rsfitsio::{NullValue, cs};
    use std::ffi::CString;

    // Test dimensions for all images
    const IMAGE_WIDTH: c_long = 50;
    const IMAGE_HEIGHT: c_long = 50;

    #[test]
    fn test_byte_img_null_read_write() {
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;
            const NULL_VALUE: c_uchar = 99;
            const SUBSTITUTE_VALUE: c_uchar = 101;

            // Create test data with nulls on diagonal
            let mut write_data: Vec<c_uchar> = vec![0; nelements];
            for j in 0..IMAGE_HEIGHT {
                for i in 0..IMAGE_WIDTH {
                    let index = (j * IMAGE_WIDTH + i) as usize;
                    if i == j {
                        write_data[index] = NULL_VALUE; // Set diagonal to null value
                    } else {
                        write_data[index] = ((i + j) % 98) as c_uchar; // Avoid 99
                    }
                }
            }

            // Create FITS file and write image
            let filename_cstr = CString::new(filename).unwrap();
            unsafe {
                fits_create_file(
                    &mut fptr,
                    cast_slice(filename_cstr.to_bytes_with_nul()),
                    &mut status,
                );
            }
            assert_eq!(status, 0, "Failed to create file");

            if let Some(ref mut fptr_box) = fptr {
                // Create image HDU
                unsafe {
                    fits_create_img(fptr_box, BYTE_IMG, naxis as c_int, &naxes, &mut status);
                }
                assert_eq!(status, 0, "Failed to create BYTE_IMG");

                // Set BLANK keyword to define null value for integer image
                unsafe {
                    fits_update_key_lng(
                        fptr_box,
                        cs!(c"BLANK"),
                        NULL_VALUE as c_long,
                        None,
                        &mut status,
                    );
                }
                assert_eq!(status, 0, "Failed to set BLANK keyword");

                // Write image data normally - pixels with value NULL_VALUE should be treated as nulls
                unsafe {
                    fits_write_img(
                        fptr_box,
                        TBYTE,
                        1,
                        nelements as LONGLONG,
                        cast_slice(&write_data),
                        &mut status,
                    );
                }
                assert_eq!(status, 0, "Failed to write BYTE_IMG data");
            }

            // Close and reopen for reading
            unsafe {
                fits_close_file(fptr.take().unwrap(), &mut status);
            }
            assert_eq!(status, 0, "Failed to close file");

            // Read back the data with null substitution
            fptr = None;
            unsafe {
                fits_open_file(
                    &mut fptr,
                    cast_slice(filename_cstr.to_bytes_with_nul()),
                    READWRITE,
                    &mut status,
                );
            }
            assert_eq!(status, 0, "Failed to open file for reading");

            if let Some(ref mut fptr_box) = fptr {
                let mut read_data: Vec<c_uchar> = vec![0; nelements];
                let mut anynull: c_int = 0;

                unsafe {
                    fits_read_img(
                        fptr_box,
                        TBYTE,
                        1,
                        nelements as LONGLONG,
                        Some(NullValue::UByte(SUBSTITUTE_VALUE)), // Don't substitute null values yet
                        cast_slice_mut(&mut read_data),
                        Some(&mut anynull),
                        &mut status,
                    );
                }
                assert_eq!(status, 0, "Failed to read BYTE_IMG data");
                assert_eq!(anynull, 1, "Should have detected null values");

                // Verify data - diagonal should be substituted, others unchanged
                for j in 0..IMAGE_HEIGHT {
                    for i in 0..IMAGE_WIDTH {
                        let index = (j * IMAGE_WIDTH + i) as usize;
                        if i == j {
                            assert_eq!(
                                read_data[index], SUBSTITUTE_VALUE,
                                "Null value not substituted at diagonal position ({}, {})",
                                i, j
                            );
                        } else {
                            assert_eq!(
                                write_data[index], read_data[index],
                                "Data mismatch at index {}",
                                index
                            );
                        }
                    }
                }
            }

            // Close file
            unsafe {
                fits_close_file(fptr.take().unwrap(), &mut status);
            }
            assert_eq!(status, 0, "Failed to close file after reading");
        });
    }

    #[test]
    fn test_short_img_null_read_write() {
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;
            const NULL_VALUE: c_short = 99;
            const SUBSTITUTE_VALUE: c_short = 101;

            // Create test data with nulls on diagonal
            let mut write_data: Vec<c_short> = vec![0; nelements];
            for j in 0..IMAGE_HEIGHT {
                for i in 0..IMAGE_WIDTH {
                    let index = (j * IMAGE_WIDTH + i) as usize;
                    if i == j {
                        write_data[index] = NULL_VALUE;
                    } else {
                        let val = ((i - 50) * (j - 50)) as c_short;
                        // Avoid NULL_VALUE in non-diagonal positions
                        write_data[index] = if val == NULL_VALUE {
                            NULL_VALUE + 1
                        } else {
                            val
                        };
                    }
                }
            }

            // Create FITS file and write image
            let filename_cstr = CString::new(filename).unwrap();
            unsafe {
                fits_create_file(
                    &mut fptr,
                    cast_slice(filename_cstr.to_bytes_with_nul()),
                    &mut status,
                );
            }
            assert_eq!(status, 0, "Failed to create file");

            if let Some(ref mut fptr_box) = fptr {
                unsafe {
                    fits_create_img(fptr_box, SHORT_IMG, naxis as c_int, &naxes, &mut status);
                }
                assert_eq!(status, 0, "Failed to create SHORT_IMG");

                // Set BLANK keyword to define null value for integer image
                unsafe {
                    fits_update_key_lng(
                        fptr_box,
                        cs!(c"BLANK"),
                        NULL_VALUE as c_long,
                        None,
                        &mut status,
                    );
                }
                assert_eq!(status, 0, "Failed to set BLANK keyword");

                unsafe {
                    fits_write_imgnull_sht(
                        fptr_box,
                        1, // group
                        1, // firstelem
                        nelements as LONGLONG,
                        &write_data,
                        NULL_VALUE,
                        &mut status,
                    );
                }
                assert_eq!(status, 0, "Failed to write SHORT_IMG data with nulls");
            }

            // Close and reopen for reading
            unsafe {
                fits_close_file(fptr.take().unwrap(), &mut status);
            }
            assert_eq!(status, 0, "Failed to close file");

            fptr = None;
            unsafe {
                fits_open_file(
                    &mut fptr,
                    cast_slice(filename_cstr.to_bytes_with_nul()),
                    READWRITE,
                    &mut status,
                );
            }
            assert_eq!(status, 0, "Failed to open file for reading");

            if let Some(ref mut fptr_box) = fptr {
                let mut read_data: Vec<c_short> = vec![0; nelements];
                let mut anynull: c_int = 0;

                unsafe {
                    fits_read_img(
                        fptr_box,
                        TSHORT,
                        1,
                        nelements as LONGLONG,
                        Some(NullValue::Short(SUBSTITUTE_VALUE)),
                        cast_slice_mut(&mut read_data),
                        Some(&mut anynull),
                        &mut status,
                    );
                }
                assert_eq!(status, 0, "Failed to read SHORT_IMG data");
                assert_eq!(anynull, 1, "Should have detected null values");

                // Verify data
                for j in 0..IMAGE_HEIGHT {
                    for i in 0..IMAGE_WIDTH {
                        let index = (j * IMAGE_WIDTH + i) as usize;
                        if i == j {
                            assert_eq!(
                                read_data[index], SUBSTITUTE_VALUE,
                                "Null value not substituted at diagonal position ({}, {})",
                                i, j
                            );
                        } else {
                            assert_eq!(
                                write_data[index], read_data[index],
                                "Data mismatch at index {}",
                                index
                            );
                        }
                    }
                }
            }

            unsafe {
                fits_close_file(fptr.take().unwrap(), &mut status);
            }
            assert_eq!(status, 0, "Failed to close file after reading");
        });
    }

    #[test]
    fn test_ushort_img_null_read_write() {
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;
            const NULL_VALUE: u16 = 99;
            const SUBSTITUTE_VALUE: u16 = 101;

            // Create test data with nulls on diagonal
            let mut write_data: Vec<u16> = vec![0; nelements];
            for j in 0..IMAGE_HEIGHT {
                for i in 0..IMAGE_WIDTH {
                    let index = (j * IMAGE_WIDTH + i) as usize;
                    if i == j {
                        write_data[index] = NULL_VALUE;
                    } else {
                        let val = ((i + j) * 200) as u16;
                        // Avoid NULL_VALUE in non-diagonal positions
                        write_data[index] = if val == NULL_VALUE {
                            NULL_VALUE + 1
                        } else {
                            val
                        };
                    }
                }
            }

            // Create FITS file and write image
            let filename_cstr = CString::new(filename).unwrap();
            unsafe {
                fits_create_file(
                    &mut fptr,
                    cast_slice(filename_cstr.to_bytes_with_nul()),
                    &mut status,
                );
            }
            assert_eq!(status, 0, "Failed to create file");

            if let Some(ref mut fptr_box) = fptr {
                unsafe {
                    fits_create_img(fptr_box, USHORT_IMG, naxis as c_int, &naxes, &mut status);
                }
                assert_eq!(status, 0, "Failed to create USHORT_IMG");

                // Set BLANK keyword to define null value for integer image
                unsafe {
                    fits_update_key_lng(
                        fptr_box,
                        cs!(c"BLANK"),
                        NULL_VALUE as c_long,
                        None,
                        &mut status,
                    );
                }
                assert_eq!(status, 0, "Failed to set BLANK keyword");

                unsafe {
                    fits_write_imgnull_usht(
                        fptr_box,
                        1, // group
                        1, // firstelem
                        nelements as LONGLONG,
                        &write_data,
                        NULL_VALUE,
                        &mut status,
                    );
                }
                assert_eq!(status, 0, "Failed to write USHORT_IMG data with nulls");
            }

            // Close and reopen for reading
            unsafe {
                fits_close_file(fptr.take().unwrap(), &mut status);
            }
            assert_eq!(status, 0, "Failed to close file");

            fptr = None;
            unsafe {
                fits_open_file(
                    &mut fptr,
                    cast_slice(filename_cstr.to_bytes_with_nul()),
                    READWRITE,
                    &mut status,
                );
            }
            assert_eq!(status, 0, "Failed to open file for reading");

            if let Some(ref mut fptr_box) = fptr {
                let mut read_data: Vec<u16> = vec![0; nelements];
                let mut anynull: c_int = 0;

                unsafe {
                    fits_read_img(
                        fptr_box,
                        TUSHORT,
                        1,
                        nelements as LONGLONG,
                        Some(NullValue::UShort(SUBSTITUTE_VALUE)),
                        cast_slice_mut(&mut read_data),
                        Some(&mut anynull),
                        &mut status,
                    );
                }
                assert_eq!(status, 0, "Failed to read USHORT_IMG data");
                assert_eq!(anynull, 1, "Should have detected null values");

                // Verify data
                for j in 0..IMAGE_HEIGHT {
                    for i in 0..IMAGE_WIDTH {
                        let index = (j * IMAGE_WIDTH + i) as usize;
                        if i == j {
                            assert_eq!(
                                read_data[index], SUBSTITUTE_VALUE,
                                "Null value not substituted at diagonal position ({}, {})",
                                i, j
                            );
                        } else {
                            assert_eq!(
                                write_data[index], read_data[index],
                                "Data mismatch at index {}",
                                index
                            );
                        }
                    }
                }
            }

            unsafe {
                fits_close_file(fptr.take().unwrap(), &mut status);
            }
            assert_eq!(status, 0, "Failed to close file after reading");
        });
    }

    #[test]
    fn test_long_img_null_read_write() {
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;
            const NULL_VALUE: c_long = 99;
            const SUBSTITUTE_VALUE: c_long = 101;

            // Create test data with nulls on diagonal
            let mut write_data: Vec<c_long> = vec![0; nelements];
            for j in 0..IMAGE_HEIGHT {
                for i in 0..IMAGE_WIDTH {
                    let index = (j * IMAGE_WIDTH + i) as usize;
                    if i == j {
                        write_data[index] = NULL_VALUE;
                    } else {
                        let val = ((i - 25) * 1000 + (j - 25) * 10) as c_long;
                        // Avoid NULL_VALUE in non-diagonal positions
                        write_data[index] = if val == NULL_VALUE {
                            NULL_VALUE + 1
                        } else {
                            val
                        };
                    }
                }
            }

            // Create FITS file and write image
            let filename_cstr = CString::new(filename).unwrap();
            unsafe {
                fits_create_file(
                    &mut fptr,
                    cast_slice(filename_cstr.to_bytes_with_nul()),
                    &mut status,
                );
            }
            assert_eq!(status, 0, "Failed to create file");

            if let Some(ref mut fptr_box) = fptr {
                unsafe {
                    fits_create_img(fptr_box, LONG_IMG, naxis as c_int, &naxes, &mut status);
                }
                assert_eq!(status, 0, "Failed to create LONG_IMG");

                // Set BLANK keyword to define null value for integer image
                unsafe {
                    fits_update_key_lng(
                        fptr_box,
                        cs!(c"BLANK"),
                        NULL_VALUE as c_long,
                        None,
                        &mut status,
                    );
                }
                assert_eq!(status, 0, "Failed to set BLANK keyword");

                unsafe {
                    fits_write_imgnull_lng(
                        fptr_box,
                        1, // group
                        1, // firstelem
                        nelements as LONGLONG,
                        &write_data,
                        NULL_VALUE,
                        &mut status,
                    );
                }
                assert_eq!(status, 0, "Failed to write LONG_IMG data with nulls");
            }

            // Close and reopen for reading
            unsafe {
                fits_close_file(fptr.take().unwrap(), &mut status);
            }
            assert_eq!(status, 0, "Failed to close file");

            fptr = None;
            unsafe {
                fits_open_file(
                    &mut fptr,
                    cast_slice(filename_cstr.to_bytes_with_nul()),
                    READWRITE,
                    &mut status,
                );
            }
            assert_eq!(status, 0, "Failed to open file for reading");

            if let Some(ref mut fptr_box) = fptr {
                let mut read_data: Vec<c_long> = vec![0; nelements];
                let mut anynull: c_int = 0;

                unsafe {
                    fits_read_img(
                        fptr_box,
                        TLONG,
                        1,
                        nelements as LONGLONG,
                        Some(NullValue::Long(SUBSTITUTE_VALUE)),
                        cast_slice_mut(&mut read_data),
                        Some(&mut anynull),
                        &mut status,
                    );
                }
                assert_eq!(status, 0, "Failed to read LONG_IMG data");
                assert_eq!(anynull, 1, "Should have detected null values");

                // Verify data
                for j in 0..IMAGE_HEIGHT {
                    for i in 0..IMAGE_WIDTH {
                        let index = (j * IMAGE_WIDTH + i) as usize;
                        if i == j {
                            assert_eq!(
                                read_data[index], SUBSTITUTE_VALUE,
                                "Null value not substituted at diagonal position ({}, {})",
                                i, j
                            );
                        } else {
                            assert_eq!(
                                write_data[index], read_data[index],
                                "Data mismatch at index {}",
                                index
                            );
                        }
                    }
                }
            }

            unsafe {
                fits_close_file(fptr.take().unwrap(), &mut status);
            }
            assert_eq!(status, 0, "Failed to close file after reading");
        });
    }

    #[test]
    fn test_longlong_img_null_read_write() {
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;
            const NULL_VALUE: c_longlong = 99;
            const SUBSTITUTE_VALUE: c_longlong = 101;

            // Create test data with nulls on diagonal
            let mut write_data: Vec<c_longlong> = vec![0; nelements];
            for j in 0..IMAGE_HEIGHT {
                for i in 0..IMAGE_WIDTH {
                    let index = (j * IMAGE_WIDTH + i) as usize;
                    if i == j {
                        write_data[index] = NULL_VALUE;
                    } else {
                        let val = ((i - 25) as c_longlong * 100000 + (j - 25) as c_longlong * 1000);
                        // Avoid NULL_VALUE in non-diagonal positions
                        write_data[index] = if val == NULL_VALUE {
                            NULL_VALUE + 1
                        } else {
                            val
                        };
                    }
                }
            }

            // Create FITS file and write image
            let filename_cstr = CString::new(filename).unwrap();
            unsafe {
                fits_create_file(
                    &mut fptr,
                    cast_slice(filename_cstr.to_bytes_with_nul()),
                    &mut status,
                );
            }
            assert_eq!(status, 0, "Failed to create file");

            if let Some(ref mut fptr_box) = fptr {
                unsafe {
                    fits_create_img(fptr_box, LONGLONG_IMG, naxis as c_int, &naxes, &mut status);
                }
                assert_eq!(status, 0, "Failed to create LONGLONG_IMG");

                // Set BLANK keyword to define null value for integer image
                unsafe {
                    fits_update_key_lng(
                        fptr_box,
                        cs!(c"BLANK"),
                        NULL_VALUE as c_long,
                        None,
                        &mut status,
                    );
                }
                assert_eq!(status, 0, "Failed to set BLANK keyword");

                unsafe {
                    fits_write_imgnull_lnglng(
                        fptr_box,
                        1, // group
                        1, // firstelem
                        nelements as LONGLONG,
                        &write_data,
                        NULL_VALUE,
                        &mut status,
                    );
                }
                assert_eq!(status, 0, "Failed to write LONGLONG_IMG data with nulls");
            }

            // Close and reopen for reading
            unsafe {
                fits_close_file(fptr.take().unwrap(), &mut status);
            }
            assert_eq!(status, 0, "Failed to close file");

            fptr = None;
            unsafe {
                fits_open_file(
                    &mut fptr,
                    cast_slice(filename_cstr.to_bytes_with_nul()),
                    READWRITE,
                    &mut status,
                );
            }
            assert_eq!(status, 0, "Failed to open file for reading");

            if let Some(ref mut fptr_box) = fptr {
                let mut read_data: Vec<c_longlong> = vec![0; nelements];
                let mut anynull: c_int = 0;

                unsafe {
                    fits_read_img(
                        fptr_box,
                        TLONGLONG,
                        1,
                        nelements as LONGLONG,
                        Some(NullValue::LONGLONG(SUBSTITUTE_VALUE)),
                        cast_slice_mut(&mut read_data),
                        Some(&mut anynull),
                        &mut status,
                    );
                }
                assert_eq!(status, 0, "Failed to read LONGLONG_IMG data");
                assert_eq!(anynull, 1, "Should have detected null values");

                // Verify data
                for j in 0..IMAGE_HEIGHT {
                    for i in 0..IMAGE_WIDTH {
                        let index = (j * IMAGE_WIDTH + i) as usize;
                        if i == j {
                            assert_eq!(
                                read_data[index], SUBSTITUTE_VALUE,
                                "Null value not substituted at diagonal position ({}, {})",
                                i, j
                            );
                        } else {
                            assert_eq!(
                                write_data[index], read_data[index],
                                "Data mismatch at index {}",
                                index
                            );
                        }
                    }
                }
            }

            unsafe {
                fits_close_file(fptr.take().unwrap(), &mut status);
            }
            assert_eq!(status, 0, "Failed to close file after reading");
        });
    }

    #[test]
    fn test_float_img_null_read_write() {
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;
            const NULL_VALUE: c_float = f32::NAN;
            const SUBSTITUTE_VALUE: c_float = 101.0;

            // Create test data with nulls on diagonal
            let mut write_data: Vec<c_float> = vec![0.0; nelements];
            for j in 0..IMAGE_HEIGHT {
                for i in 0..IMAGE_WIDTH {
                    let index = (j * IMAGE_WIDTH + i) as usize;
                    if i == j {
                        write_data[index] = f32::NAN; // Use NaN for null values
                    } else {
                        write_data[index] = ((i - 25) as f32 * 10.0 + (j - 25) as f32 * 1.0);
                    }
                }
            }

            // Create FITS file and write image
            let filename_cstr = CString::new(filename).unwrap();
            unsafe {
                fits_create_file(
                    &mut fptr,
                    cast_slice(filename_cstr.to_bytes_with_nul()),
                    &mut status,
                );
            }
            assert_eq!(status, 0, "Failed to create file");

            if let Some(ref mut fptr_box) = fptr {
                unsafe {
                    fits_create_img(fptr_box, FLOAT_IMG, naxis as c_int, &naxes, &mut status);
                }
                assert_eq!(status, 0, "Failed to create FLOAT_IMG");

                unsafe {
                    fits_write_imgnull_flt(
                        fptr_box,
                        1, // group
                        1, // firstelem
                        nelements as LONGLONG,
                        &write_data,
                        NULL_VALUE,
                        &mut status,
                    );
                }
                assert_eq!(status, 0, "Failed to write FLOAT_IMG data with nulls");
            }

            // Close and reopen for reading
            unsafe {
                fits_close_file(fptr.take().unwrap(), &mut status);
            }
            assert_eq!(status, 0, "Failed to close file");

            fptr = None;
            unsafe {
                fits_open_file(
                    &mut fptr,
                    cast_slice(filename_cstr.to_bytes_with_nul()),
                    READWRITE,
                    &mut status,
                );
            }
            assert_eq!(status, 0, "Failed to open file for reading");

            if let Some(ref mut fptr_box) = fptr {
                let mut read_data: Vec<c_float> = vec![0.0; nelements];
                let mut anynull: c_int = 0;

                unsafe {
                    fits_read_img(
                        fptr_box,
                        TFLOAT,
                        1,
                        nelements as LONGLONG,
                        Some(NullValue::Float(SUBSTITUTE_VALUE)),
                        cast_slice_mut(&mut read_data),
                        Some(&mut anynull),
                        &mut status,
                    );
                }
                assert_eq!(status, 0, "Failed to read FLOAT_IMG data");
                assert_eq!(anynull, 1, "Should have detected null values");

                // Verify data
                for j in 0..IMAGE_HEIGHT {
                    for i in 0..IMAGE_WIDTH {
                        let index = (j * IMAGE_WIDTH + i) as usize;
                        if i == j {
                            // Diagonal elements had NaN, should be substituted with 0
                            assert!(
                                floats_close_f32(read_data[index], SUBSTITUTE_VALUE),
                                "Null value not substituted at diagonal position ({}, {})",
                                i,
                                j
                            );
                        } else {
                            // Non-diagonal elements should match original values
                            assert!(
                                floats_close_f32(write_data[index], read_data[index]),
                                "Data mismatch at index {}: {} != {}",
                                index,
                                write_data[index],
                                read_data[index]
                            );
                        }
                    }
                }
            }

            unsafe {
                fits_close_file(fptr.take().unwrap(), &mut status);
            }
            assert_eq!(status, 0, "Failed to close file after reading");
        });
    }

    #[test]
    fn test_double_img_null_read_write() {
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;
            const NULL_VALUE: c_double = f64::NAN;
            const SUBSTITUTE_VALUE: c_double = 101.0;

            // Create test data with nulls on diagonal
            let mut write_data: Vec<c_double> = vec![0.0; nelements];
            for j in 0..IMAGE_HEIGHT {
                for i in 0..IMAGE_WIDTH {
                    let index = (j * IMAGE_WIDTH + i) as usize;
                    if i == j {
                        write_data[index] = f64::NAN; // Use NaN for null values
                    } else {
                        write_data[index] = ((i - 25) as f64 * 100.0 + (j - 25) as f64 * 10.0);
                    }
                }
            }

            // Create FITS file and write image
            let filename_cstr = CString::new(filename).unwrap();
            unsafe {
                fits_create_file(
                    &mut fptr,
                    cast_slice(filename_cstr.to_bytes_with_nul()),
                    &mut status,
                );
            }
            assert_eq!(status, 0, "Failed to create file");

            if let Some(ref mut fptr_box) = fptr {
                unsafe {
                    fits_create_img(fptr_box, DOUBLE_IMG, naxis as c_int, &naxes, &mut status);
                }
                assert_eq!(status, 0, "Failed to create DOUBLE_IMG");

                unsafe {
                    fits_write_imgnull_dbl(
                        fptr_box,
                        1, // group
                        1, // firstelem
                        nelements as LONGLONG,
                        &write_data,
                        NULL_VALUE,
                        &mut status,
                    );
                }
                assert_eq!(status, 0, "Failed to write DOUBLE_IMG data with nulls");
            }

            // Close and reopen for reading
            unsafe {
                fits_close_file(fptr.take().unwrap(), &mut status);
            }
            assert_eq!(status, 0, "Failed to close file");

            fptr = None;
            unsafe {
                fits_open_file(
                    &mut fptr,
                    cast_slice(filename_cstr.to_bytes_with_nul()),
                    READWRITE,
                    &mut status,
                );
            }
            assert_eq!(status, 0, "Failed to open file for reading");

            if let Some(ref mut fptr_box) = fptr {
                let mut read_data: Vec<c_double> = vec![0.0; nelements];
                let mut anynull: c_int = 0;

                unsafe {
                    fits_read_img(
                        fptr_box,
                        TDOUBLE,
                        1,
                        nelements as LONGLONG,
                        Some(NullValue::Double(SUBSTITUTE_VALUE)),
                        cast_slice_mut(&mut read_data),
                        Some(&mut anynull),
                        &mut status,
                    );
                }
                assert_eq!(status, 0, "Failed to read DOUBLE_IMG data");
                assert_eq!(anynull, 1, "Should have detected null values");

                // Verify data
                for j in 0..IMAGE_HEIGHT {
                    for i in 0..IMAGE_WIDTH {
                        let index = (j * IMAGE_WIDTH + i) as usize;
                        if i == j {
                            // Diagonal elements had NaN, should be substituted with 0
                            assert!(
                                floats_close_f64(read_data[index], SUBSTITUTE_VALUE),
                                "Null value not substituted at diagonal position ({}, {})",
                                i,
                                j
                            );
                        } else {
                            // Non-diagonal elements should match original values
                            assert!(
                                floats_close_f64(write_data[index], read_data[index]),
                                "Data mismatch at index {}: {} != {}",
                                index,
                                write_data[index],
                                read_data[index]
                            );
                        }
                    }
                }
            }

            unsafe {
                fits_close_file(fptr.take().unwrap(), &mut status);
            }
            assert_eq!(status, 0, "Failed to close file after reading");
        });
    }

    #[test]
    fn test_sbyte_img_null_read_write() {
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;
            const NULL_VALUE: i8 = 127; // Use max i8 value as null
            const SUBSTITUTE_VALUE: i8 = 101;

            // Create test data with nulls on diagonal
            let mut write_data: Vec<i8> = vec![0; nelements];
            for j in 0..IMAGE_HEIGHT {
                for i in 0..IMAGE_WIDTH {
                    let index = (j * IMAGE_WIDTH + i) as usize;
                    if i == j {
                        write_data[index] = NULL_VALUE;
                    } else {
                        // Avoid NULL_VALUE and stay within signed byte range
                        let val = ((i + j) % 98 - 49) as i8;
                        write_data[index] = if val == NULL_VALUE {
                            NULL_VALUE - 1
                        } else {
                            val
                        };
                    }
                }
            }

            // Create FITS file and write image
            let filename_cstr = CString::new(filename).unwrap();
            unsafe {
                fits_create_file(
                    &mut fptr,
                    cast_slice(filename_cstr.to_bytes_with_nul()),
                    &mut status,
                );
            }
            assert_eq!(status, 0, "Failed to create file");

            if let Some(ref mut fptr_box) = fptr {
                unsafe {
                    fits_create_img(fptr_box, SBYTE_IMG, naxis as c_int, &naxes, &mut status);
                }
                assert_eq!(status, 0, "Failed to create SBYTE_IMG");

                // Set BLANK keyword to define null value for integer image
                unsafe {
                    fits_update_key_lng(
                        fptr_box,
                        cs!(c"BLANK"),
                        (NULL_VALUE as i32 + 128) as c_long, // Convert to stored unsigned value
                        None,
                        &mut status,
                    );
                }
                assert_eq!(status, 0, "Failed to set BLANK keyword");

                unsafe {
                    fits_write_imgnull_sbyt(
                        fptr_box,
                        1, // group
                        1, // firstelem
                        nelements as LONGLONG,
                        &write_data,
                        NULL_VALUE,
                        &mut status,
                    );
                }
                assert_eq!(status, 0, "Failed to write SBYTE_IMG data with nulls");
            }

            // Close and reopen for reading
            unsafe {
                fits_close_file(fptr.take().unwrap(), &mut status);
            }
            assert_eq!(status, 0, "Failed to close file");

            fptr = None;
            unsafe {
                fits_open_file(
                    &mut fptr,
                    cast_slice(filename_cstr.to_bytes_with_nul()),
                    READWRITE,
                    &mut status,
                );
            }
            assert_eq!(status, 0, "Failed to open file for reading");

            if let Some(ref mut fptr_box) = fptr {
                let mut read_data: Vec<i8> = vec![0; nelements];
                let mut anynull: c_int = 0;

                unsafe {
                    fits_read_img(
                        fptr_box,
                        TSBYTE,
                        1,
                        nelements as LONGLONG,
                        Some(NullValue::Byte(SUBSTITUTE_VALUE)),
                        cast_slice_mut(&mut read_data),
                        Some(&mut anynull),
                        &mut status,
                    );
                }
                assert_eq!(status, 0, "Failed to read SBYTE_IMG data");
                assert_eq!(anynull, 1, "Should have detected null values");

                // Verify data
                for j in 0..IMAGE_HEIGHT {
                    for i in 0..IMAGE_WIDTH {
                        let index = (j * IMAGE_WIDTH + i) as usize;
                        if i == j {
                            assert_eq!(
                                read_data[index], SUBSTITUTE_VALUE,
                                "Null value not substituted at diagonal position ({}, {})",
                                i, j
                            );
                        } else {
                            assert_eq!(
                                write_data[index], read_data[index],
                                "Data mismatch at index {}",
                                index
                            );
                        }
                    }
                }
            }

            unsafe {
                fits_close_file(fptr.take().unwrap(), &mut status);
            }
            assert_eq!(status, 0, "Failed to close file after reading");
        });
    }

    #[test]
    fn test_ulong_img_null_read_write() {
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;
            const NULL_VALUE: c_ulong = 99;
            const SUBSTITUTE_VALUE: c_ulong = 101;

            // Create test data with nulls on diagonal
            let mut write_data: Vec<c_ulong> = vec![0; nelements];
            for j in 0..IMAGE_HEIGHT {
                for i in 0..IMAGE_WIDTH {
                    let index = (j * IMAGE_WIDTH + i) as usize;
                    if i == j {
                        write_data[index] = NULL_VALUE;
                    } else {
                        let val = ((i + j) * 1000) as c_ulong;
                        // Avoid NULL_VALUE in non-diagonal positions
                        write_data[index] = if val == NULL_VALUE {
                            NULL_VALUE + 1
                        } else {
                            val
                        };
                    }
                }
            }

            // Create FITS file and write image
            let filename_cstr = CString::new(filename).unwrap();
            unsafe {
                fits_create_file(
                    &mut fptr,
                    cast_slice(filename_cstr.to_bytes_with_nul()),
                    &mut status,
                );
            }
            assert_eq!(status, 0, "Failed to create file");

            if let Some(ref mut fptr_box) = fptr {
                unsafe {
                    fits_create_img(fptr_box, ULONG_IMG, naxis as c_int, &naxes, &mut status);
                }
                assert_eq!(status, 0, "Failed to create ULONG_IMG");

                // Set BLANK keyword to define null value for integer image
                unsafe {
                    fits_update_key_lng(
                        fptr_box,
                        cs!(c"BLANK"),
                        NULL_VALUE as c_long,
                        None,
                        &mut status,
                    );
                }
                assert_eq!(status, 0, "Failed to set BLANK keyword");

                unsafe {
                    fits_write_imgnull_ulng(
                        fptr_box,
                        1, // group
                        1, // firstelem
                        nelements as LONGLONG,
                        &write_data,
                        NULL_VALUE,
                        &mut status,
                    );
                }
                assert_eq!(status, 0, "Failed to write ULONG_IMG data with nulls");
            }

            // Close and reopen for reading
            unsafe {
                fits_close_file(fptr.take().unwrap(), &mut status);
            }
            assert_eq!(status, 0, "Failed to close file");

            fptr = None;
            unsafe {
                fits_open_file(
                    &mut fptr,
                    cast_slice(filename_cstr.to_bytes_with_nul()),
                    READWRITE,
                    &mut status,
                );
            }
            assert_eq!(status, 0, "Failed to open file for reading");

            if let Some(ref mut fptr_box) = fptr {
                let mut read_data: Vec<c_ulong> = vec![0; nelements];
                let mut anynull: c_int = 0;

                unsafe {
                    fits_read_img(
                        fptr_box,
                        TULONG,
                        1,
                        nelements as LONGLONG,
                        Some(NullValue::ULong(SUBSTITUTE_VALUE)),
                        cast_slice_mut(&mut read_data),
                        Some(&mut anynull),
                        &mut status,
                    );
                }
                assert_eq!(status, 0, "Failed to read ULONG_IMG data");
                assert_eq!(anynull, 1, "Should have detected null values");

                // Verify data
                for j in 0..IMAGE_HEIGHT {
                    for i in 0..IMAGE_WIDTH {
                        let index = (j * IMAGE_WIDTH + i) as usize;
                        if i == j {
                            assert_eq!(
                                read_data[index], SUBSTITUTE_VALUE,
                                "Null value not substituted at diagonal position ({}, {})",
                                i, j
                            );
                        } else {
                            assert_eq!(
                                write_data[index], read_data[index],
                                "Data mismatch at index {}",
                                index
                            );
                        }
                    }
                }
            }

            unsafe {
                fits_close_file(fptr.take().unwrap(), &mut status);
            }
            assert_eq!(status, 0, "Failed to close file after reading");
        });
    }

    #[test]
    fn test_ulonglong_img_null_read_write() {
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;
            const NULL_VALUE: c_ulonglong = 99;
            const SUBSTITUTE_VALUE: c_ulonglong = 101;

            // Create test data with nulls on diagonal
            let mut write_data: Vec<c_ulonglong> = vec![0; nelements];
            for j in 0..IMAGE_HEIGHT {
                for i in 0..IMAGE_WIDTH {
                    let index = (j * IMAGE_WIDTH + i) as usize;
                    if i == j {
                        write_data[index] = NULL_VALUE;
                    } else {
                        let val = ((i as c_ulonglong + j as c_ulonglong) * 100000);
                        // Avoid NULL_VALUE in non-diagonal positions
                        write_data[index] = if val == NULL_VALUE {
                            NULL_VALUE + 1
                        } else {
                            val
                        };
                    }
                }
            }

            // Create FITS file and write image
            let filename_cstr = CString::new(filename).unwrap();
            unsafe {
                fits_create_file(
                    &mut fptr,
                    cast_slice(filename_cstr.to_bytes_with_nul()),
                    &mut status,
                );
            }
            assert_eq!(status, 0, "Failed to create file");

            if let Some(ref mut fptr_box) = fptr {
                unsafe {
                    fits_create_img(fptr_box, ULONGLONG_IMG, naxis as c_int, &naxes, &mut status);
                }
                assert_eq!(status, 0, "Failed to create ULONGLONG_IMG");

                // Set BLANK keyword to define null value for integer image
                unsafe {
                    fits_update_key_lng(
                        fptr_box,
                        cs!(c"BLANK"),
                        NULL_VALUE as c_long,
                        None,
                        &mut status,
                    );
                }
                assert_eq!(status, 0, "Failed to set BLANK keyword");

                unsafe {
                    fits_write_imgnull_ulnglng(
                        fptr_box,
                        1, // group
                        1, // firstelem
                        nelements as LONGLONG,
                        &write_data,
                        NULL_VALUE,
                        &mut status,
                    );
                }
                assert_eq!(status, 0, "Failed to write ULONGLONG_IMG data with nulls");
            }

            // Close and reopen for reading
            unsafe {
                fits_close_file(fptr.take().unwrap(), &mut status);
            }
            assert_eq!(status, 0, "Failed to close file");

            fptr = None;
            unsafe {
                fits_open_file(
                    &mut fptr,
                    cast_slice(filename_cstr.to_bytes_with_nul()),
                    READWRITE,
                    &mut status,
                );
            }
            assert_eq!(status, 0, "Failed to open file for reading");

            if let Some(ref mut fptr_box) = fptr {
                let mut read_data: Vec<c_ulonglong> = vec![0; nelements];
                let mut anynull: c_int = 0;

                unsafe {
                    fits_read_img(
                        fptr_box,
                        TULONGLONG,
                        1,
                        nelements as LONGLONG,
                        Some(NullValue::ULONGLONG(SUBSTITUTE_VALUE)),
                        cast_slice_mut(&mut read_data),
                        Some(&mut anynull),
                        &mut status,
                    );
                }
                assert_eq!(status, 0, "Failed to read ULONGLONG_IMG data");
                assert_eq!(anynull, 1, "Should have detected null values");

                // Verify data
                for j in 0..IMAGE_HEIGHT {
                    for i in 0..IMAGE_WIDTH {
                        let index = (j * IMAGE_WIDTH + i) as usize;
                        if i == j {
                            assert_eq!(
                                read_data[index], SUBSTITUTE_VALUE,
                                "Null value not substituted at diagonal position ({}, {})",
                                i, j
                            );
                        } else {
                            assert_eq!(
                                write_data[index], read_data[index],
                                "Data mismatch at index {}",
                                index
                            );
                        }
                    }
                }
            }

            unsafe {
                fits_close_file(fptr.take().unwrap(), &mut status);
            }
            assert_eq!(status, 0, "Failed to close file after reading");
        });
    }
}
