#[cfg(test)]
mod tests {
    use bytemuck::{cast_slice, cast_slice_mut};
    use libc::{
        c_double, c_float, c_int, c_long, c_longlong, c_short, c_uchar, c_ulong, c_ulonglong,
    };
    use rsfitsio::aliases::rust_api::*;
    use rsfitsio::fitsio::{
        BYTE_IMG, DOUBLE_IMG, FLOAT_IMG, LONG_IMG, LONGLONG, LONGLONG_IMG, READWRITE, SBYTE_IMG,
        SHORT_IMG, TBYTE, TDOUBLE, TFLOAT, TLONG, TLONGLONG, TSBYTE, TSHORT, TULONG, TULONGLONG,
        TUSHORT, ULONG_IMG, ULONGLONG_IMG, USHORT_IMG, fitsfile,
    };
    use rsfitsio::helpers::testhelpers::{floats_close_f32, floats_close_f64, with_temp_file};
    use std::ffi::CString;

    // Test dimensions for all images
    const IMAGE_WIDTH: c_long = 50;
    const IMAGE_HEIGHT: c_long = 50;

    #[test]
    fn test_byte_img_read_write() {
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;

            // Create test data - unsigned byte values
            let mut write_data: Vec<c_uchar> = vec![0; nelements];
            for j in 0..IMAGE_HEIGHT {
                for i in 0..IMAGE_WIDTH {
                    let index = (j * IMAGE_WIDTH + i) as usize;
                    write_data[index] = ((i + j) % 256) as c_uchar;
                }
            }

            // Create FITS file and write image
            let filename_cstr = CString::new(filename).unwrap();
            fits_create_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                &mut status,
            );

            assert_eq!(status, 0, "Failed to create file");

            if let Some(ref mut fptr_box) = fptr {
                // Create image HDU
                fits_create_img(fptr_box, BYTE_IMG, naxis as c_int, &naxes, &mut status);

                assert_eq!(status, 0, "Failed to create BYTE_IMG");

                // Write image data
                fits_write_img(
                    fptr_box,
                    TBYTE,
                    1,
                    nelements as LONGLONG,
                    cast_slice(&write_data),
                    &mut status,
                );

                assert_eq!(status, 0, "Failed to write BYTE_IMG data");
            }

            // Close and reopen for reading
            fits_close_file(fptr.take().unwrap(), &mut status);
            assert_eq!(status, 0, "Failed to close file");

            // Read back the data
            fptr = None;
            fits_open_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                READWRITE,
                &mut status,
            );
            assert_eq!(status, 0, "Failed to open file for reading");

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

                assert_eq!(status, 0, "Failed to read BYTE_IMG data");
                assert_eq!(anynull, 0, "Unexpected null values");

                // Verify data
                for i in 0..nelements {
                    assert_eq!(write_data[i], read_data[i], "Data mismatch at index {}", i);
                }
            }

            // Close file
            fits_close_file(fptr.take().unwrap(), &mut status);
            assert_eq!(status, 0, "Failed to close file after reading");
        });
    }

    #[test]
    fn test_short_img_read_write() {
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;

            // Create test data - signed short values
            let mut write_data: Vec<c_short> = vec![0; nelements];
            for j in 0..IMAGE_HEIGHT {
                for i in 0..IMAGE_WIDTH {
                    let index = (j * IMAGE_WIDTH + i) as usize;
                    write_data[index] = ((i - 50) * (j - 50)) as c_short;
                }
            }

            // Create FITS file and write image
            let filename_cstr = CString::new(filename).unwrap();
            fits_create_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                &mut status,
            );

            assert_eq!(status, 0, "Failed to create file");

            if let Some(ref mut fptr_box) = fptr {
                // Create image HDU
                fits_create_img(fptr_box, SHORT_IMG, naxis as c_int, &naxes, &mut status);

                assert_eq!(status, 0, "Failed to create SHORT_IMG");

                // Write image data

                fits_write_img(
                    fptr_box,
                    TSHORT,
                    1,
                    nelements as LONGLONG,
                    cast_slice(&write_data),
                    &mut status,
                );

                assert_eq!(status, 0, "Failed to write SHORT_IMG data");
            }

            // Close and reopen for reading
            fits_close_file(fptr.take().unwrap(), &mut status);
            assert_eq!(status, 0, "Failed to close file");

            // Read back the data
            fptr = None;
            fits_open_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                READWRITE,
                &mut status,
            );
            assert_eq!(status, 0, "Failed to open file for reading");

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

                assert_eq!(status, 0, "Failed to read SHORT_IMG data");
                assert_eq!(anynull, 0, "Unexpected null values");

                // Verify data
                for i in 0..nelements {
                    assert_eq!(write_data[i], read_data[i], "Data mismatch at index {}", i);
                }
            }

            // Close file
            fits_close_file(fptr.take().unwrap(), &mut status);
            assert_eq!(status, 0, "Failed to close file after reading");
        });
    }

    #[test]
    fn test_ushort_img_read_write() {
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;

            // Create test data - unsigned short values
            let mut write_data: Vec<u16> = vec![0; nelements];
            for j in 0..IMAGE_HEIGHT {
                for i in 0..IMAGE_WIDTH {
                    let index = (j * IMAGE_WIDTH + i) as usize;
                    write_data[index] = ((i + j) * 200) as u16;
                }
            }

            // Create FITS file and write image
            let filename_cstr = CString::new(filename).unwrap();
            fits_create_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                &mut status,
            );

            assert_eq!(status, 0, "Failed to create file");

            if let Some(ref mut fptr_box) = fptr {
                // Create image HDU - this will store as BITPIX=16 with BZERO=32768
                fits_create_img(fptr_box, USHORT_IMG, naxis as c_int, &naxes, &mut status);
                assert_eq!(status, 0, "Failed to create USHORT_IMG");

                // Write image data
                fits_write_img(
                    fptr_box,
                    TUSHORT,
                    1,
                    nelements as LONGLONG,
                    cast_slice(&write_data),
                    &mut status,
                );
                assert_eq!(status, 0, "Failed to write USHORT_IMG data");
            }

            // Close and reopen for reading
            fits_close_file(fptr.take().unwrap(), &mut status);
            assert_eq!(status, 0, "Failed to close file");

            // Read back the data
            fptr = None;
            fits_open_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                READWRITE,
                &mut status,
            );
            assert_eq!(status, 0, "Failed to open file for reading");

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

                assert_eq!(status, 0, "Failed to read USHORT_IMG data");
                assert_eq!(anynull, 0, "Unexpected null values");

                // Verify data
                for i in 0..nelements {
                    assert_eq!(write_data[i], read_data[i], "Data mismatch at index {}", i);
                }
            }

            // Close file
            fits_close_file(fptr.take().unwrap(), &mut status);
            assert_eq!(status, 0, "Failed to close file after reading");
        });
    }

    #[test]
    fn test_long_img_read_write() {
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;

            // Create test data - signed long values
            let mut write_data: Vec<c_long> = vec![0; nelements];
            for j in 0..IMAGE_HEIGHT {
                for i in 0..IMAGE_WIDTH {
                    let index = (j * IMAGE_WIDTH + i) as usize;
                    write_data[index] = ((i - 25) * 1000 + (j - 25) * 10) as c_long;
                }
            }

            // Create FITS file and write image
            let filename_cstr = CString::new(filename).unwrap();
            fits_create_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                &mut status,
            );

            assert_eq!(status, 0, "Failed to create file");

            if let Some(ref mut fptr_box) = fptr {
                // Create image HDU
                fits_create_img(fptr_box, LONG_IMG, naxis as c_int, &naxes, &mut status);
                assert_eq!(status, 0, "Failed to create LONG_IMG");

                // Write image data

                fits_write_img(
                    fptr_box,
                    TLONG,
                    1,
                    nelements as LONGLONG,
                    cast_slice(&write_data),
                    &mut status,
                );

                assert_eq!(status, 0, "Failed to write LONG_IMG data");
            }

            // Close and reopen for reading
            fits_close_file(fptr.take().unwrap(), &mut status);
            assert_eq!(status, 0, "Failed to close file");

            // Read back the data
            fptr = None;
            fits_open_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                READWRITE,
                &mut status,
            );
            assert_eq!(status, 0, "Failed to open file for reading");

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

                assert_eq!(status, 0, "Failed to read LONG_IMG data");
                assert_eq!(anynull, 0, "Unexpected null values");

                // Verify data
                for i in 0..nelements {
                    assert_eq!(write_data[i], read_data[i], "Data mismatch at index {}", i);
                }
            }

            // Close file
            fits_close_file(fptr.take().unwrap(), &mut status);
            assert_eq!(status, 0, "Failed to close file after reading");
        });
    }

    #[test]
    fn test_longlong_img_read_write() {
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;

            // Create test data - signed long long values
            let mut write_data: Vec<c_longlong> = vec![0; nelements];
            for j in 0..IMAGE_HEIGHT {
                for i in 0..IMAGE_WIDTH {
                    let index = (j * IMAGE_WIDTH + i) as usize;
                    write_data[index] =
                        (i - 25) as c_longlong * 100000 + (j - 25) as c_longlong * 1000;
                }
            }

            // Create FITS file and write image
            let filename_cstr = CString::new(filename).unwrap();
            fits_create_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                &mut status,
            );

            assert_eq!(status, 0, "Failed to create file");

            if let Some(ref mut fptr_box) = fptr {
                // Create image HDU
                fits_create_img(fptr_box, LONGLONG_IMG, naxis as c_int, &naxes, &mut status);

                assert_eq!(status, 0, "Failed to create LONGLONG_IMG");

                // Write image data

                fits_write_img(
                    fptr_box,
                    TLONGLONG,
                    1,
                    nelements as LONGLONG,
                    cast_slice(&write_data),
                    &mut status,
                );

                assert_eq!(status, 0, "Failed to write LONGLONG_IMG data");
            }

            // Close and reopen for reading
            fits_close_file(fptr.take().unwrap(), &mut status);
            assert_eq!(status, 0, "Failed to close file");

            // Read back the data
            fptr = None;
            fits_open_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                READWRITE,
                &mut status,
            );
            assert_eq!(status, 0, "Failed to open file for reading");

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

                assert_eq!(status, 0, "Failed to read LONGLONG_IMG data");
                assert_eq!(anynull, 0, "Unexpected null values");

                // Verify data
                for i in 0..nelements {
                    assert_eq!(write_data[i], read_data[i], "Data mismatch at index {}", i);
                }
            }

            // Close file
            fits_close_file(fptr.take().unwrap(), &mut status);
            assert_eq!(status, 0, "Failed to close file after reading");
        });
    }

    #[test]
    fn test_float_img_read_write() {
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;

            // Create test data - float values
            let mut write_data: Vec<c_float> = vec![0.0; nelements];
            for j in 0..IMAGE_HEIGHT {
                for i in 0..IMAGE_WIDTH {
                    let index = (j * IMAGE_WIDTH + i) as usize;
                    write_data[index] = (i - 25) as f32 * 10.0 + (j - 25) as f32 * 1.0;
                }
            }

            // Create FITS file and write image
            let filename_cstr = CString::new(filename).unwrap();
            fits_create_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                &mut status,
            );

            assert_eq!(status, 0, "Failed to create file");

            if let Some(ref mut fptr_box) = fptr {
                // Create image HDU
                fits_create_img(fptr_box, FLOAT_IMG, naxis as c_int, &naxes, &mut status);

                assert_eq!(status, 0, "Failed to create FLOAT_IMG");

                // Write image data

                fits_write_img(
                    fptr_box,
                    TFLOAT,
                    1,
                    nelements as LONGLONG,
                    cast_slice(&write_data),
                    &mut status,
                );

                assert_eq!(status, 0, "Failed to write FLOAT_IMG data");
            }

            // Close and reopen for reading
            fits_close_file(fptr.take().unwrap(), &mut status);
            assert_eq!(status, 0, "Failed to close file");

            // Read back the data
            fptr = None;
            fits_open_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                READWRITE,
                &mut status,
            );
            assert_eq!(status, 0, "Failed to open file for reading");

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

                assert_eq!(status, 0, "Failed to read FLOAT_IMG data");
                assert_eq!(anynull, 0, "Unexpected null values");

                // Verify data with floating point tolerance
                for i in 0..nelements {
                    assert!(
                        floats_close_f32(write_data[i], read_data[i]),
                        "Data mismatch at index {}: {} != {}",
                        i,
                        write_data[i],
                        read_data[i]
                    );
                }
            }

            // Close file
            fits_close_file(fptr.take().unwrap(), &mut status);
            assert_eq!(status, 0, "Failed to close file after reading");
        });
    }

    #[test]
    fn test_double_img_read_write() {
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;

            // Create test data - double values
            let mut write_data: Vec<c_double> = vec![0.0; nelements];
            for j in 0..IMAGE_HEIGHT {
                for i in 0..IMAGE_WIDTH {
                    let index = (j * IMAGE_WIDTH + i) as usize;
                    write_data[index] = (i - 25) as f64 * 100.0 + (j - 25) as f64 * 10.0;
                }
            }

            // Create FITS file and write image
            let filename_cstr = CString::new(filename).unwrap();
            fits_create_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                &mut status,
            );

            assert_eq!(status, 0, "Failed to create file");

            if let Some(ref mut fptr_box) = fptr {
                // Create image HDU
                fits_create_img(fptr_box, DOUBLE_IMG, naxis as c_int, &naxes, &mut status);

                assert_eq!(status, 0, "Failed to create DOUBLE_IMG");

                // Write image data

                fits_write_img(
                    fptr_box,
                    TDOUBLE,
                    1,
                    nelements as LONGLONG,
                    cast_slice(&write_data),
                    &mut status,
                );

                assert_eq!(status, 0, "Failed to write DOUBLE_IMG data");
            }

            // Close and reopen for reading
            fits_close_file(fptr.take().unwrap(), &mut status);
            assert_eq!(status, 0, "Failed to close file");

            // Read back the data
            fptr = None;
            fits_open_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                READWRITE,
                &mut status,
            );
            assert_eq!(status, 0, "Failed to open file for reading");

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

                assert_eq!(status, 0, "Failed to read DOUBLE_IMG data");
                assert_eq!(anynull, 0, "Unexpected null values");

                // Verify data with floating point tolerance
                for i in 0..nelements {
                    assert!(
                        floats_close_f64(write_data[i], read_data[i]),
                        "Data mismatch at index {}: {} != {}",
                        i,
                        write_data[i],
                        read_data[i]
                    );
                }
            }

            // Close file
            fits_close_file(fptr.take().unwrap(), &mut status);
            assert_eq!(status, 0, "Failed to close file after reading");
        });
    }

    #[test]
    fn test_sbyte_img_read_write() {
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;

            // Create test data - signed byte values
            let mut write_data: Vec<i8> = vec![0; nelements];
            for j in 0..IMAGE_HEIGHT {
                for i in 0..IMAGE_WIDTH {
                    let index = (j * IMAGE_WIDTH + i) as usize;
                    write_data[index] = ((i + j) % 255 - 127) as i8;
                }
            }

            // Create FITS file and write image
            let filename_cstr = CString::new(filename).unwrap();
            fits_create_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                &mut status,
            );

            assert_eq!(status, 0, "Failed to create file");

            if let Some(ref mut fptr_box) = fptr {
                // Create image HDU - this will store as BITPIX=8 with BZERO=-128
                fits_create_img(fptr_box, SBYTE_IMG, naxis as c_int, &naxes, &mut status);

                assert_eq!(status, 0, "Failed to create SBYTE_IMG");

                // Write image data

                fits_write_img(
                    fptr_box,
                    TSBYTE,
                    1,
                    nelements as LONGLONG,
                    cast_slice(&write_data),
                    &mut status,
                );

                assert_eq!(status, 0, "Failed to write SBYTE_IMG data");
            }

            // Close and reopen for reading
            fits_close_file(fptr.take().unwrap(), &mut status);
            assert_eq!(status, 0, "Failed to close file");

            // Read back the data
            fptr = None;
            fits_open_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                READWRITE,
                &mut status,
            );
            assert_eq!(status, 0, "Failed to open file for reading");

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

                assert_eq!(status, 0, "Failed to read SBYTE_IMG data");
                assert_eq!(anynull, 0, "Unexpected null values");

                // Verify data
                for i in 0..nelements {
                    assert_eq!(write_data[i], read_data[i], "Data mismatch at index {}", i);
                }
            }

            // Close file
            fits_close_file(fptr.take().unwrap(), &mut status);
            assert_eq!(status, 0, "Failed to close file after reading");
        });
    }

    #[test]
    fn test_ulong_img_read_write() {
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;

            // Create test data - unsigned long values
            let mut write_data: Vec<c_ulong> = vec![0; nelements];
            for j in 0..IMAGE_HEIGHT {
                for i in 0..IMAGE_WIDTH {
                    let index = (j * IMAGE_WIDTH + i) as usize;
                    write_data[index] = ((i + j) * 1000) as c_ulong;
                }
            }

            // Create FITS file and write image
            let filename_cstr = CString::new(filename).unwrap();
            fits_create_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                &mut status,
            );

            assert_eq!(status, 0, "Failed to create file");

            if let Some(ref mut fptr_box) = fptr {
                // Create image HDU - this will store as BITPIX=32 with BZERO=2147483648
                fits_create_img(fptr_box, ULONG_IMG, naxis as c_int, &naxes, &mut status);

                assert_eq!(status, 0, "Failed to create ULONG_IMG");

                // Write image data

                fits_write_img(
                    fptr_box,
                    TULONG,
                    1,
                    nelements as LONGLONG,
                    cast_slice(&write_data),
                    &mut status,
                );

                assert_eq!(status, 0, "Failed to write ULONG_IMG data");
            }

            // Close and reopen for reading
            fits_close_file(fptr.take().unwrap(), &mut status);
            assert_eq!(status, 0, "Failed to close file");

            // Read back the data
            fptr = None;
            fits_open_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                READWRITE,
                &mut status,
            );
            assert_eq!(status, 0, "Failed to open file for reading");

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

                assert_eq!(status, 0, "Failed to read ULONG_IMG data");
                assert_eq!(anynull, 0, "Unexpected null values");

                // Verify data
                for i in 0..nelements {
                    assert_eq!(write_data[i], read_data[i], "Data mismatch at index {}", i);
                }
            }

            // Close file
            fits_close_file(fptr.take().unwrap(), &mut status);
            assert_eq!(status, 0, "Failed to close file after reading");
        });
    }

    #[test]
    fn test_ulonglong_img_read_write() {
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;
            let naxis: c_long = 2;
            let naxes: [c_long; 2] = [IMAGE_WIDTH, IMAGE_HEIGHT];
            let nelements = (IMAGE_WIDTH * IMAGE_HEIGHT) as usize;

            // Create test data - unsigned long long values
            let mut write_data: Vec<c_ulonglong> = vec![0; nelements];
            for j in 0..IMAGE_HEIGHT {
                for i in 0..IMAGE_WIDTH {
                    let index = (j * IMAGE_WIDTH + i) as usize;
                    write_data[index] = (i as c_ulonglong + j as c_ulonglong) * 100000;
                }
            }

            // Create FITS file and write image
            let filename_cstr = CString::new(filename).unwrap();
            fits_create_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                &mut status,
            );

            assert_eq!(status, 0, "Failed to create file");

            if let Some(ref mut fptr_box) = fptr {
                // Create image HDU - this will store as BITPIX=64 with BZERO=9223372036854775808
                fits_create_img(fptr_box, ULONGLONG_IMG, naxis as c_int, &naxes, &mut status);

                assert_eq!(status, 0, "Failed to create ULONGLONG_IMG");

                // Write image data

                fits_write_img(
                    fptr_box,
                    TULONGLONG,
                    1,
                    nelements as LONGLONG,
                    cast_slice(&write_data),
                    &mut status,
                );

                assert_eq!(status, 0, "Failed to write ULONGLONG_IMG data");
            }

            // Close and reopen for reading
            fits_close_file(fptr.take().unwrap(), &mut status);
            assert_eq!(status, 0, "Failed to close file");

            // Read back the data
            fptr = None;
            fits_open_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                READWRITE,
                &mut status,
            );
            assert_eq!(status, 0, "Failed to open file for reading");

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

                assert_eq!(status, 0, "Failed to read ULONGLONG_IMG data");
                assert_eq!(anynull, 0, "Unexpected null values");

                // Verify data
                for i in 0..nelements {
                    assert_eq!(write_data[i], read_data[i], "Data mismatch at index {}", i);
                }
            }

            // Close file
            fits_close_file(fptr.take().unwrap(), &mut status);
            assert_eq!(status, 0, "Failed to close file after reading");
        });
    }
}
