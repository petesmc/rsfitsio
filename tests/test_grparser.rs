// The extern "C" entry points are #[deprecated] so that internal callers reach
// for the _safe forms instead. These tests deliberately exercise the C ABI
// surface itself, which is what the attribute is warning about.
#![allow(deprecated)]

mod common;

#[cfg(test)]
mod tests {
    use crate::common::with_temp_file;
    use bytemuck::cast_slice;
    use rsfitsio::aliases::rust_api::*;
    use rsfitsio::c_types::{c_char, c_int, c_long};
    use rsfitsio::fitsio::fitsfile;
    use std::ffi::CString;
    use std::fs;
    use std::io::Write;

    /// Helper to create a temporary template file with given content
    fn create_temp_template<F>(content: &str, test_fn: F)
    where
        F: Fn(&str),
    {
        with_temp_file(|temp_path| {
            // Write template content
            {
                let mut file = fs::File::create(temp_path).unwrap();
                file.write_all(content.as_bytes()).unwrap();
            }

            // Run test
            test_fn(temp_path);
        });
    }

    #[test]
    fn test_fits_execute_template_simple_image() {
        with_temp_file(|filename| {
            create_temp_template(
                "SIMPLE  = T\nBITPIX  = 16\nNAXIS   = 2\nNAXIS1  = 10\nNAXIS2  = 10\nEND\n",
                |template_path| {
                    let mut fptr: Option<Box<fitsfile>> = None;
                    let mut status: c_int = 0;

                    // Create FITS file from template - this is the correct way!
                    let filename_cstr = CString::new(filename).unwrap();
                    let template_cstr = CString::new(template_path).unwrap();
                    fits_create_template(
                        &mut fptr,
                        cast_slice(filename_cstr.to_bytes_with_nul()),
                        cast_slice(template_cstr.to_bytes_with_nul()),
                        &mut status,
                    );
                    assert_eq!(status, 0, "Failed to create file from template");

                    // Verify the file was created correctly
                    let mut naxis: c_int = 0;
                    fits_get_img_dim(fptr.as_mut().unwrap().as_mut(), &mut naxis, &mut status);
                    assert_eq!(status, 0);
                    assert_eq!(naxis, 2);

                    let mut naxes: [c_long; 2] = [0; 2];
                    fits_get_img_size(fptr.as_mut().unwrap().as_mut(), 2, &mut naxes, &mut status);
                    assert_eq!(status, 0);
                    assert_eq!(naxes[0], 10);
                    assert_eq!(naxes[1], 10);

                    // Close file
                    fits_close_file(fptr.take().unwrap(), &mut status);
                    assert_eq!(status, 0);
                },
            );
        });
    }

    #[test]
    fn test_fits_execute_template_with_keywords() {
        with_temp_file(|filename| {
            create_temp_template(
                "SIMPLE  = T\n\
                 BITPIX  = 8\n\
                 NAXIS   = 0\n\
                 EXTEND  = T\n\
                 TELESCOP= 'TEST TELESCOPE'\n\
                 OBSERVER= 'Test Observer'\n\
                 DATE-OBS= '2024-01-01'\n\
                 END\n",
                |template_path| {
                    let mut fptr: Option<Box<fitsfile>> = None;
                    let mut status: c_int = 0;

                    // Create FITS file
                    let filename_cstr = CString::new(filename).unwrap();
                    fits_create_file(
                        &mut fptr,
                        cast_slice(filename_cstr.to_bytes_with_nul()),
                        &mut status,
                    );
                    assert_eq!(status, 0);

                    // Execute template
                    let template_cstr = CString::new(template_path).unwrap();
                    unsafe {
                        rsfitsio::grparser::fits_execute_template(
                            fptr.as_mut().unwrap().as_mut(),
                            template_cstr.as_ptr(),
                            &raw mut status,
                        );
                    }
                    assert_eq!(status, 0, "Failed to execute template");

                    // Verify keywords were written
                    let mut telescop: Vec<c_char> = vec![0; 80];
                    let telescop_key = CString::new("TELESCOP").unwrap();
                    fits_read_key_str(
                        fptr.as_mut().unwrap().as_mut(),
                        cast_slice(telescop_key.to_bytes_with_nul()),
                        &mut telescop,
                        None,
                        &mut status,
                    );
                    assert_eq!(status, 0);

                    // Close file
                    fits_close_file(fptr.take().unwrap(), &mut status);
                    assert_eq!(status, 0);
                },
            );
        });
    }

    #[test]
    fn test_fits_execute_template_binary_table() {
        with_temp_file(|filename| {
            create_temp_template(
                "XTENSION= BINTABLE\n\
                 BITPIX  = 8\n\
                 NAXIS   = 2\n\
                 NAXIS1  = 16\n\
                 NAXIS2  = 10\n\
                 PCOUNT  = 0\n\
                 GCOUNT  = 1\n\
                 TFIELDS = 2\n\
                 TTYPE1  = 'Column1'\n\
                 TFORM1  = '1J'\n\
                 TTYPE2  = 'Column2'\n\
                 TFORM2  = '1E'\n\
                 EXTNAME = 'TEST_TABLE'\n\
                 END\n",
                |template_path| {
                    let mut fptr: Option<Box<fitsfile>> = None;
                    let mut status: c_int = 0;

                    // Create FITS file with primary HDU first
                    let filename_cstr = CString::new(filename).unwrap();
                    fits_create_file(
                        &mut fptr,
                        cast_slice(filename_cstr.to_bytes_with_nul()),
                        &mut status,
                    );
                    assert_eq!(status, 0);

                    // Create a minimal primary HDU
                    let naxes: [c_long; 0] = [];
                    fits_create_img(fptr.as_mut().unwrap().as_mut(), 8, 0, &naxes, &mut status);
                    assert_eq!(status, 0);

                    // Execute template for binary table extension
                    let template_cstr = CString::new(template_path).unwrap();
                    unsafe {
                        rsfitsio::grparser::fits_execute_template(
                            fptr.as_mut().unwrap().as_mut(),
                            template_cstr.as_ptr(),
                            &raw mut status,
                        );
                    }
                    assert_eq!(status, 0, "Failed to execute template");

                    // Verify table was created
                    let mut ncols: c_int = 0;
                    fits_get_num_cols(fptr.as_mut().unwrap().as_mut(), &mut ncols, &mut status);
                    assert_eq!(status, 0);
                    assert_eq!(ncols, 2);

                    // Close file
                    fits_close_file(fptr.take().unwrap(), &mut status);
                    assert_eq!(status, 0);
                },
            );
        });
    }

    #[test]
    fn test_fits_execute_template_with_comments() {
        with_temp_file(|filename| {
            create_temp_template(
                "# This is a comment line\n\
                 SIMPLE  = T\n\
                 # Another comment\n\
                 BITPIX  = 16\n\
                 NAXIS   = 0\n\
                 # Final comment\n\
                 END\n",
                |template_path| {
                    let mut fptr: Option<Box<fitsfile>> = None;
                    let mut status: c_int = 0;

                    // Create FITS file
                    let filename_cstr = CString::new(filename).unwrap();
                    fits_create_file(
                        &mut fptr,
                        cast_slice(filename_cstr.to_bytes_with_nul()),
                        &mut status,
                    );
                    assert_eq!(status, 0);

                    // Execute template
                    let template_cstr = CString::new(template_path).unwrap();
                    unsafe {
                        rsfitsio::grparser::fits_execute_template(
                            fptr.as_mut().unwrap().as_mut(),
                            template_cstr.as_ptr(),
                            &raw mut status,
                        );
                    }
                    assert_eq!(
                        status, 0,
                        "Template with comments should execute successfully"
                    );

                    // Close file
                    fits_close_file(fptr.take().unwrap(), &mut status);
                    assert_eq!(status, 0);
                },
            );
        });
    }

    #[test]
    fn test_fits_execute_template_history_and_comment_keywords() {
        with_temp_file(|filename| {
            create_temp_template(
                "SIMPLE  = T\n\
                 BITPIX  = 8\n\
                 NAXIS   = 0\n\
                 HISTORY This file was created by test\n\
                 HISTORY Second history line\n\
                 COMMENT This is a comment\n\
                 COMMENT Another comment line\n\
                 END\n",
                |template_path| {
                    let mut fptr: Option<Box<fitsfile>> = None;
                    let mut status: c_int = 0;

                    // Create FITS file
                    let filename_cstr = CString::new(filename).unwrap();
                    fits_create_file(
                        &mut fptr,
                        cast_slice(filename_cstr.to_bytes_with_nul()),
                        &mut status,
                    );
                    assert_eq!(status, 0);

                    // Execute template
                    let template_cstr = CString::new(template_path).unwrap();
                    unsafe {
                        rsfitsio::grparser::fits_execute_template(
                            fptr.as_mut().unwrap().as_mut(),
                            template_cstr.as_ptr(),
                            &raw mut status,
                        );
                    }
                    assert_eq!(
                        status, 0,
                        "Template with HISTORY/COMMENT should execute successfully"
                    );

                    // Close file
                    fits_close_file(fptr.take().unwrap(), &mut status);
                    assert_eq!(status, 0);
                },
            );
        });
    }

    #[test]
    fn test_fits_execute_template_nonexistent_file() {
        with_temp_file(|filename| {
            let mut fptr: Option<Box<fitsfile>> = None;
            let mut status: c_int = 0;

            // Create FITS file
            let filename_cstr = CString::new(filename).unwrap();
            fits_create_file(
                &mut fptr,
                cast_slice(filename_cstr.to_bytes_with_nul()),
                &mut status,
            );
            assert_eq!(status, 0);

            // Try to execute non-existent template
            let template_cstr = CString::new("/nonexistent/template.tpl").unwrap();
            unsafe {
                rsfitsio::grparser::fits_execute_template(
                    fptr.as_mut().unwrap().as_mut(),
                    template_cstr.as_ptr(),
                    &raw mut status,
                );
            }
            assert_ne!(status, 0, "Should fail with non-existent template file");

            // Close file (ignore status since we expect an error)
            status = 0;
            fits_close_file(fptr.take().unwrap(), &mut status);
        });
    }
}
