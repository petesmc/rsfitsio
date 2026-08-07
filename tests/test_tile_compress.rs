//! Round-trip tests for tile-compressed images.

mod common;

#[cfg(test)]
mod tests {
    use crate::common::with_temp_file;
    use bytemuck::cast_slice;
    use libc::c_int;
    use rsfitsio::aliases::rust_api::*;
    use rsfitsio::fitsio::{
        DOUBLE_IMG, FLOAT_IMG, GZIP_2, LONGLONG, READONLY, SHORT_IMG, fitsfile,
    };
    use rsfitsio::imcompress::{fits_set_compression_type_safe, fits_set_quantize_level_safe};
    use std::ffi::CString;

    const NX: LONGLONG = 16;
    const NY: LONGLONG = 8;
    const NELEM: usize = (NX * NY) as usize;

    /// Regression test for https://github.com/cruzzil/rsfitsio/issues/89
    ///
    /// The GZIP_2 byte unshuffling used on read decremented its destination
    /// cursor once more than there were bytes to store, which underflowed the
    /// `usize` index and panicked in any build with overflow checks on.
    #[test]
    fn test_gzip2_float_round_trip() {
        with_temp_file(|filename| {
            let filename_cstr = CString::new(filename).unwrap();
            let name = cast_slice(filename_cstr.to_bytes_with_nul());
            let data: Vec<f32> = (0..NELEM).map(|i| i as f32 * 0.5).collect();
            let mut status: c_int = 0;

            let mut fptr: Option<Box<fitsfile>> = None;
            fits_create_diskfile(&mut fptr, name, &mut status);
            let f = fptr.as_mut().unwrap();
            fits_set_quantize_level_safe(f, 0.0, &mut status); // lossless floats
            fits_set_compression_type_safe(f, GZIP_2, &mut status);
            fits_create_imgll(f, FLOAT_IMG, 2, &[NX, NY], &mut status);
            fits_write_img_flt(f, 1, 1, NELEM as LONGLONG, &data, &mut status);
            fits_close_file(fptr.take().unwrap(), &mut status);
            assert_eq!(status, 0, "Failed to write GZIP_2 float image");

            fits_open_diskfile(&mut fptr, name, READONLY, &mut status);
            let f = fptr.as_mut().unwrap();
            fits_movabs_hdu(f, 2, None, &mut status); // the image is in the extension
            let mut out = vec![0f32; NELEM];
            fits_read_img_flt(f, 1, 1, NELEM as LONGLONG, 0.0, &mut out, None, &mut status);
            fits_close_file(fptr.take().unwrap(), &mut status);

            assert_eq!(status, 0, "Failed to read GZIP_2 float image");
            assert_eq!(out, data, "GZIP_2 float round trip differs");
        });
    }

    /// Regression test for https://github.com/cruzzil/rsfitsio/issues/89
    #[test]
    fn test_gzip2_double_round_trip() {
        with_temp_file(|filename| {
            let filename_cstr = CString::new(filename).unwrap();
            let name = cast_slice(filename_cstr.to_bytes_with_nul());
            let data: Vec<f64> = (0..NELEM).map(|i| i as f64 * 0.25).collect();
            let mut status: c_int = 0;

            let mut fptr: Option<Box<fitsfile>> = None;
            fits_create_diskfile(&mut fptr, name, &mut status);
            let f = fptr.as_mut().unwrap();
            fits_set_quantize_level_safe(f, 0.0, &mut status); // lossless doubles
            fits_set_compression_type_safe(f, GZIP_2, &mut status);
            fits_create_imgll(f, DOUBLE_IMG, 2, &[NX, NY], &mut status);
            fits_write_img_dbl(f, 1, 1, NELEM as LONGLONG, &data, &mut status);
            fits_close_file(fptr.take().unwrap(), &mut status);
            assert_eq!(status, 0, "Failed to write GZIP_2 double image");

            fits_open_diskfile(&mut fptr, name, READONLY, &mut status);
            let f = fptr.as_mut().unwrap();
            fits_movabs_hdu(f, 2, None, &mut status);
            let mut out = vec![0f64; NELEM];
            fits_read_img_dbl(f, 1, 1, NELEM as LONGLONG, 0.0, &mut out, None, &mut status);
            fits_close_file(fptr.take().unwrap(), &mut status);

            assert_eq!(status, 0, "Failed to read GZIP_2 double image");
            assert_eq!(out, data, "GZIP_2 double round trip differs");
        });
    }

    /// Regression test for https://github.com/cruzzil/rsfitsio/issues/89
    #[test]
    fn test_gzip2_short_round_trip() {
        with_temp_file(|filename| {
            let filename_cstr = CString::new(filename).unwrap();
            let name = cast_slice(filename_cstr.to_bytes_with_nul());
            let data: Vec<i16> = (0..NELEM).map(|i| (i as i16) * 3 - 100).collect();
            let mut status: c_int = 0;

            let mut fptr: Option<Box<fitsfile>> = None;
            fits_create_diskfile(&mut fptr, name, &mut status);
            let f = fptr.as_mut().unwrap();
            fits_set_compression_type_safe(f, GZIP_2, &mut status);
            fits_create_imgll(f, SHORT_IMG, 2, &[NX, NY], &mut status);
            fits_write_img_sht(f, 1, 1, NELEM as LONGLONG, &data, &mut status);
            fits_close_file(fptr.take().unwrap(), &mut status);
            assert_eq!(status, 0, "Failed to write GZIP_2 short image");

            fits_open_diskfile(&mut fptr, name, READONLY, &mut status);
            let f = fptr.as_mut().unwrap();
            fits_movabs_hdu(f, 2, None, &mut status);
            let mut out = vec![0i16; NELEM];
            fits_read_img_sht(f, 1, 1, NELEM as LONGLONG, 0, &mut out, None, &mut status);
            fits_close_file(fptr.take().unwrap(), &mut status);

            assert_eq!(status, 0, "Failed to read GZIP_2 short image");
            assert_eq!(out, data, "GZIP_2 short round trip differs");
        });
    }
}
