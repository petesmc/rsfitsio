#![allow(deprecated)]

use std::ffi::CString;

use std::{process::ExitCode, ptr};

use bytemuck::cast_slice;

use rsfitsio::aliases::fits_is_compressed_image;
use rsfitsio::c_types::{FILE, c_char};
use rsfitsio::{STDERR, cs};

#[cfg(windows)]
use rsfitsio::__acrt_iob_func;

use rsfitsio::fitsio::{
    BYTE_IMG, DOUBLE_IMG, END_OF_FILE, FLOAT_IMG, IMAGE_HDU, LONG_IMG, LONGLONG, READONLY,
    SHORT_IMG, TBYTE, TDOUBLE, TFLOAT, TINT, TSHORT, TYP_CMPRS_KEY,
};
use rsfitsio::wrappers::{strchr_safe, strcpy_safe};
use rsfitsio::{aliases::fits_open_file, fitsio::fitsfile};

use rsfitsio::cfileio::ffclos as fits_close_file;
use rsfitsio::cfileio::ffinit as fits_create_file;
use rsfitsio::cfileio::ffrprt as fits_report_error;
use rsfitsio::edithdu::ffcopy as fits_copy_hdu;
use rsfitsio::fitscore::ffghdn as fits_get_hdu_num;
use rsfitsio::fitscore::ffghdt as fits_get_hdu_type;
use rsfitsio::fitscore::ffgipr as fits_get_img_param;
use rsfitsio::fitscore::ffgkcl as fits_get_keyclass;
use rsfitsio::fitscore::ffmrhd as fits_movrel_hdu;
use rsfitsio::getcol::ffgpv as fits_read_img;
use rsfitsio::getkey::ffgcrd as fits_read_card;
use rsfitsio::getkey::ffghsp as fits_get_hdrspace;
use rsfitsio::getkey::ffgrec as fits_read_record;
use rsfitsio::putcol::ffppr as fits_write_img;
use rsfitsio::putkey::ffcrim as fits_create_img;
use rsfitsio::putkey::ffprec as fits_write_record;
use rsfitsio::scalnull::ffpscl as fits_set_bscale;

pub fn main() -> ExitCode {
    /* FITS file pointers defined in fitsio.h */
    let mut infptr: Option<Box<fitsfile>> = None;
    let mut outfptr: Option<Box<fitsfile>> = None;

    let mut status = 0;
    let mut tstatus;

    let mut iteration;
    let mut single = 0;
    let mut hdupos = 0;
    let mut hdutype = 0;
    let mut bitpix = 0;
    let mut bytepix;
    let mut naxis = 0;
    let mut nkeys = 0;
    let mut datatype;
    let mut anynul = 0;
    let mut naxes = [1, 1, 1, 1, 1, 1, 1, 1, 1];
    let mut first;
    let mut totpix = 0;
    let mut npix;
    let mut array: Vec<f64>;
    let bscale = 1.0;
    let bzero = 0.0;
    let nulval = 0.0;
    let mut card: [c_char; 81] = [0; 81];

    let args = std::env::args();

    if args.len() != 3 {
        println!();
        println!("Usage:  imcopy inputImage outputImage[compress]");
        println!();
        println!("Copy an input image to an output image, optionally compressing");
        println!("or uncompressing the image in the process.  If the [compress]");
        println!("qualifier is appended to the output file name then the input image");
        println!("will be compressed using the tile-compressed format.  In this format,");
        println!("the image is divided into rectangular tiles and each tile of pixels");
        println!("is compressed and stored in a variable-length row of a binary table.");
        println!("If the [compress] qualifier is omitted, and the input image is");
        println!("in tile-compressed format, then the output image will be uncompressed.");
        println!();
        println!("If an extension name or number is appended to the input file name, ");
        println!("enclosed in square brackets, then only that single extension will be");
        println!("copied to the output file.  Otherwise, every extension in the input file");
        println!("will be processed in turn and copied to the output file.");
        println!();
        println!("Examples:");
        println!();
        println!("1)  imcopy image.fit 'cimage.fit[compress]'");
        println!();
        println!("    This compresses the input image using the default parameters, i.e.,");
        println!("    using the Rice compression algorithm and using row by row tiles.");
        println!();
        println!("2)  imcopy cimage.fit image2.fit");
        println!();
        println!("    This uncompresses the image created in the first example.");
        println!("    image2.fit should be identical to image.fit if the image");
        println!("    has an integer datatype.  There will be small differences");
        println!("    in the pixel values if it is a floating point image.");
        println!();
        println!("3)  imcopy image.fit 'cimage.fit[compress GZIP 100,100;q 16]'");
        println!();
        println!("    This compresses the input image using the following parameters:");
        println!("         GZIP compression algorithm;");
        println!("         100 X 100 pixel compression tiles;");
        println!("         quantization level = 16 (only used with floating point images)");
        println!();
        println!("The full syntax of the compression qualifier is:");
        println!("    [compress ALGORITHM TDIM1,TDIM2,...; q QLEVEL s SCALE]");
        println!("where the allowed ALGORITHM values are:");
        println!("      Rice, HCOMPRESS, HSCOMPRESS, GZIP, or PLIO. ");
        println!("       (HSCOMPRESS is a variant of HCOMPRESS in which a small");
        println!("        amount of smoothing is applied to the uncompressed image");
        println!("        to help suppress blocky compression artifacts in the image");
        println!("        when using large values for the 'scale' parameter).");
        println!("TDIMn is the size of the compression tile in each dimension,");
        println!();
        println!("QLEVEL specifies the quantization level when converting a floating");
        println!("point image into integers, prior to compressing the image.  The");
        println!("default value = 16, which means the image will be quantized into");
        println!("integer levels that are spaced at intervals of sigma/16., where ");
        println!("sigma is the estimated noise level in background areas of the image.");
        println!("If QLEVEL is negative, this means use the absolute value for the");
        println!("quantization spacing (e.g. 'q -0.005' means quantize the floating");
        println!("point image such that the scaled integers represent steps of 0.005");
        println!("in the original image).");
        println!();
        println!("SCALE is the integer scale factor that only applies to the HCOMPRESS");
        println!("algorithm.  The default value SCALE = 0 forces the image to be");
        println!("losslessly compressed; Greater amounts of lossy compression (resulting");
        println!("in smaller compressed files) can be specified with larger SCALE values.");
        println!();
        println!();
        println!("Note that it may be necessary to enclose the file names");
        println!("in single quote characters on the Unix command line.");
        return ExitCode::from(0);
    }

    let args: Vec<String> = args.collect();
    let infile = CString::new(args[1].clone()).unwrap();
    let outfile = CString::new(args[2].clone()).unwrap();

    unsafe {
        /* Open the input file and create output file */
        fits_open_file(&mut infptr, infile.as_ptr(), READONLY, &mut status);
        fits_create_file(&mut outfptr, outfile.as_ptr(), &mut status);

        if status != 0 {
            #[cfg(windows)]
            {
                use std::os::windows::io::AsRawHandle;
                let hd = std::io::stderr().as_raw_handle();
                let stream = unsafe { &mut libc::open_osfhandle(hd as isize, 0) };
                fits_report_error(stream as *mut _ as *mut FILE, status);
            }

            #[cfg(not(windows))]
            {
                use std::os::unix::io::AsRawFd;
                let stream = &mut std::io::stderr().as_raw_fd();
                fits_report_error(stream as *mut _ as *mut FILE, status);
            }

            return ExitCode::from(status as u8);
        }

        let mut infptr = infptr.unwrap();
        let mut outfptr = outfptr.unwrap();

        fits_get_hdu_num(infptr.as_mut(), &mut hdupos); /* Get the current HDU position */

        /* Copy only a single HDU if a specific extension was given */
        if hdupos != 1
            || strchr_safe(cast_slice(infile.to_bytes_with_nul()), b'[' as c_char).is_some()
        {
            single = 1;
        }

        while status == 0 {
            /* Main loop through each extension */

            fits_get_hdu_type(infptr.as_mut(), &mut hdutype, &mut status);

            if hdutype == IMAGE_HDU {
                /* get image dimensions and total number of pixels in image */
                naxes.fill(1);

                fits_get_img_param(
                    infptr.as_mut(),
                    9,
                    &mut bitpix,
                    &mut naxis,
                    naxes.as_mut_ptr(),
                    &mut status,
                );

                totpix = naxes[0]
                    * naxes[1]
                    * naxes[2]
                    * naxes[3]
                    * naxes[4]
                    * naxes[5]
                    * naxes[6]
                    * naxes[7]
                    * naxes[8];
            }

            if hdutype != IMAGE_HDU || naxis == 0 || totpix == 0 {
                /* just copy tables and null images */
                fits_copy_hdu(infptr.as_mut(), outfptr.as_mut(), 0, &mut status);
            } else {
                /* Explicitly create new image, to support compression */
                fits_create_img(outfptr.as_mut(), bitpix, naxis, naxes.as_ptr(), &mut status);
                if status != 0 {
                    fits_report_error(STDERR!() as *mut _ as *mut FILE, status);
                    return ExitCode::from(status as u8);
                }

                if fits_is_compressed_image(outfptr.as_mut(), &mut status) != 0 {
                    /* write default EXTNAME keyword if it doesn't already exist */
                    tstatus = 0;
                    fits_read_card(
                        infptr.as_mut(),
                        c"EXTNAME".as_ptr(),
                        card.as_mut_ptr(),
                        &mut tstatus,
                    );

                    if tstatus != 0 {
                        strcpy_safe(
                            &mut card,
                            cs!(
                                c"EXTNAME = 'COMPRESSED_IMAGE'   / name of this binary table extension"
                            ),
                        );
                        fits_write_record(outfptr.as_mut(), card.as_ptr(), &mut status);
                    }
                }

                /* copy all the user keywords (not the structural keywords) */
                fits_get_hdrspace(infptr.as_mut(), &mut nkeys, ptr::null_mut(), &mut status);

                for ii in 1..=nkeys {
                    fits_read_record(infptr.as_mut(), ii, card.as_mut_ptr(), &mut status);
                    if fits_get_keyclass(card.as_mut_ptr()) > TYP_CMPRS_KEY {
                        fits_write_record(outfptr.as_mut(), card.as_ptr(), &mut status);
                    }
                }

                match bitpix {
                    BYTE_IMG => {
                        datatype = TBYTE;
                    }
                    SHORT_IMG => {
                        datatype = TSHORT;
                    }
                    LONG_IMG => {
                        datatype = TINT;
                    }
                    FLOAT_IMG => {
                        datatype = TFLOAT;
                    }
                    DOUBLE_IMG => {
                        datatype = TDOUBLE;
                    }
                    _ => {
                        datatype = TBYTE; // Assume this is correct
                    }
                }

                bytepix = (bitpix.abs()) / 8;

                npix = totpix;
                iteration = 0;

                /* try to allocate memory for the entire image */
                /* use double type to force memory alignment */
                array = Vec::new();

                /* if allocation failed, divide size by 2 and try again */
                let mut additional = ((npix as f64 * bytepix as f64) / 8.0).ceil() as usize;
                while array.try_reserve_exact(additional).is_err() && iteration < 10 {
                    iteration += 1;
                    npix /= 2;
                    additional = ((npix as f64 * bytepix as f64) / 8.0).ceil() as usize;
                }

                if array.capacity() == 0 {
                    println!("Memory allocation error\n");
                    return ExitCode::from(0);
                } else {
                    array.resize(additional, 0.0);
                }

                /* turn off any scaling so that we copy the raw pixel values */
                fits_set_bscale(infptr.as_mut(), bscale, bzero, &mut status);
                fits_set_bscale(outfptr.as_mut(), bscale, bzero, &mut status);

                first = 1;
                while totpix > 0 && status == 0 {
                    /* read all or part of image then write it back to the output file */
                    fits_read_img(
                        infptr.as_mut(),
                        datatype,
                        first,
                        npix as LONGLONG,
                        &nulval as *const _ as *const _,
                        array.as_mut_ptr() as *mut _ as *mut _,
                        &mut anynul as *mut _ as *mut _,
                        &mut status,
                    );

                    fits_write_img(
                        outfptr.as_mut(),
                        datatype,
                        first,
                        npix as LONGLONG,
                        array.as_ptr() as *const _,
                        &mut status,
                    );
                    totpix -= npix;
                    first += npix as LONGLONG;
                }
            }

            if single != 0 {
                break;
            } /* quit if only copying a single HDU */
            fits_movrel_hdu(infptr.as_mut(), 1, ptr::null_mut(), &mut status); /* try to move to next HDU */

            hdupos += 1;
        }

        if status == END_OF_FILE {
            status = 0; /* Reset after normal error */
        }

        fits_close_file(Some(outfptr), &mut status);
        fits_close_file(Some(infptr), &mut status);

        /* if error occurred, print out error message */
        if status != 0 {
            fits_report_error(STDERR!() as *mut _ as *mut FILE, status);
        }
    }

    ExitCode::from(status as u8)
}
