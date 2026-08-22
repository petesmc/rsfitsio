//! Quantization of floating point images for lossy compression.
//!
//! The tile compressors work on integers, so a `float` or `double` image is
//! first quantized: each tile's values are divided by a scale factor and
//! rounded, and the scale and zero point are recorded in the `ZSCALE` and
//! `ZZERO` columns so the values can be reconstructed. The scale is chosen from
//! an estimate of the tile's noise, so that the quantization error stays below
//! the noise already present -- which is what makes the loss acceptable.
//!
//! Rounding alone leaves the reconstructed values visibly striped at low signal
//! levels, so a dither may be added before rounding and subtracted after; see
//! `DitherType`, which is crate-private. `SubtractiveDither2` additionally
//! preserves exact zeros, which would otherwise be perturbed.
//!
//! Based on algorithms written by Richard White at STScI, made available for
//! use in CFITSIO in July 1999 and updated in January 2008.
#![warn(missing_docs)]

use core::slice;

use crate::c_types::*;

use bytemuck::{cast, cast_slice_mut};

use crate::{
    fitsio::{LONGLONG, MEMORY_ALLOCATION},
    fitsio2::N_RANDOM,
    imcompress::{FITS_RAND_VALUE, fits_init_randoms},
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
/// Which dithering scheme to apply when quantizing a floating point tile.
pub(crate) enum DitherType {
    /// Round without dithering.
    NoDither = -1,
    /// Add a pseudo-random dither before rounding and subtract it after.
    SubtractiveDither1 = 1,
    /// As [`DitherType::SubtractiveDither1`], but a value of exactly zero is
    /// preserved rather than dithered.
    SubtractiveDither2 = 2,
}

impl From<c_int> for DitherType {
    fn from(v: c_int) -> Self {
        match v {
            -1 => DitherType::NoDither,
            1 => DitherType::SubtractiveDither1,
            2 => DitherType::SubtractiveDither2,
            _ => DitherType::NoDither,
        }
    }
}

/* nearest integer function */
fn nint_f64(x: f64) -> i32 {
    if x >= 0.0 {
        (x + 0.5) as i32
    } else {
        (x - 0.5) as i32
    }
}

const NULL_VALUE: i32 = -2147483647; /* value used to represent undefined pixels */
const ZERO_VALUE: i32 = -2147483646; /* value used to represent zero-valued pixels */

/* number of reserved values, starting with
and including NULL_VALUE.  These values
may not be used to represent the quantized
and scaled floating point pixel values
If lossy Hcompression is used, and the
array contains null values, then it is also
possible for the compressed values to slightly
exceed the range of the actual (lossless) values
so we must reserve a little more space */
const N_RESERVED_VALUES: i32 = 10;

/* more than this many standard deviations from the mean is an outlier */
/// Values further than this many sigma from the mean are clipped when
/// estimating a tile's noise.
const SIGMA_CLIP: f64 = 5.0;

/// Number of sigma-clipping iterations used when estimating a tile's noise.
const NITER: i32 = 3;

/// Writes through an optional output reference, doing nothing when the caller
/// passed `None`.
///
/// The C writes to an output parameter only after testing the pointer against
/// NULL; this is that idiom, as a method.
pub trait SomeSet<T> {
    /// Store `v` if there is somewhere to store it.
    fn set_if_some(self, v: T);
}

impl<T> SomeSet<T> for Option<&mut T> {
    fn set_if_some(self, v: T) {
        let s = self;
        if let Some(s) = s {
            *s = v;
        }
    }
}

/// The function value will be one if the input fdata were copied to idata;
/// in this case the parameters bscale and bzero can be used to convert back to
/// nearly the original floating point values:  fdata ~= idata * bscale + bzero.
/// If the function value is zero, the data were not copied to idata.
/// Assumes that f32 and c_int are the same size (4 bytes).
pub(crate) fn fits_quantize_float_inplace(
    row: usize, // tile number = row number in the binary table (this is only used when dithering the quantized values)
    fdata: &mut [f32], // array of image pixels to be compressed. Modified in place after applying bzero and bscale
    nxpix: usize,      // number of pixels in each row of fdata
    nypix: usize,      // number of rows in fdata
    nullcheck: bool,   // check for nullvalues in fdata?
    in_null_value: f32, // value used to represent undefined pixels in fdata
    qlevel: f32,       // quantization level
    dither_method: DitherType, // which dithering method to use
    bscale: &mut f64,  // scale factor
    bzero: &mut f64,   // zero offset
    iminval: &mut i32, // minimum quantized value that is returned
    imaxval: &mut i32, // maximum quantized value that is returned
) -> c_int {
    let mut iseed: usize = 0;

    let mut ngood: usize = 0;
    let mut stdev: f64;
    let mut noise2: f64 = 0.0;
    let mut noise3: f64 = 0.0;
    let mut noise5: f64 = 0.0; /* MAD 2nd, 3rd, and 5th order noise values */
    let mut minval: f32 = 0.0;
    let mut maxval: f32 = 0.0; /* min & max of fdata */
    let delta: f64; /* bscale, 1 in idata = delta in fdata */
    let mut zeropt: f64; /* bzero */

    let mut nextrand: usize = 0;
    let iqfactor: i64;

    let mut status = 0;

    let nx: usize = nxpix * nypix;

    if nx <= 1 {
        *bscale = 1.;
        *bzero = 0.;
        return 0;
    }

    if FITS_RAND_VALUE.get().is_none() && fits_init_randoms() != 0 {
        return MEMORY_ALLOCATION;
    }

    let fits_rand_value = FITS_RAND_VALUE.get().unwrap();

    if qlevel >= 0.0 {
        /* estimate background noise using MAD pixel differences */
        FnNoise5_float(
            fdata,
            nxpix,
            nypix,
            nullcheck,
            in_null_value,
            Some(&mut ngood),
            Some(&mut minval),
            Some(&mut maxval),
            Some(&mut noise2),
            Some(&mut noise3),
            Some(&mut noise5),
            &mut status,
        );

        if nullcheck && ngood == 0 {
            /* special case of an image filled with Nulls */
            /* set parameters to dummy values, which are not used */
            minval = 0.0;
            maxval = 1.0;
            stdev = 1.0;
        } else {
            /* use the minimum of noise2, noise3, and noise5 as the best noise value */
            stdev = noise3;
            if noise2 != 0.0 && noise2 < stdev {
                stdev = noise2;
            }
            if noise5 != 0.0 && noise5 < stdev {
                stdev = noise5;
            }
        }

        if qlevel == 0.0 {
            delta = stdev / 4.0; /* default quantization */
        } else {
            delta = stdev / f64::from(qlevel);
        }
        if delta == 0.0 {
            return 0; /* don't quantize */
        }
    } else {
        /* negative value represents the absolute quantization level */
        delta = f64::from(-qlevel);

        /* only nned to calculate the min and max values */
        FnNoise3_float(
            fdata,
            nxpix,
            nypix,
            nullcheck,
            in_null_value,
            Some(&mut ngood),
            Some(&mut minval),
            Some(&mut maxval),
            None,
            &mut status,
        );
    }

    /* check that the range of quantized levels is not > range of int */
    if f64::from(maxval - minval) / delta > 2.0 * 2147483647.0 - f64::from(N_RESERVED_VALUES) {
        return 0; /* don't quantize */
    }

    if row > 0 {
        /* we need to dither the quantized values */
        /* initialize the index to the next random number in the list */
        iseed = (row - 1) % N_RANDOM;
        nextrand = (fits_rand_value[iseed] * 500.0) as usize;
    }

    if ngood == nx {
        /* don't have to check for nulls */
        /* return all positive values, if possible since some */
        /* compression algorithms either only work for positive integers, */
        /* or are more efficient.  */

        if dither_method == DitherType::SubtractiveDither2 {
            /* shift the range to be close to the value used to represent zeros */
            zeropt = f64::from(minval) - delta * f64::from(NULL_VALUE + N_RESERVED_VALUES);
        } else if f64::from(maxval - minval) / delta < 2147483647.0 - f64::from(N_RESERVED_VALUES) {
            zeropt = minval.into();
            /* fudge the zero point so it is an integer multiple of delta */
            /* This helps to ensure the same scaling will be performed if the */
            /* file undergoes multiple fpack/funpack cycles */
            iqfactor = (zeropt / delta + 0.5) as i64;
            zeropt = iqfactor as f64 * delta;
        } else {
            /* center the quantized levels around zero */
            zeropt = f64::from(minval + maxval) / 2.;
        }

        if row > 0 {
            /* dither the values when quantizing */
            for i in 0..nx {
                // for (i = 0;  i < nx;  i+=1) {

                if dither_method == DitherType::SubtractiveDither2 && fdata[i] == 0.0 {
                    fdata[i] = cast(ZERO_VALUE);
                } else {
                    fdata[i] = cast(nint_f64(
                        ((f64::from(fdata[i]) - zeropt) / delta)
                            + f64::from(fits_rand_value[nextrand])
                            - 0.5,
                    ));
                }

                nextrand += 1;
                if nextrand == N_RANDOM {
                    iseed += 1;
                    if iseed == N_RANDOM {
                        iseed = 0;
                    }
                    nextrand = (fits_rand_value[iseed] * 500.0) as usize;
                }
            }
        } else {
            /* do not dither the values */

            for i in 0..nx {
                // for (i = 0;  i < nx;  i+=1) {
                fdata[i] = cast(nint_f64((f64::from(fdata[i]) - zeropt) / delta));
            }
        }
    } else {
        /* data contains null values; shift the range to be */
        /* close to the value used to represent null values */
        zeropt = f64::from(minval) - delta * f64::from(NULL_VALUE + N_RESERVED_VALUES);

        if row > 0 {
            /* dither the values */
            for i in 0..nx {
                // for (i = 0;  i < nx;  i+=1) {
                if fdata[i] != in_null_value {
                    if dither_method == DitherType::SubtractiveDither2 && fdata[i] == 0.0 {
                        fdata[i] = cast(ZERO_VALUE);
                    } else {
                        fdata[i] = cast(nint_f64(
                            ((f64::from(fdata[i]) - zeropt) / delta)
                                + f64::from(fits_rand_value[nextrand])
                                - 0.5,
                        ));
                    }
                } else {
                    fdata[i] = cast(NULL_VALUE);
                }

                /* increment the random number index, regardless */
                nextrand += 1;
                if nextrand == N_RANDOM {
                    iseed += 1;
                    if iseed == N_RANDOM {
                        iseed = 0;
                    }
                    nextrand = (fits_rand_value[iseed] * 500.0) as usize;
                }
            }
        } else {
            /* do not dither the values */
            for i in 0..nx {
                // for (i = 0;  i < nx;  i+=1) {
                if fdata[i] != in_null_value {
                    fdata[i] = cast(nint_f64((f64::from(fdata[i]) - zeropt) / delta));
                } else {
                    fdata[i] = cast(NULL_VALUE);
                }
            }
        }
    }

    /* calc min and max values */
    let mut temp: f64 = (f64::from(minval) - zeropt) / delta;
    *iminval = nint_f64(temp);
    temp = (f64::from(maxval) - zeropt) / delta;
    *imaxval = nint_f64(temp);

    *bscale = delta;
    *bzero = zeropt;

    /* yes, data have been quantized */
    1
}

/// The function value will be one if the input fdata were copied to idata;
/// in this case the parameters bscale and bzero can be used to convert back to
/// nearly the original floating point values:  fdata ~= idata * bscale + bzero.
/// If the function value is zero, the data were not copied to idata.
pub(crate) fn fits_quantize_double_inplace(
    row: usize, // tile number = row number in the binary table (this is only used when dithering the quantized values)
    fdata: &mut [f64], // array of image pixels to be compressed. Modified in place after applying bzero and bscale
    nxpix: usize,      // number of pixels in each row of fdata
    nypix: usize,      // number of rows in fdata
    nullcheck: bool,   // check for nullvalues in fdata?
    in_null_value: f64, // value used to represent undefined pixels in fdata
    qlevel: f32,       // quantization level
    dither_method: DitherType, // which dithering method to use
    bscale: &mut f64,  // scale factor
    bzero: &mut f64,   // zero offset
    iminval: &mut i32, // minimum quantized value that is returned
    imaxval: &mut i32, // maximum quantized value that is returned
) -> c_int {
    let mut iseed: usize = 0;

    let mut ngood: usize = 0;
    let mut stdev: f64;
    let mut noise2: f64 = 0.0;
    let mut noise3: f64 = 0.0;
    let mut noise5: f64 = 0.0; /* MAD 2nd, 3rd, and 5th order noise values */
    let mut minval: f64 = 0.0;
    let mut maxval: f64 = 0.0; /* min & max of fdata */
    let delta: f64; /* bscale, 1 in idata = delta in fdata */
    let mut zeropt: f64; /* bzero */

    let mut nextrand: usize = 0;
    let iqfactor: i64;

    let mut status = 0;

    let nx: usize = nxpix * nypix;
    if nx <= 1 {
        *bscale = 1.;
        *bzero = 0.;
        return 0;
    }

    if FITS_RAND_VALUE.get().is_none() && fits_init_randoms() != 0 {
        return MEMORY_ALLOCATION;
    }

    let fits_rand_value = FITS_RAND_VALUE.get().unwrap();

    if qlevel >= 0.0 {
        /* estimate background noise using MAD pixel differences */
        FnNoise5_double(
            fdata,
            nxpix,
            nypix,
            nullcheck,
            in_null_value,
            Some(&mut ngood),
            Some(&mut minval),
            Some(&mut maxval),
            Some(&mut noise2),
            Some(&mut noise3),
            Some(&mut noise5),
            &mut status,
        );

        if nullcheck && ngood == 0 {
            /* special case of an image filled with Nulls */
            /* set parameters to dummy values, which are not used */
            minval = 0.0;
            maxval = 1.0;
            stdev = 1.0;
        } else {
            /* use the minimum of noise2, noise3, and noise5 as the best noise value */
            stdev = noise3;
            if noise2 != 0.0 && noise2 < stdev {
                stdev = noise2;
            }
            if noise5 != 0.0 && noise5 < stdev {
                stdev = noise5;
            }
        }

        if qlevel == 0.0 {
            delta = stdev / 4.0; /* default quantization */
        } else {
            delta = stdev / f64::from(qlevel);
        }
        if delta == 0.0 {
            return 0; /* don't quantize */
        }
    } else {
        /* negative value represents the absolute quantization level */
        delta = f64::from(-qlevel);

        /* only nned to calculate the min and max values */
        FnNoise3_double(
            fdata,
            nxpix,
            nypix,
            nullcheck,
            in_null_value,
            Some(&mut ngood),
            Some(&mut minval),
            Some(&mut maxval),
            None,
            &mut status,
        );
    }

    /* check that the range of quantized levels is not > range of int */
    if (maxval - minval) / delta > 2.0 * 2147483647.0 - f64::from(N_RESERVED_VALUES) {
        return 0; /* don't quantize */
    }

    if row > 0 {
        /* we need to dither the quantized values */
        /* initialize the index to the next random number in the list */
        iseed = (row - 1) % N_RANDOM;
        nextrand = (fits_rand_value[iseed] * 500.0) as usize;
    }

    if ngood == nx {
        /* don't have to check for nulls */
        /* return all positive values, if possible since some */
        /* compression algorithms either only work for positive integers, */
        /* or are more efficient.  */

        if dither_method == DitherType::SubtractiveDither2 {
            /* shift the range to be close to the value used to represent zeros */
            zeropt = minval - delta * f64::from(NULL_VALUE + N_RESERVED_VALUES);
        } else if (maxval - minval) / delta < 2147483647.0 - f64::from(N_RESERVED_VALUES) {
            zeropt = minval;
            /* fudge the zero point so it is an integer multiple of delta */
            /* This helps to ensure the same scaling will be performed if the */
            /* file undergoes multiple fpack/funpack cycles */
            iqfactor = (zeropt / delta + 0.5) as i64;
            zeropt = iqfactor as f64 * delta;
        } else {
            /* center the quantized levels around zero */
            zeropt = (minval + maxval) / 2.;
        }

        if row > 0 {
            /* dither the values when quantizing */
            for i in 0..nx {
                let fdata_i = fdata[i];
                let idata = cast_slice_mut(fdata);

                if dither_method == DitherType::SubtractiveDither2 && fdata_i == 0.0 {
                    idata[i] = ZERO_VALUE;
                } else {
                    idata[i] = nint_f64(
                        ((fdata_i - zeropt) / delta) + f64::from(fits_rand_value[nextrand]) - 0.5,
                    );
                }

                nextrand += 1;
                if nextrand == N_RANDOM {
                    iseed += 1;
                    if iseed == N_RANDOM {
                        iseed = 0;
                    }
                    nextrand = (fits_rand_value[iseed] * 500.0) as usize;
                }
            }
        } else {
            /* do not dither the values */

            for i in 0..nx {
                let fdata_i = fdata[i];
                let idata = cast_slice_mut(fdata);
                idata[i] = nint_f64((fdata_i - zeropt) / delta);
            }
        }
    } else {
        /* data contains null values; shift the range to be */
        /* close to the value used to represent null values */
        zeropt = minval - delta * f64::from(NULL_VALUE + N_RESERVED_VALUES);

        if row > 0 {
            /* dither the values */
            for i in 0..nx {
                let fdata_i = fdata[i];
                let idata = cast_slice_mut(fdata);
                if fdata_i != in_null_value {
                    if dither_method == DitherType::SubtractiveDither2 && fdata_i == 0.0 {
                        idata[i] = ZERO_VALUE;
                    } else {
                        idata[i] = nint_f64(
                            ((fdata_i - zeropt) / delta) + f64::from(fits_rand_value[nextrand])
                                - 0.5,
                        );
                    }
                } else {
                    idata[i] = NULL_VALUE;
                }

                /* increment the random number index, regardless */
                nextrand += 1;
                if nextrand == N_RANDOM {
                    iseed += 1;
                    if iseed == N_RANDOM {
                        iseed = 0;
                    }
                    nextrand = (fits_rand_value[iseed] * 500.0) as usize;
                }
            }
        } else {
            /* do not dither the values */
            for i in 0..nx {
                let fdata_i = fdata[i];
                let idata = cast_slice_mut(fdata);
                if fdata_i != in_null_value {
                    idata[i] = nint_f64((fdata_i - zeropt) / delta);
                } else {
                    idata[i] = NULL_VALUE;
                }
            }
        }
    }

    /* calc min and max values */
    let mut temp: f64 = (minval - zeropt) / delta;
    *iminval = nint_f64(temp);
    temp = (maxval - zeropt) / delta;
    *imaxval = nint_f64(temp);

    *bscale = delta;
    *bzero = zeropt;

    /* yes, data have been quantized */
    1
}

/// Compute statistics of the input short integer image.
///
/// # Parameters
///
/// * `array`     — 2 dimensional array of image pixels
/// * `nx`        — number of pixels in each row of the image
/// * `ny`        — number of rows in the image
/// * `nullcheck` — check for null values, if true
/// * `nullvalue` — value of null pixels, if nullcheck is true
/// * `ngoodpix`  — number of non-null pixels in the image
/// * `minvalue`  — returned minimum non-null value in the array
/// * `maxvalue`  — returned maximum non-null value in the array
/// * `mean`      — returned mean value of all non-null pixels
/// * `sigma`     — returned R.M.S. value of all non-null pixels
/// * `noise1`    — 1st order estimate of noise in image background level
/// * `noise2`    — 2nd order estimate of noise in image background level
/// * `noise3`    — 3rd order estimate of noise in image background level
/// * `noise5`    — 5th order estimate of noise in image background level
/// * `status`    — error status
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn fits_img_stats_short(
    array: *const c_short,
    nx: c_long,
    ny: c_long,
    /* (if this is a 3D image, then ny should be the */
    /* product of the no. of rows times the no. of planes) */
    nullcheck: c_int,
    nullvalue: c_short,

    /* returned parameters (if the pointer is not null)  */
    ngoodpix: *mut c_long,
    minvalue: *mut c_short,
    maxvalue: *mut c_short,
    mean: *mut f64,
    sigma: *mut f64,
    noise1: *mut f64,
    noise2: *mut f64,
    noise3: *mut f64,
    noise5: *mut f64,
    status: *mut c_int,
) -> c_int {
    // FFI WRAPPER
    unsafe {
        let ngoodpix = ngoodpix.as_mut();
        let minvalue = minvalue.as_mut();
        let maxvalue = maxvalue.as_mut();
        let mean = mean.as_mut();
        let sigma = sigma.as_mut();
        let noise1 = noise1.as_mut();
        let noise2 = noise2.as_mut();
        let noise3 = noise3.as_mut();
        let noise5 = noise5.as_mut();
        let status = status.as_mut().unwrap();

        let array = slice::from_raw_parts(array, (nx * ny) as usize);

        let nullcheck = nullcheck != 0;

        fits_img_stats_short_safe(
            array, nx, ny, nullcheck, nullvalue, ngoodpix, minvalue, maxvalue, mean, sigma, noise1,
            noise2, noise3, noise5, status,
        )
    }
}

/// Compute statistics of the input short integer image.
///
/// # Parameters
///
/// * `array`     — 2 dimensional array of image pixels
/// * `nx`        — number of pixels in each row of the image
/// * `ny`        — number of rows in the image
/// * `nullcheck` — check for null values, if true
/// * `nullvalue` — value of null pixels, if nullcheck is true
/// * `ngoodpix`  — number of non-null pixels in the image
/// * `minvalue`  — returned minimum non-null value in the array
/// * `maxvalue`  — returned maximum non-null value in the array
/// * `mean`      — returned mean value of all non-null pixels
/// * `sigma`     — returned R.M.S. value of all non-null pixels
/// * `noise1`    — 1st order estimate of noise in image background level
/// * `noise2`    — 2nd order estimate of noise in image background level
/// * `noise3`    — 3rd order estimate of noise in image background level
/// * `noise5`    — 5th order estimate of noise in image background level
/// * `status`    — error status
pub fn fits_img_stats_short_safe(
    array: &[c_short],
    nx: c_long,
    ny: c_long,
    /* (if this is a 3D image, then ny should be the */
    /* product of the no. of rows times the no. of planes) */
    nullcheck: bool,
    nullvalue: c_short,

    /* returned parameters (if the pointer is not null)  */
    mut ngoodpix: Option<&mut c_long>,
    minvalue: Option<&mut c_short>,
    maxvalue: Option<&mut c_short>,
    mean: Option<&mut f64>,
    sigma: Option<&mut f64>,
    noise1: Option<&mut f64>,
    noise2: Option<&mut f64>,
    noise3: Option<&mut f64>,
    noise5: Option<&mut f64>,
    status: &mut c_int,
) -> c_int {
    let mut ngood: usize = 0;
    let mut minval: c_short = 0;
    let mut maxval: c_short = 0;
    let mut xmean: f64 = 0.;
    let mut xsigma: f64 = 0.;
    let mut xnoise: f64 = 0.;
    let mut xnoise2: f64 = 0.;
    let mut xnoise3: f64 = 0.;
    let mut xnoise5: f64 = 0.;

    /* need to calculate mean and/or sigma and/or limits? */
    if mean.is_some() || sigma.is_some() {
        FnMeanSigma_short(
            array,
            (nx * ny) as usize,
            nullcheck,
            nullvalue,
            &mut ngood,
            &mut xmean,
            &mut xsigma,
            status,
        );

        if let Some(ngoodpix) = ngoodpix.as_deref_mut() {
            *ngoodpix = ngood as c_long;
        }

        if let Some(mean) = mean {
            *mean = xmean;
        }

        if let Some(sigma) = sigma {
            *sigma = xsigma;
        }
    }

    if let Some(noise1) = noise1 {
        FnNoise1_short(
            array,
            nx as usize,
            ny as usize,
            nullcheck,
            nullvalue,
            &mut xnoise,
            status,
        );

        *noise1 = xnoise;
    }

    if minvalue.is_some() || maxvalue.is_some() || noise3.is_some() {
        FnNoise5_short(
            array,
            nx as usize,
            ny as usize,
            nullcheck,
            nullvalue,
            Some(&mut ngood),
            Some(&mut minval),
            Some(&mut maxval),
            Some(&mut xnoise2),
            Some(&mut xnoise3),
            Some(&mut xnoise5),
            status,
        );

        if let Some(ngoodpix) = ngoodpix {
            *ngoodpix = ngood as c_long;
        }
        if let Some(minvalue) = minvalue {
            *minvalue = minval;
        }
        if let Some(maxvalue) = maxvalue {
            *maxvalue = maxval;
        }
        if let Some(noise2) = noise2 {
            *noise2 = xnoise2;
        }
        if let Some(noise3) = noise3 {
            *noise3 = xnoise3;
        }
        if let Some(noise5) = noise5 {
            *noise5 = xnoise5;
        }
    }
    *status
}

/// Compute statistics of the input integer image.
///
/// # Parameters
///
/// * `array`     — 2 dimensional array of image pixels
/// * `nx`        — number of pixels in each row of the image
/// * `ny`        — number of rows in the image
/// * `nullcheck` — check for null values, if true
/// * `nullvalue` — value of null pixels, if nullcheck is true
/// * `ngoodpix`  — number of non-null pixels in the image
/// * `minvalue`  — returned minimum non-null value in the array
/// * `maxvalue`  — returned maximum non-null value in the array
/// * `mean`      — returned mean value of all non-null pixels
/// * `sigma`     — returned R.M.S. value of all non-null pixels
/// * `noise1`    — 1st order estimate of noise in image background level
/// * `noise2`    — 2nd order estimate of noise in image background level
/// * `noise3`    — 3rd order estimate of noise in image background level
/// * `noise5`    — 5th order estimate of noise in image background level
/// * `status`    — error status
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn fits_img_stats_int(
    array: *const c_int,
    nx: c_long,
    ny: c_long,
    /* (if this is a 3D image, then ny should be the */
    /* product of the no. of rows times the no. of planes) */
    nullcheck: c_int,
    nullvalue: c_int,

    /* returned parameters (if the pointer is not null)  */
    ngoodpix: *mut c_long,
    minvalue: *mut c_int,
    maxvalue: *mut c_int,
    mean: *mut f64,
    sigma: *mut f64,
    noise1: *mut f64,
    noise2: *mut f64,
    noise3: *mut f64,
    noise5: *mut f64,
    status: *mut c_int,
) -> c_int {
    // FFI WRAPPER
    unsafe {
        let ngoodpix = ngoodpix.as_mut();
        let minvalue = minvalue.as_mut();
        let maxvalue = maxvalue.as_mut();
        let mean = mean.as_mut();
        let sigma = sigma.as_mut();
        let noise1 = noise1.as_mut();
        let noise2 = noise2.as_mut();
        let noise3 = noise3.as_mut();
        let noise5 = noise5.as_mut();
        let status = status.as_mut().unwrap();

        let array = slice::from_raw_parts(array, (nx * ny) as usize);

        let nullcheck = nullcheck != 0;

        fits_img_stats_int_safe(
            array, nx, ny, nullcheck, nullvalue, ngoodpix, minvalue, maxvalue, mean, sigma, noise1,
            noise2, noise3, noise5, status,
        )
    }
}

/// Compute statistics of the input integer image.
///
/// # Parameters
///
/// * `array`     — 2 dimensional array of image pixels
/// * `nx`        — number of pixels in each row of the image
/// * `ny`        — number of rows in the image
/// * `nullcheck` — check for null values, if true
/// * `nullvalue` — value of null pixels, if nullcheck is true
/// * `ngoodpix`  — number of non-null pixels in the image
/// * `minvalue`  — returned minimum non-null value in the array
/// * `maxvalue`  — returned maximum non-null value in the array
/// * `mean`      — returned mean value of all non-null pixels
/// * `sigma`     — returned R.M.S. value of all non-null pixels
/// * `noise1`    — 1st order estimate of noise in image background level
/// * `noise2`    — 2nd order estimate of noise in image background level
/// * `noise3`    — 3rd order estimate of noise in image background level
/// * `noise5`    — 5th order estimate of noise in image background level
/// * `status`    — error status
pub fn fits_img_stats_int_safe(
    array: &[c_int],
    nx: c_long,
    ny: c_long,
    /* (if this is a 3D image, then ny should be the */
    /* product of the no. of rows times the no. of planes) */
    nullcheck: bool,
    nullvalue: c_int,

    /* returned parameters (if the pointer is not null)  */
    mut ngoodpix: Option<&mut c_long>,
    minvalue: Option<&mut c_int>,
    maxvalue: Option<&mut c_int>,
    mean: Option<&mut f64>,
    sigma: Option<&mut f64>,
    noise1: Option<&mut f64>,
    noise2: Option<&mut f64>,
    noise3: Option<&mut f64>,
    noise5: Option<&mut f64>,
    status: &mut c_int,
) -> c_int {
    let mut ngood: usize = 0;
    let mut minval: c_int = 0;
    let mut maxval: c_int = 0;
    let mut xmean: f64 = 0.;
    let mut xsigma: f64 = 0.;
    let mut xnoise: f64 = 0.;
    let mut xnoise2: f64 = 0.;
    let mut xnoise3: f64 = 0.;
    let mut xnoise5: f64 = 0.;

    /* need to calculate mean and/or sigma and/or limits? */
    if mean.is_some() || sigma.is_some() {
        FnMeanSigma_int(
            array,
            (nx * ny) as usize,
            nullcheck,
            nullvalue,
            &mut ngood,
            &mut xmean,
            &mut xsigma,
            status,
        );

        if let Some(ngoodpix) = ngoodpix.as_deref_mut() {
            *ngoodpix = ngood as c_long;
        }
        if let Some(mean) = mean {
            *mean = xmean;
        }
        if let Some(sigma) = sigma {
            *sigma = xsigma;
        }
    }

    if let Some(noise1) = noise1 {
        FnNoise1_int(
            array,
            nx as usize,
            ny as usize,
            nullcheck,
            nullvalue,
            &mut xnoise,
            status,
        );

        *noise1 = xnoise;
    }

    if minvalue.is_some() || maxvalue.is_some() || noise3.is_some() {
        FnNoise5_int(
            array,
            nx as usize,
            ny as usize,
            nullcheck,
            nullvalue,
            Some(&mut ngood),
            Some(&mut minval),
            Some(&mut maxval),
            Some(&mut xnoise2),
            Some(&mut xnoise3),
            Some(&mut xnoise5),
            status,
        );

        if let Some(ngoodpix) = ngoodpix {
            *ngoodpix = ngood as c_long
        };
        if let Some(minvalue) = minvalue {
            *minvalue = minval
        };
        if let Some(maxvalue) = maxvalue {
            *maxvalue = maxval
        };
        if let Some(noise2) = noise2 {
            *noise2 = xnoise2
        };
        if let Some(noise3) = noise3 {
            *noise3 = xnoise3
        };
        if let Some(noise5) = noise5 {
            *noise5 = xnoise5
        };
    }
    *status
}

/// Compute statistics of the input float image.
///
/// # Parameters
///
/// * `array`     — 2 dimensional array of image pixels
/// * `nx`        — number of pixels in each row of the image
/// * `ny`        — number of rows in the image
/// * `nullcheck` — check for null values, if true
/// * `nullvalue` — value of null pixels, if nullcheck is true
/// * `ngoodpix`  — number of non-null pixels in the image
/// * `minvalue`  — returned minimum non-null value in the array
/// * `maxvalue`  — returned maximum non-null value in the array
/// * `mean`      — returned mean value of all non-null pixels
/// * `sigma`     — returned R.M.S. value of all non-null pixels
/// * `noise1`    — 1st order estimate of noise in image background level
/// * `noise2`    — 2nd order estimate of noise in image background level
/// * `noise3`    — 3rd order estimate of noise in image background level
/// * `noise5`    — 5th order estimate of noise in image background level
/// * `status`    — error status
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn fits_img_stats_float(
    array: *const f32,
    nx: c_long,
    ny: c_long,
    /* (if this is a 3D image, then ny should be the */
    /* product of the no. of rows times the no. of planes) */
    nullcheck: c_int,
    nullvalue: f32,

    /* returned parameters (if the pointer is not null)  */
    ngoodpix: *mut c_long,
    minvalue: *mut f32,
    maxvalue: *mut f32,
    mean: *mut f64,
    sigma: *mut f64,
    noise1: *mut f64,
    noise2: *mut f64,
    noise3: *mut f64,
    noise5: *mut f64,
    status: *mut c_int,
) -> c_int {
    // FFI WRAPPER
    unsafe {
        let ngoodpix = ngoodpix.as_mut();
        let minvalue = minvalue.as_mut();
        let maxvalue = maxvalue.as_mut();
        let mean = mean.as_mut();
        let sigma = sigma.as_mut();
        let noise1 = noise1.as_mut();
        let noise2 = noise2.as_mut();
        let noise3 = noise3.as_mut();
        let noise5 = noise5.as_mut();
        let status = status.as_mut().unwrap();

        let array = slice::from_raw_parts(array, (nx * ny) as usize);

        let nullcheck = nullcheck != 0;

        fits_img_stats_float_safe(
            array, nx, ny, nullcheck, nullvalue, ngoodpix, minvalue, maxvalue, mean, sigma, noise1,
            noise2, noise3, noise5, status,
        )
    }
}

/// Compute statistics of the input float image.
///
/// # Parameters
///
/// * `array`     — 2 dimensional array of image pixels
/// * `nx`        — number of pixels in each row of the image
/// * `ny`        — number of rows in the image
/// * `nullcheck` — check for null values, if true
/// * `nullvalue` — value of null pixels, if nullcheck is true
/// * `ngoodpix`  — number of non-null pixels in the image
/// * `minvalue`  — returned minimum non-null value in the array
/// * `maxvalue`  — returned maximum non-null value in the array
/// * `mean`      — returned mean value of all non-null pixels
/// * `sigma`     — returned R.M.S. value of all non-null pixels
/// * `noise1`    — 1st order estimate of noise in image background level
/// * `noise2`    — 2nd order estimate of noise in image background level
/// * `noise3`    — 3rd order estimate of noise in image background level
/// * `noise5`    — 5th order estimate of noise in image background level
/// * `status`    — error status
pub fn fits_img_stats_float_safe(
    array: &[f32],
    nx: c_long,
    ny: c_long,
    /* (if this is a 3D image, then ny should be the */
    /* product of the no. of rows times the no. of planes) */
    nullcheck: bool,
    nullvalue: f32,

    /* returned parameters (if the pointer is not null)  */
    mut ngoodpix: Option<&mut c_long>,
    minvalue: Option<&mut f32>,
    maxvalue: Option<&mut f32>,
    mean: Option<&mut f64>,
    sigma: Option<&mut f64>,
    noise1: Option<&mut f64>,
    noise2: Option<&mut f64>,
    noise3: Option<&mut f64>,
    noise5: Option<&mut f64>,
    status: &mut c_int,
) -> c_int {
    let mut ngood: usize = 0;
    let mut minval: f32 = 0.0;
    let mut maxval: f32 = 0.0;
    let mut xmean: f64 = 0.;
    let mut xsigma: f64 = 0.;
    let mut xnoise: f64 = 0.;
    let mut xnoise2: f64 = 0.;
    let mut xnoise3: f64 = 0.;
    let mut xnoise5: f64 = 0.;

    /* need to calculate mean and/or sigma and/or limits? */
    if mean.is_some() || sigma.is_some() {
        FnMeanSigma_float(
            array,
            (nx * ny) as usize,
            nullcheck,
            nullvalue,
            &mut ngood,
            &mut xmean,
            &mut xsigma,
            status,
        );

        if let Some(ngoodpix) = ngoodpix.as_deref_mut() {
            *ngoodpix = ngood as c_long;
        };
        if let Some(mean) = mean {
            *mean = xmean;
        };
        if let Some(sigma) = sigma {
            *sigma = xsigma;
        };
    }

    if let Some(noise1) = noise1 {
        FnNoise1_float(
            array,
            nx as usize,
            ny as usize,
            nullcheck,
            nullvalue,
            &mut xnoise,
            status,
        );

        *noise1 = xnoise;
    }

    if minvalue.is_some() || maxvalue.is_some() || noise3.is_some() {
        FnNoise5_float(
            array,
            nx as usize,
            ny as usize,
            nullcheck,
            nullvalue,
            Some(&mut ngood),
            Some(&mut minval),
            Some(&mut maxval),
            Some(&mut xnoise2),
            Some(&mut xnoise3),
            Some(&mut xnoise5),
            status,
        );

        if let Some(ngoodpix) = ngoodpix {
            *ngoodpix = ngood as c_long
        };
        if let Some(minvalue) = minvalue {
            *minvalue = minval
        };
        if let Some(maxvalue) = maxvalue {
            *maxvalue = maxval
        };
        if let Some(noise2) = noise2 {
            *noise2 = xnoise2
        };
        if let Some(noise3) = noise3 {
            *noise3 = xnoise3
        };
        if let Some(noise5) = noise5 {
            *noise5 = xnoise5
        };
    }
    *status
}

/// Compute mean and RMS sigma of the non-null pixels in the input array.
///
/// # Parameters
///
/// * `array`     — 2 dimensional array of image pixels
/// * `npix`      — number of pixels in the image
/// * `nullcheck` — check for null values, if true
/// * `nullvalue` — value of null pixels, if nullcheck is true
/// * `ngoodpix`  — number of non-null pixels in the image
/// * `mean`      — returned mean value of all non-null pixels
/// * `sigma`     — returned R.M.S. value of all non-null pixels
/// * `status`    — error status
fn FnMeanSigma_short(
    array: &[c_short],
    npix: usize,
    nullcheck: bool,
    nullvalue: c_short,
    ngoodpix: &mut usize,
    mean: &mut f64,
    sigma: &mut f64,
    status: &mut c_int,
) -> c_int {
    let mut ngood: usize = 0;

    let mut sum: f64 = 0.0;
    let mut sum2: f64 = 0.0;
    let mut xtemp: f64;

    let mut value: usize = 0;

    if nullcheck {
        for _ii in 0..npix {
            //for (ii = 0; ii < npix; ii+=1, value+=1) {
            if array[value] != nullvalue {
                ngood += 1;
                xtemp = array[value].into();
                sum += xtemp;
                sum2 += xtemp * xtemp;
            }
            value += 1;
        }
    } else {
        ngood = npix;
        for _ii in 0..npix {
            //for (ii = 0; ii < npix; ii+=1, value+=1) {
            xtemp = array[value].into();
            sum += xtemp;
            sum2 += xtemp * xtemp;
            value += 1;
        }
    }

    if ngood > 1 {
        *ngoodpix = ngood;
        xtemp = sum / ngood as f64;
        *mean = xtemp;
        *sigma = ((sum2 / ngood as f64) - (xtemp * xtemp)).sqrt();
    } else if ngood == 1 {
        *ngoodpix = 1;
        *mean = sum;
        *sigma = 0.0;
    } else {
        *ngoodpix = 0;
        *mean = 0.0;
        *sigma = 0.0;
    }

    *status
}

/// Compute mean and RMS sigma of the non-null pixels in the input array.
///
/// # Parameters
///
/// * `array`     — 2 dimensional array of image pixels
/// * `npix`      — number of pixels in the image
/// * `nullcheck` — check for null values, if true
/// * `nullvalue` — value of null pixels, if nullcheck is true
/// * `ngoodpix`  — number of non-null pixels in the image
/// * `mean`      — returned mean value of all non-null pixels
/// * `sigma`     — returned R.M.S. value of all non-null pixels
/// * `status`    — error status
fn FnMeanSigma_int(
    array: &[c_int],
    npix: usize,
    nullcheck: bool,
    nullvalue: c_int,
    ngoodpix: &mut usize,
    mean: &mut f64,
    sigma: &mut f64,
    status: &mut c_int,
) -> c_int {
    let mut ngood: usize = 0;

    let mut sum: f64 = 0.0;
    let mut sum2: f64 = 0.0;
    let mut xtemp: f64;

    let mut value: usize = 0;

    if nullcheck {
        for _ii in 0..npix {
            //for (ii = 0; ii < npix; ii+=1, value+=1) {
            if array[value] != nullvalue {
                ngood += 1;
                xtemp = array[value].into();
                sum += xtemp;
                sum2 += xtemp * xtemp;
            }
            value += 1;
        }
    } else {
        ngood = npix;
        for _ii in 0..npix {
            //for (ii = 0; ii < npix; ii+=1, value+=1) {
            xtemp = array[value].into();
            sum += xtemp;
            sum2 += xtemp * xtemp;
            value += 1;
        }
    }

    if ngood > 1 {
        *ngoodpix = ngood;
        xtemp = sum / ngood as f64;
        *mean = xtemp;
        *sigma = ((sum2 / ngood as f64) - (xtemp * xtemp)).sqrt();
    } else if ngood == 1 {
        *ngoodpix = 1;
        *mean = sum;
        *sigma = 0.0;
    } else {
        *ngoodpix = 0;
        *mean = 0.0;
        *sigma = 0.0;
    }

    *status
}

/// Compute mean and RMS sigma of the non-null pixels in the input array.
///
/// # Parameters
///
/// * `array`     — 2 dimensional array of image pixels
/// * `npix`      — number of pixels in the image
/// * `nullcheck` — check for null values, if true
/// * `nullvalue` — value of null pixels, if nullcheck is true
/// * `ngoodpix`  — number of non-null pixels in the image
/// * `mean`      — returned mean value of all non-null pixels
/// * `sigma`     — returned R.M.S. value of all non-null pixels
/// * `status`    — error status
fn FnMeanSigma_float(
    array: &[f32],
    npix: usize,
    nullcheck: bool,
    nullvalue: f32,
    ngoodpix: &mut usize,
    mean: &mut f64,
    sigma: &mut f64,
    status: &mut c_int,
) -> c_int {
    let mut ngood: usize = 0;

    let mut sum: f64 = 0.0;
    let mut sum2: f64 = 0.0;
    let mut xtemp: f64;

    let mut value: usize = 0;

    if nullcheck {
        for _ii in 0..npix {
            //for (ii = 0; ii < npix; ii+=1, value+=1) {
            if array[value] != nullvalue {
                ngood += 1;
                xtemp = array[value].into();
                sum += xtemp;
                sum2 += xtemp * xtemp;
            }
            value += 1;
        }
    } else {
        ngood = npix;
        for _ii in 0..npix {
            //for (ii = 0; ii < npix; ii+=1, value+=1) {
            xtemp = array[value].into();
            sum += xtemp;
            sum2 += xtemp * xtemp;
            value += 1;
        }
    }

    if ngood > 1 {
        *ngoodpix = ngood;
        xtemp = sum / ngood as f64;
        *mean = xtemp;
        *sigma = ((sum2 / ngood as f64) - (xtemp * xtemp)).sqrt();
    } else if ngood == 1 {
        *ngoodpix = 1;
        *mean = sum;
        *sigma = 0.0;
    } else {
        *ngoodpix = 0;
        *mean = 0.0;
        *sigma = 0.0;
    }

    *status
}

/// Compute mean and RMS sigma of the non-null pixels in the input array.
///
/// # Parameters
///
/// * `array`     — 2 dimensional array of image pixels
/// * `npix`      — number of pixels in the image
/// * `nullcheck` — check for null values, if true
/// * `nullvalue` — value of null pixels, if nullcheck is true
/// * `ngoodpix`  — number of non-null pixels in the image
/// * `mean`      — returned mean value of all non-null pixels
/// * `sigma`     — returned R.M.S. value of all non-null pixels
/// * `status`    — error status
fn FnMeanSigma_double(
    array: &[f64],
    npix: usize,
    nullcheck: bool,
    nullvalue: f64,
    ngoodpix: &mut usize,
    mean: &mut f64,
    sigma: &mut f64,
    status: &mut c_int,
) -> c_int {
    let mut ngood: usize = 0;

    let mut sum: f64 = 0.0;
    let mut sum2: f64 = 0.0;
    let mut xtemp: f64;

    let mut value: usize = 0;

    if nullcheck {
        for _ii in 0..npix {
            //for (ii = 0; ii < npix; ii+=1, value+=1) {
            if array[value] != nullvalue {
                ngood += 1;
                xtemp = array[value];
                sum += xtemp;
                sum2 += xtemp * xtemp;
            }
            value += 1;
        }
    } else {
        ngood = npix;
        for _ii in 0..npix {
            //for (ii = 0; ii < npix; ii+=1, value+=1) {
            xtemp = array[value];
            sum += xtemp;
            sum2 += xtemp * xtemp;
            value += 1;
        }
    }

    if ngood > 1 {
        *ngoodpix = ngood;
        xtemp = sum / ngood as f64;
        *mean = xtemp;
        *sigma = ((sum2 / ngood as f64) - (xtemp * xtemp)).sqrt();
    } else if ngood == 1 {
        *ngoodpix = 1;
        *mean = sum;
        *sigma = 0.0;
    } else {
        *ngoodpix = 0;
        *mean = 0.0;
        *sigma = 0.0;
    }

    *status
}

/// Estimate the median and background noise in the input image using 2nd, 3rd and 5th order Median Absolute Differences.
///
/// The noise in the background of the image is calculated using the MAD algorithms
/// developed for deriving the signal to noise ratio in spectra
/// (see issue #42 of the ST-ECF newsletter, <http://www.stecf.org/documents/newsletter/>)
///
/// 3rd order:  noise = 1.482602 / sqrt(6) * median (abs(2*flux(i) - flux(i-2) - flux(i+2)))
///
/// The returned estimates are the median of the values that are computed for each
/// row of the image.
///
/// # Parameters
///
/// * `array`     — 2 dimensional array of image pixels
/// * `nx`        — number of pixels in each row of the image
/// * `ny`        — number of rows in the image
/// * `nullcheck` — check for null values, if true
/// * `nullvalue` — value of null pixels, if nullcheck is true
/// * `ngood`     — number of good, non-null pixels?
/// * `minval`    — minimum non-null value
/// * `maxval`    — maximum non-null value
/// * `noise2`    — returned 2nd order MAD of all non-null pixels
/// * `noise3`    — returned 3rd order MAD of all non-null pixels
/// * `noise5`    — returned 5th order MAD of all non-null pixels
/// * `status`    — error status
fn FnNoise5_short(
    array: &[c_short],
    nx: usize,
    ny: usize,
    nullcheck: bool,
    nullvalue: c_short,
    ngood: Option<&mut usize>,
    minval: Option<&mut c_short>,
    maxval: Option<&mut c_short>,
    noise2: Option<&mut f64>,
    noise3: Option<&mut f64>,
    noise5: Option<&mut f64>,
    status: &mut c_int,
) -> c_int {
    let mut ii: usize;
    let mut _jj: usize;
    let mut nrows: usize = 0;
    let mut nrows2: usize = 0;
    let mut nvals: usize;
    let mut nvals2: usize;
    let mut ngoodpix: usize = 0;

    let mut rowpix: &[c_short];
    let mut v1: c_short;
    let mut v2: c_short;
    let mut v3: c_short;
    let mut v4: c_short;
    let mut v5: c_short;
    let mut v6: c_short;
    let mut v7: c_short;
    let mut v8: c_short;
    let mut v9: c_short;
    let mut xminval = c_short::MAX;
    let mut xmaxval = c_short::MIN;

    let xnoise2: f64;
    let xnoise3: f64;
    let xnoise5: f64;

    let mut do_range = false;

    let mut nx = nx;
    let mut ny = ny;

    if nx < 9 {
        /* treat entire array as an image with a single row */
        nx *= ny;
        ny = 1;
    }

    /* rows must have at least 9 pixels */
    if nx < 9 {
        for &item in array.iter().take(nx) {
            //for (ii = 0; ii < nx; ii+=1) {
            if nullcheck && item == nullvalue {
                continue;
            } else {
                if item < xminval {
                    xminval = item;
                }
                if item > xmaxval {
                    xmaxval = item;
                }
                ngoodpix += 1;
            }
        }

        minval.set_if_some(xminval);
        maxval.set_if_some(xmaxval);
        ngood.set_if_some(ngoodpix);
        noise2.set_if_some(0.);
        noise3.set_if_some(0.);
        noise5.set_if_some(0.);

        return *status;
    }

    /* do we need to compute the min and max value? */
    if minval.is_some() || maxval.is_some() {
        do_range = true;
    }

    /* allocate arrays used to compute the median and noise estimates */
    let mut differences2: Vec<c_int> = Vec::new();
    let mut differences3: Vec<c_int> = Vec::new();
    let mut differences5: Vec<c_int> = Vec::new();
    let mut diffs2: Vec<f64> = Vec::new();
    let mut diffs3: Vec<f64> = Vec::new();
    let mut diffs5: Vec<f64> = Vec::new();

    if differences2.try_reserve_exact(nx).is_err() {
        *status = MEMORY_ALLOCATION;
        return *status;
    } else {
        differences2.resize(nx, 0);
    }

    if differences3.try_reserve_exact(nx).is_err() {
        *status = MEMORY_ALLOCATION;
        return *status;
    } else {
        differences3.resize(nx, 0);
    }

    if differences5.try_reserve_exact(nx).is_err() {
        *status = MEMORY_ALLOCATION;
        return *status;
    } else {
        differences5.resize(nx, 0);
    }

    if diffs2.try_reserve_exact(ny).is_err() {
        *status = MEMORY_ALLOCATION;
        return *status;
    } else {
        diffs2.resize(ny, 0.0);
    }

    if diffs3.try_reserve_exact(ny).is_err() {
        *status = MEMORY_ALLOCATION;
        return *status;
    } else {
        diffs3.resize(ny, 0.0);
    }

    if diffs5.try_reserve_exact(ny).is_err() {
        *status = MEMORY_ALLOCATION;
        return *status;
    } else {
        diffs5.resize(ny, 0.0);
    }

    /* loop over each row of the image */
    for jj in 0..ny {
        //for (jj=0; jj < ny; jj+=1) {

        rowpix = &array[(jj * nx)..]; /* point to first pixel in the row */

        ii = 0;
        if nullcheck {
            while ii < nx && rowpix[ii] == nullvalue {
                ii += 1;
            }
        }

        if ii == nx {
            continue;
        } /* hit end of row */

        v1 = rowpix[ii]; /* store the good pixel value */
        ngoodpix += 1;

        if do_range {
            if v1 < xminval {
                xminval = v1;
            }
            if v1 > xmaxval {
                xmaxval = v1;
            }
        }

        ii += 1;
        if nullcheck {
            while ii < nx && rowpix[ii] == nullvalue {
                ii += 1;
            }
        }

        if ii == nx {
            continue;
        } /* hit end of row */

        v2 = rowpix[ii]; /* store the good pixel value */
        ngoodpix += 1;

        if do_range {
            if v2 < xminval {
                xminval = v2;
            }
            if v2 > xmaxval {
                xmaxval = v2;
            }
        }

        ii += 1;
        if nullcheck {
            while ii < nx && rowpix[ii] == nullvalue {
                ii += 1;
            }
        }

        if ii == nx {
            continue;
        } /* hit end of row */

        v3 = rowpix[ii]; /* store the good pixel value */
        ngoodpix += 1;

        if do_range {
            if v3 < xminval {
                xminval = v3;
            }
            if v3 > xmaxval {
                xmaxval = v3;
            }
        }

        /* find the 4nd valid pixel in row (to be skipped) */
        ii += 1;
        if nullcheck {
            while ii < nx && rowpix[ii] == nullvalue {
                ii += 1;
            }
        }

        if ii == nx {
            continue;
        } /* hit end of row */

        v4 = rowpix[ii]; /* store the good pixel value */
        ngoodpix += 1;

        if do_range {
            if v4 < xminval {
                xminval = v4;
            }
            if v4 > xmaxval {
                xmaxval = v4;
            }
        }

        /* find the 5th valid pixel in row (to be skipped) */
        ii += 1;
        if nullcheck {
            while ii < nx && rowpix[ii] == nullvalue {
                ii += 1;
            }
        }

        if ii == nx {
            continue;
        } /* hit end of row */

        v5 = rowpix[ii]; /* store the good pixel value */
        ngoodpix += 1;

        if do_range {
            if v5 < xminval {
                xminval = v5;
            }
            if v5 > xmaxval {
                xmaxval = v5;
            }
        }

        /* find the 6th valid pixel in row (to be skipped) */
        ii += 1;
        if nullcheck {
            while ii < nx && rowpix[ii] == nullvalue {
                ii += 1;
            }
        }

        if ii == nx {
            continue;
        } /* hit end of row */

        v6 = rowpix[ii]; /* store the good pixel value */
        ngoodpix += 1;

        if do_range {
            if v6 < xminval {
                xminval = v6;
            }
            if v6 > xmaxval {
                xmaxval = v6;
            }
        }

        /* find the 7th valid pixel in row (to be skipped) */
        ii += 1;
        if nullcheck {
            while ii < nx && rowpix[ii] == nullvalue {
                ii += 1;
            }
        }

        if ii == nx {
            continue;
        } /* hit end of row */

        v7 = rowpix[ii]; /* store the good pixel value */
        ngoodpix += 1;

        if do_range {
            if v7 < xminval {
                xminval = v7;
            }
            if v7 > xmaxval {
                xmaxval = v7;
            }
        }

        /* find the 8th valid pixel in row (to be skipped) */
        ii += 1;
        if nullcheck {
            while ii < nx && rowpix[ii] == nullvalue {
                ii += 1;
            }
        }

        if ii == nx {
            continue;
        } /* hit end of row */

        v8 = rowpix[ii]; /* store the good pixel value */
        ngoodpix += 1;

        if do_range {
            if v8 < xminval {
                xminval = v8;
            }
            if v8 > xmaxval {
                xmaxval = v8;
            }
        }

        /* now populate the differences arrays */
        /* for the remaining pixels in the row */
        nvals = 0;
        nvals2 = 0;
        ii += 1;
        while ii < nx {
            /* find the next valid pixel in row */
            if nullcheck {
                while ii < nx && rowpix[ii] == nullvalue {
                    ii += 1;
                }
            }

            if ii == nx {
                break;
            } /* hit end of row */

            v9 = rowpix[ii]; /* store the good pixel value */

            if do_range {
                if v9 < xminval {
                    xminval = v9;
                }
                if v9 > xmaxval {
                    xmaxval = v9;
                }
            }

            /* construct array of absolute differences */

            if !(v5 == v6 && v6 == v7) {
                differences2[nvals2] = (c_int::from(v5) - c_int::from(v7)).abs();

                nvals2 += 1;
            }

            if !(v3 == v4 && v4 == v5 && v5 == v6 && v6 == v7) {
                differences3[nvals] =
                    ((2 * c_int::from(v5)) - c_int::from(v3) - c_int::from(v7)).abs();
                differences5[nvals] =
                    ((6 * c_int::from(v5)) - (4 * c_int::from(v3)) - (4 * c_int::from(v7))
                        + c_int::from(v1)
                        + c_int::from(v9))
                    .abs();

                nvals += 1;
            } else {
                /* ignore constant background regions */
                ngoodpix += 1;
            }

            /* shift over 1 pixel */
            v1 = v2;
            v2 = v3;
            v3 = v4;
            v4 = v5;
            v5 = v6;
            v6 = v7;
            v7 = v8;
            v8 = v9;

            ii += 1;
        } /* end of loop over pixels in the row */

        /* compute the median diffs */
        /* Note that there are 8 more pixel values than there are diffs values. */
        ngoodpix += nvals;

        if nvals == 0 {
            continue; /* cannot compute medians on this row */
        } else if nvals == 1 {
            if nvals2 == 1 {
                diffs2[nrows2] = differences2[0].into();
                nrows2 += 1;
            }

            diffs3[nrows] = differences3[0].into();
            diffs5[nrows] = differences5[0].into();
        } else {
            /* quick_select returns the median MUCH faster than using qsort */
            if nvals2 > 1 {
                diffs2[nrows2] = quick_select_int(&mut differences2, nvals).into();
                nrows2 += 1;
            }

            diffs3[nrows] = quick_select_int(&mut differences3, nvals).into();
            diffs5[nrows] = quick_select_int(&mut differences5, nvals).into();
        }

        nrows += 1;
    } /* end of loop over rows */

    /* compute median of the values for each row */
    if nrows == 0 {
        xnoise3 = 0.0;
        xnoise5 = 0.0;
    } else if nrows == 1 {
        xnoise3 = diffs3[0];
        xnoise5 = diffs5[0];
    } else {
        diffs3.sort_by(f64::total_cmp);
        diffs5.sort_by(f64::total_cmp);
        xnoise3 = (diffs3[(nrows - 1) / 2] + diffs3[nrows / 2]) / 2.0;
        xnoise5 = (diffs5[(nrows - 1) / 2] + diffs5[nrows / 2]) / 2.0;
    }

    if nrows2 == 0 {
        xnoise2 = 0.0;
    } else if nrows2 == 1 {
        xnoise2 = diffs2[0];
    } else {
        diffs2.sort_by(f64::total_cmp);
        xnoise2 = (diffs2[(nrows2 - 1) / 2] + diffs2[nrows2 / 2]) / 2.0;
    }

    minval.set_if_some(xminval);
    maxval.set_if_some(xmaxval);
    ngood.set_if_some(ngoodpix);
    noise2.set_if_some(1.0483579 * xnoise2);
    noise3.set_if_some(0.6052697 * xnoise3);
    noise5.set_if_some(0.1772048 * xnoise5);

    *status
}

/// Estimate the median and background noise in the input image using 2nd, 3rd and 5th order Median Absolute Differences.
///
/// The noise in the background of the image is calculated using the MAD algorithms
/// developed for deriving the signal to noise ratio in spectra
/// (see issue #42 of the ST-ECF newsletter, <http://www.stecf.org/documents/newsletter/>)
///
/// 3rd order:  noise = 1.482602 / sqrt(6) * median (abs(2*flux(i) - flux(i-2) - flux(i+2)))
///
/// The returned estimates are the median of the values that are computed for each
/// row of the image.
///
/// # Parameters
///
/// * `array`     — 2 dimensional array of image pixels
/// * `nx`        — number of pixels in each row of the image
/// * `ny`        — number of rows in the image
/// * `nullcheck` — check for null values, if true
/// * `nullvalue` — value of null pixels, if nullcheck is true
/// * `ngood`     — number of good, non-null pixels?
/// * `minval`    — minimum non-null value
/// * `maxval`    — maximum non-null value
/// * `noise2`    — returned 2nd order MAD of all non-null pixels
/// * `noise3`    — returned 3rd order MAD of all non-null pixels
/// * `noise5`    — returned 5th order MAD of all non-null pixels
/// * `status`    — error status
fn FnNoise5_int(
    array: &[c_int],
    nx: usize,
    ny: usize,
    nullcheck: bool,
    nullvalue: c_int,
    ngood: Option<&mut usize>,
    minval: Option<&mut c_int>,
    maxval: Option<&mut c_int>,
    noise2: Option<&mut f64>,
    noise3: Option<&mut f64>,
    noise5: Option<&mut f64>,
    status: &mut c_int,
) -> c_int {
    let mut ii: usize;
    let mut _jj: usize;
    let mut nrows: usize = 0;
    let mut nrows2: usize = 0;
    let mut nvals: usize;
    let mut nvals2: usize;
    let mut ngoodpix: usize = 0;

    let mut rowpix: &[c_int];
    let mut v1: c_int;
    let mut v2: c_int;
    let mut v3: c_int;
    let mut v4: c_int;
    let mut v5: c_int;
    let mut v6: c_int;
    let mut v7: c_int;
    let mut v8: c_int;
    let mut v9: c_int;
    let mut xminval: c_int = c_int::MAX;
    let mut xmaxval = c_int::MIN;

    let xnoise2: f64;
    let xnoise3: f64;
    let xnoise5: f64;

    let mut do_range = false;

    let mut nx = nx;
    let mut ny = ny;

    if nx < 9 {
        /* treat entire array as an image with a single row */
        nx *= ny;
        ny = 1;
    }

    /* rows must have at least 9 pixels */
    if nx < 9 {
        for &item in array.iter().take(nx) {
            //for (ii = 0; ii < nx; ii+=1) {
            if nullcheck && item == nullvalue {
                continue;
            } else {
                if item < xminval {
                    xminval = item;
                }
                if item > xmaxval {
                    xmaxval = item;
                }
                ngoodpix += 1;
            }
        }

        minval.set_if_some(xminval);
        maxval.set_if_some(xmaxval);
        ngood.set_if_some(ngoodpix);
        noise2.set_if_some(0.);
        noise3.set_if_some(0.);
        noise5.set_if_some(0.);

        return *status;
    }

    /* do we need to compute the min and max value? */
    if minval.is_some() || maxval.is_some() {
        do_range = true;
    }

    /* allocate arrays used to compute the median and noise estimates */
    let mut differences2: Vec<LONGLONG> = Vec::new();
    let mut differences3: Vec<LONGLONG> = Vec::new();
    let mut differences5: Vec<LONGLONG> = Vec::new();
    let mut diffs2: Vec<f64> = Vec::new();
    let mut diffs3: Vec<f64> = Vec::new();
    let mut diffs5: Vec<f64> = Vec::new();

    if differences2.try_reserve_exact(nx).is_err() {
        *status = MEMORY_ALLOCATION;
        return *status;
    } else {
        differences2.resize(nx, 0);
    }

    if differences3.try_reserve_exact(nx).is_err() {
        *status = MEMORY_ALLOCATION;
        return *status;
    } else {
        differences3.resize(nx, 0);
    }

    if differences5.try_reserve_exact(nx).is_err() {
        *status = MEMORY_ALLOCATION;
        return *status;
    } else {
        differences5.resize(nx, 0);
    }

    if diffs2.try_reserve_exact(ny).is_err() {
        *status = MEMORY_ALLOCATION;
        return *status;
    } else {
        diffs2.resize(ny, 0.0);
    }

    if diffs3.try_reserve_exact(ny).is_err() {
        *status = MEMORY_ALLOCATION;
        return *status;
    } else {
        diffs3.resize(ny, 0.0);
    }

    if diffs5.try_reserve_exact(ny).is_err() {
        *status = MEMORY_ALLOCATION;
        return *status;
    } else {
        diffs5.resize(ny, 0.0);
    }

    /* loop over each row of the image */
    for jj in 0..ny {
        //for (jj=0; jj < ny; jj+=1) {

        rowpix = &array[(jj * nx)..]; /* point to first pixel in the row */

        ii = 0;
        if nullcheck {
            while ii < nx && rowpix[ii] == nullvalue {
                ii += 1;
            }
        }

        if ii == nx {
            continue;
        } /* hit end of row */
        v1 = rowpix[ii]; /* store the good pixel value */
        ngoodpix += 1;

        if do_range {
            if v1 < xminval {
                xminval = v1;
            }
            if v1 > xmaxval {
                xmaxval = v1;
            }
        }

        ii += 1;
        if nullcheck {
            while ii < nx && rowpix[ii] == nullvalue {
                ii += 1;
            }
        }

        if ii == nx {
            continue;
        } /* hit end of row */
        v2 = rowpix[ii]; /* store the good pixel value */
        ngoodpix += 1;

        if do_range {
            if v2 < xminval {
                xminval = v2;
            }
            if v2 > xmaxval {
                xmaxval = v2;
            }
        }

        ii += 1;
        if nullcheck {
            while ii < nx && rowpix[ii] == nullvalue {
                ii += 1;
            }
        }

        if ii == nx {
            continue;
        } /* hit end of row */
        v3 = rowpix[ii]; /* store the good pixel value */
        ngoodpix += 1;

        if do_range {
            if v3 < xminval {
                xminval = v3;
            }
            if v3 > xmaxval {
                xmaxval = v3;
            }
        }

        /* find the 4nd valid pixel in row (to be skipped) */
        ii += 1;
        if nullcheck {
            while ii < nx && rowpix[ii] == nullvalue {
                ii += 1;
            }
        }

        if ii == nx {
            continue;
        } /* hit end of row */
        v4 = rowpix[ii]; /* store the good pixel value */
        ngoodpix += 1;

        if do_range {
            if v4 < xminval {
                xminval = v4;
            }
            if v4 > xmaxval {
                xmaxval = v4;
            }
        }

        /* find the 5th valid pixel in row (to be skipped) */
        ii += 1;
        if nullcheck {
            while ii < nx && rowpix[ii] == nullvalue {
                ii += 1;
            }
        }

        if ii == nx {
            continue;
        } /* hit end of row */
        v5 = rowpix[ii]; /* store the good pixel value */
        ngoodpix += 1;

        if do_range {
            if v5 < xminval {
                xminval = v5;
            }
            if v5 > xmaxval {
                xmaxval = v5;
            }
        }

        /* find the 6th valid pixel in row (to be skipped) */
        ii += 1;
        if nullcheck {
            while ii < nx && rowpix[ii] == nullvalue {
                ii += 1;
            }
        }

        if ii == nx {
            continue;
        } /* hit end of row */
        v6 = rowpix[ii]; /* store the good pixel value */
        ngoodpix += 1;

        if do_range {
            if v6 < xminval {
                xminval = v6;
            }
            if v6 > xmaxval {
                xmaxval = v6;
            }
        }

        /* find the 7th valid pixel in row (to be skipped) */
        ii += 1;
        if nullcheck {
            while ii < nx && rowpix[ii] == nullvalue {
                ii += 1;
            }
        }

        if ii == nx {
            continue;
        } /* hit end of row */
        v7 = rowpix[ii]; /* store the good pixel value */
        ngoodpix += 1;

        if do_range {
            if v7 < xminval {
                xminval = v7;
            }
            if v7 > xmaxval {
                xmaxval = v7;
            }
        }

        /* find the 8th valid pixel in row (to be skipped) */
        ii += 1;
        if nullcheck {
            while ii < nx && rowpix[ii] == nullvalue {
                ii += 1;
            }
        }

        if ii == nx {
            continue;
        } /* hit end of row */
        v8 = rowpix[ii]; /* store the good pixel value */
        ngoodpix += 1;

        if do_range {
            if v8 < xminval {
                xminval = v8;
            }
            if v8 > xmaxval {
                xmaxval = v8;
            }
        }

        /* now populate the differences arrays */
        /* for the remaining pixels in the row */
        nvals = 0;
        nvals2 = 0;
        ii += 1;
        while ii < nx {
            /* find the next valid pixel in row */
            if nullcheck {
                while ii < nx && rowpix[ii] == nullvalue {
                    ii += 1;
                }
            }

            if ii == nx {
                break;
            } /* hit end of row */
            v9 = rowpix[ii]; /* store the good pixel value */

            if do_range {
                if v9 < xminval {
                    xminval = v9;
                }
                if v9 > xmaxval {
                    xmaxval = v9;
                }
            }

            /* construct array of absolute differences */

            if !(v5 == v6 && v6 == v7) {
                let tdiff = LONGLONG::from(v5) - LONGLONG::from(v7);
                if tdiff < 0 {
                    differences2[nvals2] = -tdiff;
                } else {
                    differences2[nvals2] = tdiff;
                }

                nvals2 += 1;
            }

            if !(v3 == v4 && v4 == v5 && v5 == v6 && v6 == v7) {
                let tdiff = (2 * LONGLONG::from(v5)) - LONGLONG::from(v3) - LONGLONG::from(v7);
                if tdiff < 0 {
                    differences3[nvals] = -tdiff;
                } else {
                    differences3[nvals] = tdiff;
                }

                let tdiff =
                    (6 * LONGLONG::from(v5)) - (4 * LONGLONG::from(v3)) - (4 * LONGLONG::from(v7))
                        + LONGLONG::from(v1)
                        + LONGLONG::from(v9);
                if tdiff < 0 {
                    differences5[nvals] = -tdiff;
                } else {
                    differences5[nvals] = tdiff;
                }

                nvals += 1;
            } else {
                /* ignore constant background regions */
                ngoodpix += 1;
            }

            /* shift over 1 pixel */
            v1 = v2;
            v2 = v3;
            v3 = v4;
            v4 = v5;
            v5 = v6;
            v6 = v7;
            v7 = v8;
            v8 = v9;

            ii += 1;
        } /* end of loop over pixels in the row */

        /* compute the median diffs */
        /* Note that there are 8 more pixel values than there are diffs values. */
        ngoodpix += nvals;

        if nvals == 0 {
            continue; /* cannot compute medians on this row */
        } else if nvals == 1 {
            if nvals2 == 1 {
                diffs2[nrows2] = differences2[0] as f64;
                nrows2 += 1;
            }

            diffs3[nrows] = differences3[0] as f64;
            diffs5[nrows] = differences5[0] as f64;
        } else {
            /* quick_select returns the median MUCH faster than using qsort */
            if nvals2 > 1 {
                diffs2[nrows2] = quick_select_longlong(&mut differences2, nvals) as f64;
                nrows2 += 1;
            }

            diffs3[nrows] = quick_select_longlong(&mut differences3, nvals) as f64;
            diffs5[nrows] = quick_select_longlong(&mut differences5, nvals) as f64;
        }

        nrows += 1;
    } /* end of loop over rows */

    /* compute median of the values for each row */
    if nrows == 0 {
        xnoise3 = 0.0;
        xnoise5 = 0.0;
    } else if nrows == 1 {
        xnoise3 = diffs3[0];
        xnoise5 = diffs5[0];
    } else {
        diffs3.sort_by(f64::total_cmp);
        diffs5.sort_by(f64::total_cmp);
        xnoise3 = (diffs3[(nrows - 1) / 2] + diffs3[nrows / 2]) / 2.0;
        xnoise5 = (diffs5[(nrows - 1) / 2] + diffs5[nrows / 2]) / 2.0;
    }

    if nrows2 == 0 {
        xnoise2 = 0.0;
    } else if nrows2 == 1 {
        xnoise2 = diffs2[0];
    } else {
        diffs2.sort_by(f64::total_cmp);
        xnoise2 = (diffs2[(nrows2 - 1) / 2] + diffs2[nrows2 / 2]) / 2.0;
    }

    minval.set_if_some(xminval);
    maxval.set_if_some(xmaxval);
    ngood.set_if_some(ngoodpix);
    noise2.set_if_some(1.0483579 * xnoise2);
    noise3.set_if_some(0.6052697 * xnoise3);
    noise5.set_if_some(0.1772048 * xnoise5);

    *status
}

/// Estimate the median and background noise in the input image using 2nd, 3rd and 5th order Median Absolute Differences.
///
/// The noise in the background of the image is calculated using the MAD algorithms
/// developed for deriving the signal to noise ratio in spectra
/// (see issue #42 of the ST-ECF newsletter, <http://www.stecf.org/documents/newsletter/>)
///
/// 3rd order:  noise = 1.482602 / sqrt(6) * median (abs(2*flux(i) - flux(i-2) - flux(i+2)))
///
/// The returned estimates are the median of the values that are computed for each
/// row of the image.
///
/// # Parameters
///
/// * `array`     — 2 dimensional array of image pixels
/// * `nx`        — number of pixels in each row of the image
/// * `ny`        — number of rows in the image
/// * `nullcheck` — check for null values, if true
/// * `nullvalue` — value of null pixels, if nullcheck is true
/// * `ngood`     — number of good, non-null pixels?
/// * `minval`    — minimum non-null value
/// * `maxval`    — maximum non-null value
/// * `noise2`    — returned 2nd order MAD of all non-null pixels
/// * `noise3`    — returned 3rd order MAD of all non-null pixels
/// * `noise5`    — returned 5th order MAD of all non-null pixels
/// * `status`    — error status
fn FnNoise5_float(
    array: &[f32],
    nx: usize,
    ny: usize,
    nullcheck: bool,
    nullvalue: f32,
    ngood: Option<&mut usize>,
    minval: Option<&mut f32>,
    maxval: Option<&mut f32>,
    noise2: Option<&mut f64>,
    noise3: Option<&mut f64>,
    noise5: Option<&mut f64>,
    status: &mut c_int,
) -> c_int {
    let mut ii: usize;
    let mut _jj: usize;
    let mut nrows: usize = 0;
    let mut nrows2: usize = 0;
    let mut nvals: usize;
    let mut nvals2: usize;
    let mut ngoodpix: usize = 0;

    let mut rowpix: &[f32];
    let mut v1: f32;
    let mut v2: f32;
    let mut v3: f32;
    let mut v4: f32;
    let mut v5: f32;
    let mut v6: f32;
    let mut v7: f32;
    let mut v8: f32;
    let mut v9: f32;
    let mut xminval: f32 = f32::MAX;
    let mut xmaxval = -f32::MAX;

    let xnoise2: f64;
    let xnoise3: f64;
    let xnoise5: f64;

    let mut do_range = false;

    let mut nx = nx;
    let mut ny = ny;

    if nx < 9 {
        /* treat entire array as an image with a single row */
        nx *= ny;
        ny = 1;
    }

    /* rows must have at least 9 pixels */
    if nx < 9 {
        for &item in array.iter().take(nx) {
            //for (ii = 0; ii < nx; ii+=1) {
            if nullcheck && item == nullvalue {
                continue;
            } else {
                if item < xminval {
                    xminval = item;
                }
                if item > xmaxval {
                    xmaxval = item;
                }
                ngoodpix += 1;
            }
        }

        minval.set_if_some(xminval);
        maxval.set_if_some(xmaxval);
        ngood.set_if_some(ngoodpix);
        noise2.set_if_some(0.);
        noise3.set_if_some(0.);
        noise5.set_if_some(0.);

        return *status;
    }

    /* do we need to compute the min and max value? */
    if minval.is_some() || maxval.is_some() {
        do_range = true;
    }

    /* allocate arrays used to compute the median and noise estimates */
    let mut differences2: Vec<f32> = Vec::new();
    let mut differences3: Vec<f32> = Vec::new();
    let mut differences5: Vec<f32> = Vec::new();
    let mut diffs2: Vec<f64> = Vec::new();
    let mut diffs3: Vec<f64> = Vec::new();
    let mut diffs5: Vec<f64> = Vec::new();

    if differences2.try_reserve_exact(nx).is_err() {
        *status = MEMORY_ALLOCATION;
        return *status;
    } else {
        differences2.resize(nx, 0.0);
    }

    if differences3.try_reserve_exact(nx).is_err() {
        *status = MEMORY_ALLOCATION;
        return *status;
    } else {
        differences3.resize(nx, 0.0);
    }

    if differences5.try_reserve_exact(nx).is_err() {
        *status = MEMORY_ALLOCATION;
        return *status;
    } else {
        differences5.resize(nx, 0.0);
    }

    if diffs2.try_reserve_exact(ny).is_err() {
        *status = MEMORY_ALLOCATION;
        return *status;
    } else {
        diffs2.resize(ny, 0.0);
    }

    if diffs3.try_reserve_exact(ny).is_err() {
        *status = MEMORY_ALLOCATION;
        return *status;
    } else {
        diffs3.resize(ny, 0.0);
    }

    if diffs5.try_reserve_exact(ny).is_err() {
        *status = MEMORY_ALLOCATION;
        return *status;
    } else {
        diffs5.resize(ny, 0.0);
    }

    /* loop over each row of the image */
    for jj in 0..ny {
        //for (jj=0; jj < ny; jj+=1) {

        rowpix = &array[(jj * nx)..]; /* point to first pixel in the row */

        ii = 0;
        if nullcheck {
            while ii < nx && rowpix[ii] == nullvalue {
                ii += 1;
            }
        }

        if ii == nx {
            continue;
        } /* hit end of row */
        v1 = rowpix[ii]; /* store the good pixel value */
        ngoodpix += 1;

        if do_range {
            if v1 < xminval {
                xminval = v1;
            }
            if v1 > xmaxval {
                xmaxval = v1;
            }
        }

        ii += 1;
        if nullcheck {
            while ii < nx && rowpix[ii] == nullvalue {
                ii += 1;
            }
        }

        if ii == nx {
            continue;
        } /* hit end of row */
        v2 = rowpix[ii]; /* store the good pixel value */
        ngoodpix += 1;

        if do_range {
            if v2 < xminval {
                xminval = v2;
            }
            if v2 > xmaxval {
                xmaxval = v2;
            }
        }

        ii += 1;
        if nullcheck {
            while ii < nx && rowpix[ii] == nullvalue {
                ii += 1;
            }
        }

        if ii == nx {
            continue;
        } /* hit end of row */
        v3 = rowpix[ii]; /* store the good pixel value */
        ngoodpix += 1;

        if do_range {
            if v3 < xminval {
                xminval = v3;
            }
            if v3 > xmaxval {
                xmaxval = v3;
            }
        }

        /* find the 4nd valid pixel in row (to be skipped) */
        ii += 1;
        if nullcheck {
            while ii < nx && rowpix[ii] == nullvalue {
                ii += 1;
            }
        }

        if ii == nx {
            continue;
        } /* hit end of row */
        v4 = rowpix[ii]; /* store the good pixel value */
        ngoodpix += 1;

        if do_range {
            if v4 < xminval {
                xminval = v4;
            }
            if v4 > xmaxval {
                xmaxval = v4;
            }
        }

        /* find the 5th valid pixel in row (to be skipped) */
        ii += 1;
        if nullcheck {
            while ii < nx && rowpix[ii] == nullvalue {
                ii += 1;
            }
        }

        if ii == nx {
            continue;
        } /* hit end of row */
        v5 = rowpix[ii]; /* store the good pixel value */
        ngoodpix += 1;

        if do_range {
            if v5 < xminval {
                xminval = v5;
            }
            if v5 > xmaxval {
                xmaxval = v5;
            }
        }

        /* find the 6th valid pixel in row (to be skipped) */
        ii += 1;
        if nullcheck {
            while ii < nx && rowpix[ii] == nullvalue {
                ii += 1;
            }
        }

        if ii == nx {
            continue;
        } /* hit end of row */
        v6 = rowpix[ii]; /* store the good pixel value */
        ngoodpix += 1;

        if do_range {
            if v6 < xminval {
                xminval = v6;
            }
            if v6 > xmaxval {
                xmaxval = v6;
            }
        }

        /* find the 7th valid pixel in row (to be skipped) */
        ii += 1;
        if nullcheck {
            while ii < nx && rowpix[ii] == nullvalue {
                ii += 1;
            }
        }

        if ii == nx {
            continue;
        } /* hit end of row */
        v7 = rowpix[ii]; /* store the good pixel value */
        ngoodpix += 1;

        if do_range {
            if v7 < xminval {
                xminval = v7;
            }
            if v7 > xmaxval {
                xmaxval = v7;
            }
        }

        /* find the 8th valid pixel in row (to be skipped) */
        ii += 1;
        if nullcheck {
            while ii < nx && rowpix[ii] == nullvalue {
                ii += 1;
            }
        }

        if ii == nx {
            continue;
        } /* hit end of row */
        v8 = rowpix[ii]; /* store the good pixel value */
        ngoodpix += 1;

        if do_range {
            if v8 < xminval {
                xminval = v8;
            }
            if v8 > xmaxval {
                xmaxval = v8;
            }
        }

        /* now populate the differences arrays */
        /* for the remaining pixels in the row */
        nvals = 0;
        nvals2 = 0;
        ii += 1;
        while ii < nx {
            /* find the next valid pixel in row */
            if nullcheck {
                while ii < nx && rowpix[ii] == nullvalue {
                    ii += 1;
                }
            }

            if ii == nx {
                break;
            } /* hit end of row */
            v9 = rowpix[ii]; /* store the good pixel value */

            if do_range {
                if v9 < xminval {
                    xminval = v9;
                }
                if v9 > xmaxval {
                    xmaxval = v9;
                }
            }

            /* construct array of absolute differences */

            if !(v5 == v6 && v6 == v7) {
                differences2[nvals2] = (v5 - v7).abs();
                nvals2 += 1;
            }

            if !(v3 == v4 && v4 == v5 && v5 == v6 && v6 == v7) {
                differences3[nvals] = ((2.0 * v5) - v3 - v7).abs();
                differences5[nvals] = ((6.0 * v5) - (4.0 * v3) - (4.0 * v7) + v1 + v9).abs();
                nvals += 1;
            } else {
                /* ignore constant background regions */
                ngoodpix += 1;
            }

            /* shift over 1 pixel */
            v1 = v2;
            v2 = v3;
            v3 = v4;
            v4 = v5;
            v5 = v6;
            v6 = v7;
            v7 = v8;
            v8 = v9;

            ii += 1;
        } /* end of loop over pixels in the row */

        /* compute the median diffs */
        /* Note that there are 8 more pixel values than there are diffs values. */
        ngoodpix += nvals;

        if nvals == 0 {
            continue; /* cannot compute medians on this row */
        } else if nvals == 1 {
            if nvals2 == 1 {
                diffs2[nrows2] = differences2[0].into();
                nrows2 += 1;
            }

            diffs3[nrows] = differences3[0].into();
            diffs5[nrows] = differences5[0].into();
        } else {
            /* quick_select returns the median MUCH faster than using qsort */
            if nvals2 > 1 {
                diffs2[nrows2] = quick_select_float(&mut differences2, nvals).into();
                nrows2 += 1;
            }

            diffs3[nrows] = quick_select_float(&mut differences3, nvals).into();
            diffs5[nrows] = quick_select_float(&mut differences5, nvals).into();
        }

        nrows += 1;
    } /* end of loop over rows */

    /* compute median of the values for each row */
    if nrows == 0 {
        xnoise3 = 0.0;
        xnoise5 = 0.0;
    } else if nrows == 1 {
        xnoise3 = diffs3[0];
        xnoise5 = diffs5[0];
    } else {
        diffs3.sort_by(f64::total_cmp);
        diffs5.sort_by(f64::total_cmp);
        xnoise3 = (diffs3[(nrows - 1) / 2] + diffs3[nrows / 2]) / 2.0;
        xnoise5 = (diffs5[(nrows - 1) / 2] + diffs5[nrows / 2]) / 2.0;
    }

    if nrows2 == 0 {
        xnoise2 = 0.0;
    } else if nrows2 == 1 {
        xnoise2 = diffs2[0];
    } else {
        diffs2.sort_by(f64::total_cmp);
        xnoise2 = (diffs2[(nrows2 - 1) / 2] + diffs2[nrows2 / 2]) / 2.0;
    }

    minval.set_if_some(xminval);
    maxval.set_if_some(xmaxval);
    ngood.set_if_some(ngoodpix);
    noise2.set_if_some(1.0483579 * xnoise2);
    noise3.set_if_some(0.6052697 * xnoise3);
    noise5.set_if_some(0.1772048 * xnoise5);

    *status
}

/// Estimate the median and background noise in the input image using 2nd, 3rd and 5th order Median Absolute Differences.
///
/// The noise in the background of the image is calculated using the MAD algorithms
/// developed for deriving the signal to noise ratio in spectra
/// (see issue #42 of the ST-ECF newsletter, <http://www.stecf.org/documents/newsletter/>)
///
/// 3rd order:  noise = 1.482602 / sqrt(6) * median (abs(2*flux(i) - flux(i-2) - flux(i+2)))
///
/// The returned estimates are the median of the values that are computed for each
/// row of the image.
///
/// # Parameters
///
/// * `array`     — 2 dimensional array of image pixels
/// * `nx`        — number of pixels in each row of the image
/// * `ny`        — number of rows in the image
/// * `nullcheck` — check for null values, if true
/// * `nullvalue` — value of null pixels, if nullcheck is true
/// * `ngood`     — number of good, non-null pixels?
/// * `minval`    — minimum non-null value
/// * `maxval`    — maximum non-null value
/// * `noise2`    — returned 2nd order MAD of all non-null pixels
/// * `noise3`    — returned 3rd order MAD of all non-null pixels
/// * `noise5`    — returned 5th order MAD of all non-null pixels
/// * `status`    — error status
fn FnNoise5_double(
    array: &[f64],
    nx: usize,
    ny: usize,
    nullcheck: bool,
    nullvalue: f64,
    ngood: Option<&mut usize>,
    minval: Option<&mut f64>,
    maxval: Option<&mut f64>,
    noise2: Option<&mut f64>,
    noise3: Option<&mut f64>,
    noise5: Option<&mut f64>,
    status: &mut c_int,
) -> c_int {
    let mut ii: usize;
    let mut _jj: usize;
    let mut nrows: usize = 0;
    let mut nrows2: usize = 0;
    let mut nvals: usize;
    let mut nvals2: usize;
    let mut ngoodpix: usize = 0;

    let mut rowpix: &[f64];
    let mut v1: f64;
    let mut v2: f64;
    let mut v3: f64;
    let mut v4: f64;
    let mut v5: f64;
    let mut v6: f64;
    let mut v7: f64;
    let mut v8: f64;
    let mut v9: f64;
    let mut xminval: f64 = f64::MAX;
    let mut xmaxval = -f64::MAX;

    let xnoise2: f64;
    let xnoise3: f64;
    let xnoise5: f64;

    let mut do_range = false;

    let mut nx = nx;
    let mut ny = ny;

    if nx < 9 {
        /* treat entire array as an image with a single row */
        nx *= ny;
        ny = 1;
    }

    /* rows must have at least 9 pixels */
    if nx < 9 {
        for &item in array.iter().take(nx) {
            //for (ii = 0; ii < nx; ii+=1) {
            if nullcheck && item == nullvalue {
                continue;
            } else {
                if item < xminval {
                    xminval = item;
                }
                if item > xmaxval {
                    xmaxval = item;
                }
                ngoodpix += 1;
            }
        }

        minval.set_if_some(xminval);
        maxval.set_if_some(xmaxval);
        ngood.set_if_some(ngoodpix);
        noise2.set_if_some(0.);
        noise3.set_if_some(0.);
        noise5.set_if_some(0.);

        return *status;
    }

    /* do we need to compute the min and max value? */
    if minval.is_some() || maxval.is_some() {
        do_range = true;
    }

    /* allocate arrays used to compute the median and noise estimates */
    let mut differences2: Vec<f64> = Vec::new();
    let mut differences3: Vec<f64> = Vec::new();
    let mut differences5: Vec<f64> = Vec::new();
    let mut diffs2: Vec<f64> = Vec::new();
    let mut diffs3: Vec<f64> = Vec::new();
    let mut diffs5: Vec<f64> = Vec::new();

    if differences2.try_reserve_exact(nx).is_err() {
        *status = MEMORY_ALLOCATION;
        return *status;
    } else {
        differences2.resize(nx, 0.0);
    }

    if differences3.try_reserve_exact(nx).is_err() {
        *status = MEMORY_ALLOCATION;
        return *status;
    } else {
        differences3.resize(nx, 0.0);
    }

    if differences5.try_reserve_exact(nx).is_err() {
        *status = MEMORY_ALLOCATION;
        return *status;
    } else {
        differences5.resize(nx, 0.0);
    }

    if diffs2.try_reserve_exact(ny).is_err() {
        *status = MEMORY_ALLOCATION;
        return *status;
    } else {
        diffs2.resize(ny, 0.0);
    }

    if diffs3.try_reserve_exact(ny).is_err() {
        *status = MEMORY_ALLOCATION;
        return *status;
    } else {
        diffs3.resize(ny, 0.0);
    }

    if diffs5.try_reserve_exact(ny).is_err() {
        *status = MEMORY_ALLOCATION;
        return *status;
    } else {
        diffs5.resize(ny, 0.0);
    }

    /* loop over each row of the image */
    for jj in 0..ny {
        //for (jj=0; jj < ny; jj+=1) {

        rowpix = &array[(jj * nx)..]; /* point to first pixel in the row */

        ii = 0;
        if nullcheck {
            while ii < nx && rowpix[ii] == nullvalue {
                ii += 1;
            }
        }

        if ii == nx {
            continue;
        } /* hit end of row */
        v1 = rowpix[ii]; /* store the good pixel value */
        ngoodpix += 1;

        if do_range {
            if v1 < xminval {
                xminval = v1;
            }
            if v1 > xmaxval {
                xmaxval = v1;
            }
        }

        ii += 1;
        if nullcheck {
            while ii < nx && rowpix[ii] == nullvalue {
                ii += 1;
            }
        }

        if ii == nx {
            continue;
        } /* hit end of row */
        v2 = rowpix[ii]; /* store the good pixel value */
        ngoodpix += 1;

        if do_range {
            if v2 < xminval {
                xminval = v2;
            }
            if v2 > xmaxval {
                xmaxval = v2;
            }
        }

        ii += 1;
        if nullcheck {
            while ii < nx && rowpix[ii] == nullvalue {
                ii += 1;
            }
        }

        if ii == nx {
            continue;
        } /* hit end of row */
        v3 = rowpix[ii]; /* store the good pixel value */
        ngoodpix += 1;

        if do_range {
            if v3 < xminval {
                xminval = v3;
            }
            if v3 > xmaxval {
                xmaxval = v3;
            }
        }

        /* find the 4nd valid pixel in row (to be skipped) */
        ii += 1;
        if nullcheck {
            while ii < nx && rowpix[ii] == nullvalue {
                ii += 1;
            }
        }

        if ii == nx {
            continue;
        } /* hit end of row */
        v4 = rowpix[ii]; /* store the good pixel value */
        ngoodpix += 1;

        if do_range {
            if v4 < xminval {
                xminval = v4;
            }
            if v4 > xmaxval {
                xmaxval = v4;
            }
        }

        /* find the 5th valid pixel in row (to be skipped) */
        ii += 1;
        if nullcheck {
            while ii < nx && rowpix[ii] == nullvalue {
                ii += 1;
            }
        }

        if ii == nx {
            continue;
        } /* hit end of row */
        v5 = rowpix[ii]; /* store the good pixel value */
        ngoodpix += 1;

        if do_range {
            if v5 < xminval {
                xminval = v5;
            }
            if v5 > xmaxval {
                xmaxval = v5;
            }
        }

        /* find the 6th valid pixel in row (to be skipped) */
        ii += 1;
        if nullcheck {
            while ii < nx && rowpix[ii] == nullvalue {
                ii += 1;
            }
        }

        if ii == nx {
            continue;
        } /* hit end of row */
        v6 = rowpix[ii]; /* store the good pixel value */
        ngoodpix += 1;

        if do_range {
            if v6 < xminval {
                xminval = v6;
            }
            if v6 > xmaxval {
                xmaxval = v6;
            }
        }

        /* find the 7th valid pixel in row (to be skipped) */
        ii += 1;
        if nullcheck {
            while ii < nx && rowpix[ii] == nullvalue {
                ii += 1;
            }
        }

        if ii == nx {
            continue;
        } /* hit end of row */
        v7 = rowpix[ii]; /* store the good pixel value */
        ngoodpix += 1;

        if do_range {
            if v7 < xminval {
                xminval = v7;
            }
            if v7 > xmaxval {
                xmaxval = v7;
            }
        }

        /* find the 8th valid pixel in row (to be skipped) */
        ii += 1;
        if nullcheck {
            while ii < nx && rowpix[ii] == nullvalue {
                ii += 1;
            }
        }

        if ii == nx {
            continue;
        } /* hit end of row */
        v8 = rowpix[ii]; /* store the good pixel value */
        ngoodpix += 1;

        if do_range {
            if v8 < xminval {
                xminval = v8;
            }
            if v8 > xmaxval {
                xmaxval = v8;
            }
        }

        /* now populate the differences arrays */
        /* for the remaining pixels in the row */
        nvals = 0;
        nvals2 = 0;
        ii += 1;
        while ii < nx {
            /* find the next valid pixel in row */
            if nullcheck {
                while ii < nx && rowpix[ii] == nullvalue {
                    ii += 1;
                }
            }

            if ii == nx {
                break;
            } /* hit end of row */
            v9 = rowpix[ii]; /* store the good pixel value */

            if do_range {
                if v9 < xminval {
                    xminval = v9;
                }
                if v9 > xmaxval {
                    xmaxval = v9;
                }
            }

            /* construct array of absolute differences */

            if !(v5 == v6 && v6 == v7) {
                differences2[nvals2] = (v5 - v7).abs();
                nvals2 += 1;
            }

            if !(v3 == v4 && v4 == v5 && v5 == v6 && v6 == v7) {
                differences3[nvals] = ((2.0 * v5) - v3 - v7).abs();
                differences5[nvals] = ((6.0 * v5) - (4.0 * v3) - (4.0 * v7) + v1 + v9).abs();
                nvals += 1;
            } else {
                /* ignore constant background regions */
                ngoodpix += 1;
            }

            /* shift over 1 pixel */
            v1 = v2;
            v2 = v3;
            v3 = v4;
            v4 = v5;
            v5 = v6;
            v6 = v7;
            v7 = v8;
            v8 = v9;

            ii += 1;
        } /* end of loop over pixels in the row */

        /* compute the median diffs */
        /* Note that there are 8 more pixel values than there are diffs values. */
        ngoodpix += nvals;

        if nvals == 0 {
            continue; /* cannot compute medians on this row */
        } else if nvals == 1 {
            if nvals2 == 1 {
                diffs2[nrows2] = differences2[0];
                nrows2 += 1;
            }

            diffs3[nrows] = differences3[0];
            diffs5[nrows] = differences5[0];
        } else {
            /* quick_select returns the median MUCH faster than using qsort */
            if nvals2 > 1 {
                diffs2[nrows2] = quick_select_double(&mut differences2, nvals);
                nrows2 += 1;
            }

            diffs3[nrows] = quick_select_double(&mut differences3, nvals);
            diffs5[nrows] = quick_select_double(&mut differences5, nvals);
        }

        nrows += 1;
    } /* end of loop over rows */

    /* compute median of the values for each row */
    if nrows == 0 {
        xnoise3 = 0.0;
        xnoise5 = 0.0;
    } else if nrows == 1 {
        xnoise3 = diffs3[0];
        xnoise5 = diffs5[0];
    } else {
        diffs3.sort_by(f64::total_cmp);
        diffs5.sort_by(f64::total_cmp);
        xnoise3 = (diffs3[(nrows - 1) / 2] + diffs3[nrows / 2]) / 2.;
        xnoise5 = (diffs5[(nrows - 1) / 2] + diffs5[nrows / 2]) / 2.;
    }

    if nrows2 == 0 {
        xnoise2 = 0.0;
    } else if nrows2 == 1 {
        xnoise2 = diffs2[0];
    } else {
        diffs2.sort_by(f64::total_cmp);
        xnoise2 = (diffs2[(nrows2 - 1) / 2] + diffs2[nrows2 / 2]) / 2.;
    }

    minval.set_if_some(xminval);
    maxval.set_if_some(xmaxval);
    ngood.set_if_some(ngoodpix);
    noise2.set_if_some(1.0483579 * xnoise2);
    noise3.set_if_some(0.6052697 * xnoise3);
    noise5.set_if_some(0.1772048 * xnoise5);

    *status
}

/// Estimate the median and background noise in the input image using 3rd order differences.
///
/// The noise in the background of the image is calculated using the 3rd order algorithm
/// developed for deriving the signal to noise ratio in spectra
/// (see issue #42 of the ST-ECF newsletter, <http://www.stecf.org/documents/newsletter/>)
///
///   noise = 1.482602 / sqrt(6) * median (abs(2*flux(i) - flux(i-2) - flux(i+2)))
///
/// The returned estimates are the median of the values that are computed for each
/// row of the image.
///
/// # Parameters
///
/// * `array`     — 2 dimensional array of image pixels
/// * `nx`        — number of pixels in each row of the image
/// * `ny`        — number of rows in the image
/// * `nullcheck` — check for null values, if true
/// * `nullvalue` — value of null pixels, if nullcheck is true
/// * `ngood`     — number of good, non-null pixels?
/// * `minval`    — minimum non-null value
/// * `maxval`    — maximum non-null value
/// * `noise`     — returned R.M.S. value of all non-null pixels
/// * `status`    — error status
fn FnNoise3_short(
    array: &[c_short],
    nx: usize,
    ny: usize,
    nullcheck: bool,
    nullvalue: c_short,
    ngood: Option<&mut usize>,
    minval: Option<&mut c_short>,
    maxval: Option<&mut c_short>,
    noise: Option<&mut f64>,
    status: &mut c_int,
) -> c_int {
    let mut ii: usize;
    let mut _jj: usize;
    let mut nrows: usize = 0;
    let mut nvals: usize;
    let mut ngoodpix: usize = 0; /* number of good, non-null pixels? */
    let mut differences: Vec<c_short> = Vec::new();
    let mut rowpix: &[c_short];
    let mut v1: c_short;
    let mut v2: c_short;
    let mut v3: c_short;
    let mut v4: c_short;
    let mut v5: c_short;
    let mut xminval: c_short = c_short::MAX; /* minimum non-null value */
    let mut xmaxval: c_short = c_short::MIN; /* maximum non-null value */
    let mut diffs: Vec<f64> = Vec::new();
    let mut xnoise: f64 = 0.0; /* returned R.M.S. value of all non-null pixels */
    let mut sigma: f64 = 0.0;

    let mut do_range = false;

    let mut nx = nx;
    let mut ny = ny;

    if nx < 5 {
        /* treat entire array as an image with a single row */
        nx *= ny;
        ny = 1;
    }

    /* rows must have at least 5 pixels */
    if nx < 5 {
        for &item in array.iter().take(nx) {
            //for (ii = 0; ii < nx; ii+=1) {
            if nullcheck && item == nullvalue {
                continue;
            } else {
                if item < xminval {
                    xminval = item;
                }
                if item > xmaxval {
                    xmaxval = item;
                }
                ngoodpix += 1;
            }
        }

        minval.set_if_some(xminval);
        maxval.set_if_some(xmaxval);
        ngood.set_if_some(ngoodpix);
        noise.set_if_some(0.0);

        return *status;
    }

    /* do we need to compute the min and max value? */
    if minval.is_some() || maxval.is_some() {
        do_range = true;
    }

    /* allocate arrays used to compute the median and noise estimates */
    if differences.try_reserve_exact(nx).is_err() {
        *status = MEMORY_ALLOCATION;
        return *status;
    } else {
        differences.resize(nx, 0);
    }

    if diffs.try_reserve_exact(ny).is_err() {
        *status = MEMORY_ALLOCATION;
        return *status;
    } else {
        diffs.resize(ny, 0.0);
    }

    /* loop over each row of the image */
    for jj in 0..ny {
        //for (jj=0; jj < ny; jj+=1) {

        rowpix = &array[(jj * nx)..]; /* point to first pixel in the row */

        ii = 0;
        if nullcheck {
            while ii < nx && rowpix[ii] == nullvalue {
                ii += 1;
            }
        }

        if ii == nx {
            continue;
        } /* hit end of row */
        v1 = rowpix[ii]; /* store the good pixel value */

        if do_range {
            if v1 < xminval {
                xminval = v1;
            }
            if v1 > xmaxval {
                xmaxval = v1;
            }
        }

        ii += 1;
        if nullcheck {
            while ii < nx && rowpix[ii] == nullvalue {
                ii += 1;
            }
        }

        if ii == nx {
            continue;
        } /* hit end of row */
        v2 = rowpix[ii]; /* store the good pixel value */

        if do_range {
            if v2 < xminval {
                xminval = v2;
            }
            if v2 > xmaxval {
                xmaxval = v2;
            }
        }

        ii += 1;
        if nullcheck {
            while ii < nx && rowpix[ii] == nullvalue {
                ii += 1;
            }
        }

        if ii == nx {
            continue;
        } /* hit end of row */
        v3 = rowpix[ii]; /* store the good pixel value */

        if do_range {
            if v3 < xminval {
                xminval = v3;
            }
            if v3 > xmaxval {
                xmaxval = v3;
            }
        }

        /* find the 4nd valid pixel in row (to be skipped) */
        ii += 1;
        if nullcheck {
            while ii < nx && rowpix[ii] == nullvalue {
                ii += 1;
            }
        }

        if ii == nx {
            continue;
        } /* hit end of row */
        v4 = rowpix[ii]; /* store the good pixel value */

        if do_range {
            if v4 < xminval {
                xminval = v4;
            }
            if v4 > xmaxval {
                xmaxval = v4;
            }
        }

        /* now populate the differences arrays */
        /* for the remaining pixels in the row */
        nvals = 0;
        ii += 1;
        while ii < nx {
            /* find the next valid pixel in row */
            if nullcheck {
                while ii < nx && rowpix[ii] == nullvalue {
                    ii += 1;
                }
            }

            if ii == nx {
                break;
            } /* hit end of row */
            v5 = rowpix[ii]; /* store the good pixel value */

            if v5 < xminval {
                xminval = v5;
            }
            if v5 > xmaxval {
                xmaxval = v5;
            }

            /* construct array of 3rd order absolute differences */
            if !(v1 == v2 && v2 == v3 && v3 == v4 && v4 == v5) {
                differences[nvals] = c_short::abs((2 * v3) - v1 - v5);
                nvals += 1;
            } else {
                /* ignore constant background regions */
                ngoodpix += 1;
            }

            /* shift over 1 pixel */
            v1 = v2;
            v2 = v3;
            v3 = v4;
            v4 = v5;

            ii += 1;
        } /* end of loop over pixels in the row */

        /* compute the 3rd order diffs */
        /* Note that there are 4 more pixel values than there are diffs values. */
        ngoodpix += nvals + 4;

        if nvals == 0 {
            continue; /* cannot compute medians on this row */
        } else if nvals == 1 {
            diffs[nrows] = differences[0].into();
        } else {
            /* quick_select returns the median MUCH faster than using qsort */
            diffs[nrows] = quick_select_short(&mut differences, nvals).into();
        }
        nrows += 1;
    } /* end of loop over rows */

    /* compute median of the values for each row */
    if nrows == 0 {
        xnoise = 0.0;
    } else if nrows == 1 {
        xnoise = diffs[0];
    } else {
        diffs.sort_by(f64::total_cmp);
        xnoise = (diffs[(nrows - 1) / 2] + diffs[nrows / 2]) / 2.;

        FnMeanSigma_double(
            &diffs,
            nrows,
            false,
            0.0,
            &mut 0,
            &mut xnoise,
            &mut sigma,
            status,
        );

        /* do a 4.5 sigma rejection of outliers */
        let mut jj = 0;
        let mut ii = 0;
        sigma *= 4.5;
        while ii < nrows {
            if f64::abs(diffs[ii] - xnoise) <= sigma {
                if jj != ii {
                    diffs[jj] = diffs[ii];
                }
                jj += 1;
            }
            ii += 1;
        }

        if ii != jj {
            FnMeanSigma_double(
                &diffs,
                jj,
                false,
                0.0,
                &mut 0,
                &mut xnoise,
                &mut sigma,
                status,
            );
        }
    }

    minval.set_if_some(xminval);
    maxval.set_if_some(xmaxval);
    ngood.set_if_some(ngoodpix);
    noise.set_if_some(0.6052697 * xnoise);

    *status
}

/// Estimate the median and background noise in the input image using 3rd order differences.
///
/// The noise in the background of the image is calculated using the 3rd order algorithm
/// developed for deriving the signal to noise ratio in spectra
/// (see issue #42 of the ST-ECF newsletter, <http://www.stecf.org/documents/newsletter/>)
///
///   noise = 1.482602 / sqrt(6) * median (abs(2*flux(i) - flux(i-2) - flux(i+2)))
///
/// The returned estimates are the median of the values that are computed for each
/// row of the image.
///
/// # Parameters
///
/// * `array`     — 2 dimensional array of image pixels
/// * `nx`        — number of pixels in each row of the image
/// * `ny`        — number of rows in the image
/// * `nullcheck` — check for null values, if true
/// * `nullvalue` — value of null pixels, if nullcheck is true
/// * `ngood`     — number of good, non-null pixels?
/// * `minval`    — minimum non-null value
/// * `maxval`    — maximum non-null value
/// * `noise`     — returned R.M.S. value of all non-null pixels
/// * `status`    — error status
fn FnNoise3_int(
    array: &[c_int],
    nx: usize,
    ny: usize,
    nullcheck: bool,
    nullvalue: c_int,
    ngood: Option<&mut usize>,
    minval: Option<&mut c_int>,
    maxval: Option<&mut c_int>,
    noise: Option<&mut f64>,
    status: &mut c_int,
) -> c_int {
    let mut ii: usize;
    let mut _jj: usize;
    let mut nrows: usize = 0;
    let mut nvals: usize;
    let mut ngoodpix: usize = 0; /* number of good, non-null pixels? */
    let mut differences: Vec<c_int> = Vec::new();
    let mut rowpix: &[c_int];
    let mut v1: c_int;
    let mut v2: c_int;
    let mut v3: c_int;
    let mut v4: c_int;
    let mut v5: c_int;
    let mut xminval: c_int = c_int::MAX; /* minimum non-null value */
    let mut xmaxval: c_int = c_int::MIN; /* maximum non-null value */
    let mut diffs: Vec<f64> = Vec::new();
    let mut xnoise: f64 = 0.0; /* returned R.M.S. value of all non-null pixels */
    let mut sigma: f64 = 0.0;

    let mut do_range = false;

    let mut nx = nx;
    let mut ny = ny;

    if nx < 5 {
        /* treat entire array as an image with a single row */
        nx *= ny;
        ny = 1;
    }

    /* rows must have at least 5 pixels */
    if nx < 5 {
        for &item in array.iter().take(nx) {
            //for (ii = 0; ii < nx; ii+=1) {
            if nullcheck && item == nullvalue {
                continue;
            } else {
                if item < xminval {
                    xminval = item;
                }
                if item > xmaxval {
                    xmaxval = item;
                }
                ngoodpix += 1;
            }
        }

        minval.set_if_some(xminval);
        maxval.set_if_some(xmaxval);
        ngood.set_if_some(ngoodpix);
        noise.set_if_some(0.0);

        return *status;
    }

    /* do we need to compute the min and max value? */
    if minval.is_some() || maxval.is_some() {
        do_range = true;
    }

    /* allocate arrays used to compute the median and noise estimates */
    if differences.try_reserve_exact(nx).is_err() {
        *status = MEMORY_ALLOCATION;
        return *status;
    } else {
        differences.resize(nx, 0);
    }

    if diffs.try_reserve_exact(ny).is_err() {
        *status = MEMORY_ALLOCATION;
        return *status;
    } else {
        diffs.resize(ny, 0.0);
    }

    /* loop over each row of the image */
    for jj in 0..ny {
        //for (jj=0; jj < ny; jj+=1) {

        rowpix = &array[(jj * nx)..]; /* point to first pixel in the row */

        ii = 0;
        if nullcheck {
            while ii < nx && rowpix[ii] == nullvalue {
                ii += 1;
            }
        }

        if ii == nx {
            continue;
        } /* hit end of row */
        v1 = rowpix[ii]; /* store the good pixel value */

        if do_range {
            if v1 < xminval {
                xminval = v1;
            }
            if v1 > xmaxval {
                xmaxval = v1;
            }
        }

        ii += 1;
        if nullcheck {
            while ii < nx && rowpix[ii] == nullvalue {
                ii += 1;
            }
        }

        if ii == nx {
            continue;
        } /* hit end of row */
        v2 = rowpix[ii]; /* store the good pixel value */

        if do_range {
            if v2 < xminval {
                xminval = v2;
            }
            if v2 > xmaxval {
                xmaxval = v2;
            }
        }

        ii += 1;
        if nullcheck {
            while ii < nx && rowpix[ii] == nullvalue {
                ii += 1;
            }
        }

        if ii == nx {
            continue;
        } /* hit end of row */
        v3 = rowpix[ii]; /* store the good pixel value */

        if do_range {
            if v3 < xminval {
                xminval = v3;
            }
            if v3 > xmaxval {
                xmaxval = v3;
            }
        }

        /* find the 4nd valid pixel in row (to be skipped) */
        ii += 1;
        if nullcheck {
            while ii < nx && rowpix[ii] == nullvalue {
                ii += 1;
            }
        }

        if ii == nx {
            continue;
        } /* hit end of row */
        v4 = rowpix[ii]; /* store the good pixel value */

        if do_range {
            if v4 < xminval {
                xminval = v4;
            }
            if v4 > xmaxval {
                xmaxval = v4;
            }
        }

        /* now populate the differences arrays */
        /* for the remaining pixels in the row */
        nvals = 0;
        ii += 1;
        while ii < nx {
            /* find the next valid pixel in row */
            if nullcheck {
                while ii < nx && rowpix[ii] == nullvalue {
                    ii += 1;
                }
            }

            if ii == nx {
                break;
            } /* hit end of row */
            v5 = rowpix[ii]; /* store the good pixel value */

            if v5 < xminval {
                xminval = v5;
            }
            if v5 > xmaxval {
                xmaxval = v5;
            }

            /* construct array of 3rd order absolute differences */
            if !(v1 == v2 && v2 == v3 && v3 == v4 && v4 == v5) {
                differences[nvals] = c_int::abs((2 * v3) - v1 - v5);
                nvals += 1;
            } else {
                /* ignore constant background regions */
                ngoodpix += 1;
            }

            /* shift over 1 pixel */
            v1 = v2;
            v2 = v3;
            v3 = v4;
            v4 = v5;

            ii += 1;
        } /* end of loop over pixels in the row */

        /* compute the 3rd order diffs */
        /* Note that there are 4 more pixel values than there are diffs values. */
        ngoodpix += nvals + 4;

        if nvals == 0 {
            continue; /* cannot compute medians on this row */
        } else if nvals == 1 {
            diffs[nrows] = differences[0].into();
        } else {
            /* quick_select returns the median MUCH faster than using qsort */
            diffs[nrows] = quick_select_int(&mut differences, nvals).into();
        }
        nrows += 1;
    } /* end of loop over rows */

    /* compute median of the values for each row */
    if nrows == 0 {
        xnoise = 0.0;
    } else if nrows == 1 {
        xnoise = diffs[0];
    } else {
        diffs.sort_by(f64::total_cmp);
        xnoise = (diffs[(nrows - 1) / 2] + diffs[nrows / 2]) / 2.;

        FnMeanSigma_double(
            &diffs,
            nrows,
            false,
            0.0,
            &mut 0,
            &mut xnoise,
            &mut sigma,
            status,
        );

        /* do a 4.5 sigma rejection of outliers */
        let mut jj = 0;
        let mut ii = 0;
        sigma *= 4.5;

        while ii < nrows {
            if f64::abs(diffs[ii] - xnoise) <= sigma {
                if jj != ii {
                    diffs[jj] = diffs[ii];
                }
                jj += 1;
            }
            ii += 1;
        }

        if ii != jj {
            FnMeanSigma_double(
                &diffs,
                jj,
                false,
                0.0,
                &mut 0,
                &mut xnoise,
                &mut sigma,
                status,
            );
        }
    }

    minval.set_if_some(xminval);
    maxval.set_if_some(xmaxval);
    ngood.set_if_some(ngoodpix);
    noise.set_if_some(0.6052697 * xnoise);

    *status
}

/// Estimate the median and background noise in the input image using 3rd order differences.
///
/// The noise in the background of the image is calculated using the 3rd order algorithm
/// developed for deriving the signal to noise ratio in spectra
/// (see issue #42 of the ST-ECF newsletter, <http://www.stecf.org/documents/newsletter/>)
///
///   noise = 1.482602 / sqrt(6) * median (abs(2*flux(i) - flux(i-2) - flux(i+2)))
///
/// The returned estimates are the median of the values that are computed for each
/// row of the image.
///
/// # Parameters
///
/// * `array`     — 2 dimensional array of image pixels
/// * `nx`        — number of pixels in each row of the image
/// * `ny`        — number of rows in the image
/// * `nullcheck` — check for null values, if true
/// * `nullvalue` — value of null pixels, if nullcheck is true
/// * `ngood`     — number of good, non-null pixels?
/// * `minval`    — minimum non-null value
/// * `maxval`    — maximum non-null value
/// * `noise`     — returned R.M.S. value of all non-null pixels
/// * `status`    — error status
fn FnNoise3_float(
    array: &[f32],
    nx: usize,
    ny: usize,
    nullcheck: bool,
    nullvalue: f32,
    ngood: Option<&mut usize>,
    minval: Option<&mut f32>,
    maxval: Option<&mut f32>,
    noise: Option<&mut f64>,
    status: &mut c_int,
) -> c_int {
    let mut ii: usize;
    let mut _jj: usize;
    let mut nrows: usize = 0;
    let mut nvals: usize;
    let mut ngoodpix: usize = 0; /* number of good, non-null pixels? */
    let mut differences: Vec<f32> = Vec::new();
    let mut rowpix: &[f32];
    let mut v1: f32;
    let mut v2: f32;
    let mut v3: f32;
    let mut v4: f32;
    let mut v5: f32;
    let mut xminval: f32 = f32::MAX; /* minimum non-null value */
    let mut xmaxval: f32 = -f32::MAX; /* maximum non-null value */
    let mut diffs: Vec<f64> = Vec::new();
    let mut xnoise: f64 = 0.0; /* returned R.M.S. value of all non-null pixels */

    let mut do_range = false;

    let mut nx = nx;
    let mut ny = ny;

    if nx < 5 {
        /* treat entire array as an image with a single row */
        nx *= ny;
        ny = 1;
    }

    /* rows must have at least 5 pixels */
    if nx < 5 {
        for &item in array.iter().take(nx) {
            //for (ii = 0; ii < nx; ii+=1) {
            if nullcheck && item == nullvalue {
                continue;
            } else {
                if item < xminval {
                    xminval = item;
                }
                if item > xmaxval {
                    xmaxval = item;
                }
                ngoodpix += 1;
            }
        }

        minval.set_if_some(xminval);
        maxval.set_if_some(xmaxval);
        ngood.set_if_some(ngoodpix);
        noise.set_if_some(0.0);

        return *status;
    }

    /* do we need to compute the min and max value? */
    if minval.is_some() || maxval.is_some() {
        do_range = true;
    }

    /* allocate arrays used to compute the median and noise estimates */
    if noise.is_some() {
        if differences.try_reserve_exact(nx).is_err() {
            *status = MEMORY_ALLOCATION;
            return *status;
        } else {
            differences.resize(nx, 0.0);
        }

        if diffs.try_reserve_exact(ny).is_err() {
            *status = MEMORY_ALLOCATION;
            return *status;
        } else {
            diffs.resize(ny, 0.0);
        }
    }

    /* loop over each row of the image */
    for jj in 0..ny {
        //for (jj=0; jj < ny; jj+=1) {

        rowpix = &array[(jj * nx)..]; /* point to first pixel in the row */

        ii = 0;
        if nullcheck {
            while ii < nx && rowpix[ii] == nullvalue {
                ii += 1;
            }
        }

        if ii == nx {
            continue;
        } /* hit end of row */
        v1 = rowpix[ii]; /* store the good pixel value */

        if do_range {
            if v1 < xminval {
                xminval = v1;
            }
            if v1 > xmaxval {
                xmaxval = v1;
            }
        }

        ii += 1;
        if nullcheck {
            while ii < nx && rowpix[ii] == nullvalue {
                ii += 1;
            }
        }

        if ii == nx {
            continue;
        } /* hit end of row */
        v2 = rowpix[ii]; /* store the good pixel value */

        if do_range {
            if v2 < xminval {
                xminval = v2;
            }
            if v2 > xmaxval {
                xmaxval = v2;
            }
        }

        ii += 1;
        if nullcheck {
            while ii < nx && rowpix[ii] == nullvalue {
                ii += 1;
            }
        }

        if ii == nx {
            continue;
        } /* hit end of row */
        v3 = rowpix[ii]; /* store the good pixel value */

        if do_range {
            if v3 < xminval {
                xminval = v3;
            }
            if v3 > xmaxval {
                xmaxval = v3;
            }
        }

        /* find the 4nd valid pixel in row (to be skipped) */
        ii += 1;
        if nullcheck {
            while ii < nx && rowpix[ii] == nullvalue {
                ii += 1;
            }
        }

        if ii == nx {
            continue;
        } /* hit end of row */
        v4 = rowpix[ii]; /* store the good pixel value */

        if do_range {
            if v4 < xminval {
                xminval = v4;
            }
            if v4 > xmaxval {
                xmaxval = v4;
            }
        }

        /* now populate the differences arrays */
        /* for the remaining pixels in the row */
        nvals = 0;
        ii += 1;
        while ii < nx {
            /* find the next valid pixel in row */
            if nullcheck {
                while ii < nx && rowpix[ii] == nullvalue {
                    ii += 1;
                }
            }

            if ii == nx {
                break;
            } /* hit end of row */
            v5 = rowpix[ii]; /* store the good pixel value */

            if v5 < xminval {
                xminval = v5;
            }
            if v5 > xmaxval {
                xmaxval = v5;
            }

            /* construct array of 3rd order absolute differences */
            if noise.is_some() {
                if !(v1 == v2 && v2 == v3 && v3 == v4 && v4 == v5) {
                    differences[nvals] = f32::abs((2.0 * v3) - v1 - v5);
                    nvals += 1;
                } else {
                    /* ignore constant background regions */
                    ngoodpix += 1;
                }
            } else {
                /* just increment the number of non-null pixels */
                ngoodpix += 1;
            }

            /* shift over 1 pixel */
            v1 = v2;
            v2 = v3;
            v3 = v4;
            v4 = v5;

            ii += 1;
        } /* end of loop over pixels in the row */

        /* compute the 3rd order diffs */
        /* Note that there are 4 more pixel values than there are diffs values. */
        ngoodpix += nvals + 4;

        if noise.is_some() {
            if nvals == 0 {
                continue; /* cannot compute medians on this row */
            } else if nvals == 1 {
                diffs[nrows] = differences[0].into();
            } else {
                /* quick_select returns the median MUCH faster than using qsort */
                diffs[nrows] = quick_select_float(&mut differences, nvals).into();
            }
        }
        nrows += 1;
    } /* end of loop over rows */

    /* compute median of the values for each row */
    if noise.is_some() {
        if nrows == 0 {
            xnoise = 0.0;
        } else if nrows == 1 {
            xnoise = diffs[0];
        } else {
            diffs.sort_by(f64::total_cmp);
            xnoise = (diffs[(nrows - 1) / 2] + diffs[nrows / 2]) / 2.;
        }
    }

    minval.set_if_some(xminval);
    maxval.set_if_some(xmaxval);
    ngood.set_if_some(ngoodpix);
    noise.set_if_some(0.6052697 * xnoise);

    *status
}

/// Estimate the median and background noise in the input image using 3rd order differences.
///
///The noise in the background of the image is calculated using the 3rd order algorithm
///developed for deriving the signal to noise ratio in spectra
///(see issue #42 of the ST-ECF newsletter, <http://www.stecf.org/documents/newsletter/>)
///
///  noise = 1.482602 / sqrt(6) * median (abs(2*flux(i) - flux(i-2) - flux(i+2)))
///
///The returned estimates are the median of the values that are computed for each
///row of the image.
///
/// # Parameters
///
/// * `array`     — 2 dimensional array of image pixels
/// * `nx`        — number of pixels in each row of the image
/// * `ny`        — number of rows in the image
/// * `nullcheck` — check for null values, if true
/// * `nullvalue` — value of null pixels, if nullcheck is true
/// * `ngood`     — number of good, non-null pixels?
/// * `minval`    — minimum non-null value
/// * `maxval`    — maximum non-null value
/// * `noise`     — returned R.M.S. value of all non-null pixels
/// * `status`    — error status
fn FnNoise3_double(
    array: &[f64],
    nx: usize,
    ny: usize,
    nullcheck: bool,
    nullvalue: f64,
    ngood: Option<&mut usize>,
    minval: Option<&mut f64>,
    maxval: Option<&mut f64>,
    noise: Option<&mut f64>,
    status: &mut c_int,
) -> c_int {
    let mut ii: usize;
    let mut _jj: usize;
    let mut nrows: usize = 0;
    let mut nvals: usize;
    let mut ngoodpix: usize = 0; /* number of good, non-null pixels? */
    let mut differences: Vec<f64> = Vec::new();
    let mut rowpix: &[f64];
    let mut v1: f64;
    let mut v2: f64;
    let mut v3: f64;
    let mut v4: f64;
    let mut v5: f64;
    let mut xminval: f64 = f64::MAX; /* minimum non-null value */
    let mut xmaxval: f64 = -f64::MAX; /* maximum non-null value */
    let mut diffs: Vec<f64> = Vec::new();
    let mut xnoise: f64 = 0.0; /* returned R.M.S. value of all non-null pixels */

    let mut do_range = false;

    let mut nx = nx;
    let mut ny = ny;

    if nx < 5 {
        /* treat entire array as an image with a single row */
        nx *= ny;
        ny = 1;
    }

    /* rows must have at least 5 pixels */
    if nx < 5 {
        for &item in array.iter().take(nx) {
            //for (ii = 0; ii < nx; ii+=1) {
            if nullcheck && item == nullvalue {
                continue;
            } else {
                if item < xminval {
                    xminval = item;
                }
                if item > xmaxval {
                    xmaxval = item;
                }
                ngoodpix += 1;
            }
        }

        minval.set_if_some(xminval);
        maxval.set_if_some(xmaxval);
        ngood.set_if_some(ngoodpix);
        noise.set_if_some(0.0);

        return *status;
    }

    /* do we need to compute the min and max value? */
    if minval.is_some() || maxval.is_some() {
        do_range = true;
    }

    /* allocate arrays used to compute the median and noise estimates */
    if noise.is_some() {
        if differences.try_reserve_exact(nx).is_err() {
            *status = MEMORY_ALLOCATION;
            return *status;
        } else {
            differences.resize(nx, 0.0);
        }

        if diffs.try_reserve_exact(ny).is_err() {
            *status = MEMORY_ALLOCATION;
            return *status;
        } else {
            diffs.resize(ny, 0.0);
        }
    }

    /* loop over each row of the image */
    for jj in 0..ny {
        //for (jj=0; jj < ny; jj+=1) {

        rowpix = &array[(jj * nx)..]; /* point to first pixel in the row */

        ii = 0;
        if nullcheck {
            while ii < nx && rowpix[ii] == nullvalue {
                ii += 1;
            }
        }

        if ii == nx {
            continue;
        } /* hit end of row */
        v1 = rowpix[ii]; /* store the good pixel value */

        if do_range {
            if v1 < xminval {
                xminval = v1;
            }
            if v1 > xmaxval {
                xmaxval = v1;
            }
        }

        ii += 1;
        if nullcheck {
            while ii < nx && rowpix[ii] == nullvalue {
                ii += 1;
            }
        }

        if ii == nx {
            continue;
        } /* hit end of row */
        v2 = rowpix[ii]; /* store the good pixel value */

        if do_range {
            if v2 < xminval {
                xminval = v2;
            }
            if v2 > xmaxval {
                xmaxval = v2;
            }
        }

        ii += 1;
        if nullcheck {
            while ii < nx && rowpix[ii] == nullvalue {
                ii += 1;
            }
        }

        if ii == nx {
            continue;
        } /* hit end of row */
        v3 = rowpix[ii]; /* store the good pixel value */

        if do_range {
            if v3 < xminval {
                xminval = v3;
            }
            if v3 > xmaxval {
                xmaxval = v3;
            }
        }

        /* find the 4nd valid pixel in row (to be skipped) */
        ii += 1;
        if nullcheck {
            while ii < nx && rowpix[ii] == nullvalue {
                ii += 1;
            }
        }

        if ii == nx {
            continue;
        } /* hit end of row */
        v4 = rowpix[ii]; /* store the good pixel value */

        if do_range {
            if v4 < xminval {
                xminval = v4;
            }
            if v4 > xmaxval {
                xmaxval = v4;
            }
        }

        /* now populate the differences arrays */
        /* for the remaining pixels in the row */
        nvals = 0;
        ii += 1;
        while ii < nx {
            /* find the next valid pixel in row */
            if nullcheck {
                while ii < nx && rowpix[ii] == nullvalue {
                    ii += 1;
                }
            }

            if ii == nx {
                break;
            } /* hit end of row */
            v5 = rowpix[ii]; /* store the good pixel value */

            if v5 < xminval {
                xminval = v5;
            }
            if v5 > xmaxval {
                xmaxval = v5;
            }

            /* construct array of 3rd order absolute differences */
            if noise.is_some() {
                if !(v1 == v2 && v2 == v3 && v3 == v4 && v4 == v5) {
                    differences[nvals] = f64::abs((2.0 * v3) - v1 - v5);
                    nvals += 1;
                } else {
                    /* ignore constant background regions */
                    ngoodpix += 1;
                }
            } else {
                /* just increment the number of non-null pixels */
                ngoodpix += 1;
            }

            /* shift over 1 pixel */
            v1 = v2;
            v2 = v3;
            v3 = v4;
            v4 = v5;

            ii += 1;
        } /* end of loop over pixels in the row */

        /* compute the 3rd order diffs */
        /* Note that there are 4 more pixel values than there are diffs values. */
        ngoodpix += nvals + 4;

        if noise.is_some() {
            if nvals == 0 {
                continue; /* cannot compute medians on this row */
            } else if nvals == 1 {
                diffs[nrows] = differences[0];
            } else {
                /* quick_select returns the median MUCH faster than using qsort */
                diffs[nrows] = quick_select_double(&mut differences, nvals);
            }
        }
        nrows += 1;
    } /* end of loop over rows */

    /* compute median of the values for each row */
    if noise.is_some() {
        if nrows == 0 {
            xnoise = 0.0;
        } else if nrows == 1 {
            xnoise = diffs[0];
        } else {
            diffs.sort_by(f64::total_cmp);
            xnoise = (diffs[(nrows - 1) / 2] + diffs[nrows / 2]) / 2.;
        }
    }

    minval.set_if_some(xminval);
    maxval.set_if_some(xmaxval);
    ngood.set_if_some(ngoodpix);
    noise.set_if_some(0.6052697 * xnoise);

    *status
}

/// Estimate the background noise in the input image using sigma of 1st order differences.
///
///   noise = 1.0 / sqrt(2) * rms of (`flux[i]` - `flux[i-1]`)
///
/// The returned estimate is the median of the values that are computed for each
/// row of the image.
///
/// # Parameters
///
/// * `array`     — 2 dimensional array of image pixels
/// * `nx`        — number of pixels in each row of the image
/// * `ny`        — number of rows in the image
/// * `nullcheck` — check for null values, if true
/// * `nullvalue` — value of null pixels, if nullcheck is true
/// * `noise`     — returned R.M.S. value of all non-null pixels
/// * `status`    — error status
fn FnNoise1_short(
    array: &[c_short],
    nx: usize,
    ny: usize,
    nullcheck: bool,
    nullvalue: c_short,
    noise: &mut f64,
    status: &mut c_int,
) -> c_int {
    let mut ii: usize;
    let mut _jj: usize;
    let mut kk: usize;
    let mut nrows: usize = 0;
    let mut nvals: usize;

    let mut differences: Vec<c_short> = Vec::new();
    let mut rowpix: &[c_short];
    let mut v1: c_short;

    let mut xnoise: f64 = 0.0;
    let mut mean: f64 = 0.0;
    let mut stdev: f64 = 0.0;

    let mut diffs: Vec<f64> = Vec::new();

    /* rows must have at least 3 pixels to estimate noise */
    if nx < 3 {
        *noise = 0.0;
        return *status;
    }

    /* allocate arrays used to compute the median and noise estimates */
    if differences.try_reserve_exact(nx).is_err() {
        *status = MEMORY_ALLOCATION;
        return *status;
    } else {
        differences.resize(nx, 0);
    }

    if diffs.try_reserve_exact(ny).is_err() {
        *status = MEMORY_ALLOCATION;
        return *status;
    } else {
        diffs.resize(ny, 0.0);
    }

    /* loop over each row of the image */
    for jj in 0..ny {
        //for (jj=0; jj < ny; jj+=1) {

        rowpix = &array[(jj * nx)..]; /* point to first pixel in the row */

        ii = 0;
        if nullcheck {
            while ii < nx && rowpix[ii] == nullvalue {
                ii += 1;
            }
        }

        if ii == nx {
            continue;
        } /* hit end of row */
        v1 = rowpix[ii]; /* store the good pixel value */

        /* now continue populating the differences arrays */
        /* for the remaining pixels in the row */
        nvals = 0;
        ii += 1;
        while ii < nx {
            /* find the next valid pixel in row */
            if nullcheck {
                while ii < nx && rowpix[ii] == nullvalue {
                    ii += 1;
                }
            }

            if ii == nx {
                break;
            } /* hit end of row */

            /* construct array of 1st order differences */
            differences[nvals] = v1 - rowpix[ii];

            nvals += 1;
            /* shift over 1 pixel */
            v1 = rowpix[ii];

            ii += 1;
        } /* end of loop over pixels in the row */

        if nvals < 2 {
            continue;
        } else {
            FnMeanSigma_short(
                &differences,
                nvals,
                false,
                0,
                &mut 0,
                &mut mean,
                &mut stdev,
                status,
            );

            if stdev > 0.0 {
                for _iter in 0..NITER {
                    //for (iter = 0;  iter < NITER;  iter+=1) {
                    kk = 0;
                    for ii in 0..nvals {
                        //for (ii = 0;  ii < nvals;  ii+=1) {
                        if (f64::from(differences[ii]) - mean).abs() < SIGMA_CLIP * stdev {
                            if kk < ii {
                                differences[kk] = differences[ii];
                            }
                            kk += 1;
                        }
                    }
                    if kk == nvals {
                        break;
                    }

                    nvals = kk;
                    FnMeanSigma_short(
                        &differences,
                        nvals,
                        false,
                        0,
                        &mut 0,
                        &mut mean,
                        &mut stdev,
                        status,
                    );
                }
            }

            diffs[nrows] = stdev;
            nrows += 1;
        }
    } /* end of loop over rows */

    /* compute median of the values for each row */
    if nrows == 0 {
        xnoise = 0.0;
    } else if nrows == 1 {
        xnoise = diffs[0];
    } else {
        diffs.sort_by(f64::total_cmp);
        xnoise = (diffs[(nrows - 1) / 2] + diffs[nrows / 2]) / 2.;
    }

    *noise = xnoise * core::f64::consts::FRAC_1_SQRT_2;

    *status
}

/// Estimate the background noise in the input image using sigma of 1st order differences.
///
///   noise = 1.0 / sqrt(2) * rms of (`flux[i]` - `flux[i-1]`)
///
/// The returned estimate is the median of the values that are computed for each
/// row of the image.
///
/// # Parameters
///
/// * `array`     — 2 dimensional array of image pixels
/// * `nx`        — number of pixels in each row of the image
/// * `ny`        — number of rows in the image
/// * `nullcheck` — check for null values, if true
/// * `nullvalue` — value of null pixels, if nullcheck is true
/// * `noise`     — returned R.M.S. value of all non-null pixels
/// * `status`    — error status
fn FnNoise1_int(
    array: &[c_int],
    nx: usize,
    ny: usize,
    nullcheck: bool,
    nullvalue: c_int,
    noise: &mut f64,
    status: &mut c_int,
) -> c_int {
    let mut ii: usize;
    let mut _jj: usize;
    let mut kk: usize;
    let mut nrows: usize = 0;
    let mut nvals: usize;

    let mut differences: Vec<c_int> = Vec::new();
    let mut rowpix: &[c_int];
    let mut v1: c_int;

    let mut xnoise: f64 = 0.0;
    let mut mean: f64 = 0.0;
    let mut stdev: f64 = 0.0;

    let mut diffs: Vec<f64> = Vec::new();

    /* rows must have at least 3 pixels to estimate noise */
    if nx < 3 {
        *noise = 0.0;
        return *status;
    }

    /* allocate arrays used to compute the median and noise estimates */
    if differences.try_reserve_exact(nx).is_err() {
        *status = MEMORY_ALLOCATION;
        return *status;
    } else {
        differences.resize(nx, 0);
    }

    if diffs.try_reserve_exact(ny).is_err() {
        *status = MEMORY_ALLOCATION;
        return *status;
    } else {
        diffs.resize(ny, 0.0);
    }

    /* loop over each row of the image */
    for jj in 0..ny {
        //for (jj=0; jj < ny; jj+=1) {

        rowpix = &array[(jj * nx)..]; /* point to first pixel in the row */

        ii = 0;
        if nullcheck {
            while ii < nx && rowpix[ii] == nullvalue {
                ii += 1;
            }
        }

        if ii == nx {
            continue;
        } /* hit end of row */
        v1 = rowpix[ii]; /* store the good pixel value */

        /* now continue populating the differences arrays */
        /* for the remaining pixels in the row */
        nvals = 0;
        ii += 1;
        while ii < nx {
            /* find the next valid pixel in row */
            if nullcheck {
                while ii < nx && rowpix[ii] == nullvalue {
                    ii += 1;
                }
            }

            if ii == nx {
                break;
            } /* hit end of row */

            /* construct array of 1st order differences */
            differences[nvals] = v1 - rowpix[ii];

            nvals += 1;
            /* shift over 1 pixel */
            v1 = rowpix[ii];

            ii += 1;
        } /* end of loop over pixels in the row */

        if nvals < 2 {
            continue;
        } else {
            FnMeanSigma_int(
                &differences,
                nvals,
                false,
                0,
                &mut 0,
                &mut mean,
                &mut stdev,
                status,
            );

            if stdev > 0.0 {
                for _iter in 0..NITER {
                    //for (iter = 0;  iter < NITER;  iter+=1) {
                    kk = 0;
                    for ii in 0..nvals {
                        //for (ii = 0;  ii < nvals;  ii+=1) {
                        if (f64::from(differences[ii]) - mean).abs() < SIGMA_CLIP * stdev {
                            if kk < ii {
                                differences[kk] = differences[ii];
                            }
                            kk += 1;
                        }
                    }
                    if kk == nvals {
                        break;
                    }

                    nvals = kk;
                    FnMeanSigma_int(
                        &differences,
                        nvals,
                        false,
                        0,
                        &mut 0,
                        &mut mean,
                        &mut stdev,
                        status,
                    );
                }
            }

            diffs[nrows] = stdev;
            nrows += 1;
        }
    } /* end of loop over rows */

    /* compute median of the values for each row */
    if nrows == 0 {
        xnoise = 0.0;
    } else if nrows == 1 {
        xnoise = diffs[0];
    } else {
        diffs.sort_by(f64::total_cmp);
        xnoise = (diffs[(nrows - 1) / 2] + diffs[nrows / 2]) / 2.;
    }

    *noise = xnoise * core::f64::consts::FRAC_1_SQRT_2;

    *status
}

/// Estimate the background noise in the input image using sigma of 1st order differences.
///
///   noise = 1.0 / sqrt(2) * rms of (`flux[i]` - `flux[i-1]`)
///
/// The returned estimate is the median of the values that are computed for each
/// row of the image.
///
/// # Parameters
///
/// * `array`     — 2 dimensional array of image pixels
/// * `nx`        — number of pixels in each row of the image
/// * `ny`        — number of rows in the image
/// * `nullcheck` — check for null values, if true
/// * `nullvalue` — value of null pixels, if nullcheck is true
/// * `noise`     — returned R.M.S. value of all non-null pixels
/// * `status`    — error status
fn FnNoise1_float(
    array: &[f32],
    nx: usize,
    ny: usize,
    nullcheck: bool,
    nullvalue: f32,
    noise: &mut f64,
    status: &mut c_int,
) -> c_int {
    let mut ii: usize;
    let mut _jj: usize;
    let mut kk: usize;
    let mut nrows: usize = 0;
    let mut nvals: usize;

    let mut differences: Vec<f32> = Vec::new();
    let mut rowpix: &[f32];
    let mut v1: f32;

    let mut xnoise: f64 = 0.0;
    let mut mean: f64 = 0.0;
    let mut stdev: f64 = 0.0;

    let mut diffs: Vec<f64> = Vec::new();

    /* rows must have at least 3 pixels to estimate noise */
    if nx < 3 {
        *noise = 0.0;
        return *status;
    }

    /* allocate arrays used to compute the median and noise estimates */
    if differences.try_reserve_exact(nx).is_err() {
        *status = MEMORY_ALLOCATION;
        return *status;
    } else {
        differences.resize(nx, 0.0);
    }

    if diffs.try_reserve_exact(ny).is_err() {
        *status = MEMORY_ALLOCATION;
        return *status;
    } else {
        diffs.resize(ny, 0.0);
    }

    /* loop over each row of the image */
    for jj in 0..ny {
        //for (jj=0; jj < ny; jj+=1) {

        rowpix = &array[(jj * nx)..]; /* point to first pixel in the row */

        ii = 0;
        if nullcheck {
            while ii < nx && rowpix[ii] == nullvalue {
                ii += 1;
            }
        }

        if ii == nx {
            continue;
        } /* hit end of row */
        v1 = rowpix[ii]; /* store the good pixel value */

        /* now continue populating the differences arrays */
        /* for the remaining pixels in the row */
        nvals = 0;
        ii += 1;
        while ii < nx {
            /* find the next valid pixel in row */
            if nullcheck {
                while ii < nx && rowpix[ii] == nullvalue {
                    ii += 1;
                }
            }

            if ii == nx {
                break;
            } /* hit end of row */

            /* construct array of 1st order differences */
            differences[nvals] = v1 - rowpix[ii];

            nvals += 1;
            /* shift over 1 pixel */
            v1 = rowpix[ii];

            ii += 1;
        } /* end of loop over pixels in the row */

        if nvals < 2 {
            continue;
        } else {
            FnMeanSigma_float(
                &differences,
                nvals,
                false,
                0.0,
                &mut 0,
                &mut mean,
                &mut stdev,
                status,
            );

            if stdev > 0.0 {
                for _iter in 0..NITER {
                    //for (iter = 0;  iter < NITER;  iter+=1) {
                    kk = 0;
                    for ii in 0..nvals {
                        //for (ii = 0;  ii < nvals;  ii+=1) {
                        if (f64::from(differences[ii]) - mean).abs() < SIGMA_CLIP * stdev {
                            if kk < ii {
                                differences[kk] = differences[ii];
                            }
                            kk += 1;
                        }
                    }
                    if kk == nvals {
                        break;
                    }

                    nvals = kk;
                    FnMeanSigma_float(
                        &differences,
                        nvals,
                        false,
                        0.0,
                        &mut 0,
                        &mut mean,
                        &mut stdev,
                        status,
                    );
                }
            }

            diffs[nrows] = stdev;
            nrows += 1;
        }
    } /* end of loop over rows */

    /* compute median of the values for each row */
    if nrows == 0 {
        xnoise = 0.0;
    } else if nrows == 1 {
        xnoise = diffs[0];
    } else {
        diffs.sort_by(f64::total_cmp);
        xnoise = (diffs[(nrows - 1) / 2] + diffs[nrows / 2]) / 2.;
    }

    *noise = xnoise * core::f64::consts::FRAC_1_SQRT_2;

    *status
}

/// Estimate the background noise in the input image using sigma of 1st order differences.
///
///   noise = 1.0 / sqrt(2) * rms of (`flux[i]` - `flux[i-1]`)
///
/// The returned estimate is the median of the values that are computed for each
/// row of the image.
///
/// # Parameters
///
/// * `array`     — 2 dimensional array of image pixels
/// * `nx`        — number of pixels in each row of the image
/// * `ny`        — number of rows in the image
/// * `nullcheck` — check for null values, if true
/// * `nullvalue` — value of null pixels, if nullcheck is true
/// * `noise`     — returned R.M.S. value of all non-null pixels
/// * `status`    — error status
fn FnNoise1_double(
    array: &[f64],
    nx: usize,
    ny: usize,
    nullcheck: bool,
    nullvalue: f64,
    noise: &mut f64,
    status: &mut c_int,
) -> c_int {
    let mut ii: usize;
    let mut _jj: usize;
    let mut kk: usize;
    let mut nrows: usize = 0;
    let mut nvals: usize;

    let mut differences: Vec<f64> = Vec::new();
    let mut rowpix: &[f64];
    let mut v1: f64;

    let mut xnoise: f64 = 0.0;
    let mut mean: f64 = 0.0;
    let mut stdev: f64 = 0.0;

    let mut diffs: Vec<f64> = Vec::new();

    /* rows must have at least 3 pixels to estimate noise */
    if nx < 3 {
        *noise = 0.0;
        return *status;
    }

    /* allocate arrays used to compute the median and noise estimates */
    if differences.try_reserve_exact(nx).is_err() {
        *status = MEMORY_ALLOCATION;
        return *status;
    } else {
        differences.resize(nx, 0.0);
    }

    if diffs.try_reserve_exact(ny).is_err() {
        *status = MEMORY_ALLOCATION;
        return *status;
    } else {
        diffs.resize(ny, 0.0);
    }

    /* loop over each row of the image */
    for jj in 0..ny {
        //for (jj=0; jj < ny; jj+=1) {

        rowpix = &array[(jj * nx)..]; /* point to first pixel in the row */

        ii = 0;
        if nullcheck {
            while ii < nx && rowpix[ii] == nullvalue {
                ii += 1;
            }
        }

        if ii == nx {
            continue;
        } /* hit end of row */
        v1 = rowpix[ii]; /* store the good pixel value */

        /* now continue populating the differences arrays */
        /* for the remaining pixels in the row */
        nvals = 0;
        ii += 1;
        while ii < nx {
            /* find the next valid pixel in row */
            if nullcheck {
                while ii < nx && rowpix[ii] == nullvalue {
                    ii += 1;
                }
            }

            if ii == nx {
                break;
            } /* hit end of row */

            /* construct array of 1st order differences */
            differences[nvals] = v1 - rowpix[ii];

            nvals += 1;
            /* shift over 1 pixel */
            v1 = rowpix[ii];

            ii += 1;
        } /* end of loop over pixels in the row */

        if nvals < 2 {
            continue;
        } else {
            FnMeanSigma_double(
                &differences,
                nvals,
                false,
                0.0,
                &mut 0,
                &mut mean,
                &mut stdev,
                status,
            );

            if stdev > 0.0 {
                for _iter in 0..NITER {
                    //for (iter = 0;  iter < NITER;  iter+=1) {
                    kk = 0;
                    for ii in 0..nvals {
                        //for (ii = 0;  ii < nvals;  ii+=1) {
                        if (differences[ii] - mean).abs() < SIGMA_CLIP * stdev {
                            if kk < ii {
                                differences[kk] = differences[ii];
                            }
                            kk += 1;
                        }
                    }
                    if kk == nvals {
                        break;
                    }

                    nvals = kk;
                    FnMeanSigma_double(
                        &differences,
                        nvals,
                        false,
                        0.0,
                        &mut 0,
                        &mut mean,
                        &mut stdev,
                        status,
                    );
                }
            }

            diffs[nrows] = stdev;
            nrows += 1;
        }
    } /* end of loop over rows */

    /* compute median of the values for each row */
    if nrows == 0 {
        xnoise = 0.0;
    } else if nrows == 1 {
        xnoise = diffs[0];
    } else {
        diffs.sort_by(f64::total_cmp);
        xnoise = (diffs[(nrows - 1) / 2] + diffs[nrows / 2]) / 2.;
    }

    *noise = xnoise * core::f64::consts::FRAC_1_SQRT_2;

    *status
}

/// Selects the `n`th smallest element of `arr`, reordering it in the process.
///
/// Based on the Quickselect algorithm described in *Numerical Recipes in C*,
/// Second Edition, Cambridge University Press, 1992, Section 8.5, ISBN
/// 0-521-43108-5. This code by Nicolas Devillard, 1998, public domain.
pub fn quick_select_float(arr: &mut [f32], n: usize) -> f32 {
    let mut middle: usize;
    let mut ll: usize;
    let mut hh: usize;

    let mut low: usize = 0;
    let mut high: usize = n - 1;
    let median: usize = (low + high) / 2;
    loop {
        //for (;;) {
        if high <= low {
            /* One element only */
            return arr[median];
        }

        if high == low + 1 {
            /* Two elements only */
            if arr[low] > arr[high] {
                arr.swap(low, high);
            }
            return arr[median];
        }

        /* Find median of low, middle and high items; swap into position low */
        middle = (low + high) / 2;
        if arr[middle] > arr[high] {
            arr.swap(middle, high);
        }
        if arr[low] > arr[high] {
            arr.swap(low, high);
        }
        if arr[middle] > arr[low] {
            arr.swap(middle, low);
        }

        /* Swap low item (now in position middle) into position (low+1) */
        arr.swap(middle, low + 1);

        /* Nibble from each end towards middle, swapping items when stuck */
        ll = low + 1;
        hh = high;
        loop {
            //for (;;) {
            loop {
                ll += 1;
                if arr[low] <= arr[ll] {
                    break;
                }
            }
            loop {
                hh -= 1;
                if arr[hh] <= arr[low] {
                    break;
                }
            }

            if hh < ll {
                break;
            }
            arr.swap(ll, hh);
        }

        /* Swap middle item (in position low) back into correct position */
        arr.swap(low, hh);

        /* Re-set active partition */
        if hh <= median {
            low = ll;
        }
        if hh >= median {
            high = hh - 1;
        }
    }
}

/// Selects the `n`th smallest element of `arr`, reordering it in the process.
///
/// Based on the Quickselect algorithm described in *Numerical Recipes in C*,
/// Second Edition, Cambridge University Press, 1992, Section 8.5, ISBN
/// 0-521-43108-5. This code by Nicolas Devillard, 1998, public domain.
pub fn quick_select_short(arr: &mut [i16], n: usize) -> i16 {
    let mut middle: usize;
    let mut ll: usize;
    let mut hh: usize;

    let mut low: usize = 0;
    let mut high: usize = n - 1;
    let median: usize = (low + high) / 2;
    loop {
        //for (;;) {
        if high <= low {
            /* One element only */
            return arr[median];
        }

        if high == low + 1 {
            /* Two elements only */
            if arr[low] > arr[high] {
                arr.swap(low, high);
            }
            return arr[median];
        }

        /* Find median of low, middle and high items; swap into position low */
        middle = (low + high) / 2;
        if arr[middle] > arr[high] {
            arr.swap(middle, high);
        }
        if arr[low] > arr[high] {
            arr.swap(low, high);
        }
        if arr[middle] > arr[low] {
            arr.swap(middle, low);
        }

        /* Swap low item (now in position middle) into position (low+1) */
        arr.swap(middle, low + 1);

        /* Nibble from each end towards middle, swapping items when stuck */
        ll = low + 1;
        hh = high;
        loop {
            //for (;;) {
            loop {
                ll += 1;
                if arr[low] <= arr[ll] {
                    break;
                }
            }
            loop {
                hh -= 1;
                if arr[hh] <= arr[low] {
                    break;
                }
            }

            if hh < ll {
                break;
            }
            arr.swap(ll, hh);
        }

        /* Swap middle item (in position low) back into correct position */
        arr.swap(low, hh);

        /* Re-set active partition */
        if hh <= median {
            low = ll;
        }
        if hh >= median {
            high = hh - 1;
        }
    }
}

/// Selects the `n`th smallest element of `arr`, reordering it in the process.
///
/// Based on the Quickselect algorithm described in *Numerical Recipes in C*,
/// Second Edition, Cambridge University Press, 1992, Section 8.5, ISBN
/// 0-521-43108-5. This code by Nicolas Devillard, 1998, public domain.
pub fn quick_select_int(arr: &mut [i32], n: usize) -> i32 {
    let mut middle: usize;
    let mut ll: usize;
    let mut hh: usize;

    let mut low: usize = 0;
    let mut high: usize = n - 1;
    let median: usize = (low + high) / 2;
    loop {
        //for (;;) {
        if high <= low {
            /* One element only */
            return arr[median];
        }

        if high == low + 1 {
            /* Two elements only */
            if arr[low] > arr[high] {
                arr.swap(low, high);
            }
            return arr[median];
        }

        /* Find median of low, middle and high items; swap into position low */
        middle = (low + high) / 2;
        if arr[middle] > arr[high] {
            arr.swap(middle, high);
        }
        if arr[low] > arr[high] {
            arr.swap(low, high);
        }
        if arr[middle] > arr[low] {
            arr.swap(middle, low);
        }

        /* Swap low item (now in position middle) into position (low+1) */
        arr.swap(middle, low + 1);

        /* Nibble from each end towards middle, swapping items when stuck */
        ll = low + 1;
        hh = high;
        loop {
            //for (;;) {
            loop {
                ll += 1;
                if arr[low] <= arr[ll] {
                    break;
                }
            }
            loop {
                hh -= 1;
                if arr[hh] <= arr[low] {
                    break;
                }
            }

            if hh < ll {
                break;
            }
            arr.swap(ll, hh);
        }

        /* Swap middle item (in position low) back into correct position */
        arr.swap(low, hh);

        /* Re-set active partition */
        if hh <= median {
            low = ll;
        }
        if hh >= median {
            high = hh - 1;
        }
    }
}

/// Selects the `n`th smallest element of `arr`, reordering it in the process.
///
/// Based on the Quickselect algorithm described in *Numerical Recipes in C*,
/// Second Edition, Cambridge University Press, 1992, Section 8.5, ISBN
/// 0-521-43108-5. This code by Nicolas Devillard, 1998, public domain.
pub fn quick_select_longlong(arr: &mut [i64], n: usize) -> i64 {
    let mut middle: usize;
    let mut ll: usize;
    let mut hh: usize;

    let mut low: usize = 0;
    let mut high: usize = n - 1;
    let median: usize = (low + high) / 2;
    loop {
        //for (;;) {
        if high <= low {
            /* One element only */
            return arr[median];
        }

        if high == low + 1 {
            /* Two elements only */
            if arr[low] > arr[high] {
                arr.swap(low, high);
            }
            return arr[median];
        }

        /* Find median of low, middle and high items; swap into position low */
        middle = (low + high) / 2;
        if arr[middle] > arr[high] {
            arr.swap(middle, high);
        }
        if arr[low] > arr[high] {
            arr.swap(low, high);
        }
        if arr[middle] > arr[low] {
            arr.swap(middle, low);
        }

        /* Swap low item (now in position middle) into position (low+1) */
        arr.swap(middle, low + 1);

        /* Nibble from each end towards middle, swapping items when stuck */
        ll = low + 1;
        hh = high;
        loop {
            //for (;;) {
            loop {
                ll += 1;
                if arr[low] <= arr[ll] {
                    break;
                }
            }
            loop {
                hh -= 1;
                if arr[hh] <= arr[low] {
                    break;
                }
            }

            if hh < ll {
                break;
            }

            arr.swap(ll, hh);
        }

        /* Swap middle item (in position low) back into correct position */
        arr.swap(low, hh);

        /* Re-set active partition */
        if hh <= median {
            low = ll;
        }
        if hh >= median {
            high = hh - 1;
        }
    }
}

/// Selects the `n`th smallest element of `arr`, reordering it in the process.
///
/// Based on the Quickselect algorithm described in *Numerical Recipes in C*,
/// Second Edition, Cambridge University Press, 1992, Section 8.5, ISBN
/// 0-521-43108-5. This code by Nicolas Devillard, 1998, public domain.
pub fn quick_select_double(arr: &mut [f64], n: usize) -> f64 {
    let mut middle: usize;
    let mut ll: usize;
    let mut hh: usize;

    let mut low: usize = 0;
    let mut high: usize = n - 1;
    let median: usize = (low + high) / 2;
    loop {
        //for (;;) {
        if high <= low {
            /* One element only */
            return arr[median];
        }

        if high == low + 1 {
            /* Two elements only */
            if arr[low] > arr[high] {
                arr.swap(low, high);
            }
            return arr[median];
        }

        /* Find median of low, middle and high items; swap into position low */
        middle = (low + high) / 2;
        if arr[middle] > arr[high] {
            arr.swap(middle, high);
        }
        if arr[low] > arr[high] {
            arr.swap(low, high);
        }
        if arr[middle] > arr[low] {
            arr.swap(middle, low);
        }

        /* Swap low item (now in position middle) into position (low+1) */
        arr.swap(middle, low + 1);

        /* Nibble from each end towards middle, swapping items when stuck */
        ll = low + 1;
        hh = high;
        loop {
            //for (;;) {
            loop {
                ll += 1;
                if arr[low] <= arr[ll] {
                    break;
                }
            }
            loop {
                hh -= 1;
                if arr[hh] <= arr[low] {
                    break;
                }
            }

            if hh < ll {
                break;
            }

            arr.swap(ll, hh);
        }

        /* Swap middle item (in position low) back into correct position */
        arr.swap(low, hh);

        /* Re-set active partition */
        if hh <= median {
            low = ll;
        }
        if hh >= median {
            high = hh - 1;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fitsio::{NO_DITHER, SUBTRACTIVE_DITHER_1, SUBTRACTIVE_DITHER_2};
    use libc::{c_int, c_long, c_short};

    fn dither(v: c_int) -> DitherType {
        match v {
            SUBTRACTIVE_DITHER_1 => DitherType::SubtractiveDither1,
            SUBTRACTIVE_DITHER_2 => DitherType::SubtractiveDither2,
            _ => DitherType::NoDither,
        }
    }

    /// Read back the int that the in-place float quantizer stored in a f32 slot.
    fn idata_f(fdata: &[f32], i: usize) -> i32 {
        fdata[i].to_bits() as i32
    }

    /// Read back the int that the in-place double quantizer stored.
    /// The double quantizer reinterprets the `&mut [f64]` slice as `&mut [i32]`
    /// (via bytemuck::cast_slice_mut) and writes idata[i] there, so index `i`
    /// maps to bytes 4*i..4*i+4 of the f64 buffer.
    fn idata_d(ddata: &[f64], i: usize) -> i32 {
        let idata: &[i32] = bytemuck::cast_slice(ddata);
        idata[i]
    }

    /*
     * Test basic float quantization with gradient data
     */
    #[test]
    fn test_quantize_float_basic() {
        let mut fdata = [0.0f32; 64];

        /* Create gradient data with some noise */
        for i in 0..64 {
            fdata[i] = (i as f64 * 10.0 + (i % 3) as f64 * 0.5) as f32;
        }
        let orig = fdata;

        let mut bscale = 0.0;
        let mut bzero = 0.0;
        let mut iminval = 0;
        let mut imaxval = 0;

        /* Quantize with no dithering, negative qlevel for absolute quant */
        let result = fits_quantize_float_inplace(
            0,
            &mut fdata,
            8,
            8,
            false,
            0.0,
            -1.0,
            dither(NO_DITHER),
            &mut bscale,
            &mut bzero,
            &mut iminval,
            &mut imaxval,
        );

        assert_eq!(result, 1); /* Should succeed */
        assert_ne!(bscale, 0.0);

        /* Verify roundtrip approximate equality */
        for i in 0..64 {
            let recovered = idata_f(&fdata, i) as f64 * bscale + bzero;
            let diff = (recovered - orig[i] as f64).abs();
            assert!(diff <= bscale); /* Error should be less than 1 quant */
        }
    }

    /*
     * Test basic double quantization
     */
    #[test]
    fn test_quantize_double_basic() {
        let mut ddata = [0.0f64; 64];

        for i in 0..64 {
            ddata[i] = i as f64 * 10.0 + (i % 3) as f64 * 0.5;
        }
        let orig = ddata;

        let mut bscale = 0.0;
        let mut bzero = 0.0;
        let mut iminval = 0;
        let mut imaxval = 0;

        let result = fits_quantize_double_inplace(
            0,
            &mut ddata,
            8,
            8,
            false,
            0.0,
            -1.0,
            dither(NO_DITHER),
            &mut bscale,
            &mut bzero,
            &mut iminval,
            &mut imaxval,
        );

        assert_eq!(result, 1);
        assert_ne!(bscale, 0.0);

        for i in 0..64 {
            let recovered = idata_d(&ddata, i) as f64 * bscale + bzero;
            let diff = (recovered - orig[i]).abs();
            assert!(diff <= bscale);
        }
    }

    /*
     * Test quantization with dithering method 1
     */
    #[test]
    fn test_quantize_dither1() {
        let mut fdata = [0.0f32; 64];
        for i in 0..64 {
            fdata[i] = (i as f64 * 5.0 + (i % 5) as f64 * 0.3) as f32;
        }

        let mut bscale = 0.0;
        let mut bzero = 0.0;
        let mut iminval = 0;
        let mut imaxval = 0;

        let result = fits_quantize_float_inplace(
            1,
            &mut fdata,
            8,
            8,
            false,
            0.0,
            -1.0,
            dither(SUBTRACTIVE_DITHER_1),
            &mut bscale,
            &mut bzero,
            &mut iminval,
            &mut imaxval,
        );

        assert_eq!(result, 1);
        assert_ne!(bscale, 0.0);
    }

    /*
     * Test quantization with dithering method 2
     */
    #[test]
    fn test_quantize_dither2() {
        let mut fdata = [0.0f32; 64];
        for i in 0..64 {
            if i % 8 == 0 {
                fdata[i] = 0.0;
            } else {
                fdata[i] = (i as f64 * 5.0 + (i % 5) as f64 * 0.3) as f32;
            }
        }

        let mut bscale = 0.0;
        let mut bzero = 0.0;
        let mut iminval = 0;
        let mut imaxval = 0;

        let result = fits_quantize_float_inplace(
            1,
            &mut fdata,
            8,
            8,
            false,
            0.0,
            -1.0,
            dither(SUBTRACTIVE_DITHER_2),
            &mut bscale,
            &mut bzero,
            &mut iminval,
            &mut imaxval,
        );

        assert_eq!(result, 1);
        assert_ne!(bscale, 0.0);
    }

    /*
     * Test quantization with null values
     */
    #[test]
    fn test_quantize_with_nulls() {
        let mut fdata = [0.0f32; 64];
        let null_value = -999.0f32;
        for i in 0..64 {
            if i % 10 == 0 {
                fdata[i] = null_value;
            } else {
                fdata[i] = (i as f64 * 10.0 + (i % 3) as f64 * 0.5) as f32;
            }
        }

        let mut bscale = 0.0;
        let mut bzero = 0.0;
        let mut iminval = 0;
        let mut imaxval = 0;

        let result = fits_quantize_float_inplace(
            0,
            &mut fdata,
            8,
            8,
            true,
            null_value,
            -1.0,
            dither(NO_DITHER),
            &mut bscale,
            &mut bzero,
            &mut iminval,
            &mut imaxval,
        );

        assert_eq!(result, 1);

        /* Null values should be preserved as NULL_VALUE (-2147483647) */
        for i in 0..64 {
            if i % 10 == 0 {
                assert_eq!(idata_f(&fdata, i), NULL_VALUE);
            }
        }
    }

    /*
     * Test uniform data (should return 0 - no quantization needed)
     */
    #[test]
    fn test_quantize_uniform() {
        let mut fdata = [42.0f32; 64];

        let mut bscale = 0.0;
        let mut bzero = 0.0;
        let mut iminval = 0;
        let mut imaxval = 0;

        /* With positive qlevel, uniform data has zero noise -> delta=0 */
        let result = fits_quantize_float_inplace(
            0,
            &mut fdata,
            8,
            8,
            false,
            0.0,
            4.0,
            dither(NO_DITHER),
            &mut bscale,
            &mut bzero,
            &mut iminval,
            &mut imaxval,
        );

        /* Should return 0 (no quantization) because delta=0 */
        assert_eq!(result, 0);
    }

    /*
     * Test quantization with negative qlevel (absolute quantization)
     */
    #[test]
    fn test_quantize_absolute() {
        let mut fdata = [0.0f32; 64];
        for i in 0..64 {
            fdata[i] = (i as f64 * 100.0) as f32;
        }

        let mut bscale = 0.0;
        let mut bzero = 0.0;
        let mut iminval = 0;
        let mut imaxval = 0;

        /* Negative qlevel = absolute quantization level */
        let result = fits_quantize_float_inplace(
            0,
            &mut fdata,
            8,
            8,
            false,
            0.0,
            -10.0,
            dither(NO_DITHER),
            &mut bscale,
            &mut bzero,
            &mut iminval,
            &mut imaxval,
        );

        assert_eq!(result, 1);
        /* bscale should be close to 10 (the absolute level) */
        assert!((bscale - 10.0_f64).abs() <= 0.001);
    }

    /*
     * Test small image (nx * ny <= 1 should return 0)
     */
    #[test]
    fn test_quantize_small() {
        let mut fdata = [42.0f32; 1];

        let mut bscale = 0.0;
        let mut bzero = 0.0;
        let mut iminval = 0;
        let mut imaxval = 0;

        let result = fits_quantize_float_inplace(
            0,
            &mut fdata,
            1,
            1,
            false,
            0.0,
            -1.0,
            dither(NO_DITHER),
            &mut bscale,
            &mut bzero,
            &mut iminval,
            &mut imaxval,
        );

        assert_eq!(result, 0);
        assert_eq!(bscale, 1.0);
        assert_eq!(bzero, 0.0);
    }

    /*
     * Test larger image with noise
     */
    #[test]
    fn test_quantize_large() {
        let mut fdata = vec![0.0f32; 256 * 256];

        for i in 0..256 * 256 {
            let x = i % 256;
            let y = i / 256;
            fdata[i] = (x as f64 + y as f64 * 100.0 + (i % 7) as f64 * 0.1) as f32;
        }

        let mut bscale = 0.0;
        let mut bzero = 0.0;
        let mut iminval = 0;
        let mut imaxval = 0;

        let result = fits_quantize_float_inplace(
            0,
            &mut fdata,
            256,
            256,
            false,
            0.0,
            -1.0,
            dither(NO_DITHER),
            &mut bscale,
            &mut bzero,
            &mut iminval,
            &mut imaxval,
        );

        assert_eq!(result, 1);
        assert_ne!(bscale, 0.0);
    }

    /*
     * Test double quantization with null values and dithering
     */
    #[test]
    fn test_quantize_double_nulls_dither() {
        let mut ddata = [0.0f64; 64];
        let null_value = -999.0f64;
        for i in 0..64 {
            if i % 8 == 0 {
                ddata[i] = null_value;
            } else {
                ddata[i] = i as f64 * 10.0 + (i % 3) as f64 * 0.5;
            }
        }

        let mut bscale = 0.0;
        let mut bzero = 0.0;
        let mut iminval = 0;
        let mut imaxval = 0;

        let result = fits_quantize_double_inplace(
            1,
            &mut ddata,
            8,
            8,
            true,
            null_value,
            -1.0,
            dither(SUBTRACTIVE_DITHER_1),
            &mut bscale,
            &mut bzero,
            &mut iminval,
            &mut imaxval,
        );

        assert_eq!(result, 1);

        for i in 0..64 {
            if i % 8 == 0 {
                assert_eq!(idata_d(&ddata, i), NULL_VALUE);
            }
        }
    }

    /*
     * Test fits_img_stats_short
     */
    #[test]
    fn test_img_stats_short_basic() {
        let mut array = [0i16; 64];
        for i in 0..64 {
            array[i] = (i * 10 + (i % 3)) as i16;
        }

        let mut status = 0;
        let mut ngoodpix: c_long = 0;
        let mut minvalue: c_short = 0;
        let mut maxvalue: c_short = 0;
        let mut mean = 0.0;
        let mut sigma = 0.0;
        let mut noise1 = 0.0;
        let mut noise2 = 0.0;
        let mut noise3 = 0.0;
        let mut noise5 = 0.0;

        fits_img_stats_short_safe(
            &array,
            8,
            8,
            false,
            0,
            Some(&mut ngoodpix),
            Some(&mut minvalue),
            Some(&mut maxvalue),
            Some(&mut mean),
            Some(&mut sigma),
            Some(&mut noise1),
            Some(&mut noise2),
            Some(&mut noise3),
            Some(&mut noise5),
            &mut status,
        );

        assert_eq!(status, 0);
        assert_eq!(ngoodpix, 64);
        assert_eq!(minvalue, 0);
        assert_eq!(maxvalue, 630); /* 63 * 10 + 0 (63 % 3 = 0) */
    }

    /*
     * Test fits_img_stats_short with null values
     */
    #[test]
    fn test_img_stats_short_nulls() {
        let mut array = [0i16; 64];
        let null_value: c_short = -999;
        for i in 0..64 {
            if i % 8 == 0 {
                array[i] = null_value;
            } else {
                array[i] = (i * 10 + (i % 3)) as i16;
            }
        }

        let mut status = 0;
        let mut ngoodpix: c_long = 0;
        let mut minvalue: c_short = 0;
        let mut maxvalue: c_short = 0;
        let mut mean = 0.0;
        let mut sigma = 0.0;
        let mut noise1 = 0.0;
        let mut noise2 = 0.0;
        let mut noise3 = 0.0;
        let mut noise5 = 0.0;

        fits_img_stats_short_safe(
            &array,
            8,
            8,
            true,
            null_value,
            Some(&mut ngoodpix),
            Some(&mut minvalue),
            Some(&mut maxvalue),
            Some(&mut mean),
            Some(&mut sigma),
            Some(&mut noise1),
            Some(&mut noise2),
            Some(&mut noise3),
            Some(&mut noise5),
            &mut status,
        );

        assert_eq!(status, 0);
        assert_eq!(ngoodpix, 56); /* 64 - 8 nulls */
    }

    /*
     * Test fits_img_stats_int
     */
    #[test]
    fn test_img_stats_int_basic() {
        let mut array = [0i32; 64];
        for i in 0..64 {
            array[i] = (i * 1000 + (i % 7)) as i32;
        }

        let mut status = 0;
        let mut ngoodpix: c_long = 0;
        let mut minvalue: c_int = 0;
        let mut maxvalue: c_int = 0;
        let mut mean = 0.0;
        let mut sigma = 0.0;
        let mut noise1 = 0.0;
        let mut noise2 = 0.0;
        let mut noise3 = 0.0;
        let mut noise5 = 0.0;

        fits_img_stats_int_safe(
            &array,
            8,
            8,
            false,
            0,
            Some(&mut ngoodpix),
            Some(&mut minvalue),
            Some(&mut maxvalue),
            Some(&mut mean),
            Some(&mut sigma),
            Some(&mut noise1),
            Some(&mut noise2),
            Some(&mut noise3),
            Some(&mut noise5),
            &mut status,
        );

        assert_eq!(status, 0);
        assert_eq!(ngoodpix, 64);
        assert_eq!(minvalue, 0);
        assert_eq!(maxvalue, 63000); /* 63 * 1000 + 0 */
    }

    /*
     * Test fits_img_stats_int with null values
     */
    #[test]
    fn test_img_stats_int_nulls() {
        let mut array = [0i32; 64];
        let null_value: c_int = -99999;
        for i in 0..64 {
            if i % 10 == 0 {
                array[i] = null_value;
            } else {
                array[i] = (i * 1000 + (i % 7)) as i32;
            }
        }

        let mut status = 0;
        let mut ngoodpix: c_long = 0;
        let mut minvalue: c_int = 0;
        let mut maxvalue: c_int = 0;
        let mut mean = 0.0;
        let mut sigma = 0.0;
        let mut noise1 = 0.0;
        let mut noise2 = 0.0;
        let mut noise3 = 0.0;
        let mut noise5 = 0.0;

        fits_img_stats_int_safe(
            &array,
            8,
            8,
            true,
            null_value,
            Some(&mut ngoodpix),
            Some(&mut minvalue),
            Some(&mut maxvalue),
            Some(&mut mean),
            Some(&mut sigma),
            Some(&mut noise1),
            Some(&mut noise2),
            Some(&mut noise3),
            Some(&mut noise5),
            &mut status,
        );

        assert_eq!(status, 0);
        assert_eq!(ngoodpix, 57); /* 64 - 7 nulls (0, 10, 20, 30, 40, 50, 60) */
    }

    /*
     * Test fits_img_stats_float
     */
    #[test]
    fn test_img_stats_float_basic() {
        let mut array = [0.0f32; 64];
        for i in 0..64 {
            array[i] = (i as f64 * 10.5 + (i % 5) as f64 * 0.3) as f32;
        }

        let mut status = 0;
        let mut ngoodpix: c_long = 0;
        let mut minvalue: f32 = 0.0;
        let mut maxvalue: f32 = 0.0;
        let mut mean = 0.0;
        let mut sigma = 0.0;
        let mut noise1 = 0.0;
        let mut noise2 = 0.0;
        let mut noise3 = 0.0;
        let mut noise5 = 0.0;

        fits_img_stats_float_safe(
            &array,
            8,
            8,
            false,
            0.0,
            Some(&mut ngoodpix),
            Some(&mut minvalue),
            Some(&mut maxvalue),
            Some(&mut mean),
            Some(&mut sigma),
            Some(&mut noise1),
            Some(&mut noise2),
            Some(&mut noise3),
            Some(&mut noise5),
            &mut status,
        );

        assert_eq!(status, 0);
        assert_eq!(ngoodpix, 64);
    }

    /*
     * Test fits_img_stats_float with null values
     */
    #[test]
    fn test_img_stats_float_nulls() {
        let mut array = [0.0f32; 64];
        let null_value = -999.0f32;
        for i in 0..64 {
            if i % 8 == 0 {
                array[i] = null_value;
            } else {
                array[i] = (i as f64 * 10.5 + (i % 5) as f64 * 0.3) as f32;
            }
        }

        let mut status = 0;
        let mut ngoodpix: c_long = 0;
        let mut minvalue: f32 = 0.0;
        let mut maxvalue: f32 = 0.0;
        let mut mean = 0.0;
        let mut sigma = 0.0;
        let mut noise1 = 0.0;
        let mut noise2 = 0.0;
        let mut noise3 = 0.0;
        let mut noise5 = 0.0;

        fits_img_stats_float_safe(
            &array,
            8,
            8,
            true,
            null_value,
            Some(&mut ngoodpix),
            Some(&mut minvalue),
            Some(&mut maxvalue),
            Some(&mut mean),
            Some(&mut sigma),
            Some(&mut noise1),
            Some(&mut noise2),
            Some(&mut noise3),
            Some(&mut noise5),
            &mut status,
        );

        assert_eq!(status, 0);
        assert_eq!(ngoodpix, 56); /* 64 - 8 nulls */
    }

    /*
     * Test quantization with positive qlevel (noise-based)
     */
    #[test]
    fn test_quantize_positive_qlevel() {
        let mut fdata = [0.0f32; 256];
        for i in 0..256 {
            fdata[i] = (i as f64 * 10.0 + (i % 7) as f64 * 0.5) as f32;
        }

        let mut bscale = 0.0;
        let mut bzero = 0.0;
        let mut iminval = 0;
        let mut imaxval = 0;

        let result = fits_quantize_float_inplace(
            0,
            &mut fdata,
            16,
            16,
            false,
            0.0,
            4.0,
            dither(NO_DITHER),
            &mut bscale,
            &mut bzero,
            &mut iminval,
            &mut imaxval,
        );

        /* May or may not quantize depending on noise estimate */
        assert!(result >= 0); /* Should not fail */
    }

    /*
     * Test double quantization with positive qlevel
     */
    #[test]
    fn test_quantize_double_positive_qlevel() {
        let mut ddata = [0.0f64; 256];
        for i in 0..256 {
            ddata[i] = i as f64 * 10.0 + (i % 7) as f64 * 0.5;
        }

        let mut bscale = 0.0;
        let mut bzero = 0.0;
        let mut iminval = 0;
        let mut imaxval = 0;

        let result = fits_quantize_double_inplace(
            0,
            &mut ddata,
            16,
            16,
            false,
            0.0,
            4.0,
            dither(NO_DITHER),
            &mut bscale,
            &mut bzero,
            &mut iminval,
            &mut imaxval,
        );

        assert!(result >= 0); /* Should not fail */
    }

    /*
     * Test double quantization with dithering
     */
    #[test]
    fn test_quantize_double_dither1() {
        let mut ddata = [0.0f64; 64];
        for i in 0..64 {
            ddata[i] = i as f64 * 5.0 + (i % 5) as f64 * 0.3;
        }

        let mut bscale = 0.0;
        let mut bzero = 0.0;
        let mut iminval = 0;
        let mut imaxval = 0;

        let result = fits_quantize_double_inplace(
            1,
            &mut ddata,
            8,
            8,
            false,
            0.0,
            -1.0,
            dither(SUBTRACTIVE_DITHER_1),
            &mut bscale,
            &mut bzero,
            &mut iminval,
            &mut imaxval,
        );

        assert_eq!(result, 1);
        assert_ne!(bscale, 0.0);
    }

    /*
     * Test double quantization with dithering method 2
     */
    #[test]
    fn test_quantize_double_dither2() {
        let mut ddata = [0.0f64; 64];
        for i in 0..64 {
            if i % 8 == 0 {
                ddata[i] = 0.0;
            } else {
                ddata[i] = i as f64 * 5.0 + (i % 5) as f64 * 0.3;
            }
        }

        let mut bscale = 0.0;
        let mut bzero = 0.0;
        let mut iminval = 0;
        let mut imaxval = 0;

        let result = fits_quantize_double_inplace(
            1,
            &mut ddata,
            8,
            8,
            false,
            0.0,
            -1.0,
            dither(SUBTRACTIVE_DITHER_2),
            &mut bscale,
            &mut bzero,
            &mut iminval,
            &mut imaxval,
        );

        assert_eq!(result, 1);
        assert_ne!(bscale, 0.0);
    }

    /*
     * Test quantization with range exceeding int
     */
    #[test]
    fn test_quantize_range_too_large() {
        let mut fdata = [0.0f32; 64];
        for i in 0..64 {
            if i < 32 {
                fdata[i] = -1e38;
            } else {
                fdata[i] = 1e38;
            }
        }

        let mut bscale = 0.0;
        let mut bzero = 0.0;
        let mut iminval = 0;
        let mut imaxval = 0;

        /* Should return 0 because range is too large */
        let result = fits_quantize_float_inplace(
            0,
            &mut fdata,
            8,
            8,
            false,
            0.0,
            -0.001,
            dither(NO_DITHER),
            &mut bscale,
            &mut bzero,
            &mut iminval,
            &mut imaxval,
        );

        assert_eq!(result, 0); /* Should not quantize */
    }

    /*
     * Test double quantization small image
     */
    #[test]
    fn test_quantize_double_small() {
        let mut ddata = [42.0f64; 1];

        let mut bscale = 0.0;
        let mut bzero = 0.0;
        let mut iminval = 0;
        let mut imaxval = 0;

        let result = fits_quantize_double_inplace(
            0,
            &mut ddata,
            1,
            1,
            false,
            0.0,
            -1.0,
            dither(NO_DITHER),
            &mut bscale,
            &mut bzero,
            &mut iminval,
            &mut imaxval,
        );

        assert_eq!(result, 0);
        assert_eq!(bscale, 1.0);
        assert_eq!(bzero, 0.0);
    }

    /*
     * Test img_stats with larger image
     */
    #[test]
    fn test_img_stats_short_large() {
        let mut array = vec![0i16; 256 * 256];
        for i in 0..256 * 256 {
            let x = i % 256;
            let y = i / 256;
            array[i] = ((x + y * 10) % 32000) as i16;
        }

        let mut status = 0;
        let mut ngoodpix: c_long = 0;
        let mut minvalue: c_short = 0;
        let mut maxvalue: c_short = 0;
        let mut mean = 0.0;
        let mut sigma = 0.0;
        let mut noise1 = 0.0;
        let mut noise2 = 0.0;
        let mut noise3 = 0.0;
        let mut noise5 = 0.0;

        fits_img_stats_short_safe(
            &array,
            256,
            256,
            false,
            0,
            Some(&mut ngoodpix),
            Some(&mut minvalue),
            Some(&mut maxvalue),
            Some(&mut mean),
            Some(&mut sigma),
            Some(&mut noise1),
            Some(&mut noise2),
            Some(&mut noise3),
            Some(&mut noise5),
            &mut status,
        );

        assert_eq!(status, 0);
        assert_eq!(ngoodpix, 256 * 256);
    }

    /*
     * Test img_stats_int with larger image
     */
    #[test]
    fn test_img_stats_int_large() {
        let mut array = vec![0i32; 256 * 256];
        for i in 0..256 * 256 {
            let x = i % 256;
            let y = i / 256;
            array[i] = (x + y * 1000) as i32;
        }

        let mut status = 0;
        let mut ngoodpix: c_long = 0;
        let mut minvalue: c_int = 0;
        let mut maxvalue: c_int = 0;
        let mut mean = 0.0;
        let mut sigma = 0.0;
        let mut noise1 = 0.0;
        let mut noise2 = 0.0;
        let mut noise3 = 0.0;
        let mut noise5 = 0.0;

        fits_img_stats_int_safe(
            &array,
            256,
            256,
            false,
            0,
            Some(&mut ngoodpix),
            Some(&mut minvalue),
            Some(&mut maxvalue),
            Some(&mut mean),
            Some(&mut sigma),
            Some(&mut noise1),
            Some(&mut noise2),
            Some(&mut noise3),
            Some(&mut noise5),
            &mut status,
        );

        assert_eq!(status, 0);
        assert_eq!(ngoodpix, 256 * 256);
    }

    /*
     * Test img_stats_float with larger image
     */
    #[test]
    fn test_img_stats_float_large() {
        let mut array = vec![0.0f32; 256 * 256];
        for i in 0..256 * 256 {
            let x = i % 256;
            let y = i / 256;
            array[i] = (x as f64 + y as f64 * 10.5 + (i % 7) as f64 * 0.1) as f32;
        }

        let mut status = 0;
        let mut ngoodpix: c_long = 0;
        let mut minvalue: f32 = 0.0;
        let mut maxvalue: f32 = 0.0;
        let mut mean = 0.0;
        let mut sigma = 0.0;
        let mut noise1 = 0.0;
        let mut noise2 = 0.0;
        let mut noise3 = 0.0;
        let mut noise5 = 0.0;

        fits_img_stats_float_safe(
            &array,
            256,
            256,
            false,
            0.0,
            Some(&mut ngoodpix),
            Some(&mut minvalue),
            Some(&mut maxvalue),
            Some(&mut mean),
            Some(&mut sigma),
            Some(&mut noise1),
            Some(&mut noise2),
            Some(&mut noise3),
            Some(&mut noise5),
            &mut status,
        );

        assert_eq!(status, 0);
        assert_eq!(ngoodpix, 256 * 256);
    }

    /*
     * Test img_stats with very small image (nx < 9)
     */
    #[test]
    fn test_img_stats_short_small() {
        let array: [i16; 4] = [10, 20, 30, 40];

        let mut status = 0;
        let mut ngoodpix: c_long = 0;
        let mut minvalue: c_short = 0;
        let mut maxvalue: c_short = 0;
        let mut mean = 0.0;
        let mut sigma = 0.0;
        let mut noise1 = 0.0;
        let mut noise2 = 0.0;
        let mut noise3 = 0.0;
        let mut noise5 = 0.0;

        fits_img_stats_short_safe(
            &array,
            4,
            1,
            false,
            0,
            Some(&mut ngoodpix),
            Some(&mut minvalue),
            Some(&mut maxvalue),
            Some(&mut mean),
            Some(&mut sigma),
            Some(&mut noise1),
            Some(&mut noise2),
            Some(&mut noise3),
            Some(&mut noise5),
            &mut status,
        );

        assert_eq!(status, 0);
        assert_eq!(ngoodpix, 4);
        assert_eq!(minvalue, 10);
        assert_eq!(maxvalue, 40);
    }

    /*
     * Test img_stats with single pixel (ngood == 1 branch)
     */
    #[test]
    fn test_img_stats_short_single() {
        let array: [i16; 1] = [42];

        let mut status = 0;
        let mut ngoodpix: c_long = 0;
        let mut minvalue: c_short = 0;
        let mut maxvalue: c_short = 0;
        let mut mean = 0.0;
        let mut sigma = 0.0;
        let mut noise1 = 0.0;
        let mut noise2 = 0.0;
        let mut noise3 = 0.0;
        let mut noise5 = 0.0;

        fits_img_stats_short_safe(
            &array,
            1,
            1,
            false,
            0,
            Some(&mut ngoodpix),
            Some(&mut minvalue),
            Some(&mut maxvalue),
            Some(&mut mean),
            Some(&mut sigma),
            Some(&mut noise1),
            Some(&mut noise2),
            Some(&mut noise3),
            Some(&mut noise5),
            &mut status,
        );

        assert_eq!(status, 0);
        assert_eq!(ngoodpix, 1);
        assert_eq!(mean, 42.0);
        assert_eq!(sigma, 0.0);
    }

    /*
     * Test img_stats_int with single pixel
     */
    #[test]
    fn test_img_stats_int_single() {
        let array: [i32; 1] = [42];

        let mut status = 0;
        let mut ngoodpix: c_long = 0;
        let mut minvalue: c_int = 0;
        let mut maxvalue: c_int = 0;
        let mut mean = 0.0;
        let mut sigma = 0.0;
        let mut noise1 = 0.0;
        let mut noise2 = 0.0;
        let mut noise3 = 0.0;
        let mut noise5 = 0.0;

        fits_img_stats_int_safe(
            &array,
            1,
            1,
            false,
            0,
            Some(&mut ngoodpix),
            Some(&mut minvalue),
            Some(&mut maxvalue),
            Some(&mut mean),
            Some(&mut sigma),
            Some(&mut noise1),
            Some(&mut noise2),
            Some(&mut noise3),
            Some(&mut noise5),
            &mut status,
        );

        assert_eq!(status, 0);
        assert_eq!(ngoodpix, 1);
        assert_eq!(mean, 42.0);
        assert_eq!(sigma, 0.0);
    }

    /*
     * Test img_stats_float with single pixel
     */
    #[test]
    fn test_img_stats_float_single() {
        let array: [f32; 1] = [42.5];

        let mut status = 0;
        let mut ngoodpix: c_long = 0;
        let mut minvalue: f32 = 0.0;
        let mut maxvalue: f32 = 0.0;
        let mut mean = 0.0;
        let mut sigma = 0.0;
        let mut noise1 = 0.0;
        let mut noise2 = 0.0;
        let mut noise3 = 0.0;
        let mut noise5 = 0.0;

        fits_img_stats_float_safe(
            &array,
            1,
            1,
            false,
            0.0,
            Some(&mut ngoodpix),
            Some(&mut minvalue),
            Some(&mut maxvalue),
            Some(&mut mean),
            Some(&mut sigma),
            Some(&mut noise1),
            Some(&mut noise2),
            Some(&mut noise3),
            Some(&mut noise5),
            &mut status,
        );

        assert_eq!(status, 0);
        assert_eq!(ngoodpix, 1);
        assert!((mean - 42.5).abs() <= 0.001);
        assert_eq!(sigma, 0.0);
    }

    /*
     * Test img_stats with all nulls (ngood == 0 branch)
     */
    #[test]
    fn test_img_stats_short_all_nulls() {
        let array: [i16; 4] = [-999, -999, -999, -999];

        let mut status = 0;
        let mut ngoodpix: c_long = 0;
        let mut minvalue: c_short = 0;
        let mut maxvalue: c_short = 0;
        let mut mean = 0.0;
        let mut sigma = 0.0;
        let mut noise1 = 0.0;
        let mut noise2 = 0.0;
        let mut noise3 = 0.0;
        let mut noise5 = 0.0;

        fits_img_stats_short_safe(
            &array,
            4,
            1,
            true,
            -999,
            Some(&mut ngoodpix),
            Some(&mut minvalue),
            Some(&mut maxvalue),
            Some(&mut mean),
            Some(&mut sigma),
            Some(&mut noise1),
            Some(&mut noise2),
            Some(&mut noise3),
            Some(&mut noise5),
            &mut status,
        );

        assert_eq!(status, 0);
        assert_eq!(ngoodpix, 0);
        assert_eq!(mean, 0.0);
        assert_eq!(sigma, 0.0);
    }

    /*
     * Test img_stats_int all nulls
     */
    #[test]
    fn test_img_stats_int_all_nulls() {
        let array: [i32; 4] = [-999, -999, -999, -999];

        let mut status = 0;
        let mut ngoodpix: c_long = 0;
        let mut minvalue: c_int = 0;
        let mut maxvalue: c_int = 0;
        let mut mean = 0.0;
        let mut sigma = 0.0;
        let mut noise1 = 0.0;
        let mut noise2 = 0.0;
        let mut noise3 = 0.0;
        let mut noise5 = 0.0;

        fits_img_stats_int_safe(
            &array,
            4,
            1,
            true,
            -999,
            Some(&mut ngoodpix),
            Some(&mut minvalue),
            Some(&mut maxvalue),
            Some(&mut mean),
            Some(&mut sigma),
            Some(&mut noise1),
            Some(&mut noise2),
            Some(&mut noise3),
            Some(&mut noise5),
            &mut status,
        );

        assert_eq!(status, 0);
        assert_eq!(ngoodpix, 0);
        assert_eq!(mean, 0.0);
        assert_eq!(sigma, 0.0);
    }

    /*
     * Test img_stats_float all nulls
     */
    #[test]
    fn test_img_stats_float_all_nulls() {
        let array: [f32; 4] = [-999.0, -999.0, -999.0, -999.0];

        let mut status = 0;
        let mut ngoodpix: c_long = 0;
        let mut minvalue: f32 = 0.0;
        let mut maxvalue: f32 = 0.0;
        let mut mean = 0.0;
        let mut sigma = 0.0;
        let mut noise1 = 0.0;
        let mut noise2 = 0.0;
        let mut noise3 = 0.0;
        let mut noise5 = 0.0;

        fits_img_stats_float_safe(
            &array,
            4,
            1,
            true,
            -999.0,
            Some(&mut ngoodpix),
            Some(&mut minvalue),
            Some(&mut maxvalue),
            Some(&mut mean),
            Some(&mut sigma),
            Some(&mut noise1),
            Some(&mut noise2),
            Some(&mut noise3),
            Some(&mut noise5),
            &mut status,
        );

        assert_eq!(status, 0);
        assert_eq!(ngoodpix, 0);
        assert_eq!(mean, 0.0);
        assert_eq!(sigma, 0.0);
    }

    /*
     * Test quantization with qlevel == 0 (default quantization)
     */
    #[test]
    fn test_quantize_qlevel_zero() {
        let mut fdata = [0.0f32; 256];
        for i in 0..256 {
            fdata[i] = (i as f64 * 10.0 + (i % 17) as f64 * 2.5) as f32;
        }

        let mut bscale = 0.0;
        let mut bzero = 0.0;
        let mut iminval = 0;
        let mut imaxval = 0;

        let result = fits_quantize_float_inplace(
            0,
            &mut fdata,
            16,
            16,
            false,
            0.0,
            0.0,
            dither(NO_DITHER),
            &mut bscale,
            &mut bzero,
            &mut iminval,
            &mut imaxval,
        );

        assert!(result >= 0);
    }

    /*
     * Test double quantization with qlevel == 0
     */
    #[test]
    fn test_quantize_double_qlevel_zero() {
        let mut ddata = [0.0f64; 256];
        for i in 0..256 {
            ddata[i] = i as f64 * 10.0 + (i % 17) as f64 * 2.5;
        }

        let mut bscale = 0.0;
        let mut bzero = 0.0;
        let mut iminval = 0;
        let mut imaxval = 0;

        let result = fits_quantize_double_inplace(
            0,
            &mut ddata,
            16,
            16,
            false,
            0.0,
            0.0,
            dither(NO_DITHER),
            &mut bscale,
            &mut bzero,
            &mut iminval,
            &mut imaxval,
        );

        assert!(result >= 0);
    }

    /*
     * Test float quantization all null values (ngood == 0)
     */
    #[test]
    fn test_quantize_float_all_nulls() {
        let null_value = -999.0f32;
        let mut fdata = [null_value; 64];

        let mut bscale = 0.0;
        let mut bzero = 0.0;
        let mut iminval = 0;
        let mut imaxval = 0;

        let result = fits_quantize_float_inplace(
            0,
            &mut fdata,
            8,
            8,
            true,
            null_value,
            4.0,
            dither(NO_DITHER),
            &mut bscale,
            &mut bzero,
            &mut iminval,
            &mut imaxval,
        );

        /* Should return 1 with dummy values for all nulls */
        assert_eq!(result, 1);
    }

    /*
     * Test double quantization all null values (ngood == 0)
     */
    #[test]
    fn test_quantize_double_all_nulls() {
        let null_value = -999.0f64;
        let mut ddata = [null_value; 64];

        let mut bscale = 0.0;
        let mut bzero = 0.0;
        let mut iminval = 0;
        let mut imaxval = 0;

        let result = fits_quantize_double_inplace(
            0,
            &mut ddata,
            8,
            8,
            true,
            null_value,
            4.0,
            dither(NO_DITHER),
            &mut bscale,
            &mut bzero,
            &mut iminval,
            &mut imaxval,
        );

        /* Should return 1 with dummy values for all nulls */
        assert_eq!(result, 1);
    }

    /*
     * Test float dithering with null values (row > 0 + nulls)
     */
    #[test]
    fn test_quantize_float_dither_nulls() {
        let mut fdata = [0.0f32; 64];
        let null_value = -999.0f32;
        for i in 0..64 {
            if i % 8 == 0 {
                fdata[i] = null_value;
            } else {
                fdata[i] = (i as f64 * 10.0 + (i % 3) as f64 * 0.5) as f32;
            }
        }

        let mut bscale = 0.0;
        let mut bzero = 0.0;
        let mut iminval = 0;
        let mut imaxval = 0;

        let result = fits_quantize_float_inplace(
            1,
            &mut fdata,
            8,
            8,
            true,
            null_value,
            -1.0,
            dither(SUBTRACTIVE_DITHER_1),
            &mut bscale,
            &mut bzero,
            &mut iminval,
            &mut imaxval,
        );

        assert_eq!(result, 1);

        for i in 0..64 {
            if i % 8 == 0 {
                assert_eq!(idata_f(&fdata, i), NULL_VALUE);
            }
        }
    }

    /*
     * Test double quantization with nulls and no dithering
     */
    #[test]
    fn test_quantize_double_nulls_nodither() {
        let mut ddata = [0.0f64; 64];
        let null_value = -999.0f64;
        for i in 0..64 {
            if i % 8 == 0 {
                ddata[i] = null_value;
            } else {
                ddata[i] = i as f64 * 10.0 + (i % 3) as f64 * 0.5;
            }
        }

        let mut bscale = 0.0;
        let mut bzero = 0.0;
        let mut iminval = 0;
        let mut imaxval = 0;

        let result = fits_quantize_double_inplace(
            0,
            &mut ddata,
            8,
            8,
            true,
            null_value,
            -1.0,
            dither(NO_DITHER),
            &mut bscale,
            &mut bzero,
            &mut iminval,
            &mut imaxval,
        );

        assert_eq!(result, 1);

        for i in 0..64 {
            if i % 8 == 0 {
                assert_eq!(idata_d(&ddata, i), NULL_VALUE);
            }
        }
    }

    /*
     * Test double range too large
     */
    #[test]
    fn test_quantize_double_range_too_large() {
        let mut ddata = [0.0f64; 64];
        for i in 0..64 {
            if i < 32 {
                ddata[i] = -1e300;
            } else {
                ddata[i] = 1e300;
            }
        }

        let mut bscale = 0.0;
        let mut bzero = 0.0;
        let mut iminval = 0;
        let mut imaxval = 0;

        let result = fits_quantize_double_inplace(
            0,
            &mut ddata,
            8,
            8,
            false,
            0.0,
            -0.001,
            dither(NO_DITHER),
            &mut bscale,
            &mut bzero,
            &mut iminval,
            &mut imaxval,
        );

        assert_eq!(result, 0); /* Should not quantize */
    }

    /*
     * Test float with very large range (zeropt centering)
     */
    #[test]
    fn test_quantize_float_large_range() {
        let mut fdata = [0.0f32; 64];
        for i in 0..64 {
            fdata[i] = ((i as f64 - 32.0) * 1e8) as f32;
        }

        let mut bscale = 0.0;
        let mut bzero = 0.0;
        let mut iminval = 0;
        let mut imaxval = 0;

        let result = fits_quantize_float_inplace(
            0,
            &mut fdata,
            8,
            8,
            false,
            0.0,
            -1e6,
            dither(NO_DITHER),
            &mut bscale,
            &mut bzero,
            &mut iminval,
            &mut imaxval,
        );

        assert!(result >= 0);
    }

    /*
     * Test img_stats short with small and null rows
     */
    #[test]
    fn test_img_stats_short_small_nulls() {
        let array: [i16; 8] = [-999, 10, 20, 30, 40, 50, 60, -999];

        let mut status = 0;
        let mut ngoodpix: c_long = 0;
        let mut minvalue: c_short = 0;
        let mut maxvalue: c_short = 0;
        let mut mean = 0.0;
        let mut sigma = 0.0;
        let mut noise1 = 0.0;
        let mut noise2 = 0.0;
        let mut noise3 = 0.0;
        let mut noise5 = 0.0;

        fits_img_stats_short_safe(
            &array,
            8,
            1,
            true,
            -999,
            Some(&mut ngoodpix),
            Some(&mut minvalue),
            Some(&mut maxvalue),
            Some(&mut mean),
            Some(&mut sigma),
            Some(&mut noise1),
            Some(&mut noise2),
            Some(&mut noise3),
            Some(&mut noise5),
            &mut status,
        );

        assert_eq!(status, 0);
        assert_eq!(ngoodpix, 6); /* 8 - 2 nulls */
    }

    /*
     * Test float dithering with null values and zeros (DITHER_2)
     */
    #[test]
    fn test_quantize_float_dither2_nulls_zeros() {
        let mut fdata = [0.0f32; 64];
        let null_value = -999.0f32;
        for i in 0..64 {
            if i % 8 == 0 {
                fdata[i] = null_value;
            } else if i % 4 == 0 {
                fdata[i] = 0.0; /* Zero value in DITHER_2 */
            } else {
                fdata[i] = (i as f64 * 10.0 + (i % 3) as f64 * 0.5) as f32;
            }
        }

        let mut bscale = 0.0;
        let mut bzero = 0.0;
        let mut iminval = 0;
        let mut imaxval = 0;

        let result = fits_quantize_float_inplace(
            1,
            &mut fdata,
            8,
            8,
            true,
            null_value,
            -1.0,
            dither(SUBTRACTIVE_DITHER_2),
            &mut bscale,
            &mut bzero,
            &mut iminval,
            &mut imaxval,
        );

        assert_eq!(result, 1);

        for i in 0..64 {
            if i % 8 == 0 {
                assert_eq!(idata_f(&fdata, i), NULL_VALUE);
            }
        }
    }

    /*
     * Test double dithering with null values and zeros (DITHER_2)
     */
    #[test]
    fn test_quantize_double_dither2_nulls_zeros() {
        let mut ddata = [0.0f64; 64];
        let null_value = -999.0f64;
        for i in 0..64 {
            if i % 8 == 0 {
                ddata[i] = null_value;
            } else if i % 4 == 0 {
                ddata[i] = 0.0; /* Zero value in DITHER_2 */
            } else {
                ddata[i] = i as f64 * 10.0 + (i % 3) as f64 * 0.5;
            }
        }

        let mut bscale = 0.0;
        let mut bzero = 0.0;
        let mut iminval = 0;
        let mut imaxval = 0;

        let result = fits_quantize_double_inplace(
            1,
            &mut ddata,
            8,
            8,
            true,
            null_value,
            -1.0,
            dither(SUBTRACTIVE_DITHER_2),
            &mut bscale,
            &mut bzero,
            &mut iminval,
            &mut imaxval,
        );

        assert_eq!(result, 1);

        for i in 0..64 {
            if i % 8 == 0 {
                assert_eq!(idata_d(&ddata, i), NULL_VALUE);
            }
        }
    }

    /*
     * Test large dithered image (to exercise N_RANDOM wraparound)
     * N_RANDOM is 10000, so we need >10000 pixels per row to trigger wraparound
     */
    #[test]
    fn test_quantize_float_large_dither() {
        let nx = 10500usize; /* > N_RANDOM (10000) to trigger wraparound */
        let ny = 10usize;
        let size = nx * ny;

        let mut fdata = vec![0.0f32; size];
        for i in 0..size {
            let x = i % nx;
            let y = i / nx;
            fdata[i] = (x as f64 + y as f64 * 100.0 + (i % 7) as f64 * 0.1) as f32;
        }

        let mut bscale = 0.0;
        let mut bzero = 0.0;
        let mut iminval = 0;
        let mut imaxval = 0;

        let result = fits_quantize_float_inplace(
            1,
            &mut fdata,
            nx,
            ny,
            false,
            0.0,
            -1.0,
            dither(SUBTRACTIVE_DITHER_1),
            &mut bscale,
            &mut bzero,
            &mut iminval,
            &mut imaxval,
        );

        assert_eq!(result, 1);
        assert_ne!(bscale, 0.0);
    }

    /*
     * Test large dithered double image
     */
    #[test]
    fn test_quantize_double_large_dither() {
        let nx = 10500usize;
        let ny = 10usize;
        let size = nx * ny;

        let mut ddata = vec![0.0f64; size];
        for i in 0..size {
            let x = i % nx;
            let y = i / nx;
            ddata[i] = x as f64 + y as f64 * 100.0 + (i % 7) as f64 * 0.1;
        }

        let mut bscale = 0.0;
        let mut bzero = 0.0;
        let mut iminval = 0;
        let mut imaxval = 0;

        let result = fits_quantize_double_inplace(
            1,
            &mut ddata,
            nx,
            ny,
            false,
            0.0,
            -1.0,
            dither(SUBTRACTIVE_DITHER_1),
            &mut bscale,
            &mut bzero,
            &mut iminval,
            &mut imaxval,
        );

        assert_eq!(result, 1);
        assert_ne!(bscale, 0.0);
    }

    /*
     * Test large dithered image with null values (to exercise N_RANDOM wraparound)
     */
    #[test]
    fn test_quantize_float_large_dither_nulls() {
        let nx = 10500usize;
        let ny = 10usize;
        let size = nx * ny;
        let null_value = -999.0f32;

        let mut fdata = vec![0.0f32; size];
        for i in 0..size {
            if i % 1000 == 0 {
                fdata[i] = null_value;
            } else {
                fdata[i] =
                    ((i % nx) as f64 + (i / nx) as f64 * 100.0 + (i % 7) as f64 * 0.1) as f32;
            }
        }

        let mut bscale = 0.0;
        let mut bzero = 0.0;
        let mut iminval = 0;
        let mut imaxval = 0;

        let result = fits_quantize_float_inplace(
            1,
            &mut fdata,
            nx,
            ny,
            true,
            null_value,
            -1.0,
            dither(SUBTRACTIVE_DITHER_1),
            &mut bscale,
            &mut bzero,
            &mut iminval,
            &mut imaxval,
        );

        assert_eq!(result, 1);
        assert_ne!(bscale, 0.0);
    }

    /*
     * Test large dithered double image with null values
     */
    #[test]
    fn test_quantize_double_large_dither_nulls() {
        let nx = 10500usize;
        let ny = 10usize;
        let size = nx * ny;
        let null_value = -999.0f64;

        let mut ddata = vec![0.0f64; size];
        for i in 0..size {
            if i % 1000 == 0 {
                ddata[i] = null_value;
            } else {
                ddata[i] = (i % nx) as f64 + (i / nx) as f64 * 100.0 + (i % 7) as f64 * 0.1;
            }
        }

        let mut bscale = 0.0;
        let mut bzero = 0.0;
        let mut iminval = 0;
        let mut imaxval = 0;

        let result = fits_quantize_double_inplace(
            1,
            &mut ddata,
            nx,
            ny,
            true,
            null_value,
            -1.0,
            dither(SUBTRACTIVE_DITHER_1),
            &mut bscale,
            &mut bzero,
            &mut iminval,
            &mut imaxval,
        );

        assert_eq!(result, 1);
        assert_ne!(bscale, 0.0);
    }

    /*
     * Test float with range requiring centering
     */
    #[test]
    fn test_quantize_float_needs_centering() {
        let mut fdata = [0.0f32; 64];
        for i in 0..64 {
            fdata[i] = (i as f32 - 32.0) * 5e7;
        }

        let mut bscale = 0.0;
        let mut bzero = 0.0;
        let mut iminval = 0;
        let mut imaxval = 0;

        let result = fits_quantize_float_inplace(
            0,
            &mut fdata,
            8,
            8,
            false,
            0.0,
            -1.0,
            dither(NO_DITHER),
            &mut bscale,
            &mut bzero,
            &mut iminval,
            &mut imaxval,
        );

        assert!(result >= 0);
    }

    /*
     * Test double with range requiring centering
     */
    #[test]
    fn test_quantize_double_needs_centering() {
        let mut ddata = [0.0f64; 64];
        for i in 0..64 {
            ddata[i] = (i as f64 - 32.0) * 5e7;
        }

        let mut bscale = 0.0;
        let mut bzero = 0.0;
        let mut iminval = 0;
        let mut imaxval = 0;

        let result = fits_quantize_double_inplace(
            0,
            &mut ddata,
            8,
            8,
            false,
            0.0,
            -1.0,
            dither(NO_DITHER),
            &mut bscale,
            &mut bzero,
            &mut iminval,
            &mut imaxval,
        );

        assert!(result >= 0);
    }
}
