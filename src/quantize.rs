/*
  The following code is based on algorithms written by Richard White at STScI and made
  available for use in CFITSIO in July 1999 and updated in January 2008.

*/

use core::slice;

use crate::c_types::*;

use bytemuck::{cast, cast_slice_mut};

use crate::{
    fitsio::{LONGLONG, MEMORY_ALLOCATION},
    fitsio2::N_RANDOM,
    imcompress::{FITS_RAND_VALUE, fits_init_randoms},
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum DitherType {
    NoDither = -1,
    SubtractiveDither1 = 1,
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
const SIGMA_CLIP: f64 = 5.0;

const NITER: i32 = 3; /* number of sigma-clipping iterations */

pub trait SomeSet<T> {
    fn is_some_set(self, v: T);
}

impl<T> SomeSet<T> for Option<&mut T> {
    fn is_some_set(self, v: T) {
        let s = self;
        if let Some(s) = s {
            *s = v;
        }
    }
}

/*---------------------------------------------------------------------------*/
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
            delta = stdev / qlevel as f64;
        }
        if delta == 0.0 {
            return 0; /* don't quantize */
        }
    } else {
        /* negative value represents the absolute quantization level */
        delta = -qlevel as f64;

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
    if (maxval - minval) as f64 / delta > 2.0 * 2147483647.0 - N_RESERVED_VALUES as f64 {
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
            zeropt = (minval as f64) - delta * (NULL_VALUE + N_RESERVED_VALUES) as f64;
        } else if (maxval - minval) as f64 / delta < 2147483647.0 - N_RESERVED_VALUES as f64 {
            zeropt = minval.into();
            /* fudge the zero point so it is an integer multiple of delta */
            /* This helps to ensure the same scaling will be performed if the */
            /* file undergoes multiple fpack/funpack cycles */
            iqfactor = (zeropt / delta + 0.5) as i64;
            zeropt = iqfactor as f64 * delta;
        } else {
            /* center the quantized levels around zero */
            zeropt = (minval + maxval) as f64 / 2.;
        }

        if row > 0 {
            /* dither the values when quantizing */
            for i in 0..nx {
                // for (i = 0;  i < nx;  i+=1) {

                if dither_method == DitherType::SubtractiveDither2 && fdata[i] == 0.0 {
                    fdata[i] = cast(ZERO_VALUE);
                } else {
                    fdata[i] = cast(nint_f64(
                        ((fdata[i] as f64 - zeropt) / delta) + fits_rand_value[nextrand] as f64
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
                fdata[i] = cast(nint_f64((fdata[i] as f64 - zeropt) / delta));
            }
        }
    } else {
        /* data contains null values; shift the range to be */
        /* close to the value used to represent null values */
        zeropt = (minval as f64) - delta * (NULL_VALUE + N_RESERVED_VALUES) as f64;

        if row > 0 {
            /* dither the values */
            for i in 0..nx {
                // for (i = 0;  i < nx;  i+=1) {
                if fdata[i] != in_null_value {
                    if dither_method == DitherType::SubtractiveDither2 && fdata[i] == 0.0 {
                        fdata[i] = cast(ZERO_VALUE);
                    } else {
                        fdata[i] = cast(nint_f64(
                            ((fdata[i] as f64 - zeropt) / delta) + fits_rand_value[nextrand] as f64
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
                    fdata[i] = cast(nint_f64((fdata[i] as f64 - zeropt) / delta));
                } else {
                    fdata[i] = cast(NULL_VALUE);
                }
            }
        }
    }

    /* calc min and max values */
    let mut temp: f64 = (minval as f64 - zeropt) / delta;
    *iminval = nint_f64(temp);
    temp = (maxval as f64 - zeropt) / delta;
    *imaxval = nint_f64(temp);

    *bscale = delta;
    *bzero = zeropt;

    /* yes, data have been quantized */
    1
}

/*---------------------------------------------------------------------------*/
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
            delta = stdev / qlevel as f64;
        }
        if delta == 0.0 {
            return 0; /* don't quantize */
        }
    } else {
        /* negative value represents the absolute quantization level */
        delta = -qlevel as f64;

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
    if (maxval - minval) / delta > 2.0 * 2147483647.0 - N_RESERVED_VALUES as f64 {
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
            zeropt = minval - delta * (NULL_VALUE + N_RESERVED_VALUES) as f64;
        } else if (maxval - minval) / delta < 2147483647.0 - N_RESERVED_VALUES as f64 {
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
                        ((fdata_i - zeropt) / delta) + fits_rand_value[nextrand] as f64 - 0.5,
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
        zeropt = minval - delta * (NULL_VALUE + N_RESERVED_VALUES) as f64;

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
                            ((fdata_i - zeropt) / delta) + fits_rand_value[nextrand] as f64 - 0.5,
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

/*--------------------------------------------------------------------------*/
/// Compute statistics of the input short integer image.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn fits_img_stats_short(
    array: *const c_short, /*  2 dimensional array of image pixels */
    nx: c_long,            /* number of pixels in each row of the image */
    ny: c_long,            /* number of rows in the image */
    /* (if this is a 3D image, then ny should be the */
    /* product of the no. of rows times the no. of planes) */
    nullcheck: c_int,   /* check for null values, if true */
    nullvalue: c_short, /* value of null pixels, if nullcheck is true */

    /* returned parameters (if the pointer is not null)  */
    ngoodpix: *mut c_long,  /* number of non-null pixels in the image */
    minvalue: *mut c_short, /* returned minimum non-null value in the array */
    maxvalue: *mut c_short, /* returned maximum non-null value in the array */
    mean: *mut f64,         /* returned mean value of all non-null pixels */
    sigma: *mut f64,        /* returned R.M.S. value of all non-null pixels */
    noise1: *mut f64,       /* 1st order estimate of noise in image background level */
    noise2: *mut f64,       /* 2nd order estimate of noise in image background level */
    noise3: *mut f64,       /* 3rd order estimate of noise in image background level */
    noise5: *mut f64,       /* 5th order estimate of noise in image background level */
    status: *mut c_int,     /* error status */
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

/*--------------------------------------------------------------------------*/
/// Compute statistics of the input short integer image.
pub(crate) fn fits_img_stats_short_safe(
    array: &[c_short], /*  2 dimensional array of image pixels */
    nx: c_long,        /* number of pixels in each row of the image */
    ny: c_long,        /* number of rows in the image */
    /* (if this is a 3D image, then ny should be the */
    /* product of the no. of rows times the no. of planes) */
    nullcheck: bool,    /* check for null values, if true */
    nullvalue: c_short, /* value of null pixels, if nullcheck is true */

    /* returned parameters (if the pointer is not null)  */
    mut ngoodpix: Option<&mut c_long>, /* number of non-null pixels in the image */
    minvalue: Option<&mut c_short>,    /* returned minimum non-null value in the array */
    maxvalue: Option<&mut c_short>,    /* returned maximum non-null value in the array */
    mean: Option<&mut f64>,            /* returned mean value of all non-null pixels */
    sigma: Option<&mut f64>,           /* returned R.M.S. value of all non-null pixels */
    noise1: Option<&mut f64>,          /* 1st order estimate of noise in image background level */
    noise2: Option<&mut f64>,          /* 2nd order estimate of noise in image background level */
    noise3: Option<&mut f64>,          /* 3rd order estimate of noise in image background level */
    noise5: Option<&mut f64>,          /* 5th order estimate of noise in image background level */
    status: &mut c_int,                /* error status */
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

/*--------------------------------------------------------------------------*/
/// Compute statistics of the input integer image.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn fits_img_stats_int(
    array: *const c_int, /*  2 dimensional array of image pixels */
    nx: c_long,          /* number of pixels in each row of the image */
    ny: c_long,          /* number of rows in the image */
    /* (if this is a 3D image, then ny should be the */
    /* product of the no. of rows times the no. of planes) */
    nullcheck: c_int, /* check for null values, if true */
    nullvalue: c_int, /* value of null pixels, if nullcheck is true */

    /* returned parameters (if the pointer is not null)  */
    ngoodpix: *mut c_long, /* number of non-null pixels in the image */
    minvalue: *mut c_int,  /* returned minimum non-null value in the array */
    maxvalue: *mut c_int,  /* returned maximum non-null value in the array */
    mean: *mut f64,        /* returned mean value of all non-null pixels */
    sigma: *mut f64,       /* returned R.M.S. value of all non-null pixels */
    noise1: *mut f64,      /* 1st order estimate of noise in image background level */
    noise2: *mut f64,      /* 2nd order estimate of noise in image background level */
    noise3: *mut f64,      /* 3rd order estimate of noise in image background level */
    noise5: *mut f64,      /* 5th order estimate of noise in image background level */
    status: *mut c_int,    /* error status */
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

/*--------------------------------------------------------------------------*/
/// Compute statistics of the input integer image.
pub(crate) fn fits_img_stats_int_safe(
    array: &[c_int], /*  2 dimensional array of image pixels */
    nx: c_long,      /* number of pixels in each row of the image */
    ny: c_long,      /* number of rows in the image */
    /* (if this is a 3D image, then ny should be the */
    /* product of the no. of rows times the no. of planes) */
    nullcheck: bool,  /* check for null values, if true */
    nullvalue: c_int, /* value of null pixels, if nullcheck is true */

    /* returned parameters (if the pointer is not null)  */
    mut ngoodpix: Option<&mut c_long>, /* number of non-null pixels in the image */
    minvalue: Option<&mut c_int>,      /* returned minimum non-null value in the array */
    maxvalue: Option<&mut c_int>,      /* returned maximum non-null value in the array */
    mean: Option<&mut f64>,            /* returned mean value of all non-null pixels */
    sigma: Option<&mut f64>,           /* returned R.M.S. value of all non-null pixels */
    noise1: Option<&mut f64>,          /* 1st order estimate of noise in image background level */
    noise2: Option<&mut f64>,          /* 2nd order estimate of noise in image background level */
    noise3: Option<&mut f64>,          /* 3rd order estimate of noise in image background level */
    noise5: Option<&mut f64>,          /* 5th order estimate of noise in image background level */
    status: &mut c_int,                /* error status */
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

/*--------------------------------------------------------------------------*/
/// Compute statistics of the input float image.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn fits_img_stats_float(
    array: *const f32, /*  2 dimensional array of image pixels */
    nx: c_long,        /* number of pixels in each row of the image */
    ny: c_long,        /* number of rows in the image */
    /* (if this is a 3D image, then ny should be the */
    /* product of the no. of rows times the no. of planes) */
    nullcheck: c_int, /* check for null values, if true */
    nullvalue: f32,   /* value of null pixels, if nullcheck is true */

    /* returned parameters (if the pointer is not null)  */
    ngoodpix: *mut c_long, /* number of non-null pixels in the image */
    minvalue: *mut f32,    /* returned minimum non-null value in the array */
    maxvalue: *mut f32,    /* returned maximum non-null value in the array */
    mean: *mut f64,        /* returned mean value of all non-null pixels */
    sigma: *mut f64,       /* returned R.M.S. value of all non-null pixels */
    noise1: *mut f64,      /* 1st order estimate of noise in image background level */
    noise2: *mut f64,      /* 2nd order estimate of noise in image background level */
    noise3: *mut f64,      /* 3rd order estimate of noise in image background level */
    noise5: *mut f64,      /* 5th order estimate of noise in image background level */
    status: *mut c_int,    /* error status */
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

/*--------------------------------------------------------------------------*/
/// Compute statistics of the input float image.
pub(crate) fn fits_img_stats_float_safe(
    array: &[f32], /*  2 dimensional array of image pixels */
    nx: c_long,    /* number of pixels in each row of the image */
    ny: c_long,    /* number of rows in the image */
    /* (if this is a 3D image, then ny should be the */
    /* product of the no. of rows times the no. of planes) */
    nullcheck: bool, /* check for null values, if true */
    nullvalue: f32,  /* value of null pixels, if nullcheck is true */

    /* returned parameters (if the pointer is not null)  */
    mut ngoodpix: Option<&mut c_long>, /* number of non-null pixels in the image */
    minvalue: Option<&mut f32>,        /* returned minimum non-null value in the array */
    maxvalue: Option<&mut f32>,        /* returned maximum non-null value in the array */
    mean: Option<&mut f64>,            /* returned mean value of all non-null pixels */
    sigma: Option<&mut f64>,           /* returned R.M.S. value of all non-null pixels */
    noise1: Option<&mut f64>,          /* 1st order estimate of noise in image background level */
    noise2: Option<&mut f64>,          /* 2nd order estimate of noise in image background level */
    noise3: Option<&mut f64>,          /* 3rd order estimate of noise in image background level */
    noise5: Option<&mut f64>,          /* 5th order estimate of noise in image background level */
    status: &mut c_int,                /* error status */
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

/*--------------------------------------------------------------------------*/
/// Compute mean and RMS sigma of the non-null pixels in the input array.
fn FnMeanSigma_short(
    array: &[c_short],    /*  2 dimensional array of image pixels */
    npix: usize,          /* number of pixels in the image */
    nullcheck: bool,      /* check for null values, if true */
    nullvalue: c_short,   /* value of null pixels, if nullcheck is true */
    ngoodpix: &mut usize, /* number of non-null pixels in the image */
    mean: &mut f64,       /* returned mean value of all non-null pixels */
    sigma: &mut f64,      /* returned R.M.S. value of all non-null pixels */
    status: &mut c_int,   /* error status */
) -> c_int {
    let mut ngood: usize = 0;

    let mut sum: f64 = 0.0;
    let mut sum2: f64 = 0.0;
    let mut xtemp: f64;

    let mut value: usize = 0;

    if nullcheck {
        for ii in 0..npix {
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
        for ii in 0..npix {
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

/*--------------------------------------------------------------------------*/
/// Compute mean and RMS sigma of the non-null pixels in the input array.
fn FnMeanSigma_int(
    array: &[c_int],      /*  2 dimensional array of image pixels */
    npix: usize,          /* number of pixels in the image */
    nullcheck: bool,      /* check for null values, if true */
    nullvalue: c_int,     /* value of null pixels, if nullcheck is true */
    ngoodpix: &mut usize, /* number of non-null pixels in the image */
    mean: &mut f64,       /* returned mean value of all non-null pixels */
    sigma: &mut f64,      /* returned R.M.S. value of all non-null pixels */
    status: &mut c_int,   /* error status */
) -> c_int {
    let mut ngood: usize = 0;

    let mut sum: f64 = 0.0;
    let mut sum2: f64 = 0.0;
    let mut xtemp: f64;

    let mut value: usize = 0;

    if nullcheck {
        for ii in 0..npix {
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
        for ii in 0..npix {
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

/*--------------------------------------------------------------------------*/
/// Compute mean and RMS sigma of the non-null pixels in the input array.
fn FnMeanSigma_float(
    array: &[f32],        /*  2 dimensional array of image pixels */
    npix: usize,          /* number of pixels in the image */
    nullcheck: bool,      /* check for null values, if true */
    nullvalue: f32,       /* value of null pixels, if nullcheck is true */
    ngoodpix: &mut usize, /* number of non-null pixels in the image */
    mean: &mut f64,       /* returned mean value of all non-null pixels */
    sigma: &mut f64,      /* returned R.M.S. value of all non-null pixels */
    status: &mut c_int,   /* error status */
) -> c_int {
    let mut ngood: usize = 0;

    let mut sum: f64 = 0.0;
    let mut sum2: f64 = 0.0;
    let mut xtemp: f64;

    let mut value: usize = 0;

    if nullcheck {
        for ii in 0..npix {
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
        for ii in 0..npix {
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

/*--------------------------------------------------------------------------*/
/// Compute mean and RMS sigma of the non-null pixels in the input array.
fn FnMeanSigma_double(
    array: &[f64],        /*  2 dimensional array of image pixels */
    npix: usize,          /* number of pixels in the image */
    nullcheck: bool,      /* check for null values, if true */
    nullvalue: f64,       /* value of null pixels, if nullcheck is true */
    ngoodpix: &mut usize, /* number of non-null pixels in the image */
    mean: &mut f64,       /* returned mean value of all non-null pixels */
    sigma: &mut f64,      /* returned R.M.S. value of all non-null pixels */
    status: &mut c_int,   /* error status */
) -> c_int {
    let mut ngood: usize = 0;

    let mut sum: f64 = 0.0;
    let mut sum2: f64 = 0.0;
    let mut xtemp: f64;

    let mut value: usize = 0;

    if nullcheck {
        for ii in 0..npix {
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
        for ii in 0..npix {
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

/*--------------------------------------------------------------------------*/
/// Estimate the median and background noise in the input image using 2nd, 3rd and 5th order Median Absolute Differences.
///
/// The noise in the background of the image is calculated using the MAD algorithms
/// developed for deriving the signal to noise ratio in spectra
/// (see issue #42 of the ST-ECF newsletter, http://www.stecf.org/documents/newsletter/)
///
/// 3rd order:  noise = 1.482602 / sqrt(6) * median (abs(2*flux(i) - flux(i-2) - flux(i+2)))
///
/// The returned estimates are the median of the values that are computed for each
/// row of the image.
fn FnNoise5_short(
    array: &[c_short],            /*  2 dimensional array of image pixels */
    nx: usize,                    /* number of pixels in each row of the image */
    ny: usize,                    /* number of rows in the image */
    nullcheck: bool,              /* check for null values, if true */
    nullvalue: c_short,           /* value of null pixels, if nullcheck is true */
    ngood: Option<&mut usize>,    /* number of good, non-null pixels? */
    minval: Option<&mut c_short>, /* minimum non-null value */
    maxval: Option<&mut c_short>, /* maximum non-null value */
    noise2: Option<&mut f64>,     /* returned 2nd order MAD of all non-null pixels */
    noise3: Option<&mut f64>,     /* returned 3rd order MAD of all non-null pixels */
    noise5: Option<&mut f64>,     /* returned 5th order MAD of all non-null pixels */
    status: &mut c_int,           /* error status */
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
    let mut xmaxval = c_short::MAX;
    let mut xminval = c_short::MIN;

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

        minval.is_some_set(xminval);
        maxval.is_some_set(xmaxval);
        ngood.is_some_set(ngoodpix);
        noise2.is_some_set(0.);
        noise3.is_some_set(0.);
        noise5.is_some_set(0.);

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

        /***** find the first valid pixel in row */
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

        /***** find the 2nd valid pixel in row (which we will skip over) */
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

        /***** find the 3rd valid pixel in row */
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
                differences2[nvals2] = ((v5 as c_int) - (v7 as c_int)).abs();

                nvals2 += 1;
            }

            if !(v3 == v4 && v4 == v5 && v5 == v6 && v6 == v7) {
                differences3[nvals] = ((2 * (v5 as c_int)) - (v3 as c_int) - (v7 as c_int)).abs();
                differences5[nvals] =
                    ((6 * (v5 as c_int)) - (4 * (v3 as c_int)) - (4 * (v7 as c_int))
                        + (v1 as c_int)
                        + (v9 as c_int))
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

    minval.is_some_set(xminval);
    maxval.is_some_set(xmaxval);
    ngood.is_some_set(ngoodpix);
    noise2.is_some_set(1.0483579 * xnoise2);
    noise3.is_some_set(0.6052697 * xnoise3);
    noise5.is_some_set(0.1772048 * xnoise5);

    *status
}

/*--------------------------------------------------------------------------*/
/// Estimate the median and background noise in the input image using 2nd, 3rd and 5th order Median Absolute Differences.
///
/// The noise in the background of the image is calculated using the MAD algorithms
/// developed for deriving the signal to noise ratio in spectra
/// (see issue #42 of the ST-ECF newsletter, http://www.stecf.org/documents/newsletter/)
///
/// 3rd order:  noise = 1.482602 / sqrt(6) * median (abs(2*flux(i) - flux(i-2) - flux(i+2)))
///
/// The returned estimates are the median of the values that are computed for each
/// row of the image.
fn FnNoise5_int(
    array: &[c_int],            /*  2 dimensional array of image pixels */
    nx: usize,                  /* number of pixels in each row of the image */
    ny: usize,                  /* number of rows in the image */
    nullcheck: bool,            /* check for null values, if true */
    nullvalue: c_int,           /* value of null pixels, if nullcheck is true */
    ngood: Option<&mut usize>,  /* number of good, non-null pixels? */
    minval: Option<&mut c_int>, /* minimum non-null value */
    maxval: Option<&mut c_int>, /* maximum non-null value */
    noise2: Option<&mut f64>,   /* returned 2nd order MAD of all non-null pixels */
    noise3: Option<&mut f64>,   /* returned 3rd order MAD of all non-null pixels */
    noise5: Option<&mut f64>,   /* returned 5th order MAD of all non-null pixels */
    status: &mut c_int,         /* error status */
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
    let mut xminval: c_int = c_int::MIN;
    let mut xmaxval = c_int::MAX;

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

        minval.is_some_set(xminval);
        maxval.is_some_set(xmaxval);
        ngood.is_some_set(ngoodpix);
        noise2.is_some_set(0.);
        noise3.is_some_set(0.);
        noise5.is_some_set(0.);

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

        /***** find the first valid pixel in row */
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

        /***** find the 2nd valid pixel in row (which we will skip over) */
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

        /***** find the 3rd valid pixel in row */
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
                let tdiff = (v5 as LONGLONG) - (v7 as LONGLONG);
                if tdiff < 0 {
                    differences2[nvals2] = -tdiff;
                } else {
                    differences2[nvals2] = tdiff;
                }

                nvals2 += 1;
            }

            if !(v3 == v4 && v4 == v5 && v5 == v6 && v6 == v7) {
                let tdiff = (2 * (v5 as LONGLONG)) - (v3 as LONGLONG) - (v7 as LONGLONG);
                if tdiff < 0 {
                    differences3[nvals] = -tdiff;
                } else {
                    differences3[nvals] = tdiff;
                }

                let tdiff =
                    (6 * (v5 as LONGLONG)) - (4 * (v3 as LONGLONG)) - (4 * (v7 as LONGLONG))
                        + (v1 as LONGLONG)
                        + (v9 as LONGLONG);
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

    minval.is_some_set(xminval);
    maxval.is_some_set(xmaxval);
    ngood.is_some_set(ngoodpix);
    noise2.is_some_set(1.0483579 * xnoise2);
    noise3.is_some_set(0.6052697 * xnoise3);
    noise5.is_some_set(0.1772048 * xnoise5);

    *status
}

/*--------------------------------------------------------------------------*/
/// Estimate the median and background noise in the input image using 2nd, 3rd and 5th order Median Absolute Differences.
///
/// The noise in the background of the image is calculated using the MAD algorithms
/// developed for deriving the signal to noise ratio in spectra
/// (see issue #42 of the ST-ECF newsletter, http://www.stecf.org/documents/newsletter/)
///
/// 3rd order:  noise = 1.482602 / sqrt(6) * median (abs(2*flux(i) - flux(i-2) - flux(i+2)))
///
/// The returned estimates are the median of the values that are computed for each
/// row of the image.
fn FnNoise5_float(
    array: &[f32],             /*  2 dimensional array of image pixels */
    nx: usize,                 /* number of pixels in each row of the image */
    ny: usize,                 /* number of rows in the image */
    nullcheck: bool,           /* check for null values, if true */
    nullvalue: f32,            /* value of null pixels, if nullcheck is true */
    ngood: Option<&mut usize>, /* number of good, non-null pixels? */
    minval: Option<&mut f32>,  /* minimum non-null value */
    maxval: Option<&mut f32>,  /* maximum non-null value */
    noise2: Option<&mut f64>,  /* returned 2nd order MAD of all non-null pixels */
    noise3: Option<&mut f64>,  /* returned 3rd order MAD of all non-null pixels */
    noise5: Option<&mut f64>,  /* returned 5th order MAD of all non-null pixels */
    status: &mut c_int,        /* error status */
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
    let mut xminval: f32 = f32::MIN;
    let mut xmaxval = f32::MAX;

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

        minval.is_some_set(xminval);
        maxval.is_some_set(xmaxval);
        ngood.is_some_set(ngoodpix);
        noise2.is_some_set(0.);
        noise3.is_some_set(0.);
        noise5.is_some_set(0.);

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

        /***** find the first valid pixel in row */
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

        /***** find the 2nd valid pixel in row (which we will skip over) */
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

        /***** find the 3rd valid pixel in row */
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

    minval.is_some_set(xminval);
    maxval.is_some_set(xmaxval);
    ngood.is_some_set(ngoodpix);
    noise2.is_some_set(1.0483579 * xnoise2);
    noise3.is_some_set(0.6052697 * xnoise3);
    noise5.is_some_set(0.1772048 * xnoise5);

    *status
}

/*--------------------------------------------------------------------------*/
/// Estimate the median and background noise in the input image using 2nd, 3rd and 5th order Median Absolute Differences.
///
/// The noise in the background of the image is calculated using the MAD algorithms
/// developed for deriving the signal to noise ratio in spectra
/// (see issue #42 of the ST-ECF newsletter, http://www.stecf.org/documents/newsletter/)
///
/// 3rd order:  noise = 1.482602 / sqrt(6) * median (abs(2*flux(i) - flux(i-2) - flux(i+2)))
///
/// The returned estimates are the median of the values that are computed for each
/// row of the image.
fn FnNoise5_double(
    array: &[f64],             /*  2 dimensional array of image pixels */
    nx: usize,                 /* number of pixels in each row of the image */
    ny: usize,                 /* number of rows in the image */
    nullcheck: bool,           /* check for null values, if true */
    nullvalue: f64,            /* value of null pixels, if nullcheck is true */
    ngood: Option<&mut usize>, /* number of good, non-null pixels? */
    minval: Option<&mut f64>,  /* minimum non-null value */
    maxval: Option<&mut f64>,  /* maximum non-null value */
    noise2: Option<&mut f64>,  /* returned 2nd order MAD of all non-null pixels */
    noise3: Option<&mut f64>,  /* returned 3rd order MAD of all non-null pixels */
    noise5: Option<&mut f64>,  /* returned 5th order MAD of all non-null pixels */
    status: &mut c_int,        /* error status */
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
    let mut xminval: f64 = f64::MIN;
    let mut xmaxval = f64::MAX;

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

        minval.is_some_set(xminval);
        maxval.is_some_set(xmaxval);
        ngood.is_some_set(ngoodpix);
        noise2.is_some_set(0.);
        noise3.is_some_set(0.);
        noise5.is_some_set(0.);

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

        /***** find the first valid pixel in row */
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

        /***** find the 2nd valid pixel in row (which we will skip over) */
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

        /***** find the 3rd valid pixel in row */
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

    minval.is_some_set(xminval);
    maxval.is_some_set(xmaxval);
    ngood.is_some_set(ngoodpix);
    noise2.is_some_set(1.0483579 * xnoise2);
    noise3.is_some_set(0.6052697 * xnoise3);
    noise5.is_some_set(0.1772048 * xnoise5);

    *status
}

/*--------------------------------------------------------------------------*/
/// Estimate the median and background noise in the input image using 3rd order differences.
///
/// The noise in the background of the image is calculated using the 3rd order algorithm
/// developed for deriving the signal to noise ratio in spectra
/// (see issue #42 of the ST-ECF newsletter, http://www.stecf.org/documents/newsletter/)
///
///   noise = 1.482602 / sqrt(6) * median (abs(2*flux(i) - flux(i-2) - flux(i+2)))
///
/// The returned estimates are the median of the values that are computed for each
/// row of the image.
fn FnNoise3_short(
    array: &[c_short],            /*  2 dimensional array of image pixels */
    nx: usize,                    /* number of pixels in each row of the image */
    ny: usize,                    /* number of rows in the image */
    nullcheck: bool,              /* check for null values, if true */
    nullvalue: c_short,           /* value of null pixels, if nullcheck is true */
    ngood: Option<&mut usize>,    /* number of good, non-null pixels? */
    minval: Option<&mut c_short>, /* minimum non-null value */
    maxval: Option<&mut c_short>, /* maximum non-null value */
    noise: Option<&mut f64>,      /* returned R.M.S. value of all non-null pixels */
    status: &mut c_int,           /* error status */
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
    let mut xminval: c_short = c_short::MIN; /* minimum non-null value */
    let mut xmaxval: c_short = -c_short::MAX; /* maximum non-null value */
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

        minval.is_some_set(xminval);
        maxval.is_some_set(xmaxval);
        ngood.is_some_set(ngoodpix);
        noise.is_some_set(0.0);

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

        /***** find the first valid pixel in row */
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

        /***** find the 2nd valid pixel in row (which we will skip over) */
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

        /***** find the 3rd valid pixel in row */
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

    minval.is_some_set(xminval);
    maxval.is_some_set(xmaxval);
    ngood.is_some_set(ngoodpix);
    noise.is_some_set(0.6052697 * xnoise);

    *status
}

/*--------------------------------------------------------------------------*/
/// Estimate the median and background noise in the input image using 3rd order differences.
///
/// The noise in the background of the image is calculated using the 3rd order algorithm
/// developed for deriving the signal to noise ratio in spectra
/// (see issue #42 of the ST-ECF newsletter, http://www.stecf.org/documents/newsletter/)
///
///   noise = 1.482602 / sqrt(6) * median (abs(2*flux(i) - flux(i-2) - flux(i+2)))
///
/// The returned estimates are the median of the values that are computed for each
/// row of the image.
fn FnNoise3_int(
    array: &[c_int],            /*  2 dimensional array of image pixels */
    nx: usize,                  /* number of pixels in each row of the image */
    ny: usize,                  /* number of rows in the image */
    nullcheck: bool,            /* check for null values, if true */
    nullvalue: c_int,           /* value of null pixels, if nullcheck is true */
    ngood: Option<&mut usize>,  /* number of good, non-null pixels? */
    minval: Option<&mut c_int>, /* minimum non-null value */
    maxval: Option<&mut c_int>, /* maximum non-null value */
    noise: Option<&mut f64>,    /* returned R.M.S. value of all non-null pixels */
    status: &mut c_int,         /* error status */
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
    let mut xminval: c_int = c_int::MIN; /* minimum non-null value */
    let mut xmaxval: c_int = -c_int::MAX; /* maximum non-null value */
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

        minval.is_some_set(xminval);
        maxval.is_some_set(xmaxval);
        ngood.is_some_set(ngoodpix);
        noise.is_some_set(0.0);

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

        /***** find the first valid pixel in row */
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

        /***** find the 2nd valid pixel in row (which we will skip over) */
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

        /***** find the 3rd valid pixel in row */
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

    minval.is_some_set(xminval);
    maxval.is_some_set(xmaxval);
    ngood.is_some_set(ngoodpix);
    noise.is_some_set(0.6052697 * xnoise);

    *status
}

/*--------------------------------------------------------------------------*/
/// Estimate the median and background noise in the input image using 3rd order differences.
///
/// The noise in the background of the image is calculated using the 3rd order algorithm
/// developed for deriving the signal to noise ratio in spectra
/// (see issue #42 of the ST-ECF newsletter, http://www.stecf.org/documents/newsletter/)
///
///   noise = 1.482602 / sqrt(6) * median (abs(2*flux(i) - flux(i-2) - flux(i+2)))
///
/// The returned estimates are the median of the values that are computed for each
/// row of the image.
fn FnNoise3_float(
    array: &[f32],             /*  2 dimensional array of image pixels */
    nx: usize,                 /* number of pixels in each row of the image */
    ny: usize,                 /* number of rows in the image */
    nullcheck: bool,           /* check for null values, if true */
    nullvalue: f32,            /* value of null pixels, if nullcheck is true */
    ngood: Option<&mut usize>, /* number of good, non-null pixels? */
    minval: Option<&mut f32>,  /* minimum non-null value */
    maxval: Option<&mut f32>,  /* maximum non-null value */
    noise: Option<&mut f64>,   /* returned R.M.S. value of all non-null pixels */
    status: &mut c_int,        /* error status */
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
    let mut xminval: f32 = f32::MIN; /* minimum non-null value */
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

        minval.is_some_set(xminval);
        maxval.is_some_set(xmaxval);
        ngood.is_some_set(ngoodpix);
        noise.is_some_set(0.0);

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

        /***** find the first valid pixel in row */
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

        /***** find the 2nd valid pixel in row (which we will skip over) */
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

        /***** find the 3rd valid pixel in row */
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

    minval.is_some_set(xminval);
    maxval.is_some_set(xmaxval);
    ngood.is_some_set(ngoodpix);
    noise.is_some_set(0.6052697 * xnoise);

    *status
}

/*--------------------------------------------------------------------------*/
/// Estimate the median and background noise in the input image using 3rd order differences.
///
///The noise in the background of the image is calculated using the 3rd order algorithm
///developed for deriving the signal to noise ratio in spectra
///(see issue #42 of the ST-ECF newsletter, http://www.stecf.org/documents/newsletter/)
///
///  noise = 1.482602 / sqrt(6) * median (abs(2*flux(i) - flux(i-2) - flux(i+2)))
///
///The returned estimates are the median of the values that are computed for each
///row of the image.
fn FnNoise3_double(
    array: &[f64],             /*  2 dimensional array of image pixels */
    nx: usize,                 /* number of pixels in each row of the image */
    ny: usize,                 /* number of rows in the image */
    nullcheck: bool,           /* check for null values, if true */
    nullvalue: f64,            /* value of null pixels, if nullcheck is true */
    ngood: Option<&mut usize>, /* number of good, non-null pixels? */
    minval: Option<&mut f64>,  /* minimum non-null value */
    maxval: Option<&mut f64>,  /* maximum non-null value */
    noise: Option<&mut f64>,   /* returned R.M.S. value of all non-null pixels */
    status: &mut c_int,        /* error status */
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
    let mut xminval: f64 = f64::MIN; /* minimum non-null value */
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

        minval.is_some_set(xminval);
        maxval.is_some_set(xmaxval);
        ngood.is_some_set(ngoodpix);
        noise.is_some_set(0.0);

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

        /***** find the first valid pixel in row */
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

        /***** find the 2nd valid pixel in row (which we will skip over) */
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

        /***** find the 3rd valid pixel in row */
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

    minval.is_some_set(xminval);
    maxval.is_some_set(xmaxval);
    ngood.is_some_set(ngoodpix);
    noise.is_some_set(0.6052697 * xnoise);

    *status
}

/*--------------------------------------------------------------------------*/
/// Estimate the background noise in the input image using sigma of 1st order differences.
///
///   noise = 1.0 / sqrt(2) * rms of (flux[i] - flux[i-1])
///
/// The returned estimate is the median of the values that are computed for each
/// row of the image.
fn FnNoise1_short(
    array: &[c_short],  /*  2 dimensional array of image pixels */
    nx: usize,          /* number of pixels in each row of the image */
    ny: usize,          /* number of rows in the image */
    nullcheck: bool,    /* check for null values, if true */
    nullvalue: c_short, /* value of null pixels, if nullcheck is true */
    noise: &mut f64,    /* returned R.M.S. value of all non-null pixels */
    status: &mut c_int, /* error status */
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

        /***** find the first valid pixel in row */
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
                        if (differences[ii] as f64 - mean).abs() < SIGMA_CLIP * stdev {
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

    *noise *= std::f64::consts::FRAC_1_SQRT_2;

    *status
}

/*--------------------------------------------------------------------------*/
/// Estimate the background noise in the input image using sigma of 1st order differences.
///
///   noise = 1.0 / sqrt(2) * rms of (flux[i] - flux[i-1])
///
/// The returned estimate is the median of the values that are computed for each
/// row of the image.
fn FnNoise1_int(
    array: &[c_int],    /*  2 dimensional array of image pixels */
    nx: usize,          /* number of pixels in each row of the image */
    ny: usize,          /* number of rows in the image */
    nullcheck: bool,    /* check for null values, if true */
    nullvalue: c_int,   /* value of null pixels, if nullcheck is true */
    noise: &mut f64,    /* returned R.M.S. value of all non-null pixels */
    status: &mut c_int, /* error status */
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

        /***** find the first valid pixel in row */
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
                        if (differences[ii] as f64 - mean).abs() < SIGMA_CLIP * stdev {
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

    *noise *= std::f64::consts::FRAC_1_SQRT_2;

    *status
}

/*--------------------------------------------------------------------------*/
/// Estimate the background noise in the input image using sigma of 1st order differences.
///
///   noise = 1.0 / sqrt(2) * rms of (flux[i] - flux[i-1])
///
/// The returned estimate is the median of the values that are computed for each
/// row of the image.
fn FnNoise1_float(
    array: &[f32],      /*  2 dimensional array of image pixels */
    nx: usize,          /* number of pixels in each row of the image */
    ny: usize,          /* number of rows in the image */
    nullcheck: bool,    /* check for null values, if true */
    nullvalue: f32,     /* value of null pixels, if nullcheck is true */
    noise: &mut f64,    /* returned R.M.S. value of all non-null pixels */
    status: &mut c_int, /* error status */
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

        /***** find the first valid pixel in row */
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
                        if (differences[ii] as f64 - mean).abs() < SIGMA_CLIP * stdev {
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

    *noise *= std::f64::consts::FRAC_1_SQRT_2;

    *status
}

/*--------------------------------------------------------------------------*/
/// Estimate the background noise in the input image using sigma of 1st order differences.
///
///   noise = 1.0 / sqrt(2) * rms of (flux[i] - flux[i-1])
///
/// The returned estimate is the median of the values that are computed for each
/// row of the image.
fn FnNoise1_double(
    array: &[f64],      /*  2 dimensional array of image pixels */
    nx: usize,          /* number of pixels in each row of the image */
    ny: usize,          /* number of rows in the image */
    nullcheck: bool,    /* check for null values, if true */
    nullvalue: f64,     /* value of null pixels, if nullcheck is true */
    noise: &mut f64,    /* returned R.M.S. value of all non-null pixels */
    status: &mut c_int, /* error status */
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

        /***** find the first valid pixel in row */
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

    *noise *= std::f64::consts::FRAC_1_SQRT_2;

    *status
}

/*--------------------------------------------------------------------------*/

/*
 *  These Quickselect routines are based on the algorithm described in
 *  "Numerical recipes in C", Second Edition,
 *  Cambridge University Press, 1992, Section 8.5, ISBN 0-521-43108-5
 *  This code by Nicolas Devillard - 1998. Public domain.
 */

/*--------------------------------------------------------------------------*/

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
                if arr[hh] <= arr[ll] {
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

/*--------------------------------------------------------------------------*/
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
                if arr[hh] <= arr[ll] {
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

/*--------------------------------------------------------------------------*/
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
                if arr[hh] <= arr[ll] {
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

/*--------------------------------------------------------------------------*/
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
                if arr[hh] <= arr[ll] {
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

/*--------------------------------------------------------------------------*/
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
                if arr[hh] <= arr[ll] {
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
