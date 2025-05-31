/*  This file, scalnull.c, contains the FITSIO routines used to define     */
/*  the starting heap address, the value scaling and the null values.      */
/*  The FITSIO software was written by William Pence at the High Energy    */
/*  Astrophysic Science Archive Research Center (HEASARC) at the NASA      */
/*  Goddard Space Flight Center.                                           */
/*--------------------------------------------------------------------------*/

use std::ffi::CStr;

use crate::c_types::*;

use bytemuck::cast_slice;

use crate::cs;
use crate::fitscore::{ffghdt_safe, ffmahd_safe, fits_is_compressed_image_safe};
use crate::fitsio::*;
use crate::modkey::ffukyj_safe;
use crate::wrappers::*;

/*--------------------------------------------------------------------------*/
/// Define the starting address for the heap for a binary table.
///
/// The default address is NAXIS1 * NAXIS2.  It is in units of
/// bytes relative to the beginning of the regular binary table data.
/// This routine also writes the appropriate THEAP keyword to the
/// FITS header.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffpthp(
    fptr: *mut fitsfile, /* I - FITS file pointer */
    theap: c_long,       /* I - starting addrss for the heap */
    status: *mut c_int,  /* IO - error status     */
) -> c_int {
    // FFI WRAPPER
    unsafe {
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);

        ffpthp_safe(fptr, theap, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Define the starting address for the heap for a binary table.
///
/// The default address is NAXIS1 * NAXIS2.  It is in units of
/// bytes relative to the beginning of the regular binary table data.
/// This routine also writes the appropriate THEAP keyword to the
/// FITS header.
pub(crate) fn ffpthp_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer */
    theap: c_long,       /* I - starting addrss for the heap */
    status: &mut c_int,  /* IO - error status     */
) -> c_int {
    if *status > 0 || theap < 1 {
        return *status;
    }

    /* reset position to the correct HDU if necessary */
    if fptr.HDUposition != fptr.Fptr.curhdu {
        ffmahd_safe(fptr, (fptr.HDUposition) + 1, None, status);
    }

    fptr.Fptr.heapstart = theap as LONGLONG;

    ffukyj_safe(
        fptr,
        cs!(c"THEAP"),
        theap as LONGLONG,
        Some(cs!(c"byte offset to heap area")),
        status,
    );
    *status
}

/*--------------------------------------------------------------------------*/
/// Define the linear scaling factor for the primary array or image extension
/// pixel values. This routine overrides the scaling values given by the
/// BSCALE and BZERO keywords if present.  Note that this routine does not
/// write or modify the BSCALE and BZERO keywords, but instead only modifies
/// the values temporarily in the internal buffer.  Thus, a subsequent call to
/// the ffrdef routine will reset the scaling back to the BSCALE and BZERO
/// keyword values (or 1. and 0. respectively if the keywords are not present).
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffpscl(
    fptr: *mut fitsfile, /* I - FITS file pointer               */
    scale: f64,          /* I - scaling factor: value of BSCALE */
    zero: f64,           /* I - zero point: value of BZERO      */
    status: *mut c_int,  /* IO - error status                   */
) -> c_int {
    // FFI WRAPPER
    unsafe {
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);

        ffpscl_safe(fptr, scale, zero, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Define the linear scaling factor for the primary array or image extension
/// pixel values. This routine overrides the scaling values given by the
/// BSCALE and BZERO keywords if present.  Note that this routine does not
/// write or modify the BSCALE and BZERO keywords, but instead only modifies
/// the values temporarily in the internal buffer.  Thus, a subsequent call to
/// the ffrdef routine will reset the scaling back to the BSCALE and BZERO
/// keyword values (or 1. and 0. respectively if the keywords are not present).
pub(crate) fn ffpscl_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer               */
    scale: f64,          /* I - scaling factor: value of BSCALE */
    zero: f64,           /* I - zero point: value of BZERO      */
    status: &mut c_int,  /* IO - error status                   */
) -> c_int {
    let mut hdutype = 0;
    if *status > 0 {
        return *status;
    }

    if scale == 0.0 {
        /* zero scale value is illegal */
        *status = ZERO_SCALE;
        return *status;
    }

    /* get HDU type */
    if ffghdt_safe(fptr, &mut hdutype, status) > 0 {
        return *status;
    }

    if hdutype != IMAGE_HDU {
        /* not proper HDU type */
        *status = NOT_IMAGE;
        return *status;
    }

    if fits_is_compressed_image_safe(fptr, status) != 0 {
        /* compressed images */
        fptr.Fptr.cn_bscale = scale;
        fptr.Fptr.cn_bzero = zero;
        return *status;
    }

    /* set pointer to the first 'column' (contains group parameters if any) */
    let c = fptr.Fptr.get_tableptr_as_mut_slice();
    let mut colptr = c[1]; /* increment to the 2nd 'column' pointer  (the image itself) */

    (colptr).tscale = scale;
    (colptr).tzero = zero;
    *status
}

/*--------------------------------------------------------------------------*/
/// Define the value used to represent undefined pixels in the primary array or
/// image extension. This only applies to integer image pixel (i.e. BITPIX > 0).
/// This routine overrides the null pixel value given by the BLANK keyword
/// if present.  Note that this routine does not write or modify the BLANK
/// keyword, but instead only modifies the value temporarily in the internal
/// buffer. Thus, a subsequent call to the ffrdef routine will reset the null
/// value back to the BLANK  keyword value (or not defined if the keyword is not
/// present).
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffpnul(
    fptr: *mut fitsfile, /* I - FITS file pointer                */
    nulvalue: LONGLONG,  /* I - null pixel value: value of BLANK */
    status: *mut c_int,  /* IO - error status                    */
) -> c_int {
    // FFI WRAPPER
    unsafe {
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);

        ffpnul_safe(fptr, nulvalue, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Define the value used to represent undefined pixels in the primary array or
/// image extension. This only applies to integer image pixel (i.e. BITPIX > 0).
/// This routine overrides the null pixel value given by the BLANK keyword
/// if present.  Note that this routine does not write or modify the BLANK
/// keyword, but instead only modifies the value temporarily in the internal
/// buffer. Thus, a subsequent call to the ffrdef routine will reset the null
/// value back to the BLANK  keyword value (or not defined if the keyword is not
/// present).
pub(crate) fn ffpnul_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer                */
    nulvalue: LONGLONG,  /* I - null pixel value: value of BLANK */
    status: &mut c_int,  /* IO - error status                    */
) -> c_int {
    let mut hdutype = 0;

    if *status > 0 {
        return *status;
    }

    /* get HDU type */
    if ffghdt_safe(fptr, &mut hdutype, status) > 0 {
        return *status;
    }

    if hdutype != IMAGE_HDU {
        /* not proper HDU type */
        *status = NOT_IMAGE;
        return *status;
    }

    if fits_is_compressed_image_safe(fptr, status) != 0 {
        /* ignore compressed images */
        return *status;
    }

    /* set pointer to the first 'column' (contains group parameters if any) */
    let c = fptr.Fptr.get_tableptr_as_mut_slice();
    let colptr = &mut c[1]; /* increment to the 2nd 'column' pointer  (the image itself) */
    colptr.tnull = nulvalue;
    *status
}

/*--------------------------------------------------------------------------*/
/// Define the linear scaling factor for the TABLE or BINTABLE extension
/// column values. This routine overrides the scaling values given by the
/// TSCALn and TZEROn keywords if present.  Note that this routine does not
/// write or modify the TSCALn and TZEROn keywords, but instead only modifies
/// the values temporarily in the internal buffer.  Thus, a subsequent call to
/// the ffrdef routine will reset the scaling back to the TSCALn and TZEROn
/// keyword values (or 1. and 0. respectively if the keywords are not present).
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn fftscl(
    fptr: *mut fitsfile, /* I - FITS file pointer                 */
    colnum: c_int,       /* I - column number to apply scaling to */
    scale: f64,          /* I - scaling factor: value of TSCALn   */
    zero: f64,           /* I - zero point: value of TZEROn       */
    status: *mut c_int,  /* IO - error status                     */
) -> c_int {
    // FFI WRAPPER
    unsafe {
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);

        fftscl_safe(fptr, colnum, scale, zero, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Define the linear scaling factor for the TABLE or BINTABLE extension
/// column values. This routine overrides the scaling values given by the
/// TSCALn and TZEROn keywords if present.  Note that this routine does not
/// write or modify the TSCALn and TZEROn keywords, but instead only modifies
/// the values temporarily in the internal buffer.  Thus, a subsequent call to
/// the ffrdef routine will reset the scaling back to the TSCALn and TZEROn
/// keyword values (or 1. and 0. respectively if the keywords are not present).
pub(crate) fn fftscl_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer                 */
    colnum: c_int,       /* I - column number to apply scaling to */
    scale: f64,          /* I - scaling factor: value of TSCALn   */
    zero: f64,           /* I - zero point: value of TZEROn       */
    status: &mut c_int,  /* IO - error status                     */
) -> c_int {
    let mut hdutype = 0;

    if *status > 0 {
        return *status;
    }

    if scale == 0.0 {
        /* zero scale value is illegal */
        *status = ZERO_SCALE;
        return *status;
    }

    /* get HDU type */
    if ffghdt_safe(fptr, &mut hdutype, status) > 0 {
        return *status;
    }

    if hdutype == IMAGE_HDU {
        /* not proper HDU type */
        *status = NOT_TABLE;
        return *status;
    }

    let c = fptr.Fptr.get_tableptr_as_mut_slice();
    let mut colptr = 0; /* set pointer to the first column */
    colptr += colnum as usize - 1; /* increment to the correct column */
    (c[colptr]).tscale = scale;
    (c[colptr]).tzero = zero;
    *status
}

/*--------------------------------------------------------------------------*/
/// Define the value used to represent undefined pixels in the BINTABLE column.
///
/// This only applies to integer datatype columns (TFORM = B, I, or J).
/// This routine overrides the null pixel value given by the TNULLn keyword
/// if present.  Note that this routine does not write or modify the TNULLn
/// keyword, but instead only modifies the value temporarily in the internal
/// buffer. Thus, a subsequent call to the ffrdef routine will reset the null
/// value back to the TNULLn  keyword value (or not defined if the keyword is not
/// present).
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn fftnul(
    fptr: *mut fitsfile, /* I - FITS file pointer                  */
    colnum: c_int,       /* I - column number to apply nulvalue to */
    nulvalue: LONGLONG,  /* I - null pixel value: value of TNULLn  */
    status: *mut c_int,  /* IO - error status                      */
) -> c_int {
    // FFI WRAPPER
    unsafe {
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);

        fftnul_safe(fptr, colnum, nulvalue, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Define the value used to represent undefined pixels in the BINTABLE column.
///
/// This only applies to integer datatype columns (TFORM = B, I, or J).
/// This routine overrides the null pixel value given by the TNULLn keyword
/// if present.  Note that this routine does not write or modify the TNULLn
/// keyword, but instead only modifies the value temporarily in the internal
/// buffer. Thus, a subsequent call to the ffrdef routine will reset the null
/// value back to the TNULLn  keyword value (or not defined if the keyword is not
/// present).
pub(crate) fn fftnul_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer                  */
    colnum: c_int,       /* I - column number to apply nulvalue to */
    nulvalue: LONGLONG,  /* I - null pixel value: value of TNULLn  */
    status: &mut c_int,  /* IO - error status                      */
) -> c_int {
    let mut hdutype = 0;
    if *status > 0 {
        return *status;
    }

    /* get HDU type */
    if ffghdt_safe(fptr, &mut hdutype, status) > 0 {
        return *status;
    }

    if hdutype != BINARY_TBL {
        /* not proper HDU type */
        *status = NOT_BTABLE;
        return *status;
    }

    let c = fptr.Fptr.get_tableptr_as_mut_slice(); /* set pointer to the first column */
    let mut colptr = 0; /* set pointer to the first column */
    colptr += colnum as usize - 1; /* increment to the correct column */
    c[colptr].tnull = nulvalue;
    *status
}

/*--------------------------------------------------------------------------*/
/// Define the string used to represent undefined pixels in the ASCII TABLE
/// column. This routine overrides the null  value given by the TNULLn keyword
/// if present.  Note that this routine does not write or modify the TNULLn
/// keyword, but instead only modifies the value temporarily in the internal
/// buffer. Thus, a subsequent call to the ffrdef routine will reset the null
/// value back to the TNULLn keyword value (or not defined if the keyword is not
/// present).
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffsnul(
    fptr: *mut fitsfile,      /* I - FITS file pointer                  */
    colnum: c_int,            /* I - column number to apply nulvalue to */
    nulstring: *const c_char, /* I - null pixel value: value of TNULLn  */
    status: *mut c_int,       /* IO - error status                      */
) -> c_int {
    // FFI WRAPPER
    unsafe {
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);

        let nulstring = CStr::from_ptr(nulstring).to_bytes_with_nul();

        ffsnul_safe(fptr, colnum, cast_slice(nulstring), status)
    }
}

/*--------------------------------------------------------------------------*/
/// Define the string used to represent undefined pixels in the ASCII TABLE
/// column. This routine overrides the null  value given by the TNULLn keyword
/// if present.  Note that this routine does not write or modify the TNULLn
/// keyword, but instead only modifies the value temporarily in the internal
/// buffer. Thus, a subsequent call to the ffrdef routine will reset the null
/// value back to the TNULLn keyword value (or not defined if the keyword is not
/// present).
pub(crate) fn ffsnul_safe(
    fptr: &mut fitsfile,  /* I - FITS file pointer                  */
    colnum: c_int,        /* I - column number to apply nulvalue to */
    nulstring: &[c_char], /* I - null pixel value: value of TNULLn  */
    status: &mut c_int,   /* IO - error status                      */
) -> c_int {
    let mut hdutype = 0;
    if *status > 0 {
        return *status;
    }

    /* get HDU type */
    if ffghdt_safe(fptr, &mut hdutype, status) > 0 {
        return *status;
    }

    if hdutype != ASCII_TBL {
        /* not proper HDU type */
        *status = NOT_ATABLE;
        return *status;
    }

    let c = fptr.Fptr.get_tableptr_as_mut_slice(); /* set pointer to the first column */
    let colptr = colnum as usize - 1; /* increment to the correct column */
    c[colptr].strnull[0] = 0;
    strncat_safe(&mut c[colptr].strnull, nulstring, 19); /* limit string to 19 chars */
    *status
}
