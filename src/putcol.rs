/*  This file, putcol.c, contains routines that write data elements to     */
/*  a FITS image or table. These are the generic routines.                 */

/*  The FITSIO software was written by William Pence at the High Energy    */
/*  Astrophysic Science Archive Research Center (HEASARC) at the NASA      */
/*  Goddard Space Flight Center.                                           */

use core::slice;
use std::ffi::CStr;

use crate::c_types::*;

use bytemuck::cast_slice;

use crate::NullValue;
use crate::calculate_subsection_length_unit;
use crate::fitscore::{ffgidm_safe, ffgisz_safe, ffgiszll_safe};
use crate::putcolb::{ffpclb_safe, ffpcnb_safe, ffppnb_safe, ffpprb_safe, ffpssb_safe};
use crate::putcold::{ffpcld_safe, ffpcnd_safe, ffppnd_safe, ffpprd_safe, ffpssd_safe};
use crate::putcole::{ffpcle_safe, ffpcne_safe, ffppne_safe, ffppre_safe, ffpsse_safe};
use crate::putcoli::{ffpcli_safe, ffpcni_safe, ffppni_safe, ffppri_safe, ffpssi_safe};
use crate::putcolj::{
    ffpclj_safe, ffpcljj_safe, ffpcnj_safe, ffpcnjj_safe, ffppnj_safe, ffppnjj_safe, ffpprj_safe,
    ffpprjj_safe, ffpssj_safe, ffpssjj_safe,
};
use crate::putcolk::{ffpclk_safe, ffpcnk_safe, ffppnk_safe, ffpprk_safe, ffpssk_safe};
use crate::putcoll::{ffpcll_safe, ffpclx_safe, ffpcnl_safe};
use crate::putcols::{ffpcls_safe, ffpcns_safe};
use crate::putcolsb::{ffpclsb_safe, ffpcnsb_safe, ffppnsb_safe, ffpprsb_safe, ffpsssb_safe};
use crate::putcolui::{ffpclui_safe, ffpcnui_safe, ffppnui_safe, ffpprui_safe, ffpssui_safe};
use crate::putcoluj::{
    ffpcluj_safe, ffpclujj_safe, ffpcnuj_safe, ffpcnujj_safe, ffppnuj_safe, ffppnujj_safe,
    ffppruj_safe, ffpprujj_safe, ffpssuj_safe, ffpssujj_safe,
};
use crate::putcoluk::{ffpcluk_safe, ffpcnuk_safe, ffppnuk_safe, ffppruk_safe, ffpssuk_safe};
use crate::raw_to_slice;
use crate::wrappers::*;
use crate::{bytes_per_datatype, fitsio::*};

/*--------------------------------------------------------------------------*/
/// Write an array of pixels to the primary array.  The datatype of the
/// input array is defined by the 2nd argument. Data conversion
/// and scaling will be performed if necessary (e.g, if the datatype of
/// the FITS array is not the same as the array being written).
///
/// This routine is simillar to ffppr, except it supports writing to
/// large images with more than 2**31 pixels.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffppx(
    fptr: *mut fitsfile,     /* I - FITS file pointer                       */
    datatype: c_int,         /* I - datatype of the value                   */
    firstpix: *const c_long, /* I - coord of  first pixel to write(1 based) */
    nelem: LONGLONG,         /* I - number of values to write               */
    array: *const c_void,    /* I - array of values that are written        */
    status: *mut c_int,      /* IO - error status                           */
) -> c_int {
    unsafe {
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);

        let bytes = match bytes_per_datatype(datatype) {
            Some(x) => x * nelem as usize,
            None => {
                *status = BAD_DATATYPE;
                return *status;
            }
        };

        let array = slice::from_raw_parts(array as *const u8, bytes);

        ffppx_safer(fptr, datatype, firstpix, nelem, array, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Write an array of pixels to the primary array.  The datatype of the
/// input array is defined by the 2nd argument. Data conversion
/// and scaling will be performed if necessary (e.g, if the datatype of
/// the FITS array is not the same as the array being written).
///
/// This routine is simillar to ffppr, except it supports writing to
/// large images with more than 2**31 pixels.
pub unsafe fn ffppx_safer(
    fptr: &mut fitsfile,     /* I - FITS file pointer                       */
    datatype: c_int,         /* I - datatype of the value                   */
    firstpix: *const c_long, /* I - coord of  first pixel to write(1 based) */
    nelem: LONGLONG,         /* I - number of values to write               */
    array: &[u8],            /* I - array of values that are written        */
    status: &mut c_int,      /* IO - error status                           */
) -> c_int {
    unsafe {
        let mut naxis: c_int = 0;
        let group: c_long = 1;
        let mut firstelem: LONGLONG = 0;
        let mut dimsize: LONGLONG = 1;
        let mut naxes: [LONGLONG; 9] = [0; 9];

        if *status > 0 {
            /* inherit input status value if > 0 */
            return *status;
        }

        /* get the size of the image */
        ffgidm_safe(fptr, &mut naxis, status);
        ffgiszll_safe(fptr, 9, &mut naxes, status);

        let firstpix = slice::from_raw_parts(firstpix, naxis as usize);

        firstelem = 0;
        for ii in 0..(naxis as usize) {
            firstelem += (firstpix[ii] as LONGLONG - 1) * dimsize;
            dimsize *= naxes[ii];
        }
        firstelem += 1;

        if datatype == TBYTE {
            let array = cast_slice(array);
            ffpprb_safe(fptr, group, firstelem, nelem, array, status);
        } else if datatype == TSBYTE {
            let array = cast_slice(array);
            ffpprsb_safe(fptr, group, firstelem, nelem, array, status);
        } else if datatype == TUSHORT {
            let array = cast_slice(array);
            ffpprui_safe(fptr, group, firstelem, nelem, array, status);
        } else if datatype == TSHORT {
            let array = cast_slice(array);
            ffppri_safe(fptr, group, firstelem, nelem, array, status);
        } else if datatype == TUINT {
            let array = cast_slice(array);
            ffppruk_safe(fptr, group, firstelem, nelem, array, status);
        } else if datatype == TINT {
            let array = cast_slice(array);
            ffpprk_safe(fptr, group, firstelem, nelem, array, status);
        } else if datatype == TULONG {
            let array = cast_slice(array);
            ffppruj_safe(fptr, group, firstelem, nelem, array, status);
        } else if datatype == TLONG {
            let array = cast_slice(array);
            ffpprj_safe(fptr, group, firstelem, nelem, array, status);
        } else if datatype == TULONGLONG {
            let array = cast_slice(array);
            ffpprujj_safe(fptr, group, firstelem, nelem, array, status);
        } else if datatype == TLONGLONG {
            let array = cast_slice(array);
            ffpprjj_safe(fptr, group, firstelem, nelem, array, status);
        } else if datatype == TFLOAT {
            let array = cast_slice(array);
            ffppre_safe(fptr, group, firstelem, nelem, array, status);
        } else if datatype == TDOUBLE {
            let array = cast_slice(array);
            ffpprd_safe(fptr, group, firstelem, nelem, array, status);
        } else {
            *status = BAD_DATATYPE;
        }

        *status
    }
}

/*--------------------------------------------------------------------------*/
/// Write an array of pixels to the primary array.  The datatype of the
/// input array is defined by the 2nd argument. Data conversion
/// and scaling will be performed if necessary (e.g, if the datatype of
/// the FITS array is not the same as the array being written).
///
/// This routine is simillar to ffppr, except it supports writing to
/// large images with more than 2**31 pixels.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffppxll(
    fptr: *mut fitsfile,       /* I - FITS file pointer                       */
    datatype: c_int,           /* I - datatype of the value                   */
    firstpix: *const LONGLONG, /* I - coord of  first pixel to write(1 based) */
    nelem: LONGLONG,           /* I - number of values to write               */
    array: *const c_void,      /* I - array of values that are written        */
    status: *mut c_int,        /* IO - error status                           */
) -> c_int {
    unsafe {
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);

        let bytes = match bytes_per_datatype(datatype) {
            Some(x) => x * nelem as usize,
            None => {
                *status = BAD_DATATYPE;
                return *status;
            }
        };

        let array = slice::from_raw_parts(array as *const u8, bytes);

        ffppxll_safer(fptr, datatype, firstpix, nelem, array, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Write an array of pixels to the primary array.  The datatype of the
/// input array is defined by the 2nd argument. Data conversion
/// and scaling will be performed if necessary (e.g, if the datatype of
/// the FITS array is not the same as the array being written).
///
/// This routine is simillar to ffppr, except it supports writing to
/// large images with more than 2**31 pixels.
pub unsafe fn ffppxll_safer(
    fptr: &mut fitsfile,       /* I - FITS file pointer                       */
    datatype: c_int,           /* I - datatype of the value                   */
    firstpix: *const LONGLONG, /* I - coord of  first pixel to write(1 based) */
    nelem: LONGLONG,           /* I - number of values to write               */
    array: &[u8],              /* I - array of values that are written        */
    status: &mut c_int,        /* IO - error status                           */
) -> c_int {
    unsafe {
        let mut naxis: c_int = 0;
        let group: c_long = 1;
        let mut firstelem: LONGLONG = 0;
        let mut dimsize = 1;
        let mut naxes: [LONGLONG; 9] = [0; 9];

        if *status > 0 {
            /* inherit input status value if > 0 */
            return *status;
        }

        /* get the size of the image */
        ffgidm_safe(fptr, &mut naxis, status);
        ffgiszll_safe(fptr, 9, &mut naxes, status);

        let firstpix = slice::from_raw_parts(firstpix, naxis as usize);

        firstelem = 0;
        for ii in 0..(naxis as usize) {
            firstelem += (firstpix[ii] - 1) * dimsize;
            dimsize *= naxes[ii];
        }
        firstelem += 1;

        if datatype == TBYTE {
            let array = cast_slice(array);
            ffpprb_safe(fptr, group, firstelem, nelem, array, status);
        } else if datatype == TSBYTE {
            let array = cast_slice(array);
            ffpprsb_safe(fptr, group, firstelem, nelem, array, status);
        } else if datatype == TUSHORT {
            let array = cast_slice(array);
            ffpprui_safe(fptr, group, firstelem, nelem, array, status);
        } else if datatype == TSHORT {
            let array = cast_slice(array);
            ffppri_safe(fptr, group, firstelem, nelem, array, status);
        } else if datatype == TUINT {
            let array = cast_slice(array);
            ffppruk_safe(fptr, group, firstelem, nelem, array, status);
        } else if datatype == TINT {
            let array = cast_slice(array);
            ffpprk_safe(fptr, group, firstelem, nelem, array, status);
        } else if datatype == TULONG {
            let array = cast_slice(array);
            ffppruj_safe(fptr, group, firstelem, nelem, array, status);
        } else if datatype == TLONG {
            let array = cast_slice(array);
            ffpprj_safe(fptr, group, firstelem, nelem, array, status);
        } else if datatype == TULONGLONG {
            let array = cast_slice(array);
            ffpprujj_safe(fptr, group, firstelem, nelem, array, status);
        } else if datatype == TLONGLONG {
            let array = cast_slice(array);
            ffpprjj_safe(fptr, group, firstelem, nelem, array, status);
        } else if datatype == TFLOAT {
            let array = cast_slice(array);
            ffppre_safe(fptr, group, firstelem, nelem, array, status);
        } else if datatype == TDOUBLE {
            let array = cast_slice(array);
            ffpprd_safe(fptr, group, firstelem, nelem, array, status);
        } else {
            *status = BAD_DATATYPE;
        }

        *status
    }
}

/*--------------------------------------------------------------------------*/
/// Write an array of values to the primary array.  The datatype of the
/// input array is defined by the 2nd argument. Data conversion
/// and scaling will be performed if necessary (e.g, if the datatype of
/// the FITS array is not the same as the array being written).
///
/// This routine supports writing to large images with
/// more than 2**31 pixels.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffppxn(
    fptr: *mut fitsfile,   /* I - FITS file pointer                       */
    datatype: c_int,       /* I - datatype of the value                   */
    firstpix: *mut c_long, /* I - first vector element to write(1 = 1st)  */
    nelem: LONGLONG,       /* I - number of values to write               */
    array: *const c_void,  /* I - array of values that are written        */
    nulval: *const c_void, /* I - pointer to the null value               */
    status: *mut c_int,    /* IO - error status                           */
) -> c_int {
    unsafe {
        let mut naxis = 0;

        let group: c_long = 1;

        let mut firstelem: LONGLONG = 0;
        let mut dimsize: LONGLONG = 1;
        let mut naxes: [LONGLONG; 9] = [0; 9];

        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        let bytes = match bytes_per_datatype(datatype) {
            Some(x) => x * nelem as usize,
            None => {
                *status = BAD_DATATYPE;
                return *status;
            }
        };

        let array = slice::from_raw_parts(array as *const u8, bytes);

        if *status > 0 {
            /* inherit input status value if > 0 */
            return *status;
        }

        if nulval.is_null() {
            /* null value not defined? */
            ffppx_safer(fptr, datatype, firstpix, nelem, array, status);
            return *status;
        }

        /* get the size of the image */
        ffgidm_safe(fptr, &mut naxis, status);
        ffgiszll_safe(fptr, 9, &mut naxes, status);

        firstelem = 0;
        let firstpix = slice::from_raw_parts(firstpix, naxis as usize);

        for ii in 0..(naxis as usize) {
            firstelem += (firstpix[ii] as LONGLONG - 1) * dimsize;
            dimsize *= naxes[ii];
        }
        firstelem += 1;

        if datatype == TBYTE {
            let array = cast_slice(array);
            ffppnb_safe(
                fptr,
                group,
                firstelem,
                nelem,
                array,
                *(nulval as *const u8),
                status,
            );
        } else if datatype == TSBYTE {
            let array = cast_slice(array);
            ffppnsb_safe(
                fptr,
                group,
                firstelem,
                nelem,
                array,
                *(nulval as *const i8),
                status,
            );
        } else if datatype == TUSHORT {
            let array = cast_slice(array);
            ffppnui_safe(
                fptr,
                group,
                firstelem,
                nelem,
                array,
                *(nulval as *const c_ushort),
                status,
            );
        } else if datatype == TSHORT {
            let array = cast_slice(array);
            ffppni_safe(
                fptr,
                group,
                firstelem,
                nelem,
                array,
                *(nulval as *const c_short),
                status,
            );
        } else if datatype == TUINT {
            let array = cast_slice(array);
            ffppnuk_safe(
                fptr,
                group,
                firstelem,
                nelem,
                array,
                *(nulval as *const c_uint),
                status,
            );
        } else if datatype == TINT {
            let array = cast_slice(array);
            ffppnk_safe(
                fptr,
                group,
                firstelem,
                nelem,
                array,
                *(nulval as *const c_int),
                status,
            );
        } else if datatype == TULONG {
            let array = cast_slice(array);
            ffppnuj_safe(
                fptr,
                group,
                firstelem,
                nelem,
                array,
                *(nulval as *const c_ulong),
                status,
            );
        } else if datatype == TLONG {
            let array = cast_slice(array);
            ffppnj_safe(
                fptr,
                group,
                firstelem,
                nelem,
                array,
                *(nulval as *const c_long),
                status,
            );
        } else if datatype == TULONGLONG {
            let array = cast_slice(array);
            ffppnujj_safe(
                fptr,
                group,
                firstelem,
                nelem,
                array,
                *(nulval as *const ULONGLONG),
                status,
            );
        } else if datatype == TLONGLONG {
            let array = cast_slice(array);
            ffppnjj_safe(
                fptr,
                group,
                firstelem,
                nelem,
                array,
                *(nulval as *const LONGLONG),
                status,
            );
        } else if datatype == TFLOAT {
            let array = cast_slice(array);
            ffppne_safe(
                fptr,
                group,
                firstelem,
                nelem,
                array,
                *(nulval as *const f32),
                status,
            );
        } else if datatype == TDOUBLE {
            let array = cast_slice(array);
            ffppnd_safe(
                fptr,
                group,
                firstelem,
                nelem,
                array,
                *(nulval as *const f64),
                status,
            );
        } else {
            *status = BAD_DATATYPE;
        }
        *status
    }
}

/*--------------------------------------------------------------------------*/
/// Write an array of values to the primary array.  The datatype of the
/// input array is defined by the 2nd argument. Data conversion
/// and scaling will be performed if necessary (e.g, if the datatype of
/// the FITS array is not the same as the array being written).
///
/// This routine supports writing to large images with
/// more than 2**31 pixels.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffppxnll(
    fptr: *mut fitsfile,       /* I - FITS file pointer                       */
    datatype: c_int,           /* I - datatype of the value                   */
    firstpix: *const LONGLONG, /* I - first vector element to write(1 = 1st)  */
    nelem: LONGLONG,           /* I - number of values to write               */
    array: *const c_void,      /* I - array of values that are written        */
    nulval: *const c_void,     /* I - pointer to the null value               */
    status: *mut c_int,        /* IO - error status                           */
) -> c_int {
    unsafe {
        let mut naxis = 0;

        let group: c_long = 1;

        let mut firstelem: LONGLONG = 0;
        let mut dimsize: LONGLONG = 1;
        let mut naxes: [LONGLONG; 9] = [0; 9];

        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        let bytes = match bytes_per_datatype(datatype) {
            Some(x) => x * nelem as usize,
            None => {
                *status = BAD_DATATYPE;
                return *status;
            }
        };

        let array = slice::from_raw_parts(array as *const u8, bytes);

        if *status > 0 {
            /* inherit input status value if > 0 */
            return *status;
        }

        if nulval.is_null() {
            /* null value not defined? */
            ffppxll_safer(fptr, datatype, firstpix, nelem, array, status);
            return *status;
        }

        /* get the size of the image */
        ffgidm_safe(fptr, &mut naxis, status);
        ffgiszll_safe(fptr, 9, &mut naxes, status);

        firstelem = 0;
        let firstpix = slice::from_raw_parts(firstpix, naxis as usize);

        for ii in 0..(naxis as usize) {
            firstelem += (firstpix[ii] - 1) * dimsize;
            dimsize *= naxes[ii];
        }
        firstelem += 1;

        if datatype == TBYTE {
            let array = cast_slice(array);
            ffppnb_safe(
                fptr,
                group,
                firstelem,
                nelem,
                array,
                *(nulval as *const u8),
                status,
            );
        } else if datatype == TSBYTE {
            let array = cast_slice(array);
            ffppnsb_safe(
                fptr,
                group,
                firstelem,
                nelem,
                array,
                *(nulval as *const i8),
                status,
            );
        } else if datatype == TUSHORT {
            let array = cast_slice(array);
            ffppnui_safe(
                fptr,
                group,
                firstelem,
                nelem,
                array,
                *(nulval as *const c_ushort),
                status,
            );
        } else if datatype == TSHORT {
            let array = cast_slice(array);
            ffppni_safe(
                fptr,
                group,
                firstelem,
                nelem,
                array,
                *(nulval as *const c_short),
                status,
            );
        } else if datatype == TUINT {
            let array = cast_slice(array);
            ffppnuk_safe(
                fptr,
                group,
                firstelem,
                nelem,
                array,
                *(nulval as *const c_uint),
                status,
            );
        } else if datatype == TINT {
            let array = cast_slice(array);
            ffppnk_safe(
                fptr,
                group,
                firstelem,
                nelem,
                array,
                *(nulval as *const c_int),
                status,
            );
        } else if datatype == TULONG {
            let array = cast_slice(array);
            ffppnuj_safe(
                fptr,
                group,
                firstelem,
                nelem,
                array,
                *(nulval as *const c_ulong),
                status,
            );
        } else if datatype == TLONG {
            let array = cast_slice(array);
            ffppnj_safe(
                fptr,
                group,
                firstelem,
                nelem,
                array,
                *(nulval as *const c_long),
                status,
            );
        } else if datatype == TULONGLONG {
            let array = cast_slice(array);
            ffppnujj_safe(
                fptr,
                group,
                firstelem,
                nelem,
                array,
                *(nulval as *const ULONGLONG),
                status,
            );
        } else if datatype == TLONGLONG {
            let array = cast_slice(array);
            ffppnjj_safe(
                fptr,
                group,
                firstelem,
                nelem,
                array,
                *(nulval as *const LONGLONG),
                status,
            );
        } else if datatype == TFLOAT {
            let array = cast_slice(array);
            ffppne_safe(
                fptr,
                group,
                firstelem,
                nelem,
                array,
                *(nulval as *const f32),
                status,
            );
        } else if datatype == TDOUBLE {
            let array = cast_slice(array);

            ffppnd_safe(
                fptr,
                group,
                firstelem,
                nelem,
                array,
                *(nulval as *const f64),
                status,
            );
        } else {
            *status = BAD_DATATYPE;
        }

        *status
    }
}

/*--------------------------------------------------------------------------*/
/// Write an array of values to the primary array.  The datatype of the
/// input array is defined by the 2nd argument. Data conversion
/// and scaling will be performed if necessary (e.g, if the datatype of
/// the FITS array is not the same as the array being written).
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffppr(
    fptr: *mut fitsfile,  /* I - FITS file pointer                       */
    datatype: c_int,      /* I - datatype of the value                   */
    firstelem: LONGLONG,  /* I - first vector element to write(1 = 1st)  */
    nelem: LONGLONG,      /* I - number of values to write               */
    array: *const c_void, /* I - array of values that are written        */
    status: *mut c_int,   /* IO - error status                           */
) -> c_int {
    unsafe {
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);

        let bytes = match bytes_per_datatype(datatype) {
            Some(x) => x * nelem as usize,
            None => {
                *status = BAD_DATATYPE;
                return *status;
            }
        };

        let array = slice::from_raw_parts(array as *const u8, bytes);

        ffppr_safe(fptr, datatype, firstelem, nelem, array, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Write an array of values to the primary array.  The datatype of the
/// input array is defined by the 2nd argument. Data conversion
/// and scaling will be performed if necessary (e.g, if the datatype of
/// the FITS array is not the same as the array being written).
pub fn ffppr_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer                       */
    datatype: c_int,     /* I - datatype of the value                   */
    firstelem: LONGLONG, /* I - first vector element to write(1 = 1st)  */
    nelem: LONGLONG,     /* I - number of values to write               */
    array: &[u8],        /* I - array of values that are written        */
    status: &mut c_int,  /* IO - error status                           */
) -> c_int {
    let group: c_long = 1;

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    if datatype == TBYTE {
        let array = cast_slice(array);
        ffpprb_safe(fptr, group, firstelem, nelem, array, status);
    } else if datatype == TSBYTE {
        let array = cast_slice(array);
        ffpprsb_safe(fptr, group, firstelem, nelem, array, status);
    } else if datatype == TUSHORT {
        let array = cast_slice(array);
        ffpprui_safe(fptr, group, firstelem, nelem, array, status);
    } else if datatype == TSHORT {
        let array = cast_slice(array);
        ffppri_safe(fptr, group, firstelem, nelem, array, status);
    } else if datatype == TUINT {
        let array = cast_slice(array);
        ffppruk_safe(fptr, group, firstelem, nelem, array, status);
    } else if datatype == TINT {
        let array = cast_slice(array);
        ffpprk_safe(fptr, group, firstelem, nelem, array, status);
    } else if datatype == TULONG {
        let array = cast_slice(array);
        ffppruj_safe(fptr, group, firstelem, nelem, array, status);
    } else if datatype == TLONG {
        let array = cast_slice(array);
        ffpprj_safe(fptr, group, firstelem, nelem, array, status);
    } else if datatype == TULONGLONG {
        let array = cast_slice(array);
        ffpprujj_safe(fptr, group, firstelem, nelem, array, status);
    } else if datatype == TLONGLONG {
        let array = cast_slice(array);
        ffpprjj_safe(fptr, group, firstelem, nelem, array, status);
    } else if datatype == TFLOAT {
        let array = cast_slice(array);
        ffppre_safe(fptr, group, firstelem, nelem, array, status);
    } else if datatype == TDOUBLE {
        let array = cast_slice(array);
        ffpprd_safe(fptr, group, firstelem, nelem, array, status);
    } else {
        *status = BAD_DATATYPE;
    }
    *status
}

/*--------------------------------------------------------------------------*/
/// Write an array of values to the primary array.  The datatype of the
/// input array is defined by the 2nd argument. Data conversion
/// and scaling will be performed if necessary (e.g, if the datatype of
/// the FITS array is not the same as the array being written).
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffppn(
    fptr: *mut fitsfile,   /* I - FITS file pointer                       */
    datatype: c_int,       /* I - datatype of the value                   */
    firstelem: LONGLONG,   /* I - first vector element to write(1 = 1st)  */
    nelem: LONGLONG,       /* I - number of values to write               */
    array: *const c_void,  /* I - array of values that are written        */
    nulval: *const c_void, /* I - pointer to the null value               */
    status: *mut c_int,    /* IO - error status                           */
) -> c_int {
    unsafe {
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);

        let bytes = match bytes_per_datatype(datatype) {
            Some(x) => x * nelem as usize,
            None => {
                *status = BAD_DATATYPE;
                return *status;
            }
        };

        let array = slice::from_raw_parts(array as *const u8, bytes);
        let nulval = NullValue::from_raw_ptr(datatype, nulval);

        ffppn_safe(fptr, datatype, firstelem, nelem, array, nulval, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Write an array of values to the primary array.  The datatype of the
/// input array is defined by the 2nd argument. Data conversion
/// and scaling will be performed if necessary (e.g, if the datatype of
/// the FITS array is not the same as the array being written).
pub fn ffppn_safe(
    fptr: &mut fitsfile,       /* I - FITS file pointer                       */
    datatype: c_int,           /* I - datatype of the value                   */
    firstelem: LONGLONG,       /* I - first vector element to write(1 = 1st)  */
    nelem: LONGLONG,           /* I - number of values to write               */
    array: &[u8],              /* I - array of values that are written        */
    nulval: Option<NullValue>, /* I - pointer to the null value               */
    status: &mut c_int,        /* IO - error status                           */
) -> c_int {
    let group: c_long = 1;

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    if nulval.is_none() {
        /* null value not defined? */
        ffppr_safe(fptr, datatype, firstelem, nelem, array, status);
        return *status;
    }

    let nulval = nulval.unwrap(); // safe to unwrap since we checked for None above

    match datatype {
        TBYTE => {
            let array = cast_slice(array);
            ffppnb_safe(
                fptr,
                group,
                firstelem,
                nelem,
                array,
                match nulval {
                    NullValue::UByte(x) => x,
                    _ => 0,
                },
                status,
            );
        }
        TSBYTE => {
            let array = cast_slice(array);
            ffppnsb_safe(
                fptr,
                group,
                firstelem,
                nelem,
                array,
                match nulval {
                    NullValue::Byte(x) => x,
                    _ => 0,
                },
                status,
            );
        }
        TUSHORT => {
            let array = cast_slice(array);
            ffppnui_safe(
                fptr,
                group,
                firstelem,
                nelem,
                array,
                match nulval {
                    NullValue::UShort(x) => x,
                    _ => 0,
                },
                status,
            );
        }
        TSHORT => {
            let array = cast_slice(array);
            ffppni_safe(
                fptr,
                group,
                firstelem,
                nelem,
                array,
                match nulval {
                    NullValue::Short(x) => x,
                    _ => 0,
                },
                status,
            );
        }
        TUINT => {
            let array = cast_slice(array);
            ffppnuk_safe(
                fptr,
                group,
                firstelem,
                nelem,
                array,
                match nulval {
                    NullValue::UInt(x) => x,
                    _ => 0,
                },
                status,
            );
        }
        TINT => {
            let array = cast_slice(array);
            ffppnk_safe(
                fptr,
                group,
                firstelem,
                nelem,
                array,
                match nulval {
                    NullValue::Int(x) => x,
                    _ => 0,
                },
                status,
            );
        }
        TULONG => {
            let array = cast_slice(array);
            ffppnuj_safe(
                fptr,
                group,
                firstelem,
                nelem,
                array,
                match nulval {
                    NullValue::ULong(x) => x,
                    _ => 0,
                },
                status,
            );
        }
        TLONG => {
            let array = cast_slice(array);
            ffppnj_safe(
                fptr,
                group,
                firstelem,
                nelem,
                array,
                match nulval {
                    NullValue::Long(x) => x,
                    _ => 0,
                },
                status,
            );
        }
        TULONGLONG => {
            let array = cast_slice(array);
            ffppnujj_safe(
                fptr,
                group,
                firstelem,
                nelem,
                array,
                match nulval {
                    NullValue::ULONGLONG(x) => x,
                    _ => 0,
                },
                status,
            );
        }
        TLONGLONG => {
            let array = cast_slice(array);
            ffppnjj_safe(
                fptr,
                group,
                firstelem,
                nelem,
                array,
                match nulval {
                    NullValue::LONGLONG(x) => x,
                    _ => 0,
                },
                status,
            );
        }
        TFLOAT => {
            let array = cast_slice(array);
            ffppne_safe(
                fptr,
                group,
                firstelem,
                nelem,
                array,
                match nulval {
                    NullValue::Float(x) => x,
                    _ => 0.0,
                },
                status,
            );
        }
        TDOUBLE => {
            let array = cast_slice(array);

            ffppnd_safe(
                fptr,
                group,
                firstelem,
                nelem,
                array,
                match nulval {
                    NullValue::Double(x) => x,
                    _ => 0.0,
                },
                status,
            );
        }
        _ => {
            *status = BAD_DATATYPE;
        }
    }

    *status
}

/*--------------------------------------------------------------------------*/
/// Write a section of values to the primary array. The datatype of the
/// input array is defined by the 2nd argument.  Data conversion
/// and scaling will be performed if necessary (e.g, if the datatype of
/// the FITS array is not the same as the array being written).
///
/// This routine supports writing to large images with
/// more than 2**31 pixels.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffpss(
    fptr: *mut fitsfile,  /* I - FITS file pointer                       */
    datatype: c_int,      /* I - datatype of the value                   */
    blc: *const c_long,   /* I - 'bottom left corner' of the subsection  */
    trc: *const c_long,   /* I - 'top right corner' of the subsection    */
    array: *const c_void, /* I - array of values that are written        */
    status: *mut c_int,   /* IO - error status                           */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        let mut naxis = 0;

        /* get the size of the image */
        ffgidm_safe(fptr, &mut naxis, status);

        let blc = slice::from_raw_parts(blc, naxis as usize);
        let trc = slice::from_raw_parts(trc, naxis as usize);

        let nelem = calculate_subsection_length_unit(blc, trc);

        let bytes = match bytes_per_datatype(datatype) {
            Some(x) => x * nelem as usize,
            None => {
                *status = BAD_DATATYPE;
                return *status;
            }
        };

        let array = slice::from_raw_parts(array as *const u8, bytes);

        ffpss_safe(fptr, datatype, blc, trc, array, status)
    }
}

/*--------------------------------------------------------------------------*/
/// Write a section of values to the primary array. The datatype of the
/// input array is defined by the 2nd argument.  Data conversion
/// and scaling will be performed if necessary (e.g, if the datatype of
/// the FITS array is not the same as the array being written).
///
/// This routine supports writing to large images with
/// more than 2**31 pixels.
pub fn ffpss_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer                       */
    datatype: c_int,     /* I - datatype of the value                   */
    blc: &[c_long],      /* I - 'bottom left corner' of the subsection  */
    trc: &[c_long],      /* I - 'top right corner' of the subsection    */
    array: &[u8],        /* I - array of values that are written        */
    status: &mut c_int,  /* IO - error status                           */
) -> c_int {
    let mut naxis: c_int = 0;
    let mut naxes: [c_long; 9] = [0; 9];

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    /* get the size of the image */
    ffgidm_safe(fptr, &mut naxis, status);
    ffgisz_safe(fptr, 9, &mut naxes, status);

    let naxis = naxis as c_long;

    if datatype == TBYTE {
        let array = cast_slice(array);
        ffpssb_safe(fptr, 1, naxis, &naxes, blc, trc, array, status);
    } else if datatype == TSBYTE {
        let array = cast_slice(array);
        ffpsssb_safe(fptr, 1, naxis, &naxes, blc, trc, array, status);
    } else if datatype == TUSHORT {
        let array = cast_slice(array);
        ffpssui_safe(fptr, 1, naxis, &naxes, blc, trc, array, status);
    } else if datatype == TSHORT {
        let array = cast_slice(array);
        ffpssi_safe(fptr, 1, naxis, &naxes, blc, trc, array, status);
    } else if datatype == TUINT {
        let array = cast_slice(array);
        ffpssuk_safe(fptr, 1, naxis, &naxes, blc, trc, array, status);
    } else if datatype == TINT {
        let array = cast_slice(array);
        ffpssk_safe(fptr, 1, naxis, &naxes, blc, trc, array, status);
    } else if datatype == TULONG {
        let array = cast_slice(array);
        ffpssuj_safe(fptr, 1, naxis, &naxes, blc, trc, array, status);
    } else if datatype == TLONG {
        let array = cast_slice(array);
        ffpssj_safe(fptr, 1, naxis, &naxes, blc, trc, array, status);
    } else if datatype == TULONGLONG {
        let array = cast_slice(array);
        ffpssujj_safe(fptr, 1, naxis, &naxes, blc, trc, array, status);
    } else if datatype == TLONGLONG {
        let array = cast_slice(array);
        ffpssjj_safe(fptr, 1, naxis, &naxes, blc, trc, array, status);
    } else if datatype == TFLOAT {
        let array = cast_slice(array);
        ffpsse_safe(fptr, 1, naxis, &naxes, blc, trc, array, status);
    } else if datatype == TDOUBLE {
        let array = cast_slice(array);
        ffpssd_safe(fptr, 1, naxis, &naxes, blc, trc, array, status);
    } else {
        *status = BAD_DATATYPE;
    }

    *status
}

/*--------------------------------------------------------------------------*/
/// Write an array of values to a table column.  The datatype of the
/// input array is defined by the 2nd argument. Data conversion
/// and scaling will be performed if necessary (e.g, if the datatype of
/// the FITS column is not the same as the array being written).
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffpcl(
    fptr: *mut fitsfile,  /* I - FITS file pointer                       */
    datatype: c_int,      /* I - datatype of the value                   */
    colnum: c_int,        /* I - number of column to write (1 = 1st col) */
    firstrow: LONGLONG,   /* I - first row to write (1 = 1st row)        */
    firstelem: LONGLONG,  /* I - first vector element to write (1 = 1st) */
    nelem: LONGLONG,      /* I - number of elements to write             */
    array: *const c_void, /* I - array of values that are written        */
    status: *mut c_int,   /* IO - error status                           */
) -> c_int {
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        let bytes = match bytes_per_datatype(datatype) {
            Some(x) => x * nelem as usize,
            None => {
                *status = BAD_DATATYPE;
                return *status;
            }
        };

        let array = slice::from_raw_parts(array as *const u8, bytes);

        ffpcl_safer(
            fptr,
            datatype,
            colnum,
            firstrow as LONGLONG,
            firstelem as LONGLONG,
            nelem as LONGLONG,
            array,
            status,
        )
    }
}

/*--------------------------------------------------------------------------*/
/// Write an array of values to a table column.  The datatype of the
/// input array is defined by the 2nd argument. Data conversion
/// and scaling will be performed if necessary (e.g, if the datatype of
/// the FITS column is not the same as the array being written).
pub unsafe fn ffpcl_safer(
    fptr: &mut fitsfile, /* I - FITS file pointer                       */
    datatype: c_int,     /* I - datatype of the value                   */
    colnum: c_int,       /* I - number of column to write (1 = 1st col) */
    firstrow: LONGLONG,  /* I - first row to write (1 = 1st row)        */
    firstelem: LONGLONG, /* I - first vector element to write (1 = 1st) */
    nelem: LONGLONG,     /* I - number of elements to write             */
    array: &[u8],        /* I - array of values that are written        */
    status: &mut c_int,  /* IO - error status                           */
) -> c_int {
    unsafe {
        if *status > 0 {
            /* inherit input status value if > 0 */
            return *status;
        }

        match datatype {
            TBIT => {
                let array = cast_slice(array);
                ffpclx_safe(
                    fptr,
                    colnum,
                    firstrow,
                    firstelem as c_long,
                    nelem as c_long,
                    array,
                    status,
                );
            }
            TBYTE => {
                let array = cast_slice(array);
                ffpclb_safe(
                    fptr,
                    colnum,
                    firstrow as LONGLONG,
                    firstelem as LONGLONG,
                    nelem as LONGLONG,
                    array,
                    status,
                );
            }
            TSBYTE => {
                let array = cast_slice(array);
                ffpclsb_safe(
                    fptr,
                    colnum,
                    firstrow as LONGLONG,
                    firstelem as LONGLONG,
                    nelem as LONGLONG,
                    array,
                    status,
                );
            }
            TUSHORT => {
                let array = cast_slice(array);
                ffpclui_safe(
                    fptr,
                    colnum,
                    firstrow as LONGLONG,
                    firstelem as LONGLONG,
                    nelem as LONGLONG,
                    array,
                    status,
                );
            }
            TSHORT => {
                let array = cast_slice(array);
                ffpcli_safe(
                    fptr,
                    colnum,
                    firstrow as LONGLONG,
                    firstelem as LONGLONG,
                    nelem as LONGLONG,
                    array,
                    status,
                );
            }
            TUINT => {
                let array = cast_slice(array);
                ffpcluk_safe(
                    fptr,
                    colnum,
                    firstrow as LONGLONG,
                    firstelem as LONGLONG,
                    nelem as LONGLONG,
                    array,
                    status,
                );
            }
            TINT => {
                let array = cast_slice(array);
                ffpclk_safe(
                    fptr,
                    colnum,
                    firstrow as LONGLONG,
                    firstelem as LONGLONG,
                    nelem as LONGLONG,
                    array,
                    status,
                );
            }
            TULONG => {
                let array = cast_slice(array);
                ffpcluj_safe(
                    fptr,
                    colnum,
                    firstrow as LONGLONG,
                    firstelem as LONGLONG,
                    nelem as LONGLONG,
                    array,
                    status,
                );
            }
            TLONG => {
                let array = cast_slice(array);
                ffpclj_safe(
                    fptr,
                    colnum,
                    firstrow as LONGLONG,
                    firstelem as LONGLONG,
                    nelem as LONGLONG,
                    array,
                    status,
                );
            }
            TULONGLONG => {
                let array = cast_slice(array);
                ffpclujj_safe(
                    fptr,
                    colnum,
                    firstrow as LONGLONG,
                    firstelem as LONGLONG,
                    nelem as LONGLONG,
                    array,
                    status,
                );
            }
            TLONGLONG => {
                let array = cast_slice(array);
                ffpcljj_safe(
                    fptr,
                    colnum,
                    firstrow as LONGLONG,
                    firstelem as LONGLONG,
                    nelem as LONGLONG,
                    array,
                    status,
                );
            }
            TFLOAT => {
                let array = cast_slice(array);
                ffpcle_safe(
                    fptr,
                    colnum,
                    firstrow as LONGLONG,
                    firstelem as LONGLONG,
                    nelem as LONGLONG,
                    array,
                    status,
                );
            }
            TDOUBLE => {
                let array = cast_slice(array);
                ffpcld_safe(
                    fptr,
                    colnum,
                    firstrow as LONGLONG,
                    firstelem as LONGLONG,
                    nelem as LONGLONG,
                    array,
                    status,
                );
            }
            TCOMPLEX => {
                let array = cast_slice(array);
                ffpcle_safe(
                    fptr,
                    colnum,
                    firstrow,
                    (firstelem - 1) * 2 + 1,
                    nelem * 2,
                    array,
                    status,
                );
            }
            TDBLCOMPLEX => {
                let array = cast_slice(array);
                ffpcld_safe(
                    fptr,
                    colnum,
                    firstrow,
                    (firstelem - 1) * 2 + 1,
                    nelem * 2,
                    array,
                    status,
                );
            }
            TLOGICAL => {
                let array = cast_slice(array);
                ffpcll_safe(
                    fptr,
                    colnum,
                    firstrow as LONGLONG,
                    firstelem as LONGLONG,
                    nelem as LONGLONG,
                    array,
                    status,
                );
            }
            TSTRING => {
                let array =
                    slice::from_raw_parts(array.as_ptr() as *const *const c_char, nelem as usize);
                let mut v_array = Vec::new();
                for item in array {
                    let array_item = slice::from_raw_parts(*item, FLEN_VALUE);
                    v_array.push(array_item);
                }

                ffpcls_safe(
                    fptr,
                    colnum,
                    firstrow as LONGLONG,
                    firstelem as LONGLONG,
                    nelem as LONGLONG,
                    &v_array,
                    status,
                );
            }
            _ => {
                *status = BAD_DATATYPE;
            }
        }
        *status
    }
}

/*--------------------------------------------------------------------------*/
/// Write an array of values to a table column.  The datatype of the
/// input array is defined by the 2nd argument. Data conversion
/// and scaling will be performed if necessary (e.g, if the datatype of
/// the FITS column is not the same as the array being written).
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffpcn(
    fptr: *mut fitsfile,   /* I - FITS file pointer                       */
    datatype: c_int,       /* I - datatype of the value                   */
    colnum: c_int,         /* I - number of column to write (1 = 1st col) */
    firstrow: LONGLONG,    /* I - first row to write (1 = 1st row)        */
    firstelem: LONGLONG,   /* I - first vector element to write (1 = 1st) */
    nelem: LONGLONG,       /* I - number of elements to write             */
    array: *const c_void,  /* I - array of values that are written        */
    nulval: *const c_void, /* I - pointer to the null value               */
    status: *mut c_int,    /* IO - error status                           */
) -> c_int {
    unsafe {
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let status = status.as_mut().expect(NULL_MSG);

        let bytes = match bytes_per_datatype(datatype) {
            Some(x) => x * nelem as usize,
            None => {
                *status = BAD_DATATYPE;
                return *status;
            }
        };

        let array: &[u8] = slice::from_raw_parts(array as *const u8, bytes);
        let nulval = NullValue::from_raw_ptr(datatype, nulval);

        ffpcn_safer(
            fptr,
            datatype,
            colnum,
            firstrow as LONGLONG,
            firstelem as LONGLONG,
            nelem as LONGLONG,
            array,
            nulval,
            status,
        )
    }
}

/*--------------------------------------------------------------------------*/
/// Write an array of values to a table column.  The datatype of the
/// input array is defined by the 2nd argument. Data conversion
/// and scaling will be performed if necessary (e.g, if the datatype of
/// the FITS column is not the same as the array being written).
pub unsafe fn ffpcn_safer(
    fptr: &mut fitsfile,       /* I - FITS file pointer                       */
    datatype: c_int,           /* I - datatype of the value                   */
    colnum: c_int,             /* I - number of column to write (1 = 1st col) */
    firstrow: LONGLONG,        /* I - first row to write (1 = 1st row)        */
    firstelem: LONGLONG,       /* I - first vector element to write (1 = 1st) */
    nelem: LONGLONG,           /* I - number of elements to write             */
    array: &[u8],              /* I - array of values that are written        */
    nulval: Option<NullValue>, /* I - pointer to the null value               */
    status: &mut c_int,        /* IO - error status                           */
) -> c_int {
    unsafe {
        if *status > 0 {
            /* inherit input status value if > 0 */
            return *status;
        }

        if nulval.is_none() {
            /* null value not defined? */
            ffpcl_safer(
                fptr,
                datatype,
                colnum,
                firstrow as LONGLONG,
                firstelem as LONGLONG,
                nelem as LONGLONG,
                array,
                status,
            );
            return *status;
        }

        let nulval = nulval.unwrap();

        if datatype == TBYTE {
            let array = cast_slice(array);
            ffpcnb_safe(
                fptr,
                colnum,
                firstrow,
                firstelem,
                nelem,
                array,
                match nulval {
                    NullValue::UByte(x) => x,
                    _ => 0,
                },
                status,
            );
        } else if datatype == TSBYTE {
            let array = cast_slice(array);
            ffpcnsb_safe(
                fptr,
                colnum,
                firstrow,
                firstelem,
                nelem,
                array,
                match nulval {
                    NullValue::Byte(x) => x,
                    _ => 0,
                },
                status,
            );
        } else if datatype == TUSHORT {
            let array = cast_slice(array);
            ffpcnui_safe(
                fptr,
                colnum,
                firstrow,
                firstelem,
                nelem,
                array,
                match nulval {
                    NullValue::UShort(x) => x,
                    _ => 0,
                },
                status,
            );
        } else if datatype == TSHORT {
            let array = cast_slice(array);
            ffpcni_safe(
                fptr,
                colnum,
                firstrow,
                firstelem,
                nelem,
                array,
                match nulval {
                    NullValue::Short(x) => x,
                    _ => 0,
                },
                status,
            );
        } else if datatype == TUINT {
            let array = cast_slice(array);
            ffpcnuk_safe(
                fptr,
                colnum,
                firstrow,
                firstelem,
                nelem,
                array,
                match nulval {
                    NullValue::UInt(x) => x,
                    _ => 0,
                },
                status,
            );
        } else if datatype == TINT {
            let array = cast_slice(array);
            ffpcnk_safe(
                fptr,
                colnum,
                firstrow,
                firstelem,
                nelem,
                array,
                match nulval {
                    NullValue::Int(x) => x,
                    _ => 0,
                },
                status,
            );
        } else if datatype == TULONG {
            let array = cast_slice(array);
            ffpcnuj_safe(
                fptr,
                colnum,
                firstrow,
                firstelem,
                nelem,
                array,
                match nulval {
                    NullValue::ULong(x) => x,
                    _ => 0,
                },
                status,
            );
        } else if datatype == TLONG {
            let array = cast_slice(array);
            ffpcnj_safe(
                fptr,
                colnum,
                firstrow,
                firstelem,
                nelem,
                array,
                match nulval {
                    NullValue::Long(x) => x,
                    _ => 0,
                },
                status,
            );
        } else if datatype == TULONGLONG {
            let array = cast_slice(array);
            ffpcnujj_safe(
                fptr,
                colnum,
                firstrow,
                firstelem,
                nelem,
                array,
                match nulval {
                    NullValue::ULONGLONG(x) => x,
                    _ => 0,
                },
                status,
            );
        } else if datatype == TLONGLONG {
            let array = cast_slice(array);
            ffpcnjj_safe(
                fptr,
                colnum,
                firstrow,
                firstelem,
                nelem,
                array,
                match nulval {
                    NullValue::LONGLONG(x) => x,
                    _ => 0,
                },
                status,
            );
        } else if datatype == TFLOAT {
            let array = cast_slice(array);
            ffpcne_safe(
                fptr,
                colnum,
                firstrow,
                firstelem,
                nelem,
                array,
                match nulval {
                    NullValue::Float(x) => x,
                    _ => 0.0,
                },
                status,
            );
        } else if datatype == TDOUBLE {
            let array = cast_slice(array);
            ffpcnd_safe(
                fptr,
                colnum,
                firstrow,
                firstelem,
                nelem,
                array,
                match nulval {
                    NullValue::Double(x) => x,
                    _ => 0.0,
                },
                status,
            );
        } else if datatype == TCOMPLEX {
            let array = cast_slice(array);
            ffpcne_safe(
                fptr,
                colnum,
                firstrow,
                (firstelem - 1) * 2 + 1,
                nelem * 2,
                array,
                match nulval {
                    NullValue::Float(x) => x,
                    _ => 0.0,
                },
                status,
            );
        } else if datatype == TDBLCOMPLEX {
            let array = cast_slice(array);
            ffpcnd_safe(
                fptr,
                colnum,
                firstrow,
                (firstelem - 1) * 2 + 1,
                nelem * 2,
                array,
                match nulval {
                    NullValue::Double(x) => x,
                    _ => 0.0,
                },
                status,
            );
        } else if datatype == TLOGICAL {
            let array = cast_slice(array);
            ffpcnl_safe(
                fptr,
                colnum,
                firstrow,
                firstelem,
                nelem,
                array,
                match nulval {
                    NullValue::Byte(x) => x as c_char,
                    NullValue::UByte(x) => x as c_char,
                    _ => 0,
                },
                status,
            );
        } else if datatype == TSTRING {
            let array: &[*const c_char] =
                slice::from_raw_parts(array.as_ptr() as *const _, nelem as usize);
            let mut v_array = Vec::new();
            for item in array {
                let array_item = slice::from_raw_parts(*item, FLEN_VALUE);
                v_array.push(array_item);
            }

            match nulval {
                // safe to unwrap since we checked for None above
                NullValue::String(x) => {
                    ffpcns_safe(
                        fptr,
                        colnum,
                        firstrow,
                        firstelem,
                        nelem,
                        &v_array,
                        cast_slice(x.as_bytes_with_nul()),
                        status,
                    );
                }
                _ => {
                    ffpcns_safe(
                        fptr,
                        colnum,
                        firstrow,
                        firstelem,
                        nelem,
                        &v_array,
                        &[0],
                        status,
                    );
                }
            }
        } else {
            *status = BAD_DATATYPE;
        }

        *status
    }
}

/*--------------------------------------------------------------------------*/
/// Write arrays of values to NCOLS table columns. This is an optimization
/// to write all columns in one pass through the table.  The datatypes of the
/// input arrays are defined by the 3rd argument.  Data conversion
/// and scaling will be performed if necessary (e.g, if the datatype of
/// the FITS array is not the same as the array being written).
/// Undefined elements for column i that are equal to *(nulval[i]) are set to
/// the defined null value, unless nulval[i]=0,
/// in which case no checking for undefined values will be performed.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffpcln(
    fptr: *mut fitsfile,          /* I - FITS file pointer                       */
    ncols: c_int,                 /* I - number of columns to write              */
    datatype: *const c_int,       /* I - datatypes of the values                 */
    colnum: *const c_int,         /* I - columns numbers to write (1 = 1st col)  */
    firstrow: LONGLONG,           /* I - first row to write (1 = 1st row)    */
    nrows: LONGLONG,              /* I - number of rows to write             */
    array: *const *const c_void,  /* I - array of pointers to values to write    */
    nulval: *const *const c_void, /* I - array of pointers to values for undefined pixels */
    status: *mut c_int,           /* IO - error status                           */
) -> c_int {
    todo!();
}

/*--------------------------------------------------------------------------*/
/// set all the parameters for an iterator column, by column name
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn fits_iter_set_by_name(
    col: *mut iteratorCol,  /* I - iterator col structure */
    fptr: *mut fitsfile,    /* I - FITS file pointer                      */
    colname: *const c_char, /* I - column name                            */
    datatype: c_int,        /* I - column datatype                        */
    iotype: c_int,          /* I - InputCol, InputOutputCol, or OutputCol */
) -> c_int {
    unsafe {
        let col = col.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        raw_to_slice!(colname);

        fits_iter_set_by_name_safe(col, fptr, colname, datatype, iotype)
    }
}

/*--------------------------------------------------------------------------*/
/// set all the parameters for an iterator column, by column name
pub fn fits_iter_set_by_name_safe(
    col: &mut iteratorCol, /* I - iterator col structure */
    fptr: &mut fitsfile,   /* I - FITS file pointer                      */
    colname: &[c_char],    /* I - column name                            */
    datatype: c_int,       /* I - column datatype                        */
    iotype: c_int,         /* I - InputCol, InputOutputCol, or OutputCol */
) -> c_int {
    col.fptr = fptr;
    strncpy_safe(&mut col.colname, colname, 69);
    col.colname[69] = 0;
    col.colnum = 0; /* set column number undefined since name is given */
    col.datatype = datatype;
    col.iotype = iotype;

    0
}

/*--------------------------------------------------------------------------*/
/// set all the parameters for an iterator column, by column number
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn fits_iter_set_by_num(
    col: *mut iteratorCol, /* I - iterator col structure */
    fptr: *mut fitsfile,   /* I - FITS file pointer                      */
    colnum: c_int,         /* I - column number                          */
    datatype: c_int,       /* I - column datatype                        */
    iotype: c_int,         /* I - InputCol, InputOutputCol, or OutputCol */
) -> c_int {
    unsafe {
        let col = col.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        fits_iter_set_by_num_safe(col, fptr, colnum, datatype, iotype)
    }
}

/*--------------------------------------------------------------------------*/
/// set all the parameters for an iterator column, by column number
pub fn fits_iter_set_by_num_safe(
    col: &mut iteratorCol, /* I - iterator col structure */
    fptr: &mut fitsfile,   /* I - FITS file pointer                      */
    colnum: c_int,         /* I - column number                          */
    datatype: c_int,       /* I - column datatype                        */
    iotype: c_int,         /* I - InputCol, InputOutputCol, or OutputCol */
) -> c_int {
    col.fptr = fptr;
    col.colnum = colnum;
    col.datatype = datatype;
    col.iotype = iotype;

    0
}

/*--------------------------------------------------------------------------*/
/// set iterator column parameter
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn fits_iter_set_file(
    col: *mut iteratorCol, /* I - iterator column structure   */
    fptr: *mut fitsfile,   /* I - FITS file pointer                      */
) -> c_int {
    unsafe {
        let col = col.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);

        fits_iter_set_file_safe(col, fptr)
    }
}

/*--------------------------------------------------------------------------*/
/// set iterator column parameter
pub fn fits_iter_set_file_safe(
    col: &mut iteratorCol, /* I - iterator column structure   */
    fptr: &mut fitsfile,   /* I - FITS file pointer                      */
) -> c_int {
    col.fptr = fptr;
    0
}

/*--------------------------------------------------------------------------*/
/// set iterator column parameter
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn fits_iter_set_colname(
    col: *mut iteratorCol,  /* I - iterator col structure  */
    colname: *const c_char, /* I - column name                            */
) -> c_int {
    unsafe {
        let col = col.as_mut().expect(NULL_MSG);
        raw_to_slice!(colname);
        fits_iter_set_colname_safe(col, colname)
    }
}

/*--------------------------------------------------------------------------*/
/// set iterator column parameter
pub fn fits_iter_set_colname_safe(
    col: &mut iteratorCol, /* I - iterator col structure  */
    colname: &[c_char],    /* I - column name                            */
) -> c_int {
    strncpy_safe(&mut col.colname, colname, 69);
    col.colname[69] = 0;
    col.colnum = 0; /* set column number undefined since name is given */
    0
}

/*--------------------------------------------------------------------------*/
/// set iterator column parameter
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn fits_iter_set_colnum(
    col: *mut iteratorCol, /* I - iterator column structure */
    colnum: c_int,         /* I - column number                          */
) -> c_int {
    unsafe {
        let col = col.as_mut().expect(NULL_MSG);
        fits_iter_set_colnum_safe(col, colnum)
    }
}

/*--------------------------------------------------------------------------*/
/// set iterator column parameter
pub fn fits_iter_set_colnum_safe(
    col: &mut iteratorCol, /* I - iterator column structure */
    colnum: c_int,         /* I - column number                          */
) -> c_int {
    col.colnum = colnum;
    0
}

/*--------------------------------------------------------------------------*/
/// set iterator column parameter
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn fits_iter_set_datatype(
    col: *mut iteratorCol, /* I - iterator col structure */
    datatype: c_int,       /* I - column datatype                        */
) -> c_int {
    unsafe {
        let col = col.as_mut().expect(NULL_MSG);
        fits_iter_set_datatype_safe(col, datatype)
    }
}

/*--------------------------------------------------------------------------*/
/// set iterator column parameter
pub fn fits_iter_set_datatype_safe(
    col: &mut iteratorCol, /* I - iterator col structure */
    datatype: c_int,       /* I - column datatype                        */
) -> c_int {
    col.datatype = datatype;
    0
}

/*--------------------------------------------------------------------------*/
/// set iterator column parameter
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn fits_iter_set_iotype(
    col: *mut iteratorCol, /* I - iterator column structure */
    iotype: c_int,         /* I - InputCol, InputOutputCol, or OutputCol */
) -> c_int {
    unsafe {
        let col = col.as_mut().expect(NULL_MSG);
        fits_iter_set_iotype_safe(col, iotype)
    }
}

/*--------------------------------------------------------------------------*/
/// set iterator column parameter
pub fn fits_iter_set_iotype_safe(
    col: &mut iteratorCol, /* I - iterator column structure */
    iotype: c_int,         /* I - InputCol, InputOutputCol, or OutputCol */
) -> c_int {
    col.iotype = iotype;
    0
}

/*--------------------------------------------------------------------------*/
/// get iterator column parameter
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn fits_iter_get_file(
    col: *mut iteratorCol, /* I -iterator col structure */
) -> *mut fitsfile {
    unsafe {
        let col = col.as_mut().expect(NULL_MSG);
        fits_iter_get_file_safe(col)
    }
}

/*--------------------------------------------------------------------------*/
/// get iterator column parameter
pub fn fits_iter_get_file_safe(
    col: &mut iteratorCol, /* I -iterator col structure */
) -> *mut fitsfile {
    col.fptr
}

/*--------------------------------------------------------------------------*/
/// get iterator column parameter
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn fits_iter_get_colname(
    col: *mut iteratorCol, /* I -iterator col structure */
) -> *const c_char {
    unsafe {
        let col = col.as_mut().expect(NULL_MSG);
        fits_iter_get_colname_safe(col).as_ptr()
    }
}

/*--------------------------------------------------------------------------*/
/// get iterator column parameter
pub fn fits_iter_get_colname_safe(
    col: &mut iteratorCol, /* I -iterator col structure */
) -> &[c_char; 70] {
    &col.colname
}

/*--------------------------------------------------------------------------*/
/// get iterator column parameter
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn fits_iter_get_colnum(
    col: *mut iteratorCol, /* I - iterator column structure */
) -> c_int {
    unsafe {
        let col = col.as_mut().expect(NULL_MSG);
        fits_iter_get_colnum_safe(col)
    }
}

/*--------------------------------------------------------------------------*/
/// get iterator column parameter
pub fn fits_iter_get_colnum_safe(
    col: &mut iteratorCol, /* I - iterator column structure */
) -> c_int {
    col.colnum
}

/*--------------------------------------------------------------------------*/
/// get iterator column parameter
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn fits_iter_get_datatype(
    col: *mut iteratorCol, /* I - iterator col structure */
) -> c_int {
    unsafe {
        let col = col.as_mut().expect(NULL_MSG);
        fits_iter_get_datatype_safe(col)
    }
}

/*--------------------------------------------------------------------------*/
/// get iterator column parameter
pub fn fits_iter_get_datatype_safe(
    col: &mut iteratorCol, /* I - iterator col structure */
) -> c_int {
    col.datatype
}

/*--------------------------------------------------------------------------*/
/// get iterator column parameter
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn fits_iter_get_iotype(
    col: *mut iteratorCol, /* I - iterator column structure */
) -> c_int {
    unsafe {
        let col = col.as_mut().expect(NULL_MSG);
        fits_iter_get_iotype_safe(col)
    }
}

/*--------------------------------------------------------------------------*/
/// get iterator column parameter
pub fn fits_iter_get_iotype_safe(
    col: &mut iteratorCol, /* I - iterator column structure */
) -> c_int {
    col.iotype
}

/*--------------------------------------------------------------------------*/
/// get iterator column parameter
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn fits_iter_get_array(
    col: *mut iteratorCol, /* I - iterator col structure */
) -> *const c_void {
    unsafe {
        let col = col.as_mut().expect(NULL_MSG);
        fits_iter_get_array_safe(col)
    }
}

/*--------------------------------------------------------------------------*/
/// get iterator column parameter
pub fn fits_iter_get_array_safe(
    col: &mut iteratorCol, /* I - iterator col structure */
) -> *const c_void {
    col.array
}

/*--------------------------------------------------------------------------*/
/// get iterator column parameter
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn fits_iter_get_tlmin(
    col: *mut iteratorCol, /* I - iterator column structure */
) -> c_long {
    unsafe {
        let col = col.as_mut().expect(NULL_MSG);
        fits_iter_get_tlmin_safe(col)
    }
}

/*--------------------------------------------------------------------------*/
/// get iterator column parameter
pub fn fits_iter_get_tlmin_safe(
    col: &mut iteratorCol, /* I - iterator column structure */
) -> c_long {
    col.tlmin
}

/*--------------------------------------------------------------------------*/
/// get iterator column parameter
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn fits_iter_get_tlmax(
    col: *mut iteratorCol, /* I - iterator column structure */
) -> c_long {
    unsafe {
        let col = col.as_mut().expect(NULL_MSG);
        fits_iter_get_tlmax_safe(col)
    }
}

/*--------------------------------------------------------------------------*/
/// get iterator column parameter
pub fn fits_iter_get_tlmax_safe(
    col: &mut iteratorCol, /* I - iterator column structure */
) -> c_long {
    col.tlmax
}

/*--------------------------------------------------------------------------*/
/// get iterator column parameter
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn fits_iter_get_repeat(
    col: *mut iteratorCol, /* I - iterator col structure */
) -> c_long {
    unsafe {
        let col = col.as_mut().expect(NULL_MSG);
        fits_iter_get_repeat_safe(col)
    }
}

/*--------------------------------------------------------------------------*/
/// get iterator column parameter
pub fn fits_iter_get_repeat_safe(
    col: &mut iteratorCol, /* I - iterator col structure */
) -> c_long {
    col.repeat
}

/*--------------------------------------------------------------------------*/
/// get iterator column parameter
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn fits_iter_get_tunit(
    col: *mut iteratorCol, /* I - iterator col structure */
) -> *const c_char {
    unsafe {
        let col = col.as_mut().expect(NULL_MSG);
        fits_iter_get_tunit_safe(col).as_ptr()
    }
}

/*--------------------------------------------------------------------------*/
/// get iterator column parameter
pub fn fits_iter_get_tunit_safe(
    col: &mut iteratorCol, /* I - iterator col structure */
) -> &[c_char; 70] {
    &col.tunit
}

/*--------------------------------------------------------------------------*/
/// get iterator column parameter
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn fits_iter_get_tdisp(
    col: *mut iteratorCol, /* I -iterator col structure   */
) -> *const c_char {
    unsafe {
        let col = col.as_mut().expect(NULL_MSG);
        fits_iter_get_tdisp_safe(col).as_ptr()
    }
}

/*--------------------------------------------------------------------------*/
/// get iterator column parameter
pub fn fits_iter_get_tdisp_safe(
    col: &mut iteratorCol, /* I -iterator col structure   */
) -> &[c_char; 70] {
    &col.tdisp
}

/*--------------------------------------------------------------------------*/
/// The iterator function.  This function will pass the specified
/// columns from a FITS table or pixels from a FITS image to the
/// user-supplied function.  Depending on the size of the table
/// or image, only a subset of the rows or pixels may be passed to the
/// function on each call, in which case the function will be called
/// multiple times until all the rows or pixels have been processed.
#[cfg_attr(not(test), unsafe(no_mangle), deprecated)]
pub unsafe extern "C" fn ffiter(
    n_cols: c_int,
    cols: *mut iteratorCol,
    offset: c_long,
    n_per_loop: c_long,
    workfn: extern "C" fn(
        total_n: c_long,
        offset: c_long,
        first_n: c_long,
        n_values: c_long,
        n_cols: c_int,
        cols: *mut iteratorCol,
        userPointer: *mut c_void,
    ) -> c_int,
    userPointer: *mut c_void,
    status: *mut c_int,
) -> c_int {
    todo!()
}
