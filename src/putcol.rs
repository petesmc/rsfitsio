/*  This file, putcol.rs, contains routines that write data elements to     */
/*  a FITS image or table. These are the generic routines.                 */

/*  The FITSIO software was written by William Pence at the High Energy    */
/*  Astrophysic Science Archive Research Center (HEASARC) at the NASA      */
/*  Goddard Space Flight Center.                                           */

use core::ffi::CStr;
use core::slice;

use crate::buffers::ffgrsz_safe;
use crate::c_types::*;
use crate::cs;
use crate::fitscore::ffc2ii;
use crate::fitscore::ffgcno_safe;
use crate::fitscore::ffgdes_safe;
use crate::fitscore::ffghdt_safe;
use crate::fitscore::ffgidm_safe;
use crate::fitscore::ffgnrwll_safe;
use crate::fitscore::ffgtcl_safe;
use crate::fitscore::ffgtclll_safe;
use crate::fitscore::ffkeyn_safe;
use crate::fitscore::ffpmsg_slice;
use crate::fitscore::ffpmsg_str;
use crate::fitsio2::INT32_MAX;
use crate::getcol::ffgcv_safe;
use crate::getcol::ffgpv_safe;
use crate::getkey::ffgkyj_safe;
use crate::getkey::ffgkys_safe;

use bytemuck::{cast_slice, cast_slice_mut};
use core::mem::size_of;
use libc::{calloc, free, memcmp, memcpy, memset};
use std::os::raw::{c_char, c_int, c_long, c_schar, c_short, c_uchar, c_uint, c_ulong, c_ushort};

use crate::NullValue;
use crate::aliases::rust_api::*;
use crate::calculate_subsection_length_unit;
use crate::fitscore::{ffgisz_safe, ffgiszll_safe};
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
use crate::{bytes_per_datatype, fitsio::*, int_snprintf};

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
    // FFI WRAPPER
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

        let array = slice::from_raw_parts(array.cast::<u8>(), bytes);

        let mut naxis: c_int = 0;
        ffgidm_safe(fptr, &mut naxis, status);

        let firstpix = slice::from_raw_parts(firstpix, naxis as usize);

        ffppx_safe(fptr, datatype, firstpix, nelem, array, status)
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
pub fn ffppx_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer                       */
    datatype: c_int,     /* I - datatype of the value                   */
    firstpix: &[c_long], /* I - coord of  first pixel to write(1 based) */
    nelem: LONGLONG,     /* I - number of values to write               */
    array: &[u8],        /* I - array of values that are written        */
    status: &mut c_int,  /* IO - error status                           */
) -> c_int {
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
    // FFI WRAPPER
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

        let array = slice::from_raw_parts(array.cast::<u8>(), bytes);

        let mut naxis: c_int = 0;
        ffgidm_safe(fptr, &mut naxis, status);
        let firstpix = slice::from_raw_parts(firstpix, naxis as usize);

        ffppxll_safe(fptr, datatype, firstpix, nelem, array, status)
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
pub fn ffppxll_safe(
    fptr: &mut fitsfile,   /* I - FITS file pointer                       */
    datatype: c_int,       /* I - datatype of the value                   */
    firstpix: &[LONGLONG], /* I - coord of  first pixel to write(1 based) */
    nelem: LONGLONG,       /* I - number of values to write               */
    array: &[u8],          /* I - array of values that are written        */
    status: &mut c_int,    /* IO - error status                           */
) -> c_int {
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
    // FFI WRAPPER
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

        let array = slice::from_raw_parts(array.cast::<u8>(), bytes);

        let mut naxis = 0;
        ffgidm_safe(fptr, &mut naxis, status);
        let firstpix = slice::from_raw_parts(firstpix, naxis as usize);

        let nulval = NullValue::from_raw_ptr(datatype, nulval);

        ffppxn_safe(fptr, datatype, firstpix, nelem, array, nulval, status)
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
pub fn ffppxn_safe(
    fptr: &mut fitsfile,       /* I - FITS file pointer                       */
    datatype: c_int,           /* I - datatype of the value                   */
    firstpix: &[c_long],       /* I - first vector element to write(1 = 1st)  */
    nelem: LONGLONG,           /* I - number of values to write               */
    array: &[u8],              /* I - array of values that are written        */
    nulval: Option<NullValue>, /* I - pointer to the null value               */
    status: &mut c_int,        /* IO - error status                           */
) -> c_int {
    let mut naxis = 0;

    let group: c_long = 1;

    let mut firstelem: LONGLONG = 0;
    let mut dimsize: LONGLONG = 1;
    let mut naxes: [LONGLONG; 9] = [0; 9];

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    if nulval.is_none() {
        /* null value not defined? */
        ffppx_safe(fptr, datatype, firstpix, nelem, array, status);
        return *status;
    }

    let nulval = nulval.unwrap();

    /* get the size of the image */
    ffgidm_safe(fptr, &mut naxis, status);
    ffgiszll_safe(fptr, 9, &mut naxes, status);

    firstelem = 0;

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
            nulval.get_value_as_f64() as u8,
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
            nulval.get_value_as_f64() as i8,
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
            nulval.get_value_as_f64() as c_ushort,
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
            nulval.get_value_as_f64() as c_short,
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
            nulval.get_value_as_f64() as c_uint,
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
            nulval.get_value_as_f64() as c_int,
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
            nulval.get_value_as_f64() as c_ulong,
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
            nulval.get_value_as_f64() as c_long,
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
            nulval.get_value_as_f64() as ULONGLONG,
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
            nulval.get_value_as_f64() as LONGLONG,
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
            nulval.get_value_as_f64() as f32,
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
            nulval.get_value_as_f64(),
            status,
        );
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
    // FFI WRAPPER
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

        let array = slice::from_raw_parts(array.cast::<u8>(), bytes);
        let mut naxis = 0;
        ffgidm_safe(fptr, &mut naxis, status);
        let firstpix = slice::from_raw_parts(firstpix, naxis as usize);

        let nulval = NullValue::from_raw_ptr(datatype, nulval);

        ffppxnll_safe(fptr, datatype, firstpix, nelem, array, nulval, status)
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
pub fn ffppxnll_safe(
    fptr: &mut fitsfile,       /* I - FITS file pointer                       */
    datatype: c_int,           /* I - datatype of the value                   */
    firstpix: &[LONGLONG],     /* I - first vector element to write(1 = 1st)  */
    nelem: LONGLONG,           /* I - number of values to write               */
    array: &[u8],              /* I - array of values that are written        */
    nulval: Option<NullValue>, /* I - pointer to the null value               */
    status: &mut c_int,        /* IO - error status                           */
) -> c_int {
    let mut naxis = 0;

    let group: c_long = 1;

    let mut firstelem: LONGLONG = 0;
    let mut dimsize: LONGLONG = 1;
    let mut naxes: [LONGLONG; 9] = [0; 9];

    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    if nulval.is_none() {
        /* null value not defined? */
        ffppxll_safe(fptr, datatype, firstpix, nelem, array, status);
        return *status;
    }

    let nulval = nulval.unwrap();

    /* get the size of the image */
    ffgidm_safe(fptr, &mut naxis, status);
    ffgiszll_safe(fptr, 9, &mut naxes, status);

    firstelem = 0;

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
            nulval.get_value_as_f64() as u8,
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
            nulval.get_value_as_f64() as i8,
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
            nulval.get_value_as_f64() as c_ushort,
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
            nulval.get_value_as_f64() as c_short,
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
            nulval.get_value_as_f64() as c_uint,
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
            nulval.get_value_as_f64() as c_int,
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
            nulval.get_value_as_f64() as c_ulong,
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
            nulval.get_value_as_f64() as c_long,
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
            nulval.get_value_as_f64() as ULONGLONG,
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
            nulval.get_value_as_f64() as LONGLONG,
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
            nulval.get_value_as_f64() as f32,
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
            nulval.get_value_as_f64(),
            status,
        );
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
pub unsafe extern "C" fn ffppr(
    fptr: *mut fitsfile,  /* I - FITS file pointer                       */
    datatype: c_int,      /* I - datatype of the value                   */
    firstelem: LONGLONG,  /* I - first vector element to write(1 = 1st)  */
    nelem: LONGLONG,      /* I - number of values to write               */
    array: *const c_void, /* I - array of values that are written        */
    status: *mut c_int,   /* IO - error status                           */
) -> c_int {
    // FFI WRAPPER
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

        let array = slice::from_raw_parts(array.cast::<u8>(), bytes);

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
    // FFI WRAPPER
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

        let array = slice::from_raw_parts(array.cast::<u8>(), bytes);
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
    // FFI WRAPPER
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

        let array = slice::from_raw_parts(array.cast::<u8>(), bytes);

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

    let naxis = c_long::from(naxis);

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
    // FFI WRAPPER
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

        let array = slice::from_raw_parts(array.cast::<u8>(), bytes);

        ffpcl_safe(
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
pub fn ffpcl_safe(
    fptr: &mut fitsfile, /* I - FITS file pointer                       */
    datatype: c_int,     /* I - datatype of the value                   */
    colnum: c_int,       /* I - number of column to write (1 = 1st col) */
    firstrow: LONGLONG,  /* I - first row to write (1 = 1st row)        */
    firstelem: LONGLONG, /* I - first vector element to write (1 = 1st) */
    nelem: LONGLONG,     /* I - number of elements to write             */
    array: &[u8],        /* I - array of values that are written        */
    status: &mut c_int,  /* IO - error status                           */
) -> c_int {
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
        // SAFETY: for TSTRING the caller's buffer holds an array of `char *`,
        // as in the C, so it carries pointer alignment even though the
        // signature types it as bytes (clippy's cast_ptr_alignment); assert
        // that rather than assuming it.
        TSTRING => unsafe {
            #[allow(clippy::cast_ptr_alignment)]
            let p = array.as_ptr().cast::<*const c_char>();
            debug_assert!(p.is_aligned());
            let array = slice::from_raw_parts(p, nelem as usize);
            let mut v_array = Vec::new();
            for item in array {
                let array_item = cast_slice(CStr::from_ptr(*item).to_bytes_with_nul());
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
        },
        _ => {
            *status = BAD_DATATYPE;
        }
    }
    *status
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
    // FFI WRAPPER
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

        let array: &[u8] = slice::from_raw_parts(array.cast::<u8>(), bytes);
        let nulval = NullValue::from_raw_ptr(datatype, nulval);

        ffpcn_safe(
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
pub fn ffpcn_safe(
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
    if *status > 0 {
        /* inherit input status value if > 0 */
        return *status;
    }

    if nulval.is_none() {
        /* null value not defined? */
        ffpcl_safe(
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
        unsafe {
            let array: &[*const c_char] =
                slice::from_raw_parts(array.as_ptr().cast(), nelem as usize);
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
        }
    } else {
        *status = BAD_DATATYPE;
    }

    *status
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
    // FFI WRAPPER
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let fptr = fptr.as_mut().expect(NULL_MSG);
        let datatype = slice::from_raw_parts(datatype, ncols as usize);
        let colnum = slice::from_raw_parts(colnum, ncols as usize);
        let array = slice::from_raw_parts(array, ncols as usize);
        let nulval = slice::from_raw_parts(nulval, ncols as usize);

        ffpcln_safe(
            fptr, ncols, datatype, colnum, firstrow, nrows, array, nulval, status,
        )
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
pub fn ffpcln_safe(
    fptr: &mut fitsfile,      /* I - FITS file pointer                       */
    ncols: c_int,             /* I - number of columns to write              */
    datatype: &[c_int],       /* I - datatypes of the values                 */
    colnum: &[c_int],         /* I - columns numbers to write (1 = 1st col)  */
    firstrow: LONGLONG,       /* I - first row to write (1 = 1st row)    */
    nrows: LONGLONG,          /* I - number of rows to write             */
    array: &[*const c_void],  /* I - array of pointers to values to write    */
    nulval: &[*const c_void], /* I - array of pointers to values for undefined pixels */
    status: &mut c_int,       /* IO - error status                           */
) -> c_int {
    let mut ntotrows: LONGLONG = 0;
    let mut nrowbuf: c_long = 0;
    let mut sizes: [usize; 255] = [0; 255];

    sizes[TBYTE as usize] = size_of::<c_char>();
    sizes[TSBYTE as usize] = size_of::<c_char>();
    sizes[TLOGICAL as usize] = size_of::<c_char>();
    sizes[TUSHORT as usize] = size_of::<c_short>();
    sizes[TSHORT as usize] = size_of::<c_short>();
    sizes[TINT as usize] = size_of::<c_int>();
    sizes[TUINT as usize] = size_of::<c_int>();
    sizes[TLONG as usize] = size_of::<c_long>();
    sizes[TULONG as usize] = size_of::<c_long>();
    sizes[TLONGLONG as usize] = size_of::<LONGLONG>();
    sizes[TULONGLONG as usize] = size_of::<LONGLONG>();
    sizes[TFLOAT as usize] = size_of::<f32>();
    sizes[TDOUBLE as usize] = size_of::<f64>();
    sizes[TDBLCOMPLEX as usize] = 2 * size_of::<f64>();

    if *status > 0 {
        return *status;
    }

    if ncols <= 0 {
        *status = 0;
        return *status;
    }

    /* C used malloc(sizeof(LONGLONG)*ncols); a Vec owns the storage and frees it
    on return, so the explicit MEMORY_ALLOCATION failure path and free() calls
    are dropped. */
    let mut repeats: Vec<LONGLONG> = vec![0; ncols as usize];

    ffgnrwll_safe(fptr, &mut ntotrows, status);
    ffgrsz_safe(fptr, &mut nrowbuf, status);

    /* Retrieve column repeats */
    let mut icol: c_int = 0;
    while icol < ncols && icol < 1000 {
        let mut typecode: c_int = 0;
        let mut repeat: LONGLONG = 0;
        let mut width: LONGLONG = 0;
        ffgtclll_safe(
            fptr,
            colnum[icol as usize],
            Some(&mut typecode),
            Some(&mut repeat),
            Some(&mut width),
            status,
        );
        repeats[icol as usize] = repeat;

        let dt = datatype[icol as usize];
        if dt == TBIT || dt == TSTRING || sizes.get(dt as usize).copied().unwrap_or(0) == 0 {
            ffpmsg_str("Cannot write to TBIT or TSTRING datatypes (ffpcln)");
            *status = BAD_DATATYPE;
        }
        if typecode < 0 {
            ffpmsg_str("Cannot write to variable-length data (ffpcln)");
            *status = BAD_DIMEN;
        }

        if *status != 0 {
            break;
        }
        icol += 1;
    }
    if *status != 0 {
        return *status;
    }

    /* Optimize for 1 column */
    if ncols == 1 {
        let dt = datatype[0];
        let elsize = sizes.get(dt as usize).copied().unwrap_or(0);
        let nelem = nrows * repeats[0];
        // SAFETY: array[0] points to a buffer of at least nelem*elsize bytes
        // supplied by the caller; ffpcn_safe expects exactly that many bytes.
        let array0: &[u8] =
            unsafe { slice::from_raw_parts(array[0].cast::<u8>(), nelem as usize * elsize) };
        let nulval0 = NullValue::from_raw_ptr(dt, nulval[0]);
        ffpcn_safe(
            fptr, dt, colnum[0], firstrow, 1, nelem, array0, nulval0, status,
        );
        return *status;
    }

    /* Scan through file, in chunks of nrowbuf */
    let mut currow = firstrow;
    let mut ndone: LONGLONG = 0;
    while ndone < nrows {
        let mut nwrite = nrows - ndone;
        if nwrite > nrowbuf as LONGLONG {
            nwrite = nrowbuf as LONGLONG;
        }

        let mut icol: c_int = 0;
        while icol < ncols {
            let ic = icol as usize;
            let dt = datatype[ic];
            let elsize = sizes.get(dt as usize).copied().unwrap_or(0);
            let nelem1: LONGLONG = nwrite * repeats[ic];
            let byteoff = (repeats[ic] * ndone) as usize * elsize;
            // SAFETY: array[ic] points to a caller-supplied buffer; the offset
            // and length mirror the C pointer arithmetic and stay within it.
            let array1: &[u8] = unsafe {
                slice::from_raw_parts(
                    array[ic].cast::<u8>().add(byteoff),
                    nelem1 as usize * elsize,
                )
            };
            let nulval1 = NullValue::from_raw_ptr(dt, nulval[ic]);

            ffpcn_safe(
                fptr,
                dt,
                colnum[ic],
                ndone + 1,
                1,
                nelem1,
                array1,
                nulval1,
                status,
            );
            if *status != 0 {
                let mut errmsg: [c_char; 100] = [0; 100];
                int_snprintf!(
                    &mut errmsg,
                    100,
                    "Failed to write column {} data rows {}-{} (ffpcln)",
                    colnum[ic],
                    currow,
                    currow + nwrite - 1
                );
                ffpmsg_slice(&errmsg);
                break;
            }
            icol += 1;
        }

        if *status != 0 {
            break;
        }
        currow += nwrite;
        ndone += nwrite;
    }

    *status
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
    // FFI WRAPPER
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
    // FFI WRAPPER
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
    fptr: *mut fitsfile,   /* I - FITS file pointer                      */
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
    // FFI WRAPPER
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
    // FFI WRAPPER
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
    // FFI WRAPPER
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
    // FFI WRAPPER
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
    // FFI WRAPPER
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
    // FFI WRAPPER
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
    // FFI WRAPPER
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
    // FFI WRAPPER
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
    // FFI WRAPPER
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
    // FFI WRAPPER
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
    // FFI WRAPPER
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
    // FFI WRAPPER
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
    // FFI WRAPPER
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
    // FFI WRAPPER
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
    // FFI WRAPPER
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
    // FFI WRAPPER
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
    // FFI WRAPPER
    unsafe {
        let status = status.as_mut().expect(NULL_MSG);
        let cols = slice::from_raw_parts_mut(cols, n_cols as usize);

        ffiter_safe(
            n_cols,
            cols,
            offset,
            n_per_loop,
            workfn,
            userPointer,
            status,
        )
    }
}

/*--------------------------------------------------------------------------*/
/// Iterate through the specified columns or rows of data in a table
/// or image, calling a user-supplied function for each group of rows.
/// The columns are passed by reference as iteratorCol structs. This
/// function provides a way to efficiently process large amounts of data
/// in a table or image by reading and processing it in chunks rather
/// than all at once.
pub fn ffiter_safe(
    n_cols: c_int,
    cols: &mut [iteratorCol],
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
    status: &mut c_int,
) -> c_int {
    let cols_len = cols.len();
    ffiter_internal(
        n_cols,
        cols.as_mut_ptr(),
        cols_len,
        offset,
        n_per_loop,
        workfn,
        userPointer,
        status,
    )
}

/*--------------------------------------------------------------------------*/
/// The body of [`ffiter_safe`], with the column array passed as a raw pointer
/// instead of as a `&mut [iteratorCol]`.
///
/// The work function is free to reach the very same `iteratorCol` array by a
/// second route -- the parser's work function does exactly that, because
/// `ParseData::colData` *is* the array handed to the iterator. A `&mut` slice
/// parameter cannot express that: it is exclusive for the whole call, so the
/// work function's own access to the array would invalidate it.
///
/// `cols_ptr` is therefore kept as a raw pointer for the lifetime of the
/// iteration and re-borrowed as a slice only between work-function calls; the
/// pointer handed to the work function is `cols_ptr` itself, not a pointer
/// derived from the local re-borrow, so that the two routes to the array are
/// siblings rather than parent and child. Every work-function call is followed
/// by a fresh re-borrow, because the call may have written to the array through
/// the other route.
///
/// # Safety
/// `cols_ptr` must point at `cols_len` initialised `iteratorCol`s that stay put
/// (and stay allocated) for the whole call, and no `&mut` derived from it may
/// be live in the caller across the call. `cols_len` must be at least
/// `max(n_cols, 1)`.
#[allow(clippy::too_many_arguments)]
pub(crate) fn ffiter_internal(
    n_cols: c_int,
    cols_ptr: *mut iteratorCol,
    cols_len: usize,
    mut offset: c_long,
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
    status: &mut c_int,
) -> c_int {
    /* Re-borrowed from `cols_ptr` again after every work-function call. */
    let mut cols: &mut [iteratorCol] = unsafe { slice::from_raw_parts_mut(cols_ptr, cols_len) };

    // Structure to store the column null value
    #[derive(Default, Clone, Copy)]
    struct ColNulls {
        nullsize: usize, // length of the null value, in bytes
        null: ColNullValue,
    }

    #[derive(Default, Clone, Copy)]
    enum ColNullValue {
        StringNull(*mut c_char),
        UCharNull(c_uchar),
        SCharNull(c_schar),
        IntNull(c_int),
        ShortNull(c_short),
        LongNull(c_long),
        UIntNull(c_uint),
        UShortNull(c_ushort),
        ULongNull(c_ulong),
        FloatNull(f32),
        DoubleNull(f64),
        LONGLONGNull(LONGLONG),
        #[default]
        None,
    }

    impl ColNullValue {
        /// The null value as the C `union` held it: the bytes of the value
        /// itself, in the column's own type.
        ///
        /// The C passes `col[jj].null` -- a `char[8]` union -- straight to the
        /// read routines as the `void*` null value, and copies it into the
        /// first array element for the caller. A Rust enum cannot stand in for
        /// that union: its first bytes are the discriminant, and where the
        /// payload lands is not something the layout promises. Casting a
        /// `*const ColNullValue` to `*const f64` therefore reads the tag, and
        /// which garbage that yields varies by platform -- `c_long` is 4 bytes
        /// here and 8 on Linux, so the two do not even agree on the size of
        /// the enum's variants.
        fn to_bytes(self) -> [u8; 8] {
            let mut bytes = [0u8; 8];
            {
                let mut put = |src: &[u8]| bytes[..src.len()].copy_from_slice(src);
                match self {
                    ColNullValue::UCharNull(v) => put(&v.to_ne_bytes()),
                    ColNullValue::SCharNull(v) => put(&v.to_ne_bytes()),
                    ColNullValue::IntNull(v) => put(&v.to_ne_bytes()),
                    ColNullValue::ShortNull(v) => put(&v.to_ne_bytes()),
                    ColNullValue::LongNull(v) => put(&v.to_ne_bytes()),
                    ColNullValue::UIntNull(v) => put(&v.to_ne_bytes()),
                    ColNullValue::UShortNull(v) => put(&v.to_ne_bytes()),
                    ColNullValue::ULongNull(v) => put(&v.to_ne_bytes()),
                    ColNullValue::FloatNull(v) => put(&v.to_ne_bytes()),
                    ColNullValue::DoubleNull(v) => put(&v.to_ne_bytes()),
                    ColNullValue::LONGLONGNull(v) => put(&v.to_ne_bytes()),
                    /* a string null is a `char*` the caller already holds, and
                    the None case has no value to hand out */
                    ColNullValue::StringNull(_) | ColNullValue::None => {}
                }
            }
            bytes
        }
    }

    let mut dataptr: *mut core::ffi::c_void = core::ptr::null_mut();
    let mut defaultnull: *mut core::ffi::c_void = core::ptr::null_mut();
    let mut col: Vec<ColNulls> = Vec::new();

    let mut tstatus: c_int = 0;
    let mut naxis: c_int = 0;
    let mut bitpix: c_int = 0;
    let mut typecode: c_int = 0;
    let mut hdutype: c_int = 0;
    let mut jtype: c_int = 0;
    let mut dtype: c_int = 0;
    let mut anynul: c_int = 0;
    let mut nfiles: c_int = 0;
    let mut nbytes: c_int = 0;
    let mut totaln: c_long = 0;
    let mut nleft: c_long = 0;
    let mut frow: c_long = 0;
    let mut felement: c_long = 0;
    let mut n_optimum: c_long = 0;
    let mut i_optimum: c_long = 0;
    let mut ntodo: c_long = 0;
    let mut rept: c_long = 0;
    let mut rowrept: c_long = 0;
    let mut width: c_long = 0;
    let mut tnull: c_long = 0;
    let mut naxes: [c_long; 9] = [1, 1, 1, 1, 1, 1, 1, 1, 1];
    let mut groups: c_long = 0;
    let zeros: f64 = 0.0;
    let mut message: [c_char; FLEN_ERRMSG] = [0; FLEN_ERRMSG];
    let mut keyname: [c_char; FLEN_KEYWORD] = [0; FLEN_KEYWORD];
    let mut nullstr: [c_char; FLEN_VALUE] = [0; FLEN_VALUE];
    let mut stringptr: *mut *mut std::os::raw::c_char = core::ptr::null_mut();
    let mut nullpointer: *mut std::os::raw::c_char = core::ptr::null_mut();

    if *status > 0 {
        return *status;
    }

    if n_cols < 0 || n_cols > 999 {
        ffpmsg_str("Illegal number of columms (ffiter)");
        *status = BAD_COL_NUM; /* negative number of columns */
        return *status;
    }

    /*------------------------------------------------------------*/
    /* Make sure column numbers and datatypes are in legal range  */
    /* and column numbers and datatypes are legal.                */
    /* Also fill in other parameters in the column structure.     */
    /*------------------------------------------------------------*/

    /* Assumptions throughout regarding the existence of cols[jj].fptr:
     1) cols[0].fptr must not be NULL
     2) Only TemporaryCols can have NULL fptr.
     3) IMAGE_HDU must not have a TemporaryCol.
    Check the first 2 here. */
    for jj in 0..n_cols as usize {
        if (jj == 0 || cols[jj].iotype != TEMPORARY_COL) && cols[jj].fptr.is_null() {
            ffpmsg_str("Iterator column is missing FITS file pointer (ffiter)");
            *status = NULL_INPUT_PTR;
            return *status;
        }
    }

    unsafe {
        ffghdt_safe(cols[0].fptr.as_mut().unwrap(), &mut hdutype, status); /* type of first HDU */

        let mut jj: usize = 0;
        while jj < n_cols as usize {
            /* check that output datatype code value is legal */
            dtype = cols[jj].datatype;

            /* Allow variable length arrays for InputCol and
            InputOutputCol columns, but not for OutputCol/TemporaryCol
            columns.  Variable length arrays have a negative type code
            value. */

            if !((cols[jj].iotype == OUTPUT_COL) || (cols[jj].iotype == TEMPORARY_COL))
                && (dtype < 0)
            {
                dtype *= -1;
            }

            /* TemporaryCol must have defined datatype and repeat */
            if cols[jj].iotype == TEMPORARY_COL && (dtype <= 0 || cols[jj].repeat <= 0) {
                int_snprintf!(
                    message,
                    FLEN_ERRMSG,
                    "TemporaryCol column must have defined datatype and repeat for column {} (ffiter)",
                    jj + 1
                );
                ffpmsg_slice(&message);
                *status = BAD_DATATYPE;
                return *status;
            }

            /* Check for variable length or illegal data types */
            if dtype != 0
                && dtype != TBYTE
                && dtype != TSBYTE
                && dtype != TLOGICAL
                && dtype != TSTRING
                && dtype != TSHORT
                && dtype != TINT
                && dtype != TLONG
                && dtype != TFLOAT
                && dtype != TDOUBLE
                && dtype != TCOMPLEX
                && dtype != TULONG
                && dtype != TUSHORT
                && dtype != TDBLCOMPLEX
                && dtype != TLONGLONG
            {
                if dtype < 0 {
                    int_snprintf!(
                        message,
                        FLEN_ERRMSG,
                        "Variable length array not allowed for output column number {} (ffiter)",
                        jj + 1
                    );
                } else {
                    int_snprintf!(
                        message,
                        FLEN_ERRMSG,
                        "Illegal datatype for column number {}: {}  (ffiter)",
                        jj + 1,
                        cols[jj].datatype
                    );
                }

                ffpmsg_slice(&message);
                *status = BAD_DATATYPE;
                return *status;
            }

            /* initialize TLMINn, TLMAXn, column name, and display format */
            cols[jj].tlmin = 0;
            cols[jj].tlmax = 0;
            cols[jj].tunit[0] = 0;
            cols[jj].tdisp[0] = 0;

            /* Determine HDU type of this table (or BINARY_TBL for TemporaryCol) */
            if cols[jj].iotype != TEMPORARY_COL {
                ffghdt_safe(cols[jj].fptr.as_mut().unwrap(), &mut jtype, status); /* get HDU type */
            } else {
                hdutype = BINARY_TBL;
            }

            if hdutype == IMAGE_HDU
            /* operating on FITS images */
            {
                if jtype != IMAGE_HDU {
                    int_snprintf!(
                        message,
                        FLEN_ERRMSG,
                        "File {} not positioned to an image extension (ffiter)",
                        jj + 1
                    );
                    *status = NOT_IMAGE;
                    return *status;
                }

                /* since this is an image, set a dummy column number = 0 */
                cols[jj].colnum = 0;
                strcpy_safe(&mut cols[jj].colname, cs!(c"IMAGE")); /* dummy name for images */

                if cols[jj].iotype == TEMPORARY_COL {
                    int_snprintf!(
                        message,
                        FLEN_ERRMSG,
                        "Column type TemporaryCol not permitted for IMAGE HDUs (ffiter)"
                    );
                    *status = BAD_DATATYPE;
                    return *status;
                }

                tstatus = 0;
                ffgkys_safe(
                    cols[jj].fptr.as_mut().unwrap(),
                    cs!(c"BUNIT"),
                    &mut cols[jj].tunit[..],
                    None,
                    &mut tstatus,
                );
            } else {
                /* operating on FITS tables */

                if jtype == IMAGE_HDU {
                    int_snprintf!(
                        message,
                        FLEN_ERRMSG,
                        "File {} not positioned to a table extension (ffiter)",
                        jj + 1
                    );
                    *status = NOT_TABLE;
                    return *status;
                }

                if cols[jj].iotype != TEMPORARY_COL {
                    if cols[jj].colnum < 1 {
                        /* find the column number for the named column */
                        if ffgcno_safe(
                            cols[jj].fptr.as_mut().unwrap(),
                            CASEINSEN as c_int,
                            &cols[jj].colname[..],
                            &mut cols[jj].colnum,
                            status,
                        ) != 0
                        {
                            int_snprintf!(
                                message,
                                FLEN_ERRMSG,
                                "Column '{}' not found for column number {}  (ffiter)",
                                CStr::from_ptr(cols[jj].colname.as_ptr()).to_string_lossy(),
                                jj + 1
                            );
                            ffpmsg_slice(&message);
                            return *status;
                        }
                    }

                    /* check that the column number is valid */
                    if cols[jj].colnum < 1
                        || cols[jj].colnum > ((cols[jj].fptr).as_mut().unwrap().Fptr).tfield
                    {
                        int_snprintf!(
                            message,
                            FLEN_ERRMSG,
                            "Column {} has illegal table position number: {}  (ffiter)",
                            jj + 1,
                            cols[jj].colnum
                        );
                        ffpmsg_slice(&message);
                        *status = BAD_COL_NUM;
                        return *status;
                    }

                    /* look for column description keywords and update structure */
                    tstatus = 0;
                    ffkeyn_safe(
                        cs!(c"TLMIN"),
                        cols[jj].colnum,
                        &mut keyname[..],
                        &mut tstatus,
                    );
                    ffgkyj_safe(
                        cols[jj].fptr.as_mut().unwrap(),
                        &keyname[..],
                        &mut cols[jj].tlmin,
                        None,
                        &mut tstatus,
                    );

                    tstatus = 0;
                    ffkeyn_safe(
                        cs!(c"TLMAX"),
                        cols[jj].colnum,
                        &mut keyname[..],
                        &mut tstatus,
                    );
                    ffgkyj_safe(
                        cols[jj].fptr.as_mut().unwrap(),
                        &keyname[..],
                        &mut cols[jj].tlmax,
                        None,
                        &mut tstatus,
                    );

                    tstatus = 0;
                    ffkeyn_safe(
                        cs!(c"TTYPE"),
                        cols[jj].colnum,
                        &mut keyname[..],
                        &mut tstatus,
                    );
                    ffgkys_safe(
                        cols[jj].fptr.as_mut().unwrap(),
                        &keyname[..],
                        &mut cols[jj].colname[..],
                        None,
                        &mut tstatus,
                    );
                    if tstatus != 0 {
                        cols[jj].colname[0] = 0;
                    }

                    tstatus = 0;
                    ffkeyn_safe(
                        cs!(c"TUNIT"),
                        cols[jj].colnum,
                        &mut keyname[..],
                        &mut tstatus,
                    );
                    ffgkys_safe(
                        cols[jj].fptr.as_mut().unwrap(),
                        &keyname[..],
                        &mut cols[jj].tunit[..],
                        None,
                        &mut tstatus,
                    );

                    tstatus = 0;
                    ffkeyn_safe(
                        cs!(c"TDISP"),
                        cols[jj].colnum,
                        &mut keyname[..],
                        &mut tstatus,
                    );
                    ffgkys_safe(
                        cols[jj].fptr.as_mut().unwrap(),
                        &keyname[..],
                        &mut cols[jj].tdisp[..],
                        None,
                        &mut tstatus,
                    );
                }
            }
            jj += 1;
        } /* end of loop over all columns */

        /*-----------------------------------------------------------------*/
        /* use the first file to set the total number of values to process */
        /*-----------------------------------------------------------------*/

        offset = offset.max(0); /* make sure offset is legal */

        felement = 0;

        /* get total number of pixels in the image */
        if hdutype == IMAGE_HDU {
            fits_get_img_dim(cols[0].fptr.as_mut().unwrap(), &mut naxis, status);
            fits_get_img_size(cols[0].fptr.as_mut().unwrap(), 9, &mut naxes[..], status);

            tstatus = 0;
            ffgkyj_safe(
                cols[0].fptr.as_mut().unwrap(),
                cs!(c"GROUPS"),
                &mut groups,
                None,
                &mut tstatus,
            );

            if tstatus == 0 && groups != 0 && (naxis > 1) && (naxes[0] == 0) {
                /* this is a random groups file, with NAXIS1 = 0 */
                /* Use GCOUNT, the number of groups, as the first multiplier  */
                /* to calculate the total number of pixels in all the groups. */
                ffgkyj_safe(
                    cols[0].fptr.as_mut().unwrap(),
                    cs!(c"GCOUNT"),
                    &mut totaln,
                    None,
                    status,
                );
            } else {
                totaln = naxes[0];
            }

            for ii in 1..naxis as usize {
                totaln *= naxes[ii];
            }

            frow = 1;
            felement = 1 + offset;
        } else {
            /* get total number or rows in the table */

            /* Note the maxvalue here is a special case to deal with
            how the calculator treats expressions that have NO
            referenced columns, just constants and other derivable
            values like #ROW.  In that case, the calculator creates
            a cols[0].fptr even though there is no column for it,
            and the iterator is not meant to allocate any space,
            etc for the column.  So the maxvalue() here assures
            that cols[0] is always checked, even if ncols==0, which
            is how the original logic worked.  This is a bit
            dangerous in the sense that, what happens if the user
            passes a non-calculator input to this iterator, and
            has NOT set fptr to a legitimate FITS handle.  Boom! */
            for jj in 0..n_cols.max(1) as usize {
                if cols[jj].iotype != TEMPORARY_COL && !cols[jj].fptr.is_null() {
                    ffgkyj_safe(
                        cols[jj].fptr.as_mut().unwrap(),
                        cs!(c"NAXIS2"),
                        &mut totaln,
                        None,
                        status,
                    );
                    frow = 1 + offset;
                    felement = 1;
                    break;
                }
            }
            if felement != 1 {
                int_snprintf!(
                    message,
                    FLEN_ERRMSG,
                    "There must be at least one input or output column in iterator (ffiter)"
                );
                ffpmsg_slice(&message);
                *status = BAD_COL_NUM;
                return *status;
            }
        }

        /*  adjust total by the input starting offset value */
        totaln -= offset;
        totaln = totaln.max(0); /* don't allow negative number */

        /*------------------------------------------------------------------*/
        /* Determine number of values to pass to work function on each loop */
        /*------------------------------------------------------------------*/

        if n_per_loop == 0 {
            /* Determine optimum number of values for each iteration.    */
            /* Look at all the fitsfile pointers to determine the number */
            /* of unique files.                                          */

            nfiles = 1;
            n_optimum = 0;
            if cols[0].iotype != TEMPORARY_COL {
                ffgrsz_safe(cols[0].fptr.as_mut().unwrap(), &mut n_optimum, status);
            }

            for jj in 1..n_cols as usize {
                if cols[jj].iotype == TEMPORARY_COL {
                    continue;
                }
                let mut ii = 0;
                while ii < jj {
                    if core::ptr::eq(cols[ii].fptr, cols[jj].fptr) {
                        break;
                    }
                    ii += 1;
                }

                /* this is a new file */
                if ii == jj {
                    nfiles += 1;
                    ffgrsz_safe(cols[jj].fptr.as_mut().unwrap(), &mut i_optimum, status);
                    if n_optimum == 0 {
                        /* If first column is TemporaryCol */
                        n_optimum = i_optimum;
                    } else {
                        n_optimum = n_optimum.min(i_optimum);
                    }
                }
            }

            /* divid n_optimum by the number of files that will be processed */
            n_optimum /= c_long::from(nfiles);
            n_optimum = n_optimum.max(1);
        } else if n_per_loop < 0 {
            /* must pass all the values at one time */
            n_optimum = totaln;
        } else {
            /* calling routine specified how many values to pass at a time */
            n_optimum = n_per_loop.min(totaln);
        }

        /*--------------------------------------*/
        /* allocate work arrays for each column */
        /* and determine the null pixel value   */
        /*--------------------------------------*/

        col = vec![ColNulls::default(); n_cols as usize]; /* memory for the null values */

        #[allow(clippy::never_loop)]
        'cleanup: loop {
            for jj in 0..n_cols as usize {
                /* get image or column datatype and vector length */
                if hdutype == IMAGE_HDU
                /* get total number of pixels in the image */
                {
                    fits_get_img_type(cols[jj].fptr.as_mut().unwrap(), &mut bitpix, status);
                    match bitpix {
                        BYTE_IMG => {
                            typecode = TBYTE;
                        }
                        SHORT_IMG => {
                            typecode = TSHORT;
                        }
                        LONG_IMG => {
                            typecode = TLONG;
                        }
                        FLOAT_IMG => {
                            typecode = TFLOAT;
                        }
                        DOUBLE_IMG => {
                            typecode = TDOUBLE;
                        }
                        LONGLONG_IMG => {
                            typecode = TLONGLONG;
                        }
                        _ => {}
                    }
                } else if cols[jj].iotype != TEMPORARY_COL {
                    if ffgtcl_safe(
                        cols[jj].fptr.as_mut().unwrap(),
                        cols[jj].colnum,
                        Some(&mut typecode),
                        Some(&mut rept),
                        Some(&mut width),
                        status,
                    ) > 0
                    {
                        break 'cleanup;
                    }

                    if typecode < 0 {
                        /* if any variable length arrays, then the */
                        n_optimum = 1; /* must process the table 1 row at a time */

                        /* Allow variable length arrays for InputCol and InputOutputCol columns,
                        but not for OutputCol columns.  Variable length arrays have a
                        negative type code value. */

                        if cols[jj].iotype == OUTPUT_COL {
                            int_snprintf!(
                                message,
                                FLEN_ERRMSG,
                                "Variable length array not allowed for output column number {} (ffiter)",
                                jj + 1
                            );
                            ffpmsg_slice(&message);
                            *status = BAD_DATATYPE;
                            return *status;
                        }
                    }
                } else {
                    /* TemporaryCol - datatype etc must be defined */
                    typecode = cols[jj].datatype;
                    if typecode <= 0 || typecode == TBIT || typecode == TSTRING {
                        int_snprintf!(
                            message,
                            FLEN_ERRMSG,
                            "Invalid typecode for temporary output column number {} (ffiter)",
                            jj + 1
                        );
                        ffpmsg_slice(&message);
                        *status = BAD_DATATYPE;
                        return *status;
                    }

                    rept = cols[jj].repeat;
                    if rept <= 0 {
                        int_snprintf!(
                            message,
                            FLEN_ERRMSG,
                            "Invalid repeat ({}) for temporary output column number {} (ffiter)",
                            rept,
                            jj + 1
                        );
                        ffpmsg_slice(&message);
                        *status = BAD_DIMEN;
                        return *status;
                    }
                }

                /* special case where size_of::<c_long>() = 8: use TINT instead of TLONG */
                if typecode.abs() == TLONG && size_of::<c_long>() == 8 && size_of::<c_int>() == 4 {
                    if typecode < 0 {
                        typecode = -TINT;
                    } else {
                        typecode = TINT;
                    }
                }

                /* Special case: interprete 'X' column as 'B' */
                if typecode.abs() == TBIT {
                    typecode = typecode / TBIT * TBYTE;
                    rept = (rept + 7) / 8;
                }

                if cols[jj].datatype == 0
                /* output datatype not specified? */
                {
                    /* special case if size_of::<c_long>() = 8: use TINT instead of TLONG */
                    if typecode.abs() == TLONG
                        && size_of::<c_long>() == 8
                        && size_of::<c_int>() == 4
                    {
                        cols[jj].datatype = TINT;
                    } else {
                        cols[jj].datatype = typecode.abs();
                    }
                }

                /* calc total number of elements to do on each iteration */
                if hdutype == IMAGE_HDU || cols[jj].datatype == TSTRING {
                    ntodo = n_optimum;
                    cols[jj].repeat = 1;
                    /* handle special case of a 0-width string column */
                    if hdutype == BINARY_TBL && rept == 0 {
                        cols[jj].repeat = 0;
                    }

                    /* get the BLANK keyword value, if it exists */
                    if typecode.abs() == TBYTE
                        || typecode.abs() == TSHORT
                        || typecode.abs() == TLONG
                        || typecode.abs() == TINT
                        || typecode.abs() == TLONGLONG
                    {
                        tstatus = 0;
                        ffgkyj_safe(
                            cols[jj].fptr.as_mut().unwrap(),
                            cs!(c"BLANK"),
                            &mut tnull,
                            None,
                            &mut tstatus,
                        );
                        if tstatus != 0 {
                            tnull = 0; /* no null values */
                        }
                    }
                } else {
                    if typecode < 0 {
                        /* get max size of the variable length vector; dont't trust the value
                        given by the TFORM keyword  */
                        rept = 1;
                        for ii in 0..totaln {
                            ffgdes_safe(
                                cols[jj].fptr.as_mut().unwrap(),
                                cols[jj].colnum,
                                (frow + ii) as LONGLONG,
                                Some(&mut rowrept),
                                None,
                                status,
                            );

                            rept = rept.max(rowrept);
                        }
                    }

                    ntodo = n_optimum * rept; /* vector columns */
                    cols[jj].repeat = rept;

                    /* get the TNULL keyword value, if it exists */
                    if typecode.abs() == TBYTE
                        || typecode.abs() == TSHORT
                        || typecode.abs() == TLONG
                        || typecode.abs() == TINT
                        || typecode.abs() == TLONGLONG
                    {
                        if cols[jj].fptr.is_null() {
                            ffpmsg_str(
                                "Iterator column for table is missing FITS file pointer (ffiter)",
                            );
                            *status = NULL_INPUT_PTR;
                            return *status;
                        }

                        tstatus = 0;
                        if hdutype == ASCII_TBL
                        /* TNULLn value is a string */
                        {
                            ffkeyn_safe(
                                cs!(c"TNULL"),
                                cols[jj].colnum,
                                &mut keyname[..],
                                &mut tstatus,
                            );
                            ffgkys_safe(
                                cols[jj].fptr.as_mut().unwrap(),
                                &keyname[..],
                                &mut nullstr[..],
                                None,
                                &mut tstatus,
                            );
                            if tstatus != 0 {
                                tnull = 0; /* keyword doesn't exist; no null values */
                            } else {
                                let mut cptr = 0; //nullstr

                                while nullstr[cptr] == b' ' as c_char {
                                    /* skip over leading blanks */
                                    cptr += 1;
                                }

                                if nullstr[cptr] == 0 {
                                    /* TNULLn is all blanks? */
                                    tnull = LONG_MIN;
                                } else {
                                    /* attempt to read TNULLn string as an integer */
                                    ffc2ii(&nullstr, &mut tnull, &mut tstatus);

                                    if tstatus != 0 {
                                        tnull = LONG_MIN; /* choose smallest value to represent nulls*/
                                    }
                                }
                            }
                        } else {
                            /* Binary table; TNULLn value is an integer */

                            ffkeyn_safe(
                                cs!(c"TNULL"),
                                cols[jj].colnum,
                                &mut keyname[..],
                                &mut tstatus,
                            );
                            ffgkyj_safe(
                                cols[jj].fptr.as_mut().unwrap(),
                                &keyname[..],
                                &mut tnull,
                                None,
                                &mut tstatus,
                            );
                            if tstatus != 0 {
                                tnull = 0; /* keyword doesn't exist; no null values */
                            } else if tnull == 0 {
                                /* worst possible case: a value of 0 is used to   */
                                /* represent nulls in the FITS file.  We have to  */
                                /* use a non-zero null value here (zero is used to */
                                /* mean there are no null values in the array) so we */
                                /* will use the smallest possible integer instead. */

                                tnull = LONG_MIN; /* choose smallest possible value */
                            }
                        }
                    }
                }

                /* Note that the data array starts with 2nd element;  */
                /* 1st element of the array gives the null data value */

                match cols[jj].datatype {
                    TBYTE => {
                        cols[jj].array = calloc((ntodo + 1) as usize, size_of::<c_char>());
                        col[jj].nullsize = size_of::<c_char>(); /* number of bytes per value */

                        if typecode.abs() == TBYTE
                            || typecode.abs() == TSHORT
                            || typecode.abs() == TLONG
                            || typecode.abs() == TINT
                            || typecode.abs() == TLONGLONG
                        {
                            tnull = tnull.min(255);
                            tnull = tnull.max(0);
                            col[jj].null = ColNullValue::UCharNull(tnull as c_uchar);
                        } else {
                            col[jj].null = ColNullValue::UCharNull(255 as c_uchar); /* use 255 as null */
                        }
                    }
                    TSBYTE => {
                        cols[jj].array = calloc((ntodo + 1) as usize, size_of::<c_char>());
                        col[jj].nullsize = size_of::<c_char>(); /* number of bytes per value */

                        if typecode.abs() == TBYTE
                            || typecode.abs() == TSHORT
                            || typecode.abs() == TLONG
                            || typecode.abs() == TINT
                            || typecode.abs() == TLONGLONG
                        {
                            tnull = tnull.min(127);
                            tnull = tnull.max(-128);
                            col[jj].null = ColNullValue::SCharNull(tnull as c_schar);
                        } else {
                            col[jj].null = ColNullValue::SCharNull(-128 as c_schar); /* use -128  null */
                        }
                    }
                    TSHORT => {
                        cols[jj].array = calloc((ntodo + 1) as usize, size_of::<c_short>());
                        col[jj].nullsize = size_of::<c_short>(); /* number of bytes per value */

                        if typecode.abs() == TBYTE
                            || typecode.abs() == TSHORT
                            || typecode.abs() == TLONG
                            || typecode.abs() == TINT
                            || typecode.abs() == TLONGLONG
                        {
                            tnull = tnull.min(c_long::from(c_short::MAX));
                            tnull = tnull.max(c_long::from(c_short::MIN));
                            col[jj].null = ColNullValue::ShortNull(tnull as c_short);
                        } else {
                            col[jj].null = ColNullValue::ShortNull(c_short::MIN); /* use minimum as null */
                        }
                    }
                    TUSHORT => {
                        cols[jj].array = calloc((ntodo + 1) as usize, size_of::<c_ushort>());
                        col[jj].nullsize = size_of::<c_ushort>(); /* bytes per value */

                        if typecode.abs() == TBYTE
                            || typecode.abs() == TSHORT
                            || typecode.abs() == TLONG
                            || typecode.abs() == TINT
                            || typecode.abs() == TLONGLONG
                        {
                            tnull = tnull.min(c_long::from(c_ushort::MAX));
                            tnull = tnull.max(0); /* don't allow negative value */
                            col[jj].null = ColNullValue::UShortNull(tnull as c_ushort);
                        } else {
                            col[jj].null = ColNullValue::UShortNull(c_ushort::MAX); /* use maximum null */
                        }
                    }
                    TINT => {
                        cols[jj].array = calloc((ntodo + 1) as usize, size_of::<c_int>());
                        col[jj].nullsize = size_of::<c_int>(); /* number of bytes per value */

                        if typecode.abs() == TBYTE
                            || typecode.abs() == TSHORT
                            || typecode.abs() == TLONG
                            || typecode.abs() == TINT
                            || typecode.abs() == TLONGLONG
                        {
                            tnull = tnull.min(c_long::from(c_int::MAX));
                            tnull = tnull.max(c_long::from(c_int::MIN));
                            col[jj].null = ColNullValue::IntNull(tnull as c_int);
                        } else {
                            col[jj].null = ColNullValue::IntNull(c_int::MIN); /* use minimum as null */
                        }
                    }
                    TUINT => {
                        cols[jj].array = calloc((ntodo + 1) as usize, size_of::<c_uint>());
                        col[jj].nullsize = size_of::<c_uint>(); /* bytes per value */

                        if typecode.abs() == TBYTE
                            || typecode.abs() == TSHORT
                            || typecode.abs() == TLONG
                            || typecode.abs() == TINT
                            || typecode.abs() == TLONGLONG
                        {
                            tnull = tnull.min(c_long::from(INT32_MAX));
                            tnull = tnull.max(0);
                            col[jj].null = ColNullValue::UIntNull(tnull as c_uint);
                        } else {
                            col[jj].null = ColNullValue::UIntNull(c_uint::MAX); /* use maximum as null */
                        }
                    }
                    TLONG => {
                        cols[jj].array = calloc((ntodo + 1) as usize, size_of::<c_long>());
                        col[jj].nullsize = size_of::<c_long>(); /* number of bytes per value */

                        if typecode.abs() == TBYTE
                            || typecode.abs() == TSHORT
                            || typecode.abs() == TLONG
                            || typecode.abs() == TINT
                            || typecode.abs() == TLONGLONG
                        {
                            col[jj].null = ColNullValue::LongNull(tnull);
                        } else {
                            col[jj].null = ColNullValue::LongNull(LONG_MIN); /* use minimum as null */
                        }
                    }
                    TULONG => {
                        cols[jj].array = calloc((ntodo + 1) as usize, size_of::<c_ulong>());
                        col[jj].nullsize = size_of::<c_ulong>(); /* bytes per value */

                        if typecode.abs() == TBYTE
                            || typecode.abs() == TSHORT
                            || typecode.abs() == TLONG
                            || typecode.abs() == TINT
                            || typecode.abs() == TLONGLONG
                        {
                            if tnull < 0 {
                                /* can't use a negative null value */
                                col[jj].null = ColNullValue::ULongNull(LONG_MAX as c_ulong);
                            } else {
                                col[jj].null = ColNullValue::ULongNull(tnull as c_ulong);
                            }
                        } else {
                            col[jj].null = ColNullValue::ULongNull(LONG_MAX as c_ulong); /* use maximum as null */
                        }
                    }
                    TFLOAT => {
                        cols[jj].array = calloc((ntodo + 1) as usize, size_of::<f32>());
                        col[jj].nullsize = size_of::<f32>(); /* number of bytes per value */

                        if typecode.abs() == TBYTE
                            || typecode.abs() == TSHORT
                            || typecode.abs() == TLONG
                            || typecode.abs() == TINT
                            || typecode.abs() == TLONGLONG
                        {
                            col[jj].null = ColNullValue::FloatNull(tnull as f32);
                        } else {
                            col[jj].null = ColNullValue::FloatNull(FLOATNULLVALUE); /* special value */
                        }
                    }
                    TCOMPLEX => {
                        cols[jj].array = calloc(((ntodo * 2) + 1) as usize, size_of::<f32>());
                        col[jj].nullsize = size_of::<f32>(); /* number of bytes per value */
                        col[jj].null = ColNullValue::FloatNull(FLOATNULLVALUE); /* special value */
                    }
                    TDOUBLE => {
                        cols[jj].array = calloc((ntodo + 1) as usize, size_of::<f64>());
                        col[jj].nullsize = size_of::<f64>(); /* number of bytes per value */

                        if typecode.abs() == TBYTE
                            || typecode.abs() == TSHORT
                            || typecode.abs() == TLONG
                            || typecode.abs() == TINT
                            || typecode.abs() == TLONGLONG
                        {
                            col[jj].null = ColNullValue::DoubleNull(tnull as f64);
                        } else {
                            col[jj].null = ColNullValue::DoubleNull(DOUBLENULLVALUE); /* special value */
                        }
                    }

                    TDBLCOMPLEX => {
                        cols[jj].array = calloc(((ntodo * 2) + 1) as usize, size_of::<f64>());
                        col[jj].nullsize = size_of::<f64>(); /* number of bytes per value */
                        col[jj].null = ColNullValue::DoubleNull(DOUBLENULLVALUE); /* special value */
                    }
                    TSTRING => {
                        /* allocate array of pointers to all the strings  */
                        if hdutype == ASCII_TBL {
                            rept = width;
                        }
                        stringptr = calloc((ntodo + 1) as usize, size_of::<*mut c_char>())
                            .cast::<*mut c_char>();
                        cols[jj].array = stringptr.cast::<c_void>();
                        col[jj].nullsize = (rept + 1) as usize; /* number of bytes per value */

                        if !stringptr.is_null() {
                            /* allocate string to store the null string value */
                            let stringnull =
                                calloc((rept + 1) as usize, size_of::<c_char>()).cast::<c_char>();
                            col[jj].null = ColNullValue::StringNull(stringnull);
                            if rept > 0
                                && let ColNullValue::StringNull(ptr) = col[jj].null
                            {
                                *ptr.add(1) = 1;
                                /* to make sure string != 0 */
                            }

                            /* allocate big block for the array of table column strings */

                            (*stringptr) =
                                calloc(((ntodo + 1) * (rept + 1)) as usize, size_of::<c_char>())
                                    .cast::<c_char>();

                            if !(*stringptr).is_null() {
                                for ii in 1..=ntodo as usize {
                                    /* pointer to each string */

                                    let prev = *stringptr.add(ii - 1);
                                    *stringptr.add(ii) = prev.add((rept + 1) as usize);
                                }

                                /* get the TNULL keyword value, if it exists */
                                tstatus = 0;
                                ffkeyn_safe(
                                    cs!(c"TNULL"),
                                    cols[jj].colnum,
                                    &mut keyname[..],
                                    &mut tstatus,
                                );
                                ffgkys_safe(
                                    cols[jj].fptr.as_mut().unwrap(),
                                    &keyname[..],
                                    &mut nullstr[..],
                                    None,
                                    &mut tstatus,
                                );
                                if tstatus == 0
                                    && let ColNullValue::StringNull(ptr) = col[jj].null
                                {
                                    // SAFETY: `ptr` is the calloc'd rept+1 byte
                                    // null-value buffer allocated just above.
                                    let nullval =
                                        slice::from_raw_parts_mut(ptr, (rept + 1) as usize);
                                    strncat_safe(nullval, &nullstr, rept as usize);
                                }
                            } else {
                                ffpmsg_str("ffiter failed to allocate memory arrays");
                                *status = MEMORY_ALLOCATION; /* memory allocation failed */
                                break 'cleanup;
                            }
                        }
                    }

                    TLOGICAL => {
                        cols[jj].array = calloc((ntodo + 1) as usize, size_of::<c_char>());
                        col[jj].nullsize = size_of::<c_char>(); /* number of bytes per value */

                        /* use value = 2 to flag null values in logical columns */
                        col[jj].null = ColNullValue::UCharNull(2);
                    }
                    TLONGLONG => {
                        cols[jj].array = calloc((ntodo + 1) as usize, size_of::<LONGLONG>());
                        col[jj].nullsize = size_of::<LONGLONG>(); /* number of bytes per value */

                        if typecode.abs() == TBYTE
                            || typecode.abs() == TSHORT
                            || typecode.abs() == TLONG
                            || typecode.abs() == TLONGLONG
                            || typecode.abs() == TINT
                        {
                            col[jj].null = ColNullValue::LONGLONGNull(tnull as LONGLONG);
                        } else {
                            col[jj].null = ColNullValue::LONGLONGNull(LONGLONG_MIN); /* use minimum as null */
                        }
                    }
                    _ => {
                        int_snprintf!(
                            message,
                            FLEN_ERRMSG,
                            "Column {} datatype currently not supported: {}:  (ffiter)",
                            jj + 1,
                            cols[jj].datatype
                        );
                        ffpmsg_slice(&message);
                        *status = BAD_DATATYPE;
                        break 'cleanup;
                    }
                } /* end of switch block */

                /* check that all the arrays were allocated successfully */
                if cols[jj].array.is_null() {
                    ffpmsg_str("ffiter failed to allocate memory arrays");
                    *status = MEMORY_ALLOCATION; /* memory allocation failed */
                    break 'cleanup;
                }
            }

            /*--------------------------------------------------*/
            /* main loop while there are values left to process */
            /*--------------------------------------------------*/

            nleft = totaln;

            while nleft > 0 {
                ntodo = nleft.min(n_optimum); /* no. of values for this loop */

                /*  read input columns from FITS file(s)  */
                for jj in 0..n_cols as usize {
                    if cols[jj].iotype != OUTPUT_COL && cols[jj].iotype != TEMPORARY_COL {
                        /* the null value's own bytes, standing in for the C
                        union that `defaultnull` points at */
                        let mut nullbytes = col[jj].null.to_bytes();
                        if cols[jj].datatype == TSTRING {
                            stringptr = cols[jj].array.cast::<*mut c_char>();
                            dataptr = stringptr.add(1).cast::<c_void>();
                            defaultnull = match col[jj].null {
                                ColNullValue::StringNull(ptr) => ptr.cast::<c_void>(),
                                _ => core::ptr::null_mut(),
                            }; /* ptr to the null value */
                        } else {
                            dataptr = cols[jj].array.add(col[jj].nullsize);
                            defaultnull = nullbytes.as_mut_ptr().cast::<c_void>(); /* ptr to the null value */
                        }

                        if hdutype == IMAGE_HDU {
                            if ffgpv_safe(
                                cols[jj].fptr.as_mut().unwrap(),
                                cols[jj].datatype,
                                felement as LONGLONG,
                                (cols[jj].repeat * ntodo) as LONGLONG,
                                NullValue::from_raw_ptr(
                                    cols[jj].datatype,
                                    defaultnull.cast_const(),
                                ),
                                slice::from_raw_parts_mut(
                                    dataptr.cast::<u8>(),
                                    bytes_per_datatype(cols[jj].datatype).unwrap()
                                        * (cols[jj].repeat * ntodo) as usize,
                                ),
                                Some(&mut anynul),
                                status,
                            ) > 0
                            {
                                break;
                            }
                        } else {
                            if ffgtcl_safe(
                                cols[jj].fptr.as_mut().unwrap(),
                                cols[jj].colnum,
                                Some(&mut typecode),
                                Some(&mut rept),
                                Some(&mut width),
                                status,
                            ) > 0
                            {
                                break 'cleanup;
                            }

                            if typecode < 0 {
                                /* get size of the variable length vector */
                                ffgdes_safe(
                                    cols[jj].fptr.as_mut().unwrap(),
                                    cols[jj].colnum,
                                    frow as LONGLONG,
                                    Some(&mut cols[jj].repeat),
                                    None,
                                    status,
                                );
                            }

                            if ffgcv_safe(
                                cols[jj].fptr.as_mut().unwrap(),
                                cols[jj].datatype,
                                cols[jj].colnum,
                                frow as LONGLONG,
                                felement as LONGLONG,
                                (cols[jj].repeat * ntodo) as LONGLONG,
                                NullValue::from_raw_ptr(
                                    cols[jj].datatype,
                                    defaultnull.cast_const(),
                                ),
                                slice::from_raw_parts_mut(
                                    dataptr.cast::<u8>(),
                                    bytes_per_datatype(cols[jj].datatype).unwrap()
                                        * (cols[jj].repeat * ntodo) as usize,
                                ),
                                Some(&mut anynul),
                                status,
                            ) > 0
                            {
                                break;
                            }
                        }

                        /* copy the appropriate null value into first array element */

                        /* are there any nulls in the data? */
                        if anynul != 0 {
                            if cols[jj].datatype == TSTRING {
                                /* cols[jj].array is the `char **' block; the null
                                value goes into the string it points at, not over
                                the pointer array itself. */
                                stringptr = cols[jj].array.cast::<*mut c_char>();
                                if let ColNullValue::StringNull(ptr) = col[jj].null {
                                    memcpy(
                                        (*stringptr).cast::<c_void>(),
                                        ptr as *const c_void,
                                        col[jj].nullsize,
                                    );
                                }
                            } else {
                                memcpy(cols[jj].array, defaultnull, col[jj].nullsize);
                            }
                        } else {
                            /* no null values so copy zero into first element */
                            if cols[jj].datatype == TSTRING {
                                stringptr = cols[jj].array.cast::<*mut c_char>();

                                memset((*stringptr).cast::<c_void>(), 0, col[jj].nullsize);
                            } else {
                                memset(cols[jj].array, 0, col[jj].nullsize);
                            }
                        }
                    }
                }

                if *status > 0 {
                    break; /* looks like an error occurred; quit immediately */
                }

                /* call work function */

                if hdutype == IMAGE_HDU {
                    *status =
                        workfn(totaln, offset, felement, ntodo, n_cols, cols_ptr, userPointer);
                } else {
                    *status = workfn(totaln, offset, frow, ntodo, n_cols, cols_ptr, userPointer);
                }

                /* The work function may have reached the column array by its
                own route, so the previous re-borrow is stale: take a fresh
                one from `cols_ptr` before touching the columns again. */
                cols = slice::from_raw_parts_mut(cols_ptr, cols_len);

                if *status > 0 || *status < -1 {
                    break; /* looks like an error occurred; quit immediately */
                }

                /*  write output columns  before quiting if status = -1 */
                tstatus = 0;
                for jj in 0..n_cols as usize {
                    if cols[jj].iotype != INPUT_COL && cols[jj].iotype != TEMPORARY_COL {
                        if cols[jj].datatype == TSTRING {
                            stringptr = cols[jj].array.cast::<*mut c_char>();
                            dataptr = stringptr.add(1).cast::<c_void>();
                            nullpointer = *stringptr;
                            nbytes = 2;
                        } else {
                            dataptr = cols[jj].array.add(col[jj].nullsize);
                            nullpointer = cols[jj].array.cast::<c_char>();
                            nbytes = col[jj].nullsize as c_int;
                        }

                        if memcmp(
                            nullpointer as *const c_void,
                            core::ptr::from_ref::<f64>(&zeros).cast::<c_void>(),
                            nbytes as usize,
                        ) != 0
                        {
                            /* null value flag not zero; must check for and write nulls */
                            if hdutype == IMAGE_HDU {
                                if ffppn_safe(
                                    cols[jj].fptr.as_mut().unwrap(),
                                    cols[jj].datatype,
                                    felement as LONGLONG,
                                    (cols[jj].repeat * ntodo) as LONGLONG,
                                    slice::from_raw_parts(
                                        dataptr as *const u8,
                                        bytes_per_datatype(cols[jj].datatype).unwrap()
                                            * (cols[jj].repeat * ntodo) as usize,
                                    ),
                                    NullValue::from_raw_ptr(
                                        cols[jj].datatype,
                                        nullpointer as *const c_void,
                                    ),
                                    &mut tstatus,
                                ) > 0
                                {
                                    break;
                                }
                            } else {
                                if ffgtcl_safe(
                                    cols[jj].fptr.as_mut().unwrap(),
                                    cols[jj].colnum,
                                    Some(&mut typecode),
                                    Some(&mut rept),
                                    Some(&mut width),
                                    status,
                                ) > 0
                                {
                                    break 'cleanup;
                                }

                                if typecode < 0
                                /* variable length array colum */
                                {
                                    ffgdes_safe(
                                        cols[jj].fptr.as_mut().unwrap(),
                                        cols[jj].colnum,
                                        frow as LONGLONG,
                                        Some(&mut cols[jj].repeat),
                                        None,
                                        status,
                                    );
                                }

                                if ffpcn_safe(
                                    cols[jj].fptr.as_mut().unwrap(),
                                    cols[jj].datatype,
                                    cols[jj].colnum,
                                    frow as LONGLONG,
                                    felement as LONGLONG,
                                    (cols[jj].repeat * ntodo) as LONGLONG,
                                    slice::from_raw_parts(
                                        dataptr as *const u8,
                                        bytes_per_datatype(cols[jj].datatype).unwrap()
                                            * (cols[jj].repeat * ntodo) as usize,
                                    ),
                                    NullValue::from_raw_ptr(
                                        cols[jj].datatype,
                                        nullpointer as *const c_void,
                                    ),
                                    &mut tstatus,
                                ) > 0
                                {
                                    break;
                                }
                            }
                        } else {
                            /* no null values; just write the array */
                            if hdutype == IMAGE_HDU {
                                if ffppr_safe(
                                    cols[jj].fptr.as_mut().unwrap(),
                                    cols[jj].datatype,
                                    felement as LONGLONG,
                                    (cols[jj].repeat * ntodo) as LONGLONG,
                                    slice::from_raw_parts(
                                        dataptr as *const u8,
                                        bytes_per_datatype(cols[jj].datatype).unwrap()
                                            * (cols[jj].repeat * ntodo) as usize,
                                    ),
                                    &mut tstatus,
                                ) > 0
                                {
                                    break;
                                }
                            } else {
                                if ffgtcl_safe(
                                    cols[jj].fptr.as_mut().unwrap(),
                                    cols[jj].colnum,
                                    Some(&mut typecode),
                                    Some(&mut rept),
                                    Some(&mut width),
                                    status,
                                ) > 0
                                {
                                    break 'cleanup;
                                }

                                if typecode < 0
                                /* variable length array column */
                                {
                                    ffgdes_safe(
                                        cols[jj].fptr.as_mut().unwrap(),
                                        cols[jj].colnum,
                                        frow as LONGLONG,
                                        Some(&mut cols[jj].repeat),
                                        None,
                                        status,
                                    );
                                }

                                if ffpcl_safe(
                                    cols[jj].fptr.as_mut().unwrap(),
                                    cols[jj].datatype,
                                    cols[jj].colnum,
                                    frow as LONGLONG,
                                    felement as LONGLONG,
                                    (cols[jj].repeat * ntodo) as LONGLONG,
                                    slice::from_raw_parts_mut(
                                        dataptr.cast::<u8>(),
                                        bytes_per_datatype(cols[jj].datatype).unwrap()
                                            * (cols[jj].repeat * ntodo) as usize,
                                    ),
                                    &mut tstatus,
                                ) > 0
                                {
                                    break;
                                }
                            }
                        }
                    }
                }

                if *status == 0 {
                    *status = tstatus; /* propagate any error status from the writes */
                }

                if *status != 0 {
                    break; /* exit on any error */
                }

                nleft -= ntodo;

                if hdutype == IMAGE_HDU {
                    felement += ntodo;
                } else {
                    frow += ntodo;
                }
            }
            break;
        }

        // cleanup:

        /*----------------------------------*/
        /* free work arrays for the columns */
        /*----------------------------------*/

        for jj in 0..n_cols as usize {
            /* The C guards these on `if (cols[jj].array)`, i.e. non-null.
            Testing is_null() instead meant the frees only ran for columns
            that had nothing to free, so every allocated work array leaked. */
            if cols[jj].datatype == TSTRING && !cols[jj].array.is_null() {
                stringptr = cols[jj].array.cast::<*mut c_char>();

                free((*stringptr).cast::<c_void>()); /* free the block of strings */
                if let ColNullValue::StringNull(ptr) = col[jj].null {
                    free(ptr.cast::<c_void>()); /* free the null string */
                }
            }
            if !cols[jj].array.is_null() {
                free(cols[jj].array);
                /* memory for the array of values from the col */
            }
        }
    }

    // col is a Vec, so it will be automatically freed
    *status
}

// Ported from test_putcol.c - generic column/pixel write functions
// (ffppx, ffppxll, ffppxnll, ffppn, ffpss, ffpcn, ffpthp). These take an
// explicit `datatype` argument. Data is written then read back to verify.

#[cfg(test)]
mod tests {
    use super::*;
    use crate::NullValue;

    use crate::fitsio::{
        BINARY_TBL, BYTE_IMG, DOUBLE_IMG, FLOAT_IMG, LONG_IMG, LONGLONG, LONGLONG_IMG, READONLY,
        SBYTE_IMG, SHORT_IMG, TBYTE, TCOMPLEX, TDBLCOMPLEX, TDOUBLE, TFLOAT, TINT, TLOGICAL, TLONG,
        TLONGLONG, TSBYTE, TSHORT, TSTRING, TUINT, TULONG, TULONGLONG, TUSHORT, ULONG_IMG,
        USHORT_IMG, fitsfile,
    };
    use crate::helpers::testhelpers::{to_buf, with_temp_file};
    use alloc::ffi::CString;
    use bytemuck::{cast_slice, cast_slice_mut};
    use libc::{c_char, c_int, c_long};

    /// Make a NUL-terminated `Vec<c_char>` from a `&str`.
    fn cc(s: &str) -> Vec<c_char> {
        let mut v: Vec<c_char> = s.bytes().map(|b| b as c_char).collect();
        v.push(0);
        v
    }

    /// Create a single-column binary table and return the open file.
    fn make_table(
        name: &[c_char],
        ttype: &str,
        tform: &str,
        nrows: LONGLONG,
        status: &mut c_int,
    ) -> Option<Box<fitsfile>> {
        let mut f: Option<Box<fitsfile>> = None;
        fits_create_file(&mut f, name, status);
        fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], status);
        let ttype_v = [Some(cc(ttype))];
        let ttype_ref: Vec<Option<&[c_char]>> = ttype_v.iter().map(|o| o.as_deref()).collect();
        let tform_v = [cc(tform)];
        let tform_ref: Vec<&[c_char]> = tform_v.iter().map(|v| v.as_slice()).collect();
        fits_create_tbl(
            f.as_deref_mut().unwrap(),
            BINARY_TBL,
            nrows,
            1,
            &ttype_ref,
            &tform_ref,
            None,
            None,
            status,
        );
        f
    }

    // ----------------------------------------------------------------------
    // ffppx tests (long *firstpix)
    // ----------------------------------------------------------------------

    #[test]
    fn test_ffppx_sbyte() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 1] = [5];
            let data: [i8; 5] = [-100, -50, 0, 50, 100];
            let firstpix: [c_long; 1] = [1];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            // SBYTE stored in short
            fits_write_imghdr(f.as_deref_mut().unwrap(), SHORT_IMG, 1, &naxes, &mut status);
            fits_write_pix(
                f.as_deref_mut().unwrap(),
                TSBYTE,
                &firstpix,
                5,
                cast_slice(&data),
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);
            assert_eq!(status, 0, "setup");

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut result = [0i8; 5];
            let mut anynull = -1;
            fits_read_img(
                f.as_deref_mut().unwrap(),
                TSBYTE,
                1,
                5,
                None,
                cast_slice_mut(&mut result),
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result, data);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_ffppx_ushort() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 1] = [5];
            let data: [u16; 5] = [0, 1000, 30000, 50000, 65535];
            let firstpix: [c_long; 1] = [1];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(
                f.as_deref_mut().unwrap(),
                USHORT_IMG,
                1,
                &naxes,
                &mut status,
            );
            fits_write_pix(
                f.as_deref_mut().unwrap(),
                TUSHORT,
                &firstpix,
                5,
                cast_slice(&data),
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut result = [0u16; 5];
            let mut anynull = -1;
            fits_read_img(
                f.as_deref_mut().unwrap(),
                TUSHORT,
                1,
                5,
                None,
                cast_slice_mut(&mut result),
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result, data);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_ffppx_uint() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 1] = [3];
            let data: [u32; 3] = [0, 1000000, 4294967295];
            let firstpix: [c_long; 1] = [1];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), ULONG_IMG, 1, &naxes, &mut status);
            fits_write_pix(
                f.as_deref_mut().unwrap(),
                TUINT,
                &firstpix,
                3,
                cast_slice(&data),
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut result = [0u32; 3];
            let mut anynull = -1;
            fits_read_img(
                f.as_deref_mut().unwrap(),
                TUINT,
                1,
                3,
                None,
                cast_slice_mut(&mut result),
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result, data);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_ffppx_ulong() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 1] = [3];
            let data: [libc::c_ulong; 3] = [0, 100000, 4294967295];
            let firstpix: [c_long; 1] = [1];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), ULONG_IMG, 1, &naxes, &mut status);
            fits_write_pix(
                f.as_deref_mut().unwrap(),
                TULONG,
                &firstpix,
                3,
                cast_slice(&data),
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut result = [0 as libc::c_ulong; 3];
            let mut anynull = -1;
            fits_read_img(
                f.as_deref_mut().unwrap(),
                TULONG,
                1,
                3,
                None,
                cast_slice_mut(&mut result),
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result, data);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_ffppx_ulonglong() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 1] = [3];
            // Use values that fit in LONGLONG_IMG range (signed 64-bit).
            let data: [u64; 3] = [0, 1000000000, 9000000000000000000];
            let firstpix: [c_long; 1] = [1];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(
                f.as_deref_mut().unwrap(),
                LONGLONG_IMG,
                1,
                &naxes,
                &mut status,
            );
            fits_write_pix(
                f.as_deref_mut().unwrap(),
                TULONGLONG,
                &firstpix,
                3,
                cast_slice(&data),
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut result = [0u64; 3];
            let mut anynull = -1;
            fits_read_img(
                f.as_deref_mut().unwrap(),
                TULONGLONG,
                1,
                3,
                None,
                cast_slice_mut(&mut result),
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result, data);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_ffppx_longlong() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 1] = [3];
            let data: [LONGLONG; 3] = [-9000000000, 0, 9000000000];
            let firstpix: [c_long; 1] = [1];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(
                f.as_deref_mut().unwrap(),
                LONGLONG_IMG,
                1,
                &naxes,
                &mut status,
            );
            fits_write_pix(
                f.as_deref_mut().unwrap(),
                TLONGLONG,
                &firstpix,
                3,
                cast_slice(&data),
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut result = [0 as LONGLONG; 3];
            let mut anynull = -1;
            fits_read_img(
                f.as_deref_mut().unwrap(),
                TLONGLONG,
                1,
                3,
                None,
                cast_slice_mut(&mut result),
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result, data);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    // ----------------------------------------------------------------------
    // ffppxll tests (LONGLONG *firstpix)
    // ----------------------------------------------------------------------

    #[test]
    fn test_ffppxll_byte() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 1] = [5];
            let data: [u8; 5] = [10, 20, 30, 40, 50];
            let firstpix: [LONGLONG; 1] = [1];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 1, &naxes, &mut status);
            fits_write_pixll(
                f.as_deref_mut().unwrap(),
                TBYTE,
                &firstpix,
                5,
                &data,
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut result = [0u8; 5];
            let mut anynull = -1;
            fits_read_img(
                f.as_deref_mut().unwrap(),
                TBYTE,
                1,
                5,
                None,
                &mut result,
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result, data);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_ffppxll_short() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 1] = [5];
            let data: [i16; 5] = [-100, 0, 100, 200, 300];
            let firstpix: [LONGLONG; 1] = [1];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), SHORT_IMG, 1, &naxes, &mut status);
            fits_write_pixll(
                f.as_deref_mut().unwrap(),
                TSHORT,
                &firstpix,
                5,
                cast_slice(&data),
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut result = [0i16; 5];
            let mut anynull = -1;
            fits_read_img(
                f.as_deref_mut().unwrap(),
                TSHORT,
                1,
                5,
                None,
                cast_slice_mut(&mut result),
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result, data);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_ffppxll_int() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 1] = [3];
            let data: [c_int; 3] = [-1000000, 0, 1000000];
            let firstpix: [LONGLONG; 1] = [1];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), LONG_IMG, 1, &naxes, &mut status);
            fits_write_pixll(
                f.as_deref_mut().unwrap(),
                TINT,
                &firstpix,
                3,
                cast_slice(&data),
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut result = [0 as c_int; 3];
            let mut anynull = -1;
            fits_read_img(
                f.as_deref_mut().unwrap(),
                TINT,
                1,
                3,
                None,
                cast_slice_mut(&mut result),
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result, data);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_ffppxll_longlong() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 1] = [3];
            let data: [LONGLONG; 3] = [-9000000000, 0, 9000000000];
            let firstpix: [LONGLONG; 1] = [1];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(
                f.as_deref_mut().unwrap(),
                LONGLONG_IMG,
                1,
                &naxes,
                &mut status,
            );
            fits_write_pixll(
                f.as_deref_mut().unwrap(),
                TLONGLONG,
                &firstpix,
                3,
                cast_slice(&data),
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut result = [0 as LONGLONG; 3];
            let mut anynull = -1;
            fits_read_img(
                f.as_deref_mut().unwrap(),
                TLONGLONG,
                1,
                3,
                None,
                cast_slice_mut(&mut result),
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result, data);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_ffppxll_float() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 1] = [4];
            let data: [f32; 4] = [1.5, 2.5, 3.5, 4.5];
            let firstpix: [LONGLONG; 1] = [1];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), FLOAT_IMG, 1, &naxes, &mut status);
            fits_write_pixll(
                f.as_deref_mut().unwrap(),
                TFLOAT,
                &firstpix,
                4,
                cast_slice(&data),
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut result = [0f32; 4];
            let mut anynull = -1;
            fits_read_img(
                f.as_deref_mut().unwrap(),
                TFLOAT,
                1,
                4,
                None,
                cast_slice_mut(&mut result),
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result, data);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_ffppxll_double() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 1] = [4];
            let data: [f64; 4] = [1.5, 2.5, 3.5, 4.5];
            let firstpix: [LONGLONG; 1] = [1];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(
                f.as_deref_mut().unwrap(),
                DOUBLE_IMG,
                1,
                &naxes,
                &mut status,
            );
            fits_write_pixll(
                f.as_deref_mut().unwrap(),
                TDOUBLE,
                &firstpix,
                4,
                cast_slice(&data),
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut result = [0f64; 4];
            let mut anynull = -1;
            fits_read_img(
                f.as_deref_mut().unwrap(),
                TDOUBLE,
                1,
                4,
                None,
                cast_slice_mut(&mut result),
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result, data);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_ffppxll_ushort() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 1] = [3];
            let data: [u16; 3] = [0, 32768, 65535];
            let firstpix: [LONGLONG; 1] = [1];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(
                f.as_deref_mut().unwrap(),
                USHORT_IMG,
                1,
                &naxes,
                &mut status,
            );
            fits_write_pixll(
                f.as_deref_mut().unwrap(),
                TUSHORT,
                &firstpix,
                3,
                cast_slice(&data),
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut result = [0u16; 3];
            let mut anynull = -1;
            fits_read_img(
                f.as_deref_mut().unwrap(),
                TUSHORT,
                1,
                3,
                None,
                cast_slice_mut(&mut result),
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result, data);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_ffppxll_ulong() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 1] = [3];
            let data: [libc::c_ulong; 3] = [0, 2147483648, 4294967295];
            let firstpix: [LONGLONG; 1] = [1];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), ULONG_IMG, 1, &naxes, &mut status);
            fits_write_pixll(
                f.as_deref_mut().unwrap(),
                TULONG,
                &firstpix,
                3,
                cast_slice(&data),
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut result = [0 as libc::c_ulong; 3];
            let mut anynull = -1;
            fits_read_img(
                f.as_deref_mut().unwrap(),
                TULONG,
                1,
                3,
                None,
                cast_slice_mut(&mut result),
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result, data);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_ffppxll_sbyte() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 1] = [5];
            let data: [i8; 5] = [-128, -64, 0, 64, 127];
            let firstpix: [LONGLONG; 1] = [1];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), SBYTE_IMG, 1, &naxes, &mut status);
            fits_write_pixll(
                f.as_deref_mut().unwrap(),
                TSBYTE,
                &firstpix,
                5,
                cast_slice(&data),
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut result = [0i8; 5];
            let mut anynull = -1;
            fits_read_img(
                f.as_deref_mut().unwrap(),
                TSBYTE,
                1,
                5,
                None,
                cast_slice_mut(&mut result),
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result, data);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_ffppxll_uint() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 1] = [3];
            let data: [u32; 3] = [0, 2000000000, 4000000000];
            let firstpix: [LONGLONG; 1] = [1];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), ULONG_IMG, 1, &naxes, &mut status);
            fits_write_pixll(
                f.as_deref_mut().unwrap(),
                TUINT,
                &firstpix,
                3,
                cast_slice(&data),
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut result = [0u32; 3];
            let mut anynull = -1;
            fits_read_img(
                f.as_deref_mut().unwrap(),
                TUINT,
                1,
                3,
                None,
                cast_slice_mut(&mut result),
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result, data);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_ffppxll_long() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 1] = [3];
            let data: [c_long; 3] = [-2000000000, 0, 2000000000];
            let firstpix: [LONGLONG; 1] = [1];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), LONG_IMG, 1, &naxes, &mut status);
            fits_write_pixll(
                f.as_deref_mut().unwrap(),
                TLONG,
                &firstpix,
                3,
                cast_slice(&data),
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut result = [0 as c_long; 3];
            let mut anynull = -1;
            fits_read_img(
                f.as_deref_mut().unwrap(),
                TLONG,
                1,
                3,
                None,
                cast_slice_mut(&mut result),
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result, data);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    // ----------------------------------------------------------------------
    // ffppn tests (write primary array with null value)
    // ----------------------------------------------------------------------

    #[test]
    fn test_ffppn_short() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 1] = [5];
            let data: [i16; 5] = [10, -999, 20, -999, 30];
            let nulval: i16 = -999;

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), SHORT_IMG, 1, &naxes, &mut status);
            fits_set_imgnull(f.as_deref_mut().unwrap(), nulval as LONGLONG, &mut status);
            fits_write_imgnull(
                f.as_deref_mut().unwrap(),
                TSHORT,
                1,
                5,
                cast_slice(&data),
                Some(NullValue::Short(nulval)),
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut result = [0i16; 5];
            let mut anynull = -1;
            fits_read_img(
                f.as_deref_mut().unwrap(),
                TSHORT,
                1,
                5,
                None,
                cast_slice_mut(&mut result),
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[0], 10);
            assert_eq!(result[2], 20);
            assert_eq!(result[4], 30);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_ffppn_int() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 1] = [5];
            let data: [c_int; 5] = [100, -99999, 200, -99999, 300];
            let nulval: c_int = -99999;

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), LONG_IMG, 1, &naxes, &mut status);
            fits_set_imgnull(f.as_deref_mut().unwrap(), nulval as LONGLONG, &mut status);
            fits_write_imgnull(
                f.as_deref_mut().unwrap(),
                TINT,
                1,
                5,
                cast_slice(&data),
                Some(NullValue::Int(nulval)),
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut result = [0 as c_int; 5];
            let mut anynull = -1;
            fits_read_img(
                f.as_deref_mut().unwrap(),
                TINT,
                1,
                5,
                None,
                cast_slice_mut(&mut result),
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[0], 100);
            assert_eq!(result[2], 200);
            assert_eq!(result[4], 300);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_ffppn_float() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 1] = [5];
            let data: [f32; 5] = [1.5, -999.0, 2.5, -999.0, 3.5];
            let nulval: f32 = -999.0;

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), FLOAT_IMG, 1, &naxes, &mut status);
            fits_write_imgnull(
                f.as_deref_mut().unwrap(),
                TFLOAT,
                1,
                5,
                cast_slice(&data),
                Some(NullValue::Float(nulval)),
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            // Use a different null value to get actual stored values
            let mut result = [0f32; 5];
            let mut anynull = -1;
            fits_read_img(
                f.as_deref_mut().unwrap(),
                TFLOAT,
                1,
                5,
                Some(NullValue::Float(-888.0)),
                cast_slice_mut(&mut result),
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[0], 1.5);
            assert_eq!(result[2], 2.5);
            assert_eq!(result[4], 3.5);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_ffppn_double() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 1] = [5];
            let data: [f64; 5] = [1.5, -999.0, 2.5, -999.0, 3.5];
            let nulval: f64 = -999.0;

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(
                f.as_deref_mut().unwrap(),
                DOUBLE_IMG,
                1,
                &naxes,
                &mut status,
            );
            fits_write_imgnull(
                f.as_deref_mut().unwrap(),
                TDOUBLE,
                1,
                5,
                cast_slice(&data),
                Some(NullValue::Double(nulval)),
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut result = [0f64; 5];
            let mut anynull = -1;
            fits_read_img(
                f.as_deref_mut().unwrap(),
                TDOUBLE,
                1,
                5,
                Some(NullValue::Double(-888.0)),
                cast_slice_mut(&mut result),
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[0], 1.5);
            assert_eq!(result[2], 2.5);
            assert_eq!(result[4], 3.5);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_ffppn_byte() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 1] = [5];
            let data: [u8; 5] = [10, 255, 20, 255, 30];
            let nulval: u8 = 255;

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 1, &naxes, &mut status);
            fits_set_imgnull(f.as_deref_mut().unwrap(), nulval as LONGLONG, &mut status);
            fits_write_imgnull(
                f.as_deref_mut().unwrap(),
                TBYTE,
                1,
                5,
                &data,
                Some(NullValue::UByte(nulval)),
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut result = [0u8; 5];
            let mut anynull = -1;
            fits_read_img(
                f.as_deref_mut().unwrap(),
                TBYTE,
                1,
                5,
                None,
                &mut result,
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[0], 10);
            assert_eq!(result[2], 20);
            assert_eq!(result[4], 30);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_ffppn_sbyte() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 1] = [5];
            let data: [i8; 5] = [10, -128, 20, -128, 30];
            let nulval: i8 = -128;

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), SBYTE_IMG, 1, &naxes, &mut status);
            fits_set_imgnull(f.as_deref_mut().unwrap(), nulval as LONGLONG, &mut status);
            fits_write_imgnull(
                f.as_deref_mut().unwrap(),
                TSBYTE,
                1,
                5,
                cast_slice(&data),
                Some(NullValue::Byte(nulval)),
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut result = [0i8; 5];
            let mut anynull = -1;
            fits_read_img(
                f.as_deref_mut().unwrap(),
                TSBYTE,
                1,
                5,
                None,
                cast_slice_mut(&mut result),
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[0], 10);
            assert_eq!(result[2], 20);
            assert_eq!(result[4], 30);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_ffppn_ushort() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 1] = [5];
            let data: [u16; 5] = [100, 65535, 200, 65535, 300];
            let nulval: u16 = 65535;

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(
                f.as_deref_mut().unwrap(),
                USHORT_IMG,
                1,
                &naxes,
                &mut status,
            );
            fits_set_imgnull(f.as_deref_mut().unwrap(), nulval as LONGLONG, &mut status);
            fits_write_imgnull(
                f.as_deref_mut().unwrap(),
                TUSHORT,
                1,
                5,
                cast_slice(&data),
                Some(NullValue::UShort(nulval)),
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut result = [0u16; 5];
            let mut anynull = -1;
            fits_read_img(
                f.as_deref_mut().unwrap(),
                TUSHORT,
                1,
                5,
                None,
                cast_slice_mut(&mut result),
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[0], 100);
            assert_eq!(result[2], 200);
            assert_eq!(result[4], 300);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_ffppn_uint() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 1] = [5];
            let data: [u32; 5] = [100, 4294967295, 200, 4294967295, 300];
            let nulval: u32 = 4294967295;

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), ULONG_IMG, 1, &naxes, &mut status);
            fits_set_imgnull(f.as_deref_mut().unwrap(), nulval as LONGLONG, &mut status);
            fits_write_imgnull(
                f.as_deref_mut().unwrap(),
                TUINT,
                1,
                5,
                cast_slice(&data),
                Some(NullValue::UInt(nulval)),
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut result = [0u32; 5];
            let mut anynull = -1;
            fits_read_img(
                f.as_deref_mut().unwrap(),
                TUINT,
                1,
                5,
                None,
                cast_slice_mut(&mut result),
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[0], 100);
            assert_eq!(result[2], 200);
            assert_eq!(result[4], 300);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_ffppn_ulong() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 1] = [5];
            let data: [libc::c_ulong; 5] = [100, 4294967295, 200, 4294967295, 300];
            let nulval: libc::c_ulong = 4294967295;

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), ULONG_IMG, 1, &naxes, &mut status);
            fits_set_imgnull(f.as_deref_mut().unwrap(), nulval as LONGLONG, &mut status);
            fits_write_imgnull(
                f.as_deref_mut().unwrap(),
                TULONG,
                1,
                5,
                cast_slice(&data),
                Some(NullValue::ULong(nulval)),
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut result = [0 as libc::c_ulong; 5];
            let mut anynull = -1;
            fits_read_img(
                f.as_deref_mut().unwrap(),
                TULONG,
                1,
                5,
                None,
                cast_slice_mut(&mut result),
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[0], 100);
            assert_eq!(result[2], 200);
            assert_eq!(result[4], 300);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_ffppn_long() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 1] = [5];
            let data: [c_long; 5] = [100, -2147483648, 200, -2147483648, 300];
            let nulval: c_long = -2147483648;

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), LONG_IMG, 1, &naxes, &mut status);
            fits_set_imgnull(f.as_deref_mut().unwrap(), nulval as LONGLONG, &mut status);
            fits_write_imgnull(
                f.as_deref_mut().unwrap(),
                TLONG,
                1,
                5,
                cast_slice(&data),
                Some(NullValue::Long(nulval)),
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut result = [0 as c_long; 5];
            let mut anynull = -1;
            fits_read_img(
                f.as_deref_mut().unwrap(),
                TLONG,
                1,
                5,
                None,
                cast_slice_mut(&mut result),
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[0], 100);
            assert_eq!(result[2], 200);
            assert_eq!(result[4], 300);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_ffppn_longlong() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 1] = [5];
            let data: [LONGLONG; 5] = [100, -9223372036854775807, 200, -9223372036854775807, 300];
            let nulval: LONGLONG = -9223372036854775807;

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(
                f.as_deref_mut().unwrap(),
                LONGLONG_IMG,
                1,
                &naxes,
                &mut status,
            );
            fits_set_imgnull(f.as_deref_mut().unwrap(), nulval, &mut status);
            fits_write_imgnull(
                f.as_deref_mut().unwrap(),
                TLONGLONG,
                1,
                5,
                cast_slice(&data),
                Some(NullValue::LONGLONG(nulval)),
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut result = [0 as LONGLONG; 5];
            let mut anynull = -1;
            fits_read_img(
                f.as_deref_mut().unwrap(),
                TLONGLONG,
                1,
                5,
                None,
                cast_slice_mut(&mut result),
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[0], 100);
            assert_eq!(result[2], 200);
            assert_eq!(result[4], 300);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    // ----------------------------------------------------------------------
    // ffpss tests (write subset)
    // ----------------------------------------------------------------------

    #[test]
    fn test_ffpss_short() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 2] = [4, 4];
            let data: [i16; 4] = [11, 12, 21, 22];
            let blc: [c_long; 2] = [2, 2];
            let trc: [c_long; 2] = [3, 3];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), SHORT_IMG, 2, &naxes, &mut status);
            fits_write_subset(
                f.as_deref_mut().unwrap(),
                TSHORT,
                &blc,
                &trc,
                cast_slice(&data),
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut result = [0i16; 16];
            let mut anynull = -1;
            fits_read_img(
                f.as_deref_mut().unwrap(),
                TSHORT,
                1,
                16,
                None,
                cast_slice_mut(&mut result),
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            // Data at (2,2) should be at index 5 (row 2, col 2 = 1*4 + 1)
            assert_eq!(result[5], 11);
            assert_eq!(result[6], 12);
            assert_eq!(result[9], 21);
            assert_eq!(result[10], 22);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_ffpss_float() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 2] = [5, 5];
            let data: [f32; 9] = [1.1, 1.2, 1.3, 2.1, 2.2, 2.3, 3.1, 3.2, 3.3];
            let blc: [c_long; 2] = [2, 2];
            let trc: [c_long; 2] = [4, 4];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), FLOAT_IMG, 2, &naxes, &mut status);
            fits_write_subset(
                f.as_deref_mut().unwrap(),
                TFLOAT,
                &blc,
                &trc,
                cast_slice(&data),
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut result = [0f32; 25];
            let mut anynull = -1;
            fits_read_img(
                f.as_deref_mut().unwrap(),
                TFLOAT,
                1,
                25,
                None,
                cast_slice_mut(&mut result),
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            // Data at (2,2) should be at index 6 (row 2, col 2 = 1*5 + 1)
            assert_eq!(result[6], 1.1);
            assert_eq!(result[7], 1.2);
            assert_eq!(result[8], 1.3);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    // ----------------------------------------------------------------------
    // ffpcn tests (write column with null)
    // ----------------------------------------------------------------------

    #[test]
    fn test_ffpcn_short() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [i16; 5] = [10, -999, 20, -999, 30];
            let nulval: i16 = -999;

            let mut f = make_table(&name, "VALUE", "1J", 5, &mut status);
            fits_set_btblnull(
                f.as_deref_mut().unwrap(),
                1,
                nulval as LONGLONG,
                &mut status,
            );
            fits_write_colnull(
                f.as_deref_mut().unwrap(),
                TSHORT,
                1,
                1,
                1,
                5,
                cast_slice(&data),
                Some(NullValue::Short(nulval)),
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0i16; 5];
            let mut anynull = -1;
            fits_read_col(
                f.as_deref_mut().unwrap(),
                TSHORT,
                1,
                1,
                1,
                5,
                Some(NullValue::Short(-888)),
                cast_slice_mut(&mut result),
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[0], 10);
            assert_eq!(result[2], 20);
            assert_eq!(result[4], 30);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_ffpcn_int() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [c_int; 5] = [100, -99999, 200, -99999, 300];
            let nulval: c_int = -99999;

            let mut f = make_table(&name, "VALUE", "1J", 5, &mut status);
            // A TNULL value must be defined for an integer column before null
            // values can be written.
            fits_set_btblnull(
                f.as_deref_mut().unwrap(),
                1,
                nulval as LONGLONG,
                &mut status,
            );
            fits_write_colnull(
                f.as_deref_mut().unwrap(),
                TINT,
                1,
                1,
                1,
                5,
                cast_slice(&data),
                Some(NullValue::Int(nulval)),
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0 as c_int; 5];
            let mut anynull = -1;
            fits_read_col(
                f.as_deref_mut().unwrap(),
                TINT,
                1,
                1,
                1,
                5,
                Some(NullValue::Int(-88888)),
                cast_slice_mut(&mut result),
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[0], 100);
            assert_eq!(result[2], 200);
            assert_eq!(result[4], 300);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_ffpcn_float() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [f32; 5] = [1.5, -999.0, 2.5, -999.0, 3.5];
            let nulval: f32 = -999.0;

            let mut f = make_table(&name, "VALUE", "1E", 5, &mut status);
            fits_write_colnull(
                f.as_deref_mut().unwrap(),
                TFLOAT,
                1,
                1,
                1,
                5,
                cast_slice(&data),
                Some(NullValue::Float(nulval)),
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0f32; 5];
            let mut anynull = -1;
            fits_read_col(
                f.as_deref_mut().unwrap(),
                TFLOAT,
                1,
                1,
                1,
                5,
                Some(NullValue::Float(-888.0)),
                cast_slice_mut(&mut result),
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[0], 1.5);
            assert_eq!(result[2], 2.5);
            assert_eq!(result[4], 3.5);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_ffpcn_double() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [f64; 5] = [1.5, -999.0, 2.5, -999.0, 3.5];
            let nulval: f64 = -999.0;

            let mut f = make_table(&name, "VALUE", "1D", 5, &mut status);
            fits_write_colnull(
                f.as_deref_mut().unwrap(),
                TDOUBLE,
                1,
                1,
                1,
                5,
                cast_slice(&data),
                Some(NullValue::Double(nulval)),
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0f64; 5];
            let mut anynull = -1;
            fits_read_col(
                f.as_deref_mut().unwrap(),
                TDOUBLE,
                1,
                1,
                1,
                5,
                Some(NullValue::Double(-888.0)),
                cast_slice_mut(&mut result),
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[0], 1.5);
            assert_eq!(result[2], 2.5);
            assert_eq!(result[4], 3.5);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_ffpcn_byte() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [u8; 5] = [10, 255, 20, 255, 30];
            let nulval: u8 = 255;

            let mut f = make_table(&name, "VALUE", "1B", 5, &mut status);
            fits_set_btblnull(
                f.as_deref_mut().unwrap(),
                1,
                nulval as LONGLONG,
                &mut status,
            );
            fits_write_colnull(
                f.as_deref_mut().unwrap(),
                TBYTE,
                1,
                1,
                1,
                5,
                &data,
                Some(NullValue::UByte(nulval)),
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0u8; 5];
            let mut anynull = -1;
            fits_read_col(
                f.as_deref_mut().unwrap(),
                TBYTE,
                1,
                1,
                1,
                5,
                Some(NullValue::UByte(254)),
                &mut result,
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[0], 10);
            assert_eq!(result[2], 20);
            assert_eq!(result[4], 30);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_ffpcn_sbyte() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [i8; 5] = [10, -128, 20, -128, 30];
            let nulval: i8 = -128;

            let mut f = make_table(&name, "VALUE", "1S", 5, &mut status);
            fits_set_btblnull(
                f.as_deref_mut().unwrap(),
                1,
                nulval as LONGLONG,
                &mut status,
            );
            fits_write_colnull(
                f.as_deref_mut().unwrap(),
                TSBYTE,
                1,
                1,
                1,
                5,
                cast_slice(&data),
                Some(NullValue::Byte(nulval)),
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0i8; 5];
            let mut anynull = -1;
            fits_read_col(
                f.as_deref_mut().unwrap(),
                TSBYTE,
                1,
                1,
                1,
                5,
                Some(NullValue::Byte(-127)),
                cast_slice_mut(&mut result),
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[0], 10);
            assert_eq!(result[2], 20);
            assert_eq!(result[4], 30);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_ffpcn_ushort() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [u16; 5] = [100, 65535, 200, 65535, 300];
            let nulval: u16 = 65535;

            let mut f = make_table(&name, "VALUE", "1U", 5, &mut status);
            fits_set_btblnull(
                f.as_deref_mut().unwrap(),
                1,
                nulval as LONGLONG,
                &mut status,
            );
            fits_write_colnull(
                f.as_deref_mut().unwrap(),
                TUSHORT,
                1,
                1,
                1,
                5,
                cast_slice(&data),
                Some(NullValue::UShort(nulval)),
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0u16; 5];
            let mut anynull = -1;
            fits_read_col(
                f.as_deref_mut().unwrap(),
                TUSHORT,
                1,
                1,
                1,
                5,
                Some(NullValue::UShort(65534)),
                cast_slice_mut(&mut result),
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[0], 100);
            assert_eq!(result[2], 200);
            assert_eq!(result[4], 300);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_ffpcn_uint() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [u32; 5] = [100, 4294967295, 200, 4294967295, 300];
            let nulval: u32 = 4294967295;

            let mut f = make_table(&name, "VALUE", "1V", 5, &mut status);
            fits_set_btblnull(
                f.as_deref_mut().unwrap(),
                1,
                nulval as LONGLONG,
                &mut status,
            );
            fits_write_colnull(
                f.as_deref_mut().unwrap(),
                TUINT,
                1,
                1,
                1,
                5,
                cast_slice(&data),
                Some(NullValue::UInt(nulval)),
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0u32; 5];
            let mut anynull = -1;
            fits_read_col(
                f.as_deref_mut().unwrap(),
                TUINT,
                1,
                1,
                1,
                5,
                Some(NullValue::UInt(4294967294)),
                cast_slice_mut(&mut result),
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[0], 100);
            assert_eq!(result[2], 200);
            assert_eq!(result[4], 300);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_ffpcn_ulong() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [libc::c_ulong; 5] = [100, 4294967295, 200, 4294967295, 300];
            let nulval: libc::c_ulong = 4294967295;

            let mut f = make_table(&name, "VALUE", "1V", 5, &mut status);
            fits_set_btblnull(
                f.as_deref_mut().unwrap(),
                1,
                nulval as LONGLONG,
                &mut status,
            );
            fits_write_colnull(
                f.as_deref_mut().unwrap(),
                TULONG,
                1,
                1,
                1,
                5,
                cast_slice(&data),
                Some(NullValue::ULong(nulval)),
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0 as libc::c_ulong; 5];
            let mut anynull = -1;
            fits_read_col(
                f.as_deref_mut().unwrap(),
                TULONG,
                1,
                1,
                1,
                5,
                Some(NullValue::ULong(4294967294)),
                cast_slice_mut(&mut result),
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[0], 100);
            assert_eq!(result[2], 200);
            assert_eq!(result[4], 300);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_ffpcn_long() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [c_long; 5] = [100, -2147483647, 200, -2147483647, 300];
            let nulval: c_long = -2147483647;

            let mut f = make_table(&name, "VALUE", "1J", 5, &mut status);
            fits_set_btblnull(
                f.as_deref_mut().unwrap(),
                1,
                nulval as LONGLONG,
                &mut status,
            );
            fits_write_colnull(
                f.as_deref_mut().unwrap(),
                TLONG,
                1,
                1,
                1,
                5,
                cast_slice(&data),
                Some(NullValue::Long(nulval)),
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0 as c_long; 5];
            let mut anynull = -1;
            fits_read_col(
                f.as_deref_mut().unwrap(),
                TLONG,
                1,
                1,
                1,
                5,
                Some(NullValue::Long(-2147483646)),
                cast_slice_mut(&mut result),
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[0], 100);
            assert_eq!(result[2], 200);
            assert_eq!(result[4], 300);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_ffpcn_longlong() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [LONGLONG; 5] = [100, -9223372036854775807, 200, -9223372036854775807, 300];
            let nulval: LONGLONG = -9223372036854775807;

            let mut f = make_table(&name, "VALUE", "1K", 5, &mut status);
            fits_set_btblnull(f.as_deref_mut().unwrap(), 1, nulval, &mut status);
            fits_write_colnull(
                f.as_deref_mut().unwrap(),
                TLONGLONG,
                1,
                1,
                1,
                5,
                cast_slice(&data),
                Some(NullValue::LONGLONG(nulval)),
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0 as LONGLONG; 5];
            let mut anynull = -1;
            fits_read_col(
                f.as_deref_mut().unwrap(),
                TLONGLONG,
                1,
                1,
                1,
                5,
                Some(NullValue::LONGLONG(-9223372036854775806)),
                cast_slice_mut(&mut result),
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[0], 100);
            assert_eq!(result[2], 200);
            assert_eq!(result[4], 300);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_ffpcn_logical() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            // The C test used { 'T', 0, 'F', 0, 'T' }, but for TLOGICAL any
            // nonzero byte is TRUE, so 'F' (70) round-trips to 1 (true), not 0.
            // The 0 entries equal the null value and are written as nulls.
            let data: [c_char; 5] = [1, 0, 1, 0, 1];
            let nulval: c_schar = 0;

            let mut f = make_table(&name, "FLAG", "1L", 5, &mut status);
            fits_write_colnull(
                f.as_deref_mut().unwrap(),
                TLOGICAL,
                1,
                1,
                1,
                5,
                cast_slice(&data),
                Some(NullValue::Byte(nulval)),
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0 as c_char; 5];
            let mut anynull = -1;
            fits_read_col(
                f.as_deref_mut().unwrap(),
                TLOGICAL,
                1,
                1,
                1,
                5,
                Some(NullValue::Byte(2)),
                cast_slice_mut(&mut result),
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[0], 1);
            assert_eq!(result[2], 1);
            assert_eq!(result[4], 1);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_ffpcn_string() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let strings = ["Hello", "NULL", "World", "NULL", "Test"];

            let mut f = make_table(&name, "NAME", "10A", 5, &mut status);

            // The byte buffer holds an array of `*const c_char` pointers to
            // per-element string buffers. The generic ffpcn TSTRING path reads
            // FLEN_VALUE bytes from each pointer, so back each string with a
            // full-width, NUL-padded buffer.
            let str_bufs: Vec<Vec<c_char>> = strings
                .iter()
                .map(|s| {
                    let mut v = vec![0 as c_char; 71];
                    for (i, b) in s.bytes().enumerate() {
                        v[i] = b as c_char;
                    }
                    v
                })
                .collect();
            let ptrs: Vec<*const c_char> = str_bufs.iter().map(|v| v.as_ptr()).collect();
            let ptr_bytes: &[u8] = unsafe {
                core::slice::from_raw_parts(
                    ptrs.as_ptr().cast::<u8>(),
                    core::mem::size_of_val(&ptrs[..]),
                )
            };
            let nulval = CString::new("NULL").unwrap();
            fits_write_colnull(
                f.as_deref_mut().unwrap(),
                TSTRING,
                1,
                1,
                1,
                5,
                ptr_bytes,
                Some(NullValue::String(nulval)),
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);

            let mut buf: [[c_char; 11]; 5] = [[0; 11]; 5];
            let mut read_ptrs: [*mut c_char; 5] = core::array::from_fn(|i| buf[i].as_mut_ptr());
            let read_ptr_bytes: &mut [u8] = unsafe {
                core::slice::from_raw_parts_mut(
                    read_ptrs.as_mut_ptr().cast::<u8>(),
                    core::mem::size_of_val(&read_ptrs),
                )
            };
            let readnul = CString::new("UNDEF").unwrap();
            let mut anynull = -1;
            fits_read_col(
                f.as_deref_mut().unwrap(),
                TSTRING,
                1,
                1,
                1,
                5,
                Some(NullValue::String(readnul)),
                read_ptr_bytes,
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            let read_str = |row: usize| -> Vec<u8> {
                buf[row]
                    .iter()
                    .take_while(|&&c| c != 0)
                    .map(|&c| c as u8)
                    .collect()
            };
            assert_eq!(read_str(0), b"Hello");
            assert_eq!(read_str(2), b"World");
            assert_eq!(read_str(4), b"Test");
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_ffpcn_complex() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [f32; 6] = [1.0, 2.0, -999.0, -999.0, 3.0, 4.0];
            let nulval: f32 = -999.0;

            let mut f = make_table(&name, "CPLX", "1C", 5, &mut status);
            fits_write_colnull(
                f.as_deref_mut().unwrap(),
                TCOMPLEX,
                1,
                1,
                1,
                3,
                cast_slice(&data),
                Some(NullValue::Float(nulval)),
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0f32; 6];
            let mut anynull = -1;
            fits_read_col(
                f.as_deref_mut().unwrap(),
                TCOMPLEX,
                1,
                1,
                1,
                3,
                Some(NullValue::Float(-888.0)),
                cast_slice_mut(&mut result),
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[0], 1.0);
            assert_eq!(result[1], 2.0);
            assert_eq!(result[4], 3.0);
            assert_eq!(result[5], 4.0);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_ffpcn_dblcomplex() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let data: [f64; 6] = [1.0, 2.0, -999.0, -999.0, 3.0, 4.0];
            let nulval: f64 = -999.0;

            let mut f = make_table(&name, "DCPLX", "1M", 5, &mut status);
            fits_write_colnull(
                f.as_deref_mut().unwrap(),
                TDBLCOMPLEX,
                1,
                1,
                1,
                3,
                cast_slice(&data),
                Some(NullValue::Double(nulval)),
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            fits_movabs_hdu(f.as_deref_mut().unwrap(), 2, None, &mut status);
            let mut result = [0f64; 6];
            let mut anynull = -1;
            fits_read_col(
                f.as_deref_mut().unwrap(),
                TDBLCOMPLEX,
                1,
                1,
                1,
                3,
                Some(NullValue::Double(-888.0)),
                cast_slice_mut(&mut result),
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[0], 1.0);
            assert_eq!(result[1], 2.0);
            assert_eq!(result[4], 3.0);
            assert_eq!(result[5], 4.0);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    // ----------------------------------------------------------------------
    // ffppxnll tests (write pixels with null value, LONGLONG *firstpix)
    // ----------------------------------------------------------------------

    #[test]
    fn test_ffppxnll_short() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 1] = [5];
            let data: [i16; 5] = [10, -999, 20, -999, 30];
            let nulval: i16 = -999;
            let firstpix: [LONGLONG; 1] = [1];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), SHORT_IMG, 1, &naxes, &mut status);
            fits_set_imgnull(f.as_deref_mut().unwrap(), nulval as LONGLONG, &mut status);
            fits_write_pixnullll(
                f.as_deref_mut().unwrap(),
                TSHORT,
                &firstpix,
                5,
                cast_slice(&data),
                Some(NullValue::Short(nulval)),
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut result = [0i16; 5];
            let mut anynull = -1;
            fits_read_img(
                f.as_deref_mut().unwrap(),
                TSHORT,
                1,
                5,
                None,
                cast_slice_mut(&mut result),
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[0], 10);
            assert_eq!(result[2], 20);
            assert_eq!(result[4], 30);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    #[test]
    fn test_ffppxnll_float() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);
            let naxes: [c_long; 1] = [5];
            let data: [f32; 5] = [1.5, -999.0, 2.5, -999.0, 3.5];
            let nulval: f32 = -999.0;
            let firstpix: [LONGLONG; 1] = [1];

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), FLOAT_IMG, 1, &naxes, &mut status);
            fits_write_pixnullll(
                f.as_deref_mut().unwrap(),
                TFLOAT,
                &firstpix,
                5,
                cast_slice(&data),
                Some(NullValue::Float(nulval)),
                &mut status,
            );
            fits_close_file(f.take().unwrap(), &mut status);

            fits_open_file(&mut f, &name, READONLY, &mut status);
            let mut result = [0f32; 5];
            let mut anynull = -1;
            fits_read_img(
                f.as_deref_mut().unwrap(),
                TFLOAT,
                1,
                5,
                Some(NullValue::Float(-888.0)),
                cast_slice_mut(&mut result),
                Some(&mut anynull),
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(result[0], 1.5);
            assert_eq!(result[2], 2.5);
            assert_eq!(result[4], 3.5);
            fits_close_file(f.take().unwrap(), &mut status);
        });
    }

    // ----------------------------------------------------------------------
    // ffpthp test (set heap pointer for binary table)
    // ----------------------------------------------------------------------

    #[test]
    fn test_ffpthp() {
        with_temp_file(|filename| {
            let mut status = 0;
            let name = to_buf(filename);

            let mut f: Option<Box<fitsfile>> = None;
            fits_create_file(&mut f, &name, &mut status);
            fits_write_imghdr(f.as_deref_mut().unwrap(), BYTE_IMG, 0, &[], &mut status);

            let ttype_v = [Some(cc("DATA"))];
            let ttype_ref: Vec<Option<&[c_char]>> = ttype_v.iter().map(|o| o.as_deref()).collect();
            // Variable length array.
            let tform_v = [cc("1PJ(10)")];
            let tform_ref: Vec<&[c_char]> = tform_v.iter().map(|v| v.as_slice()).collect();
            let extname = cc("DATA");
            fits_create_tbl(
                f.as_deref_mut().unwrap(),
                BINARY_TBL,
                10,
                1,
                &ttype_ref,
                &tform_ref,
                None,
                Some(&extname),
                &mut status,
            );

            // Set heap start at offset 1000 bytes from start of data.
            fits_write_theap(f.as_deref_mut().unwrap(), 1000, &mut status);

            // Verify THEAP keyword was written.
            let mut theap: c_long = 0;
            fits_read_key_lng(
                f.as_deref_mut().unwrap(),
                &cc("THEAP"),
                &mut theap,
                None,
                &mut status,
            );
            assert_eq!(status, 0);
            assert_eq!(theap, 1000);

            fits_close_file(f.take().unwrap(), &mut status);
        });
    }
}
